module Main where

import Options (Options(..), getArgs, genOutputBS, genOutput)
import FreqCount (mk_judy, add_count, assocs)
import Serialize (toByteString, readIndex, readIndices, getmagic, unpackPairs, packPairs, addmagic)
import Kmers (kmers_rc, unkmer)
import Filter
import Entropy
import Correlate (collect, collectSqrt, correlate, regression, corr0, regr0, merge, merge2With, mergePlus)

import Bio.Core.Sequence
import Bio.Sequence.FastQ
import Bio.Sequence.Fasta
import qualified Data.ByteString.Lazy as B
import Control.Monad (when, unless, replicateM)
import Text.Printf
import Data.List (intersperse)
import System.IO.Unsafe

import Data.Array.IO

main :: IO ()
main = do
  opts <- getArgs
  case opts of 
    Count {} -> count opts
    Hist {} -> hist opts
    Verify {} -> verify opts
    Correlate {} -> corr opts
    Dump {} -> dump opts
    Merge {} -> mergeindices opts
    Heatmap {} -> heatmap opts

-- | Build a k-mer count index
count :: Options -> IO ()
count opts = do
  let rfa f = map (unSD . seqdata) `fmap` readFasta f
      rfq f = map (unSD . seqdata) `fmap` readSangerQ f
  s <- concat `fmap` mapM (if fasta opts then rfa else rfq) (files opts)
  let kmer_seq = kmers_rc (fromIntegral $ kval opts)
  freqs <- mk_judy (2*fromIntegral (kval opts))
  let kms = case filter_bits opts of
             0 -> concatMap kmer_seq s
             _ -> Filter.runfilter opts (concatMap kmer_seq s)
  mapM_ (add_count freqs) kms
  genOutputBS opts =<< Serialize.toByteString freqs

-- filter input according to mincount and maxcount values
  
-- | Output a histogram of count frequencies, optionally grouped by k-mer complexity (entropy)
hist :: Options -> IO ()
hist opts = do
  let my_filter = (if mincount opts > 0 then filter ((>= fromIntegral (mincount opts)) . fst) else id)
                  . (if maxcount opts > 0 then filter ((<= fromIntegral (maxcount opts)) . fst) else id)
  (k,kvs) <- readIndex opts
  if complexity_classes opts > 0
    then do
      counters <- replicateM (complexity_classes opts) (mk_judy 16)
      let sel_counter (x:xs@(_:_)) e = if e <= 0.0 then x else sel_counter xs (e-step)
          sel_counter [x] _ = x
          step = max_entropy / fromIntegral (complexity_classes opts)
          ent = let k' = fromIntegral k in case complexity_mersize opts of
            1 -> return . entropy k'
            2 -> entropy2 k'
            n -> entropyN n k'
          max_entropy = let x = 4.0**fromIntegral  (complexity_mersize opts) in negate x * 1/x * log (1/x) / log 2 -- uh, -1 * log(1/x)/log 2? = log_2 1/x = -log_2 x ?
      let cnt (hash,occurs) = do
            e <- ent hash
            let c0 = sel_counter counters e
            add_count c0 (fromIntegral occurs)
      mapM_ cnt kvs
      my_assocs <- mapM assocs counters
      let format_line (ky,vs) = concat $ intersperse "\t" (show ky:map show vs)
          header = concat ("#count":[printf "\t<=%.2f" x | x <- [step,2*step..max_entropy]])
      genOutput opts $ unlines (header : [ format_line ln | ln <- my_filter (merge my_assocs)])
    else do  -- just a single histogram
      cts <- mk_judy 16
      mapM_ (add_count cts . fromIntegral . snd) kvs 
      as <- assocs cts
      genOutput opts $ unlines [show ky ++ "\t" ++ show v | (ky,v) <- my_filter as]

heatmap :: Options -> IO ()
heatmap opts =   case indices opts of
  [_,_] -> do 
    a <- newArray ((mincount opts,mincount opts),(maxcount opts,maxcount opts)) 0 :: IO (IOUArray (Int,Int) Int)
    [(k1,str1), (k2,str2)] <- readIndices opts
    when (k1 /= k2) (error "K-mer sizes don't match")
    let update :: (Int,Int) -> IO ()
        update (x,y) = do let ix = ( max (mincount opts) (min x (maxcount opts))
                                   , max (mincount opts) (min y (maxcount opts)))
                          v <- readArray a ix
                          writeArray a ix (v+1)
    mapM_ (update . snd) $ merge2With (,) str1 str2
    -- Generate output
    let format ((x,y),v) = show x++" "++show y++" "++show v
        intersperseAt n xs = case splitAt n xs of
          (hd,[]) -> hd
          (hd,tl) -> hd++"":intersperseAt n tl
    xs <- getAssocs a
    genOutput opts $ unlines $ intersperseAt (maxcount opts+1) $ map format xs
  _ -> error "heatmap requires exactly two input files."
               
-- | Check that a count index file is consistent
verify :: Options -> IO ()
verify opts = do
  (k,str) <- readIndex opts
  putStrLn ("KMX Index: k="++show k++" ("++show (k*2)++" bits)")
  let limit = 4^k-1
      check ((k1,_v1):(k2,v2):rest)
        | k1 < 0 || k1 > limit = error ("KMX Index: key out of range: "++show k1)
        | k2 <= k1             = error ("KMX Index: wrapped! "++show k1++" -> "++show k2)
        | otherwise            = check ((k2,v2):rest)
      check [(k1,_v1)]
        | k1 < 0 || k1 > limit = error ("KMX Index: key out of range: "++show k1)
        | otherwise            = putStrLn "KMX Index: OK"
      check [] = putStrLn "KMX Index: OK"
  check str

-- | Exctract k-mers and counts for specific k-mers and/or limited frequency counts
dump :: Options -> IO ()
dump opts = do
  (k',idx) <- readIndex opts
  let k = fromIntegral k'
      showPair (w,v) = ((if hashes opts then (show w ++ "\t" ++ show v) else (unkmer k w ++ "\t"++show v))
                        ++if complexity opts then "\t"++printf "%.3f\t%.3f\t%.3f" (entropy k w) (unsafePerformIO $ entropy2 k w) (unsafePerformIO $ entropyN 3 k w) else "")
      my_filter = (if mincount opts > 0 then filter ((>= fromIntegral (mincount opts)) . snd) else id) 
                  . (if maxcount opts > 0 then filter ((<= fromIntegral (maxcount opts)) . snd) else id)
      header = "# k="++show k
  genOutput opts $ unlines $ (header:) $ map showPair $ my_filter idx

-- | Compare k-mer counts between two or more indexes/files
corr :: Options -> IO ()
corr opts = case indices opts of 
  [f1,f2] -> do 
    (k1,s1) <- Serialize.getmagic `fmap` B.readFile f1
    (k2,s2) <- Serialize.getmagic `fmap` B.readFile f2
    when (k1 /= k2) $ error "Indices have different k-mer sizes"
    let z = (if sqrt_transform opts then collectSqrt else collect) (Serialize.unpackPairs s1) (Serialize.unpackPairs s2)
        c0 = corr0 z
        c  = correlate z
        (beta,alpha) = regression z
        beta0 = regr0 z
    when (sqrt_transform opts) (printf "Sqrt-transformed data:")
    printf "Correlation coefficient:     \tr=%.3f  R²=%.3f\n" c (c*c)
    printf "Linear regression parameters:\tbeta=%.2f alpha=%.2f\n" beta alpha
    printf "Zero intercept correlation:  \tr=%.3f  R²=%.3f\n" c0 (c0*c0)
    printf "Zero intercept regression:   \tbeta=%.2f\n" beta0
  _ -> error "Correlation requires exactly two index files"

mergeindices :: Options -> IO ()
mergeindices opts = do
  (ks1:kss) <- mapM (\f -> Serialize.getmagic `fmap` B.readFile f) (indices opts)
  unless (all (==fst ks1) (map fst kss)) $ error "Incorrect k-mer size"
  let ps = map (Serialize.unpackPairs . snd) (ks1:kss)
  genOutputBS opts $ addmagic (fst ks1) $ Serialize.packPairs $ mergePlus ps
