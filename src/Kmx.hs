module Main where

import Options (Options(..), getArgs, genOutputBS, genOutput)
import FreqCount (mk_judy, FreqCount(..))
import Serialize (toByteString, readIndex, readIndices, getmagic, unpackPairs, packPairs, addmagic)
import Kmers (kmers_rc, unkmer, initials)
import Filter
import Entropy
import Correlate (collect, collectSqrt, correlate, regression, corr0, regr0, merge2With, mergePlus, jaccard)
import Reseq (simpleeval,simple_illumina,paths,showpath)
import Stats
import Hist

import Bio.Core.Sequence
import Bio.Sequence.FastQ
import Bio.Sequence.Fasta
import qualified Data.ByteString.Lazy.Char8 as B
import Control.Monad (when, unless)
import Text.Printf
import Data.List (sort)
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
    Classify {} -> classify opts
    Jaccard {} -> jacc opts
    Reseq {} -> reseq opts
    Stats {} -> genstats opts

-- | Build a k-mer count index
count :: Options -> IO ()
count opts = do
  s <- concat `fmap` readSequenceData opts
  let kmer_seq = kmers_rc (kval opts)
  freqs <- mk_judy (2*fromIntegral (kval opts))
  let kms = case filter_bits opts of
             0 -> concatMap (kmer_seq . snd) s
             _ -> Filter.runfilter opts (concatMap (kmer_seq . snd) s)
  mapM_ (add_count freqs) kms
  genOutputBS opts =<< Serialize.toByteString freqs

-- Read fasta or fastq sequences specified on the command line, return sequence contents
readSequenceData :: Options -> IO [[(SeqLabel,SeqData)]]
readSequenceData opts = case files opts of
                             ["-"] -> stdinput
                             [] -> stdinput
                             fs -> mapM parseSeq fs
  where stdinput = error "Can't read sequence data from standard input - sorry!"
        parseSeq f = if fasta opts then map (\s -> (seqid s, seqdata s)) `fmap` readFasta f
                     else map (\s -> (seqid s, seqdata s)) `fmap` readSangerQ f

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
      header = "# k="++show k
  genOutput opts $ unlines $ (header:) $ map showPair idx

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

jacc :: Options -> IO ()
jacc opts = do
    [(k1,vs1), (k2,vs2)] <- readIndices opts
    when (k1 /= k2) (error "K-mer sizes don't match")
    -- printf "Jaccard probability: %.4f\n" (jaccard vs1 vs2)
    print (jaccard vs1 vs2)

mergeindices :: Options -> IO ()
mergeindices opts = do
  (ks1:kss) <- mapM (\f -> Serialize.getmagic `fmap` B.readFile f) (indices opts)
  unless (all (==fst ks1) (map fst kss)) $ error "Incorrect k-mer size"
  let ps = map (Serialize.unpackPairs . snd) (ks1:kss)
  genOutputBS opts $ addmagic (fst ks1) $ Serialize.packPairs $ mergePlus ps

-- read Fastq or Fasta sequences, classify them by k-mer spectra
classify :: Options -> IO ()
classify opts = do
  -- read index into a Judy array (mincount=2 might be a good idea?)
  (k,vs) <- readIndex opts
  idx <- mk_judy k
  mapM_ (uncurry (set_count idx)) vs

  -- output: median (quantiles), average, or full list
  let classSingle [x] = genOutputBS opts . B.unlines =<< mapM class1 x
      classSingle _ = error "can currently only classify a single input file"
      class1 (h,s) = do
        -- gen keys, lookup counts, output
        let kms = kmers_rc k s
        cs <- mapM (FreqCount.get_count idx) kms
        let avg, q25, q50, q75 :: Double
            avg = fromIntegral (sum cs)/fromIntegral (length cs)
            css = sort cs
            l   = length cs `div` 4
            head0 (x:_) = fromIntegral x
            head0 []    = 0
            q25 = head0 $ take l css
            q50 = head0 $ take (2*l) css
            q75 = head0 $ take (3*l) css
        return (B.pack $ printf "%s Avg: %.1f Q25: %.1f Q50: %.1f Q75: %.1f" (B.unpack $ unSL h) avg q25 q50 q75) -- quantiles?
  if not (paired opts) then
    classSingle =<< readSequenceData opts
    else case readSequenceData opts of
{-          [r1,r2] -> classPairs (zip r1 r2)
          [rs]    -> classPairs (unterleave rs)
            where unterleave (x1:x2:xs) = (x1,x2):unterleave xs
                  unterleave [] = []
                  unterleave _  = error "Odd number of sequences in interleaved paired file?" -}
          _ -> error "Pair classification needs either two sequence files\nor one interleaved file"

-- Resequence by guided path traversal (reconstruct sequences)
reseq :: Options -> IO ()
reseq opts = do
  (k,vs) <- readIndex opts -- TODO: support multiple indices
  putStrLn "Reading index..."
  idx <- mk_judy k
  mapM_ (uncurry (set_count idx)) vs
  putStrLn "...done"
  ss <- concat `fmap` readSequenceData opts
  mapM_ (process (maxlen opts) k idx) ss

process m k idx (h,s) = let
  ev = simpleeval simple_illumina s idx
  in case initials (fromIntegral k) s of
      Nothing -> return ()
      Just is  -> do
        putStr (">"++toString h++" ")
        putStrLn $ showpath (fromIntegral k) $ head $ paths ev is m

-- Output histograms - print stats on stdout: file name, params, goodness of fit
-- TODO: better format, run to convergence, calculate fit
-- With output option: generate histograms for the estimated distributions
genstats :: Options -> IO ()
genstats opts = do
  h <- readHistogram (histogram opts)
  -- mapM_ putStrLn $ map (\x -> showDist Nothing x h) (take 20 (calcStats h))
  let d = estimate (diploid opts) h
      mkl = if kval opts == 0 then Nothing else Just (kval opts,readlength opts)
  mapM_ putStrLn (showDist mkl (diploid opts) d h)
  when (not . null $ output opts) $ do
    let go :: Histogram -> [Histogram] -> [String]
        go ((k,v):rest) kvs = format ((k,v):map head kvs) : go rest (map tail kvs)
        go [] _ = []
        format xs = (show . fst . head $ xs) ++ concat [printf "\t%.2f" w | (_,w) <- xs ]
        ts = map sort (expect_disp (dispersion opts) d h) -- why are these always reversed?
    genOutput opts $ unlines $ go h ts
    return ()
