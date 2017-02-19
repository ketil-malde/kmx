module Hist where

import Options
import FreqCount
import Stats
import Entropy
import Correlate (merge)
import Data.List (intersperse)
import Control.Monad (replicateM)
import Serialize (readIndex)
import Text.Printf

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
          sel_counter [] _ = error "empty input to sel_counter"
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
      if stats opts -- WTF: commenting out this (up to ..else) reduces memory footprint, even if not (stats opts)!
        then do ss <- calcstats opts k cts
                genOutput opts $ unlines $ (map ("# "++) ss) ++
                  [show ky ++ "\t" ++ show v | (ky,v) <- my_filter as] -- my_filter should be part of readIndex
        else genOutput opts $ unlines [show ky ++ "\t" ++ show v | (ky,v) <- my_filter as]

-- this causes a huge residency, why?  Strictify d?
{-# NoInline calcstats #-}
calcstats :: Options -> Int -> FreqCount -> IO [String]
calcstats opts k cts = do
      as' <- assocs cts
      let as = [(fromIntegral x,fromIntegral y) | (x,y) <- as']
          hdr = "k="++show k++(case indices opts of
                                [i] ->" input="++i
                                []  ->" input=<stdin>"
                                _   ->"")
                ++if diploid opts then " (diploid statistics)" else " (haploid statistics)"
          d = estimate (diploid opts) as
          stats_out = showDist (Just k) (diploid opts) d as
      return (hdr:stats_out)
