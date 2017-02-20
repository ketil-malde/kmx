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
  if complexity_classes opts > 0
    then hist_with_complexity opts
    else if stats opts
         then hist_with_stats opts
         else hist_simple opts

mkcount :: Options -> IO (Int,FreqCount)
mkcount opts = do
      (k,kvs) <- readIndex opts
      cts <- mk_judy 16
      mapM_ (add_count cts . fromIntegral . snd) kvs
      return (k,cts)

hist_simple :: Options -> IO ()
hist_simple opts = do
      (_k,cts) <- mkcount opts
      as <- assocs cts
      genOutput opts $ unlines [show ky ++ "\t" ++ show v | (ky,v) <- as]

hist_with_stats :: Options -> IO ()
hist_with_stats opts = do
      (k,cts) <- mkcount opts
      as <- assocs cts
      ss <- calcstats opts k cts
      genOutput opts $ unlines $ (map ("# "++) ss) ++ [show ky ++ "\t" ++ show v | (ky,v) <- as]

hist_with_complexity :: Options -> IO ()
hist_with_complexity opts = do
      (k,kvs) <- readIndex opts
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
      genOutput opts $ unlines (header : [ format_line ln | ln <- merge my_assocs ])

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
