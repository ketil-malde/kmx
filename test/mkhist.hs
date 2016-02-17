module Main where

import FreqCount
import Serialize
import Kmers
import Bio.Core.Sequence
import Bio.Sequence.FastQ
import System.Environment (getArgs)

main :: IO ()
main = do
  [k',f] <- getArgs
  s <- readSangerQ f
  let k = read k'
      kmer_seq = kmers_rc k . unSD . seqdata
  freqs <- mk_judy (2*fromIntegral k)
  mapM_ (add_count freqs) $ concatMap kmer_seq s  -- concat $ parMap rpar <- no speedup
  cts <- mk_vector 16
  mapM_ (add_count cts . fromIntegral) =<< counts freqs
  writeHist cts (f++".hist")

writeHist :: FreqCount -> FilePath -> IO ()
writeHist c f = do
  es <- keys c
  ks <- counts c
  writeFile f $ unlines $ [ showKey k v | (k,v) <- zip es ks]