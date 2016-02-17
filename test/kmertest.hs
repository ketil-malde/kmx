{-# Language BangPatterns #-}

module Main where

import Prelude hiding (sum)

import FreqCount
import Kmers
import Serialize as S

import Data.Word
import Data.Int
import qualified Data.List as L
import Criterion.Main
import qualified Data.ByteString.Lazy.Char8 as B
import Bio.Core.Sequence
import Bio.Sequence.FastQ
import System.Environment (getArgs)
import Data.Time.Clock

-- notes: kmers are reasonably fast, 1.1s to sum kmers from 250K reads (100bp), i.e. 22Mb/sec - still CPU bound though
-- judy is slow...ish, 13 seconds for counting kmers in a judy array - only 2Mb/sec (but better with insertWith)
-- that means 2000 seconds for five gigs, maybe not too shabby? (jellyfish: 270s, but 8 cores)

-- Criterion tests will only run if no cmdline parameter is given?
main :: IO ()
main = do
  -- main_vector5  -- 5s, vector10 = 50 secs! 800MB
  x <- getArgs
  if null x
     then main_criterion 
     else do  
       -- the rest requires a cmdline arg
       print "Summing kmers from file"
       time main_sum       -- 1.5 secs, 3.5MB
       print "Storing kmers from file in judy array"  
       j <- time main_judy     -- 11 secs, 250MB
       print . sum =<< keys j  
       print "Writing to file"
       time $ B.writeFile "tmp.freq" =<< S.toByteString j
       print "Reading from file"
       j2 <- time (S.fromByteString mk_judy =<< B.readFile "tmp.freq")
       print . sum =<< keys j2
    
time :: IO a -> IO a
time a = do
  t1 <- getCurrentTime
  r <- a
  t2 <- getCurrentTime
  print (diffUTCTime t2 t1)
  return r

main_criterion :: IO ()
main_criterion = defaultMain 
       [ bench "kmers, 10K, k=10" $ nf (sum . kmers 10) (B.replicate 10000 'A')
       , bench "kmers, 10K, k=20" $ nf (sum . kmers 20) (B.replicate 10000 'A')
       , bench "kmers_noerr, 1M, k=20" $ nf (sum . kmers_noerr 20) (B.replicate 1000000 'A') -- 16ms, 60MB/s
       , bench "kmers,       1M, k=20" $ nf (sum . kmers 20) (B.replicate 1000000 'A')    -- 22ms, 44MB/s
       , bench "kmers_rc,    1M, k=20" $ nf (sum . kmers_rc 20) (B.replicate 1000000 'A') -- 26ms, 38MB/s
       --
       , bench "bytePack/w" $ nf (sum . unpackList . packList) ([1000000*x| x <- [1..10000::Word]]) -- 1.7ms
       , bench "bytePack/i" $ nf (sum . unpackList . packList) ([1000000*x| x <- [1..10000::Int]])  -- 1.7ms
       , bench "bytePack/n" $ nf (sum . unpackList . packList) ([1000000*x| x <- [1..10000::Integer]]) -- 7 ms
       --
       , bench "packPairs"  $ nf (sum . map fst . unpackPairs . packPairs) (zip [1000000*x | x <- [1..10000::Word]] [1..10000::Int]) -- 3.3ms
       , bench "entropy"    $ nf (L.foldl' (+) 0 . map (entropy 20)) [1..100000::Word]  -- 560ms - 436ms /10: 44ms
       ]

sum :: [Word] -> Word
sum = L.foldl' (+) 0

main_sum  :: IO ()
main_sum = print . sum . concatMap (kmers_rc 22 . unSD . seqdata) =<< readSangerQ . head =<< getArgs
  
main_judy :: IO FreqCount
main_judy = do
  [f] <- getArgs
  s <- readSangerQ f
  freqs <- mk_judy 44
  mapM_ (add_count freqs) $ concatMap (kmer_seq 22) s
  return freqs
  
main_vector5 :: IO ()
main_vector5 = do
  [f] <- getArgs
  s <- readSangerQ f
  freqs <- mk_vector 10
  mapM_ (add_count freqs) $ concatMap (kmer_seq 5) s  
  return ()

kmer_seq :: BioSeq s => Int64 -> s -> [Word]
kmer_seq k s = kmers_rc k $ unSD $ seqdata s

  