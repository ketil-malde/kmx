{-# Language BangPatterns #-}
module Kmers (kmers_noerr, kmers, kmers_rc, unkmer, val, unval, next1, next1_rc, initials) where

import Prelude hiding (null, scanl, sum)
import qualified Data.ByteString.Lazy.Char8 as B
import Data.ByteString.Lazy.Char8 (ByteString)
import Bio.Core.Sequence

-- import Data.ByteString.Lazy as BW
import Data.Word
import Data.Int
import Data.Bits
import Data.List (unfoldr)

-- | Calculate K-mers, ignoring Ns (counting them as As, so a toy solution)
-- this is insanely fast, five gigawords would take about 75 seconds on my laptop
-- but -fno-full-laziness brings it from 15 to 80 ms!
kmers_noerr :: Int64 -> ByteString -> [Word]
kmers_noerr k bs = if B.length bs < k then [] else first : go first rest
  where
    (w0,rest) = B.splitAt k bs
    first = B.foldl accum 0 w0
    go !w !s = case B.uncons s of 
      Just (!c,!s') -> let !w1 = accum w c in w1 : go w1 s'
      Nothing -> []
    accum w x = (w `shiftL` 2 .|. val x) .&. (4^k-1)
    
-- about two minutes for 5Gb
kmers :: Int64 -> ByteString -> [Word]
kmers k bs' = if B.length bs' < k then [] else go_incomplete k 0 bs'
  where
    go_incomplete 0 !w0 !bs = go_complete w0 bs
    go_incomplete i !w0 !bs = case B.uncons bs of
      Just (!c,!rs) -> if val c == 4 then go_incomplete k 0 rs
                       else go_incomplete (i-1) (accum w0 c) rs
      Nothing -> []
    go_complete !w0 !bs = w0 : case B.uncons bs of
      Just (!c,!rs) -> if val c /= 4 then go_complete (accum w0 c) rs
                       else go_incomplete k 0 rs
      Nothing -> []
    accum w x = (w `shiftL` 2 .|. val x) .&. (4^k-1)

-- both forward and reverse complement, equally fast(!)
kmers_rc :: Integral i => i -> SeqData -> [Word]
kmers_rc k' bs'' = if B.length bs' < k then [] else go_incomplete k 0 0 bs'
  where
    k = fromIntegral k'
    bs' = unSD bs''
    go_incomplete 0 !wf !wr !bs = go_complete wf wr bs
    go_incomplete i !wf !wr !bs = case B.uncons bs of
      Just (!c,!rs) -> if val c == 4 then go_incomplete k 0 0 rs
                       else go_incomplete (i-1) (accum wf c) (accum_rc wr c) rs
      Nothing -> []
    go_complete !wf !wr !bs = let !z = min wf wr in z : case B.uncons bs of
      Just (!c,!rs) -> if val c /= 4 then go_complete (accum wf c) (accum_rc wr c) rs
                       else go_incomplete k 0 0 rs
      Nothing -> []
    accum w x    = (w `shiftL` 2 .|. val x) .&. (4^k-1)
    accum_rc w x = (w `shiftR` 2 .|. (val_rc x `shiftL` (2*(fromIntegral k-1)))) .&. (4^k-1)

-- TODO: use these, but benchmark to make sure we don't lose performance
next1, next1_rc :: Integral i => i -> Word -> Char -> Word
next1 k w x = (w `shiftL` 2 .|. val x) .&. (4^k-1)
next1_rc k w x = (w `shiftR` 2 .|. (val_rc x `shiftL` (2*(fromIntegral k-1)))) .&. (4^k-1)

-- get first fwd and reverse hash from a sequence, needed for Reseq.paths 
initials :: Integral i => Int64 -> SeqData -> (i,Word,Word)
initials k' bs'' = if B.length bs' < k then error "too short input sequence" else go_incomplete 0 k 0 0 bs'
  where
    k = fromIntegral k'
    bs' = unSD bs''
    go_incomplete :: Integral i => i -> Int64 -> Word -> Word -> ByteString -> (i,Word,Word)
    go_incomplete !p 0 !wf !wr _ = (p,wf,wr)
    go_incomplete !p i !wf !wr !bs = case B.uncons bs of
      Just (!c,!rs) -> if val c == 4 then go_incomplete (p+1) k 0 0 rs
                       else go_incomplete (p+1) (i-1) (next1 k wf c) (next1_rc k wr c) rs
      Nothing -> error "ran out of sequence"

-- | Decode a k-mer to the corresponding sequence
unkmer :: Int64 -> Word -> String
unkmer k = reverse . take (fromIntegral k) . unfoldr dec1
  where dec1 w = Just (unval (w .&. 0x3), w `shiftR` 2)

unval :: Word -> Char
unval 0 = 'a'
unval 1 = 'g'
unval 2 = 't'
unval 3 = 'c'

val :: Char -> Word
val 'A' = 0    
val 'C' = 3
val 'G' = 1
val 'T' = 2
val 'a' = 0    
val 'c' = 3
val 'g' = 1
val 't' = 2
val _ = 4

val_rc :: Char -> Word
val_rc 'A' = 2    
val_rc 'C' = 1
val_rc 'G' = 3
val_rc 'T' = 0
val_rc 'a' = 2    
val_rc 'c' = 1
val_rc 'g' = 3
val_rc 't' = 0
val_rc _ = 4
