module Filter (runfilter, combineKvs) where

import Data.Bits
import Data.Word
import Options (Options(..))

runfilter :: Options -> [Word] -> [Word]
runfilter opts = let
  b = filter_bits opts :: Int
  k = kval opts
  fv = filter_value opts
  off = 6 -- shift the filter a bit to avoid too much skew in the partition sizes
  s = 2*k-off-b
  m = (2^b-1) `shiftL` s -- 000011..10000000000
  check x = (x .&. m) `shiftR` s == fv
  in if k < b+off then error "Can't use filter on very short kmers."
     else if fv >= 2^b then error "Illegal filter value, must be less than 2^bits"
     else filter check


-- | Project values down to a smaller k-value
combineKvs :: Int -> [(Word,Int)] -> [(Word,Int)]
combineKvs _ [] = []
combineKvs k ((w',c'):rest) = go (mask w') c' rest
  where go :: Word -> Int -> [(Word,Int)] -> [(Word,Int)]
        go m c [] = [(unmask m,c)] 
        go m c rs@(x:xs) = if mask (fst x) /= m then (unmask m,c) : combineKvs k rs
                           else go m (c+snd x) xs
        mask :: Word -> Word
        mask = (.&.) ((2^(k*2)-1) `shiftL` (64-2*k)) -- correct?
        unmask y = y `shiftR` (64-2*k)
