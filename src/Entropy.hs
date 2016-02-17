{-# Language BangPatterns #-}

module Entropy where

import Data.Int
import Data.Word
import Data.Bits
import FreqCount (mk_vector, add_count, counts)

-- | Calculate (1-)entropy for k-mer
entropy :: Int64 -> Word -> Double
entropy k w = let
    (as,cs,gs,ts) = go k 0 0 0 0 w
    go :: Int64 -> Int -> Int -> Int -> Int -> Word -> (Double,Double,Double,Double)
    go !0 !a !c !g !t _ = (fromIntegral a,fromIntegral c,fromIntegral g,fromIntegral t)
    go !x !a !c !g !t cur = 
      let w' = (cur `shiftR` 2)
          x' = x-1
      in x' `seq` w' `seq` case cur .&. 0x3 of 
            0 -> let r = a+1 in r `seq` go x' r c g t w'
            1 -> let r = g+1 in r `seq` go x' a c r t w'
            2 -> let r = t+1 in r `seq` go x' a c g r w'
            3 -> let r = c+1 in r `seq` go x' a r g t w'
    n = as+cs+gs+ts
  in negate $ sum [ f * log f/log 2 | f <- filter (/=0) [as/n,cs/n,gs/n,ts/n]]
-- todo: faster if storing all counters in the same word?
-- lookup-table for 16bit qtys? (four and four nucs?)

-- | Calculate dinuc entropy for k-mer
--   Currently only marginally faster than entropyN 2
entropy2 :: Int64 -> Word -> IO Double
entropy2 k w = do
    cs <- mk_vector 4
    let go :: Int64 -> Word -> IO ()
        go !1 !_ = return ()
        go !x !cur = do
          add_count cs (cur .&. 0xF)
          go (x-1) (cur `shiftR` 2)
    go k w
    vs <- map fromIntegral `fmap` counts cs
    let vtot = sum vs
    return $ negate $ sum [ (v/vtot) * log (v/vtot)/log 2 | v <- filter (/=0) vs]

-- | Calculate dinuc entropy for k-mer
entropyN :: Int -> Int64 -> Word -> IO Double
entropyN n k w = do
    cs <- mk_vector (2*n)
    let go :: Int64 -> Word -> IO ()
        go 0 _ = return ()
        go !x !cur = do
          add_count cs (cur .&. (4^n-1))
          go (x-1) (cur `shiftR` 2)
    go (k-fromIntegral n+1) w
    vs <- map fromIntegral `fmap` counts cs
    let vtot = sum vs
    return $ negate $ sum [ (v/vtot) * log (v/vtot)/log 2 | v <- filter (/=0) vs]
