{-# Language BangPatterns #-}

module Serialize where

import Filter (combineKvs)
import FreqCount
import Options

import Data.Bits
import qualified Data.ByteString.Lazy as B
import qualified Data.ByteString.Lazy.Char8 as BC
import System.IO.Unsafe
import Data.Char (ord)

addmagic :: Int ->  B.ByteString -> B.ByteString
addmagic k = B.append (B.pack (map (fromIntegral . ord) "kmx:" ++ [fromIntegral k]))

getmagic :: B.ByteString -> (Int,B.ByteString)
getmagic b = let (ks,rest) = B.splitAt 5 b
             in if B.take 4 ks == BC.pack "kmx:" then (fromIntegral $ last $ B.unpack ks ,rest)
                else error ("Wrong magic number - expected \"kmx:\", but got "++show (BC.take 4 ks)
                            ++"\nThis does not look like a valid kmx index file.")

type Index = (Int,[(Word,Int)]) -- k and a list of key/value pairs

-- | readIndex reads an index from a single input file or standard input if no file is specified
--   If a k-value is specified in the options (and the index is built with a larger k) it
--   combines k-mer counts to calculate counts (approximately) using the shorter k-value.
readIndex :: Options -> IO Index
readIndex opts = do
  str <- case indices opts of
          [] -> B.getContents
          ["-"] -> B.getContents
          [f] -> B.readFile f
          _ -> error ("Multiple indices specified, but we only want one:\n\t"++show (indices opts))
  case opts of
   Verify {} -> let (k,s) = getmagic str in return (k `div` 2,unpackPairs s) -- verify doesn't have the k option
   _ -> return (parse1idx opts str)

readIndices :: Options -> IO [Index]
readIndices opts = if length (indices opts) <= 1 then do {c <- readIndex opts; return [c] }
                   else mapM (\f -> parse1idx opts `fmap` B.readFile f) (indices opts)

parse1idx :: Options -> B.ByteString -> Index
parse1idx opts inp = do
  let (bits,str) = getmagic inp
      k = bits `div` 2
      idx = unpackPairs str
  case kval opts of
   0 -> (k,idx)
   myk -> if myk <= k && myk > 0 then (myk, combineKvs (fromIntegral myk) idx)
          else error ("Illegal kmer value '"++show myk++"', must be less than "++show k++".")

toByteString :: FreqCount -> IO B.ByteString
toByteString f = (addmagic (key_bits f)) `fmap` packPairs `fmap` assocs f

mapM' :: (a -> IO b) -> [a] -> IO [b]
mapM' fn = sequence' . map fn
  where
    sequence' [] = return []
    sequence' (x:xs) = do
      r <- x
      rs <- unsafeInterleaveIO (sequence' xs)
      return (r:rs)

fromByteString :: (Int -> IO FreqCount) -> B.ByteString -> IO FreqCount
fromByteString mf bs = do
  let (k,rest) = getmagic bs
  f <- mf k
  mapM_ (uncurry (set_count f)) (unpackPairs rest)
  return f

-- | Pack a list of pairs to a bytestring, using delta coding and
--   byte-packing

packPairs :: [(Word,Int)] -> B.ByteString
packPairs = packList . delta

-- | Unpack a delta-coded, byte-packed list of pairs
unpackPairs :: B.ByteString -> [(Word,Int)]
unpackPairs = undelta . unpackList

-- | Delta coding, store only difference to next key
delta :: [(Word,Int)] -> [Word]
delta = go 0
  where go !_ [] = []
        go !i ((!k,!v):rest) = let x = k-i in x `seq` x : fromIntegral v : go k rest

-- | Inverse delta coding.
undelta :: [Word] -> [(Word,Int)]
undelta = go 0
  where go !i (w:v:rest) = let i' = i+w in i' `seq` (i',fromIntegral v) : go i' rest
        go !_ [] = []
        go !_ [_] = error "undelta given an odd length list"


packList :: (Integral i,Bits i) => [i] -> B.ByteString
packList = B.concat . map bytePack
{-# SPECIALIZE INLINE packList :: [Word] -> B.ByteString #-}

unpackList :: (Integral i,Bits i) => B.ByteString -> [i]
unpackList s | B.null s = []
             | otherwise = let (!x,!y) = byteUnpack s in x : unpackList y
{-# SPECIALIZE INLINE unpackList :: B.ByteString -> [Word] #-}

-- | bytePack produces a variable length string, each byte packs
-- 7 bits of payload, with the MSB set to 0 in all but the last byte.
-- Byte order is big endian. See also http://wiki.tcl.tk/4334
bytePack :: (Integral i,Bits i) => i -> B.ByteString
bytePack x = B.pack $ go [first] (x `shiftR` 7)
  where !first = 128 .|. fromIntegral (x .&. 127)
        go !acc !0 = acc
        go !acc !n = let !a = fromIntegral (n .&. 127) : acc
                     in go a (n `shiftR` 7)
{-# Specialize Inline bytePack :: Integer -> B.ByteString #-}
{-# Specialize Inline bytePack :: Int -> B.ByteString #-}
{-# Specialize Inline bytePack :: Word -> B.ByteString #-}

-- | unpack a bytePacked integer, warning, may overflow!
byteUnpack :: (Integral i, Bits i) => B.ByteString -> (i,B.ByteString)
byteUnpack b = go 0 b
  where go !i !s = case B.uncons s of
          Just (!x,!rest) -> let n = (i `shiftL` 7)+(fromIntegral x .&. 127)
                             in if x .&. 128 == 0 then go n rest
                                else n `seq` rest `seq` (n,rest)
          Nothing -> error ("Failed to decode bytestring"++show (B.take 30 b))
{-# Specialize byteUnpack ::  B.ByteString -> (Integer, B.ByteString) #-}
{-# Specialize byteUnpack ::  B.ByteString -> (Int, B.ByteString) #-}
{-# Specialize byteUnpack ::  B.ByteString -> (Word, B.ByteString) #-}

-- | Dump a FreqCount as text
writeHist :: FreqCount -> FilePath -> IO ()
writeHist c f = do
  as <- assocs c
  writeFile f $ unlines $ [ showKey k v | (k,v) <- as]

showKey :: Word -> Int -> String
showKey k v = show k++" "++show v
