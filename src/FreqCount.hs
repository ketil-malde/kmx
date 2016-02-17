module FreqCount where

import qualified Data.Judy as J
import Data.Vector.Unboxed.Mutable as V
import Data.Maybe (fromMaybe)
-- import Data.Bits
-- import Data.HashTable.IO as H
import Data.Word
import System.IO.Unsafe

readv :: Unbox a => IOVector a -> Word -> IO a
readv v = V.unsafeRead v . fromIntegral -- use this for speed
-- readv v = V.read v . fromIntegral        -- use this for safety

data FreqCount = FreqCount
     { key_bits :: Int
     , add_count :: Word -> IO ()           -- ^ increase count for key by one
     , get_count      :: Word -> IO Int     -- ^ get the count for a single key
     , set_count :: Word -> Int -> IO ()    -- ^ set the count for a single key
       -- the following should be lazy (i.e. use interleaveIO)
     , keys   :: IO [Word]
     , counts :: IO [Int]
     , assocs :: IO [(Word,Int)]
     }

mk_vector :: Int -> IO FreqCount
mk_vector l = do
  -- putStrLn ("Init vector "++show l)
  v <- V.replicate (2^l) (0::Int)
  let ac k = do { x <- readv v k; V.write v (fromIntegral k) $! x+1 }
      gc k = readv v k
      sc k x = V.write v (fromIntegral k) x -- or unsafe?
      ks = go 0 where go i = if i >= 2^l then return []
                             else do 
                               x <- readv v i
                               rest <- unsafeInterleaveIO (go (i+1))
                               if x /= 0 then return (i:rest)
                                 else return rest
      es = do xs <- mapM (V.unsafeRead v . fromIntegral) =<< ks  -- NOT LAZY!
              return (filter (/=0) xs)
      -- like ks, but return key and value
      as = go 0 where go i = if i >= 2^l then return []
                             else do 
                               x <- readv v i
                               rest <- unsafeInterleaveIO (go (i+1))
                               if x /= 0 then return ((i,x):rest)
                                 else return rest
  return $ FreqCount { key_bits = l
                    , add_count = ac, get_count = gc, set_count = sc
                    , keys = ks, counts = es, assocs = as }

mk_judy :: Int -> IO FreqCount
mk_judy l = do
  -- putStrLn ("Init judy "++show l)
  j <- J.new :: IO (J.JudyL Int)
  let ac k = k `seq` J.insertWith (+) k 1 j
      gc k = k `seq` fromMaybe 0 `fmap` J.lookup k j
      sc k v = J.insert k v j
      es = J.elems j
      ks = J.keys j
      as = J.toList j
  return $ FreqCount { key_bits = l
                     , add_count = ac, get_count = gc, set_count = sc
                     , keys = ks, counts = es, assocs = as
                     }

{-
type HT = BasicHashTable Word Int -- smaller count type?

mk_hash :: IO FreqCount
mk_hash = do
  putStrLn ("Init Hash")
  h <- H.new :: IO HT
  let ac k = do { x <- H.lookup h k ; H.insert h k (1+fromMaybe 0 x) }
      gc k = k `seq` fromMaybe 0 `fmap` H.lookup h k 
      sc k v = H.insert h k v
      es = map snd `fmap` H.toList h
      ks = map fst `fmap` H.toList h
  return $ FreqCount { prefix_bits = 0, main_bits = 0
                     , add_count = ac, get_count = gc
                     , set_count = sc
                     , keys = ks
                     , counts = es
                     }
-}

{-
-- Combined data structure, a vector prefix, then Judy.
mk_vector_judy :: Int -> Int -> IO FreqCount
mk_vector_judy k l = do
  putStrLn ("Init Vector "++show k++" Judy "++show l)
  vj <- V.replicateM (2^k) (J.new :: IO (J.JudyL Int))
  let split i = (i .&. (2^k-1),i `shiftR` k)
      ac i = do 
        let (x,y) = split i
        j <- readv vj x
        J.insertWith (+) (fromIntegral y) 1 j
      gc i = do
        let (x,y) = split i
        j <- readv vj x
        fromMaybe 0 `fmap` J.lookup (fromIntegral y) j
      sc i t = do
        let (x,y) = split i
        j <- readv vj x
        J.insert (fromIntegral y) t j 
      es = undefined
      ks = undefined
  return $ FreqCount { prefix_bits = 0, main_bits = l
                     , add_count = ac, get_count = gc, set_count = sc
                     , keys = ks, counts = es
                     }
-}
