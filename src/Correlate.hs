-- | Correlate k-mer occurrences between data sets
--   Probably room for efficiency improvements!

module Correlate
       ( collect, collectSqrt, correlate, regression, corr0, regr0
       , merge, merge2With, mergeWith, mergePlus
       , jaccard ) where

import Data.List (foldl')
import Data.Word  -- needed for old GHC

-- merge key-value pairs, make a row for each minimum value
merge :: (Ord k, Num v) => [[(k,v)]] -> [(k,[v])]
merge kvs = if all null kvs then []
            else let kmin = minimum $ map (fst . head) $ filter (not . null) kvs
                     (this,rest) = go kmin [] [] kvs
                     go i ac1 ac2 (((k,v):ks):kks)
                       | k == i    = go i (v:ac1) (ks:ac2) kks
                       | otherwise = go i (0:ac1) (((k,v):ks):ac2) kks
                     go i ac1 ac2 ([]:kks) = go i (0:ac1) ([]:ac2) kks
                     go _ ac1 ac2 [] = (reverse ac1, reverse ac2)
                 in (kmin,this) : merge rest

-- Deprecated version of mergeWith, it is slower than the one below
-- mergeWith' :: (Ord k, Num a) => (a -> a -> a) -> [[(k,a)]] -> [(k,a)]
-- mergeWith' f = foldr1 (merge2With f)

-- merging in a triangle pattern, more limited than 'merge'
-- mergeWith :: (Ord k, Num a) => (a -> a -> a) -> [[(k,a)]] -> [(k,a)]
mergeWith :: (Int -> Int -> Int) -> [[(Word,Int)]] -> [(Word,Int)]
mergeWith _ [] = []
mergeWith _ [as] = as
mergeWith f (as:bs:rest) = merge2With f (merge2With f as bs) (mergeWith f rest)

-- Specialized version of the above
mergePlus :: [[(Word,Int)]] -> [(Word,Int)]
mergePlus [] = []
mergePlus [as] = as
mergePlus (as:bs:rest) = merge2With (+) (merge2With (+) as bs) (mergePlus rest)

{-# SPECIALIZE INLINE merge2With :: (Int -> Int -> Int) -> [(Word,Int)] -> [(Word,Int)] -> [(Word,Int)] #-}
{-# SPECIALIZE INLINE merge2With :: (Int -> Int -> (Int,Int)) -> [(Word,Int)] -> [(Word,Int)] -> [(Word,(Int,Int))] #-}
merge2With :: (Ord k, Num a, Num b) => (a -> b -> c) -> [(k,a)] -> [(k,b)] -> [(k,c)]
merge2With f ((k1,a):rest1) ((k2,b):rest2) 
  | k1 == k2 = (k1,f a b) : merge2With f rest1 rest2
  | k1 < k2  = (k1,f a 0) : merge2With f rest1 ((k2,b):rest2) 
  | k1 > k2  = (k2,f 0 b) : merge2With f ((k1,a):rest1) rest2
merge2With f [] rest = map (\(k,v) -> (k,f 0 v)) rest
merge2With f rest [] = map (\(k,v) -> (k,f v 0)) rest
                                         
-- correlate? assuming through origo
regr0 :: C -> Double
regr0 r = xys r / x2s r

corr0 :: C -> Double  -- pearson's r
corr0 r = xys r / sqrt (x2s r * y2s r)

regression :: C -> (Double,Double) -- alpha and beta
regression r = let
  ns = fromIntegral (n r)
  beta = (xys r - xs r*ys r/ns) / (x2s r- xs r*xs r/ns)
  alpha = (ys r-beta*xs r)/ns
  in (beta,alpha)

-- | Calculate Pearson's r (correlation coefficient)
-- Example from WP: correlate $ collect (zip [1..] [0,10,101,102]) (zip [1..] [1,100,500,2000]) = 0.7544...
correlate :: C -> Double
correlate res = (fromIntegral (n res) * xys res - xs res * ys res)
                  / sqrt ((fromIntegral (n res) * x2s res - xs res * xs res)
                     * (fromIntegral (n res) * y2s res - ys res * ys res))
        
-- | Collect data in order to calculate various statistics
--   Note that as the values (k-mer counts) are binomial sampling,         
--   stdev scales with the sqrt of the count.  We should probably normalize, 
--   in order to get a uniform stdev.
collect :: Ord k => [(k, Int)] -> [(k, Int)] -> C
collect as bs = foldl' acc1 zero $ map snd $ merge2With (,) as bs
  where zero = C 0 0 0 0 0 0

data C = C { n :: !Integer, xs, ys, x2s, xys, y2s :: !Double }

acc1 :: C -> (Int,Int) -> C
acc1 c (x',y') = C { n = n c + 1
               , xs = xs c + x
               , ys = ys c + y
               , x2s = x2s c + x*x
               , xys = xys c + x*y
               , y2s = y2s c + y*y
               }
  where x = fromIntegral x'
        y = fromIntegral y'

collectSqrt:: Ord k => [(k, Int)] -> [(k, Int)] -> C
collectSqrt as bs = foldl' accSqrt1 zero $ map snd $ merge2With (,) as bs
  where zero = C 0 0 0 0 0 0

accSqrt1 :: C -> (Int,Int) -> C
accSqrt1 c (x',y') = C { n = n c + 1
               , xs = xs c + sqrt x
               , ys = ys c + sqrt y
               , x2s = x2s c + x
               , xys = xys c + sqrt(x*y)
               , y2s = y2s c + y
               }
  where x = fromIntegral x'
        y = fromIntegral y'

-- Calculate a generalized jaccard distance, i.e. the probability of a random kmer picked
-- from a random data set being occurring in both data sets.
jaccard :: Ord k => [(k, Int)] -> [(k, Int)] -> (Int,Int,Int,Int,Double)
jaccard as bs = jac $ foldl' jcount z $ map snd $ merge2With (,) as bs
  where z = (0,0,0,0)
        jcount (a,ab,ba,b) (0,x) = let b' = b+x in b' `seq` (a,ab,ba,b')
        jcount (a,ab,ba,b) (x,0) = let a' = a+x in a' `seq` (a',ab,ba,b)
        jcount (a,ab,ba,b) (x,y) = let ab' = ab+x; ba' = ba+y
                                   in ab' `seq` ba' `seq` (a,ab',ba',b)
        jac (a,ab,ba,b) = (a,ab,ba,b
                          ,(fromIntegral ab/fromIntegral (a+ab) + fromIntegral ba/fromIntegral (ba+b))/2)

