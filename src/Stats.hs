module Stats where

import qualified Data.ByteString.Lazy.Char8 as B

type Histogram a = [(Int,a)]
data Distribution = Dist { lambda_err, lambda_dip, w_err, w_hap, w_dip, w_tetra :: Double } deriving Show  

-- invariant: sum of weights constant?

readHistogram :: FilePath -> IO (Histogram Double)
readHistogram f = do
  let parse1 :: B.ByteString -> (Int,Double)
      parse1 l = case B.words l of
        [w1,w2] -> case B.readInt w1 of
          Just (x,_) -> case B.readInt w2 of
            Just (y,_) -> (x,fromIntegral y)
            _ -> my_error ("Couldn't parse :"++B.unpack w2)
          _ -> my_error ("Couldn't parse :"++B.unpack w1)            
        _ -> my_error ""
        where my_error x = error ("Couldn't parse file '"++f++"' as a histogram.\nOffending line:\n"++B.unpack (B.take 72 l)++"\n"++x)
  map parse1 `fmap` B.lines `fmap` B.readFile f

calcStats :: Histogram Double -> [Distribution]
calcStats hist = let
  start = Dist (lambda (average (take 4 hist)) 1) (lambda (average (drop 4 hist)) 2) 0.45 0.05 0.45 0.05
  step d = maximization (expectation d hist)
  in iterate step start

average :: Histogram Double -> Double
average xs = total xs / cnt xs
  where cnt = sum . map snd
        total = sum . map (\(x,y) -> fromIntegral x*y)

-- expectation: assign data to distributions
expectation :: Distribution -> Histogram Double -> [Histogram Double] -- histograms: error, hap, dip, tetra+
expectation (Dist le ld we wh wd wr) = go [] [] [] []
  where go errs haps dips reps ((x,v):xs) = let
          ws = po_ratio2 [le,ld/2,ld,2*ld] (fromIntegral x)
          in go ((x,we*ws!!0):errs) ((x,wh*ws!!1):haps) ((x,wd*ws!!2):dips) ((x,wr*ws!!3):reps) xs
        go errs haps dips reps [] = [errs,haps,dips,reps]

-- maximization: determine parameters from assigned data
maximization :: [Histogram Double] -> Distribution
maximization [he,hh,hd,hr] = Dist e d (tot he) (tot hh) (tot hd) (tot hr)
  where e = lambda (average he) 1
        d1 = lambda (average hd) 3
        d2 = lambda (average hh) 2
        d  = (d1*tot hd+2*d2*tot hh)/(tot hd+tot hh)
        tot = sum . map snd
        
-- calculate zero-truncated poisson distribution (naively)
ztpoisson :: Double -> Double -> Double
ztpoisson lam x = lam**x * exp(-lam)/( (1-exp(-lam)) * product [2..x])

po_ratio0, po_ratio1, po_ratio2 :: [Double] -> Double -> [Double]
po_ratio0 ls x = map (/ sum prs) prs
  where prs = [ztpoisson l x | l <- ls] -- ztpoisson overflows and is slooow

po_ratio1 ls x = map (/ sum prs) prs
  where prs = map (pr1 x) ls
        pr1 k l = exp (x*log l -l - log (1-exp(-l))) -- overflows for large x

-- stable calculation of prob ratios (max likelihood) for large x
po_ratio2 ls x = [prs k | k <- exps]
  where exps = [x*log l -l - log (1-exp(-l)) | l <- ls]
        prs k = 1/sum [exp (e-k) | e <- exps]

-- estimate lambda from average of zero-truncated poisson data
-- thanks: Simpson (2014), Bioinformatics 30:9
lambda :: Double -> Double -> Double
lambda avg l0 = until' 0.0001 . iterate (\l -> avg*(1-exp(negate l))) $ l0
  where until' eps xs = if abs (xs!!0-xs!!1) < eps then xs!!1 else until' eps (drop 1 xs)
