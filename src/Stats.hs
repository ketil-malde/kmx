module Stats where

import Text.Printf
import qualified Data.ByteString.Lazy.Char8 as B

type Histogram = [(Int,Double)]
data Distribution = Dist { lambda_err, lambda_dip, w_err, w_hap, w_dip, w_tetra :: Double } deriving Show  

-- invariant: sum of weights constant?

readHistogram :: FilePath -> IO Histogram
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

calcStats :: Histogram -> [Distribution]
calcStats hist = let
  start = Dist (lambda (average (take 4 hist)) 1) (lambda (average (drop 4 hist)) 2) 0.45 0.05 0.45 0.05
  step d = maximization d (expectation d hist)
  in iterate step start

average :: Histogram -> Double
average xs = total xs/cnt
  where cnt = sum . map snd $ xs

total :: Histogram -> Double
total xs = sum . map (\(x,y) -> fromIntegral x*y) $ xs

-- expectation: assign data to distributions
expectation :: Distribution -> Histogram -> [Histogram] -- histograms: error, hap, dip, tetra+
expectation (Dist le ld we wh wd wr) = go [] [] [] []
  where go errs haps dips reps ((x,v):xs) = let
          ws = normalize [we,wh,wd,wr] $ po_ratio2 [le,ld/2,ld,2*ld] (fromIntegral x)
          normalize ps qs = let ns = zipWith (*) ps qs in map (/sum ns) ns
          in go ((x,v*ws!!0):errs) ((x,v*ws!!1):haps) ((x,v*ws!!2):dips) ((x,v*ws!!3):reps) xs
        go errs haps dips reps [] = [errs,haps,dips,reps]

-- maximization: determine parameters from assigned data
maximization :: Distribution -> [Histogram] -> Distribution
maximization (Dist e' d' _ _ _ _) [he,hh,hd,hr] = Dist e d (tot he) (tot hh) (tot hd) (tot hr)
  where e = lambda (average he) e'
        dd = lambda (average hd) d'
        dh = lambda (average hh) (d'/2)
        d  = (dd*tot hd+2*dh*tot hh)/(tot hd+tot hh)
        tot = sum . map snd
        
-- calculate zero-truncated poisson distribution (naively), expensive for large x
ztpoisson :: Double -> Int -> Double
ztpoisson lam x = lam**fromIntegral x * exp(-lam)/( (1-exp(-lam)) * product [2..fromIntegral x])

po_ratio0, po_ratio1, po_ratio2 :: [Double] -> Int -> [Double]
po_ratio0 ls x = map (/ sum prs) prs
  where prs = [ztpoisson l x | l <- ls] -- ztpoisson overflows and is slooow

po_ratio1 ls x = map (/ sum prs) prs
  where prs = map (pr1 x) ls
        pr1 k l = exp (fromIntegral x*log l -l - log (1-exp(-l))) -- overflows for large x

-- stable calculation of prob ratios (max likelihood) for large x
-- tricksy: e^a / e^a + e^b + e^c + ... = 1/(1+e^(b-a)+e^(c-a)+...)
po_ratio2 ls x = [prs k | k <- exps]
  where exps = [fromIntegral x*log l -l - log (1-exp(-l)) | l <- ls]
        prs k = 1/sum [exp (e-k) | e <- exps]

-- estimate lambda from average of zero-truncated poisson data
-- thanks: Simpson (2014), Bioinformatics 30:9
lambda :: Double -> Double -> Double
lambda avg l0 = until' 0.0001 . iterate (\l -> avg*(1-exp(negate l))) $ l0
  where until' eps xs = if abs (xs!!0-xs!!1) < eps then xs!!1 else until' eps (drop 1 xs)

-- output

showDist :: Distribution -> Histogram -> IO ()
showDist d@(Dist le ld we wh wd wr) h = do
  let hs = expectation d h
      [te,th,td,tr] = map total hs
      fit = undefined  -- pointwise diff h and hs
      ect = we*le
      hct = wh*ld/2
      dct = wd*ld
      rct = wr*ld*2 -- dubious, but..
  print d
  printf "Diffs: %.0f %.0f %.0f %.0f\n" ect hct dct rct
  printf "       %.0f %.0f %.0f %.0f\n" te th td tr
  -- todo: really take into account read lenght and k-mer size
  printf "Genome size: %d, " ((round ((hct+dct+rct)/ld))::Int)
  printf "Error rate: %.4f, " (ect/(ect+hct+dct+rct))
  printf "Heterozygosity: %.4f, " (hct/(hct+dct))
  printf "Repeats: %.4f.\n" (rct/(hct+dct+rct))
  printf "Dist: lambda_e=%.5f, lambda_d=%.4f\n" le ld
  printf "Genome size: %d, " (round ((th+td+tr)/ld)::Int)
  printf "Error rate: %.4f, " (te/(te+th+td+tr))
  printf "Heterozygosity: %.4f, "(th/(th+td))
  printf "Repeats: %.4f.\n" (tr/(th+td+tr))
