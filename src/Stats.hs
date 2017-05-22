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
      noComment l = B.head l /= '#'
      check xs | null xs = error ("Couldn't read histogram '"++f++"'- empty input file?")
               | otherwise = xs
  check `fmap` map parse1 `fmap` filter noComment `fmap` B.lines `fmap` B.readFile f

estimate :: Bool -> Histogram -> Distribution
estimate dip = until' . take 100 . calcStats dip
  where until' (d1:d2:rest) = if similar d1 d2 then d2 else until' (d2:rest)
        until' [d] = d
        until' []  = error "estimate failed inexplicably, this is as surprising to me as it is to you."
        similar a b = abs (lambda_dip a - lambda_dip b) < 0.0001 && abs (lambda_err a - lambda_err b) < 0.0001

calcStats :: Bool -> Histogram -> [Distribution]
calcStats diploid hist = let
  start = Dist (lambda (average (take 1 hist)) 1) (lambda (average (drop 1 hist)) 2) 0.45 (if diploid then 0.05 else 0.00) 0.45 0.05 -- 0.00 for haploid
  step d = maximization diploid d (expectation d hist)
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
          [pe,ph,pd,pr] = normalize [we,wh,wd,wr] $ po_ratio2 [le,ld/2,ld,2*ld] (fromIntegral x) -- weighted probabilities
          normalize ps qs = let ns = zipWith (*) ps qs in map (/sum ns) ns
          in go ((x,v*pe):errs) ((x,v*ph):haps) ((x,v*pd):dips) ((x,v*pr):reps) xs
        go errs haps dips reps [] = [errs,haps,dips,reps]

-- maximization: determine parameters from assigned data
maximization :: Bool -> Distribution -> [Histogram] -> Distribution
maximization diploid (Dist e' d' _ _ _ _) [he,hh,hd,hr] = Dist e d (tot he) (tot hh) (tot hd) (tot hr)
  where e = lambda (average he) e'
        dd = lambda (average hd) d'
        dh = lambda (average hh) (d'/2)
        d  = if diploid then (dd*tot hd+2*dh*tot hh)/(tot hd+tot hh) -- avg of dip and hap
             else dd
        tot = sum . map snd
maximization _ _ _ = error "Internal error: maximization needs four histograms."

-- calculate zero-truncated poisson distribution (naively), expensive for large x
ztpoisson :: Double -> Int -> Double
ztpoisson lam x = lam**fromIntegral x * exp(-lam)/( (1-exp(-lam)) * product [2..fromIntegral x])

po_ratio0, po_ratio1, po_ratio2 :: [Double] -> Int -> [Double]
po_ratio0 ls x = map (/ sum prs) prs
  where prs = [ztpoisson l x | l <- ls] -- ztpoisson overflows and is slooow

po_ratio1 ls x = map (/ sum prs) prs
  where prs = map pr1 ls
        pr1 l = exp (fromIntegral x*log l -l - log (1-exp(-l))) -- overflows for large x

-- stable calculation of prob ratios (max likelihood) for large x
-- tricksy: e^a / e^a + e^b + e^c + ... = 1/(1+e^(b-a)+e^(c-a)+...)
po_ratio2 ls x = [prs k | k <- exps]
  where exps = [fromIntegral x*log l -l - log (1-exp(-l)) | l <- ls]
        prs k = 1/sum [exp (e-k) | e <- exps]

-- ----------------------------------------------------------------------------
-- add a dispersion factor alpha to the first (i.e. error) distribution
-- errors are inflated by x^alpha, so 0 gives no dispersion, positive gives
-- overdispersion (biased towards larger x), negative alpha gives underdispersion.

expect_disp :: Double -> Distribution -> Histogram -> [Histogram]
expect_disp alpha (Dist le ld we wh wd wr) = go [] [] [] []
  where go errs haps dips reps ((x,v):xs) = let
          [pe,ph,pd,pr] = normalize [we,wh,wd,wr] $ po_ratio_disp alpha [le,ld/2,ld,2*ld] (fromIntegral x)
          normalize ps qs = let ns = zipWith (*) ps qs in map (/sum ns) ns
          in go ((x,v*pe):errs) ((x,v*ph):haps) ((x,v*pd):dips) ((x,v*pr):reps) xs
        go errs haps dips reps [] = [errs,haps,dips,reps]

po_ratio_disp :: Double -> [Double] -> Int -> [Double]
po_ratio_disp alpha ls x = [prs k | k <- exps]
  where exps = case [fromIntegral x*log l -l - log (1-exp(-l)) | l <- ls] of
                (e0:es) -> (e0+alpha*log (fromIntegral x):es)
        prs k = 1/sum [exp (e-k) | e <- exps]

-- estimate lambda from average of zero-truncated poisson data
-- thanks: Simpson (2014), Bioinformatics 30:9
lambda :: Double -> Double -> Double
lambda avg l0 = until' 0.0001 . iterate (\l -> avg*(1-exp(negate l))) $ l0
  where until' eps xs = if abs (xs!!0-xs!!1) < eps then xs!!1 else until' eps (drop 1 xs)

-- output
showDist :: Maybe (Int,Int,Double) -> Bool -> Distribution -> Histogram -> [String]
showDist mkl diploid d@(Dist le ld _we _wh _wd _wr) h =
  let hs = case mkl of Just (_,_,a) -> expect_disp a d h
                       _            -> expectation d h
      [te,th,td,tr] = map total hs -- number of k-mers in data assigned to each category
  -- todo: check that size calcs make sense
  in (if diploid then printf "Dist: lambda_e=%.5f, lambda_d=%.4f, errs: %.0fM hap: %.0fM dip: %.0fM rep: %.0fM" le ld (te/1e6) (th/1e6) (td/1e6) (tr/1e6)
      else printf "Dist: lambda_e=%.5f, lambda_d=%.4f, errs: %.0fM hap: %.0fM rep: %.0fM" le ld (te/1e6) (td/1e6) (tr/1e6))
     : case mkl of Just (k,l,_) -> [concat [
                                       -- ld is k-mer coverage, don't adjust for base coverage
                                       printf "Genome size: %d, " (round ((th+td+tr)/ld)::Int) 
                                       -- Each (single nucleotide) error gives k-1 erroneous kmers
                                       , printf "Error rate: %.4f, " (te/(te+th+td+tr)/(fromIntegral k-1))
                                         -- Het: divide by k to get actual rate, assuming each SNP gives rise to k-1 k-mers in th
                                         -- Surely also divide by two?  (each variant gives two sets of k-1 k-mers)
                                       , if diploid then printf "Heterozygosity: %.4f, "(th/(th+td)/(fromIntegral k-1)) else ""
                                         -- repeats are considered continuous (significantly longer than k)
                                       , printf "Repeats: %.4f, " (tr/(th+td+tr))
                                         -- Waterman adjusted
                                       , printf "Base coverage: %.2f." (ld*fromIntegral (l-k+1)/fromIntegral l)
                                       ]]
                   Nothing -> []
