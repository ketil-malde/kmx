-- time various combinations of freq counters, needs 
-- to be updated to match current status

module Main where 

import FreqCount

import System.Random
import System.Environment (getArgs)

main :: IO ()
main = do
  [a',b',c'] <- getArgs
  let abits = read a'
      bbits = read b'
      amax = 2^abits-1
      bmax = 2^bbits-1
      c = read c'
      (g1,g2) = split $ mkStdGen 4711
  putStrLn ("a="++show abits++" b="++show bbits++" c="++show c)
  fc <- mk_counter abits bbits
  mapM_ (add_count fc) $ take c $ zipWith Key (randomRs (0,amax) g1) (randomRs (0,bmax) g2)
  let pr i = print =<< mapM (get_count fc) [Key i m | m <- [0..min 20 bmax]]
  mapM_ pr [0..min 20 amax]
