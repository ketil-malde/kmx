module Reseq where
import Data.Word

{- Use one or more kmer-indices to resequence.
   Inputs can be
   - single (low-quality) reads (PacBio, 454)
   - paired-ends (filling in the gap in the middle)
   - mate-pairs (jumping libs, fosmid- or bac-ends)
   - protein sequences

Cycles: resolved by statistical distribution of fragment sizes

-}

{-
  *********** Todo: ************
  - merge Seen into Paths
  - switch Seen from (Pos -> Kmers) around to (Kmers -> Pos) (for easier paired-end analysis)
  - add a "greed" option, giving a positive bonus to scores
  - add a "trailing" option, to cut threads falling too far behind
  - support Bloom filters

************ Bugs/Issues ***********
  - handle out of sequence conditions

-}

import Kmers
import FreqCount (FreqCount(..))

import Bio.Core.Sequence
import Data.ByteString.Lazy.Char8 (index)
import qualified Data.ByteString.Lazy.Char8 as B

import System.IO.Unsafe
import Data.List (sortBy, intercalate)
import Data.Function (on)
import Data.Bits
import Data.Char (toLower)

import qualified Data.PQueue.Min as H
import qualified Data.IntMap.Strict as S
import qualified Data.IntSet as SS

import Debug.Trace

wlog :: String -> a -> a
-- log = trace
wlog s x = x

{-# ANN module "HLint: ignore Use camelCase" #-}

{- IMPORTANT: To get speed, I think we need to cache (position,kmer) info
   (dynamic programming) so that we don't calculate the same (remaining)
   path many times  ACTUALLY: this is branch-and-bound.

   Order by max achievable score (add scores assuming matches all the way to end)!
   prove that this is admissible (easy?) and consisten (hard?)

   Maybe use Giegerich - array (pos,kmer) -> path

   The scoring model needs to produce non-negative numbers, we could use
   the score in bits (i.e. positive from matches), but then we can't do a
   most-probable-first search, since an early mismatch may lead to a better
   path overall.  We can (probably?) give the best option zero score, which
   means it is immediately expanded.
-}

-- some definitions
type Prob = Int
type Pos  = Int
type Kmer = Word

data Path = Path { prob :: Prob    --  Negative logprob
                 , pos  :: Pos     -- position in reference sequence
                 , rc  :: Kmer     -- ugly, but cache revcompl of current mer
                 , path :: [Kmer]
                 } deriving Eq

type Paths = H.MinQueue Path
type Seen  = S.IntMap SS.IntSet

-- This is highly questionable with the Eq instance above, isn't it?
instance Ord Path where compare = compare `on` prob

instance Show Path where
  show (Path pr p rc ws) = "Path "++show pr++" "++show p++" "++intercalate ":" (map (unkmer 5) (take 4 ws))

-- bug: this will display the full kmer, even if that is longer than the max specified with -m
showpath k (Path pr p rc ws) = "path: "++show pr++"\n"++reverse (go ws)
  where go (x0:x:xs) = last (unkmer k x0) : go (x:xs)
        go [x]    = reverse (unkmer k x)

-- Gives probabilities for operations at various positions, used for
-- scoring.  This can typically be created from a reference sequence and
-- a kmer-index, but could be anything, really.
data Eval = Eval { subst, ins :: Pos -> Char -> Prob  -- probability of char in position
                 , del :: Pos -> Prob       -- here we don't care about actual letter, but we could
                 , ksize :: Int                      -- kmer size (we need this?)
                 , next :: Kmer -> Kmer -> [(Char,Kmer,Prob)] -- Calculating options and probabilities to move forward
                 }

data SimpleScores = SimpleScores { match, unknown, mismatch, indel :: Prob  }

-- guesstimated illumina values 
simple_illumina, simple454 :: SimpleScores
simple_illumina = SimpleScores { match    =  0 -- 97% chance = 0.14
                               , unknown  =  6 -- 25% chance
                               , mismatch = 30 --  0.1% chance
                               , indel    = 40 -- 0.01% chance
                               }

-- 454 has higher prob of indels
simple454 = SimpleScores { match = 0, unknown = 6, mismatch = 25, indel = 15 }


-- should probably recurse, and use uncons on sd to avoid index lookup
-- should cache reversecompl kmer as well
simpleeval :: SimpleScores -> SeqData -> FreqCount -> Eval
simpleeval ss sd kx = Eval
               { subst = \i x -> let s = unSD sd
                                     c = s `index` fromIntegral i  -- here, i is the position being considered
                                 in -- trace ("c="++show c) $
                                    if fromIntegral i >= B.length s || toLower c == toLower x then match ss -- no penalty beyond sequence end
                                    else if c `notElem` "ACGTacgt" then unknown ss
                                         else mismatch ss
               , ins   = \_i _x -> indel ss 
               , del   = \_i    -> indel ss
               , ksize = my_k
               , next  = \w wr -> let xs = ['a','c','g','t']
                                      ws = [next1 my_k w c | c <- xs] 
                                      cs :: [Int] -- counts for scores
                                      cs = let gc = unsafePerformIO . get_count kx -- gc x = max 1 $ min 50 $ unsafePerformIO (get_count kx x) -- pseudocount and ceiling
                                           in map gc $ zipWith min ws [next1_rc my_k wr c | c <- xs] -- kmer index stores minimums
                                      cs' = map (>=3) cs -- threshold for accepted link

                                      f = case length (filter (==True) cs') of
                                        0 -> (undefined,2) -- all equally likely
                                        1 -> (0,20) -- one taken, rest penalized (1% chance)
                                        2 -> (1,20) -- 50% for each of two, rest 1%
                                        3 -> (1,20) -- 33% and 1%
                                        4 -> (1,undefined) -- again equally likely
                                        _ -> error "wot? say it ain't so"
                                      prs :: [Prob]         -- prs = [(10*) . round . negate . logBase 10 $ (fromIntegral c / fromIntegral (sum cs)) | c <- cs]
                                      prs = map (\x -> if x then fst f else snd f) cs'
                                  in -- trace (show $ map (unkmer 5) ws) $
                                   zip3 xs ws prs
               }
  where my_k     = key_bits kx
        
-- ------------------------------------------------------------------------
-- Traversing the graph.  These should all return [Path], I think.
-- ------------------------------------------------------------------------

paths :: Eval -> (Pos,Kmer,Kmer) -> Pos -> [Path] -- iterate until position is reached, return lazy list of all paths
paths e (p0,start,startrc) endpos = go S.empty (H.singleton (Path { prob = 0, pos = p0, rc = startrc, path = [start] }))
  where go sn ps
          | H.null ps = []
          | otherwise = let p = H.findMin ps in if pos p >= endpos then p : go sn (H.deleteMin ps)
                                                else let (s,p) = step e sn ps in go s p
        
paths_pe :: a -> [Path]
paths_pe = undefined -- iterate from both ends, until met, and add penalty for distance (first may not be best?)

-- ------------------------------------------------------------------------
-- Utility functions
-- ------------------------------------------------------------------------

-- Pop the top, expand it, and merge it back into the list           
step :: Eval -> Seen -> Paths -> (Seen,Paths)
step e sn ps' | H.null ps' = error "Empty starting path for 'step'."
              | otherwise = let (p,ps) = H.deleteFindMin ps'
                                (s,new) = expand e sn p
                            in (s,H.union new ps)

expand :: Eval -> Seen -> Path -> (Seen,Paths)
expand sd sn pt@(Path p i wr ws) = wlog ("Current: "++show pt) $ -- ++show substs++show dels++show inss) $
                                   wlog ("  candidates: "++show (map (\p -> (unkmer 5 . head $ path p,pos p,prob p)) unseen)) $
  (newpos, H.fromList unseen) -- hmm. does this work?
  where
    -- substitution or match - move forward and add one letter to path
    substs = [ Path (p+pr+subst sd i x) (i+1) wr (w:ws) | (x,w,pr) <- nextmers ]
    -- deletion in reference - add one letter to path but don't progress in reference
    dels   = [ Path (p+pr+ins sd i x) i wr (w:ws) | (x,w,pr) <- nextmers ]
    nextmers = next sd (head ws) wr 
    -- insertion in reference - progress reference but don't add letter to path
    inss   = [ Path (p+del sd i) (i+1) wr ws ]

    unseen :: [Path]
    unseen = filter unmember (substs++inss++dels)

    unmember :: Path -> Bool
    unmember x = case S.lookup (pos x) sn of
                  Just kms -> not ((fromIntegral $ head $ path x) `SS.member` kms)
                  Nothing  -> True

    newpos :: Seen
    newpos = go sn unseen
      where go s (x:xs) = case S.lookup (pos x) s of
              Just kms -> go (S.insert (pos x) (SS.insert (fromIntegral $ head $ path x) kms) s) xs
              Nothing  -> go (S.insert (pos x) (SS.insert (fromIntegral $ head $ path x) SS.empty)  s) xs
            go s [] = s
