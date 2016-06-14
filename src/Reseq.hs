module Reseq where

{- Use one or more kmer-indices to resequence.
   Inputs can be
   - single (low-quality) reads (PacBio, 454)
   - paired-ends (filling in the gap in the middle)
   - mate-pairs (jumping libs, fosmid- or bac-ends)
   - protein sequences

Cycles: resolved by statistical distribution of fragment sizes

-}

import Kmers
import FreqCount (FreqCount(..))

import Bio.Core.Sequence
import Data.ByteString.Lazy.Char8 (index)

import System.IO.Unsafe
import Data.List (sortBy, intercalate)
import Data.Function (on)
import Data.Bits
import Data.Char (toLower)
import qualified Data.PQueue.Min as H

import Debug.Trace

{-# ANN module "HLint: ignore Use camelCase" #-}

{- IMPORTANT: To get speed, I think we need to cache (position,kmer) info
   (dynamic programming) so that we don't calculate the same (remaining)
   path many times

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

-- This is highly questionable with the Eq instance above, isn't it?
instance Ord Path where compare = compare `on` prob

instance Show Path where
  show (Path pr p rc ws) = "Path "++show pr++" "++show p++" "++intercalate ":" (map (unkmer 5) (take 4 ws))

-- bug: this will display the full kmer, even if that is longer than the max specified with -m
showpath k (Path pr p rc (w:ws)) = "path: "++show pr++"\n"++reverse (map (last . unkmer k) ws) ++ unkmer k w

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
               { subst = \i x -> let c = unSD sd `index` fromIntegral i
                                 in -- trace ("c="++show c) $
                                    if toLower c == toLower x then match ss
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
                                        0 -> (undefined,6) -- all equally likely
                                        1 -> (0,20) -- one taken, rest penalized (1% chance)
                                        2 -> (3,20) -- 50% for each of two, rest 1%
                                        3 -> (5,20) -- 33% and 1%
                                        4 -> (6,undefined) -- again equally likely
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

type Paths = H.MinQueue Path

paths :: Eval -> (Pos,Kmer,Kmer) -> Pos -> [Path] -- iterate until position is reached, return lazy list of all paths
paths e (p0,start,startrc) endpos = go (H.singleton (Path { prob = 0, pos = p0, rc = startrc, path = [start] }))
  where go ps
          | H.null ps = []
          | otherwise = let p = H.findMin ps in if pos p >= endpos then p : go (H.deleteMin ps) else go (step e ps)
        
paths_pe :: a -> [Path]
paths_pe = undefined -- iterate from both ends, until met, and add penalty for distance (first may not be best?)

-- ------------------------------------------------------------------------
-- Utility functions
-- ------------------------------------------------------------------------

-- Pop the top, expand it, and merge it back into the list           
step :: Eval -> Paths -> Paths
step e ps' | H.null ps' = error "Empty starting path for 'step'."
           | otherwise = let (p,ps) = H.deleteFindMin ps' in H.union (expand e p) ps

expand :: Eval -> Path -> Paths
expand sd pt@(Path p i wr ws) = -- trace ("Current: "++show pt) $ -- ++show substs++show dels++show inss) $
  H.fromList (substs++inss++dels)
  where
    -- substitution or match - move forward and add one letter to path
    substs = [ Path (p+pr+subst sd i x) (i+1) wr (w:ws) | (x,w,pr) <- nextmers ]
    -- deletion in reference - add one letter to path but don't progress in reference
    dels   = [ Path (p+pr+ins sd i x) i wr (w:ws) | (x,w,pr) <- nextmers ]
    nextmers = next sd (head ws) wr 
    -- insertion in reference - progress reference but don't add letter to path
    inss   = [ Path (p+del sd i) (i+1) wr ws ]
                                
