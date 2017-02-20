{-# Language DeriveDataTypeable #-}
{-# OPTIONS_GHC -fno-cse #-}

module Options (getArgs, Options(..), genOutputBS, genOutput) where

import System.Console.CmdArgs
import Data.Word
import qualified Data.ByteString.Lazy.Char8 as B

data Options = Count { kval :: Int, fasta :: Bool
                     , filter_bits :: Int, filter_value :: Word
                     , files :: [FilePath], output :: FilePath }
             | Hist  { indices :: [FilePath], output :: FilePath
                     , kval, mincount, maxcount :: Int
                     , diploid, stats :: Bool
                     , complexity_classes, complexity_mersize :: Int
                     }
             | Stats { output :: FilePath, histogram :: FilePath, diploid :: Bool }
             | Correlate { indices :: [FilePath]
                         , kval, mincount, maxcount :: Int
                         , sqrt_transform :: Bool
                         }
             | Verify { indices :: [FilePath], mincount, maxcount :: Int }
             | Dump { indices :: [FilePath], output :: FilePath
                    , kval, mincount, maxcount :: Int
                    , hashes :: Bool
                    , complexity :: Bool }
             | Merge { indices :: [FilePath]
                     , output :: FilePath
                     -- min/maxcount?
                     }
             | Heatmap { indices :: [FilePath]
                       , output :: FilePath
                       , kval, mincount, maxcount :: Int
                       }
             | Classify { indices :: [FilePath]
                       , output :: FilePath
                       , kval, mincount, maxcount :: Int
                       , files :: [FilePath]
                       , fasta, paired :: Bool
                       }
             | Jaccard { indices :: [FilePath]
                       , kval :: Int
                       } -- todo: kval, mincount, etc
             | Reseq   { indices :: [FilePath]
                       , kval, maxlen :: Int
                       , files :: [FilePath]
                       , fasta :: Bool
                       }
             deriving (Typeable,Data)

-- | Build a k-mer index and output to specified file name
def_count :: Options
def_count = Count { kval = 32 &= help "k-mer size to use"
                  , fasta = False &= help "input is Fasta format (default is FastQ)"
                  , output = "-" &= help "output kmer index file name" &= typFile
                  , filter_value = 0 &= help "construct a partial index, using only words matching the filter"
                  , filter_bits = 0 &= help "number of bits to use for filter"
                  , files = [] &= args &= typFile }
            &= details ["Generate the count index from an input file."]

-- | Read a k-mer index (or build it directly) and output a histogram of counts
def_hist :: Options
def_hist = Hist { -- kval = 0 &= help "k-mer size to use"
                  output = "-" &= help "output file name" &= typFile
                , kval = 0 &= help "k-mer size to reduce to"
                , mincount = 0 &= help "minimum count to include"
                , maxcount = 0 &= help "maximum count to include"
                , complexity_classes = 0 &= help "number of categories for k-mer complexity (entropy)" &= name "c"
                , complexity_mersize = 1 &= help "mersize to calculate complexity for"                 &= name "m"
                , diploid = False &= help "calculate statistics for diploid organism"
                , stats   = False &= help "calculate statistics (may use more memory)"                 &= name "s"
                , indices = [] &= typFile &= args
                }
           &= details ["Output a histogram of frequency counts."
                      ,"Each line contains a frequency and the number of distinct k-mers occurring with this frequency."
                      ,"With -s, output comment header with estimated k-mer distribution (similar to 'kmx stats')"
                      ]

def_verify :: Options
def_verify = Verify { indices = [] &= typFile &= args
                     , mincount = 0 &= help "minimum count to include"  -- useless here?
                     , maxcount = 0 &= help "maximum count to include"
                    }
             &= details ["Verify the correctness of a count index."]

def_corr :: Options
def_corr = Correlate { indices = [] &= args &= typFile
                     , kval = 0 &= help "k-mer size to reduce to"
                     , mincount = 0 &= help "minimum count to include"
                     , maxcount = 0 &= help "maximum count to include"
                     , sqrt_transform = False &= help "sqrt-transform data points"
                     }
           &= details ["Calculate the correlation coefficient (Pearson's r) between k-mer frequencies,"
                      ,"as well as the linear regression coefficients."]

def_heatmap :: Options
def_heatmap = Heatmap { indices = [] &= args &= typFile
                      , output =  "" &= typFile &= help "Output file"
                      , kval = 0 &= help "k-mer size to reduce to"
                      , mincount = 0 &= help "minimum count to include"
                      , maxcount = 200 &= help "maximum count to include"
                      }
              &= details ["Construct heatmap data for correlation of k-mers in two files."]

def_dump :: Options
def_dump = Dump { indices = [] &= args &= typFile
                , output =  "" &= typFile &= help "Output file"
                , kval = 0 &= help "k-mer size to reduce to"
                , mincount = 0 &= help "minimum count to include"
                , maxcount = 0 &= help "maximum count to include"
                , hashes = False &= help "output k-mers as raw hash values instead of text"
                , complexity = False &= help "also output k-mer complexity (entropy)"
                }
           &= details ["Dump the contents of a count index."]

def_merge :: Options
def_merge = Merge { indices = [] &= args &= typFile
                , output =  "" &= typFile &= help "Output file"
                } &= details ["Merge two or more indices."]

def_class :: Options
def_class = Classify { indices = [] &= args &= typFile
                     , output =  "" &= typFile &= help "Output file"
                     , kval = 0 &= help "k-mer size to reduce to"
                     , fasta = False &= help "input is Fasta format (default is FastQ)"
                     , files = [] &= typFile &= help "sequences to classify"
                     , paired = False &= help "read sequence pairs"
                     , mincount = 0 &= help "minimum count to include"
                     , maxcount = 0 &= help "maximum count to include"
                     } &= details ["Classify fastq or fasta sequences by k-mer count spectrum"]

def_jacc :: Options
def_jacc = Jaccard { indices = [] &= args &= typFile
                   , kval = 0 &= help "k-mer size to reduce to"
                   }  &= details ["Calculate a probability-based jaccard distance between two k-indices"]

def_reseq :: Options
def_reseq = Reseq { indices = [] &= typFile
                  , kval = 0 &= help "k-mer size to reduce to"
                  , maxlen = 500 &= help "target position length"
                  , files = [] &= typFile &= args
                  , fasta = False &= help "input is Fasta format (default is FastQ)"
                  }  &= details ["Resequence input data using kmer information"]

def_stats :: Options
def_stats = Stats { output =  "" &= typFile &= help "Output file"
                  , histogram = "" &= args &= typFile
                  , diploid = False &= help "organism is diploid"
                  } &= details ["Read histograms and estimate k-mer distribution parameters"]

getArgs :: IO Options
getArgs = checkopts `fmap` cmdArgsRun (cmdArgsMode $ modes [def_count, def_hist, def_verify, def_corr, def_dump, def_merge, def_heatmap, def_class, def_jacc, def_reseq, def_stats]
  &= summary "kmx v0.5 - tool for k-mers analysis in biological sequences.\n© Ketil Malde, 2014."
  &= program "kmx")

checkopts :: Options -> Options
checkopts opts@(Verify {})
  | otherwise        = opts
checkopts opts@(Correlate {})
  | length (indices opts) /= 2 = error "Correlate takes exactly two index files as parameters."
  | otherwise                  = opts
checkopts opts = opts

genOutput :: Options -> String -> IO ()
genOutput opts = case output opts of
  "" -> putStr
  "-" -> putStr
  f -> writeFile f

genOutputBS :: Options -> B.ByteString -> IO ()
genOutputBS opts = case output opts of
  "" -> B.putStr
  "-" -> B.putStr
  f -> B.writeFile f
