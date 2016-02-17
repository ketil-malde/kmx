# KMX - K-Mer indeXing

This is (yet another) k-mer indexing/counting utility for use on
biological sequence data.  In essence, it builds an index by reading
Fastq or Fasta files, which it then can perform various analyses on.

## Why yet another one?

So, why yet another k-mer indexing tool?  Don't we have enough with
jellyfish, khmer, turtle, kmc, and so on, and so forth?  Well, kmx has
some features that may or may not be unique.

1. I think it is the only tool to use Judy arrays to store the index.
This isn't so important in itself, but the underlying, highly tuned
C-libarary gives decent performance and low memory consumption.

2. It's exact, not probabilistic. Some of the faster and more
memory-frugal tools out there use heuristics to get that way. This is
probably a worthwhile trade-off in many cases, but exact is still nice.

3. It has a compact indexing format, using a differential coding to
save space.  A benchmark indexing 200GB of Fastq data gives an index
of 20GB.  In addition, the format is very easy to parse and generate
(see src/Serialize.hs for details).  At least jellyfish produces a bit
larger indexes in my tests.

4. It is single threaded (for now, at least), but can distribute
index-building over multiple processes.  This means you can get build
an index for almost any amount of data on almost any kind of
computer.  Those 200GB were indexed using 32 separate processes in
three and a half hours, the biggest memory footprint of a single
process was 6 gigabytes.  There's a script that will use GNU parallel, 
which wil let you distribute those processes to separate machines if
you want.

5. It has a bunch of analyses that may or may not be unique:
correlations, heatmaps, and (of course) histograms.  See '--help' for
details.  All of these (so far) run in constant (and rather little)
space.

6. Analyses can downscale the k-mer size without actually doing a
recount.  So if you build with a k-mer size of 32 (the current
maximum, limited by the available Judy data structure), you get
indexing with any k-mer size less than that for free.

7. It is written in Haskell, and is thus much more beautiful than its
competitors.  In retrospect, I should have named it Galatea.

## Installation

Clone the repository, enter the directory, and type 'cabal install'.

This presupposes that you have a fairly recent and working GHC
installation, and if you just did 'apt-get install ghc' for this, you
will probably need to 'apt-get cabal-install' and do a 'cabal update'
first.

## Usage

It's pretty simple.  To build an index for k-mers of size 31, do:

    kmx count -k 31 input.fastq more.fastq -o index.32
	
To get the archetypical histogram from this, you can then:

    kmx hist index.32 -o index.hist
	
Plotting is an optional extra, I tend to use Gnuplot, where you should
be able to do something like:

	set xrange [0:200]
    plot "index.hist" with lines
	
You'll probably want to adjust yrange as well, but it's probably a
good idea to look at the data first.

In all cases, slapping on '--help' should give you more information on
usage and options.

## Distributed/parallel indexing

If you are low on memory, you may want to split the indexing into,
say, eight sequential processes.  This is how:

    for x in {0..7}; do
       kmx count --filter-bits=3 --filter-value=$x input.fastq -o partial_index.$x
    done
	kmx merge partial_index.* -o index

Basically, --filter-bits lets you specify the granularity, and the
partial indices are then numbered from zero to 2^b-1.  Of course, if
you have enough memory but would like to utilize more CPUs, just run
all of these simultaneously.  I encourage you to look at the script
included, and modify it to suit your needs.

## More information

 - [http://biohaskell.org/Applications/kmc]()
 - [http://blog.malde.org/posts/frequency-counting.html]()
 - [http://blog.malde.org/posts/k-mer-counting.html]()
