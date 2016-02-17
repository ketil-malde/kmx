#!/bin/bash

# Script to build a k-mer index in parallel
# This builds 2^$BITS partial indices simultaneously and merges them afterwards

KMX="/usr/bin/time kmx"
BITS=4
MAX=$((2**$BITS-1))
FILES=$*
K=32

echo STARTED >> LOG
date >> LOG

# Building partial indices - if you are low on memory, you may want to run them
# sequentially.  One experiment with 225M 100bp reads (80Gb fastq) split 16 ways required
# 4-17 GB RAM for each index.  This will vary by sequencing quality, genome size,
# and probably many other parameters.

DIR=kmx$K:$BITS
mkdir -p $DIR

if type parallel > /dev/null 2>&1 ; then
    # If we have GNU parallel installed, we use it to distribute the processes
    # This can be used to distribute things to different machines as well
    for f in $(eval echo {0..$MAX}); do
	echo "$KMX count --filter-bits=$BITS --filter-value=$f -k $K $FILES -o $DIR/$f"
    done | parallel
else
    # ..else we just run everything directly in the shell
    for f in $(eval echo {0..$MAX}); do
	$KMX count --filter-bits=$BITS --filter-value=$f -k $K $FILES -o $DIR/$f &
    done
    wait
fi

date >> LOG

# Merging all the partial indices
$KMX merge $DIR/* -o kmx_index.$K
date >> LOG
echo END >> LOG
