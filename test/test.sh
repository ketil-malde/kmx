#!/bin/bash

log(){
    echo "$(tput setaf 2)$*$(tput sgr0)"
}
error(){
    echo "$(tput setaf 1)$*$(tput sgr0)"
    exit -1
}

# Run some tests for kmx functionality
cabal build || error "Build failed!"

KMX=dist/build/kmx/kmx
DIR=test_tmp.d

log "Test run started - $(date)"
mkdir -p $DIR

log "Testing count"
$KMX count -k 32 test/250K-reads.fastq -o $DIR/index.32
$KMX verify $DIR/index.32

# Build histograms
log "Testing histograms"
$KMX hist $DIR/index.32 > $DIR/32.hist
$KMX hist -k 28 $DIR/index.32 > $DIR/32-to-28.hist

# Test merge functionality
log "Testing partial builds and merge"
$KMX count -k 32 test/250K-reads.fastq --filter-bits=1 --filter-value=0 -o $DIR/index.32-0
$KMX count -k 32 test/250K-reads.fastq --filter-bits=1 --filter-value=1 -o $DIR/index.32-1
$KMX merge $DIR/index.32-0 $DIR/index.32-1 -o $DIR/index.32-merged
cmp $DIR/index.32-merged $DIR/index.32 || error "Merged partial indexes differ from directly built index"

# Dump
log "Testing dump"
$KMX dump $DIR/index.32 | head -100 > $DIR/32.dump100
$KMX dump -k 30 $DIR/index.32 | head -100 > $DIR/30.dump100

# verify k-mer lenghts
# verify that count-30 for XXXXX is count-32 for XXXXXyy 
tail -n +2 $DIR/30.dump100 | head -4 | while read kmer count; do
    total=$(grep "^$kmer" $DIR/32.dump100 | awk '{N+=$2} END {print N}')
    test "$total" -eq "$count" || error "Dump counts for $kmer don't match: $total and $count"
done

# tests for correlate, heatmap?
# check md5sums

log "Test run completed - $(date)"
