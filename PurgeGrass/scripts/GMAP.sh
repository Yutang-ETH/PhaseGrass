#!/bin/bash

# map cds or transcripts to the reference

# positional argument 1: path of the reference
# positional argument 2: path of the trnascripts
# positional argument 3: number of threads

reference=$1
transcripts=$2
threads=$3

# first build the database
gmap_build -D gmap_test -d test $reference

if [ $(grep -v '^>' ${reference} | wc -c) -ge 4294967296 ]
then
     # map transcripts
     echo 'large genome >= 2^32 bp, use gmapl'
     gmapl -D gmap_test \
     -d test \
     --no-chimeras \
     -t $threads \
     -f 2 \
     -n 1 \
     --gff3-add-separators=0 \
     --min-identity=0.90 \
     --min-trimmed-coverage=0.90 \
     --max-intronlength-middle=20000 \
     $transcripts > gmap.gff
else
    # map transcripts
    echo 'small genome < 2^32 bp, use gmap'
    gmap -D gmap_test \
     -d test \
     --no-chimeras \
     -t $threads \
     -f 2 \
     -n 1 \
     --gff3-add-separators=0 \
     --min-identity=0.90 \
     --min-trimmed-coverage=0.90 \
     --max-intronlength-middle=20000 \
     $transcripts > gmap.gff
fi


