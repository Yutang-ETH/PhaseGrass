#!/bin/bash

# trim hbr

hbr=$1
fasta=$2
output=$3

mkdir $output

cut -f1 $hbr > hbr_contig

cat hbr_contig | parallel -j 40 -k "bash bin/TRIM_from_b.sh $fasta {} $output $hbr"

