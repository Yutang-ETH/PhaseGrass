#!/bin/bash

# trim hbl

hbl=$1
fasta=$2
output=$3

mkdir $output

cut -f1 $hbl > hbl_contig

cat hbl_contig | parallel -j 20 -k "bash bin/TRIM_from_e.sh $fasta {} $output $hbl"

 



