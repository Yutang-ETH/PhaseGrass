#!/bin/bash

# trim fasta from left, keep right end

fasta=$1
name=$2
output=$3
hbr=$4

echo $name > $name.list

trim=$(grep "$name" $hbr | cut -f2 | head -1)

seqtk subseq $fasta $name.list > $name.fa
seqtk trimfq -b $trim $name.fa > $output/$name.fa

rm $name.list $name.fa

