#!/bin/bash

# trim fasta from right, keep left end

fasta=$1
name=$2
output=$3
hbl=$4

echo $name > $name.list

trim=$(grep "$name" $hbl | cut -f2 | head -1) 

seqtk subseq $fasta $name.list > $name.fa
seqtk trimfq -e $trim $name.fa > $output/$name.fa

rm $name.list $name.fa
