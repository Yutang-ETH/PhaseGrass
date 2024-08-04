#!/bin/bash

reference=$1

seqtk subseq $reference final_suspect_haplotigs.txt > final_haplotigs.fa
seqtk subseq $reference final_primary_contigs.txt > final_primary_contigs.fa

if [ -s trimmed.fa ] 
then
    cat final_primary_contigs.fa trimmed.fa > final_primary_with_trim.fa
else
    cat final_primary_contigs.fa > final_primary_with_trim.fa
fi

assembly-stats final_primary_with_trim.fa > final_primary_with_trim.fa.stats 
