#!/bin/bash

# trim haplotigs if they have large non-overlapping part with the primary contig

# positional argument 1: reference

reference=$1
scripts=$2

Rscript ${scripts}/produce_trim_bed_list.R ${scripts}

# trim contigs using seqtk
if [ $? -eq 0 ] && [ -f hbl_hbr.bed ]
then
    bedtools getfasta -fi $reference -bed hbl_hbr.bed -fo trimmed.fa
else
    echo "produce_trim_bed_list.R failed"
    exit 1
fi

