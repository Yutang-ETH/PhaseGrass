#!/bin/bash

# positional argument 1: reference

reference=$1

# first using seqtk to extract DNA sequences of primary contigs and secondary contigs
seqtk subseq $reference primary_contigs.txt > primary_contigs.fasta
seqtk subseq $reference secondary_contigs.txt > secondary_contigs.fasta

# do all by all aligment between primary_contig.fasta and secondary_contig.fasta using NUCMER
nucmer -t 48 --delta=ps.delta primary_contigs.fasta secondary_contigs.fasta

# the longest alignment chain
delta-filter -g ps.delta > ps_g.delta

# convert delta format to coords format
show-coords -HTcl ps_g.delta > ps.coords

