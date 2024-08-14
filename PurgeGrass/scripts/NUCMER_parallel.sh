#!/bin/bash

# 13.03.2022

# do nucmer in parallel

# positional argument 1: reference

reference=$1
threads=$2

# first extract primary sequences
# make the primary folder
if [ -d primary ] 
then
    rm -r primary
    mkdir primary 
else
    mkdir primary 
fi

# make the name list for each primary contig
parallel -j ${threads} echo {} '>' primary/{}.lst :::: <(cut -f1 -d'&' collinear_block_ID.txt)

# now extract each primary contig
parallel -j ${threads} seqtk subseq ${reference} primary/{}.lst '>' primary/{}.fasta :::: <(cut -f1 -d'&' collinear_block_ID.txt)

# second extract secondary sequences
# make the secodnary folder
if [ -d secondary ]
then
    rm -r secondary
    mkdir secondary
else
    mkdir secondary
fi

# make the name list for each secondary contig
parallel -j ${threads} echo {} '>' secondary/{}.lst :::: <(cut -f2 -d'&' collinear_block_ID.txt)

# now extract each primary contig
parallel -j ${threads} seqtk subseq ${reference} secondary/{}.lst '>' secondary/{}.fasta :::: <(cut -f2 -d'&' collinear_block_ID.txt)

# third do nucmer alignment between primary and secondary in parallel
# make a folder to store nucmer results
if [ -d nucmer ]
then
    rm -r nucmer
    mkdir nucmer
else
    mkdir nucmer
fi

# align
parallel -j ${threads} --link nucmer -t 1 --delta=nucmer/{1}.delta primary/{2}.fasta secondary/{3}.fasta :::: <(cat collinear_block_ID.txt) :::: <(cut -f1 -d'&' collinear_block_ID.txt) :::: <(cut -f2 -d'&' collinear_block_ID.txt)

# filter
parallel -j ${threads} delta-filter -g nucmer/{}.delta '>' nucmer/{}_g.delta :::: <(cat collinear_block_ID.txt)

# convert to coords
parallel -j ${threads} show-coords -HTcl nucmer/{}_g.delta '>' nucmer/{}.coords :::: <(cat collinear_block_ID.txt)

# finally concatenate all coords to one coords
cat nucmer/*.coords > ps.coords
