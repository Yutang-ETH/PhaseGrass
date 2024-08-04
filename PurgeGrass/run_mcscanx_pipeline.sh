#!/bin/bash

# this is a wrapper to run mcscanx pipeline

# set path to reference, transcripts and reference.fai
# actually, it is not transcripts but gene models (DNA sequences of genes)

reference=sikem.asm.p_ctg.fasta
fai=sikem.asm.p_ctg.fasta.fai
transcripts=Rabiosa_v1_maker_transcripts_Lmu01_slp.fa
haplotig_log=curated.contig_associations.log
busco=full_table.tsv

echo "Hello"

#-------------------------------------------------------------------------#

# first step, run gmap, this will output a gmap.gff in currrent directory

./bin/GMAP.sh $reference $transcripts 60

#------------------------------------------------------------------------#

# second step, convert gmap.gff to gmap_protein.fasta

if [ $? -eq 0 ] && [ -s gmap.gff ]
then 
    echo "GMAP.sh finished successfully, extract protein sequences"
    ./bin/GFFREAD.sh $reference
else
    echo "GMAP.sh failed, stop pipeline"  
    exit 1
fi

#------------------------------------------------------------------------#

# third step, make mcscanx.gff and make new header names for gmap_protein.fasta

if [ $? -eq 0 ] && [ -s gmap_protein.fasta ]
then
    echo "GFFREAD.sh finished successfully, make gff for mcscanx and make protein files for all-by-all blast"
    ./bin/MAKE_mcscanx_gff_and_new_protein_header.sh
else
    echo "GFFREAD.sh failed, stop pipeline" 
    exit 1
fi

#------------------------------------------------------------------------#

# fourth step, do all-by-all protein alignment with diamond

if [ $? -eq 0 ] && [ -s mcscanx.gff ] && [ -s gmap_protein_A.fa ] && [ -s gmap_protein_B.fa ]
then
    echo "MAKE_mcscanx_gff_and_new_protein_header.sh finished successfully, do all-by-all protein blast"
    ./bin/DIAMOND.sh
else
    echo "MAKE_mcscanx_gff_and_new_protein_header.sh failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# fifth step, do mcscanx

if [ $? -eq 0 ] && [ -s mcscanx.blast ] 
then
    echo "DIAMOND.sh finished successfully, do mcscanx"
    ./bin/MCSCANX.sh
else
    echo "DIAMOND.sh failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# sixth step, extract collinear blocks from mcscanx results

if [ $? -eq 0 ] && [ -d mymcscanx ] && [ -s mymcscanx/mcscanx.collinearity ]
then
    echo "MCSCANX.sh finished successfully, extract allelic contig pairs or micro collinear blocks"
    ./bin/EXTRACT_paris.sh
else
    echo "MCSCANX.sh failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# seventh step, visualize collinear blocks

if [ $? -eq 0 ] && [ -s mymcscanx/mcscanx_cb_position.txt ]
then
    echo "EXTRACT_paris.sh finished successfully, visualize collinear blocks"
    Rscript bin/plot_cb.R $fai
else
    echo "EXTRACT_paris.sh failed, stop pipeline" 
    exit 1
fi

#------------------------------------------------------------------------#

# check if visualization finished successfully

if [ $? -eq 0 ]
then
    echo "mcscanx pipeline finished successfully, congratulations from Yutang"
else
    echo "somehting might be wrong with the visualization step, check it please"
    exit 1
fi

#------------------------------------------------------------------------#

# do all-by-all alignment between primary-secondary contigs
# last step produces primary_contigs.txt and secondary_contigs.txt

if [ $? -eq 0 ] && [ -s primary_contigs.txt ] && [ -s secondary_contigs.txt ] && [ -s collinear_block_ID.txt ]
then
    echo "got primary and secondary contig list, now do all-by-all alignment between primary and secondary contigs using NUCMER"
    ./bin/NUCMER_parallel.sh $reference
else
    echo "something must be wrong with the last step, please check it"
    exit 1
fi

#------------------------------------------------------------------------#  

# trim haplotigs

if [ $? -eq 0 ] && [ -s ps.coords ]
then
    echo "all-by-all alignment between primary and secondary contigs finished successfully, now trim haplotigs"
    ./bin/TRIM_haplotigs.sh $reference
else
    echo "all-by-all alignment between primary and secondary contigs failed"
    exit 1
fi

#------------------------------------------------------------------------#

# get haplotigs from purge haplotigs

if [ $? -eq 0 ]
then
    echo "trim haplotigs was successful and now get haplotigs from purge haplotigs"
    awk 'NF' $haplotig_log | cut -f2 -d'>' | cut -f1 -d',' > purge_haplotigs.txt && grep 'PRIMARY' $haplotig_log | cut -f1 -d',' > purge_haplotig_primary_contigs.txt
else
    echo "trim haplotigs failed, please check"
    exit 1
fi

#------------------------------------------------------------------------#
# check with busco table and output final primary and allelic contigs

if [ $? -eq 0 ] && [ -s purge_haplotigs.txt ] && [ -s purge_haplotig_primary_contigs.txt ]
then
    echo "Got purge_haplotigs.txt and now check busco table and output primary and allelic contigs"
    Rscript ./bin/tune_busco_and_output_final_haplotig.R $fai secondary_contigs.txt primary_contigs.txt purge_haplotigs.txt $busco purge_haplotig_primary_contigs.txt
else
    echo
fi 
 
#------------------------------------------------------------------------#

# make haplotig.fa and primary.fa

if [ $? -eq 0 ] && [ -s final_suspect_haplotigs.txt ] && [ -s final_primary_contigs.txt ]
then
    echo "got final haplotig and primary contig list, now make the fasta file"
    ./bin/make_final_fasta.sh $reference
else
    echo "something is wrong with the last step, check it please"
    exit 1
fi

#-------------------------------------------------------------------------#

if [ $? -eq 0 ] && [ -s final_primary_with_trim.fa.stats ]
then
    echo "whole pipeline finished successfully, congratulations, goodbye"
else
    echo "something is wrong in the last step, please check it"
    exit 1
fi

