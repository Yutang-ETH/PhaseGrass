#!/bin/bash

# using purge_haplotig to remove redundant haplotigs in sikem pctg hifiasm assembly

# 14.04.2023

# step 1, align ONT long reads to the assembly
minimap2 -t 60 -I 5G -K 5G -ax map-hifi sikem.asm.p_ctg.fasta ../all_hifi.fastq.gz --secondary=no | samtools sort -m 3G -@ 20 -o aligned.bam

# step 2, calculate read depth coverage 
purge_haplotigs hist -b aligned.bam -g sikem.asm.p_ctg.fasta -t 48

# step 3, flag haplotigs
purge_haplotigs cov -i aligned.bam.gencov -l 5 -m 37 -h 180 -j 101 -s 80

# step 4, purge haplotigs
purge_haplotigs purge -t 60 -g sikem.asm.p_ctg.fasta -c coverage_stats.csv -a 70 -I 6G

