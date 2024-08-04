#!/bin/bash

# 14.04.2023

# busco of Sikem primary assembly

assembly="/scratch/yutang/sikem/hifiasm_ctg/sikem.asm.p_ctg.fasta"
database="/scratch/yutang/busco_database/embryophyta_odb10"

mkdir unphased_genome
busco -i ${assembly} \
      -l ${database} \
      -m genome \
      -o unphased_genome \
      -f \
      -c 60 \
      --offline
