#!/bin/bash

mkdir mykat
kat comp -t 30 -m 23 -h -H 50000000000 -o mykat/ '/home/yutachen/public/Yutangchen/Rabiosa_data/tmp/Ryegrass_SG700bp_trimmed_?P.fq' final_primary_with_trim.fa
kat plot spectra-cn -x 200 -o mykat/Rabiosa_canu_hap_kat.png mykat/-main.mx


