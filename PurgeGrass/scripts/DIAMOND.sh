#!/bin/bash

threads=$1

# do all-by-all protein alignment with diamond

sed -i 's/\.//g' gmap_protein_A.fa
sed -i 's/\.//g' gmap_protein_B.fa 

diamond makedb --in gmap_protein_A.fa -d gmap_protein_A.db
diamond blastp -d gmap_protein_A.db -q gmap_protein_B.fa --query-cover 90 -f 6 -o mcscanx.blast --very-sensitive -p ${threads} --quiet -e 0.00001 --id 90

sed -i -e 's/A_//' -e 's/B_//' mcscanx.blast
