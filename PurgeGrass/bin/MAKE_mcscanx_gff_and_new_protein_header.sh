#!/bin/bash

# make the gff file required by mcscanx based on the gmap.gff

grep -w "mRNA" gmap.gff | cut -f1 > scaffold_name
grep -w "mRNA" gmap.gff | cut -f4 > start
grep -w "mRNA" gmap.gff | cut -f5 > end
paste -d "_" scaffold_name start end > gene
paste scaffold_name gene start end > mcscanx.gff
rm scaffold_name start end gene


# make new header names for gmap_protein.fa

grep -w "mRNA" gmap.gff | cut -f9 | cut -f1 -d ";" | sed -n 's/ID=/>/p' > oldnames
cut -f2 mcscanx.gff > newnames

paste oldnames newnames > names

rm oldnames newnames

sed -i 's/>//' names

python bin/CHANGE_fa_header.py
