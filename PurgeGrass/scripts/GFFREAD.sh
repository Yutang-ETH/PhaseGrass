#!/bin/bash

# convert gff3 to protein.fa

# positional argument 1: path to the reference

reference=$1

gffread gmap.gff -M -g $reference -y gmap_protein.fasta
