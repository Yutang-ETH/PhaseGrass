#!/bin/bash

# get the list of suggested repeats 

input=$1

if [ ${input##*.} = fa ] || [ ${input##*.} = fasta ]
then
    grep -w 'suggestRepeat=yes' $input | cut -f1,5 -d" " | sed 's/>//' > assembly_repeat.list
elif [ ${input##*.} = txt ]
then
   cut -f1,5 $input | grep "Y" > assembly_repeat.list
fi
