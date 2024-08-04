#!/bin/bash

# extract collinear blocks from mcscanx resutls

grep -v '#' mymcscanx/mcscanx.collinearity | sed 's/^[ \t]*//' > content.txt

cut -f 2 content.txt | cut -d "_" -f 1 > contig1
cut -f 2 content.txt | cut -d "_" -f 2 > start1
cut -f 2 content.txt | cut -d "_" -f 3 > end1

cut -f 3 content.txt | cut -d "_" -f 1 > contig2
cut -f 3 content.txt | cut -d "_" -f 2 > start2
cut -f 3 content.txt | cut -d "_" -f 3 > end2

paste contig1 start1 end1 contig2 start2 end2 > mymcscanx/mcscanx_cb_position.txt
rm contig1 start1 end1 contig2 start2 end2 content.txt 

