#!/bin/bash

mkdir busco_primary_with_trim

busco -i final_primary_with_trim.fa \
       -l /home/yutachen/public/Yutangchen/data_from_DLF/Columbus/busco_dir/embryophyta_odb10 \
       -m genome \
       -o busco_primary_with_trim \
       -f \
       -c 20 \
       --offline


