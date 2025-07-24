#!/bin/bash

# the goal of this script is the following: 
# construct micro-synteny between allelic contigs uing MCScanX
# take 1) all-by-all alignment results from purge haplotigs, 2) BUSCO results table and 3) micro-synteny between allelic contigs to purge the shorter cllelic contig 

set -euo pipefail

# defualt thread
mythread=4

usage(){
>&2 cat << EOF
Usage:
   
   ./PurgeGrass.sh -a ./asm.fasta -b ./busco.table -f ./asm.fasta.fai -g ./transcript.fasta -m ./MCScanX -p ./curated.contig_associations.log -s ./scripts -t 24

options:
   [ -a | --assembly          required, path to the fasta assembly you want to purge ]
   [ -b | --busco             required, path to the busco table ]
   [ -f | --fai               required, path to the fai index of the assembly ]
   [ -g | --gene              required, path to a fasta transcript/cds input, the genes used to construct micro-synteny ]
   [ -m | --mcscanx           required, path to mcscanx program ] 
   [ -p | --PurgeHaplotig]    required, path to the Purge Haplotig curated.contig_associations.log ]
   [ -s | --scripts           required, path to the folder with scripts of the PurgeGrass pipeline ]
   [ -t | --thread            optional, number of threads, default 4, for each multi-thread step in the PurgeGrass pipeline ]
   [ -h | --help              show usage and exit ]
EOF
exit 1
}

# if no options are given, show usage page and exit
if [ $# -eq 0 ]
then
    echo "Please provide the valide input files and paths"
    usage
fi

args=$(getopt -o a:b:f:g:m:p:s:t:h --long assembly:,busco:,fai:,gene:,mcscanx:,PurgeHaplotig:,scripts:,thread:,help -- "$@")

if [[ $? -gt 0 ]]; then
  usage
fi

eval set -- ${args}

while :
do
  case $1 in
    -a | --assembly) 
    myasm=$2 
    shift ;;
    -b | --busco)    
    mybusco=$2
    shift ;;
    -f | --fai)
    myfai=$2
    shift ;;
    -g | --gene)
    mygene=$2
    shift ;;
    -m | --mcscanx)
    mymcscanx=$2
    shift ;;
    -p | --PurgeHaplotig)
    mypurgehaplotig=$2
    shift ;;
    -s | --scripts)
    myscripts=$2
    shift ;;
    -t | --thread)   
    mythread=$2
    shift ;;
    -h | --help)
    usage
    shift ;;
    --) shift; break ;;
    *) >&2 echo Invalida option for $(basename $0)
       usage ;;
  esac
  shift
done

# check if all required files are provided
if [ ! -s ${myasm} ] || [ ! -s ${mybusco} ] || [ ! -s ${myfai} ] || [ ! -s ${mygene} ] || [ -z ${mymcscanx} ] || [ ! -s ${mypurgehaplotig} ] || [ -z ${myscripts} ]
then
    echo 'Invalid files or paths, please check the usage'
    usage
else
    echo "Seems that all input files and paths are provided, run the pipeline"
fi

reference=${myasm}
fai=${myfai}
transcripts=${mygene}
haplotig_log=${mypurgehaplotig}
busco=${mybusco}
bin=${myscripts}
threads=${mythread}
mcscanx=${mymcscanx}


#-------------------------------------------------------------------------#

# first step, run gmap, this will output a gmap.gff in currrent directory
if [ -f GMAP.success ] && [ -s gmap.gff ]
then
    echo "Found GMAP.success and gmap.gff, skip GMAP.sh, continue"
else
    echo 'Map genes to the assembly by GMAP.sh with ${threads} threads'
    ${bin}/GMAP.sh ${reference} ${transcripts} ${threads}
fi

if [ $? -eq 0 ] && [ -s gmap.gff ]
then
    touch GMAP.success
else
    echo "GMAP.sh failed, stop pipeline"  
    exit 1
fi

#------------------------------------------------------------------------#

# second step, convert gmap.gff to gmap_protein.fasta
if [ -f GFFREAD.success ] && [ -s gmap_protein.fasta ]
then
    echo "Found GFFREAD.success and gmap_protein.fasta, skip GFFREAD, continue"
else
    echo "Convert mapped genes to proteins by GFFREAD.sh"
    ${bin}/GFFREAD.sh ${reference}
fi

if [ $? -eq 0 ] && [ -s gmap_protein.fasta ]
then 
    touch GFFREAD.success
else
    echo "GFFREAD.sh failed, stop pipeline"  
    exit 1
fi

#------------------------------------------------------------------------#

# third step, make mcscanx.gff and make new header names for gmap_protein.fasta
if [ -f NEW_header.success ] && [ -s gmap_protein_A.fa ] && [ -s gmap_protein_B.fa ]
then
    echo "Found protein files with new headers, skip reheader, continue"
else
    echo "Generate two protein files with different headers by MAKE_mcscanx_gff_and_new_protein_header.sh"
    ${bin}/MAKE_mcscanx_gff_and_new_protein_header.sh ${bin}
fi

if [ $? -eq 0 ] && [ -s gmap_protein_A.fa ] && [ -s gmap_protein_B.fa ]
then
    touch NEW_header.success
else
    echo "MAKE_mcscanx_gff_and_new_protein_header.sh failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# fourth step, do all-by-all protein alignment with diamond
if [ -f DIAMOND.success ] && [ -s mcscanx.blast ]
then
    echo "Found DIAMOND.success and mcscanx.blast, skip all-by-all protein alignment, continue"
else
    echo "Do all-by-all protein alignment by DIAMOND.sh with ${threads} threads"
    ${bin}/DIAMOND.sh ${threads}
fi

if [ $? -eq 0 ] && [ -s mcscanx.blast ]
then
    touch DIAMOND.success
else
    "DIAMOND.sh failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# fifth step, do mcscanx
if [ -f MCSCANX.success ] && [ -d mymcscanx ] && [ -s mymcscanx/mcscanx.collinearity ]
then
    echo "Found MCSCANX.success and mymcscanx/mcscanx.collinearity, skip MCScanX, continue"
else
    echo "Construct micro-synteny by MCSCNX.sh"
    ${bin}/MCSCANX.sh ${mcscanx}
fi

if [ $? -eq 0 ] && [ -d mymcscanx ] && [ -s mymcscanx/mcscanx.collinearity ]
then
    touch MCSCANX.success
else
    echo "Construct micro-synteny by MCSCNX.sh failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# sixth step, extract collinear blocks from mcscanx results
if [ -f EXTRACT_pairs.success ] && [ -s mymcscanx/mcscanx_cb_position.txt ]
then
    echo "Found EXTRACT_pairs.success and mymcscanx/mcscanx_cb_position.txt, skip extract paired contigs, continue"
else
    echo "Extract paired contigs based on micro-synteny by EXTRACT_pairs.sh"
    ${bin}/EXTRACT_pairs.sh
fi

if [ $? -eq 0 ] && [ -s mymcscanx/mcscanx_cb_position.txt ]
then
    touch EXTRACT_pairs.success
else
    echo "Extract paired contigs based on micro-synteny by EXTRACT_pairs.sh, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# seventh step, visualize collinear blocks
if [ -f PLOT_BLOCK.success ] && [ -s mymcscanx/mcscanx_linux_cb.pdf ] && [ -s primary_contigs.txt ] && [ -s secondary_contigs.txt ] && [ -s collinear_block_ID.txt ]
then
    echo "Found PLOT_BLOCK.success, skip plot collinear block, continue"
else
    echo "Plot collinear block pairs"
    Rscript ${bin}/plot_cb.R ${fai} ${bin}
fi

if [ $? -eq 0 ] && [ -s mymcscanx/mcscanx_linux_cb.pdf ] && [ -s primary_contigs.txt ] && [ -s secondary_contigs.txt ] && [ -s collinear_block_ID.txt ]
then
    touch PLOT_BLOCK.success
else
    echo "Plot collinear block pairs failed, stop pipeline" 
    exit 1
fi

#------------------------------------------------------------------------#

# do all-by-all alignment between primary-secondary contigs to find non-overlap sequences
# last step produces primary_contigs.txt and secondary_contigs.txt
if [ -f NUCMER_parallel.success ] && [ -s ps.coords ]
then
    echo "Found NUCMER_parallel.success, skip all-by-all alignment between primary and sencondary contigs, continue"
else
    echo "All-by-all alignment between primary and sencondary contigs"
    ${bin}/NUCMER_parallel.sh ${reference} ${threads}
fi

if [ $? -eq 0 ] && [ -s ps.coords ]
then
    touch NUCMER_parallel.success
else
    echo "all-by-all alignment between primary and sencondary contigs failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#  

# trim haplotigs
if [ -f TRIM_haplotigs.success ] && [ -f trimmed.fa ]
then
    echo "Found TRIM_haplotigs.success, skip trim haplotigs, continue"
else
    echo "Trim haplotigs to get non-overlap sequences by TRIM_haplotigs.sh"
    ${bin}/TRIM_haplotigs.sh ${reference} ${bin}
fi

if [ $? -eq 0 ] && [ -f trimmed.fa ]
then
    touch TRIM_haplotigs.success
else
    echo "Trim haplotigs failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#

# get haplotigs from purge haplotigs
if [ -f PARSE_purgehaplotig.success ] && [ -s purge_haplotig_primary_contigs.txt ]
then
    echo "Found PARSE_purgehaplotig.success, skip parsing purge haplotig reuslts, continue"
else
    echo "Parse purge haplotig reuslts"
    awk 'NF' $haplotig_log | cut -f2 -d'>' | cut -f1 -d',' > purge_haplotigs.txt && grep 'PRIMARY' $haplotig_log | cut -f1 -d',' > purge_haplotig_primary_contigs.txt
fi

if [ $? -eq 0 ] && [ -s purge_haplotig_primary_contigs.txt ]
then
    touch PARSE_purgehaplotig.success
else
    echo "Parse purge haplotig reuslts failed, stop pipeline"
    exit 1
fi

#------------------------------------------------------------------------#
# check with busco table and output final primary and allelic contigs
if [ -f TUNE_haplotigs.success ] && [ -s final_suspect_haplotigs.txt ] && [ -s final_primary_contigs.txt ]
then
    echo "Found TUNE_haplotigs.success, skip final tuning of haplotigs, continue"
else
    echo "Finall tuning of haplotigs with evidence from BUSCO, purge haplotigs and micro-synteny by tune_busco_and_output_final_haplotig.R"
    Rscript ${bin}/tune_busco_and_output_final_haplotig.R ${fai} secondary_contigs.txt primary_contigs.txt purge_haplotigs.txt ${busco} purge_haplotig_primary_contigs.txt ${bin}
fi

if [ $? -eq 0 ] && [ -s final_suspect_haplotigs.txt ] && [ -s final_primary_contigs.txt ]
then
    touch TUNE_haplotigs.success
else
    echo "Fina tuning of haplotigs failed, stop pipeline"
    exit 1
fi 
 
#------------------------------------------------------------------------#

# make haplotig.fa and primary.fa
if [ -f FINAL_fasta.success ] && [ -s final_primary_with_trim.fa ] && [ -s final_primary_with_trim.fa.stats ]
then
    echo "Found FINAL_fasta.success and final_primary_with_trim.fa, pipeline completed, congratulations"
else
    echo "Generate the final primary assembly"
    ${bin}/make_final_fasta.sh $reference
fi

if [ $? -eq 0 ] && [ -s final_primary_with_trim.fa ] && [ -s final_primary_with_trim.fa.stats ]
then
    touch FINAL_fasta.success && echo "Found FINAL_fasta.success and final_primary_with_trim.fa, pipeline completed, congratulations"
else
    echo "Generate the final primary assembly failed, stop pipeline"
    exit 1
fi

#-------------------------------------------------------------------------#
