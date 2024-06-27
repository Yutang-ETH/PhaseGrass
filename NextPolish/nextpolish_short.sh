#!/bin/bash

# do the short read polishing manually after the long read polishing

# Set input and parameters
round=2 # do 2 rounds of short-read polishing. Each round includes two times of polishing, each time with a different algorithm. 
threads=48
read1="path/sample_R1.fastq.gz" # path to the pair-end 1 short reads
read2="path/sample_R2.fastq.gz" # path to the pair-end 2 short reads
input=genome.nextpolish.long.fa # path to the long-read polished primary assembly
nextpolish1="/home/yutachen/anaconda3/envs/nextpolish/share/nextpolish-1.4.1/lib/nextpolish1.py" # path to the nextpolish short-read polishing script

### do not edit lines below, unless you know what you are doing

for ((i=1; i<=${round};i++)); do
#step 1:
   #index the genome file and do alignment
   bwa index ${input};
   bwa mem -t ${threads} ${input} ${read1} ${read2} | samtools view --threads ${threads} -F 0x4 -b - | samtools fixmate -m --threads ${threads}  - - | samtools sort -m 2g --threads ${threads} - | samtools markdup --threads ${threads} -r - sgs.sort.bam
   #index bam and genome files
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${input};
   #polish genome file
   python ${nextpolish1} -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
   input=genome.polishtemp.fa;
#step2:
   #index genome file and do alignment
   bwa index ${input};
   bwa mem -t ${threads} ${input} ${read1} ${read2} | samtools view --threads ${threads} -F 0x4 -b - | samtools fixmate -m --threads ${threads}  - - | samtools sort -m 2g --threads ${threads} - | samtools markdup --threads ${threads} -r - sgs.sort.bam
   #index bam and genome files
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${input};
   #polish genome file
   python ${nextpolish1} -g ${input} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
   input=genome.nextpolish.fa;
done;

if [ $? -eq 0 ] && [ -s genome.nextpolish.fa ]
then
    echo 'Finally polished genome file: genome.nextpolish.fa'
else
    echo 'program died'
fi
