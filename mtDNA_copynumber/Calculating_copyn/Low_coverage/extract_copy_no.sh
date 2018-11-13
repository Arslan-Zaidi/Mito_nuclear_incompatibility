#!/bin/bash

location=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data
iid=${1}
bam=${2}

alignment=${location}/${iid}/alignment/${bam}
index=${location}/${iid}/alignment/${bam}.bai

read_length=`samtools view ${alignment} |head -n 1000 |awk '{sum+=length($10)} END{print sum/NR}'`
#Counts DNA string/read in 10th column of SAM file in the first 1000 reads present.
echo $read_length
echo "running idxstats"

samtools idxstats ${alignment} > ${iid}.idxstats

echo "calculating depth of sequencing"

Rscript cal_mtdna_copy.R ${iid}.idxstats ${read_length}



