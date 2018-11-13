#!/bin/bash

location=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data
iid=${1}
bam=${2}
coverage=${3}
# low coverage = "alignment", high coverage = "high_coverage_alignment"

alignment=${location}/${iid}/${coverage}/${bam}
index=${location}/${iid}/${coverage}/${bam}.bai

read_length=`samtools view ${alignment} |head -n 1000 |awk '{sum+=length($10)} END{print sum/NR}'`
echo $read_length
echo "running idxstats"

samtools idxstats ${alignment} > ${coverage}/idxstats/${iid}.idxstats

echo "calculating depth of sequencing"

Rscript cal_mtdna_copy.R ${coverage}/idxstats/${iid}.idxstats ${read_length}



