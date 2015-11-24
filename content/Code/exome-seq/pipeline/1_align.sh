#! /usr/bin/env bash

#BSUB -J align[1-8]
#BSUB -e logs/align.%I.%J.err
#BSUB -o logs/align.%I.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 12

set -exo pipefail

FASTQ_DIR=data
# XXX export PICARD env from modules file
picard=/vol1/software/modules-sw/Picard/build/picard-tools-1.83/

# REF=/vol3/home/jhessel/ref/genomes/hg19/hg19.fa
REF=GRCh37/genome.fa

samples=(ll-1 ll-3 ll-5 ll-7 lll-10 lll-1 lll-2 lll-9)
sample=${samples[$(($LSB_JOBINDEX - 1))]}

R1=$FASTQ_DIR/${sample}_R1_001.fastq.gz
R2=$FASTQ_DIR/${sample}_R2_001.fastq.gz

##########
# FASTQC #
##########

mkdir -p exome-qc/$sample

# run only on 2nd read
# fastqc -o exome-qc/$sample/ $R2

#############
# Alignment #
#############
result=results/alignment
if [[ ! -d $result ]]; then
    mkdir -p $result
fi

if [[ ! -f $sample.bam ]]; then

    bwa mem -R '@RG\tID:'$sample'\tSM:'$sample \
        -M -U 40 $REF $R1 $R2 -t 12 \
        | samtools view -bS - \
        | samtools sort - $result/$sample

    samtools index $result/$sample.bam

fi

################################
# Mark/Remove (PCR) Duplicates #
################################

if [[ ! -f $result/$sample.dedup.bam ]]; then

    java -Xmx4g -jar $picard/MarkDuplicates.jar \
        I=$result/$sample.bam \
        O=$result/$sample.dedup.bam \
        M=$result/$sample.dedup.metrics.txt \
        ASSUME_SORTED=true

    samtools index $result/$sample.dedup.bam

fi

