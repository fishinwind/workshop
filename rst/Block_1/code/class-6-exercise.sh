#! /usr/bin/env bash

#BSUB -J workflow[1-4]                                                                                                              
#BSUB -e %J.%I.err                                                                                                                  
#BSUB -o %J.%I.out                                                                                                                  
#BSUB -R "select[mem>4] rusage[mem=4]"

# this section defines the sample names and grabs the appropriate                                                                   
# sample name for the current process
SAMPLENAMES=(SP1 SP2 SP3 SP4)                                                                                                       

# $LSB_JOBINDEX is set at runtime for the current process; in this                                                                  
# case you asked for 4 jobs, so it be a value between 1 and 4                                                                       
sample=${SAMPLENAMES[$(($LSB_JOBINDEX - 1))]}                                                                                       

# set up file names                                                                                                                 
data=/vol1/opt/data
bwtindex=/vol1/opt/data/hg19
fastq=$data/$sample.fq.gz
align=$sample.sam
coverage=$sample.coverage.gz

# run the workflow 
bowtie2 -x $bwtindex -U $fastq > $align

# now generate counts of each position that was aligned to; fields 3                                                                
# and 4 of the same file are chrom and pos.
grep -v '^@' $align \
    cut -f3,4 \
    | sort \
    | uniq -c \
    | awk 'BEGIN {OFS="\t"} {print $2,$3,$1}' \
    | gzip -c > $coverage
