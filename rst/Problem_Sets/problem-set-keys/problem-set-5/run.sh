#! /usr/bin/env bash

#BSUB -J ps5
#BSUB -o ps5.%J.out
#BSUB -e ps5.%J.err

# run script for Problem Set 5 key

# XXX this will be user-specific
project="$HOME/devel/bio-workshop/problem-set-keys/problem-set-5"

data="/vol1/opt/data"
fasta="$data/hg19.fa"
dnase_bed="$data/wgEncodeRegDnaseClusteredV2.bed.gz"
ctcf_bed="$data/wgEncodeUwTfbsHelas3CtcfStdPkRep1.narrowPeak.gz"
clustered_tfbs_bed="$data/wgEncodeRegTfbsClusteredV3.bed.gz"

result="$project/results/2014-03-03"
results="$result/result.tab"

if [ ! -d $result ]; then
    mkdir -p $result
fi

# Problem 1.1, 1.2
ctcf_fasta="$result/ctcf.chr22.fa"
ctcf_chr22_bed="$results.ctcf.chr22.bed"

awk '$1 == chr22' < $ctcf_bed > $ctcf_chr22_bed
bedtools getfasta -fi $fasta -bed $ctcf_chip -fo $ctcf_fasta

meme -nmotifs 5 -dna -minw 6 -maxw 20 -maxsize 100000000 \
    -o problem_1_motifs $ctcf_fasta

# Problem 2
# Identify transcription factor binding peaks that do not overlap with
# DNase I hypersensitive sites

tfbs_no_dhs="$result/tfbs_no_dhs.bed"
bedtools intersect -a $clustered_tfbs_bed -b $dnase_bed -v -sorted \
    > $tfbs_no_dhs

# What transcription factors are represented in these peaks?
tfbs_no_dhs_summary="$result/tfbs_no_dhs.summary.tab"
zcat $tfbs_no_dhs | cut -f4 | sort | uniq -c | sort -k1nr \
    > $tfbs_no_dhs_summary

# identify DNase I hypersensitive sites that do not have corresponding
# transcription factor peak calls
dhs_no_tfbs="$result/dhs_no_tfbs.bed"
bedtools intersect -a $dnase_bed -b $clustered_tfbs_bed -v -sorted \
    > $dhs_no_tfbs

# What motifs are enriched in this set of hypersensitive sites?
dhs_no_tfbs_fasta="$result/dhs_no_tfbs.fa"
bedtools getfasta -fi $fasta -bed $dhs_no_tfbs -fo $dhs_no_tfbs_fasta
meme -nmotifs 10 -dna -minw 6 -maxw 20 -maxsize 100000000 \
    -o problem_2_2_2_motifs $dhs_no_tfbs_fasta

    
