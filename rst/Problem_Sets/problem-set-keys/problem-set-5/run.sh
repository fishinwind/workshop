#! /usr/bin/env bash

#BSUB -J ps5
#BSUB -o ps5.%J.out
#BSUB -e ps5.%J.err

# run script for Problem Set 5 key

set -o nounset -o pipefail -o errexit -x

# XXX this will be user-specific
project="$HOME/devel/bio-workshop/_test/problem-set-keys/problem-set-5"

# input files
data="/vol1/opt/data"
fasta="$data/hg19.fa"
dnase_bed="$data/wgEncodeRegDnaseClusteredV2.bed.gz"
ctcf_bed="$data/wgEncodeUwTfbsHelas3CtcfStdPkRep1.narrowPeak.gz"
clustered_tfbs_bed="$data/wgEncodeRegTfbsClusteredV3.bed.gz"

result="$project/results/2014-03-03"

if [ ! -d $result ]; then
    mkdir -p $result
fi

# output files
ctcf_fasta="$result/ctcf.chr22.fa"
ctcf_chr22_bed="$result/ctcf.chr22.bed"
tfbs_no_dhs="$result/tfbs_no_dhs.bed"
tfbs_no_dhs_summary="$result/tfbs_no_dhs.summary.tab"
dhs_no_tfbs="$result/dhs_no_tfbs.bed"
dhs_no_tfbs_chr22="$result/dhs_no_tfbs.chr22.bed"
dhs_no_tfbs_fasta="$result/dhs_no_tfbs.chr22.fa"
problem_1_motifs="$result/problem_1_motifs"
problem_2_motifs="$result/problem_2_motifs"

# Problem 1.1, 1.2
zcat $ctcf_bed | awk '$1 == "chr22"' > $ctcf_chr22_bed
bedtools getfasta -fi $fasta -bed $ctcf_chr22_bed -fo $ctcf_fasta

if [[ ! -d $problem_1_motifs ]]; then
    meme -nmotifs 5 -dna -minw 6 -maxw 20 -maxsize 100000000 \
        -o $problem_1_motifs $ctcf_fasta
fi

# Problem 2
# Identify transcription factor binding peaks that do not overlap with
# DNase I hypersensitive sites

bedtools intersect -a $clustered_tfbs_bed -b $dnase_bed -v -sorted \
    > $tfbs_no_dhs

# What transcription factors are represented in these peaks?
cut -f4 $tfbs_no_dhs \
    | sort \
    | uniq -c \
    | sort -k1nr \
    > $tfbs_no_dhs_summary

# identify DNase I hypersensitive sites that do not have corresponding
# transcription factor peak calls
bedtools intersect -a $dnase_bed -b $clustered_tfbs_bed -v -sorted \
    > $dhs_no_tfbs

# What motifs are enriched in this set of hypersensitive sites?
awk '$1 == "chrY" && ($3 - $2 > 8)' < $dhs_no_tfbs > $dhs_no_tfbs_chr22
bedtools getfasta -fi $fasta -bed $dhs_no_tfbs_chr22 -fo $dhs_no_tfbs_fasta

if [[ ! -d $problem_2_motifs ]]; then
    meme -nmotifs 10 -dna -minw 6 -maxw 20 -maxsize 100000000 \
        -o $problem_2_motifs $dhs_no_tfbs_fasta
fi
