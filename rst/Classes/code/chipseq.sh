#! /usr/bin/env bash

#BSUB -J chipseq
#BSUB -o %J.out
#BSUB -e %J.err

#
# load required modules before submitting this script
#
# $ module load bowtie2
# $ module load ucsc
# $ module load bedtools
# $ module load macs2
# $ module load meme
#

data=/vol1/opt/data
fasta=$data/hg19.fa
bwtindex=$data/hg19
chromsize=$data/hg19.chrom.sizes

# downloaded from ENCODE website - H3K4me3 ChIP in Hela cells
fastqfile=$data/wgEncodeBroadHistoneHelas3H3k4me3StdRawDataRep1.fastq.gz

expname="chipseq.hela.h3k4me3"
bamfile=$expname.bam

# run the alignment
bowtie2 -x $bwtindex -U $fastqfile \
    | samtools view -ShuF4 - \
    | samtools sort -o - aln.temp -m 8G \
    > $bamfile

# call peaks
macs2 callpeak --treatment $bamfile --name $expname

# this will generate a series of files, including one in
# ``narrowPeak`` format, an extended BED fromat:
peakbed=chipseq.hela.h3k4me3_peaks.narrowPeak

# find motifs
peakfasta=$expname.peaks.fa
bedtools getfasta -fi $fasta -bed $peakbed -fo $expname.peaks.fa
meme -nmotifs 100 -minw 6 -maxw 20 $peakfasta

# make coverage plots in bedgraph and bigwig
bedgraph="$expname.bg"
bigwig="$expname.bw"
bedtools genomecov -ibam $bamfile -g $chromsize -bg > $bedgraph
bedGraphToBigWig $bedgraph $chromsize $bigwig

# this needs to be all on one line:
# track type=bigWig
# bigDataUrl="http://amc-sandbox.ucdenver.edu/~jhessel/
