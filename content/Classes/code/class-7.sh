#! /usr/bin/env bash

#BSUB -J bedtools.vignette
#BSUB -o vignette.%J.out
#BSUB -e vignette.%J.err

genes=/vol1/opt/data/encode/refGene.bed.gz
tss=tss.bed

# make annotations
zcat $genes | awk '$6 == "+"' \
    | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1}' > $tssbed
zcat $genes | awk '$6 == "-"' \
    | awk 'BEGIN {OFS="\t"} {print $1,$3,$3+1}' > $tssbed

bedSort $tssbed $tssbed
gzip $tssbed

signal=/vol1/opt/data/endcode/wgEncodeBroadHistoneHelas3Pol2bStdSig.bigWig

chromsize=/vol1/opt/data/hg19/hg19.chrom.sizes
slopbed=tss.slop.2000.bed

# make slop
bedtools slop -b 2000 -i $tssbed -g $chromsize > $slopbed

# make windows
windowbed=tss.slop.2000.5bp.windows.bed
bedtools makewindows -b $slopbed -w 5 -i srcwinnum \
    | sort -k1,1 -k2,2n \
    | tr "_" "\t" \
    > $windowbed

# map data
signalmap=signal.map.bg
bedtools map -a $windowbed \
    -b <(bigWigToBedGraph $signal stdout) \
    -c 4 -o mean -null 0 \
    > $signalmap 

# group data
sort -t$'\t' -k5,5n $signalmap \
    | bedtools groupby \
        -i - \
        -g 5 -c 6 -o sum \
        > output.tab 

