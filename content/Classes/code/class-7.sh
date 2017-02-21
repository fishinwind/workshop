#! /usr/bin/env bash

genes="data-sets/bed/genes.hg19.bed.gz"
tssbed="tss.bed"

# make TSS annotations

gzcat $genes | awk '$6 == "+"' \
 | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$4}' > $tssbed
gzcat $genes | awk '$6 == "-"' \
 | awk 'BEGIN {OFS="\t"} {print $1,$3,$3+1,$4}' >> $tssbed

# extract chr22 intervals and sort bed output
awk '($1 == "chr22")' $tssbed \
  | bedtools sort -i - > tmp.bed
mv tmp.bed $tssbed

signal="data-sets/bedtools/ctcf.hela.chr22.bg.gz"

chromsize="data-sets/bedtools/hg19.genome"
slopbed="tss.slop.2000.bed"

# make slop
bedtools slop -b 2000 -i $tssbed -g $chromsize > $slopbed

# make windows
windowbed="tss.slop.2000.5bp.windows.bed"
bedtools makewindows -b $slopbed -w 5 -i srcwinnum \
    | sort -k1,1 -k2,2n \
    | tr "_" "\t" \
    > $windowbed

# map data
signalmap="signal.map.bg"
bedtools map -a $windowbed \
    -b $signal \
    -c 4 -o mean -null 0 \
    > $signalmap 

# group data
sort -t$'\t' -k5,5n $signalmap \
    | bedtools groupby \
        -i - \
        -g 5 -c 6 -o sum \
        > output.tab 

# inspect output.tab to see how to plot the data
