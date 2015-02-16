#! /usr/bin/env bash

# run script for Problem Set 3 keys

project="$HOME/devel/bio-workshop/problem-set-keys/problem-set-3"
bedfile="$project/data/2014-02-12/lamina.bed"
fastqfile="$project/data/2014-02-17/SP1.fq"

result="$project/results/2014-02-17"
results="$project/results/2014-02-17/result.tab"

if [ ! -d $result ]; then
    mkdir -p $result
fi

# run out of current directory
python problem-1-bedfile.py $bedfile >> $results
python problem-2-fastqfile.py $fastqfile >> $results

