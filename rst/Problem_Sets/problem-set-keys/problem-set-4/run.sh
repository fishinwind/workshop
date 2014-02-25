#! /usr/bin/env bash

# run script for Problem Set 4 key

project="$HOME/devel/bio-workshop/problem-set-keys/problem-set-4"
bedfile="$project/data/2014-02-2/lamina.bed"
u

result="$project/results/2014-02-24"
results="$project/results/2014-02-24/result.tab"

if [ ! -d $result ]; then
    mkdir -p $result
fi

# run out of current directory
python problem-4.py $bedfile >> $results

