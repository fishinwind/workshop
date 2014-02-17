#! /usr/bin/env bash

project="$HOME/devel/bio-workshop/problem-set-keys"
bedfile="$project/problem-set-3/data/2014-02-12/lamina.bed"
fastqfile="$project/problem-set-3/data/2014-02-17/"

result="$project/results/2014-02-17result.tab"

if [ ! -d $result ]; then
    mkdir -p $result
fi

python problem-1-1.py $bedfile >> $result
python problem-1-2.py $bedfile >> $result

