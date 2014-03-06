#! /usr/bin/env bash

#BSUB -J jaccard
#BSUB -e %J.err 
#BSUB -o %J.out

filenames=/vol1/opt/data/wgEncodeUwDnase/*.gz
result=result.tab

python analysis.py $filenames > $result 
