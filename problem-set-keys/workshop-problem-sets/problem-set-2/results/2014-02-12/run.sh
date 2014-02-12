#! /usr/bin/env bash
#
# run script for quiz 2 key

# these are bash flags the print out variables that get set when you
# run the script.
# set -o nounset -o pipefail -o errexit -x

# You will need to change the '???' strings below.
#
# define the project variable here. this should be the full path to
# your project directory, i.e. the directory at the top of the
# results/data/doc directories.

project="$HOME/workshop-problem-sets/problem-set-2"

# fill in the date here
date="2014-02-12"

# these refer to the data file that you moved into place
data=$project/data/$date
datafile=$data/lamina.bed

# these refer to the place where you will write the results of the
# "analysis"
results=$project/results/$date
resultsfile=$results/result.tab

# if the directory doesn't exist, make it
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

# Problem 1.1
# What is the region with the largest start position (2nd column) on any
# chromosome in lamina.bed? 

echo "Answer 1.1" > results.tab

sort -k2,2nr $datafile \
    | head -n 1 >> results.tab
echo >> results.tab

# Problem 1.2
# What is the region with the largest end position on chrY in lamina.bed?
# Report this region in the format: chr12:1234-5678

echo "Answer 1.2" >> results.tab

awk '$1 == "chrY"' $datafile \
    | sort -k2,2nr | head -n 1 \
    | awk 'BEGIN {OFS=""} {print $1,":",$2,"-",$3}' \
    >> results.tab
echo >> results.tab

# Problem 2.1
# What is the longest region (end - start) in lamina.bed?

echo "Answer 2.1" >> results.tab

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $3-$2}' $datafile \
    | sort -k5,5nr \
    | head -n 1 \
    >> results.tab
echo >> results.tab

# Problem 2.2
# What is the longest region (end - start) in lamina.bed with a value (4th
# column) greater than 0.9 on chr13. Report the header (1st line) in
# lamina.bed as well as the region
echo "Answer 2.2" >> results.tab

# print header line
head -n 1 $datafile >> results.tab

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $3-$2}' $datafile \
    | awk '$1 == "chr13" && $4 > 0.9' \
    | sort -k5,5nr \
    | head -n 1 \
    >> results.tab

echo >> results.tab

# Problem 3.1
# What are the regions that overlap this interval in lamina.bed:
# chr12:5,000,000-6,000,000? Report regions so that they are ordered by
# descending value (4th column), and the columns are separated by commas
# rather than tabs

echo "Answer 3.1" >> results.tab

# print header line
awk '$1 == "chr12" && $2 >= 5e6 && $3 <= 6e6' $datafile \
    | sort -k4,4gr \
    | awk 'BEGIN {OFS=","} {print $1,$2,$3,$4}' \
    >> results.tab

echo >> results.tab

# Problem 3.2
# What is the average value (4th column) of those regions from lamina.bed
# that overlap the region (5 points): chr12:5,000,000-6,000,000?

echo "Answer 3.2" >> results.tab

# print header line
awk '$1 == "chr12" && $2 >= 5e6 && $3 <= 6e6' $datafile \
    | awk '{sum += $4} END {print sum / NR}' \
    >> results.tab

echo >> results.tab

