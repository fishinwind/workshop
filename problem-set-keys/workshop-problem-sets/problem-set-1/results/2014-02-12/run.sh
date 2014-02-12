#! /usr/bin/env bash
#
# run script for quiz 1

# these are bash flags the print out variables that get set when you
# run the script.
set -o nounset -o pipefail -o errexit -x

# You will need to change the '???' strings below.
#
# define the project variable here. this should be the full path to
# your project directory, i.e. the directory at the top of the
# results/data/doc directories.

project="$HOME/workshop-problem-sets/problem-set-1"

# fill in the date here
date="2014-02-12"

# these refer to the data file that you moved into place
data=$project/data/$date
datafile=$data/states.tab

# these refer to the place where you will write the results of the
# "analysis"
results=$project/results/$date
resultsfile=$results/result.tab

# if the directory doesn't exist, make it
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

# Note how we are using redirects here. The first ">" writes a file,
# and overwrites existing data. The following ">>" append data to the
# existing file

echo "here is our starting data ..." > $resultsfile
cat $datafile >> $resultsfile
echo >> $resultsfile

echo "here are the states sorted by population size ..." \
    >> $resultsfile
sort -k2n $datafile >> $resultsfile
echo >> $resultsfile

echo "here are the states with the highest number of murders ..." >> \
    $resultsfile

grep 'Name' $datafile | sort -k6gr $datafile | \
    head -n 10 >> $resultsfile
echo >> $resultsfile

