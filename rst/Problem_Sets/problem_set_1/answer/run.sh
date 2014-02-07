#! /usr/bin/env bash
#
# answer to quiz 1
#

set -o nounset -o pipefail -o errexit -x

project=pset_1

# make project directories
dirnames="results data doc"
for dir in $dirnames; do
    mkdir -p "$project/$dir"
done

# get current date
date=$(date "+%Y-%m-%d")

