#! /usr/bin/env bash

# analysis script for quiz 1
# take a filename on input, read it and dump some data to stdout

filename=$1
if [[ -z $1 ]]; then
    echo "Usage: `basename $0` filename" >&2
    exit 1
fi


