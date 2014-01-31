Week 2 / Class 1 : (bio)awk
===========================

Goals
-----

1. learn awk basics
2. learn bioawk basics

awk
---

http://en.wikipedia.org/wiki/AWK

*AWK is an interpreted programming language designed for text processing and typically used as a data extraction and reporting tool. It is a standard feature of most Unix-like operating systems.*

Named after authors *A* ho, *W* einberger, *K* ernighan

basic program structure
-----------------------

actions on each line that match a pattern.

    awk 'PATTERN { ACTIONS }'

columns available as *variables* $1, $2, ... $n



program structure example
-------------------------

extract out p-values (in 4th column from BED file) that are < 0.04

    $ awk '$4 < 0.04' all.pvalues.bed > some.pvalues.bed

limit to chr12 (&& is and):

    $ awk '$4 < 0.04 && $1 == "chr12"' all.pvalues.bed > chr12.pvalues.bed

awk continued
-------------

$0 contains the entire line.

multiple patterns:

    $ awk '($4 < 0.05) { print $0"\tsignificant"}($4 >= 0.05) { print $0"\tlame" }' input.bed > output.classified.bed
   

bioawk
------

Bioawk is a variant of awk that knows about common sequence formats. To count
the number of records in a fastq file:

    $ bioawk -c fastx 'END { print $NR }'

biowak (names)
--------------

You can access `name`, `seq`, `qual`, `comment`:

   $ biowak -c fastx '{ print $name, $seq, $qual}'


See the README at: https://github.com/lh3/bioawk
for more info

