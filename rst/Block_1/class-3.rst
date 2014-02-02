Class 3 : The command-line (part 2)
===================================

Goals
-----
1. learn additional linux utilities (cut, sort, zless, uniq, wget)
2. understand how to combine tools with pipes (|)

wget
----
fetch a a file from the web with ``wget``::

    $ cd /opt/bio-workshop/data/
    $ wget http://ucd-bioworkshop.github.io/_downloads/states.tab

cut
---
The ``cut`` command allows you to extract certain columns of a file::

    # cut columns 1-4 and 7-10
    $ cut -f 1-4,7-10 /opt/bio-workshop/data/states.tab

    # cut columns 1-4
    $ cut -f 1,2,3,4

    # cut first column 1
    $ cut -f 1 /opt/bio-workshop/data/lamina.bed

    # output all columns after the 1st
    $ cut -f 2- /opt/bio-workshop/data/lamina.bed

uniq
----
The ``uniq`` command  allows you to get and count unique entries::

    # remove duplicate lines
    $ cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq

    # show duplicate lines
    $ cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq -d

    # count unique entries:
    $ cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq -c

.. important::

   ``uniq`` assumes that file is sorted by the column of interest.

   Use ``sort`` to sort the data before ``uniq``-ing it.

Redirection of output
---------------------
To send the output of a command (or a file) to another file, use ">"::

    $ cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq -c > output.txt
    $ head output.txt

To **append** the output of a command (or a file) to another file, use
">>"::

    $ echo "last line" >> output.txt
    $ tail output.txt

Compressed Files
----------------
The most common way to uncompress single files is ``gunzip``::

    $ gunzip /opt/bio-workshop/data/t_R1.fastq.gz

And re-zip the file with ``gzip``:: 

    $ gzip /opt/bio-workshop/data/t_R1.fastq

But if we just want to stream the uncompressed data without changing the
file::

    $ zless /opt/bio-workshop/data/t_R1.fastq.gz

Pipes
-----
We probably want to do something with the file as we uncompress it::

    $ zless /opt/bio-workshop/data/t_R1.fastq.gz | head

We already know the head command prints the first -n lines.

Try piping the output to some other commands (tail|echo|cowsay)


Sort
----
You will often want to ``sort`` your data.

Have a look at::

    $ man sort

The main flag is `-k` to indicate which column to sort on.

You will also sometimes use `-u` to get unique entries.

Sort Questions
--------------

How do you:
   1) sort by a particular column? (-k 4)
   2) sort as a number (-k4n)
   3) sort as a general number (1e-3 < 0.05) (-k4g)
   4) change the default delimiter (-t,
   5) sort by 2 columns (-k1,1 -k2,2n)
   6) sort in reverse as a number (-k1rn)
   7) get unique entries (-u)

If you know all these, you'll know 99% of what you'll use sort for.

Sort Example
------------
BED files have columns `chrom` [tab] `start` [tab] `end` [tab] ...

Sort by chrom, then by start (a lot of tools will require this)::

    $ sort -k1,1 -k2,2n /opt/bio-workshop/data/lamina.bed > /tmp/sorted.bed

This tells it to sort the chromosome [column 1] as a character and the
start [column 2] as a number.

Question:
+++++++++

    What happens if you omit the `n` ?

Sort Example (part 2)
---------------------
What if we want to sort by Income **descending** in the 3rd column?::

    $ sort -t$'\t' -k3,3rg /opt/bio-workshop/data/states.tab > /tmp/sorted.out
    $ head /tmp/sorted.out 

Sort Exercise
-------------
Print out the 10 states (1st column, contains spaces) with the highest
income (3rd column) from states.tab using ``sort`` and piping to ``cut``.

Or, use ``cut`` and pipe to ``sort`` to do the same.

Application 1
-------------
Use pipes (|) chained together to look see which transcription factor
binding sites are the most common in a set of putative sites from ENCODE.

  + data file available from http (wget)
  + compressed BED format (zless)
  + TF name in 4th column (cut)
  + count frequency (uniq -c) after sorting (sort)
  + sort resulting frequencies so most common are first (sort -rn)
  + show top 10 (head)

Application 2
-------------
Note that we are using the variable FILE for the long file name::

    # BED format file of transcription factor binding sites
    FILE=http://bit.ly/tfbs-x

    wget --quiet -O - $FILE \
        | zless \
        | head -n 7000 \
        | cut -f 4 \
        | sort \
        | uniq -c \
        | sort -k1,1rn \
        | head -n 10

.. comments::
    FILE=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV2.bed.gz

Let's go through this line by line ...

grep
----
Use ``grep`` to identify lines in a file that match a specified pattern.

To find any instance of *chr5* in the lamina.bed file::
   
    # grep [pattern] [filename]
    $ grep chr5 /opt/bio-workshop/data/lamina.bed | head

To find all lines that start with a number sign::

    # The caret (^) matches the beginning of the line
    # FYI dollar sign ($) matches the end
    $ grep '^#' /opt/bio-workshop/data/lamina.bed

To find any line that *does not* start with "chr"::

    # the -v flag inverts the match (grep "not" [pattern])
    $ grep -v '^chr' /opt/bio-workshop/data/lamina.bed

grep (2)
--------
Beware of using ``grep`` to find patterns that might be partial matches::

    # this will match chr1, chr10, chr11 etc.
    $ grep chr1 /opt/bio-workshop/data/lamina.bed | cut -f1 | uniq

Also beware of using ``grep`` to search for numbers::

    $ grep 100 /opt/bio-workshop/data/lamina.bed | head -n 20

If you're trying to find numeric values in a file, you should use ``awk``
instead.

In Class Exercises
------------------
::

  1. To learn about piping (|), use cowsay to:

     a. show your current working directory
     b. show the number of lines in /opt/bio-workshop/data/lamina.bed
     c. show the most recently modified file/dir in $HOME

  2. write a bash script that you can run to list only the 2
     most recently modified files in a given directory (using
     what you've learned in this class)
  3. make that script executable (use google to learn how to do this).

  4. With `head`, you can see the first line of a file with head -n1.
     How can you see all of a file *except* the first line. (use google)

  5. Without using your history, how few keystrokes can you use to run
     the following command (must work from any directory)?
        ls /opt/bio-workshop/data/lamina.bed

  6. How few keystrokes can you do 5. using your history?

