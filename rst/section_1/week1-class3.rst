Week 1 / Class 3 : The command-line++
=====================================

Goals
-----

1. learn some additional linux utilities (cut, sort, zless, uniq)
2. understand pipes (|)

cut
---

cut allows you to extract certain columns of a file.


.. code-block:: bash

    # cut columns 1-4 and 7-10
    cut -f 1-4,7-10 some.file.txt

    # cut columns 1-4
    cut -f 1,2,3,4

    # cut first column 1
    cut -f 1 /opt/bio-workshop/data/lamina.bed

    # output all columns after the 1st
    cut -f 2- /opt/bio-workshop/data/lamina.bed

uniq
----

uniq allows you to get and count unique entries.

 + remove duplicate lines
 + show duplicate lines
 + count unique entries:


.. code-block:: bash

    cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq -c

.. note::

    **cut** Assumes that file is sorted by the column of interest.

(re)direction
-------------

to send the output of a command (or a file) to another file, use ">"

.. code-block:: bash

    cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq -c > output.txt
    head output.txt

to **append** the output of a command (or a file) to another file, use ">>"

.. code-block:: bash

    echo "last line" >> output.txt
    tail output.txt


Compressed Files
----------------

Most common way to compress single files is `gzip`

    $ ls /opt/bio-workshop/data/\*.gz

We can (g)unzip that file as:

    $ gunzip /opt/bio-workshop/data/t_R1.fastq.gz

And re-zip is as:

    $ gzip /opt/bio-workshop/data/t_R1.fastq

But if we just want to stream the uncompressed data without changing the file:

    $ zless /opt/bio-workshop/data/t_R1.fastq.gz

Pipes
-----

We probably want to do something with the file as we uncompress it:

    $ zless /opt/bio-workshop/data/some.gz | head

We already know the head command prints the first -n lines.

Try piping the output to some other commands (tail|echo|cowsay)


Sort
----

You will often want to `sort` your data.

Have a look at

    $ man sort

The main flag is `-k` to indicate which column to sort on.

Sort Questions
--------------

How do you:
   1) sort by a particular column? (-k 4)
   2) sort as a number (-k4n)
   3) sort as a general number (1e-3 < 0.05) (-k4g)
   4) change the default delimiter (-t,
   5) sort by 2 columns (-k1,1 -k2,2n)
   6) sort in reverse (-k1r)

If you know all these, you'll know 99% of what you'll use sort for.

Sort Example
------------

BED files have columns `chrom` [tab] `start` [tab] `end` [tab] ...

Sort by chrom, then by start (a lot of tools will require this):

    $ sort -k1,1 -k2,2n some.bed > some.sorted.bed

This tells it to sort the chromosome [1] as a character and the
start as a number.

Question:
+++++++++

    What happens if you omit the `n` ?

Sort Example (2)
----------------

What if we want to sort by p-value **descending** in the 4th column?

    $ sort -k4,4rg some.pvals.txt > some.pvals.sorted.txt


Sort Question
-------------

Compress `some.pvals.txt` with gzip. Then zless that and
pipe the result to sort by p-value and show only the rows
with the 10 lowest p-values.

Application (1)
===============

Let's use pipes (|) chained together to look see which
transcription factor binding sites are the most common
in a set of putative sites from ENCODE.

  + data file available from http (wget)
  + compressed BED format (zless)
  + TF name in 4th column (cut)
  + count frequency (uniq -c) after sorting (sort)
  + sort resulting frequencies so most common are first (sort -rn)
  + show top 10 (head)

Application (2)
===============

.. code-block:: bash

    FILE=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV2.bed.gz

    wget --quiet -O - $FILE \
        | zless \
        | head -n 7000 \
        | cut -f 4 \
        | sort \
        | uniq -c \
        | sort -k1,1rn \
        | head -n 10

Let's go through this line by line...


grep
----

We use **grep** to find stuff.


In Class Exercises
------------------

Place the answers to these in the bash script:


    1. write a bash script that you can run to list only the 2 most recently
       modified files in a given directory (using what you've learned in this class)
    2. make that script executable (use google to learn how to do this).

    3. With `head`, you can see the first line of a file with head -n1.
       How can you see all of a file *except* the first line.

    4. Without using your history, how few keystrokes can you use to run the following command (must work from any directory)?

        ls /opt/bio-workshop/data/lamina.bed

    5. How few keystrokes can you do #4 using your history?

    6. To learn about piping (|), use cowsay to:

       a. show your current working directory
       b. tell  you the number of lines in /opt/bio-workshop/data/lamina.bed
       c. tell you the most recently modified file (or directory) in $HOME
