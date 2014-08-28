***********************************
Class 3 : The command-line (part 2)
***********************************

:Class date: Friday 31 January 2014

Goals
=====
1. learn additional linux utilities (cut, sort, zless, uniq, wget)
2. understand how to combine tools with pipes (|)

wget
====
fetch a a file from the web with ``wget``

.. code-block:: bash

    $ cd /opt/bio-workshop/data/
    $ wget http://ucd-bioworkshop.github.io/_downloads/states.tab

cut
===
The ``cut`` command allows you to extract certain columns of a file

.. code-block:: bash

    # cut columns 1-4 and 7-10
    $ cut -f 1-4,7-10 /opt/bio-workshop/data/states.tab

    # cut columns 1-4
    $ cut -f 1,2,3,4

    # cut first column 1
    $ cut -f 1 /opt/bio-workshop/data/lamina.bed

    # output all columns after the 1st
    $ cut -f 2- /opt/bio-workshop/data/lamina.bed

uniq
====
The ``uniq`` command  allows you to get and count unique entries

.. code-block:: bash

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
=====================
To send the output of a command (or a file) to another file, use ">"

.. code-block:: bash

    $ cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq -c > output.txt
    $ head output.txt

To **append** the output of a command (or a file) to another file, use
">>"

.. code-block:: bash

    $ echo "last line" >> output.txt
    $ tail output.txt

Compressed Files
================
The most common way to uncompress single files is ``gunzip``

.. code-block:: bash

    $ gunzip /opt/bio-workshop/data/t_R1.fastq.gz

And re-zip the file with ``gzip``

.. code-block:: bash 

    $ gzip /opt/bio-workshop/data/t_R1.fastq

But if we just want to stream the uncompressed data without changing the
file

.. code-block:: bash

    $ zless /opt/bio-workshop/data/t_R1.fastq.gz

Pipes
=====
We probably want to do something with the file as we uncompress it

.. code-block:: bash

    $ zless /opt/bio-workshop/data/t_R1.fastq.gz | head

We already know the head command prints the first -n lines.

Try piping the output to some other commands (tail|echo|cowsay)


Sort
====
You will often want to ``sort`` your data.

Have a look at

.. code-block:: bash

    $ man sort

The main flag is `-k` to indicate which column to sort on.

You will also sometimes use `-u` to get unique entries.

Sort Questions
==============
How do you:

   #. sort by a particular column? (-k 4)
   #. sort as a number (-k4n)
   #. sort as a general number (1e-3 < 0.05) (-k4g)
   #. change the default delimiter (-t,
   #. sort by 2 columns (-k1,1 -k2,2n)
   #. sort in reverse as a number (-k1rn)
   #. get unique entries (-u)

If you know all these, you'll know 99% of what you'll use sort for.

Sort Example
============
BED files have columns `chrom` [tab] `start` [tab] `end` [tab] ...

Sort by chrom, then by start (a lot of tools will require this)

.. code-block:: bash

    $ sort -k1,1 -k2,2n /opt/bio-workshop/data/lamina.bed > /tmp/sorted.bed

This tells it to sort the chromosome [column 1] as a character and the
start [column 2] as a number.

Question:
+++++++++

    What happens if you omit the `n` ?

Sort Example (part 2)
=====================
What if we want to sort by Income **descending** in the 3rd column?

.. code-block:: bash

    $ sort -t$'\t' -k3,3rg /opt/bio-workshop/data/states.tab > /tmp/sorted.out
    $ head /tmp/sorted.out 

Sort Exercise
=============
Print out the 10 states (1st column, contains spaces) with the highest
income (3rd column) from states.tab using ``sort`` and piping to ``cut``.

Or, use ``cut`` and pipe to ``sort`` to do the same.

Application 1
=============
Use pipes (|) chained together to look see which transcription factor
binding sites are the most common in a set of putative sites from ENCODE.

  + data file available from http (wget)
  + compressed BED format (zless)
  + TF name in 4th column (cut)
  + count frequency (uniq -c) after sorting (sort)
  + sort resulting frequencies so most common are first (sort -rn)
  + show top 10 (head)

Application 2
=============
Note that we are using the variable FILE for the long file name

.. code-block:: bash

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

.. FILE=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV2.bed.gz

Let's go through this line by line ...

.. _class-3-exercises:

In Class Exercises - Class 3
============================
.. code-block:: bash

    1. To learn about piping (|), use cowsay to:

         a. show your current working directory
         b. show the number of lines in /opt/bio-workshop/data/lamina.bed
         c. show the most recently modified file/dir in $HOME

    2. write a bash script that you can run to list only the 2 most
    recently modified files in a given directory (using what you've
    learned in this class)

    3. make that script executable (use google to learn how to do this).

    4. With `head`, you can see the first line of a file with head -n1.
    How can you see all of a file *except* the first line. (use google)

    5. Without using your history, how few keystrokes can you use to run
    the following command (must work from any directory)?

.. code-block:: bash

        $ ls /opt/bio-workshop/data/lamina.bed

    6. How few keystrokes can you do 5. using your history?

.. raw:: pdf

    PageBreak
