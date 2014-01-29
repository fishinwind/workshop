Class 3: The command-line++
===========================

Goals
-----

1. understand how to apply some common linux utilities to files
2. understand pipes (|)


Compressed Files
----------------

Most common way to compress single files is `gzip`

    $ ls /opt/bio-workshop/data/\*.gz

We can (g)unzip that file as:

    $ gunzip /opt/bio-workshop/data/some.gz

And re-zip is as:

    $ gzip /opt/bio-workshop/data/some

But if we just want to stream the uncompressed data without changing the file:

    $ zless /opt/bio-workshop/data/some.gz

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


Sort Questions
--------------

How do you:
   1) sort by a particular column?
   2) sort as a number
   3) sort as a general number (1e-3 < 0.05)
   4) change the default delimiter
   5) sort by 2 columns
   6) sort in reverse

If you know all these, you'll know 99% of what you'll use sort for.

Sort Example
------------

BED files have columns `chrom` [tab] `start` [tab] `end` [tab] ...

Sort by chrom, then by start (a lot of tools will require this):

    $ sort -k1,1 -k2,2n some.bed > some.sorted.bed

This tells it to sort the chromosome [1] as a character and the
start as a number.

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

grep
----

We use **grep** to find stuff.

