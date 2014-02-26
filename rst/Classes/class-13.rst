********************
Class 13 : BEDtools
********************

:Class date: Wed 2014 Feb 26 

.. code-block:: bash

    mkdir ~/learn-bedtools/
    cd ~/learn-bedtools
    cp ~/src/bedtools2/test/intersect/a.bed .
    cp ~/src/bedtools2/test/intersect/b.bed .
    curl -O https://ucd-bioworkshop.github.io/_downloads/cpg.bed.gz
    curl -O https://ucd-bioworkshop.github.io/_downloads/genes.hg19.bed.gz


or ask if you're not using the VM.

Goals
=====

#. Introduce the BEDTools suite of tools
#. Understand why using BEDTools is needed.
#. Practice common operations on BED files with BEDTools

BEDTools Overview
=================

BEDTools will be one of the tools with the best return on investment. For
example, to extract out **all genes that overlap a CpG island**:

.. code-block:: bash

    $ bedtools intersect -wa -a genes.bed -b cpg_island.bed \
                                     > genes-in-islands.bed

:ref:`intersect <bedtools:intersect>` is a bedtools tool. It follows a
common pattern in bedtools that the query file is specified after the
``-a`` flag and the *subject* file after the ``-b`` flag

BEDTools Utility
================

Finding all overlaps between a pair of BED files naively in python would look like:

.. code-block:: python

    for a in parse_bed('a.bed'):
        for b in parse_bed('b.bed'):
            if overlaps(a, b):
                # do stuff

If *'a.bed'* has 10K entries and *'b.bed'* has 100K entries, this would involved
checking for overlaps **1 billion times**. That will be slow.

BEDTools uses an indexing scheme so that it reduces the number of tests
dramatically.

.. note::
  
  See the original BEDTools paper for more information:
  http://bioinformatics.oxfordjournals.org/content/26/6/841.full

BEDTools Utility (2)
====================

 + Fast: faster than intersect code you will write
 + Terse: syntax is terse, but readable
 + Formats: handles BED, VCF and GFF formats
 + Special Cases: handles stranded-ness, 1-base overlaps, abutted intervals,
   etc.
  - (all the things that will likely be bugs in your code should you do this manually)


BEDTools Commands
=================

To see all available BEDTools commands, type

.. code-block:: bash

    $ bedtools

The most commonly used BEDtools are:

    + :ref:`intersect <bedtools:intersect>`
    + :ref:`genomecov <bedtools:genomecov>`
    + :ref:`closest <bedtools:closest>`
    + :ref:`map <bedtools:map>`


BEDTools Documentation
======================

The BEDTools documentation is quite good and ever improving.

See the documentation for :ref:`intersect <bedtools:intersect>` with:

.. code-block:: bash

    $ bedtools intersect

The online HTML help is also good and includes pictures: 
 https://bedtools.readthedocs.org/en/latest/content/tools/intersect.html


BEDTools intersect
==================
Have a browser window open to :ref:`BEDTools intersect documentation <bedtools:intersect>`.
It will likely be the BEDTools function that you use the most. It has a lot of
options.

.. image:: http://bedtools.readthedocs.org/en/latest/_images/intersect-glyph.png

"-v" means (like grep) include all intervals from `-a` that do not overlap
intervals in `-b`

Example Files
=============

.. code-block:: bash

    $ cat a.bed 
    chr1    10  20  a1  1   +
    chr1    100 200 a2  2   -

    $ cat b.bed 
    chr1    20  30  b1  1   +
    chr1    90  101 b2  2   -
    chr1    100 110 b3  3   +
    chr1    200 210 b4  4   +

What will happen if you intersect those files?
For example, the *a.bed* region `chr1:100-200` overlaps::

    chr1:90-101 
    chr1:100-101

from *b.bed*

intersect
=========

intersect with default arguments means **extract chunks of `-a` that overlap
regions in `-b`**

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed
    chr1    100 101 a2  2   -
    chr1    100 110 a2  2   -

Here is the original interval from *a.bed*::

    chr1	100	200	a2	2	-

And the overlapping intervals from *b.bed*::

    chr1	90	101	b2	2	-
    chr1	100	110	b3	3	+

intersect -wa
=============

Often, we want the *entire interval from -a if it overlaps any interval in -b*

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -wa
    chr1    100 200 a2  2   -
    chr1    100 200 a2  2   -

We can get that uniquely with (-u)



intersect -wo
=============

We can see which intervals in *-b* are associated with *-a*

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -wo
    chr1	100	200	a2	2	-	chr1	90	101	b2	2	-	1
    chr1	100	200	a2	2	-	chr1	100	110	b3	3	+	10

intersect exercise
==================

What happens if you reverse the arguments? E.g. instead of::

  -a a.bed -b b.bed

use::

   -b a.bed -a b.bed

Try that with no extra flags, with -u, -wa, -wu.

How does it compare to the original?

intersect -c
============

We can count overlaps for each interval in *-a* with those in *-b* with

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -c
    chr1	10	20	a1	1	+	0
    chr1	100	200	a2	2	-	2

This is our original `a.bed` with an **additional column indicating number of
overlaps** with `b.bed`


intersect -v
============

Extract intervals in `a.bed` that do not overlap any interval in `b.bed`

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -v
    chr1	10	20	a1	1	+

Extract intervals in `b.bed` that do not overlap any interval in `a.bed`
.. code-block:: bash

    $ bedtools intersect -a b.bed -b a.bed -v
    chr1	20	30	b1	1	+
    chr1	200	210	b4	4	+


Exercises
=========

#. zless :download:`cpg.bed.gz <../misc/data/cpg.bed.gz>` and :download:`genes.hg19.bed.gz <../misc/data/genes.hg19.bed.gz>`
#. Extract the fragment of CpG Islands that touch any gene.
#. Extract CpG's that do not touch any gene
#. Extract (uniquely) all of each CpG Island that touches any gene.
#. Extract CpG's that are completely contained within a gene (look at the help
   for a flag to indicate that you want the fraction of overlap to be 1 (for 100 %).
#. Report genes that overlap any CpG island.
#. Report genes that overlap more than 1 CpG Island (use -c and awk).

.. important::

as you are figuring these out, make sure to pipe the output to less or head

BEDTools map()
==============
The :ref:`BEDTools map <bedtools:map>` function is useful for aggregating
data across intervals and performing math operations on that data:

.. code-block:: bash

    $ bedtools map -a lamina.bed -b peaks.bed

BEDTools example problems 
=========================
#. What are all the peaks (i.e. BED regions) in this file that overlap with
   another set of peaks? (intersect)

#. Which of these features overlap with exons / introns / transcription
   start sites / 3' UTRs (in another BED file)? (intersect)


