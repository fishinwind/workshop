********************
Class 13 : BEDtools
********************

:Class date: Wed 2014 Feb 26 

Goals
=====

#. Introduce the BEDTools suite of tools
#. Understand why using BEDTools is needed.
#. Practice common operations on BED files with BEDTools

BEDTools Overview
=================

BEDTools will be one of the tools with the largest return on investment. For
example, to extract out all genes that overlap a CpG island:

.. code-block:: bash

    $ bedtools intersect -wa -a genes.bed -b cpg_island.bed \
                                     > genes-in-islands.bed

:ref:`intersect <bedtools:intersect>` is a bedtools tool. It follows a
common pattern in bedtools that the query file is specified after the
``-a`` flag and the *subject* file after the ``-b`` flag

BEDTools Utility
================

Finding all overlaps between a pair of BED files in python would look like:

.. code-block:: python

    for a in parse_bed('a.bed'):
        for b in parse_bed('b.bed'):
            if overlaps(a, b):
                # do stuff
If 'a.bed' has 10K entries and 'b.bed' has 100K entries, this would involved
checking for overlaps **1 billion times**. That will be slow.

BEDTools uses an indexing scheme so that it reduces the number of tests
dramatically.

.. note::
  
  See the original BEDTools paper for more information: http://bioinformatics.oxfordjournals.org/content/26/6/841.full

BEDTools Utility (2)
====================

 + Fast: faster than intersect code you will write
 + Terse: syntax is terse, but readable
 + Formats: handles BED, VCF and GFF formats
 + Special Cases: handles stranded-ness, 1-base overlaps, abutted intervals,
   etc.
   - (all the things that will likely be bugs in your code should you do this
     manually)





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


BEDTools required files
=======================
Many of the BEDTools functions use a "genome file" which is a list of
chromosomes and their sizes for a given chromsome build::

    # hg19 chrom sizes file example
    # <chrom> <tab> <size>

Download these into /opt/bio-workshop/data/ from:

    **hg19:** :download:`hg19.genome <../misc/data/hg19.genome>`
    **hg18:** :download:`hg18.genome <../misc/data/hg18.genome>`


BEDTools Documentation
======================

The BEDTools documentation is quite good and ever improving.

See the documentation for :ref:`intersect <bedtools:intersect>` with:

.. code-block:: bash

    $ bedtools intersect

The online HTML help is also good and includes pictures: 
 https://bedtools.readthedocs.org/en/latest/content/tools/intersect.html



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

#. What is the total signal from a bedGraph that overlaps mRNA exons?
   (map)

