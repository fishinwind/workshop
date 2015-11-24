.. include:: /_static/substitutions.txt

**************************
Class 20 : Python Workout
**************************

:Class date: |c20-date|
:Last updated: |today|

Announcements
=============

**No problem set this week.** There will be 1 more, due Mon April 20.

**Keep practicing Python.** At this point, it's use-it-or-lose-it.

Final Project Exmaple
=====================

:Hypothesis: Cell-type-specific DNase I hypersensitive sites are less
             evolutionarily conserved than constitutive hypersensitve sites

#. Determine the similarity of DNase I hypersensitivity across all
   mapped human cell types.

#. Calculate mean conservation scores for the matrix of cell types.

#. Correlate patterns of DHS sensititvity and conservation.

In-class exercises
==================

BAM file
--------

1. Write a Python program using the `pysam`_ module to parse this BAM
   file: ``/vol1/opt/data/bams/2_8-bis.bam``.

   a. Report the number of reads that map to each chromosome

   b. Identify sequences that map to more than 1 location (i.e.,
      "multiply-mapping" reads). Report the sequence with the largest number
      of different locations.

.. _pysam: http://pysam.readthedocs.org/en/latest/
