
*********************************
  Class 17 : 
*********************************

:Class date: 2014 Mar 7 Friday

Goals
=====

 #. 

 #. 

Coverage plots for ChIP and DNase I
===================================

You will want to calculate coverage plots that are appropriate for the
type of experiment. For example, in a ChIP experiment, you want to examine
the entire region covered by sequences.

.. code-block:: bash

    $ common_args="-ibam <aln.bam> -g <chrom.size> -bg"
    $ bedtools genomecov $common_args > coverage.bg

But in a DNase I mapping experiment, you only want to know the exact
position where DNase I cut the DNA, i.e. the 5' position. Use the ``-5``
flag in bedtools to only count 5' positions.

.. code-block:: bash

    $ common_args="-ibam <aln.bam> -g <chrom.size> -bg"
    $ bedtools genomecov $common_args -5 > coverage.5p.bg


Problem Set Questions
=====================

Anybody have questions on the problem set?

