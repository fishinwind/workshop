
.. include:: /_static/substitutions.txt 

.. _problem-set-5:

*************
Problem Set 5
*************

:Due date: |pset5-due| 

Overview
========

For the following problems, generate a Rmarkdown document that explains
the workflow and generates the appropriate plots.

Problem 1
=========

Use BEDTools, dplyr and ggplot2 to analyze some ENCODE data (**25
points**). 

#. Load the peaks BED file (`peaks.chr22.bed.gz`) and find the 10 factors
   that have the largest range in peak width. Inspect a ``geom_boxplot()`` or
   ``geom_violin()`` to support your answer (also add individual points to
   the plot with ``geom_jitter()``).

   You should calcuate a `min.width` and `max.width` for each transcription
   factor in the BED file. Then the `range` is the difference between the
   `max.width` and `min.width`.

   Hint: you can ``summarize()`` to calculate the min.width and max.width,
   and then ``mutate()`` the summary to add a `range` column.

#. Repeat the above, but only examine peaks do not overlap with any other
   peak (hint: use the ``-c`` argument to ``bedtools intersect`` to find these).

Problem 2
=========

Use ggplot2 and BEDTools to analyze the following ENCODE data (**25
points**). **Only look at data on chr22 for the following problem**.

Examine ENCODE CTCF and histone ChIP-seq data to determine what histone
marks are enriched near CTCF binding sites. 

ENCODE CTCF peak calls from K562 are in::

    /vol1/opt/data/wgEncodeUwTfbsK562CtcfStdPkRep1.narrowPeak.gz

#. Download ENCODE histone (H3K4me3 and H3K27ac) ChIP data for
   the K562 cell line.

#. Determine how histone marks are distributed around CTCF binding sites.
   This is very similar to the problem you did last week, i.e. use ``bedtools
   intersect`` and other bedtools

#. Generate summary statistics (e.g., how many histone peaks are near CTCF
   peaks, etc) using ``dplyr`` and ``ggplot2``.

Problem Set Submission
======================

Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak

