
.. include:: /_static/substitutions.txt 

.. _problem-set-5:

*************
Problem Set 5
*************

:Due date: |pset5-due| 

Problem 1
=========

Use bedtools, dplyr and ggplot2 to analyze some ENCODE data (**30
points**). You can use Python to filter the data prior to analysis. 

#. Load the peaks BED file and find the 10 factors that have the largest
   range in peak width. Inspect a ``geom_boxplot()`` or ``geom_violin()``
   to support your answer (also add individual points to the plot with
   ``geom_jitter()``).

   You should calcuate a `min.width` and `max.width` for each transcription
   factor in the BED file. Then the `range` is the difference between the
   `max.width` and `min.width`.

   Hint: you can ``summarize()`` to calculate the min.width and max.width,
   and then ``mutate()`` the summary to add a `range` column.

#. Repeat the above, but only examine peaks do not overlap with any other
   peak (e.g. use bedtools to find these).

Problem 2
=========

Make up your own question that uses ggplot2 and BEDTools to analyze ENCODE
data (**20 points**).

This description should include:

#. a biological question motivating the analysis
#. how you will address the question

Example
-------

I will examine ENCODE ChIA-PET and ChIP-seq data to determine what histone
marks are enriched in high confidence CTCF-mediated chromatin loops::

    ----CTCF-----------------------------------------------CTCF----
         |                                                  
         -------------------loop region----------------------
             H3K4me3 peak?                 H3K4me3 peak?

Workflow
--------

#. download ENCODE ChIA-PET/CTCF and histone ChIP data for K562 cell line

#. Determine which histone marks are present in the regions defined by
   CTCF Chia-PET boundaries (as a control, determine which are outside)
   use ``bedtools intersect`` and other bedtools

#. Generate summary statistics for which histone marks are enriched in
   CTCF loops. use ``dplyr`` and ``ggplot2`` for this.


Problem Set Submission
======================

Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak

