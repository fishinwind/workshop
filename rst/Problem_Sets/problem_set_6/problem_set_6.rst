.. _problem-set-6:

*************
Problem Set 6
*************

:Due date: 2014 Mar 25 at 9 PM MST

Problem 1
=========

Write up a short description of the goals for your final project (**15
points**).

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

Problem 2
=========

Use bedtools, dplyr and ggplot2 to analyze some ENCODE data (**15
points**).

#. Load the peaks BED file and find the 10 factors that have the largest
   range in peak width. Inspect a geom_boxplot() or geom_viol in() to
   support your answer (also add individual points to the plot with
   geom_jitter()).

#. Repeat the above, but only examine peaks do not overlap with any other
   peak (use bedtools to find these).

Problem Set Submission
======================

Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak

