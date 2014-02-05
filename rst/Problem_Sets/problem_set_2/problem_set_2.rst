
.. _label-problem-set-2:

*************
Problem Set 2
*************

:Due date: 2014 Feb 18 at 5 PM MST

Overview
--------
For this quiz you will use the tools we learned about in the last several
classes, focusing on manipulating text files with Linux command line
tools.

.. note::

    Continue to use the organization scheme that we learned about in
    **Problem Set 1** (:ref:`label-problem-set-1`). Part of our evaluation
    will include whether you are developing good organizational habits.

In the solution for all of the problems below, print out the region(s) as
they appear in the BED file, with additional columns at the end. e.g.::

    chr12 <tab> 1234 <tab> 5678 <tab> 0.93 <tab> my-extra-info

These should be separated by tabs (a ``\t`` character), **not spaces**,
unless otherwise indicated.

Combine ``awk`` with the other utilities you have learned. Create a
``run.sh`` file that executes the commands for each problem and writes out
each result in a dated directory.

Problem 1
---------

    #. What is the region with the largest start position (2nd column) on any
       chromosome in `lamina.bed`? (**5 points**)

    #. What is the region with the largest end position on chrY in
       lamina.bed? Report this region in the format: ``chr12:1234-5678``"
       (**5 points**)

Problem 2
---------

    #. What is the longest region (end - start) in `lamina.bed`? (**5 points**)
       Report as::

        chrom <tab> start <tab> end <tab> value <tab> region_length

    #. What is the longest region (end - start) in `lamina.bed` with a value
       (4th column) greater than 0.9 on chr13. Report the header (1st line) in
       lamina.bed as well as the region (**5 points**).

Problem 3
---------

    #. What are the regions that overlap this interval in `lamina.bed`:
       ``chr12:5,000,000-6,000,000``? Report regions so that they are ordered
       by descending value (4th column), and the columns are separated by commas
       rather than tabs (**5 points**).

    #. What is the average value (4th column) of those regions from `lamina.bed`
       that overlap the region (**5 points**): ``chr12:5,000,000-6,000,000?``

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`label-problem-set-submission`).
