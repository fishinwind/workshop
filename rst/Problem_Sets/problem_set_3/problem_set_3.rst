
.. _label-problem-set-3:

*************
Problem Set 3
*************

:Due date: 2014 Feb 18 at 5 PM MST

Overview
--------
For this quiz you will write programs in Python to analyze data. 

.. note::

    Continue to use the organization scheme that we learned about in
    :ref:`label-problem-set-1`. Part of our evaluation
    will include whether you are developing good organizational habits.

Create a ``run.sh`` file that executes the commands for each problem and
writes out each result in a dated directory.

Problem 1
---------
Use Python to perform a similar task as the previous problem set, where
you used ``awk``. Write a Python program to read the lamina.bed file and
report the following:

   #. What is the region with the largest start position (2nd column) **for
      each** chromosome in `lamina.bed`? (**5 points**)

   #. What is the region with the largest end position on chrY in
      lamina.bed? (**5 points**) Report as::

        chrom <tab> start <tab> end <tab> value <tab> region_length

Problem 2
---------
Use Python to read a file in FASTQ format. FASTQ records are comprised of
4 sections; in this case each section will be on a unique line::

    @cluster_2:UMI_ATTCCG             # record name
    TTTCCGGGGCACATAATCTTCAGCCGGGCGC   # DNA sequence
    +                                 # empty line; starts with '+'
    9C;=;=<9@4868>9:67AA<9>65<=>591   # phred-scaled quality scores

The FASTQ file is on amc-tesla::
    
    /vol1/opt/data/SP1.fq

Write a Python program to parse the FASTQ records and report the
following:

   #. Which of the first 10 sequence records has the largest number of 'C'
      residues in the sequence? Report the record name (**5 points**)
    
   #. XXX (**5 points**)

Problem 3
---------

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`label-problem-set-submission`).

.. raw:: pdf

    PageBreak
