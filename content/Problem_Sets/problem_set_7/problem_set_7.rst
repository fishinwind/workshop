.. include:: /_static/substitutions.txt 

.. _problem-set-7:

***************
 Problem Set 7
***************

:Due date: |pset7-due| 

Overview
--------
For this quiz you will write programs in Python to analyze data. 

Problem 1 (BED files)
---------------------
Use Python to perform a similar tasks as Problem Set 2, where you used
``awk``. Write a Python program to read the lamina.bed file and report the
following:

- What is the region with the largest start position (2nd column) **for
  each** chromosome in `lamina.bed`? (**10 points**)

- What is the region with the largest end position on chrY in
  lamina.bed? (**10 points**) Report as::

   chrom <tab> start <tab> end <tab> value <tab> region_length

Download the bed file here: :ref:`bed-file`

.. note::

    See end of this page for the solution to question 1.


Problem 2 (FASTQ files)
-----------------------
Use Python to read a file in FASTQ format. FASTQ records are comprised of
4 sections; in this case each section will be on a unique line::

    @cluster_2:UMI_ATTCCG             # record name
    TTTCCGGGGCACATAATCTTCAGCCGGGCGC   # DNA sequence
    +                                 # empty line; starts with '+'
    9C;=;=<9@4868>9:67AA<9>65<=>591   # phred-scaled quality scores

Download the fastq file here: :ref:`fastq-file`

Write a Python program to parse the FASTQ records and report the
following:

- Which of the first 10 sequence records has the largest number of 'C'
  residues in the sequence? Report its record name (**10 points**).
    
- For each of the first 10 records, Covert each character in the
  quality score to a number, and sum the numbers. Use :py:func:`ord`
  to convert characters to numbers (**10 points**).

- Report the revese complement of each of the first 10 sequences (**10
  points**). You will have to define a ``reverse_complement()``, or find
  one from another package.

.. note::

    The reverse complement of ``5'-AGCTCGTA-3''`` is ``5'-TACGAGCT-3'``

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak
