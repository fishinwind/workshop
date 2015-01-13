
.. _problem-set-3:

*************
Problem Set 3
*************

:Due date: 2014 Feb 21 at 9 PM MST

Overview
--------
For this quiz you will write programs in Python to analyze data. 

.. note::

    Continue to use the organization scheme that we learned about in
    :ref:`problem-set-1`. Part of our evaluation
    will include whether you are developing good organizational habits.

Create a ``run.sh`` file that executes the commands for each problem and
writes out each result in a dated directory.


Problem 1 (BED files)
---------------------
Use Python to perform a similar tasks as Problem Set 2, where you used
``awk``. Write a Python program to read the lamina.bed file and report the
following:

- What is the region with the largest start position (2nd column) **for
  each** chromosome in `lamina.bed`? (**5 points**)

- What is the region with the largest end position on chrY in
  lamina.bed? (**5 points**) Report as::

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
  residues in the sequence? Report its record name (**5 points**).
    
- For each of the first 10 records, Covert each character in the
  quality score to a number, and sum the numbers. Use :py:func:`ord`
  to convert characters to numbers (**5 points**).

- Use the Python :py:class:`~collections.Counter` to count unique
  sequences in the Fastq file. Report the top ten most abundant
  sequences in the Fastq file (**5 points**).

- Report the revese complement of each of the first 10 sequences (**5
  points**)

.. note::

    The reverse complement of ``5'-AGCTCGTA-3''`` is ``5'-TACGAGCT-3'``

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).


Question 1
==========

Here is one way to answer question 1. Note that we could also take
advantage of the fact that the file is sorted so the largest start
will always be the last line for each chromosome. Instead, here we
store the largest start we've seen so far, along with the value
in a dictionary keyed by chromosome.

.. code-block:: python
    :linenos:

    import sys

    bedfile = sys.argv[1]

    largest_by_chrom = {}

    for line in open(bedfile):
        if line.startswith('#'): continue
        chrom, start, end, value = line.rstrip('\r\n').split('\t')
        start = int(start)

        # hint for question 2. if we wanted to just get chrY, we could add...
        # if chrom != "chrY": continue

        # we haven't seen it before, so it has to be the largest
        if not chrom in largest_by_chrom:
            largest_by_chrom[chrom] = (chrom, start, end, value)
        else:
            # we have to check if the current start is greater
            # than the one we've stored.
            largest_chrom, largest_start, largest_end, largest_value \
                                                   = largest_by_chrom[chrom]
            if start > largest_start:
                # if it is, then we store a new region
                largest_by_chrom[chrom] = (chrom, start, end, value)

    # see what it looks like:
    # print largest_by_chrom
    for chrom in largest_by_chrom:
        chrom, start, end, value = largest_by_chrom[chrom]
        print chrom, start, end, value

.. raw:: pdf

    PageBreak
