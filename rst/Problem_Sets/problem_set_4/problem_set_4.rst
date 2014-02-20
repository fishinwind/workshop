
.. _problem-set-4:

*************
Problem Set 4
*************

:Due date: 2014 Feb 25 at 9 PM MST

Overview
--------
For this quiz you will write programs in Python to analyze data. 

.. note::

    Continue to use the organization scheme that we learned about in
    :ref:`problem-set-1`. Part of our evaluation
    will include whether you are developing good organizational habits.

Create a ``run.sh`` file that executes the commands for each problem and
writes out each result in a dated directory.

Use the :ref:`data-sets` for these promblems.

Problems 
--------
#.  Convert the code you wrote in :ref:`problem-set-3` to parse BED and FASTQ
    files to functions (**10 points**). Instead of reading records like
    this:

    ..code-block:: python
    
        for line in file(bedfilename):
            fields = line.strip()

    Write a function that returns records from a file using
    :py:obj:`yield`:

    .. code-block:: python

        def parse_bed(bedfilename):
            # parse records

        for record in parse_bed(bedfilename):
            # use the record

#.  Use the functions to create nested data structures built from records
    in a BED file.

    #. Load a BED file and create a dict() of lists() of (start, end)
       tuples. Find the largest and smallest starts for each chromosome
       (**10 points**)

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak
