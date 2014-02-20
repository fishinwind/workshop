
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

Use the provided :ref:`data-sets` for these promblems.

Problem 1
---------
Start with the code you wrote in :ref:`problem-set-3` to parse BED and
FASTQ files. Convert the parsing logic to functions (**10 points**).

Instead of reading records like this ... :

.. code-block:: python

    for line in file(bedfilename):
        fields = line.strip()

... use this skeleton code to write a function that returns records from a
file using :ref:`yield() <python:yield>`:

.. code-block:: python
   :linenos:

    def parse_bed(bedfilename):
        ''' parse records and return each record '''

        # here, you need to write code to split the fields, assign them to
        # `chrom`, `start`, `end`, and `value`
        # be sure to coerce the values to ints and floats as needed
        
        result = {'chrom':chrom, 'start':start, 'end':end, 'value':value}
        yield result

    for record in parse_bed(bedfilename):
        # use the records to:
        #
        # 1. calculate the distance between the start of the current record
        # and the previous record.
        #
        # 2. calculate the bases covered by the intervals for each
        # chromosome. note you may have to define a structure outside of
        # this loop to keep track of that information.

Problem 2
---------
Modify the following skeleton code to create nested data structures built
from records in a BED file, using the function you created above.

Load a BED file and create a dict() of lists() of (start, end)
tuples. Find the largest and smallest starts for each chromosome using
min() and max().
(**10 points**)

.. code-block:: python
   :linenos:

    from collections import defaultdict

    # specify the bedfilename 
    bedfilename = 'XXX'
    struct = defaultdict(list)

    for record in parse_bed(bedfilename):
       
        chrom = record['chrom']
        # tuple of coords 
        coords = (start, end)

        # add the coords to the growing list. replace `whichmeth` with the
        # appropriate method call
        struct[chrom].whichmeth(coords)

    for chrom in struct:
        # 1. use max() and min() in this loop to determine biggest start
        # values?
        #
        # 2. how do you change the max() and min() calls to look at the `end`
        # value instead of the `start`? (RTM)

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak
