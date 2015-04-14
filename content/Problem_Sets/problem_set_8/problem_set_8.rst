.. include:: /_static/substitutions.txt 

.. _problem-set-8:

***************
 Problem Set 8
***************

:Due date: |pset8-due| 

Overview
--------
For this quiz you will write programs in Python to analyze data. 

Problem 1 (functions)
---------------------

Write a Python function to parse records from a BED file. Then use the
function to parse the BED file. Add an option to ignore specific
chromosomes. The function should look something like:

.. code:: python

    def parse_bed(filename, ignore=[]):
        ''' Documentation '''        
        ...
        # ignore specified chroms here
        ...

        yield(fields)

    # ignore chr1 in this code
    for region in parse_bed(filename, ...):
        ...

1. Document your function using the Google [Style]_ (**10 points**).

.. [Style] http://sphinxcontrib-napoleon.readthedocs.org/en/latest/index.html

2. Parse the ``lamina.bed`` file, find the following for each of the
   odd-numbered and even-numbered chromosomes:

   a. Largest region on each (**5 points**)
   b. Mean region size on each using ``numpy`` (**5 points**).


Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak
