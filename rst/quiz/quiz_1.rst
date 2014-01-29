Week 1 - Quiz 1
===============

Reading
-------
This quiz will test your ability to do simple tasks on the command line.
For reference, the tutorial is at: http://cli.learncodethehardway.org/book/

Reading
-------
For the first quiz, you'll need to read this paper and be able to put
a set of files in the correct places. We highly recommend adopting this
scheme for all of your projects in and out of the class.

Computational biology projects inevitably accrue a lot of files, and 
This paper is a great way to get started on organizing all of this
information:

1. Organization for computational biology projects
    `Link <http://dx.plos.org/10.1371/journal.pcbi.1000424>`_

As you read in the paper, organization of projects is important for
remembering what you did, and reanalyzing data when changes are made.

You will use the command line tools discussed in the tutorial (e.g. mkdir,
cd, ls, mv) to create the directory structure, move files into place and
check whether everything looks ok.

Problem 1. make a directory structure as outlined in the paper.
The directories should be nested under a common project directory, with
directories for data, results and documentation (doc). (10 points).

Problem 2. There are 2 files in this directory: a data table (data.tab) and a
script that reads the data and generates a small summary. You need to
write a run.sh shell script that executes the analysis.sh script and
writes data into a results directory with the current date (10 points).

For example, in your run.sh script you will need a line like this that
executes the program and writes out the result:

.. code-block:: bash

   bash analysis.sh data.tab > result.tab 

Problem 3. Finally you need to create a log of what you did in the root of the
results directory to summarize the key points of your analysis (5 points).
