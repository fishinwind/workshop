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
directories for data, results and documentation (doc). You should also
create dated folders with today's date so you know where to put the dated
data and reuslts (10 points).

Problem 2. Download the the data table from Canvas (data.tab). You need to
write a run.sh shell script that runs the following code, and
writes the output into a results directory with the current date (10 points).

You should copy the following code block into a file using gedit. You will
need to change the '???' characters in the file to correspond to the path
you want to write the results in (hint: it should include today's date).

.. code-block:: bash

    #! /usr/bin/env bash
    #
    # run script for quiz 1
    # you will need to change the '???' positions in the following:
    # 
    # define the project variable here. this should be the full path to
    # your project directory.
    project=???

    # fill in the date here
    date=???

    results=$project/$date/results

    mkdir -p $results

    # echo the data file to the results.tab
    # count up the number of cars with 6 cylinders
    echo "# cylinder counts 
    cut -f 3 | sort | uniq -c > results.tab

    # count how many of each type of car there are
    # echo "car type counts"
    cut -f 1 | cut -f1 -d' ' | sort | uniq -c >> results.tab

Then, save that file as run.sh in your results directory. To run the file,
use:

.. code-block:: bash

    $ bash run.sh

If this ran correctly, you should see a new results.tab file in the
results directory you specified in the run.sh script. If you don't see the file, double check
the path you specified, and make sure you're looking in the right spot. If
it's in a different spot than you intended. remove the results file you
wrote, update the program and run it again.

Problem 3. Finally you need to create a log of what you did in the root of the
results directory to summarize the key points of your analysis (5 points).

