
.. include:: /_static/substitutions.txt 

.. _problem-set-1:

*************
Problem Set 1
*************

:Due date: |pset1-due| 

Reading
-------
For this problem set, you'll need to read Bill Noble's paper on organizing
computational biology projects [#]_ and be able to put a set of files in
the correct places. You should adopt this scheme for all of your projects
in and out of class.

.. [#] A Quick Guide to Organizing Computational Biology Projects (2009)
        PLoS Comput. Biol. William S. Noble
        http://dx.plos.org/10.1371/journal.pcbi.1000424

Computational biology projects inevitably accrue a lot of files, and this
paper is a great way to organize all of this information.  As you read in
the paper, organization of projects is important for remembering what you
did, and reanalyzing data when changes are made.

This quiz will also test your ability to do simple tasks on the command
line. You will need to take the tutorial to learn the necessary tools
[#]_.

.. [#] The Command Line Crash Course
        http://cli.learncodethehardway.org/book/

You will use the command line tools discussed in the tutorial (e.g. mkdir,
cd, ls, mv) to create the directory structure, move files into place and
check whether everything looks ok.

Problem 1
---------
Make a directory structure as outlined in the paper.  The directories
should be nested under a common project directory, with directories for
data, results and documentation (doc). You should also create dated
folders with today's date so you know where to put the dated data and
reuslts. Finally, keep a list of the commands that you used to create the
directory (**10 points**).

.. note::
    
    Name your directories with YYYY-MM-DD format (year-month-date, e.g.
    2014-01-10). If you do this, the directories will sort chronologically
    using the directory name.
    
    **Avoid MM-DD-YY, they will not sort chronologically by name.**

Problem 2
---------
Download the the following data table: :download:`states.tab`.  Then move
the data file to the appropriate *dated* directory.

.. code-block:: bash

    $ wget http://hesselberthlab.github.io/workshop/_downloads/states.tab


You need to create a run.sh shell script that runs the following code, and
writes the output into a results directory with the current date (**10
points**).

You should copy the following code block into a file using `vim`. You will
need to change the '???' characters in the file to correspond to the path
you want to write the results in (hint: it should include today's date).

.. code-block:: bash
    :linenos:

    #! /usr/bin/env bash
    #
    # run script for quiz 1

    # these are bash flags the print out variables that get set when you
    # run the script.
    set -o nounset -o pipefail -o errexit -x

    # You will need to change the '???' strings below.
    # 
    # define the project variable here. this should be the full path to
    # your project directory, i.e. the directory at the top of the
    # results/data/doc directories.

    project=???

    # fill in the date here
    date=???

    # these refer to the data file that you moved into place
    data=$project/data/$date
    datafile=$data/states.tab

    # these refer to the place where you will write the results of the
    # "analysis"
    results=$project/results/$date
    resultsfile=$results/result.tab

    # if the directory doesn't exist, make it
    if [[ ! -d $results ]]; then
        mkdir -p $results
    fi

    # Note how we are using redirects here. The first ">" writes a file,
    # and overwrites existing data. The following ">>" append data to the
    # existing file

    echo "here is our starting data ..."
    cat $datafile > $resultsfile
    echo

    echo "here are the states sorted by population size ..."
    sort -k2n $datafile >> $resultsfile
    echo

    echo "here are the states with the highest number of murders ..."
    sort -k6n $datafile | head -n 10 >> $resultsfile

Then, save the above text in a run.sh script in your results directory. To run the file,
use:

.. code-block:: bash

    $ bash run.sh

If this ran correctly, you should see a new results.tab file in the
results directory you specified in the run.sh script. If you don't see the file, double check
the path you specified, and make sure you're looking in the right spot. If
it's in a different spot than you intended. remove the results file you
wrote, update the program and run it again.

Problem 3
---------
Finally you need to create a log of what you did in the root of the
results directory to summarize the key points of your analysis (**5
points**). For example::

    Captain's log, star date 2014-07-16
    -----------------------------------
    After examining the results.tab file, learned that Alaska has the highest
    income per person. Something must be wrong ...

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak
