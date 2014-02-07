***************************************************
Class 6 : Working in a cluster environment (part 2)
***************************************************

:Class date: Friday 7 February 2014

Goals
=====
#. Leverage cluster resources to solve bioinformatic problems
#. Exercises! (``grep`` review)
#. Practice typing for a bit.
#. Learn to switch between terminal windows with <Alt>-<Tab>

Note
====
We would like to quickly move toward having your homework be executable.
In other words, we should be able to:

.. code-block:: bash

    $ bash run.homework-1.sh

and have it execute all of the questions, writing results into files into
the right places as necessary. Getting into this habit lets you save
everything you did, and also re-run things when you need to make changes.

Overview
========
We will go over a common scenario in bioinformatics where you can
leverage cluster resources to process many files simultaneously with
identical workflows.

Exercise
========
Sequencing experiments typically generate raw sequencing data from many
different experimental conditions (i.e. treatments, controls, cell types,
time points, etc.) However, the initial steps of processing sequencing
data use common steps to process the data into an interpretable form. A
typical workflow for assessing sequencing data might begin with:

    #. Align reads to a reference genome
    #. Calculate coverage from the alignment

.. todo::

    figure out how to:

        #. move the following code block to a file
        #. have it displayed in the rendered html
        #. provide a link to the file so that one can download (instead of
            copying and pasting)

    Maybe use this?: http://sphinx-doc.org/markup/inline.html#ref-role

Exercise (2)
============
Login to amc-tesla.

.. code-block:: bash

    # the -X flag starts an X11 connection 
    $ ssh -X username@amc-tesla.ucdenver.pvt

We have 4 different FASTQ files from a sequencing experiment. They
are located in::

    /vol1/opt/data

Exercise (3)
============

.. note::

    Type this out. Do not copy / paste.

.. literalinclude:: code/class-6-exercise.sh
    :language: bash

Exercise (4)
============
Run the above script and check its status immediately:

.. code-block:: bash

    $ bsub < run.sh
    $ bjobs

You should see 4 running jobs, each with its own index i.e. workflow[1],
workflow[2] etc.    

These will run for a bit. After they are done, you will see several new
files, corresponding to the output of the same analysis applied to all the
different samples.

Questions
=========

 #. Check the ``.err`` files from the run. What information do they
    contain? What does this tell you about your starting sequences?

 #. Find out how you would modify the ``bowtie2`` command to write out the
    unaligned reads into a new file. Re-run the analysis to report those
    reads.

 #. Modify the ``awk`` command in the script to print out a valid BED4
    format::
    
    chrom <tab> start <tab> end <tab> count
    
 #. Find out how many unique UMI sequences are associted with each
    chromosomal coordinate (note: not as easy).

More exercises
==============

 #. use ``grep`` to identify lines in lamina.bed where the second field
    (start) begins with ``100``.

 #. use ``grep`` to identify lines in lamina.bed where the third field
    (end) ends with 99 .

 #. use ``grep`` with its ``-w`` flag to count the number of 'chr1'
    records in lamina.bed.

 #. use ``grep`` to count how many fastq records are in the
    /vol1/opt/data/t_R1.fastq.gz file (fastq records begin with an
    '@' symbol)

 #. login to amc-tesla. use ``grep`` to count the number of fastq records
    in /vol1/opt/data/SP1.fq.gz

