***************************************************
Class 8 : Working in a cluster environment (part 2)
***************************************************

Goals
=====
#. Leverage cluster resources to solve bioinformatic problems

Overview
========
We will go over some common scenarios in bioinformatics where you can
leverage cluster resources to process many files simultaneously with
identical workflows.

Scenario 1
==========
Sequencing experiments typically generate raw sequencing data from many
different experimental conditions (i.e. treatments, controls, cell types,
time points, etc.) However, the initial steps of processing sequencing
data use common steps to processthe data into an interpretable form. A
typical workflow for assessing sequencing data is:

    #. Calculate summary statistics on the raw reads
    #. Align reads to a reference genome
    #. Calculate coverage from the alignment

.. todo::

    figure out how to:

        #. move the following code block to a file
        #. have it displayed in the rendered html
        #. provide a link to the file so that one can download (instead of
            copying and pasting)

    Maybe use this?: http://sphinx-doc.org/markup/inline.html#ref-role

Now assume we have 4 different FASTQ files from a sequencing experiment.

.. code-block:: bash

    #! /usr/bin/env bash

    #BSUB -J workflow[1-4]
    #BSUB -e %J.%I.err
    #BUSB -o %J.%I.out

    # this section defines the sample names and grabs the appropriate
    # sample name for the current process
    SAMPLENAMES=(SP1 SP2 SP3 SP4)

    # $LSB_JOBINDEX is set at runtime for the current process; in this
    # case you asked for 4 jobs, so it be a value between 1 and 4
    sample=${SAMPLENAMES[$(($LSB_JOBINDEX - 1))]}

    # set up file names
    fastq=$sample.fq.gz
    align=$sample.sam
    coverage=$sample.coverage

    # run the workflow
    fastqc $fastq
    bowtie2 -x test -U $fastq > $align

    # now generate counts of each position that was aligned to; fields 3
    # and 4 of the same file are chrom and pos.
    cut -f3,4 $align \
        | sort \
        | uniq -c \
        | awk 'BEGIN {OFS="\t"} {print $2,$3,$1}'
        > $summary

Run the above script and check its status immediately::

    $ bsub < run.sh
    $ bjobs

You should see 4 running jobs, each with its own index i.e. workflow[1],
workflow[2] etc.    

These will run for a bit. After they are done, you will see several new
files, corresponding to the output of the same analysis applied to all the
different samples.

Scenario 2
==========
