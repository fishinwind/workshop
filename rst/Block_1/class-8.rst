Class 8 : Working in a cluster environment
===================================================

We will go over a common scenario in bioinformatics in which you can
leverage cluster resources to process many files simultaneously with
identical workflows.

Sequencing experiments typically generate raw sequencing data from many
different experimental conditions (i.e. treatments, cell types, time
points, controls, etc.) However, the initial steps of processing
sequencing data rely on common steps to get the data into a manageable
form. A typical workflow for assessing sequencing data is:

    1. Calculate summary statistics on the raw reads
    2. Align reads to a reference genome
    3. Calculate coverage from the alignment

Now assume we have 4 different FASTQ files from a sequencing experiment.

.. code-block:: bash

    # /usr/bin/env bash

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

