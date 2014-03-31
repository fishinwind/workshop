
*********************
Class 24 : RNA-Seq II
*********************

Goals
=====

 #. Check Alignments
 #. Sort Alignments
 #. Count reads
 #. Differential Expression


Alignment QC
============

Since it was a long-running job, we should check the logs:

.. code-block:: bash

    ls -lhtr logs/ | tail -n 4

And *less one of the .err* logs. The end should look something like this::

    [2014-03-28 16:38:53] Run complete: 01:26:06 elapsed
    + exit

The `.out` file should contain the text "Successfully completed."

Tophat Alignment Metrics
========================

Remember we had 4 samples. We can check the alignment metrics
from all of them::

    less results/CNTRL_002_CGATGT_L002.tophat2/align_summary.txt

We can also see the % of paired reads::

    grep "concordant" results/*/align_summary.txt

Samtools
========

The alignment files are in results/\*/accepted_hits.bam we can
count the number of paired reads

.. code-block:: bash

     samtools view -cf2 results/CNTRL_002_CGATGT_L002.tophat2/accepted_hits.bam

Samtools Sort
=============

With paired-end reads, we should sort the alignments by name so
that reads from the same pair are adjacent in the BAM.

.. code-block:: bash

    samtools sort -n some.bam some.sort


FastQC
======

We can run fastqc on the alignments. *see run.sh*


Count Reads
===========

subread/featureCounts: 

 + uses same GTF as tophat to counts number of reads per gene
 + we will adjust resulting file format with python script
 + send to DESeq for differential expression

Differential Expression
=======================

From `featureCounts` we get a matrix of differential expression. We can 
send this to DESeq2 for differential expression.

See `run.sh`
