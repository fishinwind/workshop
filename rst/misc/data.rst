*****************
Sample Data files
*****************

We will use several example data files throughout the class.

.. _bed-file:

BED format
==========
Data in BED format contains region information (e.g. single nucleotides or
megbase regions) in a simple format [#]_:

**Download a sample BED file:** :download:`lamina.bed <../../data/lamina.bed>`

.. [#] BED documentation 
       http://genome.ucsc.edu/FAQ/FAQformat.html#format1

.. _fastq-file:

FASTQ format
============
FASTQ format contains DNA sequence data with quality scores::

    @cluster_2:UMI_ATTCCG             # record name
    TTTCCGGGGCACATAATCTTCAGCCGGGCGC   # DNA sequence
    +                                 # empty line; starts with '+'
    9C;=;=<9@4868>9:67AA<9>65<=>591   # phred-scaled quality scores

**Download a sample FASTQ file:** :download:`SP1.fq <../../data/SP1.fq>`

