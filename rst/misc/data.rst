
.. _data-sets:

*****************
Sample Data files
*****************

We will use several example data files throughout the class.

.. _bed-file:

BED format
==========
Data in BED format contains region information (e.g. single nucleotides or
megbase regions) in a simple format [#]_:

**Download a sample BED file:** :download:`lamina.bed <data/lamina.bed>`

.. [#] BED documentation 
       http://genome.ucsc.edu/FAQ/FAQformat.html#format1

.. _fastq-file:

FASTQ format
============
FASTQ format contains DNA sequence data with quality scores::

    @cluster_2:UMI_ATTCCG             # record name; starts with '@'
    TTTCCGGGGCACATAATCTTCAGCCGGGCGC   # DNA sequence
    +                                 # empty line; starts with '+'
    9C;=;=<9@4868>9:67AA<9>65<=>591   # phred-scaled quality scores

**Download a sample FASTQ file:** :download:`SP1.fq <data/SP1.fq>`

.. _fasta-file:

FASTA format
============
FASTA format just contains DNA sequence data; no quality scores::

    >cluster_2:UMI_ATTCCG             # record name; starts with '>'
    TTTCCGGGGCACATAATCTTCAGCCGGGCGC   # DNA sequence

**Download a sample FASTA file:** :download:`sample.fa <data/sample.fa>`

.. _encode-data:

ENCODE data
===========
All encode data are available at
https://genome.ucsc.edu/ENCODE/downloads.html

For Problem Set 5, you will need these files on the amc-tesla cluster,
available in::

    /vol1/opt/data

.. list-table::
    :header-rows: 1
    
    * - Experiment
      - Target
      - Cell line
      - Replicate
      - File Type
      - File name
    * - ChIP-seq
      - Histone H3 Lysine 4 trimethyl (H3K4me3)
      - Hela-S3
      - 1
      - FASTQ
      - wgEncodeBroadHistoneHelas3H3k4me3StdRawDataRep1.fastq.gz
    * - ChIP-seq
      - CTCF
      - Hela-S3
      - 1
      - narrowPeak
      - wgEncodeUwTfbsHelas3CtcfStdPkRep1.narrowPeak.gz
    * - Merged TFBS ChIP-seq
      - all
      - all
      - n/a
      - BED
      - wgEncodeRegTfbsClusteredV3.bed.gz
    * - Merged DNase I hypersensitive sites
      - all
      - all
      - n/a
      - BED
      - wgEncodeRegDnaseClusteredV2.bed.gz

