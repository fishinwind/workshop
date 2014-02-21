********************
Class XXX : BEDtools
********************

:Class date: XXX 

Goals
=====
#. Learn to perform operations on BED files with BEDTools

BEDTools required files
=======================
Many of the BEDTools functions use a "genome file" which is a list of
chromosomes and their sizes for a given chromsome build::

    # hg18 chrom sizes file example
    # <chrom> <tab> <size>

UCSC provides a tool to fetch these files:

.. code-block:: bash

    $ fetchChromSizes hg19 > hg19.chrom.sizes

BEDTools map()
==============
The :ref:`BEDTools map <bedtools:map>` function is useful for aggregating
data across intervals and performing math operations on that data:

.. code-block:: bash

    $ bedtools map -a lamina.bed -b peaks.bed


