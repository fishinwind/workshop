********************
Class XXX : BEDtools
********************

:Class date: XXX 

Goals
=====

#. perform operations on BED files with BEDTools

BEDTools Overview
=================

BEDTools will be one of the tools with the largest return on investment. For
example, to extract out all genes that overlap a CpG island:

.. code-block:: bash

    $ bedtools intersect -a genes.bed -b cpg_island.bed \
                                     > genes-in-islands.bed

:ref:`intersect <bedtools:intersect>` is a bedtools tool. It follows a common pattern in bedtools
that the query file is specified after the ``-a`` flag and the *subject* file
after the ``-b`` flag

BEDTools Commands
=================

To see all available BEDTools commands, type

.. code-block:: bash

    $ bedtools

BEDTools required files
=======================
Many of the BEDTools functions use a "genome file" which is a list of
chromosomes and their sizes for a given chromsome build::

    # hg19 chrom sizes file example
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


