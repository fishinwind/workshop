
********************
Class 14 : ChIP-seq
********************

:Class date: 2014 Feb 28 Friday

Goals
=====

#. Learn the workflow for analyzing ChIP-seq data

Chromatin Immunoprecipitation Overview
======================================

Chromatin Immunoprecipitation is used to determine where a protein of
interest binds on a chromatin template [Park_Chipseq]_.

.. [Park_Chipseq] http://www.nature.com/nrg/journal/v10/n10/full/nrg2641.html

.. image:: ../_static/images/chip-workflow.png

.. nextslide::

.. image:: ../_static/images/chip-data.png

ChIP-seq analysis workflow
==========================

A general workflow for visualizing ChIP-seq data (and many other types of
data) is:

.. list-table::
    :widths: 40 40
    :header-rows: 1

    * - Operation
      - File formats
    * - Align reads to reference genome
      - ``FASTQ ~~> BAM``
    * - Generate coverage plots
      - ``BAM ~~> bedGraph``
    * - Call peaks 
      - ``BAM ~~> BED``
    * - Make binary files for UCSC display
      - ``bedGraph ~~> bigWig``, ``BED ~~> bigBed``
    * - Identify motifs
      - ``BED ~~> FASTA ~~> TXT / HTML``

ChIP-seq data
=============

Look at some human ChIP-seq data [#]_.

.. [#] Genome Browser Session
       http://goo.gl/WfJxcM

Short read alignment
====================

There are several short read alignment packages available. We will use
bowtie2 [#]_ because it is easy to use and relatively fast.

.. code-block:: bash

    # minimal bowtie2 command
    $ bowtie2 -x <index> -U <read1.fq.gz> [options] > output.sam

.. [#] Bowtie2
       http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

.. nextslide::

Normally, you would also pipe the SAM-format alignment through the samtools
suite to discard unaligned reads, sort the alignment and store it in
binary format (bam).

SAM format: http://samtools.sourceforge.net/SAMv1.pdf

.. code-block:: bash

    # sort sam output with samtools and write bam output
    $ bowtie2 -x <index> -U <read1.fq.gz> [options] \
        | samtools view -ShuF4 - \ 
        | samtools sort -o - aln.temp -m 8G \
        > aln.bam

Generate and visualize coverage plots
=====================================

Once alignment is complete, you can create coverage plots from your aligned
data, so that you can visualize your data.

.. tip::

    You will need a "chromsizes" flle for many BEDtools commands. This file
    contains the sizes of each chromosome in an assembly. UCSC provides a
    tool to retrieve this information:

    .. code-block:: bash

        # retrieve chrom sizes for the hg19 assembly and write them to a file
        # inspect this file so you know what it looks like
        $ fetchChromSizes hg19 > hg19.chrom.sizes

Coverage plots with BEDtools
----------------------------

To generate coverage plots, we will use BEDtools. Here, we'll
use the :ref:`genomecov <bedtools:genomecov>` tool.

.. code-block:: bash

    # -bg : write ouptput in bedGraph format
    $ bedtools genomecov -i <aln.bam> -g <chrom.sizes> -bg > coverage.bg

This command writes a bedGraph format file called ``coverage.bg``. Use
``less`` to examine this file.

.. nextslide::
    :increment:


Words to live by: **If you make a BED file, sort the BED file**

Many strange things can happen if you use unsorted BED files for
analysis. Once you create a BED file, sort it with one of these:

.. code-block:: bash

   # using UNIX sort
   $ sort -k1,1 -k2,2 unsorted.bed > sorted.bed

   # UCSC tool. same filename twice, overwrites original file
   $ bedSort file.bed file.bed

   # or you can use bedtools; writes additional file
   $ bedtools sort -i - < unsorted.bed > sorted.bed


Coverage plots split by strand
------------------------------

For some experiments, you will analyze the data relative to each strand of
the reference genome. For example, RNA is transcribed in single-stranded
form and derives from one or the other strand.

During alignment, reads from an RNA-based experiment will map to either
the positive ('+' or ``pos``) or negative ('-' or ``neg``) strand. You can
generate signal plots for ``pos`` and ``neg`` strands separately with
``bedtools``:

.. code-block:: bash

    $ common_args="-ibam <aln.bam> -g <chrom.size> -bg"
    $ bedtools genomecov $common_args -strand + > coverage.pos.bg
    $ bedtools genomecov $common_args -strand - > coverage.neg.bg

You would then create bigWigs for each of these display the stranded data
in the Genome Browser.


Peak calling
============

There are several available software packages for identifying regions
enriched in your IP experiment (i.e. peaks). We will use macs2 here.

.. code-block:: bash

    # minimal macs2 command 
    $ macs2 callpeak --treatment <aln.bam> --name <exp.name> [options]

Identify sequence motifs in enriched regions
============================================

You can use meme [#]_ to identify over-represented motifs in groups of
sequences (e.g. sequences covered by ChIP peaks).

Use the :ref:`bedtools getfasta <bedtools:getfasta>` command to fetch
fasta sequences.

Note: meme looks at both strands of a DNA sequence by default.

.. [#] MEME 
       http://meme.nbcr.net/meme/

.. code-block:: bash

    $ bedtools getfasta -fi <ref.fa> -bed <peaks.bed> -fo peaks.fa
    $ meme -nmotifs 5 -minw 6 -maxw 20 -dna <peaks.fa>

Putting it all together
=======================
Here is a script that combines the above in a single workflow:

.. literalinclude:: code/chipseq.sh
   :language: bash
   :linenos:

