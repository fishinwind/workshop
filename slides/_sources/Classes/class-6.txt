*****************************************
Class 6 : Genomic Data Vignette: ChIP-seq
*****************************************

:Class date: Thurs 12 Feb 2015

Goals
=====

#. Learn workflows for analyzing ChIP-seq data 

#. Derive DNA sequence motifs from ChIP-seq peaks 

ENCODE
======
 
The Human Genome Project was finished, giving us a list of human genes and their 
locations. Unfortunately, we still had no idea how they were regulated. If only 
there was an `ENCyclopedia Of Dna Elements 
<http://www.sciencemag.org.hsl-ezproxy.ucdenver.edu/content/306/5696/636.full>`_…

Advantages: massive amounts of information on key cell lines, reproducible 
experiments, public data access, technology development.

ENCODE Project Cell Lines
=========================

Tier 1: GM12878 (EBV-transformed lymphoblast), K562 (CML lymphoblast), H1-hESC

Tier 2: HeLa-S3 (cervical cancer), HepG2 (liver carcinoma), HUVEC (umbilical vein)

Tier 2.5: SKNSH (neuroblastoma), IMR90 (lung fibroblast), A549 (lung carcinoma), 
MCF7 (breast carcinoma), LHCN (myoblast), CD14+, CD20+
 
`link <http://genome.ucsc.edu/ENCODE/cellTypes.html>`_ (this page also has very useful
links to cell culture protocols)

Experiments
===========

#. ChIP-seq: Histone marks, transcription factors

#. Chromatin structure: DNaseI-seq, FAIRE, 5C/Hi-C

#. RNA expression: mRNA-seq, GENCODE gene predictions

#. Data Integration: Segway / ChromHMM integration of functional data

Common File Formats
===================

+ FASTQ: Raw sequencing data. `[link] <http://maq.sourceforge.net/fastq.shtml>`

+ SAM/BAM: Aligned sequence data `[link] <http://samtools.github.io/hts-specs/SAMv1.pdf>`

+ Bed/bigBed: List of genomic regions `[link] <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`

+ Bedgraph/Wig/bigWig: Continuous signal `[link] <http://genome.ucsc.edu/goldenPath/help/bedgraph.html>` 

Many other formats are described on this `page <http://genome.ucsc.edu/FAQ/FAQformat.html>`_

References
==========

Completion of the entire project, and a ton of papers: 
`Nature <http://www.nature.com/nature/journal/v489/n7414/index.html>`_, 
`Genome Research <http://genome.cshlp.org/content/22/9.toc>`_, 
`Genome Biology <http://genomebiology.com/content/13/9>`_, 

How to Access ENCODE Data
=========================

The `ENCODE project page <https://www.encodeproject.org/>`_ is the portal
to all of the ENCODE data.


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

(We'll talk more in depth about ChIP-Seq workflows in the future,
but for now, just a brief introdcution to a few commands you can use to
work on ChIP-Seq data (and pset3)).


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

