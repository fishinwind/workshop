*************************************
 Class 5 : Working with genomic data 
*************************************

:Class date: Tues 10 Feb 2015

Goals
=====

#. What is the ENCODE project?
 
#. What kinds of data did the ENCODE project produce? 
 
#. Where can I find these data on the Internet? 

#. What can I do with these data sets?
 
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


.. raw:: pdf

    PageBreak
