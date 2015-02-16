
********************************
Class 23 : Exome-Seq - Alignment
********************************

Goals
=====

#. Not kill Tesla
#. Exome Sequencing
#. Alignment
#. QC Alignments

Exome Sequencing
================

+ captures coding regions of genome
+ deeper coverage and more samples at same cost 
+ goal is to call variants (SNPs/deletions/insertions)
+ higher concentration of interesting variants
  relative to non-coding genome

Project Design
==============
 
.. image:: ../_static/images/phenotype.png


Project Design II
=================

8 sequenced samples from single family

+ paired-end 100 base reads
+ 5 affected
+ 3 unaffected

Analysis
========

+ How likely are we to find something?

+ Naive probability that 5 individuals have a variant and 3 do not (without
  considering pedigree)?

+ What if we find 1 million SNPs total?

+ How can we filter?

Alignment
=========

+ bwa mem to align the reads. 
    
+ bwa mem very good aligner for long reads

http://arxiv.org/pdf/1303.3997v2.pdf

BWA mem
=======

.. image:: ../_static/images/bwa-mem-roc.png

Mark Duplicates
===============

Marking/removing duplicates is to remove the influence of PCR duplication on
variant calling.

+ Marking/Removing duplicates makes sense unless you have extremely deep coverage,
  especially on paired-end reads.

+ Uses position of both reads of a pair to determine duplication--if many pairs
  have the same start and sequence, all but 1 pair is marked as a duplicate.

+ Recommended by the variant caller that we will use.

+ To mark or not to mark? http://seqanswers.com/forums/showthread.php?t=6854

Read Groups
===========

Your alignments should have read-groups.

These tell downstream tools what library and sample your data came from

Used for calling duplicates (repeats must be from same library)

With read-groups, can have multiple samples in a single BAM.


Variant Calling
===============

At each site in the genome, given:

+ base calls
+ base qualities
+ mapping qualities
+ surrounding sequence
+ position in read
+ strand bias
+ all of the above for other samples

What is the genotype of each sample and the confidence in that genotype?

Variant Calling Software
========================

The most commonly used software for variant calling is GATK:
Genome Analysis Toolkit

It has a number of steps including base recalibration, indel realignment,
etc.

We will use `freebayes` which has been shown to perform as well/better
(and it is much easier to use): https://github.com/ekg/freebayes

http://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/

Freebayes
=========

https://github.com/ekg/freebayes

*FreeBayes is haplotype-based, in the sense that it calls variants based on the
literal sequences of reads aligned to a particular target, not their precise
alignment.*

*FreeBayes uses short-read alignments (BAM files) for any number of individuals
from a population and a reference genome to determine the most-likely combination
of genotypes for the population at each position in the reference. It reports
positions which it finds putatively polymorphic in variant call file (VCF) format.*


VCF Format
==========

File format for SNP/indel calls:
http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
::

    CHROM      1
    POS        2345544
    ID         rs1234..
    REF        C
    ALT        TA
    QUAL       55.2
    FILTER     PASS
    INFO       NS=3;DP=14;AF=0.5;DB
    FORMAT     GT:GQ:DP:HQ
    sample_1   0|0:48:1:51,51
    ...
    sample_n

":" delimited chunks from samples match the format in `FORMAT`


Pipeline
========

Walk through .sh scripts on *amc-tesla*

