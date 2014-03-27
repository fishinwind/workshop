
********************
Class 22 : RNA-Seq I
********************

Goals
=====

 #. Why RNA-Seq
 #. QC.
 #. Align Reads

RNA-Seq
=======

 + We use RNA-Seq primarily to measure gene expression
 + We can also find transcript use (transcript-switching/alternative-splicing)
 + Fusion Transcripts
 + Possible to call variants (SNPs) on RNA-Seq data

Library Prep
============

Jay ??

Example Data
============

We will look at 10 million reads from 4 samples:

 + 2 x controls
 + 2 x RTCB

Jay RTCB ??

75 base paired-end reads

Analysis Plan
=============

 #. QC reads with Fastqc
 #. Trim reads by quality with sickle
 #. Re-check QC on trimmed reads
 #. Align reads with Tophat
    #. give a GTF to help it to find known transcripts
    #. also allow novel junctions
    #. fusion transcripts

 #. Count reads by gene with featureCounts (from subread)
 #. Differential Expression statistics with DESeq2

QC Reads
========

We will use fastqc to check the quality of the reads and look for adapter
contamination

