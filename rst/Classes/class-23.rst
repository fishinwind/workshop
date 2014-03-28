
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

Libraries were prepared by a standard RNA-seq protocol including:

 #. Oligo-dT purification of polyadenylated mRNA
 #. Reverse transcription with random nonamers
 #. Second strand sythesis to create double stranded DNA (dsDNA)
 #. Shearing of dsDNA to 300-600 bp
 #. End polishing, adaptor ligation, PCR with indexed primers
 #. Quantify libraries, Mix, sequence.

Example Data
============

We will look at ~10 million 75 base paired-end reads from 4 RNA-seq
libraries prepared from these `C. elegans` strains:

.. list-table::
    :header-rows: 1

 * - Strain
   - Treatment
 * - RtcB wild-type
   - No treatment
 * - RtcB wild-type
   - Tunicamycin (induces UPR)
 * - RtcB mutant
   - No treatment
 * - RtcB mutant
   - Tunicamycin (induces UPR)

.. nextslide::
    :increment:

RtcB is an RNA ligase. RtcB wild-type and mutant strains were treated with
tunicamycin, which inhibits protein glycosylation in the endoplasmic
reticulum and triggers the "unfolded protein reponse"::

    ER stress ~~>
        Ire1 multimerization ~~>
            Ire1 ribonuclease domain cuts
            xbp-1 transcript to release intron ~~>
                RtcB ligates xbp-1 exons together ~~>
                    Xbp-1 protein is made, which activates UPR-related genes
                    in the nucleus.

The expectation is that RtcB mutants will be unable to ligate specific RNA
substrates together, thus the mRNA library in the mutant might contain
altered levels of RNA fragments relative to the wild-type. In particular,
xbp-1 [#]_ fragments should accumulate in the mutant.

.. [#] http://www.wormbase.org/species/c_elegans/gene/WBGene00006959?query=xbp-1#0-9e-3

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

