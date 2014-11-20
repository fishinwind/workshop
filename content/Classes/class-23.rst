
********************
Class 23 : RNA-Seq I
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

Experimental Design
===================

.. important::
   
   Please do not do RNA-Seq on fewer than 6 samples (3 vs 3).

   Even in a pilot study.

See: http://bioinformatics.oxfordjournals.org/content/30/3/301

Titled: **RNA-seq differential expression studies: more sequence or more
replication?**

Answer: **more replication**

Your bioinformatician will thank you.

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
#. Differential Transcript Use and fusion transcripts (if time)

QC Reads
========

We will use fastqc to check the quality of the reads and look for adapter
contamination

    $ fastqc -o ~/public_html/fqc-$sample $data/${sample}_R2.fastq.gz

Quality Trim Reads
==================

We will trim low-quality bases from reads using sickle:
https://github.com/najoshi/sickle

Lenient trimming has been shown to improve mapping rates in
RNA-Seq. See, e.g. http://genomebio.org/is-trimming-is-beneficial-in-rna-seq/

QC Quality Trimmed Reads
========================

Compare before/after. Makes most difference in bad datasets.

Align with Tophat2
==================

Many RNA-Seq (splice-aware) aligners:

+ Tophat/Tophat2
+ GSNAP
+ STAR
+ RUM
+ MapSplice
+ etc.

We will use tophat2.

IGenomes
=========

If you work with a model organism, you can get bowtie (tophat2) genome
indexes and feature annotation in normalized format from Illumina.

http://tophat.cbcb.umd.edu/igenomes.shtml

We will use the *C. elegans* data downloaded to `amc-tesla:~brentp/data/ce/`

+ Using this will save you a lot of trouble
+ Tophat2 uses known transcripts and attempts to align to those as well as to
  novel transcripts

Spliced Alignment
=================

From tophat paper:

.. image:: http://bioinformatics.oxfordjournals.org/content/25/9/1105/F1.large.jpg

Spliced Alignment
=================

From tophat2 paper:

.. image:: http://genomebiology.com/content/figures/gb-2013-14-4-r36-1-l.jpg

Tophat2 Invocation
==================

.. code-block:: bash

    tophat2 -o $out/results/$sample $reference $fq_1 $fq2 \
        --fusion-search -p 6 --transcriptome_index $TINDEX \
        --GTF $ANNOTATION_GTF

Output will be in accepted_hits.bam

Just Do It
==========

On amc-tesla

.. code-block:: bash

    mkdir -p ~/class-23/src/
    cp ~brentp/class-23/run.sh ~/class-23/
    cp ~brentp/class-23/src/clean-counts.py ~/class-23/src/
    cd ~/class-23/src/
    module load sickle
    module load subread
    module load bowtie2

