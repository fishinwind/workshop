
.. include:: /_static/substitutions.txt 

.. _problem-set-3:

*************
Problem Set 3
*************

:Due date: |pset3-due| 

Overview
--------

For this quiz you will analyze ChIP-seq data.

Use the provided :ref:`encode-data` for these problems. All of the data
are available on the amc-tesla cluster at::

    /vol1/opt/data

Note that the FASTA sequence for hg19 is available in the same directory
``hg19.fa``. That is a big file, so do not copy it anywhere; just use it
out of that directory.

Problem 1
---------

Transcription factors bind specific DNA sequences. In this exercise, you
will determine the specific DNA sequence bound by an abundant human
transcription factor called CTCF.

Use the CTCF peak calls in the ENCODE data to derive a binding motif for
the CTCF transcription factor. (**10 points**) You will need to:

#. Select out the peaks on chr22 (otherwise meme will take a long time
   to run)

#. create FASTA sequence from the peak calls from the hg19 genome.

#. use MEME to identify motifs from these FASTA sequences. The motifs
   in the meme output will be ranked according to their significance
   (i.e. how often this motif would occur in random sequence).

**Report the top 5 high scoring motifs, and compare your most significant
motif with what is already known about CTCF binding sites.**

Problem 2
---------

Use BEDtools to intersect peaks calls from clustered transcription factor
binding sites (TFBS) with clustered DNase I peaks. (**20 points**)

- Identify transcription factor binding sites that do not overlap with
  DNase I hypersensitive sites.
    
  * What transcription factors are represented in these regions? Analyze
    the BED file output to get this answer.

- Do the converse: identify DNase I hypersensitive sites that do not
  have corresponding transcription factor peak calls.
    
  * What motifs are enriched in this set of hypersensitive sites on chr22?

  * How would you use existing tools to identify similar motifs, ideally
    ones that are previously associated with a transcription factor?

**Report the factors in the peaks and the top 5 high scoring motifs from
each meme analysis.**

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak

