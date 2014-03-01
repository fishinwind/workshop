.. _problem-set-5:

*************
Problem Set 5
*************

:Due date: 2014 Mar 4 at 9 PM MST

Overview
--------

For this quiz you will analyze ChIP-seq data.

Use the provided :ref:`encode-data` for these promblems. ALl of the data
are available on the amc-tesla cluster at::

    /vol1/opt/data

Note that the FASTA sequence for hg19 is available in the same directory
``hg19.fa``. That is a big file, so do not copy it anywhere; just use it
out of that directory.

Problem 1
---------

Transcription factors (usually) bind specific DNA sequences. In this
exercise, you will determine the specific DNA sequence bound by an
abundant human transcription factor called CTCF.

Use the CTCF peak calls in the ENCODE data to derive a binding motif for
the CTCF transcription factor. (**10 points**) You will need to:

  #. create FASTA sequence from CTCF ChIP-seq peak calls from the hg19
     genome.

  #. Select out the peaks on chr22 (otherwise meme will take a long time
     to run)

  #. use MEME to identify motifs from these FASTA sequences. The motifs
     in the meme output will be ranked according to their significance
     (i.e. how often this motif would occur in random sequence).

**Report the top 5 high scoring motifs, and compare your most significant
motif with what is already known about CTCF binding sites.**

Problem 2
---------

Use BEDtools to intersect peaks calls from the clustered transcription factor
binding regions with the clustered DNase I peaks. (**20 points**)

 #. Identify transcription factor binding peaks that do not overlap with
    DNase I hypersensitive sites.
    
    - What transcription factors are represented in these peaks? Analyze
      the BED file output to get this answer.

    - What motifs are found in this set of peaks? How would you use
      existing tools to identify similar motifs, ideally ones that are
      previously associated with a transcription factor?

 #. Do the converse: identify DNase I hypersensitive sites that do not
    have corresponding transcription factor peak calls.
    
    - What motifs are enriched in this set of hypersensitive sites?

**Report the factors in the peaks and the top 5 high scoring motifs from
each meme analysis.**

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak

