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

Problem 1
---------

Use the peak calls from to derive a binding motif for the CTCF
transcription factor. (**10 points**) You will need to:

  - create FASTA sequence from CTCF ChIP-seq peak calls from the hg19
    genome.
  - use meme to identify motifs from these FASTA sequences.
  - write a run script to drive this analysis.

Keep the meme report and submit it with your homework.

Problem 2
---------

Use BEDtools to intersect peaks calls from the clustered transcription factor
binding regions with the clustered DNase I peaks. (**20 points**)

 #. Identify transcription factor binding peaks that were never identified
    as DNase I hypersensitive sites. What transcription factor binding sites
    are represented in these peaks?

 #. Do the converse: identify DNase I hypersensitive sites that do not
    have corresponding transcription factor peak calls. What motifs are
    enriched in this set of hypersensitive sites?

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak

