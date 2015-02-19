
.. include:: /_static/substitutions.txt 

.. _problem-set-4:

*************
Problem Set 4
*************

:Due date: |pset4-due| 

Overview
--------

For this quiz you will use BEDTools and command line tools to make some
simple summary files and will use ggplot2 to make some plots.

Problem 1
---------

#. Examine Pol II distribution in a 5 kb window around transcription
   termination sites (TTSs) at a resolution of 50 bp?  Make a plot of the
   result and save it as a PDF. (**15 points**)

#. Pick one other factor you might expect to find at a TTS and analyze its
   distribution. You can find bigWig files of ChIP-seq experiemnts done for a
   variety of factors at UCSC [#]_ and the ENCODE [#]_ website. Make a plot of the
   result and save it as a PDF. (**15 points**)

.. [#] ENCODE Project at UCSC Encyclopedia of DNA Elements at UCSC 
       http://genome.ucsc.edu/ENCODE/

.. [#] Current ENCODE data
       https://www.encodeproject.org/

Problem 2
---------

1000 genomes VCF file is at::

    /vol1/opt/data/hg19/1000genomes/phase3.vcf.gz

#. Examine the relationship between SNP density and gene structure. Use
   1000Genomes SNP data in VCF format, and assume that each
   SNP in the file is the same in terms of quality (**20 points**).

    #. Calculate SNP density (SNPs per bp) for first exons of refGene
       annotations. (hint: use the ``-n`` option of 
       :ref:`bed12tobed6 <bedtools:bed12tobed6>`)

    #. Do the same for the introns of each gene. (hint: use
       :ref:`subtract <bedtools:subtract>` with BED files for genes and
       exons) 

    #. Make a histogram of the results (use ``geom_hist()``).
    #. Add text annotations to the plot with ``geom_text()`` and
       ``ggtitle``. E.g. how many SNPs were examined? How many exon /
       intron regions?

Problem Set Submission
----------------------
Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf
    PageBreak

