
.. include:: /_static/substitutions.txt 

.. _problem-set-4:

*************
Problem Set 4
*************

:Due date: |pset4-due| 

Overview
--------

For this quiz you will use BEDTools and command line tools to make some
simple summary files and plots.

Problem 1
---------

#. Examine Pol II distribution in a 5 kb window around transcription
   termination sites (TTSs)?  Make a plot of the result and
   save it as a PDF. (**15 points**)

#. Pick one other factor you might expect to find at a TTS iand analyze its
   distribution. You can find bigWig files of ChIP-seq experiemnts done for a
   variety of factors at UCSC [#]_ and the ENCODE [#]_ website. Make a plot of the
   result and save it as a PDF. (**15 points**)

.. [#] ENCODE Project at UCSC Encyclopedia of DNA Elements at UCSC 
       http://genome.ucsc.edu/ENCODE/

.. [#] Current ENCODE data
       https://www.encodeproject.org/

Problem 2
---------

#. Examine the relationship between SNP density and gene structure. Use
   1000Genomes SNP data in VCF format, and assume that each
   SNP in the file is the same in terms of quality (**20 points**).

        #. Calculate SNP density (SNPs / bp) for first exons of refGene
           annotations. Pay attention to strand.
          
        #. Do the same for the first intron of each gene. 

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

