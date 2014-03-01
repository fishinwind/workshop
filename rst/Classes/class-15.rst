
*********************************
  Class 15 : ChIP-seq / DNase I
*********************************

:Class date: 2014 Mar 3 Monday

Goals
=====

 #. Post ChIP-seq data for the UCSC browser to see

 #. Learn to process DNase I data

ChIP-seq / DNase I data 
=======================

Look at some human ChIP-seq data [#]_.

.. [#] Genome Browser Session
       http://goo.gl/WfJxcM

Coverage plots for ChIP and DNase I
-----------------------------------

You will want to calculate coverage plots that are appropriate for the
type of experiment. For example, in a ChIP experiment, you want to examine
the entire region covered by sequences.

.. code-block:: bash

    $ common_args="-ibam <aln.bam> -g <chrom.size> -bg"
    $ bedtools genomecov $common_args > coverage.bg

But in a DNase I mapping experiment, you only want to know the exact
position where DNase I cut the DNA, i.e. the 5' position. Use the ``-5``
flag in bedtools to only count 5' positions.

.. code-block:: bash

    $ common_args="-ibam <aln.bam> -g <chrom.size> -bg"
    $ bedtools genomecov $common_args -5 > coverage.5p.bg

Plot coverage with the Genome Browser
-------------------------------------

Use the UCSC Genome Browser to plot your data. Files in bedGraph format
can be large, so UCSC created a facility for posting binary format data in
a web-accessible directory that the browser can read.

.. code-block:: bash

    # NOTE: bedGraph and BED files must be sorted first

    # convert bedGraph to binary format (bigWig) 
    $ bedGraphToBigWig <coverage.bg> <chrom.sizes> <coverage.bw> 

    # convert BED to binary format (bigBed)
    $ bedToBigBed <peaks.bed> <chrom.sizes> <peaks.bb>

Posting your data
-----------------

We need to set up your accounts so that you will be able to see data
posted in your account. When this is setup, you should be able to put
files in::

    $HOME/public_html/file.html

and be able to see them like::

    http://amc-sandbox.ucdenver.edu/~username/file.html

Writing tracklines
------------------

You can now write "tracklines" [#]_ to tell where UCSC to find your data::

    # URL = http://amc-sandbox.ucdenver.edu/~username/path-to-binaryfile
    track type=bigWig bigDataUrl=<URL> name='coverage' color=r,g,b
    track type=bigBed bigDataUrl=<URL> name='peaks' color=r,g,b

.. [#] UCSC Track configuration
       https://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK

.. nextslide::

.. tip::

    Don't pick colors yourself, they will be ugly. **Use Colorbrewer**
    http://colorbrewer2.org.
    
    RGB colors in the ``Dark2`` and ``Set1`` qualitative palettes work
    well for UCSC display.

There are a large number of additional options you can use in tracklines
to change their display.

