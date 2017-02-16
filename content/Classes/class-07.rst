
.. include:: /_static/substitutions.txt

=============================
Class 7 : BEDTools vignette
=============================

:Class date: |c7-date|
:Last updated: |today|

Goals
-----
#. Go over BEDTools analysis vignette to show how to combine different
   tools

There are several more vignettes `here.
<https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md>`_
Looking over these may be helpful to understand some of the questions on
the problem sets. 

Overview
--------
A common situation in genomic analysis is to examine how signals (e.g.,
coverage in a ChIP-seq experiment) relate to genome annotations (e.g.,
transcription start sites).

This simple process can be broken down into several individual steps:

#. Obtain genome annotations (we will use TSS)

#. Obtain signals in bedGraph format (we will use CTCF ChIP-seq)

#. Use :ref:`slop <bedtools:slop>` to create regions to examine [optional]

#. Break the regions into windows with :ref:`makewindows
   <bedtools:makewindows>`.

#. :ref:`map <bedtools:map>` data onto the windows.

#. Calculate a summary statistic for each window with :ref:`groupby
   <bedtools:groupby>`.

Annotations
-----------
First we need some annotations. Let's get some transcription start sites
(TSSs) from a BED file of gene annotations. For this analysis extract TSSs
from chr22. 

.. code-block:: bash

    $ genes="data-sets/bed/genes.hg19.bed.gz"
    $ tssbed="tss.bed"

    # need to grab the start or end based on the strand field
    $ gzcat $genes | awk '$6 == "+"' \
      | awk 'BEGIN {OFS="\t"} {print $1, $2, $2+1, $4}' > $tssbed
    $ gzcat $genes | awk '$6 == "-"' \
      | awk 'BEGIN {OFS="\t"} {print $1, $3, $3+1, $4}' >> $tssbed
    
    # extract chr22 intervals and sort bed output
    $ awk '($1 == "chr22")' $tssbed \
        | bedtools sort -i - > tmp.bed
    
    #rename tmp.bed to tss.bed
    $ mv tmp.bed $tssbed

Signal
------
The signal we will use is human CTCF ChIP-seq signal. The file is::

    signal="data-sets/bedtools/ctcf.hela.chr22.bg.gz"

Another common file format for ChIP/Clip-Seq signals is the `bigWig
<https://genome.ucsc.edu/goldenPath/help/bigWig.html>`_ format.

bigWig is a compressed form of bedGraph. You can either convert to bedGraph
and compress::

    bigWigToBedGraph <wig file> <bedgraph>
    gzip <bedgraph>

Or you can use it in a stream with::

    <(bigWigToBedGraph <wigfile> stdout)

Additional Encode data can be downloaded `here
<https://www.encodeproject.org/matrix/?type=Experiment>`_

Slop
----
Each of the TSS regions we made is 1 base in size. We need to make these
bigger so that we can examine coverage within a larger region. Let's make them
2 kb with the :ref:`slop <bedtools:slop>` command:

.. code-block:: bash

    $ chromsize="data-sets/bedtools/hg19.genome"
    $ slopbed="tss.slop.2000.bed"
    $ bedtools slop -b 2000 -i $tssbed -g $chromsize > $slopbed

Inspect the slop'd file and make sure you understand how it is different
from the input BED file. Check the new widths of each region using
``awk``.

Make windows
------------
Now that we have a 2 kb window around each TSS, we can break the regions up
into `windows` and examine signal in each. The idea here is that we would
like to see a peak of signal at the middle window, and it should descrease
as we move away from TSS.

.. code-block:: bash

    $ windowbed="tss.slop.2000.5bp.windows.bed"

    $ bedtools makewindows -b $slopbed -w 5 -i srcwinnum \
        | bedtools sort -i - \
        | tr "_" "\t" \
        > $windowbed 
    
Note that each window now has a unique number, and the relative positions
of each window of a given number are the same. For example, all the
windows labeled "1" are 2000 bases upstream of the original annotated TSS.

What would the signal look like if we didn't do this?

Map signal to the annotations
-----------------------------
Now we can map the signal to each of the windows that we made.

.. code-block:: bash

    $ signalmap="signal.map.bg"
    $ bedtools map -a $windowbed \
        -b $signal \
        -c 4 -o mean -null 0 \
        > $signalmap

- the ``-c`` flag specifies the values to operate on. 
- the ``-o`` flag specifices the operation to perform (e.g. mean)
- What does ``-null 0`` do? Why is it important?

Inspect this file and make sure you know what it looks like.

Grouping the data
-----------------
Now we can calculate summary statistics on the mapped data with :ref:`groupby
<bedtools:groupby>`:

.. code-block:: bash

    # Note that the mapped data has to be sorted by window number
    $ sort -t$'\t' -k5,5n $signalmap \
        | bedtools groupby \
            -i - \
            -g 5 -c 6 -o sum \
            > output.tab

- the ``-g`` flag specifies the column to group on
- the ``-c`` flag specifies the values to aggregate. 
- the ``-o`` flag specifices the operation to perform (e.g. sum)

Inspect the output so you know what it looks like.

Code
----

.. literalinclude:: /Classes/code/class-7.sh
    :language: bash
    :linenos:

Plotting
--------
Now we'll make a plot with R. Navigate to::

    http://amc-tesla.ucdenver.pvt/rstudio

and login with your tesla credentials. You should see R Studio open up.
Navigate to the directory where you did you work with the `Files` menu on
the lower right, and then set the working directory with `More` --> `Set
as Working Directory`.

Now make a plot!

.. code-block:: r
    
    > install.packages(ggplot2)
    > library(ggplot2)
    > col.names <- c('window','signal')
    > df <- read.table('output.tab', col.names=col.names)

    # coerce to numbers
    > df$signal <- as.double(df$signal)
    > df$window <- seq(-2000, 2000, 5)
    > qplot(window, signal, data = df)

Exercises
---------
Try to generalize this approach to other types of data. How would you:

#. Look at CTCF distribution near transcription termination sites
   (TTSs)? Write this program and run it. (easy)

#. Analyze the distribution of H3K36me3 ChIP-seq signal relative to coding
   exons?  Download some H3K36me3 data and :ref:`bed12tobed6
   <bedtools:bed12tobed6>`.

#. Analyze the distribution of H3K36me3 ChIP-seq relative to introns?
   Use the ``-n`` option in :ref:`makewindows <bedtools:makewindows>`.

#. Determine levels of evolutionary conservation surrounding different
   types of non-coding variants? (medium). Use dbSnp and phyloP scores.

