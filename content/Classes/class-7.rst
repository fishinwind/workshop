=================
BEDTools vignette
=================

:Class date: T Feb 17

Goals
-----
#. Go over BEDTools analysis vignette to show how to combine different
   tools

There are several vignettes `here.
<https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md>_`

Overview
--------
A common situation in genomic analysis is to examine how signals (e.g.,
coverage in a ChIP-seq experiment) relate to genome annotations (e.g.,
transcription start sites).

This simple process can be broken down into several individual steps:

#. Obtain genome annotations (we will use TSS)

#. Obtain signals in bedGraph format (we will use Pol II ChIP-seq)

#. Use `slop <bedtools:slop>` to create regions to examine [optional]

#. Break the regions into windows with `makewindows
   <bedtools:makewindows>`.

#. `map <bedtools:map>` data onto the windows.

#. Calculate a summary statistic for each window with `groupby
   <bedtools:groupby>`.

Annotations
-----------
First we need some annotations. Let's get some transcription start sites
(TSSs) from a BED file of gene annotations.

.. code-block:: bash

    $ genes=/vol1/opt/data/refGene.bed.gz
    $ tssbed=tss.bed

    # need to grab the start or end based on the strand field
    $ zcat $genes | awk '$6 == "+" | awk '{print $1,$2,$2+1}' > $tssbed
    $ zcat $genes | awk '$6 == "-" | awk '{print $1,$3,$3+1}' > $tssbed

    # we made a bed file, so now we sort it.
    $ bedSort $tssbed $tssbed
    $ gzip $tssbed

Signal
------
The signal we will use is human Pol II ChIP-seq signal. The file is::

    signal=/vol1/opt/data/endcode/wgEncodeBroadHistoneHelas3Pol2bStdSig.bigWig

bigWig is a compressed form of bedGraph. You can either convert to bedGraph
and compress::

    bigWigToBedGraph <wig file> <bedgraph>
    gzip <bedgraph>

Or you can use it in a stream with::

    <(bigWigToBedGraph <wigfile> stdout)

Slop
----
Each of the TSS regions we made is 1 base in size. We need to make these
bigger so that we can examine coverage in a larger region. Let's make them
2 kb with the `slop <bedtools:slop` command:

.. code-block:: bash

    $ chromsize=/vol1/opt/data/hg19.chrom.sizes
    $ slopbed=tss.slop.2000.bed
    $ bedtools slop -b 2000 -i $tssbed -g $chromsize > $slopbed

Inspect the slop'd file and make sure you understand how it is different
from the input BED file.

Make windows
------------
Now that we have a 2 kb window around each TSS, we can break the regions up
into `windows` and examine signal in each. The idea here is that we would
like to see a peak of signal at the middle window, and it should descrease
as we move away from TSS.

.. code-block:: bash

    $ windowbed=tss.slop.2000.5bp.windows.bed

    $ bedtools makewindows -b $slopbed -w 5 -i srcwinnum \
        sort -k1,1 -k2,2n \
        tr "_" "\t" \
        > $windowbed 
    
What would the signal look like if we didn't do this?

Map signal to the annotations
-----------------------------
Now we can map the signal to each of the windows that we made.

.. code-block:: bash

    $ signalmap=signal.map.bg
    $ bedtools map -a $windowbed \
        -b <(bigWigToBedGraph $signal stdout) \
        -c 4 -o mean -null 0 \
        > $signalmap

Inspect this file and make sure you know what it looks like.

Grouping the data
-----------------
Now we can calculate summary statistics on the mapped data with `groupby
<bedtools:groupby>`:

.. code-block:: bash

    # Note that the mapped data has to be sorted by window number
    $ sort -t$'\t' -k5,5n $signalmap \
        bedtools groupby \
            -i - \
            -g 5 -c 6 -o sum \
            > output.tab

Inspect the output so you know what it looks like.

Plotting
--------
Now we'll make a plot with R. Navigate to::

    http://amc-tesla.ucdenver.pvt/rstudio

and login with your tesla credentials. You should see R Studio open up.
Navigate to the directory where you did you work with the `Files` menu on
the lower right.

.. code-block:: r

    > library(ggplot2)
    > col.names <- c('window','signal')
    > df <- read.table('output.tab')
    > qplot(data = df,

