********************
Class XXX : ChIP-seq
********************

Goals
=====

 #. Learn workflow for analyzing ChIP-seq data
 #. Visualize data in the UCSC Genome Browser

.. _coverage-workflow:

Chromatin Immunoprecipitation Overview
--------------------------------------

A general workflow for visulazling ChIP-seq data (and many other types of
data) is:

.. list-table::
    :widths: 40 40
    :header-rows: 1

    * - Operation
      - File formats
    * - Align reads to reference genome
      - ``FASTQ -> BAM``
    * - Generate coverage plots
      - ``BAM -> bedGraph``
    * - Call peaks 
      - ``BAM -> BED``
    * - Make binary files for UCSC
      - ``bedGraph -> bigWig``, ``BED -> bigBed``
    * - Identify motifs
      - ``BED -> FASTA -> TXT``

.. _short-read-alignment:

Short read alignment
--------------------
There are several short read alignment packages available. We will mainly
use bowtie2 [#]_ because it is easy to use and is relatively fast.

.. code-block:: bash

    # minimal bowtie2 command
    $ bowtie2 -x <index> -U <read1.fq.gz> [options] > output.sam

.. [#] Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Normally, you would also pipe the SAM-format [#]_ alignment through the samtools
suite to discard unaligned reads, sort the alignment and store it in
binary format (bam).

.. code-block:: bash

    # sort sam output with samtools and write bam output
    $ bowtie2 -x <index> -U <read1.fq.gz> [options] \
        | samtools view -ShuF4 - \ 
        | samtools sort -o - aln.temp -m 8G \
        > aln.bam

.. [#] SAM format http://samtools.sourceforge.net/SAMv1.pdf

.. _coverage-plots:

Generate and visualize coverage plots
-------------------------------------
Once alignment is complete, you can create coverage plots from your aligned
data, so that you can visualize your data.

.. note::

    You will need a "chromsizes" flle for many BEDtools commands. This file
    contains the sizes of each chromosome in an assembly. UCSC provides a
    tool to retrieve this information:

    .. code-block:: bash

        # retrieve chrom sizes for the hg19 assembly and write them to a file
        # inspect this file so you know what it looks like
        $ fetchChromSizes hg19 > hg19.chrom.sizes

To generate coverage plots, we will use ``bedtools`` [#]_, a suite of tools
that we will cover in more detail later. Here, we'll use the :ref:`genomecov
<bedtools:genomecov>` tool.

.. code-block:: bash

    # -bg : write ouptput in bedGraph format
    $ bedtools genomecov -i <aln.bam> -g <chrom.sizes> -bg > coverage.bg

.. [#] BEDtools http://bedtools.readthedocs.org/en/latest/

This command writes a bedGraph format file called ``coverage.bg``. Use
``less`` to examine this file.

.. note::

    Words to live by: **If you make a BED file, sort the BED file**

    Many strange things can happen if you use unsorted BED files for
    analysis..  Once you create a BED file, sort it with one of these:

    .. code-block:: bash

        # same filename twice, overwrites original file
        $ bedSort file.bed file.bed

        # or you can use bedtools; writes additional file
        $ bedtools sort -i - < unsorted.bed > sorted.bed

.. _peak-calling:

Peak calling
------------
There are several available software packages for identying regions
encriched in your IP experiment (i.e. peaks). We will use macs2 here.

.. code-block:: bash

    # minimal macs2 command 
    $ macs2 callpeaks --treatment <aln.bam> --name <exp.name> [options]

.. _genome-browser-display:

Plot coverage with the Genome Browser
-------------------------------------

Use the UCSC Genome Browser to plot your data. Files in bedGraph format
can be large, so UCSC created a facility for posting binary format data in
a web-accessible directory that the browser can read.

.. code-block:: bash

    # convert bedGraph to binary format (bigWig) 
    $ bedGraphToBigWig <coverage.bg> <chrom.sizes> <coverage.bw> 

    # convert BED to binary format (bigBed)
    $ bedToBigBed <peaks.bed> <chrom.sizes> <peaks.bb>

XXX Post data to public_html directory

You can now write "tracklines" to tell where UCSC to find your data::

    # URL = http://amc-sandbox.ucdenver.edu/~username/path-to-binaryfile
    track type=bigWig bigDataUrl=<URL> name='coverage' color=r,g,b
    track type=bigBed bigDataUrl=<URL> name='peaks' color=r,g,b

.. note::

    Don't pick colors yourself, they will be ugly. **Use Colorbrewer**
    http://colorbrewer2.org.
    
    RGB colors in the ``Dark2`` and ``Set1`` qualitative palettes work
    well for UCSC display.

There are a large number of additional options you can use in tracklines
to change their display [#]_.

.. [#] UCSC Track configuration
       https://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK

.. _motif-identification:

Identify sequence motifs in enriched regions
--------------------------------------------

.. code-block:: bash

    # use ``qlogin`` to move to a compute node on tesla
    $ bedtools getfasta -fi <ref.fa> -bed <peaks.bed> -fo peaks.fa
    $ meme -nmotifs 100 -minw 6 -maxw 20 <peaks.fa>

.. _stranded-signals:

Split coverage by strand
------------------------
For some experiments, you will analyze the data relative to each strand of
the reference genome. For example, RNA is transcribed in single-stranded
form and derives from one or the other strand.

During alignment, reads from an RNA-based experiment will map to either
the positive ('+' or ``pos``) or negative ('-' or ``neg``) strand. You can
generate signal plots for ``pos`` and ``neg`` strands separately with
``bedtools``:

.. code-block:: bash

    $ common_args="-ibam <aln.bam> -g <chrom.size> -bg"
    $ bedtools genomecov $common_args -strand + > coverage.pos.bg
    $ bedtools genomecov $common_args -strand - > coverage.neg.bg

You would then create bigWigs for each of these display the stranded data
in the Genome Browser.

