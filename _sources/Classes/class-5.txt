***************************************
  Week 3 : Working with genomic data 
***************************************

:Class date: Tues 10 Feb 2015
:Class date: Thurs 12 Feb 2015

Goals
=====

#. Learn to run scripts on the cluster via the queuing system

#. Learn about genomic data types and where to get data 
 
#. Start to use ``bedtools`` to analyze genomic data

grep
====
Use :linuxman:`grep(1)` to identify lines in a file that match a specified pattern.

To find any instance of *chr5* in the lamina.bed file

.. code-block:: bash

    # grep [pattern] [filename]
    $ grep chr5 /vol1/opt/data/lamina.bed | head

To find all lines that start with a number sign:

.. code-block:: bash

    # The caret (^) matches the beginning of the line
    # FYI dollar sign ($) matches the end
    $ grep '^#' /vol1/opt/data/lamina.bed

.. nextslide::
    :increment:

To find any line that *does not* start with "chr":

.. code-block:: bash

    # the -v flag inverts the match (grep "not" [pattern])
    $ grep -v '^chr' /vol1/opt/data/lamina.bed

Beware of using ``grep`` to find patterns that might be partial matches:

.. code-block:: bash

    # this will match chr1, chr10, chr11 etc.
    $ grep chr1 /vol1/opt/data/lamina.bed | cut -f1 | uniq

You can find exact matches that are split on words with the ``-w`` flag:

.. code-block:: bash

    # this will only match chr1
    $ grep -w chr1 /vol1/opt/data/lamina.bed | cut -f1 | uniq

.. nextslide::
    :increment:

Beware of using ``grep`` to search for numbers:

.. code-block:: bash

    # finds all strings that match `100`
    $ grep 100 /vol1/opt/data/lamina.bed | head -n 20

    # better, but doesn't look at numeric value
    $ grep -w 100 /vol1/opt/data/lamina.bed | head -n 20

.. tip::

    If you're trying to find numeric values in a file, use ``awk``
    instead::

        $ awk '$2 == 500' /vol1/opt/data/lamina.bed

Cluster access
==============
We have set up accounts for the class on our departmental cluster. We will
set up your accounts at the end of class and reset your passwords:

.. code-block:: bash

    # the -X flag starts an X11 connection 
    $ ssh -X username@amc-tesla.ucdenver.pvt

    ...

    # once you are logged in, text your X11 connection with
    $ xeyes

Cluster etiquette
=================
There are some specific rules you need to know when you're operating in a
cluster environment.

.. graphviz::

    digraph cluster {
        "YOU" [shape=box];
        "amc-tesla" [shape=box];
        "filesystem" [shape=box];
        "compute nodes" [shape=box];
        "YOU" -> "amc-tesla";
        "amc-tesla" -> "filesystem";
        "amc-tesla" -> "compute nodes";
    }

.. important::

  **DO NOT** run jobs on the head node (amc-tesla). The head node is the
  brains of the cluster and it can easily be overextended. Use ``qlogin``
  instead.

Example commands on the cluster
===============================
Find the size of the file system:

.. code-block:: bash

    $ df -h

Find how much space you have allocated:

.. code-block:: bash

    $ quota -h

The queueing system
===================
First you will grab a single CPU from the queueing system so that you can
work without affecting the head node. We use ``qlogin`` for this:

.. code-block:: bash

    jhessel@amc-tesla ~
    $ qlogin 

    Job <492536> is submitted to queue <interactive>.
    <<ssh X11 forwarding job>>
    <<Waiting for dispatch ...>>
    <<Starting on compute00>>

    jhessel@compute00 ~
    $ 

.. note:: 

    The host in the prompt changed from ``amc-tesla`` to ``compute00``.
    
You can now execute long-running processes without worry of affecting the
cluster. Type ``exit`` to return back to your head node login.

.. nextslide::
    :increment: 

The cluster uses a queueing system that will run jobs that you submit to
it. You can write a small test script to see how the system works. First,
write this into a run.sh file:

.. code-block:: bash

    #!/usr/bin/env bash

    #BSUB -J sleeper
    #BSUB -e %J.err
    #BSUB -o %J.out

    sleep 20

.. nextslide::
    :increment: 

The ``#BSUB`` lines are comments, but are read by the ``bsub`` program to
identify features associated with your job. 

- ``-J`` sets the job's name
- ``%J`` is a unique job ID that is set when you run the job.
- ``-e`` and ``-o`` set the filenames for stderr and stdout from the job

.. nextslide::
    :increment: 

Now you can submit the script to the queuing system. As soon as you submit
it, you can check on its progress:

.. code-block:: bash

    $ bsub < run.sh
    $ bjobs

After the job finishes, you should see two new files that end
`.out` and `.err`; these stdout and stderr from the running job.
Look at the contents of those files so you know what is in
each one.

Killing jobs
============
Sometimes you need to kill your jobs. You can kill specific jobs using
their job ID numbers, obtained from checking ``bjobs``:

.. code-block:: bash

    $ bkill <jobid> 

You can also kill **all** of your jobs at once:

.. code-block:: bash

    $ bkill 0 

.. warning::

    ``bkill 0`` is dangerous – it will wipe out all of your jobs. If
    you have long-running jobs that you forgot about, you will kill them
    too if you are not careful!

Other cluster-specific commands
===============================
.. code-block:: bash

    $ bhosts     # hosts in the cluster
    $ man bhosts # bsub man page
    $ bqueues    # available queues
    $ lsload     # check load values for all hosts

Questions
=========

#. Check the ``.err`` files from the run. What information do they
   contain? What does this tell you about your starting sequences?

#. Find out how you would modify the ``bowtie2`` command to write out the
   unaligned reads into a new file. Re-run the analysis to report those
   reads.

#. Modify the ``awk`` command in the script to print out a valid BED4
   format::
    
        chrom <tab> start <tab> end <tab> count
    
#. Find out how many unique UMI sequences are associted with each
   chromosomal coordinate (note: not as easy).

More exercises
==============

#. use ``grep`` to identify lines in lamina.bed where the second field
   (start) begins with ``100``.

#. use ``grep`` to identify lines in lamina.bed where the third field
   (end) ends with 99 .

#. use ``grep`` with its ``-w`` flag to count the number of 'chr1'
   records in lamina.bed.

#. use ``grep`` to count how many fastq records are in the
   /vol1/opt/data/t_R1.fastq.gz file (fastq records begin with an
   '@' symbol)

#. login to amc-tesla. use ``grep`` to count the number of fastq records
   in /vol1/opt/data/SP1.fq.gz


.. raw:: pdf

    PageBreak

********************
      BEDTools
********************

Goals
=====

#. Introduce the BEDTools suite of tools.
#. Understand why using BEDTools is needed.
#. Practice common operations on BED files with BEDTools.

BEDTools Overview
=================

BEDTools will be one of the tools with the best return on investment. For
example, to extract out **all genes that overlap a CpG island**:

.. code-block:: bash

    $ bedtools intersect -u -a genes.hg19.bed.gz -b cpg.bed.gz \
                                     > genes-in-islands.bed

:ref:`intersect <bedtools:intersect>` is a bedtools tool. It follows a
common pattern in bedtools that the query file is specified after the
``-a`` flag and the *subject* file after the ``-b`` flag

BEDTools Utility
================

Finding all overlaps between a pair of BED files naively in python would look like:

.. code-block:: python

    for a in parse_bed('a.bed'):
        for b in parse_bed('b.bed'):
            if overlaps(a, b):
                # do stuff

If *'a.bed'* has 10K entries and *'b.bed'* has 100K entries, this would
involved checking for overlaps **1 billion times**. That will be slow.

BEDTools uses an indexing scheme that reduces the number of tests
dramatically.

.. note::
  
  See the original BEDTools paper for more information:
  http://bioinformatics.oxfordjournals.org/content/26/6/841.full

.. nextslide::
   :increment:

+ Fast: faster than intersect code you will write
+ Terse: syntax is terse, but readable
+ Formats: handles BED, VCF and GFF formats (gzip'ed or not)
+ Special Cases: handles stranded-ness, 1-base overlaps, abutted intervals,
  etc. (likely to be bugs if you do code in manually)

BEDTools Commands
=================

To see all available BEDTools commands, type

.. code-block:: bash

    $ bedtools

The most commonly used BEDtools are:

+ :ref:`intersect <bedtools:intersect>`
+ :ref:`genomecov <bedtools:genomecov>`
+ :ref:`closest <bedtools:closest>`
+ :ref:`map <bedtools:map>`

BEDTools Documentation
======================

The BEDTools documentation is quite good and ever improving.

See the documentation for :ref:`intersect <bedtools:intersect>` with:

.. code-block:: bash

    $ bedtools intersect

The online HTML help is also good and includes pictures: 
 https://bedtools.readthedocs.org/en/latest/content/tools/intersect.html

BEDTools intersect
==================
Have a browser window open to :ref:`BEDTools intersect documentation <bedtools:intersect>`.
It will likely be the BEDTools function that you use the most. It has a lot of
options.

.. image:: http://bedtools.readthedocs.org/en/latest/_images/intersect-glyph.png

"-v" means (like grep) include all intervals from `-a` that do not overlap
intervals in `-b`

Example Files
=============

.. code-block:: bash

    $ cat a.bed 
    chr1    10  20  a1  1   +
    chr1    100 200 a2  2   -

    $ cat b.bed 
    chr1    20  30  b1  1   +
    chr1    90  101 b2  2   -
    chr1    100 110 b3  3   +
    chr1    200 210 b4  4   +

What will happen if you intersect those files?
For example, the *a.bed* region `chr1:100-200` overlaps::

    chr1:90-101 
    chr1:100-110

from *b.bed*

intersect
=========

intersect with default arguments means **extract chunks of `-a` that overlap
regions in `-b`**

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed
    chr1    100 101 a2  2   -
    chr1    100 110 a2  2   -

Here is the original interval from *a.bed*::

    chr1	100	200	a2	2	-

And the overlapping intervals from *b.bed*::

    chr1	90	101	b2	2	-
    chr1	100	110	b3	3	+

intersect -wa
=============

Often, we want the *entire interval from -a if it overlaps any interval in -b*

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -wa
    chr1    100 200 a2  2   -
    chr1    100 200 a2  2   -

We can get that uniquely with (-u)

intersect -wo
=============

We can see which intervals in *-b* are associated with *-a*

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -wo
    chr1  100  200  a2  2  -  chr1  90  101  b2  2  -  1
    chr1  100  200  a2  2  -  chr1  100  110  b3  3  +  10

intersect exercise
==================

What happens if you reverse the arguments? E.g. instead of::

  -a a.bed -b b.bed

use::

   -b a.bed -a b.bed

Try that with no extra flags, with -u, -wa, -wo.

How does it compare to the original?

intersect -c
============

We can count overlaps for each interval in *-a* with those in *-b* with

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -c
    chr1	10	20	a1	1	+	0
    chr1	100	200	a2	2	-	2

This is our original `a.bed` with an **additional column indicating number of
overlaps** with `b.bed`

intersect -v
============

Extract intervals in `a.bed` that do not overlap any interval in `b.bed`

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -v
    chr1	10	20	a1	1	+

Extract intervals in `b.bed` that do not overlap any interval in `a.bed`

.. code-block:: bash

    $ bedtools intersect -a b.bed -b a.bed -v
    chr1	20	30	b1	1	+
    chr1	200	210	b4	4	+

Intersect Summary
=================

+ fragments of `a` that overlap `b`:
  `intersect -a a.bed -b b.bed`
+ complete regions of `a` that overlap `b`:
  `intersect -a a.bed -b b.bed -u`
+ intervals of `b` as well as `a`:
  `intersect -a a.bed -b b.bed -wo`
+ number of times each `a` overlaps `b`:
  `intersect -a a.bed -b b.bed -c`
+ intervals of `a` that do not overlap `b`:
  `intersect -a a.bed -b b.bed -v`

Exercises (Or Other Tools)
==========================

#. zless :download:`cpg.bed.gz <../misc/data/cpg.bed.gz>` and :download:`genes.hg19.bed.gz <../misc/data/genes.hg19.bed.gz>`
#. Extract the CpG islands that touch any gene [**24611**]
#. Extract CpG islands that do not touch any gene [**7012**]
#. Extract (uniquely) all of each CpG Island that touches any gene [**21679**]
#. Extract CpG's that are completely contained within a gene (look at the help
   for a flag to indicate that you want the fraction of overlap to be 1 (for 100 %). [**10714**]
#. Report genes that overlap any CpG island. [**16908**]
#. Report genes that overlap more than 1 CpG Island (use -c and awk). [**3703**].

.. note::

    as you are figuring these out, make sure to pipe the output to less or head

Other Reading
=============

+ Check out the online `documentation <https://bedtools.readthedocs.org/en/latest/content/tools/intersect.html>`_.
+ A `tutorial <http://quinlanlab.org/tutorials/cshl2013/bedtools.html>`_ by the author of BEDTools

Intersect Bam
=============

We have seen that `intersect <bedtools:intersect>` takes `-a` and `-b`
arguments. It can also intersect against an alignment BAM file by using `-abam`
in place of `-a`

e.g:

.. code-block:: bash

    $ bedtools intersect \
        -abam experiment.bam \
        -b target-regions.bed \
        > on-target.bam

Intersect Strand
================

From the `help <https://bedtools.readthedocs.org/en/latest/content/tools/intersect.html>`_ ,
one can see that intersect can consider strand. For example if both files have a
strand field then

.. code-block:: bash

    $ bedtools intersect -a a.bed -b b.bed -s

Will only consider as overlapping those intervals in `a.bed` that have the same
strand as `b.bed`.

Closest
=======

with :ref:`intersect <bedtools:intersect>` we can only get overlapping
intervals. :ref:`closest <bedtools:closest>` reports the nearest interval even
if it's not overlapping. 

Example: report the nearest CpG to each gene as long as it is within 5KB.

.. code-block:: bash

    bedtools closest \
        -a genes.hg19.bed.gz \
        -b cpg.bed.gz -d \
        | awk '$NF <= 5000'

Map
===

For each CpG print the sum of the values (4th column) of overlapping intervals from
lamina.bed (and filter out those with no overlap using awk)

.. code-block:: bash

    $ bedtools map \
        -a cpg.bed.gz \
        -b /vol1/opt/data/lamina.bed \
        -c 4 -o sum \
        | awk '$5 != "."'

Other *-o* perations include **min**, **max**, **mean**, **median**, **concat**

Sorted
======

When you start dealing with larger data-files. Look at the `-sorted` flag.
For example in :ref:`intersect <bedtools:intersect>`.

+ Uses less memory
+ Faster

Takes advantage of sorted chromosome, positions in both files so it doesn't have
to create an index.

.. image:: http://bedtools.readthedocs.org/en/latest/_images/speed-comparo.png

Genomecov
=========

Get coverage of intervals in BED by BAM 

.. image:: https://bedtools.readthedocs.org/en/latest/_images/genomecov-glyph.png

Usually want the last option `-bg -split`

ENCODE
======
 
The Human Genome Project was finished, giving us a list of human genes and their 
locations. Unfortunately, we still had no idea how they were regulated. If only 
there was an `ENCyclopedia Of Dna Elements 
<http://www.sciencemag.org.hsl-ezproxy.ucdenver.edu/content/306/5696/636.full>`_…

Advantages: massive amounts of information on key cell lines, reproducible 
experiments, public data access, technology development.

ENCODE Project Cell Lines
=========================

Tier 1: GM12878 (EBV-transformed lymphoblast), K562 (CML lymphoblast), H1-hESC

Tier 2: HeLa-S3 (cervical cancer), HepG2 (liver carcinoma), HUVEC (umbilical vein)

Tier 2.5: SKNSH (neuroblastoma), IMR90 (lung fibroblast), A549 (lung carcinoma), 
MCF7 (breast carcinoma), LHCN (myoblast), CD14+, CD20+
 
`link <http://genome.ucsc.edu/ENCODE/cellTypes.html>`_ (this page also has very useful
links to cell culture protocols)

Experiments
===========

#. ChIP-seq: Histone marks, transcription factors

#. Chromatin structure: DNaseI-seq, FAIRE, 5C/Hi-C

#. RNA expression: mRNA-seq, GENCODE gene predictions

#. Data Integration: Segway / ChromHMM integration of functional data

Common File Formats
===================

+ FASTQ: Raw sequencing data. `[link] <http://maq.sourceforge.net/fastq.shtml>`
+ SAM/BAM: Aligned sequence data `[link] <http://samtools.github.io/hts-specs/SAMv1.pdf>`
+ Bed/bigBed: List of genomic regions `[link] <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`
+ Bedgraph/Wig/bigWig: Continuous signal `[link] <http://genome.ucsc.edu/goldenPath/help/bedgraph.html>` 

Many other formats are described on this `page <http://genome.ucsc.edu/FAQ/FAQformat.html>`_

References
==========

Completion of the entire project, and a ton of papers: 
`Nature <http://www.nature.com/nature/journal/v489/n7414/index.html>`_, 
`Genome Research <http://genome.cshlp.org/content/22/9.toc>`_, 
`Genome Biology <http://genomebiology.com/content/13/9>`_, 

How to Access ENCODE Data
=========================

The `ENCODE project page <https://www.encodeproject.org/>`_ is the portal
to all of the ENCODE data.

**************************************
Chromatin Immunoprecipitation Overview
**************************************

Chromatin Immunoprecipitation is used to determine where a protein of
interest binds on a chromatin template [Park_Chipseq]_.

.. [Park_Chipseq] http://www.nature.com/nrg/journal/v10/n10/full/nrg2641.html

.. image:: ../_static/images/chip-workflow.png

.. nextslide::

.. image:: ../_static/images/chip-data.png

ChIP-seq analysis workflow
==========================

A general workflow for visualizing ChIP-seq data (and many other types of
data) is:

.. list-table::
    :widths: 40 40
    :header-rows: 1

    * - Operation
      - File formats
    * - Align reads to reference genome
      - ``FASTQ ~~> BAM``
    * - Generate coverage plots
      - ``BAM ~~> bedGraph``
    * - Call peaks 
      - ``BAM ~~> BED``
    * - Make binary files for UCSC display
      - ``bedGraph ~~> bigWig``, ``BED ~~> bigBed``
    * - Identify motifs
      - ``BED ~~> FASTA ~~> TXT / HTML``

ChIP-seq data
=============

Look at some human ChIP-seq data [#]_.

.. [#] Genome Browser Session
       http://goo.gl/WfJxcM

(We'll talk more in depth about ChIP-Seq workflows in the future,
but for now, just a brief introdcution to a few commands you can use to
work on ChIP-Seq data (and pset3)).


Peak calling
============

There are several available software packages for identifying regions
enriched in your IP experiment (i.e. peaks). We will use macs2 here.

.. code-block:: bash

    # minimal macs2 command 
    $ macs2 callpeak --treatment <aln.bam> --name <exp.name> [options]

Identify sequence motifs in enriched regions
============================================

You can use meme [#]_ to identify over-represented motifs in groups of
sequences (e.g. sequences covered by ChIP peaks).

Use the :ref:`bedtools getfasta <bedtools:getfasta>` command to fetch
fasta sequences.

Note: meme looks at both strands of a DNA sequence by default.

.. [#] MEME 
       http://meme.nbcr.net/meme/

.. code-block:: bash

    $ bedtools getfasta -fi <ref.fa> -bed <peaks.bed> -fo peaks.fa
    $ meme -nmotifs 5 -minw 6 -maxw 20 -dna <peaks.fa>


