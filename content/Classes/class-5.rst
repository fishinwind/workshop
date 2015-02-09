***************************************
 Class 5/6 : Working with genomic data 
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


.. raw:: pdf

    PageBreak
