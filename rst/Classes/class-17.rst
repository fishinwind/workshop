
*************************
  Class 17 : Bash and Bam
*************************

:Class date: 2014 Mar 7 Friday

Log in to amc-tesla and do this:

    cd ~
    cp /vol1/opt/data/class-17.tar.gz .
    tar xzvf class-17.tar.gz
    cd class-17

Goals
=====

 #. Debugging and Scripting in bash

 #. QC Sequence and Alignments


Pybedtools
==========

The pybedtools documentation is very good: https://pythonhosted.org/pybedtools/

Citable from: http://www.ncbi.nlm.nih.gov/pubmed/21949271

installable as:

.. code-block:: python

    pip install pybedtools

iterable:

.. code-block:: python

    for feature in BedTool('my.bed'):
        if feature.strand == "+":
            print feature.chrom, feature.start, feature.end


Bashing
=======

Some times you will have long bash scripts and you will misspell variables

.. code-block:: bash

    expname="chipseq.hela"
    bamfile=$exnpame.bam
   
The above code will run without error.
But, if you add:

.. code-block:: bash

    set -o nounset

To the **top of the script**. Accessing an undefined variable will raise an error.


Bashing (2)
===========

In a long bash script, you may have a series of commands:

.. code-block:: bash

    eco "hello world" > useful-file.txt
    cowsay < useful-file.txt

In this case, the first line will show an error, but the second will still run.
To make it **stop on the first error**, add the following to the **top of the
script**

.. code-block:: bash

    set -e

Bashing (3)
===========

Sometimes some of you were getting confused about what you were doing after all
of the variables. You can force `bash` to echo the expanded commands it is
running (including setting variable names) with

.. code-block:: bash

    set -x

Bashing Summary
===============

Do this at the top of every script:

.. code-block:: bash

   set -eo nounset -o pipefail
   set -x # this can sometimes be removed

Pipefail gives more useful error messages when piping (|) commands.


FASTQ
=====

remember fastq is [(name, seq, +, qual), ...]::

    @cluster_2:UMI_ATTCCG
    TTTCCGGGGCACATAATCTTCAGCCGGGCGC
    +
    9C;=;=<9@4868>9:67AA<9>65<=>591
    @cluster_8:UMI_CTTTGA
    TATCCTTGCAATACTCTCCGAACGGGAGAGC
    +
    1/04.72,(003,-2-22+00-12./.-.4-

We often want to see how quality scores degrade over the read,
check for adaptors, and see some info about our sequences...

FASTQC
======

fastqc is run as:

.. code-block:: bash

    fastqc /path/to/your/your.fastq

and it creates an output directory containing html, e.g.:

    http://amc-sandbox.ucdenver.edu/~brentp/fastqc/real_R1_fastqc/fastqc_report.html

FASTQC SP1
==========

#. Add a comment '#' before the line "<<FASTQC" in class-17/run.sh
#. Reason about what that block will do
#. Save and exit and run
#. Open the printed path in the browser


BAM
===

A BAM is **B** inary **A** lignment **F** ormat. It is the binary
version of SAM format. 
All of the alignments from high-throughput data you are likely to encounter will
be in BAM format.

You can easily transfer between Binary BAM and text SAM using samtools view:

    samtools view a.bam | python process-sam-text.py > processes.sam

http://samtools.sourceforge.net/samtools.shtml

Example Data
============

There are 4 small example BAM files in `/vol1/opt/data/bams/`

.. code-block:: bash

    ls -lh /vol1/opt/data/bams/*.bam

Since they are in binary format, you'll need to use samtools to `view` them


.. code-block:: bash

    # view the header:
    $ samtools view -H /vol1/opt/data/bams/2_8-bwa.bam | less
    # view the alignments:
    $ samtools view /vol1/opt/data/bams/2_8-bwa.bam | less

+ Alignments contain a lot of information!
+ Look at the format and read on the samtools site for more info.


picard
======

http://picard.sourceforge.net/

Picard has a number of tools for manipulating alignment files.

We will look at alignment metrics.

The 4 example bams are from a targetted sequencing project so we will
examine the percent on and off-target along with the coverage.

picard metrics
==============

 explain from terminal

picard output
=============

The output from picard is 1 sample per file with a bunch of extra lines.
We will parse them into a single, useful file with class-17/src/merge-metrics.py

Projects
========

come up with an idea for your projects.
