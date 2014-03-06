
*************************
  Class 17 : Bash and Bam
*************************

:Class date: 2014 Mar 7 Friday

Log in to amc-tesla.

Goals
=====

 #. Debugging in bash

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

BAM
===

A BAM is **B** inary **A** lignment **F** ormat. It is the binary
version of SAM format. 
All of the alignments from high-throughput data you are likely to encounter will
be in BAM format.

Example Data
============

There are 4 example BAM files in ?????



picard
======



metrics

samtools
========

view alignments...


Projects
========

come up with an idea for your projects.
