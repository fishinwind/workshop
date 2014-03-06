
*********************************
  Class 17 : 
*********************************

:Class date: 2014 Mar 7 Friday

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

(How it began)


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

fastq is ...

FASTQC
======

BAM
===

bam is



picard
======

metrics

samtools
========

view alignments...


Projects
========

come up with an idea for your projects.
