*****************************
Class 9 : Intermediate Python 
*****************************

Goals
=====

 #. review
 #. intermediate concepts
 #. start reading useful python programs
 #. start writing useful python programs 
 #.

Review
======

....

Intermediate Concepts: Streaming
================================

One of the reasons why python is so useful is that faciliates
**iteration** over a file or *iterable* without reading the entire 
dataset into computer memory.

This is similar to streaming data in the Linux tools we've discussed.
For example

.. code-block:: bash

    zless /opt/bio-workshop/data/t_R1.fastq.gz | wc -l

never holds the file in memory, it just streams the data.

We can do this in python.

Intermediate : Streaming
========================

.. warning:: 

    DO NOT DO THIS!!

.. code-block:: python

    data = list(gzip.open('/opt/bio-workshop/data/t_R1.fastq.gz'))
    lines = len(data)

.. important:: 

    DO THIS

.. code-block:: python

    lines = sum(1 for line in gzip.open('opt/bio-workshop/data/t_R1.fastq.gz'))
    # or:
    lines = 0
    for line in gzip.open('/opt/bio-workshop/data/t_R1.fastq.gz'):
        lines += 1


Intermediate : Streaming with generators
========================================

.. code-block:: python

    def bed_generator(bed_file):
        for line in open(bed_file):
            chrom, start, end, value = line.split("\t")[:4]
            start, end = int(start), int(end)
            yield dict(chrom=chrom, start=start, end=end, value=value)

Then use it:

.. code-block:: python

    for bed in bed_generator('/opt/bio-workshop/data/lamina.bed'):
        print bed


Useful python modules
=====================
There are several modules in the standard library you will use all the
time:

    - :py:mod:`sys`: :py:obj:`sys.argv` has all the arguments from the command
      line

    - :py:mod:`collections`: espcially :py:class:`~collections.defaultdict`
      and :py:class:`~collections.Counter`

    - :py:mod:`itertools`: tools for efficient aggregation and iteration

    - :py:mod:`argparse`: command line option parsing


In Class Exercise
=================

 #. foo

.. raw:: pdf

    PageBreak
