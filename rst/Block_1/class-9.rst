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

One of the reasons why python is so useful is that is allows one to
**iterate** over a file or *iterable* without reading the entire 
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

.. warning:: 

    DO THIS

.. code-block:: python

    lines = sum(1 for line in gzip.open('opt/bio-workshop/data/t_R1.fastq.gz'))


Intermediate : Streaming with generators
========================================


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

Debugging Python code
=====================
The :py:mod:`pdb` is the Python Debugger. You can use it to debug programs by
dropping you into a shell that allows you to step through the program, line by
line.

.. ipython::
    :verbatim:

    In [6]: import pdb

    # this will drop you into a shell. find the value of ``i`` at the (Pdb)
    # prompt
    In [7]: for i in range(100):
       ...:     if i == 50:
       ...:         pdb.set_trace()
       ...:         


In Class Exercise
=================

 #. foo

.. raw:: pdf

    PageBreak
