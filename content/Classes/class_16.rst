
.. include:: /_static/substitutions.txt

*****************************
Class 16 : Python : Parsing Fastqs
*****************************

:Class date: |c16-date|
:Last updated: |today|

Goals
=====
#. Finish concepts from Tuesday (list slices)
#. Generate executable scripts
#. Parsing Fastq exercises 


Python Ecosystem
================

- IPython is good for interactive work, but not for scripting

- Jupyter Notebooks is the python equivlant to RStudio and RMarkdown:
  Can run >40 languages, built first as IPython 
  https://jupyter.org/

- ``Pandas`` package for manipulating data.frames in python, ``dplyr`` was
  inspired by pandas

- ``Matplotlib`` for plotting in python

- ``numpy``, ``scipy``, and ``scikit-learn`` for math, science, and machine
  learning methods

  
Why Python
==========

- Reusable configurable scripts (i.e. write script to get 
  distribution of sequence lengths from ``FASTQ``, ``FASTA``, ``BAM``, ``BED``, etc)

- Parsing strange file formats (``FASTQ``, ``FASTA``, ``GTF``) that are not
  columnar or easily imported into ``R``

- Files too big to load into memory. ``R`` will generally load all of a
  file into your memory. What to do if file is 100Gb? Parse it line by
  line. 

Slicing of lists
================

.. code-block:: python

    nums = range(30)

    # try this again with implicit start and end
    nums[5:10]

    reversed(nums)

    # two ways to examine the contents of the iterator,
    # same principle for sorted()
    [i for i in reversed(nums)]

    list(reversed(nums))

    # skip reversed()
    range(30,0,-1)


Slicing of Strings
==================

.. code-block:: python

    dna = "ATCGATGCT"

    len(dna)
    
    # strings are iterable in python
    for base in dna:
        print(base)
        print(base.lower())
    
    # slicing works on strings

    dna[::2]
    
    # reverse the string
    dna[::-1]



Equality and Logic
==================
Use ``if``:``elif``:``else`` statements to test conditions and act on the
result. The ``==`` and ``!=`` operators test for equality and inequality, and
work on many object comparisons.

.. code-block:: python

    # animal colors
    cat = 'white'

    dog = 'black'

    if cat == dog: 
        print "same color" 
    elif cat != dog:
        print "different color" 
    else:
        print "not going to happen"      


Importing modules
=================
There are a number of modules with objects and functions in the standard
library, and there are a also a huge number of Python modules on the web
(check github).

To be able to access the contents of a module, you need to import it into
your `namespace`:

.. code-block:: python

    import math

    math.log10(1000)
    
    import math as m

    m.log10(1000)

    from math import log10

    log10(1000)


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


Hello World program
===================

Open up vim and save this script as hello-world.py

.. code-block:: python
    
    #! /usr/bin/env python3

    print("hello world, it's great to be alive")

.. code-block:: bash
    
    python3 hello-world.py

    #or make directly executable
    chmod u+x hello-world.py

    ./hello-world

In Class Exercises
==================

#.  Make a python program that takes a command line argument and prints
    the argument out as a string
#.  Build a parser for a fastq file (located here) ``data-sets/fastq/SP1.fq``
#.  Trim the first 10 nucleotides from the fastq file and return a fastq
    formatted file
#.  Remove fastq records whose sequences end with ``AAA``
#.  Count the number of bases found in the entire fastq (or alternatively,
    per record)
#.  Count the number of Unique Molecular Identifiers (``UMI_ATCGAT``) and
    the number of sequences and determine the sequence duplication rate (UMI
    counts / Total Sequences)

.. raw:: pdf

    PageBreak
