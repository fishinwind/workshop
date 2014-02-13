************************
Class 9 : Applied Python
************************

Run this in **bash**

.. code-block:: bash

       python -c "import toolshed"

If you see an error get help to install toolshed

Goal
====

Take the basic concepts we've learned and do something useful.

toolshed
========

`toolshed <https://pypi.python.org/pypi/toolshed>`_ is a python module
that simplifies common file/text-processing tasks.
For example, it assumes the first line of a file is the header
and gives a python dictionary for each line keyed by the header.

.. code-block:: python

    from toolshed import reader

    for region in reader('/opt/bio-workshop/data/lamina.bed'):
        # the first line in lamina.bed is: '#chrom  start  end  value'

        if region['chrom'] != "chr12": continue
        if float(region['value']) < 0.90: continue
        print region['chrom'], region['start'], region['end']

toolshed
========

The toolshed reader function can also take gzipped files, files
over http, bash commands, and (some) xls files.

It can also accept a python class, that, for example
converts start and end to int's.

Mostly we will use it as:

.. code-block:: python

    from toolshed import reader
    for region in reader('/opt/bio-workshop/data/lamina.bed'):
        # do something with region
        print region['chrom']


Application: Setup
==================

We have 3 sets of data:

#. a set of paired-end FASTQ sequence files
#. a file that maps the FASTQ file name to a sample-id
#. a file that maps a sample-id to a phenotype.

We need to integrate these 3 so that we know, for example which
FASTQ files are associated with which phenotype.

Application: Desired Output
===========================

The output will be a tab-delimited file with columns for

#. sample-id
#. phenotype
#. R1 fastq name
#. R2 fastq name
#. other clinical or lab information ...


Looping Review
==============
Before we move on to more advanced topics, let's do some exercises 
reviewing how loops work. First, use wget to download 
`hamlet.txt <http://www.cs.uni.edu/~schafer/1140/assignments/pa11/hamlet.txt>`_. 

In Class Exercise
=================

From hamlet.txt: 

#. Print the first word of each line.

#. Print only lines that are not indented.

#. Count the number of times that the word "therefore" appears.


Counters
========
`Counters <http://docs.python.org/2/library/collections.html>`_ are a 
specialized version of Python dict objects where the keys are objects, and 
the values are counts of those objects. Counter objects have several specific
methods that take advantage of this, like most_common(). 

.. code-block:: python

    from collections import Counter
    c = Counter()
    for line in open('hamlet.txt'):
        words = line.strip().split(' ')
        for word in words:
            c[word] += 1
    print c.most_common(5)


Looping: Reading Multiple Lines at a Time
=========================================
There are lots of biological data files that have information for one thing
spread over multiple lines. A common example is a FASTA file, which is used to 
store sequences. Each sequence has a line with '>' and some information (like a name)
followed by another line of sequence data. For example: 

.. code-block:: fasta

    >Sequence name
    AGCATCGTAGCTAGTCGTACGTAGCTATCGATCGTAGCTA

In Class Exercise
=================

#. Open sample.fas and make a dictionary with four items corresponding to the sequences 
   from the file

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

    DO NOT DO THIS!! It reads everything into memory.

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


Streaming with yield
===================================

Make a bed reader that returns a useful dict

.. code-block:: python

    def bed_generator(bed_file):
        for line in open(bed_file):
            if line.startswith('#'): continue
            chrom, start, end, value = line.split("\t")[:4]
            start, end = int(start), int(end)
            yield dict(chrom=chrom, start=start, end=end, value=value)

Then use it:

.. code-block:: python

    for bed in bed_generator('/opt/bio-workshop/data/lamina.bed'):
        print bed # bed is a useful, usable thing. with numeric start and end.

Note that only ever have 1 (**) line in memory at a time.

In Class Exercise
=================

#. Modify the `bed_generator` code from the previous slide so that it
   turns value into a :py:obj:`float` before yielding
#. In the code that calls bed_generator, print out the value
#. In the code that calls bed_generator, append value to a list.

In Class Exercise (Answer)
==========================

.. code-block:: python

    def bed_generator(bed_file):
        if line.startswith('#'): continue
        for line in open(bed_file):
            chrom, start, end, value = line.split("\t")[:4]
            start, end = int(start), int(end)
            yield {'chrom': chrom, 'start': start, 'end': end,
                   'value': float(value))}

    vals = []
    for bed in bed_generator('/opt/bio-workshop/data/lamina.bed'):
        print bed['value']
        vals.append(bed['value'])

    print vals[:10]
    print sum(vals)


.. raw:: pdf

    PageBreak
