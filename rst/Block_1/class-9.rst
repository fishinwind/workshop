*************************
Class 9 : Python Concepts
*************************

:Class date: Friday 14 February 2014

.. code-block:: bash

    $ cowsay -f dragon-and-cow "Happy Valentine's Day"

Goals
=====
#. Review Python concepts

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

(hint: the :py:meth:`continue` statement will skip to the next loop
iteration, and is usually found in an if statement)

In Class Exercise Code
======================

Print the first word of each line.

.. code-block:: python

    for line in open('hamlet.txt'):
        words = line.strip().split(' ')
        print words[0]

Print only lines that are not indented.

.. code-block:: python

    for line in open('hamlet.txt'):
        if line[0] != ' ':
            print line.strip()

Count the number of times that the word "therefore" appears.

.. code-block:: python

    num = 0
    for line in open('hamlet.txt'):
        words = line.strip().split()
        for word in words:
            if word == 'therefore':
                num += 1
    print num

Counters
========
:py:class:`~collections.Counter` objects are a type of Python dict in
which the values are counts of the keys. Counter objects have several
methods to query the counts, like
:py:meth:`~collections.Counter.most_common()`. 

.. code-block:: python

    from collections import Counter

    word_counts = Counter()

    for line in open('hamlet.txt'):
        words = line.strip().split(' ')

        for word in words:
            word_counts[word] += 1

    print word_counts.most_common(5)
    
There is more than one way to do this. It's usually a good idea to look at the
Python documentation to see if there's a method that does what you're trying to do.

.. code-block:: python

    from collections import Counter

    word_counts = Counter()

    for line in open('hamlet.txt'):
        words = line.strip().split(' ')
        word_counts.update(words)

    print word_counts.most_common(5)

Looping: Reading Multiple Lines at a Time
=========================================
There are lots of biological data files that have information spread over
multiple lines. For example, a FASTA file is used to store sequences. Each
record has a line with '>' and some information (like a name) followed by
another line of sequence data. For example::

    >Sequence name
    AGCATCGTAGCTAGTCGTACGTAGCTATCGATCGTAGCTA

**Download the sample FASTA file:** :ref:`fasta-file`

In Class Exercise
=================

#. Open sample.fa and make a dictionary with four items corresponding to
   the sequences from the file
   
.. code-block:: python

    records  = dict()
    count = 0
    name = ""
    seq = ""

    for line in open('sample.fa'):

        if (count == 0):
            name = line.strip('>\n\r')
            count = 1
        else:
            seq = line.strip()
            count = 0
            records[name] = seq

    print records.items()


Intermediate Concepts: Streaming
================================
One of the reasons why python is so useful is that faciliates
**iteration** over a file without reading the entire dataset into computer
memory.

This is similar to streaming data in the Linux tools we've discussed.
For example:

.. code-block:: bash

    zless /opt/bio-workshop/data/t_R1.fastq.gz | wc -l

never holds the file in memory, it just streams the data.

We can do this in python.

Intermediate : Streaming
========================

.. warning:: 

    DO NOT DO THIS!! It reads everything into memory.

.. code-block:: python

    import gzip
    fastq_filename = '/opt/bio-workshop/data/t_R1.fastq.gz'

    data = list(gzip.open(fastq_filename))
    lines = len(data)

.. important:: 

    DO THIS

.. code-block:: python

    import gzip
    fastq_filename = '/opt/bio-workshop/data/t_R1.fastq.gz'

    lines = 0
    for line in gzip.open(fastq_filename):
        lines += 1

    # or:

    lines = sum(1 for line in gzip.open(fastq_filename))

Streaming with yield
===================================

Make a bed reader that returns a useful dict:

.. code-block:: python

    def bed_generator(bedfilename):
        for line in open(bedfilename):
            if line.startswith('#'): continue
            chrom, start, end, value = line.split("\t")[:4]
            start, end = int(start), int(end)
            yield dict(chrom=chrom, start=start, end=end, value=value)

Then use it:

.. code-block:: python

    bedfilename = '/opt/bio-workshop/data/lamina.bed'
    for bed in bed_generator(bedfilename):
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
    for bed in bed_generator(bedfilename):
        print bed['value']
        vals.append(bed['value'])

    print vals[:10]
    print sum(vals)

Goal
====

Take the basic concepts we've learned and do something useful.

toolshed
========

`toolshed <https://pypi.python.org/pypi/toolshed>`_ is a python module
that simplifies common file/text-processing tasks.  For example, it
assumes the first line of a file is the header and gives a python
dictionary for each line keyed by the header.

.. code-block:: bash

    $ python -c "import toolshed"

    # If you see an error get help to install toolshed:
    $ pip install toolshed

.. code-block:: python

    from toolshed import reader

    bedfilename = '/opt/bio-workshop/data/lamina.bed'

    for region in reader(bedfilename):
        # the first line in lamina.bed is: '#chrom  start  end  value'
        # reader uses these names as keys in a dict

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

    bedfilename = '/opt/bio-workshop/data/lamina.bed'

    for region in reader(bedfilename):
        # do something with region
        print region['chrom']


.. raw:: pdf

    PageBreak
