************************
Class 11 : Python Idioms
************************

Goals
=====

 #. Learn some python tricks to write more concise code
 #. Nested data structures

FASTQ parsing
=============

We have seen the problem of parsing a FASTQ file.
How do we describe the format in English?

   FASTQ records occur in groups of 4 in the order: name, seq, plus, qual

How can we do this in python?


enumerate (1)
=============

We have seen that we often want to know the index (or line number)
that goes with an iterable. For example, in a for loop, if we know
that we are on the first line, we can skip the header.

We can know the index of an iterable with enumerate:

.. code-block:: python

    names = ('fred', 'sally', 'harry', 'jack', 'texan')
    for index, name in enumerate(names):
        print index, name

enumerate (2)
=============

So we can wrap any iterable in enumerate and we will get a tuple of
`index, thing`. Where `thing` was the element of the list.

We can skip the header in a file like this:

.. code-block:: python

    for i, line in enumerate(open('/opt/bio-workshop/data/lamina.bed')):
        # skip the header
        if i == 0: continue
        fields = line.rstrip().split("\t")
        # ... do stuff with fields


modulo
======

Modulo is the remainder operation.

So 12 modulo 4 is 0. 13 modulo 4 is 1

.. ipython::

    In [1]: 12 % 4
    Out[1]: 0

    In [2]: 13 % 4
    Out[2]: 1

modulo and enumerate
====================

.. ipython::

    In [1]: for i in range(12):
       ...:     print i, i % 4
       ...:     
    0 0
    1 1
    2 2
    3 3
    4 0
    5 1
    6 2
    7 3
    8 0
    9 1
    10 2
    11 3

How does this relate to our FASTQ?

modulo, enumerate, fastq
========================

.. ipython

    In [1]: for i, line in enumerate(open('misc/data/SP1.fq')):
       ...:     print i, i % 4, line.strip()
       ...:     if i > 8: break
       ...:     
    0 0 @cluster_2:UMI_ATTCCG
    1 1 TTTCCGGGGCACATAATCTTCAGCCGGGCGC
    2 2 +
    3 3 9C;=;=<9@4868>9:67AA<9>65<=>591
    4 0 @cluster_8:UMI_CTTTGA
    5 1 TATCCTTGCAATACTCTCCGAACGGGAGAGC
    6 2 +
    7 3 1/04.72,(003,-2-22+00-12./.-.4-
    8 0 @cluster_12:UMI_GGTCAA
    9 1 GCAGTTTAAGATCATTTTATTGAAGAGCAAG


modulo, enumerate, fastq: parse
===============================

Parse a fastq!!

.. code-block:: python

    for i, line in enumerate(open('/opt/bio-workshop/data/SP1.fq')):
        if i % 4 == 0:
            name == line
        elif i % 4 == 1:
            seq == line
        elif i % 4 == 3:
            qual == line
            # here have name, seq, qual from a single record

zip
===

zip is another python function. It merges items from multiple lists:

.. ipython:: 

    In [2]: a = range(5)

    In [3]: b = "abcde"

    In [4]: zip(a, b)
    Out[4]: [(0, 'a'), (1, 'b'), (2, 'c'), (3, 'd'), (4, 'e')]

    In [5]: c = [dict(), [], None, "hello", "world"]

    In [6]: zip(a, b, c)
    Out[6]: 
    [(0, 'a', {}),
     (1, 'b', []),
     (2, 'c', None),
     (3, 'd', 'hello'),
     (4, 'e', 'world')]

    
izip
====

izip is a lazy version of zip. It doesn't consume or return elements until you
 ask for them.

.. ipython::

    In [10]: from itertools import izip

    In [11]: izip(a, b, c)
    Out[11]: <itertools.izip at 0x2799d88>

    In [12]: for item_a, item_b, item_c in izip(a, b, c):
       ....:     print item_a, item_b, item_c
       ....:     
    0 a {}
    1 b []
    2 c None
    3 d hello
    4 e world



list comprehensions(1)
======================

In one problem you had to sum the ord()'s of the quality line.
The common way to do that was this:

.. code-block:: python

    qual_sum = 0
    for q in qual:
        qual_sum += ord(q)

Once could get the quals instead as:

.. code-block:: python

    integer_quals = [ord(q) for q in qual]

So the sum can be shortened to:

.. code-block:: python

    qual_sum = sum(ord(q) for q in qual)

parsing fastq
=============

what if we could get:

.. code-block:: python

    for name, seq, plus, qual in ????:
        # do stuff

Then we could use enumerate to count records:

.. code-block:: python

    for rec_no, (name, seq, plus, qual) in enumerate(????):
        if rec_no == 10: break
        # do stuff

Remember how we saw that zip and izip could merge iterables?

parsing fastq (filehandles)
===========================

when you open a file, you get a python file-handle object


.. code-block:: python

    fh = open('/opt/bio-workshop/data/lamina.bed')
    name, seq, plus, qual = fh.readline(), fh.readline(), \
                              fh.readline(), fh.readline()

But how to make that happen continuously?

.. code-block:: python

    fh = open('/opt/bio-workshop/data/lamina.bed')
    from itertools import izip
    for name, seq, plus, qual in izip(fh, fh, fh, fh):
        print (name, seq, plus, qual)
        print # add a new-line

izip *zips* iterables together and here, we zip for iterables together
that happen to be the same file handle.

.. raw:: pdf

    PageBreak
