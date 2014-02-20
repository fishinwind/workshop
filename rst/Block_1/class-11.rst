************************
Class 11 : Python Idioms
************************

Goals
=====

 #. Learn some python tricks to write more concise code
 #. Nested data structures

enumerate
=========

We have seen that we often want to know the index (or line number)
that goes with an iterable. For example, in a for loop, if we know
that we are on the first line, we can skip the header.

We can know the index of an iterable with enumerate:

.. code-block:: python

    names = ('fred', 'sally', 'harry', 'jack', 'texan')
    for index, name in enumerate(names):
        print index, name

enumerate
=========

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

Parse a fastq!!

.. code-block:: python

    for i, line in enumerate(open('/opt/bio-workshop/data/SP1.fq')):
        if i % 4 == 0:
            name == line
        elif i % 4 == 1:
            seq == line
        elif i % 4 == 3:
            qual == line
            # now we have name, seq, qual from a single record

list comprehensions
===================

In one problem you had to sum the ord()'s of the quality line.
The common way to do that was this:

.. code-block:: python

    qual_sum = 0
    for q in qual:
        qual_sum += ord(q)

This can be shortened to:

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

parsing fastq (filehandles)
===========================

when you open a file, you get a python file-handle object


.. code-block:: python

    fh = open('/opt/bio-workshop/data/lamina.bed')
    name, seq, plus, qual = fh.readline(), fh.readline(), \
                              fh.readline(), fh.readline()

But how to make that happen continuously?

.. code-block:: python

    from itertools import izip
    for name, seq, plus, qual in izip(fh, fh, fh, fh):
        print name, seq, plus, qual

izip *zips* iterables together and here, we zip for iterables together
that happen to be the same file handle.

Explore zip in ipython by zipping lists of things together.



.. raw:: pdf

    PageBreak
