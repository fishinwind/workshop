
.. include:: /_static/substitutions.txt

************************
Class 19 : Python Idioms
************************

:Class date: |c19-date|
:Last updated: |today|

Goals
=====

#. python tricks to write more concise code
#. Laziness
#. Data structures

Variable Names Vs. Builtins
===========================

you decide what to call a variable in a for loop or the name of
a list or dictionary.

There are python keywords and builtin functions

list of keywords:

.. code-block:: python

    >>> import keyword
    >>> print " ".join(keyword.kwlist)
    and as assert break class continue def del elif else except exec
    finally for from global if import in is lambda not or pass print
    raise return try while with yield

list of builtin functions:

    http://docs.python.org/2/library/functions.html

(includes range, open, file, list, dict, etc)


FASTQ parsing
=============

We have seen the problem of parsing a FASTQ file.
How do we describe the format in English?

   FASTQ records occur in groups of 4 in the order: name, seq, plus, qual

We will learn some tools to simplify this.

enumerate
=========

We have seen that we often want to know the index (or line number)
that goes with an iterable. For example, if we know
that we are on the first line, we can skip the header.

We can know the index of an iterable with enumerate:

.. code-block:: python

    names = ('fred', 'sally', 'harry', 'jack', 'texan')

    for index, name in enumerate(names):
        print index, name

    # To illustrate the freedom you have in choosing variable names
    # the above is identical to
    for really_bad_var_name, x42 in enumerate(names):
        print really_bad_var_name, x42

.. nextslide::
    :increment:

Try this in the ipython terminal:

.. code-block:: python

    enumerate("abcdefg")

Enumerate is lazy, meaning it won't consume an iterable until we ask it to

.. code-block:: python

    list(enumerate("abcdefg"))
    # OR 
    for i, letter in enumerate("abcdefg"):
        print i, letter
    # remember we have our choice of variable names. The above is identical to
    for xx, hello in enumerate("abcdefg"):
        print xx, hello

Generally we name index variables as *i* and give the other variables names that
make sense.

.. nextslide::
    :increment:

When we wrap any iterable in enumerate and we get a tuple of
`index, thing`. Where `thing` was the element of the original list.

We can skip the header in a file like this:

.. code-block:: python

    for i, line in enumerate(open('/vol1/opt/data/lamina.bed')):
        # skip the header
        if i == 0: continue
        fields = line.rstrip().split("\t")
        # or we can get the variables directly since
        # we know there are 4 cols
        chrom, start, end, val = line.rstrip().split("\t")


Using enumerate like this is safer than manually incrementing a variable
as sometimes you will forget to increment or you will *continue* before
incrementing.

modulo
======

Modulo is the remainder operation.

+ 12 modulo 4 is 0
+ 13 modulo 4 is 1

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

.. ipython::

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

    for i, line in enumerate(open('/vol1/opt/data/SP1.fq')):
        if i % 4 == 0:
            name = line
        elif i % 4 == 1:
            seq = line
        elif i % 4 == 3:
            qual = line
            # here have name, seq, qual from a single record

note how this fairly closely matches our english explanation of the fastq
format.

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
    Out[6]: [(0, 'a', {}),
     (1, 'b', []),
     (2, 'c', None),
     (3, 'd', 'hello'),
     (4, 'e', 'world')]

    
izip
====

 izip is a lazy version of zip. It doesn't consume or return elements until you
 ask for them.

.. ipython::

    In [1]: from itertools import izip

    In [2]: seq = "TTTCCGGGGCACATAATCTTCAGCCGGGCGC"

    In [3]: qual = "9C;=;=<9@4868>9:67AA<9>65<=>591"

    In [4]: izip(seq, qual)
    Out[4]: <itertools.izip at 0x2419368>

    In [5]: for base, base_qual in izip(seq, qual):
       ...:     print base, base_qual
    T 9
    T C
    T ;
    C =
    ...


izip laziness
=============

Laziness is important, if for example we are zipping over a file. If we use
**zip** it will consume the entire file immediately and read it into memory.
**izip** will only consume the file as we request the zipped items.

Note that in the previous slide, we associated each base with it's base-quality.
That's useful...

list comprehensions
===================

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

Why use functions?
==================
Functions are useful for encapsulating reusable chunks of code. For
exmaple, you don't want to write messy code for parsing a fastq file every
time you need to parse a fastq file. Instead, you define a function:

.. code-block:: python
    
    def parse_fastq(filename):
        # parse records
   
Once that is defined, you can put it in a file in your PYTHONPATH called
``mytools.py`` and use it:

.. code-block:: python

    # look ma! no messy parsing code!
    from mytools import parse_fastq

    for record in parse_fastq(filename):
        # use the record
 
in-class exercise
=================
calculate mean base-quality by base.

zip quality with sequence. append quality for each base in a dict of lists

.. code-block:: python

    # append all quality scores for A base to quals_by_base['A'] list.
    quals_by_base = {'A': [], 'C': [], 'T': [], 'G': []}
    for i, line in enumerate(open('/vol1/opt/data/SP1.fq')):
        if i % 4 == 0:
            name = line
        elif i % 4 == 1:
            seq = line
        elif i % 4 == 3:
            qual = line
            # update quals_by_base here since we have seq and qual
            # use zip/izip
            ...
    # outside the loop calculate the avg base quality:
    for base, integer_quals in quals_by_base.items():
        mean_quals = XXX_FIX_ME_XXX # remember to float()
        print base, mean_quals

exercises
=========

+ do previous exercise without a list. instead storing running sum and count of
  quals and using that at the end.
+ look at xrange, the lazy version of range
+ how can you implement your own version of enumerate using izip and xrange?
+ clean up some of your homeworks using the simpler fastq parsing.
+ look at the itertools module (http://docs.python.org/2/library/itertools.html)

Resources
=========

+ idiomatic python: http://python.net/~goodger/projects/pycon/2007/idiomatic/handout.html
+ itertools: http://naiquevin.github.io/a-look-at-some-of-pythons-useful-itertools.html

.. raw:: pdf

    PageBreak
