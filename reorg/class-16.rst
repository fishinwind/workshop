*****************************
Class 16 : Python : Basics (2)
*****************************

:Class date: Wednesday 12 February 2014

Goals
=====
#. Finish concepts from Monday (list slices)
#. New concepts: equality, undefined, functions
#. Importing and debugging
#. Exercises

Review of list concepts
=======================
We covered slicing of lists; i.e. retrieving list parts:

.. ipython::
    :verbatim:

    In [11]: nums = range(30)

    # try this again with implicit start and end
    In [12]: nums[5:10]

    In [2]: reversed(nums)

    # two ways to examine the contents of the iterator,
    # same principle for sorted()
    In [3]: [i for i in reversed(nums)]

    In [3]: list(reversed(nums))

    # skip reversed()
    In [5]: range(30,0,-1)

    # edit a file directly from ipython
    In [4]: edit filename.py

Sets
====
A :py:class:`set` is another type in python that lets you store a non-redundant
list of items. They support logical operations:

.. ipython::

    In [11]: skiiers = set(['Tom','Harry','Gurf'])

    In [12]: snowboarders = set(['Lucy','Brian','Gurf'])

    # intersection
    In [13]: skiiers & snowboarders

    # union
    In [14]: skiiers | snowboarders

    # difference 
    In [14]: skiiers - snowboarders

Equality and Logic
==================
Use ``if``:``elif``:``else`` statements to test conditions and act on the
result. The ``==`` and ``!=`` operators test for equality and inequality, and
work on many object comparisons.

.. ipython::
    :verbatim:

    # animal colors
    In [3]: cat = 'white'

    In [4]: dog = 'black'

    In [5]: if cat == dog: 
       ...:     print "same color"
       ...: elif cat != dog:
       ...:     print "different color"
       ...: else:
       ...:     print "not going to happen"
       ...:     

Undefined values in Python 
==========================
.. ipython::
    :verbatim:

    In [1]: this = None

    In [4]: bool(this)

    In [2]: not this

    # the following if statements are equivalent:
    In [6]: if this is None:
       ...:     print 'foo'
       ...:     

    In [5]: if not this:
       ...:     print 'foo'
       ...:     

    # set the following and test with ``not this``
    In [7]: this = 0

    In [9]: this = ''

Defining functions
==================
You can define functions that encapsulate work flows

.. ipython::
    :verbatim:

    In [19]: def square(numlist):
       ....:     result = []
       ....:     for i in numlist:
       ....:         sq = i * i 
       ....:         result.append(sq)
       ....:         return result

    # or replace loop with: return [i*i for i in numlist]

    In [20]: square(nums)

Importing modules
=================
There are a number of modules with objects and functions in the standard
library, and there are a also a huge number of Python modules on the web
(check github).

To be able to access the contents of a module, you need to import it into
your `namespace`:

.. ipython::

    In [1]: import math

    In [2]: math.log10(1000)

    In [3]: import sys

.. Regular Expressions
.. ===================
.. Python provides a regular expression module for pattern matching. We'll
.. cover some basics of writing regular expressions:

.. .. ipython::
..     :verbatim:

..    In [1]: phrase = 'how now brown cow'

..    In [2]: import re

..    In [3]: regex = re.compile('brown')

..    In [6]: regex.findall(phrase) 

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

In Class Exercises
==================

#. Create a :py:obj:`list` that contains multiple redundant entries.
   Covert the list to a :py:class:`set` with set(list). What happened to
   the redundant entries?

#. Open lamina.bed and print the start position of each entry

#. Print the total coverage of entries in lamina.bed

#. Convert each row in lamina.bed into a :py:obj:`list`. Then, print
   each entry in the list. 

#. Find the average value of entries in lamina.bed (watch out for
   int / float type issues).

#. Find the median value of entries in lamina.bed. Then find the mode.

Out of Class Exercises 
======================

#. Use a python :py:class:`dict` object to count the number of entries
   on each chromosome in lamina.bed. 

#. Do the same thing as the previous exercise, but using a
   :py:class:`~collections.Counter` object. Then, use the Counter()
   methods to find out which chromosomes have the largest and smallest
   number of entries. 

#. Create a python script that takes a chromosome number, and finds all entries
   in lamina.bed that are on that chromosome. 

#. Modify the previous script to use :py:mod:`argparse`, so that it
   will find entries on every chromosome by default unless given an
   argument to look on a particular chromosome (advanced)

.. raw:: pdf

    PageBreak
