*****************************
Class 8 : Python : Basics (2)
*****************************

Goals
=====
#. foo
#. NOTE for Jay to cover: slicing lists [:5], [2:4] (like linux cut, but 0-based)

Sets
====
A :py:class:`set` is another type in python that let you store a non-redundant
lists of items. They support logical operations:

.. ipython::

    In [11]: skiiers = set(['Tom','Dick','Harry','Gurf'])

    In [12]: snowboarders = set(['Lucy','Steve','Brian','Gurf'])

    # intersection
    In [13]: skiiers & snowboarders

    # union
    In [14]: skiiers | snowboarders

    # difference 
    In [14]: skiiers - snowboarders

Regular Expressions
===================
Python provides a regular expression module for pattern matching. We'll
cover some basics of writing regular expressions:

.. ipython::
    :verbatim:

    In [1]: phrase = 'how now brown cow'

    In [2]: import re

    In [3]: regex = re.compile('brown')

    In [6]: regex.findall(phrase) 

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

.. raw:: pdf

    PageBreak
