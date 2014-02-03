.. todo::
    - useful ipython directive page for decorator syntax
        http://matplotlib.org/sampledoc/ipython_directive.html

*************************
Class 5 : Python : Basics
*************************

Goals
=====
#. Learn the basics of python syntax
#. Learn to start ipython
#. Learn basic types in python

Overview
========
Python is a popular programming language that is commonly used for
bioinformatics. We will use it to process and filter files. When you can't
write a simple script in ``awk``, it is better to use python.

The Python documentation [#]_ is very helpful, with lots of examples. You
should read it to become familiar with the language and refer to it when
you get stuck.

.. [#] Python 2.x docs http://docs.python.org/2/

.. important::

    **You should begin the Python tutorial** [#]_. We will cover some language
    specifics of Python, but will quickly move to using Python  in the
    context of bioinformatic applications.

.. [#] Learn Python the Hard Way
        http://learnpythonthehardway.org/book/

Ipython
=======
Ipython is an (I)nteractive python terminal that lets you
type in python expressions and see the results immediately

.. code-block:: bash

    $ ipython

This command puts you in a shell that accepts python commands, much like
the login terminal accepts ``bash`` commands.

Python indentation
==================
Many older languages (e.g. bash, PERL, C) use curly brackets to delineate
blocks of code. In contrast, Python depends on proper indentation of your
code. This does not work:

.. code-block:: python

    for i in (1, 2, 3):
    print i

Instead you have to indent the "print" statement to nest it in the for
loop:

.. code-block:: python

    for i in (1, 2, 3):
        print i

Python indentation(2)
=====================

.. important::

    **Always** use spaces – and not *tab* characters – to indent
    your code. Change your gedit preferences so that when you press the
    tab key, four spaces are inserted. Look in `Edit -> Preferences`,
    check the `Insert spaces instead of tabs` box.

    IPython will auto-indent code blocks for you. But when you move to
    writing standalone programs, you will need to indent blocks yourself.

Python Help
===========
You can find more about any python type or function using ``pydoc``::

    # learn about the python `string` type
    $ pydoc str

At the ``ipython`` prompt, you can also use:

.. ipython::
    :verbatim:

    In [1]: str?

Finally, ask Google (e.g. python list slice).

For Loops (characters)
======================
Lots of things in python are `iterable`, meaning we can write loops
over them. For instance, a string is iterable:

.. ipython::
    :verbatim:

    In [1]: sentence = 'i LOVE programming'

    In [1]: for char in sentence:
       ...:     print char

For Loops (range)
=================
You can also automate repetitive tasks with a for loop:

.. ipython::
    :verbatim:

    # Print "hello" 5 times:
    In [1]: for i in range(5):
       ...:     print "hello"

    # now print the numbers
    In [1]: for i in range(5):
       ...:     print i

where :py:func:`range` is a python function that generates the numbers
`0, 1, 2, 3, 4`.

Python Types
============
There are several core types in Python that you will use a lot.

    - :py:obj:`str` is a collection of characters (words and sentences).
    - :py:obj:`int` and :py:obj:`float` are numbers.
    - :py:obj:`list` is a group of other objects.
    - :py:class:`dict` contains key:value mappings.

Strings
=======
Strings are collections of characters.

.. ipython::
    :verbatim:

    In [2]: phrase = 'this that other'

    In [3]: phrase 

    # uppercase
    In [3]: phrase.upper()

    # number of characters (including spaces) in phrase
    In [3]: len(phrase)

Numbers (Ints and math)
=========================
Python has an integer number representation (:py:obj:`int`) and a floating point
representation (:py:obj:`float`). Most math operations work within and across
both types:

.. ipython::
    :verbatim:

    # set up some ints
    In [6]: x = 10

    In [7]: y = 100

    In [8]: type(x)

    # add
    In [9]: x + y

    # subtract
    In [10]: x - y

    # x * y
    In [11]: x * y

Numbers (Float division)
========================
Division is a case where you need to pay attention to ``type``:

.. ipython::
    :verbatim:

    # try to divide the ints ...
    In [12]: x / y

    # need float conversion!
    In [14]: float(x) / float(y)

    # make floats directly and divide
    In [15]: x = 10.0

    In [16]: y = 100.0

    In [16]: type(x)

    In [17]: x / y

.. note:: This changed in Python 3, where 5 / 2 will return 2.5. If you
    want that behaviour, you need to add this to your code::

        from __future__ import division
    
Lists
=====
A :py:obj:`list` is a collection of other objects. You can create lists
directly using brackets (``[ ]``), or they can be created from other
objects.

Lists are *subscriptable*, meaning that you can access items in a list by
position.

.. ipython::
    :verbatim:

    In [2]: phrase = 'this that other'

    # convert to list
    In [3]: words = phrase.split()

    # number of items in list
    In [3]: len(words)

Lists (2)
=========

.. ipython::
    :verbatim:

    # two ways to add new words
    In [3]: words.append('foo')

    In [3]: words.extend(['bar','baz'])

    # first item only, zero-based
    In [3]: words[0]

    # first through third
    In [3]: words[:3]

In Class Exercises (1)
======================
Here are a few exercises:

    #. Use :py:func:`range` to count from 0 to 100 by 10.  How do you get 100 in the
      result?

    #. Get every other value of ``words`` (hint: use a slice)

    #. Use :py:func:`enumerate` on a list (hint: convert the
       result with list(result))

    #. Use :py:func:`sorted` and :py:func:`reversed` on a list.

Dictionaries (dicts)
====================
A :py:class:`dict` contains key:value mappings. 

.. ipython::
    :verbatim:

    # set up new dicts with {}
    In [3]: produce  = {'apple':'red', 'banana':'yellow', 'lettuce':'green'}

    In [5]: produce.keys()

    In [7]: produce.values()

    # sorted by keys
    In [8]: sorted(produce.items())

    # test for membership
    In [9]: 'apple' in produce

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

Useful python modules
=====================
There are several modules in the standard library that we use all the time
for bioinformatics.

    - :py:mod:`collections`: espcially :py:class:`collections.defaultdict`
      and :py:class:`collections.Counter`
    - :py:mod:`itertools`: tools for efficient aggregation and iteration

In Class Exercises (2)
======================
Here are a few exercises:

    #. Create a :py:obj:`dict` that contains several key:value pairs. 

    #. Create a :py:obj:`list` that contains multiple redundant entries.
       Covert the list to a :py:class:`set` with set(list). What happened to
       the redundant entries?

