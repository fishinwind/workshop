.. useful ipython directive page for decorator syntax
   http://matplotlib.org/sampledoc/ipython_directive.html

*************************
Class 7 : Python : Basics
*************************

:Class date: Monday 10 February 2014

Goals
=====
#. Learn the basics of python syntax
#. Learn to start ipython
#. Learn basic types in python
#. Write simple python scripts

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

IPython
=======
Ipython is an (I)nteractive python terminal that lets you
type in python expressions and see the results immediately

.. code-block:: bash

    $ ipython

This command puts you in a shell that accepts python commands, much like
the login terminal accepts ``bash`` commands.

Python indentation
==================
Python depends on proper indentation of your code. This does not work:

.. code-block:: python

    for i in (1, 2, 3):
    print i

Instead you have to indent the ``print`` statement to nest it in the for
loop:

.. code-block:: python

    for i in (1, 2, 3):
        print i

.. note::

    Use spaces instead of *tab* characters for indentation.

    Update your gedit preferences: `Edit -> Preferences`,
    check the `Insert spaces instead of tabs` box.

Python Help
===========
You can find more about any python type or function using ``pydoc``:

.. code-block:: bash

    # learn about the python ``string`` type
    $ pydoc str

At the ``ipython`` prompt, you can also use:

.. ipython::
    :verbatim:

    In [1]: str?

In these slides, links will take you to the python docs: :py:obj:`str`

Finally, ask Google (e.g. python string split).

For Loops (iteration)
======================
Many things in python are **iterable**, meaning we can write loops over
them. For example, a string is iterable:

.. ipython::
    :verbatim:

    In [1]: sentence = 'i LOVE programming'

    In [1]: for char in sentence:
       ...:     print char

For Loops (range)
=================
Automate repetitive tasks with a for loop:

.. ipython::
    :verbatim:

    # Print "hello" 5 times:
    In [1]: for i in range(5):
       ...:     print "hello"

    # now print the numbers
    In [1]: for i in range(5):
       ...:     print i

where the :py:func:`range` function generates the numbers `0, 1, 2, 3, 4`.

Python Types
============
There are several core types in Python that you will use a lot.

    - :py:obj:`str` is a collection of characters (words and sentences).
    - :py:obj:`int` and :py:obj:`float` are numbers.
    - :py:obj:`list` is a group of other objects.
    - :py:class:`dict` contains key:value mappings.

Strings
=======
A :py:obj:`str` is a collection of characters. You can make strings with
single, double and triple quotes.

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
Python has integer numbers (:py:obj:`int`) and floating point numbers
(:py:obj:`float`). Math operations work within and across both types:

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
For division you need to pay attention to ``type``:

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

Lists
=====
A :py:obj:`list` is a collection of other objects. You create lists
directly using brackets (``[ ]``), or they can be created from other
objects.

Lists are *subscriptable*, meaning that you can access items in a list by
position.

.. ipython::
    :verbatim:

    # convert to list, str.split() defaults to space
    In [3]: words = phrase.split()

    # number of items in list
    In [3]: len(words)

    # two ways to add new words
    In [3]: words.append('foo')

    In [3]: words.extend(['bar','baz'])

Lists (2)
=========

.. ipython::
    :verbatim:

    # first item only, zero-based
    In [3]: words[0]

    # first through third, start is implicit
    In [3]: words[:3]

    # iterate over the list
    In [7]: for word in words:
       ...:     print word.capitalize()
       ...:     

    # mix types in lists
    In [1]: words.extend([1,2,3])

    # a "list comprehension"
    In [2]: [type(i) for i in words]

Python Exceptions
=================
When you're learning to program in Python, you will see lots of errors.
Examples: :py:class:`~exceptions.ValueError`,
:py:class:`~exceptions.IndexError` and :py:class:`~exceptions.KeyError`

.. ipython::
    :verbatim:

    In [6]: int('blah')

    In [8]: words[100]

    In [7]: parts = {'hip':'thigh'}

    In [9]: parts['nose']

    # Catch errors, print useful debugging messages
    In [12]: try:
       ....:     nums[100]
       ....: except IndexError:
       ....:     print "error: not enough nums"
       ....:     

Reading data from a file
========================
Now we'll read some data from a file and operate on each line:

.. ipython::
    :verbatim:

    In [3]: filename = '/opt/bio-workshop/data/lamina.bed'

    # What is the BUG in this block?
    In [4]: for line in open(filename):
       ...:     fields = line.strip().split('\t')
       ...:     start = fields[1]
       ...:     if start > 5000:
       ...:         print fields 
    
In Class Exercises (1)
======================

    #. Use :py:func:`range` to count from 0 to 100 **by 10**. How do you get
       100 in the result?

    #. Get **every other** value of ``words`` (hint: use a slice)

    #. Use :py:func:`enumerate` on a list (hint: convert the
       result with list(result))

    #. Use :py:func:`sorted` and :py:func:`reversed` on a list.

    #. Do type conversion on each of the fields in the lamina.bed file

Dictionaries (dicts)
====================
A :py:class:`dict` contains key:value mappings. 

.. ipython::
    :verbatim:

    # set up new dicts with {}
    In [14]: produce  = {'lettuce':'green', 'apple':'red',
       ....: 'banana':'yellow'}

    In [5]: produce.keys()

    In [7]: produce.values()

    In [7]: produce.items()

    # sorted by keys
    In [8]: sorted(produce.items())

    # test for membership
    In [9]: 'apple' in produce

    In [10]: not 'orange' in produce

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

Simple Python script
====================
You can also write scripts as files and execute them.
Save this in file called ``run.py``:

.. code-block:: python

    #! /usr/bin/env python

    import sys

    # sys.argv[1] has the filename. the will open the file and read line
    # by line
    for line in open(sys.argv[1]):
        # strip() removes the newline char (\n)
        fields = line.strip().split()

        print fields[0]
        
And run it:

.. code-block:: bash

    $ python run.py /opt/bio-workshop/data/lamina.bed

In Class Exercises (2)
======================

    #. Create a :py:obj:`list` that contains multiple redundant entries.
       Covert the list to a :py:class:`set` with set(list). What happened to
       the redundant entries?

    #. Open lamina.bed and print the start position of each entry

    #. Print the total coverage of entries in lamina.bed

    #. Convert lamina.bed into a :py:obj:`list`. Then, print each entry in the list. 

    #. Find the average value of entries in lamina.bed (watch out for
       int/float issues).

    #. Find the median value of entries in lamina.bed.

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
       argument to look on a particular chromosome. 

    #. Construct regular expressions with :py:mod:`re` to find all states
       in states.tab hose names include 'c' and 'o' in that order, with or
       without other letters in between.

.. raw:: pdf

    PageBreak
