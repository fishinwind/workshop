.. useful ipython directive page for decorator syntax
   http://matplotlib.org/sampledoc/ipython_directive.html

*************************
Class 15 : Python : Basics
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

.. nextslide::
    :increment::

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

    In [3]: filename = '/vol1/opt/data/lamina.bed'

    # What is the BUG in this block?
    In [4]: for line in open(filename):
       ...:     fields = line.strip().split('\t')
       ...:     start = fields[1]
       ...:     if start > 5000:
       ...:         print fields 
    
In Class Exercises
==================

#. Use :py:func:`range` to count from 0 to 100 **by 10**. How do you get
   100 in the result?

#. Get **every other** value of ``words`` (hint: use a slice)

#. Use :py:func:`enumerate` on a list (hint: convert the
   result with list(result))

#. Use :py:func:`sorted` and :py:func:`reversed` on a list.

#. Do type conversion on each of the fields in the lamina.bed file

.. raw:: pdf

    PageBreak
