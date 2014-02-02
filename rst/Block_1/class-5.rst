.. todo::
    - useful ipython directive page for decorator syntax
        http://matplotlib.org/sampledoc/ipython_directive.html

Class 5 : Python : Basics
=========================

Goals
-----
1. Learn to start ipython
2. Learn the basics of python syntax
3. Learn basic types in python

Python overview
---------------
Python is a popular programming language that is commonly used for
bioinformatics. We will use it to process and filter files. When you can't
write a simple script in ``awk``, it is better to use python.

The Python documentation [#]_ is very helpful, with lots of examples. You
should read it to become familiar with the language and refer to it when
you get stuck.

.. [#] Python 2.x docs http://docs.python.org/2/

.. note::

    You should start the Python tutorial [#]_. We will cover some language
    specifics of Python, but will quickly move to using Python  in the
    context of bioinformatic applications.

.. [#] Learn Python the Hard Way
        http://learnpythonthehardway.org/book/

Ipython
-------
Ipython is an (I)nteractive python terminal that lets you
type in python expressions and see the results immediately::

    $ ipython

This command puts you in a shell that accepts python commands, much like
the login terminal accepts ``bash`` commands.

Python indentation
------------------
Many older languages (e.g. PERL, C) use curly brackets to delineate blocks of
code. Python depends on proper indentation of your code. This does not
work:

.. code-block:: python

    for i in (1, 2, 3):
    print i

Instead you have to indent the "print" statement to nest it in the for
loop. 

.. code-block:: python

    for i in (1, 2, 3):
        print i

.. important::

    **Always** use spaces (not tab characters) to indent your code. Change
    your gedit preferences so that when you press the tab key, four spaces are
    inserted. Look in `Edit -> Preferences`, check the `Insert spaces instead
    of tabs` box.

.. note:: 

    IPython will auto-indent code blocks for you. But when you move to
    writing standalone programs, you will need to indent blocks yourself.

For Loops (characters)
----------------------
Lots of things in python are `iterable`, meaning we can write loops
over them. For instance, a string is iterable:

.. ipython::
    :verbatim:

    In [1]: sentence = 'i LOVE programming'

    In [1]: for char in sentence:
       ...:     print char

For Loops (range)
-----------------
We saw an example of a ``for`` loop in the previous example. A ``for``
loop is how you can automaite repetitive tasks.

For example. Print "hello" 5 times:

.. ipython::
    :verbatim:

    In [1]: for i in range(5):
       ...:     print i

where ``range`` is a python function that generates the numbers
`0, 1, 2, 3, 4`.

Try executing the ``range`` function alone at the ``ipython`` prompt.

Python Types
------------
There are several core types in Python that you will use a lot.

- ``Strings`` are collections of characters (words and sentences).
- ``Ints`` and ``Floats`` are numbers.
- ``Lists`` are groups of other objects.
- ``Dictionaries`` contain key:value mappings.

Python Help
-----------
You can find more about any python type or function using ``pydoc``::

    # learn about the python `string` type
    $ pydoc str

At the ``ipython`` prompt, you can also use:

.. ipython::
    :verbatim:

    In [1]: str?

Strings
-------
Strings are collections of characters.

.. ipython::
    :verbatim:

    In [2]: phrase = 'this that other'

    In [3]: phrase 

    # uppercase
    In [3]: phrase.upper()

    # number of characters (including spaces) in phrase
    In [3]: len(phrase)

Lists
-----
Lists are collection of other objects. You can create lists directly,
using brackets (``[ ]``), or they can be created from other objects. Lists
are *subscriptable*, meaning that you can access items in a list by
position.

.. note::

    ``Lists are 0-based and half-open``. This means that:
        - the first item is at index 0
        - you need to specify 

.. ipython::
    :verbatim:

    # convert to list
    In [3]: words = phrase.split()

    # number of items in list
    In [3]: len(words)

    # first item only
    In [3]: phrase[0]

    # add a new word
    In [3]: words.append('foo')

Exercises (1)
-------------
::
    - Try using the ``enumerate`` function on a list.

    - Try using the ``sorted`` and ``reversed`` functions on a list.

Dictionaries
------------

Sets
----
Sets are another type in python that let you store a non-redundant
lists of items. They support logical operations:

.. ipython::
    :verbatim:

    In [11]: skiiers = set(['Tom','Dick','Harry','Gurf'])

    In [12]: snowboarders = set(['Lucy','Steve','Brian','Gurf'])

    # intersection
    In [13]: skiiers & snowboarders

    # union
    In [14]: skiiers | snowboarders

    # difference 
    In [14]: skiiers - snowboarders

Importing modules
-----------------
There are a number of modules with objects and functions in the standard
library, and there are a also a huge number of Python modules on the web
(check github).

To be able to access the contents of a module, you need to import it into
your `namespace`:

.. ipython::

    In [1]: import math

    In [2]: math.log10(1000)


Useful python modules
---------------------
There are several modules in the standard library that we use all the time
for bioinformatics.

    - ``collections``

In Class Exercise
------------------
::
    - stub
