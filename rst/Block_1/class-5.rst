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
bioinformatics. 

We will use it to process and filter files. When you can't write a simple
script in ``awk``, it is better to use python.

The Python documentation [#]_ is very helpful, with lots of examples. You
should read it to become familiar with the language.

.. [#] Python 2.x docs http://docs.python.org/2/

.. note::

    You should be going through the Python tutorial at this point. We will
    cover some language specifics of Python in the context of
    bioinformatic applications.

Ipython
-------
Ipython is an (I)nteractive python terminal that lets you
type in python expressions and see the results immediately::

    $ ipython

This command puts you in a shell that accepts python commands, much like
the login terminal accepts `bash` commands.

Formatting
----------
Python depends on the alignment of your code to understand. This does not
work:

.. code-block:: python

    for i in (1, 2, 3):
    print i

Instead you have to indent the "print" statement to nest it in the for
loop. 

.. code-block:: python

    for i in (1, 2, 3):
        print i

**Always** use spaces (not tab characters) to indent your code. We will
set up your gedit preferences so that when you press the tab key, four
spaces are inserted.

.. note:: 

    IPython will auto-indent code blocks for you. But when you move to
    writing standalone programs, you will need to indent blocks yourself.

For Loops (range)
-----------------
We have seen an example of a ``for`` loop in the previous
example. A ``for`` loop is how you can automaite repetitive
tasks.

For example. Print "hello" 5 times:

.. ipython::
    :verbatim:

    In [1]: for i in range(5):
       ...:     print i

where ``range`` is a python function that generates the numbers
`0, 1, 2, 3, 4`. Try executing the ``range`` function alone at the ipython
prompt.

Try this yourself in ``ipython``.

For Loops (characters)
----------------------
Lots of things in python are `iterable`, meaning we can write loops
over them. For instance, a string is iterable:

.. ipython::
    :verbatim:

    In [1]: for char in 'i LOVE programming':
       ...:     print char

Python Types
------------
There are several core types in Python that you will use a lot.

- ``Strings`` are collections of characters (words and sentences).
- ``Ints`` and ``Floats`` are numbers.
- ``Lists`` are groups of other objects.
- ``Dictionaries`` contain key:value mappings.

Strings
-------
Strings are collections of characters.

.. ipython::
    :verbatim:

    In [2]: words = 'this that other'

    In [3]: words

    In [3]: words.upper()

    # convert to list
    In [3]: words.split()

    # convert to list
    In [4]: words.split().capitalize()

Everything in Python is an object. In practice this means that there is an
expected presentation of everything, but everything has additional methods
that can be called.

.. code-block:: python


In Class Exercise
------------------
::

