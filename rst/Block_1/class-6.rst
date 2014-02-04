*****************************
Class 5 : Python : Basics (2)
*****************************

Goals
=====
#. foo

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

    - :py:mod:`collections`: espcially :py:class:`collections.defaultdict`
      and :py:class:`collections.Counter`

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
