*****************************
Class 5 : Python : Basics (2)
*****************************

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

In Class Exercise
=================
