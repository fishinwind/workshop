Week 2 / Class 2 : Python Files
===============================

Goals
-----
1. start with python
2. ipython
3. python programming

Python overview
---------------
Python is a popular programming language that is commonly used for bioinformatics.

We will use it to process and filter files. When you can't write a simple script
in `awk` it is better to use python.

The Python documentation is very helpful, with lots of examples. You should
read it to become familiar with the language: http://docs.python.org/2/

.. note::
    We will cover all of the Python programming in Python 2.x. We will *not* be using Python 3.

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

Ipython
-------
Ipython is an (I)nteractive python terminal that lets you
type in python expressions and see the results immediately.

Start it by typing `ipython` in your terminal.

Once you do so, you are in shell that accepts python commands
much like the normal terminal accepts `bash` commands.


For Loops (range)
-----------------
We have seen an example of a `for` loop in the previous
example. A `for` loop is how you can automaite repetitive
tasks.

For example. Print "hello" 5 times:

.. code-block:: python

    for i in range(5):
        print "hello"

where `range` is a python function that generates the numbers
`0, 1, 2, 3, 4`. Try executing the `range` function alone at the ipython
prompt.


For Loops (characters)
----------------------
Lots of things in python are `iterable`, meaning we can write loops
over them. For instance, a string is iterable:

.. code-block:: python

    for char in "bioinformatics programming":
        print char

Try this in `ipython` to see what happens.

Cluster access
--------------
We have set up accounts for the class on our departmental cluster. We will
set up your accounts at the end of class and reset your passwords.

.. code-block:: bash
    
    $ ssh -X username@amc-tesla.ucdenver.pvt

