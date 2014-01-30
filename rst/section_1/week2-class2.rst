Week 2 / Class 2: Python Files
==============================

Goals
-----

1. start with python
2. ipython
3. python programming

python
------

Python is a popular programming language that is commonly used for bioinformatics.

We will use it to process and filter files. When you can't write a simple script
in `awk` it is better to use python.


Note
----

.. note:
    
    Python depends on the alignment of your code to understand.  

.. code-block:: python

    for i in (1, 2, 3):
    print i

will not work. While:


.. code-block:: python

    for i in (1, 2, 3):
        print i

Will work. 

Please always use spaces (not tabs) to indent your code.

Ipython
-------

Ipython is an (I)nteractive python terminal that let's one 
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

where `range` is a python construct that generates the numbers
`0, 1, 2, 3, 4`.


For Loops (characters)
----------------------

Lots of things in python are `iterable` meaning we can write loops
over them. For instance, a string is iterable:

.. code-block:: python

    for char in "bioinformatics programming":
        print char

Try this in `ipython` to see what happens.



awk continued
-------------

$0 contains the entire line.

multiple patterns:

    $ awk '($4 < 0.05) { print $0"\tsignificant"}($4 >= 0.05) { print $0"\tlame" }' input.bed > output.classified.bed
   

bioawk
------

Bioawk is a variant of awk that knows about common sequence formats. To count
the number of records in a fastq file:

    $ bioawk -c fastx 'END { print $NR }'

biowak (names)
--------------

You can access `name`, `seq`, `qual`, `comment`:

   $ biowak -c fastx '{ print $name, $seq, $qual}'


See the README at: https://github.com/lh3/bioawk
for more info
