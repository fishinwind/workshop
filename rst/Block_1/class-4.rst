Class 4 : awk
=============

Monday: 3 February 2014

Goals
-----

 1. learn awk basics to filter and manipulate text


awk
---

http://en.wikipedia.org/wiki/AWK

``AWK`` is an interpreted programming language designed for text processing
and typically used as a data extraction and reporting tool. It is a
standard feature of most Unix-like operating systems.

Named after authors **A** ho, **W** einberger & **K** ernighan

basic program structure
-----------------------

actions on each line that match a pattern

.. code-block:: bash

    awk 'PATTERN { ACTIONS }'

columns available as *variables* $1, $2, ... $n

program structure example
-------------------------

extract p-values (in 4th column from BED file) that are < 0.04


.. code-block:: bash

    awk '$4 < 0.04' all.pvalues.bed > some.pvalues.bed

limit to chr12 (&& means "and")

.. code-block:: bash

    $ awk '$4 < 0.04 && $1 == "chr12"' \
        all.pvalues.bed \
        > chr12.pvalues.bed

awk continued
-------------

The ``$0`` variable contains the entire line.

multiple patterns

.. code-block:: bash

    $ awk '($4 < 0.05) { print $0"\tsignificant"}($4 >= 0.05) \
        { print $0"\tlame" }' \
        input.bed \
        > output.classified.bed
   
bioawk
------

Bioawk [#]_ is a variant of awk that knows about common sequence formats. To
count the number of records in a fastq file

.. code-block:: bash

    $ bioawk -c fastx 'END { print $NR }'

.. [#] Bioawk README https://github.com/lh3/bioawk

biowak (names)
--------------

You can access `name`, `seq`, `qual`, `comment`

.. code-block:: bash

    $ biowak -c fastx '{ print $name, $seq, $qual}'

