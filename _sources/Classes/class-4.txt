************************
 Class 4 : grep and awk
************************

:Class date: Thurs 5 Feb 2015

Goals
=====
#. Review

#. remember BED format (chr, start, end)

#. learn grep and awk basics to filter and manipulate text

grep
====
Use :linuxman:`grep(1)` to identify lines in a file that match a specified pattern.

To find any instance of *chr5* in the lamina.bed file

.. code-block:: bash

    # grep [pattern] [filename]
    $ grep chr5 /vol1/opt/data/lamina.bed | head

To find all lines that start with a number sign:

.. code-block:: bash

    # The caret (^) matches the beginning of the line
    # FYI dollar sign ($) matches the end
    $ grep '^#' /vol1/opt/data/lamina.bed

.. nextslide::
    :increment:

To find any line that *does not* start with "chr":

.. code-block:: bash

    # the -v flag inverts the match (grep "not" [pattern])
    $ grep -v '^chr' /vol1/opt/data/lamina.bed

Beware of using ``grep`` to find patterns that might be partial matches:

.. code-block:: bash

    # this will match chr1, chr10, chr11 etc.
    $ grep chr1 /vol1/opt/data/lamina.bed | cut -f1 | uniq

You can find exact matches that are split on words with the ``-w`` flag:

.. code-block:: bash

    # this will only match chr1
    $ grep -w chr1 /vol1/opt/data/lamina.bed | cut -f1 | uniq

.. nextslide::
    :increment:

Beware of using ``grep`` to search for numbers:

.. code-block:: bash

    # finds all strings that match `100`
    $ grep 100 /vol1/opt/data/lamina.bed | head -n 20

    # better, but doesn't look at numeric value
    $ grep -w 100 /vol1/opt/data/lamina.bed | head -n 20

.. tip::

    If you're trying to find numeric values in a file, use ``awk``
    instead::

        $ awk '$2 == 500' /vol1/opt/data/lamina.bed

Exercises
=========

#. use ``grep`` to identify lines in lamina.bed where the second field
   (start) begins with ``100``.

#. use ``grep`` to identify lines in lamina.bed where the third field
   (end) ends with 99 .

#. use ``grep`` with its ``-w`` flag to count the number of 'chr1'
   records in lamina.bed.

#. use ``grep`` to count how many fastq records are in the
   /vol1/opt/data/t_R1.fastq.gz file (fastq records begin with an
   '@' symbol)

#. use ``grep`` to count the number of fastq records in
   /vol1/opt/data/SP1.fq.gz

awk
===
http://en.wikipedia.org/wiki/AWK

``AWK`` is an interpreted **programming language** designed for text
processing and typically used as a data extraction and reporting tool. It
is a standard feature of most Unix-like operating systems.

Named after authors **A** ho, **W** einberger & **K** ernighan

**This is programming**

basic principles
================
#. awk operates on each line of a text file

#. in an awk program, $1 is an alias for the 1st column, $2 for the 2nd, etc. 

#. awk can filter lines by a pattern

awk program structure
=====================

- **BEGIN** runs before the program starts

- **END** runs after the program runs through all lines in the file

- **PATTERN** and **ACTIONS** check and execute on each line.

.. code-block:: bash

    awk 'BEGIN {} (PATTERN) { ACTIONS } END {}' some.file.txt

awk BEGIN
=========

You can use **BEGIN** without a file. We just do one thing then exit:

.. code-block:: bash

   awk 'BEGIN { print 12 * 12 }'

Same with **END**:

.. code-block:: bash

   awk 'END { print 12 * 13 }'
   # then type ctrl+d so it knows it's not getting more input.

filtering
=========
A simple and powerful use of awk is lines that match a pattern or meet set
of criteria. Here, we match (and implicitly print) only lines where the
first column is chr12:

.. code-block:: bash

    awk '($1 == "chr12")' /vol1/opt/data/lamina.bed

We can also filter on start position using '&&' which means 'and':

.. code-block:: bash

    awk '($1 == "chr12" && $2 < 9599990)' /vol1/opt/data/lamina.bed

.. important::

    ``=`` and ``==`` are not the same thing, and are frequently mixed up.

    ``=`` is the assignment operator 
    ``==`` tests for equality 
    ``!=`` tests for inequality.

program structure
=================

.. code-block:: bash

    awk '($1 == "chr12" && $2 < 9599990)' /vol1/opt/data/lamina.bed

.. important::

    + when we are checking as a character ("chr12") we need the quotes.
    + when we are checking as a number (9599990) can not use quotes.
    + can't use commas (e.g. 9,599,990) in numbers

in-class exercise
=================

we will do the first of these together.

#. how many regions (lines) in lamina.bed have a start less than 1,234,567 on any chromosome?

#. how many regions in lamina.bed have a start less than 1,234,567 on chromosome 8?

#. how many regions (lines) in lamina.bed have a start between 50,000 and 951,000

#. how many regions in lamina.bed overlap the interval **chr12:5,000,000-6,000,000** ?

.. important::

    the last question is not trivial and understanding it will be useful

awk program structure (actions)
===============================

print total bases covered on chromosome 13:

.. code-block:: bash

    awk '($1 == "chr13") { coverage = coverage + $3 - $2 }
         END { print coverage }' /vol1/opt/data/lamina.bed

.. important::
    
    1. the entire awk program must be wrapped in quotes. Nearly always best to use
        single quotes (') on the outside.

    2. *coverage* is a variable that stores values; we don't use
        a $ to access it like we do in bash or like we do for the $1,
        $2, ... columns

in-class exercise
=================

below is how we find coverage for chr13. 

.. code-block:: bash

    awk '($1 == "chr13") { coverage += $3 - $2 }
         END{ print coverage }' /vol1/opt/data/lamina.bed

how can we find the total coverage for all chromsomes **except** 13?

awk continued
=============

The ``$0`` variable contains the entire line.

multiple patterns

.. code-block:: bash

      awk '$3 >= 5000 { print $0"\tGREATER" }
           $3  < 5000   { print $0"\tLESS" }' \
            /vol1/opt/data/states.tab

remember we can simply filter to the lines > 5000 with:

.. code-block:: bash

      awk '$3 >= 5000' /vol1/opt/data/states.tab

awk special variables
=====================
#. we know *$1*, *$2*, ... for the column numbers

#. NR is a special variable that holds the line number

#. NF is a special variable that holds the number of fields in the line

#. FS and OFS are the (F)ield and (O)output (F)ield (S)eparators
   i.e. the delimiters (default is any space character)

using awk to count lines with NR
================================

.. code-block:: bash

    $ wc -l /vol1/opt/data/lamina.bed

    $ awk 'END { print NR }' /vol1/opt/data/lamina.bed


using FS and OFS
================
Let's convert lamina.bed to comma-delimited but only for chr12

remember FS is the input separator and OFS is the output delimiter

.. code-block:: bash

    $ awk 'BEGIN{FS="\t"; OFS=","}
        ($1 == "chr12"){ print $1,$2,$3 }' /vol1/opt/data/lamina.bed

regular expressions
===================
we won't cover these in detail, but you can match on *regular expressions*.

The following finds lines containing chr2 (chr2, chr20, chr21) in the first column

.. code-block:: bash

   $ awk '$1 ~ /chr2/' /vol1/opt/data/lamina.bed

Often we can get by without *regular expressions* but they are extremeley powerful
and available in nearly all programming languages.

advanced awk
============
You can do a lot more with awk, here are several resources:

- http://www.hcs.harvard.edu/~dholland/computers/awk.html

- http://doc.infosnel.nl/quickawk.html

- http://www.catonmat.net/download/awk.cheat.sheet.pdf

.. _class-4-exercises:

In Class Exercises - Class 4
============================
we will do the first 2 of these together.

#. use NR to print each line of `lamina.bed` *preceded* by it's line number

#. use NF to see how many columns are in each row of `states.tab`

review
======
+ $1, $2, $3 (default sep is space)
+ adjust sep with: OFS="\t"; FS=","
+ $0 # entire line

.. code-block:: awk

   BEGIN {} 
   (match) { coverage += $3 - $2 } 
   END { print coverage }

+ NR is line number; NF is number of fields;
+ BEGIN {} filter { action } END { }

Exercises
=========

#. are there any regions in `lamina.bed` with start > end?

#. what is the total coverage [sum of (end - start)] of regions on chr13 in `lamina.bed`?

#. what is the mean value (4th column) on chromome 3 of `lamina.bed`

#. print out only the header and the entry for colorado in `states.tab`

#. what is the (single-number) sum of all the incomes for `states.tab`
   with illiteracy rate less than 0.1? greater than 2?

#. use NR to filter out the header from `lamina.bed` (hint: what is NR for the header?)

.. raw:: pdf

    PageBreak
