*************
Class 4 : awk
*************

:Class date: Monday 3 February 2014

Goals
=====
#. Review
#. remember BED format (chr, start, end)
#. learn awk basics to filter and manipulate text

Review
======

 what we've learned

awk
---
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

 + **BEGIN** runs before the program starts
 + **END** runs after the program runs through all lines in the file
 + **PATTERN** and **ACTIONS** check and execute on each line.

.. code-block:: bash

    awk 'BEGIN {} (PATTERN) { ACTIONS } END {}' some.file.txt

awk BEGIN
=========

Using begin, we don't need a file. We just do one thing then exit:

.. code-block:: bash

   awk 'BEGIN { print 12 * 12 }'

Same with **END**:

.. code-block:: bash

   awk 'END { print 12 * 13 }' # then type ctrl+d so it knows it's not getting more input.
 
filtering
=========
A simple and powerful use of awk is lines that match a pattern or meet set
of criteria. Here, we match (and implicitly print) only lines where the
first column is chr12:

.. code-block:: bash

    awk '($1 == "chr12")' /opt/bio-workshop/data/lamina.bed

We can also filter on start position using '&&' which means 'and':

.. code-block:: bash

    awk '($1 == "chr12" && $2 < 9599990)' /opt/bio-workshop/data/lamina.bed

.. important::

    ``=`` and ``==`` are not the same thing, and are frequently mixed up.

    ``=`` is the assignment operator 
    ``==`` tests for equality 
    ``!=`` tests for inequality.

program structure
=================
.. code-block:: bash

    awk '($1 == "chr12" && $2 < 9599990)' /opt/bio-workshop/data/lamina.bed

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
         END{ print coverage }' /opt/bio-workshop/data/lamina.bed

.. important::
    
 #. the entire awk program must be wrapped in quotes. Nearly always best to use
    single quotes (') on the outside.
 #. *coverage* is a variable that stores values; we don't use
    a $ to access it like we do in bash or like we do for the $1,
    $2, ... columns


in-class exercise
=================

below is how we find coverage for chr13. 

.. code-block:: bash

    awk '($1 == "chr13") { coverage += $3 - $2 }
         END{ print coverage }' /opt/bio-workshop/data/lamina.bed

how can we find the total coverage for all chromsomes **except** 13?

awk continued
=============

The ``$0`` variable contains the entire line.

multiple patterns

.. code-block:: bash

      awk '$3 >= 5000 { print $0"\tGREATER" }
           $3  < 5000   { print $0"\tLESS" }' \
            /opt/bio-workshop/data/states.tab

remember we can simply filter to the lines > 5000 with:

.. code-block:: bash

      awk '$3 >= 5000' /opt/bio-workshop/data/states.tab

awk special variables
=====================
 #. we know *$1*, *$2*, ... for the column numbers
 #. NR is a special variable that holds the line number
 #. NF is a special variable that holds the number of fields in the line

 #. FS and OFS are the (F)ield and (O)output (F)ield (S)eparators
    --meaning the delimiters (default is any space character)

using awk to count lines with NR
================================
.. code-block:: bash

    $ wc -l /opt/bio-workshop/data/lamina.bed

    $ awk 'END { print NR }' /opt/bio-workshop/data/lamina.bed


using FS and OFS
================
Let's convert lamina.bed to comma-delimited but only for chr12

remember FS is the input separator and OFS is the output delimiter

.. code-block:: bash

    $ awk 'BEGIN{FS="\t"; OFS=","}
        ($1 == "chr12"){ print $1,$2,$3 }' /opt/bio-workshop/data/lamina.bed

regular expressions
===================
we won't cover these in detail, but you can match on *regular expressions*.

The following finds lines containing chr2 (chr2, chr20, chr21) in the first column

.. code-block:: bash

   $ awk '$1 ~ /chr2/' /opt/bio-workshop/data/lamina.bed

Often we can get by without *regular expressions* but they are extremeley powerful
and available in nearly all programming languages.

advanced awk
============
You can do a lot more with awk, here are several resources:

    - http://www.hcs.harvard.edu/~dholland/computers/awk.html

    - http://doc.infosnel.nl/quickawk.html

    - http://www.catonmat.net/download/awk.cheat.sheet.pdf

.. _label-class-4-exercises:

In Class Exercises - Class 4
============================
we will do the first 2. of these together

1. use NR to print each line of `lamina.bed` *preceded* by it's line number

  a. do the above, but only for regions on chromosome 12

2. use NF to see how many columns are in each row of `states.tab`

  a. use sort and uniq -c to see uniq column counts.
  b. why are there 2 numbers?
  c. can you adjust the file separator so that awk thinks all rows have
     the same number of columns?

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

In Class Exercises - Class 4 (2)
================================

 #. are there any regions in `lamina.bed` with start > end?

 #. what is the total coverage [sum of (end - start)] of regions on chr13 in `lamina.bed`?

 #. what is the mean value (4th column) on chromome 3 of `lamina.bed`

 #. print out only the header and the entry for colorado in `states.tab`

 #. what is the (single-number) sum of all the incomes for `states.tab` with illiteracy rate:
    a. less than 0.1?
    b. greater than 2?

 #. use NR to filter out the header from `lamina.bed` (hint: what is NR for the header?)

.. raw:: pdf

    PageBreak
