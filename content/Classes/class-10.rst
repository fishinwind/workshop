
.. include:: /_static/substitutions.txt

************************************
Class 10 : R : DataFrames & Plotting
************************************

:Class date: |c10-date|
:Last updated: |today|

Goals
=====

#. Understand data.frames in more detail
#. Practice manipulating data.frames
#. more dplyr and ggplot

data.frames
===========

Read in a data.frame, view the first few lines and
extract a single column

.. code-block:: r

    covs = read.delim('expr-geno-covs.txt')
    head(covs)
    head(covs$genotype)

    # or
    head(covs[,"genotype"])

To extract multiple selected columns:

.. code-block:: r

    head(covs[,c("genotype", "expression")])


Extracting From data.frames
---------------------------

The order is rows, columns. By having nothing before the comma:

.. code-block:: r

    head(covs[,c("genotype", "expression")])

We extract all rows. We can manually extract the first 5 rows and our 2 columns
as:

.. code-block:: r

    covs[1:5, c("genotype", "expression")]

Named rows in data.frames
-------------------------

In `expr-geno-covs.txt`, there is a sample column, we may wish to use
that as the row id:

.. code-block:: r

    covs = read.delim('expr-geno-covs.txt', row.names='sample_id')
    head(covs)

    covs[c('sample_1', 'sample_14'),]

We can still extract them with numbers:

    covs[c(1, 4),]

.. warning::

    Be careful with data.frames where the row.names are
    from an integer column


Boolean Operations on Columns
=============================

.. code-block:: r
    
    is_old = covs$age > 65

now `is_old` is a list of TRUE / FALSE values. We can extract only those > 65 as:

.. code-block:: r

    old = covs[is_old,]
    # same as:
    old = covs[covs$age > 65,]

We can combine selections with '&' for and and '|' for or

.. code-block:: r
    
    old_with_disease = (covs$age > 65) & (covs$condition == "case")
    owd = covs[old_with_disease,]

.. nextslide::
   :increment:

.. code-block:: r

   # with subset
   old <- subset(covs, age > 65)

   # select on membership
   genos <- c('AC', 'CA')
   hets  <- subset(covs, genotype %in% genos)

Excercises
==========

Remember for combining expressions, you can create a variable for each, `is_AA`,
`is_CC` and then combine after.

#. How many people have genotype 'CC'
#. How many people have genotype 'CC' or 'AA'?
#. How many people have genotype of 'CC' or 'AA' and are under 65 years old.
#. How many males have genotype of 'CC' or 'AA' and are under 65 years old.

Introduction to tidyr 
=====================

The ``tidyr`` package provides two important functions called
``gather()`` and ``separate()``.

These functions allow you to manipulate the shape of data frames. One
common operation is to convert data tables from `wide` format to `long`
format and back.

There are useful examples in the article describing the reshape package
[#]_. Check out the ``french fries`` case study.

.. [#] http://www.jstatsoft.org/v21/i12/paper

Wide format (or `unstacked`)
----------------------------
Values for each variable are in a separate column.

.. list-table::
    :header-rows: 1

    * - Person
      - Age
      - Weight
    * - Bob
      - 32
      - 128
    * - Alice
      - 24
      - 86
    * - Steve
      - 64
      - 95

Long format (or `stacked`)
--------------------------

One column contains the variables, one column contains the values.

.. list-table::
    :header-rows: 1

    * - Person
      - Variable
      - Value
    * - Bob
      - Age
      - 32
    * - Bob
      - Weight
      - 128
    * - Alice
      - Age
      - 24
    * - Alice
      - Weight
      - 86
    * - Steve
      - Age
      - 64
    * - Steve
      - Weight
      - 95

Why is ``tidyr`` useful?
------------------------

``ggplot2`` expects data in `long` format, where individual points are
categorized.

**Question:** Look at the ``summary`` data.frame. Is it in ``wide`` or
``long`` format?

.. nextslide::
   :increment:

The data.frame from dplyr is in ``wide`` format. 

.. code-block:: r

    > library(tidyr)
    > library(ggplot2)

    # covert to long format
    > long.summary <- gather(summary, chrom)

    > gp <- ggplot(long.summary, aes(x=chrom, y=value, fill=variable))
    > gp <- gp + geom_bar(stat='identity', position='dodge')
    > gp


