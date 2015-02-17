**************************************
Class 10 : R : Manipulation & Plotting
**************************************

Goals
=====

#. Begin to manipulate data with ``tidyr`` and ``dplyr``

Introduction to tidyr 
=======================

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

How is ``tidyr`` useful?
==========================

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

Exercises
=========

#. Figure out how to move overlapping points so categorical data is
   viewable (hint: look at geom_jitter() or the `position` argument to
   geom_point()) 

#. Load a BED file (e.g. ``lamina.bed``) and calculate the mean length of
   regions on each chromosome in the BED file with dplyr.  Plot the result as
   a bar plot with ggplot2.

#. Worth through the ``dplyr`` vignette.
   http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html

.. raw:: pdf

    PageBreak

