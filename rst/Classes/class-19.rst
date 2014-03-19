**************************************
Class 19 : R : Manipulation & Plotting
**************************************

Goals
=====
 #. Learn to load data with ``read.delim``
 #. Begin to manipulate data with ``reshape`` and ``plyr``

Loading data into R
===================
The main function for loading data is ``read.delim()``:

.. code-block:: r

    # dfx is a data.frame. look at it with ``summary`` and ``head``
    > dfx <- read.delim('lamina.bed')
    > is.data.frame(dfx)

If the column names are specified in a header line (begins with ``#``),
then you can load them as the column names with:

.. code-block:: r

    > dfx <- read.delim('lamina.bed', header=TRUE)

You can also specify headers explicitly with:

.. code-block:: r

    > bedfilename <- '/opt/bio-workshop/data/lamina.bed'
    > colnames <- c('chrom','start','end','name','score','strand')
    > dfx <- read.delim(bedfilename, col.names=colnames)

.. nextslide::
    :increment:

Recall that you access the column data with the ``$`` character:

.. code-block:: r

    > dfx$chrom 

You can also ``attach`` and ``detach`` data.frames if you will be
repeatedly accessing the data:

.. code-block:: r

    > attach(dfx)
    # now, columns can be referred to directly
    > chrom
    > detach(dfx)
    # back to original style
    > dfx$chrom

.. nextslide::
    :increment:

There are several data sets that are built-in to R, including:

.. code-block:: r

    > mtcars   # Motor Trend Cars Road Tests
    > baseball # in ``libarary(plyr)``

    # see all built-in data sets
    > library(help = "datasets")

Introduction to plyr
====================

The ``plyr`` package provides several methods for flexible manipulation of
data. Think "pliers".

``plyr`` uses the `split-apply-combine` strategy of data analysis.

The function we will use initially is called ``ddply``. The first ``d``
says that the input is a ``data.frame`` and the second ``d`` says that the
output is a ``data.frame``. 

Break down a ddply() function:

.. code-block:: r

   # ddply(data, variables, function)
   > ddply(mtcars, rownames(mtcars), summarize, mpg_g = mpg / gear)

``plyr`` provides several helper functions:

.. code-block:: r

   # see the bottom of this page
   > ?plyr

.. nextslide::
   :increment:

Exercise: calculate some statistics on the regions from each chromosome in
a BED file.

.. code-block:: r

    > library(plyr)

    # load data
    > colnames <- c('chrom','start','end')
    > dfx <- read.delim('/opt/bio-workshop/data/lamina.bed',
                        col.names=colnames)

    # measure lengths for each entry
    # mutate() is a plyr helper function that adds or modifies existing
    # columns
    > dfx <- mutate(dfx, length = end - start)

    > summary <- ddply(dfx, "chrom", summarize,
                       mean.len = mean(length),
                       median.len = median(length))

.. nextslide::
   :increment:


Introduction to reshape
=======================
The ``reshape`` package provides two important functions called
``melt()`` and ``cast()``.

These functions allow you to manipulate the shape of data frames. One
common operation is to convert data tables from `wide` format to `long`
format and back.

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

How is ``reshape`` useful?
==========================

The ``ggplot2`` expects data in ``long`` format, where points are
categorized.

**Question:** Look at the ``summary`` data.frame. Is it in ``wide`` or
``long`` format?

.. nextslide::
   :increment:

The data.frame from plyr is in ``wide`` format. 

.. code-block:: r

    > library(reshape2)
    > library(ggplot2)

    # covert to long format
    > long.summary <- melt(summary, id=c('chrom'))

    > gp <- ggplot(long.summary, aes(x=chrom, y=value, fill=variable))
    > gp + geom_bar(stat='identity', position='dodge')

Exercises
=========

#. Figure out how to move overlapping points so categorical data is
   viewable (hint: look at geom_jitter() or the `position` argument to
   geom_point()) 

#. Load a BED file (e.g. ``lamina.bed``) and calculate the mean length of
   regions on each chromosome in the BED file with plyr.  Plot the result as
   a bar plot with ggplot2.

#. Install ``dplyr`` and work through the vignette.

   http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html

.. raw:: pdf

    PageBreak

