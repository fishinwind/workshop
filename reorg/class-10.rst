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

dplyr
=====

``dplyr`` provides these simple methods:

#. ``summarise()``
#. ``filter()``
#. ``select()``
#. ``mutate()``
#. ``arrange()``
#. ``group_by()``

``dplyr`` also provides an operator called ``%>%`` that allows you to
string manipulations together:

.. code-block:: r

    > summary <- df %>% select() %>% group_by() %>% summarize()

dplyr example 
=============

.. code-block:: r

    dfx <- read.table('misc/data/expr-geno-covs.txt', header=TRUE)

    # calculate some simple stats with dplyr
    grouped <- group_by(dfx, condition, genotype)
    summarize(grouped, count = n(), mean.age = mean(age))

    # even better
    dfx %>% 
        group_by(condition, genotype) %>%
        summarize(count = n(), mean.age = mean(age))

dplyr example
=============

Fetch the peaks.bed.gz file <http://amc-sandbox.ucdenver.edu/~jhessel/outbox/2014/peaks.bed.gz>

.. nextslide::
   :increment:

Summarize transcription factor binding site peaks:

.. code-block:: r

    > colnames = c('chrom','start','end','name')
    > bedfilename = 'peaks.bed.gz'
    # use ``gzfile`` to load gzipped data
    > peaks <- read.table(gzfile(bedfilename), col.names=colnames)

    > peaks %>% 
        group_by(name) %>%
        mutate(peak.width = end - start) %>%
        filter(peak.width > 500 ) %>%
        summarize(count = n(), mean.width = mean(peak.width)) %>%
        arrange(desc(count))

+ ``n()`` is a special function for counting observations
+ assign the whole thing to a new data.frame

Exercises
=========

#. Melt the `expr-geno-covs.txt` data table. Recast it with ``tidyr:gather()``
   and calculate the mean for each variable conditioned on gender. Plot
   the result.

#. Use ``dplyr`` to calculate the mean age of smokers grouped by gender
   and smoking status. Plot the result.

#. Make a plot of age by expression faceted by genotype. Fit a linear
   model through these curves (use ``geom_smooth()``) on the plot.

#. Load the peaks BED file and find the 10 factors that have the largest range
   in peak width. Inspect a ``geom_boxplot()`` or ``geom_violin()`` to support
   your answer (also add individual points to the plot with ``geom_jitter()``).

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

