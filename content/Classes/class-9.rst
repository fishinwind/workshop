
.. include:: /_static/substitutions.txt

**************************************
Class 9 : R : Manipulation & Plotting
**************************************

:Class date: |c9-date|


Goals
=====

#. Begin to manipulate data with and ``dplyr`` and ``tidyr``

dplyr
=====

``dplyr`` provides several methods:

- ``summarise()``
- ``filter()``
- ``select()``
- ``mutate()``
- ``arrange()``
- ``group_by()``

``dplyr`` also provides an operator called ``%>%`` that allows you to
string manipulations together.

By combining these functions together, you can quickly and easily
manipulate data.frames and generate useful data summaries.

.. code-block:: r

    > summary <- df %>% select() %>% group_by() %>% summarize()

dplyr example 
-------------

.. code-block:: r

    dfx <- read.table('expr-geno-covs.txt', header=TRUE)

    # calculate some simple stats with dplyr
    grouped <- group_by(dfx, condition, genotype)
    summarize(grouped, count = n(), mean.age = mean(age))

    # even better
    dfx %>% 
        group_by(condition, genotype) %>%
        summarize(count = n(), mean.age = mean(age))

.. nextslide::
    :increment:

Summarize transcription factor binding site peaks in::
  
    # columns 1-4 from merged peak calls
    /vol1/opt/data/encode/peaks.bed.gz

.. code-block:: r

    > library(dplyr)
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
+ assign the summary to a new data.frame

Exercises with dplyr
--------------------

#. Use ``dplyr`` to calculate the mean age of smokers grouped by gender
   and smoking status. Plot the result.

#. Make a plot of age by expression faceted by genotype. Fit a linear
   model through these curves (use geom_smooth) on the plot.

#. Load the peaks BED file and find the 10 factors that have the largest range
   in peak width. Inspect a ``geom_boxplot()`` or ``geom_violin()`` to support
   your answer (also add individual points to the plot with ``geom_jitter()``).

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

Exercises
---------

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

