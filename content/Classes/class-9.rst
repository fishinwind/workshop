
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

Combining with ggplot
=====================

We had mean expression by condition and genotype as:

.. code-block:: r

    dfx %>% group_by(condition, genotype) \
        %>% summarize(count=n(), mean.expr=mean(expression)) \

We can add to that expression (after typing 'library(ggplot2)')

.. code-block:: r

        %>% ggplot(aes(x=genotype, y=expression) \
            + geom_histogram(stat='identity')

how can we change the color of all the bars to 'red'? [Hint, it's not
**color** ='red']

ggplot histograms
=================

Since `expr-geno-covs.txt` is already in long format, we can use it directly in
ggplot:

.. code-block:: r

    ggplot(covs, aes(x=expression)) + 
           geom_histogram() +
           scale_x_log10()

Exercise
========
Adjust this:

.. code-block:: r

    ggplot(covs, aes(x=expression)) + 
           geom_histogram() +
           scale_x_log10()

#. to color by genotype

#. and to split plots (facet_wrap) by condition (case/control)

#. to color by age > 60 vs. <= 60 (use row selection stuff from start of class to
   make a new column named, e.g. `is_old`)

Exercises
---------

#. Figure out how to move overlapping points so categorical data is
   viewable (hint: look at geom_jitter() or the `position` argument to
   geom_point()) 

#. Load a BED file (e.g. ``lamina.bed``) and calculate the mean length of
   regions on each chromosome in the BED file with dplyr.  Plot the result as
   a bar plot with ggplot2.

#. Work through the ``dplyr`` vignette.
   http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html

.. raw:: pdf

    PageBreak

