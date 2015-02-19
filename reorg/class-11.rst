**************************************
Class 11 : R : Manipulation & Plotting
**************************************

Goals
=====

#. Learn to keep a running script in RStudio to save your work.
#. more reshape and dplyr
#. Exercises 

RStudio practice
================

#. open a new R script with File -> New -> R script

#. As you run interactive analyses, paste the commands that work into the
   script.

reshape2 review
===============

There are two functions in the reshape2 package: ``melt()`` and
``dcast()``. The ``d`` means it returns a data.frame.

``dcast()`` takes a data.frame, a forumula and an optional function

.. code-block:: r

   > library(reshape2)
   > library(dplyr)
   > mycars <- mutate(mtcars, car.name = rownames(mtcars))
   > mycarsm <- melt(mycars, id=c('car.name','gear'))
   # examine mycarsm. is it wide or long?
   > dcast(mycarsm, gear ~ variable, mean)

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

#. Melt the `expr-geno-covs.txt` data table. Recast it with ``dcast()``
   and calculate the mean for each variable conditioned on gender. Plot
   the result.

#. Use ``dplyr`` to calculate the mean age of smokers grouped by gender
   and smoking status. Plot the result.

#. Make a plot of age by expression faceted by genotype. Fit a linear
   model through these curves (use geom_smooth) on the plot.

#. Load the peaks BED file and find the 10 factors that have the largest range
   in peak width. Inspect a ``geom_boxplot()`` or ``geom_violin()`` to support
   your answer (also add individual points to the plot with ``geom_jitter()``).

