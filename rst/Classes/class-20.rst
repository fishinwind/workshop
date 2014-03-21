**************************************
Class 20 : R : Manipulation & Plotting
**************************************

Goals
=====

 #. Learn to keep a running script in RStudio to save your work.

 #. Exercises 

<<<<<<< HEAD
RStudio practice
================

#. open a new R script with File -> New -> R script

#. As you run interactive analyses, paste the commands that work into the
   script.

reshape2 review
===============

There are two functions in the reshape2 package: ``melt()`` and
``dcast()``. The ``d`` means it returns a data.frame.

dplyr
=====

The ``dplyr`` package is uses syntax that is much more appealing that
``plyr``. It provides these simple methods:

    #. ``summarise()``
    #. ``filter()``
    #. ``select()``
    #. ``mutate()``
    #. ``arrange()``
    #. ``group_by()``

``dplyr`` also provides an operator called ``%.%`` that allows you to
string manipulations together:

.. code-block:: r

    > summary <- dfx %.% select() %.% group_by() %.% summarize()

.. nextslide::
    :increment:

Summarize transcription factor binding site peaks:

.. code-block:: r

    > colnames = c('chrom','start','end','name')
    > bedfilename = 'peaks.bed.gz'
    # use ``gzfile`` to load gzipped data
    > peaks <- read.table(gzfile(bedfilename), col.names=colnames)

    > peaks %.% 
        group_by(name) %.%
        mutate(peak.width = end - start) %.%
        filter(peak.width > 500 ) %.%
        summarize(count = n(), mean.width = mean(peak.width)) %.%
        arrange(desc(count))

+ ``group_by()`` takes the place of the variables in ``ddply``
+ ``n()`` is a special function for counting observations
+ assign the whole thing to a new data.frame
=======
plyr / dplyr comparison
=======================

.. code-block:: r

    dfx <- read.table('misc/data/expr-geno-covs.txt', header=TRUE)

    # calculate some simple stats with plyr and dplyr

    # with plyr
    ddply(dfx, .(condition, genotype), summarise, mean.age = mean(age), count = n())

    # with dplyr
    grouped <- group_by(dfx, condition, genotype)
    summarize(grouped, count = n(), mean.age = mean(age))
>>>>>>> FETCH_HEAD

Exercises
=========

<<<<<<< HEAD
#. Melt the `expr-geno-covs.txt` data table. Recast it with ``dcast()``
   and calculate the mean for each variable conditioned on gender. Plot
   the result.

=======
>>>>>>> FETCH_HEAD
#. Load the `expr-geno-covs.txt <../misc/data/expr-geno-covs.txt>` data.
   Use ``dplyr`` to calculate the mean age of smokers grouped by gender
   and smoking status. Plot the result.

#. Make a plot of age by expression faceted by genotype. Fit a linear
   model through these curves (use geom_smooth) on the plot.

#. Load a BED file (e.g. ``lamina.bed``) and calculate the mean length of
<<<<<<< HEAD
   regions on each chromosome in the BED file with dplyr. Plot the result as
   a bar plot.
=======
   regions on each chromosome in the BED file with plyr.  Plot the result as
   a bar plot with ggplot2.
>>>>>>> FETCH_HEAD

