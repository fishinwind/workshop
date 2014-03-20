**************************************
Class 20 : R : Manipulation & Plotting
**************************************

Goals
=====

 #. Learn to keep a running script in RStudio to save your work.

 #. Exercises 

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

Exercises
=========

#. Load the `expr-geno-covs.txt <../misc/data/expr-geno-covs.txt>` data.
   Use ``ddply`` to calculate the mean age of smokers grouped by gender
   and smoking status. Plot the result.

#. Make a plot of age by expression faceted by genotype. Fit a linear
   model through these curves (use geom_smooth) on the plot.

#. Load a BED file (e.g. ``lamina.bed``) and calculate the mean length of
   regions on each chromosome in the BED file with plyr.  Plot the result as
   a bar plot with ggplot2.

