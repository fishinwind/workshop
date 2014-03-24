************************************
Class 20 : R : DataFrames & Plotting
************************************

Goals
=====

 #. Understand data.frames in more detail
 #. Practice manipulating data.frames
 #. dplyr, reshape, ggplot

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
===========================

The order is rows, columns. By having nothing before the comma:

.. code-block:: r

    head(covs[,c("genotype", "expression")])

We extract all rows. We can manually extract the first 5 rows and our 2 columns
as:

.. code-block:: r

    covs[1:5, c("genotype", "expression")]



    

reshape2 review
===============

There are two functions in the reshape2 package: ``melt()`` and
``dcast()``. The ``d`` means it returns a data.frame.

``dcast()`` takes a data.frame, a forumula and an optional function

.. code-block:: r

   > library(reshape2)
   > library(dplyr)
   > mut.cars <- mutate(mtcars, car.name = rownames(mtcars))
   > melt.cars <- melt(mut.cars, id=c('car.name','gear'))
   # examine mycarsm. is it wide or long?
   > dcast(melt.cars, gear ~ variable, mean)

dplyr review
============

``dplyr`` provides these simple methods:

    #. ``summarise()``
    #. ``filter()``
    #. ``select()``
    #. ``mutate()``
    #. ``arrange()``
    #. ``group_by()``

``dplyr`` also provides an operator called ``%.%`` that allows you to
chain manipulations together:

.. code-block:: r

    dfx %.% 
        group_by(condition, genotype) %.%
        summarize(count = n(), mean.age = mean(age))

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

