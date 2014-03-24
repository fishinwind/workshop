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

Named rows in data.frames
=========================

In `expr-geno-covs.txt`, there is a sample column, we may wish to use
that as the row id:

.. code-block:: r

    covs = read.delim('expr-geno-covs.txt', row.names='sample_id')
    head(covs)

    covs[c('sample_1', 'sample_14'),]

Boolean Operations on Colums
============================

.. code-block:: r
    
    is_old = covs$age > 65

now `is_old` is a list of true falses. We can extract only those > 65 as:

.. code-block:: r

    old = covs[is_old,]
    # same as:
    old = covs[covs$age > 65,]

We can combine selections with '&' for and and '|' for or

.. code-block:: r
    
    old_with_disease = (covs$age > 65) & (covs$condition == "case")
    owd = covs[old_with_disease,]

Excercises
==========

Remember for combining expressions, you can create a variable for each, `is_AA`,
`is_CC` and then combine after.

#. How many people have genotype 'CC'
#. How many people have genotype 'CC' or 'AA'?
#. How many people have genotype of 'CC' or 'AA' and are under 65 years old.
#. How many males have genotype of 'CC' or 'AA' and are under 65 years old.


reshape2 review
===============

There are two functions in the reshape2 package: ``melt()`` and
``dcast()``. The ``d`` means it returns a data.frame.


.. code-block:: r

   > library(reshape2)
   > library(dplyr)
   > mut.cars <- mutate(mtcars, car.name = rownames(mtcars))
   > melt.cars <- melt(mut.cars, id=c('car.name','gear'), measure.vars="mpg")


+ Q: what is melt.cars?
+ A: for every car:gear combination, it gives, the mpg

dcast
=====

``dcast()`` takes a data.frame, a forumula and an optional function

To get the mean mpg by gear (across all cars), we can do:

.. code-block:: r

   dcast(melt.cars, gear ~ variable, mean)


Or mean mpg by car (across all gears for that car):

.. code-block:: r

   dcast(melt.cars, car.name ~ variable, mean)


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

