************************************
Class 20 : R : DataFrames & Plotting
************************************

Goals
=====

 #. Understand data.frames in more detail
 #. Practice manipulating data.frames
 #. melt, dplyr, ggplot

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

dplyr review
============

``dplyr`` provides these simple methods:

    #. ``summarise()``
    #. ``filter()``
    #. ``select()``
    #. ``mutate()``
    #. ``arrange()``
    #. ``group_by()``


dplyr
=====
``dplyr`` also provides an operator called ``%.%`` that allows you to
chain manipulations together.

To get mean expression level by condition (case/control)

.. code-block:: r

    covs %.% group_by(condition) \
         %.% summarize(count=n(), mean.expr=mean(expression))

Mean expression by condition and genotype

.. code-block:: r

    covs %.% group_by(condition, genotype) \
         %.% summarize(count=n(), mean.expr=mean(expression))


Exercise
========

#. What is the difference in mean age between cases and controls?
#. What is the difference in mean age by genotype?


ggplot
======

Above, we had mean expression by condition and genotype as:

.. code-block:: r

    covs %.% group_by(condition, genotype) \
         %.% summarize(count=n(), mean.expr=mean(expression)) \

We can add to that expression (after typing 'library(ggplot2)')

.. code-block:: r

        %.% ggplot(aes(x=genotype, y=expression) \
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
#. to color by age (> 50 vs < 50)

