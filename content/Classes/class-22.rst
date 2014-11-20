
*******************************************
Class 22 : R : Data Manipulation & Plotting
*******************************************

Goals
=====

#. Simple statistics 
#. ggplot2 manipulations

Statistics in R
===============

R provides a number of builtin statistics:

- ``t.test()``
- ``fisher.test()``
- ``wilcox.test()``
- ``ks.set()``

Each of these functions takes 2 vectors on input, and return a result
object:

.. code-block:: r

    > this <- rnorm(100)
    > that <- rnorm(100)
    > result <- t.test(this, that)

    > result$p.value

Excercises
----------

#. Use ``t.test()`` determine whether there are significant expression
   differences between in `expr-geno-covs.txt`:

- conditions
- gender

ggplot manipulations
====================

Once you make a plot that you like, you can save it with:

.. code-block:: r

   # ggsave uses the last plot by default and learns format from the file
   # suffix
   > ggsave('myplot.pdf')

.. nextslide::
    :increment:

Often you will have observations of two variables that are both 
continuous data. How can you examine the relationships between these?

Use ``geom_boxplot()`` with the ``group`` aesthetic:

.. code-block:: r

    # ``round_any()`` is provided by ``plyr``
    > library(plyr)
    > gp <- ggplot(covs, aes(x = age, y = expression, 
                   group = round_any(age, 2)))

    > gp + geom_boxplot()

.. nextslide::
    :increment:

Or you can plot the data and fit a curve or linear model to examine the
overall relationship:

.. code-block:: r

    # plot a smoothed, possibly curvy, line 
    > gp + geom_point() + geom_smooth()


    # fit and plot a straight line
    > gp + stat_smooth(method='lm')

R Model Syntax
==============

You can also fit and examine a linear model:

.. code-block:: r

    > model <- lm(expression ~ genotype + condition + gender , data = covs)
    > summary(model)

    # look at diagnostic plots
    > plot(model)

..  XXX is this model plottable?
..    # use ggplot to plot the data
..    > gp + geom_point()
..    > gp + geom_abline(intercept = coef(model)[1], 
..                       slope = coef(model)[2])

Exercises
---------

#. Does adding age to the existing model (expression ~ genotype + condition +
   gender) change the signficance of the other variables? 

#. How does removing condition from the model affect the significance of
   genotype and vice-versa?

