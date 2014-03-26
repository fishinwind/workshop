
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

 #. 

ggplot manipulations
====================

Often you will have observations of two variables that are both 
continuous data. How can you examine the relationships between these?

Use ``geom_boxplot()`` with the ``group`` aesthetic:

.. code-block:: r

    # ``round_any()`` is provided by ``plyr``
    > gp <- ggplot(covs, aes(x = expr, y = age,
                             group = round_any(expr, 100)))
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

You can also fit a linear model separately and plot the data:

.. code-block:: r

    > model <- lm(expression ~ genotype + condition + gender , data = covs)
    > summary(model)


    > gp + geom_abline(intercept = coef(model)[1], 
                       slope = coef(model)[2])

Exercises
---------

#. Does adding age to the existing model (expression ~ genotype + condition +
   gender) change the signficance of the other variables? 

#. How does removing condition from the model affect the significance of
   genotype and vice-versa?



