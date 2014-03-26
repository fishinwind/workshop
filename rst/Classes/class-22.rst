
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
    differences between:

    - conditions
    - conditions grouped by gender

 #. 

ggplot maninpulations
=====================

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

.. nextslide::
    :increment:

You can also fit a linear model separately and plot the data:

.. code-block:: r

    # intercept = coef(model)[1]
    # slope = coef(model)[2]

    > model <- lm(age ~ expr, data = covs)
    > gp + geom_abline(intercept = coef(model)[1], 
                       slope = coef(model)[2])

Exercises
---------

foo


