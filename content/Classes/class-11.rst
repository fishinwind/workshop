
.. include:: /_static/substitutions.txt

********************************************
Class 11 : R : Simple statistics and ggplot2 by Charlotte 
********************************************

:Class date: |c11-date|
:Last updated: |today|

Goals
=====

#. Simple statistics 
#. ggplot2 manipulations

Statistics in R
===============

R provides a number of builtin statistics. Let's go over some using expr-geno-covs.txt. First, let's load the dataframe.

.. codeblock:: r
> df <- read.table(file = "expr-geno-covs.txt", sep = "\t", header = TRUE)

Student t-test
--------------

Student t-test is used to determine if two sets of values are significantly different. We assume data has a normal distribution and is continuous.
- R function is ``t.test()``

.. codeblock:: r
> exp_male <- df$expression[df$gender == "Male"]
> exp_female <- df$expression[df$gender == "Female"]
> ttest <- t.test(exp_male, exp_female)
> pVal <- ttest$p.value

Fisher's Exact Test
-------------------

Contigency tables are matrices that describe the frequency of variables. They can be used to determine if a variable occurs at a higher frequency compared to another variable. Fisher's Exact Test is used to determine the statistical significance.

From our data, let's see if there is difference in frequency of former smokers between males and females.

.. codeblock:: r

# generate contingency table
> mf <- length(which((df$gender == "Male" & df$smoking == "former") == TRUE))
> mNf <- length(which((df$gender == "Male" & df$smoking != "former") == TRUE))
> cf <- length(which((df$gender == "Female" & df$smoking == "former") == TRUE))
> cNf <- length(which((df$gender == "Female" & df$smoking != "former") == TRUE))
> con <- matrix(c(mf,mNf,cf,cNf), nrow = 2, ncol = 2, byrow = FALSE)
> print(con)
> fisherTest <- fisher.test(con)
> pVal <- fisherTest$p.value

Exercises
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
    > gp <- ggplot(covs, 
                   aes(x = age, y = expression, 
                   group = plyr::round_any(age, 2)))

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

