
.. include:: /_static/substitutions.txt

********************************************
Class 11 : R : Simple statistics and ggplot2
********************************************

:Class date: |c11-date|
:Last updated: |today|

Goals
=====

#. Simple statistics 
#. ggplot2 manipulations

Statistics in R
===============

R provides a number of builtin statistics. Let's go over some using expr-geno-covs.txt. First, load the dataframe.

.. codeblock:: r
	> df <- read.table(file = "expr-geno-covs.txt", sep = "\t", header = TRUE)

.. nextslide::
    :increment:

Student t-test
--------------

Student t-test is used to determine if two sets of values are significantly different.
- Assumptions are data has a normal distribution and is continuous.

.. codeblock:: r
	> exp_male <- df$expression[df$gender == "Male"]
	> exp_female <- df$expression[df$gender == "Female"]
	> ttest <- t.test(exp_male, exp_female)
	> pVal <- ttest$p.value

.. nextslide::
    :increment:

Fisher's Exact Test
-------------------

Contingency tables are matrices that describe the frequency of variables. 
- Determine if variable "a" occurs at a higher frequency compared to variable "b".
- Fisher's Exact Test is used to determine the statistical significance.

From our data, let's see if there is difference in frequency of former smokers between males and females.

.. codeblock:: r

	# generate contingency table
	> mf <- length(which((df$gender == "Male" & df$smoking == "former") == TRUE))
	> mNf <- length(which((df$gender == "Male" & df$smoking != "former") == TRUE))
	> cf <- length(which((df$gender == "Female" & df$smoking == "former") == TRUE))
	> cNf <- length(which((df$gender == "Female" & df$smoking != "former") == TRUE))
	> con <- matrix(c(mf,mNf,cf,cNf), nrow = 2, ncol = 2, byrow = FALSE)
	
	# perform Fisher's exact test on contingency table
	> fisherTest <- fisher.test(con)
	> pVal <- fisherTest$p.value

.. nextslide::
    :increment:

Wilcoxon's Test
---------------

This test is used on two groups of data where the samples are somehow related. Unfortunately, there is no example in df so we will have to use other data. 
- Use ``immer``; describes the yield of barley of 30 different locations in two different years.

.. codeblock:: r

	> wilcTest <- wilcoxon.test(immer$Y1,immer$Y2)
	> wilcTest$p.value

.. nextslide::
    :increment:

Kolmogorov-Smirnov
------------------

This test determines if two distributions of data are similar.
- Data does not have to have a normal distribution.
- Data can be used continuous or discrete (count data).

.. codeblock:: r

	> ksTest <- ks.test(exp_male, exp_female)
	> ksTest$p.value

.. nextslide::
    :increment:

Plot the distribution!
----------------------

Does the distribution make sense with our results??

.. codeblock:: r

	> gp <- ggplot(df, aes(x = expression, fill = gender))
	> gp + geom_bar()

.. nextslide::
    :increment:

Exercises
---------

#. Use ``t.test()`` determine whether there are significant expression
   differences between in `expr-geno-covs.txt` using conditions instead of gender.

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
> gp <- ggplot(df, aes(x = age, y = expression)
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

R Model Syntax
==============

You can also fit and examine a linear model:

.. code-block:: r

    > model <- lm(expression ~ genotype + condition + gender, data = df)
    > summary(model)

    # look at diagnostic plots
    > plot(model)

What assumptions can we make from the diagnostics?
What do you suggest we do to fix it? 
- don't look ahead it's cheating >:(

	.. nextslide::
	    :increment:

Let's Normalize
---------------

If we normalize the data, the diagnostics might fit.
- Log that data.

But wait, there are 0 values - get values that do not exist. What on earth do we do??!
- Filter out values.
- After normalizing, replace again with 0.

.. code-block:: r
	# example if we filter out DNE values
	> df2 <- df
	> df2$expression <- log2(df2$expression)
	> filterOut <- which(is.infinite(df2$expression))
	> df2 <- df2[-filterOut,]
	> model <- lm(expression ~ genotype + condition + gender , data = df2)
	> plot(model)

How would you replace values with 0?

- MAJOR NOTE: you don't need to know this stuff for the homework, I'm just having fun.
- Also! knowing some of these tricks may make your life easier with your own work.

.. nextslide::
    :increment:

Exercises
---------

#. Does adding age to the existing model (expression ~ genotype + condition +
   gender) change the significance of the other variables?

#. How does removing condition from the model affect the significance of
   genotype and vice-versa?

