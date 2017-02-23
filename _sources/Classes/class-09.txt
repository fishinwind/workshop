
.. include:: /_static/substitutions.txt

**************************************
Class 9 : R : Manipulation & Plotting
**************************************

:Class date: |c9-date|
:Last updated: |today|

Goals
=====

#. Begin to make plots with ggplot2 and combine with dplyr 

ggplot2 syntax
==============

.. code-block:: r

    library(tidyverse)
    df <- read_tsv('expr-geno-covs.txt')

    ggplot(df, aes(x=genotype, y=expression))

    ggplot(df, aes(x=genotype, y=expression)) +
        geom_point()

    ggplot(df, aes(x=genotype, y=expression)) +
        geom_point(color='red', size = 3)

.. nextslide::
    :increment:

``aes()`` stands for **aesthetics**, which specifies the the coordinates,
colors, size, etc from these columns in the data.frame.

.. code-block:: r

    aes(x=genotype, y=expression, color=gender)

``geom_point()`` means plot these as points, could be ``geom_line()`` or 
a number of other `geoms`.

ggplot2 docs
=====================

The ggplot2 docs are very good: http://docs.ggplot2.org/current/

Use them find how to change the y-scale on this plot to log10:

.. code-block:: r

    ggplot(dfx, aes(x=genotype, y=expression)) +
            geom_point()

Look at the `geom_point()` documentation and change the color
of the plot above so that males and females are colored differently.

Histograms
==========
One of the simplest things to do in R, without ggplot is to look at 
a histogram of your data:

.. code-block:: r

    hist(df$expression)
    # or
    hist(log(df$expression))

You can make these look a lot nicer with ggplot2.

``hist()`` does not work with ggplot, you'll have to use the
ggplot2 machinery for that:

.. nextslide::
    :increment:

.. code-block:: r

    ggplot(df, aes(x=log10(expression))
      + geom_histogram()

Try to figure out how to overlay a density plot with ``geom_density()``.

.. nextslide::
    :increment:

You can also make effective separate visualizations with ``facet_grid()``
and ``facet_wrap()``. You need to specify a formula to determine how to
separate:

.. code-block:: r

    gp <- ggplot(df, aes(x=genotype, y=expression)) + geom_point()

    # separate by gender
    gp + facet_grid(~ gender)

    # separate by gender and condition
    gp + facet_grid(condition ~ gender)

Exercises
=========

#. Make a histogram using ggplot and separate cases from controls by
   changing fill color and symbol type.

#. Make a histogram using ggplot and separate cases from controls
   using ``facet_grid()``, ``facet_wrap()``.

Exercises with dplyr
--------------------

Use the ``expr-geno-covs.txt`` file for the following exercises.

#. Use ``dplyr`` to calculate the mean age of smokers grouped by gender
   and smoking status. Plot the result.

#. Make a plot of age by expression faceted by genotype. Fit a linear
   model through these curves (use geom_smooth) on the plot.

Combining with ggplot
=====================

We had mean expression by condition and genotype as:

.. code-block:: r

    df %>% group_by(condition, genotype) \
        %>% summarize(count = n(),
                      expr.mean = mean(expression))

We can continue the thought by adding a plot with the pipe:

.. code-block:: r

        %>% ggplot(aes(x=genotype, y=expression)) \
            + geom_histogram(stat='identity')

how can we change the color of all the bars to 'red'? [Hint, it's not
**color** ='red']

ggplot histograms
=================

Since `expr-geno-covs.txt` is already in long format, we can use it directly in
ggplot:

.. code-block:: r

    ggplot(df, aes(x=expression)) + 
           geom_histogram() +
           scale_x_log10()

Exercises
---------
#. Adjust this:

.. code-block:: r

    ggplot(df, aes(x=expression)) + 
           geom_histogram() +
           scale_x_log10()

- to color by genotype

- and to split plots (facet_wrap) by condition (case/control)

- to color by age > 60 vs. <= 60 (use row selection stuff from start of class to
  make a new column named, e.g. `is_old`)

#. Figure out how to move overlapping points so categorical data is
   viewable (hint: look at geom_jitter() or the `position` argument to
   geom_point()) 

#. Load a BED file (e.g. ``lamina.bed``) and calculate the mean length of
   regions on each chromosome in the BED file with dplyr.  Plot the result as
   a bar plot with ggplot2.

Practice
--------
Start the Datacamp dplyr and ggplot2 courses if you haven't already.

.. raw:: pdf

    PageBreak

