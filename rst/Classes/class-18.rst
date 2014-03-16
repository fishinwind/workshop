******************************
Class 18 : R : Getting Started
******************************

Goals
=====

 #. Startup RStudio
 #. Learn to navigate in RStudio
 #. load the ggplot2 library
 #. Learn to load external data

RStudio
=======

RStudio is a useful interface for R, encapsulates the R prompt, data views
and plotting in a single window.

.. code-block:: bash

    $ rstudio

ggplot2
=======
ggplot2 is a powerful plotting library.

.. code-block:: r

    # examine an internal data set
    > mtcars
    > summary(mtcars)
    # make a simple plot
    > library(ggplot2)
    > qplot(mpg, cyl, data=mtcars)

    # change colors and sizes
    > qplot(mpg, carb, size=factor(cyl), color=factor(gear), data=mtcars)

.. nextslide::
    :increment:

.. code-block:: r

    # XXX doesn't work
    > plot <- ggplot(aes(x = cyl, y = gear), data=mtcars)
    > plot <- plot + geom_bar()

Loading data in R
=================

.. code-block:: r

    # ``df`` is short for `data frame`, the basic R data structure
    # ``<-`` is the R assignment operator, can also use ``=``
    > df <- read.delim('data.tab')
    > summary(df)

In Class Exercise
=================

    #. Make a plot of XXX by YYY.

.. raw:: pdf

    PageBreak
