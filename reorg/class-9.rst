******************************
Class 9 : R : Getting Started
******************************

Goals
=====

#. Startup RStudio
#. Learn to navigate in RStudio
#. load the ggplot2 library
#. Learn to load external data

Data
====

Download the data.

**Download Link:** :download:`expr-geno-covs.txt <../misc/data/expr-geno-covs.txt>`


.. code-block:: bash

    cd /vol1/opt/data/
    wget https://ucd-bioworkshop.github.io/_downloads/expr-geno-covs.txt

RStudio
=======

RStudio is a useful interface for R, encapsulates the R prompt, data views
and plotting in a single User Interface.

.. code-block:: bash

    $ rstudio

Loading data in R
=================

 ``df`` is short for `data frame`, the basic R data structure
 ``<-`` is the R assignment operator, can also use ``=``

.. code-block:: r

    > df <- read.delim('expr-geno-covs.txt')
    > summary(df)
    > head(df)
    > head(df$expression)

.. nextslide::
    :increment:

The most common ways to read in data in R are:

   
.. code-block:: r

    read.csv('some.csv')
    read.delim('some.tab.txt')

These take common arguments. You can get help on a function in R
with:

.. code-block:: r

   ?read.delim
   ?head

print
=====

.. code-block:: r

    print("hello world")
    

libraries
=========

.. code-block:: r

    library(ggplot2)
    # or
    library('ggplot2')

R paths
=======

get/set working directory

.. code-block:: r

    getwd() # print
    setwd('C:\whatever\path\') # on windows
    setwd('/vol1/opt/data/') # on linux

ggplot2
=======
We will learn to use it to create plots like this

.. image:: ../_static/images/ggplot-ex.png

ggplot2 syntax
==============

.. code-block:: r

    library(ggplot2)
    df = read.delim('expr-geno-covs.txt')

    ggplot(df, aes(x=genotype, y=expression)) +
        geom_point()

.. nextslide::
    :increment:

`aes()` stands for **aesthetics**, means pull the coordiantes/colors/size/etc
from these columns in the data.frame.

.. code-block:: r

    aes(x=genotype, y=expression, color=gender)

`geom_point()` means plot these as points, could be geom_line() or 
   a number of other geom_ things.

googling with ggplot2
=====================

Use google to find how to change the y-scale on this plot to log10

.. code-block:: r

    library(ggplot2)
    df = read.delim('expr-geno-covs.txt')

    ggplot(df, aes(x=genotype, y=expression)) +
            geom_point()

answer
======

.. code-block:: r

    library(ggplot2)
    df = read.delim('expr-geno-covs.txt')

    ggplot(df, aes(x=genotype, y=expression)) +
            geom_point() +
            scale_y_log10()

You can find a lot of info for ggplot2 with some googling.

ggplot2 documentation
=====================

The ggplot2 docs are very good: http://docs.ggplot2.org/current/

Look at the `geom_point()` documentation and change the color
of the plot above so that males and females are color'ed differently.


DataFrame
=========
In a data.frame, we read everything into memory

+ R figures out if it is int/character/numeric
+ each column of the data.frame is accessed by `$`  e.g df$genotype

Hist
====
One of the simplest things to do in R, without ggplot is to look at 
a histogram of your data:

.. code-block:: r

    df = read.delim('expr-geno-covs.txt')
    hist(df$expression)
    # or
    hist(log(df$expression))

You can make these look a lot nicer with ggplot2.

**Hist** does not work with ggplot, you'll have to use the
ggplot2 machinery for that.

Exercise
========

Make a histogram using ggplot and separate cases from controls
either by **facet** or by **fill**. 

.. raw:: pdf

    PageBreak
