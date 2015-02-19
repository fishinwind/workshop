
.. include:: /_static/substitutions.txt

******************************
Class 8 : R : Getting Started
******************************

:Class date: |c8-date|


Setting up web access to tesla
==============================

You can access html and images via the web if you:

#. Create a `public_html` directory in your $HOME::
    
    $ mkdir $HOME/public_html

#. Copy content into the `public_html` directory.

#. Make sure that permissions of your $HOME are accessible by all::

    $ chmod 755 $HOME

#. Access tesla on campus, or via VPN at the URL::

    http://amc-tesla/~username

    http://amc-tesla.ucdenver.pvt/~username

Goals
=====

#. Startup RStudio
#. Learn to navigate in RStudio
#. load the ggplot2 library
#. Learn to load external data

Data
====

The data file we will use is::

    /vol1/opt/data/expr-geno-covs.txt

Copy this to a working directory e.g. ``$HOME/class/class-8``

RStudio
=======

RStudio is a useful interface for R, encapsulates the R prompt, data views
and plotting in a single User Interface.

Navigate to http://amc-sandbox.ucdenver.pvt/rstudio and enter your tesla
credentials.

Note you need to use VPN if you are off campus.

Loading data in R
=================

The main function for loading data is ``read.delim()``:

.. code-block:: r

    # dfx is a data.frame. look at it with ``summary`` and ``head``
    > dfx <- read.delim('lamina.bed')
    > is.data.frame(dfx)

If the column names are specified in a header line (begins with ``#``),
then you can load them as the column names with:

.. code-block:: r

    # this uses the names in the comment line as header names
    > dfx <- read.delim('lamina.bed', header=TRUE)

You can also specify headers explicitly with:

.. code-block:: r

    > bedfilename <- '/vol1/opt/data/lamina.bed'
    > colnames <- c('chrom','start','end','name','score','strand')
    > dfx <- read.delim(bedfilename, col.names=colnames)

.. nextslide::
    :increment:

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

    read.csv('file.csv')
    read.delim('file.txt')

These take common arguments. You can get help on a function in R
with:

.. code-block:: r

   ?read.delim
   ?head

DataFrame
=========
The ``read()`` functions load a ``data.frame``. In a data.frame,
everything is read into memory.

+ R figures out the types of each column e.g. numeric/character/factor
+ each column of the data.frame is accessed by `$`  e.g df$genotype

Access the column data with the ``$`` character:

.. code-block:: r

    > dfx$chrom

And some functions allow you to refer to columns by name:

.. code-block:: r

    > subset(dfx, chrom == 'chr1')

.. nextslide::
    :increment:

There are several data sets that are built-in to R, including:

.. code-block:: r

    > mtcars   # Motor Trend Cars Road Tests
    > baseball # in ``libarary(plyr)``

    # see all built-in data sets
    > library(help = "datasets")

print
=====

.. code-block:: r

    > print("hello world")
    
libraries
=========

.. code-block:: r

    > library(ggplot2)
    # or
    > library('ggplot2')

R paths
=======

get/set working directory

.. code-block:: r

    > getwd()                    # print
    > setwd('C:\whatever\path\') # on windows
    > setwd('/vol1/opt/data/')   # on linux

Writing and running R scripts in RStudio
========================================
Open up a new R Script with `File --> New File --> R Script`. You should
see a new editor window open in the top left.

You can write a simple R script like:

.. code-block:: r

    plot(rnorm(1000))

And then save the file as `plot.R`.

You can then highlight all of the lines, or just select the `plot` line,
and run that portion of the program by presseing the `Run` button in the
window, or with `<Cmd>-<Return>` on Macs and `<Ctrl>-<Return>` on Windows.

ggplot2
=======
We will learn to use it to create plots like this

.. image:: /_static/images/ggplot-ex.png

ggplot2 syntax
==============

.. code-block:: r

    library(ggplot2)
    dfx <- read.delim('expr-geno-covs.txt')

    ggplot(dfx, aes(x=genotype, y=expression)) +
        geom_point()

.. nextslide::
    :increment:

``aes()`` stands for **aesthetics**, which specifies the the coordinates,
colors, size, etc from these columns in the data.frame.

.. code-block:: r

    aes(x=genotype, y=expression, color=gender)

``geom_point()`` means plot these as points, could be ``geom_line()`` or 
a number of other `geoms`.

googling with ggplot2
=====================

Use google to find how to change the y-scale on this plot to log10

.. code-block:: r

    library(ggplot2)
    dfx = read.delim('expr-geno-covs.txt')

    ggplot(dfx, aes(x=genotype, y=expression)) +
            geom_point()

answer
======

.. code-block:: r
    :emphasize-lines: 6

    library(ggplot2)
    dfx <- read.delim('expr-geno-covs.txt')

    ggplot(dfx, aes(x=genotype, y=expression)) +
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

Histograms
==========
One of the simplest things to do in R, without ggplot is to look at 
a histogram of your data:

.. code-block:: r

    hist(dfx$expression)
    # or
    hist(log(dfx$expression))

You can make these look a lot nicer with ggplot2.

``hist()`` does not work with ggplot, you'll have to use the
ggplot2 machinery for that:

.. nextslide::
    :increment:

.. code-block:: r

    gp <- ggplot(dfx, aex(x=log10(expression))
    gp + geom_histgram()

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

.. raw:: pdf

    PageBreak
