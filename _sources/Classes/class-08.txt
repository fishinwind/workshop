
.. include:: /_static/substitutions.txt

******************************
Class 8 : R : Getting Started
******************************

:Class date: |c8-date|
:Last updated: |today|

Goals
=====

#. Startup RStudio
#. Learn to navigate in RStudio
#. Load the tidyverse and explore

RStudio
=======

RStudio is a useful interface for R, encapsulates the R prompt, data views
and plotting in a single User Interface.

Useful panes in RStudio include:

- editor
- environment
- console
- plots
- help

DataFrames and tibbles
======================

One of the base R data structures is the ``data.frame``. 

The ``data.frame`` is a rectangular data structure with columns and rows.

There are a few problems with base ``data.frame``, but the ``tibble``
package deals with these nicely.

.. code-block:: r

    > iris
    > class(iris)

    > as_data_frame(iris)
    > class(as_data_frame(iris)

Access the column data with the ``$`` character:

.. code-block:: r

    > iris$Species 

.. nextslide::
    :increment:

There are several data sets built-in to R, including:

.. code-block:: r

    # Motor Trend Cars Road Tests
    > mtcars 

    # see all built-in data sets
    > library(help = "datasets")
   
libraries
=========

.. code-block:: r

    > library(ggplot2)
    # or
    > library('ggplot2')

Paths in RStudio
================

In RStudio, you can use the Files tab to get and set the working
directory.

Alternatively you can use the "Projects" facility in RStudio, which
manages the working directory automatically. Your problem sets would each
be a project.
Goals

dplyr
=====

``dplyr`` is one of the most powerful data manipulation libraries out
there. It is also very easy to learn, once you master the core verbs

``dplyr`` verbs:

- ``filter()``
- ``select()``
- ``mutate()``
- ``arrange()``
- ``group_by()``
- ``summarise()``

``dplyr`` also provides an operator called ``%>%`` that allows you to
string manipulations together.

.. code-block:: r

    > summary <- df %>% select() %>% group_by() %>% summarize()

dplyr example 
-------------

.. code-block:: r

    library(tidyverse)

    dfx <- read_tsv('expr-geno-covs.txt')

    # calculate some simple stats with dplyr
    grouped <- group_by(dfx, condition, genotype)
    summarize(grouped, count = n(), mean.age = mean(age))

    # even better
    dfx %>% 
        group_by(condition, genotype) %>%
        summarize(count = n(), mean.age = mean(age))

.. nextslide::
    :increment:

Summarize transcription factor binding site peaks in
``encode.tfbs.chr22.bed.gz``.
  
.. code-block:: r

    # don't forget to reset the working directory
    > colnames <- c('chrom', 'start', 'end', 'name')
    > bed <- 'encode.tfbs.chr22.bed.gz'

    > peaks <- read_tsv(bed, col_names = colnames)

    > peaks %>% 
        group_by(name) %>%
        mutate(peak.width = end - start) %>%
        filter(peak.width > 500) %>%
        summarize(count = n(), mean.width = mean(peak.width)) %>%
        arrange(desc(count))

``n()`` is a special function for counting observations.

Tidying data
============

Raw data must often be restructured to enable downstream analysis. A
useful format to start with is "tidy" data.

Tidy data has a few key principles:

-  Each variable forms a column
-  Each observation forms a row
-  Each type of observational unit forms a table

Common problems with messy data
-------------------------------

-  Column headers are values, not variable names.
-  Multiple variables are stored in one column.
-  Variables are stored in both rows and columns.

The ``anscombe`` data sets illustrate a key point of data analysis:
**visualize your data**.

DataCamp and Swirl
==================

You should begin learning ``dplyr`` and ``tidyr`` using DataCamp and
Swirl.

.. raw:: pdf

    PageBreak

