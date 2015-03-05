
.. include:: /_static/substitutions.txt

************************************
Class 10 : R : DataFrames & Plotting
************************************

:Class date: |c10-date|
:Last updated: |today|

Goals
=====

#. Look at data.frames in more detail
#. Practice manipulating data.frames
#. more dplyr and ggplot

data.frames
===========

Read in a data.frame, view the first few lines and
extract a single column

.. code-block:: r

    # assumes workshop repo is in $HOME/devel
    setwd("~/devel/workshop/content/misc/data")

    # load expr data
    expr.df <- read.delim('expr-geno-covs.txt', header=TRUE)

    head(expr.df)
    head(expr.df$genotype)

    # or
    head(expr.df[,"genotype"])

    # using dplyr
    expr.df %>% select(genotype) %>% head()

To extract multiple selected columns:

.. code-block:: r

    head(expr.df[,c("genotype", "expression")])

    # using dplyr
    expr.df %>% select(genotype, expression) %>% head()

Extracting From data.frames
---------------------------

The order is rows, columns. By having nothing before the comma:

.. code-block:: r

    head(expr.df[,c("genotype", "expression")])

We extract all rows. We can manually extract the first 5 rows and our 2 columns
as:

.. code-block:: r

    expr.df[1:5, c("genotype", "expression")]

    # using dplyr
    expr.df %>% select(genotype, expression) %>% head(5)

Named rows in data.frames
-------------------------

In `expr-geno-covs.txt`, there is a sample column, we may wish to use
that as the row id:

.. code-block:: r

    expr.df = read.delim('expr-geno-expr.df.txt', row.names='sample_id')
    head(expr.df)

    expr.df[c('sample1', 'sample14'),]

    # or with row numbers:
    covs[c(1, 4),]

Boolean Operations on Columns
-----------------------------

.. code-block:: r
    
    is_old <- covs$age > 65

now `is_old` is a list of TRUE / FALSE values. We can extract only those > 65 as:

.. code-block:: r

    old <- covs[is_old,]
    # same as:
    old <- covs[covs$age > 65,]

    # using dplyr
    is_old <- covs %>% filter(age > 65)

We can combine selections with '&' for and and '|' for or

.. code-block:: r
    
    old_with_disease = (covs$age > 65) & (covs$condition == "case")
    owd = covs[old_with_disease,]

.. nextslide::
   :increment:

.. code-block:: r

   # with subset
   old <- subset(covs, age > 65)

   # select on membership
   genos <- c('AC', 'CA')
   hets  <- subset(covs, genotype %in% genos)

Excercises for data.frames
--------------------------

Remember for combining expressions, you can create a variable for each, `is_AA`,
`is_CC` and then combine after. **Do these using both dplyr and "old" R**.

#. How many people have genotype 'CC'

#. How many people have genotype 'CC' or 'AA'?

#. How many people have genotype of 'CC' or 'AA' and are under 65 years old.

#. How many *females* have genotype of 'CC' or 'AA' and are under 65 years old.

Rmarkdown
=========

We'll go over a Rmarkdown, a simple markup language that can be used
within RStudio to generate nice looking documents combining plots and
text. 

Github
------
Some may also find it useful to link their RStudio projects with Github
accounts.
