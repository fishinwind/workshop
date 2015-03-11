
.. include:: /_static/substitutions.txt 

.. _problem-set-6:

***************
 Problem Set 6
***************

:Due date: |pset6-due| 

Overview
========

**For the following problems, generate a Rmarkdown document that explains
the workflow and generates the appropriate plots.**

In this problem set you will do some exploratory data analysis with the
:url:`nihexporter <https://github.com/jayhesselberth/nihexporter>` R data
package.

Problems
========

#. Perform `left_join()` operations on all pairs of tables that have
   matching column names. Use `select()` to pull one column (in addition
   to the matching one) from each of the two tables, and print the first 5
   rows of the combined table (**10 points**).

#. Compare the costs of R01 and P01 grants over time (i.e. across fiscal
   years) for the `GM` and `CA` institutes. Plot the result using
   `geom_boxplot()` and `facet_grid()` (**10 points**).

#. Calculate "productivity" (i.e. number of publications per dollar of
   `total.cost` for R01 grants across `study.section` (**10 points**).
   Make a summary table of the result.

#. Identify grant supplements (i.e. `suffix` entries that start with "S")
   that have the highest total cost over the shorted duration (**10
   points**).

#. Use dplyr `slice()` to determine the highest producing grants (pubs / dollar)
   of any type from *each* institute (**10 points**). 

Problem Set Submission
======================

Submit your problem set as a tar file to Canvas
(:ref:`problem-set-submission`).

.. raw:: pdf

    PageBreak

