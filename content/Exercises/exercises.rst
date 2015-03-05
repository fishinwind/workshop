
.. include:: /_static/substitutions.txt 

.. _exercises:

*************
  Exercises
*************

:Last updated: |today|

Exercises for the Genome Analysis Workshop (MOLB 7621).

.. _unix-exercises:

UNIX exercises
==============

#. Use ``wc`` to determine how many *lines* are in `lamina.bed`.

#. Use ``wc`` to determine how many *characters* are in `lamina.bed`.

#. Print out the 10 states (1st column, contains spaces) with the highest
   income (3rd column) from `states.tab`.

.. code-block:: bash

    $ wc -l lamina.bed # lines

    $ wc -c lamina.bed # characters

    $ sort -t'\t' -k3nr states.tab | head -n 10

.. _awk-exercises:

AWK exercises
=============

.. _r-exercises:

R exercises
===========

R exercises in Rmarkdown format (and rendered HTML) are here:

- Rmarkdown Example :download:`[Rmd] </RStudio/Rmarkdown_example.Rmd>`
  :download:`[HTML] </RStudio/Rmarkdown_example.html>`

- Rmarkdown Excercises :download:`[Rmd] </RStudio/Exercises.Rmd>`
  :download:`[HTML] </RStudio/Exercises.html>`

.. _python-exercises:

Python excercises
=================
