
========================
Genome Analysis Workshop 
========================

:Title: Genome Analysis Workshop
:Course Number: MOLB 7621
:Semester: Spring 2015
:Homepage: http://hesselberthlab.github.io/workshop
:Author: Jay R. Hesselberth
:Organization: University of Colorado School of Medicine
:Address: Department of Biochemistry and Molecular Genetics
:Copyright: 2013-2015 Jay R. Hesselberth

Making content 
--------------
These are the ``make`` targets for content:

.. code-block:: makefile

  pdf        to make standalone PDF
  html       to make standalone HTML files
  slides     to make standalone HTML5 slides
  all        to test make html, slides and PDF
  publish    push website to github gh-pages
  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter

Content will be in ``_build/<target>`` e.g. ``_build/html``.

Making webiste
--------------
Push content to the github.io website with:

.. code-block:: bash

    $ make publish

