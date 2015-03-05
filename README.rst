
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

Organization
------------

`content/` has all of the active content.

`reorg/` is used for oragnizing and renaming files. It is outside of the
sphinx root and will not be used during build.

Updating new dates
------------------

Dates are stored in `_static/substituions.txt` and are referred to via
sphinx substitutions in the content. Update these and the changes will
propagate site-wide.

Making changes to the material
------------------------------

Everybody in the `workshop` group in Github has read access, but cannot
write directly to the `workshop` repo. Changes will be done via pull
requests.

Workflow
~~~~~~~~

- Fork the repository on the Github website

- Clone the repo into your home directory

- Create a new branch for your work::
    
    $ git branch new-class
    $ git checkout new-class

- Make changes, commit them to the new branch.

- Push your changes on the new branch::

    $ git push origin new-class

- Initiate a pull request on the Github website.


