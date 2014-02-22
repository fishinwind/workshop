***********
Class Notes
***********

Notes for Spring 2015
=====================

Merge awk class into Linux command line classes, cover grep much earlier,
and trim this to squeeze in more Python

Students struggled with problem set 3 (the first python problem set).
we failed to go into detail on nested data structures (e.g. dict of
lists of tuples). we also probably should have gone over multi-line
data formats a bit more prior to the problem set.

Next year, should focus less on trivial examples (e.g. fruit examples) and
spend more time going over biologic data.

     #. create nested data types during loops
     #. logic for determining state during looping
     #. functions to act on data structures once created:

        - min(), max()
        - Counter.most_common()

In class 10 Dont included the extra lines in `sample-seq-info.csv` it
was too confusing to use the skip_while arg to toolshed.reader to pass a
function

Name expansion seems to be confusing. one slide should have:

.. code-block:: python

    # these are all the same ...

    # explicit but cluttered
    for tuple in zip(list1, list2):
        item1 = tuple[0]
        item2 = tuple[1]
        ...

    # better, parens not needed
    for (item1, item2) in zip(list1, list2):
        ...

    # best
    for item1, item2 in zip(list1, list2):
        ...

Sphinx docs for for Genome Informatics Workshop
===============================================

Making linked slides and HTML
-----------------------------
Make linked slides with:
    
.. code-block:: bash

    $ make html slides

Making slides
-------------
Make slides with:

.. code-block:: bash

    $ make slides

Slides will be in ``_build/slides``

Making webiste
--------------
Make HTML only with:

.. code-block:: bash

    $ make html

.. todo::

   - fix the git add -f hack for posting new content
    
Move HTML to bioworkshop.github.io with:

.. code-block:: bash

    $ websitedir=$HOME/devel/UCD-BioWorkshop.github.io
    $ cp -r _build/html/* ../../bio-workshop-io/ $websitedir
    $ cp -r _build/slides/ ../../bio-workshop-io/ $websitedir
    $ cd $websitedir

    $ git pull origin master
    # add all new html
    $ git add . 

    # use git status to check all is added.
    # you will have to explicity add (and use git add -f) stuff in
    # Block_1, so that we don't post draft classes early.
    git commit -m "update html"
    git push origin master

.. note::

    this is codified in src/website.sh and can be run with:

    $ make website

Check the new website http://ucd-bioworkshop.github.io/

Pages are cached at this site, so it will take a few minutes to update.

