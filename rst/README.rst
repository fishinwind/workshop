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

Slides will be in:
     
.. code-block:: bash

    _build/slides

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
    $ cp -r _build/html/* $websitedir
    $ cp -r _build/slides/* $websitedir
    $ cd $websitedir

    $ git pull origin master
    # add all new html
    $ git add . 

    # use git status to check all is added.
    # you will have to explicity add (and use git add -f) stuff in
    # Block_1, so that we don't post draft classes early.
    git commit -m "update html"
    git push origin master

Check the new website http://ucd-bioworkshop.github.io/

Pages are cached at this site, so it will take a few minutes to update.

