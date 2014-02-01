Sphinx docs for for Genome Informatics Workshop
===============================================

Making slides
-------------

Make slides with::

    make slides

Slides will be in::
    
    $ _build/slides

Making webiste
--------------

Make HTML with::

    make html

.. todo::

   - fix the git add -f hack for posting new content
    
Move HTML to bioworkshop.github.io with::

    websitedir=$HOME/devel/UCD-BioWorkshop.github.io
    cp -r _build/html/* $websitedir
    cd $websitedir

    git pull origin master
    git add . # add all new html

    # use git status to check all is added.
    # you will have to explicity add (and use git add -f) stuff in
    # section_1, so that we don't post draft classes early.
    git commit -m "update html"
    git push origin master

Check the new website http://ucd-bioworkshop.github.io/

Pages are cached at this site, so it will take a few minutes to update.

