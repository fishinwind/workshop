README for Genome Informatics Workshop Sphinx docs
==================================================

Making HTML and slides
----------------------

Make slides with::

    make slides

Make HTML with::

    make html

Move HTML to bioworkshop.github.io with::

    websitedir=$HOME/devel/bio-workshop-io # change if necessary
    make html
    cp -r _build/html/* $websitedir
    cd $websitedir
    git pull origin master
    git add . # add all new html

    # use git status to check all is added.
    # you will have to explicity add (and use git add -f) stuff in
    # section_1, so that we don't post draft classes early.
    git commit -m "update html"
    git push origin master

Navigate to http://ucd-bioworkshop.github.io/ Pages are cached at this
site, so it will take a few minutes to update

