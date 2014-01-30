makeslides with::

   rst2slide input.rst output.html

to get new html to bioworkshop.io. From this rst directory:

    make html
    cp -r _build/html/*  ../../bio-workshop-io/
    cd ../../bio-workshop-io/
    git pull origin master
    git add . # add all new html

    # use git status to check all is added.
    # you will have to explicity add (and use git add -f) stuff in
    # section_1, so that we don't post draft classes early.
    git commit -m "update html"
    git push origin master


    # then go to: http://ucd-bioworkshop.github.io/ . It will take a few minutes to update
     
