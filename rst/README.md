makeslides with::

   rst2slide input.rst output.html

to get new html to bioworkshop.io:

    make html
    cp -r _build/html/  ../../bio-workshop-io/
    cd ../../bio-workshop-io/
    git pull origin master
    git add . # add all new html
    git commit -m "update html"
    git push origin master

    # then go to: http://ucd-bioworkshop.github.io/ . It will take a few minutes to update
     
