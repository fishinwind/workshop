#! /usr/bin/env bash
#
# Script to install script to push workshop html
# and slides to github.io site.
# 
# Contact: jay hesselberth at gmail.com

set -o nounset
set -o pipefail
set -o errexit
set -x

WEBSITEDIR=$HOME/devel/hesselberthlab.github.io/workshop
SRCDIR=$HOME/devel/genomics-workshop
RSTDIR=$SRCDIR/content
BUILDDIR=$RSTDIR/_build

if [[ ! -d WEBSITEDIR ]]; then
    mkdir -p $WEBSITEDIR
fi

echo ">> pulling content ..."
cd $WEBSITEDIR
git pull origin master
echo

echo ">> copying files to $WEBSITEDIR ..."
cd $RSTDIR
cp -r $BUILDDIR/html/* $WEBSITEDIR
cp -r $BUILDDIR/slides/ $WEBSITEDIR
echo

echo ">> pushing content ..."
# sync the directory and add new content
cd $WEBSITEDIR
git add .
git status
echo

# use git status to check all is added.
date=$(date "+%Y-%m-%d at %T")
echo "

Website pushed. check and finish with:
cd $WEBSITEDIR
git commit -m \"html and slide update on $date\"
git push origin master

"

