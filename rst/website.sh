#! /usr/bin/env bash
#
# install script to push workshop html and slides to github.io site
# 
# XXX: confirm this is working
#
set -o nounset -o pipefail -o errexit -x

WEBSITEDIR=$HOME/devel/UCD-BioWorkshop.github.io
TARGETIDR=../../bio-workshop-io/
BUILDDIR=_build

echo ">> copying files to $WEBSITEDIR ..."
cp -r $BUILDDIR/html/* $TARGETDIR $WEBSITEDIR
cp -r $BUILDDIR/slides/* $TARGETDIR $WEBSITEDIR

echo ">> syncing files with git ..."
# sync the directory and add new content
cd $WEBSITEDIR
git pull origin master
git add .
git status

# use git status to check all is added.
date=$(date "+%Y-%m-%d at %T")
git commit -m "html and slide update on $date"
git push origin master

