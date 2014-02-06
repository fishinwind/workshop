#! /usr/bin/env bash
#
# install script to push workshop html and slides to github.io site
# 
# XXX: confirm this is working
#
set -o nounset -o pipefail -o errexit -x

WEBSITEDIR=$HOME/devel/bio-workshop-io
SRCDIR=$HOME/devel/bio-workshop
RSTDIR=$SRCDIR/rst
BUILDDIR=$RSTDIR/_build

cd $RSTDIR

echo ">> updating content ..."
cd $WEBSITEDIR
git pull origin master
cd $RSTDIR
echo

echo ">> copying files to $WEBSITEDIR ..."
cp -r $BUILDDIR/html/* $WEBSITEDIR
cp -r $BUILDDIR/slides/ $WEBSITEDIR
echo

echo ">> syncing files with git ..."
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

