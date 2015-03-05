# Makefile for Sphinx documentation
#

# the `branch` and `push` commands wipe the gh-pages branch, then the
# `ghp-import` command creates an entirely new copy
#
# $ git branch -D gh-pages
# $ git push origin :gh-pages

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don't have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo
	@echo "  rmarkdown  to make Rmarkdown documents"
	@echo "  pdf        to make standalone PDF"
	@echo "  html       to make standalone HTML files"
	@echo "  slides     to make standalone HTML5 slides"
	@echo "  all        to test make html, slides and PDF"
	@echo "  publish    push website to github gh-pages"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"

clean:
	rm -rf $(BUILDDIR)/*
	cd RStudio && $(MAKE) clean

html: 
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
	@echo

slides:
	$(SPHINXBUILD) -b slides $(ALLSPHINXOPTS) $(BUILDDIR)/html/slides
	@echo
	@echo "Build finished. The HTML slides are in $(BUILDDIR)/html/slides."
	@echo

pdf:
	$(SPHINXBUILD) -b pdf $(ALLSPHINXOPTS) $(BUILDDIR)/html/_downloads
	@echo
	@echo "Build finished. The PDF is in $(BUILDDIR)/html/_downloads."
	@echo

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
          "run these through (pdf)latex."

rmarkdown:
	cd RStudio && $(MAKE) clean && $(MAKE) all

all: rmarkdown html slides pdf 

publish: all 
	ghp-import -n -p ./_build/html
	@echo
	@echo "HTML output published on github-pages"

