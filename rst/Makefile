# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = 
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build
SRCDIR        = ../src
KEYS_DIR      = problem-set-keys
PROBSET_DIR   = Problem_Sets
KEYS_TAR      = $(PROBSET_DIR)/problem-set-keys.tar

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
	@echo "  all        to create all class files and post to website"
	@echo "  pdf        to make standalone PDF"
	@echo "  html       to make standalone HTML files"
	@echo "  slides     to make standalone HTML5 slides"
	@echo "  website    to copy and push website to github"
	@echo "  keys       to create tarball of problem set keys"

clean:
	rm -rf $(BUILDDIR)/*

all: clean pdf html slides website

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
	@echo

slides:
	$(SPHINXBUILD) -b slides $(ALLSPHINXOPTS) $(BUILDDIR)/slides
	@echo
	@echo "Build finished. The HTML slides are in $(BUILDDIR)/slides."
	@echo

website:
	@echo
	. $(SRCDIR)/website.sh; 

pdf:
	$(SPHINXBUILD) -b pdf $(ALLSPHINXOPTS) $(BUILDDIR)/pdf
	@echo
	@echo "Build finished. The PDF is in $(BUILDDIR)/pdf."
	@echo

keys:
	@cd .. && tar cvf $(KEYS_TAR) $(KEYS_DIR)
	@echo "Problem Set Keys are available in $(KEYS_TAR)"
	@echo

