# Minimal makefile for Sphinx documentation.
#
# -W makes warnings fatal. The build is clean as of the docstring-markup and
# absolute-link commits, so any new warning is a real regression: a malformed
# docstring section, a relative link to a source file, or a page missing from a
# toctree. --keep-going reports every one instead of stopping at the first.
#
# Deliberately no -n (nitpicky): many functions document their types as prose
# ("dict containing ..."), which would bury the real signal in unresolvable
# type references.
#
# To triage without failing, override on the command line:
#     make html SPHINXOPTS=

SPHINXOPTS    ?= -W --keep-going
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
