# Minimal makefile for Sphinx documentation.
#
# Once the docstring cleanup is finished, add `-W --keep-going` to SPHINXOPTS to
# make warnings fatal. Not before: today that would fail on the malformed
# docstring section headers in calculations/.

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
