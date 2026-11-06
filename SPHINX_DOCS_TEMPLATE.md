# Markdown → Sphinx → GitHub Pages: a reusable template

A generic recipe, extracted from setting this up on the MHI permeation model. Nothing
here is specific to that project — placeholders are in `<ANGLE BRACKETS>`.

**What it gives you.** Existing `.md` files render as a documentation site *without
moving them* out of the repo root where people already read them on GitHub. Docstrings
become an API reference. Warnings are fatal, so the build doubles as a test.

The traps at the end are the part worth reading twice. Every one of them cost real time.

---

## Layout

```
<repo>/
├─ README.md, DESIGN.md, …          your existing .md — they stay where they are
├─ pyproject.toml                   docs extra lives here
├─ .gitignore                       docs/_build/, docs/.mplcache/
└─ docs/
   ├─ conf.py
   ├─ Makefile
   ├─ index.md                      landing page + toctrees
   ├─ DESIGN.md                     3-line stub pointing at ../DESIGN.md
   ├─ _static/                      images; add a .gitkeep
   ├─ _templates/
   ├─ api/
   │  ├─ index.md                   toctree of module pages
   │  └─ <module>.md                one automodule stub per module
   └─ _tools/                       excluded from the build
```

---

## 1. The stub pattern — the key move

Do not copy or move your root `.md` files. Give each one a three-line stub in `docs/`:

````markdown
```{include} ../DESIGN.md
```
````

One source of truth. GitHub renders the original; Sphinx renders the same bytes. Add
`DESIGN` to a toctree in `index.md` and it appears in the sidebar.

---

## 2. `docs/conf.py`

```python
"""Sphinx configuration for <PROJECT>.

Build with:  cd docs && make html
Requires:    pip install -e ".[docs]"
"""
import os
from importlib.metadata import version as _dist_version

# --- import-time side effects, neutralised BEFORE autodoc imports anything ----
# conf.py runs first, which is what makes this work. setdefault, not assignment,
# so an explicit env var on the command line still wins.
os.environ.setdefault("MPLBACKEND", "Agg")          # if any module imports pyplot
os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(os.path.dirname(__file__), ".mplcache"))
# Any env var your package needs in order to import at all goes here too:
# os.environ.setdefault("<APP>_CONFIG", "<default>")

project   = "<PROJECT>"
author    = "<AUTHOR>"
copyright = "<YEAR>, <AUTHOR>"
release   = _dist_version("<DISTRIBUTION-NAME>")   # from pyproject; loud if not installed
version   = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",      # docstrings -> reference
    "sphinx.ext.napoleon",     # NumPy / Google docstring styles
    "sphinx.ext.viewcode",     # [source] link per function
    "sphinx.ext.intersphinx",  # numpy.ndarray etc. become links
    "sphinx.ext.mathjax",      # $$ math in your .md
    "myst_parser",             # REQUIRED: lets Sphinx read .md at all
]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "_tools"]
templates_path = ["_templates"]

# --- MyST ---------------------------------------------------------------------
myst_enable_extensions = [
    "dollarmath",    # $...$ and $$...$$
    "colon_fence",   # :::{note} — avoids nested-backtick pain inside stubs
    # "amsmath",     # only if you use \begin{align} at block level
]
myst_heading_anchors = 3        # stable #slug per heading -> deep links survive
myst_dmath_double_inline = True
myst_footnote_transition = False

# --- napoleon -----------------------------------------------------------------
napoleon_google_docstring = False   # pick exactly ONE; see trap 3
napoleon_numpy_docstring  = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True

# Any docstring section header that is NOT a napoleon built-in must be listed
# here, or docutils reads it as a section title inside a directive and errors.
# This list is load-bearing, not decorative.
napoleon_custom_sections = [
    ("Theory", "notes"),
    ("Limitations", "notes"),
    # ("<Your Header>", "notes"),
]

# --- autodoc ------------------------------------------------------------------
autodoc_member_order = "bysource"   # alphabetical destroys a deliberate order
autodoc_typehints = "none"          # "signature" if you annotate
autodoc_preserve_defaults = True    # keeps seed=42 readable instead of repr-mangled
autodoc_inherit_docstrings = False
autodoc_default_options = {"members": True, "show-inheritance": False}
# Do NOT enable undoc-members globally — see trap 4.

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy":  ("https://numpy.org/doc/stable/", None),
}
intersphinx_timeout = 5             # fail fast when building offline

nitpicky = False                    # True drowns you if types are documented as prose

html_theme = "furo"                 # persistent sidebar, search, dark mode
html_static_path = ["_static"]
html_title = f"{project} {release}"
```

---

## 3. `docs/Makefile`

```make
# -W makes warnings fatal: the build becomes a test.
# --keep-going reports every warning instead of stopping at the first.
# To triage without failing:  make html SPHINXOPTS=
SPHINXOPTS    ?= -W --keep-going
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

.PHONY: help Makefile
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)
```

---

## 4. `docs/index.md`

````markdown
# <PROJECT>

<One paragraph: what it does and who it is for.>

## Where to start

- **New here** — {doc}`getting-started`
- **Reading the code** — {doc}`api/index`

```{toctree}
:maxdepth: 2
:caption: Guide

getting-started
DESIGN
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/index
```
````

---

## 5. `docs/api/` — one stub per module

`docs/api/<module>.md`:

````markdown
# `<module>` — <one line on what it covers>

```{eval-rst}
.. automodule:: <package>.<module>
   :members:
   :member-order: bysource
```
````

`docs/api/index.md` is a toctree listing them in whatever order makes sense. Concept
order usually beats alphabetical.

---

## 6. Dependencies — `pyproject.toml`

```toml
[project.optional-dependencies]
docs = [
    "sphinx>=7.2,<8",      # <8 ONLY if you still support Python 3.9
    "myst-parser>=2,<3",   # <3 likewise
    "furo>=2024.1",
]
```

Those ceilings are **Python-version ceilings**, not caprice: Sphinx 8 and myst-parser 3
both require Python >= 3.10. On 3.10+, drop the upper bounds.

---

## 7. `.gitignore`

```gitignore
docs/_build/
docs/.mplcache/

# If you blanket-ignore images, re-admit the docs ones:
# !docs/_static/**/*.png
# !docs/_static/**/*.svg
```

---

## 8. Publish to GitHub Pages

`gh-pages` holds the **built HTML at its root**. Build from a committed state, in a
worktree, so your working tree is never involved:

```bash
git worktree add --detach /tmp/src <COMMIT>
(cd /tmp/src/docs && make html)

git worktree add /tmp/pages gh-pages          # orphan branch, created once
rsync -a --delete --exclude='.nojekyll' --exclude='.git' \
      /tmp/src/docs/_build/html/ /tmp/pages/

cd /tmp/pages && git add -A && git commit -m "Publish docs at <COMMIT>" && git push
git worktree remove /tmp/src && git worktree remove /tmp/pages
```

Creating the orphan branch, once:

```bash
git checkout --orphan gh-pages
git rm -rf .
touch .nojekyll
git add .nojekyll && git commit -m "Initialise gh-pages"
git push -u origin gh-pages
git checkout main
```

**`.nojekyll` at the gh-pages root is mandatory, and Sphinx does not emit it.** Without
it, GitHub Pages ignores every underscore-prefixed directory: `_static/` disappears and
the site renders with no CSS at all. Create it once, and `--exclude` it from every sync
afterwards. This is the easiest way to ship a silently broken docs site.

---

## The traps, in the order they bite

**1. `make html` is the only thing that catches docstring RST.** A custom checker that
validates identifiers and return keys will pass happily while a malformed docstring
breaks the build. Mine did exactly that, and the docs were unbuildable for a week
because I ran the checker and not the build. **Put `make html` in CI.**

**2. Wrapped bullets need a blank line after them.** This is invalid RST:

```
- an option whose description
  wraps onto a second line
Default: 'x'
```

Put a blank line before `Default:`, or lead with a paragraph and put the list below it.
Single-line bullets tolerate the missing blank line; wrapped ones do not. The error
surfaces as `Bullet list ends without a blank line; unexpected unindent`.

**3. Enable exactly one napoleon style.** Leaving both `napoleon_google_docstring` and
`napoleon_numpy_docstring` on means the Google parser claims malformed NumPy blocks and
emits a parameter literally named `-----------` with the real type as its description.
Plausible-looking HTML that is wrong. A loud docutils error beats a quiet lie.

**4. Do not enable `undoc-members` globally.** If any module does
`globals().update(...)` at import, that page becomes a wall of undocumented data plus a
stray `np`.

**5. Relative links to source files break.** `[x](../src/thing.py)` is correct on GitHub
and 404s in the built site. Use absolute `https://github.com/<org>/<repo>/blob/main/...`
URLs, which work in both renderers.

**6. Neutralise import-time side effects in `conf.py`, not in the package.** `conf.py`
executes before autodoc imports anything, which is the only place an `MPLBACKEND`
setdefault lands in time.

**7. Keep caches out of `_build/`.** `make clean` deletes `_build/`, so a font cache
living there is rebuilt on every clean build.

**8. A raw docstring is required for LaTeX.** In a non-raw string Python turns `\frac`
into a form feed and `\right` into a carriage return before Sphinx ever sees them. The
build can still succeed with silently mangled maths. Use `r"""`.

---

## Optional: a doc verifier

Worth adding once the docs quote identifiers. Roughly 150 lines that walk every `.md`
and check:

- every `` `name(` `` resolves to a real function in the package
- every `` `pkg.module` `` exists
- no link points at a source file by a relative path

Exit 1 on findings. Give it an escape hatch for deliberate mentions of things that do
not exist:

```markdown
<!-- verify-docs: allow old_function_name -->
```

**It is not a substitute for `make html`** — it checks names, not markup. Run both.

---

## Bootstrap order

1. `mkdir -p docs/{_static,_templates,api,_tools}`; add `docs/_static/.gitkeep`
2. Drop in `conf.py` and `Makefile`; fill the placeholders
3. `pip install -e ".[docs]"`
4. Write `index.md` with a single toctree
5. `make html` — fix until clean
6. Add stubs for root `.md` files **one at a time**, rebuilding after each
7. Add `api/` stubs **module by module** — this is where napoleon errors surface
8. Create the orphan `gh-pages` branch with its `.nojekyll`, then publish
