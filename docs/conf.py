"""Sphinx configuration for the hierarchical hydrogen permeation model.

Build with:

    conda activate mace_env
    cd docs && make html

Requires ``pip install -e ".[docs]"`` from the repo root. There is deliberately no
``sys.path.insert`` here — PACKAGING_1.0.0.md records several incidents caused by
sys.path manipulation, and this file will not be the next one.
"""

import os
from importlib.metadata import version as _dist_version

# --- import-time side effects of the package, neutralised before autodoc runs ---

# calculations/sensitivity.py imports matplotlib.pyplot at MODULE level, so every
# autodoc run instantiates a backend. On macOS that resolves to the interactive
# `macosx` backend, which wants a window server. conf.py is executed before any
# autodoc import, which is what makes this work. setdefault, not assignment, so an
# explicit `MPLBACKEND=... make html` still wins.
os.environ.setdefault("MPLBACKEND", "Agg")

# The same import triggers matplotlib's font-cache build. Keep it under docs/ so
# the docs never fail on an unwritable $HOME (same reasoning as
# MPLCONFIGDIR=/tmp/.matplotlib in container/Dockerfile) — but deliberately NOT
# inside _build/, which `make clean` deletes: that would rebuild the font cache
# on every clean build, costing ~15 s and emitting matplotlib log records.
os.environ.setdefault(
    "MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), ".mplcache")
)

# --- project metadata ---------------------------------------------------------

project = "Hierarchical Hydrogen Permeation Model"
author = "Azeez Akinyemi"
copyright = "2026, Azeez Akinyemi"

# Single source of truth is pyproject.toml's `version`, read back from installed
# metadata. PackageNotFoundError here means `pip install -e .` was skipped — a
# loud failure is correct, because without the install autodoc imports nothing.
release = _dist_version("continuum-h-permeation")
version = ".".join(release.split(".")[:2])

# --- extensions ---------------------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",      # the whole point: ~89 public functions, docstrings already written
    "sphinx.ext.napoleon",     # those docstrings are NumPy-style; without this each is one blob
    "sphinx.ext.viewcode",     # [source] link per function — the physics lives in the code
    "sphinx.ext.intersphinx",  # makes numpy.ndarray / scipy.optimize.brentq clickable
    "sphinx.ext.mathjax",      # renders the $$-blocks in TRAPPING_VALIDATION.md and STEADY_STATE.md
    "myst_parser",             # without it Sphinx cannot read the GitHub-flavoured .md files at all
]

# --- source discovery ---------------------------------------------------------

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
exclude_patterns = [
    "_build", "Thumbs.db", ".DS_Store",
    # Tooling, not content: the doc verifier and the worked-values generator.
    "_tools",
    # Legacy documents awaiting absorption into theory/ and how-to/. They are
    # kept in-tree so the rewrite can be diffed against them, but they are NOT
    # part of the site — they contain the stale material the rewrite replaces.
    # This directory should be empty (and this entry removed) when the rewrite
    # is finished.
    "_source_material",
]
templates_path = ["_templates"]

# --- MyST ---------------------------------------------------------------------

myst_enable_extensions = [
    "dollarmath",   # required: $...$ / $$...$$ is how all the root .md files write math
    "colon_fence",  # lets the stub pages use :::{warning} without nested-backtick pain
]
# TRAPPING_VALIDATION.md uses \begin{cases} inside $$, which MathJax handles
# natively; there is no block-level \begin{align}, so `amsmath` is not needed.

myst_heading_anchors = 3         # stable #slug per heading, so cross-file deep links keep working
myst_dmath_double_inline = True  # tolerate a $$...$$ that sits on one line
myst_footnote_transition = False

# --- napoleon -----------------------------------------------------------------

# FALSE on purpose. Every docstring here is NumPy-INTENT. Leaving the Google pass
# on makes it claim the malformed `Parameters:`-with-colon blocks and emit a
# parameter literally named "-----------" with the real type string as its
# description — plausible-looking HTML that is wrong. A loud docutils error is
# better than a quiet lie in a scientific reference.
napoleon_google_docstring = False
napoleon_numpy_docstring = True

napoleon_include_init_with_doc = False     # no classes in the package; explicit for future-proofing
napoleon_include_private_with_doc = False  # keeps sensitivity.py's _helpers out of the reference
napoleon_include_special_with_doc = False

napoleon_use_param = True         # one :param:/:type: per argument -> a definition list with units
napoleon_use_rtype = True         # these functions return dicts; a separate line reads better
napoleon_use_ivar = False         # no classes
napoleon_preprocess_types = True  # normalises bare `float`/`dict` so intersphinx can resolve them
napoleon_attr_annotations = True

napoleon_use_admonition_for_examples = True    # boxes the `>>>` blocks so they stop running into prose
napoleon_use_admonition_for_notes = True       # ditto the physics Notes sections
napoleon_use_admonition_for_references = True  # nearly every function cites DOIs

# Non-standard sections this codebase actually uses, mapped onto rendering styles.
# Registering them here means the docstring cleanup is a HEADER rename only — the
# content never has to be rewritten into a `Notes` section.
# Anything NOT registered here and not one of napoleon's built-in section names
# is left unparsed, and docutils then reads the header + underline as an RST
# section title inside a directive body — a CRITICAL "Unexpected section title".
# So this list is load-bearing, not decorative.
napoleon_custom_sections = [
    ("Theory", "notes"),
    ("Mathematical Derivation", "notes"),
    ("Mathematical Model", "notes"),
    ("Mathematical Framework", "notes"),
    ("Physics Note", "notes"),
    ("Limitations", "notes"),
    ("Typical ranges", "notes"),
    ("Module Structure", "notes"),
    ("Physical Parameters", "notes"),
    ("Usage", "notes"),
]

# --- autodoc ------------------------------------------------------------------

autodoc_member_order = "bysource"  # the modules are ordered by model Level; alphabetical would destroy that
autodoc_typehints = "none"         # there are no annotations in the package; "signature" would add noise
autodoc_preserve_defaults = True   # keeps `seed=42`, `figsize=(9, 6.5)` readable instead of repr-mangled
autodoc_inherit_docstrings = False

# Do NOT enable undoc-members or imported-members globally. model_config.py does
# globals().update(_study_ns) at import, injecting ~15 dicts AND numpy as `np`;
# with undoc-members that page becomes a wall of undocumented data plus a stray
# `np`. Default autodoc skips both (no docstring, foreign __module__).
autodoc_default_options = {
    "members": True,
    "show-inheritance": False,
}
autodoc_mock_imports = []  # everything importable is a real dependency; mocking would hide a broken install

# --- intersphinx --------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}
intersphinx_timeout = 5  # fail fast when building offline; the build still succeeds, with warnings

# --- warnings -----------------------------------------------------------------

# No suppression. The root .md files used to link into the source tree with
# GitHub-relative paths, which are correct on GitHub and unresolvable from
# docs/; they are now absolute github.com URLs, which work in both renderers.
# Keep it that way — a relative link to a source file will now fail the build.

nitpicky = False  # many functions document types as prose; nitpicky would drown the real signal

# --- HTML ---------------------------------------------------------------------

html_theme = "furo"  # good default for a function-heavy reference: persistent sidebar, search, dark mode
html_static_path = ["_static"]
html_title = f"{project} {release}"
html_show_sourcelink = True
