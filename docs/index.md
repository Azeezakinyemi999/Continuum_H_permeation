# Hierarchical Hydrogen Permeation Model

A steady-state analytical model for hydrogen permeation through an oxide-coated
structural alloy wall, built as a hierarchy of increasing physical complexity
(Levels 1–6). Given a material, an operating temperature and a pressure
difference, it predicts the hydrogen flux through the wall and — more usefully —
says **which physical process is limiting that flux**.

Every level is closed-form and steady-state. There is no time axis; see
{doc}`STEADY_STATE` for why that is sufficient here and what would
invalidate it.

## Where to start

- **New to the model** — read {doc}`getting-started`, then the theory chapters
  in order. They build up one physical effect at a time.
- **Looking up a parameter or a returned quantity** — {doc}`PARAMETERS`
  and {doc}`OUTPUTS`.
- **Running an analysis** — the sensitivity-analysis how-to (in progress) for the
  regime-stratified sensitivity workflow, {doc}`how-to/switch-study` to change
  material system.
- **Reading the code** — {doc}`api/index`, generated from the
  docstrings.

```{toctree}
:maxdepth: 2
:caption: Getting started

getting-started
```

```{toctree}
:maxdepth: 2
:caption: Theory

theory/foundations
theory/equilibrium-models
STEADY_STATE
theory/level1-metal
theory/level2-oxide
theory/level3-defective-oxide
theory/level4-microstructure
```

```{toctree}
:maxdepth: 2
:caption: How-to

how-to/switch-study
CONTAINER
```

```{toctree}
:maxdepth: 2
:caption: Reference

PARAMETERS
OUTPUTS
api/index
```

```{toctree}
:maxdepth: 2
:caption: Validation and release

TRAPPING_VALIDATION
PACKAGING_1.0.0
```

```{note}
This documentation is mid-rewrite. The theory chapters and how-to guides are
being authored from the project's earlier design notes, with every claim checked
against the current code; the pages still carrying their original uppercase
filenames are the ones not yet absorbed. Source material awaiting absorption sits
in `docs/_source_material/` and is excluded from the build.
```

## How this documentation is kept honest

Documentation about a model that has moved on is worse than none, so the factual
claims here are machine-checked rather than trusted:

- `docs/_tools/verify_docs.py` fails if a page names a function, module or
  parameter that does not exist, or links to a source file with a
  repo-relative path.
- `docs/_tools/worked_values.py` regenerates every number quoted in the theory
  chapters by running the current model. `--check` diffs against the stored copy,
  so a model change surfaces as a diff rather than as quietly wrong prose.
- The Sphinx build runs with `-W`, so a malformed docstring section or a broken
  cross-reference fails it.

Numbers in the theory chapters are therefore reproducible by construction. Where
a chapter quotes a flux or a coverage, it is the value the code returns today for
the active study — not a figure carried over from an earlier configuration.
