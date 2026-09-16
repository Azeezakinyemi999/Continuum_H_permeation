# Hierarchical Hydrogen Permeation Model

Steady-state analytical models for hydrogen permeation through an oxide-coated
structural alloy wall, built as a hierarchy of increasing physical complexity
(Levels 1–6). Which material system is modelled is selected by a single line in
{py:mod}`calculations.config.model_config`.

```{toctree}
:maxdepth: 2
:caption: Overview

README
STEADY_STATE
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
:caption: Validation

TRAPPING_VALIDATION
```

```{toctree}
:maxdepth: 1
:caption: Release and distribution

CONTAINER
PACKAGING_1.0.0
```

## Scope and conventions

- **This build is strict.** `docs/Makefile` passes `-W --keep-going`, so a
  malformed docstring section or a repo-relative link to a source file fails it.
  Links from these pages into the source tree are absolute `github.com` URLs on
  the `main` branch, which resolve both here and on GitHub. Keep them that way.
- **Links point at `main`.** Work merged to `main` after a page was written will
  be reflected; work still on a feature branch will not.
- Not included by design: the notebooks (`Application/*.ipynb` and the three
  study directories, whose outputs are committed and large) and the LaTeX
  derivation `latex/Model_Equations.tex`.
- Partially documented: `surface_kinetics` writes some parameter blocks in a
  compressed one-line form that napoleon reads as a type rather than a
  description, so those entries render less well than the rest.
