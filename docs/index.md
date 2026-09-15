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

## Known gaps

- Links from the included pages into the source tree (for example
  `calculations/defective_metal.py#L280`) are GitHub-relative and do not resolve
  in this HTML build. They are suppressed via `suppress_warnings` in `conf.py`;
  the fix is to rewrite them as absolute `github.com` URLs, which work in both
  renderers.
- The notebooks (`Application/*.ipynb` and the three study directories) and the
  LaTeX derivation `latex/Model_Equations.tex` are intentionally not included
  here.
- `README.md` and Appendix C of `TRAPPING_VALIDATION.md` still reference the
  removed `data/` package; each carries a warning banner until they are rewritten.
- The docstring section headers in `oxide_permeation`, `defective_metal`,
  `interface_solver` and `parallel_oxide_defect_paths` use a non-standard
  `Parameters:` form; those four pages are pending cleanup.
