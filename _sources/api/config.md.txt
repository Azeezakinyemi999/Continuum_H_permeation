# `config` — study selection

## How a study is selected

Every module in `calculations/` imports its material data from
`calculations.config.model_config`, which re-exports exactly one module from
`calculations/config/studies/`. Switching study is the single line below, and no
other file changes:

```{literalinclude} ../../calculations/config/model_config.py
:language: python
:lines: 40-46
:caption: calculations/config/model_config.py
```

The re-export is done by hand rather than with `from studies.X import *`, because
a star import's module path is resolved at compile time and cannot be a variable.
`_load_active_study` reproduces star-import semantics exactly, then validates the
result against a 15-name `REQUIRED_EXPORTS` frozenset — raising `ValueError`
naming both the study and the gap, rather than letting a missing export surface as
an `ImportError` from deep inside `calculations.sensitivity`.

After switching, run `calculations.sensitivity.check_against_config()`: preset
yields, sweep ranges and draw counts were tuned for a specific material and are
not automatically valid for another one.

```{note}
`model_config` populates its own namespace with `globals().update()` at import
time, so `METALS`, `OXIDES`, `CONDITIONS`, `DEFAULT_PARAMS_LEVEL5L6` and the rest
are **not statically defined there** and therefore do not appear below. For their
values see {doc}`../PARAMETERS`, which is generated against the active study.
```

## `calculations.config.model_config`

```{eval-rst}
.. automodule:: calculations.config.model_config
   :members:
   :member-order: bysource
```

## `calculations.config.studies`

```{eval-rst}
.. automodule:: calculations.config.studies
```

### Available studies

| Module | System |
|---|---|
| `Guo_etal_2025_316L` | 316L / Cr₂O₃ — **active** |
| `incoloy802_cr2o3` | Incoloy 802 (X40 NiCrAlTi 31/19) / Cr₂O₃ |
| `fuerst_etal_2024_model_config` | Hastelloy N |

The three modules are structurally identical and define the same interface, so
only the active one is documented here. Note that a per-directive `:members:`
list would be silently replaced by `autodoc_default_options["members"] = True`
in `conf.py`, so this documents all public functions rather than a subset.

```{eval-rst}
.. automodule:: calculations.config.studies.Guo_etal_2025_316L
   :members:
   :member-order: bysource
```
