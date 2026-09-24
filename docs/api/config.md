# `config` — study selection

## How a study is selected

Every module in `calculations/` imports its material data from
`calculations.config.model_config`, which loads exactly one study. Which one is
decided by **location**, first hit wins:

1. `$MHI_STUDY` — a path to a study file, resolved from the repository root.
2. `study_config.py` in the working directory, then each parent up to the
   repository root. This is how a notebook gets the study sitting in its own
   folder.
3. Nothing else — there is no default. A run that named no study raises, listing the
   study folders it could have meant.

Each study folder's `study_config.py` **is** the study: the material properties,
operating conditions, sweep ranges and presets live there, not in the package.

```{literalinclude} ../../calculations/config/model_config.py
:language: python
:start-at: "# WHERE STUDIES LIVE"
:end-at: "STUDY_ENV_VAR = 'MHI_STUDY'"
:caption: calculations/config/model_config.py
```

The study file is executed by path and its public names copied across by hand,
because a star import's module path is resolved at compile time and cannot be a
variable. `_load_active_study` reproduces star-import semantics exactly, then
validates the result against a 15-name `REQUIRED_EXPORTS` frozenset — raising
`ValueError` naming both the study and the gap, rather than letting a missing export
surface as an `ImportError` from deep inside `calculations.sensitivity`. `load_study`
applies the same check, so reading a second study cannot enforce a weaker contract.

The choice is announced on stderr at import and recorded in `ACTIVE_STUDY_NAME`,
`ACTIVE_STUDY_PATH` and `ACTIVE_STUDY_SOURCE`. A notebook that must have its own
folder's study calls `assert_study_is_local` in its first cell.

After changing study, run `calculations.sensitivity.check_against_config()`: preset
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

### Available studies

| Folder | `STUDY_NAME` | System |
|---|---|---|
| `Guo_etal_2025/` | `Guo_etal_2025_316L` | 316L / Cr₂O₃ |
| `Fuerst_et_al_2024/` | `fuerst_etal_2024_model_config` | Hastelloy N / Cr₂O₃ |
| `incoloy802_cr2o3/` | `incoloy802_cr2o3` | Incoloy 802 (X40 NiCrAlTi 31/19) / Cr₂O₃ |
| `Application/` | `Guo_etal_2025_316L` | 316L / Cr₂O₃ — the delivery copy, and what the container gate tests |

The files are structurally identical and define the same interface; they differ only
in their values. `STUDY_NAME` is declared by the file rather than taken from the folder,
so renaming a folder cannot silently relabel a result — and `Application/` reports the
316L study it holds a copy of, which is why its reference set resolves.

Because a study is a file outside the package, it is not importable by a dotted name and
is not documented by autodoc here. For its structure and the meaning of every block, see
{doc}`../PARAMETERS`; to read one without making it active, use `load_study` above.
