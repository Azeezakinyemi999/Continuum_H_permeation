# How to switch material study

All material data, operating conditions, microstructure and defect parameters
resolve through one loader, so exactly one study is active at a time. Which one
is decided by **where you are running from**, not by a switch you edit.

## The normal case: run the notebook in its study folder

Each study folder holds a `study_config.py`, and that file **is** the study — every
material property, operating condition, sweep range and preset it runs on. Running a
notebook with that folder as its working directory selects it; nothing to edit, and the
study follows the same rule the outputs already follow (`RESULTS_DIR = "sa_results"` is
resolved the same way).

| Notebook folder | Study |
|---|---|
| `Guo_etal_2025/` | 316L stainless / Cr₂O₃ |
| `Fuerst_et_al_2024/` | Hastelloy N / Cr₂O₃ |
| `incoloy802_cr2o3/` | Incoloy 802 (X40 NiCrAlTi 31/19) / Cr₂O₃ |
| `Application/` | 316L stainless / Cr₂O₃ (the delivery copy) |

So switching study means opening a different folder's notebook, and changing a
material means editing that folder's own file. Every module in `calculations/` still
imports from `model_config`; only its resolution changed.

```{important}
**In VS Code, set the kernel's directory.** VS Code starts Jupyter kernels in the
*workspace root*, not the notebook's folder, where no study applies and the import
fails. `.vscode/settings.json` in this repo sets

    "jupyter.notebookFileRoot": "${fileDirname}"

Jupyter Lab already uses the notebook's own directory. **Restart the kernel** after
changing anything about which study applies — the binding happens at import time, so
`importlib.reload` on a notebook module is not enough.
```

## Resolution order

```{list-table}
:header-rows: 1
:widths: 26 44 30

* - Order
  - Source
  - Used by
* - 1
  - `$MHI_STUDY` — a path to a study file, relative to the repository root
  - the container gate, the docs tools, CI, any deliberate override
* - 2
  - `study_config.py` in the working directory, then each parent up to the
    repository root
  - **notebooks**
* - 3
  - nothing — with no study selected the import raises, listing the folders it
    could have meant
  - —
```

The search never looks above the repository root, so a stray `study_config.py`
elsewhere on the machine can never be picked up.

There is deliberately no default. A fallback would be taken by exactly the runs that
never considered which material they wanted, and taken silently. Tools that legitimately
run from the root name their study instead: `docs/conf.py` and `docs/_tools/*` use
`Guo_etal_2025/study_config.py`, and `container/model_smoke_test.py` uses
`Application/study_config.py`, the folder the release ships.

## Know which study you got

The loader announces itself on stderr at import:

```text
study: fuerst_etal_2024_model_config  (Fuerst_et_al_2024/study_config.py, via folder)
```

and exposes `ACTIVE_STUDY_NAME`, `ACTIVE_STUDY_PATH` and `ACTIVE_STUDY_SOURCE`
(`'env'` or `'folder'`). A notebook that must have its *own* folder's
study says so in its first cell:

```python
from calculations.config import assert_study_is_local
print(assert_study_is_local())
```

It raises unless the active study came from a `study_config.py` in the working
directory — turning a kernel started in the wrong place into a stack trace instead
of a plausible-looking number in the wrong folder's `sa_results/`.

The sensitivity runs also write a `study.txt` beside their CSVs, and warn if a
results directory that already held one study's output starts receiving another's.

## Reading one study from anywhere

```bash
MHI_STUDY=Fuerst_et_al_2024/study_config.py python your_script.py
```

The env var beats the folder file, and a relative path is resolved from the repository
root — not from the working directory, so the same command works from anywhere.

## Then validate

Preset yields and sweep ranges were tuned for a specific material and are not
automatically valid for another one. Run:

```python
from calculations.sensitivity import check_against_config
check_against_config()
```

It verifies that every swept parameter has a default, that defaults are finite
and scalar, that each default lies inside its own sweep range, and that the L5
and L5L6 parameter sets agree where they overlap.

```{warning}
On the 316L study (`Guo_etal_2025/`, `Application/`) this currently reports **four
pre-existing problems** — `D_ref`, `k_diss_metal_ref` and `E_diss_metal` defaults lying
outside their own sweep ranges. They are not caused by switching; they are present in
that study's file as shipped, and the Incoloy and Hastelloy N studies report clean.

A default outside its sweep range means the sensitivity analysis never evaluates the
nominal case, so either the default or the range is wrong. Resolve before trusting a
sensitivity result on those parameters.
```

## Do not hardcode material keys

The `METALS` keys are deliberately different in each study:

```text
316L_Guo_2025
Hastelloy_N_fuerst_2024
metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985
```

So a notebook containing `METALS['316L_Guo_2025']` breaks the moment it is run in
another study's folder. Resolve the active names instead:

```python
from calculations.config.model_config import build_simulation_config
SIM = build_simulation_config()
SIM['metal_name'], SIM['oxide_name']
```

`OXIDES` is uniform — all three studies hold only `Cr2O3_sample4` — so oxide
lookups are safe either way, but resolving through `build_simulation_config` is
still the supported route.

## Reading a non-active study without switching

The property getters take an explicit registry, so you can read another study's
data without changing directory or setting anything:

```python
from calculations.config.model_config import load_study
from calculations.oxide_permeation import get_metal_properties_at_T

other = load_study('Fuerst_et_al_2024/study_config.py')
props = get_metal_properties_at_T('Hastelloy_N_fuerst_2024', 900.0,
                                  metals=other.METALS)
```

`load_study` validates the file against `REQUIRED_EXPORTS` exactly as the active study
is validated, and resolves a relative path from the repository root rather than the
working directory.

This is the supported way to compare studies side by side. Note that it reads
material properties only — the solvers still take their defaults from whichever
study was resolved at import.

## When to switch

The 316L study is metal-limited at its default operating point *under Level 5*: the
1 µm Cr₂O₃ layer contributes `frac_oxide = 7.8e-3` of the resistance and `PRF ≈
1.004`. If you are investigating oxide or defect physics with Level 5, that operating
point will show you little. Level 5L6 is oxide-limited on the same wall.

To make the oxide matter, either increase `L_oxide`, lower the temperature toward
the validated [473, 773] K window where the oxide is less permeable, or work in a
folder whose study has an oxide that is a genuine barrier.
