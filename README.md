# MHI Permeation Model

Analytical model for hydrogen permeation through oxide-coated structural alloys, built as a hierarchy of increasing physical complexity (Levels 1–6). Every level is steady-state and closed-form — there is no time axis; see [STEADY_STATE.md](STEADY_STATE.md) for why that is sufficient and what would break the argument.

Three material systems are configured, one per folder. Exactly one is loaded per process, decided by which folder you run in — see [Configuration](#configuration).

---

## Model Hierarchy

```text
Level 1  — Perfect metal (bulk diffusion + Sieverts solubility)
Level 2a — Perfect oxide layer on metal
Level 2b — Oxide + metal in series
Level 3  — Defective oxide (pinholes, cracks, grain boundaries)
Level 4  — Defective metal (grain boundary enhancement, trapping)
Level 5  — Full system: defective oxide + defective metal
Level 6  — Surface kinetics: dissociation / recombination (Langmuir coverage)
```

Level 5L6 = Level 5 + Level 6 combined.

---

## Quick Start

```bash
conda activate mace_env          # or any Python >= 3.9 environment
pip install -e .                 # makes `import calculations` work from anywhere
```

The editable install is what lets the notebooks and scripts resolve
`from calculations.config.model_config import ...` without a `sys.path` hack.
Dependency bounds are in [pyproject.toml](https://github.com/Azeezakinyemi999/MHI_permeation/blob/main/pyproject.toml); exact pins for the
reproducible container live in `container/requirements.lock.txt`.

Entry points are the notebooks in `Application/`:

| Notebook | What it does |
|---|---|
| `Proposal.ipynb` | Main walkthrough — Levels 1 through 5 (bulk transport, defective oxide, trapping) |
| `Surface_proposal.ipynb` | Adds Level 6 surface kinetics (dissociation, Langmuir coverage) |
| `sensitivity_regime_L5L6.ipynb` | Regime-stratified SA, Level 5L6 — regimes surface / oxide / metal |
| `sensitivity_regime_L5.ipynb` | Regime-stratified SA, Level 5 — regimes oxide / metal / defect |
| `regime_parallel_coords.ipynb` | Parallel-coordinates views of the L5L6 regime clusters |
| `regime_parallel_coords_L5.ipynb` | Same, for the Level 5 clusters |

`Application/` is the working copy. `Guo_etal_2025/`, `Fuerst_et_al_2024/` and
`incoloy802_cr2o3/` hold the same six notebooks as run against each study, with
outputs committed. Each of the four also holds the `study_config.py` its notebooks
run on.

The sensitivity notebooks write CSVs into `sa_results*/`, which is gitignored —
the parallel-coordinates notebooks read those files and cannot regenerate them,
so run the SA notebooks first.

---

## Repository Structure

```text
calculations/           Physics modules
  permeation_calc.py      Sieverts/diffusion core (Levels 1–2)
  oxide_permeation.py     Oxide molecular transport (Levels 2–3)
  parallel_oxide_defect_paths.py  Parallel-path assembly (Level 3)
  defective_metal.py      Grain boundary enhancement + Oriani trapping (Level 4)
  interface_solver.py     Oxide–metal interface pressure solve
  surface_kinetics.py     Dissociation kinetics (Level 6)
  classify_regime.py      Hierarchical regime labels (Levels 1–4)
  sensitivity.py          Regime-stratified given-data SA (PAWN + Borgonovo delta)
  config/
    model_config.py       Loads the one active study — see Configuration below
    __init__.py           Where studies live, and the add-a-study recipe

Application/            Notebooks + study_config.py (working copy, and what ships)
Guo_etal_2025/          One folder per study: its notebooks AND its input file
Fuerst_et_al_2024/
incoloy802_cr2o3/

container/              Offline Docker deliverable + its acceptance gates
docs/                   Sphinx documentation source
latex/                  Model_Equations.tex — full equation catalogue with code locations
```

---

## Configuration

All material data and operating conditions resolve through
[`calculations/config/model_config.py`](https://github.com/Azeezakinyemi999/MHI_permeation/blob/main/calculations/config/model_config.py),
which loads exactly one study. **A study is a folder**, and its
`study_config.py` is the input file — material properties, operating conditions,
sweep ranges and presets all live there. Running a notebook in that folder uses
that folder's study. There is no switch to edit.

| Study folder | Input file | System |
|---|---|---|
| `Guo_etal_2025/` | `study_config.py` | 316L / Cr₂O₃ |
| `Fuerst_et_al_2024/` | `study_config.py` | Hastelloy N / Cr₂O₃ |
| `incoloy802_cr2o3/` | `study_config.py` | Incoloy 802 (X40 NiCrAlTi 31/19) / Cr₂O₃ |
| `Application/` | `study_config.py` | 316L / Cr₂O₃ — the delivery copy, and what the container gate tests |

Resolution order is `$MHI_STUDY` (a path, relative to the repo root), then a
`study_config.py` in the working directory or a parent. **There is no default** —
anything run from the repository root must name a study or it raises, listing the
folders it could have meant. The choice is printed on stderr at import, and the
notebooks call `assert_study_is_local()` in their first cell so a kernel started in
the wrong directory fails instead of quietly computing another study.

```bash
MHI_STUDY=incoloy802_cr2o3/study_config.py python your_script.py
```

In VS Code this needs `"jupyter.notebookFileRoot": "${fileDirname}"`, which
`.vscode/settings.json` sets — otherwise kernels start at the workspace root.
Restart any running kernel after changing study, since the binding happens at
import time.

**After switching, run:**

```python
from calculations.sensitivity import check_against_config
check_against_config()
```

The sweep ranges and regime presets were tuned for a specific material and are
not automatically valid for another one. Adding a study — make a folder, put the
notebooks and a `study_config.py` in it — is documented in
[`calculations/config/__init__.py`](https://github.com/Azeezakinyemi999/MHI_permeation/blob/main/calculations/config/__init__.py).

---

## Documentation

| Document | Contents |
|---|---|
| [PARAMETERS.md](PARAMETERS.md) | Every model input, by physical category, with units and provenance |
| [OUTPUTS.md](OUTPUTS.md) | Every returned quantity, function by function |
| [STEADY_STATE.md](STEADY_STATE.md) | Why the model needs no time axis, and what would invalidate that |
| [TRAPPING_VALIDATION.md](TRAPPING_VALIDATION.md) | Zero-free-parameter TDS validation of the trapping block, plus equation appendices |
| [PACKAGING_1.0.0.md](PACKAGING_1.0.0.md) | Scope and acceptance record for the container release |
| `latex/Model_Equations.tex` | Every equation transcribed, with its `file → function` location |

A browsable HTML version, combining the above with an API reference generated
from the docstrings:

```bash
pip install -e ".[docs]"
cd docs && make html      # then open _build/html/index.html
```

The build runs with warnings-as-errors, so a malformed docstring or a
repo-relative link to a source file will fail it.

---

## Tests

```bash
pip install -e ".[test]"
pytest
```

`tests/test_study_resolution.py` covers how a study is selected: that each folder
resolves its own `study_config.py`, that the repository root raises instead of
guessing, that `$MHI_STUDY` overrides and resolves relative paths from the repo root,
that `assert_study_is_local()` rejects an inherited or env-selected study, and that
every study satisfies `REQUIRED_EXPORTS` and has a reference set in
`container/model_smoke_test.py`. Study folders are discovered, not listed, so a new
study is covered automatically.

---

## Offline container

`container/` builds `hydrogen-model:1.1.0`, a self-contained JupyterLab image
for running the model where packages cannot be installed. See
[container/README.md](https://github.com/Azeezakinyemi999/MHI_permeation/blob/main/container/README.md) for end-user instructions and
[PACKAGING_1.0.0.md](PACKAGING_1.0.0.md) for what is pinned and why. The image
bind-mounts the code rather than installing it, so `pyproject.toml` does not
affect the release.

`container/model_smoke_test.py` is the numerical regression gate: it asserts
pinned reference values through the real solvers and is study-aware.

---

## Key References

| Source | Used for |
|---|---|
| Guo et al. 2025 | 316L diffusivity, solubility and permeability (active study) |
| Fuerst et al. 2024 | Hastelloy N permeation data |
| Schmidt et al. 1985 | Incoloy 802 diffusivity and solubility |
| Nemanic et al. 2023 | Cr₂O₃ transport properties |
| Stover 1986 | Cr₂O₃ activation energies |
| Grant et al. 1988 | Surface dissociation rate constants |
| Lu et al. 2022 | Trap binding energies (TDS) |
| Zhu et al. 2021 | Grain size and dislocation density (EBSD) |
| Young et al. 1997 | M₆C carbide trap density |
| Oriani 1970 | Local-equilibrium trapping model |
| Strehlow & Savage 1974 | Parallel-path transport through defective coatings |

Full bibliography in `latex/references.bib` (32 entries).
