# Hierarchical Hydrogen Permeation Model

Analytical model for hydrogen permeation through oxide-coated structural alloys, built as a hierarchy of increasing physical complexity (Levels 1–6). Every level is steady-state and closed-form — there is no time axis; see [STEADY_STATE.md](STEADY_STATE.md) for why that is sufficient and what would break the argument.

Three material systems are configured. Exactly one is active at a time, currently **316L stainless steel with a Cr₂O₃ surface oxide** (`ACTIVE_STUDY = 'Guo_etal_2025_316L'`).

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
Dependency bounds are in [pyproject.toml](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/pyproject.toml); exact pins for the
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
outputs committed.

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
    model_config.py       The ACTIVE_STUDY switch — see Configuration below
    studies/              One module per material system

Application/            Notebooks (working copy)
Guo_etal_2025/          Notebooks as run for each study
Fuerst_et_al_2024/
incoloy802_cr2o3/

container/              Offline Docker deliverable + its acceptance gates
docs/                   Sphinx documentation source
latex/                  Model_Equations.tex — full equation catalogue with code locations
```

---

## Configuration

All material data and operating conditions resolve through a single switch in
[`calculations/config/model_config.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/config/model_config.py):

```python
ACTIVE_STUDY = 'Guo_etal_2025_316L'
```

| Study module | System |
|---|---|
| `Guo_etal_2025_316L` | 316L / Cr₂O₃ — **active** |
| `fuerst_etal_2024_model_config` | Hastelloy N / Cr₂O₃ |
| `incoloy802_cr2o3` | Incoloy 802 (X40 NiCrAlTi 31/19) / Cr₂O₃ |

Changing that one line switches the whole model; no other file changes. Restart
any running kernel afterwards, since the binding happens at import time.

**After switching, run:**

```python
from calculations.sensitivity import check_against_config
check_against_config()
```

The sweep ranges and regime presets were tuned for a specific material and are
not automatically valid for another one. Adding a study is documented in
[`calculations/config/studies/__init__.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/config/studies/__init__.py).

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

## Offline container

`container/` builds `hydrogen-model:1.1.0`, a self-contained JupyterLab image
for running the model where packages cannot be installed. See
[container/README.md](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/container/README.md) for end-user instructions and
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
