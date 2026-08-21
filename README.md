# Hierarchical Hydrogen Permeation Model

Analytical model for hydrogen permeation through oxide-coated structural alloys, built as a hierarchy of increasing physical complexity (Levels 1–6). The active system is **Incoloy 802 (X40 NiCrAlTi 31/19)** with a **Cr₂O₃** surface oxide.

---

## Quick Start

All entry points are Jupyter notebooks in `Application/`. Open them in order:

| Notebook | What it does |
|---|---|
| `Application/Proposal.ipynb` | Main model walkthrough — Levels 1 through 5 (bulk transport, defective oxide, trapping) |
| `Application/Surface_proposal.ipynb` | Adds Level 6 surface kinetics (chemosorption, Langmuir coverage) |
| `Application/complete_model.ipynb` | Single end-to-end Level 5L6 run showing all resistance paths |
| `Application/model_analysis.ipynb` | Parametric sweeps and regime classification |
| `Application/sensitivity_regime_L5L6.ipynb` | Regime-stratified SA, Level 5L6 — regimes surface / oxide / metal |
| `Application/sensitivity_regime_L5.ipynb` | Regime-stratified SA, Level 5 (no surface kinetics) — regimes oxide / metal / defect |
| `Application/regime_parallel_coords.ipynb` | Parallel-coordinates views of the regime clusters and their sensitivities |

---

## Model Hierarchy

```
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

## Repository Structure

```
calculations/           Physics modules
  permeation_calc.py      Sieverts/diffusion core (Levels 1–2)
  oxide_permeation.py     Oxide transport (Level 2–3)
  parallel_oxide_defect_paths.py  Defect path assembly (Level 3)
  defective_metal.py      Grain boundary + trapping (Level 4)
  interface_solver.py     Oxide–metal interface coupling
  surface_kinetics.py     Chemosorption kinetics (Level 6)
  classify_regime.py      Regime classifier (surface- vs diffusion-limited)
  sensitivity.py          Regime-stratified given-data SA (PAWN + Borgonovo delta)
  config/
    model_config.py       Active material parameters (Incoloy 802 / Cr2O3_sample4)
    model_config1.py      Alternative config for comparison
    model_analysis_config.py  Parameter ranges for analysis workflows

data/                   Material property databases
  material_data.py        76 alloys from 8 literature sources (METALSssssssss dict)
  oxide_properties.py     Cr2O3 and Al2O3 transport properties
  literature_permeation_data.csv  CSV export of material_data.py
  sensitivity_parameters.csv      40 sensitivity parameters with ranges and references
  permeation_literature_database_v3 (3).xlsx  Raw literature data (Excel)

plot&docs/              Documentation, design notes, and archived plots

Application/            Notebooks and figures (main working directory)
```

---

## Configuration

Active material is set in `calculations/config/model_config.py`.  
To change the metal or oxide, edit the `build_simulation_config()` call at the bottom of that file:

```python
CONFIG = build_simulation_config(
    metal='metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985',
    oxide='Cr2O3_sample4',
)
```

Available metals and oxides are defined in the dicts above that line in the same file.

---

## Dependencies

- Python ≥ 3.9
- `numpy`, `scipy`, `matplotlib`
- `SALib` (for sensitivity analysis notebooks only)
- `jupyter` or VS Code with Jupyter extension

---

## Key References

| Source | Used for |
|---|---|
| Schmidt et al. 1985 | Incoloy 802 diffusivity and solubility |
| Nemanic et al. 2023 | Cr₂O₃ transport properties |
| Stover 1986 | Cr₂O₃ activation energies |
| Grant et al. 1988 | Surface dissociation rate constants |
| Lu et al. 2022 | Hastelloy N trap binding energies (TDS) |
| Zhu et al. 2021 | Hastelloy N grain/dislocation density (EBSD) |
| Young et al. 1997 | M6C carbide trap density |
