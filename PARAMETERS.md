# Complete Parameter Reference

Every input required to run the model, extracted from [`calculations/`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/tree/main/calculations/) and the active study config. Generated against `ACTIVE_STUDY = 'Guo_etal_2025_316L'`; values change per study, names and units do not.

**The authoritative flat set is `DEFAULT_PARAMS_LEVEL5L6` — 35 parameters.** `DEFAULT_PARAMS_LEVEL5` is the same set minus the 10 surface-kinetics entries (25 parameters).

Counts: **35 total** = 27 with a range in the L5L6 SA + 8 never sampled.

There are no trap or lattice-site-density parameters: hydrogen trapping does not change the steady-state flux and was removed from the model (see [STEADY_STATE.md](STEADY_STATE.md)).

---

## A. Metal bulk transport — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `D_ref` | 2.8e-11 | m²/s | lattice diffusivity at `T_ref_metal` |
| `E_D` | 52102 | J/mol | activation energy for diffusion |
| `K_s_ref` | 0.039 | mol m⁻³ Pa⁻⁰·⁵ | Sieverts constant at `T_ref_metal` |
| `H_s` | 9648 | J/mol | heat of solution |
| `T_ref_metal` | 673 | K | measurement reference temperature — **fixed** |

## B. Oxide bulk transport — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `D_ox_ref` | 7.8e-19 | m²/s | oxide diffusivity at `T_ref_oxide` |
| `E_D_ox` | 70434 | J/mol | oxide diffusion activation energy |
| `K_ox_ref` | 0.35417 | mol m⁻³ Pa⁻⁰·⁵ | oxide solubility constant |
| `H_sol_ox` | 163566 | J/mol | oxide heat of solution |
| `T_ref_oxide` | 673 | K | measurement reference — **fixed** |

## C. Geometry — 2

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `metal_thickness` | 1e-3 | m | wall thickness — from `CONDITIONS['L_metal']` |
| `oxide_thickness` | 1e-6 | m | oxide layer — from `CONDITIONS['L_oxide']` |

> Both thicknesses now come from `CONDITIONS`, which is the single description of the wall being simulated. The oxide entry keeps a `sample_thickness` (4.8e-8 m) recording the specimen its transport data were measured on; nothing reads it. Until recently that value silently won, and every model call ran a geometry 20× thinner than this table claimed.

## D. Operating conditions — 3

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `P_upstream` | 1e5 | Pa | upstream H₂ partial pressure |
| `P_downstream` | 0 | Pa | downstream pressure — **fixed** |
| `temperature` | 873 | K | operating temperature (pinned in SA) |

## E. Oxide defects (Level 3) — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `f_pinhole` | 0.01 | – | pinhole area fraction |
| `f_crack` | 0.005 | – | crack area fraction |
| `f_gb_defect` | 0.005 | – | oxide grain-boundary area fraction |
| `crack_thickness_factor` | 0.1 | – | $L_{crack} = f \times L_{oxide}$ |
| `gb_diffusivity_factor` | 10 | – | $D_{gb} = f \times D_{oxide}$ |

## F. Metal grain structure (Level 4) — 4

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `grain_size` | 1e-4 | m | mean grain diameter (100 µm) |
| `gb_thickness` | 5e-10 | m | grain-boundary width |
| `grain_shape` | `'equiaxed'` | – | stereology model — **fixed** |
| `gb_type` | `'LAGB'` | – | selects GB enhancement data — **fixed** |

## G. Model options — 1

A switch, not a physical quantity. **Fixed** (never sampled).

| Parameter | Value | Meaning |
|---|---|---|
| `enable_gb_enhancement` | `True` | enable GB short-circuit diffusion path |

## H. Oxide surface kinetics (Level 6) — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `k_diss_ref` | 9.487e-8 | mol m⁻² s⁻¹ Pa⁻¹ | dissociation rate constant at `T_ref_surface` |
| `E_diss` | 57950 | J/mol | dissociation activation energy |
| `K_eq_ref` | 1e-4 | Pa⁻¹ | Langmuir adsorption equilibrium constant |
| `H_eq` | 20000 | J/mol | adsorption enthalpy |
| `T_ref_surface` | 1623 | K | reference temperature — **fixed** |

## I. Metal surface kinetics (Level 6) — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `k_diss_metal_ref` | 2.6287e-11 | mol m⁻² s⁻¹ Pa⁻¹ | metal dissociation rate constant |
| `E_diss_metal` | 54996 | J/mol | metal dissociation activation energy |
| `K_eq_metal_ref` | 1.52e-3 | Pa⁻¹ | metal Langmuir constant |
| `H_eq_metal` | 19297 | J/mol | metal adsorption enthalpy |
| `T_ref_surface_metal` | 673 | K | reference temperature — **fixed** |

---

## J. Never sampled in the SA — 8

`T_ref_metal`, `T_ref_oxide`, `P_downstream`, `grain_shape`, `gb_type`, `enable_gb_enhancement`, `T_ref_surface`, `T_ref_surface_metal`

Reference temperatures are measurement anchors (varying them is meaningless), `P_downstream` is 0 by construction, and the rest are switches or categorical choices.

`temperature` **is** in `SUGGESTED_RANGES_*` but is removed before sampling by `presets_without(REGIME_PRESETS, 'temperature')` and pinned via `fixed_params` — see [`sensitivity.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/sensitivity.py). So the L5L6 SA samples **26** parameters, not 27 (and the L5 SA 18, not 19).

---

## K. Additional inputs not in the flat dict

These live in structured config dicts or as function defaults, and never appear in `DEFAULT_PARAMS_*`.

### Structured-dict only

| Parameter | Value | Location | Meaning |
|---|---|---|---|
| `use_sieverts_pinhole` | `False` | `OXIDE_DEFECTS` | pinhole boundary-condition mode |
| `area_fraction` | 0.02 | `OXIDE_DEFECTS` | total defect fraction (= sum of components) |
| `type` | `'mixed'` | `OXIDE_DEFECTS` | defect-population label |

### Function defaults

| Parameter | Default | Function |
|---|---|---|
| `n_profile_points` | 100 | [`calculate_concentration_profile`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/interface_solver.py) |
| `rule`, `threshold` | `'argmax'`, 0.5 | [`assign_regime`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/sensitivity.py) |
| `temperature_unit` | `'K'` | [`gb_enhancement_factor`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/defective_metal.py) |
| `data_source` | `'default'` | [`gb_enhancement_factor`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/defective_metal.py) |
| `method` | `'brentq'` | [`solve_interface_pressure`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/interface_solver.py) |

### Sweep and SA controls

| Parameter | Value | Location |
|---|---|---|
| `T_range`, `n_T_points` | (623, 1200), 20 | `CONDITIONS` |
| `P_range`, `n_P_points` | (1e-4, 1e12), 40 | `CONDITIONS` |
| `L_metal_range`, `n_L_points` | (1e-4, 5e-3), 20 | `CONDITIONS` |
| `L_oxide_range` | (1e-7, 1e-5) | `CONDITIONS` |
| `DEFAULT_SWEEP_TEMPERATURES` | (773, 1073, 1273) K | study config |
| `SEED` | 42 | notebook cell 3 |
| `TARGET_CLUSTER_SIZE` | — | `sensitivity.py` |

---

## L. Outputs — not inputs

Do not supply these; the model computes them.

$D_{eff}$ (the grain-boundary-enhanced metal diffusivity), `enhancement_ratio`, $\theta$ (surface coverage), $\Phi$, `PRF`, `flux`, `permeability`, `C_up`/`C_down`/`C_array`, `P_interface`, `regime` / `regime_hierarchy`, `frac_surface`/`frac_oxide`/`frac_metal`/`frac_defect`.

Also derived, not independent: `D_0`, `K_s0`, `Phi_0`, `Phi_ref`, `Q_p` in `METALS`, and `D_ox_0`, `K_ox_0`, `Phi_ox_ref`, `Q_p_ox_J_per_mol` in `OXIDES`. These are alternate Arrhenius parameterisations of the same transport data — the model reads the `*_ref` + `E_*`/`H_*` form.

---

## M. Minimum set per level

| Level | Needs | Count |
|---|---|---|
| L1 — perfect metal | A + C(metal) + D | 9 |
| L2a — perfect oxide | B + C(oxide) + D | 9 |
| L2b — oxide + metal series | A + B + C + D | 15 |
| L3 — defective oxide | + E | 20 |
| L4 — defective metal | A + C(metal) + D + F + G | 14 |
| L5 — full system | A–G | 25 |
| L5L6 — + surface kinetics | A–I | 35 |
