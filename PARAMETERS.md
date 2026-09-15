# Complete Parameter Reference

Every input required to run the model, extracted from [`calculations/`](calculations/) and the active study config. Generated against `ACTIVE_STUDY = 'Guo_etal_2025_316L'`; values change per study, names and units do not.

**The authoritative flat set is `DEFAULT_PARAMS_LEVEL5L6` — 46 parameters.** `DEFAULT_PARAMS_LEVEL5` is the same set minus the 10 surface-kinetics entries (36 parameters).

Counts: **46 total** = 36 sampled in the L5L6 SA + 10 never sampled.

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
| `oxide_thickness` | 4.8e-8 | m | oxide layer — from `OXIDES[...]['thickness']` |

> ⚠️ `CONDITIONS['L_oxide']` is **1e-6 m**, not 4.8e-8 m. The flat parameter dict takes oxide thickness from the *oxide* dict, so `L_oxide` in `CONDITIONS` is only used by thickness sweeps. Two different oxide thicknesses coexist in the config — know which one your call path reads.

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

## F. Metal grain structure (Level 4) — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `grain_size` | 1e-4 | m | mean grain diameter (100 µm) |
| `gb_thickness` | 5e-10 | m | grain-boundary width |
| `grain_shape` | `'equiaxed'` | – | stereology model — **fixed** |
| `gb_type` | `'LAGB'` | – | selects GB enhancement data — **fixed** |
| `lattice_density` | 8.774e28 | m⁻³ | interstitial site density $N_L$ |

## G. Traps (Level 4) — 8

Four populations × (binding energy, density).

| Parameter | Value | Units | Source |
|---|---|---|---|
| `trap_vacancy_E_b` | 41489 | J/mol | Ni estimate (0.43 eV) — **no alloy data** |
| `trap_vacancy_N_T` | 1e26 | m⁻³ | model default — **no measurement** |
| `trap_dislocation_E_b` | 19297 | J/mol | Lu 2022 TDS Peak 1 (0.20 eV) |
| `trap_dislocation_N_T` | 8.16e12 | m⁻³ | Zhu 2021 GND — **see units warning** |
| `trap_gb_E_b` | 26051 | J/mol | Lu 2022 TDS Peak 2 (0.27 eV) |
| `trap_gb_N_T` | 6e14 | m⁻³ | geometric estimate — **see units warning** |
| `trap_carbide_E_b` | 26051 | J/mol | Lu 2022 M6C TDS (0.27 eV) |
| `trap_carbide_N_T` | 2e25 | m⁻³ | Young 1997 upper bound |

> ⚠️ **Two densities look dimensionally inconsistent.** `grain_boundary_density(grain_size=1e-4)` in [`defective_metal.py:280`](calculations/defective_metal.py#L280) returns **3.0e23 m⁻³**, while the config hardcodes 6e14 — a factor of 5e8. And EBSD GND density is reported in m⁻² (line length per volume); converting to trap sites needs $\rho/b \approx 3.3\times10^{22}\,\mathrm{m^{-3}}$, not the raw 8.16e12 used as m⁻³. Consequence: `vacancies` carries 97–99 % of the trapping term at every temperature, so the two TDS-derived traps currently contribute nothing.

## H. Model options — 3

Switches, not physical quantities. All **fixed** (never sampled).

| Parameter | Value | Meaning |
|---|---|---|
| `include_gb_enhancement` | `True` | enable GB short-circuit diffusion path |
| `include_trapping` | `True` | enable Oriani trapping reduction |
| `D_eff_method` | `'average'` | how $D_{eff}(x)$ is collapsed to a scalar |

## I. Oxide surface kinetics (Level 6) — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `k_diss_ref` | 9.487e-8 | mol m⁻² s⁻¹ Pa⁻¹ | dissociation rate constant at `T_ref_surface` |
| `E_diss` | 57950 | J/mol | dissociation activation energy |
| `K_eq_ref` | 1e-4 | Pa⁻¹ | Langmuir adsorption equilibrium constant |
| `H_eq` | 20000 | J/mol | adsorption enthalpy |
| `T_ref_surface` | 1623 | K | reference temperature — **fixed** |

## J. Metal surface kinetics (Level 6) — 5

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| `k_diss_metal_ref` | 2.6287e-11 | mol m⁻² s⁻¹ Pa⁻¹ | metal dissociation rate constant |
| `E_diss_metal` | 54996 | J/mol | metal dissociation activation energy |
| `K_eq_metal_ref` | 1.52e-3 | Pa⁻¹ | metal Langmuir constant |
| `H_eq_metal` | 19297 | J/mol | metal adsorption enthalpy |
| `T_ref_surface_metal` | 673 | K | reference temperature — **fixed** |

---

## K. Never sampled in the SA — 10

`T_ref_metal`, `T_ref_oxide`, `P_downstream`, `grain_shape`, `gb_type`, `include_gb_enhancement`, `include_trapping`, `D_eff_method`, `T_ref_surface`, `T_ref_surface_metal`

Reference temperatures are measurement anchors (varying them is meaningless), `P_downstream` is 0 by construction, and the rest are switches or categorical choices.

`temperature` **is** in `SUGGESTED_RANGES_*` but is removed before sampling by `presets_without(REGIME_PRESETS, 'temperature')` and pinned via `fixed_params` — see [`sensitivity.py:75-90`](calculations/sensitivity.py#L75). So the SA samples **35** parameters, not 36.

---

## L. Additional inputs not in the flat dict

These live in structured config dicts or as function defaults, and never appear in `DEFAULT_PARAMS_*`.

### Structured-dict only

| Parameter | Value | Location | Meaning |
|---|---|---|---|
| `use_sieverts_pinhole` | `False` | `OXIDE_DEFECTS` | pinhole boundary-condition mode |
| `include_gb_trapping` | `False` | `MICROSTRUCTURE` | GB **trapping** — distinct from `include_gb_enhancement` |
| `area_fraction` | 0.02 | `OXIDE_DEFECTS` | total defect fraction (= sum of components) |
| `type` | `'mixed'` | `OXIDE_DEFECTS` | defect-population label |
| `N_L` | 8.774e28 | `MICROSTRUCTURE` | alias of `lattice_density` |

### Function defaults

| Parameter | Default | Function |
|---|---|---|
| `sites_per_area` | 1e19 m⁻² | [`grain_boundary_density`](calculations/defective_metal.py#L280) |
| `n_points` | 100 | [`calculate_concentration_profile`](calculations/interface_solver.py#L303) |
| `threshold_traps` | 0.5 | [`classify_regime_level4_metal`](calculations/classify_regime.py#L45) |
| `rule`, `threshold` | `'argmax'`, 0.5 | [`assign_regime`](calculations/sensitivity.py#L511) |
| `temperature_unit` | `'K'` | [`gb_enhancement_factor`](calculations/defective_metal.py#L470) |
| `data_source` | `'default'` | [`gb_enhancement_factor`](calculations/defective_metal.py#L470) |
| `material` | `'Incoloy800'` | [`vacancy_concentration`](calculations/defective_metal.py#L680) |
| `condition` | `'equilibrium'` | [`vacancy_concentration`](calculations/defective_metal.py#L680) |
| `quench_temperature` | `None` | [`vacancy_concentration`](calculations/defective_metal.py#L680) |
| `method` | `'brentq'` | [`solve_interface_pressure`](calculations/interface_solver.py#L89) |

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

## M. Outputs — not inputs

Do not supply these; the model computes them.

$D_{eff}$, `reduction_factor`, `mobile_fraction`, `trapping_term`, $\theta_i$ (trap occupancy), `dominant_trap`, $\Phi$, `PRF`, `flux`, `permeability`, `C_up`/`C_down`/`C_array`, `P_interface`, `regime` / `regime_hierarchy`, `frac_surface`/`frac_oxide`/`frac_metal`/`frac_defect`.

Also derived, not independent: `D_0`, `K_s0`, `Phi_0`, `Phi_ref`, `Q_p` in `METALS`, and `D_ox_0`, `K_ox_0`, `Phi_ox_ref`, `Q_p_ox_J_per_mol` in `OXIDES`. These are alternate Arrhenius parameterisations of the same transport data — the model reads the `*_ref` + `E_*`/`H_*` form.

---

## N. Minimum set per level

| Level | Needs | Count |
|---|---|---|
| L1 — perfect metal | A + C(metal) + D | 9 |
| L2a — perfect oxide | B + C(oxide) + D | 9 |
| L2b — oxide + metal series | A + B + C + D | 15 |
| L3 — defective oxide | + E | 20 |
| L4 — defective metal | A + C(metal) + D + F + G + H | 25 |
| L5 — full system | A–H | 36 |
| L5L6 — + surface kinetics | A–J | 46 |
