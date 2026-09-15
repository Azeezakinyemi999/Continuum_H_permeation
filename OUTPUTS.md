# Complete Output Reference

Every quantity the model returns, function by function, extracted by AST-parsing [`calculations/`](calculations/) and verified by execution where possible. Companion to [PARAMETERS.md](PARAMETERS.md).

Generated against `ACTIVE_STUDY = 'Guo_etal_2025_316L'`.

---

## 1. The headline outputs — what a full run gives you

These are the two top-level entry points the sensitivity analysis calls. **Verified by execution**, not just parsed.

### `level5_model_wrapper(params_dict)` → 20 keys

[`sensitivity.py:188`](calculations/sensitivity.py#L188)

| Key | Example | Units | Meaning |
|---|---|---|---|
| `flux` | 3.197e-06 | mol m⁻² s⁻¹ | **total permeation flux** |
| `permeability` | 4.0e-15 | mol m⁻¹ s⁻¹ Pa⁻⁰·⁵ | effective system permeability |
| `PRF` | 1.0 | – | permeation reduction factor |
| `D_eff` | 1.747e-10 | m²/s | effective metal diffusivity after GB + trapping |
| `D_modification` | 0.7389 | – | $D_{eff}/D_{lattice}$ |
| `modification_factor` | 0.7389 | – | duplicate of `D_modification` |
| `P_interface` | 9.996e+04 | Pa | pressure at the oxide/metal interface |
| `flux_intact` | 3.133e-06 | mol m⁻² s⁻¹ | flux through intact oxide |
| `flux_defect` | 6.395e-08 | mol m⁻² s⁻¹ | flux through defect paths |
| `flux_bare_metal` | 3.197e-06 | mol m⁻² s⁻¹ | reference flux, no oxide |
| `defect_enhancement` | 1.0 | – | defect path enhancement |
| `frac_oxide` | 3.758e-04 | – | fraction of total resistance in oxide |
| `frac_metal` | 0.9796 | – | fraction in metal |
| `frac_defect` | 0.02 | – | fraction in defect paths |
| `regime` | `'metal'` | – | **argmax regime label** |
| `D_metal` | 2.364e-10 | m²/s | lattice diffusivity at T |
| `K_s_metal` | 0.05789 | mol m⁻³ Pa⁻⁰·⁵ | Sieverts constant at T |
| `D_ox` | 1.395e-17 | m²/s | oxide diffusivity at T |
| `K_ox` | 286.9 | mol m⁻³ Pa⁻⁰·⁵ | oxide solubility at T |
| `temperature` | 873 | K | echo of input |

With `return_full_record=True`, adds 2 more: `regime_hierarchy`, `dominant_path`.

### `level5L6_model_wrapper(params_dict)` → 22 keys

[`sensitivity.py:548`](calculations/sensitivity.py#L548)

Same as above except: **no** `flux_bare_metal`, `frac_defect`, `defect_enhancement`, `modification_factor`; **adds** these 6:

| Key | Example | Units | Meaning |
|---|---|---|---|
| `frac_surface` | 2.924e-04 | – | fraction of resistance in surface step |
| `theta` | 0.6163 | – | **Langmuir surface coverage** |
| `k_diss` | 2.37e-09 | mol m⁻² s⁻¹ Pa⁻¹ | oxide dissociation rate at T |
| `K_eq` | 2.799e-05 | Pa⁻¹ | oxide adsorption equilibrium at T |
| `k_diss_metal` | 2.498e-10 | mol m⁻² s⁻¹ Pa⁻¹ | metal dissociation rate at T |
| `K_eq_metal` | 3.349e-03 | Pa⁻¹ | metal adsorption equilibrium at T |

With `return_full_record=True`, adds 2 more: `system_rate_limiting`, `dominant_path`.

> ⚠️ Two things the verification run surfaced. `PRF` returns **NaN** from the L5L6 wrapper at default parameters (it is finite from L5). And `T_operating = 873 K` sits just below the GB-enhancement data range `[873.1, 1273.2] K`, so every default L5/L5L6 call emits an extrapolation warning from [`defective_metal.py:619`](calculations/defective_metal.py#L619).

---

## 2. `permeation_calc.py` — bulk transport (Levels 1–4)

| Function | Returns |
|---|---|
| `sieverts_concentration()` | scalar: `concentration` [mol/m³] |
| `fick_flux()` | scalar: `flux` [mol m⁻² s⁻¹] |
| `calculate_simple_metal_flux()` | **7 keys:** `flux`, `C_up`, `C_down`, `permeability`, `Diffusivity`, `solubility`, `units` |
| `calculate_defective_metal_flux()` | **14 keys:** `flux`, `C_up`, `C_down`, `permeability`, `D_eff`, `D_lattice`, `modification_factor`, `regime_classification`, `regime`, `regime_base`, `regime_detail`, `microstructure_details`, `profiles`, `units` |
| `calculate_defective_metal_flux_sieverts()` | scalar flux (expression) |

`profiles` is a nested dict holding the spatial arrays — `x_array`, `C_array`, `D_array`, `theta_array`, `gb_factor_array`.

## 3. `defective_metal.py` — microstructure (Level 4)

| Function | Returns |
|---|---|
| `trap_occupancy()` | **6 keys:** `theta`, `K_equilibrium`, `approximation_used`, `theta_lattice`, `saturation_warning`, `trap_concentration` |
| `grain_boundary_density()` | **7 keys:** `trap_density`, `volume_fraction`, `surface_per_volume`, `mean_intercept`, `shape_factor`, `dimensionality`, `warnings` |
| `gb_enhancement_factor()` | **7 keys:** `enhancement_factor`, `uncertainty`, `interpolation_method`, `temperature_K`, `temperature_C`, `gb_type_factor`, `data_range_K` |
| `vacancy_concentration()` | **10 keys:** `concentration`, `site_fraction`, `formation_energy`, `lattice_sites`, `is_equilibrium`, `effective_temperature`, `minimum_concentration`, `material`, `melting_point`, `warnings` |
| `calculate_effective_diffusivity_trapping()` | **10 keys:** `D_eff`, `D_lattice`, `trapping_term`, `theta_total`, `trap_contributions`, `reduction_factor`, `dominant_trap`, `mobile_fraction`, `saturation_warnings`, `temperature` |
| `calculate_gb_enhanced_diffusivity()` | **14 keys:** `D_eff`, `D_bulk`, `D_gb`, `f_gb`, `f_bulk`, `enhancement_ratio`, `gb_enhancement_factor`, `regime`, `model_used`, `percolation_warning`, `grain_size`, `gb_thickness`, `gb_type`, `temperature` |
| `combined_microstructure_model()` | **16 keys:** `D_eff`, `D_lattice`, `D_gb_enhanced`, `gb_enhancement`, `trapping`, `overall_factor`, `dominant_effect`, `regime`, `mode`, `calculation_sequence`, `parameters`, `warnings`, `theta_total`, `trap_details`, `gb_enhancement_factor`, `trapping_reduction_factor` |

`trap_contributions` is a **list of per-trap dicts**, each with 7 keys: `name`, `theta`, `binding_energy`, `density`, `trapped_concentration`, `K_equilibrium`, `trapping_contribution`.

## 4. `oxide_permeation.py` — oxide transport (Levels 2–3)

| Function | Returns |
|---|---|
| `molecular_diffusion_flux()` | scalar: `flux` |
| `calculate_oxide_resistance()` | scalar: `resistance` |
| `calculate_metal_resistance()` | scalar: `resistance` |
| `get_oxide_properties_at_T()` | **3 keys:** `D_ox`, `K_ox`, `thickness` |
| `get_metal_properties_at_T()` | **2 keys:** `D_metal`, `K_s_metal` |
| `compare_resistances()` | **5 keys:** `R_oxide`, `R_metal`, `ratio`, `limiting_mechanism`, `P_interface` |
| `calculate_transition_pressure()` | scalar: `P_transition` |
| `pressure_dependence_analysis()` | `results` — array/dict over the pressure sweep |

## 5. `interface_solver.py` — oxide/metal coupling

| Function | Returns |
|---|---|
| `calculate_metal_flux_sieverts()` | scalar: `flux` |
| `solve_interface_pressure()` | **7 keys:** `P_interface`, `P_upstream`, `P_downstream`, `flux`, `flux_error`, `converged`, `P_interface_normalized` |
| `calculate_concentration_profile()` | **8 keys:** `x_oxide`, `C_oxide`, `x_metal`, `C_metal`, `x_all`, `C_all`, `P_interface`, `C_discontinuity` |
| `solve_interface_pressure_defective_metal()` | **14 keys:** the 7 above **+** `D_eff`, `D_lattice`, `modification_factor`, `level4_iterations`, `level4_converged`, `microstructure_details`, `temperature` |
| `calculate_oxide_metal_system()` | `solution` dict |
| `calculate_oxide_defective_metal_system()` | `solution` dict |

## 6. `parallel_oxide_defect_paths.py` — defect paths (Level 3)

| Function | Returns |
|---|---|
| `calculate_defect_path_flux()` | scalar: `flux_defect` |
| `calculate_parallel_path_flux()` | **18 keys:** `flux_total`, `flux_intact_contribution`, `flux_defect_contribution`, `flux_intact_per_area`, `flux_defect_per_area`, `dominant_path`, `defect_enhancement_factor`, `area_fraction_defect`, `D_eff_metal`, `modification_factor`, `level4_converged`, `P_interface_intact`, `regime_intact`, `regime_classification`, `regime`, `regime_base`, `regime_detail`, `microstructure_details` |
| `calculate_parallel_path_flux_defective_metal()` | same **18 keys** |
| `calculate_PRF()` | **8 keys:** `PRF`, `PRF_perfect`, `efficiency`, `regime`, `test_pressure`, `flux_bare_metal`, `flux_with_oxide`, `flux_reduction_factor` |
| `calculate_PRF_defective_metal()` | **13 keys:** the 8 above **+** `D_eff_bare_metal`, `modification_factor_bare`, `D_eff_with_oxide`, `modification_factor_oxide`, `microstructure_details` |

## 7. `surface_kinetics.py` — surface steps (Level 6)

| Function | Returns |
|---|---|
| `get_all_properties()` | **8 keys:** `k_diss`, `k_recomb`, `K_eq`, `D_ox`, `K_ox`, `L_ox`, `D_m`, `K_s_m` |
| `compute_permeances()` | tuple: `alpha`, `beta` |
| `g_theta()`, `sqrt_P_int_from_theta()`, `smooth_surface_resistance()` | scalars |
| `surface_flux()`, `oxide_flux()`, `metal_flux()` | scalar fluxes |
| `solve_steady_state_flux_L1L6()` | **6 keys:** `theta`, `P_int`, `J_ss`, `beta`, `rate_limiting`, `resistances` |
| `solve_steady_state_flux_L2aL6()` | **6 keys:** `theta`, `P_surf`, `J_ss`, `alpha`, `rate_limiting`, `resistances` |
| `solve_steady_state_flux()` | **10 keys:** `theta`, `P_int`, `J_ss`, `J_surface`, `J_oxide`, `J_metal`, `alpha`, `beta`, `rate_limiting`, `resistances` |
| `solve_steady_state_flux_direct()` | same **10 keys** |
| `calculate_defective_metal_flux_L6()` | **22 keys:** `flux`, `J_surface`, `J_oxide`, `J_metal`, `theta_surface`, `P_int`, `C_up`, `C_down`, `alpha`, `beta_lattice`, `beta_eff`, `D_eff`, `D_lattice`, `modification_factor`, `microstructure_details`, `flux_balance`, `profiles`, `rate_limiting`, `resistances`, `convergence`, `error`, `convergence_history` |
| `calculate_path_flux_L6()` | **11 keys:** `flux`, `theta`, `P_int`, `path_type`, `alpha`, `beta`, `kinetics_used`, `flux_balance`, `resistances`, `rate_limiting`, `error` |
| `calculate_parallel_path_flux_L6()` | **16 keys:** `J_total`, `enhancement_factor`, `dominant_path`, `fraction_intact`, `fraction_defect`, `flux_from_intact`, `flux_from_defect`, `fraction_from_defect`, `intact_path`, `defect_path`, `alpha_intact`, `alpha_defect`, `alpha_ratio`, `defect_type`, `thickness_factor`, `diffusivity_factor` |
| `calculate_mixed_defect_flux_L6()` | **14 keys:** `J_total`, `enhancement_factor`, `dominant_path`, `dominant_fraction`, `flux_breakdown`, `fraction_intact`, `total_defect_fraction`, `intact_path`, `defect_paths`, `system_rate_limiting`, `dominant_resistances`, `alpha_intact`, `defect_config`, `error` |
| `calculate_path_flux_L346_v2()` | **19 keys:** `flux`, `theta`, `P_int`, `path_type`, `kinetics_used`, `alpha`, `beta_lattice`, `beta_eff`, `D_eff`, `D_lattice`, `modification_factor`, `microstructure`, `flux_balance`, `resistances`, `rate_limiting`, `convergence`, `profiles`, `error`, `iteration` |
| `calculate_full_model_flux_L346_v2()` | **18 keys:** `J_total`, `enhancement_vs_intact`, `dominant_path`, `dominant_fraction`, `flux_breakdown`, `fraction_intact`, `total_defect_fraction`, `intact_path`, `defect_paths`, `D_eff_avg`, `D_lattice`, `overall_modification_factor`, `system_rate_limiting`, `system_resistances`, `alpha_intact`, `defect_config`, `microstructure_params`, `flux_weighted_resistances` |

## 8. `classify_regime.py` — regime labels

| Function | Returns |
|---|---|
| `classify_regime_level2()` | **6 keys:** `model_level`, `base_regime`, `regime_hierarchy`, `regime_detail`, `regime_subdetail`, `classification_depth` |
| `classify_regime_level3()` | the 6 above **+** `flux_ratio_defect_to_intact` |
| `classify_regime_level4_metal()` | **4 keys:** `metal_regime_detail`, `modification_factor`, `trapping_significant`, `trapping_reduction_percent` |
| `classify_regime_level14()` | **9 keys** — the 6 base **+** `modification_factor`, `trapping_significant`, `trapping_reduction_percent` |
| `classify_regime_level24()` | same **9 keys** |
| `classify_regime_level34()` | **10 keys** — the 9 above **+** `flux_ratio_defect_to_intact` |

**Regime label vocabulary.** `base_regime`: `surface_limited`, `oxide_limited`, `metal_limited`. `regime_detail`: `defect_limited`, `regime_intact_oxide`, `traps_defect_limited`, `lattice_limited`. `regime` (argmax, used by the SA): `surface`, `oxide`, `metal`, `defect`.

## 9. `sensitivity.py` — SA outputs

### `givendata_sensitivity_by_regime()` → nested `results[regime][metric]`

Each entry: `n` (cluster size), `log` (whether log₁₀ was applied), `pawn`, `delta`, `floor`, and `skipped` if the cluster was too small.

### `summarize_givendata()` → tidy DataFrame, 9 columns

`regime`, `metric`, `parameter`, `pawn_median`, `delta`, `S1_givendata`, `n`, `floor`, `delta_over_floor`

`delta_over_floor` ≤ 1 means **not resolved above the dummy-parameter noise floor** — the go/no-go column for reading any ranking.

Metrics available: **L5** `['flux', 'permeability', 'PRF']`, **L5L6** `['flux', 'permeability', 'theta']`.

### Other returns

| Function | Returns |
|---|---|
| `run_global_lhs_scan()` | tuple of 2 — samples DataFrame, records DataFrame |
| `run_targeted_regime_scans()` | tuple of 4 |
| `load_regime_scans()` | tuple of 3 |
| `partition_by_regime()` | `partition` dict, keyed by regime label |
| `size_draws_for_target()` | `sized` dict — draws per preset |
| `regime_comparison_matrix()` | `df`, `mat` |
| `check_against_config()` | `problems` list |
| `top_drivers()` | ranked parameter list |
| `plot_*()` (6 functions) | matplotlib `fig` |

---

## 10. Quick index — where to get a given quantity

| You want | Read |
|---|---|
| Permeation flux | `flux` / `J_ss` / `J_total` |
| Permeability | `permeability` |
| Barrier effectiveness | `PRF`, `PRF_perfect`, `efficiency`, `flux_reduction_factor` |
| Which step limits transport | `regime`, `rate_limiting`, `system_rate_limiting`, `limiting_mechanism` |
| Resistance split | `frac_surface`, `frac_oxide`, `frac_metal`, `frac_defect`, `resistances` |
| Effective diffusivity | `D_eff`, `D_modification`, `reduction_factor`, `mobile_fraction` |
| Trapping detail | `trapping_term`, `theta_total`, `trap_contributions`, `dominant_trap` |
| Surface coverage | `theta`, `theta_surface` |
| Spatial profiles | `profiles`, `x_all`, `C_all`, `C_array`, `D_array`, `theta_array` |
| Interface pressure | `P_interface`, `P_int`, `C_discontinuity` |
| Intact vs defect split | `flux_intact`, `flux_defect`, `dominant_path`, `flux_breakdown` |
| Numerical health | `converged`, `convergence`, `flux_error`, `flux_balance`, `error`, `warnings` |
| Sensitivity ranking | `delta`, `pawn_median`, `delta_over_floor` |
