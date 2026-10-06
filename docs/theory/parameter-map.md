# Parameter map

One quantity, many names. The metal's diffusivity is called `D` at Level 1, `D_m`
at Level 1+L6, `D_lattice` at Level 4, `D_metal` inside a properties dict, and
`D_ref` in the study configuration — five spellings of one number, and that is
before the sensitivity analysis renames the wall thicknesses again.

That is survivable. What is not survivable is the reverse: a single name carrying
more than one meaning. `thickness` is the metal in {doc}`level1-metal` and the
oxide in {doc}`level2-oxide`. `f_gb` is a volume fraction inside the metal in one
place and an area fraction of oxide defects in another, two quantities that differ
by more than two orders of magnitude. Nothing errors when those are confused.

This page inventories both, and proposes one canonical name per quantity.

```{warning}
**The canonical names below are a PROPOSAL. None of them exist in the code.**

Every name in a "canonical" column is what a quantity *would* be called under a
unified scheme. The names in the "replaces" columns are the ones that are real
today and that you must use. Nothing here has been renamed, and the tables are a
decision document, not a description of the codebase.
```

## How to read it

**Three naming layers.** The same physical quantity is named independently at
each, and a rename has to cross all three to be complete:

```{list-table}
:header-rows: 1
:widths: 22 78

* - layer
  - where it lives
* - function argument
  - the positional parameters of each level's entry point
* - dict key
  - inside `oxide_props`, `metal_props`, `microstructure_params`, `trap_list`
    entries and the defect configuration
* - study configuration
  - `CONDITIONS`, `METALS`, `OXIDES`, `MICROSTRUCTURE`, `OXIDE_DEFECTS`, and the
    flat `DEFAULT_PARAMS_LEVEL5` / `DEFAULT_PARAMS_LEVEL5L6` the sensitivity
    analysis varies
```

**Three suffixes, used consistently in the proposal.** A bare name is the value
*at the operating temperature*; `_ref` is the value *at the reference
temperature*; `_energy` or `_enthalpy` is the Arrhenius exponent that connects
them. The current code blurs this — `D_ref` and `D_lattice` are related but
genuinely different numbers, and no present naming makes that visible.

---

## A. Metal bulk transport

```{list-table}
:header-rows: 1
:widths: 34 30 36

* - canonical
  - quantity (units)
  - replaces
* - `metal_diffusivity`
  - lattice diffusivity at T, m²/s
  - `D` · `D_m` · `D_lattice` · `D_metal`
* - `metal_diffusivity_ref`
  - same, at reference T
  - `D_ref`
* - `metal_diffusivity_prefactor`
  - Arrhenius D₀, derived
  - `D_0`
* - `metal_diffusion_energy`
  - activation energy, J/mol
  - `E_D`
* - `metal_solubility`
  - Sieverts constant at T
  - `K_s` · `K_s_m` · `K_s_metal`
* - `metal_solubility_ref`
  - same, at reference T
  - `K_s_ref`
* - `metal_solubility_prefactor`
  - derived
  - `K_s0`
* - `metal_solution_enthalpy`
  - enthalpy of solution, J/mol
  - `H_s`
* - `metal_permeability_ref`
  - Φ = D·K_s at reference T
  - `Phi_ref`
* - `metal_permeability_prefactor`
  - derived
  - `Phi_0`
* - `metal_permeability_energy`
  - E_D + ΔH_s, derived
  - `Q_p`
* - `metal_thickness`
  - wall thickness, m
  - `thickness` · `L_m` · `L_metal` · `metal_thickness`
* - `metal_reference_temperature`
  - Arrhenius anchor, K
  - `T_ref` (in METALS) · `T_ref_metal`
```

## B. Oxide bulk transport

```{list-table}
:header-rows: 1
:widths: 34 30 36

* - canonical
  - quantity (units)
  - replaces
* - `oxide_diffusivity`
  - H₂ (Model 1) or H (Model 2) diffusivity at T, m²/s
  - `D_ox` — already uniform
* - `oxide_diffusivity_ref`
  - same, at reference T
  - `D_ox_ref`
* - `oxide_diffusivity_prefactor`
  - derived
  - `D_ox_0`
* - `oxide_diffusion_energy`
  - activation energy
  - `E_D_ox`
* - `oxide_solubility`
  - Henry (Model 1) or Sieverts (Model 2) constant at T
  - `K_ox` — already uniform
* - `oxide_solubility_ref`
  - same, at reference T
  - `K_ox_ref`
* - `oxide_solubility_prefactor`
  - derived
  - `K_ox_0`
* - `oxide_solution_enthalpy`
  - enthalpy of solution
  - `H_sol_ox`
* - `oxide_permeability_ref`
  - Φ_ox at reference T
  - `Phi_ox_ref`
* - `oxide_permeability_energy`
  - its activation energy
  - `Q_p_ox_J_per_mol` — the only name carrying its units
* - `oxide_thickness`
  - layer thickness, m
  - `thickness` · `L_ox` · `L_oxide` · `oxide_thickness`
* - `oxide_reference_temperature`
  - Arrhenius anchor, K
  - `T_ref` (in OXIDES) · `T_ref_oxide`
```

Which model the oxide is in decides the *units* of its solubility, not its name —
see {doc}`two-models`. A single symbol serving both laws is an unresolved problem
in its own right, recorded in {doc}`../caveats`.

## C. Operating point and sweep ranges

```{list-table}
:header-rows: 1
:widths: 34 30 36

* - canonical
  - quantity
  - replaces
* - `upstream_pressure`
  - charging-face H₂ pressure, Pa
  - `P_up` · `P_upstream`
* - `downstream_pressure`
  - exit-face pressure, Pa
  - `P_down` · `P_downstream`
* - `operating_temperature`
  - wall temperature, K
  - `temperature` · `T_K` · `temperature_K` · `T_operating`
* - `sweep_pressure_range`
  - SA span for upstream pressure
  - `P_range`
* - `sweep_temperature_range`
  - SA span for temperature
  - `T_range`
* - `sweep_metal_thickness_range`
  - SA span
  - `L_metal_range`
* - `sweep_oxide_thickness_range`
  - SA span
  - `L_oxide_range`
* - `sweep_defect_fraction_range`
  - SA span
  - `defect_fraction_range`
* - `sweep_dissociation_rate_range`
  - SA span
  - `k_diss_sweep_range`
* - `n_sweep_samples`
  - samples per SA axis
  - `n_points` (in CONDITIONS)
* - `defect_fractions_compared`
  - the discrete fractions a figure contrasts
  - `defect_fractions_compare`
* - `defect_limit_check`
  - the limit case the notebooks assert
  - unchanged
* - `max_gb_volume_fraction`
  - cap on the metal grain-boundary fraction
  - `f_gb_max`
```

## D. Surface kinetics — Model 2 only

```{list-table}
:header-rows: 1
:widths: 40 26 34

* - canonical
  - quantity
  - replaces
* - `oxide_dissociation_rate`
  - rate constant, gas/oxide surface
  - `k_diss` · `k_diss_ref`
* - `oxide_dissociation_energy`
  - its activation energy
  - `E_diss`
* - `oxide_surface_equilibrium`
  - adsorption/absorption equilibrium constant
  - `K_eq` · `K_eq_ref`
* - `oxide_surface_enthalpy`
  - its reaction enthalpy
  - `H_eq`
* - `oxide_surface_reference_temperature`
  - Arrhenius anchor
  - `T_ref` (oxide surface block) · `T_ref_surface`
* - `metal_dissociation_rate`
  - rate constant, bare metal at a pinhole
  - `k_diss_metal` · `k_diss_metal_ref`
* - `metal_dissociation_energy`
  - its activation energy
  - `E_diss_metal`
* - `metal_surface_equilibrium`
  - equilibrium constant, bare metal
  - `K_eq_metal` · `K_eq_metal_ref`
* - `metal_surface_enthalpy`
  - its reaction enthalpy
  - `H_eq_metal`
* - `metal_surface_reference_temperature`
  - Arrhenius anchor
  - `T_ref` (metal surface block) · `T_ref_surface_metal`
```

```{note}
This is the most systematic block in the configuration as it stands. The
oxide/metal pairing is consistent throughout — one surface per exposed material,
five parameters each — and only the `_ref` suffix wobbles between the nested and
flattened forms. At the configuration layer the two models are cleanly nested:
`DEFAULT_PARAMS_LEVEL5L6` is a strict superset of `DEFAULT_PARAMS_LEVEL5`, adding
exactly these ten names and removing none.
```

## E. Metal microstructure and traps

```{list-table}
:header-rows: 1
:widths: 34 30 36

* - canonical
  - quantity
  - replaces
* - `grain_size`
  - mean grain diameter, m
  - unchanged
* - `gb_width`
  - grain-boundary thickness, m
  - `gb_thickness`
* - `gb_character`
  - HAGB / LAGB / twin / special
  - `gb_type`
* - `grain_morphology`
  - equiaxed / columnar / planar
  - `grain_shape`
* - `lattice_site_density`
  - H sites per m³
  - `lattice_density` · `N_L` (a deliberate alias)
* - `metal_gb_volume_fraction`
  - 3δ/d, derived
  - `f_gb`
* - `gb_diffusivity_ratio`
  - D_gb / D_lattice
  - `alpha` (local) · `enhancement_factor`
* - `trap_vacancy_density`
  - sites per m³
  - `trap_vacancy_N_T`
* - `trap_carbide_density`
  - sites per m³
  - `trap_carbide_N_T`
* - `trap_gb_density`
  - sites per m³
  - `trap_gb_N_T`
* - `trap_dislocation_density`
  - sites per m³
  - `trap_dislocation_N_T`
* - `trap_<name>_binding_energy`
  - trap depth, J/mol
  - `trap_*_E_b`
* - `trap_density` · `trap_binding_energy` · `trap_name`
  - inside a `trap_list` entry
  - `density` · `binding_energy` · `name`
```

## F. Oxide defects

```{list-table}
:header-rows: 1
:widths: 34 30 36

* - canonical
  - quantity
  - replaces
* - `pinhole_area_fraction`
  - area share with no oxide
  - `f_pinhole` · `components['pinholes']`
* - `crack_area_fraction`
  - area share with thinned oxide
  - `f_crack` · `components['cracks']`
* - `oxide_gb_area_fraction`
  - area share through oxide boundaries
  - `f_gb_defect` · `components['grain_boundaries']`
* - `total_defect_area_fraction`
  - their sum
  - `area_fraction` · `defect_area_fraction`
* - `crack_thickness_ratio`
  - L_crack / L_oxide
  - `thickness_factor` · `crack_thickness_factor`
* - `oxide_gb_diffusivity_ratio`
  - D_gb / D_oxide
  - `diffusivity_factor` · `gb_diffusivity_factor` · `beta` (in `calculations.permeability`)
* - `defect_population`
  - which defect type a call models
  - `type` · `defect_type`
* - `pinhole_uses_sieverts`
  - whether a pinhole takes the Sieverts limit
  - `use_sieverts_pinhole`
```

## G. Model switches and numerics

```{list-table}
:header-rows: 1
:widths: 34 30 36

* - canonical
  - quantity
  - replaces
* - `enable_gb_enhancement`
  - grain-boundary **fast path** on/off
  - `include_gb_enhancement`
* - `enable_trapping`
  - trapping on/off
  - `include_trapping`
* - `enable_gb_as_trap`
  - auto-add grain boundaries as a **trap** population
  - `include_gb_trapping`
* - `microstructure_mode`
  - the composed switch the solvers take
  - `mode`
* - `d_eff_averaging`
  - average / harmonic / inlet / outlet
  - `method` · `D_eff_method`
* - `n_profile_points`
  - wall discretisation
  - `n_points` (function argument)
* - `max_iterations` · `tolerance`
  - outer-loop convergence controls
  - unchanged
```

## H. Provenance — not physics

These describe where a number came from. Mistaking one for a model input is a
recorded error in this project: a sample thickness is the specimen the data were
measured on, never the wall being simulated.

```{list-table}
:header-rows: 1
:widths: 34 30 36

* - canonical
  - quantity
  - replaces
* - `data_sample_thickness`
  - specimen thickness the source data used
  - `sample_thickness`
* - `data_temperature_range`
  - validity window of the source data
  - `temperature_range`
* - `data_pressure_range`
  - pressure validity of the source data
  - `pressure_range` (in METALS)
* - `data_pressure`
  - the pressure the source measured at
  - `pressure`
* - `gas_species`
  - H₂ / D₂ / T₂
  - `gas`
* - `source_reference`
  - citation
  - `reference`
* - `notes`
  - free text
  - `Notes`
```

---

## I. Collisions — one name, more than one meaning

This is the half that is a defect rather than a tidiness problem. An alias costs a
reader a moment; a collision silently produces the wrong number.

```{list-table}
:header-rows: 1
:widths: 18 40 42

* - current
  - splits into
  - why it must split
* - `thickness`
  - `metal_thickness` · `oxide_thickness` · `layer_thickness`
  - the metal in `calculate_simple_metal_flux` and `calculate_defective_metal_flux`,
    the oxide in `molecular_diffusion_flux` and `calculate_oxide_resistance`, and
    **genuinely generic** in `fick_flux`, which is applied to both layers
* - `f_gb`
  - `metal_gb_volume_fraction` · `oxide_gb_area_fraction`
  - a volume fraction inside the metal against an area fraction of oxide defects.
    For the shipped studies they differ by more than two orders of magnitude
* - `n_points`
  - `n_sweep_samples` · `n_profile_points`
  - the sensitivity analysis's samples per axis against the wall discretisation,
    whose default differs per level
* - `beta`
  - `metal_permeance` · `oxide_gb_diffusivity_ratio`
  - a permeance carrying units in `calculations.surface_kinetics`, a dimensionless
    ratio in `calculations.permeability`
* - `alpha`
  - `oxide_permeance` · `gb_diffusivity_ratio`
  - the same problem, for the other layer
* - `T_ref`
  - four scoped reference temperatures
  - one name in four dicts. The flattened sensitivity layer already splits it
    correctly; the nested configuration relies on scope
* - `include_gb_*`
  - `enable_gb_enhancement` · `enable_gb_as_trap`
  - **near-identical names, opposite physics.** Enhancement speeds diffusion up
    through a fast path; trapping slows it down by holding hydrogen still
* - pressure ranges
  - `sweep_pressure_range` · `data_pressure_range`
  - `P_range` is the span the sensitivity analysis sweeps; `pressure_range` is the
    window the source data are valid over. Unrelated, nearly homonymous
```

## Totals

**117 distinct names in use today → 86 under the proposal.** Thirty-one
eliminated as aliases, eight split because they carry more than one meaning.

---

## What a rename would cost

Measured rather than estimated, against the codebase as of 2026-10-06.

### Where the work is

`thickness` is the most-feared case and the most tractable. Classifying every
occurrence by its syntactic role shows that nearly all of it self-identifies:

```{list-table}
:header-rows: 1
:widths: 46 14 40

* - role
  - count
  - attributable by
* - `metal_props['thickness']`
  - 25
  - the dict — metal
* - `oxide_props['thickness']`
  - 19
  - the dict — oxide
* - dict-literal key
  - 19
  - which dict is being built
* - bare name use
  - 25
  - the enclosing function's parameter
* - parameter declaration
  - 10
  - the function
* - keyword argument
  - 12
  - the callee
* - `crack_oxide_props['thickness']`
  - 4
  - the dict — oxide, a thinned copy
```

Forty-eight sites resolve from the dict they are indexed on, and the ten
declarations split cleanly into seven metal and two oxide. **The residue needing
human judgement is one function**: `fick_flux` takes a `thickness` that is
layer-agnostic by design, and should become `layer_thickness` rather than either
canonical name.

### Which layer is expensive

The notebooks invert the intuition. Configuration names are everywhere; solver
argument names are not.

```{list-table}
:header-rows: 1
:widths: 30 24 46

* - name
  - layer
  - notebooks containing it
* - `lattice_density`
  - configuration
  - 25 of 25
* - `D_ref`
  - configuration
  - 25 of 25
* - `P_upstream`
  - configuration
  - 24 of 25
* - `thickness`
  - solver
  - 10 of 25
* - `D_lattice`
  - solver
  - 9 of 25
* - `L_m`
  - solver
  - 9 of 25
```

So renaming function arguments is comparatively cheap, and renaming configuration
keys touches every notebook in the repository.

### One rename path fails silently

```{danger}
A model wrapper absorbs an unrecognised parameter name without complaint.
`level5_model_wrapper` copies the defaults and applies `.update(params_dict)`, so
a key that no longer exists is stored and then never read. Verified: passing a
stale name alongside no other change returns a flux bit-identical to the baseline,
with no error and no warning.

`check_against_config` does **not** close this. It demands a default for every
swept range, and it cross-compares the Level 5 and Level 5+L6 preset blocks
against each other — but it never checks that a preset key exists in
`DEFAULT_PARAMS`. A rename applied consistently to both preset dicts and not to
the defaults passes clean.

The exposure is the regime presets: 36 distinct parameter names pinned across the
three Level 5+L6 regimes, 28 across the Level 5 regimes. None are stale today. A
preset that silently stops pinning its parameter does not fail a test — it quietly
degrades a targeted regime scan into an untargeted one, which is the whole basis
of the regime-stratified workflow.

Eight defaults have no sweep range at all, so even the working half of the check
does not cover them: `D_eff_method`, `P_downstream`, `T_ref_metal`,
`T_ref_oxide`, `gb_type`, `grain_shape`, `include_gb_enhancement`,
`include_trapping`.
```

### Suggested order

1. **Add a preset-key guard to `check_against_config`** — assert every regime
   preset key exists in the corresponding defaults. A few lines, passes on the
   configuration as shipped, and converts the only silent failure mode into a loud
   one. Worth doing whether or not a rename follows, because a *typo* in a preset
   name is indistinguishable from a stale one today.
2. **Fix the collisions in section I.** Small, self-contained, and the only part
   that can silently produce a wrong number. It does not require touching the
   hundreds of alias sites.
3. **The alias rename in sections A–H, if at all.** Cheapest at the solver layer,
   most expensive at the configuration layer, and of no scientific consequence —
   this is readability, weighed against churn in 25 notebooks.
