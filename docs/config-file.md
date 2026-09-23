# The study configuration file

Every number the model runs on comes from one file. This page explains that
file's structure: what each block holds, which blocks you are meant to edit,
which are derived from the ones above them, and which entries the code never
reads at all.

It describes the *shape* of the file, not its values — for the values of the
active study see {doc}`PARAMETERS`. To change which study is active, see
{doc}`how-to/switch-study`; for the loading mechanism, {doc}`api/config`.

## One study, one file

A study is a complete physical system — a metal, an oxide, a microstructure, a
set of defects and an operating point — and each one lives in a single module
under
[`calculations/config/studies/`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/tree/main/calculations/config/studies).
Nothing in `calculations/` imports a study directly:

```text
calculations/*.py  ->  calculations.config.model_config  ->  studies/<the active one>.py
```

`model_config` is a switch, not a config. It holds one line naming the active
study and re-exports that module's public names into its own namespace, so every
module and notebook reads the same material through the same import.

The study modules are **structurally identical**. They define the same names in
the same order and differ only in their values, which is what lets one switch
line replace an entire material system.

## Three kinds of entry

The single most useful thing to know about the file is that its blocks are not
all the same kind of thing. Reading them the same way is how the config drifted
before.

```{list-table}
:header-rows: 1
:widths: 20 40 40

* - Kind
  - Blocks
  - How to treat it
* - **Measured data**
  - `METALS`, `OXIDES`, `MICROSTRUCTURE`, `OXIDE_DEFECTS`, `CONDITIONS`
  - Authored by hand from the literature, each with a `reference` string. This
    is what you edit when you change material.
* - **Derived**
  - `DEFAULT_PARAMS_LEVEL5`, `DEFAULT_PARAMS_LEVEL5L6`
  - Computed *from the blocks above* by reference, never retyped. Flat
    name → value dicts, which is the form the model wrappers want.
* - **Study design**
  - `SUGGESTED_RANGES_*`, `REGIME_PRESETS*`, `DEFAULT_SWEEP_TEMPERATURES`
  - Deliberately literal. A range that suppresses `f_pinhole` to force an
    oxide-limited cluster expresses intent about the analysis, not a property
    of the material, so it is written out explicitly.
```

```{warning}
**Never hardcode a number into `DEFAULT_PARAMS_*`.** Those dicts exist to be
derived; a literal there is a second copy of a value that already lives in
`METALS`/`OXIDES`/`CONDITIONS`, and the two will diverge. They already did once:
when the defaults were hand-copied, the defect area fractions, the operating
temperature and both surface equilibrium constants disagreed with the material
blocks above them by up to six orders of magnitude, silently. Change the source
dict instead.
```

## The blocks, in file order

```{list-table}
:header-rows: 1
:widths: 26 40 34

* - Block
  - What it describes
  - Read by
* - `R`, `F`
  - gas constant; Faraday constant, used to write trap energies in eV
  - the file itself
* - `METALS`
  - bulk metal transport + metal-surface kinetics, per alloy
  - `calculations.oxide_permeation`, `calculations.surface_kinetics`
* - `OXIDES`
  - bulk oxide transport + gas/oxide-surface kinetics, per oxide
  - `calculations.oxide_permeation`, `calculations.surface_kinetics`
* - `MICROSTRUCTURE`
  - grains, grain boundaries and trap populations (Level 4)
  - `calculations.sensitivity`, notebooks
* - `OXIDE_DEFECTS`
  - pinhole / crack / grain-boundary bypass paths (Level 3)
  - `calculations.sensitivity`, notebooks
* - `CONDITIONS`
  - the operating point and the sweep ranges around it
  - `calculations.sensitivity`, notebooks
* - `DEFAULT_PARAMS_LEVEL5`, `DEFAULT_PARAMS_LEVEL5L6`
  - the same data flattened into one name → value dict per model
  - `calculations.sensitivity`
* - `SUGGESTED_RANGES_LEVEL5`, `SUGGESTED_RANGES_LEVEL5L6`
  - the interval each parameter is sampled over
  - `calculations.sensitivity`
* - `REGIME_PRESETS`, `REGIME_PRESETS_L5`
  - ranges overridden to target one rate-limiting regime
  - `calculations.sensitivity`
* - `PLOT_STYLE`, `COLORS`, `CURVE_STYLES`
  - figure styling, one entry per model level
  - notebooks, via `apply_style`
* - helper functions
  - `get_surface_kinetics_from_config`, `build_simulation_config`, `apply_style`
  - notebooks and `calculations.surface_kinetics`
* - `VALIDATION`
  - per-level test grids for the validation notebooks
  - notebooks
```

### `METALS`

A registry keyed by alloy name. Each entry carries transport data at a reference
temperature, a nested surface-kinetics block, and provenance.

```python
METALS = {
    '<alloy_key>': {
        'T_ref':   1223,        # K — the temperature the *_ref values were measured at
        'D_ref':   5.285e-09,   # m²/s          lattice diffusivity at T_ref
        'E_D':     7220,        # J/mol         diffusion activation energy
        'K_s_ref': 0.05093,     # mol/m³/Pa^0.5 Sieverts constant at T_ref
        'H_s':     7200,        # J/mol         heat of solution
        ...
        'surface_kinetics': {...},
        'reference': '...',
    }
}
```

Properties are evaluated by `get_metal_properties_at_T` in the **reference-point
Arrhenius** form, which is why a `T_ref` accompanies every `*_ref` value:

$$X(T) = X_{\text{ref}}\exp\!\left[\frac{-E}{R}\left(\frac{1}{T} -
\frac{1}{T_{\text{ref}}}\right)\right]$$

Only five keys are read: `T_ref`, `D_ref`, `E_D`, `K_s_ref`, `H_s`. The rest —
`D_0`, `K_s0`, `Phi_ref`, `Phi_0`, `Q_p` — are **recorded, not used**. They are
the same physics in the absolute-Arrhenius form, kept for traceability to the
source paper. Editing `Phi_ref` changes nothing the model computes.

The remaining fields are metadata: `reference`, `temp_range`, `pressure_range`,
`pressure`, `gas`, `Notes`, and `metal_thickness`.

```{warning}
**`METALS[...]['metal_thickness']` is a list, and is not the wall thickness.**
It records the specimen thicknesses of the source experiment, e.g. `[5e-3]`. The
wall thickness the model uses comes from `CONDITIONS['L_metal']`, and
`DEFAULT_PARAMS_LEVEL5` reads it from there for exactly this reason. Wiring the
`METALS` entry into a model call passes a *list* where a float is expected and
poisons every result; `check_against_config` exists partly to catch that class
of mistake.
```

### `OXIDES`

Same shape as `METALS`, one entry per oxide, with `D_ox_ref` / `E_D_ox` /
`K_ox_ref` / `H_sol_ox` in place of the metal quartet, plus two extra fields the
metal does not have:

- `thickness` — the oxide layer thickness the model uses, and
- `temperature_range` — the validity window of the source data, which
  `get_oxide_properties_at_T` checks and warns about on extrapolation.

As with the metal, `D_ox_0`, `K_ox_0`, `Phi_ox_ref` and `Q_p_ox_J_per_mol` are
recorded but never read.

```{warning}
**Two oxide thicknesses coexist in the file.** `OXIDES[...]['thickness']` and
`CONDITIONS['L_oxide']` are independent entries and in the shipped studies they
differ by more than an order of magnitude. The flat parameter dicts take the
oxide thickness from the **oxide** entry; `CONDITIONS['L_oxide']` is used by
thickness sweeps. Know which one your call path reads before quoting a barrier
thickness.
```

```{note}
`K_ox_ref` means different things in the two models — a Henry constant under
Model 1 and a Sieverts constant under the Level 6 family — and the two cannot
share units. The config carries one entry for both. See
{doc}`theory/equilibrium-models`.
```

### `surface_kinetics` (nested in both)

Each material entry carries its own dissociation/recombination block, evaluated
in the same reference-point Arrhenius form, with `k_recomb = k_diss / K_eq`
derived rather than stored. **The two use different key names:**

| Oxide entry (gas/oxide surface) | Metal entry (pinhole paths) |
|---|---|
| `k_diss_ref` | `k_diss_metal_ref` |
| `E_diss` | `E_diss_metal` |
| `K_eq_ref` | `K_eq_metal_ref` |
| `H_eq` | `H_eq_metal` |

`get_surface_kinetics_from_config` discriminates on the presence of
`k_diss_ref`: found, and it reads the oxide names; absent, it falls through to
the metal names. A typo in an oxide key therefore does not fall back silently —
it routes to the metal branch and raises `KeyError` on a name that is not there.

Check the comments before trusting the equilibrium constants: the oxide's
`K_eq_ref` and `H_eq` are marked `PLACEHOLDER` in every shipped study, and the
metal's `K_eq_metal_ref` and `H_eq_metal` are placeholders in all but the 316L
study. They are order-of-magnitude stand-ins, not measurements.

### `MICROSTRUCTURE`

Level 4 inputs: `grain_size`, `gb_thickness`, `grain_shape`, `gb_type`,
`include_gb_trapping`, `lattice_density`, and `trap_list` — a *list* of
`{name, binding_energy, density}` dicts, one per trap population, with binding
energies written as `0.27 * F` so the eV value from the source paper stays
readable.

Two consequences of it being a list rather than a dict: order does not matter,
and the flattening step looks traps up by name (case-insensitively, since the
carbide population is capitalised in some studies). A renamed trap raises a
`KeyError` naming the populations that do exist.

Immediately after the block, `MICROSTRUCTURE['N_L']` is aliased to
`lattice_density` — the same quantity under the symbol the derivations and
notebooks use, aliased rather than written twice so the two cannot drift.

### `OXIDE_DEFECTS`

The Level 3 bypass paths. `components` holds the three area fractions
(`pinholes`, `cracks`, `grain_boundaries`); `thickness_factor` sets the crack
short-circuit, $L_{\text{crack}} = \gamma L_{\text{oxide}}$; and
`diffusivity_factor` sets the grain-boundary enhancement,
$D_{\text{gb}} = \delta D_{\text{oxide}}$.

`use_sieverts_pinhole` is a model *switch*, not data: it chooses whether a
pinhole assumes fast metal kinetics or uses the metal's own surface-kinetics
block. It is deliberately absent from the flat parameter dicts — the notebooks
read it from the assembled simulation config, so `OXIDE_DEFECTS` is its only
home.

### `CONDITIONS`

The operating point (`T_operating`, `P_upstream`, `P_downstream`, `L_metal`,
`L_oxide`) together with the sweep ranges and point counts the notebooks plot
over (`T_range`, `P_range`, `L_metal_range`, `L_oxide_range`, `n_*_points`).
Single values and sweep bounds sit in one block, distinguished only by the
`_range` suffix.

### The active-material shorthand

```python
ACTIVE_METAL = _first_key(METALS, 'METALS')
ACTIVE_OXIDE = _first_key(OXIDES, 'OXIDES')
```

Everything below this line refers to the material **positionally** — the first
entry of each registry — so adding or renaming a material needs no edit further
down, and a key name cannot fall out of sync with the dict it indexes.

```{warning}
**The first entry wins, so add new materials at the bottom.** Dicts preserve
insertion order, so `ACTIVE_METAL` is the first literal in the `METALS` block.
Inserting a second alloy *above* the current one silently re-points every
derived default, every sweep range and every model call at the new entry while
the rest of the file looks unchanged.
```

### `DEFAULT_PARAMS_LEVEL5` and `DEFAULT_PARAMS_LEVEL5L6`

The flattening step. The nested, physically-organised blocks above become one
flat `name → scalar` dict, because that is what the model wrappers and the
sensitivity sampler consume. Every value is a *reference* to a block above —
`_SA_M['D_ref']`, `_SA_DC['pinholes']`, `CONDITIONS['L_metal']` — so the flat
view cannot disagree with the structured one.

`DEFAULT_PARAMS_LEVEL5L6` starts with `**DEFAULT_PARAMS_LEVEL5` and adds the ten
surface-kinetics entries, which is what keeps the two models sharing one value
for everything they have in common. The split mirrors the model fork described
in {doc}`theory/two-models`.

A few entries have no config source because they are not physics: the flags
`include_gb_enhancement`, `include_trapping` and `D_eff_method` say how to run
the model, and are written literally.

### `SUGGESTED_RANGES_LEVEL5` / `..._LEVEL5L6`

One `[lo, hi]` per sampled parameter. The convention is roughly ±2 decades
around the reference value for prefactors, and physically motivated bounds for
activation energies. `SUGGESTED_RANGES_LEVEL5L6` inherits all of the Level 5
ranges and adds the surface block.

The `T_ref_*` entries are deliberately **not** sampled: they are measurement
metadata, redundant with the `*_ref` prefactors they accompany, and varying them
injects an artificial prefactor sweep of several decades.

### `REGIME_PRESETS` and `REGIME_PRESETS_L5`

Each preset is the base ranges with only the regime-controlling parameters
overridden, composed from named building blocks (`_SA_SLOW_OXIDE`,
`_SA_FAST_METAL`, `_SA_DEFECTS_OFF`, `_SA_DEFECTS_LARGE`, `_SA_FAST_SURFACE`,
`_SA_SLOW_SURFACE`) so that the parts both models share are defined exactly
once and cannot drift apart.

Targeting is not a convenience. Over the default ranges the full model yields
essentially no surface-limited samples and is overwhelmingly metal-limited —
one regime, nothing to stratify. The presets are what create the study. See
{doc}`how-to/sensitivity-analysis`.

### `DEFAULT_SWEEP_TEMPERATURES`

The temperatures of the isothermal sweep. Temperature is pinned rather than
sampled because, allowed to vary, it monopolises the variance of the log flux
and pushes every other parameter to the noise floor.

### `PLOT_STYLE`, `COLORS`, `CURVE_STYLES`

Presentation only, but part of the export contract, so a study file must define
all three. `CURVE_STYLES` gives every level and mode a combination of colour,
linestyle and marker that is distinguishable by any one of the three alone;
`apply_style` applies one by name.

### Helper functions

- `get_surface_kinetics_from_config` — evaluates either surface block at a
  temperature and returns `k_diss`, `k_recomb`, `K_eq`.
- `build_simulation_config` — assembles one solver-ready dict for a metal/oxide
  pair, merging the material entries with `CONDITIONS` and translating the
  `OXIDE_DEFECTS` fractions into the per-path `defect_config` the Level 3
  solvers expect. A path with a zero area fraction is *omitted*, not
  zero-weighted, so setting a fraction to zero removes that transport path.
- `apply_style` — plots a curve with a named entry from `CURVE_STYLES`.

```{note}
`build_simulation_config` applies its overrides with `or`, so a falsy `0` takes
the default instead. Pass a small positive number rather than zero.
```

### `VALIDATION`

Test grids for the validation notebooks, keyed by level (`L1`, `L2a`, `L2b`,
`L3`, `L4`, `L5`, `L6`): the pressures, temperatures and sweep arrays each
level's checks run over, including the fast-kinetics limit `L6` uses to verify
that each Level 6 model collapses onto its non-surface counterpart. Nothing in
`calculations/` reads this block.

## The export contract

A study file is valid only if it defines all fifteen names in
`REQUIRED_EXPORTS`:

```text
DEFAULT_PARAMS_LEVEL5     DEFAULT_PARAMS_LEVEL5L6
SUGGESTED_RANGES_LEVEL5   SUGGESTED_RANGES_LEVEL5L6
REGIME_PRESETS            REGIME_PRESETS_L5
METALS   OXIDES   MICROSTRUCTURE   OXIDE_DEFECTS   CONDITIONS
VALIDATION   build_simulation_config
PLOT_STYLE   COLORS   CURVE_STYLES
```

The list is derived from what `calculations/` actually imports, not from what
one study happens to define. The first six are pulled at module level by
`calculations.sensitivity`, so a study missing one of them would otherwise fail
at `import calculations.sensitivity`, before any error message could help. The
loader checks the set at import and raises a `ValueError` naming both the study
and the gap.

## After editing a study file

```python
from calculations.sensitivity import check_against_config
check_against_config()
```

Derivation guarantees that the flat defaults match the material blocks. This
checks the four things it does not: that every derived value is a finite scalar
(catching an entry that has quietly become a list or a range), that every
sampled parameter has a default, that every default lies strictly *inside* its
own sweep range — a default sitting exactly on a bound is only half-sampled —
and that the two models still agree wherever their defaults, ranges and presets
overlap.

Restart the kernel afterwards. The study is bound into `model_config` at import
time, so a reload of a notebook module is not enough.
