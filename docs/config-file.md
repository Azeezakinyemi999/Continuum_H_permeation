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
set of defects and an operating point — and each one lives in a `study_config.py`
**inside its own study folder**, beside the notebooks that read it. That file is
the real input file, not a pointer to one. Nothing in `calculations/` imports a
study directly:

```text
calculations/*.py  ->  calculations.config.model_config  ->  <study folder>/study_config.py
```

There is no default study and no switch to edit. `model_config` resolves one, in
this order:

1. `$MHI_STUDY`, if set, naming a study file explicitly;
2. otherwise a `study_config.py` in the working directory, or any parent up to the
   repository root;
3. otherwise it raises, listing the studies it can see.

So a notebook picks up its own folder's file by virtue of running there, and
running from the repository root — which belongs to no study — is an error rather
than a silent default. Each study is self-contained: `Guo_etal_2025/` has nothing
to do with `Fuerst_et_al_2024/`, and adding a study means copying a folder, not
registering it anywhere.

The study files are **structurally identical**. They define the same names in the
same three sections and differ only in their values, which is what lets one folder
replace an entire material system.

## The three sections

The file is ordered by dependency, and the section banners say so. Python executes
top to bottom, so nothing may reference a block below it.

| Section | Holds |
|---|---|
| 1. Identity and constants | `STUDY_NAME`, `R`, `F` |
| 2. What is being studied | `CONDITIONS`, `MODEL_OPTIONS`, `OXIDE_DEFECTS`, the `SUGGESTED_RANGES_*`, the presets, `DEFAULT_SWEEP_TEMPERATURES` |
| 3. What the materials are | `METALS`, `OXIDES`, `MICROSTRUCTURE` |

The split between 2 and 3 is the one to internalise. Section 3 says what a
substance *does* at a temperature; section 2 says what wall is being built from it
and at what operating point. A material entry therefore carries no thickness and
no operating point — see `sample_thickness` below for the one place that reads
like an exception and is not.

```{note}
There used to be two more sections — **4. Derived** (the `DEFAULT_PARAMS_*`,
`build_simulation_config`, `get_surface_kinetics_from_config`) and
**5. Presentation** (`PLOT_STYLE`, `COLORS`, `CURVE_STYLES`, `apply_style`). They
were byte-identical in all four studies, 408 lines duplicated four ways, and they
have moved into
[`calculations/config/model_config.py`](https://github.com/Azeezakinyemi999/MHI_permeation/blob/main/calculations/config/model_config.py).

You still import those names from exactly the same place — `model_config`
re-exports everything — so nothing you write changes. What changed is that they are
no longer yours to copy: a new study inherits them, and "a literal in the derived
block is a bug" is now true by construction instead of by four separate checks on
four separate copies. The copies had already drifted; the Level 5 header comment
named Incoloy 802 in the 316L study.
```

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
  - `METALS`, `OXIDES`, `MICROSTRUCTURE`, `OXIDE_DEFECTS`
  - Authored by hand from the literature, each with a `reference` string. This
    is what you edit when you change material.
* - **The scenario**
  - `CONDITIONS`, `MODEL_OPTIONS`
  - The operating point and the switches. Not a property of any material — this
    is what you edit to ask a different question of the same materials.
* - **Derived** — *no longer in this file*
  - `DEFAULT_PARAMS_LEVEL5`, `DEFAULT_PARAMS_LEVEL5L6`
  - Computed *from the blocks above* by reference, never retyped. Flat
    name → value dicts, which is the form the model wrappers want. Now held once
    in `model_config`, and listed here because they are still what the two blocks
    above feed.
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
* - `STUDY_NAME`
  - the study's own name, reported by the loader on import
  - `calculations.config`, the sensitivity stamps, the container gate
* - `R`, `F`
  - gas constant; Faraday constant, used to write trap energies in eV
  - the file itself
* - `MODEL_OPTIONS`
  - how the model is run — five switches, not data about anything
  - `MICROSTRUCTURE`, `OXIDE_DEFECTS` and `DEFAULT_PARAMS_LEVEL5`, by reference
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
* - `SUGGESTED_RANGES_LEVEL5`, `SUGGESTED_RANGES_LEVEL5L6`
  - the interval each parameter is sampled over
  - `calculations.sensitivity`
* - `REGIME_PRESETS`, `REGIME_PRESETS_L5`
  - ranges overridden to target one rate-limiting regime
  - `calculations.sensitivity`
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

The remaining fields are metadata: `reference`, `temperature_range`,
`pressure_range`, `pressure`, `gas`, `Notes`, and `sample_thickness`.

`temperature_range` is the window the source data were measured over.
`get_metal_properties_at_T` checks it and warns on extrapolation — as does
`get_oxide_properties_at_T` for the oxide's. Both materials spell it the same way
for that reason: the metal's used to be `temp_range`, which nothing read, so a
metal extrapolation was silent. `pressure_range` is recorded but read by nothing.

`sample_thickness` is provenance — the specimens the measurements were taken on,
as a *list*, e.g. `[5e-3]` — and is not the wall. The wall comes from
`CONDITIONS['L_metal']`; see {doc}`caveats`.

### `OXIDES`

Same shape as `METALS`, one entry per oxide, with `D_ox_ref` / `E_D_ox` /
`K_ox_ref` / `H_sol_ox` in place of the metal quartet, and the same
`temperature_range` validity window.

Its `sample_thickness` is a scalar rather than a list — the single coupon the
oxide data came from, e.g. `4.8e-08` — and, as for the metal, it is provenance and
not the layer being simulated. The oxide layer the model builds is
`CONDITIONS['L_oxide']`, which `build_simulation_config` injects into the assembled
oxide dict as `thickness`. There is exactly one authored oxide thickness; an
earlier revision carried two, differing by more than an order of magnitude.

As with the metal, `D_ox_0`, `K_ox_0`, `Phi_ox_ref` and `Q_p_ox_J_per_mol` are
recorded but never read.

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
block. The entry here is a **reference** to `MODEL_OPTIONS`, not a second value —
the solvers read a defect dict and look inside it, so the key has to be present,
but it is authored in one place only.

### `MODEL_OPTIONS`

The five switches that say how the model is run, rather than describing anything
physical: `include_gb_trapping`, `use_sieverts_pinhole`, `include_gb_enhancement`,
`include_trapping` and `D_eff_method`.

They were previously scattered — one inside `MICROSTRUCTURE`, one inside
`OXIDE_DEFECTS`, three inside `DEFAULT_PARAMS_LEVEL5` — so answering "what does
this study switch on?" meant knowing all three places. The consuming dicts still
carry the keys their solvers look up, but each is a reference to this block.

### `CONDITIONS`

**The authority for the operating point.** Nothing below it may restate a value it
owns; where a material entry and `CONDITIONS` disagree about the wall,
`CONDITIONS` wins by construction, because the derived blocks read it rather than
retyping it.

It holds the operating point (`T_operating`, `P_upstream`, `P_downstream`,
`L_metal`, `L_oxide`), the sweep ranges the notebooks plot over (`T_range`,
`P_range`, `L_metal_range`, `L_oxide_range`, `defect_fraction_range`,
`k_diss_sweep_range`), the comparison and limit values the level notebooks need
(`defect_fractions_compare`, `defect_limit_check`, `f_gb_max`), and a **single**
`n_points` count shared by every sweep. Single values and sweep bounds sit in one
block, distinguished only by the `_range` suffix.

Sweep bounds live here, not with the materials, because a sweep is part of the
question being asked rather than a property of the substance. Grids are built at
the point of use from a range and `n_points`:

```python
C = CONDITIONS
temps     = np.linspace(*C['T_range'], C['n_points'])
pressures = np.logspace(*np.log10(C['P_range']), C['n_points'])
```

### The active-material shorthand

*Defined in `calculations/config/model_config.py`, not in the study file — shared by every study. Imported from `model_config` exactly as before.*

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

*Defined in `calculations/config/model_config.py`, not in the study file — shared by every study. Imported from `model_config` exactly as before.*

The flattening step. The nested, physically-organised blocks above become one
flat `name → scalar` dict, because that is what the model wrappers and the
sensitivity sampler consume. Every value is a *reference* to a block above —
`_SA_M['D_ref']`, `_SA_DC['pinholes']`, `CONDITIONS['L_metal']` — so the flat
view cannot disagree with the structured one.

`DEFAULT_PARAMS_LEVEL5L6` starts with `**DEFAULT_PARAMS_LEVEL5` and adds the ten
surface-kinetics entries, which is what keeps the two models sharing one value
for everything they have in common. The split mirrors the model fork described
in {doc}`theory/two-models`.

A few entries are not physics: the flags `include_gb_enhancement`,
`include_trapping` and `D_eff_method` say how to run the model. They are read from
`MODEL_OPTIONS` like everything else, so even these are references rather than
literals.

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

*Defined in `calculations/config/model_config.py`, not in the study file — shared by every study. Imported from `model_config` exactly as before.*

Presentation only, but part of the export contract, so a study file must define
all three. `CURVE_STYLES` gives every level and mode a combination of colour,
linestyle and marker that is distinguishable by any one of the three alone;
`apply_style` applies one by name.

### Helper functions

*Defined in `calculations/config/model_config.py`, not in the study file — shared by every study. Imported from `model_config` exactly as before.*

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

## The export contract

A study file is valid only if it defines all thirteen names in `REQUIRED_EXPORTS`:

```text
CONDITIONS   METALS   MICROSTRUCTURE   MODEL_OPTIONS   OXIDES   OXIDE_DEFECTS   R
SUGGESTED_RANGES_LEVEL5   SUGGESTED_RANGES_LEVEL5L6
REGIME_PRESETS            REGIME_PRESETS_L5
F   DEFAULT_SWEEP_TEMPERATURES
```

A name is on that list for one of two reasons, and the difference matters:

- **The first seven are read by the derived block in `model_config` itself.** A study
  missing one of these does not merely disappoint a later import — it breaks the
  derivation, so without the check you would get a bare `NameError` raised from
  inside the loader.
- **The rest are imported by something.** Four are pulled at module level by
  `calculations.sensitivity`, so a study missing one would otherwise fail at
  `import calculations.sensitivity`, before any error message could help.

The list is derived from what the codebase actually imports and what the derived
block actually reads, not from what one study happens to define — it is recomputed
by an AST scan rather than maintained by hand. It grew from fifteen names to twenty
when the scan first ran, then fell to thirteen when the derived and presentation
sections moved into `model_config`: eight names stopped being a study's
responsibility, and `MODEL_OPTIONS` became one, since the derivation reads it
although nothing imports it.

The loader checks the set at import and raises a `ValueError` naming both the study
and the gap.

There is no `VALIDATION` block. Earlier revisions carried one — per-level test
grids for the level notebooks — but it was a second copy of the operating point,
level-prefixed, and it drifted from `CONDITIONS`. Its live entries moved into
`CONDITIONS`, its grids are now built at the point of use, and its dead keys were
deleted.

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
