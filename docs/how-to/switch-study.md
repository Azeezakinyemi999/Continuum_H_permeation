# How to switch material study

All material data, operating conditions, microstructure and defect parameters
resolve through a single switch, so exactly one study is active at a time.

## The switch

Edit one line in `calculations/config/model_config.py`:

```python
ACTIVE_STUDY = 'Guo_etal_2025_316L'
```

| Value | System |
|---|---|
| `Guo_etal_2025_316L` | 316L stainless / Cr₂O₃ |
| `fuerst_etal_2024_model_config` | Hastelloy N / Cr₂O₃ |
| `incoloy802_cr2o3` | Incoloy 802 (X40 NiCrAlTi 31/19) / Cr₂O₃ |

Nothing else changes. Every module in `calculations/` imports its data from
`model_config`, which re-exports whichever study is named.

```{important}
**Restart the kernel.** The binding happens at import time — `model_config`
imports the study module and copies its public names into its own namespace, so
`importlib.reload` on a notebook module is not enough. A running kernel will keep
serving the old study's data.
```

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
On the active study this currently reports **four pre-existing problems** —
`D_ref`, `k_diss_metal_ref` and `E_diss_metal` defaults lying outside their own
sweep ranges. They are not caused by switching; they are present in the
configuration as shipped. A default outside its sweep range means the sensitivity
analysis never evaluates the nominal case, so either the default or the range is
wrong. Resolve before trusting a sensitivity result on those parameters.
```

## Do not hardcode material keys

The `METALS` keys are deliberately different in each study:

```text
316L_Guo_2025
Hastelloy_N_fuerst_2024
metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985
```

So a notebook containing `METALS['316L_Guo_2025']` breaks the moment you switch.
Resolve the active names instead:

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
data without touching `ACTIVE_STUDY`:

```python
import importlib
from calculations.oxide_permeation import get_metal_properties_at_T

other = importlib.import_module(
    'calculations.config.studies.fuerst_etal_2024_model_config')
props = get_metal_properties_at_T('Hastelloy_N_fuerst_2024', 900.0,
                                  metals=other.METALS)
```

This is the supported way to compare studies side by side. Note that it reads
material properties only — the solvers still take their defaults from whichever
study is active.

## When to switch

The active 316L study is metal-limited at its default operating point: the 48 nm
Cr₂O₃ layer contributes `frac_oxide = 3.8e-4` of the resistance and `PRF ≈
1.0002`. If you are investigating oxide or defect physics, that operating point
will show you almost nothing.

To make the oxide matter, either increase `L_oxide`, lower the temperature toward
the validated [473, 773] K window where the oxide is less permeable, or switch to
a study whose oxide is a genuine barrier. See {doc}`../theory/equilibrium-models`
for why the oxide's permeability depends so strongly on which level you are
running.
