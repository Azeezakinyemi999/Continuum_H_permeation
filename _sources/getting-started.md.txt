# Getting started

## What the model computes

Given a wall — a structural alloy with an oxide layer on its upstream face — and
an operating point (temperature, upstream and downstream hydrogen pressure), the
model returns the steady-state hydrogen flux through that wall, together with a
breakdown of **where the resistance to that flux actually sits**.

That breakdown is the point. A flux on its own is hard to act on; knowing that
94% of the resistance is in the metal and 0.04% in the oxide tells you that
improving the coating is wasted effort at this operating point.

## Install

```bash
conda activate mace_env          # or any Python >= 3.9 environment
pip install -e .                 # from the repository root
```

The editable install is what makes `from calculations... import ...` resolve from
any working directory. Dependency bounds are deliberately ranges rather than
pins, so this resolves to a no-op in an environment that already has numpy,
scipy, pandas, matplotlib and SALib. Exact pins for the reproducible container
live in `container/requirements.lock.txt`.

## The model hierarchy

Each level adds one physical effect to the one before it. The levels are
cumulative, not alternatives.

| Level | Adds | Chapter |
|---|---|---|
| 1 | Bulk metal: Fickian diffusion with Sieverts' law boundaries | Level 1 |
| 2a | A perfect oxide layer, treated as molecular (Henry's law) transport | Level 2 |
| 2b | Oxide and metal in series, coupled through an interface pressure | Level 2 |
| 3 | Defects short-circuiting the oxide — pinholes, cracks, grain boundaries | Level 3 |
| 4 | Metal microstructure: grain-boundary fast paths and trapping | Level 4 |
| 5 | The full wall — defective oxide over defective metal | Level 5 |
| 6 | Finite-rate surface dissociation, with Langmuir coverage | Level 6 |

Level 6 composes with the others rather than sitting on top of them, so you will
see combinations written `L1+L6`, `L2a+L6` and `L5L6`. `L5L6` — the full wall
with surface kinetics — is the most complete configuration.

## First run

The one-call entry point is the Level 5 wrapper, which assembles the whole wall
from the active study's configuration:

```python
import warnings
warnings.simplefilter("ignore")

from calculations.config.model_config import ACTIVE_STUDY, DEFAULT_PARAMS_LEVEL5
from calculations.sensitivity import level5_model_wrapper

print(ACTIVE_STUDY)
r = level5_model_wrapper(DEFAULT_PARAMS_LEVEL5)
for k in ("flux", "regime", "PRF", "frac_oxide", "frac_metal", "frac_defect"):
    print(f"{k:<12} {r[k]}")
```

For the active study this prints:

```text
Guo_etal_2025_316L
flux         3.196893...e-06
regime       metal
PRF          1.000188...
frac_oxide   0.000375...
frac_metal   0.979620...
frac_defect  0.020003...
```

Read that as: 3.2 µmol/m²/s of hydrogen crosses the wall; the metal carries 98%
of the resistance; and the coating reduces the flux by a factor of 1.0002, which
is to say **not at all**.

## Read these caveats before trusting a number

Four behaviours will mislead you if you meet them without warning. None is a bug;
all are consequences of how the model is defined.

```{note}
An earlier version of this page described `permeability` as a harmonic mean of
the two bulk permeabilities, disagreeing with `flux` by construction and to be
treated as a bulk material property only. That was an accurate description of the
old metric, which was a material-constant formula applied to a composite wall. It
no longer describes what is computed.
```

## Which material system is active

All material data resolve through a single switch, so exactly one study is active
at a time:

| Study | System | Status |
|---|---|---|
| `Guo_etal_2025_316L` | 316L stainless / Cr₂O₃ | active |
| `fuerst_etal_2024_model_config` | Hastelloy N / Cr₂O₃ | available |
| `incoloy802_cr2o3` | Incoloy 802 / Cr₂O₃ | available |

Changing it is a one-line edit followed by a validation step — see
{doc}`how-to/switch-study`.

## Where to go next

{doc}`theory/foundations` starts the physics from chemical potential and builds
to the driving force the code actually uses. If you would rather see the wall
assembled first and the thermodynamics second, the level chapters stand on their
own.
