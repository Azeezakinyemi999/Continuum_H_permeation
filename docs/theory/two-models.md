# Two models, not six levels

The chapters that follow are numbered Level 1 to Level 6, which makes them look
like six refinements of one model. They are not. Levels 1–5 and the Level 6
family rest on **different physical hypotheses about where H₂ dissociates**, and
that choice propagates into the oxide's sorption law, the units of its
permeability, and which results may legitimately be compared with which.

Read this page before the level chapters. It says only which model you are in
and what that costs you; the physics is derived in
{doc}`equilibrium-models` and the thermodynamics in {doc}`foundations`.

## The fork

Dissociation is the branch point, because hydrogen crosses the oxide as whatever
species the model leaves it as.

```{list-table}
:header-rows: 1
:widths: 22 39 39

* -
  - **Model 1** — Levels 1–5
  - **Model 2** — the Level 6 family
* - dissociation happens at
  - the oxide/metal interface
  - the gas/oxide surface, at a finite rate
* - the oxide therefore carries
  - intact H₂ molecules
  - atomic H
* - so the oxide obeys
  - Henry's law, $C = K_{ox}P$
  - Sieverts' law, $C = K_{ox}\sqrt{P}$
* - the metal obeys
  - Sieverts
  - Sieverts
* - fluxes matched
  - 2, at one interface
  - 3, at two interfaces
* - oxide flux function
  - `molecular_diffusion_flux`
  - `oxide_flux`
```

Neither is a mistake. A model with no dissociation step has no atomic hydrogen to
transport, and one with a dissociation step has no intact molecules left. See
{doc}`equilibrium-models` for the full argument and for the unresolved `K_ox`
calibration that follows from it.

## Which model each level belongs to

| Level | Wall | Model |
|---|---|---|
| 1 | pristine metal | 1 (no oxide, so the fork does not bite) |
| 2a | pristine oxide alone | 1 |
| 2b | pristine oxide + pristine metal | 1 |
| 3 | defective oxide + metal | 1 |
| 4 | defective metal, bare | 1 (no oxide) |
| 5 | defective oxide + defective metal | 1 |
| 1+L6 | pristine metal + surface | 2 (no oxide) |
| 2a+L6 | pristine oxide + surface | 2 |
| 2b+L6 | pristine bilayer + surface | 2 |
| 3+L6 | defective oxide + surface | 2 |
| 4+L6 | defective metal + oxide + surface | 2 |
| 5+L6 | full system + surface | 2 |

Levels 1 and 4 model bare metal with no oxide at all, so they sit in Model 1 by
numbering only — nothing in them depends on the oxide's sorption law.

## What may be compared with what

This is the practical consequence, and the easiest thing to get wrong.

```{list-table}
:header-rows: 1
:widths: 26 18 56

* - Comparison
  - Valid?
  - Why
* - within Model 1 (L1…L5)
  - yes
  - one set of laws throughout
* - within Model 2 (L1L6…L5L6)
  - yes
  - one set of laws throughout
* - **L1 vs L1+L6**
  - **yes**
  - neither has an oxide, so the only difference is the surface term. This is the
    clean way to isolate the cost of finite dissociation.
* - L5 vs L5+L6
  - **no**
  - the oxide's sorption law changes as well as the surface term. Any difference
    mixes the two and neither flux nor permeability separates them.
* - L2a vs L2a+L6
  - **no**
  - same problem, with the oxide alone
* - L4 vs L4+L6
  - **no**
  - L4 is bare metal; L4+L6 adds a pristine oxide *and* the surface
```

```{warning}
The 15% flux difference between Level 5 and Level 5L6, and the factor-of-284 jump
in `frac_oxide` between them, are **not** measurements of what surface kinetics
cost. They are a change of surface model *and* a change of oxide sorption law,
reported together. {doc}`level6-surface-kinetics` works through the `frac_oxide`
jump; the size of it is inflated by the `K_ox` calibration issue recorded in
{doc}`equilibrium-models`.

If you want the cost of finite dissociation on its own, compare L1 against L1+L6.
```

## Units, in one place

The two models give the oxide's constants different dimensions. Nothing may be
combined across that boundary — it is what makes a single "system permeability"
impossible in Model 1, as {doc}`permeability` sets out.

| Quantity | Symbol | Model 1 | Model 2 |
|---|---|---|---|
| metal solubility | $K_s$ | mol m⁻³ Pa⁻⁰·⁵ | mol m⁻³ Pa⁻⁰·⁵ |
| oxide solubility | $K_{ox}$ | mol m⁻³ **Pa⁻¹** | mol m⁻³ **Pa⁻⁰·⁵** |
| metal permeability | $\Phi_m = D_mK_s$ | mol m⁻¹ s⁻¹ Pa⁻⁰·⁵ | mol m⁻¹ s⁻¹ Pa⁻⁰·⁵ |
| oxide permeability | $\Phi_{ox} = D_{ox}K_{ox}$ | mol m⁻¹ s⁻¹ **Pa⁻¹** | mol m⁻¹ s⁻¹ **Pa⁻⁰·⁵** |
| oxide permeance | $\alpha = \Phi_{ox}/L_{ox}$ | mol m⁻² s⁻¹ **Pa⁻¹** | mol m⁻² s⁻¹ **Pa⁻⁰·⁵** |
| metal permeance | $\beta = \Phi_m/L_m$ | mol m⁻² s⁻¹ Pa⁻⁰·⁵ | mol m⁻² s⁻¹ Pa⁻⁰·⁵ |
| dissociation rate | $k_{\text{diss}}$ | — | mol m⁻² s⁻¹ Pa⁻¹ |
| surface equilibrium | $K_{eq}$ | — | Pa⁻¹ |
| flux | $J$ | mol m⁻² s⁻¹ | mol m⁻² s⁻¹ |
| permeance | $\Pi = J/\Delta\sqrt{P}$ | mol m⁻² s⁻¹ Pa⁻⁰·⁵ | mol m⁻² s⁻¹ Pa⁻⁰·⁵ |

```{important}
$\Phi_{ox}$ and $\Phi_m$ are commensurable in Model 2 and **not** in Model 1.
Any expression that adds, averages or harmonically combines them is therefore
valid at most in Model 2 — and only with thickness weighting. This is exactly the
defect that made the old `permeability` column unreadable; see {doc}`permeability`.

The configuration declares `K_ox_ref` in mol m⁻³ Pa⁻⁰·⁵, which is right for
Model 2 and wrong for Model 1. That is a known, recorded issue with a stated fix
— see the warning in {doc}`equilibrium-models`.
```

## Where to go next

- {doc}`foundations` — chemical potential, coverage, and why a $\sqrt{P}$
  driving force is legitimate
- {doc}`equilibrium-models` — Sieverts, Henry and Langmuir, where each is
  enforced, and the `K_ox` calibration problem
- {doc}`permeability` — what a permeability means in each model, and the four
  tiers that follow
