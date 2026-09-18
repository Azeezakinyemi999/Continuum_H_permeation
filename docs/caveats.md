# Warnings and caveats

Every warning in this documentation, collected in one place. Each is
reproduced verbatim from the page it belongs to, and each section links back
to that page for the surrounding argument.

```{note}
These warnings were written to sit alongside the material they qualify. Read
them here for an overview of what can go wrong, but follow the link back
before acting on one — several are procedural and depend on where you are in
a workflow.
```

## Getting started

Source: {doc}`getting-started`

### The active study is metal-limited, so the oxide barely matters

```{warning}
**The active study is metal-limited, so the oxide barely matters.** With a 48 nm
Cr₂O₃ layer on 316L at 873 K, the oxide contributes `frac_oxide = 3.8e-4` of the
resistance and `PRF ≈ 1.0002`. If you are trying to see oxide or defect physics,
this operating point will show you almost nothing — switch to a study or
thickness where the oxide is actually a barrier. See {doc}`how-to/switch-study`.
```

### Oxide properties are extrapolated at the default operating point

```{warning}
**The oxide is extrapolated at the default operating point, and the sensitivity
analysis extrapolates much further — silently.** The Cr₂O₃ data carry
`temperature_range = [473, 773]` K. The active study runs at 873 K, and the
sensitivity sweeps run to **1273 K**, 500 K past the validated ceiling.

Whether you are told depends on which entry point you use:

- `get_oxide_properties_at_T` checks the range and prints
  `Warning: Temperature 873K outside validated range [473, 773]K`. It uses
  `print()`, not `warnings.warn`, so `simplefilter("ignore")` does **not**
  silence it.
- `level5_model_wrapper` and `level5L6_model_wrapper` — the canonical entry
  points, and the ones the sensitivity analysis calls — compute $D_{ox}$ and
  $K_{ox}$ from raw parameters via `arrhenius` and never consult the range.
  They are **silent**.

So the workflow that extrapolates hardest is the one that warns least. Treat
oxide properties above 773 K as an Arrhenius extrapolation you have chosen, not
as validated data. 873 K also sits just below the grain-boundary enhancement
data range of [873.1, 1273.2] K.
```

### `permeability` is apparent and needs its operating point

```{warning}
**`permeability` is an *apparent* permeability and must be quoted with its
operating point.** It is derived from the flux the level actually solved,

$$\Phi_{\text{app}} = \frac{J\,L_{\text{tot}}}{\sqrt{P_{\text{up}}} - \sqrt{P_{\text{down}}}}$$

so it carries the defect paths, the metal microstructure and the surface
kinetics, and it can no longer disagree with `flux`, `PRF` or `regime` — all of
them are now functions of the same flux.

What it is *not* is a material constant. It depends on the layer thicknesses, and
under Model 2 it depends on pressure as well. Report it with
$(T, P_{\text{up}}, P_{\text{down}})$ or not at all.

Do not try to reconstruct it from $D$ and $K$ products. In Model 1 there is no
stack permeability to reconstruct: the Henry oxide and the Sieverts metal are
dimensionally incommensurable, so no weighting combines them. See
{doc}`theory/permeability` for the four tiers of permeability this model
supports, and {doc}`theory/two-models` for why the two halves of the level ladder
differ.
```

### `PRF` is `nan` from the Level 5L6 wrapper

```{warning}
**`PRF` is `nan` from the Level 5L6 wrapper** at default parameters, though it is
finite from Level 5. Use the Level 5 value when you need a coating-effectiveness
number.
```

## References

Source: {doc}`references`

### Quoted parameter values belong to the active study

```{warning}
Parameter values quoted in the theory chapters are the **active study's** real
configured values, regenerated from the running code. They are not illustrative
placeholders. But they are also specific to `Guo_etal_2025_316L` at its default
operating point — switching study changes them, and
`docs/_tools/worked_values.py` is what regenerates them. See
{doc}`how-to/switch-study`.
```

## Equilibrium models

Source: {doc}`theory/equilibrium-models`

### `K_ox` serves both the Henry and Sieverts forms, which cannot share units

```{warning}
**One numerical consequence is unresolved.** Both forms currently read the same
configured constant, `K_ox_ref`, but a Henry constant and a Sieverts constant
cannot share units:

| Form | Required units of $K_{ox}$ |
|---|---|
| Henry, $C = K_{ox}P$ | mol m⁻³ Pa⁻¹ |
| Sieverts, $C = K_{ox}\sqrt{P}$ | mol m⁻³ Pa⁻⁰·⁵ |

The configuration declares `K_ox_ref = 0.35417` as **mol m⁻³ Pa⁻⁰·⁵**, and
`H_sol_ox` is derived as $Q_p - E_D$, which is the Sieverts decomposition
$\Phi = D K_s$. Both point at the Sieverts calibration. The Henry path therefore
uses a Sieverts-calibrated number, which rescales its flux by a factor of
$\sqrt{P}$ — at 1 bar, a factor of 316.

At the default operating point:

| Oxide flux at 873 K, 1 bar → 0 | value [mol/m²/s] |
|---|---|
| Henry form, as the code computes it | 8.336797e-03 |
| Sieverts form, same constant | 2.636327e-05 |

For the active study this does not change any conclusion — the 48 nm oxide is
non-limiting either way, so `regime` is `metal` and `PRF ≈ 1.0002`. It *would*
matter for a thicker or less permeable oxide, and it means a Level 3 defect
result is not directly comparable against a Level 5L6 one.

**The fix is to split the symbol** into two separately calibrated constants —
one Henry, one Sieverts — rather than to change either flux law. That requires a
calibration decision (what is the molecular solubility of H₂ in Cr₂O₃?) which is
not recoverable from the current configuration, so it is recorded here rather
than guessed at.
```

## Level 1 — the perfect metal

Source: {doc}`theory/level1-metal`

### The stored `Phi_ref` is not the permeability the model uses

```{warning}
**The stored `Phi_ref` is not the permeability the model uses.** The
configuration asserts the invariant in its own notes field — `Phi_ref = D_ref *
Ks_ref` — but the values do not satisfy it:

| | stored $\Phi_{\text{ref}}$ | $D_{\text{ref}}K_{\text{ref}}$ | disagreement |
|---|---|---|---|
| 316L metal | 1.0000e-12 | 1.0920e-12 | 9.2% |
| Cr₂O₃ oxide | 3.4000e-19 | 2.7625e-19 | 18.7% |

The activation energies, by contrast, satisfy their derived relations *exactly*:
$Q_p = E_D + \Delta H_s$ (61750 = 52102 + 9648) and
$Q_{p,ox} = E_{D,ox} + \Delta H_{sol,ox}$ (234000). The energies were obtained by
arithmetic; the pre-factors were taken from separately reported and independently
rounded measurements.

No computed result is affected, because `Phi_ref`, `Phi_ox_ref`, `Q_p` and
`Q_p_ox_J_per_mol` are **read by no solver** — the model always works from the
`*_ref` + activation-energy pair. But anyone quoting `Phi_ref` from the
configuration is reporting a permeability that differs from the model's by the
amounts above, so cite $D_{\text{ref}}K_{s,\text{ref}}$ instead.
```

## Level 2 — the oxide layer

Source: {doc}`theory/level2-oxide`

### Always check the `converged` flag

```{warning}
**Always check the `converged` flag.** `solve_interface_pressure` has three
fallback paths that return a plausible-looking result without solving anything:

- a degenerate early return setting $P_{\text{int}} = P_{\text{down}} + \epsilon$
- invalid bounds ($P_{\min} \geq P_{\max}$), which falls back to the geometric
  mean $\sqrt{P_{\text{up}}P_{\text{down}}}$
- $f$ having the same sign at both ends, which sets $P_{\text{int}} = P_{\min}$

Each returns `'converged': False`, and the returned `flux` is computed from that
unsolved interface pressure. Nothing raises. A caller that ignores the flag can
silently use a fabricated answer, so treat `converged` as part of the result, not
as diagnostics.
```

## Level 4 — metal microstructure

Source: {doc}`theory/level4-microstructure`

### Two similar-sounding diffusivities differ by 32%

```{warning}
**Two similar-sounding quantities are not the same, and confusing them is a 32%
error.** The code distinguishes:

| Key | Meaning | Value here |
|---|---|---|
| `trapping_term` | $\sum_i N_{T,i}K_i/N_L$ — the denominator | 0.354429 |
| `theta_total` | $\sum_i \theta_i$, sum of trap *occupancies* | 0.026371 |

`reduction_factor` $= 1/(1 + \texttt{trapping\_term})$, verified exact. The
occupancy sum is carried "for info" and appears in **no** formula.

The project's earlier design note wrote the denominator as $1 + \theta_{\text{total}}$,
using the name the code gives the occupancy sum. Substituting the occupancy sum
into the formula gives $D_{\text{eff}} = 2.304828\times10^{-10}$ instead of the
correct $1.746573\times10^{-10}$ — **32% high**. If you have inherited that
formula, check which quantity you fed it.
```

## Level 6 — surface kinetics

Source: {doc}`theory/level6-surface-kinetics`

### `rate_limiting` and `fraction_surface` do not measure what the surface cost

```{warning}
**`rate_limiting` and `fraction_surface` do not measure what the surface cost
you.** They are computed at the *solved* operating point, so they answer "given
this coverage, where does the remaining driving force drop?" — not "how much flux
was lost relative to an infinitely fast surface?"

The surface acts by **depressing the effective upstream pressure**, and that loss
happens before the series decomposition begins. A surface holding 2.78% of the
resistance can still be responsible for a ninefold flux reduction.

This is why `calculations.sensitivity` documents `flux` as the primary metric —
its comment notes it is *the only metric that responds to surface kinetics*.
`permeability` is now derived from that same flux, so it does move with the
surface — but only as the flux moves, and it remains an apparent value that means
nothing without its operating point ({doc}`theory/permeability`). `frac_surface` sees
only the residual series share. If you want to know whether surface kinetics
matter, compare fluxes with and without them, or inspect $\theta$ against
$\theta_{eq}$.
```

## Two models, not six levels

Source: {doc}`theory/two-models`

### L5-vs-L5L6 differences mix surface kinetics with a change of oxide law

```{warning}
The 15% flux difference between Level 5 and Level 5L6, and the factor-of-284 jump
in `frac_oxide` between them, are **not** measurements of what surface kinetics
cost. They are a change of surface model *and* a change of oxide sorption law,
reported together. {doc}`theory/level6-surface-kinetics` works through the `frac_oxide`
jump; the size of it is inflated by the `K_ox` calibration issue recorded in
{doc}`theory/equilibrium-models`.

If you want the cost of finite dissociation on its own, compare L1 against L1+L6.
```

## How-to: sensitivity analysis

Source: {doc}`how-to/sensitivity-analysis`

### Presets suppress the very parameters they pin

```{warning}
**Presets suppress the very parameters they pin.** A preset narrows or fixes
parameter ranges to steer the sampling. A parameter that has been narrowed cannot
show sensitivity, because it barely varies — so a low score for a
preset-suppressed parameter carries no information at all.

Always check a parameter's range in the preset before concluding it does not
matter. `presets_without(presets, 'name')` removes a parameter from the sampled
set so it can be pinned through `fixed_params` instead, which is the supported
way to take something out of competition deliberately.
```

### The default `n_dummy=3` is too small and inflates every margin

```{warning}
**The default `n_dummy=3` is too small, and it inflates every margin.** With
only three dummies the maximum underestimates the true noise ceiling.

Measured on the L5 clusters, raising `n_dummy` from 3 to 20 lifted the floor by
9–34% and cut the number of parameters scoring above 1.0× the floor **from 21 to
8** in the defect regime. Two thirds of the apparently significant parameters
were noise.

Pass `n_dummy >= 20` for anything that will be published. The shipped default is
convenient for exploration, not defensible for a result.
```

### The sweep range extrapolates past the oxide's validated window

```{warning}
That operating range extrapolates well past the oxide's validated window of
[473, 773] K, and the wrappers used by this pipeline do not warn about it — see
{doc}`getting-started`.
```

### Run the SA notebooks before the parallel-coordinates ones

```{warning}
**Run the SA notebooks before the parallel-coordinates ones.** The SA notebooks
write scan CSVs into `sa_results*/`, which is gitignored; the
parallel-coordinates notebooks only *read* those files and cannot regenerate
them. On a fresh clone the visualisation notebooks will fail until the scans have
been produced.
```

## How-to: switch study

Source: {doc}`how-to/switch-study`

### The verifier reports four pre-existing problems on the active study

```{warning}
On the active study this currently reports **four pre-existing problems** —
`D_ref`, `k_diss_metal_ref` and `E_diss_metal` defaults lying outside their own
sweep ranges. They are not caused by switching; they are present in the
configuration as shipped. A default outside its sweep range means the sensitivity
analysis never evaluates the nominal case, so either the default or the range is
wrong. Resolve before trusting a sensitivity result on those parameters.
```
