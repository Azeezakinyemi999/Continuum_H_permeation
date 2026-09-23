# Level 3 — the defective oxide

```{note}
**Model 1.** Dissociation sits at the oxide/metal interface, so the oxide carries
intact H₂ under Henry's law. The Level 6 family makes a different choice — see
{doc}`two-models` before comparing results across the two.
```

Real oxide layers are not continuous. They contain pinholes where the coating
never formed or spalled off, cracks from thermal cycling and growth stress, and
grain boundaries that conduct faster than the lattice. Level 3 adds those paths.

The standard claim about coating defects is that a tiny area fraction dominates
transport. That claim is **conditional**, and the condition is the subject of this
chapter.

## The parallel-path model

Divide the wall's area into intact and defective regions,
$f_{\text{intact}} + f_{\text{defect}} = 1$. Each region passes its own flux per
unit area, and the total flux per unit wall area is the area-weighted sum:

$$\boxed{\;J_{\text{total}} = j_{\text{intact}}f_{\text{intact}} + j_{\text{defect}}f_{\text{defect}}\;}$$

```{important}
This is an **additive** combination, not a reciprocal one. Parallel electrical
resistors combine as $1/R = \sum 1/R_i$; parallel permeation paths combine as
$J = \sum j_i f_i$. The difference is that flux is extensive — two paths side by
side carry the sum of what each carries — whereas resistance is intensive.

Getting this backwards is a natural mistake for anyone carrying the resistance
analogy over from Level 2, where the layers really are in series.
```

Rearranging for a single defect type makes the leverage explicit:

$$J_{\text{total}} = j_{\text{intact}} + \left(j_{\text{defect}} - j_{\text{intact}}\right)f_{\text{defect}}$$

$$\eta \equiv \frac{J_{\text{total}}}{j_{\text{intact}}}
= 1 + \left(\frac{j_{\text{defect}}}{j_{\text{intact}}} - 1\right)f_{\text{defect}}$$

So the enhancement depends on the *contrast* $j_{\text{defect}}/j_{\text{intact}}$
multiplied by the area fraction. A small fraction matters only when the contrast
is large.

## The three defect types

Each is modelled by re-running the Level 2 series solve with one property
modified, rather than by a separate physical model:

| Type | Physical picture | Implementation | Configured |
|---|---|---|---|
| Pinhole | oxide absent; gas meets bare metal | Level 1 metal-only flux | $f = 0.01$ |
| Crack | oxide locally thinner | Level 2b with $L_{\text{crack}} = \phi L_{ox}$ | $f = 0.005$, $\phi = 0.1$ |
| Grain boundary | faster diffusion path through oxide | Level 2b with $D_{gb} = \psi D_{ox}$ | $f = 0.005$, $\psi = 10$ |

The three fractions sum to the declared total, `area_fraction = 0.02` — verified
to hold exactly for the active study. They are independent populations, so a wall
can carry all three at once, and `calculate_mixed_defect_flux_L6` exists for that
case at Level 6.

## When do defects actually matter?

Here is the conditional claim made concrete. Holding the defect population fixed
at the configured 2% and varying only the intact oxide thickness:

| $L_{ox}$ | $j_{\text{defect}}/j_{\text{intact}}$ | defect share of flux | $\eta$ | PRF | PRF if perfect |
|---|---|---|---|---|---|
| 48 nm (as configured) | 1.000 | 2.000% | 1.000005 | 1.0003 | 1.0003 |
| 48 µm (×10³) | 1.276 | 2.538% | 1.0055 | 1.2856 | 1.2927 |
| 4.8 mm (×10⁵) | 30.79 | 38.59% | 1.5959 | 32.54 | 51.93 |

Read the third column against the fact that the defects occupy **2% of the area
in every row**. At the configured thickness they carry exactly 2% of the flux —
no leverage whatsoever. At 4.8 mm they carry 38.6%, nineteen times their area
share.

```{important}
**Defect leverage is a property of the barrier, not of the defects.** The reason
2% of the area can carry 39% of the flux is that the intact path is bad at its
job; the defect path wins by comparison. When the intact oxide is already more
permeable than the metal underneath it — as the active study's 48 nm layer is, by
a factor of 1926 — there is nothing for a defect to short-circuit, and
$j_{\text{defect}} \approx j_{\text{intact}}$.

So for the active study, Level 3 changes essentially nothing:
$\eta = 1.000005$ and PRF moves from 1.0003 to 1.0003. Any sensitivity analysis
on this study will correctly report the defect parameters as unimportant, and
that is a statement about this operating point, not about defects in general.
```

The last two columns separate two different questions. `PRF_perfect` is what the
coating *could* deliver if it were flawless; `PRF` is what it actually delivers.
At 4.8 mm the defects cost 37% of the achievable protection (32.5 against 51.9) —
which is the regime where defect control is worth engineering effort.

## Permeation Reduction Factor

$$\mathrm{PRF} = \frac{J_{\text{bare metal}}}{J_{\text{coated}}}$$

PRF is the coating-effectiveness metric: how many times less hydrogen crosses
because the coating is there. PRF = 1 means no benefit, PRF = 100 means a
hundredfold reduction. `calculate_PRF` returns both the achieved `PRF` and
`PRF_perfect` for the same oxide without defects, so the gap between them
attributes the shortfall to the defects specifically.

## Verified limit checks

| Check | Expected | Result |
|---|---|---|
| $f_{\text{defect}} \to 0$ | $J_{\text{total}} = j_{\text{intact}}$ | exact — `4.326475e-06` both |
| component fractions sum to the declared total | $0.01 + 0.005 + 0.005 = 0.02$ | holds exactly |
| $\eta \geq 1$ whenever $j_{\text{defect}} \geq j_{\text{intact}}$ | monotone | holds across all three thicknesses |

## Values at the configured operating point

```text
flux_total                  4.326496e-06  mol/m²/s
flux_intact_contribution    4.239945e-06     (area-weighted)
flux_defect_contribution    8.655084e-08     (area-weighted)
flux_intact_per_area        4.326475e-06
flux_defect_per_area        4.327542e-06
defect_enhancement_factor   1.000005
dominant_path               intact_oxide
PRF                         1.000255
regime                      metal_limited/lattice_limited
```

```{note}
**The returned key is `defect_enhancement_factor`, not `enhancement_factor`.**
There is no key by the shorter name, so `result.get('enhancement_factor', ...)`
silently returns your default — which is how a `nan` found its way into a draft
of this chapter.

Note also the two flux pairs. `flux_*_per_area` is the local flux density in each
region; `flux_*_contribution` is that density multiplied by the region's area
fraction. Only the contributions sum to `flux_total`. Comparing a per-area value
against a contribution understates the defect path by a factor of $1/f$.
```

The `regime` string is hierarchical — `metal_limited/lattice_limited`. Because
tier 1 resolved to `metal_limited` rather than `oxide_limited`, the
defect-versus-intact comparison was **never applied**: that test only subdivides
an oxide-limited wall. This is intended behaviour, and
{py:mod}`calculations.classify_regime` explains why, but it means a reader
looking for a defect verdict in `regime` will not find one for this study. The
answer is in `dominant_path`, which reports `intact_oxide`.

## Level 4 metal underneath

`calculate_parallel_path_flux_defective_metal` and
`calculate_PRF_defective_metal` are the Level 3+4 versions, which replace the
perfect-lattice metal beneath each path with the microstructural model of
{doc}`level4-microstructure`. Their `D_eff_metal`, `modification_factor` and
`level4_converged` keys are `None` in the pure Level 3 functions — a signal that
no microstructural calculation ran, not that it failed.

## What Level 3 cannot tell you

The metal beneath every path is still a perfect lattice: no traps, no
grain-boundary enhancement, a diffusivity that is a pure material constant. For
the active study that is the dominant remaining approximation, precisely because
the wall is metal-limited — the layer being modelled most crudely is the one
carrying 98% of the resistance. {doc}`level4-microstructure` addresses it.

Dissociation is also still infinitely fast, which {doc}`level6-surface-kinetics`
addresses.
