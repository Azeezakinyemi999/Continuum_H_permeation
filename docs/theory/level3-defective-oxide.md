# Level 3 — the defective oxide

```{note}
**Model 1.** Dissociation sits at the oxide/metal interface, so the oxide carries
intact H₂ under Henry's law. The Level 6 family makes a different choice — see
{doc}`two-models` before comparing results across the two.
```

## Purpose

Real oxide layers are not continuous. They contain pinholes where the coating never
formed or spalled off, cracks from thermal cycling and growth stress, and grain
boundaries that conduct faster than the lattice. Level 3 adds those paths.

The standard claim about coating defects is that a tiny area fraction dominates
transport. That claim is **conditional**, and the condition is the subject of this
chapter.

**What it cannot tell you.** The metal beneath every path is still a perfect lattice:
no traps, no grain-boundary enhancement, a diffusivity that is a pure material constant.
For the active study that is the dominant remaining approximation, precisely because the
wall is metal-limited — the layer being modelled most crudely is the one carrying 97% of
the resistance. {doc}`level4-microstructure` addresses it. Dissociation is also still
infinitely fast, which {doc}`level6-surface-kinetics` addresses.

## At a glance

| | |
|---|---|
| **Inputs** | Level 2b inputs, plus `area_fraction`, defect `type`, and per-type modifiers |
| **Outputs** | `flux_total`, `flux_*_contribution`, `flux_*_per_area`, `defect_enhancement_factor`, `dominant_path`, `PRF`, `PRF_perfect` |
| **Code** | [`parallel_oxide_defect_paths`](../api/parallel_oxide_defect_paths) |
| **Solver** | One Level 2b solve per path; no new solver |

## Concepts

### The parallel-path model

Divide the wall's area into intact and defective regions,
$f_{\text{intact}} + f_{\text{defect}} = 1$. Each region passes its own flux per unit
area, and the total flux per unit wall area is the area-weighted sum:

$$\boxed{\;J_{\text{total}} = j_{\text{intact}}f_{\text{intact}} + j_{\text{defect}}f_{\text{defect}}\;}$$

```{important}
This is an **additive** combination, not a reciprocal one. Parallel electrical resistors
combine as $1/R = \sum 1/R_i$; parallel permeation paths combine as $J = \sum j_i f_i$.
The difference is that flux is extensive — two paths side by side carry the sum of what
each carries — whereas resistance is intensive.

Getting this backwards is a natural mistake for anyone carrying the resistance analogy
over from Level 2, where the layers really are in series. Figure 3 in {doc}`../overview`
draws the distinction.
```

Rearranging for a single defect type makes the leverage explicit:

$$J_{\text{total}} = j_{\text{intact}} + \left(j_{\text{defect}} - j_{\text{intact}}\right)f_{\text{defect}}$$

$$\eta \equiv \frac{J_{\text{total}}}{j_{\text{intact}}}
= 1 + \left(\frac{j_{\text{defect}}}{j_{\text{intact}}} - 1\right)f_{\text{defect}}$$

So the enhancement depends on the *contrast* $j_{\text{defect}}/j_{\text{intact}}$
multiplied by the area fraction. A small fraction matters only when the contrast is
large.

### The three defect types

Each is modelled by re-running the Level 2 series solve with one property modified,
rather than by a separate physical model:

| Type | Physical picture | Implementation | Configured |
|---|---|---|---|
| Pinhole | oxide absent; gas meets bare metal | Level 1 metal-only flux | $f = 0.01$ |
| Crack | oxide locally thinner | Level 2b with $L_{\text{crack}} = \phi L_{ox}$ | $f = 0.005$, $\phi = 0.1$ |
| Grain boundary | faster diffusion path through oxide | Level 2b with $D_{gb} = \psi D_{ox}$ | $f = 0.005$, $\psi = 10$ |

The three fractions sum to the declared total, `area_fraction = 0.02` — verified to hold
exactly for the active study. They are independent populations, so a wall can carry all
three at once, and `calculate_mixed_defect_flux_L6` exists for that case at Level 6.

### When do defects actually matter?

```{figure} ../figures/fig-11-level3-threshold.svg
:alt: The defect takeover fraction moving three decades between two oxide qualities.
:width: 90%

**Figure 11.** Look at how far the two crossings are apart. The takeover point — where
the defect path carries as much flux as the whole intact area — moves by orders of
magnitude with the quality of the coating. **A threshold quoted without its oxide is
meaningless**, and the better the oxide, the smaller the defect fraction that ruins it.
```

Here is the conditional claim made concrete. Holding the defect population fixed at the
configured 2% and varying only the intact oxide thickness:

| $L_{ox}$ | $j_{\text{defect}}/j_{\text{intact}}$ | defect share of flux | $\eta$ | PRF | PRF if perfect |
|---|---|---|---|---|---|
| 1 µm (as configured) | 1.005 | 2.010% | 1.000103 | 1.0053 | 1.0054 |
| 48 µm (×48) | 1.276 | 2.538% | 1.0055 | 1.2856 | 1.2927 |
| 4.8 mm (×4800) | 30.79 | 38.59% | 1.5959 | 32.54 | 51.93 |

Read the third column against the fact that the defects occupy **2% of the area in every
row**. At the configured thickness they carry 2.01% of the flux — no leverage
whatsoever. At 4.8 mm they carry 38.6%, nineteen times their area share.

```{important}
**Defect leverage is a property of the barrier, not of the defects.** The reason 2% of
the area can carry 39% of the flux is that the intact path is bad at its job; the defect
path wins by comparison. When the intact oxide is already more permeable than the metal
underneath it — as the active study's 1 µm layer is, by a factor of 92 — there is
nothing for a defect to short-circuit, and $j_{\text{defect}} \approx j_{\text{intact}}$.

So for the active study, Level 3 changes essentially nothing: $\eta = 1.000103$, and PRF
moves only from a perfect-coating 1.0054 to an achieved 1.0053. Any sensitivity analysis
on this study will correctly report the defect parameters as unimportant, and **that is
a statement about this operating point, not about defects in general.**
```

The last two columns separate two different questions. `PRF_perfect` is what the coating
*could* deliver if it were flawless; `PRF` is what it actually delivers. At 4.8 mm the
defects cost 37% of the achievable protection (32.5 against 51.9) — which is the regime
where defect control is worth engineering effort.

### Permeation Reduction Factor

$$\mathrm{PRF} = \frac{J_{\text{bare metal}}}{J_{\text{coated}}}$$

PRF is the coating-effectiveness metric: how many times less hydrogen crosses because
the coating is there. PRF = 1 means no benefit, PRF = 100 means a hundredfold reduction.
`calculate_PRF` returns both the achieved `PRF` and `PRF_perfect` for the same oxide
without defects, so the gap between them attributes the shortfall to the defects
specifically.

## Procedure

In the order the code runs it:

1. Read `area_fraction` as $f_{\text{defect}}$, defaulting to 0.01, and set
   $f_{\text{intact}} = 1 - f_{\text{defect}}$.
2. Validate $0 \le f_{\text{defect}} \le 1$ and raise otherwise.
3. **Intact path** — run the full Level 2b series solve, giving $j_{\text{intact}}$.
4. **Defect path** — `calculate_defect_path_flux`, branching on `type`: a pinhole is the
   bare-metal Level 1 flux, a crack is Level 2b with $L_{ox}$ scaled by
   `thickness_factor`, a grain boundary is Level 2b with $D_{ox}$ scaled.
5. Weight each by its area fraction to get the two `*_contribution` values.
6. Sum them into `flux_total`.
7. Compare the contributions to set `dominant_path`.
8. Form `defect_enhancement_factor` $= J_{\text{total}}/j_{\text{intact}}$.
9. Classify hierarchically with `classify_regime_level34`.

## Design choices

**Why defect types are re-runs rather than separate models.** A pinhole, a crack and a
fast grain boundary differ only in which property of the same transport problem is
modified. Writing three physical models would triple the surface area for a bug and
guarantee they drift apart; re-running the Level 2 solve with one property changed keeps
exactly one series solver in the codebase.

**Why the regime label is hierarchical.** `regime` reads
`metal_limited/lattice_limited`: tier 1 resolved to `metal_limited` rather than
`oxide_limited`, so the defect-versus-intact comparison was **never applied** — that
test only subdivides an oxide-limited wall. This is intended, and
{py:mod}`calculations.classify_regime` explains why, but it means a reader looking for a
defect verdict in `regime` will not find one for this study. The answer is in
`dominant_path`, which reports `intact_oxide`.

**Why both `PRF` and `PRF_perfect` are returned.** One number cannot answer both "is
this coating worth having" and "how much is the defect population costing me". Returning
the pair makes the attribution explicit instead of leaving it to be reconstructed.

## Checks and failure modes

| Check | Expected | Where |
|---|---|---|
| $J_{\text{total}}$ equals the area-weighted sum | exact | `check_equations.py` — *L3: J_total = sum* |
| contributions sum to the total | exact | *L3: contributions sum* |
| $f_{\text{defect}} \to 0$ recovers Level 2b | `4.304261e-06` both | *L3 → L2b* |
| component fractions sum to the declared total | $0.01 + 0.005 + 0.005 = 0.02$ | `worked_values.py` |
| $\eta \ge 1$ whenever $j_{\text{defect}} \ge j_{\text{intact}}$ | monotone | holds across all three thicknesses |

```{important}
**The returned key is `defect_enhancement_factor`, not `enhancement_factor`.** There is
no key by the shorter name, so `result.get('enhancement_factor', ...)` silently returns
your default — which is how a `nan` found its way into a draft of this chapter.

Note also the two flux pairs. `flux_*_per_area` is the local flux density in each
region; `flux_*_contribution` is that density multiplied by the region's area fraction.
Only the contributions sum to `flux_total`. Comparing a per-area value against a
contribution understates the defect path by a factor of $1/f$.
```

| Symptom | Meaning |
|---|---|
| `ValueError`: area fraction must be 0–1 | A fraction outside $[0,1]$ reached the solver |
| `defect_enhancement_factor` is `nan` | Almost certainly the key-name trap above |
| No defect verdict in `regime` | Expected when tier 1 is not `oxide_limited`; read `dominant_path` |
| `D_eff_metal`, `modification_factor`, `level4_converged` are `None` | A signal that no microstructural calculation ran, not that it failed — these are the Level 3+4 keys |
| Defects report as unimportant in a sensitivity run | Check the oxide first. On a non-limiting coating this is correct and uninformative. |

```{note}
**Worked values.** Active study, 1 µm coating, defects at the configured 2%.
Regenerated by `docs/_tools/worked_values.py`.
```

```text
flux_total                  4.304704e-06  mol/m²/s
flux_intact_contribution    4.218176e-06     (area-weighted)
flux_defect_contribution    8.652857e-08     (area-weighted)
flux_intact_per_area        4.304261e-06
flux_defect_per_area        4.326428e-06
defect_enhancement_factor   1.000103
dominant_path               intact_oxide
PRF                         1.005318
PRF_perfect                 1.005422
regime                      metal_limited/lattice_limited
```

## Level 4 metal underneath

`calculate_parallel_path_flux_defective_metal` and `calculate_PRF_defective_metal` are
the Level 3+4 versions, which replace the perfect-lattice metal beneath each path with
the microstructural model of {doc}`level4-microstructure`.
