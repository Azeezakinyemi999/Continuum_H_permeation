# Permeability — two models, four tiers

```{note}
Every number in this chapter comes from `docs/_tools/worked_values.py`, under
`[permeability tiers]`. Regenerate it after any model change and the diff shows
up here rather than quietly rotting in the prose.
```

## The mechanism decides everything: where does H₂ dissociate?

Hydrogen changes chemical identity on its way through a coated wall, and the
permeability question is settled entirely by **where** that change happens. This
codebase contains two deliberate answers, and they are different physical
hypotheses, not a discrepancy.

### Model 1 — dissociation at the oxide/metal interface

H₂ dissolves in Cr₂O₃ **without breaking**. The molecule diffuses through the
oxide intact, and only on reaching the metal does it dissociate, `H₂ → 2H`, to
enter the lattice as atomic hydrogen.

$$C_{ox} = K_{ox}P \quad\text{(Henry)} \qquad\qquad C_m = K_s\sqrt{P} \quad\text{(Sieverts)}$$

Two fluxes, matched at one internal interface. This is the L1–L5 family:
`oxide_permeation.py` → `interface_solver.py` → `parallel_oxide_defect_paths.py`.

### Model 2 — dissociation at the gas/oxide surface

Dissociative adsorption happens on the **outer** surface, at a finite rate,
described by a Langmuir–Hinshelwood coverage $\theta$. Everything downstream of
that surface therefore carries **atomic** hydrogen — including the oxide. So the
oxide obeys Sieverts too.

$$J_{\text{surf}} = k_{\text{diss}}P(1-\theta)^2 - k_{\text{rec}}\theta^2$$

Three fluxes — surface, oxide, metal — matched at two interfaces. This is the L6
family, `surface_kinetics.py`.

The sorption law of the oxide is thus a *consequence* of where dissociation is
placed, not an independent modelling choice. Molecular transport through the
oxide implies Henry; atomic transport implies Sieverts.

## Why the sorption law decides whether permeability is intrinsic

Permeability is a constant only when flux is proportional to a fixed function of
pressure. Writing the driving force as $\Delta(P^{n})$,

$$\Phi = \frac{J\,L}{\Delta(P^{n})}$$

is constant only if $n$ is fixed. A Henry layer has $n=1$; a Sieverts layer has
$n=1/2$; a surface rate contributes $n=1$ through $\theta$.

**Model 1 mixes exponents across the two layers.** **Model 2 shares one exponent
across both layers** and confines the other to the surface. That single
difference decides everything below.

## Four tiers

Adding a thickness is a real demotion — $\Phi = DK$ is intrinsic, and anything
containing $L_{ox}/L_m$ is a property of one particular stack. But it is not the
same demotion as needing a pressure. Hence four tiers, not two:

| tier | depends on | meaning |
|---|---|---|
| **1 — intrinsic** | nothing | material constant |
| **2 — effective layer** | nothing | microstructure homogenised; sorption law preserved |
| **3 — stack** | $L$ | geometry-dependent, still pressure-independent |
| **4 — apparent** | $L$ and $P$ | property of the wall *and* the experiment |

Tiers 1–3 are reproducible from parameters alone. Tier 4 must be quoted with its
operating point.

## Model 1: no stack permeability exists

### Tier 1 — per layer, and irreducibly incommensurable

| | formula | units | value at 873 K |
|---|---|---|---|
| oxide | $D_{ox}K_{ox}$ | mol/m/s/**Pa** | 4.0017e-15 |
| metal | $D_mK_s$ | mol/m/s/**Pa^0.5** | 1.3685e-11 |

These cannot be combined arithmetically. Not by a harmonic mean, not by a
thickness-weighted harmonic mean, not at all — the units differ. **There is no
$\Phi_{\text{transport}}$ for Model 1.** This is not a defect to be fixed; it is
what a molecular-transport barrier on an atomic-transport substrate costs you.

### Tier 4 — the wall, and it genuinely drifts with pressure

Sweeping $P_{up}$ over five decades on a pristine bilayer:

| control | $J L/\Delta\sqrt{P}$ | $J L/\Delta P$ |
|---|---|---|
| metal-limited | 1.16e-11 → 1.37e-11 (18%) | varies 3 decades |
| oxide-limited | 4.0e-15 → 1.26e-12 (**2.5 decades**) | **4.0e-16, constant** |

Neither normalisation works for both, because the wall's effective $n$ slides
between 1/2 and 1 as control passes between the layers.

**A physical corollary worth stating on its own.** Because the oxide is linear in
$P$ and the metal goes as $\sqrt{P}$, **the rate-limiting layer is
pressure-dependent**. The oxide carries ~16% of the resistance at 100 Pa and
essentially none at 10 MPa, on identical material. An oxide that is an excellent
barrier at low pressure can be irrelevant at high pressure.

### The bilayer is still analytic

Tier 4 does not mean numerical-only. With $\alpha = \Phi_{ox}/L_{ox}$ and
$\beta = \Phi_m/L_m$, matching a Henry layer to a Sieverts layer gives a
quadratic in $\sqrt{P_{int}}$:

$$\alpha\left(P_{up}-P_{int}\right) = \beta\left(\sqrt{P_{int}}-\sqrt{P_{down}}\right)$$

$$\sqrt{P_{int}} = \frac{-\beta+\sqrt{\beta^{2}+4\alpha\left(\alpha P_{up}+\beta\sqrt{P_{down}}\right)}}{2\alpha},
\qquad J = \beta\left(\sqrt{P_{int}}-\sqrt{P_{down}}\right)$$

Verified against the `brentq` solver to below 1e-9 relative across oxide
thickness, oxide diffusivity and pressure. It replaces a root-find with an
expression, and gives the interface solver a free regression test.

## Model 2: the stack permeability exists and is exact

Both layers share $n=1/2$, so $\alpha$ and $\beta$ are commensurable permeances
and the interface pressure mixes linearly:

$$\sqrt{P_{int}} = \frac{\alpha\,g(\theta)+\beta\sqrt{P_{down}}}{\alpha+\beta}$$

Substituting into either flux collapses the whole transport stack to a single
permeance:

$$\boxed{\;J = \underbrace{\frac{1}{1/\alpha+1/\beta}}_{\Pi_{\text{transport}}}
\Big(g(\theta)-\sqrt{P_{down}}\Big)\;}$$

Checked against `solve_steady_state_flux_direct` over five decades of pressure:
**ratio 1.000000 at every point.**

### Tier 3 — a genuine stack permeability

$$\Phi_{\text{transport}} = \frac{L_{ox}+L_m}{L_{ox}/\Phi_{ox}+L_m/\Phi_m}
= 1.1756\times10^{-11}\ \text{mol/m/s/Pa}^{0.5}$$

Exactly pressure-independent. Note it is **86% of the bare metal's**
$\Phi_m = 1.3685\times10^{-11}$ — the oxide adds only 14% of the resistance.
This wall is decisively metal-controlled, in agreement with the regime
classifier.

### The surface is a driving force, not a resistance

The surface contributes **nothing** to $\Pi_{\text{transport}}$. It acts only by
replacing $\sqrt{P_{up}}$ with $g(\theta)$. That gives an exact decomposition:

$$\Phi_{\text{app}} = \Phi_{\text{transport}}\cdot\eta_{\text{surf}},
\qquad \eta_{\text{surf}} = \frac{g(\theta)-\sqrt{P_{down}}}{\sqrt{P_{up}}-\sqrt{P_{down}}}\in(0,1]$$

Swept on the pristine bilayer, `L2b+L6`:

| $P_{up}$ | $\eta_{\text{surf}}$ | $n = \mathrm{d}\ln J/\mathrm{d}\ln P$ |
|---|---|---|
| 100 Pa | 0.767 | 0.617 |
| 1 bar | 0.949 | 0.494 |
| 100 bar | 0.838 | 0.443 |

$\eta_{\text{surf}}$ is non-monotonic and $n$ leaves the $[0.5,1]$ interval in
*both* directions: it rises toward 1 at low pressure where the surface rate
limits, and falls below 0.5 at high pressure where coverage saturation makes
$g(\theta)$ grow more slowly than $\sqrt{P}$.

This formalises an observation already recorded in
{doc}`level6-surface-kinetics`: *"the surface suppresses pressure, not flux
directly"*, where $P_{up}/P_{\text{virtual}} = 83.99$ and $\sqrt{83.99} = 9.1644$
reproduces the measured flux ratio exactly. $\eta_{\text{surf}}$ is that same
ratio, written as a multiplicative factor on permeability rather than as a
pressure ratio.

So Model 2 pays a different price from Model 1. Its stack permeability is clean;
what makes the *measured* permeability pressure-dependent is the surface alone —
and that dependence is isolated in a single dimensionless factor.

### The pinhole path carries the metal's surface kinetics

One Model 2 subtlety that matters for any defect row: a pinhole exposes **bare
metal** to the gas, so dissociation there happens on the metal surface, not the
oxide's. The code passes `k_diss_metal` and `K_eq_metal` for that path
(`surface_kinetics.calculate_path_flux_L6`), with a
`use_sieverts_pinhole` switch falling back to the instantaneous-dissociation
limit ($\theta = 0$, $\sqrt{P_{int}} = \sqrt{P_{up}}$).

So the Model 2 pinhole branch is not merely $\alpha = \infty$ — it is different
surface chemistry, with its own $\theta$. Any $\eta_{\text{surf}}$ reported for a
wall containing pinholes is a flux-weighted blend of two different surface
kinetics, and should be labelled as such.

## Homogenising oxide defects

Two questions hide here and they have different answers: what is the effective
permeability of the **oxide alone**, and may that value be substituted into the
**coupled wall**?

### The oxide alone: a closed form for cracks and grain boundaries

A free-standing defective oxide is pure parallel paths under a single sorption
law, so conductances add by area and no solver is needed:

$$\frac{\Phi_{ox,\text{eff}}}{\Phi_{ox}}
= \left(1-f_{cr}-f_{gb}\right) + \frac{f_{cr}}{\gamma} + f_{gb}\,\beta$$

with $\gamma = L_{crack}/L_{ox}$ and $\beta = D_{gb}/D_{ox}$. Checked against
`calculate_parallel_path_flux` with the metal made effectively free — **ratio
1.00000 at every point**:

| $f_{cr}$ | $f_{gb}$ | $\Phi_{\text{eff}}/\Phi_{ox}$ | flux |
|---|---|---|---|
| 0 | 0 | 1.000 | 4.0017e-04 |
| 1e-3 | 0 | 1.009 | 4.0377e-04 |
| 0 | 1e-3 | 1.009 | 4.0377e-04 |
| 1e-2 | 1e-2 | 1.180 | 4.7220e-04 |
| 5e-2 | 5e-2 | 1.900 | 7.6032e-04 |

This is a proper tier-2 quantity: pressure-independent, and dependent only on the
*ratio* $\gamma$ rather than on $L_{ox}$ itself. It is also not a new model — it
is the code's own per-path $\alpha$ ($\alpha_{cr} = \alpha_{int}/\gamma$,
$\alpha_{gb} = \beta\,\alpha_{int}$) summed by area instead of solved path by path.

### The oxide alone: pinholes are undefined, not merely awkward

A pinhole in a free-standing oxide is an open hole. There is no oxide there, so
there is no resistance, and the permeance is infinite. What crosses is not
permeation at all but gas flow through an aperture — viscous or Knudsen — which
this model does not contain and should not.

The implementation states this directly: `alpha_defect = np.inf` and
`R_oxide = 0.0` on a pinhole path. Hence:

```{important}
A pinhole's resistance is supplied entirely by whatever sits behind it.
Free-standing, it is infinite. With metal behind it, it is finite and set by the
metal. Its contribution therefore belongs to the **assembly**, never to the oxide.
```

### The coupled wall: pinholes are conditional

Once a metal is behind it the pinhole flux is finite, and the question becomes
whether a single lumped oxide can stand in for the branch structure. The obstacle
is that each parallel branch has **its own interface pressure** — the pinhole
branch sits at $P_{int}\approx P_{up}$, the intact branch far below it — while one
lumped oxide forces a single interface pressure on the whole wall.

That makes the validity condition physical and law-independent: lumping works
whenever the **metal** dominates, because then all branches share nearly the same
interface pressure.

| $R_{ox}/R_m$ | lumped / true |
|---|---|
| 3.4 | **1.005** |
| 3.4e2 | 1.68 |
| 3.4e4 | ~1e2 |
| 3.4e6 | up to 1.07e4 |

At the default operating point — which is metal-limited — lumping a pinhole is
accurate to 0.5%. It degrades only as the oxide takes over, which is exactly
where oxide defects become the interesting lever. So the coupled case ships with
its validity condition attached, separately from the oxide-only closed form above.

### Contrast, not area fraction, is the lever

Both results above are governed by one quantity. For a single defect type,

$$\eta = \frac{J_{\text{total}}}{j_{\text{intact}}}
= 1 + \left(\text{contrast} - 1\right)f,
\qquad \text{contrast} = \frac{j_{\text{defect}}}{j_{\text{intact}}}$$

Cracks give contrast $1/\gamma = 10$ and oxide grain boundaries $\beta = 10$, so
$f = 10^{-3}$ buys 0.9% — negligible. A pinhole backed by metal has contrast
$\sim\!10^{4}$, three decades more leverage per unit area fraction.

This is the explicit condition behind the claim {doc}`level3-defective-oxide`
opens with — that a tiny area fraction dominates transport *conditionally*. The
condition is contrast. It is also why `f_pinhole` carries ~2.2 decades in the
sensitivity analysis while `f_crack` and `f_gb_defect` barely register: they are
not weak because their area fractions are small, they are weak because their
contrast is 10.

## What is actually wrong with the reported metric

```python
Phi_oxide = D_ox * K_ox
Phi_metal = D_eff * K_s_metal
permeability = 1.0 / (1.0/Phi_oxide + 1.0/Phi_metal)
```

**The category error.** One number is asked to be two different things: a
material permeability (per layer, no $L$, intrinsic) and a wall permeability (for
the stack, necessarily $L$-dependent). No single quantity can serve both, and
this one is shaped like the first while being read as the second.

**The thickness weighting.** Permeation resistance is $L/\Phi$, not $1/\Phi$. For
Model 2, where the harmonic-mean *structure* is correct, this is the whole bug:
weighting by thickness turns the expression into exactly
$\Phi_{\text{transport}}$. Unweighted it returns **4.0005e-15** where the correct
value is **1.1756e-11** — understating by a factor of **2,939**, and pinning the
result to $\Phi_{ox}$ so the wall reads oxide-limited when it is metal-limited by
a factor of seven. That is the reported contradiction, quantified.

**The missing defect branch.** `f_pinhole`, `f_crack`, `f_gb_defect` and every
surface-kinetics parameter are absent from those three lines, even though at L5
the defect branch can carry most of the flux.

**The consequence for L6.** A formula assembled from $DK$ products has no slot
for a rate constant, so it cannot see the surface at all — which is why it
returned identical values for L5 and L5L6 on cases whose fluxes differ by 15%.

## Tier assignment

### Model 1 — Henry oxide, Sieverts metal, no surface

| level | wall | tier |
|---|---|---|
| L1 | pristine metal | 1 |
| L2a | pristine oxide | 1 (Pa units) |
| L4 | defective metal, bare | 2 |
| OXeff | defective oxide, cracks/GB | 2 |
| L2b | pristine oxide + pristine metal | 4 |
| L3 | defective oxide + pristine metal | 4 |
| L5 | defective oxide + defective metal | 4 |

No tier 3 exists. The layers are incommensurable.

### Model 2 — surface, Sieverts oxide, Sieverts metal

| level | wall | tier |
|---|---|---|
| L1L6 | pristine metal + surface | 3 ⊗ $\eta_{\text{surf}}$ → 4 |
| L2aL6 | pristine oxide + surface | 3 ⊗ $\eta_{\text{surf}}$ → 4 |
| L2bL6 | pristine bilayer + surface | 3 ⊗ $\eta_{\text{surf}}$ → 4 |
| L3L6 | defective oxide + surface | 4 (branch topology) |
| L4L6 | defective metal + surface | 3 ⊗ $\eta_{\text{surf}}$ → 4 |
| L5L6 | defective oxide + defective metal + surface | 4 |

Every Model 2 wall has a tier-3 transport permeability; the surface then demotes
the *measured* value to tier 4 through one dimensionless factor. Parallel defect
branches break tier 3, because each branch has its own $\theta$ and its own
interface pressure.

## One comparison to make with care

L5 and L5L6 are not the same wall with and without surface kinetics. They place
dissociation in different locations and therefore give the oxide different
sorption laws. Any difference between them mixes the surface effect with a change
of oxide physics, and neither flux nor permeability separates the two. Stating
that explicitly matters wherever the two levels are compared.

## Implementation check

Nothing in this document requires a change to how any flux is computed. The
defective oxide is already implemented consistently across all four levels, in
both the module and the notebooks, as per-path solves that are then area-weighted:

| path | L3 | L5 | L3L6 | L5L6 |
|---|---|---|---|---|
| intact | oxide+metal | oxide + `D_eff` metal | $\alpha_{int}=\Phi_{ox}/L_{ox}$ | + `D_eff` iteration |
| pinhole | bare metal | bare defective metal | $\alpha=\infty$, metal kinetics | same + `D_eff` |
| crack | oxide $\times\gamma$ | $\times\gamma$ + `D_eff` metal | $\alpha_{int}/\gamma$ | same + `D_eff` |
| oxide GB | $D_{ox}\times\beta$ | $\times\beta$ + `D_eff` metal | $\beta\,\alpha_{int}$ | same + `D_eff` |
| combine | $(1-f)j_{int}+f j_{def}$ | same | $f_{int}j_{int}+\sum f_i j_i$ | same |

Note that **every path includes the metal**. No function anywhere computes a bare
defective oxide, which is why the oxide-only closed form above is new rather than
a refactor — and why the pinhole case never had to be excluded: the model was
never in a position to ask for it.

## Where this lives in the code

`calculations.permeability`, with one function per level and the tier recorded on
every result. Three refusals are deliberate and are covered by tests:

- `transport_permeability` raises for Model 1. There is no stack permeability to
  return, and inventing one would reintroduce the original fault.
- `oxide_only_permeability` raises on a pinhole fraction. A free-standing oxide
  with a hole in it has infinite permeance.
- `combine_series` raises when the two operands carry different units, which is
  every Model 1 oxide/metal pair.

`lumped_oxide_permeability` does not refuse — it returns a validity flag and the
measured $R_{ox}/R_m$, because the coupled case is quantitatively fine while the
metal dominates.
