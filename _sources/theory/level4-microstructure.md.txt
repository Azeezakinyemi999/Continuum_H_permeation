# Level 4 — metal microstructure

```{note}
**Model 1** by numbering, though nothing here turns on it: this level has no
oxide, so the Henry/Sieverts fork of {doc}`two-models` does not bite.
```

## Purpose

Levels 1 to 3 treat the metal as a perfect lattice with a single diffusivity.
Real alloys have grain boundaries, which conduct hydrogen *faster* than the
lattice, and traps — vacancies, dislocations, carbides — which hold it and slow
its net advance. Level 4 puts both in, and they pull in opposite directions.

For the active study this is the most consequential level, because the wall is
metal-limited: the layer that Levels 1–3 model most crudely is the one carrying
97% of the resistance.

**What it cannot tell you.** See the closing section — in short, it is a
local-equilibrium model with no concentration dependence in $D_{\text{eff}}$, and it
says nothing about the oxide or the surface.

## At a glance

| | |
|---|---|
| **Inputs** | $D_L$, $T$, `MICROSTRUCTURE` (grain size, GB thickness, trap list), `lattice_concentration`, `lattice_density` |
| **Outputs** | `D_eff`, `overall_factor`, `gb_enhancement`, `trapping`, `dominant_effect`, `regime` |
| **Code** | [`defective_metal`](../api/defective_metal) — `combined_microstructure_model`, `trap_occupancy`, `calculate_effective_diffusivity_trapping`, `calculate_gb_enhanced_diffusivity` |
| **Solver** | None. Closed form. |

## Concepts

### Two effects, applied in sequence

```{figure} ../figures/fig-12-level4-competition.svg
:alt: Grain-boundary enhancement near one across realistic grain sizes, against trapping spanning orders of magnitude.
:width: 100%

**Figure 12.** Look at the vertical scales. Grain-boundary enhancement stays within a
few per cent of 1 for any grain size a structural alloy actually has — 1.0016× at
20 µm, and only 1.66× even at 50 nm — while trapping changes the diffusivity by orders
of magnitude over the same window. For an engineering wall, **Level 4 is a trapping
model with a grain-boundary term that rounds away.** Worth knowing before spending
effort measuring grain size.
```

```{important}
The two effects are often described as competing. At realistic grain sizes they do not
compete: in panel (B) the combined curve and the trapping-only curve are
indistinguishable. Grain-boundary enhancement becomes a real multiplier only in the
ultrafine range, which is a different material problem.
```

$$D_{\text{eff}} = \underbrace{\left[(1-f_{gb})D_L + f_{gb}\alpha D_L\right]}_{\text{grain-boundary enhancement}}
\times \underbrace{\frac{1}{1 + \sum_i N_{T,i}K_i/N_L}}_{\text{trapping reduction}}$$

and the headline output is the ratio

$$\eta = \frac{D_{\text{eff}}}{D_L}$$

which is above 1 when grain boundaries win and below 1 when trapping wins. In
the code this is `overall_factor`, and it factorises exactly:

$$\eta = \underbrace{f_{\text{gb}}^{\text{factor}}}_{\texttt{gb\_enhancement['factor']}}
\times \underbrace{f_{\text{red}}}_{\texttt{trapping['reduction\_factor']}}$$

Verified for the active study: $1.000011492 \times 0.738318563 = 0.738327048$,
matching `overall_factor` to machine precision.

### Grain-boundary enhancement

#### How much boundary is there?

Grain boundaries are treated as a parallel fast path occupying a volume
fraction set by stereology. For equiaxed grains the boundary area per unit
volume is $S_v = 3/d$, so a boundary of thickness $\delta$ occupies

$$f_{gb} = \frac{3\delta}{d}$$

With the configured $\delta = 5\times10^{-10}$ m and $d = 10^{-4}$ m this gives
$f_{gb} = 1.5\times10^{-5}$ — verified exactly against the formula.

#### The parallel-path diffusivity

$$D_{gb\text{-enh}} = (1-f_{gb})D_L + f_{gb}D_{gb}, \qquad D_{gb} = \alpha D_L$$

The enhancement factor $\alpha$ is itself Arrhenius, and it *decreases* with
temperature: boundary diffusion has a lower activation energy than lattice
diffusion, so the two converge as thermal energy makes the lattice route
competitive. For the configured LAGB type at 873 K it comes out at
$D_{gb} = 4.174906\times10^{-10}$ against $D_L = 2.363871\times10^{-10}$, so
$\alpha = 1.766$.

```{important}
**Two small factors multiply, and the result is negligible.** Grain boundaries here
are only 1.77× faster than the lattice, *and* they occupy fifteen parts per million of
the volume. The enhancement factor is
$1 + f_{gb}(\alpha - 1) = 1 + 1.5\times10^{-5} \times 0.766 = 1.0000115$ — a change of
**0.001%**.

To make grain boundaries matter you need $f_{gb}(\alpha-1) \gtrsim 1$, which at
$\alpha = 1.766$ would take $f_{gb} \gtrsim 0.57$: grain boundaries occupying more than
half the volume. That is not a material. **Through this route, on this alloy, grain
boundaries cannot matter at any grain size.**
```

```{note}
$\alpha$ used to be 50, read from a table of austenitic steel, and the enhancement was
correspondingly 0.07% rather than 0.001%. Commit `366b70d` made it come from the
metal's own data instead, and the measured ratio is 1.77. Any note quoting
$\alpha = 50$, an enhancement of 1.000735, or an `overall_factor` of 0.738861227
predates that change. Figure 12 is drawn on the Example Wall, whose generic
grain-boundary correlation is more optimistic than this alloy's own data — which is why
its panel (A) still rises at ultrafine grain sizes while the paragraph above says this
alloy's cannot.
```

### Trapping — the Oriani model

#### Local equilibrium between traps and lattice

Oriani's assumption is that hydrogen redistributes between lattice sites and
trap sites *fast* compared with the time it takes to cross the wall, so the two
populations are in local equilibrium everywhere. Each trap type has an
equilibrium constant set by its binding energy:

$$K_i = \exp\!\left(\frac{E_{b,i}}{RT}\right)$$

Only lattice hydrogen is mobile — trapped atoms contribute concentration but not
flux. The effective diffusivity is therefore the lattice diffusivity scaled by
the mobile fraction:

$$D_{\text{eff}} = D_L\frac{C_{\text{mobile}}}{C_{\text{total}}}
= \frac{D_L}{1 + \sum_i N_{T,i}K_i/N_L}$$

#### Where the trapping actually comes from

The configured trap population, with each type's contribution to
`trapping_term`:

| Trap | $N_T$ [m⁻³] | $E_b$ [J/mol] | $K_i$ | $\theta_i$ | contribution |
|---|---|---|---|---|---|
| vacancies | 1.000e+26 | 41489 | 303.73 | 0.020421 | **0.346176** |
| Carbides | 2.000e+25 | 26051 | 36.21 | 0.002485 | 0.008253 |
| grain_boundaries | 6.000e+14 | 26051 | 36.21 | 0.002485 | 0.000000 |
| dislocations | 8.160e+12 | 19297 | 14.28 | 0.000980 | 0.000000 |

The contributions sum to 0.354429, matching `trapping_term` exactly.

Two things are worth noticing. **Vacancies supply 97.7% of the entire trapping
effect** — the model's trapping behaviour is, for practical purposes, a
one-parameter story in the vacancy density. And **two of the four trap types
contribute exactly nothing**: their densities are 10–12 orders of magnitude below
the vacancy density, so their $N_T K/N_L$ terms vanish at the printed precision
regardless of their binding energies. A sensitivity analysis will correctly find
`trap_dislocation_N_T` and `trap_gb_N_T` to be irrelevant, and that is a
consequence of the configured densities, not of the physics being wrong.

Note also that the grain-boundary trap and the carbide trap share a binding
energy of 26051 J/mol, so their $K$ and $\theta$ are identical and only their
densities separate them.

## Procedure

### Composing the two

```text
D_lattice            2.363871e-10  m²/s
D_gb_enhanced        2.363898e-10  m²/s   (+0.001%)
D_eff                1.745310e-10  m²/s   (−26%)
overall_factor       0.738327048
trapping_term        0.354428901
reduction_factor     0.738318563
theta_total          0.026370842
f_gb                 1.500000e-05
D_gb/D_L             1.766
dominant_trap        vacancies
dominant_effect      balanced
regime               competitive
```

Trapping reduces the diffusivity by 26%; grain boundaries add back 0.001%. The net
is a 26% reduction.

:::{warning}
**`dominant_effect = 'balanced'` does not mean the two effects are comparable.**
It is the fallthrough branch of

```python
if   gb_factor > 2.0 and reduction_factor > 0.5:  'gb_enhancement'
elif gb_factor < 2.0 and reduction_factor < 0.5:  'trapping'
else:                                             'balanced'
```

Here `gb_factor = 1.0000115` fails the first test and `reduction_factor = 0.7383`
fails the second, so the label falls through — even though trapping outweighs
grain-boundary enhancement by a factor of roughly **23,000**. Read `overall_factor`,
`reduction_factor` and `gb_enhancement['factor']` to find out what actually
happened; `dominant_effect` only tells you which threshold pair fired.

A separate warning does fire when the effects genuinely near-cancel, at
$0.8 < \eta < 1.25$. With $\eta = 0.7383$ this study sits just outside it.
:::

### Modes

`combined_microstructure_model` takes a `mode` argument that isolates the
effects, which is the cleanest way to check either one:

| `mode` | $D_{\text{eff}}$ [m²/s] | $\eta$ | `dominant_effect` | `regime` |
|---|---|---|---|---|
| `'both'` | 1.745310e-10 | 0.738327048 | `balanced` | `competitive` |
| `'gb_only'` | 2.363898e-10 | 1.000011492 | `gb_enhancement` | `bulk_dominated` |
| `'trapping_only'` | 1.745290e-10 | 0.738318563 | `trapping` | `weak_trapping` |
| `'none'` | 2.363871e-10 | 1.000000000 | `none` | `perfect_lattice` |

```{important}
Compare the `'both'` and `'trapping_only'` rows: they agree to five significant
figures. Grain-boundary enhancement contributes **+0.001%** here, so at this grain size
Level 4 is a trapping model and nothing else — see Figure 12.

These numbers changed when the grain-boundary enhancement stopped being read from a
generic austenitic-steel table and started coming from the metal's own data. The ratio
$D_{gb}/D_L$ went from a tabulated 50 to the metal's 1.77, which is what collapsed the
grain-boundary contribution from +0.07% to +0.001%. Any note quoting the older figures
predates that change.
```

All four measured. Note `'both'` is the product of the two single-effect factors,
not the product of their diffusivities, so $D_{\text{eff}}(\text{both})$ exceeds
$D_{\text{eff}}(\text{trapping\_only})$ by exactly the grain-boundary factor.

`'none'` is the Level 4 → Level 1 limit and should reproduce `D_lattice`
exactly; it is the first thing to check after touching this module.

## Design choices

### Concentration dependence, and a caveat

`combined_microstructure_model` takes a `lattice_concentration` argument, and the
values above use `10.0` mol/m³ — the same value the container smoke test pins, so
they are directly comparable against it.

```{note}
The Oriani form as implemented has **no $C_L$ dependence** in $D_{\text{eff}}$:
the denominator $1 + \sum N_{T,i}K_i/N_L$ contains only trap densities, binding
energies, the lattice site density and temperature. The occupancies $\theta_i$
do depend on concentration, and they are reported, but they do not feed back
into the diffusivity.

The practical consequence is that a concentration profile computed across the
wall comes out with a *flat* $D_{\text{eff}}$ while $\theta$ varies — which is
correct behaviour in the dilute-trap limit, and is the discrepancy recorded in
`STEADY_STATE.md` between the implementation and a docstring that claims
concentration dependence.
```

Oriani also assumes low occupancy. The code warns above $\theta > 0.9$, where
local equilibrium stops being a good approximation and a kinetic trapping model
such as McNabb–Foster would be needed instead. At $\theta_{\max} = 0.0204$ here
there is plenty of margin.

## Checks and failure modes

### Where the validation evidence lives

The trapping block is the one part of the model validated against independent
experiment without fitting: `TRAPPING_VALIDATION.md` works through a
zero-free-parameter comparison against thermal desorption spectroscopy, using
trap densities and binding energies measured elsewhere. That document is the
reason to trust the numbers in this chapter, and it is worth reading alongside
it.

## What Level 4 cannot tell you

It still assumes a perfect oxide. Combining a defective oxide with a defective
metal is Level 5, and it is not a matter of multiplying the two corrections —
each defect path through the oxide sits above its own patch of metal and must be
solved with the microstructural diffusivity in place.

Dissociation remains infinitely fast, which {doc}`level6-surface-kinetics` addresses.
