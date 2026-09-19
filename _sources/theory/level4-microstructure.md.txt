# Level 4 — metal microstructure

```{note}
**Model 1** by numbering, though nothing here turns on it: this level has no
oxide, so the Henry/Sieverts fork of {doc}`two-models` does not bite.
```

Levels 1 to 3 treat the metal as a perfect lattice with a single diffusivity.
Real alloys have grain boundaries, which conduct hydrogen *faster* than the
lattice, and traps — vacancies, dislocations, carbides — which hold it and slow
its net advance. Level 4 puts both in, and they pull in opposite directions.

For the active study this is the most consequential level, because the wall is
metal-limited: the layer that Levels 1–3 model most crudely is the one carrying
98% of the resistance.

## Two effects, applied in sequence

$$D_{\text{eff}} = \underbrace{\left[(1-f_{gb})D_L + f_{gb}\alpha D_L\right]}_{\text{grain-boundary enhancement}}
\times \underbrace{\frac{1}{1 + \sum_i N_{T,i}K_i/N_L}}_{\text{trapping reduction}}$$

and the headline output is the ratio

$$\eta = \frac{D_{\text{eff}}}{D_L}$$

which is above 1 when grain boundaries win and below 1 when trapping wins. In
the code this is `overall_factor`, and it factorises exactly:

$$\eta = \underbrace{f_{\text{gb}}^{\text{factor}}}_{\texttt{gb\_enhancement['factor']}}
\times \underbrace{f_{\text{red}}}_{\texttt{trapping['reduction\_factor']}}$$

Verified for the active study: $1.000735 \times 0.738318563 = 0.738861227$,
matching `overall_factor` to machine precision.

## Grain-boundary enhancement

### How much boundary is there?

Grain boundaries are treated as a parallel fast path occupying a volume
fraction set by stereology. For equiaxed grains the boundary area per unit
volume is $S_v = 3/d$, so a boundary of thickness $\delta$ occupies

$$f_{gb} = \frac{3\delta}{d}$$

With the configured $\delta = 5\times10^{-10}$ m and $d = 10^{-4}$ m this gives
$f_{gb} = 1.5\times10^{-5}$ — verified exactly against the formula.

### The parallel-path diffusivity

$$D_{gb\text{-enh}} = (1-f_{gb})D_L + f_{gb}D_{gb}, \qquad D_{gb} = \alpha D_L$$

The enhancement factor $\alpha$ is itself Arrhenius, and it *decreases* with
temperature: boundary diffusion has a lower activation energy than lattice
diffusion, so the two converge as thermal energy makes the lattice route
competitive. For the configured LAGB type at 873 K, $\alpha = 50.0$ exactly.

```{important}
**A large $\alpha$ does not imply a large effect.** Grain boundaries are 50×
faster than the lattice here, and they still change the diffusivity by
**0.07%** — because they occupy fifteen parts per million of the volume. The
enhancement factor is $1 + f_{gb}(\alpha - 1) = 1 + 1.5\times10^{-5} \times 49 =
1.000735$.

To make grain boundaries matter you need $f_{gb}\alpha \gtrsim 1$, which for
$\alpha = 50$ means $f_{gb} \gtrsim 0.02$, i.e. a grain size of order 75 nm.
Nanocrystalline material qualifies; the configured 100 µm grain does not, by
three orders of magnitude.
```

## Trapping — the Oriani model

### Local equilibrium between traps and lattice

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

### Where the trapping actually comes from

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

## Composing the two

```text
D_lattice            2.363871e-10  m²/s
D_gb_enhanced        2.365609e-10  m²/s   (+0.07%)
D_eff                1.746573e-10  m²/s   (−26%)
overall_factor       0.738861227
trapping_term        0.354428901
reduction_factor     0.738318563
theta_total          0.026370842
f_gb                 1.500000e-05
alpha (D_gb/D_L)     50.0
dominant_trap        vacancies
dominant_effect      balanced
regime               competitive
```

Trapping reduces the diffusivity by 26%; grain boundaries add back 0.07%. The net
is a 26% reduction.

:::{warning}
**`dominant_effect = 'balanced'` does not mean the two effects are comparable.**
It is the fallthrough branch of

```python
if   gb_factor > 2.0 and reduction_factor > 0.5:  'gb_enhancement'
elif gb_factor < 2.0 and reduction_factor < 0.5:  'trapping'
else:                                             'balanced'
```

Here `gb_factor = 1.0007` fails the first test and `reduction_factor = 0.7383`
fails the second, so the label falls through — even though trapping outweighs
grain-boundary enhancement by a factor of roughly 350. Read `overall_factor`,
`reduction_factor` and `gb_enhancement['factor']` to find out what actually
happened; `dominant_effect` only tells you which threshold pair fired.

A separate warning does fire when the effects genuinely near-cancel, at
$0.8 < \eta < 1.25$. With $\eta = 0.7389$ this study sits just outside it.
:::

## Modes

`combined_microstructure_model` takes a `mode` argument that isolates the
effects, which is the cleanest way to check either one:

| `mode` | $D_{\text{eff}}$ [m²/s] | $\eta$ | `dominant_effect` | `regime` |
|---|---|---|---|---|
| `'both'` | 1.746573e-10 | 0.738861227 | `balanced` | `competitive` |
| `'gb_only'` | 2.365609e-10 | 1.000735000 | `gb_enhancement` | `bulk_dominated` |
| `'trapping_only'` | 1.745290e-10 | 0.738318563 | `trapping` | `weak_trapping` |
| `'none'` | 2.363871e-10 | 1.000000000 | `none` | `perfect_lattice` |

All four measured. Note `'both'` is the product of the two single-effect factors,
not the product of their diffusivities, so $D_{\text{eff}}(\text{both})$ exceeds
$D_{\text{eff}}(\text{trapping\_only})$ by exactly the grain-boundary factor.

`'none'` is the Level 4 → Level 1 limit and should reproduce `D_lattice`
exactly; it is the first thing to check after touching this module.

## Concentration dependence, and a caveat

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

## Where the validation evidence lives

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
