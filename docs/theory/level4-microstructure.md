# Level 4 — metal microstructure

```{note}
**Model 1** by numbering, though nothing here turns on it: this level has no
oxide, so the Henry/Sieverts fork of {doc}`two-models` does not bite.
```

## Purpose

Levels 1 to 3 treat the metal as a perfect lattice with a single diffusivity.
Real alloys have grain boundaries, which conduct hydrogen *faster* than the
lattice. Level 4 puts that fast path in: it replaces the lattice diffusivity
$D_L$ with a grain-boundary-enhanced metal diffusivity $D_{metal} \ge D_L$, and
then runs Level 1 with it.

For the active study this matters more than any other level in principle, because
the wall is metal-limited: the layer that Levels 1–3 model most crudely is the one
carrying 97% of the resistance. In practice, as the worked values below show, the
correction on this alloy is a thousandth of a per cent.

**What it cannot tell you.** See the closing section — in short, it is a
parallel-path model with no concentration dependence, it has no trapping, and it
says nothing about the oxide or the surface.

## At a glance

| | |
|---|---|
| **Inputs** | $D_L$, $T$, `MICROSTRUCTURE` (`grain_size`, `gb_thickness`, `gb_type`), `MODEL_OPTIONS['enable_gb_enhancement']` |
| **Outputs** | `D_eff`, `enhancement_ratio`, `gb_enhancement_factor`, `f_gb`, `regime` |
| **Code** | [`defective_metal`](../api/defective_metal) — `metal_diffusivity`, `calculate_gb_enhanced_diffusivity`, `gb_enhancement_factor`; the flux is `permeation_calc.calculate_defective_metal_flux` |
| **Solver** | None. Closed form. |

## Concepts

### One effect: grain boundaries as fast paths

```{figure} ../figures/fig-12-level4-grain-boundaries.svg
:alt: Grain-boundary enhancement of the metal diffusivity against grain size, near one for realistic grains.
:width: 100%

**Figure 12.** Level 4 replaces the lattice diffusivity with a grain-boundary-enhanced
one. The enhancement stays within a few per cent of 1 for any grain size a structural
alloy actually has, and only becomes a real multiplier in the ultrafine range. Trapping
is not drawn because it is not in the model: at steady state traps exchange no net
hydrogen, so they do not change the flux.
```

$$D_{metal} = (1-f_{gb})D_L + f_{gb}\alpha D_L$$

and the headline output is the ratio

$$\eta = \frac{D_{metal}}{D_L} = 1 + f_{gb}(\alpha - 1) \;\ge\; 1$$

In the code $D_{metal}$ is `D_eff` and $\eta$ is `enhancement_ratio` (the Level 4
and Level 5 wrappers report it as `modification_factor` or `D_modification`).
Because it depends on temperature, grain size and boundary type but not on the
hydrogen concentration, a wall has **one** metal diffusivity, and every Level 4
or Level 5 solver is its Level 1–3 counterpart run with $D_{metal}$ in place of
$D_L$.

### How much boundary is there?

Grain boundaries are treated as a parallel fast path occupying a volume
fraction set by stereology. For equiaxed grains the boundary area per unit
volume is $S_v = 3/d$, so a boundary of thickness $\delta$ occupies

$$f_{gb} = \frac{3\delta}{d}$$

With the configured $\delta = 5\times10^{-10}$ m and $d = 10^{-4}$ m this gives
$f_{gb} = 1.5\times10^{-5}$ — verified exactly against the formula.

### How fast is a boundary?

$$D_{gb} = \alpha D_L$$

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
$\alpha = 50$ or an enhancement of 1.000735 predates that change. Figure 12 is drawn on
the Example Wall, whose generic grain-boundary correlation is more optimistic than this
alloy's own data — which is why it still rises at ultrafine grain sizes while the
paragraph above says this alloy's cannot.
```

### Why there is no trapping

Grain boundaries, vacancies, dislocations and carbides also *trap* hydrogen, and
earlier versions of the model folded an Oriani trap-reduced diffusivity into this
level. That was removed. At steady state the traps are full to their equilibrium
occupancy and exchange no net hydrogen, so only lattice hydrogen carries flux and
$J = D\,\Delta C_L/L$ whatever the trap density, binding energy or saturation.
Trapping sets the time lag and the hydrogen inventory, which a steady-state
permeation model does not predict. The argument, and what the old double-counting
did to the numbers, is in {doc}`../STEADY_STATE`.

## Procedure

1. Compute $f_{gb} = 3\delta/d$ (capped at 0.5, with a warning above 0.1).
2. Compute $\alpha = D_{gb}/D_L$ for the configured `gb_type` at $T$
   (`gb_enhancement_factor`).
3. Combine them by the parallel-path rule (or the Hart equation, if `model='hart'`)
   into $D_{metal}$ (`calculate_gb_enhanced_diffusivity`).
4. Hand $D_{metal}$ to the Level 1 flux (`calculate_defective_metal_flux`) or, at
   Level 5, to the Level 2b/3 solvers. `metal_diffusivity` does steps 1–3 and returns
   $D_L$ unchanged when `enable_gb_enhancement` is false.

**Worked values** (active study, 873 K):

```text
D_bulk                   2.363871e-10  m²/s
D_eff                    2.363898e-10  m²/s   (+0.001%)
enhancement_ratio        1.000011
gb_enhancement_factor    1.766131
f_gb                     1.500000e-05
regime                   bulk_dominated
```

`regime` is `bulk_dominated` when $f_{gb}\alpha < 0.1$, `gb_dominated` above 10, and
`mixed` between.

## Design choices

### Parallel path, not Hart, by default

The parallel-path rule treats bulk and boundary as independent conductors in
proportion to their volume. The Hart equation accounts for boundary connectivity and
is the better choice once boundaries percolate ($f_{gb} > 0.1$), which is where the
code warns. At structural grain sizes $f_{gb}$ is of order $10^{-5}$ and the two give
the same answer.

### $\alpha$ from the metal's own data

$\alpha$ is computed from the active metal's lattice activation energy, with the
boundary activation energy taken as a fraction of it, rather than from a generic
table. See the note above for what that changed.

### Trapping removed, not switched off

There is no trapping switch to turn back on. Keeping one would invite comparing
"with" and "without" fluxes that differ only by an error. The removal and its reasons
are recorded in {doc}`../STEADY_STATE`.

## Checks and failure modes

| Symptom | Meaning |
|---|---|
| `enable_gb_enhancement=False` does not reproduce `D_lattice` exactly | Level 4 → Level 1 limit broken; the first thing to check after touching this module |
| Level 4 flux differs from Level 1 run with `D_eff` | the wrapper is doing something other than substituting the diffusivity |
| `enhancement_ratio < 1` | impossible by construction; an input is wrong |
| warning: GB volume fraction > 10% | ultrafine grains; consider `model='hart'` |
| warning: enhancement > 100× | check grain size, boundary type and temperature |

## What Level 4 cannot tell you

It still assumes a perfect oxide. Combining a defective oxide with a defective
metal is Level 5, and it is not a matter of multiplying the two corrections —
each defect path through the oxide sits above its own patch of metal and must be
solved with the microstructural diffusivity in place.

It says nothing about how much hydrogen the metal holds, or how long it takes to
get there. Those need trapping and a time axis; see {doc}`../STEADY_STATE`.

Dissociation remains infinitely fast, which {doc}`level6-surface-kinetics` addresses.
