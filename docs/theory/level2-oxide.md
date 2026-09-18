# Level 2 — the oxide layer

```{note}
**Model 1.** Dissociation sits at the oxide/metal interface, so the oxide carries
intact H₂ under Henry's law. The Level 6 family makes a different choice — see
{doc}`two-models` before comparing results across the two.
```

Level 2 comes in two parts. **Level 2a** is the oxide alone, which is mostly a
vehicle for understanding how oxide transport differs from metal transport.
**Level 2b** is the useful one: oxide and metal in series, coupled through an
interface pressure that has to be solved for.

## Level 2a — the perfect oxide

### Henry's law, not Sieverts

Below Level 6 the model contains no dissociation step, so hydrogen is taken to
dissolve and diffuse through the oxide **as intact H₂ molecules**. A molecule
that stays whole has no factor-of-two stoichiometry to account for, so
concentration is *linear* in pressure:

$$C = K_{ox}P$$

This is Henry's law, and it is the single most important difference between the
oxide and the metal. Everything downstream follows from it. (Once Level 6 models
dissociation explicitly the species entering the oxide becomes atomic and this
changes — see {doc}`equilibrium-models`, which also records an unresolved
consequence for the shared `K_ox` constant.)

The bulk treatment is identical to Level 1 — steady state, so $d^2C/dz^2 = 0$, so
a linear profile and a constant gradient. Substituting the Henry boundary
conditions:

$$\boxed{\;J_{ox} = \frac{D_{ox}K_{ox}}{L_{ox}}\left(P_{\text{up}} - P_{\text{down}}\right)\;}$$

### What linearity buys you

**A different diagnostic slope.** With $P_{\text{down}} = 0$, $J \propto P$, so a
log–log plot of flux against pressure has **slope 1.0**, against 0.5 for the
metal. Measuring that slope tells you which layer is in control without needing
to know either material's properties.

**A pressure-independent resistance.** Defining

$$R_{ox} = \frac{L_{ox}}{D_{ox}K_{ox}}$$

gives $J = \Delta P / R_{ox}$ with $R_{ox}$ a genuine constant. The metal's
resistance is not — it depends on the pressure at which you evaluate it, because
$\sqrt{P}$ is nonlinear. That asymmetry matters as soon as the two are put in
series.

### Verified values

For Cr₂O₃ at 873 K, 48 nm thick, 1 bar → vacuum:

```text
D_ox        1.394877e-17  m²/s
K_ox        2.868829e+02  (see the units warning in equilibrium-models)
R_oxide     1.199501e+07  Pa·s·m²/mol
flux        8.336797e-03  mol/m²/s
```

Compare that flux against Level 1's `4.327598e-06` for the same conditions: the
oxide, taken alone, is roughly **1900× more permeable than the metal**. A 48 nm
Cr₂O₃ layer at this temperature is not a barrier. This is the fact that makes the
active study metal-limited at every subsequent level, and it is worth absorbing
early because it explains why coating parameters barely move the answer.

## Level 2b — oxide and metal in series

### The coupling

Two layers in series must pass the same flux at steady state, because hydrogen
cannot accumulate at the interface between them:

$$J_{ox} = J_{\text{metal}}$$

Writing each in terms of the unknown interface pressure $P_{\text{int}}$, and
defining the two permeances

$$\alpha = \frac{D_{ox}K_{ox}}{L_{ox}}, \qquad \beta = \frac{D_mK_{s,m}}{L_m}$$

flux continuity becomes

$$\alpha\left(P_{\text{up}} - P_{\text{int}}\right)
= \beta\left(\sqrt{P_{\text{int}}} - \sqrt{P_{\text{down}}}\right)$$

```{important}
**$\alpha$ and $\beta$ do not have the same units**, because one multiplies a
pressure and the other a square-root pressure:

$$[\alpha] = \mathrm{mol\,m^{-2}s^{-1}Pa^{-1}}, \qquad
[\beta] = \mathrm{mol\,m^{-2}s^{-1}Pa^{-0.5}}$$

So you cannot compare them directly, and you cannot add $1/\alpha$ and $1/\beta$
as series resistances the way you would two electrical resistors. The mixed
currency is exactly why this level needs a numerical solve instead of an
algebraic one.
```

### Why there is no closed form

The equation is linear in $P_{\text{int}}$ on the left and square-root on the
right. Substituting $u = \sqrt{P_{\text{int}}}$ turns it into a quadratic in $u$,
which *is* solvable for this particular pair of laws — but the same structure
recurs at Levels 3, 4 and 6 with additional nonlinearities that are not
reducible, so the project solves all of them the same numerical way rather than
special-casing this one.

### Solved as a root-finding problem

Define the residual

$$f(P_{\text{int}}) = J_{ox}(P_{\text{int}}) - J_{\text{metal}}(P_{\text{int}})$$

and find its zero. The solution is guaranteed to exist and be unique inside
$[P_{\text{down}}, P_{\text{up}}]$ by a sign argument:

| At | $J_{ox}$ | $J_{\text{metal}}$ | $f$ |
|---|---|---|---|
| $P_{\text{int}} \to P_{\text{down}}$ | large — full ΔP across the oxide | ≈ 0 — no driving force | $> 0$ |
| $P_{\text{int}} \to P_{\text{up}}$ | ≈ 0 — no driving force | large — full ΔP across the metal | $< 0$ |

$f$ is continuous and monotonically decreasing between those, so it crosses zero
exactly once. `solve_interface_pressure` brackets with

```python
P_min = max(P_downstream + P_upstream * 1e-10, min_pressure)
P_max = P_upstream * (1 - 1e-10)
```

and calls `scipy.optimize.brentq`, which combines bisection's guaranteed
convergence with faster interpolation when the function is well behaved.

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

### Verified behaviour at the default point

```text
alpha (oxide permeance)   8.336797e-08  mol/m²/s/Pa
beta  (metal permeance)   1.368507e-08  mol/m²/s/Pa^0.5
converged                 True
P_interface               9.994810e+04  Pa   — 99.948% of P_up
flux                      4.326475e-06  mol/m²/s
flux_error                2.608e-11
resistance_ratio          2.596154e-04
regime                    metal_limited
```

The interface pressure sits at 99.95% of upstream, meaning the oxide consumes
0.05% of the available driving force and the metal takes the rest. The composite
flux `4.326475e-06` is within 0.03% of the bare-metal Level 1 value — the coating
is very nearly invisible.

Three checks, all run against the current code:

| Check | Result |
|---|---|
| flux continuity, $\lvert J_{ox} - J_{\text{metal}}\rvert$ | `2.608e-11` |
| $P_{\text{int}}$ increases monotonically with $P_{\text{up}}$ | holds over $10^3$–$10^6$ Pa, `converged` throughout |
| log–log slope of flux vs $P_{\text{up}}$ | `0.5002`, `0.5001` |

That last row is the most informative line in this chapter. The composite system
reproduces the **metal's** slope of 0.5, not the oxide's 1.0, which is an
independent confirmation that the metal is in control — arrived at without
consulting a single resistance value.

### Reading the regime

`compare_resistances` forms the ratio $R_{ox}/R_{\text{metal}}$ and thresholds it:

| Ratio | `regime` |
|---|---|
| $> 10$ | `oxide_limited` |
| $< 0.5$ | `metal_limited` |
| otherwise | `transition` |

```{note}
The two resistances are not the same kind of object. $R_{ox} = L_{ox}/(D_{ox}K_{ox})$
is a true constant, whereas the metal's is a **linearised differential**
resistance evaluated at the interface pressure,

$$R_{\text{metal}} = \frac{2L_m\sqrt{P_{\text{int}}}}{D_mK_{s,m}}$$

which is the derivative of the $\sqrt{P}$ flux law at that operating point. So
the ratio is a local diagnostic, valid near $P_{\text{int}}$, not a global
material comparison. Verified at the default point: $R_{ox} = 1.199501\mathrm{e}{+7}$,
$R_{\text{metal}} = 4.620302\mathrm{e}{+10}$, ratio $2.596154\mathrm{e}{-4}$ →
`metal_limited`.
```

### The limiting cases

- $\alpha \gg \beta$ — the oxide is transparent, $P_{\text{int}} \to P_{\text{up}}$,
  and the result collapses to Level 1. This is the active study's situation.
- $\alpha \ll \beta$ — the oxide dominates, $P_{\text{int}} \to P_{\text{down}}$,
  and the result collapses to Level 2a.

Both limits are worth testing after any change to the solver, because they are
the two cases where a wrong answer still looks physically plausible.

## What Level 2 cannot tell you

It assumes the oxide is **perfect** — continuous, uniform, unbroken. Real oxides
crack, contain pinholes, and have grain boundaries, and because the intact oxide
is such a good barrier in most systems, those defects tend to dominate transport
even at tiny area fractions. That is Level 3.

It also assumes the metal is a **perfect lattice**, with no traps and no
grain-boundary fast paths. That is Level 4.

Finally, it assumes dissociation at the gas-facing surface is infinitely fast, so
that the surface imposes no resistance at all. That is Level 6.
