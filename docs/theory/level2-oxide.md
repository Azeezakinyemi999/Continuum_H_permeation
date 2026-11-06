# Level 2 — the oxide layer

```{note}
**Model 1.** Dissociation sits at the oxide/metal interface, so the oxide carries
intact H₂ under Henry's law. The Level 6 family makes a different choice — see
{doc}`two-models` before comparing results across the two.
```

## Purpose

Level 2 comes in two parts. **Level 2a** is the oxide alone, which is mostly a vehicle
for understanding how oxide transport differs from metal transport. **Level 2b** is the
useful one: oxide and metal in series, coupled through an interface pressure that has to
be solved for.

**What it cannot tell you.** It assumes the oxide is **perfect** — continuous, uniform,
unbroken. Real oxides crack, contain pinholes and have grain boundaries, and because the
intact oxide is such a good barrier in most systems, those defects tend to dominate
transport even at tiny area fractions. That is Level 3. It also assumes the metal is a
**perfect lattice**, with no traps and no grain-boundary fast paths — Level 4. Finally
it assumes dissociation at the gas-facing surface is infinitely fast, so the surface
imposes no resistance at all — Level 6.

## At a glance

| | |
|---|---|
| **Inputs** | $D_{ox}$, $K_{ox}$, $L_{ox}$ and, for 2b, the Level 1 metal inputs |
| **Outputs** | 2a: `flux`, `R_oxide`. 2b: `flux`, `P_interface`, `resistance_ratio`, `regime`, `converged` |
| **Code** | [`oxide_permeation`](../api/oxide_permeation), [`interface_solver`](../api/interface_solver) |
| **Solver** | 2a closed form; 2b `scipy.optimize.brentq` on a bracketed residual |

## Concepts

### Henry's law, not Sieverts

Below Level 6 the model contains no dissociation step, so hydrogen is taken to dissolve
and diffuse through the oxide **as intact H₂ molecules**. A molecule that stays whole
has no factor-of-two stoichiometry to account for, so concentration is *linear* in
pressure:

$$C = K_{ox}P$$

This is Henry's law, and it is the single most important difference between the oxide
and the metal. Everything downstream follows from it. (Once Level 6 models dissociation
explicitly the species entering the oxide becomes atomic and this changes — see
{doc}`equilibrium-models`, which also records an unresolved consequence for the shared
`K_ox` constant.)

The bulk treatment is identical to Level 1 — steady state, so $d^2C/dz^2 = 0$, so a
linear profile and a constant gradient. Substituting the Henry boundary conditions:

$$\boxed{\;J_{ox} = \frac{D_{ox}K_{ox}}{L_{ox}}\left(P_{\text{up}} - P_{\text{down}}\right)\;}$$

### What linearity buys you

**A different diagnostic slope.** With $P_{\text{down}} = 0$, $J \propto P$, so a
log–log plot of flux against pressure has **slope 1.0**, against 0.5 for the metal.
Measuring that slope tells you which layer is in control without needing to know either
material's properties — see Figure 4 in {doc}`two-models`, and {doc}`permeability` for
the exponent as a general diagnostic.

**A pressure-independent resistance.** Defining

$$R_{ox} = \frac{L_{ox}}{D_{ox}K_{ox}}$$

gives $J = \Delta P / R_{ox}$ with $R_{ox}$ a genuine constant. The metal's resistance
is not — it depends on the pressure at which you evaluate it, because $\sqrt{P}$ is
nonlinear. That asymmetry matters as soon as the two are put in series.

### The coupling

```{figure} ../figures/fig-06-interface-root.svg
:alt: Oxide flux falling and metal flux rising with interface pressure, crossing at the solution.
:width: 90%

**Figure 6.** Look at the single crossing. The oxide flux falls as the interface
pressure rises and the metal flux rises, so they meet at exactly one interface
pressure — and that crossing is the answer. The curve shapes differ because the laws
differ: the oxide line is straight in $P$, the metal curve goes as $\sqrt{P}$. The
shaded band is the bracket `brentq` searches, and the root is guaranteed to lie in it.
```

Two layers in series must pass the same flux at steady state, because hydrogen cannot
accumulate at the interface between them:

$$J_{ox} = J_{\text{metal}}$$

Writing each in terms of the unknown interface pressure $P_{\text{int}}$, and defining
the two permeances

$$\alpha = \frac{D_{ox}K_{ox}}{L_{ox}}, \qquad \beta = \frac{D_mK_{s,m}}{L_m}$$

flux continuity becomes

$$\alpha\left(P_{\text{up}} - P_{\text{int}}\right)
= \beta\left(\sqrt{P_{\text{int}}} - \sqrt{P_{\text{down}}}\right)$$

```{important}
**$\alpha$ and $\beta$ do not have the same units**, because one multiplies a pressure
and the other a square-root pressure:

$$[\alpha] = \mathrm{mol\,m^{-2}s^{-1}Pa^{-1}}, \qquad
[\beta] = \mathrm{mol\,m^{-2}s^{-1}Pa^{-0.5}}$$

So you cannot compare them directly, and you cannot add $1/\alpha$ and $1/\beta$ as
series resistances the way you would two electrical resistors. The mixed currency is
exactly why this level needs a numerical solve instead of an algebraic one, and why
Model 1 has no single system permeability at all — {doc}`permeability`.
```

### Reading the regime

`compare_resistances` forms the ratio $R_{ox}/R_{\text{metal}}$ and thresholds it:

| Ratio | `regime` |
|---|---|
| $> 10$ | `oxide_limited` |
| $< 0.5$ | `metal_limited` |
| otherwise | `transition` |

```{important}
The two resistances are not the same kind of object. $R_{ox} = L_{ox}/(D_{ox}K_{ox})$ is
a true constant, whereas the metal's is a **linearised differential** resistance
evaluated at the interface pressure,

$$R_{\text{metal}} = \frac{2L_m\sqrt{P_{\text{int}}}}{D_mK_{s,m}}$$

which is the derivative of the $\sqrt{P}$ flux law at that operating point. So the ratio
is a local diagnostic, valid near $P_{\text{int}}$, not a global material comparison.
```

### The limiting cases

- $\alpha \gg \beta$ — the oxide is transparent, $P_{\text{int}} \to P_{\text{up}}$, and
  the result collapses to Level 1. This is the active study's situation.
- $\alpha \ll \beta$ — the oxide dominates, $P_{\text{int}} \to P_{\text{down}}$, and
  the result collapses to Level 2a.

Both limits are worth testing after any change to the solver, because they are the two
cases where a wrong answer still looks physically plausible.

## Procedure

Level 2b, in the order the code runs it:

1. Evaluate both materials at the operating temperature —
   `get_oxide_properties_at_T`, `get_metal_properties_at_T`. Both warn if $T$ falls
   outside the measured window.
2. Short-circuit if $P_{\text{up}} \le$ `min_pressure` (1e-20 Pa): set
   $P_{\text{int}} = P_{\text{down}} + \texttt{min\_pressure}$, return with
   `converged=False`. The result is an approximation, not a solve.
3. Bracket the root:

   ```python
   P_min = max(P_downstream + P_upstream * 1e-10, min_pressure)
   P_max = P_upstream * (1 - 1e-10)
   ```

4. Define the residual
   $f(P_{\text{int}}) = J_{ox}(P_{\text{int}}) - J_{\text{metal}}(P_{\text{int}})$ and
   call `brentq` on it.
5. Evaluate the flux at the root, from either side — they agree by construction.
6. Form $R_{ox}$ and the linearised $R_{\text{metal}}$ at $P_{\text{int}}$, take their
   ratio, and threshold it into `regime`.

The root is guaranteed to exist and be unique inside
$[P_{\text{down}}, P_{\text{up}}]$ by a sign argument:

| At | $J_{ox}$ | $J_{\text{metal}}$ | $f$ |
|---|---|---|---|
| $P_{\text{int}} \to P_{\text{down}}$ | large — full ΔP across the oxide | ≈ 0 — no driving force | $> 0$ |
| $P_{\text{int}} \to P_{\text{up}}$ | ≈ 0 — no driving force | large — full ΔP across the metal | $< 0$ |

$f$ is continuous and monotonically decreasing between those, so it crosses zero exactly
once.

## Design choices

**Why a numerical solve when a closed form exists here.** The equation is linear in
$P_{\text{int}}$ on the left and square-root on the right. Substituting
$u = \sqrt{P_{\text{int}}}$ turns it into a quadratic in $u$, which *is* solvable for
this particular pair of laws. But the same structure recurs at Levels 3, 4 and 6 with
additional nonlinearities that are not reducible, so the project solves all of them the
same numerical way rather than special-casing this one.

That closed form is not wasted: the documentation's own `check_equations.py` uses it as
an independent target, testing `brentq` against the exact quadratic root on the Example
Wall rather than against itself. They agree to about $2\times10^{-16}$.

**Why `brentq` rather than Newton.** The bracket is known and the sign change is
guaranteed, so a bracketing method cannot diverge. `brentq` combines bisection's
guaranteed convergence with faster interpolation when the function is well behaved,
which it is here.

**Why the regime ratio is linearised rather than exact.** There is no single resistance
for a $\sqrt{P}$ law, only a slope at a point. Reporting the derivative at
$P_{\text{int}}$ makes the comparison local and honest; a "global" metal resistance
would have to pick a pressure anyway and would hide that it had.

## Checks and failure modes

```{figure} ../figures/fig-10-level2-limit.svg
:alt: Series flux rising to the bare-metal value as the oxide is thinned to nothing.
:width: 90%

**Figure 10.** The limit check, in the notebooks' panel (D) idiom. Look at where the
curve flattens: as the oxide is thinned the series wall must climb to the bare-metal
flux, and it does. A solver that failed to recover the level below would flatten at the
wrong value — which is why this is drawn rather than asserted.
```

| Check | Expected | Where |
|---|---|---|
| $P_{\text{int}}$ matches the closed-form quadratic root | exact | `check_equations.py` — *L2b: P_int matches closed form* |
| flux continuity, $\lvert J_{ox} - J_{\text{metal}}\rvert$ | ≈ 0 | *L2b: flux balance* |
| series flux below **both** single-layer fluxes | always | *L2b: J_series < min(...)* |
| $L_{ox} \to 0$ recovers Level 1 | exact | *L2b → L1* |
| $J \propto P$ for the oxide alone | linearity | *L2a: J(P/2) = J(P)/2* |

| Symptom | Meaning |
|---|---|
| `converged: False`, $P_{\text{int}} \approx P_{\text{down}}$ | $P_{\text{up}}$ at or below `min_pressure`; approximated, not solved |
| Series flux above either single-layer flux | Impossible. A layer in series can only add resistance. |
| $L_{ox} \to 0$ does not recover Level 1 | The series composition is wrong, not the oxide data |
| Regime label flipping under a small change | Expected near a threshold; read the ratio, not only the label |
| log–log slope near 1.0 rather than 0.5 | The oxide is in control — a different study, or a much thicker coating |

:::{note}
**Worked values.** Real numbers from the active study — 316L / Cr₂O₃ at 873 K, 1 bar →
vacuum, 1 mm wall, 1 µm coating. Regenerated by `docs/_tools/worked_values.py`.

Level 2a, the oxide alone:

```text
D_ox        1.394877e-17  m²/s
K_ox        2.868829e+02  (see the units warning in equilibrium-models)
R_oxide     2.498961e+08  Pa·s·m²/mol
flux        4.001663e-04  mol/m²/s
```

Against Level 1's `4.327598e-06` for the same conditions, the oxide taken alone is
roughly **92× more permeable than the metal**. A 1 µm Cr₂O₃ layer at this temperature is
still not a barrier under Model 1, and that is what makes the study metal-limited
through Level 5. Worth absorbing early: it explains why coating parameters barely move a
Level 5 answer — and why they move a Level 5L6 answer a great deal, the oxide there
carrying atomic H under Sieverts.

Level 2b, the series wall:

```text
alpha (oxide permeance)   4.001663e-09  mol/m²/s/Pa
beta  (metal permeance)   1.368507e-08  mol/m²/s/Pa^0.5
converged                 True
P_interface               9.892438e+04  Pa   — 98.924% of P_up
flux                      4.304261e-06  mol/m²/s
flux_error                2.287e-19
R_ox                      2.498961e+08
R_metal                   4.596579e+10
resistance_ratio          5.436567e-03
regime                    metal_limited
```

The interface pressure sits at 98.9% of upstream, so the oxide consumes about 1.1% of
the available driving force and the metal takes the rest. The composite flux is 0.5%
below the bare-metal Level 1 value — the coating is very nearly invisible.

Log–log slope of flux against $P_{\text{up}}$, measured at $10^4$, $10^5$ and $10^6$ Pa:
`0.5085`, `0.5027`, `0.5009`. **This is the most informative line in the chapter.** The
composite system reproduces the *metal's* slope of 0.5, not the oxide's 1.0 — an
independent confirmation that the metal is in control, arrived at without consulting a
single resistance value. The slight excess over 0.5, shrinking with pressure, is the
oxide's small linear contribution.
:::

```{note}
An earlier version of this block was computed at a **48 nm** coating, where $\alpha$ was
`8.336797e-08` and the interface sat at 99.95% of upstream. The configured thickness is
now 1 µm. If you have inherited those figures, they are from that older geometry.
```
