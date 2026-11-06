# Level 1 — the perfect metal

```{note}
**Model 1** by numbering, though nothing here turns on it: this level has no oxide, so
the Henry/Sieverts fork of {doc}`two-models` does not bite.
```

## Purpose

A bare metal wall of thickness $L$, with hydrogen gas at $P_{\text{up}}$ on one side and
$P_{\text{down}}$ on the other. No oxide, no defects, no microstructure, no surface
resistance. Hydrogen dissolves into the metal at both faces, diffuses down its
concentration gradient, and leaves.

This is the floor of the hierarchy. Every later level is a correction to it, and the
ratio of a later flux to the Level 1 flux is the standard way this project measures how
much a coating or a microstructure actually did.

**What it cannot tell you.** By construction it has no oxide, so it cannot describe a
coating; no defects, so it cannot describe a breached one; no traps, so it predicts a
diffusivity that is a pure lattice property; and no surface resistance, so it assumes
dissociation is infinitely fast. Each of those becomes a later level.

## At a glance

| | |
|---|---|
| **Inputs** | $D_m$, $K_s$, $L_m$, $P_{\text{up}}$, $P_{\text{down}}$ |
| **Outputs** | `flux`, `C_up`, `C_down`, `permeability` |
| **Code** | [`permeation_calc`](../api/permeation_calc) — `sieverts_concentration`, `fick_flux`, `calculate_simple_metal_flux` |
| **Solver** | None. Closed form. |

## Concepts

### Surface equilibrium: Sieverts' law

Hydrogen enters the metal as atoms, not molecules. The two-step equilibrium

$$\mathrm{H_2(gas)} \rightleftharpoons \mathrm{2H(ads)}
\rightleftharpoons \mathrm{2H(dissolved)}$$

gives a concentration proportional to the *square root* of pressure, because one
molecule yields two independent dissolved atoms. Combining the adsorption equilibrium
$\theta^2 \propto P$ with dissolution $C \propto \theta$:

$$C = K_s\sqrt{P}$$

with $K_s$ the Sieverts solubility constant [mol m⁻³ Pa⁻⁰·⁵]. Applied at both faces:

$$C(0) = K_s\sqrt{P_{\text{up}}}, \qquad C(L) = K_s\sqrt{P_{\text{down}}}$$

The square root is not a fitting choice. It is a consequence of the molecule splitting,
and it propagates into every flux law in the project.

### Steady-state diffusion

In the bulk, hydrogen obeys Fick's second law,
$\partial C/\partial t = D\,\partial^2 C/\partial z^2$. At steady state the time
derivative vanishes,

$$\frac{d^2C}{dz^2} = 0$$

whose solution is a straight line between the two boundary concentrations:

$$C(z) = C(0) + \left[C(L) - C(0)\right]\frac{z}{L}$$

A linear profile is worth pausing on: it means the concentration gradient is the same
everywhere in the wall, so there is no accumulation anywhere, which is what
{doc}`../STEADY_STATE` argues is adequate for the questions this model answers.

### The flux law

Fick's first law with a constant gradient gives

$$J = -D\frac{dC}{dz} = \frac{D}{L}\left[C(0) - C(L)\right]$$

and substituting the Sieverts boundary conditions yields the Level 1 result:

$$\boxed{\;J = \frac{D K_s}{L}\left(\sqrt{P_{\text{up}}} - \sqrt{P_{\text{down}}}\right)\;}$$

### The two signatures

```{figure} ../figures/fig-09-level1-sqrt-arrhenius.svg
:alt: Flux straight against the square root of pressure, and an Arrhenius plot whose slope recovers the permeability activation energy.
:width: 100%

**Figure 9.** Look at the two straight lines. On the left, flux against $\sqrt{P}$ is
straight **through the origin** — the Sieverts signature. On the right, the fitted
Arrhenius slope recovers $Q = E_D + H_s$ exactly, which is how a measured permeability
is decomposed. Both are limit checks: if either bends, the metal block is wrong.
```

With $P_{\text{down}} = 0$ the law collapses to $J \propto \sqrt{P_{\text{up}}}$, so on
log–log axes the flux against pressure has **slope 0.5**. That slope is a diagnostic:
measure it and you learn which mechanism controls. Slope 0.5 means Sieverts-limited bulk
diffusion; slope 1.0 means molecular transport under Henry's law; anything else means
something more interesting is happening, which is what Levels 3 and 6 exist to describe.
See {doc}`permeability` for the exponent as a general diagnostic.

### Permeability is the material property

Diffusivity and solubility never appear separately in the answer, only as their product:

$$\Phi \equiv D K_s \qquad [\mathrm{mol\,m^{-1}s^{-1}Pa^{-0.5}}]$$

A metal that dissolves a lot of hydrogen but holds it tightly can have the same
permeability as one that dissolves little but lets it move freely. This is why
literature tabulates $\Phi$, and why the configuration stores `Phi_ref` alongside
`D_ref` and `K_s_ref` as an alternate parameterisation of the same data.

### The resistance analogy, and where it breaks

Rearranging to

$$J = \frac{\sqrt{P_{\text{up}}} - \sqrt{P_{\text{down}}}}{R_{\text{metal}}},
\qquad R_{\text{metal}} = \frac{L}{\Phi}$$

makes the wall an Ohm's-law resistor, which is what licenses the series-and-parallel
resistance reasoning used from Level 2 onward.

```{important}
The analogy holds in $\sqrt{P}$, not in $P$. The driving force is $\Delta\sqrt{P}$, so
resistances may be added in series only when every layer is expressed in the same
$\sqrt{P}$ variable. This is exactly why the oxide needs care: below Level 6 it is
linear in $P$, not $\sqrt{P}$, so it is not a resistor in the same currency as the
metal. See {doc}`equilibrium-models`.
```

## Procedure

Three functions, each one step of the derivation, composed in this order:

1. `sieverts_concentration(K_s, P_up)` → $C(0) = K_s\sqrt{P_{\text{up}}}$. Raises on a
   negative pressure.
2. `sieverts_concentration(K_s, P_down)` → $C(L)$.
3. `fick_flux(D, C_up, C_down, L)` → $J = D(C_{\text{up}} - C_{\text{down}})/L$. Raises
   on a non-positive thickness or negative $D$.
4. `calculate_simple_metal_flux` composes the three, computes $\Phi = DK_s$, and returns
   both with a units dictionary. It warns rather than raises if
   $P_{\text{down}} > P_{\text{up}}$, returning a negative flux — backward flow is a
   legitimate question, not an error.

Temperature enters earlier: `get_metal_properties_at_T` evaluates $D$ and $K_s$ at the
operating temperature before any of this runs.

```{important}
`get_metal_properties_at_T` returns these as `D_metal` and `K_s_metal`, not `D` and
`K_s`. The naming difference has bitten callers before.
```

## Design choices

**Why the reference-temperature Arrhenius form.** Both $D$ and $K_s$ are
Arrhenius-activated, and the project uses

$$k(T) = k_{\text{ref}}\exp\left[\frac{-E}{R}\left(\frac{1}{T} - \frac{1}{T_{\text{ref}}}\right)\right]$$

implemented once in `arrhenius`. This is deliberately *not* the classical
$k(T) = k_0\exp(-E/RT)$. The two are algebraically equivalent, with
$k_0 = k_{\text{ref}}\exp(E/RT_{\text{ref}})$, but they differ in what they ask you to
trust. The reference form is anchored on a value measured at a stated temperature; the
classical form is anchored on a pre-exponential extrapolated to infinite temperature,
which is far from any measurement and correspondingly less certain.

The configuration stores both — `D_ref` with `E_D`, and `D_0` — but the model reads the
`*_ref` pair. `D_0`, `K_s0`, `Phi_0` and `Q_p` are derived conveniences, not independent
inputs.

Because $\Phi = D K_s$, permeability carries the **sum** of the two activation energies,
$Q_p = E_D + \Delta H_s$, which is why a permeability Arrhenius plot has a steeper slope
than a diffusivity one for the same material.

**Why a warning rather than an error on reversed pressures.** Backward flow is
physically meaningful and occasionally what the caller wants. Raising would make the
legitimate case unreachable.

## Checks and failure modes

Each of these is a property the closed form must have, and each is checked against the
running code rather than asserted.

| Check | Expected | Where |
|---|---|---|
| $P_{\text{down}} = 0$ matches $\Phi\sqrt{P_{\text{up}}}/L$ | exact | `check_equations.py` — *L1: flux* |
| $\Phi$ recovered from the flux is pressure-independent | constant | *L1: permeability* |
| $J(P/4) = J(P)/2$ | the $\sqrt{P}$ signature | *L1: J(P/4) = J(P)/2* |
| $P_{\text{up}} = P_{\text{down}}$ gives no driving force | $J = 0$ | `worked_values.py` |
| $J \propto 1/L$ | ×10 per decade | `worked_values.py` |

The pressure-independence check is the sharpest: recovering $\Phi = JL/\Delta\sqrt{P}$
at five pressures and getting the same number to sixteen digits confirms the $\sqrt{P}$
exponent is exactly $1/2$ in the implementation, not merely close.

| Symptom | Meaning |
|---|---|
| `ValueError: Pressure cannot be negative` | A pressure went below zero upstream of the call. |
| `ValueError: Thickness must be positive` | A thickness of zero reached `fick_flux`; the $L \to 0$ limit must be taken on the flux, not by passing zero. |
| Printed warning about backward flow | $P_{\text{down}} > P_{\text{up}}$. The negative flux is correct. |
| Flux an order of magnitude from expectation | Check whether `D_metal`/`K_s_metal` were read under the wrong names. |
| Slope on log–log is not 0.5 | Not Level 1. Something else is controlling — Levels 2, 3 or 6. |

```{note}
**Worked values.** For the active study — 316L at 873 K, 1 bar upstream, vacuum
downstream, 1 mm wall — this level gives $D_m$ = 2.363871e-10 m²/s,
$K_s$ = 5.789261e-02 mol m⁻³ Pa⁻⁰·⁵, $\Phi$ = 1.368507e-11, $C_{up}$ = 1.830725e+01
mol m⁻³ and a flux of 4.327598e-06 mol m⁻² s⁻¹. These are study-specific and
regenerated by `docs/_tools/worked_values.py`; the figures above use the Example Wall
instead and are study-independent.
```

## Where to go next

The immediate next question — what happens when you put an oxide on the surface — is
{doc}`equilibrium-models` for the physics and {doc}`level2-oxide` for the assembled
result.
