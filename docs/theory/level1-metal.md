# Level 1 — the perfect metal

```{note}
**Model 1** by numbering, though nothing here turns on it: this level has no
oxide, so the Henry/Sieverts fork of {doc}`two-models` does not bite.
```

## What it models

A bare metal wall of thickness $L$, with hydrogen gas at $P_{\text{up}}$ on one
side and $P_{\text{down}}$ on the other. No oxide, no defects, no
microstructure, no surface resistance. Hydrogen dissolves into the metal at both
faces, diffuses down its concentration gradient, and leaves.

This is the floor of the hierarchy. Every later level is a correction to it, and
the ratio of a later flux to the Level 1 flux is the standard way this project
measures how much a coating or a microstructure actually did.

## Derivation

### Surface equilibrium: Sieverts' law

Hydrogen enters the metal as atoms, not molecules. The two-step equilibrium

$$\mathrm{H_2(gas)} \rightleftharpoons \mathrm{2H(ads)}
\rightleftharpoons \mathrm{2H(dissolved)}$$

gives a concentration proportional to the *square root* of pressure, because one
molecule yields two independent dissolved atoms. Combining the adsorption
equilibrium $\theta^2 \propto P$ with dissolution $C \propto \theta$:

$$C = K_s\sqrt{P}$$

with $K_s$ the Sieverts solubility constant [mol m⁻³ Pa⁻⁰·⁵]. Applied at both
faces:

$$C(0) = K_s\sqrt{P_{\text{up}}}, \qquad C(L) = K_s\sqrt{P_{\text{down}}}$$

The square root is not a fitting choice. It is a consequence of the molecule
splitting, and it propagates into every flux law in the project.

### Steady-state diffusion

In the bulk, hydrogen obeys Fick's second law,
$\partial C/\partial t = D\,\partial^2 C/\partial z^2$. At steady state the time
derivative vanishes,

$$\frac{d^2C}{dz^2} = 0$$

whose solution is a straight line between the two boundary concentrations:

$$C(z) = C(0) + \left[C(L) - C(0)\right]\frac{z}{L}$$

A linear profile is worth pausing on: it means the concentration gradient is the
same everywhere in the wall, so there is no accumulation anywhere, which is what
{doc}`../STEADY_STATE` argues is adequate for the questions this model answers.

### The flux law

Fick's first law with a constant gradient gives

$$J = -D\frac{dC}{dz} = \frac{D}{L}\left[C(0) - C(L)\right]$$

and substituting the Sieverts boundary conditions yields the Level 1 result:

$$\boxed{\;J = \frac{D K_s}{L}\left(\sqrt{P_{\text{up}}} - \sqrt{P_{\text{down}}}\right)\;}$$

## What the result tells you

### The √P signature

With $P_{\text{down}} = 0$ the law collapses to $J \propto \sqrt{P_{\text{up}}}$,
so on log–log axes the flux against pressure has **slope 0.5**. That slope is a
diagnostic: measure it and you learn which mechanism controls. Slope 0.5 means
Sieverts-limited bulk diffusion; slope 1.0 means molecular transport (the oxide
under Henry's law — see {doc}`equilibrium-models`); anything else means something
more interesting is happening, which is what Levels 3 and 6 exist to describe.

### Permeability is the material property

Diffusivity and solubility never appear separately in the answer, only as their
product:

$$\Phi \equiv D K_s \qquad [\mathrm{mol\,m^{-1}s^{-1}Pa^{-0.5}}]$$

A metal that dissolves a lot of hydrogen but holds it tightly can have the same
permeability as one that dissolves little but lets it move freely. This is why
literature tabulates $\Phi$, and why the configuration stores `Phi_ref`
alongside `D_ref` and `K_s_ref` as an alternate parameterisation of the same
data.

### The resistance analogy, and where it breaks

Rearranging to

$$J = \frac{\sqrt{P_{\text{up}}} - \sqrt{P_{\text{down}}}}{R_{\text{metal}}},
\qquad R_{\text{metal}} = \frac{L}{\Phi}$$

makes the wall an Ohm's-law resistor, which is what licenses the
series-and-parallel resistance reasoning used from Level 2 onward.

```{important}
The analogy holds in $\sqrt{P}$, not in $P$. The driving force is
$\Delta\sqrt{P}$, so resistances may be added in series only when every layer is
expressed in the same $\sqrt{P}$ variable. This is exactly why the oxide needs
care: below Level 6 it is linear in $P$, not $\sqrt{P}$, so it is not a resistor
in the same currency as the metal. See {doc}`equilibrium-models`.
```

## In the code

Three functions, each one step of the derivation:

| Function | Step |
|---|---|
| `sieverts_concentration` | $C = K_s\sqrt{P}$ — the boundary condition |
| `fick_flux` | $J = D(C_{\text{up}} - C_{\text{down}})/L$ — the transport step |
| `calculate_simple_metal_flux` | the two composed, returning flux and diagnostics |

For the active study (316L at 873 K, 1 bar upstream, vacuum downstream, 1 mm
wall):

```text
D_metal        2.363871e-10  m²/s
K_s_metal      5.789261e-02  mol/m³/Pa^0.5
Φ = D·K_s      1.368507e-11  mol/m/s/Pa^0.5
C_up           1.830725e+01  mol/m³
flux           4.327598e-06  mol/m²/s
```

Note `get_metal_properties_at_T` returns these as `D_metal` and `K_s_metal`, not
`D` and `K_s` — a naming difference that has bitten callers before.

## Verified limit checks

Each of these is a property the closed form must have, and each is checked
against the running code rather than asserted. Reproduce them with
`docs/_tools/worked_values.py`.

| Check | Expected | Measured |
|---|---|---|
| $P_{\text{down}} = 0$ matches $\Phi\sqrt{P_{\text{up}}}/L$ | exact | relative error `0.00e+00` |
| $P_{\text{up}} = P_{\text{down}}$ gives no driving force | $J = 0$ | `0.0` exactly |
| $\Phi$ independent of pressure over $10^3$–$10^7$ Pa | constant | spread `1.18e-16` |
| $J \propto 1/L$ | ×10 per decade | `10.0000`, `1.0000`, `0.1000` |

The third is the sharpest of the four: recovering $\Phi = JL/\sqrt{P}$ from the
flux at five pressures and getting the same number to sixteen digits confirms the
$\sqrt{P}$ exponent is exactly $1/2$ in the implementation, not merely close.

## Temperature dependence

Both $D$ and $K_s$ are Arrhenius-activated. The project uses the
**reference-temperature** form throughout:

$$k(T) = k_{\text{ref}}\exp\left[\frac{-E}{R}\left(\frac{1}{T} - \frac{1}{T_{\text{ref}}}\right)\right]$$

implemented once in `arrhenius` and used for every temperature-dependent
property.

```{note}
This is deliberately *not* the classical $k(T) = k_0\exp(-E/RT)$ form. The two
are algebraically equivalent, with $k_0 = k_{\text{ref}}\exp(E/RT_{\text{ref}})$,
but they differ in what they ask you to trust. The reference form is anchored on
a value measured at a stated temperature; the classical form is anchored on a
pre-exponential extrapolated to infinite temperature, which is far from any
measurement and correspondingly less certain.

The configuration stores both — `D_ref` with `E_D`, and `D_0` — but the model
reads the `*_ref` pair. `D_0`, `K_s0`, `Phi_0` and `Q_p` are derived
conveniences, not independent inputs.
```

Because $\Phi = D K_s$, permeability carries the *sum* of the two activation
energies, $Q_p = E_D + \Delta H_s$, which is why a permeability Arrhenius plot
has a steeper slope than a diffusivity one for the same material.

## What Level 1 cannot tell you

By construction it has no oxide, so it cannot describe a coating; no defects, so
it cannot describe a breached one; no traps, so it predicts a diffusivity that is
a pure lattice property; and no surface resistance, so it assumes dissociation is
infinitely fast.

Each of those becomes a later level. The immediate next question — what happens
when you put an oxide on the surface — is {doc}`equilibrium-models` for the
physics and Level 2 for the assembled result.
