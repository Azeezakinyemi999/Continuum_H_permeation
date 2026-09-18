# Level 6 — surface kinetics

```{note}
**Model 2.** Dissociation moves to the gas/oxide surface, so what crosses the
oxide is now *atomic* H under Sieverts' law. This is a different oxide model from
Levels 1–5, not merely an added resistance, and it is why L5 and L5L6 results are
not directly comparable — see {doc}`two-models`.
```

Every level so far has assumed that hydrogen arriving at the wall splits and
dissolves instantly, so the gas-facing surface imposes no resistance. Level 6
removes that assumption: dissociative adsorption proceeds at a finite rate, the
surface coverage $\theta$ becomes an unknown to be solved for, and the wall sees
not the gas pressure but an *effective* pressure set by that coverage.

Level 6 composes with the other levels rather than replacing them, so it appears
as `L1+L6`, `L2a+L6`, and — combined with everything — `L5L6`.

## Three fluxes, one unknown

With surface kinetics active there are three transport steps in series, and at
steady state all three carry the same flux:

$$J_{\text{surface}} = k_{\text{diss}}P_{\text{up}}(1-\theta)^2 - k_{\text{recomb}}\theta^2$$

$$J_{\text{oxide}} = \alpha\left[g(\theta) - \sqrt{P_{\text{int}}}\right]$$

$$J_{\text{metal}} = \beta\left[\sqrt{P_{\text{int}}} - \sqrt{P_{\text{down}}}\right]$$

with $g(\theta) = \sqrt{P_{\text{virtual}}}$ from {doc}`foundations` and the
permeances $\alpha, \beta$ from {doc}`level2-oxide`.

That looks like two equations in two unknowns ($\theta$ and $P_{\text{int}}$),
but the oxide-equals-metal condition can be solved for $P_{\text{int}}$ in closed
form:

$$\alpha g(\theta) + \beta\sqrt{P_{\text{down}}} = (\alpha+\beta)\sqrt{P_{\text{int}}}
\;\Longrightarrow\;
\sqrt{P_{\text{int}}}(\theta) = \frac{\alpha g(\theta) + \beta\sqrt{P_{\text{down}}}}{\alpha + \beta}$$

```{important}
$\sqrt{P_{\text{int}}}$ is an **analytical function of $\theta$**. Substituting it
into the surface equation leaves a single residual in a single unknown, bracketed
on $\theta \in (0,1)$ — a physically bounded interval requiring no initial guess.

This is why Level 6 is cheap despite being the most coupled level in the
hierarchy. `sqrt_P_int_from_theta` performs the substitution and
`surface_flux_residual` is the residual that `brentq` drives to zero.
```

## Flux matching is the correctness check

Because the solve enforces only one residual, the other equalities are
predictions and can be checked independently. Measured at the default operating
point:

| Configuration | Fluxes compared | Worst mismatch |
|---|---|---|
| L1+L6 | surface vs metal | `1.23e-15` |
| L2+L6 | surface vs oxide vs metal | `6.09e-16` |

Machine precision in both. This is the check to run after any change to the
surface module, because a wrong $\theta$ breaks flux equality immediately and
visibly, whereas it perturbs the reported flux only subtly.

## The result that matters: the surface suppresses pressure, not flux directly

Take the bare metal with surface kinetics, `L1+L6`, at 873 K and 1 bar:

```text
theta                 0.573614764
P_int = P_virtual     1190.673636   Pa
J_ss                  4.722188e-07  mol/m²/s
R_surface             2.0924e+06
R_metal               7.3072e+07
fraction_surface      2.78%
fraction_metal        97.22%
rate_limiting         metal
```

Level 1 without surface kinetics gave `4.327598e-06` — so switching the surface
on **cut the flux by a factor of 9.16**. Yet the resistance decomposition
attributes only 2.78% of the resistance to the surface, and `rate_limiting`
reports `metal`.

Both are correct, and reconciling them is the key to reading Level 6 output.

If the surface were infinitely fast, $\theta$ would equilibrate with the gas at
$\theta_{eq} = 0.924975$. Instead it sits at $0.573615$ — the surface is
**starved**, because dissociation cannot keep up with what diffusion removes. A
depleted coverage means a low virtual pressure, and the wall downstream sees that
virtual pressure rather than $P_{\text{up}}$:

$$\frac{P_{\text{up}}}{P_{\text{virtual}}} = 83.9861
\qquad\Longrightarrow\qquad
\sqrt{83.9861} = 9.1644$$

which is exactly the measured flux ratio, because flux scales as $\sqrt{P}$.

```{warning}
**`rate_limiting` and `fraction_surface` do not measure what the surface cost
you.** They are computed at the *solved* operating point, so they answer "given
this coverage, where does the remaining driving force drop?" — not "how much flux
was lost relative to an infinitely fast surface?"

The surface acts by **depressing the effective upstream pressure**, and that loss
happens before the series decomposition begins. A surface holding 2.78% of the
resistance can still be responsible for a ninefold flux reduction.

This is why `calculations.sensitivity` documents `flux` as the primary metric —
its comment notes it is *the only metric that responds to surface kinetics*.
`permeability` is now derived from that same flux, so it does move with the
surface — but only as the flux moves, and it remains an apparent value that means
nothing without its operating point ({doc}`permeability`). `frac_surface` sees
only the residual series share. If you want to know whether surface kinetics
matter, compare fluxes with and without them, or inspect $\theta$ against
$\theta_{eq}$.
```

## Reading coverage

$\theta$ is the most diagnostic single number Level 6 produces:

| $\theta$ relative to $\theta_{eq}$ | Interpretation |
|---|---|
| $\theta \approx \theta_{eq}$ | surface fast; Level 6 reduces to the level below it |
| $\theta \ll \theta_{eq}$ | surface starved; dissociation is throttling the system |
| $\theta \to 1$ | saturated; $g(\theta) \to \infty$ and the isotherm stiffens |

The $\theta \to 1$ limit is worth respecting: `g_theta` returns `inf` at
$\theta \geq 1$ by design, which keeps the bracketing well-posed rather than
producing a silent overflow.

## The full system: L5L6

Combining surface kinetics with the defective oxide and defective metal gives
the most complete configuration available:

```text
flux               2.720783e-06  mol/m²/s
theta              0.616276
frac_surface       2.923902e-04
frac_oxide         1.067605e-01
frac_metal         8.929471e-01
sum                1.000000000
regime             metal
PRF                nan
k_diss  (oxide)    2.370190e-09
k_diss_metal       2.498430e-10
K_eq    (oxide)    2.798918e-05
K_eq_metal         3.349488e-03
```

Against Level 5's `3.196893e-06`, surface kinetics reduce the flux by 15% here.
Note the two separate kinetic pairs: the oxide surface and the bare metal exposed
at pinholes have their own constants, as {doc}`equilibrium-models` describes.

```{important}
**`frac_oxide` jumps by a factor of 284 between Level 5 and Level 5L6** — from
`3.757957e-04` to `1.067605e-01` — even though the oxide's thickness and material
properties are unchanged.

This is the Henry/Sieverts split of {doc}`equilibrium-models` showing up in the
numbers. Without Level 6 the oxide transports intact molecules under Henry's law
and, with the Sieverts-calibrated `K_ox`, comes out extremely permeable. With
Level 6 the species entering the oxide is atomic, transport becomes Sieverts, and
the oxide's share of the resistance rises by more than two orders of magnitude.

The jump is therefore expected given the model structure, but its *size* is
inflated by the unresolved `K_ox` units issue. Do not read it as a physical
finding about Cr₂O₃ until that constant is split.
```

`PRF` is `nan` from this wrapper at default parameters. Use the Level 5 value
when a coating-effectiveness number is needed.

## Entry points

| Function | Configuration |
|---|---|
| `solve_steady_state_flux_L1L6` | bare metal + surface |
| `solve_steady_state_flux_L2aL6` | oxide only + surface |
| `solve_steady_state_flux` | oxide + metal + surface, resolving properties by name |
| `solve_steady_state_flux_direct` | same, with explicit property arguments |
| `calculate_full_model_flux_L346_v2` | the full L5L6 assembly |
| `level5L6_model_wrapper` | canonical — parameter dict, used by the SA |

`solve_steady_state_flux_L1L6` returns `J_ss`, not `flux` — there is no `flux`
key on it, which has caught callers before.

## What the hierarchy still does not include

With Level 6 the resistance chain is complete: surface, oxide, defect paths,
metal microstructure. What remains outside the model is time dependence — every
level is steady-state, and `STEADY_STATE.md` sets out why that is sufficient for
the questions asked here and what would invalidate it.

The next practical step is not another level but a sensitivity analysis: with
roughly thirty parameters and a regime that shifts as they vary, the useful
question is which of them the answer actually depends on. That is the
regime-stratified workflow.
