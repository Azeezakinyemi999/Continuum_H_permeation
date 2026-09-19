# Level 5 — the full system

```{note}
**Model 1.** Dissociation sits at the oxide/metal interface, so the oxide carries
intact H₂ under Henry's law. The Level 6 family makes a different choice — see
{doc}`two-models` before comparing results across the two.
```

Level 5 is a defective oxide over a defective metal: Level 3's parallel paths,
with Level 4's microstructural metal underneath every one of them. It introduces
no new physics. What it introduces is **coupling**, and the interesting question
is whether the two corrections simply multiply.

They do not, quite, and the size of the discrepancy is a measure of how real the
coupling is.

## The assembly

Each path through the oxide — intact, pinhole, crack, grain boundary — sits above
its own patch of metal. For each one, solve the two-layer problem for its own
interface pressure:

$$J_{ox}^{(i)}(P_{\text{up}} \to P_{\text{int}}^{(i)})
= J_{\text{metal}}(P_{\text{int}}^{(i)} \to P_{\text{down}};\ D_{\text{eff}})$$

then combine by area fraction as at Level 3:

$$J_{\text{total}} = \sum_i j^{(i)}f^{(i)}, \qquad \sum_i f^{(i)} = 1$$

The only change from Level 3 is the metal's diffusivity: $D_L$ becomes
$D_{\text{eff}}$ from the microstructural model. The same `brentq` solve runs,
with the microstructure threaded through to the metal flux function.

<!-- verify-docs: allow calculate_full_system_flux -->

```{note}
The source design note presented Level 5 through a function called
`calculate_full_system_flux()`. **There is no such function.** The real entry
points are:

| Function | Use |
|---|---|
| `level5_model_wrapper` | canonical — takes a parameter dict, used by the SA |
| `calculate_parallel_path_flux_defective_metal` | direct call with explicit props |
| `calculate_PRF_defective_metal` | coating effectiveness with microstructure |

The wrapper and the direct call agree exactly (`3.196893e-06` from both), so use
whichever is convenient; the wrapper is preferable because it resolves every
input from the active study rather than requiring you to assemble props dicts
correctly.
```

## The corrections do not simply multiply

Level 3 gave `4.326496e-06`. Level 4 gave a diffusivity reduction factor
$\eta = 0.738861227$. If the two were independent, Level 5 would be their
product.

| Quantity | Value |
|---|---|
| Level 5 flux / Level 3 flux | 0.738910356 |
| Level 4 `overall_factor` | 0.738861227 |
| ratio | 1.000066493 |

The 6.6e-5 discrepancy is small but it is not noise, and its origin is
instructive: **reducing the metal's diffusivity changes where the interface
pressure settles.** A slower metal backs hydrogen up, raising
$P_{\text{int}}$ slightly, which shifts a little more of the driving force onto
the oxide and slightly alters the split. Level 3's interface pressure was
`9.994810e+04` Pa; Level 5's is `9.996165e+04` Pa.

So the levels are genuinely coupled through the interface solve, not merely
composed. For this study the coupling is worth 0.007% and you would be entitled
to ignore it — but in an oxide-limited system, where $P_{\text{int}}$ has room to
move, the same mechanism is much stronger.

## Verified recovery limits

A composite model should reduce to each of its parents when the corresponding
physics is switched off. All three hold exactly:

| Limit | Should recover | Result |
|---|---|---|
| `mode='none'` (no microstructure) | Level 3 | `4.326496e-06` — exact |
| defect fractions → 0 **and** `mode='none'` | Level 2b | `4.326475e-06` — exact |
| defect fractions → 0, microstructure on | Level 2b + Level 4 | `3.196881e-06` |
| `level5_model_wrapper` vs direct call | each other | agree to machine precision |

These are the checks worth re-running after touching either parent module,
because they fail loudly and locally, whereas a wrong composite flux looks
plausible.

## Reading the resistance fractions

Level 5 reports three fractions which sum to exactly 1:

```text
frac_oxide    3.757957e-04
frac_metal    0.979621
frac_defect   0.020004
sum           1.000000000
```

```{important}
These are not three of a kind. `frac_oxide` and `frac_metal` describe a **series**
split — how the driving pressure divides between the two layers — while
`frac_defect` describes a **parallel** split: the share of flux that bypasses the
intact oxide entirely.

The code makes them commensurable deliberately, by scaling the series split by
the share of flux that actually travels the series path. That is what allows all
three to sum to 1 and lets the parallel axis compete directly against the two
series terms when `assign_regime_L5` picks a winner.

Worth knowing: these fractions are taken straight from the solved interface
pressure, and so are **exact**. They deliberately avoid
`calculate_metal_resistance`, whose linearised $2L\sqrt{P}/(DK_s)$ form is only
valid for small $\Delta P$ — and these sweeps run with $P_{\text{down}} = 0$,
which is the worst case for that approximation. Level 2b's `resistance_ratio`
*does* use the linearised form, so the two diagnostics are computed differently
on purpose and should not be expected to agree numerically.
```

For the active study `frac_metal = 0.98`, which is the whole story of this
configuration in one number: the metal carries essentially all the resistance,
so `regime = metal`.

## PRF at Level 5 isolates the coating

`flux_bare_metal` is `3.197494e-06`, which is exactly the Level 1 flux multiplied
by $\eta$ — the bare-metal reference **carries the microstructure too**.

$$4.327598\times10^{-6} \times 0.738861227 = 3.197494\times10^{-6}\;\checkmark$$

This matters for interpretation. Because both sides of

$$\mathrm{PRF} = \frac{J_{\text{bare, microstructured}}}{J_{\text{coated, microstructured}}}$$

include trapping, the ratio measures **the coating alone**, not coating plus
microstructure. A PRF of 1.000188 therefore says the oxide does essentially
nothing — it does not say the microstructure does nothing. The microstructure's
effect is the 26% reduction hiding in both numerator and denominator, visible
only in `modification_factor`.

## Values at the default operating point

```text
flux                  3.196893e-06  mol/m²/s
flux_intact           3.132944e-06
flux_defect           6.394927e-08     (intact + defect = flux, exact)
flux_bare_metal       3.197494e-06
PRF                   1.000188
P_interface           9.996165e+04  Pa
D_metal               2.363871e-10  m²/s
D_eff                 1.746573e-10  m²/s
modification_factor   0.738861227
defect_enhancement    1.000004
permeability          1.010995e-11
permeance             1.010946e-08
regime                metal
```

`permeability` is the *apparent* permeability $J L_{\text{tot}}/\Delta\sqrt{P}$,
backed out of the flux above, so it agrees with `flux` and `regime` instead of
contradicting them. Note where it lands: $1.010995\times10^{-11}$ against a bare
metal $\Phi_m = 1.368507\times10^{-11}$, a ratio of 0.739 — exactly
`modification_factor`. The 48 nm oxide is non-limiting, so the whole wall reduces
to $D_{\text{eff}}K_s$, and the number reads metal-controlled just as `regime`
does. It is still not a material constant; see {doc}`permeability`.

`PRF` is finite here but becomes `nan` at Level 5L6, so prefer the Level 5 value
when you need a coating number.

## What Level 5 still cannot tell you

Dissociation at the gas-facing surface is still infinitely fast. Every level so
far has assumed that the surface imposes no resistance — that hydrogen arrives,
splits and dissolves without delay. On an oxide with a high dissociation barrier
that assumption can be the worst one in the model, and it is the last one the
hierarchy removes.

That is Level 6, and combining it with everything here gives L5L6 — the most
complete configuration the project offers.
