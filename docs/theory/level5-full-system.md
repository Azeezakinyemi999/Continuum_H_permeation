# Level 5 — the full system

```{note}
**Model 1.** Dissociation sits at the oxide/metal interface, so the oxide carries
intact H₂ under Henry's law. The Level 6 family makes a different choice — see
{doc}`two-models` before comparing results across the two.
```

## Purpose

Level 5 is a defective oxide over a defective metal: Level 3's parallel paths,
with Level 4's microstructural metal underneath every one of them. It introduces
no new physics. What it introduces is **coupling**, and the interesting question
is whether the two corrections simply multiply.

They do not, quite, and the size of the discrepancy is a measure of how real the
coupling is.

**What it cannot tell you.** Dissociation at the gas-facing surface is still infinitely
fast — the last assumption the hierarchy removes, and the subject of
{doc}`level6-surface-kinetics`.

## At a glance

| | |
|---|---|
| **Inputs** | Every Level 3 and Level 4 input, over one operating point |
| **Outputs** | `flux`, `frac_oxide`, `frac_metal`, `frac_defect`, `regime`, `PRF`, `D_eff`, `P_interface` |
| **Code** | `calculate_parallel_path_flux_defective_metal`, `level5_model_wrapper` |
| **Solver** | One `brentq` interface solve **per path**, each with the microstructural metal |

## Concepts

### The assembly

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

The wrapper and the direct call agree exactly (`3.182694e-06` from both), so use
whichever is convenient; the wrapper is preferable because it resolves every
input from the active study rather than requiring you to assemble props dicts
correctly.
```

### The corrections do not simply multiply

```{figure} ../figures/fig-13-level5-not-multiplicative.svg
:alt: Bars comparing the baseline, each correction alone, their naive product, and the true Level 5 result.
:width: 90%

**Figure 13.** Look at the last two bars. Having measured the oxide-defect correction
and the metal-microstructure correction separately, it is natural to multiply them —
the hatched bar. The true Level 5 result is the red one. On the Example Wall the
product is off by **+8.2%**, because both corrections act through the interface
pressure and are therefore coupled, not independent.
```

Level 3 gave `4.304704e-06`. Level 4 gave a diffusivity reduction factor
$\eta = 0.738327048$. If the two were independent, Level 5 would be their
product.

| Quantity | Value |
|---|---|
| Level 5 flux / Level 3 flux | 0.739352578 |
| Level 4 `overall_factor` | 0.738327048 |
| ratio | 1.001388992 |

The 1.4e-3 discrepancy is small but it is not noise, and its origin is
instructive: **reducing the metal's diffusivity changes where the interface
pressure settles.** A slower metal backs hydrogen up, raising
$P_{\text{int}}$ slightly, which shifts a little more of the driving force onto
the oxide and slightly alters the split. Level 3's interface pressure was
`9.892438e+04` Pa; Level 5's is `9.920472e+04` Pa.

So the levels are genuinely coupled through the interface solve, not merely
composed. For this study the coupling is worth 0.007% and you would be entitled
to ignore it — but in an oxide-limited system, where $P_{\text{int}}$ has room to
move, the same mechanism is much stronger. It is already twenty times larger here
than it was at the older 48 nm coating, for exactly that reason.

### Reading the resistance fractions

Level 5 reports three fractions which sum to exactly 1:

```text
frac_oxide    7.793176e-03
frac_metal    0.972132
frac_defect   0.020074
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

For the active study `frac_metal = 0.97`, which is the whole story of this
configuration in one number: the metal carries essentially all the resistance,
so `regime = metal`.

### PRF at Level 5 isolates the coating

`flux_bare_metal` is `3.195183e-06`, which is exactly the Level 1 flux multiplied
by $\eta$ — the bare-metal reference **carries the microstructure too**.

$$4.327598\times10^{-6} \times 0.738327048 = 3.195183\times10^{-6}\;\checkmark$$

This matters for interpretation. Because both sides of

$$\mathrm{PRF} = \frac{J_{\text{bare, microstructured}}}{J_{\text{coated, microstructured}}}$$

include trapping, the ratio measures **the coating alone**, not coating plus
microstructure. A PRF of 1.003924 therefore says the oxide does almost nothing — it
does not say the microstructure does nothing. The microstructure's effect is the 26%
reduction hiding in both numerator and denominator, visible only in
`modification_factor`.

## Checks and failure modes

### Verified recovery limits

A composite model should reduce to each of its parents when the corresponding
physics is switched off. All three hold exactly:

| Limit | Should recover | Result |
|---|---|---|
| `enable_trapping=False`, `enable_gb_enhancement=False` | Level 3 | `4.304704e-06` — exact |
| the above **and** defect fractions → 0 | Level 2b | `4.304261e-06` — exact |
| defect fractions → 0, microstructure on | Level 2b + Level 4 | `3.182452e-06` |
| `level5_model_wrapper` vs direct call | each other | agree to machine precision |

These are the checks worth re-running after touching either parent module,
because they fail loudly and locally, whereas a wrong composite flux looks
plausible.

### Worked values at the default operating point

```{note}
**Worked values.** Real numbers from the active study — 316L / Cr₂O₃ at 873 K, 1 bar
upstream, vacuum downstream, a 1 mm wall and a 1 µm coating. They are specific to that
study and are regenerated by `docs/_tools/worked_values.py`, whose `--check` mode
diffs them so a model change surfaces here rather than as quietly wrong prose. The
figures in this chapter use the Example Wall instead and are study-independent.
```

```text
flux                  3.182694e-06  mol/m²/s
flux_intact           3.118803e-06
flux_defect           6.389090e-08     (intact + defect = flux, exact)
flux_bare_metal       3.195183e-06
PRF                   1.003924
P_interface           9.920472e+04  Pa
D_metal               2.363871e-10  m²/s
D_eff                 1.745310e-10  m²/s
modification_factor   0.738327048
defect_enhancement    1.000076
permeability          1.007463e-11
permeance             1.006456e-08
regime                metal
```

`permeability` is the *apparent* permeability $J L_{\text{tot}}/\Delta\sqrt{P}$,
backed out of the flux above, so it agrees with `flux` and `regime` instead of
contradicting them. Note where it lands: $1.007463\times10^{-11}$ against a bare
metal $\Phi_m = 1.368507\times10^{-11}$, a ratio of $0.736$ — close to
`modification_factor` $= 0.738$, but not equal to it. The gap is the coating: at
1 µm the oxide carries `frac_oxide` $= 7.8\times10^{-3}$ of the resistance, so the
wall does not quite reduce to $D_{\text{eff}}K_s$. The number still reads
metal-controlled, just as `regime` does. It is not a material constant; see
{doc}`permeability`.

```{note}
An earlier version of this block was computed at a **48 nm** coating, where the oxide
really was negligible and the ratio matched `modification_factor` to four digits. The
configured thickness is now 1 µm. If you have inherited a PRF of 1.000188, it is from
that older geometry.
```

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
