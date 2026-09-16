# How to run a regime-stratified sensitivity analysis

## What this answers, and why it is stratified

A sensitivity analysis asks: if each input is wiggled, how much does the output
wiggle? Which parameters actually matter?

The twist here is that **the answer depends on which regime the system is in**.
The wall is three resistances in series:

```text
   gas ──[ R_surface ]──[ R_oxide ]──[ R_metal ]──> downstream
          dissociation     oxide        metal
                         diffusion    diffusion
```

Whichever is largest is the bottleneck, and that is the regime. Surface-kinetics
parameters matter in the surface regime; metal-transport parameters matter in the
metal regime. A conventional SA averages over all regimes and reports a blend
that describes none of them. This pipeline computes sensitivity **separately
within each regime**.

| Model | Regimes |
|---|---|
| L5L6 | `surface`, `oxide`, `metal` |
| L5 (no surface kinetics) | `oxide`, `metal`, `defect` |

## The pipeline

```text
STAGE 1 — targeted sampling, once per regime preset
  ┌──────────────────────────────────────────────────────┐
  │  LHS draw from the preset's ranges                   │
  │        ↓                                             │
  │  run the model (calculate_full_model_flux_L346_v2)   │
  │        ↓                                             │
  │  label the regime  (assign_regime: argmax of the     │
  │                     flux-weighted fractions)         │
  │        ↓                                             │
  │  keep the row only if label == target regime         │
  └──────────────────────────────────────────────────────┘
        ↓  three clusters: metal / surface / oxide
STAGE 2 — given-data SA (PAWN + Borgonovo δ) on each cluster
        ↓
STAGE 3 — compare regimes side by side
```

```{important}
**Regime labelling is not a stage that runs first.** It is one line invoked once
per draw, immediately after that draw's model run, nested inside Stage 1's loop.
It cannot run earlier: labelling needs the flux-weighted resistance fractions,
which exist only after a model run, which exists only after a parameter draw.
```

### Why the sampling has to be aimed

Uniform random sampling almost never lands in the surface-limited regime — **0
hits in 1500 random runs**. Without targeting there is simply no surface cluster
to analyse. Each regime therefore gets a *preset*: a set of parameter ranges
biased toward producing that regime, with every draw still labelled and filtered
on its actual computed regime, so no row is assumed into its cluster.

```{warning}
**Presets suppress the very parameters they pin.** A preset narrows or fixes
parameter ranges to steer the sampling. A parameter that has been narrowed cannot
show sensitivity, because it barely varies — so a low score for a
preset-suppressed parameter carries no information at all.

Always check a parameter's range in the preset before concluding it does not
matter. `presets_without(presets, 'name')` removes a parameter from the sampled
set so it can be pinned through `fixed_params` instead, which is the supported
way to take something out of competition deliberately.
```

### Why not Morris or Sobol

Both need a *structured* sample — Morris needs intact trajectories, Sobol needs
paired Saltelli matrices — and filtering such a sample by regime destroys the
design that makes the estimator valid. A regime-conditioned Sobol would have to
resample a rectangular sub-box and accept leakage across a curved regime
boundary. Given-data estimators need no design at all, so they are strictly
better suited here. That path was removed; this pipeline is given-data only, and
consequently produces **no variance decomposition** — there are no S1/ST shares
and no pairwise interaction terms.

## The two estimators

| Method | Question it answers | Mechanism |
|---|---|---|
| **PAWN** | How much does the output *distribution shape* change when this input is pinned? | Kolmogorov–Smirnov distance between conditioned and unconditioned CDFs |
| **Borgonovo δ** | How much does the output *density* shift when this input is pinned? | Area between conditioned and unconditioned densities, δ ∈ [0,1] |

Both are moment-independent and valid on arbitrary scattered points, which is
exactly what a regime-filtered cluster is.

```{note}
**Flux is analysed as log₁₀(flux).** Flux spans roughly ten orders of magnitude,
and density-based estimators get hijacked by a few enormous values — δ came out
as 46 and −93389 in early runs, impossible for an index defined on [0,1]. Taking
the logarithm asks what controls the *order of magnitude* of flux, which is the
question worth answering anyway. `theta` is already on [0,1] and stays linear.
Controlled by `LOG_METRICS_DEFAULT`, currently `('flux', 'permeability')`.
```

## Sizing the draws

Cluster size is chosen adaptively rather than fixed. `size_draws_for_target`
runs a short probe (`probe_n=250`) per preset to measure its yield — what
fraction of draws actually land in the target regime — then scales up to reach
`target_cluster=1500` rows, with a `safety=1.15` margin and an `n_max=200000`
cap. `run_targeted_regime_scans(N_per_regime=None)` invokes this automatically
and enforces `min_cluster=300`.

```{note}
Earlier versions of this workflow referred to a `DEFAULT_N_PER_REGIME` constant.
It no longer exists — the adaptive sizing above replaced it.
```

## Reading the results honestly

This is the part to read before interpreting any table.

### The dummy parameter sets the significance threshold

Neither δ nor PAWN returns zero for an input that does nothing; both are
finite-sample estimates and their noise does not average to exactly zero. The
standard remedy, from the PAWN authors themselves, is to include a **dummy
parameter** the model never sees and treat its index as the noise floor. Inputs
scoring appreciably above the dummy are resolved as sensitive; the rest are not
resolved.

`givendata_sensitivity_by_regime(..., n_dummy=N)` implements this, reporting
`floor` and `delta_over_floor`. Because `floor = max(dummy indices)`, it is an
order statistic over `n_dummy` samples.

```{warning}
**The default `n_dummy=3` is too small, and it inflates every margin.** With
only three dummies the maximum underestimates the true noise ceiling.

Measured on the L5 clusters, raising `n_dummy` from 3 to 20 lifted the floor by
9–34% and cut the number of parameters scoring above 1.0× the floor **from 21 to
8** in the defect regime. Two thirds of the apparently significant parameters
were noise.

Pass `n_dummy >= 20` for anything that will be published. The shipped default is
convenient for exploration, not defensible for a result.
```

### Below the floor does not mean "no physical effect"

This is the caveat most easily got wrong. Puy, Lo Piano & Saltelli (2020) ran a
sensitivity analysis *of PAWN* and found that on the Morris test function, inputs
that are genuinely influential but act **purely through interactions** cannot be
distinguished from a dummy:

> the PAWN index might be incapable to differentiate between non-influential
> model inputs and influential model inputs whose effect in the model output is
> fully through interactions

So a below-floor score means *undetected by this estimator on this sample*, not
*no effect*.

**A worked example from this model.** `f_pinhole` scores at or below the floor in
the L5 metal cluster. Concluding that pinholes do not matter would be wrong. A
pinhole exposes bare metal, so hydrogen enters as though the oxide were absent —
which is exactly how the defect path is constructed. The effect is real but
*interactive*: how much a bypass buys you depends on how much resistance the
oxide was contributing. In the metal-limited regime the oxide was never the
bottleneck, so the marginal gain is small. An interaction-driven input with a
small marginal effect is precisely the documented blind spot.

The same parameter is worth about 2.2 decades of flux in the oxide regime. The
regime stratification is what makes that visible; a blended SA would have
averaged it away.

### High dimensionality is a flagged risk

Puy et al. also report that PAWN "especially underperformed" on test functions
with 8 and 20 inputs, and raise this as a red flag for models with tens of
parameters. This pipeline varies **36**. Treat PAWN rankings as corroborating
evidence alongside δ, not as a standalone verdict.

## What is varied, and what is not

| | L5 | L5L6 |
|---|---|---|
| parameters with defaults | 36 | 46 |
| parameters actually swept | 28 | **36** |

The four reference temperatures — `T_ref_metal`, `T_ref_oxide`,
`T_ref_surface`, `T_ref_surface_metal` — are held **fixed** and appear in no
swept range.

```{note}
This is deliberate, not an oversight. In
$X(T) = X_{\text{ref}}\exp[-E/R\,(1/T - 1/T_{\text{ref}})]$, the reference
temperature is *measurement metadata* — the temperature at which the reference
value was measured — and it is redundant with $X_{\text{ref}}$ through the
prefactor $A = X_{\text{ref}}\exp(E/RT_{\text{ref}})$. Varying it independently
injects an artificial prefactor sweep of three to four orders of magnitude.

Property uncertainty is captured through the `*_ref` values and the activation
energies instead. The *operating* temperature is varied, over 573–1273 K.
```

```{warning}
That operating range extrapolates well past the oxide's validated window of
[473, 773] K, and the wrappers used by this pipeline do not warn about it — see
{doc}`../getting-started`.
```

## Validate the configuration first

```python
from calculations.sensitivity import check_against_config
check_against_config()
```

On the active study this currently reports **four pre-existing problems**:
`D_ref`, `k_diss_metal_ref` and `E_diss_metal` defaults lying outside their own
sweep ranges. A default outside its sweep range means the nominal case is never
evaluated, so resolve these before trusting a sensitivity result on those
parameters. See {doc}`switch-study`.

## Stage-to-code map

| Stage | Functions in `calculations.sensitivity` |
|---|---|
| label a draw | `assign_regime`, `assign_regime_L5`, `level5L6_model_wrapper(return_full_record=True)` |
| targeted scans | `REGIME_PRESETS`, `size_draws_for_target`, `run_global_lhs_scan`, `run_targeted_regime_scans`, `load_regime_scans`, `partition_by_regime`, `plot_regime_exploration` |
| given-data SA | `givendata_sensitivity_by_regime`, `summarize_givendata`, `plot_givendata_results` |
| geometry | `plot_regime_geometry` |
| compare regimes | `regime_comparison_matrix`, `plot_regime_comparison_heatmap` |
| parallel coordinates | `parallel_coordinates_samples`, `parallel_coordinates_sensitivity`, `top_drivers` |
| utilities | `presets_without`, `check_against_config` |

All verified to exist.

## Notebooks and outputs

The analysis runs in `sensitivity_regime_L5L6.ipynb` and its L5 counterpart; the
parallel-coordinates views live in `regime_parallel_coords.ipynb` and
`regime_parallel_coords_L5.ipynb`.

```{warning}
**Run the SA notebooks before the parallel-coordinates ones.** The SA notebooks
write scan CSVs into `sa_results*/`, which is gitignored; the
parallel-coordinates notebooks only *read* those files and cannot regenerate
them. On a fresh clone the visualisation notebooks will fail until the scans have
been produced.
```

Scans are cached on disk so a re-run does not repeat the model evaluations —
`load_regime_scans` picks up an existing set. Given-data SA adds no model
evaluations at all: it runs on the rows Stage 1 already produced.
