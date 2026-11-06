# Stage 5 — Indices, and reading them

## Purpose

Compute PAWN and δ for every parameter inside every regime, rank them, and present the
rankings so they can be compared across regimes without being compared across physics
they do not share.

## At a glance

| | |
|---|---|
| **Inputs** | The regime clusters, `param_names`, `REGIME_SA_METRICS` |
| **Outputs** | Per-regime index tables; `compare_<index>_<metric>.csv` matrices; plots |
| **Code** | [`givendata_sensitivity_by_regime`](../api/sensitivity), `summarize_givendata`, `regime_comparison_matrix` |
| **Cost** | Hours on the largest cluster; roughly linear in the number of metrics |

## Concepts

**PAWN** measures the distance between the output distribution conditional on one
parameter and the unconditional one. **δ** is Borgonovo's moment-independent importance
measure. Neither decomposes variance, which is why both survive an output spanning ten
decades.

**The dummy parameter** is one the model never reads. Whatever index it earns is noise,
so it defines the floor.

Metrics analysed: `flux`, `permeability`, `theta`. The first two are analysed on
`log10`, so the density estimators are not dominated by a few huge values; `theta`,
bounded in [0, 1], stays linear.

## Procedure

1. Build a given-data problem from each cluster's sampled inputs.
2. For each regime, each metric in `REGIME_SA_METRICS`, and each parameter, compute PAWN
   and δ on the cluster.
3. Apply the `log10` transform to the metrics in `LOG_METRICS_DEFAULT`.
4. Rank per regime; `summarize_givendata` assembles the tables.
5. `regime_comparison_matrix` column-normalises into one row per parameter and one
   column per regime, which is what the heatmap and the parallel-coordinates plot read.
6. Plotters return their figure, or `None` for a plot they declined to draw — a skipped
   regime, an empty matrix — so wrapping them in `save_figure` cannot turn a skipped
   plot into a crash.

## Design choices

**Why two estimators rather than one.** They fail differently. Agreement between PAWN
and δ is evidence; disagreement is a flag to look at the cluster rather than a tie to
break.

**Why the regime-defining fractions are excluded.** `frac_surface`, `frac_oxide` and
`frac_metal` are degenerate within their own cluster, because the cluster was defined by
them. An index for them would measure the partition.

**Why `eta_surf` is not analysed by default**, despite being the sharpest surface metric
available — $\Phi_{app} = \Phi_{transport}\,\eta_{surf}$ exactly, and
$\Phi_{transport}$ is pressure-independent, so `eta_surf` carries the whole surface
effect in one dimensionless number. Analysing a fourth metric costs roughly 33% more
PAWN and δ work, and the largest cluster already runs for hours: the 773 K oxide preset
draws over 100k samples. That is too steep a standing price for something wanted
occasionally. It is still written to every scan CSV, so the data is there — add it to
the list, or analyse the column directly, when a run is specifically about surface
kinetics.

**Why column normalisation in the comparison matrix.** Raw indices are not comparable
across regimes, because the clusters differ in size and in output spread. Normalising
per column compares *rankings*, which is the question actually being asked.

## Checks and failure modes

| Symptom | Meaning |
|---|---|
| A parameter at or below the dummy floor | Indistinguishable from no effect **in this regime, under this preset**. Not a general statement. |
| A pinned parameter below the floor | Says nothing about its leverage — the preset removed its variance. Figure 20. |
| PAWN and δ disagree sharply | Look at the cluster before believing either. Often a small or skewed cluster. |
| Every parameter looks significant | The dummy floor has not been drawn. |
| An empty or skipped plot | The plotter returned `None`; the regime was skipped or the matrix was empty. Not a crash. |
| Indices shift between runs at the same `N` | Estimator noise. Compare against the dummy, not against the previous run. |
| The run takes hours | Expected. Do not add a fourth metric casually. |

---

## Where to go next

- The worked walkthrough, kept from the original how-to:
  {doc}`../how-to/sensitivity-analysis`
- Every warning in one place: {doc}`../caveats`
- What each output column means: {doc}`../appendices` C.1
