# Stage 3 — Partition by regime

## Purpose

Split the scan into clusters by which physical process carries most of the resistance,
so that every index computed afterwards is conditional on one physics rather than
averaged across several.

## At a glance

| | |
|---|---|
| **Inputs** | The scan DataFrame, `regime_col='regime'`, `min_cluster=300` |
| **Outputs** | `{regime: sub_df}`, plus a reported count of dropped rows |
| **Code** | [`partition_by_regime`](../api/sensitivity), `assign_regime`, `assign_regime_L5` |
| **Cost** | Negligible — no model evaluations |

## Concepts

A **regime** is the label for which step carries most of the resistance. It is assigned
from the resistance fractions, not from the flux:

| Branch | Regimes |
|---|---|
| L5L6 (surface kinetics) | `surface`, `oxide`, `metal` |
| L5 (no surface kinetics) | `oxide`, `metal`, `defect` |

The default rule is `argmax` — the largest fraction wins — with a `threshold` variant
that requires a fraction to exceed 0.5 before it is named, and reports `mixed`
otherwise.

```{figure} ../figures/fig-07-resistance-fractions.svg
:alt: Resistance fractions handing over smoothly as the oxide thickens, and the Example Wall's own split.
:width: 100%

**Figure 7.** Look at the 0.5 line in panel (A). The handover between regimes is
smooth, and the label puts a hard boundary across it — so a wall near the crossing will
flip label under a small parameter change without anything physical having changed much.
```

## Procedure

1. Count rows per regime and report the sizes as percentages.
2. Count `undefined` rows — failed solves — and drop them from the clusters by default,
   reporting the number separately rather than silently.
3. Group the remainder by regime label into `{regime: sub_df}`, re-indexing each.
4. Warn for any cluster below `min_cluster` (300) that its given-data δ and PAWN will be
   noisy, and suggest a targeted top-up.

## Design choices

**Why drop `undefined` rather than treat it as a fourth regime.** A failed solve is not
a physical state; pooling failures would make "what drives failure" look like a
sensitivity result when it is a numerical one. They are reported separately so the
count is visible.

**Why `min_cluster = 300` is a warning rather than an error.** Below 300 the indices
are noisy, not meaningless, and refusing to compute them would hide a run that is
nearly adequate. The warning names the cluster and the count so the decision is the
reader's.

**Why the fractions are excluded from the metrics.** `frac_surface`, `frac_oxide` and
`frac_metal` define the cluster, so within it they are degenerate: an index for them
would measure the partition, not the physics. See {doc}`../appendices` C.1.

**Why `argmax` rather than a threshold by default.** A threshold produces a `mixed`
class that is large near every handover, and `mixed` cannot be interpreted — it is not
a physics, it is the absence of a winner. `argmax` always names something, at the cost
of being unstable near a crossing, which Figure 7 makes visible.

## Checks and failure modes

| Symptom | Meaning |
|---|---|
| One cluster holds >90% of rows | The ranges do not reach the other regimes. Go to stage 4; this is the normal case. |
| A cluster below 300 | Indices will be noisy. Top it up with a targeted scan. |
| A cluster of zero | That regime is unreachable on these ranges. A preset is required, not a larger `N`. |
| Many `undefined` rows | Failed solves. Check whether the ranges reach unphysical combinations. |
| A label flipping under a tiny parameter change | Expected near a 0.5 crossing. Read the fractions, not only the label. |
| `mixed` dominating | The threshold rule is in force and most of the sample sits near a handover. |
