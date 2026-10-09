# Stage 2 — The global scan

## Purpose

One space-filling pass over the whole parameter space, evaluated once with the full
record kept. This single pass feeds **both** the regime partition and the indices
computed on it, so no further model evaluations are needed downstream.

## At a glance

| | |
|---|---|
| **Inputs** | `param_ranges` (e.g. `SUGGESTED_RANGES_LEVEL5L6`), `N`, `seed=42` |
| **Outputs** | A tidy DataFrame, one row per draw; and the SALib problem dict for reuse |
| **Code** | [`run_global_lhs_scan`](../api/sensitivity), `level5_model_wrapper`, `level5L6_model_wrapper` |
| **Cost** | One model solve per draw; the dominant cost of the workflow |

## Concepts

**Latin hypercube sampling** divides each parameter's range into `N` equal-probability
strata and draws once from each, so the sample covers the space more evenly than
independent random draws at the same `N`.

**Given-data estimation** means the indices are computed from whatever sample exists,
rather than from a structured design built for a particular estimator. It is why one
scan can serve both purposes.

## Procedure

1. Resolve the wrapper — `level5L6_model_wrapper` by default, `level5_model_wrapper`
   for the no-surface branch.
2. Refuse any parameter that appears in both `fixed_params` and `param_ranges`; the two
   would silently contradict each other.
3. Build the SALib problem dict from the ranges, and draw `N` points with `latin.sample`
   under `seed`.
4. Call the wrapper once per draw with `return_full_record=True`.
5. Collect `SCAN_OUTPUT_FIELDS` as numeric columns and `SCAN_LABEL_FIELDS` as strings
   beside the sampled inputs. Columns are listed in {doc}`../appendices` C.1.
6. If `save_path` is given, write the CSV and stamp the folder with `study.txt`.

## Design choices

**Why one scan rather than two.** Partitioning and ranking were separable, and running
a design per regime first looked natural. But the regime is only knowable *after* the
model has run, so a design built per regime would have to guess the partition it was
about to discover. One pass that records everything avoids the circularity.

**Why the full record.** Keeping only the flux would have halved the storage and made
every later question — which path dominated, what the interface pressure did, whether
the solve converged — require a re-run.

**Why the study stamp.** Scans are reloaded from disk weeks later. Without a stamp, a
folder of CSVs is indistinguishable from another material system's, and the analysis
would silently proceed.

**Why `seed=42` is a default rather than a constant.** Reproducibility matters more
than the particular value, but a fixed seed across every run also hides sampling
variability, so it is an argument.

## Checks and failure modes

| Symptom | Meaning |
|---|---|
| A parameter in both `fixed_params` and `param_ranges` | Refused. The two would contradict each other. |
| Many `undefined` rows | Failed solves. Dropped at stage 3 and reported separately; a large count means the ranges reach unphysical combinations. |
| One regime holds almost every row | Expected on the default ranges. Not a fault — go to stage 4. |
| `flux` spanning ten decades | Expected. It is why the indices are computed on `log10` — `LOG_METRICS_DEFAULT`. |
| Silent extrapolation beyond the measured temperature window | The wrappers do not re-raise it. The sensitivity ranges extrapolate much further than a single run; {doc}`../caveats`. |
| Reloaded scans disagree with the current study | Check `study.txt` before trusting anything downstream. |
