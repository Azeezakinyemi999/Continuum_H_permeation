# Regime-Stratified Sensitivity Analysis — How It Works

A teaching guide to the regime-stratified sensitivity analysis (SA) for the **L5L6**
permeation model.

- **Engine:** [`validation/sensitivity_level1.py`](../validation/sensitivity_level1.py)
- **Analysis notebook:** [`Application/sensitivity_regime_L5L6.ipynb`](../Application/sensitivity_regime_L5L6.ipynb)
- **Visualization notebook:** [`Application/regime_parallel_coords.ipynb`](../Application/regime_parallel_coords.ipynb)
  (parallel-coordinates plots — see §6)
- **Run with:** the `mace_env` conda kernel
  (`/Users/akinyemi.az/anaconda3/envs/mace_env/bin/python`). The base anaconda env is broken
  (numpy 2.3.4 vs numpy-1.x-compiled deps).
- **Parameters varied: 36.** The four reference temperatures (`T_ref_*`) are held **fixed** — see the scope note below.

---

## 0. The one-sentence idea

> **Sensitivity analysis** asks *"if I wiggle each input parameter, how much does the output
> (flux) wiggle?"* — i.e. **which parameters matter most.** We added a twist: **the answer
> depends on which physical regime the system is in**, so we compute the sensitivity
> *separately for each regime.*

### The basic analogy

The hydrogen barrier has three steps in series, like three resistors in a wire:

```
   gas ──[ R_surface ]──[ R_oxide ]──[ R_metal ]──> downstream
          dissociation    oxide        metal
                          diffusion    diffusion
```

The **biggest resistor is the bottleneck**. Which one is biggest = the **regime**:

- `R_surface` biggest → **surface-limited**
- `R_oxide` biggest → **oxide-limited**
- `R_metal` biggest → **metal-limited**

The parameters that matter for flux are *different in each regime* (surface-kinetics params in
the surface regime, metal-transport params in the metal regime, …). A normal SA blends all
regimes and hides this. **We separate them.**

> **Scope note — 36 parameters (T_ref fixed).** In the model's `X(T)=X_ref·exp(-E/R·(1/T-1/T_ref))`,
> the reference temperature `T_ref` is *measurement metadata* (the temperature the reference value was
> taken at) and is redundant with `X_ref` via the prefactor `A=X_ref·exp(E/(R·T_ref))`. Varying it
> independently injects an artificial ~10³–10⁴ prefactor sweep, so the four `T_ref_*` are held **fixed**
> and the SA varies **36** parameters. Property uncertainty is captured through the `*_ref` and
> activation-energy (`E_*`) ranges instead. (The operating `temperature`, 573–1273 K, *is* varied.)

---

## 1. The ENTIRE workflow

```mermaid
flowchart TD
    M["L5L6 model<br/>calculate_full_model_flux_L346_v2"] --> S0
    S0["STAGE 0 — label by regime<br/>argmax(frac_surface, frac_oxide, frac_metal)"] --> S1
    S1["STAGE 1 — targeted sampling<br/>one LHS scan per regime preset<br/>keep in-regime rows"] --> RB
    S1 --> RA
    S1 -. feeds top-k .-> SCREEN["Global Morris screen<br/>rank → top-k"]
    SCREEN -.-> RA
    RB["STAGE 2 Route B (PRIMARY)<br/>given-data PAWN + delta<br/>on the TRUE cluster"] --> S3
    RA["STAGE 2 Route A (secondary)<br/>Sobol over top-k sub-box<br/>accept-and-report"] --> S3
    S3["STAGE 3 — compare regimes<br/>heatmaps + contamination table"]
```

Same thing as a plain-text schema:

```
        ┌────────────────────────────────────────────────────────────────┐
        │                THE L5L6 PERMEATION MODEL                         │
        │      calculate_full_model_flux_L346_v2(params) → flux, fractions │
        └────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
   ╔══════════════════════════════════════════════════════════════════════╗
   ║ STAGE 0 — LABEL each model run with its REGIME                         ║
   ║   regime = argmax(frac_surface, frac_oxide, frac_metal)                ║
   ╚══════════════════════════════════════════════════════════════════════╝
                                       │
                                       ▼
   ╔══════════════════════════════════════════════════════════════════════╗
   ║ STAGE 1 — TARGETED SAMPLING  (because regimes are rare!)               ║
   ║   metal-preset   ─LHS─► run model ─► keep rows labeled "metal"         ║
   ║   surface-preset ─LHS─► run model ─► keep rows labeled "surface"       ║
   ║   oxide-preset   ─LHS─► run model ─► keep rows labeled "oxide"         ║
   ║   Output: 3 "clusters" = piles of (inputs, outputs) per regime         ║
   ╚══════════════════════════════════════════════════════════════════════╝
              │                                              │
              │  (also, once, on the FULL space:)            │
              │  ┌────────────────────────────────────┐     │
              │  │ GLOBAL MORRIS SCREEN → rank → top-k │     │
              │  │ (cheap "which knobs matter at all") │     │
              │  └────────────────────────────────────┘     │
              ▼                                              ▼
   ╔═══════════════════════════════╗      ╔════════════════════════════════╗
   ║ STAGE 2 — ROUTE B  (PRIMARY)  ║      ║ STAGE 2 — ROUTE A  (SECONDARY) ║
   ║ given-data SA on TRUE cluster ║      ║ Sobol over top-k SUB-BOX       ║
   ║   • PAWN   (CDF shift)        ║      ║   • S1 / ST variance shares    ║
   ║   • Borgonovo δ (density)     ║      ║   • "accept & report" leakage  ║
   ║ no box, no contamination      ║      ║ structured, but box leaks      ║
   ╚═══════════════════════════════╝      ╚════════════════════════════════╝
              │                                              │
              └──────────────────────┬───────────────────────┘
                                     ▼
   ╔══════════════════════════════════════════════════════════════════════╗
   ║ STAGE 3 — COMPARE regimes (parameter × regime heatmaps, tables)        ║
   ╚══════════════════════════════════════════════════════════════════════╝
```

**Plain walk-through:**

1. **Stage 0** — every model run is tagged with its regime (which resistor won).
2. **Stage 1** — random sampling almost never produces surface-limited cases (0 in 1500 random
   runs), so we **aim** the sampling at each regime with a "preset" (tuned input ranges that
   produce that regime), then keep the runs that actually landed in it → three clean data piles.
3. **Stage 2** — for each pile we ask "which parameters matter?" two ways:
   - **Route B** uses methods (PAWN, δ) that work on *any* pile of points → the honest, primary answer.
   - **Route A** uses classic Sobol, which needs a tidy rectangular box → we box the regime,
     accept that the box leaks into neighbours, and **report** how much it leaked.
4. **Stage 3** — line the regimes up and see how the important parameters differ.

---

## 2. Each stage as its own mini-workflow

### STAGE 0 — Regime labeling

```
 params ─► calculate_full_model_flux_L346_v2 ─► flux_weighted_resistances
                                                  │  fraction_surface
                                                  │  fraction_oxide
                                                  │  fraction_metal   (sum = 1.0)
                                                  ▼
                                   assign_regime(fs, fo, fm)
                                                  │
                            ┌─────────────────────┼─────────────────────┐
                          fs biggest           fo biggest            fm biggest
                            ▼                     ▼                     ▼
                        "surface"              "oxide"               "metal"
                            (NaN inputs → "undefined", thrown away)
```

**Why:** turns a continuous result into one clean label so we can group runs. `argmax` = "pick the
biggest of the three" (no fuzzy "mixed" bucket). Code: `assign_regime()`, attached inside
`level5L6_model_wrapper(..., return_full_record=True)`.

---

### STAGE 1 — Targeted sampling (per regime)

```
 for each regime R in {metal, surface, oxide}:

   REGIME_PRESETS[R]            ← input ranges tuned to PRODUCE regime R
        │
        ▼
   latin.sample (LHS)           ← space-filling draws over those ranges
        │   N draws
        ▼
   run model on each draw       ← run_global_lhs_scan(...)
        │   (tag each with regime via Stage 0)
        ▼
   keep rows where regime == R  ← the "cluster" for R
        │
        ▼
   warn if cluster < ~300-500   ← need enough points for stable stats
```

**Why a preset per regime?** A plain random sweep gave **metal 97%, oxide 3%, surface 0%** —
surface-limited is so rare you'd never get a usable pile. Each preset nudges inputs toward the
corner that produces that regime (surface = *low pressure + slow surface kinetics*; oxide =
*suppressed defects + slow/thick oxide + fast surface & metal*), while leaving every other
parameter at full range so the sensitivity result stays honest. Code: `run_targeted_regime_scans()`.
The production run (36-param) gave in-regime clusters of **metal 813, surface 714, oxide 3725** — all
comfortably above the ~300–500 floor.

> **LHS (Latin Hypercube Sampling):** a smart way to scatter sample points so they evenly cover
> every input's range (better than pure random). *Source: McKay, Beckman & Conover (1979),
> Technometrics 21(2):239–245.*

---

### GLOBAL MORRIS SCREEN (side-branch that feeds top-k)

```
 full 36-param ranges ─► Morris trajectory sample ─► run model (all metrics, ONE pass)
                                                          │
                                                          ▼
                             for each metric: μ* (mean |elementary effect|)
                                                          │
                                            normalize per metric → average
                                                          ▼
                                              rank → pick TOP-k parameters
```

**Why:** Morris is a **cheap first filter** — "of 36 knobs, which ~10 do anything at all?" Those
top-k define the Route A box (so Sobol isn't wasted on 36 dimensions). Code:
`morris_screen_global()` + `rank_and_select_top_k()`.

---

### STAGE 2 — ROUTE B (the primary, trustworthy answer)

```
 cluster_R (true pile of in-regime points)
        │
        │  X = the input columns,  Y = an output (flux → log10, theta linear)
        ▼
   ┌──────────────────────────────┬───────────────────────────────┐
   │ PAWN.analyze(X, Y)            │ delta.analyze(X, Y)           │
   │ → "median" KS distance        │ → "delta" (∈[0,1]) + S1        │
   │   per parameter               │   per parameter                │
   └──────────────────────────────┴───────────────────────────────┘
        │
        ▼
   rank parameters per regime → top drivers of flux in THAT regime
```

**Why these two methods?** Classic Morris/Sobol need a *special structured sample*; you cannot feed
them an arbitrary filtered pile. **PAWN and δ are "given-data" methods — they work on any cloud of
(input, output) points.** That's exactly what an in-regime cluster is. No box, no leakage. Code:
`givendata_sensitivity_by_regime()`.

---

### STAGE 2 — ROUTE A (the comparison answer)

```
 cluster_R + top-k params
        │
        ▼
 extract_topk_subbox            ← rectangular box = 5th–95th percentile of the
        │                          cluster on each top-k parameter
        ▼
 Saltelli sample inside the box ← (other params held at the cluster median)
        │
        ▼
 run model on each ──► S1, ST   ← Sobol variance shares
        │
        ▼
 re-label each sample's regime  ← how many fell OUTSIDE regime R?
        │
        ▼
 report "contamination %"       ← keep all points (deleting them would break
                                   Sobol's math) and just report the leakage
```

**Why "accept and report"?** Sobol's estimator needs its rows kept in a strict paired arrangement;
deleting the leaked points breaks the math. So we run it over the box, get standard S1/ST, and
**report** that e.g. 35 % of the surface box leaked. That's why Route A is the *secondary* check,
not the primary. Code: `sobol_regime_subbox()`.

---

### STAGE 3 — Compare

```
 Route B results (all regimes, all metrics)
        │
        ▼
 regime_comparison_matrix  →  parameter × regime grid of δ (or PAWN)
        │                       (union of each regime's top drivers)
        ▼
 plot_regime_comparison_heatmap   +   contamination_summary table
```

**Why:** the payoff — one picture showing *"temperature matters in all regimes, but `f_crack` only
matters in oxide, `k_diss_metal` only in surface, `gb_diffusivity` only in metal."* Code:
`regime_comparison_matrix()` / `plot_regime_comparison_heatmap()`.

---

## 3. The four SA methods in plain words

| Method | The question it answers | How | Needs special sample? |
|---|---|---|---|
| **Morris** (screening) | "Which knobs do *anything*?" | Nudge one input at a time, measure output change (μ\* = avg size, σ = nonlinearity) | Yes (trajectories) |
| **Sobol** (variance) | "What % of output *variance* does each input cause?" (S1 alone, ST with interactions) | Special A/B matrices | Yes (Saltelli) |
| **PAWN** (distribution) | "How much does the output's *distribution shape* change when I pin this input?" | Kolmogorov–Smirnov distance between conditioned vs free output CDFs | **No** (any data) |
| **Borgonovo δ** (density) | "How much does the output's *probability density* shift when I pin this input?" δ∈[0,1] | Area between conditioned vs unconditioned densities | **No** (any data) |

> **Why we log₁₀ the flux first:** flux spans ~10 orders of magnitude. Variance- and density-based
> methods get hijacked by a few giant values (δ came out as 46 and −93389 — impossible, δ must be in
> [0,1]). Taking log₁₀ fixes it: we measure what controls the *order of magnitude* of flux, which is
> what you actually care about. `theta` (already 0–1) stays linear. Controlled by `LOG_METRICS_DEFAULT`.

---

## 4. Stage → code map

| Stage | Functions in `validation/sensitivity_level1.py` |
|---|---|
| 0 — regime label | `assign_regime`, `level5L6_model_wrapper(return_full_record=True)` |
| 1 — targeted scans | `REGIME_PRESETS`, `DEFAULT_N_PER_REGIME`, `run_global_lhs_scan`, `run_targeted_regime_scans`, `load_regime_scans`, `partition_by_regime`, `plot_regime_exploration` |
| 1 — Morris screen | `morris_screen_global`, `rank_and_select_top_k` |
| 2 — Route B | `givendata_sensitivity_by_regime`, `summarize_givendata`, `plot_givendata_results` |
| 2 — Route A | `extract_topk_subbox`, `sobol_regime_subbox`, `estimate_subbox_purity` |
| 3 — compare | `regime_comparison_matrix`, `plot_regime_comparison_heatmap`, `contamination_summary`, `summarize_route_a` |
| 6 — parallel coords | `parallel_coordinates_samples`, `parallel_coordinates_sensitivity`, `top_drivers` (notebook `regime_parallel_coords.ipynb`) |
| helpers | `_make_problem`, `_prep_Y_full`, `_givendata_problem`, `_save_scan`, `_pcp_axis` |

---

## 5. Saved data & caching (avoid re-running)

Everything is written to CSV under **`Application/sa_results/`** (set by `RESULTS_DIR` in the
notebook config cell). The notebook **caches the model evaluations**: on the first run it computes
and saves; on later runs it **reloads the CSVs and skips the expensive model calls** (the cheap
SA analysis is recomputed). To recompute from scratch, set `FORCE_RECOMPUTE = True` or delete the
folder.

| File | What it holds | Reused to skip… |
| --- | --- | --- |
| `scans/scan_metal.csv`, `scan_surface.csv`, `scan_oxide.csv` | every LHS draw per regime (36 params + outputs + regime) | the targeted scans (biggest cost) |
| `master_clusters.csv` | all in-regime rows combined | — (convenience) |
| `morris_evals.csv` | Morris design X + Y per metric | the global Morris model runs |
| `morris_ranking.csv` | μ\* ranking + top-k | — (result table) |
| `routeB_givendata.csv` | tidy PAWN median, δ, given-data S1 per regime/metric/param | — (Route B is already eval-free) |
| `routeA_evals_<regime>.csv` | Saltelli X + Y + regime_hit per regime | the Route A model runs |
| `routeA_sobol.csv` | tidy S1/ST/conf per regime/metric/param | — (result table) |
| `routeA_contamination.csv` | off-regime % per regime | — (result table) |
| `compare_<index>_<metric>.csv` | parameter × regime comparison matrices | — (result table) |

**Caching is keyed by file existence + sample size/seed**, implemented via the `cache_csv`
argument of `morris_screen_global` / `sobol_regime_subbox` and `load_regime_scans` for the scans.
If you change `N_PER_REGIME`, `MORRIS_TRAJ`, `SOBOL_N`, or `SEED`, delete `sa_results/`
(or set `FORCE_RECOMPUTE = True`) so the caches are rebuilt to match.

---

## 6. Parallel-coordinates plots (Plotly)

Interactive views of the results, built **from the saved CSVs** (no recompute) in the standalone
notebook [`Application/regime_parallel_coords.ipynb`](../Application/regime_parallel_coords.ipynb).
Each is written to `sa_results/figures/*.html` and also renders inline. Functions:
`parallel_coordinates_samples`, `parallel_coordinates_sensitivity`, `top_drivers`.

| View | Shows | Source CSV | Axes / colour |
| --- | --- | --- | --- |
| **A** `pcp_A_regime.html` | regime manifold in parameter space | `master_clusters.csv` | union of regimes' top-δ drivers (+ log10 flux); colour = **regime** |
| **B** `pcp_B_<regime>.html` | within-regime structure (one per regime) | `master_clusters.csv` (one regime) | that regime's top-6 drivers; colour = **log10(flux)** |
| **C** `pcp_C_sensitivity.html` | how each parameter's δ shifts across regimes | `compare_delta_flux.csv` | leftmost = **parameter name**, then metal/oxide/surface; colour = mean importance |
| **D** `pcp_D_outputs.html` | regime separation in output space | `master_clusters.csv` | flux, theta, frac_\*, permeability; colour = **regime** |

Notes:
- Heavy-tailed columns (P, k_diss, thicknesses, flux, …) are shown on **log10** axes; regimes are
  coloured discretely (surface/oxide/metal).
- **Identifying lines:** `go.Parcoords` has no per-line legend or hover, so view **C** adds a leftmost
  categorical **parameter** axis (sorted by importance, highest on top) — every line traces back to its
  parameter. `parallel_coordinates_sensitivity(..., top_n=N)` trims to the N most important.
- **Layout:** plots with many axes (view A, 11) get a pinned wide width so titles don't overlap;
  few-axis plots stay responsive (fill the cell). Tunable via `width`/`labelangle` args.
- Interactivity: drag along any axis to **brush/filter**. Static PNG export needs `pip install kaleido`.

---

## 7. Sources

**Borgonovo δ (moment-independent importance):**

- Borgonovo, E. (2007). *A new uncertainty importance measure.* **Reliability Engineering & System
  Safety**, 92(6), 771–784. doi:10.1016/j.ress.2006.04.015 — the original δ measure.
- Plischke, E., Borgonovo, E., & Smith, C. L. (2013). *Global sensitivity measures from given data.*
  **European Journal of Operational Research**, 226(3), 536–550. doi:10.1016/j.ejor.2012.11.047 —
  the "given-data" estimator that **SALib actually implements**.

**PAWN (CDF / distribution-based):**

- Pianosi, F., & Wagener, T. (2015). *A simple and efficient method for global sensitivity analysis
  based on cumulative distribution functions.* **Environmental Modelling & Software**, 67, 1–11.
  doi:10.1016/j.envsoft.2015.01.004 — the original PAWN method.
- Pianosi, F., & Wagener, T. (2018). *Distribution-based sensitivity analysis from a generic
  input–output sample.* **Environmental Modelling & Software**, 108, 197–207.
  doi:10.1016/j.envsoft.2018.07.019 — the given-data version SALib uses.

**Supporting methods also used:**

- Morris, M. D. (1991). *Factorial sampling plans for preliminary computational experiments.*
  **Technometrics**, 33(2), 161–174. (+ Campolongo, Cariboni & Saltelli, 2007, EMS
  22(10):1509–1518, for the μ\* improvement.)
- Sobol', I. M. (2001). *Global sensitivity indices for nonlinear mathematical models and their
  Monte Carlo estimates.* **Mathematics and Computers in Simulation**, 55(1–3), 271–280.
  (+ Saltelli et al., 2010, Computer Physics Communications 181(2):259–270, for the sampling scheme.)
- McKay, M. D., Beckman, R. J., & Conover, W. J. (1979). **Technometrics**, 21(2), 239–245 —
  Latin Hypercube Sampling.
- Herman, J., & Usher, W. (2017). *SALib: an open-source Python library for sensitivity analysis.*
  **Journal of Open Source Software**, 2(9), 97. doi:10.21105/joss.00097 — the library implementing
  all of the above.
