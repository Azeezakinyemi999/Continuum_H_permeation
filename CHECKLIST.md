# Project checklist

Outstanding work, checked off in the commit that does it. Each item is meant to be one
commit: fix, test, commit, next.

Items are listed with enough context to act on without re-deriving the finding.

---

## Operating-point authority

`CONDITIONS` should define the scenario the simulation is built on — state and geometry —
and `METALS`/`OXIDES` should describe only what the material *is*. That rule is currently
enforced nowhere: nothing in `calculations/` reads `CONDITIONS`, so the operating point is
restated in four places with four different answers, and the material dict silently
overrides the geometry.

- [x] **Fix 0 — a way to regenerate the gate's reference values.** `--emit` mode on
      `container/model_smoke_test.py` that prints a ready-to-paste `REFERENCE[<study>]`
      block. Needed because two fixes below invalidate 13 pinned values per study and
      there is no generator — `container/revalidate_reference.py` is a cross-interpreter
      diff script, not this.
- [ ] **Fix 1 — unify `VALIDATION` with `CONDITIONS`.** `1e5` appears eight times across
      `P_ref`/`P_fixed*`, `P_down = 0` three times. Unify the *operating point*
      unconditionally; keep the deliberately wider *sweep grids* (`L2b.delta_ox_sweep`
      spans 1e-4→1e-12 m to show the thin-film limit) with a comment saying why.
- [ ] **Fix 2 — the gate takes its operating point from `CONDITIONS`.** ⚠ re-baseline.
      `container/model_smoke_test.py` hardcodes `P_DOWN = 1.0e2` where `CONDITIONS` says
      **0**, plus `L_METAL`. `P_down = 0` is safe for the solver but shifts `J_ss` by
      3.3 %, so the references must be re-measured in the same commit.
- [ ] **Fix 3 — `CONDITIONS` owns the oxide geometry.** ⚠ re-baseline, the main one.
      `CONDITIONS['L_oxide'] = 1e-6` is dead code; `OXIDES[...]['thickness'] = 4.8e-8`
      wins in both readers, so every model call runs a geometry 20× thinner than the
      config states. Switching to 1e-6 moves L5 by 0.4 % but **drops L5L6 flux 3× and
      flips its regime from metal-limited to oxide-limited**, so the gate references, the
      generated doc numbers and the "the oxide barely matters" prose all change with it.
- [ ] **Fix 4 — report sweep-range divergence.** `SUGGESTED_RANGES_*` restate the same
      quantities with different bounds (temperature `[573, 1273]` vs `T_range (623,
      1200)`, and three more), and `DEFAULT_SWEEP_TEMPERATURES` includes 1273 K. Report
      it in `check_against_config()`; do not enforce — a sweep may legitimately explore
      wider than the nominal window.
- [ ] **Fix 5 — name the provenance keys as provenance.** `METALS[...]['metal_thickness']`
      is a *list* of the source experiment's specimen thicknesses sitting where a reader
      looks for wall thickness.
- [ ] **Fix 6 — one shared `CONDITIONS`.** *(deferred, needs its own decision)* It is
      byte-identical in all four studies. Hoisting it would weaken "a study folder is
      self-contained", which was the point of an earlier change.

## Known defects, not yet scheduled

- [ ] **Four SA defaults lie outside their own sweep ranges** — `D_ref`,
      `k_diss_metal_ref`, `E_diss_metal`, in the 316L study only. The sensitivity analysis
      therefore never evaluates the nominal case for those parameters. Either the default
      or the range is wrong; `check_against_config()` has reported this for some time.
- [ ] **`REQUIRED_EXPORTS` is incomplete.** It omits `get_surface_kinetics_from_config`,
      which `calculations/surface_kinetics.py` imports at module level, plus
      `DEFAULT_SWEEP_TEMPERATURES`, `apply_style`, `R` and `F`. A study missing any of
      them passes the contract check and then fails with a raw `ImportError` — the exact
      failure the check exists to prevent.
- [ ] **`Proposal.ipynb` has an invalid escape sequence** (`\P`, almost certainly a
      `$\Phi$` label in a non-raw string). A `DeprecationWarning` today, a
      `SyntaxWarning` in a later Python.
- [ ] **Oxide properties are extrapolated ~500 K past their validated ceiling, silently.**
      The Cr₂O₃ data carry `temperature_range = [473, 773]` K; the SA sweeps to 1273 K.
      `get_oxide_properties_at_T` warns, but the level-5 wrappers — the ones the SA
      actually calls — compute the Arrhenius directly and never consult the range.

## Test coverage

- [ ] **Physics test suite.** The model wrappers are pinned only by the 13 values in the
      container gate. `level5_model_wrapper` / `level5L6_model_wrapper` are what the whole
      sensitivity analysis calls and nothing else covers them.
- [ ] **No test executes a notebook end to end.** The figure and study-resolution suites
      check notebook *source*; a run is verified only by hand.

## Repository hygiene

- [ ] **`docs/two-models` is unmerged**, 8 commits ahead of `main`, and edits two of the
      same documentation pages as the current branch — `docs/api/config.md` and
      `docs/how-to/switch-study.md`. Whichever merges second needs reconciling.
- [ ] **`container/README.md:157`** still describes the deleted `ACTIVE_STUDY` switch and
      `calculations/config/studies/`. It is accurate about the shipped 1.1.0 image and
      wrong about the repository, so it needs a "since 1.1.0" note rather than a rewrite.
