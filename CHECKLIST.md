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
- [x] **Fix 1 — unify `VALIDATION` with `CONDITIONS`.** `1e5` appears eight times across
      `P_ref`/`P_fixed*`, `P_down = 0` three times. Unify the *operating point*
      unconditionally; keep the deliberately wider *sweep grids* (`L2b.delta_ox_sweep`
      spans 1e-4→1e-12 m to show the thin-film limit) with a comment saying why.
- [x] **Fix 2 — the gate takes its operating point from `CONDITIONS`.** ⚠ re-baseline.
      `container/model_smoke_test.py` hardcodes `P_DOWN = 1.0e2` where `CONDITIONS` says
      **0**, plus `L_METAL`. `P_down = 0` is safe for the solver but shifts `J_ss` by
      3.3 %, so the references must be re-measured in the same commit.
- [x] **Fix 3 — `CONDITIONS` owns the oxide geometry.** ⚠ re-baseline, the main one.
      `CONDITIONS['L_oxide'] = 1e-6` is dead code; `OXIDES[...]['thickness'] = 4.8e-8`
      wins in both readers, so every model call runs a geometry 20× thinner than the
      config states. Switching to 1e-6 moves L5 by 0.4 % but **drops L5L6 flux 3× and
      flips its regime from metal-limited to oxide-limited**, so the gate references, the
      generated doc numbers and the "the oxide barely matters" prose all change with it.
- [ ] **Fix 3b — re-section the study file.** Five labelled sections in dependency
      order: identity, what is being studied (CONDITIONS, OXIDE_DEFECTS, MODEL_OPTIONS,
      ranges, presets, VALIDATION grids), what the materials are, derived, presentation.
      Gathers the five model switches — `include_gb_trapping`, `use_sieverts_pinhole`,
      `include_gb_enhancement`, `include_trapping`, `D_eff_method` — which are scattered
      across three dicts and are not data. A pure move: verified by diffing the module
      namespace before and after, which must be identical.
- [ ] **Fix 4 — `CONDITIONS` owns the swept operating ranges.** The range of conditions
      being studied *is* an operating condition, but it is stated twice and the two
      disagree — and each set has a nonsense bound: `P_range` reaches 1e12 Pa (10 million
      bar) and the SA's `oxide_thickness` floor is 1e-10 m (0.1 Å, sub-atomic).
      `CONDITIONS` takes defensible bounds — T `(573, 1273)`, P `(1e-7, 1e7)`, L_metal
      `(1e-4, 5e-3)`, L_oxide `(1e-9, 1e-5)` — and the four operating entries in
      `SUGGESTED_RANGES_*` reference them. Material-parameter ranges stay put. The
      divergence check narrows to its real job: reporting *preset* overrides that leave
      the window, which are deliberate and should be visible.
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
- [ ] **The preset yields quoted in `study_config.py` are stale.** The comment above
      `REGIME_PRESETS_L5` claims `metal 97.0%, oxide 77.7%, defect 87.3%`; measured today
      they are 34%, 78%, 100% (150 draws, seed 42). Unchanged by the geometry fix —
      verified by probing before and after — so the comment simply predates other work.
      `size_draws_for_target()` measures yields at run time, so nothing depends on the
      comment being right, but it reads as authoritative.
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
