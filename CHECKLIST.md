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
- [x] **Fix 3b — re-section the study file.** Five labelled sections in dependency
      order: identity, what is being studied (CONDITIONS, OXIDE_DEFECTS, MODEL_OPTIONS,
      ranges, presets, VALIDATION grids), what the materials are, derived, presentation.
      Gathers the five model switches — `enable_gb_as_trap`, `use_sieverts_pinhole`,
      `enable_gb_enhancement`, `enable_trapping`, `D_eff_method` — which are scattered
      across three dicts and are not data. A pure move: verified by diffing the module
      namespace before and after, which must be identical.
- [x] **Fix 4 — `CONDITIONS` owns the swept operating ranges.** The range of conditions
      being studied *is* an operating condition, but it is stated twice and the two
      disagree — and each set has a nonsense bound: `P_range` reaches 1e12 Pa (10 million
      bar) and the SA's `oxide_thickness` floor is 1e-10 m (0.1 Å, sub-atomic).
      `CONDITIONS` takes defensible bounds — T `(573, 1273)`, P `(1e-7, 1e7)`, L_metal
      `(1e-4, 5e-3)`, L_oxide `(1e-9, 1e-5)` — and the four operating entries in
      `SUGGESTED_RANGES_*` reference them. Material-parameter ranges stay put. The
      divergence check narrows to its real job: reporting *preset* overrides that leave
      the window, which are deliberate and should be visible.
- [x] **Fix 5 — name the provenance keys as provenance.** `METALS[...]['metal_thickness']`
      is a *list* of the source experiment's specimen thicknesses sitting where a reader
      looks for wall thickness.
- [x] **Fix 7 — dissolve `VALIDATION`.** It was a fourth statement of the operating
      point wearing level prefixes. Gone: the entries that were conditions moved into
      `CONDITIONS`, the sweep grids are built in the notebooks from a `CONDITIONS` range,
      eight point-count names collapsed to `n_points`, and eight dead or discarded keys
      were deleted. `REQUIRED_EXPORTS` no longer asks a study for it.
- [x] **Fix 6 — one shared `CONDITIONS`. REJECTED, 2026-09-25.** Each study is an
      independent experiment; Guo 2025 has nothing to do with Fuerst or Incoloy 802.
      The four identical blocks are therefore not duplication to be removed — they are
      four scenarios that happen to agree today, and nothing should force them to keep
      agreeing. No shared definition, and deliberately no test asserting they match.

## Known defects, not yet scheduled

- [ ] **Four SA defaults lie outside their own sweep ranges** — `D_ref`,
      `k_diss_metal_ref`, `E_diss_metal`, in the 316L study only. The sensitivity analysis
      therefore never evaluates the nominal case for those parameters. Either the default
      or the range is wrong; `check_against_config()` has reported this for some time.
- [x] **`REQUIRED_EXPORTS` was incomplete — five names, not one.** It omitted
      `get_surface_kinetics_from_config` (imported at module level by
      `surface_kinetics.py`), plus `R`, `F`, `apply_style` and
      `DEFAULT_SWEEP_TEMPERATURES`, all imported only by notebooks. A study missing any
      of them passed the contract and failed later anyway. The set is now derived rather
      than maintained: a test re-scans every import of `model_config` — across
      `calculations/`, `container/`, `docs/_tools/` and all 24 notebooks — and fails if
      the contract is short, or if it demands a name nobody imports.
- [x] **The notebook tests covered 8 of 24 files.** `test_operating_conditions.py`
      globbed `*/[A-Z]*.ipynb`, which matches only `Proposal` and `Surface_proposal`;
      `sensitivity_regime_*` and `regime_parallel_coords*` start lowercase and were never
      scanned. That is how the contract gap survived. Widened to `[A-Za-z]*`; the four
      previously-unscanned notebooks turned out clean.
- [ ] **`Proposal.ipynb` has an invalid escape sequence** (`\P`, almost certainly a
      `$\Phi$` label in a non-raw string). A `DeprecationWarning` today, a
      `SyntaxWarning` in a later Python.
- [ ] **Every study runs at 873 K because it was copied, not chosen — and for two of
      them that is outside their own metal's validated window, downward.** Guo (316L)
      operates 50 K above its metal's [468, 823] K ceiling; Fuerst (Hastelloy N) and
      Incoloy 802 operate **150 K below** their metals' [1023, 1223] K floor, and their
      `T_range` starts 450 K below it. Extrapolating an Arrhenius fit downward is the
      least forgiving direction. Now that each study is explicitly independent, each
      needs an operating point chosen for its own material — or an explicit note that
      the extrapolation is intended.
- [x] **The L6 activation-energy extraction was contaminated by borderline points.**
      Diagnosed and fixed: the Arrhenius fit accepted any point classified "pure
      surface" at a 90% threshold, but the remaining 10% belongs to a path whose
      activation energy is ~234 kJ/mol against the surface's ~58, biasing the fit upward
      by ~10 kJ/mol. Raising the threshold to 0.99 returns 59.2 against a configured
      57.95, and repairs the oxide fit too (217 -> 231 against 234). Ticks per
      Surface_proposal run: 8 -> 12. The threshold was always too loose; widening
      `T_range` merely produced enough borderline points to expose it. One panel still
      warns in each study (85.1 for 316L, and two for Incoloy 802, which runs 150 K
      below its metal's validated window) — a genuinely mixed regime, not a fit artefact.
- [ ] **The L5 sensitivity panel now runs at 873 K, outside the oxide's validated
      window.** It used to use 673 K, chosen to sit inside [473, 773] K. Deliberate, per
      the decision to reference `T_operating` directly; recorded so the count of places
      that extrapolate stays honest.
- [x] **The measurement-validity windows are named inconsistently.** The metal carried
      `temp_range` / `pressure_range`, the oxide `temperature_range`. Both mean the same
      thing — where the source data are valid — and only the oxide's was read, by
      `get_oxide_properties_at_T`, to warn on extrapolation. The metal's was read by
      nothing, so a metal extrapolation was silent.
      **Fixed:** the metal's window is now `temperature_range` in all four studies, and
      both readers call one shared `_warn_if_extrapolating`. The warning this exposed is
      not cosmetic — every study now warns at its *own* operating point, and two of them
      run 150 K **below** their metal's measured floor, the least forgiving direction for
      an Arrhenius extrapolation:

      | study | metal window | runs at |
      |---|---|---|
      | Guo_etal_2025 | [468, 823] K | 873 K — 50 K above |
      | Fuerst_et_al_2024 | [1023, 1223] K | 873 K — 150 K below |
      | incoloy802_cr2o3 | [1023, 1223] K | 873 K — 150 K below |

      Changes no number: the gate passes unchanged on all four studies and
      `worked_values --check` is clean. `pressure_range` is left alone — nothing reads it
      for either material, so renaming it would be churn without a consumer.
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

- [x] **Physics test suite.** The model wrappers are pinned only by the 13 values in the
      container gate. `level5_model_wrapper` / `level5L6_model_wrapper` are what the whole
      sensitivity analysis calls and nothing else covers them.
      **Done:** `tests/test_wrapper_physics.py`, 20 tests stating invariants rather than
      values — flux positive and finite; monotonic in pressure, both thicknesses and
      temperature; `frac_*` partitioning to 1; `regime` equal to the argmax of those
      fractions; `PRF >= 1`; sealing the defects never raising the flux; and Sieverts'
      half-power in pressure. Most apply to both wrappers, so they hold for any study.

      Deliberately absent: any assertion that L5L6 collapses onto L5 as the surface
      becomes fast. They are different models, not refinements — Henry vs Sieverts in
      the oxide — so they do not converge, and asserting it would pin a bug in place.

      Two findings from building it, both from mutation-checking rather than from the
      tests passing. The Sieverts test originally *skipped* itself unless the nominal
      point was metal-limited, and every mutation that broke the square root also moved
      the regime — so the guard skipped the test instead of failing it. It now forces a
      bare-metal wall (`oxide_thickness=1e-9`, defects sealed) and asserts the log-log
      slope is 0.5, which does catch a dropped root in `permeation_calc.py`. And only
      one of the three `np.sqrt` sites is on the L5 path at all; the two in
      `interface_solver.py` can be broken without any wrapper test noticing.
- [ ] **`level5L6_model_wrapper` returns `nan` at zero driving force.** With
      `P_upstream == P_downstream` Level 5 correctly gives 0.0 and Level 5L6 gives `nan`
      — the surface solve cannot bracket a root when there is nothing to drive it. Found
      while writing the physics suite, where it is recorded as a strict `xfail`, so
      fixing it will announce itself by turning that test into an unexpected pass.
- [x] **No test executes a notebook end to end.** The figure and study-resolution suites
      check notebook *source*; a run is verified only by hand.
      **Done:** `tests/test_notebook_execution.py` runs `Guo_etal_2025/Proposal.ipynb`
      through `nbconvert` in its own folder, asserts no cell produced an `error` output,
      and asserts all six figures were written *by that run* rather than left over from
      an earlier one. ~60 s, marked `slow`, deselectable with `-m "not slow"` but on by
      default. The executed copy goes to `tmp_path`, so running the tests never rewrites
      the committed notebook.

      One notebook, not 24, and deliberately: the `sensitivity_regime_*` notebooks run
      full sampling scans, and the `regime_parallel_coords*` ones only read scan CSVs
      that are gitignored and cannot be regenerated without running the scans first —
      they would have to be skipped on a fresh clone.

      Mutation-checked in both directions: a cell that raises is caught by the error
      assertion, and a `save_figure` turned into a no-op — which leaves the notebook
      running perfectly cleanly — is caught by the figure count.
- [ ] **The notebooks declare kernel `"name": "python3"`, which is not necessarily the
      interpreter that has the package.** On this machine that kernelspec points at a
      different conda environment and the first cell dies on `import calculations`; the
      environment with the editable install is registered as `mace_env`. Nothing in the
      repository pins the association, so whether a notebook runs depends on what the
      reader happens to have registered under the name `python3`. The new execution test
      works around it by finding the kernelspec whose `argv[0]` matches `sys.executable`,
      which is a fix for the test and not for the notebooks.

## Repository hygiene

- [ ] **The Guo SA notebooks' stored outputs predate the geometry re-baseline.** Merging
      `fix/eta-surf-not-default` brought re-executed copies of `sensitivity_regime_L5`,
      `sensitivity_regime_L5L6` and both `regime_parallel_coords*`, but that branch forked
      at `0c264dd`, before `acfe540` made `CONDITIONS` own the oxide geometry — so they
      were computed on the 48 nm layer, where L5L6 flux is 3x higher and the regime reads
      metal rather than oxide. Its outputs were therefore dropped and only its
      `sensitivity.py` change taken; the notebook *source* is identical either way.
      The committed outputs are now stale in both copies and a correct set needs a real
      re-execution, which is hours per notebook — the L5L6 773 K oxide preset alone draws
      >100k samples.

- [x] **`docs/two-models` is unmerged**, 8 commits ahead of `main`, and edits two of the
      same documentation pages as the current branch — `docs/api/config.md` and
      `docs/how-to/switch-study.md`. Whichever merges second needs reconciling.
      **Merged.** It was more than docs: it brought `calculations/permeability.py`,
      `tests/test_permeability.py` (35 tests, 290 -> 325) and the four-tier permeability
      metric. Six conflicts, all from the branch predating this one:

      | file | resolution |
      |---|---|
      | `calculations/config/model_config.py` | ours — the branch only flipped the deleted `ACTIVE_STUDY` literal |
      | `docs/getting-started.md` | hand-merged; `--theirs` would have reverted the page to the switch era |
      | `docs/how-to/sensitivity-analysis.md` | kept our figure-folder paragraph, took their warning relocation |
      | `docs/how-to/switch-study.md`, `docs/theory/equilibrium-models.md` | theirs — warnings relocated to `caveats.md` |
      | `docs/_tools/worked_values.txt` | regenerated rather than hand-merged |

      What the merge actually fixes: `permeability` was a harmonic mean of two bulk
      constants, so L5 and L5L6 both reported `4.000079e-15` despite fluxes differing
      3.5x. They now read 1.008e-11 and 2.878e-12, each backed out of its own flux. The
      old number understated by **775x**.

      The two new pages arrived stale, written against the pre-restructuring config, and
      were corrected against the running code rather than merged as-is: `caveats.md` had
      a 48 nm wall, `frac_oxide = 3.8e-4` and `PRF ~ 1.0002` (now 1 um, 7.8e-3, 1.004),
      and `config-file.md` documented `calculations/config/studies/`, a `VALIDATION`
      block, `metal_thickness` and a fifteen-name export contract — all retired.

      Note `verify_docs.py` passed on that stale page: it checks that referenced config
      *names* resolve, so it never noticed `VALIDATION` being documented after deletion.
      A prose claim about a name that no longer exists is outside what it can see.
- [x] **`container/README.md:157`** still describes the deleted `ACTIVE_STUDY` switch and
      `calculations/config/studies/`. It is accurate about the shipped 1.1.0 image and
      wrong about the repository, so it needs a "since 1.1.0" note rather than a rewrite.
      **Done that way.** The changelog entry stays true to the image it documents and is
      marked superseded; a new "Since 1.1.0" section states what a reader comparing a
      1.1.0 result against a current checkout needs to know — chiefly that the numbers
      moved, because `CONDITIONS` now owns the wall (1 um oxide, not the 48 nm coupon)
      and `permeability` is apparent rather than a bulk harmonic mean.
      `test_container_readme_mentions_the_switch_only_as_history` keeps the name in
      changelog territory: it may be described as history, never as instruction.
