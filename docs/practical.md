# Practical

What to do when something fails, how to point the model at a new material system, what
has to be true for a result to be reproducible, and where each piece of code lives.

---

## 1. Troubleshooting

Symptoms are listed as they actually appear. Most of these are the model behaving as
defined rather than breaking; the fix is usually to change what you asked for.

### 1.1 Starting up

| Symptom | Cause | Fix |
|---|---|---|
| `RuntimeError` listing available studies, on any import | No study is in scope. The repository root belongs to no study, and there is no default. | Run from a study folder, or set `PERMEATION_STUDY=Guo_etal_2025/study_config.py`. |
| `ModuleNotFoundError: No module named 'numpy'` | Wrong interpreter — the system Python, not the project environment. | `conda activate mace_env`. |
| `ImportError` on `from calculations... import ...` outside the repository root | The package is not installed in the environment. | `pip install -e .` from the repository root. |
| `SALib` import fails on `np.trapezoid` | `SALib >= 1.5.1` requires NumPy 2, but torch pins NumPy below 2. | Pin `SALib==1.5.0`. |
| A notebook writes into the wrong study's results folder | The kernel started somewhere other than the notebook's own folder. | Notebooks call `assert_study_is_local()` in the first cell; restart the kernel from the study folder. |

### 1.2 Results that look wrong but are not

| Symptom | Cause | Fix |
|---|---|---|
| The coating appears to do nothing — `PRF ≈ 1`, `frac_oxide ≈ 0.008` | At Level 5 this wall genuinely is metal-limited. | Expected. Under the L6 family the same wall is oxide-limited. See {doc}`theory/two-models`, and do not carry a conclusion across the fork. |
| `PRF` is `nan` from the Level 5L6 wrapper but finite from Level 5 | Known, recorded in {doc}`caveats`. | Use the Level 5 value, or compute the ratio you actually want from the fluxes. |
| Two permeabilities that should match, differ | They are different **tiers**. There is no single system permeability under Model 1. | Read {doc}`theory/permeability`; compare like with like, Figure 15. |
| `UnitMismatchError` from `combine_series` | Model 1's oxide and metal permeabilities are not commensurable. | This is the physics refusing, not a bug. Figure 5. |
| Multiplying the oxide-defect correction by the metal correction does not give Level 5 | The two corrections are coupled through $P_{int}$. | Expected — Figure 13 quantifies the error. Solve the full level. |
| A tiny change in a parameter flips the regime label | The label puts a hard boundary on a smooth curve. | Expected near the 0.5 line — Figure 7A. Read the fractions, not only the label. |
| Flux changes when you switch studies, by a lot | Each study is a different material system. | Expected. `METALS` keys are disjoint across studies, so nothing carries over. |

### 1.3 Numerical and solver

| Symptom | Cause | Fix |
|---|---|---|
| `converged: False` with $P_{int} \approx P_{down}$ | $P_{up}$ is at or below `min_pressure` (`1e-20` Pa); the solve was approximated, not performed. | Raise $P_{up}$, or treat the result as a limit rather than a solution. |
| `ValueError: The truth value of an array ... is ambiguous` from `arrhenius` | `arrhenius()` validates with a bare `if temperature <= 0`, so it takes **scalars only**. | Map over the array. Not documented in its docstring — see `docs/CHECKPOINT.md`, finding F2. |
| `warnings.warn` about trap occupancy θ > 0.9 | Oriani's local-equilibrium assumption needs low occupancy. | Lower the trap density or raise the temperature, or move to a McNabb–Foster kinetic model. |
| `warn`: total trap occupancy Σθ > 10 | Very strong trapping; the dilute-lattice simplification is past its range. | Same as above. The number is still returned — decide whether you trust it. |
| `warn`: reduction factor < 0.01 | Diffusion essentially stopped. | Check the trap density and binding energy are not in the wrong units. |
| Matplotlib refuses to save: image above 2¹⁶ px | `bbox_inches="tight"` sizes from the **artists**, so one annotation on a near-flat log axis can demand thousands of inches. It has reached 7415. | Widen the axis limits to a readable minimum, or drop the annotation. |
| Extrapolation beyond the measured temperature window, silently | The canonical wrappers do not re-raise the warning. | Check the data range yourself before trusting a number; {doc}`caveats`. |

### 1.4 Sensitivity analysis

| Symptom | Cause | Fix |
|---|---|---|
| A regime has almost no samples | An unaimed LHS over the default ranges lands almost entirely in one regime. | Use the presets. Figure 18 — this is what makes the sampling aimed. |
| A regime is skipped entirely | Its cluster is below `min_cluster` (300). | Raise `N`, or target that regime with `run_targeted_regime_scans`. |
| A parameter you believe matters ranks below the floor | It may be **pinned by the preset** that created the regime. A pinned parameter has no variance, so no index. | Check the preset first. A below-floor result for a pinned parameter says nothing about its leverage — Figure 20. |
| Every index looks significant | The dummy floor has not been drawn. | Carry the dummy parameter and read the ranking against it — Figure 19. |
| A preset that pins nothing | Previously silent; now fails loudly. | Fix the preset. |
| Reloaded scans disagree with the current study | `study.txt` records which study produced them. | Check the stamp; re-run rather than mixing material systems. |
| The run takes hours | The largest cluster draws >100k samples, and each metric costs PAWN and δ work. | Expected. Do not add a fourth metric casually — see Appendix B.2 on `eta_surf`. |

---

## 2. Adapting the model to a new material system

Every number the model runs on comes from one file. A new system means a new study
folder, not an edit to `calculations/`.

1. **Copy a study folder.** `cp -r Guo_etal_2025 MyAlloy_2026`. The folder is the unit;
   it holds `study_config.py` and the notebooks that read it.
2. **Edit `study_config.py` only.** The study files are structurally identical and
   define the same names in the same five sections, in dependency order. Change values,
   not names. {doc}`config-file` describes each block.
3. **Set `METALS` and `OXIDES`.** One entry each. The *first* key is taken as active, so
   the name is free but the order is not.
4. **Give every material its measured temperature window.** This is what makes the
   extrapolation warning possible. Omitting it does not disable the warning, it
   disables the protection.
5. **Add `surface_kinetics` to both the metal and the oxide** if you want the L6 family.
   A study without it supports Levels 1–5 only.
6. **Set `CONDITIONS`.** The operating point and the swept ranges both live here —
   state and geometry, not material identity.
7. **Set `OXIDE_DEFECTS`** including the `components` split, if defects matter.
8. **Set the regime presets** — `REGIME_PRESETS` and `REGIME_PRESETS_L5`. Without these
   the sensitivity analysis sees one regime. They are not an optimisation; they are what
   creates the study.
9. **Validate before trusting anything.** `python container/model_smoke_test.py`, then
   `check_against_config()` in the sensitivity workflow.
10. **Regenerate the reference values.** `--emit` on the container gate prints a
    ready-to-paste `REFERENCE[<study>]` block. Do not hand-edit pinned numbers.
11. **Re-run the notebooks from inside the new folder**, so the results land beside them.

```{important}
Do not hardcode material keys anywhere in `calculations/`. The keys are disjoint across
studies by design, and a hardcoded name is how a switch silently keeps reading the old
material. {doc}`how-to/switch-study` has the detail.
```

---

## 3. Reproducibility checklist

A result is reproducible when every line below can be answered.

- [ ] **Which study?** `ACTIVE_STUDY` printed, and `study.txt` present beside any scan.
- [ ] **Which commit?** The model is not frozen; a flux is only meaningful with a revision.
- [ ] **Which environment?** `mace_env`, or the container. Dependency bounds in
      `pyproject.toml` are ranges; the exact pins are in
      `container/requirements.lock.txt`.
- [ ] **Which operating point?** Temperature, both pressures, both thicknesses. A tier-4
      permeability without its operating point is not a number.
- [ ] **Which model?** Level 5 and Level 5L6 are different physics on the same wall.
- [ ] **Which seed?** `seed=42` by default; record it if changed.
- [ ] **Which preset?** A sensitivity ranking cannot be read without the preset that
      produced it.
- [ ] **Inside the measured temperature window?** If not, say so — the wrappers will not.
- [ ] **Gate passing?** `python container/model_smoke_test.py`.
- [ ] **Documentation checks passing?** `check_equations.py`, `make_figures.py`,
      `verify_docs.py`, `worked_values.py --check`, and `make html` under `-W`.

---

## 4. Code map

### 4.1 The model

| File | Role |
|---|---|
| [`calculations/config/model_config.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/config/model_config.py) | Resolves the active study and derives everything from it. The authoritative source for every value. |
| [`calculations/utils.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/utils.py) | Reference-temperature Arrhenius, and the property getters built on it. |
| [`calculations/permeation_calc.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/permeation_calc.py) | Level 1: Sieverts, Fick, bare-metal flux. Also the defective-metal flux entry points. |
| [`calculations/oxide_permeation.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/oxide_permeation.py) | Level 2a: Henry transport through the oxide, layer resistances, property lookup with the extrapolation warning. |
| [`calculations/interface_solver.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/interface_solver.py) | Level 2b: the flux balance and the `brentq` solve for $P_{int}$. The heart of the model. |
| [`calculations/parallel_oxide_defect_paths.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/parallel_oxide_defect_paths.py) | Level 3: area-weighted parallel paths, defect types, PRF. |
| [`calculations/defective_metal.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/defective_metal.py) | Level 4: Oriani trapping, grain-boundary enhancement, and their composition. The largest module. |
| [`calculations/surface_kinetics.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/surface_kinetics.py) | Level 6: Langmuir coverage, the surface rate law, and every L6 combination. |
| [`calculations/permeability.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/permeability.py) | Tiered permeability reporting, the unit refusal, and the pressure exponent. Changes no flux. |
| [`calculations/classify_regime.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/classify_regime.py) | Turns resistance fractions into regime labels, hierarchically by level. |
| [`calculations/sensitivity.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/sensitivity.py) | The whole sensitivity workflow, plus the model wrappers the sampler calls. |
| [`calculations/figures.py`](https://github.com/Azeezakinyemi999/Continuum_H_permeation/blob/main/calculations/figures.py) | Where a notebook's figures go, and the house text and frame style. |

### 4.2 The documentation's own tooling

| File | Role |
|---|---|
| `docs/figures/example_wall.py` | The fictitious wall every figure and check runs on. Also vendors the palette. |
| `docs/figures/make_figures.py` | Regenerates all 20 figures and `captions.md`. |
| `docs/figures/check_equations.py` | Restates every Part II equation independently and compares. 44 checks. |
| `docs/_tools/verify_docs.py` | Fails if a page names a function, module or parameter that does not exist. |
| `docs/_tools/worked_values.py` | Regenerates every real number quoted in the theory chapters; `--check` diffs. |
| `docs/_tools/verify_propagation.py` | Checks that a change propagates where it should. |
| `docs/CHECKPOINT.md` | Working state of the documentation rebuild. Not part of the built site. |

### 4.3 Studies and release

| Path | Role |
|---|---|
| `Guo_etal_2025/` | 316L / Cr₂O₃ — the active study. |
| `Fuerst_et_al_2024/` | Hastelloy N / Cr₂O₃. |
| `incoloy802_cr2o3/` | Incoloy 802 / Cr₂O₃. |
| `container/` | The reproducible image and the smoke-test gate. |
| `tests/` | Study resolution, operating conditions, notebook conventions. |

---

## 5. Planned improvements

```{note}
**Planned:** none of this is implemented. Each item is recorded with what it would
cost, so it can be picked up without re-deriving the finding. `CHECKLIST.md` at the
repository root is the live list; this section is the documentation-facing summary.
```

| Item | Why | Cost |
|---|---|---|
| Split `K_ox` into model-specific names | One symbol carries two incommensurable quantities. The configuration declares `K_ox_ref` in Sieverts units, which is right for Model 2 and wrong for Model 1. | A rename across the oxide path, plus re-pinning the container gate's reference values. |
| Document `arrhenius()` as scalar-only, or vectorise it | It raises on array input, and its docstring does not say so. | One docstring line, or a `np.asarray` guard. |
| Drop `Application/` from the study-resolution listing | It is offered as a valid study although it is excluded from releases. | One line, plus the tests that reference it. |
| Reconcile `METALS` keys across studies | Disjoint keys block a clean study switch. | Touches every study file; invalidates pinned values. |
| Add `eta_surf` to the regime metrics on demand | It carries the whole surface effect in one dimensionless number, but costs ~33% more analysis time. | A flag rather than a constant. |

---

## Where to go next

- The physics, from chemical potential up: {doc}`theory/foundations`
- The fork that most results depend on: {doc}`theory/two-models`
- Every warning in one place: {doc}`caveats`
- Symbols, settings and output columns: {doc}`appendices`
