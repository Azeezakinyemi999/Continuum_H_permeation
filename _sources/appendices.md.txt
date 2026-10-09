# Appendices

Reference tables. Nothing here introduces physics; every entry points back at the
chapter that derives it.

```{note}
Symbols are **generic** throughout the theory chapters — $D_m$, $K_s$, $\alpha$. The
project's own identifiers (`D_metal`, `K_s_metal`, `f_pinhole`) appear only in
Appendix B, Appendix C and the code map. Appendix A is the bridge between the two.
```

---

## A. Notation

Units follow the configuration, which is SI throughout. Where a symbol's units depend
on which model is in force, both are given and the difference is the subject of
{doc}`theory/two-models`.

### A.1 Transport

| Symbol | Meaning | Units | Derived in |
|---|---|---|---|
| $J$ | hydrogen flux through the wall | mol m⁻² s⁻¹ | {doc}`theory/level1-metal` |
| $D_m$ | metal lattice diffusivity | m² s⁻¹ | {doc}`theory/foundations` |
| $D_{ox}$ | oxide diffusivity | m² s⁻¹ | {doc}`theory/level2-oxide` |
| $D_{eff}$ | metal diffusivity after microstructure | m² s⁻¹ | {doc}`theory/level4-microstructure` |
| $K_s$ | metal solubility, Sieverts | mol m⁻³ Pa⁻⁰·⁵ | {doc}`theory/foundations` |
| $K_{ox}$ | oxide solubility | Model 1: mol m⁻³ Pa⁻¹ · Model 2: mol m⁻³ Pa⁻⁰·⁵ | {doc}`theory/two-models` |
| $L_m$ | metal thickness | m | — |
| $L_{ox}$ | oxide thickness | m | — |
| $\Phi_m = D_mK_s$ | intrinsic metal permeability | mol m⁻¹ s⁻¹ Pa⁻⁰·⁵ | {doc}`theory/permeability` |
| $\Phi_{ox} = D_{ox}K_{ox}$ | intrinsic oxide permeability | units follow $K_{ox}$ | {doc}`theory/permeability` |
| $\alpha = \Phi_{ox}/L_{ox}$ | oxide permeance | mol m⁻² s⁻¹ Pa⁻¹ (Model 1) | {doc}`theory/level2-oxide` |
| $\beta = \Phi_m/L_m$ | metal permeance | mol m⁻² s⁻¹ Pa⁻⁰·⁵ | {doc}`theory/level2-oxide` |
| $\Pi = J/\Delta\sqrt{P}$ | apparent permeance | mol m⁻² s⁻¹ Pa⁻⁰·⁵ | {doc}`theory/permeability` |
| $R_{ox}, R_m$ | layer resistances | Pa s m² mol⁻¹ | {doc}`theory/level2-oxide` |
| $n = \mathrm{d}\ln J/\mathrm{d}\ln P$ | pressure exponent | — | {doc}`theory/permeability` |

### A.2 Pressures and the operating point

| Symbol | Meaning | Units |
|---|---|---|
| $P_{up}$ | upstream hydrogen pressure | Pa |
| $P_{int}$ | oxide/metal interface pressure, **solved** | Pa |
| $P_{down}$ | downstream (sweep) pressure | Pa |
| $P_{virtual} = g(\theta)^2$ | pressure the coverage is equivalent to | Pa |
| $T$ | absolute temperature | K |

```{important}
$P_{int}$ is not an input. It is the root of a flux balance, and almost every
surprising result in this model is a statement about what $P_{int}$ did — see
Figure 6 and Figure 14.
```

### A.3 Surface kinetics — Model 2 only

| Symbol | Meaning | Units |
|---|---|---|
| $\theta$ | fractional surface coverage | — (0 to 1) |
| $k_{diss}$ | dissociation rate constant | mol m⁻² s⁻¹ Pa⁻¹ |
| $k_{recomb} = k_{diss}/K_{eq}$ | recombination rate constant | mol m⁻² s⁻¹ |
| $K_{eq}$ | surface equilibrium constant | Pa⁻¹ |
| $g(\theta) = \dfrac{\theta}{(1-\theta)\sqrt{K_{eq}}}$ | coverage as an effective $\sqrt{P}$ | Pa⁰·⁵ |
| $\eta_{surf}$ | surface efficiency, $\Phi_{app}=\Phi_{transport}\eta_{surf}$ | — |

### A.4 Microstructure and defects

| Symbol | Meaning | Units |
|---|---|---|
| $N_T$ | trap density | m⁻³ |
| $N_L$ | interstitial lattice site density | m⁻³ |
| $E_b$ | trap binding energy | J mol⁻¹ |
| $K_i = \exp(E_b/RT)$ | trap equilibrium constant | — |
| $d$ | grain size | m |
| $\delta_{gb}$ | grain-boundary thickness | m |
| $f_{gb}$ | grain-boundary volume fraction | — |
| $s = D_{gb}/D_m$ | grain-boundary enhancement factor | — |
| $f_d$ | defect area fraction | — |
| PRF | permeation reduction factor | — |

### A.5 Thermodynamics

| Symbol | Meaning | Units |
|---|---|---|
| $R$ | gas constant, `8.314` in `calculations.utils` | J mol⁻¹ K⁻¹ |
| $E_D$ | activation energy for diffusion | J mol⁻¹ |
| $H_s$ | heat of solution | J mol⁻¹ |
| $Q_{perm} = E_D + H_s$ | permeability activation energy | J mol⁻¹ |
| $T_{ref}$ | temperature a reference value was measured at | K |

```{note}
The code uses the **reference-temperature** Arrhenius form,
$k(T)=k_{ref}\exp\left[-\frac{E}{R}\left(\frac{1}{T}-\frac{1}{T_{ref}}\right)\right]$,
not $k_0\exp(-E/RT)$. It is anchored on a value measured at a stated temperature
rather than on a pre-exponential extrapolated to infinite temperature. Every
`*_ref` setting in Appendix B carries an implied $T_{ref}$.
```

---

## B. Settings reference

Values are those in force at the documented commit. Study-specific material values are
**not** listed here — they change with the active study and live in {doc}`PARAMETERS`.
These are the settings that are part of the *method*.

### B.1 Study resolution

| Setting | Value | Where it acts | Meaning |
|---|---|---|---|
| `STUDY_FILENAME` | `study_config.py` | `config.model_config` | the filename searched for |
| `STUDY_ENV_VAR` | `MHI_STUDY` | `config.model_config` | names a study explicitly; takes precedence |
| `ACTIVE_STUDY` | resolved, not set | `config.model_config` | the study actually in force |

There is no default study. Resolution order is `$MHI_STUDY`, then a `study_config.py`
in the working directory or any parent up to the repository root, then an error that
lists what it can see. See {doc}`config-file`.

### B.2 Sensitivity workflow

| Setting | Value | Where it acts | Meaning |
|---|---|---|---|
| `REGIME_LABELS` | `('surface', 'oxide', 'metal')` | L5L6 | regimes under surface kinetics |
| `REGIME_LABELS_L5` | `('oxide', 'metal', 'defect')` | L5 | regimes without surface kinetics |
| `REGIME_SA_METRICS` | `['flux', 'permeability', 'theta']` | stage 6 | metrics an index is computed for |
| `LOG_METRICS_DEFAULT` | `('flux', 'permeability')` | stage 6 | analysed on a log₁₀ scale |
| `TARGET_CLUSTER_SIZE` | `1500` | stage 4 | rows wanted in each regime cluster |
| `min_cluster` | `300` | stage 3 | below this a regime is not analysed |
| `seed` | `42` | stages 2, 5 | sampling seed |
| `STUDY_STAMP_FILE` | `study.txt` | stages 2, 5 | records which study produced a scan |

```{important}
`eta_surf` is deliberately **absent** from `REGIME_SA_METRICS`, although it is the
sharpest surface metric available: $\Phi_{app}=\Phi_{transport}\,\eta_{surf}$ exactly.
Analysing a fourth metric costs roughly 33% more PAWN and δ work, and the largest
cluster already runs for hours. It is still written to every scan CSV, so the data is
there — add it to the list when a run is specifically about surface kinetics.
```

### B.3 Figures

| Setting | Value | Where it acts | Meaning |
|---|---|---|---|
| `DPI` | `450` | `calculations.figures` | raster resolution for notebook output |
| `VECTOR_FORMAT` | `"pdf"` | `calculations.figures` | vector copy written beside each PNG |
| `RESULTS_SUFFIX` | `"-results"` | `calculations.figures` | appended to a notebook's name |
| `INLINE_FORMAT` | `"retina"` | `calculations.figures` | inline display density in Jupyter |

```{important}
Matplotlib refuses an image above 2¹⁶ px in either direction, which at 450 dpi is
145.6 inches. Because `bbox_inches="tight"` sizes the output from the **artists** and
not from `figsize`, one annotation on a near-flat log axis can ask for thousands of
inches. It has happened: a Level 4 grain-size panel demanded 7415.
```

### B.4 Numerics

| Setting | Value | Where it acts | Meaning |
|---|---|---|---|
| solver | `brentq` | `interface_solver`, `surface_kinetics` | bracketed root finder |
| bracket | $[P_{down}, P_{up}]$ | `solve_interface_pressure` | exactly one root, guaranteed |
| coverage bracket | $[10^{-10},\,1-10^{-10}]$ | `solve_steady_state_flux_*` | $\theta$ is bounded by construction |
| `min_pressure` | `1e-20` Pa | `interface_solver` | below this the solve is approximated, `converged=False` |

---

## C. Output files

### C.1 Scan CSVs — `scan_<regime>.csv`

Written by `run_targeted_regime_scans(save_dir=...)`, one per regime, reloaded by
`load_regime_scans`. One row per sample.

| Column group | Columns | Meaning |
|---|---|---|
| inputs | one per swept parameter | the sampled values |
| flux | `flux` | wall flux, mol m⁻² s⁻¹ |
| permeability | `permeability`, `permeance`, `Phi_transport` | tiered, see {doc}`theory/permeability` |
| surface | `theta`, `eta_surf` | coverage and surface efficiency |
| split | `frac_surface`, `frac_oxide`, `frac_metal` | resistance fractions; these define the regime |
| wall | `PRF`, `D_eff`, `P_interface` | reduction factor, effective diffusivity, solved interface pressure |
| paths | `flux_intact`, `flux_defect` | the two parallel contributions |
| labels | `regime`, `system_rate_limiting`, `dominant_path` | string diagnostics |

```{note}
The `frac_*` columns are **excluded** from the sensitivity metrics: within a regime
cluster they are degenerate, because the cluster was defined by them.
```

### C.2 Study stamp — `study.txt`

One line recording which study produced the scans in that folder. Its purpose is to
make a mismatched reload fail loudly instead of silently analysing another material
system's samples.

### C.3 Comparison matrices

`compare_delta_flux.csv` and siblings: one row per parameter, one column per regime,
holding the index for that metric. This is what the regime comparison heatmap and the
parallel-coordinates plot read.

### C.4 Figures

Each notebook writes into `<Notebook>-results/` beside itself — a PNG at 450 dpi and a
PDF of the same stem. The folder is resolved against the **working directory**, which
for a correctly launched notebook is its own study folder.

---

## D. Unit conversions

| From | To | Multiply by |
|---|---|---|
| bar | Pa | $10^{5}$ |
| MPa | Pa | $10^{6}$ |
| torr | Pa | $133.322$ |
| atm | Pa | $101\,325$ |
| µm | m | $10^{-6}$ |
| nm | m | $10^{-9}$ |
| kJ mol⁻¹ | J mol⁻¹ | $10^{3}$ |
| eV | J mol⁻¹ | $96\,485$ |
| °C | K | add $273.15$ |
| mol m⁻² s⁻¹ | µmol m⁻² s⁻¹ | $10^{6}$ |

Pressure exponents convert the same way only when the exponent is the same. A
permeability in mol m⁻¹ s⁻¹ Pa⁻⁰·⁵ and one in mol m⁻¹ s⁻¹ Pa⁻¹ are **not** related by
a unit factor; see Figure 5.

```{important}
Whether a quoted flux counts moles of H or moles of H₂ changes it by a factor of two,
and that convention is not recorded in the configuration. The model is internally
consistent — $K_s$, $C$ and $J$ all use the same species — so every ratio, fraction
and PRF in this documentation is unaffected. Fix the convention against your own
source data before comparing an **absolute** flux with an experiment.
```

---

## E. Checking the equations

`docs/figures/check_equations.py` restates every equation in Part II independently, in
plain arithmetic on the Example Wall, and compares it with what `calculations` returns.
A check that merely called the same function twice would pass forever and prove
nothing, so no expected value is borrowed from the module under test.

Run it from the repository root:

```bash
python docs/figures/check_equations.py
```

### E.1 What it checks

| Group | Checks |
|---|---|
| sorption and Fick | Sieverts $C=K_s\sqrt{P}$, $C(0)=0$, Henry $C=K_{ox}P$, Fick |
| Level 1 | flux, $\Phi=D_mK_s$, $C_{up}$, and the $\sqrt{P}$ signature |
| Level 2a | flux, $R_{ox}$, linearity in $P$, and $R_{metal}$ |
| Level 2b | $P_{int}$ against the **closed-form quadratic root**, flux, flux balance, and that the series flux is below both single-layer fluxes |
| Level 2b limit | $L_{ox}\to0$ recovers Level 1 |
| Level 3 | the area-weighted **sum**, and that contributions sum to the total |
| Level 3 limit | $f_d\to0$ recovers Level 2b |
| Level 4 | $K=\exp(E_b/RT)$, the trapping term, $D_{eff}$, and that trapping retards |
| Level 4 limit | an empty trap list recovers $D_{lattice}$ exactly |
| Arrhenius | the reference-temperature form, and $k(T_{ref})=k_{ref}$ |
| Level 6 | $g(\theta)$, the closed-form $\sqrt{P_{int}}$, the rate law, flux balance at the solution, $0<\theta<1$, and that a finite surface rate costs flux |
| permeability tiers | tier-1 values, the units each model assigns, that Model 1 **refuses** to combine, and the Model 2 series value |
| pressure exponent | $n=\tfrac12$ for Sieverts, $n=1$ for Henry, and strictly between for the mixed wall |
| documentation | the vendored palette matches `model_config.COLORS` |

Where no independent form exists — the Level 2b and Level 6 solves — the check is an
**invariant** instead: flux balance must hold at the solution, and a limit must
collapse to the level below.

### E.2 Expected output

The last line is the one to read:

```text
44 / 44 checks passed
```

Exit code is 0 when every check passes and 1 otherwise. A `FAIL` line names the check,
what was computed, what the documented equation gives, and the relative difference.

```{note}
Two tolerances are used. Equations that should agree to machine precision are held to
`1e-9`; anything that went through a root finder is held to `1e-6`. The Level 2b
closed-form check currently agrees to about `2e-16`.
```

---

## F. Glossary

**Apparent permeability** — a permeability backed out of a solved flux. Tier 4: a
property of the wall *and* the experiment, so it is always reported with its operating
point.

**Dummy parameter** — a parameter the model never reads, carried through the whole
sensitivity analysis so that whatever index it earns defines the noise floor.

**Example Wall** — the small fictitious system in `docs/figures/example_wall.py` that
every figure and equation check runs on. Not a real material.

**Henry's law** — concentration linear in pressure, $C=K_{ox}P$. Applies when
molecular H₂ dissolves intact, which in this model is the oxide under Model 1.

**Intact path** — the part of the wall area where the oxide is unbroken, as opposed to
the defect path. The two are in parallel.

**Model 1** — dissociation at the oxide/metal interface. Henry oxide, Sieverts metal.
Levels 1–5.

**Model 2** — dissociation at the gas/oxide surface. Sieverts oxide, Sieverts metal.
The L6 family.

**Oriani model** — local-equilibrium trapping, in which trapped and lattice hydrogen
are in equilibrium at every point. Valid while trap occupancy stays low; the
alternative when it does not is McNabb–Foster.

**PAWN** — a density-based global sensitivity index, computed from the distance between
conditional and unconditional output distributions. Needs no variance decomposition.

**Permeance** — a permeability divided by a thickness: what a *specific* wall passes,
rather than what the material passes.

**PRF** — permeation reduction factor. How much less flux the coated wall passes than
the bare one. PRF = 1 means the coating does nothing.

**Pressure exponent** $n$ — $\mathrm{d}\ln J/\mathrm{d}\ln P$. A model-free diagnostic
of which step controls the wall: $\tfrac12$ for a Sieverts step, $1$ for a Henry layer.

**Regime** — which physical process carries most of the resistance. Assigned from the
resistance fractions, and the basis on which the sensitivity analysis is stratified.

**Sieverts' law** — concentration in the square root of pressure, $C=K_s\sqrt{P}$.
Applies when H₂ dissociates before dissolving, which is always true of the metal.

**Study** — a complete physical system: a metal, an oxide, a microstructure, a defect
set and an operating point, living in one `study_config.py` inside its own folder.

**δ (delta)** — Borgonovo's moment-independent importance measure. Like PAWN, it needs
no variance decomposition, which is why both survive here.

---

## G. References

The source list lives in {doc}`references`, which is still current and is kept as the
single bibliography rather than restated here. It covers foundational theory, the
material-property sources behind each study, and the sensitivity-analysis methodology
— including the documented failure modes of PAWN.

Material values themselves are not stored in this documentation. They live in each
study's `study_config.py`; see {doc}`config-file`.
