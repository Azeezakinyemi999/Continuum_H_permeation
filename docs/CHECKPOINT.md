# Documentation checkpoint

Working state of the documentation rebuild: every decision already made, so none is
re-litigated; what is done; what is still open.

Read this first if you are resuming the work. Update it whenever a count or a scope
changes, not only at the start.

- **Baseline commit:** `1fcd605` (2026-10-09)
- **Status:** outline and figure inventory approved. **Scaffold complete** — the Example
  Wall, 20 figures and 44 equation checks all build and pass. Prose parts not started.
- **Environment:** `mace_env` (Python 3.9.23). Nothing is to be installed into it.

Everything must stay green. From the repository root:

```bash
python docs/figures/check_equations.py    # 44 / 44 checks passed
python docs/figures/make_figures.py       # 20 / 20 figures written
python docs/_tools/verify_docs.py         # clean
python docs/_tools/worked_values.py --check
cd docs && make html                      # builds with -W, so warnings are errors
```

---

## Decisions settled

Each is settled. Do not reopen one without recording the reason here.

### D1 — Adapt the existing chapters; do not rewrite from zero

The template says *"where an older document is still correct, adapt it and note its
provenance rather than rewriting it."* The theory chapters, the cross-cut chapters and
the verification tooling are the strongest part of this repository and in three places
already exceed what the template asks. A clean-sheet rewrite would destroy them.

Rebuild is therefore: **restructure to the template skeleton, adapt good prose into it,
write the missing apparatus new.**

### D2 — Keep Sphinx and MyST; translate the template's numbering into directories

The template is written for flat GitHub markdown (`docs/README.md`, `01-overview.md`,
`[!NOTE]`). This project has a working Sphinx + MyST site on GitHub Pages, built with
`-W`, with an API reference generated from docstrings.

Translation, applied consistently:

| Template | Here | Why |
|---|---|---|
| `docs/README.md` | `docs/index.md` | Sphinx landing page and toctree root |
| `01-…` `02-…` numbering | directory + toctree order | the sidebar carries the order; numeric prefixes would duplicate it and break existing links |
| `[!NOTE]` / `[!IMPORTANT]` | ```` ```{note} ```` / ```` ```{important} ```` | MyST directives; GitHub alerts do not render in Sphinx |
| relative links to code | `{doc}` roles and repo-relative paths | `verify_docs.py` already enforces these |

### D3 — There are two stage families, not one pipeline

The template assumes a single linear pipeline. This repository has two things:

1. **The model hierarchy, Levels 1–6** — cumulative physics, not a pipeline. Documented
   under `theory/`, one file per level.
2. **The regime-stratified sensitivity workflow** — a genuine pipeline with ordered
   stages. Documented under `sensitivity/`, one file per stage.

Both use the identical six-slot stage layout (D4). Part I's pipeline diagram shows the
flux path through the wall; the sensitivity part carries its own flow diagram.

### D4 — The uniform stage layout

Every file under `theory/level*.md` and `sensitivity/*.md` uses these headings, in this
order, with no others at `##` level:

1. **Purpose** — one paragraph
2. **At a glance** — table: inputs, outputs, code
3. **Concepts** — links into Part II, never restating an equation
4. **Procedure** — numbered, in the order the code runs
5. **Design choices** — why it is done this way; what was tried and rejected
6. **Checks and failure modes** — table: symptom, meaning

The existing chapters each end with *"What Level N cannot tell you"*, which is good and
survives as the closing paragraph of **Purpose** or of **Checks and failure modes**.

### D9 — Each equation is derived once, in its level chapter; `equations.md` is the catalogue

The template says Part II holds every equation and the stage sections link to it rather
than repeating them. Applied literally here it would gut the level chapters, because in
this project the levels **are** the theory — the "stages" are physics, not pipeline
steps.

The template's intent is that no equation is stated in two places. That is met by:

- each equation derived **exactly once**, in the level chapter that motivates it,
  inside its **Concepts** section;
- {doc}`theory/equations` kept as the consolidated catalogue and index, each entry
  tagged with the file and function that implements it;
- cross-cutting relations that belong to no single level — the sorption fork, the
  permeability tiers — living in their own chapters, which the levels link to.

### D5 — The fictitious example drives figures and checks; real numbers stay, but confined

The template says report nothing from real runs. This repository deliberately does the
opposite: `docs/_tools/worked_values.py` regenerates every quoted number by running the
current model, and `--check` diffs it, so a model change surfaces as a diff rather than
as quietly wrong prose. That is a real correctness mechanism and is **kept**.

Resolution, which honours both:

- A small fictitious wall — `docs/figures/example_wall.py`, round-number properties —
  drives **every figure** and **every equation check**. Figures therefore carry no study's
  results and stay valid when the active study changes.
- Real worked values stay, but move out of the derivations into clearly marked
  **Design choices** and **Worked values** blocks, which is where the template permits
  real-run numbers as evidence.

This also fixes a real coupling: the theory chapters currently read as though their
numbers were general, when they are specific to `ACTIVE_STUDY = Guo_etal_2025_316L`.

### D6 — Figures reuse the production visual language

`calculations/config/model_config.py` defines `COLORS` and `CURVE_STYLES`, chosen against
a colour-vision validator rather than eyeballed, with the rationale in comments at
`model_config.py:643-689`. Documentation figures **import those dictionaries** rather than
choosing new hues, so the docs and the real notebook output read as one system.

The production figures follow a fixed four-panel idiom, visible across every
`Proposal.ipynb`:

> **(A)** physics sweep · **(B)** temperature / Arrhenius · **(C)** geometry sweep ·
> **(D)** limit check: this level collapses to the one below

Panel **(D)** is the house style for showing a limit, and the documentation's limit-check
figures reuse it.

### D7 — Root documents keep the include-stub pattern

`docs/*.md` stubs `{include}` the root `OUTPUTS.md`, `PARAMETERS.md`, `STEADY_STATE.md`,
`TRAPPING_VALIDATION.md`, `PACKAGING_1.0.0.md`, `CONTAINER.md`. One source of truth: GitHub
renders the original, Sphinx renders the same bytes.

This is **better** than the template's "leave a header pointing at `docs/`" and is kept.
Superseded *content* — anything genuinely replaced by a new chapter — gets the header and
stays in place, per the template.

### D8 — `verify_docs.py` must stay clean at every commit

Current baseline: `clean — checked 37 pages against 151 functions, 14 modules, 70 config
names`. A part is not finished until `verify_docs.py`, `worked_values.py --check`,
`check_equations.py` and `make_figures.py` all pass and the `-W` Sphinx build succeeds.

---

## Open questions

Listed for the maintainer; the template says to ask rather than guess.

| # | Question | Blocks |
|---|---|---|
| Q1 | `K_ox` carries different units under the Henry and Sieverts models, and one symbol serves both laws. [parameter-map.md:160](theory/parameter-map.md) calls this an unresolved problem. Rename, or document and keep? | Appendix A notation, Part II |
| Q2 | `Application/` is excluded from releases but still on disk, still referenced by `tests/conftest.py:77` and three other tests, and holds the only tracked figure PNGs. Do the docs mention it at all? | Code map, Appendix C |
| Q3 | The tracked figures under `Application/Proposal-results/` are the only committed figures and come from a dropped folder. Retire them, or re-track under the live study? | Figure inventory |
| Q4 | `METALS` keys are still disjoint across study configs, which blocks a clean Fuerst switch. Document as a known limitation, or fix first? | `how-to/switch-study` |
| Q5 | The study-resolution error message offers `PERMEATION_STUDY=Application/study_config.py` as a valid choice, although `Application/` is excluded from releases. Drop it from the listing? | `how-to/switch-study` |

Answered, 2026-10-09:

- **Real-run numbers** — dual track (D5). Fictitious wall drives figures and checks;
  real worked values stay, confined to *Worked values* and *Design choices* blocks.
- **Q1, `K_ox` units** — document only, no code change. Figure 5 plus an IMPORTANT box
  and an Appendix A entry.
- **Q3, the `Application/` figure PNGs** — ignore them; documentation figures come only
  from `make_figures.py`. The old PNGs stay in place as history.
- **Sequencing** — scaffold first, then the prose parts. Scaffold is now done.

---

## Progress

| Part | State |
|---|---|
| CHECKPOINT | current |
| Outline + figure inventory | approved |
| `figures/example_wall.py` | **done** |
| `figures/check_equations.py` | **done** — 44 / 44 passing |
| `figures/make_figures.py` | **done** — 20 / 20 figures, captions generated |
| `index.md` — scope box, reading order, stage layout, units | **done** |
| Part I — `overview.md` | **done**, Figures 1–3 |
| Part II — theory, cross-cuts | **adapted**, Figures 4, 5, 15, 16 placed |
| Part III — Levels 1–6 | **done**, all six in the six-slot layout, Figures 6, 8–14 |
| Part IV — `sensitivity/` | **done**, index + 5 stages, Figures 17–20 |
| Part V — `practical.md` | **done** |
| Appendices A–G | **done** |

All 20 figures are placed and rendered into the build.

## Findings, from reading the code to write the checks

The template says to fix what documenting uncovers rather than document the wrong
behaviour. Each of these came out of writing `check_equations.py` and `make_figures.py`.

| # | Finding | State |
|---|---|---|
| F1 | `trap_occupancy` documented its return key as `'theta_T'`; it returns `'theta'`, and its own caller at `defective_metal.py:968` reads `'theta'`. The published API reference named a key that does not exist. | **fixed**, `defective_metal.py:205` |
| F2 | `arrhenius()` validates with a bare `if temperature <= 0`, so it raises on array input. Scalar-only, and the docstring does not say so. | documented at the call site in `make_figures.py`; candidate CHECKLIST item |
| F3 | Grain-boundary enhancement is negligible at structural grain sizes — 1.0016× at 20 µm, 1.66× even at 50 nm. Level 4 is a trapping model with a GB term that rounds away. This was not stated anywhere. | now Figure 12; belongs in the Level 4 chapter |
| F4 | Permeability is reported in **four** tiers, not three (intrinsic, effective-layer, stack, apparent), per commit `5d3e94a`. Earlier notes said three. | corrected; Figure 15 |
| F5 | The study-resolution error lists `Application/` as a valid study although it is excluded from releases. | open, Q5 |
| F6 | **`worked_values.py --check` was failing on `main`.** Commit `366b70d` ("Get the GB enhancement from the metal, not from a table of austenitic steel") changed `D_eff`, and the stored copy was never regenerated — so the published docs quoted stale numbers in `getting-started.md`, `theory/level4-microstructure.md`, `theory/level5-full-system.md` and `caveats.md`. | **fixed** — regenerated and every prose quote corrected |
| F7 | `theory/level5-full-system.md`'s values block was computed at a **48 nm** coating while the configured thickness is **1 µm**, and the block is not covered by `worked_values.py`, so it could not be caught. | **fixed** — replaced with the canonical wrapper values, which are generated and diffed |
| F8 | The same 48 nm staleness ran through `level2-oxide.md` (α, $P_{int}$, flux, resistances) and `level6-surface-kinetics.md`, where the L5L6 block still reported **`regime = metal`**. At the configured 1 µm the wall is **`regime = oxide`** — `frac_oxide` is 0.696, not 0.107. The chapter asserted the opposite of what the model now returns. | **fixed**, all recomputed and verified |

### The blocks that no generator covers

F6, F7 and F8 are one root cause: **prose blocks of real numbers that `worked_values.py`
does not emit cannot be caught when the model moves.** `worked_values.py` covers the
canonical wrappers; it does not cover $\alpha$ and $\beta$, the layer resistances, the
Level 4 mode table, or the Level 5 recovery limits — and every one of those had drifted.

Each is now corrected against a computed value, and the chapters say so where an older
figure is likely to have been inherited. The durable fix is to extend
`worked_values.py` to emit them, which is **not done** and is the single most valuable
follow-up in this file.

### What F6 confirmed

The same commit that made the worked values stale is the reason Figure 12 says what it
says. Reading the GB enhancement from the metal rather than from a generic
austenitic-steel table dropped $D_{gb}/D_L$ from a tabulated **50** to the metal's
**1.77**, which collapsed the grain-boundary contribution from +0.07% to **+0.001%**.
The `'both'` and `'trapping_only'` rows of the Level 4 mode table now agree to five
significant figures. Level 4 is a trapping model.

## Corrections already made to earlier notes

- `GATE_STUDY` no longer points at the dropped `Application/` folder. Commit `a49cc31`
  ("Stop the gate defaulting to a study that no longer ships") set it to
  `Guo_etal_2025/study_config.py`. Verified at `container/model_smoke_test.py:59`.
- Permeability has four tiers, not three. See F4.
