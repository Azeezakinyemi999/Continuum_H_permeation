# MHI Permeation Model

Reference documentation for a steady-state analytical model of hydrogen permeation
through an oxide-coated structural alloy wall: what each level does, why it is done that
way, and the equations behind it.

:::{note}
**Scope.** This documents the **method**, not any result.

- It applies to **any** metal/oxide wall the configuration can describe, not only the
  study that happens to be active.
- Every **figure** is drawn from a small fictitious system — the *Example Wall* — or is
  a labelled schematic. No figure reports a real run.
- Real numbers do appear, confined to blocks marked **Worked values** or **Design
  choices**, where a measured number is the only honest justification for a threshold.
  Those are regenerated from the live model and diffed, so they cannot go quietly stale.

**Checked against** Python 3.9.23, NumPy 1.26.4, SciPy 1.13.1, matplotlib 3.9.4,
pandas 2.3.3, SALib 1.5.0, Sphinx 7.4.7 — 2026-10-09, at commit `1fcd605`.

Re-check it all:

```bash
python docs/figures/check_equations.py     # 44 / 44 checks passed
python docs/figures/make_figures.py        # 20 / 20 figures written
python docs/_tools/verify_docs.py          # clean
python docs/_tools/worked_values.py --check
cd docs && make html                       # -W: warnings are errors
```
:::

---

## How to use this document

### Reading order

Nothing here needs to be read front to back. Pick the entry that matches what you came
for.

| If you want to | Start at | Then |
|---|---|---|
| Understand what the model is for | {doc}`overview` | {doc}`theory/two-models` |
| Run something today | {doc}`getting-started` | {doc}`how-to/switch-study` |
| Understand the physics properly | {doc}`theory/foundations` | the level chapters, in order |
| Look up one equation | {doc}`theory/equations` | the level chapter that derives it |
| Know which parameters matter | {doc}`sensitivity/index` | stages 1–5 in order |
| Fix something that broke | {doc}`practical` | {doc}`caveats` |
| Look up a symbol, setting or column | {doc}`appendices` | {doc}`PARAMETERS`, {doc}`OUTPUTS` |

```{important}
Read {doc}`theory/two-models` before comparing any two levels. Levels 1–5 and the L6
family describe the same wall under different physics, and a conclusion does not carry
across that fork. It is the most common way to misread this model.
```

### The layout every stage section follows

Each level chapter and each sensitivity stage uses the same six headings, in this
order. Once you have read one, you know where to look in all of them.

| Section | Holds |
|---|---|
| **Purpose** | What this stage is for, in a paragraph |
| **At a glance** | Inputs, outputs, code, cost |
| **Concepts** | The ideas needed, linking to Part II rather than restating equations |
| **Procedure** | Numbered, in the order the code runs it |
| **Design choices** | Why it is done this way; what was tried and rejected |
| **Checks and failure modes** | Symptom and meaning, as a table |

### The running example

Every figure and every equation check runs on the **Example Wall**: a fictitious system
with round-number properties, defined in `docs/figures/example_wall.py`. It is not a
real material. It was chosen so that each intermediate quantity can be checked by hand,
the oxide and metal carry comparable resistance, and the Level 2b interface pressure has
an exact closed form to test the solver against.

Symbols in the text are **generic** — $D_m$, $K_s$, $\alpha$. The project's own
identifiers (`D_metal`, `f_pinhole`) appear only in {doc}`appendices` B and C, the code
map, and examples.

### Units

SI throughout, matching the configuration. Full table in {doc}`appendices` A; the
model-dependent ones in {doc}`theory/two-models`.

| Quantity | Units |
|---|---|
| flux $J$ | mol m⁻² s⁻¹ |
| pressure $P$ | Pa |
| diffusivity $D$ | m² s⁻¹ |
| metal solubility $K_s$ | mol m⁻³ Pa⁻⁰·⁵ |
| oxide solubility $K_{ox}$ | mol m⁻³ Pa⁻¹ (Model 1) · mol m⁻³ Pa⁻⁰·⁵ (Model 2) |
| energy | J mol⁻¹ |
| temperature | K |

### The two boxes

```{important}
Behaviour that is easy to get wrong. Every one of these has cost someone real time.
```

```{note}
**Planned:** described but not implemented. Nothing in a box like this exists yet;
see {doc}`practical` §5.
```

---

## Contents

```{toctree}
:maxdepth: 2
:caption: Overview

overview
getting-started
```

```{toctree}
:maxdepth: 2
:caption: How-to

how-to/switch-study
how-to/sensitivity-analysis
CONTAINER
```

```{toctree}
:maxdepth: 2
:caption: Part II — Theory

theory/two-models
theory/foundations
theory/equilibrium-models
STEADY_STATE
theory/equations
```

```{toctree}
:maxdepth: 2
:caption: Part III — The levels

theory/level1-metal
theory/level2-oxide
theory/level3-defective-oxide
theory/level4-microstructure
theory/level5-full-system
theory/level6-surface-kinetics
```

```{toctree}
:maxdepth: 2
:caption: Cross-cutting

theory/permeability
theory/parameter-map
caveats
```

```{toctree}
:maxdepth: 2
:caption: Part IV — Sensitivity

sensitivity/index
```

```{toctree}
:maxdepth: 2
:caption: Part V — Practical

practical
appendices
```

```{toctree}
:maxdepth: 2
:caption: Reference

config-file
PARAMETERS
OUTPUTS
api/index
```

```{toctree}
:maxdepth: 2
:caption: Validation and release

TRAPPING_VALIDATION
PACKAGING_1.0.0
references
```

---

## Provenance

```{note}
The theory chapters and how-to guides were authored from the project's earlier design
notes, with every factual claim re-checked against the current code. That absorption is
complete — the original notes have been retired.

The pages with uppercase filenames ({doc}`PARAMETERS`, {doc}`OUTPUTS`,
{doc}`STEADY_STATE`, {doc}`TRAPPING_VALIDATION`, {doc}`PACKAGING_1.0.0`,
{doc}`CONTAINER`) are generated or release-record documents. They stay at the repository
root so they remain readable on GitHub, and are included here by reference rather than
copied, so there is one source of truth for each.

`docs/CHECKPOINT.md` holds the working state of the documentation itself — decisions
settled, progress, what is still open. It is deliberately not part of this site.
```
