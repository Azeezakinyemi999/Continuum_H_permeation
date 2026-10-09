# Stage 1 — Validate the configuration

## Purpose

Catch a configuration that disagrees with itself before spending hours sampling it. The
material values cannot drift — they are derived from `METALS`, `OXIDES`,
`MICROSTRUCTURE`, `OXIDE_DEFECTS` and `CONDITIONS`, so drift is impossible by
construction. What this stage checks is everything derivation does *not* guarantee.

## At a glance

| | |
|---|---|
| **Inputs** | `DEFAULT_PARAMS_LEVEL5`, `DEFAULT_PARAMS_LEVEL5L6`, `SUGGESTED_RANGES_*`, the presets |
| **Outputs** | A list of problem strings, empty when clean |
| **Code** | [`check_against_config`](../api/sensitivity), `preset_range_divergence` |
| **Cost** | Seconds |

## Concepts

A **default** is the single value a parameter takes when it is not being swept. A
**range** is the interval it is swept over. A **preset** is a set of range overrides
that aims the sampler at one regime. The three must agree, and nothing in the
derivation makes them.

Units and symbols: {doc}`../appendices` A. Settings: {doc}`../appendices` B.2.

## Procedure

In the order the code runs them:

1. **Every derived value is a finite scalar.** Catches a config entry that has become a
   list or a sweep range. The materials' `sample_thickness` entries are lists, which is
   one reason geometry reads from `CONDITIONS` rather than from the material.
2. **Every swept parameter has a default.** A parameter with a range but no default
   cannot be held fixed when another is swept.
3. **Every default lies strictly inside its own range.** Strictly, not loosely — a
   default sitting exactly *on* a bound is only half-sampled. That is how
   `K_eq_metal_ref` once slipped through.
4. **L5 and L5L6 agree on everything they share** — the inherited defaults and ranges,
   and the preset overrides both levels have in common.
5. Problems are returned as strings; `raise_on_fail=True` turns them into an exception.

## Design choices

**Why check at all, when the values are derived.** Derivation guarantees that a value
matches its material source. It does not guarantee that the value is a scalar, that it
has a range, that the range brackets it, or that two levels agree. Each of those has
failed in this repository.

**Why "strictly inside" rather than "within".** A default on a bound is sampled from
one side only, so half its neighbourhood is never explored and the resulting index
understates it. The strict test was added after `K_eq_metal_ref` passed a loose one.

**Why a preset that pins nothing now fails loudly.** It used to pass silently, which
meant a preset could target a regime while constraining nothing, producing a scan
indistinguishable from an unaimed one. Silence made that look like a result.

## Checks and failure modes

| Symptom | Meaning |
|---|---|
| `X is a list — expected a scalar` | A config entry is a sweep range or a sample list. Check its `model_config` source. |
| `X is swept but has no default` | A range exists with nothing to hold the parameter at. |
| `X default sits on a bound` | Half-sampled. Move the default inside, or widen the range. |
| `L5 and L5L6 disagree on X` | The two levels inherited different values for a shared quantity. |
| A preset that pins nothing | The preset constrains no parameter, so it does not aim anything. |
| Clean, but the scan still yields one regime | Not a configuration fault. The ranges are self-consistent and simply do not reach the other regimes — stage 4. |
