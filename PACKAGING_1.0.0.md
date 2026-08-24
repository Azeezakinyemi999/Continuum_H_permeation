# Packaging decisions — `hydrogen-model:1.0.0` offline container

Scope record for the offline Docker deliverable. Captures decisions that are not
recoverable from the code, so a later reader does not have to re-derive them.

Status: **Phase A complete** (repo hygiene). Not yet built.

## Target environment

| Decision | Value | Why |
|---|---|---|
| Python | **3.12** | 3.9 reached end of life in Oct 2025 and the company is expected to vulnerability-scan the image. The numpy<2 ceiling that forced 3.9 in `mace_env` comes from **torch**, which this project never imports, so the constraint does not apply here. Requires a numerical revalidation gate (below). |
| Platform | `linux/amd64` | Set at build time (`docker build --platform linux/amd64`), not in `FROM`, so the image definition stays reusable. Dev host is already x86_64/amd64. |
| Base image | `python:3.12-slim` | |
| Release tag | `hydrogen-model:1.0.0` | Promoted from `1.0.0-rc*` only after every acceptance gate passes. Never `latest`. |
| Container user | non-root, arbitrary-UID tolerant | See "Bind-mount ownership". |

### Numerical revalidation gate (blocks 1.0.0)

Moving 3.9 → 3.12 changes the interpreter under the solvers, so reference values
must be re-established rather than assumed. Baseline measured on `mace_env`
(Python 3.9.23, numpy 1.26.4, scipy 1.13.1), study `incoloy802_cr2o3`, at
`T = 700.0 K`, `P_up = 1e5 Pa`, `P_down = 1e2 Pa`, `L_metal = 1e-3 m`,
`L_oxide = 4.8e-8 m`:

| quantity | 3.9 baseline |
|---|---|
| `arrhenius(1e-11, 50000, 800, 700)` | `2.926830409475082e-11` |
| L1 `calculate_simple_metal_flux(...)['flux']` | `1.955109516956673e-07` |
| L2b `flux` (brentq) | `1.888003770552079e-07` |
| L2b `P_interface` | `93462.901558976` |
| L2b `resistance_ratio` | `0.036154225782795736` |
| `classify_regime_level14(0.2)['regime_hierarchy']` | `metal_limited/traps_defect_limited` |
| `classify_regime_level14(0.9)['regime_hierarchy']` | `metal_limited/lattice_limited` |

Compare with `math.isclose(rtol=1e-9)`, not equality — BLAS differs between the
macOS dev host and linux/amd64, which can move the last few bits without any
scientific meaning. A disagreement larger than `rtol=1e-9` is a real finding and
blocks the release.

Two API details for whoever writes `model_smoke_test.py`:
`get_metal_properties_at_T` returns `D_metal` / `K_s_metal` (not `D` / `K_s`),
and `calculate_oxide_metal_system` needs a `thickness` key injected into **both**
props dicts — the T-evaluated getters do not supply it for the metal.

## What ships

```
workspace/
├── calculations/     # whole package: 14 files incl. config/ and config/studies/
└── Application/
    ├── the six delivery notebooks
    └── sa_results/   # RUNTIME INPUT, not output — see below
```

`Application/sa_results/*.csv` are **inputs**, not disposable artifacts. Both
`regime_parallel_coords*.ipynb` only read them (`master_clusters.csv`,
`routeB_givendata.csv`, `compare_delta_flux.csv`) and cannot regenerate them —
regeneration requires the expensive scans in `sensitivity_regime_L5*.ipynb`. They
must travel with the workspace. This is why scanning `import` statements is not
sufficient to establish an offline package: file dependencies matter too.

Preserve `calculations/config/` **without** an `__init__.py`. It is a working
implicit namespace package; adding one during packaging would be a gratuitous
change to a validated tree.

## What does not ship: `data/`

`data/` is excluded from the deliverable. It is legacy: commits `de9f695`
("Read material data from model_config, not the data/ shims") and `a6f1bb9`
("Centralise sensitivity parameters in model_config.py") moved the codebase off
it. Verified by grep over the live tree — the only importers are notebooks that
are **not** in the delivery set (`plot&docs/analysis.ipynb`, `docs/old_ipynb/*`,
`Application/metal_Al2O3/Surface_chemistry copy.ipynb`). No module in
`calculations/` imports it (only two commented-out lines in `interface_solver.py`),
and nothing reads `data/sensitivity_parameters.csv` or
`data/literature_permeation_data.csv`.

Excluding it also frees the name `inputs/` if company-supplied files are ever
needed, since `data/` here is a Python package rather than a data directory.

## Output paths: unchanged on purpose

`RESULTS_DIR` is a bare relative name (`"sa_results"`, `"sa_results_L5"`) and
`Proposal.ipynb` calls `savefig('Level1_Perfect_Metal.png')` with no directory, so
outputs land in each notebook's own folder. For 1.0.0 this behaviour is preserved
and documented rather than refactored — an env-var-configurable `RESULTS_DIR`
would be a code change requiring its own revalidation, and centralising outputs
is not a stated company requirement. No empty `outputs/` directory is shipped,
because nothing would write to it.

## Bind-mount ownership

A baked-in fixed-UID user is **worse** than root for bind mounts: container UID
1000 writing to a host directory owned by UID 1001 gets permission denied, so
notebooks cannot save at all, whereas root always writes (files are merely
root-owned). But `--user "$(id -u):$(id -g)"` alone also fails — verified:

```
$ docker run --rm --user 4242:4242 python:3.9-slim sh -c 'echo HOME=$HOME; touch $HOME/x'
whoami: cannot find name for user ID 4242
HOME=/
touch: cannot touch '//x': Permission denied
```

Jupyter needs a writable HOME for its runtime/config/data dirs. The image must
therefore be arbitrary-UID tolerant: create the user with `--gid 0`, `chmod -R g=u`
its home, point `JUPYTER_*_DIR` at `/tmp`, and run `USER 1000:0`. `start.sh` then
passes `--user "$(id -u):0"` on Linux. `--allow-root` is dropped, since root is.

Note this cannot be tested on the macOS dev host at all: Docker Desktop's
virtiofs layer rewrites ownership so a bind-mount permission test always appears
to pass. Arbitrary-UID tolerance removes the need to test it per host.

## Phase A changes (done)

1. Removed `from turtle import mode` from `calculations/defective_metal.py`.
   Confirmed dead — every `mode` reference in that file resolves to the local
   parameter `mode='both'` of `combined_microstructure_model`.
2. Removed the shadowing second import in `data/material_data.py` and
   `data/oxide_properties.py`. Each had
   `from ...model_config import METALS` on line 1 immediately overwritten by
   `from ...fuerst_etal_2024_model_config import METALS` on line 2, so
   `MATERIALS` resolved to Fuerst `Hastelloy_N` instead of the active study's
   `metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985`.

   **This never affected results.** Nothing in the delivery path imported the
   shims, so `Application/sa_results/*.csv` is not contaminated. It was a latent
   trap for the next caller, not an active defect.
3. Added `.ruff.toml` as a lint gate — see below.

All three changes verified behaviour-neutral: the ten `calculations` modules
still import, and every value in the baseline table above is bit-identical
before and after.

## Lint gate

`ruff check calculations data Application/*.ipynb` (ruff 0.16.4).

Deliberately narrow — only bug classes that can silently change which parameters
the model resolves. Style rules are excluded so they cannot bury that signal.

- `F811` redefinition — the exact shim-shadowing bug. Regression-tested: the gate
  exits 1 on a reintroduced duplicate `METALS` import.
- `TID251` banned imports (`turtle`, `tkinter`) — GUI modules absent from
  `python:*-slim`, so an accidental autocomplete import becomes an ImportError
  only *inside* the container, after handover. Regression-tested.
- `F821` undefined name, `F402` import shadowed by loop var, `F406` star-import
  outside module level, `F632` `is` against a literal. All currently clean.

`F401` (unused import) is **not** selected: 27 findings, nearly all notebook
re-imports and the intentional `noqa`'d star re-export in `model_config.py`.
`TID251` covers the one case that mattered without the noise.

`*.ipynb` is exempted from `F811` only. The delivery notebooks are written so
each cell runs standalone, so most cells re-import `np` / `plt` / `brentq` /
`apply_style`; ruff parses a notebook as a single module and reports all 19 as
redefinitions. That pattern is intentional.

Not gated, for later cleanup if wanted: `F841` unused local variable, 7 findings.

## Repo gotcha found while writing this

`.gitignore:19` ignores `docs/`, so this record lives at the repo root instead —
anything placed under `docs/` is untracked and would be lost on a fresh clone.
`data/` is likewise gitignored, though its files predate the rule and remain
tracked, which is why the Phase A shim fixes commit normally. `plot&docs/` is
also ignored.
