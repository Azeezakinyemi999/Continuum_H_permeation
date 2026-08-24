# Packaging decisions — `hydrogen-model:1.0.0` offline container

Scope record for the offline Docker deliverable. Captures decisions that are not
recoverable from the code, so a later reader does not have to re-derive them.

Status: **Phase A, C, D, E complete.** `hydrogen-model:1.0.0-rc1` builds and both
gates pass offline as non-root. Not yet exported / notebook-executed / promoted.

## Target environment

| Decision | Value | Why |
|---|---|---|
| Python | **3.12** (`3.12.14`) | 3.9 reached end of life in Oct 2025 and the company is expected to vulnerability-scan the image. The numpy<2 ceiling that forced 3.9 in `mace_env` comes from **torch**, which this project never imports, so the constraint does not apply here. Requires a numerical revalidation gate (below). |
| Platform | `linux/amd64` | Set at build time (`docker build --platform linux/amd64`), not in `FROM`, so the image definition stays reusable. Dev host is already x86_64/amd64. |
| Base image | `python:3.12-slim` | |
| Release tag | `hydrogen-model:1.0.0` | Promoted from `1.0.0-rc*` only after every acceptance gate passes. Never `latest`. |
| Container user | non-root, arbitrary-UID tolerant | See "Bind-mount ownership". |

### Numerical revalidation gate — PASSED 2026-08-24

Moving 3.9 -> 3.12 changes the interpreter under the solvers, so reference values
were re-measured rather than assumed. `container/revalidate_reference.py` was run
under both interpreters and the output diffed.

| | |
|---|---|
| reference | Python 3.9.23, macOS, x86_64 (`mace_env`) |
| candidate | Python 3.12.14, linux/amd64, from `container/requirements.lock.txt` |
| result | **bit-identical on every model output** |

Covered: `arrhenius`; L1 `calculate_simple_metal_flux`; L2b
`calculate_oxide_metal_system` (brentq); L1L6 `solve_steady_state_flux_L1L6`
(brentq on a different equation); `combined_microstructure_model` (4 trap types);
seeded SALib Latin-hypercube -> PAWN + Borgonovo delta; pandas CSV round-trip;
and `top_drivers` / `parallel_coordinates_*` against the real
`Application/sa_results/` CSVs.

Reference values, identical under both:

| quantity | value |
|---|---|
| `arrhenius(1e-11, 50000, 800, 700)` | `2.926830409475082e-11` |
| L1 flux | `1.955109516956673e-07` |
| L2b flux | `1.888003770552079e-07` |
| L2b `P_interface` | `93462.901558976` |
| L2b `resistance_ratio` | `0.036154225782795736` |
| L1L6 `J_ss` | `7.10792976386876e-07` |
| L1L6 `theta` | `0.7597385771009892` |
| L1L6 `P_int` | `99990.85191948501` |
| micro `D_eff` | `5.084388686889583e-11` |
| micro `overall_factor` | `0.7388612272380938` |
| micro `theta_total` | `0.026370841516212183` |
| SALib `pawn_median` | `[0.17518028846153846, 0.14708533653846154, 0.6875]` |
| SALib `delta` | `[0.10706087473576083, 0.06905348258240432, 0.6947500046510651]` |
| `classify_regime_level14(0.2)` | `metal_limited/traps_defect_limited` |
| `classify_regime_level14(0.9)` | `metal_limited/lattice_limited` |

`top_drivers` rankings are identical for all three regimes on both.

**The one difference found**, and why it does not block: the plotly figure from
`parallel_coordinates_samples` differs in 23 of 1495 values of its `line.color`
channel, by a maximum **relative** difference of `2.1e-16` — one ULP of a
float64 (worst case `-4.914129070710655` vs `-4.914129070710656`). That channel
is `log10(flux)` computed for colour mapping only; the parcoords dimension values
are bit-identical. This is macOS libm vs glibc rounding in `np.log10`, has no
scientific content, and is ~7 orders of magnitude inside the `rtol=1e-9`
tolerance the smoke test should use. Do not compare figure HTML byte-for-byte:
it also differs by 13-23 bytes purely from the embedded plotly version string.

Two API details for whoever writes `model_smoke_test.py`:
`get_metal_properties_at_T` returns `D_metal` / `K_s_metal` (not `D` / `K_s`),
and `calculate_oxide_metal_system` needs a `thickness` key injected into **both**
props dicts — the T-evaluated getters do not supply it for the metal.
`solve_steady_state_flux_L1L6` returns `J_ss` / `P_int` / `theta` / `beta` /
`rate_limiting` — there is no `flux` key.

## Dependencies

`container/requirements.in` is human-maintained intent; `container/requirements.lock.txt`
(106 packages) is what the build installs. Regenerate with `container/lock.sh`,
which resolves *inside* `python:3.12-slim` on `linux/amd64` — resolving on the
macOS dev host would bake in the wrong platform wheels.

Pinned in `requirements.in` because results depend on them: `numpy==1.26.4`,
`scipy==1.13.1`, `pandas==2.3.3`, `SALib==1.5.0`. SALib must stay at 1.5.0
because 1.5.1+ calls `np.trapezoid`, which needs numpy>=2.

Left floating, then frozen by the lock, because they are presentation-only. These
had been pinned to the last Python 3.9-compatible releases, so on 3.12 the old
pins were arbitrary:

| package | was (3.9) | now (3.12 lock) |
|---|---|---|
| matplotlib | 3.9.4 | 3.11.1 |
| plotly | 6.3.1 | 6.9.0 |
| ipython | 8.18.1 | 9.16.1 |
| jupyterlab | 4.4.9 | 4.6.3 |

All four verified on 3.12: headless Agg `savefig` with `Rectangle` patches,
`go.Parcoords().to_html()`, and `from IPython.display import display, HTML`.
IPython 8 -> 9 is a major bump but the only IPython API used is
`IPython.display.display`, which is unchanged. Scanned the delivery notebooks for
matplotlib APIs removed in 3.10/3.11 (`cm.get_cmap`, `register_cmap`, `normed=`)
— none present.

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


## Gates (Phase D/E) — both passing on 1.0.0-rc1

`container/verify_environment.py` (gate 1) and `container/model_smoke_test.py`
(gate 2), both baked into `/opt` of the image.

Gate 1 checks the **whole** locked environment, not a hand-maintained subset: it
reads the same `requirements.lock.txt` the image was built from and compares all
106 pins against `importlib.metadata.version`, so the gate cannot drift out of
step with the build. It also asserts Python 3.12.x, runs functional checks
(brentq, SALib imports, headless Agg savefig with Rectangle patches, plotly
`Parcoords.to_html`, `IPython.display`, exact pandas CSV round-trip), and asserts
`turtle` and `tkinter` are **absent** — the runtime mirror of the ruff TID251 rule.
Confirmed: they are absent, so the old `from turtle import mode` would indeed have
failed only inside the container.

Gate 2 asserts the reference-value table above through the real solvers at
`rtol=1e-9`, plus `top_drivers` rankings against the shipped
`Application/sa_results/` CSVs. It does not assert plotly HTML byte counts, only
that figures build — the counts move with the embedded plotly version string.

Verified on rc1, all with `--network none` and `--user "$(id -u):0"`:

| check | result |
|---|---|
| gate 1 offline | PASSED, 106/106 versions matched |
| gate 2 offline against assembled workspace | PASSED |
| HOME and JUPYTER_RUNTIME_DIR writable as unknown UID | yes |
| write PNG + CSV into the bind mount | yes, files owned by the host user |
| JupyterLab starts offline, prints tokenised URL | yes |
| authenticated `api/status` and `api/contents` | HTTP 200, all six notebooks served |

### The sys.path trap, third occurrence

`python /opt/model_smoke_test.py -w /workspace` failed with
`ModuleNotFoundError: No module named 'calculations'`, because running a script
puts the **script's** directory on `sys.path[0]`, not the cwd. This is the same
defect as the original exercise's `python calculations/sensitivity.py`, and it has
now bitten three times. Fixed systemically with `ENV PYTHONPATH=/workspace` in the
Dockerfile, plus a defensive `sys.path` bootstrap in the smoke test. The notebooks
are unaffected — they do their own `sys.path.insert`.

### Workspace size — needs a decision

The assembled `release/workspace` is **121 MB**:

| | size | is it needed? |
|---|---|---|
| `calculations/` | 508 KB | yes |
| `sa_results*/scans/` | 21 MB | yes — cached model evaluations, `load_regime_scans` reads them, regenerating means re-running the expensive scans |
| `sa_results*/figures/` | **72 MB** | probably not — generated HTML/PNG outputs. The notebooks *write* here and glob it; they do not need pre-existing contents. |
| `sa_results*/*.csv` | ~28 MB | yes — `master_clusters.csv`, `routeB_givendata.csv`, `compare_*` are read by the parallel-coords notebooks and cannot be regenerated without re-running the scans |

Dropping `figures/` (keeping the empty directories, since the notebooks glob
them) would take the workspace to ~49 MB. Not done yet — pending confirmation
that no notebook depends on a pre-existing figure.

### Not yet verified

Notebook **execution** end-to-end (only import/serve is proven), `docker save` /
checksum / `docker load` round-trip, and true Linux bind-mount ownership. The
last cannot be tested on macOS at all: Docker Desktop's virtiofs rewrites
ownership, so the mount test always appears to pass. The arbitrary-UID design is
what makes it work; the HOME-writability half of it *is* proven, since that is
inside the container filesystem rather than the mount.
