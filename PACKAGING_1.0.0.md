# Packaging decisions — `hydrogen-model:1.0.0` offline container

Scope record for the offline Docker deliverable. Captures decisions that are not
recoverable from the code, so a later reader does not have to re-derive them.

Status: **`hydrogen-model:1.0.0` released 2026-08-24.** Built, both gates passing
offline as non-root, notebooks executed, exported, checksummed, and verified by
deleting the image and restoring it from the archive alone.

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

### Workspace contents — decided: code and notebooks only

All generated artefacts are stripped from the shipped workspace. Decided
2026-08-24: drop every `.html`, `.csv` and `.png`, keeping the directory
skeleton (with `.gitkeep`) because the notebooks glob `FIG_DIR` and write into
`RESULTS_DIR`.

| | before | after |
|---|---|---|
| workspace total | 121 MB | **13 MB** |
| `.html` figures | 15 files, 71.8 MB | 0 |
| `.csv` results + scans | 43 files, 36.9 MB | 0 |
| `.py` (calculations) | 14 files, 0.5 MB | unchanged |
| `.ipynb` | 6 files, 12.3 MB | unchanged |

There were no `.png` files in the workspace to begin with — the figures are
plotly HTML.

**Consequence, accepted deliberately.** `regime_parallel_coords.ipynb` and
`regime_parallel_coords_L5.ipynb` read `master_clusters.csv`,
`routeB_givendata.csv` and `compare_delta_*.csv`, and cannot regenerate them.
They are therefore not runnable out of the box; the company must run
`sensitivity_regime_L5.ipynb` / `sensitivity_regime_L5L6.ipynb` first to produce
the scans and cluster tables. That ordering must be stated in the README. The
upside is that no possibly-stale results ship, and the deliverable is code plus
notebooks only.

Gate 2's "shipped SA artefacts" section detects the absent CSVs and reports
SKIP rather than failing, so it passes against both the development tree (where
the CSVs exist and rankings are asserted) and the shipped workspace.

Note the remaining 13 MB is almost entirely embedded notebook output — the six
`.ipynb` files are 12.3 MB of base64 matplotlib PNGs. Stripping outputs would
take the workspace under 1 MB, but the outputs are useful as documentation of
expected results, so they are kept.

### Notebook execution — measured 2026-08-24

All six executed headless in rc1 via `jupyter nbconvert --execute`, offline
(`--network none`), non-root (`--user "$(id -u):0"`), against the stripped
workspace, with a 600 s per-cell limit.

| notebook | result |
|---|---|
| `Proposal.ipynb` | **passed**, 41 s — wrote all six `Level*.png` into the bind mount |
| `Surface_proposal.ipynb` | **passed**, 43 s |
| `regime_parallel_coords.ipynb` | failed, 10 s — `FileNotFoundError` on the stripped CSVs (expected) |
| `regime_parallel_coords_L5.ipynb` | failed, 10 s — same |
| `sensitivity_regime_L5.ipynb` | `CellTimeoutError` at 600 s, in the scan-regeneration cell |
| `sensitivity_regime_L5L6.ipynb` | `CellTimeoutError` at 600 s, in the `givendata_sensitivity_by_regime` cell |

Neither timeout is an environment fault — both are the 600 s per-cell limit
being unrealistic for these notebooks. But they fail in **different** cells, and
the distinction matters for the shipping decision:

- `sensitivity_regime_L5` timed out inside the cell that runs `cache_is_valid`
  and writes `master_clusters.csv`, i.e. **scan regeneration**. Shipping
  `sa_results_L5/scans/` (21 MB) would let its cache check short-circuit this.
- `sensitivity_regime_L5L6` got **past** scan regeneration — it successfully wrote
  `scans/scan_{metal,oxide,surface}.csv` and `master_clusters.csv` — and then timed
  out in the PAWN/delta analysis itself. Shipping cached scans would **not** help
  it; that cost is the sensitivity computation, not the scans.

So reinstating `scans/` would speed up one of the two, not both. Left stripped as
instructed; revisit only if first-run wall time becomes a company complaint.

Practical consequence: do not use a 600 s per-cell limit for any automated
notebook gate on the sensitivity notebooks. `Proposal` and `Surface_proposal`
are fast (~40 s) and are the sensible smoke-level notebook gate.

### Archive round-trip — PASSED 2026-08-24

The full company-side restore, performed after **deleting the image locally**.
Both tags (`1.0.0` and `1.0.0-rc1`) were removed, not just one — untagging a
single tag leaves the layers resident and makes `docker load` a no-op, so the
test would prove nothing. `docker rmi` reported `Deleted: sha256:f393ab...`,
confirming the layers were actually freed.

| step | result |
|---|---|
| `docker save` | 228 MB tar (image is 1.1 GB uncompressed) |
| `shasum -a 256 -c` | `image/hydrogen-model-1.0.0.tar: OK` |
| `docker rmi` both tags | image gone |
| `docker load -i` | `Loaded image: hydrogen-model:1.0.0` |
| image ID before vs after | `sha256:f393ab4bd5e6...` **identical** |
| `scripts/verify.sh` from inside the package | both gates PASSED, offline |
| `scripts/start.sh`, host port | 403 unauthenticated, 200 with token, all six notebooks served, bound `127.0.0.1:8888` only |
| `Proposal.ipynb` executed from the restored image | passed, 43 s, wrote six PNGs owned by the host user |

Image identity: `sha256:f393ab4bd5e6cf7af1516b3ae2371b6151651aa1d1a53d63c99c08fc1cfb93b6`,
`amd64/linux`, `User=1000:0`.

`installed-packages.txt` (pip freeze from inside the built image) is
byte-identical to `requirements.lock.txt`. That equality is the audit evidence
that the build installed exactly what was locked.

Note the `.sha256` file records the path `image/hydrogen-model-1.0.0.tar`
relative to the package root, so `shasum -c` must be run from there — which is
what the README instructs.

### Two script bugs found by testing, both fixed

1. `verify.sh` mounted a non-existent `workspace/`, and **Docker silently creates
   a missing bind-mount source as an empty directory**. A broken package layout
   therefore surfaced as a confusing gate-2 import failure. It now checks for
   `workspace/calculations` first and says the package is incomplete.
2. `start.sh` used `docker run -it` unconditionally and died with `the input
   device is not a TTY` whenever stdin/stdout were redirected
   (`./scripts/start.sh > jupyter.log`, CI, a GUI launcher). Now guarded by
   `[ -t 0 ] && [ -t 1 ]`. The first attempt at that fix used a bash array and
   broke on **macOS's bash 3.2**, where expanding an empty array under `set -u`
   aborts with `TTY_FLAGS[@]: unbound variable`; it uses a scalar instead.

### Still open

- **True Linux bind-mount ownership.** Cannot be tested on macOS at all: Docker
  Desktop's virtiofs rewrites ownership, so the mount test always appears to
  pass. The arbitrary-UID design is what makes it work, and the half that *can*
  be proven here — that an unknown UID gets a writable HOME — is proven.
- **`gb_enhancement_factor` is a dead parameter.** `calculate_gb_enhanced_diffusivity`
  takes no such argument and derives alpha internally from a hardcoded
  austenitic-steel table, so the `'gb_enhancement_factor': 100` that
  `sensitivity.py` puts into `microstructure_params` is never read (verified:
  alpha of 1, 100 and 10000 give an identical `D_eff`). Harmless today — it is
  not among the 35 parameters the shipped SA actually varies — but if it is ever
  added to the SA ranges it will report delta ~ 0 and read as physics rather
  than as a wiring bug.
- **Deployment questions for the company**, none of which are Python questions:
  container runtime and version, whether Linux containers are permitted, whether
  bind mounts and localhost port binds are allowed, whether images are scanned
  before import, and whether containers must run as a specific UID.
