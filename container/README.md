# Hydrogen Permeation Model — Offline Analysis Environment

A self-contained JupyterLab environment for the MHI hydrogen-permeation model.
Everything needed to run the analysis is inside the container image. **No Python
packages are installed, and no internet access is required, at any point.**

Verified to run with networking fully disabled.

---

## What you received

```
Hydrogen_Model_Environment_1.2.0/
├── image/
│   ├── hydrogen-model-1.2.0.tar      the container image
│   └── hydrogen-model-1.2.0.sha256   checksum for the above
├── scripts/
│   ├── verify.sh                     run the two acceptance gates
│   └── start.sh                      launch JupyterLab
├── workspace/                         your working files (see below)
├── requirements.in                    what we deliberately depend on
├── requirements.lock.txt              the exact 106 packages the image installs
├── installed-packages.txt             pip freeze taken from inside the built
│                                      image; byte-identical to the lock, which
│                                      is the evidence the build honoured it
├── Dockerfile                         how the image was built
└── README.md
```

## Requirements

- A Docker-compatible container runtime (Docker Engine, Docker Desktop, or
  Podman in Docker-compatible mode) able to run **Linux `amd64`** containers.
- ~3 GB free disk for the loaded image.
- Permission to bind-mount a local directory and to bind a port on `127.0.0.1`.

The container runs as a **non-root** user and needs no elevated privileges.

## First-time setup

**1. Check the archive is intact.** From the package root:

```bash
shasum -a 256 -c image/hydrogen-model-1.2.0.sha256
```
Expect `image/hydrogen-model-1.2.0.tar: OK`. On Linux without `shasum`, use
`sha256sum -c`. If this fails, the transfer was corrupted — do not proceed.

**2. Load the image.**

```bash
docker load -i image/hydrogen-model-1.2.0.tar
docker images hydrogen-model
```

**3. Verify it.**

```bash
./scripts/verify.sh
```
This runs with `--network none` and must end with:

```
VERIFICATION PASSED — hydrogen-model:1.2.0 is the validated environment.
```

It checks two separate things: that all 106 pinned packages are exactly the
validated versions, and that the model still reproduces known reference values
through its real solvers. If either fails, do not use the image for analysis.

## Daily use

```bash
./scripts/start.sh
```

Then **copy the full URL printed in the terminal, including the `?token=...`
part**, into your browser. Without the token the page will not load. The token
is new each time the server starts.

Press `Ctrl-C` in the terminal to stop.

## The workspace

```
workspace/
├── calculations/          the model. Python package; do not rename.
├── Guo_etal_2025/         316L / Cr2O3 (Guo et al. 2025)
├── Fuerst_et_al_2024/     Hastelloy N / Cr2O3 (Fuerst et al. 2024)
└── incoloy802_cr2o3/      Incoloy 802 / Cr2O3
```

Each study folder is self-contained and holds the same six notebooks plus its own
input file:

```
<study>/
├── study_config.py                 THE INPUT FILE — every number this study runs on
├── Proposal.ipynb                  levels 1-5, self-contained
├── Surface_proposal.ipynb          surface kinetics (level 6), self-contained
├── sensitivity_regime_L5.ipynb     regime sensitivity analysis
├── sensitivity_regime_L5L6.ipynb   regime sensitivity incl. surface
├── regime_parallel_coords.ipynb    visualises sensitivity_regime_L5L6 output
├── regime_parallel_coords_L5.ipynb visualises sensitivity_regime_L5 output
├── sa_results/                     scan CSVs land here (starts empty)
└── sa_results_L5/                  scan CSVs land here (starts empty)
```

A notebook picks up the `study_config.py` in its own folder by virtue of running
there — there is no switch to set. Edit that file and you change that study alone.
To add a study, copy a folder.

Anything you put in `workspace/` is visible inside JupyterLab, and anything the
notebooks write appears on your own filesystem. Figures are written into a folder
per notebook beside it — `Proposal.ipynb` writes `Proposal-results/` — not to a
central output directory.

### Run order matters

`Proposal.ipynb` and `Surface_proposal.ipynb` are self-contained — run them in
any order.

The other four are a pipeline. The two `regime_parallel_coords*` notebooks only
*plot* results; they read CSVs produced by the sensitivity notebooks and cannot
regenerate them. **No precomputed results are shipped**, so run the producer
first:

```
sensitivity_regime_L5.ipynb    ──>  regime_parallel_coords_L5.ipynb
sensitivity_regime_L5L6.ipynb  ──>  regime_parallel_coords.ipynb
```

Running a `parallel_coords` notebook first fails with
`FileNotFoundError: sa_results_L5/master_clusters.csv`. That is expected, not a
fault — run its producer.

The sensitivity notebooks run a large parameter scan and take a while. Once they
have run, their CSVs are cached in `sa_results*/scans/`, and re-running the plots
is fast.

## Adding a Python package

Don't `pip install` inside a running container — the change is lost when the
container stops, and it silently breaks the guarantee that `verify.sh` checks.

Instead, tell us what you need. We add it to `requirements.in`, regenerate the
lock, rebuild, re-run both gates, and issue a new version. That is the whole
point of the pinned lock: the environment that produced a result can be
reconstructed exactly.

## Reporting a problem

Please include:

- the output of `./scripts/verify.sh` in full
- `docker --version` and your OS
- `docker image inspect hydrogen-model:1.2.0 --format '{{.Id}} {{.Architecture}}'`
- which notebook and which cell, with the error text

## Version

`hydrogen-model:1.2.0` — Python 3.12, linux/amd64. Cite this tag alongside any
results, so they can be tied back to a known environment.

### Changes since 1.1.0

The environment is unchanged again — same Dockerfile, same 106-package lock — so a
result computed under 1.1.0 runs in the same interpreter here. What changed is the
model, and two of the changes move numbers:

- **There is no `ACTIVE_STUDY` switch any more, and no `calculations/config/studies/`.**
  Each study folder holds its own `study_config.py` — the real input file — and it is
  selected by running in that folder, or by naming one in `$MHI_STUDY`. Running from
  the workspace root, which belongs to no study, is an error rather than a silent
  default.
- **1.1.0 results will not reproduce under 1.2.0.** Two changes are responsible. The
  operating point in `CONDITIONS` is now the authority for the wall being simulated,
  so the oxide layer is its 1 µm `L_oxide` rather than the 48 nm coupon thickness
  recorded on the material — which drops the Level 5L6 flux about threefold and
  flips that level from metal-limited to oxide-limited. And `permeability` is now an
  *apparent* permeability derived from the flux each level actually solved, rather
  than a harmonic mean of two bulk material constants: a factor of 775 at the default
  operating point, and it no longer reports the same value for two levels whose
  fluxes differ.
- **Figures are written into a folder per notebook**, `<notebook>-results/`, rather
  than beside the notebook.
- **Three studies ship, one folder each.** `Guo_etal_2025/`, `Fuerst_et_al_2024/`
  and `incoloy802_cr2o3/`, each self-contained. A fourth folder, `Application/`,
  was withdrawn on 2026-09-28: its `study_config.py` was a delivery copy of the Guo
  study — identical bar the docstring, same `STUDY_NAME` — so it advertised a fourth
  study that reproduced a third one's numbers. Use `Guo_etal_2025/` for 316L / Cr2O3.
- **The derived parameter block lives in one place.** It used to be duplicated
  verbatim in every study file; it is now in `calculations/config/model_config.py`,
  so a study file holds only its own data.

Cite the image tag with any result computed from it, which is what makes this
distinction recoverable later.

### Changes since 1.0.1

The environment is byte-for-byte the same specification — identical Dockerfile,
identical 106-package lock, so any result computed under 1.0.1 reproduces here.
What changed is the model and the workspace:

- **Multiple studies, one switch.** `calculations/config/studies/` now holds three
  studies — Incoloy 802 / Cr2O3, Hastelloy N (Fuerst 2024), and 316L (Guo 2025).
  `ACTIVE_STUDY` in `calculations/config/model_config.py` selects one, and every
  module and notebook follows it. Nothing else needs editing.
  *(True of the 1.1.0 image. Superseded in 1.2.0 — see "Changes since 1.1.0" above.)*
- **A study contract.** Selecting a study that is missing a required export now
  fails immediately, naming the study and the missing names, instead of raising
  deep inside an import.
- **Sensitivity draw counts are measured, not cached.** The analysis probes each
  regime preset and sizes the run from the yield it observes, so it stays correct
  when ranges, presets, material or temperature change.
- **`scripts/verify.sh` gate 2 is study-aware.** It holds reference values per
  study and checks the active one. A study with no measured reference set fails
  rather than passing silently.

### Changes since 1.0.0

1.0.0 required `--user "$(id -u):0"`; running it as `--user "$(id -u):$(id -g)"`
— the more natural form — failed with
`PermissionError: [Errno 13] Permission denied: '/home/appuser/.ipython'`.
1.0.1 redirects every per-user state directory (IPython, matplotlib, XDG) to
`/tmp`, so any `uid:gid` works and your files are written with your own group
rather than group 0. `scripts/start.sh` uses the natural form. Nothing about the
scientific environment changed — the 106 locked package versions and every
reference value are identical.
