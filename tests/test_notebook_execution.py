"""Execute a notebook end to end, the way a reader actually runs one.

Every other notebook test in this suite reads notebook *source* — that the first cell
calls the locality guard, that no cell subscripts a retired config name, that figures are
saved through ``calculations.figures``. None of them runs a cell, so a notebook could
import a module that no longer exists, or divide by a value that has become zero, and the
whole suite would stay green. Both of those have happened: a quote-nesting fix produced a
``SyntaxError`` that only executing the notebook revealed.

This is deliberately one notebook, not all 24. ``Proposal.ipynb`` covers levels 1-5 and
is self-contained; the ``sensitivity_regime_*`` notebooks run full sampling scans, and the
``regime_parallel_coords*`` notebooks only read scan CSVs that are gitignored and cannot be
regenerated without running the scans first. One honest end-to-end run is worth more than
four that have to be skipped on a fresh clone.

It costs about 75 seconds. Marked ``slow`` so it can be deselected with
``pytest -m "not slow"``, but it runs by default — a test that closes this gap only when
someone remembers a flag does not close it.
"""
import json
import os
import pathlib
import subprocess
import sys
import time

import pytest

from conftest import ROOT

NOTEBOOK = "Proposal.ipynb"
STUDY = "Guo_etal_2025"
# Level 1 through Level 5, one figure per level with 2a and 2b separate.
EXPECTED_FIGURES = 6


def kernel_matching_this_interpreter():
    """The kernelspec whose interpreter is the one running pytest, or None.

    The notebooks declare ``"name": "python3"``, which is whatever kernel happens to be
    registered under that name on the machine — on this one it is a different conda
    environment without the package installed, so executing with the notebook's own
    kernelspec fails at ``import calculations``. The kernel has to be chosen to match
    the interpreter that has the package, not taken from the file.
    """
    try:
        out = subprocess.run(
            [sys.executable, "-m", "jupyter", "kernelspec", "list", "--json"],
            capture_output=True, text=True, timeout=120,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None

    want = pathlib.Path(sys.executable).resolve()
    for name, spec in json.loads(out.stdout).get("kernelspecs", {}).items():
        argv = spec.get("spec", {}).get("argv") or []
        if argv and pathlib.Path(argv[0]).resolve() == want:
            return name
    return None


@pytest.mark.slow
def test_proposal_notebook_runs_and_draws_its_figures(tmp_path):
    kernel = kernel_matching_this_interpreter()
    if kernel is None:
        pytest.skip(
            f"no Jupyter kernel registered for {sys.executable}; "
            "register one with `python -m ipykernel install --user --name <env>`"
        )

    study_dir = ROOT / STUDY
    results = study_dir / f"{pathlib.Path(NOTEBOOK).stem}-results"
    started = time.time()

    # conftest sets MHI_STUDY for the in-process tests, but it must not leak into this
    # one. The notebooks open with assert_study_is_local(), which exists to insist the
    # study came from the notebook's own folder rather than from the environment — so
    # inheriting the variable makes the run fail the very guard it should satisfy.
    # Dropping it is also what makes this a real test of resolution-by-location.
    env = {k: v for k, v in os.environ.items() if k != "MHI_STUDY"}
    env["MPLBACKEND"] = "Agg"

    # --output-dir keeps the executed copy in tmp_path, so the notebook committed to
    # the repository is never rewritten by running the tests.
    proc = subprocess.run(
        [sys.executable, "-m", "jupyter", "nbconvert",
         "--execute", "--to", "notebook", f"--output-dir={tmp_path}",
         f"--ExecutePreprocessor.kernel_name={kernel}",
         "--ExecutePreprocessor.timeout=1800",
         NOTEBOOK],
        cwd=study_dir, capture_output=True, text=True, timeout=3600, env=env,
    )
    assert proc.returncode == 0, (
        f"{STUDY}/{NOTEBOOK} failed to execute:\n{proc.stderr[-4000:]}"
    )

    executed = json.loads((tmp_path / NOTEBOOK).read_text())
    errors = [
        out for cell in executed["cells"]
        for out in cell.get("outputs", [])
        if out.get("output_type") == "error"
    ]
    assert not errors, "\n".join(
        f"{e.get('ename')}: {e.get('evalue')}" for e in errors
    )

    # The run has to have produced figures, and they have to be from *this* run — the
    # folder is gitignored but survives between runs, so stale PNGs would otherwise
    # satisfy a bare existence check.
    pngs = sorted(results.glob("*.png"))
    assert len(pngs) == EXPECTED_FIGURES, f"expected {EXPECTED_FIGURES}, got {[p.name for p in pngs]}"
    stale = [p.name for p in pngs if p.stat().st_mtime < started]
    assert not stale, f"figures not rewritten by this run: {stale}"
