"""Shared fixtures for the test suite.

pytest runs from the repository root, which is not a study folder, and there is no
default study — so importing anything under ``calculations`` would raise before a test
could run. The session names a study explicitly, which is exactly what every other
root-run tool does (``docs/conf.py``, ``docs/_tools/*``, ``container/model_smoke_test.py``).

Tests that need a *different* study, or that exercise the resolution itself, must use a
subprocess: the study is bound into ``model_config`` at import time, so it cannot be
changed inside a running interpreter. ``run_py`` below is that subprocess.

Study folders are DISCOVERED, not listed, so a study added tomorrow is covered by every
parametrised test — including the one asserting the release gate knows about it.
``KNOWN_STUDIES`` is the deliberate exception: hardcoded expectations for the studies
that exist today, so that a silent change to their material data fails a test.
"""
import ast
import os
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
STUDY_FILENAME = "study_config.py"

# Set before any test imports calculations.*
os.environ.setdefault("MHI_STUDY", "Guo_etal_2025/study_config.py")
os.environ.setdefault("MHI_STUDY_QUIET", "1")


def _literal(tree, name):
    """The value assigned to a module-level ``name``, without importing the file."""
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return node.value
    return None


def _discover():
    """(folder, declared STUDY_NAME, first METALS key) for every study folder.

    Read by AST rather than by import: a study file is ~900 lines of module-level
    dict literals, and importing four of them into the collector is both slow and a
    way to make collection fail on an unrelated syntax error.
    """
    out = []
    for path in sorted(ROOT.glob(f"*/{STUDY_FILENAME}")):
        tree = ast.parse(path.read_text())

        name_node = _literal(tree, "STUDY_NAME")
        assert name_node is not None, f"{path} declares no STUDY_NAME"
        study_name = ast.literal_eval(name_node)

        metals = _literal(tree, "METALS")
        assert isinstance(metals, ast.Dict) and metals.keys, f"{path} has no METALS"
        metal = ast.literal_eval(metals.keys[0])

        out.append((path.parent.name, study_name, metal))
    return out


#: Discovered from the tree — grows automatically with a new study folder.
STUDIES = _discover()

#: Pinned expectations for the studies that exist today. Discovery alone would happily
#: agree with a file whose material data had been swapped; this would not.
KNOWN_STUDIES = {
    "Guo_etal_2025":     ("Guo_etal_2025_316L", "316L_Guo_2025"),
    "Fuerst_et_al_2024": ("fuerst_etal_2024_model_config", "Hastelloy_N_fuerst_2024"),
    "incoloy802_cr2o3":  ("incoloy802_cr2o3",
                          "metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985"),
    # Application holds its own copy of the 316L data and declares the same
    # STUDY_NAME, which is what lets the container gate find its reference set.
    "Application":       ("Guo_etal_2025_316L", "316L_Guo_2025"),
}


def run_py(code, cwd=None, env=None):
    """Run ``code`` in a fresh interpreter, and return the CompletedProcess.

    The parent's MHI_STUDY is removed first: the session sets one for in-process
    imports, and inheriting it would defeat every test of the resolution order.
    """
    e = dict(os.environ)
    e.pop("MHI_STUDY", None)
    e["MHI_STUDY_QUIET"] = "1"
    e["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT), e.get("PYTHONPATH", "")]).rstrip(os.pathsep)
    e.update(env or {})
    return subprocess.run([sys.executable, "-c", code],
                          cwd=str(cwd or ROOT), env=e,
                          capture_output=True, text=True)


@pytest.fixture(scope="session")
def root():
    return ROOT
