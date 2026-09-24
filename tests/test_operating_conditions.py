"""`CONDITIONS` is the operating point; nothing else may restate it.

`CONDITIONS` defines the scenario the simulation is built on — state and geometry —
and `METALS`/`OXIDES` describe only what the material is. Nothing in `calculations/`
reads `CONDITIONS`, so that rule is a convention with no enforcement behind it; these
tests are the enforcement.

Each test corresponds to one fix in CHECKLIST.md and is added by it.
"""
import ast
import json
import re

import pytest

from conftest import ROOT, STUDIES, run_py

STUDY_IDS = [s[0] for s in STUDIES]

# Keys in the flat VALIDATION that must equal a CONDITIONS entry, not a literal.
VALIDATION_SHARED = {
    "P_up":   "P_upstream",
    "P_down": "P_downstream",
    "T":      "T_operating",
}

REPORT = """
import json
from calculations.config.model_config import CONDITIONS, VALIDATION
print(json.dumps({
    'conditions': {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in CONDITIONS.items()},
    'validation': {k: v for k, v in VALIDATION.items()
                   if isinstance(v, (int, float, str))},
    'validation_keys': sorted(VALIDATION),
    'validation_nested': sorted(k for k, v in VALIDATION.items() if isinstance(v, dict)),
}))
"""


def config_of(folder):
    proc = run_py(REPORT, cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


# --------------------------------------------------------------------------
# Fix 1 — VALIDATION is flat, and its operating point comes from CONDITIONS
# --------------------------------------------------------------------------

@pytest.mark.parametrize("folder", STUDY_IDS)
def test_validation_operating_point_comes_from_conditions(folder):
    cfg = config_of(folder)
    for val_key, cond_key in VALIDATION_SHARED.items():
        assert cfg["validation"][val_key] == cfg["conditions"][cond_key], (
            f"VALIDATION[{val_key!r}] restates the operating point instead of reading "
            f"CONDITIONS[{cond_key!r}]")
    assert cfg["validation"]["T_min"] == cfg["conditions"]["T_range"][0]
    assert cfg["validation"]["T_max"] == cfg["conditions"]["T_range"][1]


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_validation_is_flat(folder):
    """No per-level sub-dicts: one namespace, level-prefixed where a key is specific.

    The nesting is what allowed the same pressure to be called P_ref in one level,
    P_fixed in four others and P_sensitivity in a sixth.
    """
    nested = config_of(folder)["validation_nested"]
    assert not nested, f"VALIDATION still has per-level sub-dicts: {nested}"


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_validation_has_no_duplicate_names_for_one_quantity(folder):
    """The retired names must not come back under a level prefix either."""
    retired = {"P_ref", "P_fixed", "P_fixed_metal_sweep", "P_fixed_oxide_sweep",
               "P_sensitivity", "T_mode_comparison", "T_grain_size", "T_comparison",
               "T_limit"}
    offenders = [k for k in config_of(folder)["validation_keys"]
                 if k in retired or k.split("_", 1)[-1] in retired]
    assert not offenders, f"retired VALIDATION names present: {offenders}"


# --------------------------------------------------------------------------
# the notebooks that read it
# --------------------------------------------------------------------------

NOTEBOOKS = sorted(ROOT.glob("*/Proposal.ipynb")) + sorted(ROOT.glob("*/Surface_proposal.ipynb"))


def validation_keys_read(nb_path):
    """Every literal key a notebook subscripts VALIDATION with, via the AST.

    Text inside docstrings mentions `VALIDATION['L3_*']` and must not count — that is
    prose about the config, not a read of it.
    """
    nb = json.loads(nb_path.read_text())
    keys = []
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        for node in ast.walk(ast.parse("".join(cell["source"]))):
            if (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
                    and node.value.id == "VALIDATION"
                    and isinstance(node.slice, ast.Constant)):
                keys.append(node.slice.value)
    return keys


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_notebooks_only_read_keys_that_exist(nb_path):
    """A wrong key is a KeyError minutes into a run; this finds it in milliseconds."""
    defined = set(config_of(nb_path.parent.name)["validation_keys"])
    unknown = sorted(set(validation_keys_read(nb_path)) - defined)
    assert not unknown, f"reads VALIDATION keys that do not exist: {unknown}"


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_notebooks_do_not_reach_through_a_level_alias(nb_path):
    """`VAL_L4 = VALIDATION['L4']` kept the old shape alive one indirection away."""
    nb = json.loads(nb_path.read_text())
    code = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    for alias in ("VAL_L3", "VAL_L4"):
        assert not re.search(rf"\b{alias}\b", code), f"{alias} alias still present"


def validation_value_nodes(folder):
    """The AST node assigned to each VALIDATION key in a study file.

    Comparing *values* cannot tell a CONDITIONS reference from a literal that happens
    to equal it — `'P_up': 1e5` and `'P_up': CONDITIONS['P_upstream']` are both 1e5
    today, and the first silently stops tracking the config the moment it changes.
    So the derivation has to be checked in the source, not in the result.
    """
    tree = ast.parse((ROOT / folder / "study_config.py").read_text())
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == "VALIDATION"
                        for t in node.targets)):
            return {k.value: v for k, v in zip(node.value.keys, node.value.values)
                    if isinstance(k, ast.Constant)}
    raise AssertionError(f"{folder}/study_config.py defines no VALIDATION")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_operating_point_is_written_as_a_reference_not_a_literal(folder):
    nodes = validation_value_nodes(folder)
    for key in list(VALIDATION_SHARED) + ["T_min", "T_max"]:
        value = nodes[key]
        refs_conditions = any(isinstance(n, ast.Name) and n.id == "CONDITIONS"
                              for n in ast.walk(value))
        assert refs_conditions, (
            f"VALIDATION[{key!r}] is written as a literal. It must read CONDITIONS, or "
            f"it stops tracking the operating point the moment CONDITIONS changes.")
