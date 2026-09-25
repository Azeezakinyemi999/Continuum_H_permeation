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


# --------------------------------------------------------------------------
# Fix 2 — the container gate runs at the study's operating point
# --------------------------------------------------------------------------

GATE = ROOT / "container" / "model_smoke_test.py"


def gate_assignments():
    """Module-level names the gate assigns, as AST nodes."""
    tree = ast.parse(GATE.read_text())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    out[target.id] = node.value
    return out


def test_gate_takes_its_operating_point_from_conditions():
    """It used to hardcode P_DOWN = 1.0e2 while CONDITIONS said 0.

    A gate that pins values measured at its own operating point certifies a scenario
    nobody is studying — and the disagreement was invisible, because the pinned numbers
    were measured with the same wrong pressure they were then compared against.
    """
    assigns = gate_assignments()
    for name in ("P_UP", "P_DOWN", "L_METAL"):
        assert name in assigns, f"{name} is no longer assigned at module level"
        refs = any(isinstance(n, ast.Name) and n.id == "CONDITIONS"
                   for n in ast.walk(assigns[name]))
        assert refs, f"{name} is a literal; it must read CONDITIONS"


def test_gate_temperature_stays_an_explicit_test_choice():
    """T is deliberately NOT the study's T_operating: 700 K sits inside the oxide's
    validated [473, 773] K window, so the gate is not asserting on an extrapolation."""
    value = gate_assignments()["T"]
    assert isinstance(value, ast.Constant) and value.value == 700.0
    assert "TEST choice" in GATE.read_text()


# --------------------------------------------------------------------------
# Fix 3 — CONDITIONS owns the geometry
# --------------------------------------------------------------------------

PARAMS_REPORT = """
import json
from calculations.config.model_config import (
    CONDITIONS, OXIDES, DEFAULT_PARAMS_LEVEL5, DEFAULT_PARAMS_LEVEL5L6,
    build_simulation_config)
SIM = build_simulation_config()
print(json.dumps({
    'conditions':  {k: (list(v) if isinstance(v, tuple) else v)
                    for k, v in CONDITIONS.items()},
    'l5':    {k: DEFAULT_PARAMS_LEVEL5[k] for k in
              ('metal_thickness', 'oxide_thickness', 'temperature',
               'P_upstream', 'P_downstream')},
    'l5l6':  {k: DEFAULT_PARAMS_LEVEL5L6[k] for k in
              ('metal_thickness', 'oxide_thickness', 'temperature',
               'P_upstream', 'P_downstream')},
    'sim':   {'L_oxide': SIM['L_oxide'], 'L_metal': SIM['L_metal']},
    'oxide_keys': sorted(next(iter(OXIDES.values()))),
}))
"""

GEOMETRY = {
    'metal_thickness': 'L_metal',
    'oxide_thickness': 'L_oxide',
    'temperature':     'T_operating',
    'P_upstream':      'P_upstream',
    'P_downstream':    'P_downstream',
}


def params_of(folder):
    proc = run_py(PARAMS_REPORT, cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_geometry_and_operating_point_come_from_conditions(folder):
    """Both layers follow the same rule. They did not: metal thickness read CONDITIONS
    while oxide thickness read the material dict, so the wall was described in two
    places and the material won."""
    cfg = params_of(folder)
    for level in ("l5", "l5l6"):
        for param, cond_key in GEOMETRY.items():
            assert cfg[level][param] == cfg["conditions"][cond_key], (
                f"DEFAULT_PARAMS_{level.upper()}[{param!r}] != CONDITIONS[{cond_key!r}]")
    assert cfg["sim"]["L_oxide"] == cfg["conditions"]["L_oxide"]
    assert cfg["sim"]["L_metal"] == cfg["conditions"]["L_metal"]


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_oxide_carries_no_live_thickness_key(folder):
    """`thickness` is renamed to `sample_thickness` so nothing can fall back to it.

    The old key was read by four places — the flat params, build_simulation_config,
    oxide_permeation and surface_kinetics — each of which looked like it was asking the
    material how thick the wall is.
    """
    keys = params_of(folder)["oxide_keys"]
    assert "thickness" not in keys, "OXIDES still has a live 'thickness' key"
    assert "sample_thickness" in keys, "the sample's measured thickness was lost"


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_oxide_thickness_is_written_as_a_conditions_reference(folder):
    """Value equality would pass a literal 1e-6; the derivation has to be in the source."""
    tree = ast.parse((ROOT / folder / "study_config.py").read_text())
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == "DEFAULT_PARAMS_LEVEL5"
                        for t in node.targets)):
            entries = {k.value: v for k, v in zip(node.value.keys, node.value.values)
                       if isinstance(k, ast.Constant)}
            value = entries["oxide_thickness"]
            assert any(isinstance(n, ast.Name) and n.id == "CONDITIONS"
                       for n in ast.walk(value)), \
                "oxide_thickness is a literal; it must read CONDITIONS"
            return
    raise AssertionError("DEFAULT_PARAMS_LEVEL5 not found")


# --------------------------------------------------------------------------
# Fix 3b — the study file is sectioned, in dependency order
# --------------------------------------------------------------------------

SECTION_TITLES = [
    "1. IDENTITY AND CONSTANTS",
    "2. WHAT IS BEING STUDIED",
    "3. WHAT THE MATERIALS ARE",
    "4. DERIVED",
    "5. PRESENTATION",
]


def top_level_line(folder, name):
    """Line number where a top-level statement binds `name`."""
    tree = ast.parse((ROOT / folder / "study_config.py").read_text())
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node.lineno
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return node.lineno
    raise AssertionError(f"{folder}: {name} is not defined at top level")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_study_file_has_its_five_sections_in_order(folder):
    text = (ROOT / folder / "study_config.py").read_text()
    positions = []
    for title in SECTION_TITLES:
        assert text.count(f"# {title}") == 1, f"section banner missing: {title}"
        positions.append(text.index(f"# {title}"))
    assert positions == sorted(positions), "sections are out of order"


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_derived_block_comes_after_what_it_reads(folder):
    """The ordering is load-bearing, not cosmetic: Python executes top to bottom, so a
    derived dict placed above its sources is a NameError at import."""
    derived = top_level_line(folder, "DEFAULT_PARAMS_LEVEL5")
    for source in ("CONDITIONS", "MODEL_OPTIONS", "OXIDE_DEFECTS", "METALS",
                   "OXIDES", "MICROSTRUCTURE"):
        assert top_level_line(folder, source) < derived, (
            f"{source} is defined after DEFAULT_PARAMS_LEVEL5, which reads it")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_model_switches_reference_model_options(folder):
    """The five switches had three homes. Now they have one, and the dicts that deliver
    them to the solvers carry a reference rather than a second copy."""
    tree = ast.parse((ROOT / folder / "study_config.py").read_text())
    wanted = {
        "MICROSTRUCTURE":        ["include_gb_trapping"],
        "OXIDE_DEFECTS":         ["use_sieverts_pinhole"],
        "DEFAULT_PARAMS_LEVEL5": ["include_gb_enhancement", "include_trapping",
                                  "D_eff_method"],
    }
    seen = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in wanted:
                seen[target.id] = {k.value: v for k, v in
                                   zip(node.value.keys, node.value.values)
                                   if isinstance(k, ast.Constant)}
    for dict_name, keys in wanted.items():
        assert dict_name in seen, f"{dict_name} not found"
        for key in keys:
            value = seen[dict_name][key]
            assert any(isinstance(n, ast.Name) and n.id == "MODEL_OPTIONS"
                       for n in ast.walk(value)), (
                f"{dict_name}[{key!r}] restates the switch instead of referencing "
                f"MODEL_OPTIONS")


# --------------------------------------------------------------------------
# Fix 4 — CONDITIONS owns the swept operating ranges
# --------------------------------------------------------------------------

SWEPT_FROM_CONDITIONS = {
    "temperature":     "T_range",
    "P_upstream":      "P_range",
    "metal_thickness": "L_metal_range",
    "oxide_thickness": "L_oxide_range",
}


def suggested_range_nodes(folder, dict_name="SUGGESTED_RANGES_LEVEL5"):
    tree = ast.parse((ROOT / folder / "study_config.py").read_text())
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == dict_name
                        for t in node.targets)):
            return {k.value: v for k, v in zip(node.value.keys, node.value.values)
                    if isinstance(k, ast.Constant)}
    raise AssertionError(f"{folder}: {dict_name} not found")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_swept_operating_ranges_read_conditions(folder):
    """The range of conditions under study is itself an operating condition.

    It used to be stated twice, and the two disagreed — including a 1e12 Pa ceiling
    here and a sub-atomic 1e-10 m floor there.
    """
    entries = suggested_range_nodes(folder)
    for param in SWEPT_FROM_CONDITIONS:
        assert any(isinstance(n, ast.Name) and n.id == "CONDITIONS"
                   for n in ast.walk(entries[param])), (
            f"SUGGESTED_RANGES_LEVEL5[{param!r}] is a literal; it must read CONDITIONS")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_swept_ranges_equal_their_conditions_window(folder):
    cfg = params_of(folder)
    proc = run_py(
        "import json\n"
        "from calculations.config.model_config import "
        "SUGGESTED_RANGES_LEVEL5 as S, SUGGESTED_RANGES_LEVEL5L6 as S6\n"
        "print(json.dumps({'l5': {k: list(S[k]) for k in "
        f"{list(SWEPT_FROM_CONDITIONS)!r}" "},"
        " 'l5l6': {k: list(S6[k]) for k in "
        f"{list(SWEPT_FROM_CONDITIONS)!r}" "}}))",
        cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    got = json.loads(proc.stdout)
    for param, cond_key in SWEPT_FROM_CONDITIONS.items():
        window = list(cfg["conditions"][cond_key])
        assert got["l5"][param] == window, param
        assert got["l5l6"][param] == window, param


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_material_parameter_ranges_stay_literal(folder):
    """Only the operating ranges move. D_ref and friends describe materials, and a
    material has no counterpart in CONDITIONS to be read from."""
    entries = suggested_range_nodes(folder)
    for param in ("D_ref", "E_D", "K_s_ref", "trap_gb_N_T"):
        assert not any(isinstance(n, ast.Name) and n.id == "CONDITIONS"
                       for n in ast.walk(entries[param])), (
            f"{param} is a material range and must not derive from CONDITIONS")


def test_preset_divergence_reports_the_deliberate_overrides():
    """The 'oxide' presets sweep metal_thickness below the floor on purpose — thin metal
    is how they stop the metal limiting. Reported, never failed."""
    proc = run_py(
        "import json\n"
        "from calculations.sensitivity import preset_range_divergence\n"
        "print(json.dumps(preset_range_divergence()))",
        cwd=ROOT / "Guo_etal_2025")
    assert proc.returncode == 0, proc.stderr
    notes = json.loads(proc.stdout)
    assert len(notes) == 2, notes
    assert all("metal_thickness" in n and "oxide" in n for n in notes), notes


def test_preset_divergence_catches_a_new_departure():
    """A preset that quietly widened past the studied conditions must show up."""
    proc = run_py(
        "import json\n"
        "from calculations.sensitivity import preset_range_divergence\n"
        "fake = {'made_up': {'temperature': [1, 5000]}}\n"
        "print(json.dumps(preset_range_divergence(presets=fake)))",
        cwd=ROOT / "Guo_etal_2025")
    assert proc.returncode == 0, proc.stderr
    notes = json.loads(proc.stdout)
    assert len(notes) == 1 and "temperature" in notes[0], notes
