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

REPORT = """
import json
from calculations.config.model_config import CONDITIONS
print(json.dumps({
    'conditions': {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in CONDITIONS.items()},
    'conditions_keys': sorted(CONDITIONS),
}))
"""


def config_of(folder):
    proc = run_py(REPORT, cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


# --------------------------------------------------------------------------
# Fix 7 — VALIDATION is gone; CONDITIONS carries what it held
# --------------------------------------------------------------------------

#: What VALIDATION used to state separately, now named once in CONDITIONS.
ABSORBED_ENTRIES = ("defect_fractions_compare", "defect_fraction_range",
                    "defect_limit_check", "f_gb_max", "k_diss_sweep_range",
                    "n_points")

#: The eight names that meant "how many points".
RETIRED_COUNTS = ("n_T_points", "n_P_points", "n_L_points")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_no_study_defines_validation(folder):
    """It was a fourth statement of the operating point, wearing level prefixes."""
    proc = run_py("import calculations.config.model_config as m\n"
                  "print(hasattr(m, 'VALIDATION'))", cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "False", "VALIDATION still exists"


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_conditions_absorbed_what_validation_held(folder):
    cfg = config_of(folder)
    missing = [k for k in ABSORBED_ENTRIES if k not in cfg["conditions_keys"]]
    assert not missing, f"CONDITIONS is missing {missing}"


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_one_point_count_not_eight(folder):
    """Resolution is one idea. It had eight names, two of which controlled nothing."""
    keys = config_of(folder)["conditions_keys"]
    survivors = [k for k in RETIRED_COUNTS if k in keys]
    assert not survivors, f"retired count names still present: {survivors}"
    assert "n_points" in keys


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_kinetics_sweep_derives_from_the_sensitivity_range(folder):
    """The floor carries a measured factor; the ceiling needs none.

    At the sensitivity range's ceiling the flux is already within 0.04% of the
    Sieverts limit, so the fast end needs no extension. At its floor the system is
    still metal-limited, so the slow end does — x1e-5, where the surface becomes
    limiting. x1e-20 puts brentq beyond bracketing and the notebook raises.
    """
    text = (ROOT / folder / "study_config.py").read_text()
    assert "CONDITIONS['k_diss_sweep_range']" in text
    assert "SUGGESTED_RANGES_LEVEL5L6['k_diss_ref'][0] * 1e-5" in text
    assert "SUGGESTED_RANGES_LEVEL5L6['k_diss_ref'][1]" in text


def test_validation_is_no_longer_part_of_the_study_contract():
    from calculations.config.model_config import REQUIRED_EXPORTS
    assert "VALIDATION" not in REQUIRED_EXPORTS


# [A-Za-z], not [A-Z]: the sensitivity_regime_* and regime_parallel_coords*
# notebooks start lowercase, so a capital-only glob silently covered 8 of 24.
NOTEBOOKS_ALL = sorted(ROOT.glob("*/[A-Za-z]*.ipynb"))


@pytest.mark.parametrize("nb_path", NOTEBOOKS_ALL, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_no_notebook_subscripts_validation(nb_path):
    """Grids are built where they are plotted, from a CONDITIONS range."""
    nb = json.loads(nb_path.read_text())
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        for node in ast.walk(ast.parse("".join(cell["source"]))):
            if (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
                    and node.value.id == "VALIDATION"):
                raise AssertionError(f"cell {i} still subscripts VALIDATION")


@pytest.mark.parametrize("nb_path", NOTEBOOKS_ALL, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_no_notebook_reads_a_retired_point_count(nb_path):
    code = "\n".join("".join(c["source"]) for c in json.loads(nb_path.read_text())["cells"]
                      if c["cell_type"] == "code")
    for name in RETIRED_COUNTS:
        assert not re.search(rf"\[['\"]{name}['\"]\]", code), \
            f"still reads a retired count: {name}"


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


# --------------------------------------------------------------------------
# Fix 5 — provenance is named as provenance
# --------------------------------------------------------------------------

MATERIAL_REPORT = """
import json
from calculations.config.model_config import METALS, OXIDES
print(json.dumps({
    'metal_keys': sorted(next(iter(METALS.values()))),
    'oxide_keys': sorted(next(iter(OXIDES.values()))),
}))
"""

#: Keys that read as the geometry of the wall being simulated. That is CONDITIONS' job;
#: a material dict holding one invites exactly the confusion this programme removed.
GEOMETRY_LOOKALIKES = {"thickness", "metal_thickness", "oxide_thickness",
                       "L_metal", "L_oxide"}


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_materials_hold_no_geometry_lookalike_keys(folder):
    proc = run_py(MATERIAL_REPORT, cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    keys = json.loads(proc.stdout)
    for which in ("metal_keys", "oxide_keys"):
        offenders = GEOMETRY_LOOKALIKES & set(keys[which])
        assert not offenders, (
            f"{which} includes {sorted(offenders)} — geometry belongs to CONDITIONS; "
            f"a measured specimen is 'sample_thickness'")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_both_materials_record_their_sample_thickness(folder):
    """The provenance must survive the rename: it is what ties the transport data to the
    specimen they were measured on."""
    keys = json.loads(run_py(MATERIAL_REPORT, cwd=ROOT / folder).stdout)
    assert "sample_thickness" in keys["metal_keys"]
    assert "sample_thickness" in keys["oxide_keys"]


# --------------------------------------------------------------------------
# The study contract, re-derived rather than trusted
# --------------------------------------------------------------------------

def model_config_imports():
    """Every name the codebase imports from model_config, and where from.

    Scans calculations/, container/, docs/_tools/ and ALL notebooks. The notebooks
    matter: the contract drifted precisely because the names only they imported —
    R, F, apply_style, DEFAULT_SWEEP_TEMPERATURES — were never counted.
    """
    found = {}
    sources = (list((ROOT / "calculations").rglob("*.py"))
               + list((ROOT / "container").glob("*.py"))
               + list((ROOT / "docs" / "_tools").glob("*.py")))
    blocks = [(str(p.relative_to(ROOT)), p.read_text()) for p in sources]

    for nb in sorted(ROOT.glob("*/[A-Za-z]*.ipynb")):
        for cell in json.loads(nb.read_text())["cells"]:
            if cell["cell_type"] == "code":
                blocks.append((f"{nb.parent.name}/{nb.name}", "".join(cell["source"])))

    for where, src in blocks:
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if (isinstance(node, ast.ImportFrom) and node.module
                    and "model_config" in node.module):
                for alias in node.names:
                    found.setdefault(alias.name, set()).add(where)
    return found


def study_namespace(folder):
    proc = run_py("import json, importlib.util\n"
                  f"spec = importlib.util.spec_from_file_location('s', r'{ROOT / folder}/study_config.py')\n"
                  "m = importlib.util.module_from_spec(spec)\n"
                  "spec.loader.exec_module(m)\n"
                  "print(json.dumps(sorted(n for n in dir(m) if not n.startswith('_'))))",
                  cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    return set(json.loads(proc.stdout))


def test_required_exports_covers_everything_actually_imported():
    """The contract must name every study-provided thing the codebase imports.

    A name missing here passes the study check and then fails later as a raw
    ImportError or KeyError — the exact failure the check exists to prevent. Five
    names were missing when this test was written.
    """
    from calculations.config.model_config import REQUIRED_EXPORTS

    imported = model_config_imports()
    provided = study_namespace("Guo_etal_2025")
    needed = {n for n in imported if n in provided}
    missing = sorted(needed - set(REQUIRED_EXPORTS))
    assert not missing, (
        "REQUIRED_EXPORTS is short of: "
        + ", ".join(f"{n} (imported by {sorted(imported[n])[0]})" for n in missing))


def test_required_exports_names_nothing_unused():
    """The other direction: a name nobody imports is a demand on studies for nothing."""
    from calculations.config.model_config import REQUIRED_EXPORTS

    imported = set(model_config_imports())
    unused = sorted(set(REQUIRED_EXPORTS) - imported)
    assert not unused, f"REQUIRED_EXPORTS demands unimported names: {unused}"


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_every_study_satisfies_the_widened_contract(folder):
    from calculations.config.model_config import REQUIRED_EXPORTS

    missing = sorted(set(REQUIRED_EXPORTS) - study_namespace(folder))
    assert not missing, f"{folder} does not provide {missing}"
