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
import textwrap

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


MODEL_CONFIG = ROOT / "calculations" / "config" / "model_config.py"


def derived_dict(name):
    """A top-level dict literal in model_config.py, as {key: value AST node}.

    These checks used to run against each study_config.py, parametrised over the four
    folders. The derived block now lives once in model_config.py, so they run once — and
    the property holds for every study by construction rather than by four copies each
    being policed separately.
    """
    tree = ast.parse(MODEL_CONFIG.read_text())
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == name for t in node.targets)):
            return {k.value: v for k, v in zip(node.value.keys, node.value.values)
                    if isinstance(k, ast.Constant)}
    raise AssertionError(f"{name} is not defined at top level in model_config.py")


def test_oxide_thickness_is_written_as_a_conditions_reference():
    """Value equality would pass a literal 1e-6; the derivation has to be in the source."""
    value = derived_dict("DEFAULT_PARAMS_LEVEL5")["oxide_thickness"]
    assert any(isinstance(n, ast.Name) and n.id == "CONDITIONS"
               for n in ast.walk(value)), \
        "oxide_thickness is a literal; it must read CONDITIONS"


# --------------------------------------------------------------------------
# Fix 3b — the study file is sectioned, in dependency order
# --------------------------------------------------------------------------

SECTION_TITLES = [
    "1. IDENTITY AND CONSTANTS",
    "2. WHAT IS BEING STUDIED",
    "3. WHAT THE MATERIALS ARE",
]

# Sections 4 (DERIVED) and 5 (PRESENTATION) were byte-identical in all four studies, so
# they moved into model_config.py. A study file that grows one back has started keeping a
# private copy of the derivation, which is how the four copies drifted the first time.
MOVED_SECTIONS = ["4. DERIVED", "5. PRESENTATION"]

MOVED_NAMES = [
    "DEFAULT_PARAMS_LEVEL5", "DEFAULT_PARAMS_LEVEL5L6",
    "get_surface_kinetics_from_config", "build_simulation_config",
    "PLOT_STYLE", "COLORS", "CURVE_STYLES", "apply_style",
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
def test_study_file_has_its_three_sections_in_order(folder):
    text = (ROOT / folder / "study_config.py").read_text()
    positions = []
    for title in SECTION_TITLES:
        assert text.count(f"# {title}") == 1, f"section banner missing: {title}"
        positions.append(text.index(f"# {title}"))
    assert positions == sorted(positions), "sections are out of order"


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_the_derived_and_presentation_sections_stay_out_of_the_study_file(folder):
    """A study is its own data. The shared derivation is not part of it."""
    text = (ROOT / folder / "study_config.py").read_text()
    for title in MOVED_SECTIONS:
        assert f"# {title}" not in text, (
            f"{folder} has grown a '{title}' section back; it belongs in model_config.py")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_no_study_redefines_a_moved_name(folder):
    """The duplication must not creep back one study at a time.

    A study defining its own DEFAULT_PARAMS_LEVEL5 would shadow the shared one through
    globals().update() and silently take effect for that study alone — the failure this
    move exists to make impossible.
    """
    tree = ast.parse((ROOT / folder / "study_config.py").read_text())
    defined = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            defined.add(node.name)
        elif isinstance(node, ast.Assign):
            defined |= {t.id for t in node.targets if isinstance(t, ast.Name)}
    clashes = sorted(defined & set(MOVED_NAMES))
    assert not clashes, f"{folder} redefines names that moved to model_config: {clashes}"


def test_derived_block_comes_after_the_study_namespace_is_spliced_in():
    """The ordering is load-bearing, not cosmetic.

    The derived block reads CONDITIONS, METALS and the rest as module globals, and they
    exist only because of the globals().update(_study_ns) above it. Moved above that
    line, every one of them is a NameError at import — and because model_config is
    imported by everything, nothing in the repository would run.
    """
    text = MODEL_CONFIG.read_text()
    splice = text.index("globals().update(_study_ns)")
    tree = ast.parse(text)
    lines = text.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))

    for name in ("ACTIVE_METAL", "DEFAULT_PARAMS_LEVEL5", "DEFAULT_PARAMS_LEVEL5L6"):
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == name for t in node.targets):
                assert offsets[node.lineno - 1] > splice, (
                    f"{name} is derived above globals().update(_study_ns), "
                    f"so the names it reads do not exist yet")
                break
        else:
            raise AssertionError(f"{name} is not defined at top level in model_config.py")


def test_model_switches_reference_model_options():
    """The five switches had three homes. Now they have one, and the dicts that deliver
    them to the solvers carry a reference rather than a second copy.

    MICROSTRUCTURE and OXIDE_DEFECTS are study data and keep their references in each
    study file; DEFAULT_PARAMS_LEVEL5 is derived and is checked once, in model_config.
    """
    for key in ("enable_gb_enhancement", "enable_trapping", "D_eff_method"):
        value = derived_dict("DEFAULT_PARAMS_LEVEL5")[key]
        assert any(isinstance(n, ast.Name) and n.id == "MODEL_OPTIONS"
                   for n in ast.walk(value)), (
            f"DEFAULT_PARAMS_LEVEL5[{key!r}] restates the switch instead of "
            f"referencing MODEL_OPTIONS")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_study_dicts_reference_model_options(folder):
    """The two switch-carrying dicts that are still study data."""
    tree = ast.parse((ROOT / folder / "study_config.py").read_text())
    wanted = {
        "MICROSTRUCTURE": ["enable_gb_as_trap"],
        "OXIDE_DEFECTS":  ["use_sieverts_pinhole"],
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


def derived_block_reads():
    """Names the derived block in model_config.py takes from the study namespace.

    Everything it uses that it does not define itself, minus builtins and its own
    imports. These are required of a study just as firmly as an imported name is — more
    so, since a study missing one breaks model_config's own import — but no import
    statement mentions them, so the check below has to know about them separately.
    """
    import builtins

    text = MODEL_CONFIG.read_text()
    # Bounded at both ends. Running to EOF would sweep in the provenance section below,
    # whose ACTIVE_STUDY_SOURCE and STUDY_FILENAME are model_config's own names and have
    # no business being demanded of a study.
    start = text.index("# 4. DERIVED")
    end = text.index("# Provenance.")
    assert start < end, "the derived block no longer precedes the provenance section"
    tree = ast.parse(textwrap.dedent(text[start:end]))

    bound = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            bound.add(node.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            bound.add(node.id)
        elif isinstance(node, ast.arg):
            bound.add(node.arg)
    used = {n.id for n in ast.walk(tree)
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    return {n for n in used - bound - set(dir(builtins)) if not n.startswith("_")} - {"np"}


def test_required_exports_names_nothing_unused():
    """The other direction: a name nothing needs is a demand on studies for nothing.

    'Needs' has two senses since the derived block moved here. Most names are required
    because something imports them; CONDITIONS, METALS, MODEL_OPTIONS and the rest are
    required because model_config's own derivation reads them out of the study. Only
    checking imports would call MODEL_OPTIONS unused and delete it from the contract.
    """
    from calculations.config.model_config import REQUIRED_EXPORTS

    justified = set(model_config_imports()) | derived_block_reads()
    unused = sorted(set(REQUIRED_EXPORTS) - justified)
    assert not unused, f"REQUIRED_EXPORTS demands names nothing needs: {unused}"


def test_the_derived_block_can_only_read_names_the_contract_guarantees():
    """The inverse, and the one that actually bites.

    A new line in the derived block reading, say, SUGGESTED_RANGES_LEVEL5 without adding
    it to REQUIRED_EXPORTS would work on every study that happens to define it and fail
    with a bare NameError on one that does not — raised from inside the loader, which is
    the failure REQUIRED_EXPORTS exists to replace with a sentence.
    """
    from calculations.config.model_config import REQUIRED_EXPORTS

    ungurantee = sorted(derived_block_reads() - set(REQUIRED_EXPORTS))
    assert not ungurantee, (
        "the derived block reads study names the contract does not require: "
        f"{ungurantee}")


@pytest.mark.parametrize("folder", STUDY_IDS)
def test_every_study_satisfies_the_widened_contract(folder):
    from calculations.config.model_config import REQUIRED_EXPORTS

    missing = sorted(set(REQUIRED_EXPORTS) - study_namespace(folder))
    assert not missing, f"{folder} does not provide {missing}"


# --------------------------------------------------------------------------
# Measurement-validity windows: named alike, and both actually checked
# --------------------------------------------------------------------------

@pytest.mark.parametrize("folder", STUDY_IDS)
def test_both_materials_name_their_validity_window_the_same_way(folder):
    """The metal carried 'temp_range', the oxide 'temperature_range'.

    Two names for one idea, and only the oxide's was ever read — so extrapolating the
    oxide printed a warning while extrapolating the metal said nothing.
    """
    proc = run_py(MATERIAL_REPORT, cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    keys = json.loads(proc.stdout)
    for which in ("metal_keys", "oxide_keys"):
        assert "temperature_range" in keys[which], f"{which} has no validity window"
        assert "temp_range" not in keys[which], f"{which} still uses the old name"


def test_metal_extrapolation_is_no_longer_silent(capsys):
    """Guo's metal data stop at 823 K and the study runs at 873 K."""
    from calculations.config.model_config import build_simulation_config
    from calculations.oxide_permeation import get_metal_properties_at_T

    name = build_simulation_config()["metal_name"]
    get_metal_properties_at_T(name, 873.0)
    assert "outside validated range" in capsys.readouterr().out


def test_no_warning_inside_the_metal_window(capsys):
    """A warning on every call would be noise, and noise gets filtered."""
    from calculations.config.model_config import build_simulation_config
    from calculations.oxide_permeation import get_metal_properties_at_T

    name = build_simulation_config()["metal_name"]
    get_metal_properties_at_T(name, 700.0)          # inside [468, 823]
    assert "outside validated range" not in capsys.readouterr().out


def test_oxide_extrapolation_still_warns(capsys):
    """Regression: the oxide check moved into a shared helper."""
    from calculations.config.model_config import build_simulation_config
    from calculations.oxide_permeation import get_oxide_properties_at_T

    name = build_simulation_config()["oxide_name"]
    get_oxide_properties_at_T(name, 873.0)          # outside [473, 773]
    assert "outside validated range" in capsys.readouterr().out


# --------------------------------------------------------------------------
# check_against_config: a preset must pin a parameter that exists
#
# The model wrappers build their parameters as DEFAULT_PARAMS.copy() then
# .update(params_dict), so a preset key that is misspelled, or left behind by a
# rename, is stored and never read -- no error, and a flux identical to the
# unpinned case. These two tests are the only thing standing between that and a
# regime scan that has quietly stopped being targeted.
# --------------------------------------------------------------------------

def test_the_shipped_config_has_no_preset_pinning_a_missing_parameter():
    """Scoped to THIS check, not to check_against_config as a whole.

    Guo_etal_2025 and Application report four problems from a different check --
    D_ref, k_diss_metal_ref and E_diss_metal each sit outside their own sweep range,
    so the SA never samples them at their nominal value. Those predate this guard
    (identical on main) and are a range decision, not a naming one, so they are not
    asserted away here.
    """
    from calculations.sensitivity import check_against_config

    problems = check_against_config(verbose=False, raise_on_fail=False)
    stale = [p for p in problems if "silently ignored" in p]
    assert stale == [], f"a preset pins a parameter no wrapper will read: {stale}"


def test_a_preset_pinning_an_unknown_parameter_is_caught(monkeypatch):
    """The failure this guard exists for: a name no wrapper will ever read.

    monkeypatch.setitem restores the preset afterwards, so the module-level dict
    other tests share is left as it was found.
    """
    import calculations.sensitivity as sa

    regime = next(iter(sa.REGIME_PRESETS))
    monkeypatch.setitem(sa.REGIME_PRESETS[regime], "f_pinholeee", (1e-6, 1e-5))

    problems = sa.check_against_config(verbose=False, raise_on_fail=False)
    assert any("f_pinholeee" in p and "silently ignored" in p for p in problems), (
        f"a preset pinning a nonexistent parameter was not reported: {problems}"
    )

    with pytest.raises(ValueError, match="silently ignored"):
        sa.check_against_config(verbose=False)
