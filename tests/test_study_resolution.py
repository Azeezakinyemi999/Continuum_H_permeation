"""How the active study is chosen, and what happens when nothing chooses one.

A study is a folder holding its own ``study_config.py``. That file is selected by the
working directory, or named through ``$MHI_STUDY``, and by nothing else. Because the
choice is made at import time, most of these tests run a fresh interpreter (``run_py``)
rather than asserting against the one pytest is running in.
"""
import json
import re

import pytest

from conftest import KNOWN_STUDIES, ROOT, STUDIES, run_py

REPORT = (
    "from calculations.config.model_config import ("
    "    ACTIVE_STUDY_NAME, ACTIVE_STUDY_PATH, ACTIVE_STUDY_SOURCE, METALS)\n"
    "import json, pathlib\n"
    "print(json.dumps({'name': ACTIVE_STUDY_NAME,"
    "                  'path': str(pathlib.Path(ACTIVE_STUDY_PATH)),"
    "                  'source': ACTIVE_STUDY_SOURCE,"
    "                  'metal': sorted(METALS)[0]}))"
)


def resolved(cwd=None, env=None):
    """The study a fresh interpreter resolves from ``cwd`` / ``env``."""
    proc = run_py(REPORT, cwd=cwd, env=env)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


# --------------------------------------------------------------------------
# the normal case: a folder selects its own study
# --------------------------------------------------------------------------

@pytest.mark.parametrize("folder,study_name,metal", STUDIES)
def test_folder_resolves_its_own_study(folder, study_name, metal):
    """The loaded study is the file in THIS folder — its name and its own data."""
    got = resolved(cwd=ROOT / folder)
    assert got["source"] == "folder"
    assert got["path"] == str(ROOT / folder / "study_config.py")
    assert got["name"] == study_name
    assert got["metal"] == metal


@pytest.mark.parametrize("folder,study_name,metal", STUDIES)
def test_known_studies_still_hold_the_material_they_should(folder, study_name, metal):
    """Discovery agrees with the file whatever it says; this pins what it should say."""
    if folder not in KNOWN_STUDIES:
        pytest.skip(f"{folder} is a new study with no pinned expectation")
    assert (study_name, metal) == KNOWN_STUDIES[folder]


def test_study_name_is_declared_by_the_file_not_the_folder():
    """Application/ holds a copy of the 316L study and must report the 316L name.

    Deriving the label from the folder would name it 'Application', and the container
    gate — which selects its reference set on ACTIVE_STUDY_NAME — would fail to find
    one for the study it actually ships.
    """
    assert resolved(cwd=ROOT / "Application")["name"] == "Guo_etal_2025_316L"


def test_a_subdirectory_inherits_the_parent_folders_study():
    """The search walks up, so a scratch dir inside a study folder still resolves."""
    sub = ROOT / "Guo_etal_2025" / "sa_results_L5"
    if not sub.is_dir():
        pytest.skip("no output subdirectory present to walk up from")
    got = resolved(cwd=sub)
    assert got["name"] == "Guo_etal_2025_316L"
    assert got["path"] == str(ROOT / "Guo_etal_2025" / "study_config.py")


# --------------------------------------------------------------------------
# no default: naming nothing is an error, not a guess
# --------------------------------------------------------------------------

def test_repo_root_raises_and_lists_the_candidates():
    proc = run_py("import calculations.config.model_config")
    assert proc.returncode != 0
    assert "no study selected" in proc.stderr
    for folder, _, _ in STUDIES:
        assert f"{folder}/study_config.py" in proc.stderr


def test_missing_study_file_lists_what_does_exist():
    proc = run_py("import calculations.config.model_config",
                  env={"MHI_STUDY": "Nope/study_config.py"})
    assert proc.returncode != 0
    assert "does not exist" in proc.stderr
    assert "Guo_etal_2025/study_config.py" in proc.stderr


def test_incomplete_study_names_every_missing_export(tmp_path):
    """The names checked here are ones the DERIVED block in model_config.py reads.

    That is the case worth pinning. A study missing CONDITIONS or MODEL_OPTIONS does not
    merely disappoint a later import — it breaks model_config's own derivation, and
    without this check the reader gets a bare NameError raised from inside the loader
    instead of a sentence naming the study and the gap.
    """
    broken = tmp_path / "study_config.py"
    broken.write_text("STUDY_NAME = 'broken'\nMETALS = {}\nOXIDES = {}\n")
    proc = run_py("import calculations.config.model_config",
                  env={"MHI_STUDY": str(broken)})
    assert proc.returncode != 0
    assert "missing" in proc.stderr and "required export" in proc.stderr
    for name in ("CONDITIONS", "MODEL_OPTIONS", "MICROSTRUCTURE", "REGIME_PRESETS"):
        assert name in proc.stderr


# --------------------------------------------------------------------------
# the override
# --------------------------------------------------------------------------

def test_env_var_beats_the_folder():
    got = resolved(cwd=ROOT / "Guo_etal_2025",
                   env={"MHI_STUDY": "incoloy802_cr2o3/study_config.py"})
    assert got["source"] == "env"
    assert got["name"] == "incoloy802_cr2o3"


def test_env_var_relative_path_resolves_from_the_repo_root():
    """Not from the working directory — the same value must mean the same study.

    Run from inside Fuerst_et_al_2024/, a cwd-relative reading of this path would look
    for Fuerst_et_al_2024/incoloy802_cr2o3/study_config.py and fail.
    """
    got = resolved(cwd=ROOT / "Fuerst_et_al_2024",
                   env={"MHI_STUDY": "incoloy802_cr2o3/study_config.py"})
    assert got["name"] == "incoloy802_cr2o3"
    assert got["path"] == str(ROOT / "incoloy802_cr2o3" / "study_config.py")


# --------------------------------------------------------------------------
# the guard notebooks call in their first cell
# --------------------------------------------------------------------------

GUARD = ("from calculations.config import assert_study_is_local\n"
         "print(assert_study_is_local())")


@pytest.mark.parametrize("folder,_name,_metal", STUDIES)
def test_guard_passes_inside_a_study_folder(folder, _name, _metal):
    proc = run_py(GUARD, cwd=ROOT / folder)
    assert proc.returncode == 0, proc.stderr
    assert "via folder" in proc.stdout


def test_guard_rejects_a_study_inherited_from_a_parent():
    """Walk-up is fine for resolution but wrong for a notebook.

    A notebook in a subdirectory that silently used its parent's study would write
    results into a folder whose own config it never read.
    """
    sub = ROOT / "Guo_etal_2025" / "sa_results_L5"
    if not sub.is_dir():
        pytest.skip("no output subdirectory present to walk up from")
    proc = run_py(GUARD, cwd=sub)
    assert proc.returncode != 0
    assert "not from" in proc.stderr or "came from" in proc.stderr


def test_guard_rejects_the_env_var_route():
    proc = run_py(GUARD, cwd=ROOT,
                  env={"MHI_STUDY": "Guo_etal_2025/study_config.py"})
    assert proc.returncode != 0
    assert "resolved via 'env'" in proc.stderr


# --------------------------------------------------------------------------
# reading a second study without making it active
# --------------------------------------------------------------------------

def test_load_study_does_not_change_the_active_one():
    proc = run_py(
        "from calculations.config.model_config import load_study, ACTIVE_STUDY_NAME, METALS\n"
        "other = load_study('Fuerst_et_al_2024/study_config.py')\n"
        "print(ACTIVE_STUDY_NAME, other.STUDY_NAME, sorted(METALS)[0], sorted(other.METALS)[0])",
        cwd=ROOT / "Guo_etal_2025")
    assert proc.returncode == 0, proc.stderr
    active, other, active_metal, other_metal = proc.stdout.split()
    assert (active, active_metal) == ("Guo_etal_2025_316L", "316L_Guo_2025")
    assert (other, other_metal) == ("fuerst_etal_2024_model_config",
                                    "Hastelloy_N_fuerst_2024")


def test_load_study_enforces_the_same_contract(tmp_path):
    broken = tmp_path / "study_config.py"
    broken.write_text("METALS = {}\n")
    proc = run_py(
        f"from calculations.config.model_config import load_study\n"
        f"load_study({str(broken)!r})",
        cwd=ROOT / "Guo_etal_2025")
    assert proc.returncode != 0
    assert "required export" in proc.stderr


# --------------------------------------------------------------------------
# the contract, and the things that key off it
# --------------------------------------------------------------------------

@pytest.mark.parametrize("folder,_name,_metal", STUDIES)
def test_every_study_satisfies_required_exports(folder, _name, _metal):
    from calculations.config.model_config import REQUIRED_EXPORTS, load_study

    study = load_study(f"{folder}/study_config.py")
    missing = REQUIRED_EXPORTS - set(dir(study))
    assert not missing, f"{folder} is missing {sorted(missing)}"


def test_the_release_gate_has_a_reference_set_for_every_study():
    """Otherwise the gate exits 1 on a study it is supposed to certify.

    Parsed rather than imported: importing the smoke test runs it.
    """
    source = (ROOT / "container" / "model_smoke_test.py").read_text()
    block = source[source.index("REFERENCE = {"):source.index("failures:")]
    keys = set(re.findall(r"^    '([^']+)': \{", block, re.M))
    assert keys, "could not parse any REFERENCE keys — has the file been restructured?"
    missing = {study_name: folder for folder, study_name, _ in STUDIES
               if study_name not in keys}
    assert not missing, (
        f"these studies have no REFERENCE entry in container/model_smoke_test.py, "
        f"so the release gate would exit 1 on them: {missing}")


def test_the_deleted_studies_package_stays_deleted():
    """Studies live in their folders. A module here would be a second home."""
    assert not (ROOT / "calculations" / "config" / "studies").exists()
    proc = run_py("import calculations.config.studies",
                  env={"MHI_STUDY": "Guo_etal_2025/study_config.py"})
    assert proc.returncode != 0
    assert "No module named" in proc.stderr


# --------------------------------------------------------------------------
# The container README documents a shipped image, so it is allowed to describe
# the switch as history — but only as history.
# --------------------------------------------------------------------------

CONTAINER_README = ROOT / "container" / "README.md"
SUPERSEDED_HEADING = "### Since 1.1.0"


def test_container_readme_mentions_the_switch_only_as_history():
    """`hydrogen-model:1.1.0` really does have an ACTIVE_STUDY switch.

    Its changelog must stay true to the image it describes, so the name cannot simply
    be banned here. What must not happen is a *live* instruction telling a reader to
    edit a switch that no longer exists — so every mention has to sit below the
    heading that marks the repository as having moved on.
    """
    text = CONTAINER_README.read_text()
    cut = text.find(SUPERSEDED_HEADING)
    assert cut != -1, f"{CONTAINER_README.name} has no {SUPERSEDED_HEADING!r} section"

    stale = [
        line.strip()
        for line in text[:cut].splitlines()
        if "ACTIVE_STUDY" in line or "config/studies" in line
    ]
    assert not stale, (
        "the deleted switch is described above the 'Since 1.1.0' heading, where it "
        f"reads as current instruction: {stale}"
    )
