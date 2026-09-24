"""Where a notebook's figures are written, and that they are written at all.

`calculations.figures` is small, but it sits between every notebook and its output:
if `save_figure` raised on a skipped plot, or quietly wrote to the working directory
instead of the results folder, a whole notebook run would produce nothing where it was
expected. The last test guards the convention itself against drifting back.
"""
import ast
import json
import re

import matplotlib
matplotlib.use("Agg")                     # no window server in a test run
import matplotlib.pyplot as plt           # noqa: E402
import pytest                             # noqa: E402

from calculations.figures import RESULTS_SUFFIX, results_dir, save_figure  # noqa: E402
from conftest import ROOT                 # noqa: E402

NOTEBOOKS = sorted(ROOT.glob("*/[A-Za-z]*.ipynb"))


# --------------------------------------------------------------------------
# results_dir
# --------------------------------------------------------------------------

def test_results_dir_names_and_creates_the_folder(tmp_path):
    folder = results_dir("Proposal", base=tmp_path)
    assert folder == tmp_path / f"Proposal{RESULTS_SUFFIX}"
    assert folder.is_dir()


def test_results_dir_accepts_the_notebook_filename(tmp_path):
    """`results_dir("Proposal.ipynb")` and `results_dir("Proposal")` are the same folder."""
    assert (results_dir("Proposal.ipynb", base=tmp_path)
            == results_dir("Proposal", base=tmp_path))


def test_results_dir_is_idempotent_and_keeps_its_contents(tmp_path):
    first = results_dir("Proposal", base=tmp_path)
    (first / "already_here.png").write_bytes(b"")
    again = results_dir("Proposal", base=tmp_path)
    assert again == first
    assert (again / "already_here.png").exists(), "must not clear the folder"


def test_results_dir_defaults_to_the_working_directory(tmp_path, monkeypatch):
    """A notebook passes no base; the folder must land beside it, not at the repo root."""
    monkeypatch.chdir(tmp_path)
    assert results_dir("Surface_proposal") == tmp_path / f"Surface_proposal{RESULTS_SUFFIX}"


def test_results_dir_rejects_a_nameless_notebook(tmp_path):
    with pytest.raises(ValueError):
        results_dir("", base=tmp_path)


# --------------------------------------------------------------------------
# save_figure
# --------------------------------------------------------------------------

def test_save_figure_writes_a_png_and_returns_its_path(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    path = save_figure(fig, "L1L6_core_validation", tmp_path)
    assert path == tmp_path / "L1L6_core_validation.png"
    assert path.stat().st_size > 0
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", "not a real PNG"


def test_save_figure_closes_the_figure_by_default(tmp_path):
    """Nothing else closes them, and a notebook run saves dozens."""
    fig, _ = plt.subplots()
    save_figure(fig, "closed", tmp_path)
    assert fig.number not in plt.get_fignums()


def test_save_figure_can_keep_the_figure_open(tmp_path):
    """Proposal.ipynb saves BEFORE plt.show(); closing there would show nothing."""
    fig, _ = plt.subplots()
    save_figure(fig, "kept", tmp_path, close=False)
    assert fig.number in plt.get_fignums()
    plt.close(fig)


def test_save_figure_tolerates_no_figure(tmp_path):
    """plot_givendata_results returns None for a skipped regime; that is not an error."""
    assert save_figure(None, "skipped", tmp_path) is None
    assert list(tmp_path.iterdir()) == []


def test_save_figure_does_not_double_the_extension(tmp_path):
    fig, _ = plt.subplots()
    assert save_figure(fig, "already.png", tmp_path).name == "already.png"


def test_save_figure_creates_a_missing_folder(tmp_path):
    fig, _ = plt.subplots()
    path = save_figure(fig, "fresh", tmp_path / "Proposal-results")
    assert path.exists()


# --------------------------------------------------------------------------
# the convention, across every notebook
# --------------------------------------------------------------------------

def _code(nb_path):
    nb = json.loads(nb_path.read_text())
    return "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_no_notebook_writes_a_figure_to_a_bare_filename(nb_path):
    """Every figure goes through save_figure into a results folder.

    A `plt.savefig('Level1.png')` writes to the working directory, which is how the
    six loose PNGs ended up beside the notebooks in the first place.
    """
    stray = re.findall(r"""plt\.savefig\(\s*['"][^'"]+['"]""", _code(nb_path))
    assert not stray, f"bare savefig still present: {stray}"


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_every_notebook_declares_its_results_folder(nb_path):
    code = _code(nb_path)
    assert "from calculations.figures import" in code
    match = re.search(r'results_dir\(\s*"([^"]+)"\s*\)', code)
    assert match, "no results_dir(...) call"
    assert match.group(1) == nb_path.stem, (
        f"{nb_path.name} writes into {match.group(1)}{RESULTS_SUFFIX}/ — a notebook's "
        f"folder must be named after the notebook")


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_every_notebook_cell_still_parses(nb_path):
    """The figure saves were inserted by script; none of them may break a cell."""
    nb = json.loads(nb_path.read_text())
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] == "code":
            ast.parse("".join(cell["source"]))       # raises SyntaxError on failure


# --------------------------------------------------------------------------
# the figure NAMES, which are only exercised when a notebook actually runs
# --------------------------------------------------------------------------

def _save_sites(nb_path):
    """(cell index, line number, name expression, names it interpolates, scope)."""
    nb = json.loads(nb_path.read_text())
    out = []
    for ci, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"])
        if "save_figure(" not in src:
            continue
        tree = ast.parse(src)
        assigned = {t.id for n in ast.walk(tree) if isinstance(n, ast.Assign)
                    for t in n.targets if isinstance(t, ast.Name)}
        for ln, line in enumerate(src.splitlines(), start=1):
            m = re.search(r'save_figure\([^,]+,\s*f?"([^"]+)"', line)
            if not m:
                continue
            loops = {n.target.id for n in ast.walk(tree)
                     if isinstance(n, ast.For) and isinstance(n.target, ast.Name)
                     and n.lineno <= ln <= (n.end_lineno or n.lineno)}
            out.append((ci, ln, m.group(1),
                        set(re.findall(r"\{(\w+)", m.group(1))), loops | assigned))
    return out


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_figure_names_only_use_variables_in_scope(nb_path):
    """An f-string naming a loop variable the cell does not have is a NameError at run time.

    The notebooks are several generations deep: the pre-isothermal copies have no
    `for T in TEMPERATURES` loop, so a name templated on T raises there and nowhere else
    — and only when that cell is reached, minutes into a run.
    """
    for ci, ln, name, needed, scope in _save_sites(nb_path):
        missing = needed - scope
        assert not missing, (
            f"cell {ci} line {ln}: figure name {name!r} interpolates {sorted(missing)}, "
            f"which is not bound by an enclosing loop or an assignment in that cell")


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_figure_names_are_unique_within_a_notebook(nb_path):
    """Two saves with one name means the second silently overwrites the first."""
    fixed = [name for *_, name, needed, _ in
             ((c, l, n, nd, s) for c, l, n, nd, s in _save_sites(nb_path)) if not needed]
    duplicates = {n for n in fixed if fixed.count(n) > 1}
    assert not duplicates, f"figure names used more than once: {sorted(duplicates)}"
