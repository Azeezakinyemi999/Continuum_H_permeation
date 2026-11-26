"""Where a notebook's figures are written, and that they are written at all.

`calculations.figures` is small, but it sits between every notebook and its output:
if `save_figure` raised on a skipped plot, or quietly wrote to the working directory
instead of the results folder, a whole notebook run would produce nothing where it was
expected. The last test guards the convention itself against drifting back.
"""
import ast
import json
import re
import subprocess

import matplotlib
matplotlib.use("Agg")                     # no window server in a test run
import matplotlib.pyplot as plt           # noqa: E402
import pytest                             # noqa: E402

from calculations.figures import RESULTS_SUFFIX, Panel, results_dir, save_figure, save_panels  # noqa: E402
from conftest import ROOT                 # noqa: E402

def _not_ignored(paths):
    """Drop paths .gitignore excludes.

    The glob walks the filesystem, so a scratch folder beside the studies -- a
    one-off calculation, a copy kept for a paper -- was collected and failed the
    convention tests, which exist to police the STUDY notebooks. Asking git rather
    than hardcoding a folder name means a new scratch directory is handled by adding
    one line to .gitignore, where it belongs.

    `git check-ignore` exits 1 when nothing matches, which is the normal case, so the
    return code is not an error here. If git is unavailable the glob is used as-is:
    over-collecting is the safer failure, since it only ever adds coverage.
    """
    if not paths:
        return paths
    try:
        done = subprocess.run(
            ["git", "check-ignore", "--stdin"],
            input="\n".join(str(p) for p in paths),
            capture_output=True, text=True, cwd=ROOT, timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return paths
    ignored = {line.strip() for line in done.stdout.splitlines() if line.strip()}
    return [p for p in paths if str(p) not in ignored]


NOTEBOOKS = _not_ignored(sorted(ROOT.glob("*/[A-Za-z]*.ipynb")))


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
# save_panels: one figure per panel, in the house style
# --------------------------------------------------------------------------

def test_save_panels_writes_each_panel_and_clears_what_an_earlier_run_left(tmp_path):
    """Each panel becomes <name>-<panel>.pdf, .png (180 mm wide at 600 dpi) and _caption.tex; the
    old one-file figure and panels no longer drawn go; other figures stay."""
    from PIL import Image
    from calculations import house_style as hs
    for name in ("L3_arrhenius.png", "L3_arrhenius.pdf", "L3_arrhenius-gone.png", "L3_arrhenius-gone_caption.tex",
                 "L3_arrhenius_3P.png", "L5_arrhenius-flux.png"):
        (tmp_path / name).write_text("")
    panels = []
    for name in ("flux", "1e-7Pa-phi"):
        fig, ax = hs.new_figure()
        ax.plot([1, 2], [1, 3])
        panels.append(Panel(name, fig, ax, "The caption, 600~K.", [1, 2], [1, 3]))
    printed = []
    written = save_panels(panels, "L3_arrhenius", tmp_path, module="tests", log=printed.append)
    assert [p.name for p in written] == ["L3_arrhenius-flux.png", "L3_arrhenius-1e-7Pa-phi.png"]
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "L3_arrhenius-1e-7Pa-phi.pdf", "L3_arrhenius-1e-7Pa-phi.png", "L3_arrhenius-1e-7Pa-phi_caption.tex",
        "L3_arrhenius-flux.pdf", "L3_arrhenius-flux.png", "L3_arrhenius-flux_caption.tex",
        "L3_arrhenius_3P.png", "L5_arrhenius-flux.png"]
    caption = (tmp_path / "L3_arrhenius-flux_caption.tex").read_text()
    assert r"\includegraphics[width=180mm]{L3_arrhenius-flux.pdf}" in caption and "600~K" in caption
    with Image.open(tmp_path / "L3_arrhenius-flux.png") as im:
        assert abs(im.size[0] - round(180 / 25.4 * 600)) <= 2
    assert not [p for p in printed if "outside" in p]
    assert not plt.get_fignums()                            # every panel is closed


def test_save_panels_tolerates_a_plot_that_drew_nothing(tmp_path):
    assert save_panels([], "skipped", tmp_path) == []


def test_house_colours_keep_one_meaning():
    """A +L6 level keeps its level's hue; the limiting-step colours of wall_labels are the house ones."""
    from calculations import house_style as hs
    from calculations.wall_labels import STEP_STYLE
    assert hs.level_style("L3L6")["color"] == hs.level_style("L3")["color"] and hs.level_style("L3L6")["ls"] == "--"
    for step, colour in hs.STEP_COLOURS.items():
        assert STEP_STYLE[step][1].lower() == colour.lower(), step


def _two_curves():
    import numpy as np
    from matplotlib.lines import Line2D
    from calculations import house_style as hs
    hs.use_style()
    fig, ax = hs.new_figure()
    x = np.linspace(0, 1, 40)
    ax.plot(x, 0.2 + 0.1 * x, color=hs.STEP_COLOURS["metal diffusion"], lw=3)
    ax.plot(x, 0.8 - 0.1 * x, color=hs.STEP_COLOURS["oxide diffusion"], lw=3, ls="--")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    handles = [Line2D([], [], color=hs.STEP_COLOURS["metal diffusion"], lw=3),
               Line2D([], [], color=hs.STEP_COLOURS["oxide diffusion"], lw=3, ls="--")]
    return hs, fig, ax, handles


def test_labels_go_inside_the_plot_next_to_their_data():
    """Each entry is written next to its own curve, in its colour, with no key below the plot."""
    hs, fig, ax, handles = _two_curves()
    kind, artists = hs.labels_or_key(fig, ax, handles, ["Metal diffusion", "Oxide diffusion"],
                                     colours=[h.get_color() for h in handles])
    assert kind == "inside" and not fig.legends
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    low, high = (hs.bbox_frac(ax, t, renderer) for t in artists[:2])
    assert low.y1 < 0.55 < high.y0, "each label sits on its own curve's side of the plot"
    plt.close(fig)


def test_labels_fall_back_to_a_key_when_the_plot_is_full():
    """With no free space (a dense cloud of points over the whole plot) the key goes below the plot."""
    import numpy as np
    hs, fig, ax, handles = _two_curves()
    g = np.linspace(0.01, 0.99, 120)
    gx, gy = np.meshgrid(g, g)
    ax.plot(gx.ravel(), gy.ravel(), ls="none", marker=".", ms=2, color="0.5")
    kind, legend = hs.labels_or_key(fig, ax, handles, ["Metal diffusion", "Oxide diffusion"])
    assert kind == "key" and fig.legends and not [t for t in ax.texts if t.get_text()]
    plt.close(fig)


def test_report_groups_a_figures_panels_in_drawing_order(tmp_path):
    """The study report shows each panel of a figure, in the order save_panels wrote them."""
    import os
    import time
    from types import SimpleNamespace
    from calculations import report, report_latex
    folder = tmp_path / "effect-of-temperature-results"
    folder.mkdir()
    for k, name in enumerate(("L3_arrhenius-phi", "L3_arrhenius-flux", "L3_arrhenius-1e-7Pa-flux")):
        for ext in ("pdf", "png"):
            path = folder / f"{name}.{ext}"
            path.write_text("")
            t = time.time() - 100 + k
            os.utime(path, (t, t))
    (folder / "L1_arrhenius.png").write_text("")
    rep = SimpleNamespace(out=str(tmp_path))
    rep.files = lambda f, pattern="*": report.StudyReport.files(rep, f, pattern)
    figs = report.StudyReport.figures(rep, folder.name)
    assert figs["L1_arrhenius"] == "effect-of-temperature-results/L1_arrhenius.png"
    assert [os.path.basename(f) for f in figs["L3_arrhenius"]] == [
        "L3_arrhenius-phi.pdf", "L3_arrhenius-flux.pdf", "L3_arrhenius-1e-7Pa-flux.pdf"]
    assert report_latex.panel_note("x/L3_arrhenius-1e-7Pa-flux.pdf") == "Panel: 1e-7 Pa, flux."
    assert report_latex.panel_note("x/L1_arrhenius.png") == ""


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
def test_no_notebook_sets_its_own_style(nb_path):
    """Every figure takes the house style from calculations/house_style.py: no notebook brings back
    its own rcParams, palettes or the old PLOT_STYLE / COLORS / CURVE_STYLES settings."""
    code = _code(nb_path)
    old = re.findall(r"\b(PLOT_STYLE|CURVE_STYLES|apply_style|rcParams\.update|INK, GREY)\b|\bCOLORS\[", code)
    assert not old, f"old style settings still in the notebook: {sorted(set(old))}"


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
    """(cell index, line number, name expression, names it interpolates, scope) of every
    save_figure(...) and save_panels(...) call whose name is a string or an f-string."""
    nb = json.loads(nb_path.read_text())
    out = []
    for ci, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"])
        if "save_figure(" not in src and "save_panels(" not in src:
            continue
        tree = ast.parse(src)
        assigned = {t.id for n in ast.walk(tree) if isinstance(n, ast.Assign)
                    for t in n.targets if isinstance(t, ast.Name)}
        for call in ast.walk(tree):
            if not (isinstance(call, ast.Call) and getattr(call.func, "id", None) in ("save_figure", "save_panels")
                    and len(call.args) >= 2):
                continue
            arg = call.args[1]
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                name = arg.value
            elif isinstance(arg, ast.JoinedStr):
                name = ast.get_source_segment(src, arg)[2:-1]
            else:
                continue
            ln = call.lineno
            loops = {n.target.id for n in ast.walk(tree)
                     if isinstance(n, ast.For) and isinstance(n.target, ast.Name)
                     and n.lineno <= ln <= (n.end_lineno or n.lineno)}
            # a save inside a function may name that function's parameters: they are bound when it runs
            params = {a.arg for n in ast.walk(tree)
                      if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
                      and n.lineno <= ln <= (n.end_lineno or n.lineno)
                      for a in (*n.args.posonlyargs, *n.args.args, *n.args.kwonlyargs)}
            out.append((ci, ln, name, set(re.findall(r"\{(\w+)", name)), loops | assigned | params))
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
            f"which is not bound by an enclosing loop, an assignment in that cell, or a "
            f"parameter of the enclosing function")


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_figure_names_are_unique_within_a_notebook(nb_path):
    """Two saves with one name means the second silently overwrites the first."""
    fixed = [name for *_, name, needed, _ in
             ((c, l, n, nd, s) for c, l, n, nd, s in _save_sites(nb_path)) if not needed]
    duplicates = {n for n in fixed if fixed.count(n) > 1}
    assert not duplicates, f"figure names used more than once: {sorted(duplicates)}"
