"""Where a notebook's figures go.

Each notebook writes into one folder beside it, named after it::

    Guo_etal_2025/
    ├── Proposal.ipynb
    └── Proposal-results/        Level1_Perfect_Metal.png, ...

The folder is resolved against the WORKING DIRECTORY, which for a correctly launched
notebook is its own study folder — the same rule ``RESULTS_DIR = "sa_results"`` and the
``study_config.py`` lookup already follow. Notebooks call
``calculations.config.assert_study_is_local()`` in their first cell, so a kernel started
in the wrong place fails there, before anything is written to the wrong study.

Usage, in the notebook's first cell::

    from calculations.figures import results_dir, save_figure
    FIGS = results_dir("Proposal")

and then, wherever a figure's panels are drawn (each a ``Panel`` from a house-style
plot function, calculations/house_style.py)::

    save_panels(level1_panels(...), "Level1_Perfect_Metal", FIGS)

Each panel is written as ``<name>-<panel>.pdf``, a 600 dpi ``.png`` and
``<name>-<panel>_caption.tex``, a LaTeX caption with the numbers filled in. What an earlier
run wrote for the same figure (the old one-file figure, panels no longer drawn) is removed.
A single finished figure goes through ``save_figure``.
"""

import os
import pathlib
import re

import matplotlib.pyplot as plt

#: Appended to a notebook's name to make its results folder.
RESULTS_SUFFIX = "-results"

#: Raster resolution of the PNG: 600 dpi at the figure's final printed size (180 mm wide), as the
#: house style asks. The figures are drawn at that size, so they are never scaled later.
DPI = 600

#: Written beside the PNG, sharp at any zoom. PDF rather than SVG because these figures
#: carry scatter and parallel-coordinate plots with thousands of elements, which SVG
#: renders slowly and stores verbosely. Set to None to write the raster alone.
VECTOR_FORMAT = "pdf"


#: Inline display density. Jupyter embeds a figure at matplotlib's default 100 dpi,
#: roughly a quarter of the pixel density of a current display -- which is why figures
#: look soft IN THE NOTEBOOK while the saved PNG, written at DPI above, is sharp. The
#: two are separate settings and only this one affects what you see while working.
INLINE_FORMAT = "retina"


def _sharpen_inline_display(fmt=INLINE_FORMAT):
    """Ask Jupyter for a high-density inline render. No-op outside IPython.

    Called once on import, which is deliberate. Every notebook already imports this
    module -- the convention test in tests/test_figures.py requires it -- so doing it
    here needs no per-notebook edit and cannot drift between the four studies.

    Not done with a ``%config`` line in the notebooks, which would be the obvious
    way: the test suite runs ``ast.parse`` over every notebook cell, and a magic is
    not valid Python, so it would fail two tests. No notebook currently contains a
    magic and this keeps that true.

    Failure is swallowed on purpose. This is a display nicety; a notebook must not
    fail to start because a version of matplotlib_inline renamed something.
    """
    if not fmt:
        return
    try:
        from IPython import get_ipython
    except ImportError:
        return
    ip = get_ipython()
    if ip is None:                      # plain interpreter, a script, or pytest
        return
    try:
        ip.run_line_magic("config", f"InlineBackend.figure_format = {fmt!r}")
    except Exception:                   # noqa: BLE001 - never block a notebook's first cell
        pass


_sharpen_inline_display()


def results_dir(notebook, base=None):
    """The results folder for ``notebook``, created if it does not exist.

    Parameters
    ----------
    notebook : str
        The notebook's name, with or without the ``.ipynb`` extension. Case is kept, so
        the folder sorts next to the notebook it belongs to.
    base : str or pathlib.Path, optional
        Where to put it. Defaults to the working directory — the study folder.

    Returns
    -------
    pathlib.Path
        The folder. Safe to call repeatedly; it does not clear what is already there.
    """
    stem = pathlib.Path(str(notebook)).stem
    if not stem:
        raise ValueError(f"cannot build a results folder from notebook={notebook!r}")
    folder = pathlib.Path(base or pathlib.Path.cwd()) / f"{stem}{RESULTS_SUFFIX}"
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def save_figure(fig, name, folder=None, dpi=DPI, close=True, verbose=False,
                vector=VECTOR_FORMAT, caption=None, module="a notebook"):
    """Write one finished figure into ``folder`` as a PNG and a vector copy, and close it.

    For a figure drawn in the house style at its final size (``house_style.new_figure``); a
    figure with several panels is written panel by panel with :func:`save_panels` instead.

    Parameters
    ----------
    fig : matplotlib.figure.Figure or None
        ``None`` is accepted and does nothing: a plotter that declined to draw returns it.
    name : str
        Filename without extension. ``.png`` is appended unless already present.
    folder : str or pathlib.Path, optional
        Defaults to the working directory, but pass the ``results_dir`` you made in the
        first cell.
    dpi : int
        Raster resolution, defaulting to :data:`DPI`.
    vector : str or None
        Extension of a vector copy written beside the PNG (:data:`VECTOR_FORMAT`).
    close : bool
        Close the figure once written (a notebook run saves dozens of figures).
    caption : str, optional
        LaTeX caption text; written beside the figure as ``<name>_caption.tex``.

    Returns
    -------
    pathlib.Path or None
        Where the PNG was written, or ``None`` if there was no figure.
    """
    if fig is None:
        return None

    folder = pathlib.Path(folder) if folder is not None else pathlib.Path.cwd()
    folder.mkdir(parents=True, exist_ok=True)

    name = str(name)
    path = folder / (name if name.endswith(".png") else f"{name}.png")
    stem = str(path)[:-4]                      # not with_suffix: '1e-7Pa' in a name is not an extension
    fig.savefig(path, dpi=dpi, facecolor="white")
    if vector:
        fig.savefig(f"{stem}.{vector.lstrip('.')}", facecolor="white")
    if caption is not None:
        from calculations.house_style import caption_tex
        with open(f"{stem}_caption.tex", "w") as fh:
            fh.write(caption_tex(os.path.basename(stem), module, caption))
    if close:
        plt.close(fig)
    if verbose:
        print(f"  saved figure -> {path}")
    return path


class Panel:
    """One figure of what used to be a multi-panel figure, written as <name>-<panel>.pdf, .png and
    _caption.tex. name: the panel's part of the file name ('flux', '600C-flux'); caption: LaTeX text
    with the numbers filled in; x, y: the data drawn, checked to lie inside the axes; keyed: the
    figure has a key below the plot (expected, not a problem)."""

    def __init__(self, name, fig, ax, caption, x=None, y=None, keyed=False):
        self.name, self.fig, self.ax, self.caption = name, fig, ax, caption
        self.x, self.y, self.keyed = x, y, keyed


def panel_files_pattern(base):
    """What an earlier run may have written for figure `base`: base.png/.pdf, and its panels
    base-<panel>.png/.pdf/_caption.tex (group 1: the panel's stem)."""
    return re.compile("(" + re.escape(base) + "|" + re.escape(base) + r"-.+?)(?:\.(?:png|pdf)|_caption\.tex)")


def save_panels(panels, name, folder=None, module="a notebook", log=print):
    """Write each panel of the figure ``name`` into ``folder`` and close it. Files an earlier run
    wrote for this figure that are not written now (the old one-file figure, panels no longer
    drawn) are removed. Problems found by the checks are printed, never raised. Returns the PNG
    paths, in drawing order."""
    from calculations.house_style import finish
    panels = [p for p in (panels or []) if p is not None]
    folder = pathlib.Path(folder) if folder is not None else pathlib.Path.cwd()
    folder.mkdir(parents=True, exist_ok=True)
    base = pathlib.Path(str(name)).name
    base = base[:-4] if base.endswith((".png", ".pdf")) else base
    keep = {f"{base}-{p.name}" for p in panels}
    old = panel_files_pattern(base)
    for f in os.listdir(folder):
        m = old.fullmatch(f)
        if m and m.group(1) not in keep:
            os.remove(folder / f)
    written = []
    for p in panels:
        files, problems = finish(p.fig, p.ax, str(folder / f"{base}-{p.name}"), module, p.caption,
                                 p.x, p.y, p.keyed)
        written.append(pathlib.Path(files[1]))
        for problem in problems:
            log(f"  ! {base}-{p.name}: {problem}")
    return written
