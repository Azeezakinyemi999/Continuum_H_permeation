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

and then, wherever a figure is finished::

    save_figure(fig, "Level1_Perfect_Metal", FIGS)

The sensitivity plotters in ``calculations.sensitivity`` all return their figure, so
they wrap directly::

    save_figure(plot_givendata_results(gd, reg, 'flux', PARAMS), f"givendata_{reg}", FIGS)
"""

import pathlib

import matplotlib.pyplot as plt

#: Appended to a notebook's name to make its results folder.
RESULTS_SUFFIX = "-results"

#: Raster resolution. Was 300, copied from the plt.savefig calls this module replaced.
#: Raised to 450 for sharper figures on a high-density display.
#:
#: There is a ceiling. Matplotlib refuses an image above 2**16 px in either direction,
#: which at 450 dpi is 145.6 inches -- down from 218 at 300 dpi. That is normally vast
#: next to a 14-inch figure, but `bbox_inches="tight"` sizes the output from the
#: ARTISTS, not from figsize, so one annotation on a near-flat log axis can ask for
#: thousands of inches. It has happened: the Level 4 grain-size panel demanded 7415.
DPI = 450

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
                vector=VECTOR_FORMAT):
    """Write ``fig`` into ``folder`` as a PNG and a vector copy, and close it.

    Parameters
    ----------
    fig : matplotlib.figure.Figure or None
        ``None`` is accepted and does nothing. The sensitivity plotters return it for a
        plot they declined to draw — ``plot_givendata_results`` when a regime was
        skipped, ``plot_regime_comparison_heatmap`` on an empty matrix — and wrapping
        those calls must not turn a skipped plot into a crash.
    name : str
        Filename without extension. ``.png`` is appended unless already present.
    folder : str or pathlib.Path, optional
        Defaults to the working directory, but pass the ``results_dir`` you made in the
        first cell — that is the point of this module.
    dpi : int
        Raster resolution, defaulting to :data:`DPI`.
    vector : str or None
        Extension of a vector copy written beside the PNG, defaulting to
        :data:`VECTOR_FORMAT`. ``None`` writes the raster alone. The vector file
        carries no dpi: it is resolution-independent, which is the point.
    close : bool
        Close the figure once written. On by default: a notebook run saves dozens of
        figures and nothing else closes them, so they would all stay open for the life
        of the kernel. The figure has normally been displayed by the plotter's own
        ``plt.show()`` before this is called, so closing costs nothing visible.

    Returns
    -------
    pathlib.Path or None
        Where the PNG was written, or ``None`` if there was no figure. The vector
        copy sits beside it with the same stem; the PNG path is returned because
        callers and tests treat it as the figure's identity.
    """
    if fig is None:
        return None

    folder = pathlib.Path(folder) if folder is not None else pathlib.Path.cwd()
    folder.mkdir(parents=True, exist_ok=True)

    name = str(name)
    path = folder / (name if name.endswith(".png") else f"{name}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    if vector:
        fig.savefig(path.with_suffix(f".{vector.lstrip('.')}"),
                    bbox_inches="tight", facecolor="white")
    if close:
        plt.close(fig)
    if verbose:
        print(f"  saved figure -> {path}")
    return path
