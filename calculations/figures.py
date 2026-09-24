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

#: Defaults copied from the existing plt.savefig calls in Proposal.ipynb, so figures
#: that already existed come out unchanged.
DPI = 300


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


def save_figure(fig, name, folder=None, dpi=DPI, close=True, verbose=False):
    """Write ``fig`` into ``folder`` as a PNG, and close it.

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
        Resolution. 300 matches the figures this repo already produced.
    close : bool
        Close the figure once written. On by default: a notebook run saves dozens of
        figures and nothing else closes them, so they would all stay open for the life
        of the kernel. The figure has normally been displayed by the plotter's own
        ``plt.show()`` before this is called, so closing costs nothing visible.

    Returns
    -------
    pathlib.Path or None
        Where it was written, or ``None`` if there was no figure.
    """
    if fig is None:
        return None

    folder = pathlib.Path(folder) if folder is not None else pathlib.Path.cwd()
    folder.mkdir(parents=True, exist_ok=True)

    name = str(name)
    path = folder / (name if name.endswith(".png") else f"{name}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    if close:
        plt.close(fig)
    if verbose:
        print(f"  saved figure -> {path}")
    return path
