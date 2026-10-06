"""Study configuration: how the active study is chosen.

No study lives in this package. **A study is a folder** at the repository root: its
notebooks and its input file, ``study_config.py``, together. That file is the real
thing — material properties, operating conditions, microstructure, defects, sweep
ranges and regime presets — so editing it changes that study and nothing else::

    Guo_etal_2025/study_config.py        316L stainless / Cr2O3
    Fuerst_et_al_2024/study_config.py    Hastelloy N / Cr2O3
    incoloy802_cr2o3/study_config.py     Incoloy 802 (X40 NiCrAlTi 31/19) / Cr2O3

``model_config`` resolves exactly one of them per process, from ``$PERMEATION_STUDY``, then
from a ``study_config.py`` in the working directory or a parent — and from nothing
else, so a run that named no study raises rather than guessing.
``assert_study_is_local`` below is the guard a notebook uses to insist it got the one
in its OWN folder, not a parent's.

To add a study
--------------
1. Make a folder at the repository root and put the notebooks in it.
2. Copy an existing ``study_config.py`` into it and edit METALS / OXIDES / CONDITIONS /
   MICROSTRUCTURE / OXIDE_DEFECTS, and the SENSITIVITY ANALYSIS PARAMETERS section if
   the ranges or presets should differ.
3. Set ``STUDY_NAME`` to something that identifies the study. It labels figures and
   CSVs and selects the reference set in ``container/model_smoke_test.py``, and it is
   deliberately independent of the folder name, so renaming a folder cannot silently
   relabel results.
4. Run the notebooks from inside that folder. ``assert_study_is_local()`` in their
   first cell fails if the kernel started anywhere else.
5. Run ``calculations.sensitivity.check_against_config()`` before trusting results —
   the preset yields and sweep ranges were tuned for a specific material and will not
   automatically be right for another one.
6. Measure a reference set and add it to ``REFERENCE`` in
   ``container/model_smoke_test.py``, or the container release gate fails on the new
   study — deliberately, since passing without one would check nothing.
7. If the study is to ship in a container release, add its folder to the assembled
   workspace: a study folder outside the assembled set will not travel with the image.
   As of 1.2.0 that set is ``calculations/`` plus ``Guo_etal_2025/``,
   ``Fuerst_et_al_2024/`` and ``incoloy802_cr2o3/``. (PACKAGING_1.0.0.md describes the
   1.0.0 set, ``calculations/`` plus ``Application/``; ``Application/`` was a delivery
   copy of the Guo study and was withdrawn on 2026-09-28.)

Why the file is named study_config.py
-------------------------------------
NOT ``model_config.py``. ``from calculations.config.model_config import ...`` is an
ABSOLUTE package import, but a notebook kernel also has its own folder on sys.path, so
a folder-local ``model_config.py`` would shadow a bare ``import model_config`` and give
two different modules the same name — the notebook reporting one material while the
model computed another. A distinct filename cannot collide, and ``model_config`` stays
the only loader.
"""

import pathlib as _pathlib


def assert_study_is_local(directory=None):
    """Fail unless the active study came from a ``study_config.py`` in ``directory``.

    Call this in a notebook's first cell. Without it, a kernel started in the wrong
    directory — the VS Code default is the workspace root, not the notebook's folder —
    either fails to find any study at all, or picks up a parent folder's one and computes
    a different material than the folder implies, writing the result into that folder's
    ``sa_results*/``. With it, the same mistake is a stack trace before any physics runs.

    Parameters
    ----------
    directory : str or pathlib.Path, optional
        The folder whose config must be the active one. Defaults to the working
        directory, which for a correctly launched notebook is the notebook's folder.

    Returns
    -------
    str
        The provenance line, so the cell that calls this also displays it.

    Raises
    ------
    RuntimeError
        If the study was not resolved from a folder config, or was resolved from a
        different folder than ``directory``.
    """
    from calculations.config.model_config import (
        ACTIVE_STUDY_NAME, ACTIVE_STUDY_PATH, ACTIVE_STUDY_SOURCE,
        STUDY_FILENAME, study_provenance,
    )

    want = _pathlib.Path(directory).resolve() if directory else _pathlib.Path.cwd().resolve()

    if ACTIVE_STUDY_SOURCE != 'folder':
        raise RuntimeError(
            f"expected the study in {want / STUDY_FILENAME}, but the "
            f"active study {ACTIVE_STUDY_NAME!r} was resolved via "
            f"{ACTIVE_STUDY_SOURCE!r} ({ACTIVE_STUDY_PATH}).\n"
            f"The working directory is {_pathlib.Path.cwd()}. In VS Code set "
            f'"jupyter.notebookFileRoot": "${{fileDirname}}" so the kernel starts in '
            f"the notebook's own folder, then restart the kernel — reloading is not "
            f"enough, the study is bound at import."
        )

    if ACTIVE_STUDY_PATH.parent != want:
        raise RuntimeError(
            f"the active study came from {ACTIVE_STUDY_PATH}, not from {want}. A parent "
            f"directory's {STUDY_FILENAME} was picked up because this folder has none. "
            f"Add one, or start the kernel in the folder you meant."
        )

    return study_provenance()
