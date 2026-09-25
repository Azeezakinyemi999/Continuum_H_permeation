"""Active study configuration.

This module is a LOADER, not a config. It loads exactly one study and re-exports its
names, and everything else in the codebase imports from here::

    calculations/*.py  ->  model_config  ->  the one active study

Each study IS a folder. ``<Folder>/study_config.py`` holds that study's real input
data — material properties, operating conditions, sensitivity ranges and presets — not
a pointer to data kept elsewhere. Which study loads is decided by LOCATION:

1. ``$MHI_STUDY`` — a path to a study file. For the container gate, the docs tools, CI,
   and any run that is not inside a study folder.
2. ``study_config.py`` in the current working directory, then in each parent up to the
   repository root. **This is the normal case for notebooks**: a study folder holds its
   own config, and running a notebook in that folder uses it. No switch to flip.
3. Nothing else. There is no default — a run that named no study raises, listing the
   folders it could have meant. A fallback would be taken by exactly the runs that
   never considered which material they wanted.

The working directory therefore decides the physics, and that is worth one moment of
care: a notebook run from the wrong directory computes a different material than its
folder implies. Three things guard against it — the one-line announcement on stderr at
import, ``ACTIVE_STUDY_SOURCE`` recording HOW the study was chosen, and
``calculations.config.assert_study_is_local()``, which notebooks call in their first
cell to turn a wrong directory into a stack trace instead of a wrong number.

In VS Code, set ``"jupyter.notebookFileRoot": "${fileDirname}"`` or every notebook runs
from the workspace root, where no study applies and the import fails.
``.vscode/settings.json`` in this repo does that. Jupyter Lab already uses the
notebook's own directory.

Restart any running kernel after changing which study applies: the binding happens at
import time, so ``importlib.reload`` on a notebook module is not enough — the study's
dicts are already bound into calculations.sensitivity.

To add a study: make a folder, put its notebooks and a ``study_config.py`` in it.
See calculations/config/__init__.py for the full recipe.

EXACTLY ONE study is ever loaded. This used to be a literal ``from studies.X import *``
with ``ACTIVE_STUDY`` as a separate label describing it, which let the two disagree —
and they did: the label read one study while the import loaded another, and at one
point BOTH studies were star-imported at once, so Python applied them in order and
produced a hybrid. That state had Hastelloy N's material dicts with Incoloy's derived
sensitivity parameters, and it moved L1 flux by a factor of ten while ``ACTIVE_STUDY``
still claimed Incoloy. The resolution order above has exactly one winner, and the
winner is reported, so neither failure can recur.

After changing study, run:

    from calculations.sensitivity import check_against_config
    check_against_config()

Preset yields, sweep ranges and draw counts were tuned for a specific material and
are not automatically valid for another one.
"""

import importlib.util as _importlib_util
import os as _os
import pathlib as _pathlib
import sys as _sys

# =============================================================================
# WHERE STUDIES LIVE
# =============================================================================
# There is no default study and no switch to edit. Each study folder holds its own
# study_config.py — the real input file, not a pointer — and that file is selected by
# running in its folder, or named explicitly through MHI_STUDY. Anything run from the
# repository root must name one; see _no_study_found below for why silence is worse.

# The filename a study folder must use. Deliberately NOT 'model_config.py': a notebook
# kernel has its own folder on sys.path, so a folder-local model_config.py would shadow
# a bare `import model_config` and give two different modules the same name.
STUDY_FILENAME = 'study_config.py'

# Override variable. Read once, at import.
STUDY_ENV_VAR = 'MHI_STUDY'

# calculations/config/model_config.py -> config -> calculations -> repo root.
# The folder search never looks above this, so a stray study_config.py in a home
# directory can never be picked up.
_REPO_ROOT = _pathlib.Path(__file__).resolve().parents[2]

# =============================================================================


def _available_folder_studies():
    """Repo-relative paths of every ``study_config.py`` in a repo-root subdirectory.

    Only one level down — that is where the study folders live, and a full walk would
    descend into _build/, __pycache__/ and the sa_results*/ output trees.
    """
    found = []
    for entry in sorted(_REPO_ROOT.iterdir()):
        if entry.is_dir() and not entry.name.startswith(('.', '_')):
            candidate = entry / STUDY_FILENAME
            if candidate.is_file():
                found.append(str(candidate.relative_to(_REPO_ROOT)))
    return found


# Names any study must provide, derived from what `calculations/` actually imports
# rather than from whatever one study happens to define. The first seven are pulled
# at MODULE level by calculations/sensitivity.py, so a study missing one of them
# fails at `import calculations.sensitivity` — before any error message of ours can
# help. Checking here turns that into a sentence naming the study and the gap.
REQUIRED_EXPORTS = frozenset({
    # sensitivity.py module-level imports
    'DEFAULT_PARAMS_LEVEL5', 'DEFAULT_PARAMS_LEVEL5L6',
    'SUGGESTED_RANGES_LEVEL5', 'SUGGESTED_RANGES_LEVEL5L6',
    'REGIME_PRESETS', 'REGIME_PRESETS_L5',
    # material / condition dicts the model and the notebooks read
    'METALS', 'OXIDES', 'MICROSTRUCTURE', 'OXIDE_DEFECTS', 'CONDITIONS',
    'build_simulation_config',
    # presentation
    'PLOT_STYLE', 'COLORS', 'CURVE_STYLES',
})


def _resolve_path(value):
    """A study path as given, resolved from the REPO ROOT when relative.

    Not from the working directory: the same MHI_STUDY value has to mean the same
    study whether it was set for a notebook in a study folder or for a tool at the
    root, and a cwd-relative path would not.
    """
    path = _pathlib.Path(value).expanduser()
    if not path.is_absolute():
        path = _REPO_ROOT / path
    return path.resolve()


def _no_study_found():
    """Raise, rather than guess, when nothing named a study.

    There is deliberately no fallback. A default would be loaded by exactly the runs
    that did not think about which material they wanted — scripts and tools started
    from the repository root — and it would be loaded silently, producing numbers for
    one study under the name of whatever the caller assumed.
    """
    folders = _available_folder_studies()
    raise ValueError(
        f"no study selected. Nothing was found for the working directory "
        f"{_pathlib.Path.cwd()}, and there is no default.\n"
        f"Either run from a study folder (one that contains a {STUDY_FILENAME}), "
        f"or name one explicitly:\n"
        + "".join(f"    {STUDY_ENV_VAR}={f} python ...\n" for f in folders or
                  [f'<folder>/{STUDY_FILENAME}'])
    )


def _resolve_study():
    """Decide which study to load, and say where the decision came from.

    Returns ``(path, source)`` where ``source`` is 'env' or 'folder' and ``path`` is a
    ``pathlib.Path`` to a study file. Raises when nothing selected a study.
    """
    override = _os.environ.get(STUDY_ENV_VAR, '').strip()
    if override:
        return _resolve_path(override), 'env'

    # Walk up from the working directory, but only while still inside the repo.
    try:
        here = _pathlib.Path.cwd().resolve()
    except OSError:                      # deleted cwd — rare, but do not crash the import
        here = None
    if here is not None and (here == _REPO_ROOT or _REPO_ROOT in here.parents):
        for directory in (here, *here.parents):
            candidate = directory / STUDY_FILENAME
            if candidate.is_file():
                return candidate, 'folder'
            if directory == _REPO_ROOT:
                break

    _no_study_found()


def _exported_names(mod):
    """Star-import semantics: every non-underscore name, honouring __all__.

A star import cannot be written with a variable — it is a static statement whose
    module path is resolved at compile time — so the re-export is done by hand,
    reproducing those semantics exactly.

    Note that includes ``np``, which the studies import for their VALIDATION arrays.
    The original star import leaked it into this namespace too; it is kept so that
    nothing observable changed when this replaced it.
    """
    exported = getattr(mod, '__all__', None)
    if exported is None:
        exported = [n for n in dir(mod) if not n.startswith('_')]
    return list(exported)


def _check_required(label, exported):
    """Raise unless the study defines every name `calculations/` imports from it.

    Shared by the active study and by load_study, so reading a second study cannot
    end up enforcing a weaker contract than loading one.
    """
    missing = REQUIRED_EXPORTS - set(exported)
    if missing:
        raise ValueError(
            f"study {label!r} is missing {len(missing)} required export(s): "
            f"{sorted(missing)}. Every study must define these — copy them from an "
            f"existing {STUDY_FILENAME} and adapt the values "
            f"to the material. Derived dicts (DEFAULT_PARAMS_*) should be copied as "
            f"CODE so they recompute from that study's own METALS/OXIDES; only the "
            f"study-design literals (ranges, presets) are written out by hand."
        )


def _import_study_file(path):
    """Execute a study file and return it as a module.

    Registered in sys.modules under a private name so repeated imports and anything
    that pickles a value defined here stay coherent. Never named 'model_config'.
    """
    if not path.is_file():
        raise ValueError(
            f"study file {str(path)!r} does not exist. "
            f"Available studies: {_available_folder_studies()}."
        )

    mod_name = f'{__package__}._file_study_{path.stem}'
    spec = _importlib_util.spec_from_file_location(mod_name, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"cannot load a Python module from {str(path)!r}")

    mod = _importlib_util.module_from_spec(spec)
    _sys.modules[mod_name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        # A half-executed study must not stay in sys.modules pretending to be loaded.
        _sys.modules.pop(mod_name, None)
        raise
    return mod


def _load_active_study(path):
    """Load the study file at ``path`` and return ``(module, namespace)``."""
    mod = _import_study_file(path)
    exported = _exported_names(mod)
    _check_required(str(path), exported)
    return mod, {n: getattr(mod, n) for n in exported}


def load_study(spec):
    """Load any study WITHOUT making it active, and return it as a module.

    The supported way to read a second study — for a side-by-side comparison, or to
    pass another material registry to a property getter:

    >>> other = load_study('Fuerst_et_al_2024/study_config.py')   # doctest: +SKIP
    >>> get_metal_properties_at_T('Hastelloy_N_fuerst_2024', 900.0,
    ...                           metals=other.METALS)            # doctest: +SKIP

    Parameters
    ----------
    spec : str or pathlib.Path
        Path to a study file. A relative path resolves from the repository root, not
        from the working directory, so the same call works from anywhere.

    Returns
    -------
    module
        The study, validated against REQUIRED_EXPORTS exactly as the active one is.

    Notes
    -----
    This does not change what the solvers use. They bound their defaults at import
    time from whichever study was resolved then.
    """
    mod, _ = _load_active_study(_resolve_path(spec))
    return mod


_study_path, ACTIVE_STUDY_SOURCE = _resolve_study()
_study_module, _study_ns = _load_active_study(_study_path)
globals().update(_study_ns)

# ---------------------------------------------------------------------------
# Provenance. A consequence of the load above, not a parallel claim about it —
# anything labelling a figure, a CSV or a reference set should read these.
# ---------------------------------------------------------------------------

#: Absolute path of the study file that was loaded.
ACTIVE_STUDY_PATH = _pathlib.Path(_study_module.__file__).resolve()


def _study_label(mod, path):
    """The short name of a study, in order of how much it can be trusted.

    1. ``STUDY_NAME`` declared by the study itself. A folder shim star-imports the
       shared study module, so it inherits that module's name — which is what keeps
       reference-value lookups (container/model_smoke_test.py) working no matter
       which folder the run was launched from.
    2. The folder name, when the file is a ``study_config.py`` with no declared name.
       The stem would be the literal 'study_config' for every folder study.
    3. The file stem, which for a module in studies/ is its module name.
    """
    declared = str(getattr(mod, 'STUDY_NAME', '') or '').strip()
    if declared:
        return declared
    if path.name == STUDY_FILENAME:
        return path.parent.name
    return path.stem


#: Short label for figures, CSV provenance and per-study reference sets.
ACTIVE_STUDY_NAME = _study_label(_study_module, ACTIVE_STUDY_PATH)

#: Backwards compatibility. Notebooks, docs/_tools/worked_values.py and
#: container/model_smoke_test.py read ACTIVE_STUDY; for a studies/ module its value is
#: unchanged. A folder config reports the stem of its own file.
ACTIVE_STUDY = ACTIVE_STUDY_NAME

#: Dotted name of the loaded module (private for file-based studies).
ACTIVE_STUDY_MODULE = _study_module.__name__


def study_provenance():
    """One line naming the active study, where it came from, and how it was chosen."""
    try:
        where = ACTIVE_STUDY_PATH.relative_to(_REPO_ROOT)
    except ValueError:
        where = ACTIVE_STUDY_PATH
    return f"study: {ACTIVE_STUDY_NAME}  ({where}, via {ACTIVE_STUDY_SOURCE})"


# Announce the choice. On STDERR deliberately: stdout is piped to a file by
# docs/_tools/worked_values.py, and a banner there would end up inside the generated
# reference output. Jupyter shows stderr in the cell, which is where it is wanted.
if not _os.environ.get('MHI_STUDY_QUIET'):
    print(study_provenance(), file=_sys.stderr)
