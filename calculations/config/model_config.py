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

import numpy as np

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


# Names any study must provide. Derived from what the codebase actually imports from
# this module — not from whatever one study happens to define — so the check fails with
# a sentence naming the study and the gap, rather than an ImportError from deep inside
# another module.
#
# It had drifted: five names were imported somewhere but absent here, so a study missing
# one of them passed the contract and then failed later anyway. The set is regenerated by
# tests/test_operating_conditions.py, which re-derives it from the imports and fails if
# this list is short again — including for notebooks, which is how the gap arose.
# There are two reasons a name lands here, and they are not the same reason. Most are
# required because something imports them from this module. The first group is required
# because the DERIVED block at the bottom of THIS file reads them out of the study
# namespace — nothing imports MODEL_OPTIONS anywhere, but a study without it would fail
# as a bare NameError raised from inside the loader, which is precisely the failure this
# check exists to replace.
REQUIRED_EXPORTS = frozenset({
    # read by the DERIVED block below, out of the loaded study's namespace
    'CONDITIONS', 'METALS', 'MICROSTRUCTURE', 'MODEL_OPTIONS', 'OXIDES', 'OXIDE_DEFECTS',
    'R',
    # calculations/sensitivity.py, at module level
    'SUGGESTED_RANGES_LEVEL5', 'SUGGESTED_RANGES_LEVEL5L6',
    'REGIME_PRESETS', 'REGIME_PRESETS_L5',
    # the notebooks import these from this module
    'F',
    'DEFAULT_SWEEP_TEMPERATURES',
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

    The studies used to import numpy, and it leaked through here into this namespace.
    They no longer need it: the only code that used it was the derivation, which now
    lives below. ``np`` is still importable from this module because it is imported
    here directly, so nothing observable changed.
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
            f"existing {STUDY_FILENAME} and adapt the values to the material. "
            f"A study file holds its own data only: the material registries, the "
            f"operating point, the switches and the study-design literals (ranges, "
            f"presets). The derived parameter sets, build_simulation_config and the "
            f"plot styling are NOT copied — they live once in model_config.py and "
            f"recompute from whichever study is loaded."
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

# =============================================================================
# EVERYTHING BELOW IS DERIVED FROM THE STUDY LOADED ABOVE.
#
# It used to be sections 4 and 5 of every study_config.py, byte-identical in all four
# — 408 lines duplicated four ways, and a copy that can drift is a copy that will. The
# Level 5 header below still named Incoloy 802 in the 316L study because of exactly
# that. Held once, "a literal here is a bug" is true by construction.
#
# THE ORDER MATTERS. This reads CONDITIONS, METALS, MICROSTRUCTURE, MODEL_OPTIONS,
# OXIDES, OXIDE_DEFECTS and R as module globals, and they exist only because of the
# globals().update() immediately above. Moving this block above that line makes every
# one of them a NameError.
# =============================================================================

# ==============================================================================
# 4. DERIVED
#
# Assembled from sections 2 and 3 by reference, never by retyping. This is why the section
# comes last: Python executes top to bottom, so everything these read must already exist.
#
# A literal here is a bug. It is a second copy of a number that lives above, and the two
# will diverge — they already did once, by up to six orders of magnitude.
# ==============================================================================

# The flat parameter sets the model wrappers and the sensitivity analysis consume,
# assembled from sections 2 and 3. Run calculations.sensitivity.check_against_config()
# after editing anything they read.

# Which METALS / OXIDES entry the sensitivity analysis characterises: THE FIRST
# ENTRY of each dict. Taken positionally rather than spelled out, so adding or
# renaming a material needs no edit here and the name cannot fall out of sync with
# the dict it indexes.
#
# The convention this encodes: the first entry is the active one. Dicts preserve
# insertion order (Python 3.7+), so that means the first literal in the file. If you
# add a second material, put it BELOW the active one — inserting above silently
# switches every derived default.


def _first_key(d, what):
    """First key of ``d``, with a clear error instead of a bare StopIteration."""
    try:
        return next(iter(d))
    except StopIteration:
        raise ValueError(f"{what} is empty — cannot determine the active material") from None

ACTIVE_METAL = _first_key(METALS, 'METALS')

ACTIVE_OXIDE = _first_key(OXIDES, 'OXIDES')

_SA_M  = METALS[ACTIVE_METAL]

_SA_O  = OXIDES[ACTIVE_OXIDE]

_SA_MS = _SA_M['surface_kinetics']          # metal surface kinetics

_SA_OS = _SA_O['surface_kinetics']          # oxide surface kinetics

_SA_DC = OXIDE_DEFECTS['components']

def _sa_trap(name):
    """MICROSTRUCTURE trap entry by name (case-insensitive — config uses 'Carbides')."""
    for t in MICROSTRUCTURE['trap_list']:
        if t['name'].lower() == name.lower():
            return t
    raise KeyError(f"no trap named {name!r} in MICROSTRUCTURE['trap_list']; "
                   f"have {[t['name'] for t in MICROSTRUCTURE['trap_list']]}")

# -----------------------------------------------------------------------------
# LEVEL 5 defaults — defective oxide (L3) + defective metal microstructure (L4)
#
# Reference-point Arrhenius: X(T) = X_ref * exp(-E/R * (1/T - 1/T_ref))
#
# T_ref_metal and T_ref_oxide are whatever the loaded study measured its own *_ref
# values at — they are read from the material entries below, not stated here. This
# comment used to name Incoloy 802 and Cr2O3_sample4 as though they were the
# materials, which was true of one study and copied verbatim into the other three.
# -----------------------------------------------------------------------------

DEFAULT_PARAMS_LEVEL5 = {
    # Metal transport  (METALS[ACTIVE_METAL])
    'D_ref':       _SA_M['D_ref'],      # m²/s           D_metal at T_ref_metal
    'E_D':         _SA_M['E_D'],        # J/mol          activation energy for diffusion
    'K_s_ref':     _SA_M['K_s_ref'],    # mol/m³/Pa^0.5  K_s at T_ref_metal
    'H_s':         _SA_M['H_s'],        # J/mol          heat of solution
    'T_ref_metal': _SA_M['T_ref'],      # K              measurement reference

    # Oxide transport  (OXIDES[ACTIVE_OXIDE])
    'D_ox_ref':    _SA_O['D_ox_ref'],   # m²/s
    'E_D_ox':      _SA_O['E_D_ox'],     # J/mol
    'K_ox_ref':    _SA_O['K_ox_ref'],   # mol/m³/Pa^0.5
    'H_sol_ox':    _SA_O['H_sol_ox'],   # J/mol
    'T_ref_oxide': _SA_O['T_ref'],      # K

    # Geometry — both layers from CONDITIONS, which describes the wall being simulated.
    # The materials carry only a sample_thickness apiece, recording the specimens their
    # data were measured on.
    'metal_thickness': CONDITIONS['L_metal'],   # m
    'oxide_thickness': CONDITIONS['L_oxide'],   # m

    # Operating conditions  (CONDITIONS)
    'P_upstream':   CONDITIONS['P_upstream'],   # Pa
    'P_downstream': CONDITIONS['P_downstream'], # Pa
    'temperature':  CONDITIONS['T_operating'],  # K

    # Level 3: oxide defects  (OXIDE_DEFECTS)
    'f_pinhole':              _SA_DC['pinholes'],          # area fraction
    'f_crack':                _SA_DC['cracks'],            # area fraction
    'f_gb_defect':            _SA_DC['grain_boundaries'],  # area fraction
    'crack_thickness_factor': OXIDE_DEFECTS['thickness_factor'],    # L_crack = f × L_oxide
    'gb_diffusivity_factor':  OXIDE_DEFECTS['diffusivity_factor'],  # D_gb = f × D_oxide
    # 'use_sieverts_pinhole' is NOT here: neither model wrapper reads it from the
    # parameter dict. Surface_proposal.ipynb takes it from SIM['oxide_defects'],
    # i.e. straight from OXIDE_DEFECTS, so that is where it lives.

    # Level 4: grain structure  (MICROSTRUCTURE; Zhu 2021 baseline)
    'grain_size':      MICROSTRUCTURE['grain_size'],       # m
    'grain_shape':     MICROSTRUCTURE['grain_shape'],
    'gb_type':         MICROSTRUCTURE['gb_type'],
    'gb_thickness':    MICROSTRUCTURE['gb_thickness'],     # m
    'lattice_density': MICROSTRUCTURE['lattice_density'],  # m⁻³

    # Level 4: traps  (MICROSTRUCTURE['trap_list'])
    # Lu 2022 binding energies + Young 1997 carbide density
    'trap_dislocation_E_b': _sa_trap('dislocations')['binding_energy'],      # J/mol
    'trap_dislocation_N_T': _sa_trap('dislocations')['density'],             # m⁻³
    'trap_gb_E_b':          _sa_trap('grain_boundaries')['binding_energy'],  # J/mol
    'trap_gb_N_T':          _sa_trap('grain_boundaries')['density'],         # m⁻³
    'trap_vacancy_E_b':     _sa_trap('vacancies')['binding_energy'],         # J/mol
    'trap_vacancy_N_T':     _sa_trap('vacancies')['density'],                # m⁻³
    'trap_carbide_E_b':     _sa_trap('carbides')['binding_energy'],          # J/mol
    'trap_carbide_N_T':     _sa_trap('carbides')['density'],                 # m⁻³

    # Model options — how to run the model, not what it is; no config source
    'enable_gb_enhancement': MODEL_OPTIONS['enable_gb_enhancement'],
    'enable_trapping':       MODEL_OPTIONS['enable_trapping'],
    'D_eff_method':           MODEL_OPTIONS['D_eff_method'],
}

# -----------------------------------------------------------------------------
# LEVEL 5L6 defaults — adds oxide + metal surface dissociation/recombination.
# Inherits every Level 5 parameter via **DEFAULT_PARAMS_LEVEL5, which is what
# keeps the two levels sharing one set of values for everything they have in common.
# -----------------------------------------------------------------------------

DEFAULT_PARAMS_LEVEL5L6 = {
    **DEFAULT_PARAMS_LEVEL5,
    # Oxide surface kinetics — OXIDES[ACTIVE_OXIDE]['surface_kinetics'] (Grant 1988)
    'k_diss_ref':    _SA_OS['k_diss_ref'],   # mol/m²/s/Pa
    'E_diss':        _SA_OS['E_diss'],       # J/mol
    'K_eq_ref':      _SA_OS['K_eq_ref'],     # Pa⁻¹
    'H_eq':          _SA_OS['H_eq'],         # J/mol
    'T_ref_surface': _SA_OS['T_ref'],        # K

    # Metal surface kinetics — METALS[ACTIVE_METAL]['surface_kinetics'] (Grant 1988)
    'k_diss_metal_ref':    _SA_MS['k_diss_metal_ref'],  # mol/m²/s/Pa
    'E_diss_metal':        _SA_MS['E_diss_metal'],      # J/mol
    'K_eq_metal_ref':      _SA_MS['K_eq_metal_ref'],    # Pa⁻¹
    'H_eq_metal':          _SA_MS['H_eq_metal'],        # J/mol
    'T_ref_surface_metal': _SA_MS['T_ref'],             # K
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_surface_kinetics_from_config(material_key, temperature_K, material_dict):
    r"""
    Compute temperature-dependent surface kinetics from a METALS or OXIDES entry.

    Uses the 'surface_kinetics' sub-dict within the material entry.
    Arrhenius form:

    .. math::

        k_{diss}(T) &= k_{diss,ref}\exp\!\left[\frac{-E_{diss}}{R}
            \left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right] \\
        K_{eq}(T) &= K_{eq,ref}\exp\!\left[\frac{-H_{eq}}{R}
            \left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right] \\
        k_{recomb} &= \frac{k_{diss}}{K_{eq}}

    Parameters
    ----------
    material_key  : str   — key in METALS or OXIDES
    temperature_K : float — temperature [K]
    material_dict : dict  — either METALS or OXIDES

    Returns
    -------
    dict : k_diss, k_recomb, K_eq
    """
    entry = material_dict[material_key]
    sk    = entry['surface_kinetics']
    T_ref = sk['T_ref']

    # Oxide surface kinetics keys
    if 'k_diss_ref' in sk:
        k_diss_ref = sk['k_diss_ref']
        E_diss     = sk['E_diss']
        K_eq_ref   = sk['K_eq_ref']
        H_eq       = sk['H_eq']
    # Metal surface kinetics keys
    else:
        k_diss_ref = sk['k_diss_metal_ref']
        E_diss     = sk['E_diss_metal']
        K_eq_ref   = sk['K_eq_metal_ref']
        H_eq       = sk['H_eq_metal']

    inv_T_diff = 1.0 / temperature_K - 1.0 / T_ref
    k_diss  = k_diss_ref * np.exp((-E_diss / R) * inv_T_diff)
    K_eq    = K_eq_ref   * np.exp((-H_eq   / R) * inv_T_diff)
    k_recomb = k_diss / K_eq

    return {'k_diss': k_diss, 'k_recomb': k_recomb, 'K_eq': K_eq}

def build_simulation_config(
    metal=ACTIVE_METAL,   # first entry of METALS — see ACTIVE_METAL above
    oxide=ACTIVE_OXIDE,   # first entry of OXIDES
    T_operating=None,
    P_upstream=None,
    L_metal=None,
    L_oxide=None,
    microstructure_overrides=None,
    defect_overrides=None,
):
    """Assemble a solver-ready configuration dict for one metal/oxide pair.

    Pulls the chosen entries out of METALS and OXIDES and merges them with the
    operating conditions in CONDITIONS, so a caller hands the model a single dict
    instead of threading a dozen parameters through by hand.

    The oxide defect fractions in OXIDE_DEFECTS are also translated here into the
    per-path ``defect_config`` format the Level 3 solvers expect. A pinhole, crack
    or grain-boundary path is included only when its area fraction is non-zero, so
    setting a fraction to zero removes that transport path entirely rather than
    leaving a zero-weighted one in the parallel sum.

    Parameters
    ----------
    metal, oxide : str
        Keys into METALS and OXIDES. Default to ACTIVE_METAL / ACTIVE_OXIDE, which
        are the first entry of each dict.
    T_operating, P_upstream, L_metal, L_oxide : float, optional
        Override the corresponding CONDITIONS value. CONDITIONS is the only source
        for both thicknesses; the materials carry a sample_thickness recording the
        specimen their data came from, which is not the wall and is never read here.
        Note these use ``or``, so a falsy 0 also takes the default — pass a small
        positive number instead.
    microstructure_overrides, defect_overrides : dict, optional
        Shallow-merged onto MICROSTRUCTURE and OXIDE_DEFECTS respectively. The
        defect overrides are applied BEFORE ``defect_config`` is derived, so they
        do change which paths exist.

    Returns
    -------
    dict
        Material properties, operating conditions and sweep ranges, plus both the
        raw 'oxide_defects' block and the derived 'defect_config'.
    """
    od = OXIDE_DEFECTS.copy()
    if defect_overrides:
        od.update(defect_overrides)

    # Build solver-ready defect config from OXIDE_DEFECTS structure
    defect_config = {}
    f_pin = od['components'].get('pinholes', 0.0)
    f_cra = od['components'].get('cracks', 0.0)
    f_gb  = od['components'].get('grain_boundaries', 0.0)
    if f_pin > 0:
        defect_config['pinhole'] = {'area_fraction': f_pin}
    if f_cra > 0:
        defect_config['crack'] = {
            'area_fraction':    f_cra,
            'thickness_factor': od.get('thickness_factor', 0.1),
        }
    if f_gb > 0:
        defect_config['grain_boundary'] = {
            'area_fraction':      f_gb,
            'diffusivity_factor': od.get('diffusivity_factor', 100.0),
        }

    config = {
        'metal_name':  metal,
        'oxide_name':  oxide,
        'metal_props': METALS[metal].copy(),
        'oxide_props': OXIDES[oxide].copy(),

        'T_operating':  T_operating or CONDITIONS['T_operating'],
        'P_upstream':   P_upstream  or CONDITIONS['P_upstream'],
        'P_downstream': CONDITIONS['P_downstream'],

        'L_metal': L_metal or CONDITIONS['L_metal'],
        'L_oxide': L_oxide or CONDITIONS['L_oxide'],

        'microstructure': MICROSTRUCTURE.copy(),
        'oxide_defects':  od,              # raw OXIDE_DEFECTS (with any overrides)
        'defect_config':  defect_config,   # ← solver-ready format, built automatically

        'T_range':       CONDITIONS['T_range'],
        'P_range':       CONDITIONS['P_range'],
        'L_metal_range': CONDITIONS['L_metal_range'],
        'L_oxide_range': CONDITIONS['L_oxide_range'],
        'n_points':      CONDITIONS['n_points'],
    }

    if microstructure_overrides:
        config['microstructure'].update(microstructure_overrides)

    return config


# ==============================================================================
# 5. PRESENTATION
#
# Figure styling. Affects no number.
# ==============================================================================

# =============================================================================
# PLOTTING STYLE
# =============================================================================
PLOT_STYLE = {
    'figsize':              (14, 12),
    # Panel titles used to be 10 against 14 for axis labels -- a panel's title read
    # SMALLER than the label beneath it, which no amount of resolution fixes. Sizes
    # raised throughout and the ordering corrected: suptitle > title >= axis > tick.
    'fontsize_title':       15,
    'fontsize_suptitle':    20,
    'fontsize_axis':        15,
    'fontsize_tick':        13,
    'fontsize_legend':      13,
    'fontsize_annotation':  12,
    # Marks stay as they were. The ask was for text and colour; thickening the curves
    # too would crowd a four-panel figure, and a heavy line hides where two levels
    # coincide -- which is exactly what the limit-check panels exist to show.
    'linewidth':            2.5,
    'markersize':           9,
    'grid_alpha':           0.25,
}

# Hues chosen against the data-visualisation palette and CHECKED, not eyeballed: the
# validator reports lightness band, chroma floor, contrast, and pairwise separation
# under normal vision and under protan/deutan/tritan colour-vision deficiency.
#
# What the old set got wrong. 'cyan' is #00FFFF -- almost invisible on white. 'crimson'
# and 'red' sat on L4_both and L5, which share a panel, at a normal-vision separation
# far under the floor. And 'orange' against 'green', the L4 trapping and GB curves on
# one axes, came to deltaE 3.2 for a protan reader: indistinguishable.
#
# Checked per PANEL rather than globally, because what must be separable is what is
# drawn together. The binding case is the Level 5 summary, six hues at once, which
# clears the hard normal-vision floor at 15.6 with worst CVD 6.9 -- legal because every
# level also carries its own linestyle and marker in CURVE_STYLES below. L2a may take
# orange, whose only conflicts are with red and green, because it appears beside
# neither. L1 stays black: it is the baseline every other curve is measured against,
# and reads as reference ink rather than as one series among many.
COLORS = {
    'L1':       'black',        # baseline reference
    'L2a':      '#eb6834',      # orange
    'L2b':      '#4a3aa7',      # violet   (was purple)
    'L3':       '#1baf7a',      # aqua     (was cyan -- the worst offender)
    'L4_gb':    '#008300',      # green
    'L4_trap':  '#eda100',      # yellow   (was orange, which clashed with green)
    'L4_both':  '#2a78d6',      # blue     (was crimson, which clashed with L5 red)
    'L5':       '#e34948',      # red
}

# =============================================================================
# CURVE STYLE SYSTEM
# Each level/mode gets a unique combination of color + linestyle + marker
# Distinguishable by color alone, linestyle alone, or marker alone
# =============================================================================
CURVE_STYLES = {
    # Model levels
    'L1':       {'color': 'black',   'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L1_L6':    {'color': 'black',   'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L2a':      {'color': '#eb6834', 'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},  # blue
    'L2a_L6':   {'color': '#eb6834', 'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L2b':      {'color': '#4a3aa7',  'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L2_L6':    {'color': '#4a3aa7',  'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L3':       {'color': '#1baf7a',    'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L3_L6':    {'color': '#1baf7a',    'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L4_gb':    {'color': '#008300',   'ls': '-',    'marker': '^',  'lw': 2.5, 'ms': 6},
    'L4_trap':  {'color': '#eda100',  'ls': '-',    'marker': 'v',  'lw': 2.5, 'ms': 6},
    'L4_both':  {'color': '#2a78d6', 'ls': '-',    'marker': 'D',  'lw': 2.5, 'ms': 6},
    'L5':       {'color': '#e34948',     'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L3_L4_L6': {'color': '#e34948',     'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},

    # Defect fractions (used in L3 sweep plots)
    'f_0':      {'color': 'black',   'ls': '-',    'marker': None, 'lw': 2.5, 'ms': 0},
    'f_001':    {'color': '#2a78d6', 'ls': '--',   'marker': None, 'lw': 2.0, 'ms': 0},  # 0.1%
    'f_01':     {'color': '#eda100', 'ls': '-.',   'marker': None, 'lw': 2.0, 'ms': 0},  # 1%
    'f_10':     {'color': '#e34948', 'ls': ':',    'marker': None, 'lw': 2.0, 'ms': 0},  # 10%

    # Defect types
    'intact':   {'color': 'black',   'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 5},
    'pinhole':  {'color': '#e34948', 'ls': '--',   'marker': '^',  'lw': 2.0, 'ms': 5},  # red
    'crack':    {'color': '#eda100', 'ls': '-.',   'marker': 's',  'lw': 2.0, 'ms': 5},  # orange
    'gb':       {'color': '#008300', 'ls': ':',    'marker': 'D',  'lw': 2.0, 'ms': 5},  # green

    # Rate-limiting regions (shading)
    'surface_region': {'color': '#e34948', 'alpha': 0.9},  # red
    'oxide_region':   {'color': '#eda100', 'alpha': 0.9},  # orange
    'metal_region':   {'color': '#2a78d6', 'alpha': 0.9},  # blue
    'mixed_region':   {'color': '#1baf7a', 'alpha': 0.9},  # green

    # Reference lines (slopes, theory)
    'slope_1':  {'color': '#e34948', 'ls': '--',   'lw': 1.5, 'alpha': 0.6},
    'slope_05': {'color': '#008300', 'ls': '--',   'lw': 1.5, 'alpha': 0.6},
    'parity':   {'color': '#e34948',     'ls': '--',   'lw': 2.0, 'alpha': 0.8},
    'threshold':{'color': 'gray',    'ls': '--',   'lw': 1.5, 'alpha': 0.5},
}

# =============================================================================
# PLOT HELPER — apply a curve style in one call
# =============================================================================
def apply_style(ax, x, y, style_key, label=None, log=False):
    """
    Plot a curve using a named style from CURVE_STYLES.

    Parameters
    ----------
    ax        : matplotlib Axes
    x, y      : array-like
    style_key : str — key in CURVE_STYLES
    label     : str — legend label (optional)
    log       : bool — if True use loglog, else plot

    Returns
    -------
    matplotlib Line2D
    """
    s = CURVE_STYLES[style_key]
    kwargs = dict(
        color=s['color'],
        linestyle=s['ls'],
        linewidth=s['lw'],
        alpha=s.get('alpha', 1.0),
    )
    if s.get('marker') and s.get('ms', 0) > 0:
        kwargs['marker']     = s['marker']
        kwargs['markersize'] = s['ms']
        kwargs['markevery']  = 5   # don't crowd markers
    if label:
        kwargs['label'] = label

    if log:
        return ax.loglog(x, y, **kwargs)[0]
    else:
        return ax.plot(x, y, **kwargs)[0]

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
