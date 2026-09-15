"""Active study configuration.

This module is a SWITCH, not a config. It re-exports one module from
`calculations/config/studies/`, and everything else in the codebase imports from
here — so selecting a study is the single line below and no other file changes.

    calculations/*.py  ->  model_config  ->  studies/<the active one>.py

To switch study: change ACTIVE_STUDY below, then restart any running kernel (the
binding happens at import time, so `importlib.reload` on a notebook module is not
enough — the study's dicts are already bound into calculations.sensitivity).

To add a study: copy a module in studies/ , edit it, set ACTIVE_STUDY to its
module name. See studies/__init__.py for the full recipe and for why the studies
live there rather than beside the notebooks.

ACTIVE_STUDY is the ONLY switch. It used to be a literal `from studies.X import *`
with ACTIVE_STUDY as a separate label describing it, which let the two disagree —
and they did: the label read one study while the import loaded another, and at one
point BOTH studies were star-imported at once, so Python applied them in order and
produced a hybrid. That state had Hastelloy N's material dicts with Incoloy's
derived sensitivity parameters, and it moved L1 flux by a factor of ten while
`ACTIVE_STUDY` still claimed Incoloy. Exactly one study is imported here now, named
by exactly one variable, so neither failure can recur.

After switching, run:

    from calculations.sensitivity import check_against_config
    check_against_config()

Preset yields, sweep ranges and draw counts were tuned for a specific material and
are not automatically valid for another one.
"""

import importlib as _importlib
import pkgutil as _pkgutil

from calculations.config import studies as _studies_pkg

# =============================================================================
# ACTIVE STUDY — change this one line to switch
# =============================================================================
# ACTIVE_STUDY = 'incoloy802_cr2o3'
# ACTIVE_STUDY = 'fuerst_etal_2024_model_config'
ACTIVE_STUDY = 'Guo_etal_2025_316L'


# =============================================================================


def _available_studies():
    """Module names in calculations/config/studies/, excluding the package init."""
    return sorted(m.name for m in _pkgutil.iter_modules(_studies_pkg.__path__)
                  if not m.name.startswith('_'))


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
    'VALIDATION', 'build_simulation_config',
    # presentation
    'PLOT_STYLE', 'COLORS', 'CURVE_STYLES',
})


def _load_active_study(name):
    """Import one study module and copy its public names into this namespace.

    `from studies.<name> import *` cannot be written with a variable — a star import
    is a static statement whose module path is resolved at compile time. So the
    re-export is done by hand, reproducing star-import semantics exactly: every
    non-underscore name, honouring __all__ if the study defines one.

    Note that includes `np`, which the studies import for their VALIDATION arrays.
    The old star import leaked it into this namespace too; it is kept so switching
    to this mechanism changes nothing observable.
    """
    try:
        mod = _importlib.import_module(f'{__package__}.studies.{name}')
    except ModuleNotFoundError as exc:
        # Only mask a missing STUDY; a broken import inside a study that does exist
        # must surface as itself, or a typo'd dependency looks like a typo'd study.
        if getattr(exc, 'name', None) != f'{__package__}.studies.{name}':
            raise
        raise ValueError(
            f"ACTIVE_STUDY = {name!r} does not name a module in "
            f"calculations/config/studies/. Available: {_available_studies()}. "
            f"Use the module name without the '.py' extension."
        ) from None

    exported = getattr(mod, '__all__', None)
    if exported is None:
        exported = [n for n in dir(mod) if not n.startswith('_')]

    missing = REQUIRED_EXPORTS - set(exported)
    if missing:
        raise ValueError(
            f"study {name!r} is missing {len(missing)} required export(s): "
            f"{sorted(missing)}. Every study must define these — copy them from an "
            f"existing module in calculations/config/studies/ and adapt the values "
            f"to the material. Derived dicts (DEFAULT_PARAMS_*) should be copied as "
            f"CODE so they recompute from that study's own METALS/OXIDES; only the "
            f"study-design literals (ranges, presets) are written out by hand."
        )

    return mod, {n: getattr(mod, n) for n in exported}


_study_module, _study_ns = _load_active_study(ACTIVE_STUDY)
globals().update(_study_ns)

# Provenance for saved results and for anything labelling a figure or a CSV. This
# is now a consequence of the import above, not a parallel claim about it.
ACTIVE_STUDY_MODULE = _study_module.__name__
