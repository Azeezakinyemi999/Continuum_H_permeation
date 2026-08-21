"""Active study configuration.

This module is a SWITCH, not a config. It re-exports one module from
`calculations/config/studies/`, and everything else in the codebase imports from
here — so selecting a study is the single line below and no other file changes.

    calculations/*.py  ->  model_config  ->  studies/<the active one>.py

To switch study: change the import below, then restart any running kernel (the
binding happens at import time, so `importlib.reload` on a notebook module is not
enough — the study's dicts are already bound into calculations.sensitivity).

To add a study: copy a module in studies/ , edit it, point this line at it. See
studies/__init__.py for the full recipe and for why the studies live there rather
than beside the notebooks.

After switching, run:

    from calculations.sensitivity import check_against_config
    check_against_config()

Preset yields, sweep ranges and draw counts were tuned for a specific material and
are not automatically valid for another one.
"""

# =============================================================================
# ACTIVE STUDY — change this one line to switch
# =============================================================================
from calculations.config.studies.incoloy802_cr2o3 import *   # noqa: F401,F403

# Record which study is active, for provenance in saved results and for anything
# that wants to label a figure or a CSV with its source.
ACTIVE_STUDY = 'incoloy802_cr2o3'
