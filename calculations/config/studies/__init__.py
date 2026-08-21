"""Per-study configurations.

Each module in this package is one complete, self-contained study configuration:
material data, operating conditions, microstructure, defects, and the sensitivity
analysis parameters derived from them.

**Exactly one study is active at a time.** `calculations/config/model_config.py` is a
three-line switch that re-exports whichever one you name there, and every module in
`calculations/` imports from `model_config` — so switching studies is a one-line edit
and nothing else in the codebase changes.

To add a study
--------------
1. Copy an existing module here (e.g. `incoloy802_cr2o3.py`) to a new name.
2. Edit its METALS / OXIDES / CONDITIONS / MICROSTRUCTURE / OXIDE_DEFECTS, and the
   SENSITIVITY ANALYSIS PARAMETERS section if the ranges or presets should differ.
3. Point the import in `model_config.py` at it.
4. Run `calculations.sensitivity.check_against_config()` before trusting results —
   the preset yields and sweep ranges were tuned for a specific material and will
   not automatically be right for another one.

Why not a study-local model_config.py per folder
------------------------------------------------
Because `from calculations.config.model_config import ...` is an ABSOLUTE package
import: it resolves through sys.path, not the current directory. A `model_config.py`
sitting next to a notebook would be a *different module* from the one every file in
`calculations/` reads — the notebook would report one material while the model
computed another. Keeping the studies here, behind one switch, is what makes that
impossible.
"""
