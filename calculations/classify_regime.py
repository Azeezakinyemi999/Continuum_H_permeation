"""
calculations/classify_regime.py

Hierarchical regime classification utilities used by Level 1-4 models.

This module centralizes the various `classify_regime_*` helper functions so
that other calculation modules can import and reuse a single implementation.

What "regime" means here
------------------------
Hydrogen crossing an oxide-coated wall meets several resistances in series and
in parallel: dissociation at the gas-facing surface, molecular diffusion through
the intact oxide, short-circuit transport through oxide defects, and atomic
diffusion through the metal where traps compete with grain-boundary fast paths.
Only one of these usually sets the flux. Naming that one is what a "regime"
label does, and it is the first thing to look at when a flux moves unexpectedly:
a change in a parameter the regime does not depend on should do nothing.

These functions compute NO physics. Every number they inspect was already
computed upstream; they only interpret it and assemble the labels. Keeping that
interpretation in one place is why the module exists — the thresholds below are
the project's definition of a regime boundary, and the plot legends, regime
clusters and sensitivity tables all key off the strings produced here.

The label vocabulary
--------------------
Tier 1, ``base_regime``, comes from the oxide/metal resistance ratio in
:func:`calculations.interface_solver` (``R_oxide/R_metal`` > 10 gives
``'oxide_limited'``, < 0.5 gives ``'metal_limited'``, otherwise
``'transition'``); callers pass ``'unknown'`` when no solve was done.

Tier 2, ``regime_detail``, subdivides tier 1 and is the reason this module is
hierarchical rather than flat:

- under ``'oxide_limited'`` — ``'defect_limited'`` if the defect paths carry more
  flux than the intact oxide, else ``'regime_intact_oxide'``
- under ``'metal_limited'`` — ``'traps_defect_limited'`` if trapping has cut the
  effective diffusivity by more than the threshold, else ``'lattice_limited'``

``regime_hierarchy`` joins the filled tiers with ``/`` (e.g.
``'metal_limited/traps_defect_limited'``) and is the string to use for grouping
and labelling. ``classification_depth`` says how many tiers were resolved.
``regime_subdetail`` is always ``None`` — a third tier that was never defined.

Which of these are live
-----------------------
:func:`classify_regime_level14` (from
:func:`calculations.permeation_calc.calculate_defective_metal_flux`) and
:func:`classify_regime_level34` (from
:mod:`calculations.parallel_oxide_defect_paths`) are the two the model calls.
:func:`classify_regime_level4_metal` is their shared tier-2 helper.

:func:`classify_regime_level2`, :func:`classify_regime_level3` and
:func:`classify_regime_level24` have no callers in the current source tree. They
were reachable from the archived ``plot&docs/analysis.ipynb``, and the 1.0.1
release imported all six names into ``parallel_oxide_defect_paths`` while
calling only ``level34`` — those unused imports were trimmed for 1.1.0. They are
kept because they complete the Level 1-4 matrix and their return shapes are
published in ``OUTPUTS.md``; treat them as API, not as dead code to delete
silently.

Note that Level 5 and Level 5+L6 do **not** use this module. They assign a
single flat label from resistance fractions via ``assign_regime_L5`` /
``assign_regime`` in :mod:`calculations.sensitivity` — ``oxide`` / ``metal`` /
``defect`` without surface kinetics, ``surface`` / ``oxide`` / ``metal`` with
them. So the hierarchical ``a/b`` strings built here and the flat L5 labels are
two separate vocabularies; do not mix them when grouping results.
"""

def classify_regime_level2(base_regime):
    """Wrap a Level 2 regime label in the standard classification dict.

    Level 2 is a perfect oxide on a perfect metal: two resistances in series and
    nothing inside either one to subdivide. So this is a pass-through that
    records the tier-1 label and reports depth 1, existing only so Level 2
    results carry the same dict shape as the deeper levels.

    Parameters
    ----------
    base_regime : str
        Tier-1 label from the interface solve: 'oxide_limited', 'metal_limited',
        'transition' or 'unknown'.

    Returns
    -------
    dict
        Keys ``model_level``, ``base_regime``, ``regime_hierarchy``,
        ``regime_detail`` (None), ``regime_subdetail`` (None) and
        ``classification_depth`` (1). ``regime_hierarchy`` equals
        ``base_regime``, with no ``/`` since there is no second tier.

    Notes
    -----
    No caller in the current source tree — see the module docstring.
    """
    return {
        'model_level': 'Level 2',
        'base_regime': base_regime,
        'regime_hierarchy': base_regime,
        'regime_detail': None,
        'regime_subdetail': None,
        'classification_depth': 1
    }


def classify_regime_level3(base_regime, flux_intact_contribution, flux_defect_contribution):
    """Classify a Level 3 result, subdividing an oxide-limited wall by defect path.

    Level 3 adds pinholes, cracks and oxide grain boundaries, which carry
    hydrogen around the intact oxide rather than through it. When the oxide is
    what limits the flux, the question becomes which oxide path dominates: the
    defects are a short circuit whose area fraction is tiny but whose local
    permeability is large, so either side can win depending on defect density.
    Whichever contribution is larger names the tier-2 regime.

    The comparison is on total contributions (area-weighted), not on
    per-area flux, so a defect population can dominate transport while still
    occupying a negligible fraction of the wall.

    If the wall is not oxide-limited, the defects cannot be controlling and no
    second tier is assigned.

    Parameters
    ----------
    base_regime : str
        Tier-1 label: 'oxide_limited', 'metal_limited', 'transition' or
        'unknown'. Only 'oxide_limited' triggers subdivision.
    flux_intact_contribution : float
        Flux carried by the intact oxide, area-weighted [mol/m²/s].
    flux_defect_contribution : float
        Flux carried by all defect paths combined, area-weighted [mol/m²/s].

    Returns
    -------
    dict
        The standard classification keys plus
        ``flux_ratio_defect_to_intact``, the defect/intact contribution ratio.
        Under 'oxide_limited', ``regime_detail`` is 'defect_limited' or
        'regime_intact_oxide' and depth is 2; otherwise ``regime_detail`` is
        None and depth is 1.

    Notes
    -----
    ``flux_ratio_defect_to_intact`` is ``inf`` when the intact contribution is
    zero — a fully short-circuited oxide, not an error. It is reported
    regardless of ``base_regime``, so it stays readable even when no second tier
    was assigned.

    No caller in the current source tree — see the module docstring. The live
    Level 3 path is :func:`classify_regime_level34`, which applies this same
    test.
    """
    if base_regime == 'oxide_limited':
        if flux_defect_contribution > flux_intact_contribution:
            regime_detail = 'defect_limited'
        else:
            regime_detail = 'regime_intact_oxide'
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
    else:
        regime_detail = None
        regime_hierarchy = base_regime
        classification_depth = 1

    return {
        'model_level': 'Level 3',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': classification_depth,
        'flux_ratio_defect_to_intact': flux_defect_contribution / flux_intact_contribution if flux_intact_contribution > 0 else float('inf')
    }


def classify_regime_level4_metal(modification_factor, threshold_traps=0.5):
    """Decide whether traps or the lattice control diffusion through the metal.

    Two microstructural effects pull in opposite directions inside the metal.
    Traps — dislocations, grain boundaries, vacancies, carbides — hold hydrogen
    at binding energies of tens of kJ/mol and slow its net advance, while grain
    boundaries also offer fast paths that speed it up. Their combined outcome is
    the ratio of effective to lattice diffusivity, and the sign of the deviation
    from 1 says which effect won.

    This is the shared tier-2 helper for every level combination that includes
    Level 4. It returns the detail fields only, not a full classification dict.

    Parameters
    ----------
    modification_factor : float
        ``D_eff / D_lattice`` [-], as produced by
        :mod:`calculations.defective_metal` and carried through the interface
        solve. Below 1 means trapping dominates, above 1 means grain-boundary
        enhancement dominates, 1 means they cancel or neither is active.
    threshold_traps : float, optional
        Boundary below which trapping is called controlling. Default 0.5, i.e.
        trapping must have at least halved the effective diffusivity. This is a
        project convention, not a physical constant — it defines the regime
        boundary rather than describing one.

    Returns
    -------
    dict
        ``metal_regime_detail`` ('traps_defect_limited' or 'lattice_limited'),
        ``modification_factor`` echoed back, ``trapping_significant`` (bool) and
        ``trapping_reduction_percent``.

    Notes
    -----
    ``trapping_reduction_percent`` is ``(1 - modification_factor) * 100``, so it
    goes **negative** when grain-boundary enhancement outweighs trapping and
    ``modification_factor`` exceeds 1. Read it as a signed change in
    diffusivity, not as a reduction that is always positive.

    Not called directly outside this module; reached via
    :func:`classify_regime_level14`, :func:`classify_regime_level24` and
    :func:`classify_regime_level34`.
    """
    if modification_factor < threshold_traps:
        metal_regime_detail = 'traps_defect_limited'
        trapping_significant = True
    else:
        metal_regime_detail = 'lattice_limited'
        trapping_significant = False

    return {
        'metal_regime_detail': metal_regime_detail,
        'modification_factor': modification_factor,
        'trapping_significant': trapping_significant,
        'trapping_reduction_percent': (1 - modification_factor) * 100
    }


def classify_regime_level14(modification_factor, threshold_traps=0.5):
    """Classify a bare defective metal (Level 1+4) — always metal-limited.

    With no oxide on the wall there is no other resistance to compete with, so
    tier 1 is fixed at 'metal_limited' and is not inferred from anything. The
    only open question is the one inside the metal: traps or lattice. That makes
    this the simplest of the classifiers and the only one whose tier-1 label is
    a constant.

    Parameters
    ----------
    modification_factor : float
        ``D_eff / D_lattice`` [-]. See :func:`classify_regime_level4_metal`.
    threshold_traps : float, optional
        Trapping-significance boundary, default 0.5.

    Returns
    -------
    dict
        The standard classification keys with ``base_regime`` always
        'metal_limited' and ``classification_depth`` always 2, plus
        ``modification_factor``, ``trapping_significant`` and
        ``trapping_reduction_percent``.

    Examples
    --------
    >>> classify_regime_level14(0.2)['regime_hierarchy']
    'metal_limited/traps_defect_limited'
    >>> classify_regime_level14(0.9)['regime_hierarchy']
    'metal_limited/lattice_limited'

    Notes
    -----
    Called by
    :func:`calculations.permeation_calc.calculate_defective_metal_flux`. Both
    labels above are pinned in ``container/model_smoke_test.py``, so changing
    the threshold or the label strings will fail the release gate — which is the
    intent, since downstream plots and regime clusters match on these strings.
    """
    metal_detail = classify_regime_level4_metal(modification_factor, threshold_traps)

    base_regime = 'metal_limited'
    regime_detail = metal_detail['metal_regime_detail']
    regime_hierarchy = f"{base_regime}/{regime_detail}"

    return {
        'model_level': 'Level 1,4',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': 2,
        'modification_factor': modification_factor,
        'trapping_significant': metal_detail['trapping_significant'],
        'trapping_reduction_percent': metal_detail['trapping_reduction_percent']
    }


def classify_regime_level24(base_regime, modification_factor, threshold_traps=0.5):
    """Classify a perfect oxide over a defective metal (Level 2+4).

    Here the oxide is intact but the metal has traps and grain boundaries, so
    unlike Level 1+4 the wall may well be oxide-limited. Subdividing by
    trapping is only meaningful when the metal is the bottleneck: if the oxide
    controls, hydrogen barely reaches the traps in quantity and the trap
    population cannot set the flux no matter how dense it is.

    Parameters
    ----------
    base_regime : str
        Tier-1 label: 'oxide_limited', 'metal_limited', 'transition' or
        'unknown'. Only 'metal_limited' triggers subdivision.
    modification_factor : float
        ``D_eff / D_lattice`` [-]. See :func:`classify_regime_level4_metal`.
    threshold_traps : float, optional
        Trapping-significance boundary, default 0.5.

    Returns
    -------
    dict
        The standard classification keys plus ``modification_factor``,
        ``trapping_significant`` and ``trapping_reduction_percent``.

    Notes
    -----
    When ``base_regime`` is not 'metal_limited', ``trapping_significant`` is
    forced to False and ``trapping_reduction_percent`` to 0.0 even if the traps
    really did change the diffusivity. Those two fields describe whether
    trapping is *controlling the flux*, not whether it is present, and
    ``modification_factor`` is still echoed back unmodified if you need the
    underlying number.

    No caller in the current source tree — see the module docstring.
    """
    if base_regime == 'metal_limited':
        metal_detail = classify_regime_level4_metal(modification_factor, threshold_traps)
        regime_detail = metal_detail['metal_regime_detail']
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
        trapping_significant = metal_detail['trapping_significant']
        trapping_reduction = metal_detail['trapping_reduction_percent']
    else:
        regime_detail = None
        regime_hierarchy = base_regime
        classification_depth = 1
        trapping_significant = False
        trapping_reduction = 0.0

    return {
        'model_level': 'Level 2,4',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': classification_depth,
        'modification_factor': modification_factor,
        'trapping_significant': trapping_significant,
        'trapping_reduction_percent': trapping_reduction
    }


def classify_regime_level34(base_regime, flux_intact_contribution, flux_defect_contribution,
                            modification_factor, threshold_traps=0.5):
    """Classify a defective oxide over a defective metal (Level 3+4).

    This is the full picture below Level 5: defects short-circuit the oxide and
    traps impede the metal at the same time. Both mechanisms are available, but
    they sit on opposite sides of the wall, so only the one on the controlling
    side can be the explanation. Tier 1 therefore selects which tier-2 test is
    even asked — the defect/intact split if the oxide controls, the trap/lattice
    split if the metal does, and neither in the transition band where the two
    resistances are comparable.

    Parameters
    ----------
    base_regime : str
        Tier-1 label: 'oxide_limited', 'metal_limited', 'transition' or
        'unknown'.
    flux_intact_contribution, flux_defect_contribution : float
        Area-weighted flux through the intact oxide and through all defect
        paths [mol/m²/s]. Used only when ``base_regime`` is 'oxide_limited'.
    modification_factor : float
        ``D_eff / D_lattice`` [-]. Used only when ``base_regime`` is
        'metal_limited'.
    threshold_traps : float, optional
        Trapping-significance boundary, default 0.5.

    Returns
    -------
    dict
        The standard classification keys plus ``modification_factor``,
        ``trapping_significant``, ``trapping_reduction_percent`` and
        ``flux_ratio_defect_to_intact``.

    Notes
    -----
    Only the fields belonging to the branch that ran are populated, so the
    unused diagnostic is blanked rather than merely stale:

    - oxide-limited: ``flux_ratio_defect_to_intact`` is set (``inf`` if the
      intact contribution is zero); ``trapping_significant`` is forced False and
      ``trapping_reduction_percent`` to 0.0.
    - metal-limited: the trapping fields are set;
      ``flux_ratio_defect_to_intact`` is None.
    - transition or unknown: depth 1, no ``regime_detail``, trapping fields
      blanked and ``flux_ratio_defect_to_intact`` None.

    ``modification_factor`` is echoed back in every branch, so a suppressed
    ``trapping_significant`` can always be re-derived. Unlike
    :func:`classify_regime_level3`, the flux ratio is **not** computed outside
    the oxide branch — do not read a None there as "no defects".

    Called twice in :mod:`calculations.parallel_oxide_defect_paths`, for the
    perfect-metal and defective-metal parallel-path assemblies.
    """
    if base_regime == 'oxide_limited':
        if flux_defect_contribution > flux_intact_contribution:
            regime_detail = 'defect_limited'
        else:
            regime_detail = 'regime_intact_oxide'
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
        trapping_significant = False
        trapping_reduction = 0.0
        flux_ratio = flux_defect_contribution / flux_intact_contribution if flux_intact_contribution > 0 else float('inf')

    elif base_regime == 'metal_limited':
        metal_detail = classify_regime_level4_metal(modification_factor, threshold_traps)
        regime_detail = metal_detail['metal_regime_detail']
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
        trapping_significant = metal_detail['trapping_significant']
        trapping_reduction = metal_detail['trapping_reduction_percent']
        flux_ratio = None

    else:
        regime_detail = None
        regime_hierarchy = base_regime
        classification_depth = 1
        trapping_significant = False
        trapping_reduction = 0.0
        flux_ratio = None

    return {
        'model_level': 'Level 3,4',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': classification_depth,
        'modification_factor': modification_factor,
        'trapping_significant': trapping_significant,
        'trapping_reduction_percent': trapping_reduction,
        'flux_ratio_defect_to_intact': flux_ratio
    }
