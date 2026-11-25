"""
calculations/classify_regime.py

Hierarchical regime classification utilities used by the Level 2-4 models.

What "regime" means here
------------------------
Hydrogen crossing an oxide-coated wall meets several resistances in series and
in parallel: dissociation at the gas-facing surface, molecular diffusion through
the intact oxide, short-circuit transport through oxide defects, and atomic
diffusion through the metal. Only one of these usually sets the flux. Naming
that one is what a "regime" label does, and it is the first thing to look at
when a flux moves unexpectedly: a change in a parameter the regime does not
depend on should do nothing.

These functions compute NO physics. Every number they inspect was already
computed upstream; they only interpret it and assemble the labels. The plot
legends, regime clusters and sensitivity tables all key off the strings
produced here.

The label vocabulary
--------------------
Tier 1, ``base_regime``, comes from the oxide/metal resistance ratio in
:mod:`calculations.interface_solver` (``R_oxide/R_metal`` > 10 gives
``'oxide_limited'``, < 0.5 gives ``'metal_limited'``, otherwise
``'transition'``); callers pass ``'unknown'`` when no solve was done.

Tier 2, ``regime_detail``, subdivides an ``'oxide_limited'`` wall:
``'defect_limited'`` if the defect paths carry more flux than the intact oxide,
else ``'regime_intact_oxide'``. A metal-limited wall has no second tier. (It used
to be split into ``'traps_defect_limited'`` / ``'lattice_limited'`` by a
trap-reduced diffusivity; trapping does not change steady-state flux and was
removed, so that split no longer exists.)

``regime_hierarchy`` joins the filled tiers with ``/`` (e.g.
``'oxide_limited/defect_limited'``) and is the string to use for grouping and
labelling. ``classification_depth`` says how many tiers were resolved.

:func:`classify_regime_level3` is called by
:mod:`calculations.parallel_oxide_defect_paths` for every Level 3 and Level 5
parallel-path assembly; :func:`classify_regime_level2` completes the matrix.

Note that Level 5 and Level 5+L6 regime *clusters* do **not** use this module.
They assign a single flat label from resistance fractions via
``assign_regime_L5`` / ``assign_regime`` in :mod:`calculations.sensitivity` --
``oxide`` / ``metal`` / ``defect`` without surface kinetics, ``surface`` /
``oxide`` / ``metal`` with them. Do not mix the two vocabularies when grouping.
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

    Called by every parallel-path assembly in
    :mod:`calculations.parallel_oxide_defect_paths`.
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
