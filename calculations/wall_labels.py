# calculations/wall_labels.py
"""
Regime labels that say which path a wall's hydrogen takes and which step limits it.

A wall with a defective oxide is several paths side by side: the intact oxide, pinholes,
cracks and oxide grain boundaries. They exchange no hydrogen, so their fluxes add, and two
questions describe a point:

* **Which path dominates?** The one carrying more than half the wall's flux
  (:data:`DOMINANCE`). Parallel paths add, so the dominant path is the fastest route, not the
  slowest. If no path carries more than half, the point is ``mixed``.
* **Which step limits it?** The slowest of that path's steps in series, as the solver itself
  classifies that path: the gas-facing surface (Model 2 only), oxide diffusion or metal
  diffusion, or ``mixed`` when no step holds the majority.

A surface step is named for the surface the gas meets on that path: the oxide on the intact
area, under cracks and over oxide grain boundaries; the bare metal under pinholes and on a wall
with no oxide.

Labels read ``<path>/<step>[/<metal detail>]`` for a wall with defects, e.g.
``pinholes/metal surface`` or ``intact/metal diffusion/traps``, and ``<step>[/<metal detail>]``
for a wall that is a single path, e.g. ``oxide surface``. The metal detail, ``traps`` or
``lattice``, appears only on the levels with a defective metal (L4, L5 and their +L6 forms).

Model 1 has no surface step, because its surface is at equilibrium. Its oxide/metal call is the
interface solver's own: oxide when R_oxide/R_metal > 10, metal when < 0.5, otherwise mixed.
Model 2's is the solver's ``rate_limiting``: the step holding more than half the path's
resistance, otherwise mixed. Nothing here re-derives either rule; this module only names what
the solvers return and decides which path dominates.
"""

import numpy as np

from calculations.classify_regime import classify_regime_level4_metal
from calculations.permeation_calc import calculate_simple_metal_flux, calculate_defective_metal_flux
# calculations.interface_solver is imported inside model1_paths(): it loads the active study's config,
# and describe()/fitted() must import without one -- the study comparison runs outside every study.

#: A path dominates when it carries more than this share of the wall's flux.
DOMINANCE = 0.5

#: Solver path key -> the path's name in a label.
PATH_NAMES = {'intact': 'intact', 'pinhole': 'pinholes', 'crack': 'cracks', 'grain_boundary': 'oxide GBs'}

#: The surface the gas meets on each path of a wall with defects.
SURFACE_OF = {'intact': 'oxide surface', 'pinhole': 'metal surface',
              'crack': 'oxide surface', 'grain_boundary': 'oxide surface'}

#: The surface of each Model 2 level that is a single path.
SINGLE_PATH_SURFACE = {'L1L6': 'metal surface', 'L2aL6': 'oxide surface',
                       'L2bL6': 'oxide surface', 'L4L6': 'oxide surface'}

#: Levels whose metal is defective, so a metal-diffusion step says whether traps control it.
DEFECTIVE_METAL = {'L4', 'L5', 'L4L6', 'L5L6'}

_COMPONENT_KEY = {'pinhole': 'pinholes', 'crack': 'cracks', 'grain_boundary': 'grain_boundaries'}
_MODEL2_STEP   = {'oxide': 'oxide diffusion', 'metal': 'metal diffusion', 'mixed': 'mixed'}


# =============================================================================
# Building a label
# =============================================================================

def compose(step, detail='', path=None):
    """``'pinholes/metal surface'`` from its parts; ``path`` is a solver key or None for one path."""
    return '/'.join(([PATH_NAMES[path]] if path else []) + [step] + ([detail] if detail else []))


def parse(label):
    """``(path name or None, step, metal detail)`` from a label: the inverse of :func:`compose`."""
    parts = label.split('/')
    path = parts.pop(0) if parts[0] in PATH_NAMES.values() else None
    return path, parts[0], (parts[1] if len(parts) > 1 else '')


def trap_detail(modification_factor):
    """``'traps'`` when trapping slows the metal (D_eff/D_lattice below the Level 4 cut), else ``'lattice'``."""
    return 'traps' if classify_regime_level4_metal(modification_factor)['trapping_significant'] else 'lattice'


def model1_step(regime):
    """A Model 1 regime string as ``(step, metal detail)``.

    ``'oxide_limited'`` -> ``('oxide diffusion', '')``; ``'metal_limited/traps_defect_limited'``
    -> ``('metal diffusion', 'traps')``; ``'transition'`` -> ``('mixed', '')``.
    """
    base, *rest = regime.split('/')
    step = {'oxide_limited': 'oxide diffusion', 'metal_limited': 'metal diffusion'}.get(base, 'mixed')
    if step != 'metal diffusion' or not rest:
        return step, ''
    return step, ('traps' if rest[0].startswith('traps') else 'lattice')


def label_from_paths(paths, J_total):
    """Label of a wall with defects: the dominant path, then the limiting step inside it.

    Parameters
    ----------
    paths : dict
        ``{solver path key: {'contribution': ..., 'step': ..., 'detail': ...}}``, as returned by
        :func:`model1_paths` or :func:`model2_paths`.
    J_total : float
        The wall's flux, which the contributions add up to.
    """
    if not (np.isfinite(J_total) and J_total > 0):
        return 'undefined'
    share = {p: d['contribution'] / J_total for p, d in paths.items()}
    top = max(share, key=share.get)
    if not share[top] > DOMINANCE:
        return 'mixed'
    return compose(paths[top]['step'], paths[top]['detail'], top)


# =============================================================================
# Model 1 (L1-L5)
# =============================================================================

def defect_fractions(defects):
    """Area fraction of each defect type, read the way the Level 3 solvers read an OXIDE_DEFECTS dict.

    A single type carries the whole ``area_fraction``. ``'mixed'`` splits it in proportion to
    ``components``; with no components the solvers treat the whole area as pinholes.
    """
    f = defects.get('area_fraction', 0.01)
    kind = defects.get('type', 'pinhole')
    if kind != 'mixed':
        return {kind: f}
    components = defects.get('components', {})
    total = sum(components.values())
    if total == 0:
        return {'pinhole': f}
    return {t: f * components[key] / total for t, key in _COMPONENT_KEY.items() if components.get(key, 0) > 0}


def model1_paths(level, P_up, P_down, oxide_props, metal_props, defects, temperature=None, micro=None):
    """Every path of a Model 1 wall with defects (L3 or L5), each solved on its own.

    The parallel-path solvers return the wall's flux and a label for the intact path only, so
    each path is re-solved here with the solver those modules use for it, with the same defaults:
    the intact, crack and grain-boundary paths with the oxide + metal system (thin oxide for a
    crack, faster oxide diffusion for a grain boundary), and a pinhole with the metal alone.

    Parameters
    ----------
    level : {'L3', 'L5'}
    P_up, P_down : float
        Upstream and downstream pressure [Pa].
    oxide_props, metal_props : dict
        The layer properties the wall was solved with, ``thickness`` included.
    defects : dict
        The OXIDE_DEFECTS-format dict the wall was solved with.
    temperature : float, optional
        [K]. L5 only.
    micro : dict, optional
        The keyword arguments the L5 solver was given: ``microstructure_params``,
        ``lattice_density`` and ``mode``. L5 only.

    Returns
    -------
    dict
        ``{path key: {'contribution', 'flux', 'step', 'detail', 'fraction_surface',
        'fraction_oxide', 'fraction_metal'}}``. Contributions are area-weighted, so they add up
        to the wall's flux.
    """
    from calculations.interface_solver import (calculate_oxide_metal_system,
                                               calculate_oxide_defective_metal_system)
    defective = level == 'L5'

    def two_layer(ox):
        if defective:
            return calculate_oxide_defective_metal_system(P_up, P_down, ox, metal_props, temperature, **micro)
        return calculate_oxide_metal_system(P_up, P_down, ox, metal_props)

    def described(r, area):
        step, _ = model1_step(r['regime'])
        detail = trap_detail(r['modification_factor']) if (defective and step == 'metal diffusion') else ''
        R_ox, R_m = r['R_oxide'], r['R_metal']
        return {'contribution': area * r['flux'], 'flux': r['flux'], 'step': step, 'detail': detail,
                'fraction_surface': 0.0,                          # the surface is at equilibrium
                'fraction_oxide': R_ox / (R_ox + R_m), 'fraction_metal': R_m / (R_ox + R_m)}

    fractions = defect_fractions(defects)
    paths = {'intact': described(two_layer(oxide_props), 1.0 - sum(fractions.values()))}
    for kind, area in fractions.items():
        if kind == 'pinhole':                                     # no oxide on this path: the metal alone
            if defective:
                r = calculate_defective_metal_flux(D_lattice=metal_props['D_metal'], K_s=metal_props['K_s_metal'],
                                                   metal_thickness=metal_props['thickness'], P_up=P_up,
                                                   P_down=P_down, temperature=temperature, **micro)
                step, detail = model1_step(r['regime'])
            else:
                r = calculate_simple_metal_flux(metal_props['D_metal'], metal_props['K_s_metal'],
                                                metal_props['thickness'], P_up, P_down)
                step, detail = 'metal diffusion', ''
            paths[kind] = {'contribution': area * r['flux'], 'flux': r['flux'], 'step': step, 'detail': detail,
                           'fraction_surface': 0.0, 'fraction_oxide': 0.0, 'fraction_metal': 1.0}
            continue
        ox = dict(oxide_props)
        if kind == 'crack':
            ox['thickness'] = oxide_props['thickness'] * defects.get('thickness_factor', 0.1)
        else:
            ox['D_ox'] = oxide_props['D_ox'] * defects.get('diffusivity_factor', 10)
        paths[kind] = described(two_layer(ox), area)
    return paths


def label_model1(level, result=None, P_up=None, P_down=None, oxide_props=None, metal_props=None,
                 defects=None, temperature=None, micro=None):
    """Label for one Model 1 point, from the level's own solver result.

    L1 and L2a need no result: a perfect metal is limited by metal diffusion and an oxide alone
    by oxide diffusion. L2b and L4 are one path, labelled from their solver's regime. L3 and L5
    re-solve each path (:func:`model1_paths`, which takes the remaining arguments) and raise if
    the paths do not add up to the wall flux the parallel-path solver returned in ``result``.
    """
    if level == 'L1':
        return 'metal diffusion'
    if level == 'L2a':
        return 'oxide diffusion'
    if level in ('L2b', 'L4'):
        step, detail = model1_step(result['regime'])
        return compose(step, detail if level in DEFECTIVE_METAL else '')
    paths = model1_paths(level, P_up, P_down, oxide_props, metal_props, defects, temperature, micro)
    J_total = result['flux_total']
    J_paths = sum(d['contribution'] for d in paths.values())
    if not np.isclose(J_paths, J_total, rtol=1e-6, atol=0):
        raise ValueError(f'{level}: the paths add up to {J_paths:.6e}, the wall solver gave {J_total:.6e}')
    return label_from_paths(paths, J_total)


# =============================================================================
# Model 2 (the L6 family)
# =============================================================================

def model2_step(path_result, surface):
    """The solver's ``rate_limiting`` for one path, with 'surface' named for that path's surface."""
    raw = path_result.get('rate_limiting', 'undefined')
    return surface if raw == 'surface' else _MODEL2_STEP.get(raw, 'undefined')


def model2_paths(level, result):
    """Every path of a Model 2 wall with defects (L3+L6 or L5+L6), from the solver's own breakdown.

    Returns the same structure as :func:`model1_paths`. The surface/oxide/metal fractions are the
    solver's resistance split for each path.
    """
    def described(path_result, surface, contribution):
        step = model2_step(path_result, surface)
        detail = (trap_detail(path_result['D_eff'] / path_result['D_lattice'])
                  if (level == 'L5L6' and step == 'metal diffusion') else '')
        res = path_result.get('resistances') or {}
        return {'contribution': contribution, 'flux': path_result.get('flux', np.nan), 'step': step,
                'detail': detail, **{k: res.get(k, np.nan)
                                     for k in ('fraction_surface', 'fraction_oxide', 'fraction_metal')}}

    paths = {'intact': described(result['intact_path'], SURFACE_OF['intact'],
                                 result['flux_breakdown']['intact']['contribution'])}
    for kind, d in result.get('defect_paths', {}).items():
        paths[kind] = described(d['path_result'], SURFACE_OF[kind], d['flux_contribution'])
    return paths


def label_model2(level, result, defects=None):
    """Label for one Model 2 point, from the level's own solver result.

    For L3+L6 and L5+L6 pass the ``defect_config`` the wall was solved with: the L3+L6 solver
    drops a defect path whose solve failed, which would make the point look defect-free, so a
    point missing a requested path is ``undefined`` instead.
    """
    if level in ('L3L6', 'L5L6'):
        J_total = result.get('J_total', np.nan)
        requested = [k for k, cfg in (defects or {}).items() if cfg.get('area_fraction', 0) > 0]
        if any(k not in result.get('defect_paths', {}) for k in requested) or not np.isfinite(J_total):
            return 'undefined'
        return label_from_paths(model2_paths(level, result), J_total)
    step = model2_step(result, SINGLE_PATH_SURFACE[level])
    detail = (trap_detail(result['D_eff'] / result['D_lattice'])
              if (level in DEFECTIVE_METAL and step == 'metal diffusion') else '')
    return compose(step, detail)


# =============================================================================
# Drawing a label
# =============================================================================

#: Limiting step -> (words, colour, marker). Colour means the step on every figure.
STEP_STYLE = {
    'oxide surface':   ('oxide-surface limited (H$_2$ dissociation on the oxide)',      '#1baf7a', 'P'),
    'metal surface':   ('metal-surface limited (H$_2$ dissociation on the bare metal)', '#8e5cc4', 'P'),
    'oxide diffusion': ('oxide-diffusion limited',                                     '#eb6834', 's'),
    'metal diffusion': ('metal-diffusion limited',                                     '#2a78d6', 'o'),
    'mixed':           ('no single limiting step (not fitted)',                        '#898781', 'X'),
    'undefined':       ('undefined (solve failed, not fitted)',                        '#898781', 'x'),
}
_DOMINANT_WORDS = {'intact': 'Intact oxide dominant', 'pinholes': 'Pinholes dominant',
                   'cracks': 'Cracks dominant', 'oxide GBs': 'Oxide GBs dominant'}
_DETAIL_WORDS   = {'lattice': 'lattice', 'traps': 'trap-controlled'}


def describe(label):
    """``(legend name, colour, marker, fill)`` for one label.

    Colour is the limiting step, a diamond marks trap control, and a hollow marker means a defect
    path (pinholes, cracks or oxide grain boundaries) carries most of the flux.
    """
    path, step, detail = parse(label)
    if path is None and step == 'mixed':                          # no path carries more than half
        return ('No path carries more than half the flux (not fitted)', '#898781', 'X', 'filled')
    words, colour, marker = STEP_STYLE.get(step, (step, '#52514e', 'o'))
    if detail:
        words += f', {_DETAIL_WORDS[detail]}'
        marker = 'D' if detail == 'traps' else marker
    name = f'{_DOMINANT_WORDS[path]}: {words}' if path else words[0].upper() + words[1:]
    return name, colour, marker, ('hollow' if path not in (None, 'intact') else 'filled')


def fitted(label):
    """Points with no limiting step, or whose solve failed, are drawn but never fitted."""
    return parse(label)[1] not in ('mixed', 'undefined')
