"""
Sweep one study for the cross-study comparison: every level, over temperature, pressure
and thickness, at the operating point in that study's study_config.py; and every level with a
defective oxide over the pinhole area fraction, with cracks and oxide grain boundaries set to zero.

model_config resolves exactly one study per process, and the grain-boundary model reads the
active metal from it, so each study runs in its own process:

    PERMEATION_STUDY=<study>/study_config.py  python compare_worker.py  <out.csv>

study-comparison.ipynb launches one of these per study and overlays the results. Every value --
material properties, operating point, ranges, wall, microstructure, defects -- comes from the
study's own config, even where the studies currently agree, so editing one config changes only that
study's curves. The CSV holds the flux and regime label of every point, plus what the study was run
with; Φ_app is derived from the flux in the notebook. Labels come from calculations.wall_labels, the
same module the sweep notebooks use.
"""
import contextlib
import hashlib
import io
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings('ignore')        # extrapolation notices: the inputs are the study's own

from calculations.config.model_config import (CONDITIONS, METALS, OXIDES, MICROSTRUCTURE, OXIDE_DEFECTS,
                                              STUDY_NAME, build_simulation_config,
                                              get_surface_kinetics_from_config)
from calculations.oxide_permeation import (get_oxide_properties_at_T, get_metal_properties_at_T,
                                           molecular_diffusion_flux)
from calculations.permeation_calc import calculate_simple_metal_flux, calculate_defective_metal_flux
from calculations.interface_solver import calculate_oxide_metal_system
from calculations.parallel_oxide_defect_paths import (calculate_parallel_path_flux,
                                                      calculate_parallel_path_flux_defective_metal)
from calculations.surface_kinetics import (get_all_properties, solve_steady_state_flux_L1L6,
                                           solve_steady_state_flux_L2aL6, solve_steady_state_flux_direct,
                                           calculate_mixed_defect_flux_L6, calculate_defective_metal_flux_L6,
                                           calculate_full_model_flux_L346_v2)
from calculations.wall_labels import label_model1, label_model2

SIM    = build_simulation_config()
T_op, P_up, P_down = CONDITIONS['T_operating'], CONDITIONS['P_upstream'], SIM['P_downstream']
L_ox0, L_m0, N     = SIM['L_oxide'], SIM['L_metal'], CONDITIONS['n_points']
PINHOLE_LEVELS     = ('L3', 'L5', 'L3L6', 'L5L6')         # the levels with a defective oxide


def fingerprint(block):
    """Short hash of a config block, so the notebook can tell whether two studies really share it."""
    return hashlib.sha1(json.dumps(block, sort_keys=True, default=str).encode()).hexdigest()[:10]


# What this study was run with, written onto every row so the notebook never has to assume it.
META = dict(study=STUDY_NAME, metal=SIM['metal_name'], oxide=SIM['oxide_name'],
            T_op=T_op, P_up=P_up, L_ox0=L_ox0, L_m0=L_m0,
            oxide_fp=fingerprint(OXIDES[SIM['oxide_name']]), micro_fp=fingerprint(MICROSTRUCTURE),
            defects_fp=fingerprint(OXIDE_DEFECTS), conditions_fp=fingerprint(CONDITIONS))

LAYERS = {'L1': {'metal'}, 'L2a': {'oxide'}, 'L2b': {'oxide', 'metal'}, 'L3': {'oxide', 'metal'},
          'L4': {'metal'}, 'L5': {'oxide', 'metal'},
          'L1L6': {'metal'}, 'L2aL6': {'oxide'}, 'L2bL6': {'oxide', 'metal'}, 'L3L6': {'oxide', 'metal'},
          'L4L6': {'oxide', 'metal'}, 'L5L6': {'oxide', 'metal'}}

_props = {}

def props_at(T):
    """Every temperature-dependent input, computed once per temperature."""
    if T not in _props:
        with contextlib.redirect_stdout(io.StringIO()):       # the property lookups print warnings
            _props[T] = dict(ox=get_oxide_properties_at_T(SIM['oxide_name'], T),
                             me=get_metal_properties_at_T(SIM['metal_name'], T),
                             m2=get_all_properties(SIM['oxide_name'], SIM['metal_name'], T),
                             sm=get_surface_kinetics_from_config(SIM['metal_name'], T, METALS))
    return _props[T]


def model1(level, T, P, Lox, Lm, f=None):
    """Flux and label; ``f`` replaces the study's defects with pinholes alone at that area fraction."""
    p = props_at(T)
    ox, me = dict(p['ox'], thickness=Lox), dict(p['me'], thickness=Lm)
    D, K  = me['D_metal'], me['K_s_metal']
    micro = dict(microstructure_params=MICROSTRUCTURE, lattice_density=MICROSTRUCTURE['N_L'], mode='both')
    # the single-type branch for pinholes alone, never 'mixed', which would re-weight rather than add area
    defects = OXIDE_DEFECTS if f is None else {'type': 'pinhole', 'area_fraction': f}
    if level == 'L1':  return calculate_simple_metal_flux(D, K, Lm, P, P_down)['flux'], label_model1('L1')
    if level == 'L2a': return molecular_diffusion_flux(ox['D_ox'], ox['K_ox'], Lox, P, P_down), label_model1('L2a')
    if level == 'L2b': r = calculate_oxide_metal_system(P, P_down, ox, me)
    if level == 'L3':  r = calculate_parallel_path_flux(P, P_down, ox, me, defects)
    if level == 'L4':  r = calculate_defective_metal_flux(D_lattice=D, K_s=K, metal_thickness=Lm, P_up=P,
                                                          P_down=P_down, temperature=T, **micro)
    if level == 'L5':  r = calculate_parallel_path_flux_defective_metal(P, P_down, ox, me, defects, T, **micro)
    return r.get('flux_total', r.get('flux')), label_model1(level, r, P, P_down, ox, me, defects, T, micro)


def model2(level, T, P, Lox, Lm, f=None):
    """Flux and label, exactly as in the sweep notebooks; ``f`` as in model1()."""
    p, micro = props_at(T), SIM['microstructure']
    defects  = SIM['defect_config'] if f is None else {'pinhole': {'area_fraction': f}}
    k, Keq, Dox, Kox, Dm, Ksm = (p['m2'][key] for key in ('k_diss', 'K_eq', 'D_ox', 'K_ox', 'D_m', 'K_s_m'))
    sm   = p['sm']
    pin  = {} if SIM['oxide_defects'].get('use_sieverts_pinhole', False) else dict(k_diss_metal=sm['k_diss'],
                                                                                    K_eq_metal=sm['K_eq'])
    mode = micro.get('mode', 'both')
    if level == 'L1L6':  r = solve_steady_state_flux_L1L6(P, P_down, Lm, sm['k_diss'], sm['K_eq'], Dm, Ksm)
    if level == 'L2aL6': r = solve_steady_state_flux_L2aL6(P, P_down, k, Keq, Dox, Kox, Lox)
    if level == 'L2bL6': r = solve_steady_state_flux_direct(P, P_down, Lm, k, Keq, Dox, Kox, Lox, Dm, Ksm)
    if level == 'L3L6':  r = calculate_mixed_defect_flux_L6(P, P_down, Lm, k, Keq, Dox, Kox, Lox, Dm, Ksm,
                                                            defects, **pin)
    if level == 'L4L6':  r = calculate_defective_metal_flux_L6(P, P_down, Lm, T, k, Keq, Dox, Kox, Lox, Dm, Ksm,
                                                               microstructure_params=micro, mode=mode)
    if level == 'L5L6':  r = calculate_full_model_flux_L346_v2(P, P_down, Lm, T, k, Keq, Dox, Kox, Lox, Dm, Ksm,
                                                               micro, defects, mode=mode, **pin)
    return r.get('J_total', r.get('J_ss', r.get('flux'))), label_model2(level, r, defects)


def points():
    """(sweep, T, P, L_ox, L_m, f) for every point of every sweep, around the study's operating point.

    f is the pinhole area fraction of the 'pinholes' sweep, and None (the study's own defects) elsewhere.
    """
    T_grid = np.linspace(*CONDITIONS['T_range'], N)
    P_grid = np.logspace(*np.log10(CONDITIONS['P_range']), N)
    P_grid = P_grid[P_grid > P_down]
    for T in T_grid:
        yield 'temperature', T, P_up, L_ox0, L_m0, None
    for P in P_grid:
        yield 'pressure', T_op, P, L_ox0, L_m0, None
    for Lox in np.logspace(*np.log10(CONDITIONS['L_oxide_range']), N):
        yield 'oxide', T_op, P_up, Lox, L_m0, None
    for Lm in np.logspace(*np.log10(CONDITIONS['L_metal_range']), N):
        yield 'metal', T_op, P_up, L_ox0, Lm, None
    for s in np.logspace(*np.log10(np.array(CONDITIONS['L_metal_range']) / L_m0), N):
        yield 'total', T_op, P_up, s * L_ox0, s * L_m0, None  # both layers x s, ratio fixed
    for f in np.logspace(*np.log10(CONDITIONS['defect_fraction_range']), N):
        yield 'pinholes', T_op, P_up, L_ox0, L_m0, f


def main(out):
    rows = []
    for level, layers in LAYERS.items():
        solve = model2 if level.endswith('L6') else model1
        for sweep, T, P, Lox, Lm, f in points():
            if sweep in ('oxide', 'metal') and sweep not in layers:
                continue                                       # that layer is not in this level
            if sweep == 'pinholes' and level not in PINHOLE_LEVELS:
                continue                                       # no oxide defects in this level
            with contextlib.redirect_stdout(io.StringIO()):
                flux, regime = solve(level, T, P, Lox, Lm, f)
            rows.append(dict(**META, level=level, sweep=sweep,
                             T=T, P=P, P_down=P_down, L_ox=Lox, L_m=Lm,
                             L_total=Lox * ('oxide' in layers) + Lm * ('metal' in layers),
                             f=np.nan if f is None else f, flux=flux, regime=regime))
        print(f'{STUDY_NAME}: {level} done', flush=True)
    tmp = f'{out}.part'
    pd.DataFrame(rows).to_csv(tmp, index=False)
    os.replace(tmp, out)                                       # never leave a half-written file behind


if __name__ == '__main__':
    main(sys.argv[1])
