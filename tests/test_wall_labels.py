"""The regime labels in calculations.wall_labels.

A label makes two claims about a point: which path carries most of the wall's flux, and which
step limits that path. These tests pin the bookkeeping the module adds on top of the solvers
rather than the solvers' own physics: that a path dominates only above half the flux, that each
surface is named for the surface its path actually has, that the Model 1 paths re-solved one by
one add up to the flux the wall solver returned, and that a dropped defect path is never read
as a defect-free wall. Which step the solvers call limiting is theirs, and is not re-tested here.
"""
import math

import numpy as np
import pytest

from calculations.config.model_config import (CONDITIONS, METALS, MICROSTRUCTURE, OXIDE_DEFECTS,
                                              build_simulation_config, get_surface_kinetics_from_config)
from calculations.oxide_permeation import get_oxide_properties_at_T, get_metal_properties_at_T
from calculations.parallel_oxide_defect_paths import (calculate_parallel_path_flux,
                                                      calculate_parallel_path_flux_defective_metal)
from calculations.surface_kinetics import (get_all_properties, calculate_mixed_defect_flux_L6,
                                           calculate_full_model_flux_L346_v2)
from calculations.wall_labels import (DOMINANCE, SURFACE_OF, compose, defect_fractions, describe, fitted,
                                      label_from_paths, label_model1, label_model2, model1_paths,
                                      model1_step, model2_paths, parse)

SIM = build_simulation_config()
MICRO_M1 = dict(microstructure_params=MICROSTRUCTURE, lattice_density=MICROSTRUCTURE['N_L'], mode='both')

# The operating point, and the low-pressure end where the oxide (Model 1) or the surface
# (Model 2) holds hydrogen back, so both kinds of label are exercised.
POINTS = [(CONDITIONS['T_operating'], CONDITIONS['P_upstream']),
          (CONDITIONS['T_operating'], CONDITIONS['P_range'][0])]


def layers(T):
    return (dict(get_oxide_properties_at_T(SIM['oxide_name'], T), thickness=SIM['L_oxide']),
            dict(get_metal_properties_at_T(SIM['metal_name'], T), thickness=SIM['L_metal']))


def model2_wall(level, T, P, defects):
    p  = get_all_properties(SIM['oxide_name'], SIM['metal_name'], T)
    sm = get_surface_kinetics_from_config(SIM['metal_name'], T, METALS)
    args = (P, SIM['P_downstream'], SIM['L_metal'])
    kin  = dict(k_diss_metal=sm['k_diss'], K_eq_metal=sm['K_eq'])
    if level == 'L3L6':
        return calculate_mixed_defect_flux_L6(*args, p['k_diss'], p['K_eq'], p['D_ox'], p['K_ox'], SIM['L_oxide'],
                                              p['D_m'], p['K_s_m'], defects, **kin)
    return calculate_full_model_flux_L346_v2(*args, T, p['k_diss'], p['K_eq'], p['D_ox'], p['K_ox'], SIM['L_oxide'],
                                             p['D_m'], p['K_s_m'], SIM['microstructure'], defects, mode='both', **kin)


@pytest.mark.parametrize('regime,expected', [
    ('oxide_limited',                      ('oxide diffusion', '')),
    ('oxide_limited/defect_limited',       ('oxide diffusion', '')),
    ('metal_limited',                      ('metal diffusion', '')),
    ('metal_limited/lattice_limited',      ('metal diffusion', 'lattice')),
    ('metal_limited/traps_defect_limited', ('metal diffusion', 'traps')),
    ('transition',                         ('mixed', '')),
])
def test_model1_step_names_the_solver_label(regime, expected):
    assert model1_step(regime) == expected


@pytest.mark.parametrize('step,detail,path', [('oxide surface', '', None), ('metal diffusion', 'traps', 'pinhole'),
                                              ('mixed', '', 'intact'), ('oxide diffusion', '', 'grain_boundary')])
def test_compose_and_parse_are_inverses(step, detail, path):
    from calculations.wall_labels import PATH_NAMES
    assert parse(compose(step, detail, path)) == ((PATH_NAMES[path] if path else None), step, detail)


def test_a_path_dominates_only_above_half_the_flux():
    """At exactly half neither path dominates; a hair above, the larger one does."""
    def paths(share):
        return {'intact':  {'contribution': 1 - share, 'step': 'oxide diffusion', 'detail': ''},
                'pinhole': {'contribution': share,     'step': 'metal surface',   'detail': ''}}
    assert label_from_paths(paths(DOMINANCE), 1.0) == 'mixed'
    assert label_from_paths(paths(DOMINANCE + 1e-9), 1.0) == 'pinholes/metal surface'
    assert label_from_paths(paths(DOMINANCE - 1e-9), 1.0) == 'intact/oxide diffusion'
    assert label_from_paths(paths(0.3), math.nan) == 'undefined'


def test_mixed_defects_split_the_area_like_the_solver():
    fractions = defect_fractions(OXIDE_DEFECTS)
    assert math.isclose(sum(fractions.values()), OXIDE_DEFECTS['area_fraction'])
    components = OXIDE_DEFECTS['components']
    assert math.isclose(fractions['pinhole'] / fractions['crack'], components['pinholes'] / components['cracks'])
    assert defect_fractions({'type': 'crack', 'area_fraction': 0.2}) == {'crack': 0.2}
    assert defect_fractions({'type': 'mixed', 'area_fraction': 0.2, 'components': {}}) == {'pinhole': 0.2}


@pytest.mark.parametrize('T,P', POINTS)
@pytest.mark.parametrize('level', ['L3', 'L5'])
def test_model1_paths_add_up_to_the_wall_flux(level, T, P):
    """Re-solving each path alone must reproduce the wall the parallel-path solver assembled."""
    ox, me = layers(T)
    if level == 'L3':
        wall  = calculate_parallel_path_flux(P, SIM['P_downstream'], ox, me, OXIDE_DEFECTS)
        paths = model1_paths(level, P, SIM['P_downstream'], ox, me, OXIDE_DEFECTS)
    else:
        wall  = calculate_parallel_path_flux_defective_metal(P, SIM['P_downstream'], ox, me, OXIDE_DEFECTS, T, **MICRO_M1)
        paths = model1_paths(level, P, SIM['P_downstream'], ox, me, OXIDE_DEFECTS, T, MICRO_M1)
    assert set(paths) == {'intact', 'pinhole', 'crack', 'grain_boundary'}
    assert np.isclose(sum(d['contribution'] for d in paths.values()), wall['flux_total'], rtol=1e-9)
    assert np.isclose(paths['intact']['flux'], wall['flux_intact_per_area'], rtol=1e-9)


@pytest.mark.parametrize('T,P', POINTS)
@pytest.mark.parametrize('level', ['L3L6', 'L5L6'])
def test_each_surface_is_named_for_its_own_path(level, T, P):
    """Where the solver says 'surface', the label names the surface that path really has."""
    r = model2_wall(level, T, P, SIM['defect_config'])
    for kind, d in model2_paths(level, r).items():
        raw = r['intact_path'] if kind == 'intact' else r['defect_paths'][kind]['path_result']
        if raw['rate_limiting'] == 'surface':
            assert d['step'] == SURFACE_OF[kind]
    assert np.isclose(sum(d['contribution'] for d in model2_paths(level, r).values()), r['J_total'], rtol=1e-9)
    assert SURFACE_OF['pinhole'] == 'metal surface' and SURFACE_OF['intact'] == 'oxide surface'


def test_a_dropped_defect_path_is_undefined_not_defect_free():
    """calculate_mixed_defect_flux_L6 skips a path whose solve failed; the label must not hide that."""
    r = model2_wall('L3L6', *POINTS[0], {'pinhole': {'area_fraction': 0.01}})
    assert label_model2('L3L6', r, {'pinhole': {'area_fraction': 0.01}}) != 'undefined'
    r['defect_paths'].pop('pinhole')
    assert label_model2('L3L6', r, {'pinhole': {'area_fraction': 0.01}}) == 'undefined'


def test_model1_wall_label_refuses_paths_that_do_not_add_up():
    T, P = POINTS[0]
    ox, me = layers(T)
    wall = dict(calculate_parallel_path_flux(P, SIM['P_downstream'], ox, me, OXIDE_DEFECTS))
    wall['flux_total'] *= 1.01
    with pytest.raises(ValueError, match='add up'):
        label_model1('L3', wall, P, SIM['P_downstream'], ox, me, OXIDE_DEFECTS)


def test_drawing_marks_defect_paths_hollow_and_mixed_unfitted():
    assert describe('pinholes/metal surface')[3] == 'hollow'
    assert describe('intact/oxide surface')[3] == 'filled'
    assert describe('oxide surface')[3] == 'filled'
    assert describe('intact/metal diffusion/traps')[2] == 'D'
    assert describe('pinholes/metal surface')[1] != describe('intact/oxide surface')[1]
    assert not fitted('mixed') and not fitted('intact/mixed') and not fitted('undefined')
    assert fitted('cracks/oxide diffusion')
