"""Acceptance tests for :mod:`calculations.permeability`.

Numbered to match the criteria in the issue that motivated the module. The
criteria worth reading twice are 8 (Model 1 has no tier 3, and asking for one
must fail loudly), 11 (L5 and L5L6 permeability must now differ wherever their
fluxes differ -- the regression the module exists to fix) and 14 (no flux
changes).
"""

import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations import permeability as pm
from calculations.config.model_config import (
    DEFAULT_PARAMS_LEVEL5 as P5, DEFAULT_PARAMS_LEVEL5L6 as P6)

R_GAS = 8.314


def _arr(v, E, T, T_ref):
    return v * np.exp(-E / R_GAS * (1.0 / T - 1.0 / T_ref))


@pytest.fixture(scope='module')
def m1():
    """Model 1 properties at the active study's default operating point."""
    T = P5['temperature']
    return {
        'T': T,
        'D_m':  _arr(P5['D_ref'],    P5['E_D'],      T, P5['T_ref_metal']),
        'K_s':  _arr(P5['K_s_ref'],  P5['H_s'],      T, P5['T_ref_metal']),
        'D_ox': _arr(P5['D_ox_ref'], P5['E_D_ox'],   T, P5['T_ref_oxide']),
        'K_ox': _arr(P5['K_ox_ref'], P5['H_sol_ox'], T, P5['T_ref_oxide']),
        'L_m':  P5['metal_thickness'], 'L_ox': P5['oxide_thickness'],
        'P_up': P5['P_upstream'], 'P_down': P5['P_downstream'],
    }


@pytest.fixture(scope='module')
def m2():
    """Model 2 properties, including the surface kinetics."""
    T = P6['temperature']
    return {
        'T': T,
        'D_m':   _arr(P6['D_ref'],       P6['E_D'],      T, P6['T_ref_metal']),
        'K_s':   _arr(P6['K_s_ref'],     P6['H_s'],      T, P6['T_ref_metal']),
        'D_ox':  _arr(P6['D_ox_ref'],    P6['E_D_ox'],   T, P6['T_ref_oxide']),
        'K_ox':  _arr(P6['K_ox_ref'],    P6['H_sol_ox'], T, P6['T_ref_oxide']),
        'k_diss': _arr(P6['k_diss_ref'], P6['E_diss'],   T, P6['T_ref_surface']),
        'K_eq':  _arr(P6['K_eq_ref'],    P6['H_eq'],     T, P6['T_ref_surface']),
        'L_m':   P6['metal_thickness'], 'L_ox': P6['oxide_thickness'],
        'P_up':  P6['P_upstream'], 'P_down': P6['P_downstream'],
    }


@pytest.fixture(scope='module')
def micro():
    return {'grain_size': P5['grain_size'], 'grain_shape': P5['grain_shape'],
            'gb_type': P5['gb_type'], 'gb_thickness': P5['gb_thickness'],
            'trap_list': []}


# --- 1: intrinsic reproduces D*K ---------------------------------------------

def test_01_intrinsic(m1):
    assert pm.intrinsic_metal_permeability(m1['D_m'], m1['K_s'])['permeability'] \
        == m1['D_m'] * m1['K_s']
    ox1 = pm.intrinsic_oxide_permeability(m1['D_ox'], m1['K_ox'], pm.MODEL_1)
    ox2 = pm.intrinsic_oxide_permeability(m1['D_ox'], m1['K_ox'], pm.MODEL_2)
    assert ox1['permeability'] == ox2['permeability'] == m1['D_ox'] * m1['K_ox']
    assert ox1['units'] == pm.UNITS_HENRY        # molecular H2 -> per Pa
    assert ox2['units'] == pm.UNITS_SIEVERTS     # atomic H     -> per Pa^0.5
    # L1's apparent value collapses onto the intrinsic one: a single Sieverts layer
    r = pm.level1(m1['D_m'], m1['K_s'], m1['L_m'], m1['P_up'], m1['P_down'])
    assert np.isclose(r['permeability'], m1['D_m'] * m1['K_s'], rtol=1e-12)


# --- 2: effective metal reproduces D_eff*K_s ---------------------------------

def test_02_effective_metal(m1, micro):
    r = pm.level4(m1['D_m'], m1['K_s'], m1['L_m'], m1['P_up'], m1['P_down'],
                  m1['T'], micro)
    assert np.isclose(r['permeability'], r['D_eff'] * m1['K_s'], rtol=1e-12)
    assert pm.effective_metal_permeability(r['D_eff'], m1['K_s'])['permeability'] \
        == r['D_eff'] * m1['K_s']


# --- 3 / 3b: oxide-only closed form, and the coupled lumping ----------------

def test_03_oxide_only_closed_form(m1):
    Phi_ox, gamma, beta = m1['D_ox'] * m1['K_ox'], 0.1, 10.0
    for fc, fg in [(0, 0), (1e-3, 0), (0, 1e-3), (1e-2, 1e-2), (5e-2, 5e-2)]:
        got = pm.oxide_only_permeability(Phi_ox, f_crack=fc, gamma=gamma,
                                         f_gb=fg, beta=beta)
        want = (1 - fc - fg) + fc / gamma + fg * beta
        assert np.isclose(got['enhancement_ratio'], want, rtol=0, atol=0)


def test_03_oxide_only_matches_solver_with_free_metal(m1):
    """The closed form must reproduce the parallel-path solver when the metal is free."""
    from calculations.parallel_oxide_defect_paths import calculate_parallel_path_flux
    Phi_ox, gamma, beta = m1['D_ox'] * m1['K_ox'], 0.1, 10.0
    L_ox = 1e-6
    ox = {'D_ox': m1['D_ox'], 'K_ox': m1['K_ox'], 'thickness': L_ox}
    me = {'D_metal': m1['D_m'] * 1e12, 'K_s_metal': m1['K_s'], 'thickness': m1['L_m']}
    for fc, fg in [(0, 0), (1e-3, 0), (0, 1e-3), (1e-2, 1e-2), (5e-2, 5e-2)]:
        dp = {'area_fraction': fc + fg, 'type': 'mixed',
              'components': {'pinholes': 0, 'cracks': fc, 'grain_boundaries': fg},
              'thickness_factor': gamma, 'diffusivity_factor': beta}
        J_code = calculate_parallel_path_flux(m1['P_up'], m1['P_down'], ox, me, dp)['flux_total']
        Phi_eff = pm.oxide_only_permeability(Phi_ox, f_crack=fc, gamma=gamma,
                                            f_gb=fg, beta=beta)['permeability']
        J_closed = Phi_eff * (m1['P_up'] - m1['P_down']) / L_ox   # Henry
        assert abs(J_closed / J_code - 1.0) < 1e-6


def test_03_oxide_only_rejects_pinhole(m1):
    """A free-standing oxide with a pinhole is aperture flow, not permeation."""
    with pytest.raises(ValueError, match='aperture flow'):
        pm.oxide_only_permeability(m1['D_ox'] * m1['K_ox'], f_pinhole=1e-6)


def test_03b_lumped_validity_flag(m1):
    Phi_ox, Phi_m = m1['D_ox'] * m1['K_ox'], m1['D_m'] * m1['K_s']
    # metal-dominated -> lumping a pinhole is allowed
    ok = pm.lumped_oxide_permeability(Phi_ox, 1e-6, Phi_m, m1['L_m'], f_pinhole=1e-6)
    assert ok['R_ox_over_R_m'] < 10.0          # 3.42 -> measured error 0.5%
    # oxide-dominated -> refused
    bad = pm.lumped_oxide_permeability(Phi_ox * 1e-6, 1e-6, Phi_m, m1['L_m'], f_pinhole=1e-6)
    assert ok['valid'] is True and bad['valid'] is False
    assert bad['R_ox_over_R_m'] > ok['R_ox_over_R_m']
    assert np.isinf(ok['permeability'])          # pinhole -> zero oxide resistance


# --- 4: Model 1 closed-form bilayer vs brentq --------------------------------

@pytest.mark.parametrize('L_ox', [4.8e-8, 1e-6])
@pytest.mark.parametrize('scale', [1.0, 1e-4])
@pytest.mark.parametrize('P_up', [1e3, 1e5, 1e7])
def test_04_model1_closed_form(m1, L_ox, scale, P_up):
    from calculations.interface_solver import calculate_oxide_metal_system
    ox = {'D_ox': m1['D_ox'] * scale, 'K_ox': m1['K_ox'], 'thickness': L_ox}
    me = {'D_metal': m1['D_m'], 'K_s_metal': m1['K_s'], 'thickness': m1['L_m']}
    J_code = calculate_oxide_metal_system(P_up, 0.0, ox, me)['flux']
    J_cf = pm.model1_bilayer_closed_form(
        ox['D_ox'] * ox['K_ox'], L_ox, m1['D_m'] * m1['K_s'], m1['L_m'], P_up, 0.0)['flux']
    assert abs(J_cf / J_code - 1.0) < 1e-9


# --- 5-7: Model 2 decomposition ----------------------------------------------

@pytest.mark.parametrize('P_up', [1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
def test_05_model2_transport_times_driving_force(m2, P_up):
    """J == Pi_transport * (g(theta) - sqrt(P_down)): the surface adds no resistance."""
    from calculations.surface_kinetics import solve_steady_state_flux_direct, g_theta
    r = solve_steady_state_flux_direct(P_up, m2['P_down'], m2['L_m'], m2['k_diss'],
                                       m2['K_eq'], m2['D_ox'], m2['K_ox'],
                                       m2['L_ox'], m2['D_m'], m2['K_s'])
    t = pm.transport_permeability(m2['D_ox'] * m2['K_ox'], m2['L_ox'],
                                  m2['D_m'] * m2['K_s'], m2['L_m'], pm.MODEL_2)
    pred = t['permeance'] * (g_theta(r['theta'], m2['K_eq']) - np.sqrt(m2['P_down']))
    assert abs(pred / r['J_ss'] - 1.0) < 1e-9


def test_06_transport_permeability_is_pressure_invariant(m2):
    vals = [pm.transport_permeability(m2['D_ox'] * m2['K_ox'], m2['L_ox'],
                                      m2['D_m'] * m2['K_s'], m2['L_m'],
                                      pm.MODEL_2)['permeability']
            for _ in range(3)]
    assert max(vals) / min(vals) - 1.0 < 1e-12
    rows = [pm.level2b_L6(p, m2['P_down'], m2['L_m'], m2['k_diss'], m2['K_eq'],
                          m2['D_ox'], m2['K_ox'], m2['L_ox'], m2['D_m'], m2['K_s'])
            for p in (1e2, 1e5, 1e7)]
    pt = [r['Phi_transport'] for r in rows]
    assert max(pt) / min(pt) - 1.0 < 1e-12       # invariant...
    pa = [r['permeability'] for r in rows]
    assert max(pa) / min(pa) > 1.01              # ...while the apparent value is not


@pytest.mark.parametrize('P_up', [1e2, 1e5, 1e7])
def test_07_phi_app_equals_transport_times_eta(m2, P_up):
    r = pm.level2b_L6(P_up, m2['P_down'], m2['L_m'], m2['k_diss'], m2['K_eq'],
                      m2['D_ox'], m2['K_ox'], m2['L_ox'], m2['D_m'], m2['K_s'])
    assert 0.0 < r['eta_surf'] <= 1.0
    assert np.isclose(r['permeability'], r['Phi_transport'] * r['eta_surf'], rtol=1e-12)


# --- 8: Model 1 has no tier 3 -------------------------------------------------

def test_08_model1_has_no_tier_three(m1):
    with pytest.raises(pm.TierError, match='tier 3 does not exist for Model 1'):
        pm.transport_permeability(m1['D_ox'] * m1['K_ox'], m1['L_ox'],
                                  m1['D_m'] * m1['K_s'], m1['L_m'], pm.MODEL_1)


# --- 9: L5 with no defects collapses to L2b ----------------------------------

def test_09_l5_no_defects_equals_l2b(m1, micro):
    ox = {'D_ox': m1['D_ox'], 'K_ox': m1['K_ox'], 'thickness': m1['L_ox']}
    me = {'D_metal': m1['D_m'], 'K_s_metal': m1['K_s'], 'thickness': m1['L_m']}
    dp = {'area_fraction': 0.0, 'type': 'mixed',
          'components': {'pinholes': 0.0, 'cracks': 0.0, 'grain_boundaries': 0.0},
          'thickness_factor': 0.1, 'diffusivity_factor': 10.0}
    a = pm.level5(ox, me, dp, m1['P_up'], m1['P_down'], m1['T'], micro, mode='none')
    b = pm.level2b(ox, me, m1['P_up'], m1['P_down'], m1['T'])
    assert abs(a['permeability'] / b['permeability'] - 1.0) < 1e-6


# --- 10: PRF is the permeance ratio ------------------------------------------

def test_10_prf_is_permeance_ratio(m1, micro):
    from calculations.permeation_calc import calculate_defective_metal_flux
    ox = {'D_ox': m1['D_ox'], 'K_ox': m1['K_ox'], 'thickness': m1['L_ox']}
    me = {'D_metal': m1['D_m'], 'K_s_metal': m1['K_s'], 'thickness': m1['L_m']}
    dp = {'area_fraction': 1e-4, 'type': 'mixed',
          'components': {'pinholes': 1e-4, 'cracks': 0.0, 'grain_boundaries': 0.0},
          'thickness_factor': 0.1, 'diffusivity_factor': 10.0}
    coated = pm.level5(ox, me, dp, m1['P_up'], m1['P_down'], m1['T'], micro)
    bare_flux = calculate_defective_metal_flux(
        D_lattice=m1['D_m'], K_s=m1['K_s'], thickness=m1['L_m'],
        P_up=m1['P_up'], P_down=m1['P_down'], temperature=m1['T'],
        microstructure_params=micro,
        lattice_density=P5['lattice_density'])['flux']
    PRF = bare_flux / coated['flux']
    ratio = pm.permeance(bare_flux, m1['P_up'], m1['P_down']) / coated['permeance']
    assert np.isclose(PRF, ratio, rtol=1e-12)


# --- 11: the regression -- L5 and L5L6 must now differ -----------------------

def test_11_l5_and_l5l6_permeability_differ(m1, m2, micro):
    """The old metric returned identical values for these two. It must not now."""
    ox = {'D_ox': m1['D_ox'], 'K_ox': m1['K_ox'], 'thickness': m1['L_ox']}
    me = {'D_metal': m1['D_m'], 'K_s_metal': m1['K_s'], 'thickness': m1['L_m']}
    dp = {'area_fraction': 1e-4, 'type': 'mixed',
          'components': {'pinholes': 1e-4, 'cracks': 0.0, 'grain_boundaries': 0.0},
          'thickness_factor': 0.1, 'diffusivity_factor': 10.0}
    a = pm.level5(ox, me, dp, m1['P_up'], m1['P_down'], m1['T'], micro)
    b = pm.level5_L6(m2['P_up'], m2['P_down'], m2['L_m'], m2['T'], m2['k_diss'],
                     m2['K_eq'], m2['D_ox'], m2['K_ox'], m2['L_ox'], m2['D_m'],
                     m2['K_s'], micro,
                     {'pinhole': {'area_fraction': 1e-4}})
    # old behaviour: identical. new: they track their own fluxes.
    assert not np.isclose(a['permeability'], b['permeability'], rtol=1e-6, atol=0.0)
    assert np.isclose(a['permeability'] / b['permeability'],
                      a['flux'] / b['flux'], rtol=1e-9)


# --- 12: tier hygiene ---------------------------------------------------------

def test_12_tier_hygiene(m1, m2, micro):
    t12 = [pm.level1(m1['D_m'], m1['K_s'], m1['L_m'], m1['P_up'], m1['P_down']),
           pm.level4(m1['D_m'], m1['K_s'], m1['L_m'], m1['P_up'], m1['P_down'],
                     m1['T'], micro)]
    for r in t12:
        assert r['tier'] in (1, 2)
        assert r['pressure_exponent'] is None
    ox = {'D_ox': m1['D_ox'], 'K_ox': m1['K_ox'], 'thickness': m1['L_ox']}
    me = {'D_metal': m1['D_m'], 'K_s_metal': m1['K_s'], 'thickness': m1['L_m']}
    r4 = pm.level2b(ox, me, m1['P_up'], m1['P_down'], m1['T'])
    assert r4['tier'] == 4
    for k in ('P_up', 'P_down', 'temperature'):
        assert r4[k] is not None
    assert r4['pressure_exponent'] is not None


# --- 13: refuse to mix sorption laws -----------------------------------------

def test_13_unit_mismatch_raises(m1):
    ox = pm.intrinsic_oxide_permeability(m1['D_ox'], m1['K_ox'], pm.MODEL_1)
    me = pm.intrinsic_metal_permeability(m1['D_m'], m1['K_s'])
    with pytest.raises(pm.UnitMismatchError, match='not commensurable'):
        pm.combine_series(ox, me, m1['L_ox'], m1['L_m'])
    ox2 = pm.intrinsic_oxide_permeability(m1['D_ox'], m1['K_ox'], pm.MODEL_2)
    out = pm.combine_series(ox2, me, m1['L_ox'], m1['L_m'])          # Model 2: allowed
    assert out['tier'] == 3


# --- 14: no flux changed ------------------------------------------------------

def test_14_fluxes_are_untouched(m1, micro):
    """Every level's reported flux must be the solver's own, bit-for-bit."""
    from calculations.permeation_calc import calculate_simple_metal_flux
    from calculations.interface_solver import calculate_oxide_metal_system
    ox = {'D_ox': m1['D_ox'], 'K_ox': m1['K_ox'], 'thickness': m1['L_ox']}
    me = {'D_metal': m1['D_m'], 'K_s_metal': m1['K_s'], 'thickness': m1['L_m']}
    assert pm.level1(m1['D_m'], m1['K_s'], m1['L_m'], m1['P_up'], m1['P_down'])['flux'] \
        == calculate_simple_metal_flux(m1['D_m'], m1['K_s'], m1['L_m'],
                                       m1['P_up'], m1['P_down'])['flux']
    assert pm.level2b(ox, me, m1['P_up'], m1['P_down'])['flux'] \
        == calculate_oxide_metal_system(m1['P_up'], m1['P_down'], ox, me)['flux']
