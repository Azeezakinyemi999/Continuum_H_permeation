"""The metal's steady-state diffusivity: lattice, raised by grain-boundary fast paths.

Hydrogen trapping was removed from the model because it cannot change a steady-state
flux: with the traps in equilibrium and no net exchange, only lattice hydrogen moves and
J = D·ΔC_L/L whatever the trap density (see STEADY_STATE.md). Grain-boundary enhancement
stays, and it does not depend on concentration, so every Level 4/5 solver reduces to its
Level 1-3 counterpart run with one diffusivity, D_gb. These tests pin that reduction.

The study defaults use low-angle boundaries at 100 µm, where the enhancement is ~1e-5.
A fine-grained, high-angle microstructure is used here so the enhancement is large
enough that a dropped or doubled factor would show.
"""
import math
import pathlib
import re

import numpy as np
import pytest

from calculations import permeability as pm
from calculations.config.model_config import DEFAULT_PARAMS_LEVEL5 as P5
from calculations.defective_metal import metal_diffusivity
from calculations.parallel_oxide_defect_paths import calculate_parallel_path_flux
from calculations.permeation_calc import (
    calculate_defective_metal_flux, calculate_simple_metal_flux,
)
from calculations.sensitivity import level5_model_wrapper
from calculations.surface_kinetics import solve_steady_state_flux_direct

ROOT = pathlib.Path(__file__).resolve().parents[1]

FINE_HAGB = {'grain_size': 1e-6, 'gb_thickness': 0.5e-9,
             'grain_shape': 'equiaxed', 'gb_type': 'HAGB'}

T = 700.0
D_L, K_S, L_M = 1e-10, 0.05, 1e-3
P_UP, P_DOWN = 1e5, 10.0


def test_fine_hagb_really_enhances():
    """The premise of the tests below, asserted rather than assumed."""
    assert metal_diffusivity(D_L, T, FINE_HAGB) / D_L > 1.01


def test_disabled_enhancement_returns_the_lattice_value():
    assert metal_diffusivity(D_L, T, FINE_HAGB, enable_gb_enhancement=False) == D_L


def test_level4_flux_is_fick_with_the_enhanced_diffusivity():
    D_gb = metal_diffusivity(D_L, T, FINE_HAGB)
    r = calculate_defective_metal_flux(D_L, K_S, L_M, P_UP, P_DOWN, T, FINE_HAGB)
    expected = D_gb * K_S * (math.sqrt(P_UP) - math.sqrt(P_DOWN)) / L_M
    assert r['flux'] == pytest.approx(expected, rel=1e-12)
    assert r['modification_factor'] == pytest.approx(D_gb / D_L, rel=1e-12)


def test_level4_without_enhancement_is_level1():
    a = calculate_defective_metal_flux(D_L, K_S, L_M, P_UP, P_DOWN, T, FINE_HAGB,
                                       enable_gb_enhancement=False)['flux']
    b = calculate_simple_metal_flux(D_L, K_S, L_M, P_UP, P_DOWN)['flux']
    assert a == b


def test_level5_wrapper_without_enhancement_is_level3():
    """The full L5 wrapper, GB off, is exactly the L3 parallel-path assembly."""
    p = dict(P5, enable_gb_enhancement=False)
    T5 = p['temperature']

    def arr(v, E, T_ref):
        return v * np.exp(-E / 8.314 * (1.0 / T5 - 1.0 / T_ref))

    oxide = {'D_ox': arr(p['D_ox_ref'], p['E_D_ox'], p['T_ref_oxide']),
             'K_ox': arr(p['K_ox_ref'], p['H_sol_ox'], p['T_ref_oxide']),
             'thickness': p['oxide_thickness']}
    metal = {'D_metal': arr(p['D_ref'], p['E_D'], p['T_ref_metal']),
             'K_s_metal': arr(p['K_s_ref'], p['H_s'], p['T_ref_metal']),
             'thickness': p['metal_thickness']}
    defects = {'area_fraction': p['f_pinhole'] + p['f_crack'] + p['f_gb_defect'],
               'type': 'mixed',
               'components': {'pinholes': p['f_pinhole'], 'cracks': p['f_crack'],
                              'grain_boundaries': p['f_gb_defect']},
               'thickness_factor': p['crack_thickness_factor'],
               'diffusivity_factor': p['gb_diffusivity_factor']}
    l3 = calculate_parallel_path_flux(p['P_upstream'], p['P_downstream'],
                                      oxide, metal, defects)['flux_total']
    l5 = level5_model_wrapper(p)
    assert l5['flux'] == pytest.approx(l3, rel=1e-12)
    assert l5['D_modification'] == 1.0


def test_level4_l6_is_the_l2_l6_solver_with_the_enhanced_diffusivity():
    k_diss, K_eq = 1e-9, 1e-4
    D_ox, K_ox, L_ox = 1e-18, 0.1, 1e-6
    D_gb = metal_diffusivity(D_L, T, FINE_HAGB)
    row = pm.level4_L6(P_UP, P_DOWN, L_M, T, k_diss, K_eq, D_ox, K_ox, L_ox,
                       D_L, K_S, FINE_HAGB)
    direct = solve_steady_state_flux_direct(P_UP, P_DOWN, L_M, k_diss, K_eq,
                                            D_ox, K_ox, L_ox, D_gb, K_S)
    assert row['flux'] == pytest.approx(direct['J_ss'], rel=1e-12)


def test_d_modification_is_never_below_one():
    """Grain boundaries can only add fast paths; nothing in the model slows the metal."""
    for gs in (1e-7, 1e-6, 1e-5, 1e-4):
        for gb in ('HAGB', 'LAGB', 'twin', 'special'):
            ms = dict(FINE_HAGB, grain_size=gs, gb_type=gb)
            assert metal_diffusivity(D_L, T, ms) >= D_L


TRAPPING_NAMES = re.compile(
    r"trap_list|trap_\w+_(?:N_T|E_b)|lattice_density|theta_total|D_eff_method"
    r"|enable_trapping|enable_gb_as_trap|combined_microstructure_model"
    r"|calculate_effective_diffusivity_trapping|trap_occupancy")


@pytest.mark.parametrize("path", sorted(
    [*ROOT.glob("calculations/**/*.py"), *ROOT.glob("*/study_config.py")]),
    ids=lambda p: str(p.relative_to(ROOT)))
def test_no_trapping_code_remains(path):
    """Trapping does not change steady-state flux and was removed; keep it out.

    If it is ever wanted again — for inventory or time lag, which are transient
    quantities — it belongs in a separate post-processing step, not in the flux path.
    """
    hits = [f"{i}: {line.strip()}" for i, line in
            enumerate(path.read_text().splitlines(), 1) if TRAPPING_NAMES.search(line)]
    assert not hits, "\n".join(hits)
