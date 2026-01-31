"""
Unit Tests for Level 6: Surface Kinetics Module (Dissociation Only)

Tests cover:
1. Basic function behavior and input validation
2. Physical limiting cases
3. Langmuir isotherm consistency
4. Damköhler number regime classification
5. Main function with material/explicit kinetics
6. Comparison with Level 1 (Sieverts' Law) in fast kinetics limit

Run with: python -m pytest validation/test_surface_kinetics.py -v
Or directly: python validation/test_surface_kinetics.py
"""

import numpy as np
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations.surface_kinetics import (
    surface_equilibrium_coverage,
    dissociation_flux,
    calculate_damkohler_number,
    calculate_surface_limited_flux
)
from calculations.permeation_calc import calculate_simple_metal_flux


# =============================================================================
# TEST PARAMETERS
# =============================================================================

# Standard test conditions
T_TEST = 1073  # K (800°C)
P_UP_TEST = 1e5  # Pa (1 bar)
P_DOWN_TEST = 0  # Pa
THICKNESS_TEST = 1e-3  # m (1 mm)

# Bulk properties (typical for Incoloy800 at 800°C)
D_TEST = 1e-10  # m²/s
K_S_TEST = 0.5  # mol/m³/Pa^0.5

# Surface kinetics
K_DISS_TEST = 1e-2  # dissociation rate constant
K_RECOMB_TEST = 1e-7  # recombination rate constant (for K_eq)


# =============================================================================
# TEST 1: surface_equilibrium_coverage()
# =============================================================================

def test_equilibrium_coverage_zero_pressure():
    """θ should be 0 when P = 0."""
    result = surface_equilibrium_coverage(pressure=0, k_diss=K_DISS_TEST, k_recomb=K_RECOMB_TEST)
    assert result['theta'] == 0.0, f"Expected θ=0 at P=0, got {result['theta']}"
    assert result['regime'] == 'low_coverage'
    print("✓ test_equilibrium_coverage_zero_pressure passed")


def test_equilibrium_coverage_langmuir_limit():
    """At high K_eq × P, θ should approach 1."""
    k_diss_high = 1e10
    k_recomb_low = 1e-10
    result = surface_equilibrium_coverage(pressure=1e5, k_diss=k_diss_high, k_recomb=k_recomb_low)
    assert result['theta'] > 0.99, f"Expected θ→1 at high K_eq×P, got {result['theta']}"
    assert result['regime'] == 'high_coverage'
    print("✓ test_equilibrium_coverage_langmuir_limit passed")


def test_equilibrium_coverage_low_limit():
    """At low K_eq × P, θ should approach √(K_eq × P)."""
    k_diss_low = 1e-10
    k_recomb_high = 1e-5
    P_low = 1e-3
    result = surface_equilibrium_coverage(pressure=P_low, k_diss=k_diss_low, k_recomb=k_recomb_high)
    
    K_eq = k_diss_low / k_recomb_high
    expected_theta = np.sqrt(K_eq * P_low)
    
    assert abs(result['theta'] - expected_theta) < 0.01, \
        f"Expected θ≈{expected_theta:.4f} at low coverage, got {result['theta']:.4f}"
    assert result['regime'] == 'low_coverage'
    print("✓ test_equilibrium_coverage_low_limit passed")


def test_equilibrium_coverage_input_validation():
    """Should raise ValueError for invalid inputs."""
    try:
        surface_equilibrium_coverage(pressure=-1, k_diss=K_DISS_TEST, k_recomb=K_RECOMB_TEST)
        assert False, "Should raise ValueError for negative pressure"
    except ValueError:
        pass
    
    try:
        surface_equilibrium_coverage(pressure=1e5, k_diss=-1, k_recomb=K_RECOMB_TEST)
        assert False, "Should raise ValueError for negative k_diss"
    except ValueError:
        pass
    
    print("✓ test_equilibrium_coverage_input_validation passed")


# =============================================================================
# TEST 2: dissociation_flux()
# =============================================================================

def test_dissociation_flux_zero_coverage():
    """At θ=0, J_diss = k_diss × P (no blocking)."""
    flux = dissociation_flux(pressure=P_UP_TEST, theta=0.0, k_diss=K_DISS_TEST)
    expected = K_DISS_TEST * P_UP_TEST
    assert abs(flux - expected) < 1e-10, f"Expected {expected}, got {flux}"
    print("✓ test_dissociation_flux_zero_coverage passed")


def test_dissociation_flux_full_coverage():
    """At θ=1, J_diss = 0 (complete blocking)."""
    flux = dissociation_flux(pressure=P_UP_TEST, theta=1.0, k_diss=K_DISS_TEST)
    assert flux == 0.0, f"Expected 0 at θ=1, got {flux}"
    print("✓ test_dissociation_flux_full_coverage passed")


def test_dissociation_flux_half_coverage():
    """At θ=0.5, blocking factor = (1-0.5)² = 0.25."""
    flux = dissociation_flux(pressure=P_UP_TEST, theta=0.5, k_diss=K_DISS_TEST)
    expected = K_DISS_TEST * P_UP_TEST * 0.25
    assert abs(flux - expected) < 1e-10, f"Expected {expected}, got {flux}"
    print("✓ test_dissociation_flux_half_coverage passed")


# =============================================================================
# TEST 3: calculate_damkohler_number()
# =============================================================================

def test_damkohler_diffusion_limited():
    """High Da should give diffusion-limited regime."""
    result = calculate_damkohler_number(
        k_diss=1e5, D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST
    )
    assert result['Da'] > 10, f"Expected Da > 10, got {result['Da']}"
    assert result['regime'] == 'diffusion-limited', f"Expected diffusion-limited, got {result['regime']}"
    print("✓ test_damkohler_diffusion_limited passed")


def test_damkohler_dissociation_limited():
    """Low Da should give dissociation-limited regime."""
    result = calculate_damkohler_number(
        k_diss=1e-15, D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST
    )
    assert result['Da'] < 0.1, f"Expected Da < 0.1, got {result['Da']}"
    assert result['regime'] == 'dissociation-limited', f"Expected dissociation-limited, got {result['regime']}"
    print("✓ test_damkohler_dissociation_limited passed")


def test_damkohler_mixed_regime():
    """Intermediate Da should give mixed regime."""
    k_diss_mixed = D_TEST / (K_S_TEST * THICKNESS_TEST)  # Da = 1
    result = calculate_damkohler_number(
        k_diss=k_diss_mixed, D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST
    )
    assert 0.1 <= result['Da'] <= 10, f"Expected 0.1 ≤ Da ≤ 10, got {result['Da']}"
    assert result['regime'] == 'mixed', f"Expected mixed, got {result['regime']}"
    print("✓ test_damkohler_mixed_regime passed")


# =============================================================================
# TEST 4: calculate_surface_limited_flux() - Main function
# =============================================================================

def test_main_function_with_material():
    """Main function should work with material name."""
    result = calculate_surface_limited_flux(
        D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        material_name='Incoloy800'
    )
    assert 'flux' in result
    assert 'theta' in result
    assert 'damkohler' in result
    assert result['coverage_mode'] == 'calculated'
    print(f"✓ test_main_function_with_material passed (flux={result['flux']:.2e})")


def test_main_function_with_explicit_kinetics():
    """Main function should work with explicit kinetics parameters."""
    result = calculate_surface_limited_flux(
        D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        k_diss=K_DISS_TEST, k_recomb=K_RECOMB_TEST
    )
    assert 'flux' in result
    assert result['coverage_mode'] == 'calculated'
    print(f"✓ test_main_function_with_explicit_kinetics passed (flux={result['flux']:.2e})")


def test_main_function_fixed_coverage():
    """Main function should work with fixed coverage input."""
    result = calculate_surface_limited_flux(
        D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        k_diss=K_DISS_TEST, k_recomb=K_RECOMB_TEST,
        theta_fixed=0.5
    )
    assert result['coverage_mode'] == 'fixed'
    assert result['theta'] == 0.5
    assert result['converged'] == True
    print(f"✓ test_main_function_fixed_coverage passed (flux={result['flux']:.2e})")


# =============================================================================
# TEST 5: Comparison with Level 1 (fast kinetics limit)
# =============================================================================

def test_fast_kinetics_approaches_sieverts():
    """With very fast kinetics, Level 6 should approach Level 1 (Sieverts' Law)."""
    k_diss_fast = 1e10
    k_recomb_fast = 1e5
    
    result_L6 = calculate_surface_limited_flux(
        D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        k_diss=k_diss_fast, k_recomb=k_recomb_fast
    )
    
    result_L1 = calculate_simple_metal_flux(
        D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST
    )
    
    flux_L6 = result_L6['flux']
    flux_L1 = result_L1['flux']
    SRF = result_L6['surface_reduction_factor']
    
    print(f"  Level 6 flux: {flux_L6:.4e} mol/m²/s")
    print(f"  Level 1 flux: {flux_L1:.4e} mol/m²/s")
    print(f"  Surface reduction factor: {SRF:.3f}")
    print(f"  θ: {result_L6['theta']:.4f}")
    print(f"  Damköhler regime: {result_L6['damkohler']['regime']}")
    
    assert SRF == 1.0, f"Expected SRF = 1.0 for fast kinetics, got {SRF}"
    assert flux_L6 == flux_L1, f"Expected L6 flux = L1 flux for fast kinetics"
    print("✓ test_fast_kinetics_approaches_sieverts passed")


def test_fast_kinetics_approaches_level4():
    """With very fast kinetics, Level 6 should approach Level 4."""
    from calculations.permeation_calc import calculate_defective_metal_flux
    
    k_diss_fast = 1e10
    k_recomb_fast = 1e5
    
    microstructure_params = {
        'grain_size': 50e-6,
        'grain_shape': 'equiaxed',
        'gb_type': 'HAGB',
        'trap_list': [
            {'name': 'dislocations', 'density': 1e14, 'binding_energy': 27e3},
        ]
    }
    
    result_L4 = calculate_defective_metal_flux(
        D_lattice=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        microstructure_params=microstructure_params
    )
    
    D_eff_L4 = result_L4['D_eff']
    
    result_L6 = calculate_surface_limited_flux(
        D=D_eff_L4, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        k_diss=k_diss_fast, k_recomb=k_recomb_fast
    )
    
    flux_L4 = result_L4['flux']
    flux_L6 = result_L6['flux']
    
    print(f"  Level 4 flux: {flux_L4:.4e} mol/m²/s (D_eff={D_eff_L4:.2e})")
    print(f"  Level 6 flux: {flux_L6:.4e} mol/m²/s")
    print(f"  Level 6 SRF: {result_L6['surface_reduction_factor']:.3f}")
    
    assert result_L6['surface_reduction_factor'] == 1.0
    assert flux_L6 == flux_L4, f"Expected L6 = L4 for fast kinetics"
    
    print("✓ test_fast_kinetics_approaches_level4 passed")


def test_slow_kinetics_reduces_flux():
    """With slow kinetics, Level 6 flux should be less than Level 1 and Level 4."""
    from calculations.permeation_calc import calculate_defective_metal_flux
    
    # Need J_diss < J_Sieverts for surface limitation
    # J_Sieverts ~ D*K_s*√P/L ~ 1e-10 * 0.5 * 316 / 1e-3 ~ 1.6e-5 mol/m²/s
    # J_diss = k_diss × P × (1-θ)²
    # For J_diss ~ 1e-6 with P=1e5, (1-θ)²~0.25: k_diss ~ 4e-11
    k_diss_slow = 1e-11
    k_recomb_slow = 1e-3  # K_eq = 1e-8 (low coverage)
    
    K_eq = k_diss_slow / k_recomb_slow
    sqrt_K_eq_P = np.sqrt(K_eq * P_UP_TEST)
    theta_expected = sqrt_K_eq_P / (1 + sqrt_K_eq_P)
    
    print(f"  Design: K_eq = {K_eq:.2e}")
    print(f"  √(K_eq×P) = {sqrt_K_eq_P:.3f}")
    print(f"  Expected θ = {theta_expected:.4f}")
    
    # Level 1: Simple metal (Sieverts)
    result_L1 = calculate_simple_metal_flux(
        D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST
    )
    
    # Level 4: Defective metal
    microstructure_params = {
        'grain_size': 50e-6,
        'grain_shape': 'equiaxed',
        'gb_type': 'HAGB',
        'trap_list': [
            {'name': 'dislocations', 'density': 1e14, 'binding_energy': 27e3},
        ]
    }
    
    result_L4 = calculate_defective_metal_flux(
        D_lattice=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        microstructure_params=microstructure_params
    )
    
    D_eff_L4 = result_L4['D_eff']
    
    # Level 6 with slow kinetics (using Level 1 D for comparison)
    result_L6_vs_L1 = calculate_surface_limited_flux(
        D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        k_diss=k_diss_slow, k_recomb=k_recomb_slow
    )
    
    # Level 6 with slow kinetics (using Level 4 D_eff)
    result_L6_vs_L4 = calculate_surface_limited_flux(
        D=D_eff_L4, K_s=K_S_TEST, thickness=THICKNESS_TEST,
        P_up=P_UP_TEST, P_down=P_DOWN_TEST,
        temperature=T_TEST,
        k_diss=k_diss_slow, k_recomb=k_recomb_slow
    )
    
    flux_L1 = result_L1['flux']
    flux_L4 = result_L4['flux']
    flux_L6_vs_L1 = result_L6_vs_L1['flux']
    flux_L6_vs_L4 = result_L6_vs_L4['flux']
    
    SRF_vs_L1 = result_L6_vs_L1['surface_reduction_factor']
    SRF_vs_L4 = result_L6_vs_L4['surface_reduction_factor']
    J_diss = result_L6_vs_L1['flux_dissociation']
    
    reduction_vs_L1 = (1 - flux_L6_vs_L1 / flux_L1) * 100 if flux_L1 > 0 else 0
    reduction_vs_L4 = (1 - flux_L6_vs_L4 / flux_L4) * 100 if flux_L4 > 0 else 0
    
    print(f"\n  --- Comparison with Level 1 (Sieverts) ---")
    print(f"  Level 1 flux: {flux_L1:.4e} mol/m²/s")
    print(f"  Level 6 flux: {flux_L6_vs_L1:.4e} mol/m²/s")
    print(f"  J_diss: {J_diss:.4e} mol/m²/s")
    print(f"  Flux reduction vs L1: {reduction_vs_L1:.1f}%")
    print(f"  SRF vs L1: {SRF_vs_L1:.6f}")
    
    print(f"\n  --- Comparison with Level 4 (Defective Metal) ---")
    print(f"  Level 4 flux: {flux_L4:.4e} mol/m²/s (D_eff={D_eff_L4:.2e})")
    print(f"  Level 6 flux: {flux_L6_vs_L4:.4e} mol/m²/s")
    print(f"  Flux reduction vs L4: {reduction_vs_L4:.1f}%")
    print(f"  SRF vs L4: {SRF_vs_L4:.6f}")
    
    print(f"\n  --- Surface State ---")
    print(f"  θ: {result_L6_vs_L1['theta']:.4f}")
    print(f"  Damköhler regime: {result_L6_vs_L1['damkohler']['regime']}")
    print(f"  Da: {result_L6_vs_L1['damkohler']['Da']:.2e}")
    
    # Assertions
    assert flux_L6_vs_L1 < flux_L1, "Expected Level 6 flux < Level 1 flux"
    assert flux_L6_vs_L4 < flux_L4, "Expected Level 6 flux < Level 4 flux"
    assert SRF_vs_L1 < 1.0, f"Expected SRF < 1 vs L1, got {SRF_vs_L1}"
    assert SRF_vs_L4 < 1.0, f"Expected SRF < 1 vs L4, got {SRF_vs_L4}"
    assert reduction_vs_L1 > 10, f"Expected >10% reduction vs L1, got {reduction_vs_L1:.1f}%"
    assert reduction_vs_L4 > 10, f"Expected >10% reduction vs L4, got {reduction_vs_L4:.1f}%"
    
    print("✓ test_slow_kinetics_reduces_flux passed")


# =============================================================================
# TEST 6: Input validation
# =============================================================================

def test_input_validation_negative_pressure():
    """Should raise error for negative pressure."""
    try:
        calculate_surface_limited_flux(
            D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
            P_up=-1, P_down=P_DOWN_TEST,
            temperature=T_TEST,
            material_name='Incoloy800'
        )
        assert False, "Should raise ValueError for negative pressure"
    except ValueError:
        pass
    print("✓ test_input_validation_negative_pressure passed")


def test_input_validation_missing_material():
    """Should raise error if material_name missing and kinetics not provided."""
    try:
        calculate_surface_limited_flux(
            D=D_TEST, K_s=K_S_TEST, thickness=THICKNESS_TEST,
            P_up=P_UP_TEST, P_down=P_DOWN_TEST,
            temperature=T_TEST
        )
        assert False, "Should raise ValueError for missing parameters"
    except ValueError:
        pass
    print("✓ test_input_validation_missing_material passed")


# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all_tests():
    """Run all unit tests."""
    print("=" * 70)
    print("LEVEL 6 SURFACE KINETICS UNIT TESTS (Dissociation Only)")
    print("=" * 70)
    
    print("\n--- Test 1: surface_equilibrium_coverage() ---")
    test_equilibrium_coverage_zero_pressure()
    test_equilibrium_coverage_langmuir_limit()
    test_equilibrium_coverage_low_limit()
    test_equilibrium_coverage_input_validation()
    
    print("\n--- Test 2: dissociation_flux() ---")
    test_dissociation_flux_zero_coverage()
    test_dissociation_flux_full_coverage()
    test_dissociation_flux_half_coverage()
    
    print("\n--- Test 3: calculate_damkohler_number() ---")
    test_damkohler_diffusion_limited()
    test_damkohler_dissociation_limited()
    test_damkohler_mixed_regime()
    
    print("\n--- Test 4: calculate_surface_limited_flux() ---")
    test_main_function_with_material()
    test_main_function_with_explicit_kinetics()
    test_main_function_fixed_coverage()
    
    print("\n--- Test 5: Comparison with Level 1 AND Level 4 ---")
    test_fast_kinetics_approaches_sieverts()
    test_fast_kinetics_approaches_level4()
    test_slow_kinetics_reduces_flux()
    
    print("\n--- Test 6: Input validation ---")
    test_input_validation_negative_pressure()
    test_input_validation_missing_material()
    
    print("\n" + "=" * 70)
    print("ALL TESTS PASSED!")
    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
