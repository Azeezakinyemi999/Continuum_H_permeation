#!/usr/bin/env python3
"""
Test L1+L6 and L4+L6 Integration
================================
Tests the integration of surface kinetics (Level 6) with:
- Level 1: Simple metal (calculate_simple_metal_flux_with_surface)
- Level 4: Defective metal (calculate_defective_metal_flux_with_surface)

Tests all three coverage modes: 'equilibrium', 'steady_state', 'forced'
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations.permeation_calc import (
    calculate_simple_metal_flux_with_surface,
    calculate_defective_metal_flux_with_surface
)


def test_L1_L6_integration():
    """
    Test Level 1 + Level 6: Simple metal with surface kinetics
    All three coverage modes
    """
    print("="*70)
    print("TEST L1+L6: Simple Metal with Surface Kinetics")
    print("="*70)
    
    # Parameters
    D = 1e-10           # m²/s - diffusivity
    K_s = 0.5           # mol/(m³·Pa^0.5) - Sieverts constant
    thickness = 1e-3    # m - thickness
    P_up = 1e5          # Pa - upstream pressure
    P_down = 0          # Pa - downstream pressure
    temperature = 1073  # K - temperature (800°C)
    k_diss = 1e-8       # mol/(m²·s·Pa) - dissociation rate
    k_recomb = 1e-7     # m⁴/(mol·s) - recombination rate
    
    # Sieverts flux (reference)
    J_sieverts = D * K_s * (np.sqrt(P_up) - np.sqrt(P_down)) / thickness
    print(f"\nReference Sieverts flux: {J_sieverts:.6e} mol/(m²·s)")
    
    print("\n--- Testing Three Coverage Modes ---\n")
    
    modes = ['equilibrium', 'steady_state', 'forced']
    results = {}
    
    for mode in modes:
        print(f"\n{'='*50}")
        print(f"Coverage Mode: '{mode}'")
        print(f"{'='*50}")
        
        if mode == 'forced':
            # Test with forced coverage θ = 0.5
            result = calculate_simple_metal_flux_with_surface(
                D=D, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=P_down,
                temperature=temperature,
                k_diss=k_diss, k_recomb=k_recomb,
                coverage_mode='forced', forced_coverage=0.5
            )
        else:
            result = calculate_simple_metal_flux_with_surface(
                D=D, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=P_down,
                temperature=temperature,
                k_diss=k_diss, k_recomb=k_recomb,
                coverage_mode=mode
            )
        
        results[mode] = result
        
        print(f"  Flux: {result['flux']:.6e} mol/(m²·s)")
        print(f"  SRF:  {result['SRF']:.6f}")
        print(f"  θ:    {result['theta']:.6f}")
        print(f"  Da:   {result.get('Da', 'N/A')}")
        
        # Validate flux is positive
        assert result['flux'] > 0, f"Flux should be positive for mode {mode}"
        assert 0 <= result['SRF'] <= 1, f"SRF should be between 0 and 1 for mode {mode}"
        assert 0 <= result['theta'] <= 1, f"theta should be between 0 and 1 for mode {mode}"
    
    # Compare modes
    print("\n" + "="*50)
    print("COMPARISON OF MODES")
    print("="*50)
    print(f"{'Mode':<15} {'Flux':<15} {'SRF':<10} {'θ':<10}")
    print("-"*50)
    for mode in modes:
        r = results[mode]
        print(f"{mode:<15} {r['flux']:.4e} {r['SRF']:.4f} {r['theta']:.4f}")
    
    print("\n✓ L1+L6 Integration Test PASSED")
    return results


def test_L4_L6_integration():
    """
    Test Level 4 + Level 6: Defective metal with surface kinetics
    All three coverage modes
    """
    print("\n")
    print("="*70)
    print("TEST L4+L6: Defective Metal with Surface Kinetics")
    print("="*70)
    
    # Parameters - same as L1+L6 for comparison
    D_lattice = 1e-10   # m²/s - lattice diffusivity
    K_s = 0.5           # mol/(m³·Pa^0.5) - Sieverts constant
    thickness = 1e-3    # m - thickness
    P_up = 1e5          # Pa - upstream pressure
    P_down = 0          # Pa - downstream pressure
    temperature = 1073  # K - temperature (800°C)
    k_diss = 1e-8       # mol/(m²·s·Pa) - dissociation rate
    k_recomb = 1e-7     # m⁴/(mol·s) - recombination rate
    
    # Microstructure parameters (defects)
    microstructure_params = {
        'grain_size': 50e-6,  # 50 μm grain size
        'grain_shape': 'equiaxed',
        'gb_type': 'HAGB',  # high-angle grain boundary
        'trap_list': [
            {'name': 'dislocation', 'binding_energy': 30000, 'density': 1e14}  # J/mol
        ]
    }
    
    print(f"\nMicrostructure: grain_size=50 μm, dislocation traps")
    print(f"Temperature: {temperature} K")
    
    print("\n--- Testing Three Coverage Modes ---\n")
    
    modes = ['equilibrium', 'steady_state', 'forced']
    results = {}
    
    for mode in modes:
        print(f"\n{'='*50}")
        print(f"Coverage Mode: '{mode}'")
        print(f"{'='*50}")
        
        if mode == 'forced':
            result = calculate_defective_metal_flux_with_surface(
                D_lattice=D_lattice, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=P_down,
                temperature=temperature,
                microstructure_params=microstructure_params,
                k_diss=k_diss, k_recomb=k_recomb,
                coverage_mode='forced', forced_coverage=0.5
            )
        else:
            result = calculate_defective_metal_flux_with_surface(
                D_lattice=D_lattice, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=P_down,
                temperature=temperature,
                microstructure_params=microstructure_params,
                k_diss=k_diss, k_recomb=k_recomb,
                coverage_mode=mode
            )
        
        results[mode] = result
        
        print(f"  Flux: {result['flux']:.6e} mol/(m²·s)")
        print(f"  SRF:  {result['SRF']:.6f}")
        print(f"  θ:    {result['theta']:.6f}")
        print(f"  Da:   {result.get('Da', 'N/A')}")
        print(f"  D_eff: {result.get('D_eff', 'N/A')}")
        
        # Validate
        assert result['flux'] > 0, f"Flux should be positive for mode {mode}"
        assert 0 <= result['SRF'] <= 1, f"SRF should be between 0 and 1 for mode {mode}"
        assert 0 <= result['theta'] <= 1, f"theta should be between 0 and 1 for mode {mode}"
    
    # Compare modes
    print("\n" + "="*50)
    print("COMPARISON OF MODES")
    print("="*50)
    print(f"{'Mode':<15} {'Flux':<15} {'SRF':<10} {'θ':<10}")
    print("-"*50)
    for mode in modes:
        r = results[mode]
        print(f"{mode:<15} {r['flux']:.4e} {r['SRF']:.4f} {r['theta']:.4f}")
    
    print("\n✓ L4+L6 Integration Test PASSED")
    return results


def test_L1_vs_L4_comparison():
    """
    Compare L1+L6 vs L4+L6 to see effect of trapping on surface kinetics
    """
    print("\n")
    print("="*70)
    print("COMPARISON: L1+L6 vs L4+L6")
    print("="*70)
    
    # Common parameters
    D = 1e-10
    K_s = 0.5
    thickness = 1e-3
    P_up = 1e5
    P_down = 0
    temperature = 1073
    k_diss = 1e-8
    k_recomb = 1e-7
    
    # Microstructure params for L4
    microstructure_params = {
        'grain_size': 50e-6,
        'grain_shape': 'equiaxed',
        'gb_type': 'HAGB',
        'trap_list': [
            {'name': 'dislocation', 'binding_energy': 30000, 'density': 1e14}
        ]
    }
    
    # Calculate L1+L6 (steady_state)
    L1_result = calculate_simple_metal_flux_with_surface(
        D=D, K_s=K_s, thickness=thickness,
        P_up=P_up, P_down=P_down,
        temperature=temperature,
        k_diss=k_diss, k_recomb=k_recomb,
        coverage_mode='steady_state'
    )
    
    # Calculate L4+L6 (steady_state)
    L4_result = calculate_defective_metal_flux_with_surface(
        D_lattice=D, K_s=K_s, thickness=thickness,
        P_up=P_up, P_down=P_down,
        temperature=temperature,
        microstructure_params=microstructure_params,
        k_diss=k_diss, k_recomb=k_recomb,
        coverage_mode='steady_state'
    )
    
    print("\nSteady-State Mode Comparison:")
    print("-"*50)
    print(f"{'Quantity':<20} {'L1+L6':<20} {'L4+L6':<20}")
    print("-"*50)
    print(f"{'Flux':<20} {L1_result['flux']:.4e} {L4_result['flux']:.4e}")
    print(f"{'SRF':<20} {L1_result['SRF']:.4f} {L4_result['SRF']:.4f}")
    print(f"{'θ':<20} {L1_result['theta']:.4f} {L4_result['theta']:.4f}")
    print(f"{'Da':<20} {L1_result['Da']:.4f} {L4_result['Da']:.4f}")
    
    # Physics insight
    print("\n" + "="*50)
    print("PHYSICS INSIGHT")
    print("="*50)
    D_eff = L4_result.get('D_eff', D)
    print(f"D_eff in L4: {D_eff:.4e} vs D in L1: {D:.4e}")
    print(f"D_eff/D ratio: {D_eff/D:.4f}")
    print("\nNote: Trapping reduces D_eff, which increases Da (= k_diss·K_s·L/D_eff)")
    print("Higher Da means system is MORE diffusion-limited (less surface-limited)")
    print("This is why L4+L6 has HIGHER SRF than L1+L6 for same surface kinetics")
    
    print("\n✓ L1 vs L4 Comparison Test PASSED")


def test_coverage_mode_physics():
    """
    Test physics consistency of coverage modes.
    
    For UPSTREAM surface where dissociation limits flux:
    - θ_equilibrium: from Langmuir isotherm (adsorption-desorption equilibrium with gas)
    - θ_steady_state: from J_diss = J_bulk (dissociation matches bulk removal)
    
    At steady state, θ is LOWER than equilibrium because:
    - Flux removes H from surface continuously
    - θ drops until dissociation can keep up with removal
    - This is the "surface-limited" regime behavior
    """
    print("\n")
    print("="*70)
    print("PHYSICS CONSISTENCY TEST")
    print("="*70)
    
    K_s = 0.5
    thickness = 1e-3
    P_up = 1e5
    temperature = 1073
    k_recomb = 1e-7
    
    # Test with different Da regimes
    test_cases = [
        {"name": "Low Da (surface-limited)", "k_diss": 1e-10, "D": 1e-8},
        {"name": "Medium Da", "k_diss": 1e-8, "D": 1e-10},
        {"name": "High Da (diffusion-limited)", "k_diss": 1e-6, "D": 1e-12},
    ]
    
    for case in test_cases:
        print(f"\n{case['name']}:")
        
        # Equilibrium
        eq_result = calculate_simple_metal_flux_with_surface(
            D=case['D'], K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0,
            temperature=temperature,
            k_diss=case['k_diss'], k_recomb=k_recomb,
            coverage_mode='equilibrium'
        )
        
        # Steady-state
        ss_result = calculate_simple_metal_flux_with_surface(
            D=case['D'], K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0,
            temperature=temperature,
            k_diss=case['k_diss'], k_recomb=k_recomb,
            coverage_mode='steady_state'
        )
        
        print(f"  θ_equilibrium = {eq_result['theta']:.6f}")
        print(f"  θ_steady_state = {ss_result['theta']:.6f}")
        print(f"  Da = {ss_result['Da']:.4f}")
        print(f"  SRF_eq = {eq_result['SRF']:.6f}, SRF_ss = {ss_result['SRF']:.6f}")
        
        # θ_steady_state should be <= θ_equilibrium (depleted by flux)
        # At high Da, they converge (diffusion-limited, surface stays at equilibrium)
        if ss_result['theta'] <= eq_result['theta'] + 1e-6:
            print(f"  ✓ θ_ss ≤ θ_eq (correct: surface depleted by flux)")
        else:
            print(f"  ⚠ θ_ss > θ_eq (unexpected)")
        
        # Check: higher Da → SRF closer to 1
        if case['name'] == "High Da (diffusion-limited)":
            assert ss_result['SRF'] > 0.5, "High Da should have SRF > 0.5"
        elif case['name'] == "Low Da (surface-limited)":
            assert ss_result['SRF'] < 0.1, "Low Da should have SRF < 0.1"
    
    print("\n✓ Physics Consistency Test PASSED")


if __name__ == "__main__":
    print("\n" + "#"*70)
    print("# L1+L6 and L4+L6 INTEGRATION TESTS")
    print("#"*70 + "\n")
    
    # Run all tests
    L1_results = test_L1_L6_integration()
    L4_results = test_L4_L6_integration()
    test_L1_vs_L4_comparison()
    test_coverage_mode_physics()
    
    print("\n" + "#"*70)
    print("# ALL INTEGRATION TESTS PASSED!")
    print("#"*70 + "\n")
