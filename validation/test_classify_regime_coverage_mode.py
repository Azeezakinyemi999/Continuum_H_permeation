#!/usr/bin/env python3
"""
Test classify_regime.py updates for coverage_mode support.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations.classify_regime import (
    classify_regime_level6, classify_regime_level16, 
    classify_regime_level46, classify_regime_level56
)


def test_level6_coverage_modes():
    """Test classify_regime_level6 with all coverage modes."""
    print("="*60)
    print("Testing classify_regime_level6 with coverage_mode")
    print("="*60)
    
    # Test equilibrium mode
    result = classify_regime_level6(Da=0.05, theta=0.95, coverage_mode='equilibrium')
    print(f"\nEquilibrium mode:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  coverage_mode: {result['coverage_mode']}")
    print(f"  coverage_mode_info: {result['coverage_mode_info']}")
    assert result['coverage_mode'] == 'equilibrium'
    assert 'equilibrium' in result['regime_hierarchy']
    
    # Test steady_state mode with theta_equilibrium
    result = classify_regime_level6(
        Da=0.05, theta=0.85, 
        coverage_mode='steady_state', 
        theta_equilibrium=0.95
    )
    print(f"\nSteady-state mode:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  coverage_mode: {result['coverage_mode']}")
    print(f"  theta: {result['theta']}, theta_eq: {result['theta_equilibrium']}")
    print(f"  coverage_depletion: {result['coverage_depletion']:.1%}")
    assert result['coverage_mode'] == 'steady_state'
    assert result['coverage_depletion'] > 0  # Should be depleted
    
    # Test forced mode
    result = classify_regime_level6(Da=0.05, theta=0.5, coverage_mode='forced')
    print(f"\nForced mode:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  coverage_mode_info: {result['coverage_mode_info']}")
    assert result['coverage_mode'] == 'forced'
    assert 'forced' in result['regime_hierarchy']
    
    print("\n✓ Level 6 coverage mode tests PASSED")


def test_level16_coverage_modes():
    """Test classify_regime_level16 with coverage_mode."""
    print("\n" + "="*60)
    print("Testing classify_regime_level16")
    print("="*60)
    
    # Equilibrium mode (default)
    result = classify_regime_level16(Da=0.05, theta=0.95, SRF=0.2)
    print(f"\nL1+L6 equilibrium (default):")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    assert 'steady_state' not in result['regime_hierarchy']
    
    # Steady-state mode
    result = classify_regime_level16(
        Da=0.05, theta=0.85, SRF=0.2, 
        coverage_mode='steady_state', 
        theta_equilibrium=0.95
    )
    print(f"\nL1+L6 steady_state:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  coverage_mode: {result['coverage_mode']}")
    print(f"  coverage_depletion: {result['coverage_depletion']:.1%}")
    assert result['coverage_mode'] == 'steady_state'
    assert 'steady_state' in result['regime_hierarchy']
    
    print("\n✓ Level 1+6 coverage mode tests PASSED")


def test_level46_coverage_modes():
    """Test classify_regime_level46 with coverage_mode."""
    print("\n" + "="*60)
    print("Testing classify_regime_level46")
    print("="*60)
    
    # Equilibrium mode
    result = classify_regime_level46(
        Da=0.05, theta=0.95, SRF=0.2, 
        modification_factor=0.8,
        coverage_mode='equilibrium'
    )
    print(f"\nL4+L6 equilibrium:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  dominant_limitation: {result['dominant_limitation']}")
    
    # Steady-state mode
    result = classify_regime_level46(
        Da=0.05, theta=0.85, SRF=0.2, 
        modification_factor=0.8,
        coverage_mode='steady_state', 
        theta_equilibrium=0.95
    )
    print(f"\nL4+L6 steady_state:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  dominant_limitation: {result['dominant_limitation']}")
    print(f"  coverage_mode: {result['coverage_mode']}")
    print(f"  coverage_depletion: {result['coverage_depletion']:.1%}")
    assert result['coverage_mode'] == 'steady_state'
    assert 'steady_state' in result['regime_hierarchy']
    
    print("\n✓ Level 4+6 coverage mode tests PASSED")


def test_level56_coverage_modes():
    """Test classify_regime_level56 with coverage_mode."""
    print("\n" + "="*60)
    print("Testing classify_regime_level56")
    print("="*60)
    
    # Metal-limited with steady-state coverage
    result = classify_regime_level56(
        'metal_limited', 
        Da=0.05, theta=0.85, SRF=0.2, 
        modification_factor=0.8, 
        coverage_mode='steady_state',
        theta_equilibrium=0.95
    )
    print(f"\nL5+L6 metal_limited, steady_state:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  dominant_limitation: {result['dominant_limitation']}")
    print(f"  coverage_mode: {result['coverage_mode']}")
    assert result['coverage_mode'] == 'steady_state'
    
    # Oxide-limited with equilibrium
    result = classify_regime_level56(
        'oxide_limited', 
        Da=5.0, theta=0.99, SRF=0.9, 
        modification_factor=0.8,
        flux_intact_contribution=1e-6,
        flux_defect_contribution=5e-6,
        coverage_mode='equilibrium'
    )
    print(f"\nL5+L6 oxide_limited, equilibrium:")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    print(f"  dominant_limitation: {result['dominant_limitation']}")
    
    print("\n✓ Level 5+6 coverage mode tests PASSED")


def test_regime_hierarchy_format():
    """Test that regime hierarchy includes coverage_mode correctly."""
    print("\n" + "="*60)
    print("Testing regime_hierarchy format")
    print("="*60)
    
    # Low Da (surface-significant) with steady_state
    result = classify_regime_level16(
        Da=0.01, theta=0.5, SRF=0.1, 
        coverage_mode='steady_state'
    )
    print(f"\nLow Da (surface-limited), steady_state:")
    print(f"  surface_significant: {result['surface_significant']}")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    assert 'dissociation-limited' in result['regime_hierarchy']
    assert 'steady_state' in result['regime_hierarchy']
    
    # High Da (diffusion-limited) - surface not significant
    result = classify_regime_level16(
        Da=100, theta=0.99, SRF=0.95, 
        coverage_mode='steady_state'
    )
    print(f"\nHigh Da (diffusion-limited), steady_state:")
    print(f"  surface_significant: {result['surface_significant']}")
    print(f"  regime_hierarchy: {result['regime_hierarchy']}")
    assert 'sieverts_equilibrium' in result['regime_hierarchy']
    
    print("\n✓ Regime hierarchy format tests PASSED")


if __name__ == "__main__":
    print("\n" + "#"*60)
    print("# CLASSIFY_REGIME COVERAGE_MODE TESTS")
    print("#"*60)
    
    test_level6_coverage_modes()
    test_level16_coverage_modes()
    test_level46_coverage_modes()
    test_level56_coverage_modes()
    test_regime_hierarchy_format()
    
    print("\n" + "#"*60)
    print("# ALL CLASSIFY_REGIME TESTS PASSED!")
    print("#"*60 + "\n")
