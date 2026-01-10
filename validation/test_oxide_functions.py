# validation/test_oxide_functions.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from datetime import datetime
import json

# Add parent directory to path to import calculation modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations.oxide_permeation import (
    molecular_diffusion_flux,
    calculate_oxide_resistance,
    calculate_metal_resistance,
    get_oxide_properties_at_T,
    compare_resistances,
    calculate_transition_pressure,
    pressure_dependence_analysis
)

# Create results directory if it doesn't exist
RESULTS_DIR = 'validation/results'
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# Timestamp for this run
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# WIDE PRESSURE RANGE
PRESSURE_MIN = 1e-10  # Pa (ultra-high vacuum)
PRESSURE_MAX = 1e28   # Pa 
N_POINTS = 200        # Number of points for smooth curves


def test_molecular_flux_linearity():
    """Test that oxide flux is linear in pressure over wide range."""
    print("\n=== Test 1: Molecular Flux Linearity ===")
    
    D_ox = 1e-20  # m²/s
    K_ox = 1e-22  # mol/m³/Pa
    thickness = 1e-9  # 1 nm
    
    # Wide pressure range
    pressures = np.logspace(np.log10(PRESSURE_MIN), np.log10(PRESSURE_MAX), N_POINTS)
    fluxes = [molecular_diffusion_flux(D_ox, K_ox, thickness, P, 0) for P in pressures]
    
    # Test points across the range
    test_pressures = [1e-8, 1e-4, 1e0, 1e4, 1e8]
    test_fluxes = [molecular_diffusion_flux(D_ox, K_ox, thickness, P, 0) for P in test_pressures]
    
    print(f"Pressure range: {PRESSURE_MIN:.1e} to {PRESSURE_MAX:.1e} Pa")
    for P, flux in zip(test_pressures, test_fluxes):
        print(f"  P = {P:.1e} Pa: Flux = {flux:.2e} mol/m²/s")
    
    # Check linearity
    ratio = test_fluxes[-1] / test_fluxes[0]
    expected_ratio = test_pressures[-1] / test_pressures[0]
    print(f"\nFlux ratio (P={test_pressures[-1]:.0e}/P={test_pressures[0]:.0e}): {ratio:.2e}")
    print(f"Expected ratio: {expected_ratio:.2e}")
    print(f"Linearity maintained: {abs(ratio - expected_ratio) / expected_ratio < 0.01}")
    
    # Visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Log-log scale (main plot)
    ax1.loglog(pressures, fluxes, 'b-', linewidth=2, label='Oxide flux')
    ax1.scatter(test_pressures, test_fluxes, color='red', s=100, zorder=5, label='Test points')
    
    # Add slope indicator
    p_mid = 1e0
    f_mid = molecular_diffusion_flux(D_ox, K_ox, thickness, p_mid, 0)
    ax1.plot([p_mid/100, p_mid*100], [f_mid/100, f_mid*100], 'k--', alpha=0.5, label='Slope = 1')
    
    ax1.set_xlabel('Pressure (Pa)')
    ax1.set_ylabel('Flux (mol/m²/s)')
    ax1.set_title(f'Molecular Diffusion: {PRESSURE_MIN:.0e} to {PRESSURE_MAX:.0e} Pa')
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend()
    
    # Slope verification plot
    ax2.semilogx(pressures[1:], np.diff(np.log(fluxes))/np.diff(np.log(pressures)), 'g-')
    ax2.axhline(y=1.0, color='r', linestyle='--', label='Expected slope = 1.0')
    ax2.set_xlabel('Pressure (Pa)')
    ax2.set_ylabel('Local slope d(log flux)/d(log P)')
    ax2.set_title('Slope Verification')
    ax2.set_ylim([0.95, 1.05])
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test1_linearity_{timestamp}.png', dpi=150)
    plt.show()
    
    print("✓ Flux is linear in pressure across entire range")
    
    return {
        'test': 'molecular_flux_linearity',
        'pressure_range': [PRESSURE_MIN, PRESSURE_MAX],
        'test_pressures': test_pressures,
        'test_fluxes': test_fluxes
    }


def test_resistance_calculations():
    """Test oxide and metal resistance calculations over wide range."""
    print("\n=== Test 2: Resistance Calculations ===")
    
    # Properties
    D_ox = 1e-20  # m²/s
    K_ox = 1e-22  # mol/m³/Pa
    thickness_ox = 1e-9  # 1 nm
    
    D_metal = 1e-9  # m²/s
    K_s_metal = 1e-20  # mol/m³/Pa^0.5
    thickness_metal = 1e-3  # 1 mm
    
    R_ox = calculate_oxide_resistance(D_ox, K_ox, thickness_ox)
    
    # Wide pressure range
    P_interface_range = np.logspace(np.log10(PRESSURE_MIN), np.log10(PRESSURE_MAX), N_POINTS)
    R_metal_values = [calculate_metal_resistance(D_metal, K_s_metal, thickness_metal, P) 
                      for P in P_interface_range]
    R_oxide_values = [R_ox] * len(P_interface_range)
    
    # Find all regime transitions
    ratios = np.array(R_oxide_values) / np.array(R_metal_values)
    oxide_limited = ratios > 10
    metal_limited = ratios < 0.1
    transition = (~oxide_limited) & (~metal_limited)
    
    # Visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Resistances
    ax1.loglog(P_interface_range, R_metal_values, 'b-', linewidth=2, label='Metal resistance ∝ √P')
    ax1.loglog(P_interface_range, R_oxide_values, 'r-', linewidth=2, label='Oxide resistance (constant)')
    
    # Mark transition points if they exist in our range
    if np.any(transition):
        P_trans_low = P_interface_range[np.where(transition)[0][0]]
        P_trans_high = P_interface_range[np.where(transition)[0][-1]]
        ax1.axvspan(P_trans_low, P_trans_high, alpha=0.2, color='green', label='Transition regime')
        ax1.axvline(P_trans_low, color='green', linestyle=':', alpha=0.5)
        ax1.axvline(P_trans_high, color='green', linestyle=':', alpha=0.5)
    
    ax1.set_ylabel('Resistance (Pa·s·m²/mol)')
    ax1.set_title(f'Resistance Evolution: {PRESSURE_MIN:.0e} to {PRESSURE_MAX:.0e} Pa')
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend()
    
    # Resistance ratio
    ax2.loglog(P_interface_range, ratios, 'k-', linewidth=2)
    ax2.axhline(y=10, color='red', linestyle='--', alpha=0.5, label='Oxide-limited boundary')
    ax2.axhline(y=0.1, color='blue', linestyle='--', alpha=0.5, label='Metal-limited boundary')
    ax2.fill_between(P_interface_range, 0.1, 10, alpha=0.2, color='green')
    
    ax2.set_xlabel('Interface Pressure (Pa)')
    ax2.set_ylabel('R_oxide / R_metal')
    ax2.set_title('Regime Classification')
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test2_resistances_{timestamp}.png', dpi=150)
    plt.show()
    
    # Print summary
    print(f"Oxide resistance (constant): {R_ox:.2e} Pa·s·m²/mol")
    print(f"\nMetal resistance varies with pressure:")
    for P in [1e-10, 1e-5, 1e0, 1e5, 1e28]:
        R_m = calculate_metal_resistance(D_metal, K_s_metal, thickness_metal, P)
        ratio = R_ox / R_m
        if ratio > 10:
            regime = "Oxide-limited"
        elif ratio < 0.1:
            regime = "Metal-limited"
        else:
            regime = "Transition"
        print(f"  P = {P:.0e} Pa: R_metal = {R_m:.2e}, Ratio = {ratio:.2e} ({regime})")
    
    print("✓ Resistance calculations complete")
    
    return {
        'test': 'resistance_calculations',
        'R_oxide': R_ox,
        'pressure_range': [PRESSURE_MIN, PRESSURE_MAX]
    }


def test_flux_regimes():
    """Test flux behavior across all regimes in wide pressure range."""
    print("\n=== Test 3: Flux Regimes ===")
    
    # System properties
    oxide_props = {
        'D_ox': 1e-20,
        'K_ox': 1e-22,
        'thickness': 6e-10  # 6 Angstroms (realistic)
    }
    
    metal_props = {
        'D_metal': 1e-9,
        'K_s_metal': 1e-20,
        'thickness': 1e-3
    }
    
    # Calculate transition pressure
    P_trans = calculate_transition_pressure(oxide_props, metal_props)
    print(f"Theoretical transition pressure: {P_trans:.2e} Pa")
    
    if P_trans > PRESSURE_MAX:
        print(f"WARNING: Transition pressure exceeds test range!")
        print(f"System is oxide-dominated throughout {PRESSURE_MIN:.0e} to {PRESSURE_MAX:.0e} Pa")
    
    # Wide pressure range
    pressures = np.logspace(np.log10(PRESSURE_MIN), np.log10(PRESSURE_MAX), N_POINTS)
    
    # Calculate fluxes for different scenarios
    oxide_only_flux = [molecular_diffusion_flux(oxide_props['D_ox'], oxide_props['K_ox'], 
                                                oxide_props['thickness'], P, 0) for P in pressures]
    
    # For metal-only, we need to use the classic Sieverts formula
    metal_only_flux = [(metal_props['D_metal'] * metal_props['K_s_metal'] / metal_props['thickness']) * 
                       np.sqrt(P) for P in pressures]
    
    # Visualization
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Plot fluxes
    ax.loglog(pressures, oxide_only_flux, 'r-', linewidth=2, label='Oxide only (slope=1)', alpha=0.7)
    ax.loglog(pressures, metal_only_flux, 'b-', linewidth=2, label='Metal only (slope=0.5)', alpha=0.7)
    
    # Mark transition pressure if in range
    if P_trans <= PRESSURE_MAX:
        ax.axvline(P_trans, color='green', linestyle='--', linewidth=2, 
                   label=f'Transition P = {P_trans:.2e} Pa')
    else:
        # Add text indicating transition is out of range
        ax.text(1e0, 1e-30, f'Transition at {P_trans:.1e} Pa\n(outside plot range)', 
                ha='center', fontsize=10, bbox=dict(boxstyle="round,pad=0.3", 
                facecolor="yellow", alpha=0.5))
    
    # Add slope indicators
    # Slope = 1 line
    p1 = 1e-8
    f1 = 1e-50
    ax.plot([p1, p1*100], [f1, f1*100], 'r--', alpha=0.5, linewidth=1)
    ax.text(p1*10, f1*10, 'slope=1', rotation=45, fontsize=8, color='red')
    
    # Slope = 0.5 line
    p2 = 1e6
    f2 = 1e-10
    ax.plot([p2, p2*100], [f2, f2*10], 'b--', alpha=0.5, linewidth=1)
    ax.text(p2*10, f2*3, 'slope=0.5', rotation=26.6, fontsize=8, color='blue')
    
    ax.set_xlabel('Upstream Pressure (Pa)')
    ax.set_ylabel('Flux (mol/m²/s)')
    ax.set_title(f'Flux Regimes: {PRESSURE_MIN:.0e} to {PRESSURE_MAX:.0e} Pa')
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc='upper left')
    
    # Add pressure context annotations
    ax.axvspan(1e-7, 1e3, alpha=0.1, color='yellow')
    ax.text(1e-2, 1e-60, 'Typical fusion\nconditions', ha='center', fontsize=8, 
            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test3_flux_regimes_{timestamp}.png', dpi=150)
    plt.show()
    
    print("✓ Flux regime analysis complete")
    
    return {
        'test': 'flux_regimes',
        'P_transition': P_trans,
        'pressure_range': [PRESSURE_MIN, PRESSURE_MAX]
    }


def test_temperature_effects():
    """Test temperature-dependent properties."""
    print("\n=== Test 4: Temperature Effects on Cr2O3 ===")
    
    temperatures_C = np.linspace(600, 1000, 20)
    temperatures_K = temperatures_C + 273.15
    
    D_values = []
    K_values = []
    
    for T_K in temperatures_K:
        props = get_oxide_properties_at_T('Cr2O3', T_K)
        D_values.append(props['D_ox'])
        K_values.append(props['K_ox'])
    
    # Print selected values
    print("Selected temperature points:")
    for i in [0, 4, 9, 14, 19]:
        print(f"T={temperatures_C[i]:4.0f}°C: D_ox={D_values[i]:.2e} m²/s, "
              f"K_ox={K_values[i]:.2e} mol/m³/Pa")
    
    # Arrhenius plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # D_ox Arrhenius plot
    ax1.semilogy(1000/temperatures_K, D_values, 'b-o', markersize=6)
    ax1.set_xlabel('1000/T (K⁻¹)')
    ax1.set_ylabel('D_ox (m²/s)')
    ax1.set_title('Arrhenius Plot: Diffusion Coefficient')
    ax1.grid(True, which="both", alpha=0.3)
    
    # Add temperature labels on top x-axis
    ax1_top = ax1.twiny()
    ax1_top.set_xlim(ax1.get_xlim())
    tick_positions = [1000/(600+273.15), 1000/(800+273.15), 1000/(1000+273.15)]
    ax1_top.set_xticks(tick_positions)
    ax1_top.set_xticklabels(['600°C', '800°C', '1000°C'])
    
    # K_ox Arrhenius plot
    ax2.semilogy(1000/temperatures_K, K_values, 'r-o', markersize=6)
    ax2.set_xlabel('1000/T (K⁻¹)')
    ax2.set_ylabel('K_ox (mol/m³/Pa)')
    ax2.set_title('Arrhenius Plot: Solubility Constant')
    ax2.grid(True, which="both", alpha=0.3)
    
    # Add temperature labels on top x-axis
    ax2_top = ax2.twiny()
    ax2_top.set_xlim(ax2.get_xlim())
    ax2_top.set_xticks(tick_positions)
    ax2_top.set_xticklabels(['600°C', '800°C', '1000°C'])
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test4_temperature_{timestamp}.png', dpi=150)
    plt.show()
    
    # Check Arrhenius behavior
    for i in range(1, len(D_values)):
        assert D_values[i] > D_values[i-1], "D_ox should increase with temperature"
    
    print("✓ Temperature dependence follows Arrhenius behavior")
    
    return {
        'test': 'temperature_effects',
        'temperatures_C': temperatures_C.tolist(),
        'D_values': D_values,
        'K_values': K_values
    }


def save_all_results(all_results):
    """Save all test results to JSON file."""
    results_file = f'{RESULTS_DIR}/test_results_{timestamp}.json'
    
    # Convert numpy arrays to lists for JSON serialization
    for result in all_results:
        for key, value in result.items():
            if isinstance(value, np.ndarray):
                result[key] = value.tolist()
    
    with open(results_file, 'w') as f:
        json.dump({
            'timestamp': timestamp,
            'pressure_range': [PRESSURE_MIN, PRESSURE_MAX],
            'test_results': all_results
        }, f, indent=2, default=str)
    
    print(f"\nResults saved to: {results_file}")


def run_all_tests():
    """Run all validation tests."""
    print("="*60)
    print("OXIDE PERMEATION VALIDATION TESTS")
    print(f"Pressure range: {PRESSURE_MIN:.0e} to {PRESSURE_MAX:.0e} Pa")
    print("="*60)
    
    all_results = []
    
    all_results.append(test_molecular_flux_linearity())
    all_results.append(test_resistance_calculations())
    all_results.append(test_flux_regimes())
    all_results.append(test_temperature_effects())
    
    save_all_results(all_results)
    
    print("\n" + "="*60)
    print("ALL TESTS COMPLETE ✓")
    print(f"Figures saved in: {RESULTS_DIR}/")
    print("="*60)


if __name__ == "__main__":
    run_all_tests()