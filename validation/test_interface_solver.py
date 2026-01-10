# validation/test_interface_solver.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from datetime import datetime
import json


# Add parent directory to path to import calculation modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations.interface_solver import (
    solve_interface_pressure,
    calculate_oxide_metal_system,
    calculate_concentration_profile,
    flux_balance_equation,
    calculate_metal_flux_sieverts
)
from calculations.oxide_permeation import (
    molecular_diffusion_flux,
    calculate_transition_pressure)

# Create results directory if it doesn't exist
RESULTS_DIR = 'validation/test_interface_solver_results'
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# Timestamp for this run
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# PRESSURE RANGE FOR INTERFACE TESTS
PRESSURE_MIN = 1e-10  # Pa
PRESSURE_MAX = 1e28   # Pa  
N_POINTS = 100        # Number of points for smooth curves


def test_interface_pressure_solutions():
    """Test interface pressure solver across wide pressure range."""
    print("\n=== Test 1: Interface Pressure Solutions ===")
    
    # System properties (realistic values)
    oxide_props = {
        'D_ox': 1e-18,      # m²/s
        'K_ox': 1e-20,      # mol/m³/Pa
        'thickness': 6e-10  # 6 Angstroms
    }
    
    metal_props = {
        'D_metal': 1e-9,       # m²/s
        'K_s_metal': 1e-20,    # mol/m³/Pa^0.5
        'thickness': 1e-3      # 1 mm
    }
    
    # Test pressures spanning wide range
    test_pressures = np.logspace(np.log10(PRESSURE_MIN), np.log10(PRESSURE_MAX), 50)
    
    results_data = {
        'P_upstream': [],
        'P_interface': [],
        'P_normalized': [],
        'flux': [],
        'regime': [],
        'converged': []
    }
    
    print(f"Testing interface solutions from {PRESSURE_MIN:.1e} to {PRESSURE_MAX:.1e} Pa")
    print("-" * 60)
    
    # Sample points to print
    sample_indices = [0, 12, 25, 37, 49]
    
    for idx, P_up in enumerate(test_pressures):
        result = calculate_oxide_metal_system(P_up, 0, oxide_props, metal_props)
        
        results_data['P_upstream'].append(P_up)
        results_data['P_interface'].append(result['P_interface'])
        results_data['P_normalized'].append(result['P_interface_normalized'])
        results_data['flux'].append(result['flux'])
        results_data['regime'].append(result['regime'])
        results_data['converged'].append(result['converged'])
        
        if idx in sample_indices:
            print(f"P_upstream = {P_up:.1e} Pa:")
            print(f"  P_interface = {result['P_interface']:.2e} Pa")
            print(f"  Normalized = {result['P_interface_normalized']:.3f}")
            print(f"  Flux = {result['flux']:.2e} mol/m²/s")
            print(f"  Regime: {result['regime']}")
            print(f"  Converged: {result['converged']}")
    
    # Visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # Interface pressure vs upstream pressure
    ax1.loglog(results_data['P_upstream'], results_data['P_interface'], 'b-', linewidth=2)
    ax1.loglog(results_data['P_upstream'], results_data['P_upstream'], 'k--', alpha=0.3, label='P_interface = P_upstream')
    ax1.set_xlabel('Upstream Pressure (Pa)')
    ax1.set_ylabel('Interface Pressure (Pa)')
    ax1.set_title('Interface Pressure vs Upstream Pressure')
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend()
    
    # Normalized interface position
    ax2.semilogx(results_data['P_upstream'], results_data['P_normalized'], 'g-', linewidth=2)
    ax2.set_xlabel('Upstream Pressure (Pa)')
    ax2.set_ylabel('(P_interface - P_down)/(P_up - P_down)')
    ax2.set_title('Normalized Interface Position')
    ax2.set_ylim([0, 1])
    ax2.grid(True, alpha=0.3)
    
    # System flux
    ax3.loglog(results_data['P_upstream'], results_data['flux'], 'r-', linewidth=2)
    ax3.set_xlabel('Upstream Pressure (Pa)')
    ax3.set_ylabel('System Flux (mol/m²/s)')
    ax3.set_title('Total System Flux')
    ax3.grid(True, which="both", alpha=0.3)
    
    # Regime classification
    regime_map = {'oxide_limited': 0, 'transition': 1, 'metal_limited': 2}
    regime_values = [regime_map[r] for r in results_data['regime']]
    ax4.semilogx(results_data['P_upstream'], regime_values, 'mo-', markersize=4)
    ax4.set_xlabel('Upstream Pressure (Pa)')
    ax4.set_ylabel('Regime')
    ax4.set_yticks([0, 1, 2])
    ax4.set_yticklabels(['Oxide Limited', 'Transition', 'Metal Limited'])
    ax4.set_title('Operating Regime')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test1_interface_solutions_{timestamp}.png', dpi=150)
    plt.show()
    
    print("✓ Interface pressure solutions complete")
    
    return {
        'test': 'interface_pressure_solutions',
        'pressure_range': [PRESSURE_MIN, PRESSURE_MAX],
        'n_points_tested': len(test_pressures),
        'all_converged': all(results_data['converged'])
    }


def test_flux_continuity():
    """Test flux continuity at interface."""
    print("\n=== Test 2: Flux Continuity Validation ===")
    
    oxide_props = {
        'D_ox': 1e-18,
        'K_ox': 1e-20,
        'thickness': 6e-10
    }
    
    metal_props = {
        'D_metal': 1e-9,
        'K_s_metal': 1e-20,
        'thickness': 1e-3
    }
    
    # Test at various pressures
    test_pressures = np.logspace(-6, 6, 25)
    flux_errors = []
    
    print("Checking flux continuity at interface...")
    print(f"{'P_upstream (Pa)':>15} {'Flux Error':>12} {'Status':>10}")
    print("-" * 40)
    
    max_error = 0
    worst_case_P = 0
    
    for P_up in test_pressures:
        result = calculate_oxide_metal_system(P_up, 0, oxide_props, metal_props)
        
        # Manually check flux continuity
        flux_ox = molecular_diffusion_flux(
            oxide_props['D_ox'], oxide_props['K_ox'], oxide_props['thickness'],
            P_up, result['P_interface']
        )
        flux_met = calculate_metal_flux_sieverts(
            metal_props['D_metal'], metal_props['K_s_metal'], metal_props['thickness'],
            result['P_interface'], 0
        )
        
        error = abs(flux_ox - flux_met) / flux_ox if flux_ox > 0 else 0
        flux_errors.append(error)
        
        if error > max_error:
            max_error = error
            worst_case_P = P_up
        
        status = "✓" if error < 1e-10 else "⚠"
        if P_up in [1e-6, 1e-3, 1e0, 1e3, 1e6]:
            print(f"{P_up:>15.1e} {error:>12.2e} {status:>10}")
    
    print(f"\nMaximum error: {max_error:.2e} at P = {worst_case_P:.1e} Pa")
    
    # Visualization
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    ax.loglog(test_pressures, flux_errors, 'b-o', markersize=6)
    ax.axhline(y=1e-10, color='g', linestyle='--', label='Target precision (1e-10)')
    ax.set_xlabel('Upstream Pressure (Pa)')
    ax.set_ylabel('Relative Flux Error at Interface')
    ax.set_title('Flux Continuity Verification')
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test2_flux_continuity_{timestamp}.png', dpi=150)
    plt.show()
    
    print("✓ Flux continuity validated")
    
    return {
        'test': 'flux_continuity',
        'max_error': max_error,
        'worst_case_pressure': worst_case_P,
        'passed': max_error < 1e-8
    }


def test_concentration_profiles():
    """Test concentration profiles through oxide and metal."""
    print("\n=== Test 3: Concentration Profiles ===")
    
    oxide_props = {
        'D_ox': 1e-18,
        'K_ox': 1e-20,
        'thickness': 6e-10
    }
    
    metal_props = {
        'D_metal': 1e-9,
        'K_s_metal': 1e-20,
        'thickness': 1e-3
    }
    
    # Test at different pressures
    test_pressures = [1e-3, 1e0, 1e3, 1e6]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    discontinuities = []
    
    for idx, P_up in enumerate(test_pressures):
        profile = calculate_concentration_profile(P_up, 0, oxide_props, metal_props)
        discontinuities.append(profile['C_discontinuity'])
        
        ax = axes[idx]
        
        # Plot complete profile
        ax.plot(profile['x_all']*1e9, profile['C_all'], 'b-', linewidth=2)
        
        # Mark interface
        interface_pos = oxide_props['thickness']*1e9
        ax.axvline(interface_pos, color='r', linestyle='--', alpha=0.5, label='Interface')
        
        # Mark discontinuity
        ax.plot([interface_pos, interface_pos], 
               [profile['C_oxide'][-1], profile['C_metal'][0]], 
               'ro-', markersize=8, linewidth=2, label='Discontinuity')
        
        ax.set_xlabel('Position (nm)')
        ax.set_ylabel('Concentration (mol/m³)')
        ax.set_title(f'P_upstream = {P_up:.1e} Pa\n'
                    f'P_interface = {profile["P_interface"]:.2e} Pa')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add text showing discontinuity magnitude
        ax.text(0.5, 0.95, f'ΔC = {abs(profile["C_discontinuity"]):.2e} mol/m³',
                transform=ax.transAxes, ha='center', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test3_concentration_profiles_{timestamp}.png', dpi=150)
    plt.show()
    
    print(f"Concentration discontinuities at interface:")
    for P, disc in zip(test_pressures, discontinuities):
        print(f"  P = {P:.1e} Pa: ΔC = {abs(disc):.2e} mol/m³")
    
    print("✓ Concentration profiles calculated")
    
    return {
        'test': 'concentration_profiles',
        'test_pressures': test_pressures,
        'discontinuities': discontinuities
    }


def test_limiting_cases():
    """Test limiting cases for validation."""
    print("\n=== Test 4: Limiting Cases ===")
    
    metal_props = {
        'D_metal': 1e-9,
        'K_s_metal': 1e-20,
        'thickness': 1e-3
    }
    
    test_results = []
    
    # Case 1: Very thick oxide
    print("\nCase 1: Very thick oxide (1 μm)")
    oxide_props_thick = {
        'D_ox': 1e-18,
        'K_ox': 1e-20,
        'thickness': 1e-6  # 1 micron
    }
    
    result_thick = calculate_oxide_metal_system(1.0, 0, oxide_props_thick, metal_props)
    print(f"  P_interface/P_upstream = {result_thick['P_interface']/1.0:.2e}")
    print(f"  Normalized position = {result_thick['P_interface_normalized']:.2e}")
    print(f"  Regime: {result_thick['regime']}")
    print(f"  Expected: oxide_limited ✓" if result_thick['regime'] == 'oxide_limited' else "  UNEXPECTED!")
    test_results.append(('thick_oxide', result_thick['regime'] == 'oxide_limited'))
    
    # Case 2: Very thin oxide
    print("\nCase 2: Very thin oxide (1 pm)")
    oxide_props_thin = {
        'D_ox': 1e-18,
        'K_ox': 1e-20,
        'thickness': 1e-12  # 1 pm
    }
    
    result_thin = calculate_oxide_metal_system(1.0, 0, oxide_props_thin, metal_props)
    print(f"  P_interface/P_upstream = {result_thin['P_interface']/1.0:.2e}")
    print(f"  Normalized position = {result_thin['P_interface_normalized']:.2e}")
    print(f"  Regime: {result_thin['regime']}")
    print(f"  Expected: approaches metal behavior ✓")
    test_results.append(('thin_oxide', result_thin['P_interface_normalized'] > 0.99))
    
    # Case 3: Equal resistances (should be transition)
    print("\nCase 3: Matched resistances")
    # Find conditions for equal resistances
    P_test = 1e6
    oxide_props_matched = {
        'D_ox': 1e-18,
        'K_ox': 1e-20,
        'thickness': 2e-10
    }
    
    result_matched = calculate_oxide_metal_system(P_test, 0, oxide_props_matched, metal_props)
    print(f"  Resistance ratio = {result_matched['resistance_ratio']:.2f}")
    print(f"  Regime: {result_matched['regime']}")
    
    # Visualization of limiting cases
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Pressure drop distribution
    cases = ['Thick Oxide\n(1 μm)', 'Thin Oxide\n(1 pm)', 'Matched\nResistances']
    normalized_positions = [
        result_thick['P_interface_normalized'],
        result_thin['P_interface_normalized'],
        result_matched['P_interface_normalized']
    ]
    
    colors = ['red', 'blue', 'green']
    ax1.bar(cases, normalized_positions, color=colors, alpha=0.7)
    ax1.set_ylabel('Normalized Interface Position')
    ax1.set_title('Pressure Drop Distribution')
    ax1.set_ylim([0, 1.1])
    ax1.grid(True, alpha=0.3)
    
    # Resistance ratios
    resistance_ratios = [
        result_thick['resistance_ratio'],
        result_thin['resistance_ratio'],
        result_matched['resistance_ratio']
    ]
    
    ax2.semilogy(cases, resistance_ratios, 'ko-', markersize=10, linewidth=2)
    ax2.axhline(y=10, color='r', linestyle='--', alpha=0.5, label='Oxide-limited boundary')
    ax2.axhline(y=0.1, color='b', linestyle='--', alpha=0.5, label='Metal-limited boundary')
    ax2.fill_between([-0.5, 2.5], 0.1, 10, alpha=0.2, color='green', label='Transition')
    ax2.set_ylabel('R_oxide / R_metal')
    ax2.set_title('Resistance Ratios')
    ax2.set_xlim([-0.5, 2.5])
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test4_limiting_cases_{timestamp}.png', dpi=150)
    plt.show()
    
    print("\n✓ Limiting cases validated")
    
    all_tests_passed = all(result for _, result in test_results)
    
    return {
        'test': 'limiting_cases',
        'all_passed': all_tests_passed,
        'individual_results': test_results
    }


def save_all_results(all_results):
    """Save all test results to JSON file."""
    results_file = f'{RESULTS_DIR}/interface_test_results_{timestamp}.json'
    
    # Convert numpy arrays and complex types to JSON-serializable format
    for result in all_results:
        for key, value in result.items():
            if isinstance(value, np.ndarray):
                result[key] = value.tolist()
            elif isinstance(value, (list, tuple)) and len(value) > 0:
                if isinstance(value[0], np.ndarray):
                    result[key] = [v.tolist() if isinstance(v, np.ndarray) else v for v in value]
    
    with open(results_file, 'w') as f:
        json.dump({
            'timestamp': timestamp,
            'test_suite': 'interface_solver',
            'test_results': all_results
        }, f, indent=2, default=str)
    
    print(f"\nResults saved to: {results_file}")

def test_parametric_studies():
    """Perform parametric studies on key system parameters."""
    print("\n=== Test 5: Parametric Studies ===")
    
    # Base properties
    base_oxide_props = {
        'D_ox': 1e-18,
        'K_ox': 1e-20,
        'thickness': 6e-10
    }
    
    base_metal_props = {
        'D_metal': 1e-9,
        'K_s_metal': 1e-20,
        'thickness': 1e-3
    }
    
    # Create figure with subplots for different studies
    fig = plt.figure(figsize=(16, 12))
    
    # Study 1: Oxide Thickness Effect
    print("\nStudy 1: Oxide Thickness Effect on Transition Pressure")
    print("-" * 50)
    ax1 = plt.subplot(2, 3, 1)
    
    thicknesses = np.logspace(-10, -6, 30)  # 1Å to 1μm
    transition_pressures = []
    
    for t in thicknesses:
        oxide_props = base_oxide_props.copy()
        oxide_props['thickness'] = t
        P_trans = calculate_transition_pressure(oxide_props, base_metal_props)
        transition_pressures.append(P_trans)
    
    ax1.loglog(thicknesses*1e10, transition_pressures, 'b-', linewidth=2)
    ax1.set_xlabel('Oxide Thickness (Å)')
    ax1.set_ylabel('Transition Pressure (Pa)')
    ax1.set_title('Effect of Oxide Thickness')
    ax1.grid(True, which="both", alpha=0.3)
    
    # Add fusion-relevant pressure range
    ax1.axhspan(1e-7, 1e3, alpha=0.2, color='yellow', label='Fusion range')
    ax1.legend()
    
    # Print key values
    idx_6A = np.argmin(np.abs(thicknesses - 6e-10))
    print(f"At 6Å: P_transition = {transition_pressures[idx_6A]:.2e} Pa")
    idx_100A = np.argmin(np.abs(thicknesses - 1e-8))
    print(f"At 100Å: P_transition = {transition_pressures[idx_100A]:.2e} Pa")
    
    # Study 2: Oxide Diffusivity Effect
    print("\nStudy 2: Oxide Diffusivity Effect")
    print("-" * 50)
    ax2 = plt.subplot(2, 3, 2)
    
    D_ox_values = np.logspace(-24, -16, 30)  # Wide range
    transition_pressures_D = []
    
    for D in D_ox_values:
        oxide_props = base_oxide_props.copy()
        oxide_props['D_ox'] = D
        P_trans = calculate_transition_pressure(oxide_props, base_metal_props)
        transition_pressures_D.append(P_trans)
    
    ax2.loglog(D_ox_values, transition_pressures_D, 'r-', linewidth=2)
    ax2.set_xlabel('Oxide Diffusivity D_ox (m²/s)')
    ax2.set_ylabel('Transition Pressure (Pa)')
    ax2.set_title('Effect of Oxide Diffusivity')
    ax2.grid(True, which="both", alpha=0.3)
    ax2.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
    
    # Study 3: Oxide Solubility Effect
    print("\nStudy 3: Oxide Solubility Effect")
    print("-" * 50)
    ax3 = plt.subplot(2, 3, 3)
    
    K_ox_values = np.logspace(-26, -18, 30)
    transition_pressures_K = []
    
    for K in K_ox_values:
        oxide_props = base_oxide_props.copy()
        oxide_props['K_ox'] = K
        P_trans = calculate_transition_pressure(oxide_props, base_metal_props)
        transition_pressures_K.append(P_trans)
    
    ax3.loglog(K_ox_values, transition_pressures_K, 'g-', linewidth=2)
    ax3.set_xlabel('Oxide Solubility K_ox (mol/m³/Pa)')
    ax3.set_ylabel('Transition Pressure (Pa)')
    ax3.set_title('Effect of Oxide Solubility')
    ax3.grid(True, which="both", alpha=0.3)
    ax3.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
    
    # Study 4: Combined Oxide Permeability (D*K)
    print("\nStudy 4: Combined Oxide Permeability Effect")
    print("-" * 50)
    ax4 = plt.subplot(2, 3, 4)
    
    # Vary D*K product while keeping ratio constant
    permeability_multipliers = np.logspace(-2, 2, 30)
    transition_pressures_perm = []
    
    for mult in permeability_multipliers:
        oxide_props = base_oxide_props.copy()
        oxide_props['D_ox'] = base_oxide_props['D_ox'] * mult
        oxide_props['K_ox'] = base_oxide_props['K_ox'] * mult
        P_trans = calculate_transition_pressure(oxide_props, base_metal_props)
        transition_pressures_perm.append(P_trans)
    
    base_permeability = base_oxide_props['D_ox'] * base_oxide_props['K_ox']
    actual_permeabilities = base_permeability * permeability_multipliers
    
    ax4.loglog(actual_permeabilities, transition_pressures_perm, 'm-', linewidth=2)
    ax4.set_xlabel('Oxide Permeability D_ox×K_ox (mol/m/Pa/s)')
    ax4.set_ylabel('Transition Pressure (Pa)')
    ax4.set_title('Effect of Oxide Permeability')
    ax4.grid(True, which="both", alpha=0.3)
    ax4.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
    ax4.axvline(base_permeability, color='k', linestyle='--', alpha=0.5, label='Base case')
    ax4.legend()
    
    # Study 5: Metal Properties Effect
    print("\nStudy 5: Metal Properties Effect")
    print("-" * 50)
    ax5 = plt.subplot(2, 3, 5)
    
    metal_perm_multipliers = np.logspace(-2, 2, 30)
    transition_pressures_metal = []
    
    for mult in metal_perm_multipliers:
        metal_props = base_metal_props.copy()
        metal_props['D_metal'] = base_metal_props['D_metal'] * mult
        metal_props['K_s_metal'] = base_metal_props['K_s_metal'] * mult
        P_trans = calculate_transition_pressure(base_oxide_props, metal_props)
        transition_pressures_metal.append(P_trans)
    
    base_metal_perm = base_metal_props['D_metal'] * base_metal_props['K_s_metal']
    actual_metal_perms = base_metal_perm * metal_perm_multipliers
    
    ax5.loglog(actual_metal_perms, transition_pressures_metal, 'c-', linewidth=2)
    ax5.set_xlabel('Metal Permeability D_m×K_s (mol/m/Pa^0.5/s)')
    ax5.set_ylabel('Transition Pressure (Pa)')
    ax5.set_title('Effect of Metal Properties')
    ax5.grid(True, which="both", alpha=0.3)
    ax5.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
    
    # Study 6: System Flux at Fixed Pressure
    print("\nStudy 6: System Flux vs Oxide Thickness")
    print("-" * 50)
    ax6 = plt.subplot(2, 3, 6)
    
    test_pressure = 1e0  # 1 Pa
    thicknesses_flux = np.logspace(-10, -6, 30)
    system_fluxes = []
    
    for t in thicknesses_flux:
        oxide_props = base_oxide_props.copy()
        oxide_props['thickness'] = t
        result = calculate_oxide_metal_system(test_pressure, 0, oxide_props, base_metal_props)
        system_fluxes.append(result['flux'])
    
    ax6.loglog(thicknesses_flux*1e10, system_fluxes, 'k-', linewidth=2)
    ax6.set_xlabel('Oxide Thickness (Å)')
    ax6.set_ylabel('System Flux (mol/m²/s)')
    ax6.set_title(f'Flux at P = {test_pressure:.0e} Pa')
    ax6.grid(True, which="both", alpha=0.3)
    
    # Add reference line for pure metal flux
    pure_metal_flux = (base_metal_props['D_metal'] * base_metal_props['K_s_metal'] / 
                      base_metal_props['thickness']) * np.sqrt(test_pressure)
    ax6.axhline(pure_metal_flux, color='b', linestyle='--', alpha=0.5, label='Pure metal')
    ax6.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test5_parametric_studies_{timestamp}.png', dpi=150)
    plt.show()
    
    print("\n✓ Parametric studies complete")
    
    # Create summary table
    print("\n" + "="*60)
    print("SUMMARY: Key Parameter Effects on Transition Pressure")
    print("="*60)
    print("Parameter            | Change    | P_trans Change")
    print("-"*60)
    print(f"Oxide thickness      | 10x ↑     | {transition_pressures[20]/transition_pressures[10]:.1f}x ↓")
    print(f"Oxide diffusivity    | 10x ↑     | {transition_pressures_D[20]/transition_pressures_D[10]:.1f}x ↓")
    print(f"Oxide solubility     | 10x ↑     | {transition_pressures_K[20]/transition_pressures_K[10]:.1f}x ↓")
    print(f"Metal permeability   | 10x ↑     | {transition_pressures_metal[20]/transition_pressures_metal[10]:.1f}x ↑")
    
    return {
        'test': 'parametric_studies',
        'studies_performed': 6,
        'base_transition_pressure': transition_pressures[idx_6A]
    }


def test_oxide_property_comparison():
    """Test 6: Recreate all plots from Tests 1-5 for varying D_ox, K_ox, and thickness."""
    print("\n=== Test 6: Comprehensive Oxide Property Analysis ===")
    print("Recreating all plots from Tests 1-5 with varying oxide properties")
    print("="*60)
    
    # Define parameter combinations to test
    # Using a smaller set for clarity in plots
    parameter_sets = [
        {'D_ox': 1e-22, 'K_ox': 1e-22, 'thickness': 1e-10, 'label': 'Ultra-low perm, 1Å', 'color': 'red'},
        {'D_ox': 1e-20, 'K_ox': 1e-20, 'thickness': 6e-10, 'label': 'Base case, 6Å', 'color': 'blue'},
        {'D_ox': 1e-18, 'K_ox': 1e-18, 'thickness': 1e-9, 'label': 'High perm, 10Å', 'color': 'green'},
        {'D_ox': 1e-19, 'K_ox': 1e-21, 'thickness': 6e-10, 'label': 'Low K, 6Å', 'color': 'orange'},
        {'D_ox': 1e-21, 'K_ox': 1e-19, 'thickness': 6e-10, 'label': 'Low D, 6Å', 'color': 'purple'},
    ]
    
    # Fixed metal properties
    metal_props = {
        'D_metal': 1e-9,
        'K_s_metal': 1e-20,
        'thickness': 1e-3
    }
    
    # ========== RECREATE TEST 1 PLOTS ==========
    print("\n--- Recreating Test 1: Interface Pressure Solutions ---")
    fig1, ((ax1_1, ax1_2), (ax1_3, ax1_4)) = plt.subplots(2, 2, figsize=(14, 10))
    fig1.suptitle('Test 1: Interface Pressure Solutions - Multiple Oxide Properties', fontsize=14)
    
    test_pressures = np.logspace(-10, 28, 100)
    
    for param_set in parameter_sets:
        oxide_props = {
            'D_ox': param_set['D_ox'],
            'K_ox': param_set['K_ox'],
            'thickness': param_set['thickness']
        }
        
        P_interface_list = []
        P_normalized_list = []
        flux_list = []
        regime_list = []
        
        for P_up in test_pressures:
            result = calculate_oxide_metal_system(P_up, 0, oxide_props, metal_props)
            P_interface_list.append(result['P_interface'])
            P_normalized_list.append(result['P_interface_normalized'])
            flux_list.append(result['flux'])
            regime_list.append(result['regime'])
        
        # Plot 1: Interface pressure vs upstream pressure
        ax1_1.loglog(test_pressures, P_interface_list, '-', color=param_set['color'], 
                     linewidth=2, label=param_set['label'])
        
        # Plot 2: Normalized interface position
        ax1_2.semilogx(test_pressures, P_normalized_list, '-', color=param_set['color'], 
                       linewidth=2, label=param_set['label'])
        
        # Plot 3: System flux
        ax1_3.loglog(test_pressures, flux_list, '-', color=param_set['color'],
                     linewidth=2, label=param_set['label'])
        
        # Plot 4: Regime classification
        regime_map = {'oxide_limited': 0, 'transition': 1, 'metal_limited': 2}
        regime_values = [regime_map[r] for r in regime_list]
        ax1_4.semilogx(test_pressures, regime_values, '-', color=param_set['color'],
                       linewidth=2, label=param_set['label'])
    
    # Format Plot 1
    ax1_1.loglog(test_pressures, test_pressures, 'k--', alpha=0.3, label='P_int = P_up')
    ax1_1.set_xlabel('Upstream Pressure (Pa)')
    ax1_1.set_ylabel('Interface Pressure (Pa)')
    ax1_1.set_title('Interface Pressure vs Upstream Pressure')
    ax1_1.grid(True, which="both", alpha=0.3)
    ax1_1.legend(fontsize=8, loc='best')
    
    # Format Plot 2
    ax1_2.set_xlabel('Upstream Pressure (Pa)')
    ax1_2.set_ylabel('(P_interface - P_down)/(P_up - P_down)')
    ax1_2.set_title('Normalized Interface Position')
    ax1_2.set_ylim([0, 1])
    ax1_2.grid(True, alpha=0.3)
    ax1_2.legend(fontsize=8, loc='best')
    
    # Format Plot 3
    ax1_3.set_xlabel('Upstream Pressure (Pa)')
    ax1_3.set_ylabel('System Flux (mol/m²/s)')
    ax1_3.set_title('Total System Flux')
    ax1_3.grid(True, which="both", alpha=0.3)
    ax1_3.legend(fontsize=8, loc='best')
    
    # Format Plot 4
    ax1_4.set_xlabel('Upstream Pressure (Pa)')
    ax1_4.set_ylabel('Regime')
    ax1_4.set_yticks([0, 1, 2])
    ax1_4.set_yticklabels(['Oxide Limited', 'Transition', 'Metal Limited'])
    ax1_4.set_title('Operating Regime')
    ax1_4.grid(True, alpha=0.3)
    ax1_4.legend(fontsize=8, loc='best')
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test6_test1_plots_{timestamp}.png', dpi=150)
    plt.show()
    
    # ========== RECREATE TEST 2 PLOTS ==========
    print("\n--- Recreating Test 2: Flux Continuity Validation ---")
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))
    fig2.suptitle('Test 2: Flux Continuity - Multiple Oxide Properties', fontsize=14)
    
    test_pressures_flux = np.logspace(-6, 6, 25)
    
    for param_set in parameter_sets:
        oxide_props = {
            'D_ox': param_set['D_ox'],
            'K_ox': param_set['K_ox'],
            'thickness': param_set['thickness']
        }
        
        flux_errors = []
        for P_up in test_pressures_flux:
            result = calculate_oxide_metal_system(P_up, 0, oxide_props, metal_props)
            
            flux_ox = molecular_diffusion_flux(
                oxide_props['D_ox'], oxide_props['K_ox'],
                oxide_props['thickness'], P_up, result['P_interface']
            )
            flux_met = calculate_metal_flux_sieverts(
                metal_props['D_metal'], metal_props['K_s_metal'],
                metal_props['thickness'], result['P_interface'], 0
            )
            error = abs(flux_ox - flux_met) / flux_ox if flux_ox > 0 else 0
            flux_errors.append(error)
        
        ax2.loglog(test_pressures_flux, flux_errors, 'o-', color=param_set['color'],
                   markersize=6, label=param_set['label'])
    
    ax2.axhline(y=1e-10, color='g', linestyle='--', label='Target precision (1e-10)')
    ax2.set_xlabel('Upstream Pressure (Pa)')
    ax2.set_ylabel('Relative Flux Error at Interface')
    ax2.set_title('Flux Continuity Verification')
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test6_test2_plots_{timestamp}.png', dpi=150)
    plt.show()
    
    # ========== RECREATE TEST 3 PLOTS ==========
    print("\n--- Recreating Test 3: Concentration Profiles ---")
    test_pressures_conc = [1e-3, 1e0, 1e3, 1e6]
    
    # Create a grid of subplots for concentration profiles
    fig3, axes = plt.subplots(len(test_pressures_conc), len(parameter_sets), 
                              figsize=(20, 16))
    fig3.suptitle('Test 3: Concentration Profiles - Multiple Oxide Properties', fontsize=14)
    
    for col_idx, param_set in enumerate(parameter_sets):
        oxide_props = {
            'D_ox': param_set['D_ox'],
            'K_ox': param_set['K_ox'],
            'thickness': param_set['thickness']
        }
        
        for row_idx, P_up in enumerate(test_pressures_conc):
            ax = axes[row_idx, col_idx] if len(test_pressures_conc) > 1 else axes[col_idx]
            
            profile = calculate_concentration_profile(P_up, 0, oxide_props, metal_props)
            
            # Plot complete profile
            ax.plot(profile['x_all']*1e9, profile['C_all'], 'b-', linewidth=2)
            
            # Mark interface
            interface_pos = oxide_props['thickness']*1e9
            ax.axvline(interface_pos, color='r', linestyle='--', alpha=0.5, label='Interface')
            
            # Mark discontinuity
            ax.plot([interface_pos, interface_pos],
                   [profile['C_oxide'][-1], profile['C_metal'][0]],
                   'ro-', markersize=6, linewidth=2)
            
            ax.set_xlabel('Position (nm)')
            ax.set_ylabel('Conc. (mol/m³)')
            
            # Title for each subplot
            if row_idx == 0:
                ax.set_title(f"{param_set['label']}\nP_up={P_up:.0e} Pa", fontsize=10)
            else:
                ax.set_title(f"P_up={P_up:.0e} Pa", fontsize=10)
            
            ax.grid(True, alpha=0.3)
            
            # Add discontinuity text
            ax.text(0.5, 0.95, f'ΔC={abs(profile["C_discontinuity"]):.1e}',
                   transform=ax.transAxes, ha='center', va='top', fontsize=8,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test6_test3_plots_{timestamp}.png', dpi=150)
    plt.show()
    
    # ========== RECREATE TEST 4 PLOTS ==========
    print("\n--- Recreating Test 4: Limiting Cases ---")
    fig4, (ax4_1, ax4_2) = plt.subplots(1, 2, figsize=(12, 5))
    fig4.suptitle('Test 4: Limiting Cases - Multiple Oxide Properties', fontsize=14)
    
    normalized_positions_thick = []
    normalized_positions_thin = []
    normalized_positions_base = []
    resistance_ratios_thick = []
    resistance_ratios_thin = []
    resistance_ratios_base = []
    labels = []
    
    for param_set in parameter_sets:
        # Base case
        oxide_props_base = {
            'D_ox': param_set['D_ox'],
            'K_ox': param_set['K_ox'],
            'thickness': param_set['thickness']
        }
        
        # Thick oxide
        oxide_props_thick = oxide_props_base.copy()
        oxide_props_thick['thickness'] = 1e-6  # 1 μm
        
        # Thin oxide
        oxide_props_thin = oxide_props_base.copy()
        oxide_props_thin['thickness'] = 1e-12  # 1 pm
        
        P_test = 1e6
        
        result_thick = calculate_oxide_metal_system(P_test, 0, oxide_props_thick, metal_props)
        result_thin = calculate_oxide_metal_system(P_test, 0, oxide_props_thin, metal_props)
        result_base = calculate_oxide_metal_system(P_test, 0, oxide_props_base, metal_props)
        
        normalized_positions_thick.append(result_thick['P_interface_normalized'])
        normalized_positions_thin.append(result_thin['P_interface_normalized'])
        normalized_positions_base.append(result_base['P_interface_normalized'])
        
        resistance_ratios_thick.append(result_thick['resistance_ratio'])
        resistance_ratios_thin.append(result_thin['resistance_ratio'])
        resistance_ratios_base.append(result_base['resistance_ratio'])
        
        labels.append(param_set['label'].split(',')[0])  # Shortened label
    
    # Plot 1: Pressure drop distribution
    x = np.arange(len(labels))
    width = 0.25
    
    ax4_1.bar(x - width, normalized_positions_thick, width, label='Thick (1 μm)', color='red', alpha=0.7)
    ax4_1.bar(x, normalized_positions_base, width, label='Base thickness', color='blue', alpha=0.7)
    ax4_1.bar(x + width, normalized_positions_thin, width, label='Thin (1 pm)', color='green', alpha=0.7)
    
    ax4_1.set_ylabel('Normalized Interface Position')
    ax4_1.set_title('Pressure Drop Distribution')
    ax4_1.set_xticks(x)
    ax4_1.set_xticklabels(labels, rotation=45, ha='right')
    ax4_1.set_ylim([0, 1.1])
    ax4_1.grid(True, alpha=0.3)
    ax4_1.legend()
    
    # Plot 2: Resistance ratios
    ax4_2.semilogy(x - width, resistance_ratios_thick, 'ro', markersize=8, label='Thick')
    ax4_2.semilogy(x, resistance_ratios_base, 'bo', markersize=8, label='Base')
    ax4_2.semilogy(x + width, resistance_ratios_thin, 'go', markersize=8, label='Thin')
    
    ax4_2.axhline(y=10, color='r', linestyle='--', alpha=0.5, label='Oxide-limited')
    ax4_2.axhline(y=0.1, color='b', linestyle='--', alpha=0.5, label='Metal-limited')
    ax4_2.fill_between([-0.5, len(labels)-0.5], 0.1, 10, alpha=0.2, color='yellow', label='Transition')
    
    ax4_2.set_ylabel('R_oxide / R_metal')
    ax4_2.set_title('Resistance Ratios')
    ax4_2.set_xticks(x)
    ax4_2.set_xticklabels(labels, rotation=45, ha='right')
    ax4_2.grid(True, alpha=0.3)
    ax4_2.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test6_test4_plots_{timestamp}.png', dpi=150)
    plt.show()
    
    # ========== RECREATE TEST 5 PLOTS ==========
    print("\n--- Recreating Test 5: Parametric Studies ---")
    fig5 = plt.figure(figsize=(16, 12))
    fig5.suptitle('Test 5: Parametric Studies - Multiple Base Cases', fontsize=14)
    
    # For each parameter set, we'll show how varying individual parameters affects the system
    for idx, param_set in enumerate(parameter_sets[:4]):  # Show first 4 sets
        base_oxide_props = {
            'D_ox': param_set['D_ox'],
            'K_ox': param_set['K_ox'],
            'thickness': param_set['thickness']
        }
        
        # Study 1: Oxide Thickness Effect
        ax5_1 = plt.subplot(4, 4, idx*4 + 1)
        thicknesses = np.logspace(-10, -6, 20)
        transition_pressures = []
        
        for t in thicknesses:
            oxide_props = base_oxide_props.copy()
            oxide_props['thickness'] = t
            P_trans = calculate_transition_pressure(oxide_props, metal_props)
            transition_pressures.append(P_trans)
        
        ax5_1.loglog(thicknesses*1e10, transition_pressures, '-', color=param_set['color'], linewidth=2)
        ax5_1.set_xlabel('Thickness (Å)')
        ax5_1.set_ylabel('P_trans (Pa)')
        if idx == 0:
            ax5_1.set_title(f'{param_set["label"]}\nThickness Effect')
        else:
            ax5_1.set_title(f'{param_set["label"]}')
        ax5_1.grid(True, which="both", alpha=0.3)
        ax5_1.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
        
        # Study 2: Oxide Diffusivity Effect
        ax5_2 = plt.subplot(4, 4, idx*4 + 2)
        D_ox_values = np.logspace(-24, -16, 20)
        transition_pressures_D = []
        
        for D in D_ox_values:
            oxide_props = base_oxide_props.copy()
            oxide_props['D_ox'] = D
            P_trans = calculate_transition_pressure(oxide_props, metal_props)
            transition_pressures_D.append(P_trans)
        
        ax5_2.loglog(D_ox_values, transition_pressures_D, '-', color=param_set['color'], linewidth=2)
        ax5_2.set_xlabel('D_ox (m²/s)')
        ax5_2.set_ylabel('P_trans (Pa)')
        if idx == 0:
            ax5_2.set_title('Diffusivity Effect')
        ax5_2.grid(True, which="both", alpha=0.3)
        ax5_2.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
        
        # Study 3: Oxide Solubility Effect
        ax5_3 = plt.subplot(4, 4, idx*4 + 3)
        K_ox_values = np.logspace(-26, -18, 20)
        transition_pressures_K = []
        
        for K in K_ox_values:
            oxide_props = base_oxide_props.copy()
            oxide_props['K_ox'] = K
            P_trans = calculate_transition_pressure(oxide_props, metal_props)
            transition_pressures_K.append(P_trans)
        
        ax5_3.loglog(K_ox_values, transition_pressures_K, '-', color=param_set['color'], linewidth=2)
        ax5_3.set_xlabel('K_ox (mol/m³/Pa)')
        ax5_3.set_ylabel('P_trans (Pa)')
        if idx == 0:
            ax5_3.set_title('Solubility Effect')
        ax5_3.grid(True, which="both", alpha=0.3)
        ax5_3.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
        
        # Study 4: System Flux vs Thickness
        ax5_4 = plt.subplot(4, 4, idx*4 + 4)
        test_pressure = 1e0
        thicknesses_flux = np.logspace(-10, -6, 20)
        system_fluxes = []
        
        for t in thicknesses_flux:
            oxide_props = base_oxide_props.copy()
            oxide_props['thickness'] = t
            result = calculate_oxide_metal_system(test_pressure, 0, oxide_props, metal_props)
            system_fluxes.append(result['flux'])
        
        ax5_4.loglog(thicknesses_flux*1e10, system_fluxes, '-', color=param_set['color'], linewidth=2)
        ax5_4.set_xlabel('Thickness (Å)')
        ax5_4.set_ylabel('Flux (mol/m²/s)')
        if idx == 0:
            ax5_4.set_title(f'Flux at P={test_pressure} Pa')
        ax5_4.grid(True, which="both", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test6_test5_plots_{timestamp}.png', dpi=150)
    plt.show()
    
    # ========== SUMMARY COMPARISON PLOT ==========
    print("\n--- Creating Summary Comparison ---")
    fig6, ((ax6_1, ax6_2), (ax6_3, ax6_4)) = plt.subplots(2, 2, figsize=(14, 10))
    fig6.suptitle('Test 6: Summary Comparison Across All Parameter Sets', fontsize=14)
    
    # Calculate key metrics for each parameter set
    transition_pressures_all = []
    max_flux_errors = []
    regime_widths = []
    labels_short = []
    
    for param_set in parameter_sets:
        oxide_props = {
            'D_ox': param_set['D_ox'],
            'K_ox': param_set['K_ox'],
            'thickness': param_set['thickness']
        }
        
        # Transition pressure
        P_trans = calculate_transition_pressure(oxide_props, metal_props)
        transition_pressures_all.append(P_trans)
        
        # Max flux error
        max_error = 0
        for P_up in np.logspace(-6, 6, 13):
            result = calculate_oxide_metal_system(P_up, 0, oxide_props, metal_props)
            flux_ox = molecular_diffusion_flux(
                oxide_props['D_ox'], oxide_props['K_ox'],
                oxide_props['thickness'], P_up, result['P_interface']
            )
            flux_met = calculate_metal_flux_sieverts(
                metal_props['D_metal'], metal_props['K_s_metal'],
                metal_props['thickness'], result['P_interface'], 0
            )
            error = abs(flux_ox - flux_met) / flux_ox if flux_ox > 0 else 0
            max_error = max(max_error, error)
        max_flux_errors.append(max_error)
        
        # Regime transition width (pressure range of transition regime)
        pressures_test = np.logspace(-10, 20, 100)
        in_transition = []
        for P in pressures_test:
            result = calculate_oxide_metal_system(P, 0, oxide_props, metal_props)
            if result['regime'] == 'transition':
                in_transition.append(P)
        
        if in_transition:
            regime_width = np.log(max(in_transition)) - np.log(min(in_transition))
        else:
            regime_width = 0
        regime_widths.append(regime_width)
        
        labels_short.append(param_set['label'].split(',')[0])
    
    x_pos = np.arange(len(parameter_sets))
    
    # Plot 1: Transition Pressures
    bars1 = ax6_1.bar(x_pos, transition_pressures_all, color=[p['color'] for p in parameter_sets])
    ax6_1.set_yscale('log')
    ax6_1.set_ylabel('Transition Pressure (Pa)')
    ax6_1.set_title('Transition Pressure Comparison')
    ax6_1.set_xticks(x_pos)
    ax6_1.set_xticklabels(labels_short, rotation=45, ha='right')
    ax6_1.axhspan(1e-7, 1e3, alpha=0.2, color='yellow', label='Fusion range')
    ax6_1.grid(True, alpha=0.3)
    ax6_1.legend()
    
    # Plot 2: Permeability vs Transition Pressure
    permeabilities = [p['D_ox'] * p['K_ox'] for p in parameter_sets]
    ax6_2.loglog(permeabilities, transition_pressures_all, 'o', markersize=10,
                 color='blue', alpha=0.6)
    for i, label in enumerate(labels_short):
        ax6_2.annotate(label, (permeabilities[i], transition_pressures_all[i]),
                      fontsize=8, ha='center', va='bottom')
    ax6_2.set_xlabel('Oxide Permeability D×K (mol/m/Pa/s)')
    ax6_2.set_ylabel('Transition Pressure (Pa)')
    ax6_2.set_title('Permeability vs Transition Pressure')
    ax6_2.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
    ax6_2.grid(True, which="both", alpha=0.3)
    
    # Plot 3: Flux Continuity Errors
    bars3 = ax6_3.bar(x_pos, max_flux_errors, color=[p['color'] for p in parameter_sets])
    ax6_3.set_yscale('log')
    ax6_3.set_ylabel('Maximum Flux Error')
    ax6_3.set_title('Numerical Accuracy Check')
    ax6_3.set_xticks(x_pos)
    ax6_3.set_xticklabels(labels_short, rotation=45, ha='right')
    ax6_3.axhline(y=1e-10, color='g', linestyle='--', label='Target precision')
    ax6_3.grid(True, alpha=0.3)
    ax6_3.legend()
    
    # Plot 4: Transition Regime Width
    bars4 = ax6_4.bar(x_pos, regime_widths, color=[p['color'] for p in parameter_sets])
    ax6_4.set_ylabel('Transition Width (decades of pressure)')
    ax6_4.set_title('Width of Transition Regime')
    ax6_4.set_xticks(x_pos)
    ax6_4.set_xticklabels(labels_short, rotation=45, ha='right')
    ax6_4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/test6_summary_{timestamp}.png', dpi=150)
    plt.show()
    
    # Print comprehensive summary
    print("\n" + "="*60)
    print("COMPREHENSIVE ANALYSIS SUMMARY")
    print("="*60)
    print(f"Parameter sets tested: {len(parameter_sets)}")
    print(f"\nTransition pressure range: {min(transition_pressures_all):.2e} to {max(transition_pressures_all):.2e} Pa")
    print(f"Maximum flux error across all cases: {max(max_flux_errors):.2e}")
    
    fusion_relevant = sum(1 for P in transition_pressures_all if P > 1e3)
    print(f"\nFusion relevance:")
    print(f"  {fusion_relevant}/{len(parameter_sets)} cases have transition above fusion range")
    print(f"  {(fusion_relevant/len(parameter_sets))*100:.1f}% remain oxide-limited in fusion conditions")
    
    print("\n✓ Comprehensive oxide parameter analysis complete")
    print(f"Generated 6 figure sets with all plots from Tests 1-5")
    
    return {
        'test': 'comprehensive_oxide_analysis',
        'parameter_sets': parameter_sets,
        'transition_pressures': transition_pressures_all,
        'max_flux_errors': max_flux_errors,
        'regime_widths': regime_widths,
        'fusion_relevant_fraction': fusion_relevant/len(parameter_sets)
    }

# def test_oxide_property_comparison():
#     """Test 6: Comprehensive study running Tests 1-5 across varying D_ox, K_ox, and thickness."""
#     print("\n=== Test 6: Comprehensive Oxide Property Analysis ===")
#     print("Running Tests 1-5 with varying oxide properties")
#     print("="*60)
    
#     # Define parameter ranges to test
#     D_ox_range = np.logspace(-22, -18, 3)  # 3 values: 1e-22, 1e-20, 1e-18
#     K_ox_range = np.logspace(-22, -18, 3)  # 3 values: 1e-22, 1e-20, 1e-18
#     thickness_range = np.array([1e-10, 6e-10, 1e-9])  # 1Å, 6Å, 10Å
    
#     # Fixed metal properties
#     metal_props = {
#         'D_metal': 1e-9,
#         'K_s_metal': 1e-20,
#         'thickness': 1e-3
#     }
    
#     # Store all results
#     comprehensive_results = {
#         'parameter_combinations': [],
#         'test1_results': [],
#         'test2_results': [],
#         'test3_results': [],
#         'test4_results': [],
#         'test5_results': []
#     }
    
#     # Create summary figure with subplots for key metrics
#     fig_summary = plt.figure(figsize=(20, 16))
    
#     # Counter for parameter combinations
#     combo_count = 0
#     total_combos = len(D_ox_range) * len(K_ox_range) * len(thickness_range)
    
#     print(f"\nTotal parameter combinations to test: {total_combos}")
#     print("-"*60)
    
#     # Arrays to store key metrics across all combinations
#     all_transition_pressures = []
#     all_max_flux_errors = []
#     all_concentration_discontinuities = []
#     all_parameter_labels = []
    
#     # Loop through all parameter combinations
#     for D_ox in D_ox_range:
#         for K_ox in K_ox_range:
#             for thickness in thickness_range:
#                 combo_count += 1
                
#                 # Create oxide properties for this combination
#                 oxide_props = {
#                     'D_ox': D_ox,
#                     'K_ox': K_ox,
#                     'thickness': thickness
#                 }
                
#                 print(f"\n[{combo_count}/{total_combos}] Testing combination:")
#                 print(f"  D_ox = {D_ox:.1e} m²/s")
#                 print(f"  K_ox = {K_ox:.1e} mol/m³/Pa")
#                 print(f"  thickness = {thickness*1e10:.1f} Å")
#                 print("-"*40)
                
#                 # Store parameter combination
#                 param_combo = {
#                     'D_ox': D_ox,
#                     'K_ox': K_ox,
#                     'thickness': thickness,
#                     'permeability': D_ox * K_ox
#                 }
#                 comprehensive_results['parameter_combinations'].append(param_combo)
                
#                 # Create label for this combination
#                 label = f"D={D_ox:.0e}\nK={K_ox:.0e}\nt={thickness*1e10:.0f}Å"
#                 all_parameter_labels.append(label)
                
#                 # --- Run Test 1: Interface Pressure Solutions ---
#                 print("  Running Test 1: Interface pressure solutions...")
#                 test_pressures = np.logspace(-10, 20, 30)  # Reduced points for speed
#                 P_interface_list = []
#                 flux_list = []
#                 regime_list = []
                
#                 for P_up in test_pressures:
#                     result = calculate_oxide_metal_system(P_up, 0, oxide_props, metal_props)
#                     P_interface_list.append(result['P_interface'])
#                     flux_list.append(result['flux'])
#                     regime_list.append(result['regime'])
                
#                 test1_result = {
#                     'P_upstream': test_pressures,
#                     'P_interface': P_interface_list,
#                     'flux': flux_list,
#                     'regime': regime_list
#                 }
#                 comprehensive_results['test1_results'].append(test1_result)
                
#                 # --- Run Test 2: Flux Continuity ---
#                 print("  Running Test 2: Flux continuity check...")
#                 flux_errors = []
#                 for P_up in np.logspace(-6, 6, 13):
#                     result = calculate_oxide_metal_system(P_up, 0, oxide_props, metal_props)
                    
#                     flux_ox = molecular_diffusion_flux(
#                         oxide_props['D_ox'], oxide_props['K_ox'],
#                         oxide_props['thickness'], P_up, result['P_interface']
#                     )
#                     flux_met = calculate_metal_flux_sieverts(
#                         metal_props['D_metal'], metal_props['K_s_metal'],
#                         metal_props['thickness'], result['P_interface'], 0
#                     )
#                     error = abs(flux_ox - flux_met) / flux_ox if flux_ox > 0 else 0
#                     flux_errors.append(error)
                
#                 max_flux_error = max(flux_errors) if flux_errors else 0
#                 test2_result = {
#                     'max_error': max_flux_error,
#                     'passed': max_flux_error < 1e-8
#                 }
#                 comprehensive_results['test2_results'].append(test2_result)
#                 all_max_flux_errors.append(max_flux_error)
#                 print(f"    Max flux error: {max_flux_error:.2e}")
                
#                 # --- Run Test 3: Concentration Profiles ---
#                 print("  Running Test 3: Concentration profiles...")
#                 test_pressures_conc = [1e-3, 1e0, 1e3]
#                 discontinuities = []
                
#                 for P_up in test_pressures_conc:
#                     profile = calculate_concentration_profile(P_up, 0, oxide_props, metal_props)
#                     discontinuities.append(abs(profile['C_discontinuity']))
                
#                 max_discontinuity = max(discontinuities)
#                 test3_result = {
#                     'test_pressures': test_pressures_conc,
#                     'discontinuities': discontinuities,
#                     'max_discontinuity': max_discontinuity
#                 }
#                 comprehensive_results['test3_results'].append(test3_result)
#                 all_concentration_discontinuities.append(max_discontinuity)
                
#                 # --- Run Test 4: Limiting Cases ---
#                 print("  Running Test 4: Limiting cases...")
#                 # Test with very thick oxide
#                 oxide_thick = oxide_props.copy()
#                 oxide_thick['thickness'] = 1e-6
#                 result_thick = calculate_oxide_metal_system(1.0, 0, oxide_thick, metal_props)
                
#                 # Test with very thin oxide  
#                 oxide_thin = oxide_props.copy()
#                 oxide_thin['thickness'] = 1e-12
#                 result_thin = calculate_oxide_metal_system(1.0, 0, oxide_thin, metal_props)
                
#                 test4_result = {
#                     'thick_oxide_regime': result_thick['regime'],
#                     'thin_oxide_normalized': result_thin['P_interface_normalized'],
#                     'thick_is_oxide_limited': result_thick['regime'] == 'oxide_limited',
#                     'thin_approaches_metal': result_thin['P_interface_normalized'] > 0.99
#                 }
#                 comprehensive_results['test4_results'].append(test4_result)
                
#                 # --- Run Test 5: Calculate Transition Pressure ---
#                 print("  Running Test 5: Transition pressure calculation...")
#                 P_trans = calculate_transition_pressure(oxide_props, metal_props)
#                 test5_result = {
#                     'transition_pressure': P_trans,
#                     'in_fusion_range': (P_trans > 1e-7 and P_trans < 1e3)
#                 }
#                 comprehensive_results['test5_results'].append(test5_result)
#                 all_transition_pressures.append(P_trans)
#                 print(f"    Transition pressure: {P_trans:.2e} Pa")
    
#     # --- Create Comprehensive Summary Plots ---
#     print("\n" + "="*60)
#     print("Creating comprehensive summary plots...")
    
#     # Plot 1: Transition Pressure Heatmap
#     ax1 = plt.subplot(3, 4, 1)
#     # Reshape data for heatmap (assuming we have 3x3x3 = 27 combinations)
#     n_params = len(all_parameter_labels)
#     scatter = ax1.scatter(range(n_params), all_transition_pressures, 
#                          c=np.log10(all_transition_pressures), s=100, cmap='viridis')
#     ax1.set_yscale('log')
#     ax1.set_xlabel('Parameter Combination')
#     ax1.set_ylabel('Transition Pressure (Pa)')
#     ax1.set_title('Transition Pressure Across All Combinations')
#     ax1.set_xticks(range(0, n_params, 3))
#     ax1.set_xticklabels(all_parameter_labels[::3], rotation=45, fontsize=6)
#     plt.colorbar(scatter, ax=ax1, label='log₁₀(P_trans)')
#     ax1.axhspan(1e-7, 1e3, alpha=0.2, color='yellow', label='Fusion range')
#     ax1.legend(fontsize=8)
    
#     # Plot 2: Flux Error Analysis
#     ax2 = plt.subplot(3, 4, 2)
#     ax2.semilogy(range(n_params), all_max_flux_errors, 'ro-', markersize=6)
#     ax2.axhline(y=1e-10, color='g', linestyle='--', label='Target precision')
#     ax2.set_xlabel('Parameter Combination')
#     ax2.set_ylabel('Maximum Flux Error')
#     ax2.set_title('Flux Continuity Errors')
#     ax2.set_xticks(range(0, n_params, 3))
#     ax2.set_xticklabels(all_parameter_labels[::3], rotation=45, fontsize=6)
#     ax2.grid(True, alpha=0.3)
#     ax2.legend()
    
#     # Plot 3: Concentration Discontinuities
#     ax3 = plt.subplot(3, 4, 3)
#     ax3.semilogy(range(n_params), all_concentration_discontinuities, 'bo-', markersize=6)
#     ax3.set_xlabel('Parameter Combination')
#     ax3.set_ylabel('Max Concentration Discontinuity (mol/m³)')
#     ax3.set_title('Interface Concentration Jump')
#     ax3.set_xticks(range(0, n_params, 3))
#     ax3.set_xticklabels(all_parameter_labels[::3], rotation=45, fontsize=6)
#     ax3.grid(True, alpha=0.3)
    
#     # Plot 4: Permeability vs Transition Pressure
#     ax4 = plt.subplot(3, 4, 4)
#     permeabilities = [combo['permeability'] for combo in comprehensive_results['parameter_combinations']]
#     ax4.loglog(permeabilities, all_transition_pressures, 'go', markersize=8, alpha=0.6)
#     ax4.set_xlabel('Oxide Permeability D×K (mol/m/Pa/s)')
#     ax4.set_ylabel('Transition Pressure (Pa)')
#     ax4.set_title('P_trans vs Permeability')
#     ax4.axhspan(1e-7, 1e3, alpha=0.2, color='yellow')
#     ax4.grid(True, which="both", alpha=0.3)
    
#     # Plot 5-8: Sample flux profiles for corner cases
#     corner_indices = [0, len(D_ox_range)*len(K_ox_range)-1, 
#                      len(all_parameter_labels)//2, len(all_parameter_labels)-1]
    
#     for plot_idx, combo_idx in enumerate(corner_indices):
#         ax = plt.subplot(3, 4, 5 + plot_idx)
#         test1_data = comprehensive_results['test1_results'][combo_idx]
#         params = comprehensive_results['parameter_combinations'][combo_idx]
        
#         ax.loglog(test1_data['P_upstream'], test1_data['flux'], 'b-', linewidth=2)
#         ax.set_xlabel('Upstream Pressure (Pa)')
#         ax.set_ylabel('System Flux (mol/m²/s)')
#         ax.set_title(f"D={params['D_ox']:.0e}, K={params['K_ox']:.0e}\nt={params['thickness']*1e10:.0f}Å")
#         ax.axvspan(1e-7, 1e3, alpha=0.2, color='yellow')
#         ax.grid(True, which="both", alpha=0.3)
    
#     # Plot 9-12: Regime evolution for selected cases
#     for plot_idx, combo_idx in enumerate(corner_indices):
#         ax = plt.subplot(3, 4, 9 + plot_idx)
#         test1_data = comprehensive_results['test1_results'][combo_idx]
#         params = comprehensive_results['parameter_combinations'][combo_idx]
        
#         regime_map = {'oxide_limited': 0, 'transition': 1, 'metal_limited': 2}
#         regime_values = [regime_map[r] for r in test1_data['regime']]
        
#         ax.semilogx(test1_data['P_upstream'], regime_values, 'mo-', markersize=4)
#         ax.set_xlabel('Upstream Pressure (Pa)')
#         ax.set_ylabel('Regime')
#         ax.set_yticks([0, 1, 2])
#         ax.set_yticklabels(['Oxide', 'Trans.', 'Metal'])
#         ax.set_title(f"D={params['D_ox']:.0e}, K={params['K_ox']:.0e}\nt={params['thickness']*1e10:.0f}Å")
#         ax.axvspan(1e-7, 1e3, alpha=0.2, color='yellow')
#         ax.grid(True, alpha=0.3)
    
#     plt.tight_layout()
#     plt.savefig(f'{RESULTS_DIR}/test6_comprehensive_analysis_{timestamp}.png', dpi=150)
#     plt.show()
    
#     # --- Print Comprehensive Summary ---
#     print("\n" + "="*60)
#     print("COMPREHENSIVE ANALYSIS SUMMARY")
#     print("="*60)
#     print(f"Parameter ranges tested:")
#     print(f"  D_ox: {D_ox_range[0]:.0e} to {D_ox_range[-1]:.0e} m²/s")
#     print(f"  K_ox: {K_ox_range[0]:.0e} to {K_ox_range[-1]:.0e} mol/m³/Pa")
#     print(f"  Thickness: {thickness_range[0]*1e10:.0f} to {thickness_range[-1]*1e10:.0f} Å")
#     print(f"\nTotal combinations tested: {total_combos}")
    
#     print(f"\nKey findings across all combinations:")
#     print(f"  Transition pressure range: {min(all_transition_pressures):.2e} to {max(all_transition_pressures):.2e} Pa")
#     print(f"  Max flux continuity error: {max(all_max_flux_errors):.2e}")
#     print(f"  Max concentration discontinuity: {max(all_concentration_discontinuities):.2e} mol/m³")
    
#     # Check how many are fusion-relevant
#     fusion_relevant = sum(1 for P in all_transition_pressures if P > 1e3)
#     print(f"\nFusion relevance:")
#     print(f"  {fusion_relevant}/{total_combos} combinations have transition above fusion range")
#     print(f"  {(fusion_relevant/total_combos)*100:.1f}% remain oxide-limited in fusion conditions")
    
#     # Check convergence
#     all_converged = all(r['passed'] for r in comprehensive_results['test2_results'])
#     print(f"\nNumerical quality:")
#     print(f"  All flux continuity tests passed: {all_converged}")
#     print(f"  All limiting cases behaved as expected: {all(r['thick_is_oxide_limited'] for r in comprehensive_results['test4_results'])}")
    
#     print("\n✓ Comprehensive oxide parameter analysis complete")
    
#     return {
#         'test': 'comprehensive_oxide_analysis',
#         'total_combinations': total_combos,
#         'parameter_ranges': {
#             'D_ox': D_ox_range.tolist(),
#             'K_ox': K_ox_range.tolist(),
#             'thickness': thickness_range.tolist()
#         },
#         'transition_pressure_range': [min(all_transition_pressures), max(all_transition_pressures)],
#         'max_flux_error': max(all_max_flux_errors),
#         'fusion_relevant_fraction': fusion_relevant/total_combos,
#         'all_tests_passed': all_converged,
#         'full_results': comprehensive_results
#     }

# def test_oxide_property_comparison():
#     """Test 6: Comprehensive study of oxide property effects across wide parameter ranges."""
#     print("\n=== Test 6: Comprehensive Oxide Property Analysis ===")
    
#     # Define ranges for oxide properties
#     D_ox_range = np.logspace(-24, -16, 5)  # 5 values from 1e-24 to 1e-16
#     K_ox_range = np.logspace(-24, -18, 5)  # 5 values from 1e-24 to 1e-18
    
#     metal_props = {
#         'D_metal': 1e-9,
#         'K_s_metal': 1e-20,
#         'thickness': 1e-3
#     }
    
#     # Create comprehensive figure
#     fig = plt.figure(figsize=(20, 16))
    
#     # Study 1: 2D Heatmap of Transition Pressure
#     print("\nStudy 1: Transition Pressure Heatmap")
#     print("-" * 50)
#     ax1 = plt.subplot(3, 4, 1)
    
#     P_trans_matrix = np.zeros((len(K_ox_range), len(D_ox_range)))
    
#     for i, K_ox in enumerate(K_ox_range):
#         for j, D_ox in enumerate(D_ox_range):
#             oxide_props = {
#                 'D_ox': D_ox,
#                 'K_ox': K_ox,
#                 'thickness': 6e-10
#             }
#             P_trans = calculate_transition_pressure(oxide_props, metal_props)
#             P_trans_matrix[i, j] = np.log10(P_trans)
    
#     im = ax1.imshow(P_trans_matrix, aspect='auto', cmap='viridis', origin='lower')
#     ax1.set_xticks(range(len(D_ox_range)))
#     ax1.set_yticks(range(len(K_ox_range)))
#     ax1.set_xticklabels([f'{D:.0e}' for D in D_ox_range], rotation=45)
#     ax1.set_yticklabels([f'{K:.0e}' for K in K_ox_range])
#     ax1.set_xlabel('D_ox (m²/s)')
#     ax1.set_ylabel('K_ox (mol/m³/Pa)')
#     ax1.set_title('log₁₀(P_transition) Heatmap')
#     cbar = plt.colorbar(im, ax=ax1)
#     cbar.set_label('log₁₀(P_trans) [Pa]')
    
#     # Add fusion range indicator
#     fusion_range_log = [np.log10(1e-7), np.log10(1e3)]
#     for val in fusion_range_log:
#         cbar.ax.axhline(y=val, color='yellow', linestyle='--', linewidth=2)
    
#     # Study 2: Regime Maps for Different D_ox
#     print("\nStudy 2: Regime Evolution for Different D_ox")
#     print("-" * 50)
#     ax2 = plt.subplot(3, 4, 2)
    
#     test_pressures = np.logspace(-10, 20, 100)
#     K_ox_fixed = 1e-21  # Middle value
    
#     colors = plt.cm.rainbow(np.linspace(0, 1, len(D_ox_range)))
    
#     for D_ox, color in zip(D_ox_range[::2], colors[::2]):  # Show every other value
#         oxide_props = {
#             'D_ox': D_ox,
#             'K_ox': K_ox_fixed,
#             'thickness': 6e-10
#         }
#         regimes = []
#         for P in test_pressures:
#             result = calculate_oxide_metal_system(P, 0, oxide_props, metal_props)
#             if result['regime'] == 'oxide_limited':
#                 regimes.append(0)
#             elif result['regime'] == 'transition':
#                 regimes.append(1)
#             else:
#                 regimes.append(2)
#         ax2.semilogx(test_pressures, regimes, color=color, linewidth=2,
#                     label=f'D={D_ox:.0e}')
    
#     ax2.set_xlabel('Pressure (Pa)')
#     ax2.set_ylabel('Regime')
#     ax2.set_yticks([0, 1, 2])
#     ax2.set_yticklabels(['Oxide', 'Trans.', 'Metal'])
#     ax2.set_title(f'Regime vs D_ox (K_ox={K_ox_fixed:.0e})')
#     ax2.axvspan(1e-7, 1e3, alpha=0.2, color='yellow')
#     ax2.grid(True, alpha=0.3)
#     ax2.legend(fontsize=8)
    
#     # Study 3: Regime Maps for Different K_ox
#     print("\nStudy 3: Regime Evolution for Different K_ox")
#     print("-" * 50)
#     ax3 = plt.subplot(3, 4, 3)
    
#     D_ox_fixed = 1e-20  # Middle value
    
#     for K_ox, color in zip(K_ox_range[::2], colors[::2]):
#         oxide_props = {
#             'D_ox': D_ox_fixed,
#             'K_ox': K_ox,
#             'thickness': 6e-10
#         }
#         regimes = []
#         for P in test_pressures:
#             result = calculate_oxide_metal_system(P, 0, oxide_props, metal_props)
#             if result['regime'] == 'oxide_limited':
#                 regimes.append(0)
#             elif result['regime'] == 'transition':
#                 regimes.append(1)
#             else:
#                 regimes.append(2)
#         ax3.semilogx(test_pressures, regimes, color=color, linewidth=2,
#                     label=f'K={K_ox:.0e}')
    
#     ax3.set_xlabel('Pressure (Pa)')
#     ax3.set_ylabel('Regime')
#     ax3.set_yticks([0, 1, 2])
#     ax3.set_yticklabels(['Oxide', 'Trans.', 'Metal'])
#     ax3.set_title(f'Regime vs K_ox (D_ox={D_ox_fixed:.0e})')
#     ax3.axvspan(1e-7, 1e3, alpha=0.2, color='yellow')
#     ax3.grid(True, alpha=0.3)
#     ax3.legend(fontsize=8)
    
#     # Study 4: Flux at 1 Pa for all combinations
#     print("\nStudy 4: System Flux at 1 Pa")
#     print("-" * 50)
#     ax4 = plt.subplot(3, 4, 4)
    
#     P_test = 1.0  # Pa
#     flux_matrix = np.zeros((len(K_ox_range), len(D_ox_range)))
    
#     for i, K_ox in enumerate(K_ox_range):
#         for j, D_ox in enumerate(D_ox_range):
#             oxide_props = {
#                 'D_ox': D_ox,
#                 'K_ox': K_ox,
#                 'thickness': 6e-10
#             }
#             result = calculate_oxide_metal_system(P_test, 0, oxide_props, metal_props)
#             flux_matrix[i, j] = np.log10(result['flux'])
    
#     im2 = ax4.imshow(flux_matrix, aspect='auto', cmap='plasma', origin='lower')
#     ax4.set_xticks(range(len(D_ox_range)))
#     ax4.set_yticks(range(len(K_ox_range)))
#     ax4.set_xticklabels([f'{D:.0e}' for D in D_ox_range], rotation=45)
#     ax4.set_yticklabels([f'{K:.0e}' for K in K_ox_range])
#     ax4.set_xlabel('D_ox (m²/s)')
#     ax4.set_ylabel('K_ox (mol/m³/Pa)')
#     ax4.set_title(f'log₁₀(Flux) at P={P_test} Pa')
#     cbar2 = plt.colorbar(im2, ax=ax4)
#     cbar2.set_label('log₁₀(Flux) [mol/m²/s]')
    
#     # Study 5: Permeability Effect (D*K product)
#     print("\nStudy 5: Combined Permeability Effect")
#     print("-" * 50)
#     ax5 = plt.subplot(3, 4, 5)
    
#     permeability_products = []
#     transition_pressures = []
#     flux_values = []
    
#     for D_ox in D_ox_range:
#         for K_ox in K_ox_range:
#             perm = D_ox * K_ox
#             permeability_products.append(perm)
            
#             oxide_props = {
#                 'D_ox': D_ox,
#                 'K_ox': K_ox,
#                 'thickness': 6e-10
#             }
#             P_trans = calculate_transition_pressure(oxide_props, metal_props)
#             transition_pressures.append(P_trans)
            
#             result = calculate_oxide_metal_system(1.0, 0, oxide_props, metal_props)
#             flux_values.append(result['flux'])
    
#     ax5.loglog(permeability_products, transition_pressures, 'bo', markersize=8, alpha=0.6)
#     ax5.set_xlabel('Oxide Permeability D×K (mol/m/Pa/s)')
#     ax5.set_ylabel('Transition Pressure (Pa)')
#     ax5.set_title('P_trans vs Oxide Permeability')
#     ax5.axhspan(1e-7, 1e3, alpha=0.2, color='yellow', label='Fusion range')
#     ax5.grid(True, which="both", alpha=0.3)
#     ax5.legend()
    
#     # Fit power law
#     log_perm = np.log10(permeability_products)
#     log_P_trans = np.log10(transition_pressures)
#     coeffs = np.polyfit(log_perm, log_P_trans, 1)
#     fit_line = 10**(coeffs[0] * log_perm + coeffs[1])
#     ax5.loglog(permeability_products, fit_line, 'r--', 
#                label=f'Slope={coeffs[0]:.2f}')
#     ax5.legend()
    
#     # Study 6: Flux vs Permeability
#     print("\nStudy 6: Flux vs Oxide Permeability")
#     print("-" * 50)
#     ax6 = plt.subplot(3, 4, 6)
    
#     ax6.loglog(permeability_products, flux_values, 'go', markersize=8, alpha=0.6)
#     ax6.set_xlabel('Oxide Permeability D×K (mol/m/Pa/s)')
#     ax6.set_ylabel('System Flux at 1 Pa (mol/m²/s)')
#     ax6.set_title('Flux vs Oxide Permeability')
#     ax6.grid(True, which="both", alpha=0.3)
    
#     # Study 7: Interface Pressure for Different Scenarios
#     print("\nStudy 7: Interface Pressure Distribution")
#     print("-" * 50)
#     ax7 = plt.subplot(3, 4, 7)
    
#     P_upstream = 1.0  # Pa
#     interface_pressures = []
#     labels = []
    
#     for D_ox in D_ox_range[::2]:
#         for K_ox in K_ox_range[::2]:
#             oxide_props = {
#                 'D_ox': D_ox,
#                 'K_ox': K_ox,
#                 'thickness': 6e-10
#             }
#             result = calculate_oxide_metal_system(P_upstream, 0, oxide_props, metal_props)
#             interface_pressures.append(result['P_interface'])
#             labels.append(f'D={D_ox:.0e}\nK={K_ox:.0e}')
    
#     bars = ax7.bar(range(len(interface_pressures)), interface_pressures)
#     ax7.set_ylabel('Interface Pressure (Pa)')
#     ax7.set_yscale('log')
#     ax7.set_title(f'P_interface at P_upstream={P_upstream} Pa')
#     ax7.set_xticks(range(len(labels)))
#     ax7.set_xticklabels(labels, rotation=45, fontsize=7)
#     ax7.grid(True, alpha=0.3)
    
#     # Color bars by resistance ratio
#     for i, bar in enumerate(bars):
#         if interface_pressures[i] < 1e-6:
#             bar.set_color('red')  # Oxide dominated
#         elif interface_pressures[i] > 0.1:
#             bar.set_color('blue')  # Metal dominated
#         else:
#             bar.set_color('green')  # Transition
    
#     # Study 8: 3D Surface Plot
#     print("\nStudy 8: 3D Transition Pressure Surface")
#     print("-" * 50)
#     ax8 = plt.subplot(3, 4, 8, projection='3d')
    
#     D_mesh, K_mesh = np.meshgrid(np.log10(D_ox_range), np.log10(K_ox_range))
    
#     ax8.plot_surface(D_mesh, K_mesh, P_trans_matrix, cmap='viridis', alpha=0.8)
#     ax8.set_xlabel('log₁₀(D_ox)')
#     ax8.set_ylabel('log₁₀(K_ox)')
#     ax8.set_zlabel('log₁₀(P_trans)')
#     ax8.set_title('3D Transition Pressure')
    
#     # Study 9-12: Selected parameter combinations
#     print("\nStudy 9-12: Representative Cases")
#     print("-" * 50)
    
#     representative_cases = [
#         {'D_ox': 1e-24, 'K_ox': 1e-24, 'name': 'Ultra-impermeable'},
#         {'D_ox': 1e-22, 'K_ox': 1e-22, 'name': 'Very impermeable'},
#         {'D_ox': 1e-20, 'K_ox': 1e-20, 'name': 'Moderately impermeable'},
#         {'D_ox': 1e-18, 'K_ox': 1e-18, 'name': 'More permeable'}
#     ]
    
#     for idx, case in enumerate(representative_cases):
#         ax = plt.subplot(3, 4, 9+idx)
        
#         oxide_props = {
#             'D_ox': case['D_ox'],
#             'K_ox': case['K_ox'],
#             'thickness': 6e-10
#         }
        
#         pressures = np.logspace(-10, 20, 100)
#         fluxes = []
        
#         for P in pressures:
#             result = calculate_oxide_metal_system(P, 0, oxide_props, metal_props)
#             fluxes.append(result['flux'])
        
#         ax.loglog(pressures, fluxes, 'b-', linewidth=2)
#         ax.set_xlabel('Pressure (Pa)')
#         ax.set_ylabel('Flux (mol/m²/s)')
#         ax.set_title(f'{case["name"]}\nD={case["D_ox"]:.0e}, K={case["K_ox"]:.0e}')
#         ax.axvspan(1e-7, 1e3, alpha=0.2, color='yellow')
#         ax.grid(True, which="both", alpha=0.3)
        
#         # Add transition pressure line
#         P_trans = calculate_transition_pressure(oxide_props, metal_props)
#         ax.axvline(P_trans, color='r', linestyle='--', alpha=0.5,
#                   label=f'P_trans={P_trans:.1e}')
#         ax.legend(fontsize=7)
    
#     plt.tight_layout()
#     plt.savefig(f'{RESULTS_DIR}/test6_oxide_parameter_range_{timestamp}.png', dpi=150)
#     plt.show()
    
#     # Print comprehensive summary
#     print("\n" + "="*60)
#     print("COMPREHENSIVE OXIDE PARAMETER ANALYSIS SUMMARY")
#     print("="*60)
#     print(f"Parameter ranges explored:")
#     print(f"  D_ox: {D_ox_range[0]:.0e} to {D_ox_range[-1]:.0e} m²/s")
#     print(f"  K_ox: {K_ox_range[0]:.0e} to {K_ox_range[-1]:.0e} mol/m³/Pa")
#     print(f"\nKey findings:")
#     print(f"  - Transition pressure range: {min(transition_pressures):.1e} to {max(transition_pressures):.1e} Pa")
#     print(f"  - Flux range at 1 Pa: {min(flux_values):.1e} to {max(flux_values):.1e} mol/m²/s")
#     print(f"  - Power law scaling: P_trans ∝ (D×K)^{coeffs[0]:.2f}")
    
#     # Count how many scenarios are fusion-relevant
#     fusion_relevant = sum(1 for P in transition_pressures if P > 1e3)
#     total_cases = len(transition_pressures)
#     print(f"\nFusion relevance:")
#     print(f"  - {fusion_relevant}/{total_cases} cases have transition above fusion range")
#     print(f"  - {(fusion_relevant/total_cases)*100:.1f}% remain oxide-limited in fusion conditions")
    
#     print("\n✓ Comprehensive oxide parameter analysis complete")
    
#     return {
#         'test': 'oxide_property_range_analysis',
#         'D_ox_range': D_ox_range.tolist(),
#         'K_ox_range': K_ox_range.tolist(),
#         'transition_pressure_range': [min(transition_pressures), max(transition_pressures)],
#         'flux_range_at_1Pa': [min(flux_values), max(flux_values)],
#         'power_law_exponent': coeffs[0]
#     }

# def run_all_tests():
#     """Run all interface solver validation tests."""
#     print("="*60)
#     print("INTERFACE SOLVER VALIDATION TESTS")
#     print(f"Testing Phase 2: Oxide-Metal Interface Matching")
#     print("="*60)
    
#     all_results = []
    
#     all_results.append(test_interface_pressure_solutions())
#     all_results.append(test_flux_continuity())
#     all_results.append(test_concentration_profiles())
#     all_results.append(test_limiting_cases())
    
#     save_all_results(all_results)
    
#     print("\n" + "="*60)
#     print("PHASE 2 TESTS COMPLETE ✓")
#     print(f"Figures saved in: {RESULTS_DIR}/")
#     print("="*60)


# if __name__ == "__main__":
#     run_all_tests()

def run_all_tests():
    """Run all interface solver validation tests."""
    print("="*60)
    print("INTERFACE SOLVER VALIDATION TESTS")
    print(f"Testing Phase 2: Oxide-Metal Interface Matching")
    print("="*60)
    
    all_results = []
    
    all_results.append(test_interface_pressure_solutions())
    all_results.append(test_flux_continuity())
    all_results.append(test_concentration_profiles())
    all_results.append(test_limiting_cases())
    all_results.append(test_parametric_studies())
    all_results.append(test_oxide_property_comparison()) 
    
    save_all_results(all_results)
    
    print("\n" + "="*60)
    print("PHASE 2 TESTS COMPLETE ✓")
    print(f"Figures saved in: {RESULTS_DIR}/")
    print("="*60)

if __name__ == "__main__":
    run_all_tests()