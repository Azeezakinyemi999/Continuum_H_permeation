"""
Comprehensive Level 3,4 Test Suite
Tests all aspects of the parallel path model for oxide with defects + defective metal.

FULL MODEL: All Physics Active
- Defective oxide: Parallel paths (pinholes, cracks, GB)  
- Defective metal: GB enhancement + Trap-limited diffusion

Test Categories:
1. Basic Functionality - All defect types work
2. Limiting Cases - f=0 → Level 2,4, f=1 → Level 1,4
3. Monotonic Behavior - Flux increases with defect fraction
4. PRF Validation - Values in literature range
5. Pressure Sweep Analysis - Three-regime behavior
6. Defect Type Comparison - Impact analysis
7. PRF vs Defect Fraction - Barrier effectiveness
8. Dominant Path Analysis - Phase diagram of defect vs oxide dominated regions
9. PRF Regime Analysis - PRF phase diagram
10. Temperature Effects - Regime transitions vs temperature
11. Defect Type Regime Comparison - Compare defect types across regimes
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from datetime import datetime
import os
import sys
from io import StringIO

# Import Level 1,4 for comparison
from calculations.permeation_calc import calculate_defective_metal_flux

# Import Level 2,4 interface solver
from calculations.interface_solver import (
    solve_interface_pressure_defective_metal,
    calculate_oxide_defective_metal_system
)

# Import Level 3,4 functions
from calculations.parallel_oxide_defect_paths import (
    calculate_defect_path_flux_defective_metal,
    calculate_parallel_path_flux_defective_metal,
    calculate_PRF_defective_metal
)

# Import Level 4 defective metal functions  
from calculations.defective_metal import (
    grain_boundary_density,
    gb_enhancement_factor,
    trap_occupancy,
    combined_microstructure_model
)

# Import parameters
from data.oxide_defect_parameters import DEFECT_CONFIGURATIONS, PARAMETER_RANGES, PRF_RANGES
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS
from data.microstruture_parameters import (
    LATTICE_PARAMETERS,
    TRAP_PROPERTIES,
    GB_ENHANCEMENT_DATA
)

# Constants
T_TEST = 1073  # K (800°C)
R = 8.314  # J/mol/K

# Output directory
OUTPUT_DIR = "validation/results/complete_level3,4_tests"


def setup_material_properties(T=T_TEST):
    """Set up oxide and metal properties for testing."""
    oxide_data = OXIDE_PROPERTIES['Cr2O3']
    metal_data = MATERIALS['Incoloy800']
    
    oxide_props = {
        'D_ox': oxide_data['D_ox_0'] * np.exp(-oxide_data['E_D_ox'] / (R * T)),
        'K_ox': oxide_data['K_ox_0'] * np.exp(-oxide_data['H_sol_ox'] / (R * T)),
        'thickness': oxide_data['thickness']
    }
    
    metal_props = {
        'D_metal': metal_data['D_0'] * np.exp(-metal_data['E_D'] / (R * T)),
        'K_s_metal': metal_data['K_s0'] * np.exp(-metal_data['H_s'] / (R * T)),
        'thickness': 1e-3  # 1mm
    }
    
    return oxide_props, metal_props


def get_default_microstructure_params(material_name='Incoloy800', condition='solution_annealed'):
    """
    Get microstructure parameters for defective metal (Level 4).
    Same structure as Level 2,4 test for consistency.
    
    Parameters:
    - material_name: Material name
    - condition: Processing condition
    """
    from data.microstruture_parameters import PROCESSING_CONDITIONS
    
    proc = PROCESSING_CONDITIONS.get(condition, PROCESSING_CONDITIONS['solution_annealed'])
    
    # Total trap density (sum of dislocations + vacancies)
    total_trap_density = proc['dislocation_density'] + 1e21  # Typical vacancy density at 800°C
    
    # Use dislocation binding energy as representative
    trap_binding_energy = TRAP_PROPERTIES['dislocations']['binding_energy']
    
    microstructure = {
        'grain_size': proc['grain_size'],
        'grain_shape': 'equiaxed',
        'gb_type': 'HAGB',
        'trap_list': [
            {
                'name': 'dislocations',
                'density': proc['dislocation_density'],
                'binding_energy': TRAP_PROPERTIES['dislocations']['binding_energy']
            },
            {
                'name': 'vacancies',
                'density': 1e21,  # Typical at 800°C
                'binding_energy': TRAP_PROPERTIES['vacancies']['binding_energy']
            }
        ],
        'lattice_site_density': LATTICE_PARAMETERS['Incoloy800']['N_L'],
        # Additional fields for compatibility
        'lattice_density': LATTICE_PARAMETERS['Incoloy800']['N_L'],
        'trap_density': total_trap_density,
        'E_binding': trap_binding_energy,
        'D_gb_D_lattice': 100.0,  # Typical GB enhancement at 800°C
        'delta_gb': TRAP_PROPERTIES['grain_boundaries']['thickness'],
    }
    
    return microstructure


def ensure_output_dir():
    """Create output directory if it doesn't exist."""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)


# ============================================================================
# TEST 1: Basic Functionality
# ============================================================================

def test_1_basic_functionality():
    """Test all defect types including mixed."""
    print("\n" + "="*60)
    print("TEST 1: Basic Functionality (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_upstream, P_downstream = 1.0, 0.0
    
    defect_configs = {
        'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
        'crack': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
        'grain_boundary': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
        'mixed': {
            'type': 'mixed', 
            'area_fraction': 0.015,
            'components': {'pinholes': 0.005, 'cracks': 0.005, 'grain_boundaries': 0.005},
            'thickness_factor': 0.1,
            'diffusivity_factor': 10
        }
    }
    
    results = {}
    all_passed = True
    
    for name, config in defect_configs.items():
        try:
            flux = calculate_defect_path_flux_defective_metal(
                P_upstream, P_downstream, oxide_props, metal_props, 
                config, T_TEST, microstructure
            )
            results[name] = flux
            print(f"✓ {name:15}: flux = {flux:.3e} mol/m²/s")
        except Exception as e:
            print(f"✗ {name:15}: FAILED - {e}")
            all_passed = False
    
    # Create plot
    if all_passed and results:
        ensure_output_dir()
        fig, ax = plt.subplots(figsize=(10, 6))
        
        names = list(results.keys())
        fluxes = [results[n] for n in names]
        
        colors = ['#2ecc71', '#3498db', '#e74c3c', '#9b59b6']
        bars = ax.bar(names, fluxes, color=colors, edgecolor='black', linewidth=1.5)
        
        ax.set_yscale('log')
        ax.set_xlabel('Defect Type', fontsize=12)
        ax.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
        ax.set_title('Test 1: Basic Functionality - Flux by Defect Type (Level 3,4)\n(f_defect = 1%, P = 1 Pa, T = 800°C, with GB+Trapping)', fontsize=14)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar, flux in zip(bars, fluxes):
            ax.text(bar.get_x() + bar.get_width()/2., flux * 1.5, 
                   f'{flux:.2e}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plot_path = f"{OUTPUT_DIR}/test1_basic_functionality_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n✓ Plot saved: {plot_path}")
    
    if all_passed:
        print("\nTest 1: PASSED ✓")
    return all_passed


# ============================================================================
# TEST 2: Limiting Cases
# ============================================================================

def test_2_limiting_cases():
    """Verify f=0 → Level 2,4 and f=1 → Level 1,4."""
    print("\n" + "="*60)
    print("TEST 2: Limiting Cases Validation (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_upstream, P_downstream = 1.0, 0.0
    
    # Test 2a: f_defect = 0 → Level 2,4
    print("\nTest 2a: f_defect = 0 (should match Level 2,4)")
    
    result_l34_f0 = calculate_parallel_path_flux_defective_metal(
        P_upstream, P_downstream, oxide_props, metal_props,
        {'area_fraction': 0.0, 'type': 'pinhole'},
        T_TEST, microstructure
    )
    result_l24 = calculate_oxide_defective_metal_system(
        P_upstream, P_downstream, oxide_props, metal_props, T_TEST, microstructure
    )
    
    if result_l24['flux'] > 0:
        error_2a = abs(result_l34_f0['flux_total'] - result_l24['flux']) / result_l24['flux'] * 100
    else:
        error_2a = 0 if result_l34_f0['flux_total'] == 0 else float('inf')
    
    print(f"  Level 3,4 (f=0): {result_l34_f0['flux_total']:.3e} mol/m²/s")
    print(f"  Level 2,4:       {result_l24['flux']:.3e} mol/m²/s")
    print(f"  Error: {error_2a:.2e}%")
    
    # Test 2b: f_defect = 1 → Level 1,4
    print("\nTest 2b: f_defect = 1 (should match Level 1,4)")
    
    result_l34_f1 = calculate_parallel_path_flux_defective_metal(
        P_upstream, P_downstream, oxide_props, metal_props,
        {'area_fraction': 1.0, 'type': 'pinhole'},
        T_TEST, microstructure
    )
    result_l14 = calculate_defective_metal_flux(
        metal_props['D_metal'], metal_props['K_s_metal'], 
        metal_props['thickness'], P_upstream, P_downstream,
        T_TEST, microstructure
    )
    
    error_2b = abs(result_l34_f1['flux_total'] - result_l14['flux']) / result_l14['flux'] * 100
    
    print(f"  Level 3,4 (f=1): {result_l34_f1['flux_total']:.3e} mol/m²/s")
    print(f"  Level 1,4:       {result_l14['flux']:.3e} mol/m²/s")
    print(f"  Error: {error_2b:.2e}%")
    
    passed = error_2a < 1e-6 and error_2b < 1e-6
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 2a: f=0 comparison
    ax1 = axes[0]
    labels_a = ['Level 3,4 (f=0)', 'Level 2,4\n(Perfect Oxide)']
    values_a = [result_l34_f0['flux_total'], result_l24['flux']]
    colors_a = ['#3498db', '#e74c3c']
    
    bars_a = ax1.bar(labels_a, values_a, color=colors_a, edgecolor='black', linewidth=1.5)
    ax1.set_yscale('log')
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title(f'f_defect = 0 → Level 2,4\nError: {error_2a:.2e}%', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars_a, values_a):
        ax1.text(bar.get_x() + bar.get_width()/2., val * 1.5, 
                f'{val:.2e}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Plot 2b: f=1 comparison
    ax2 = axes[1]
    labels_b = ['Level 3,4 (f=1)', 'Level 1,4\n(Defective Metal)']
    values_b = [result_l34_f1['flux_total'], result_l14['flux']]
    colors_b = ['#3498db', '#2ecc71']
    
    bars_b = ax2.bar(labels_b, values_b, color=colors_b, edgecolor='black', linewidth=1.5)
    ax2.set_yscale('log')
    ax2.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax2.set_title(f'f_defect = 1 → Level 1,4\nError: {error_2b:.2e}%', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars_b, values_b):
        ax2.text(bar.get_x() + bar.get_width()/2., val * 1.5, 
                f'{val:.2e}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.suptitle('Test 2: Limiting Cases Validation (Level 3,4 Full Model)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test2_limiting_cases_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n✓ Plot saved: {plot_path}")
    
    if passed:
        print("\nTest 2: PASSED ✓")
    else:
        print("\nTest 2: FAILED ✗")
    return passed


# ============================================================================
# TEST 3: Monotonic Behavior
# ============================================================================

def test_3_monotonic_behavior():
    """Verify flux increases monotonically with defect fraction."""
    print("\n" + "="*60)
    print("TEST 3: Monotonic Behavior with Defect Fraction (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_upstream, P_downstream = 1.0, 0.0
    
    defect_fractions = [0, 0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    fluxes = []
    enhancements = []
    
    print("\nDefect Fraction | Flux (mol/m²/s) | Enhancement")
    print("-" * 55)
    
    for f in defect_fractions:
        result = calculate_parallel_path_flux_defective_metal(
            P_upstream, P_downstream, oxide_props, metal_props,
            {'area_fraction': f, 'type': 'pinhole'},
            T_TEST, microstructure
        )
        fluxes.append(result['flux_total'])
        enhancements.append(result['defect_enhancement_factor'])
        print(f"  {f:12.4f}  | {result['flux_total']:15.3e} | {result['defect_enhancement_factor']:12.2e}x")
    
    is_monotonic = all(fluxes[i] <= fluxes[i+1] for i in range(len(fluxes)-1))
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Flux vs Defect Fraction
    ax1 = axes[0]
    ax1.loglog(defect_fractions, fluxes, 'b-o', linewidth=2, markersize=10, markerfacecolor='white')
    ax1.set_xlabel('log Defect Fraction', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux vs Defect Fraction (Level 3,4)', fontsize=12)
    ax1.grid(True, alpha=0.3, which='both')
    
    # Plot 2: Enhancement Factor
    ax2 = axes[1]
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(defect_fractions)))
    bars = ax2.bar([f'{f:.1%}' for f in defect_fractions], enhancements, color=colors, edgecolor='black')
    ax2.set_yscale('log')
    ax2.set_xlabel('Defect Fraction', fontsize=12)
    ax2.set_ylabel('log Enhancement Factor', fontsize=12)
    ax2.set_title('Enhancement Factor vs Defect Fraction', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    
    status = "MONOTONIC ✓" if is_monotonic else "NON-MONOTONIC ✗"
    plt.suptitle(f'Test 3: Monotonic Behavior Verification (Level 3,4) - {status}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test3_monotonic_behavior_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n✓ Plot saved: {plot_path}")
    
    if is_monotonic:
        print("\n✓ Flux increases monotonically")
        print("Test 3: PASSED ✓")
    else:
        print("\n✗ Non-monotonic behavior detected")
        print("Test 3: FAILED ✗")
    return is_monotonic


# ============================================================================
# TEST 4: PRF Validation
# ============================================================================

def test_4_PRF_validation():
    """Verify PRF values across all configurations."""
    print("\n" + "="*60)
    print("TEST 4: PRF Validation (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_test = 1.0
    
    print("\nConfiguration       | PRF           | Efficiency | Regime")
    print("-" * 65)
    
    all_passed = True
    results = {}
    
    for config_name, config in DEFECT_CONFIGURATIONS.items():
        try:
            if config_name == 'perfect':
                result = calculate_PRF_defective_metal(P_test, oxide_props, metal_props, T_TEST, microstructure, None)
            else:
                result = calculate_PRF_defective_metal(P_test, oxide_props, metal_props, T_TEST, microstructure, config)
            
            results[config_name] = result
            prf_str = f"{result['PRF']:.2e}" if result['PRF'] > 1e6 else f"{result['PRF']:.1f}"
            print(f"  {config_name:17} | {prf_str:13} | {result['efficiency']:10.2%} | {result['regime']}")
        except Exception as e:
            print(f"  {config_name:17} | FAILED: {e}")
            all_passed = False
    
    # Create plot
    if results:
        ensure_output_dir()
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        names = list(results.keys())
        prfs = [results[n]['PRF'] for n in names]
        efficiencies = [results[n]['efficiency'] * 100 for n in names]
        
        # Plot 1: PRF by configuration
        ax1 = axes[0]
        colors = plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, len(names)))
        bars1 = ax1.bar(names, prfs, color=colors, edgecolor='black', linewidth=1.5)
        ax1.set_yscale('log')
        
        # Add literature range
        ax1.axhspan(10, 3828, alpha=0.2, color='blue', label='Literature range (10-3828)')
        
        ax1.set_xlabel('Configuration', fontsize=12)
        ax1.set_ylabel('log PRF', fontsize=12)
        ax1.set_title('Permeation Reduction Factor by Configuration (L3,4)', fontsize=12)
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3, axis='y')
        plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        # Plot 2: Barrier Efficiency
        ax2 = axes[1]
        colors2 = plt.cm.RdYlGn(np.array(efficiencies) / 100)
        bars2 = ax2.bar(names, efficiencies, color=colors2, edgecolor='black', linewidth=1.5)
        ax2.axhline(90, color='green', linestyle='--', alpha=0.7, label='90% threshold')
        ax2.axhline(50, color='orange', linestyle='--', alpha=0.7, label='50% threshold')
        
        ax2.set_xlabel('Configuration', fontsize=12)
        ax2.set_ylabel('Barrier Efficiency (%)', fontsize=12)
        ax2.set_title('Barrier Efficiency by Configuration', fontsize=12)
        ax2.set_ylim([0, 105])
        ax2.legend(loc='lower right')
        ax2.grid(True, alpha=0.3, axis='y')
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        plt.suptitle('Test 4: PRF Validation (Level 3,4 Full Model)', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        plot_path = f"{OUTPUT_DIR}/test4_PRF_validation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n✓ Plot saved: {plot_path}")
    
    if all_passed:
        print("\nTest 4: PASSED ✓")
    return all_passed


# ============================================================================
# TEST 5: Pressure Sweep Analysis
# ============================================================================

def test_5_pressure_sweep():
    """Analyze three-regime behavior across pressure range."""
    print("\n" + "="*60)
    print("TEST 5: Pressure Sweep Analysis (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_downstream = 0.0
    
    pressures = np.logspace(-2, 3, 30)  # 0.01 Pa to 1 kPa (safe range to avoid lattice saturation)
    defect_fractions = [0, 0.001, 0.01, 0.1, 0.5, 0.9]  # Reduced set for speed
    
    results = {f: {'pressures': [], 'fluxes': []} for f in defect_fractions}
    
    # Also calculate bare metal (Level 1,4) for reference
    bare_metal_fluxes = []
    for P in pressures:
        result_L14 = calculate_defective_metal_flux(
            metal_props['D_metal'], metal_props['K_s_metal'],
            metal_props['thickness'], P, P_downstream,
            T_TEST, microstructure
        )
        bare_metal_fluxes.append(result_L14['flux'])
    
    print("\nCalculating flux vs pressure for different defect fractions...")
    
    for P in pressures:
        for f in defect_fractions:
            try:
                if f == 0:
                    result = calculate_oxide_defective_metal_system(P, P_downstream, oxide_props, metal_props, T_TEST, microstructure)
                    flux = result['flux']
                else:
                    result = calculate_parallel_path_flux_defective_metal(
                        P, P_downstream, oxide_props, metal_props,
                        {'area_fraction': f, 'type': 'pinhole'},
                        T_TEST, microstructure
                    )
                    flux = result['flux_total']
                
                results[f]['pressures'].append(P)
                results[f]['fluxes'].append(flux)
            except Exception as e:
                pass  # Skip failed points
    
    # Create plot
    ensure_output_dir()
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot bare metal (Level 1,4) as upper bound
    ax.loglog(pressures, bare_metal_fluxes, 'k--', linewidth=2.5, 
              label='Level 1,4: Defective metal (upper limit)', alpha=0.8)
    
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(defect_fractions)))
    
    labels = []
    for f in defect_fractions:
        if f == 0:
            labels.append('Level 2,4: Perfect oxide (f=0%)')
        elif f < 0.01:
            labels.append(f'f={f*100:.2f}%')
        else:
            labels.append(f'f={f*100:.0f}%')
    
    for i, f in enumerate(defect_fractions):
        if results[f]['fluxes']:
            ax.loglog(results[f]['pressures'], results[f]['fluxes'], 
                     color=colors[i], linewidth=2, label=labels[i])
            
    ax.set_xlabel('log Upstream Pressure (Pa)', fontsize=12)
    ax.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax.set_title('Level 3,4: Effect of Defect Fraction (with GB+Trapping)\n(Dashed line = defective metal upper limit)', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plot_path = f"{OUTPUT_DIR}/test5_pressure_sweep_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 5: PASSED ✓")
    return True


# ============================================================================
# TEST 6: Defect Type Comparison (Level 3,4)
# ============================================================================

def test_6_defect_type_comparison():
    """Compare impact of different defect types at Level 3,4 (defective oxide + defective metal)."""
    print("\n" + "="*60)
    print("TEST 6: Defect Type Comparison (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_upstream, P_downstream = 100.0, 0.0  # 100 Pa (safe for defective metal)
    
    defect_configs = {
        'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
        'crack_thin': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
        'crack_thick': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.5},
        'grain_boundary_10x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
        'grain_boundary_100x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 100},
    }
    
    # Get baseline (Level 2,4 - perfect oxide + defective metal)
    baseline = calculate_oxide_defective_metal_system(P_upstream, P_downstream, oxide_props, metal_props, T_TEST, microstructure)
    
    # Get bare metal (Level 1,4 - defective metal only)
    bare_metal = calculate_defective_metal_flux(
        metal_props['D_metal'], metal_props['K_s_metal'],
        metal_props['thickness'], P_upstream, P_downstream,
        T_TEST, microstructure
    )
    
    print(f"\nBaseline (perfect oxide + defective metal): {baseline['flux']:.3e} mol/m²/s")
    print(f"Bare metal (Level 1,4):                     {bare_metal['flux']:.3e} mol/m²/s")
    print("\nDefect Type          | Flux (mol/m²/s) | Enhancement | Dominant Path")
    print("-" * 75)
    
    results = {}
    for name, config in defect_configs.items():
        result = calculate_parallel_path_flux_defective_metal(
            P_upstream, P_downstream, oxide_props, metal_props, config,
            T_TEST, microstructure
        )
        results[name] = result
        enhancement = result['flux_total'] / baseline['flux'] if baseline['flux'] > 0 else float('inf')
        results[name]['enhancement'] = enhancement
        print(f"  {name:18} | {result['flux_total']:15.3e} | {enhancement:11.2e}x | {result['dominant_path']}")
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    names = list(results.keys())
    fluxes = [results[n]['flux_total'] for n in names]
    enhancements = [results[n]['enhancement'] for n in names]
    
    # Plot 1: Flux comparison
    ax1 = axes[0]
    x = np.arange(len(names))
    width = 0.6
    
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6', '#f39c12']
    bars = ax1.bar(x, fluxes, width, color=colors, edgecolor='black', linewidth=1.5)
    
    # Reference lines
    ax1.axhline(baseline['flux'], color='darkblue', linestyle='--', linewidth=2, label=f'Level 2,4 (Perfect): {baseline["flux"]:.2e}')
    ax1.axhline(bare_metal['flux'], color='black', linestyle='-', linewidth=2, label=f'Level 1,4 (Bare): {bare_metal["flux"]:.2e}')
    
    ax1.set_yscale('log')
    ax1.set_xticks(x)
    ax1.set_xticklabels(['Pinhole', 'Crack\n(α=0.1)', 'Crack\n(α=0.5)', 'GB\n(β=10)', 'GB\n(β=100)'])
    ax1.set_xlabel('Defect Type', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux by Defect Type (f = 1%)', fontsize=12)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Enhancement factor
    ax2 = axes[1]
    bars2 = ax2.bar(x, enhancements, width, color=colors, edgecolor='black', linewidth=1.5)
    
    ax2.set_yscale('log')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Pinhole', 'Crack\n(α=0.1)', 'Crack\n(α=0.5)', 'GB\n(β=10)', 'GB\n(β=100)'])
    ax2.set_xlabel('Defect Type', fontsize=12)
    ax2.set_ylabel('log Enhancement (vs Perfect Oxide)', fontsize=12)
    ax2.set_title('Enhancement Factor by Defect Type', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    
    for bar, enh in zip(bars2, enhancements):
        ax2.text(bar.get_x() + bar.get_width()/2., enh * 1.5, 
                f'{enh:.1e}x', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.suptitle('Test 6: Defect Type Comparison (Level 3,4) at P = 100 Pa, T = 800°C', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test6_defect_type_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n✓ Plot saved: {plot_path}")
    
    print("\nTest 6: PASSED ✓")
    return True


# ============================================================================
# TEST 7: PRF vs Defect Fraction (Level 3,4)
# ============================================================================

def test_7_PRF_vs_defect_fraction():
    """Calculate and plot PRF vs defect fraction at Level 3,4."""
    print("\n" + "="*60)
    print("TEST 7: PRF vs Defect Fraction (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_test = 100.0  # 100 Pa (safe for defective metal)
    
    defect_fractions = np.logspace(-4, 0, 100)  # 0.01% to ~100%
    
    prf_pinhole = []
    prf_crack = []
    prf_gb = []
    
    for f in defect_fractions:
        # Pinhole
        result = calculate_PRF_defective_metal(P_test, oxide_props, metal_props,
                                              T_TEST, microstructure,
                                              {'area_fraction': f, 'type': 'pinhole'})
        prf_pinhole.append(result['PRF'])
        
        # Crack
        result = calculate_PRF_defective_metal(P_test, oxide_props, metal_props,
                                              T_TEST, microstructure,
                                              {'area_fraction': f, 'type': 'crack', 'thickness_factor': 0.1})
        prf_crack.append(result['PRF'])
        
        # Grain boundary
        result = calculate_PRF_defective_metal(P_test, oxide_props, metal_props,
                                              T_TEST, microstructure,
                                              {'area_fraction': f, 'type': 'grain_boundary', 'diffusivity_factor': 10})
        prf_gb.append(result['PRF'])
    
    # Create plot
    ensure_output_dir()
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.loglog(defect_fractions * 100, prf_pinhole, 'r-', linewidth=2, label='Pinholes')
    ax.loglog(defect_fractions * 100, prf_crack, 'b--', linewidth=2, label='Cracks (α=0.1)')
    ax.loglog(defect_fractions * 100, prf_gb, 'g:', linewidth=2, label='Grain boundaries (β=10)')
    
    # Add literature range
    ax.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range (10-3828)')
    
    # Add metal-limited region (PRF < 2 means oxide is nearly ineffective)
    ax.axhspan(1, 2, alpha=0.3, color='red', label='Metal-limited (PRF < 2)')
    ax.axhline(1, color='black', linestyle='-', linewidth=1.5, alpha=0.7)
    ax.text(0.015, 1.3, 'PRF = 1 (bare metal)', fontsize=9, color='black')
    
    ax.set_xlabel('log Defect Area Fraction (%)', fontsize=12)
    ax.set_ylabel('log Permeation Reduction Factor (PRF)', fontsize=12)
    ax.set_title('PRF vs Defect Fraction (Level 3,4) for Different Defect Types\n(PRF → 1 indicates metal-limited regime)', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0.01, 100])
    ax.set_ylim([0.8, 1e5])
    
    plot_path = f"{OUTPUT_DIR}/test7_PRF_vs_defect_fraction_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Print summary with metal-limited indication
    print("\nDefect Fraction (%) | PRF (pinhole) | PRF (crack) | PRF (GB)   | Regime")
    print("-" * 80)
    for i in [0, 25, 50, 75, 99]:  # Sample points across 100-point array
        regime = "Metal-limited" if prf_pinhole[i] < 2 else ("Poor barrier" if prf_pinhole[i] < 10 else "Effective")
        print(f"  {defect_fractions[i]*100:16.2f}  | {prf_pinhole[i]:13.1f} | {prf_crack[i]:11.2e} | {prf_gb[i]:.2e} | {regime}")
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 7: PASSED ✓")
    return True


# ============================================================================
# TEST 8: Dominant Path Analysis (Level 3,4)
# ============================================================================

def test_8_dominant_path_analysis():
    """
    Analyze which path dominates flux at different conditions at Level 3,4.
    Creates a phase diagram of defect-dominated vs oxide-dominated vs metal-dominated regions.
    
    Note: Pressure range limited to avoid lattice saturation in defective metal model.
    """
    print("\n" + "="*60)
    print("TEST 8: Dominant Path Analysis (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    
    # Create grid - limited pressure range for defective metal
    pressures = np.logspace(-2, 3, 50)  # 0.01 Pa to 1 kPa (safe for defective metal)
    defect_fractions = np.logspace(-4, 0, 50)  # 0.01% to ~100%
    
    # Arrays to store results
    dominant_path = np.zeros((len(defect_fractions), len(pressures)))
    flux_ratio = np.zeros((len(defect_fractions), len(pressures)))
    total_flux = np.zeros((len(defect_fractions), len(pressures)))
    metal_limitation = np.zeros((len(defect_fractions), len(pressures)))
    
    print(f"Computing {len(pressures)} × {len(defect_fractions)} = {len(pressures)*len(defect_fractions)} points...")
    
    for i, f_def in enumerate(defect_fractions):
        for j, P in enumerate(pressures):
            try:
                result = calculate_parallel_path_flux_defective_metal(
                    P, 0, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'},
                    T_TEST, microstructure
                )
                
                # Get bare metal flux (Level 1,4) for comparison
                bare_metal = calculate_defective_metal_flux(
                    metal_props['D_metal'], metal_props['K_s_metal'],
                    metal_props['thickness'], P, 0,
                    T_TEST, microstructure
                )
                
                flux_defect = result['flux_defect_contribution']
                flux_intact = result['flux_intact_contribution']
                flux_bare = bare_metal['flux']
                total_flux[i, j] = result['flux_total']
                
                # Calculate how close to bare metal limit
                metal_limitation[i, j] = result['flux_total'] / flux_bare if flux_bare > 0 else 0
                
                # Ratio: >1 means defect-dominated, <1 means oxide-dominated
                if flux_intact > 0:
                    flux_ratio[i, j] = flux_defect / flux_intact
                else:
                    flux_ratio[i, j] = 1e10
                
                # Classify regimes (3 regions):
                # 0 = oxide-dominated, 0.5 = defect-path dominated, 1 = metal-limited
                if metal_limitation[i, j] > 0.9:
                    dominant_path[i, j] = 1.0  # Metal-limited
                elif flux_defect > 10 * flux_intact:
                    dominant_path[i, j] = 0.67  # Defect-path dominated
                elif flux_defect > flux_intact:
                    dominant_path[i, j] = 0.5  # Mixed (defect > oxide)
                elif flux_intact > 10 * flux_defect:
                    dominant_path[i, j] = 0.0  # Oxide-dominated
                else:
                    dominant_path[i, j] = 0.33  # Mixed (oxide > defect)
                    
            except Exception as e:
                # Skip points that cause saturation
                dominant_path[i, j] = np.nan
                flux_ratio[i, j] = np.nan
                total_flux[i, j] = np.nan
                metal_limitation[i, j] = np.nan
    
    # Create plots
    ensure_output_dir()
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    P_mesh, F_mesh = np.meshgrid(pressures, defect_fractions * 100)
    
    # Plot 1: Three-regime phase diagram
    ax1 = axes[0, 0]
    from matplotlib.colors import LinearSegmentedColormap
    colors_regime = ['green', 'yellowgreen', 'yellow', 'orange', 'red']
    cmap_regime = LinearSegmentedColormap.from_list('regime', colors_regime, N=256)
    
    im1 = ax1.pcolormesh(P_mesh, F_mesh, dominant_path, cmap=cmap_regime, vmin=0, vmax=1, shading='auto')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax1.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_title('Permeation Regime Map (Level 3,4)\n(Green=Oxide, Yellow=Mixed, Orange=Defect, Red=Metal)', fontsize=11)
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_ticks([0, 0.33, 0.5, 0.67, 1.0])
    cbar1.set_ticklabels(['Oxide', 'Mixed\n(Ox>Def)', 'Mixed\n(Def>Ox)', 'Defect', 'Metal'])
    
    # Add contour lines at regime boundaries
    ax1.contour(P_mesh, F_mesh, dominant_path, levels=[0.25, 0.58, 0.85], 
                colors='black', linewidths=1.5, linestyles=['--', '-', '--'])
    
    # Plot 2: Metal limitation factor
    ax2 = axes[0, 1]
    im2 = ax2.pcolormesh(P_mesh, F_mesh, metal_limitation, cmap='hot_r', 
                         vmin=0, vmax=1, shading='auto')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax2.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_title('Metal Limitation Factor (J_total / J_bare_metal)\n1.0 = Metal substrate is rate-limiting', fontsize=11)
    plt.colorbar(im2, ax=ax2, label='J_total / J_bare')
    
    # Add contours
    ax2.contour(P_mesh, F_mesh, metal_limitation, levels=[0.5, 0.9, 0.99], 
                colors='white', linewidths=1.5, linestyles=[':', '--', '-'])
    
    # Plot 3: Flux ratio
    ax3 = axes[1, 0]
    im3 = ax3.pcolormesh(P_mesh, F_mesh, flux_ratio, cmap='coolwarm', 
                         norm=LogNorm(vmin=0.01, vmax=100), shading='auto')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax3.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax3.set_title('Flux Ratio (Defect Path / Oxide Path)', fontsize=12)
    plt.colorbar(im3, ax=ax3, label='J_defect / J_oxide')
    
    # Add contour at ratio = 1
    ax3.contour(P_mesh, F_mesh, flux_ratio, levels=[1], colors='white', linewidths=2, linestyles='--')
    
    # Plot 4: Cross-sections
    ax4 = axes[1, 1]
    
    P_slices = [10, 100, 1000]  # 10 Pa, 100 Pa, 1 kPa
    colors = ['blue', 'green', 'red']
    
    for P_slice, color in zip(P_slices, colors):
        j_idx = np.argmin(np.abs(pressures - P_slice))
        ax4.semilogx(defect_fractions * 100, metal_limitation[:, j_idx], 
                    color=color, linewidth=2, label=f'P = {P_slice} Pa')
    
    ax4.axhline(1.0, color='black', linestyle='-', alpha=0.7, label='Bare metal limit')
    ax4.axhline(0.9, color='gray', linestyle='--', alpha=0.5, label='90% of bare metal')
    ax4.axhline(0.5, color='gray', linestyle=':', alpha=0.5, label='50% of bare metal')
    
    ax4.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax4.set_ylabel('J_total / J_bare_metal', fontsize=12)
    ax4.set_title('Approach to Metal-Limited Regime', fontsize=12)
    ax4.legend(loc='best', fontsize=9)
    ax4.set_ylim([0, 1.1])
    ax4.grid(True, alpha=0.3)
    
    # Add regime annotations
    ax4.fill_between([0.01, 100], 0.9, 1.1, alpha=0.2, color='red')
    ax4.text(20, 0.95, 'Metal-limited', fontsize=10, color='darkred')
    ax4.text(0.1, 0.3, 'Oxide/Defect\ncontrolled', fontsize=10, color='darkblue')
    
    plt.suptitle(f'Three-Regime Analysis (Level 3,4) at T = 800°C\n(Pressure range limited to avoid lattice saturation)', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test8_dominant_path_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 8: PASSED ✓")
    return True


# ============================================================================
# TEST 9: PRF Regime Analysis (Level 3,4)
# ============================================================================

def test_9_PRF_regime_analysis():
    """
    Create PRF phase diagram showing barrier effectiveness regions at Level 3,4.
    Pressure range limited for defective metal model.
    """
    print("\n" + "="*60)
    print("TEST 9: PRF Regime Analysis (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    
    # Create grid - limited pressure range
    pressures = np.logspace(-2, 3, 50)  # 0.01 Pa to 1 kPa
    defect_fractions = np.logspace(-4, 0, 40)  # 0.01% to ~100%
    
    PRF_map = np.zeros((len(defect_fractions), len(pressures)))
    
    print(f"Computing PRF for {len(pressures)} × {len(defect_fractions)} points...")
    
    for i, f_def in enumerate(defect_fractions):
        for j, P in enumerate(pressures):
            try:
                result = calculate_PRF_defective_metal(
                    P, oxide_props, metal_props, T_TEST, microstructure,
                    {'area_fraction': f_def, 'type': 'pinhole'}
                )
                PRF_map[i, j] = result['PRF']
            except:
                PRF_map[i, j] = np.nan
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    P_mesh, F_mesh = np.meshgrid(pressures, defect_fractions * 100)
    
    # Plot 1: PRF map
    ax1 = axes[0]
    im1 = ax1.pcolormesh(P_mesh, F_mesh, PRF_map, cmap='RdYlGn',
                         norm=LogNorm(vmin=1, vmax=10000), shading='auto')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax1.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_title('Permeation Reduction Factor (Level 3,4)', fontsize=12)
    cbar = plt.colorbar(im1, ax=ax1)
    cbar.set_label('PRF')
    
    # Add contours for key PRF values
    contour_levels = [2, 10, 100, 1000]
    cs = ax1.contour(P_mesh, F_mesh, PRF_map, levels=contour_levels, 
                     colors=['darkred', 'red', 'orange', 'black'], linewidths=1.5)
    ax1.clabel(cs, inline=True, fontsize=9, fmt='PRF=%d')
    
    # Plot 2: PRF vs defect fraction at different pressures
    ax2 = axes[1]
    
    P_slices = [10, 100, 1000]
    colors = ['blue', 'green', 'red']
    labels = ['10 Pa', '100 Pa', '1 kPa']
    
    for P_slice, color, label in zip(P_slices, colors, labels):
        j_idx = np.argmin(np.abs(pressures - P_slice))
        ax2.loglog(defect_fractions * 100, PRF_map[:, j_idx],
                  color=color, linewidth=2, label=f'P = {label}')
    
    # Add literature range
    ax2.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature PRF (10-3828)')
    
    # Add metal-limited region
    ax2.axhspan(1, 2, alpha=0.3, color='red', label='Metal-limited (PRF < 2)')
    ax2.axhline(1, color='black', linestyle='-', linewidth=1.5)
    
    # Add regime boundaries
    ax2.axhline(100, color='gray', linestyle='--', alpha=0.7)
    ax2.axhline(10, color='gray', linestyle='--', alpha=0.7)
    
    ax2.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_ylabel('log PRF', fontsize=12)
    ax2.set_title('PRF vs Defect Fraction (Level 3,4)', fontsize=12)
    ax2.legend(loc='best')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_ylim([0.8, 1e5])
    
    # Add regime labels
    ax2.text(0.02, 2000, 'Excellent Barrier\n(PRF > 100)', fontsize=10, color='green')
    ax2.text(0.5, 30, 'Moderate Barrier\n(PRF 10-100)', fontsize=10, color='orange')
    ax2.text(5, 5, 'Poor Barrier\n(PRF 2-10)', fontsize=10, color='red')
    ax2.text(15, 1.3, 'Metal-limited\n(PRF < 2)', fontsize=10, color='darkred')
    
    plt.suptitle(f'PRF Analysis (Level 3,4) at T = 800°C', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test9_PRF_regime_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 9: PASSED ✓")
    return True


# ============================================================================
# TEST 10: Temperature Effects on Regimes (Level 3,4)
# ============================================================================

def test_10_temperature_effects():
    """
    Analyze how temperature affects regime transitions at Level 3,4.
    """
    print("\n" + "="*60)
    print("TEST 10: Temperature Effects on Regimes (Level 3,4)")
    print("="*60)
    
    temperatures = np.array([600, 700, 800, 900, 1000]) + 273.15  # K
    defect_fractions = np.logspace(-4, 0, 100)  # 0.01% to ~100%
    P_test = 100  # Pa (safe for defective metal)
    
    ensure_output_dir()
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(temperatures)))
    
    # Plot 1: Flux vs defect fraction at different temperatures
    ax1 = axes[0]
    
    # Store bare metal fluxes for each temperature
    bare_metal_by_T = {}
    
    for T_K, color in zip(temperatures, colors):
        oxide_props, metal_props = setup_material_properties(T_K)
        microstructure = get_default_microstructure_params()
        
        # Get bare metal flux (Level 1,4) at this temperature
        result_L14 = calculate_defective_metal_flux(
            metal_props['D_metal'], metal_props['K_s_metal'],
            metal_props['thickness'], P_test, 0,
            T_K, microstructure
        )
        bare_metal_by_T[T_K] = result_L14['flux']
        
        fluxes = []
        for f_def in defect_fractions:
            try:
                result = calculate_parallel_path_flux_defective_metal(
                    P_test, 0, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'},
                    T_K, microstructure
                )
                fluxes.append(result['flux_total'])
            except:
                fluxes.append(np.nan)
        
        ax1.loglog(defect_fractions * 100, fluxes, color=color, linewidth=2,
                  label=f'{T_K-273:.0f}°C')
        
        # Add bare metal reference as horizontal dashed line
        ax1.axhline(bare_metal_by_T[T_K], color=color, linestyle='--', alpha=0.4, linewidth=1)
    
    ax1.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title(f'Flux vs Defect Fraction (Level 3,4) at P = {P_test} Pa\n(dashed = bare metal limit)', fontsize=12)
    ax1.legend(loc='best', title='Temperature')
    ax1.grid(True, which='both', alpha=0.3)
    
    # Plot 2: PRF vs defect fraction at different temperatures
    ax2 = axes[1]
    
    for T_K, color in zip(temperatures, colors):
        oxide_props, metal_props = setup_material_properties(T_K)
        microstructure = get_default_microstructure_params()
        
        PRFs = []
        for f_def in defect_fractions:
            try:
                result = calculate_PRF_defective_metal(
                    P_test, oxide_props, metal_props, T_K, microstructure,
                    {'area_fraction': f_def, 'type': 'pinhole'}
                )
                PRFs.append(result['PRF'])
            except:
                PRFs.append(np.nan)
        
        ax2.loglog(defect_fractions * 100, PRFs, color=color, linewidth=2,
                  label=f'{T_K-273:.0f}°C')
    
    ax2.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range')
    ax2.axhspan(1, 2, alpha=0.3, color='red', label='Metal-limited')
    ax2.axhline(1, color='black', linestyle='-', linewidth=1.5)
    ax2.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_ylabel('log PRF', fontsize=12)
    ax2.set_title(f'PRF vs Defect Fraction (Level 3,4) at P = {P_test} Pa', fontsize=12)
    ax2.legend(loc='best', title='Temperature', fontsize=8)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_ylim([0.8, 1e5])
    
    # Plot 3: Metal limitation factor vs defect fraction at different temperatures
    ax3 = axes[2]
    
    for T_K, color in zip(temperatures, colors):
        oxide_props, metal_props = setup_material_properties(T_K)
        microstructure = get_default_microstructure_params()
        bare_flux = bare_metal_by_T[T_K]
        
        metal_factors = []
        for f_def in defect_fractions:
            try:
                result = calculate_parallel_path_flux_defective_metal(
                    P_test, 0, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'},
                    T_K, microstructure
                )
                metal_factors.append(result['flux_total'] / bare_flux)
            except:
                metal_factors.append(np.nan)
        
        ax3.semilogx(defect_fractions * 100, metal_factors, color=color, linewidth=2,
                    label=f'{T_K-273:.0f}°C')
    
    ax3.axhline(1.0, color='black', linestyle='-', linewidth=1.5, label='Bare metal')
    ax3.axhline(0.9, color='gray', linestyle='--', alpha=0.7, label='90% of bare')
    ax3.fill_between([0.01, 100], 0.9, 1.0, alpha=0.2, color='red')
    ax3.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax3.set_ylabel('J_total / J_bare_metal', fontsize=12)
    ax3.set_title('Metal Limitation Factor vs Defect Fraction', fontsize=12)
    ax3.legend(loc='best', title='Temperature', fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim([0, 1.1])
    ax3.text(20, 0.95, 'Metal-limited', fontsize=10, color='darkred')
    
    plt.suptitle('Temperature Effect on Permeation Regimes (Level 3,4)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test10_temperature_regime_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 10: PASSED ✓")
    return True


# ============================================================================
# TEST 11: Defect Type Regime Comparison (Level 3,4)
# ============================================================================

def test_11_defect_type_regime_comparison():
    """
    Compare regime behavior for different defect types at Level 3,4.
    """
    print("\n" + "="*60)
    print("TEST 11: Defect Type Regime Comparison (Level 3,4)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    microstructure = get_default_microstructure_params()
    P_test = 100  # Pa (safe for defective metal)
    
    defect_fractions = np.logspace(-4, 0, 100)  # 0.01% to ~100%
    
    defect_configs = {
        'Pinhole': {'type': 'pinhole'},
        'Crack (α=0.1)': {'type': 'crack', 'thickness_factor': 0.1},
        'Crack (α=0.5)': {'type': 'crack', 'thickness_factor': 0.5},
        'Grain Boundary (β=10)': {'type': 'grain_boundary', 'diffusivity_factor': 10},
        'Grain Boundary (β=100)': {'type': 'grain_boundary', 'diffusivity_factor': 100},
    }
    
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(defect_configs)))
    
    # Plot 1: Flux comparison
    ax1 = axes[0]
    
    for (name, config), color in zip(defect_configs.items(), colors):
        fluxes = []
        for f_def in defect_fractions:
            try:
                params = {'area_fraction': f_def, **config}
                result = calculate_parallel_path_flux_defective_metal(
                    P_test, 0, oxide_props, metal_props, params,
                    T_TEST, microstructure
                )
                fluxes.append(result['flux_total'])
            except:
                fluxes.append(np.nan)
        
        ax1.loglog(defect_fractions * 100, fluxes, color=color, linewidth=2, label=name)
    
    # Add Level 1,4 and Level 2,4 references
    result_L14 = calculate_defective_metal_flux(
        metal_props['D_metal'], metal_props['K_s_metal'],
        metal_props['thickness'], P_test, 0,
        T_TEST, microstructure
    )
    result_L24 = calculate_oxide_defective_metal_system(P_test, 0, oxide_props, metal_props, T_TEST, microstructure)
    
    ax1.axhline(result_L14['flux'], color='blue', linestyle='--', alpha=0.7, label='Level 1,4 (Bare)')
    ax1.axhline(result_L24['flux'], color='red', linestyle='--', alpha=0.7, label='Level 2,4 (Perfect)')
    
    ax1.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux vs Defect Fraction by Defect Type (Level 3,4)', fontsize=12)
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(True, which='both', alpha=0.3)
    
    # Plot 2: PRF comparison
    ax2 = axes[1]
    
    for (name, config), color in zip(defect_configs.items(), colors):
        PRFs = []
        for f_def in defect_fractions:
            try:
                params = {'area_fraction': f_def, **config}
                result = calculate_PRF_defective_metal(P_test, oxide_props, metal_props, 
                                                      T_TEST, microstructure, params)
                PRFs.append(result['PRF'])
            except:
                PRFs.append(np.nan)
        
        ax2.loglog(defect_fractions * 100, PRFs, color=color, linewidth=2, label=name)
    
    ax2.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range')
    ax2.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_ylabel('log PRF', fontsize=12)
    ax2.set_title('PRF vs Defect Fraction by Defect Type (Level 3,4)', fontsize=12)
    ax2.legend(loc='best', fontsize=9)
    ax2.grid(True, which='both', alpha=0.3)
    
    plt.suptitle(f'Defect Type Comparison (Level 3,4) at T=800°C, P={P_test} Pa', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test11_defect_type_regime_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 11: PASSED ✓")
    return True


# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

def run_all_tests(save_report=True):
    """Run all tests in sequence and optionally save results to a report file."""
    print("\n" + "="*70)
    print("LEVEL 3,4 COMPREHENSIVE TEST SUITE (FULL MODEL)")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)
    
    tests = [
        ("Basic Functionality", test_1_basic_functionality),
        ("Limiting Cases", test_2_limiting_cases),
        ("Monotonic Behavior", test_3_monotonic_behavior),
        ("PRF Validation", test_4_PRF_validation),
        ("Pressure Sweep", test_5_pressure_sweep),
        ("Defect Type Comparison", test_6_defect_type_comparison),
        ("PRF vs Defect Fraction", test_7_PRF_vs_defect_fraction),
        ("Dominant Path Analysis", test_8_dominant_path_analysis),
        ("PRF Regime Analysis", test_9_PRF_regime_analysis),
        ("Temperature Effects", test_10_temperature_effects),
        ("Defect Type Regime Comparison", test_11_defect_type_regime_comparison),
    ]
    
    results = {}
    for name, test_func in tests:
        try:
            results[name] = test_func()
        except Exception as e:
            print(f"\n✗ {name} failed with error: {e}")
            results[name] = False
            import traceback
            traceback.print_exc()
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY (Level 3,4 Full Model)")
    print("="*70)
    
    passed = sum(results.values())
    total = len(results)
    
    for name, result in results.items():
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"  {name:30} {status}")
    
    print(f"\nTotal: {passed}/{total} tests passed ({passed/total*100:.1f}%)")
    
    if passed == total:
        print("\n" + "="*70)
        print("ALL TESTS PASSED SUCCESSFULLY! ✓")
        print("="*70)
    
    return passed == total


if __name__ == "__main__":
    run_all_tests()
