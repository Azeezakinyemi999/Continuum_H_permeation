"""
Complete Test Suite for Level 4: Defective Metal + Comparison with Level 1
==========================================================================
Tests for hydrogen permeation through metal with microstructure effects:
- Grain boundary enhancement (fast diffusion paths)
- Trapping reduction (dislocations, vacancies, precipitates)

Physics: Modified Fick's Law with D_eff
- D_eff = D_gb_enhanced / (1 + Σθᵢ)
- D_gb_enhanced = (1-f_gb)×D_bulk + f_gb×D_gb
- θᵢ = trap occupancy (Oriani equilibrium)

Comparisons with Level 1 (clean metal) to understand:
- When microstructure enhances vs reduces permeation
- Temperature and grain size effects
- Trade-off between GB enhancement and trapping
"""

import numpy as np
import matplotlib.pyplot as plt
import json
from scipy import stats
from datetime import datetime

# Import core calculation functions
from calculations.permeation_calc import calculate_simple_metal_flux, calculate_defective_metal_flux
from calculations.defective_metal import (
    trap_occupancy,
    grain_boundary_density,
    gb_enhancement_factor,
    calculate_effective_diffusivity_trapping,
    calculate_gb_enhanced_diffusivity,
    combined_microstructure_model
)
from calculations.utils import get_diffusivity, get_solubility, get_permeability
from data.material_data import MATERIALS
from data.microstruture_parameters import (
    LATTICE_PARAMETERS,
    TRAP_PROPERTIES,
    PROCESSING_CONDITIONS,
    GB_ENHANCEMENT_DATA
)


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def save_test_results(test_name, results, results_dir='validation/results/level4'):
    """Save test results to JSON file."""
    import os
    os.makedirs(results_dir, exist_ok=True)
    
    filename = f"{results_dir}/{test_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    def convert_to_serializable(obj):
        """Recursively convert numpy types to JSON-serializable Python types."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.bool_, np.generic)):
            return obj.item()
        elif isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [convert_to_serializable(item) for item in obj]
        elif isinstance(obj, bool):
            return obj
        elif isinstance(obj, (int, float, str, type(None))):
            return obj
        else:
            return str(obj)
    
    json_results = convert_to_serializable(results)
    
    with open(filename, 'w') as f:
        json.dump(json_results, f, indent=2)
    
    print(f"  ✓ Results saved to: {filename}")


def print_test_header(test_name):
    """Print formatted test header."""
    print(f"\n{'='*80}")
    print(f"TEST: {test_name}")
    print(f"{'='*80}")


def print_test_result(passed, message=""):
    """Print test result."""
    status = "✅ PASSED" if passed else "❌ FAILED"
    print(f"{status}: {message}")


def get_default_microstructure(condition='solution_annealed'):
    """
    Get default microstructure parameters for testing.
    
    Parameters
    ----------
    condition : str
        Processing condition: 'solution_annealed', 'cold_worked_20pct', 
        'aged', 'nanocrystalline'
    """
    proc = PROCESSING_CONDITIONS.get(condition, PROCESSING_CONDITIONS['solution_annealed'])
    
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
        'gb_thickness': TRAP_PROPERTIES['grain_boundaries']['thickness'],
        'include_gb_trapping': True
    }
    
    return microstructure


# ============================================================================
# GROUP 1: BASIC DEFECTIVE METAL FLUX TESTS
# ============================================================================

def test_defective_metal_basic_flux(material_name='Incoloy800', temperature_C=800,
                                    P_up=1.0, P_down=0.0, thickness=0.001,
                                    condition='solution_annealed'):
    """
    Test basic flux calculation for defective metal.
    
    Verifies:
    - Flux is positive when P_up > P_down
    - D_eff is calculated and differs from D_lattice
    - All output fields are present
    """
    print_test_header("Defective Metal: Basic Flux Calculation")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D_lattice = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    microstructure = get_default_microstructure(condition)
    
    print(f"\nConditions:")
    print(f"  Material: {material_name}")
    print(f"  Temperature: {temperature_C}°C ({T_K:.1f} K)")
    print(f"  Processing: {condition}")
    print(f"  Grain size: {microstructure['grain_size']*1e6:.1f} μm")
    print(f"  Thickness: {thickness*1000:.2f} mm")
    print(f"  P_up: {P_up} Pa, P_down: {P_down} Pa")
    
    # Calculate defective metal flux
    result = calculate_defective_metal_flux(
        D_lattice=D_lattice,
        K_s=K_s,
        thickness=thickness,
        P_up=P_up,
        P_down=P_down,
        temperature=T_K,
        microstructure_params=microstructure
    )
    
    # Verification checks
    checks = {
        'flux_positive': result['flux'] > 0,
        'D_eff_positive': result['D_eff'] > 0,
        'D_eff_differs_from_lattice': abs(result['D_eff'] - D_lattice) / D_lattice > 1e-6,
        'modification_factor_calculated': result['modification_factor'] > 0,
        'concentrations_valid': result['C_up'] > 0 and result['C_down'] >= 0
    }
    
    all_passed = all(checks.values())
    
    print(f"\nInput Properties:")
    print(f"  D_lattice = {D_lattice:.3e} m²/s")
    print(f"  K_s = {K_s:.3e} mol/m³/Pa^0.5")
    
    print(f"\nDefective Metal Results:")
    print(f"  D_eff = {result['D_eff']:.3e} m²/s")
    print(f"  Modification factor (D_eff/D_lattice) = {result['modification_factor']:.4f}")
    print(f"  Flux = {result['flux']:.3e} mol/m²/s")
    print(f"  C_up = {result['C_up']:.3e} mol/m³")
    print(f"  C_down = {result['C_down']:.3e} mol/m³")
    
    if result['modification_factor'] > 1:
        print(f"\n  → GB enhancement dominates (D_eff > D_lattice)")
    else:
        print(f"\n  → Trapping dominates (D_eff < D_lattice)")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    results_dict = {
        'test_name': 'defective_metal_basic_flux',
        'material': material_name,
        'temperature_C': temperature_C,
        'condition': condition,
        'D_lattice': D_lattice,
        'D_eff': result['D_eff'],
        'modification_factor': result['modification_factor'],
        'flux': result['flux'],
        'checks': checks,
        'passed': all_passed
    }
    
    save_test_results('defective_basic_flux', results_dict)
    return all_passed, results_dict


def test_comparison_level1_vs_level4(material_name='Incoloy800', temperature_C=800,
                                     P_up=1.0, P_down=0.0, thickness=0.001,
                                     condition='solution_annealed'):
    """
    Compare Level 1 (clean metal) vs Level 4 (defective metal) at same conditions.
    
    Verifies:
    - Both models run successfully
    - Flux ratio is physically reasonable (0.01 to 100×)
    - Results can be explained by microstructure
    """
    print_test_header("Comparison: Level 1 (Clean) vs Level 4 (Defective)")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D_lattice = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    # Level 1: Clean metal
    result_L1 = calculate_simple_metal_flux(D_lattice, K_s, thickness, P_up, P_down)
    
    # Level 4: Defective metal
    microstructure = get_default_microstructure(condition)
    result_L4 = calculate_defective_metal_flux(
        D_lattice=D_lattice,
        K_s=K_s,
        thickness=thickness,
        P_up=P_up,
        P_down=P_down,
        temperature=T_K,
        microstructure_params=microstructure
    )
    
    # Calculate ratios
    flux_ratio = result_L4['flux'] / result_L1['flux']
    D_ratio = result_L4['D_eff'] / D_lattice
    
    # Verification
    checks = {
        'L1_flux_positive': result_L1['flux'] > 0,
        'L4_flux_positive': result_L4['flux'] > 0,
        'flux_ratio_reasonable': 0.01 < flux_ratio < 100,
        'D_ratio_matches_flux_ratio': abs(D_ratio - flux_ratio) / flux_ratio < 0.01
    }
    
    all_passed = all(checks.values())
    
    print(f"\nConditions:")
    print(f"  Material: {material_name}")
    print(f"  Temperature: {temperature_C}°C")
    print(f"  Processing: {condition}")
    print(f"  Grain size: {microstructure['grain_size']*1e6:.1f} μm")
    
    print(f"\nLevel 1 (Clean Metal):")
    print(f"  D = {D_lattice:.3e} m²/s")
    print(f"  Flux = {result_L1['flux']:.3e} mol/m²/s")
    
    print(f"\nLevel 4 (Defective Metal):")
    print(f"  D_eff = {result_L4['D_eff']:.3e} m²/s")
    print(f"  Flux = {result_L4['flux']:.3e} mol/m²/s")
    
    print(f"\nComparison:")
    print(f"  Flux ratio (L4/L1) = {flux_ratio:.4f}")
    print(f"  D ratio (D_eff/D_lattice) = {D_ratio:.4f}")
    
    if flux_ratio > 1:
        print(f"\n  → Level 4 flux is {flux_ratio:.2f}× HIGHER (GB enhancement dominates)")
    else:
        print(f"\n  → Level 4 flux is {1/flux_ratio:.2f}× LOWER (trapping dominates)")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    results_dict = {
        'test_name': 'comparison_L1_vs_L4',
        'material': material_name,
        'temperature_C': temperature_C,
        'condition': condition,
        'L1_flux': result_L1['flux'],
        'L4_flux': result_L4['flux'],
        'flux_ratio': flux_ratio,
        'D_ratio': D_ratio,
        'checks': checks,
        'passed': all_passed
    }
    
    save_test_results('comparison_L1_L4', results_dict)
    return all_passed, results_dict


# ============================================================================
# GROUP 2: TRAP OCCUPANCY TESTS
# ============================================================================

def test_trap_occupancy_temperature_dependence(binding_energy=27e3, trap_density=1e15,
                                               lattice_density=1.06e29,
                                               lattice_concentration=1e-3):
    """
    Test trap occupancy varies correctly with temperature.
    
    Oriani equilibrium: θ = K×θ_L / (1 + K×θ_L) where K = exp(E_b/RT)
    As T increases, θ should decrease (traps release hydrogen)
    """
    print_test_header("Trap Occupancy: Temperature Dependence")
    
    R = 8.314  # J/mol/K
    
    temperatures_C = np.array([400, 500, 600, 700, 800, 900, 1000])
    temperatures_K = temperatures_C + 273.15
    
    occupancies = []
    K_values = []
    
    for T_K in temperatures_K:
        result = trap_occupancy(
            temperature=T_K,
            binding_energy=binding_energy,
            trap_density=trap_density,
            lattice_density=lattice_density,
            lattice_concentration=lattice_concentration
        )
        occupancies.append(result['theta'])
        K_values.append(result['K_equilibrium'])
    
    occupancies = np.array(occupancies)
    K_values = np.array(K_values)
    
    # Check monotonically decreasing
    is_monotonic = np.all(np.diff(occupancies) < 0)
    
    # Check Arrhenius behavior of K
    log_K = np.log(K_values)
    inv_T = 1000 / temperatures_K
    slope, intercept, r_value, _, _ = stats.linregress(inv_T, log_K)
    
    E_b_extracted = slope * R * 1000  # Convert back to J/mol
    E_b_error = abs(E_b_extracted - binding_energy) / binding_energy
    
    checks = {
        'monotonic_decrease': is_monotonic,
        'arrhenius_K': r_value**2 > 0.999,
        'E_b_extracted_correct': E_b_error < 0.01
    }
    
    all_passed = all(checks.values())
    
    print(f"\nTrap Parameters:")
    print(f"  Binding energy: {binding_energy/1000:.1f} kJ/mol")
    print(f"  Trap density: {trap_density:.1e} m⁻³")
    print(f"  Lattice concentration: {lattice_concentration:.1e} mol/m³")
    
    print(f"\nOccupancy vs Temperature:")
    print(f"{'T (°C)':>8} {'θ':>12} {'K_eq':>12}")
    print(f"{'-'*35}")
    for i, T_C in enumerate(temperatures_C):
        print(f"{T_C:8.0f} {occupancies[i]:12.4e} {K_values[i]:12.2e}")
    
    print(f"\nArrhenius Analysis:")
    print(f"  Extracted E_b = {E_b_extracted/1000:.2f} kJ/mol (input: {binding_energy/1000:.1f})")
    print(f"  Error: {E_b_error*100:.2f}%")
    print(f"  R² = {r_value**2:.6f}")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    # Create plot - both panels in Arrhenius format
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Trap Occupancy Analysis', fontsize=14, fontweight='bold')
    
    ax1 = axes[0]
    ax1.semilogy(inv_T, occupancies, 'o-', markersize=8, color='blue')
    ax1.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax1.set_ylabel('log Trap Occupancy θ', fontsize=11)
    ax1.set_title(f'E_b = {binding_energy/1000:.0f} kJ/mol')
    ax1.grid(True, which='both', alpha=0.3)
    
    # Add secondary x-axis with temperature in °C
    ax1_top = ax1.twiny()
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(inv_T[::2])  # Show every other temperature to avoid crowding
    ax1_top.set_xticklabels([f"{int(temperatures_C[i])}" for i in range(0, len(temperatures_C), 2)])
    ax1_top.set_xlabel('Temperature (°C)', fontsize=11)
    
    ax2 = axes[1]
    ax2.semilogy(inv_T, K_values, 'o', markersize=8, color='red', label='Data')
    fit_K = np.exp(slope * inv_T + intercept)
    ax2.semilogy(inv_T, fit_K, '--', linewidth=2, color='darkred', label='Fit')
    ax2.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax2.set_ylabel('log K (equilibrium constant)', fontsize=11)
    ax2.set_title(f'E_b extracted = {E_b_extracted/1000:.1f} kJ/mol')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.3)
    
    plt.tight_layout()
    
    import os
    os.makedirs('validation/results/level4', exist_ok=True)
    plot_file = f"validation/results/level4/trap_occupancy_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    return all_passed, {
        'temperatures_C': temperatures_C,
        'occupancies': occupancies,
        'E_b_extracted': E_b_extracted,
        'checks': checks,
        'passed': all_passed
    }


def test_trap_occupancy_binding_energy_effect(temperature=1073.15,
                                              lattice_density=1.06e29,
                                              lattice_concentration=1e-3):
    """
    Test trap occupancy increases with binding energy at fixed temperature.
    Higher binding energy → more stable trap → higher occupancy
    """
    print_test_header("Trap Occupancy: Binding Energy Effect")
    
    trap_density = 1e15  # m⁻³
    
    trap_types = {
        'dislocations': 27e3,
        'vacancies': 41e3,
        'grain_boundaries': 48e3,
        'carbides': 72e3
    }
    
    occupancies = {}
    for name, E_b in trap_types.items():
        result = trap_occupancy(
            temperature=temperature,
            binding_energy=E_b,
            trap_density=trap_density,
            lattice_density=lattice_density,
            lattice_concentration=lattice_concentration
        )
        occupancies[name] = result['theta']
    
    E_b_sorted = sorted(trap_types.items(), key=lambda x: x[1])
    theta_sorted = [occupancies[name] for name, _ in E_b_sorted]
    
    is_monotonic = all(theta_sorted[i] < theta_sorted[i+1] 
                      for i in range(len(theta_sorted)-1))
    
    checks = {
        'higher_Eb_higher_theta': is_monotonic,
        'all_theta_positive': all(θ > 0 for θ in occupancies.values()),
        'all_theta_less_than_1': all(θ < 1 for θ in occupancies.values())
    }
    
    all_passed = all(checks.values())
    
    T_C = temperature - 273.15
    print(f"\nConditions:")
    print(f"  Temperature: {T_C:.0f}°C ({temperature:.1f} K)")
    print(f"  Trap density: {trap_density:.1e} m⁻³")
    
    print(f"\nOccupancy by Trap Type:")
    print(f"{'Trap Type':<20} {'E_b (kJ/mol)':>15} {'θ':>15}")
    print(f"{'-'*55}")
    for name, E_b in sorted(trap_types.items(), key=lambda x: x[1]):
        print(f"{name:<20} {E_b/1000:15.1f} {occupancies[name]:15.4e}")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    return all_passed, {
        'trap_types': trap_types,
        'occupancies': occupancies,
        'checks': checks,
        'passed': all_passed
    }


# ============================================================================
# GROUP 3: GRAIN BOUNDARY ENHANCEMENT TESTS
# ============================================================================

# def test_gb_enhancement_factor(material_name='Incoloy800'):
#     """
#     Test GB enhancement factor decreases with temperature.
#     Theory: D_gb/D_bulk = A × exp[(E_bulk - E_gb)/RT]
#     """
#     print_test_header("GB Enhancement Factor vs Temperature")
    
#     temperatures_C = np.array([500, 600, 700, 800, 900, 1000])
#     temperatures_K = temperatures_C + 273.15
    
#     material = MATERIALS[material_name]
    
#     enhancement_factors = []
#     D_lattice_values = []
#     D_gb_values = []
    
#     for T_K in temperatures_K:
#         D_lattice = get_diffusivity(T_K, material)
        
#         result = gb_enhancement_factor(
#             temperature=T_K,
#             D_bulk=D_lattice,
#             activation_ratio=0.6
#         )
        
#         enhancement_factors.append(result['enhancement_ratio'])
#         D_lattice_values.append(D_lattice)
#         D_gb_values.append(result['D_gb'])
    
#     enhancement_factors = np.array(enhancement_factors)
    
#     is_monotonic = np.all(np.diff(enhancement_factors) < 0)
#     all_reasonable = np.all(enhancement_factors > 1) and np.all(enhancement_factors < 1e6)
    
#     checks = {
#         'monotonic_decrease': is_monotonic,
#         'values_reasonable': all_reasonable,
#         'always_greater_than_1': np.all(enhancement_factors > 1)
#     }
    
#     all_passed = all(checks.values())
    
#     print(f"\nGB Enhancement vs Temperature:")
#     print(f"{'T (°C)':>8} {'D_lattice (m²/s)':>18} {'D_gb (m²/s)':>15} {'α = D_gb/D_bulk':>18}")
#     print(f"{'-'*65}")
#     for i, T_C in enumerate(temperatures_C):
#         print(f"{T_C:8.0f} {D_lattice_values[i]:18.3e} {D_gb_values[i]:15.3e} {enhancement_factors[i]:18.1f}")
    
#     print(f"\nVerification:")
#     for check, passed in checks.items():
#         print_test_result(passed, check)
    
#     # Create plot
#     fig, ax = plt.subplots(figsize=(10, 6))
#     ax.semilogy(temperatures_C, enhancement_factors, 'o-', markersize=10, 
#                 color='green', linewidth=2)
#     ax.set_xlabel('Temperature (°C)', fontsize=12)
#     ax.set_ylabel('log Enhancement Factor (D_gb/D_bulk)', fontsize=12)
#     ax.set_title('Grain Boundary Diffusion Enhancement', fontsize=14)
#     ax.grid(True, which='both', alpha=0.3)
    
#     plt.tight_layout()
    
#     import os
#     os.makedirs('validation/results/level4', exist_ok=True)
#     plot_file = f"validation/results/level4/gb_enhancement_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
#     plt.savefig(plot_file, dpi=300, bbox_inches='tight')
#     print(f"  ✓ Plot saved: {plot_file}")
#     plt.close()
    
#     return all_passed, {
#         'temperatures_C': temperatures_C,
#         'enhancement_factors': enhancement_factors,
#         'checks': checks,
#         'passed': all_passed
#     }

def test_gb_enhancement_factor(material_name='Incoloy800'):
    """
    Test GB enhancement factor decreases with temperature.
    Theory: D_gb/D_bulk = A × exp[(E_bulk - E_gb)/RT]
    """
    print_test_header("GB Enhancement Factor vs Temperature")
    
    temperatures_C = np.array([600, 700, 800, 900, 1000])
    temperatures_K = temperatures_C + 273.15
    
    material = MATERIALS[material_name]
    
    enhancement_factors = []
    
    for T_K in temperatures_K:
        result = gb_enhancement_factor(
            temperature=T_K,
            temperature_unit='K',
            gb_type='HAGB'
        )
        
        enhancement_factors.append(result['enhancement_factor'])
    
    enhancement_factors = np.array(enhancement_factors)
    
    is_monotonic = np.all(np.diff(enhancement_factors) < 0)
    all_reasonable = np.all(enhancement_factors > 1) and np.all(enhancement_factors < 1e6)
    
    checks = {
        'monotonic_decrease': is_monotonic,
        'values_reasonable': all_reasonable,
        'always_greater_than_1': np.all(enhancement_factors > 1)
    }
    
    all_passed = all(checks.values())
    
    print(f"\nGB Enhancement vs Temperature:")
    print(f"{'T (°C)':>8} {'Enhancement Factor (D_gb/D_bulk)':>35}")
    print(f"{'-'*45}")
    for i, T_C in enumerate(temperatures_C):
        print(f"{T_C:8.0f} {enhancement_factors[i]:35.1f}")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    # Create plot - Arrhenius format (1000/T on x-axis)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    inv_T = 1000 / temperatures_K
    ax.semilogy(inv_T, enhancement_factors, 'o-', markersize=10, 
                color='green', linewidth=2)
    ax.set_xlabel('1000/T (K⁻¹)', fontsize=12)
    ax.set_ylabel('log Enhancement Factor (D_gb/D_bulk)', fontsize=12)
    ax.set_title('Grain Boundary Diffusion Enhancement (Arrhenius)', fontsize=14)
    ax.grid(True, which='both', alpha=0.3)
    
    # Add secondary x-axis with temperature in °C
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(1000/(temperatures_K))
    ax_top.set_xticklabels([f"{int(t)}" for t in temperatures_C])
    ax_top.set_xlabel('Temperature (°C)', fontsize=12)
    
    plt.tight_layout()
    
    import os
    os.makedirs('validation/results/level4', exist_ok=True)
    plot_file = f"validation/results/level4/gb_enhancement_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    return all_passed, {
        'temperatures_C': temperatures_C,
        'enhancement_factors': enhancement_factors,
        'checks': checks,
        'passed': all_passed
    }


def test_grain_size_effect_on_flux(material_name='Incoloy800', temperature_C=800,
                                   P_up=1.0, P_down=0.0, thickness=0.001):
    """
    Test effect of grain size on permeation flux.
    Smaller grains → more GB area → more fast paths → higher flux
    (if GB enhancement dominates over increased trapping)
    """
    print_test_header("Grain Size Effect on Flux")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D_lattice = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    result_L1 = calculate_simple_metal_flux(D_lattice, K_s, thickness, P_up, P_down)
    flux_L1 = result_L1['flux']
    
    grain_sizes = np.array([1e-6, 5e-6, 10e-6, 25e-6, 50e-6, 100e-6, 500e-6])
    
    fluxes = []
    D_effs = []
    flux_ratios = []
    
    for d_grain in grain_sizes:
        microstructure = get_default_microstructure()
        microstructure['grain_size'] = d_grain
        
        result = calculate_defective_metal_flux(
            D_lattice=D_lattice,
            K_s=K_s,
            thickness=thickness,
            P_up=P_up,
            P_down=P_down,
            temperature=T_K,
            microstructure_params=microstructure
        )
        
        fluxes.append(result['flux'])
        D_effs.append(result['D_eff'])
        flux_ratios.append(result['flux'] / flux_L1)
    
    fluxes = np.array(fluxes)
    D_effs = np.array(D_effs)
    flux_ratios = np.array(flux_ratios)
    
    is_monotonic = (np.all(np.diff(fluxes) > 0) or np.all(np.diff(fluxes) < 0))
    
    checks = {
        'trend_monotonic': is_monotonic,
        'all_fluxes_positive': np.all(fluxes > 0),
        'significant_variation': np.max(fluxes) / np.min(fluxes) > 1.1
    }
    
    all_passed = all(checks.values())
    
    print(f"\nConditions:")
    print(f"  Temperature: {temperature_C}°C")
    print(f"  Level 1 flux: {flux_L1:.3e} mol/m²/s")
    
    print(f"\nGrain Size Effect:")
    print(f"{'d_grain (μm)':>15} {'D_eff (m²/s)':>15} {'Flux (mol/m²/s)':>18} {'Ratio (L4/L1)':>15}")
    print(f"{'-'*68}")
    for i, d in enumerate(grain_sizes):
        print(f"{d*1e6:15.1f} {D_effs[i]:15.3e} {fluxes[i]:18.3e} {flux_ratios[i]:15.4f}")
    
    if fluxes[0] > fluxes[-1]:
        print(f"\n  → Smaller grains increase flux (GB enhancement dominates)")
    else:
        print(f"\n  → Smaller grains decrease flux (GB trapping dominates)")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f'Grain Size Effect at {temperature_C}°C', fontsize=14, fontweight='bold')
    
    ax1 = axes[0]
    ax1.loglog(grain_sizes*1e6, fluxes, 'o-', markersize=8, color='blue', label='Level 4')
    ax1.axhline(y=flux_L1, color='red', linestyle='--', linewidth=2, label='Level 1 (clean)')
    ax1.set_xlabel('Grain Size (μm)', fontsize=11)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=11)
    ax1.set_title('Permeation Flux vs Grain Size')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)
    
    ax2 = axes[1]
    ax2.semilogx(grain_sizes*1e6, flux_ratios, 's-', markersize=8, color='green')
    ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Grain Size (μm)', fontsize=11)
    ax2.set_ylabel('Flux Ratio (Level 4 / Level 1)', fontsize=11)
    ax2.set_title('Microstructure Effect')
    ax2.fill_between(grain_sizes*1e6, 1, flux_ratios, where=flux_ratios > 1,
                     alpha=0.3, color='green', label='Enhancement')
    ax2.fill_between(grain_sizes*1e6, flux_ratios, 1, where=flux_ratios < 1,
                     alpha=0.3, color='red', label='Reduction')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    import os
    os.makedirs('validation/results/level4', exist_ok=True)
    plot_file = f"validation/results/level4/grain_size_effect_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    return all_passed, {
        'grain_sizes': grain_sizes,
        'fluxes': fluxes,
        'flux_ratios': flux_ratios,
        'checks': checks,
        'passed': all_passed
    }


# ============================================================================
# GROUP 4: TEMPERATURE DEPENDENCE COMPARISON (L1 vs L4)
# ============================================================================

def test_temperature_dependence_comparison(material_name='Incoloy800',
                                           T_min_C=500, T_max_C=1000, n_points=11,
                                           P_upstream=1.0, thickness=0.001,
                                           condition='solution_annealed'):
    """
    Compare temperature dependence of Level 1 vs Level 4.
    Both should show Arrhenius behavior, but Level 4 may have different E_app.
    """
    print_test_header("Temperature Dependence: Level 1 vs Level 4 Comparison")
    
    material = MATERIALS[material_name]
    microstructure = get_default_microstructure(condition)
    
    temperatures_C = np.linspace(T_min_C, T_max_C, n_points)
    temperatures_K = temperatures_C + 273.15
    
    flux_L1 = []
    flux_L4 = []
    D_lattice_arr = []
    D_eff_arr = []
    
    for T_K in temperatures_K:
        D_lattice = get_diffusivity(T_K, material)
        K_s = get_solubility(T_K, material)
        
        result_L1 = calculate_simple_metal_flux(D_lattice, K_s, thickness, P_upstream, 0)
        flux_L1.append(result_L1['flux'])
        D_lattice_arr.append(D_lattice)
        
        result_L4 = calculate_defective_metal_flux(
            D_lattice=D_lattice,
            K_s=K_s,
            thickness=thickness,
            P_up=P_upstream,
            P_down=0,
            temperature=T_K,
            microstructure_params=microstructure
        )
        flux_L4.append(result_L4['flux'])
        D_eff_arr.append(result_L4['D_eff'])
    
    flux_L1 = np.array(flux_L1)
    flux_L4 = np.array(flux_L4)
    D_lattice_arr = np.array(D_lattice_arr)
    D_eff_arr = np.array(D_eff_arr)
    
    flux_ratio = flux_L4 / flux_L1
    
    R = 8.314
    x_data = 1000 / temperatures_K
    
    slope_L1, intercept_L1, r_L1, _, _ = stats.linregress(x_data, np.log(flux_L1))
    E_app_L1 = -slope_L1 * R * 1000
    
    slope_L4, intercept_L4, r_L4, _, _ = stats.linregress(x_data, np.log(flux_L4))
    E_app_L4 = -slope_L4 * R * 1000
    
    checks = {
        'L1_arrhenius': r_L1**2 > 0.99,
        'L4_arrhenius': r_L4**2 > 0.99,
        'E_app_differs': abs(E_app_L4 - E_app_L1) / E_app_L1 > 0.01,
        'ratio_varies_with_T': np.std(flux_ratio) / np.mean(flux_ratio) > 0.01
    }
    
    all_passed = all(checks.values())
    
    print(f"\nConditions:")
    print(f"  Material: {material_name}")
    print(f"  Temperature range: {T_min_C}-{T_max_C}°C")
    print(f"  Processing: {condition}")
    
    print(f"\nApparent Activation Energies:")
    print(f"  Level 1: E_app = {E_app_L1/1000:.2f} kJ/mol (R² = {r_L1**2:.6f})")
    print(f"  Level 4: E_app = {E_app_L4/1000:.2f} kJ/mol (R² = {r_L4**2:.6f})")
    print(f"  Difference: {(E_app_L4 - E_app_L1)/1000:.2f} kJ/mol")
    
    print(f"\nFlux Ratio Range:")
    print(f"  Min (L4/L1): {flux_ratio.min():.4f} at {temperatures_C[np.argmin(flux_ratio)]:.0f}°C")
    print(f"  Max (L4/L1): {flux_ratio.max():.4f} at {temperatures_C[np.argmax(flux_ratio)]:.0f}°C")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Temperature Comparison: Level 1 vs Level 4 ({material_name})', 
                 fontsize=14, fontweight='bold')
    
    ax1 = axes[0, 0]
    ax1.semilogy(x_data, flux_L1, 'o-', markersize=6, color='blue', label='Level 1 (Clean)')
    ax1.semilogy(x_data, flux_L4, 's-', markersize=6, color='red', label='Level 4 (Defective)')
    ax1.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=11)
    ax1.set_title('Flux Arrhenius Plot')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)
    
    ax1_top = ax1.twiny()
    temp_ticks = np.array([500, 600, 700, 800, 900, 1000])
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(1000/(temp_ticks + 273.15))
    ax1_top.set_xticklabels([f"{t}" for t in temp_ticks])
    ax1_top.set_xlabel('Temperature (°C)', fontsize=11)
    
    ax2 = axes[0, 1]
    ax2.plot(temperatures_C, flux_ratio, 'go-', markersize=8, linewidth=2)
    ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='L1 = L4')
    ax2.set_xlabel('Temperature (°C)', fontsize=11)
    ax2.set_ylabel('Flux Ratio (Level 4 / Level 1)', fontsize=11)
    ax2.set_title('Microstructure Effect vs Temperature')
    ax2.fill_between(temperatures_C, 1, flux_ratio, where=flux_ratio > 1,
                     alpha=0.3, color='green', label='Enhancement')
    ax2.fill_between(temperatures_C, flux_ratio, 1, where=flux_ratio < 1,
                     alpha=0.3, color='red', label='Reduction')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    ax3 = axes[1, 0]
    ax3.semilogy(x_data, D_lattice_arr, 'o-', markersize=6, color='blue', label='D_lattice')
    ax3.semilogy(x_data, D_eff_arr, 's-', markersize=6, color='red', label='D_eff')
    ax3.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax3.set_ylabel('log Diffusivity (m²/s)', fontsize=11)
    ax3.set_title('Diffusivity Comparison')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)
    
    ax4 = axes[1, 1]
    mod_factor = D_eff_arr / D_lattice_arr
    ax4.plot(temperatures_C, mod_factor, 'mo-', markersize=8, linewidth=2)
    ax4.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Temperature (°C)', fontsize=11)
    ax4.set_ylabel('D_eff / D_lattice', fontsize=11)
    ax4.set_title('Diffusivity Modification Factor')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    import os
    os.makedirs('validation/results/level4', exist_ok=True)
    plot_file = f"validation/results/level4/temperature_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    results_dict = {
        'test_name': 'temperature_comparison_L1_L4',
        'material': material_name,
        'temperatures_C': temperatures_C,
        'flux_L1': flux_L1,
        'flux_L4': flux_L4,
        'flux_ratio': flux_ratio,
        'E_app_L1': E_app_L1,
        'E_app_L4': E_app_L4,
        'checks': checks,
        'passed': all_passed
    }
    
    save_test_results('temperature_comparison', results_dict)
    return all_passed, results_dict


# ============================================================================
# GROUP 5: PROCESSING CONDITION COMPARISON
# ============================================================================

def test_processing_conditions_comparison(material_name='Incoloy800', temperature_C=800,
                                          P_up=1.0, P_down=0.0, thickness=0.001):
    """
    Compare different processing conditions (annealed, cold-worked, etc.).
    """
    print_test_header("Processing Condition Comparison")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D_lattice = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    result_L1 = calculate_simple_metal_flux(D_lattice, K_s, thickness, P_up, P_down)
    flux_L1 = result_L1['flux']
    
    conditions = ['solution_annealed', 'cold_worked_20pct', 'aged', 'nanocrystalline']
    
    results_by_condition = {}
    
    for condition in conditions:
        microstructure = get_default_microstructure(condition)
        
        result = calculate_defective_metal_flux(
            D_lattice=D_lattice,
            K_s=K_s,
            thickness=thickness,
            P_up=P_up,
            P_down=P_down,
            temperature=T_K,
            microstructure_params=microstructure
        )
        
        results_by_condition[condition] = {
            'flux': result['flux'],
            'D_eff': result['D_eff'],
            'modification_factor': result['modification_factor'],
            'flux_ratio': result['flux'] / flux_L1,
            'grain_size': microstructure['grain_size']
        }
    
    checks = {
        'all_positive_flux': all(r['flux'] > 0 for r in results_by_condition.values()),
        'variation_exists': (max(r['flux'] for r in results_by_condition.values()) / 
                            min(r['flux'] for r in results_by_condition.values()) > 1.1),
    }
    
    all_passed = all(checks.values())
    
    print(f"\nConditions:")
    print(f"  Material: {material_name}")
    print(f"  Temperature: {temperature_C}°C")
    print(f"  Level 1 flux: {flux_L1:.3e} mol/m²/s")
    
    print(f"\nResults by Processing Condition:")
    print(f"{'Condition':<20} {'d_grain (μm)':>12} {'D_eff/D_lat':>12} {'Flux (mol/m²/s)':>18} {'L4/L1':>10}")
    print(f"{'-'*75}")
    
    for condition in conditions:
        r = results_by_condition[condition]
        print(f"{condition:<20} {r['grain_size']*1e6:12.1f} {r['modification_factor']:12.4f} "
              f"{r['flux']:18.3e} {r['flux_ratio']:10.4f}")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    # Create bar chart
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f'Processing Condition Effects at {temperature_C}°C', fontsize=14, fontweight='bold')
    
    ax1 = axes[0]
    flux_values = [results_by_condition[c]['flux'] for c in conditions]
    colors = ['blue', 'orange', 'green', 'red']
    bars = ax1.bar(conditions, flux_values, color=colors, alpha=0.7)
    ax1.axhline(y=flux_L1, color='black', linestyle='--', linewidth=2, label='Level 1 (clean)')
    ax1.set_ylabel('Flux (mol/m²/s)', fontsize=11)
    ax1.set_title('Permeation Flux')
    ax1.legend()
    ax1.tick_params(axis='x', rotation=45)
    
    ax2 = axes[1]
    ratio_values = [results_by_condition[c]['flux_ratio'] for c in conditions]
    bars = ax2.bar(conditions, ratio_values, color=colors, alpha=0.7)
    ax2.axhline(y=1, color='black', linestyle='--', linewidth=2)
    ax2.set_ylabel('Flux Ratio (Level 4 / Level 1)', fontsize=11)
    ax2.set_title('Microstructure Effect')
    ax2.tick_params(axis='x', rotation=45)
    
    for bar, val in zip(bars, ratio_values):
        ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    import os
    os.makedirs('validation/results/level4', exist_ok=True)
    plot_file = f"validation/results/level4/processing_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    return all_passed, {
        'conditions': conditions,
        'results': results_by_condition,
        'flux_L1': flux_L1,
        'checks': checks,
        'passed': all_passed
    }


# ============================================================================
# GROUP 6: SUMMARY REPORT
# ============================================================================

def generate_level4_summary_report(test_results, material_name='Incoloy800'):
    """Generate comprehensive summary report for Level 4 tests."""
    print_test_header("Generate Level 4 Summary Report")
    
    import os
    os.makedirs('validation/results/level4', exist_ok=True)
    
    material = MATERIALS[material_name]
    
    report_file = f"validation/results/level4/Level4_Complete_Summary_{material_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("LEVEL 4 MODEL: COMPLETE TEST REPORT\n")
        f.write("Defective Metal - Microstructure Effects\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Date/Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Material: {material_name}\n")
        f.write(f"Test Suite: complete_test_level1_4_defective_metal.py\n\n")
        
        f.write("-"*80 + "\n")
        f.write("1. PHYSICS SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write("Level 4 incorporates two competing microstructure effects:\n\n")
        f.write("a) Grain Boundary Enhancement:\n")
        f.write("   D_gb_enhanced = (1-f_gb)×D_bulk + f_gb×D_gb\n")
        f.write("   where f_gb = 3δ/d (GB volume fraction)\n\n")
        f.write("b) Trapping Reduction:\n")
        f.write("   D_eff = D_gb_enhanced / (1 + Σθᵢ)\n")
        f.write("   θᵢ = trap occupancy from Oriani equilibrium\n\n")
        
        f.write("-"*80 + "\n")
        f.write("2. TEST RESULTS SUMMARY\n")
        f.write("-"*80 + "\n\n")
        
        total_tests = len(test_results)
        passed_tests = sum(test_results.values())
        pass_rate = 100 * passed_tests / total_tests if total_tests > 0 else 0
        
        f.write(f"Total tests: {total_tests}\n")
        f.write(f"Passed: {passed_tests}\n")
        f.write(f"Pass rate: {pass_rate:.1f}%\n\n")
        
        for test_name, passed in test_results.items():
            status = "PASS" if passed else "FAIL"
            f.write(f"  [{status}] {test_name}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"  ✓ Report saved to: {report_file}")
    
    return True, report_file


# ============================================================================
# MAIN RUNNER
# ============================================================================

def run_all_level4_tests():
    """Run all Level 4 (defective metal) tests."""
    print("\n" + "="*80)
    print("LEVEL 4 TEST SUITE: DEFECTIVE METAL")
    print("Microstructure Effects: GB Enhancement + Trapping")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    results = {}
    
    # Group 1
    print("\n" + "▶"*40)
    print("GROUP 1: BASIC DEFECTIVE METAL TESTS")
    print("▶"*40)
    
    passed, data = test_defective_metal_basic_flux()
    results['defective_basic_flux'] = passed
    
    passed, data = test_comparison_level1_vs_level4()
    results['L1_vs_L4_comparison'] = passed
    
    # Group 2
    print("\n" + "▶"*40)
    print("GROUP 2: TRAP OCCUPANCY TESTS")
    print("▶"*40)
    
    passed, data = test_trap_occupancy_temperature_dependence()
    results['trap_occupancy_T'] = passed
    
    passed, data = test_trap_occupancy_binding_energy_effect()
    results['trap_occupancy_Eb'] = passed
    
    # Group 3
    print("\n" + "▶"*40)
    print("GROUP 3: GRAIN BOUNDARY ENHANCEMENT TESTS")
    print("▶"*40)
    
    passed, data = test_gb_enhancement_factor()
    results['gb_enhancement'] = passed
    
    passed, data = test_grain_size_effect_on_flux()
    results['grain_size_effect'] = passed
    
    # Group 4
    print("\n" + "▶"*40)
    print("GROUP 4: TEMPERATURE DEPENDENCE COMPARISON")
    print("▶"*40)
    
    passed, data = test_temperature_dependence_comparison()
    results['temperature_comparison'] = passed
    
    # Group 5
    print("\n" + "▶"*40)
    print("GROUP 5: PROCESSING CONDITION COMPARISON")
    print("▶"*40)
    
    passed, data = test_processing_conditions_comparison()
    results['processing_comparison'] = passed
    
    # Group 6
    print("\n" + "▶"*40)
    print("GROUP 6: SUMMARY REPORT")
    print("▶"*40)
    
    passed, report_file = generate_level4_summary_report(results)
    results['summary_report'] = passed
    
    # Final Summary
    print("\n" + "="*80)
    print("FINAL TEST SUMMARY - LEVEL 4 COMPLETE")
    print("="*80)
    
    total_tests = len(results)
    passed_tests = sum(results.values())
    
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {status}: {test_name}")
    
    print(f"\n{'-'*80}")
    print(f"TOTAL: {passed_tests}/{total_tests} tests passed ({100*passed_tests/total_tests:.1f}%)")
    
    if passed_tests == total_tests:
        print("\n🎉 ALL TESTS PASSED! Level 4 model fully validated.")
        print(f"📄 Full report: {report_file}")
    
    print("="*80)
    
    return results


if __name__ == "__main__":
    results = run_all_level4_tests()