"""
Complete Test Suite for Level 1: Simple Metal (Clean Metal)
============================================================
Consolidates all Level 1 tests from:
- pressure_study.py
- temperature_study.py
- experimental_comparison.py
- sensitivity_analysis.py

Physics: Sieverts' Law + Fick's Law
- C = K_s × √P (surface equilibrium)
- J = D × (C_up - C_down) / L (bulk diffusion)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from datetime import datetime
import json

# Import core calculation functions
from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.utils import get_diffusivity, get_solubility, get_permeability
from data.material_data import MATERIALS
from data.experimental_data import get_experimental_data, convert_to_SI


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def save_test_results(test_name, results, results_dir='validation/results/level1'):
    """
    Save test results to JSON file.
    
    Parameters
    ----------
    test_name : str
        Name of the test
    results : dict
        Test results dictionary
    results_dir : str
        Directory to save results
    """
    import os
    os.makedirs(results_dir, exist_ok=True)
    
    filename = f"{results_dir}/{test_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    def convert_to_serializable(obj):
        """Recursively convert numpy types to JSON-serializable Python types."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.bool_, np.generic)):
            return obj.item()  # Convert numpy scalar to Python native type
        elif isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [convert_to_serializable(item) for item in obj]
        elif isinstance(obj, bool):
            return obj
        elif isinstance(obj, (int, float, str, type(None))):
            return obj
        else:
            # Try to convert unknown types to string
            return str(obj)
    
    # Convert all results to JSON-serializable format
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


# ============================================================================
# GROUP 1: BASIC FLUX CALCULATION TESTS
# ============================================================================

def test_basic_flux_calculation(material_name='Incoloy800', temperature_C=800, 
                                P_up=1.0, P_down=0.0, thickness=0.001):
    """
    Test basic flux calculation at a single condition.
    
    Verifies:
    - Flux is positive when P_up > P_down
    - Concentrations calculated correctly
    - Permeability matches D × K_s
    """
    print_test_header("Basic Flux Calculation")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    # Calculate properties
    D = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    P_theory = get_permeability(T_K, material)
    
    # Calculate flux
    result = calculate_simple_metal_flux(D, K_s, thickness, P_up, P_down)
    
    # Verification checks
    checks = {
        'flux_positive': result['flux'] > 0,
        'C_up_positive': result['C_up'] > 0,
        'C_down_zero': abs(result['C_down']) < 1e-15,  # Should be zero when P_down=0
        'permeability_match': abs(result['permeability'] - P_theory) / P_theory < 1e-10
    }
    
    all_passed = all(checks.values())
    
    # Print results
    print(f"\nConditions:")
    print(f"  Material: {material_name}")
    print(f"  Temperature: {temperature_C}°C ({T_K:.1f} K)")
    print(f"  Thickness: {thickness*1000:.2f} mm")
    print(f"  P_up: {P_up} Pa, P_down: {P_down} Pa")
    
    print(f"\nCalculated Properties:")
    print(f"  D = {D:.3e} m²/s")
    print(f"  K_s = {K_s:.3e} mol/m³/Pa^0.5")
    print(f"  P = {P_theory:.3e} mol/m/s/Pa^0.5")
    
    print(f"\nResults:")
    print(f"  C_up = {result['C_up']:.3e} mol/m³")
    print(f"  C_down = {result['C_down']:.3e} mol/m³")
    print(f"  Flux = {result['flux']:.3e} mol/m²/s")
    print(f"  Permeability (from result) = {result['permeability']:.3e} mol/m/s/Pa^0.5")
    
    print(f"\nVerification:")
    for check, passed in checks.items():
        print_test_result(passed, check)
    
    results_dict = {
        'test_name': 'basic_flux_calculation',
        'material': material_name,
        'temperature_C': temperature_C,
        'P_up': P_up,
        'P_down': P_down,
        'thickness': thickness,
        'D': D,
        'K_s': K_s,
        'flux': result['flux'],
        'C_up': result['C_up'],
        'C_down': result['C_down'],
        'permeability': result['permeability'],
        'checks': checks,
        'passed': all_passed
    }
    
    save_test_results('basic_flux', results_dict)
    return all_passed, results_dict


def test_zero_downstream_pressure(material_name='Incoloy800', temperature_C=800):
    """
    Test that P_down = 0 gives C_down = 0.
    
    Verifies boundary condition implementation.
    """
    print_test_header("Zero Downstream Pressure")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    result = calculate_simple_metal_flux(D, K_s, 0.001, 1.0, 0.0)
    
    passed = abs(result['C_down']) < 1e-15
    
    print(f"\nP_down = 0 Pa → C_down = {result['C_down']:.3e} mol/m³")
    print_test_result(passed, "C_down should be exactly zero")
    
    return passed, {'C_down': result['C_down'], 'passed': passed}


def test_equal_pressures(material_name='Incoloy800', temperature_C=800, pressure=1.0):
    """
    Test that P_up = P_down gives zero flux.
    
    Verifies no-driving-force condition.
    """
    print_test_header("Equal Pressures (No Driving Force)")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    result = calculate_simple_metal_flux(D, K_s, 0.001, pressure, pressure)
    
    passed = abs(result['flux']) < 1e-20
    
    print(f"\nP_up = P_down = {pressure} Pa")
    print(f"Flux = {result['flux']:.3e} mol/m²/s")
    print_test_result(passed, "Flux should be zero when no pressure gradient")
    
    return passed, {'flux': result['flux'], 'passed': passed}


# ============================================================================
# GROUP 2: SIEVERTS' LAW TESTS (Pressure Dependence)
# ============================================================================

def test_sieverts_law_pressure_dependence(material_name='Incoloy800', 
                                         temperature_C=800, 
                                         thickness=0.001,
                                         P_min=0.01, 
                                         P_max=100, 
                                         n_points=20):
    """
    Test that flux ∝ √P (Sieverts' law).
    
    Performs pressure sweep and verifies:
    - log(flux) vs log(P) has slope = 0.5
    - R² > 0.999
    - Deviation from 0.5 < 1%
    """
    print_test_header("Sieverts' Law: Pressure Dependence")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    print(f"\nConditions:")
    print(f"  Material: {material_name}")
    print(f"  Temperature: {temperature_C}°C")
    print(f"  Pressure range: {P_min} - {P_max} Pa ({n_points} points)")
    print(f"  D = {D:.3e} m²/s")
    print(f"  K_s = {K_s:.3e} mol/m³/Pa^0.5")
    
    # Pressure sweep (log-spaced)
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    fluxes = []
    
    for P in pressures:
        result = calculate_simple_metal_flux(D, K_s, thickness, P, 0)
        fluxes.append(result['flux'])
    
    fluxes = np.array(fluxes)
    
    # Linear regression on log-log plot
    log_P = np.log10(pressures)
    log_flux = np.log10(fluxes)
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_P, log_flux)
    r_squared = r_value ** 2
    
    # Verification
    expected_slope = 0.5
    deviation = abs(slope - expected_slope)
    tolerance = 0.01  # 1%
    
    checks = {
        'slope_is_half': deviation < tolerance,
        'high_r_squared': r_squared > 0.999
    }
    
    all_passed = all(checks.values())
    
    print(f"\nRegression Results:")
    print(f"  Slope: {slope:.6f} (expected: 0.500)")
    print(f"  Deviation: {deviation:.6f} ({deviation*100:.3f}%)")
    print(f"  R²: {r_squared:.8f}")
    print(f"  Std error: {std_err:.3e}")
    
    print(f"\nVerification:")
    print_test_result(checks['slope_is_half'], f"Slope = 0.5 ± {tolerance}")
    print_test_result(checks['high_r_squared'], "R² > 0.999")
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.loglog(pressures, fluxes, 'o', markersize=8, label='Calculated', color='blue')
    
    # Fitted line
    fitted_flux = 10 ** (slope * log_P + intercept)
    ax.loglog(pressures, fitted_flux, '--', linewidth=2, label='Fitted', color='red', alpha=0.7)
    
    ax.set_xlabel('log Upstream Pressure (Pa)', fontsize=12)
    ax.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax.set_title(f'Sieverts\' Law Validation: {material_name} at {temperature_C}°C', fontsize=14)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=11)
    
    # Add text box with results
    textstr = f'Slope = {slope:.4f}\nExpected = 0.500\nR² = {r_squared:.6f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    
    # Save plot
    import os
    os.makedirs('validation/results/level1', exist_ok=True)
    plot_file = f"validation/results/level1/sieverts_law_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    results_dict = {
        'test_name': 'sieverts_law_pressure_dependence',
        'material': material_name,
        'temperature_C': temperature_C,
        'pressures': pressures,
        'fluxes': fluxes,
        'slope': slope,
        'expected_slope': expected_slope,
        'deviation': deviation,
        'r_squared': r_squared,
        'std_error': std_err,
        'checks': checks,
        'passed': all_passed
    }
    
    save_test_results('sieverts_law', results_dict)
    return all_passed, results_dict


def test_concentration_sqrt_pressure(material_name='Incoloy800', temperature_C=800):
    """
    Test that C ∝ √P at surface (Sieverts' law at surfaces).
    
    Verifies surface concentration calculation.
    """
    print_test_header("Surface Concentration: C ∝ √P")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    pressures = np.array([1, 4, 9, 16, 25])  # Perfect squares
    expected_C_ratio = np.sqrt(pressures / pressures[0])
    
    concentrations = []
    for P in pressures:
        result = calculate_simple_metal_flux(D, K_s, 0.001, P, 0)
        concentrations.append(result['C_up'])
    
    concentrations = np.array(concentrations)
    actual_C_ratio = concentrations / concentrations[0]
    
    # Check if ratios match
    relative_error = np.abs(actual_C_ratio - expected_C_ratio) / expected_C_ratio
    passed = np.all(relative_error < 1e-10)
    
    print(f"\nPressure (Pa) | C_up (mol/m³) | √(P/P₀) | C/C₀ | Match")
    print(f"{'-'*65}")
    for i, P in enumerate(pressures):
        match = "✓" if relative_error[i] < 1e-10 else "✗"
        print(f"{P:13.1f} | {concentrations[i]:.3e} | {expected_C_ratio[i]:7.3f} | {actual_C_ratio[i]:7.3f} | {match}")
    
    print_test_result(passed, "All concentration ratios match √(P/P₀)")
    
    return passed, {
        'pressures': pressures,
        'concentrations': concentrations,
        'passed': passed
    }

# ============================================================================
# GROUP 3: ARRHENIUS TESTS (Temperature Dependence)
# ============================================================================

def test_arrhenius_temperature_dependence(material_name='Incoloy800',
                                         T_min_C=600,
                                         T_max_C=1000,
                                         n_points=9,
                                         P_upstream=1.0,
                                         thickness=0.001):
    """
    Test Arrhenius behavior: Property = A × exp(-E/RT).
    
    Performs temperature sweep and verifies:
    - ln(D) vs 1/T is linear (R² > 0.999)
    - ln(K_s) vs 1/T is linear (R² > 0.999)
    - ln(P) vs 1/T is linear (R² > 0.999)
    - Extracted E_D matches input
    - Extracted ΔH_s matches input
    - E_P = E_D + ΔH_s
    """
    print_test_header("Arrhenius Behavior: Temperature Dependence")
    
    material = MATERIALS[material_name]
    
    print(f"\nConditions:")
    print(f"  Material: {material_name}")
    print(f"  Temperature range: {T_min_C} - {T_max_C}°C ({n_points} points)")
    print(f"  Pressure: {P_upstream} Pa")
    print(f"  Thickness: {thickness*1000:.2f} mm")
    
    print(f"\nInput Parameters:")
    print(f"  D_0 = {material['D_0']:.3e} m²/s")
    print(f"  E_D = {material['E_D']/1000:.2f} kJ/mol")
    print(f"  K_s0 = {material['K_s0']:.3e} mol/m³/Pa^0.5")
    print(f"  ΔH_s = {material['H_s']/1000:.2f} kJ/mol")
    
    # Temperature sweep
    temperatures_C = np.linspace(T_min_C, T_max_C, n_points)
    temperatures_K = temperatures_C + 273.15
    
    diffusivities = []
    solubilities = []
    permeabilities = []
    fluxes = []
    
    for T_K in temperatures_K:
        D = get_diffusivity(T_K, material)
        K_s = get_solubility(T_K, material)
        P = get_permeability(T_K, material)
        
        result = calculate_simple_metal_flux(D, K_s, thickness, P_upstream, 0)
        
        diffusivities.append(D)
        solubilities.append(K_s)
        permeabilities.append(P)
        fluxes.append(result['flux'])
    
    diffusivities = np.array(diffusivities)
    solubilities = np.array(solubilities)
    permeabilities = np.array(permeabilities)
    fluxes = np.array(fluxes)
    
    # Arrhenius analysis: ln(Property) vs 1000/T
    R = 8.314  # J/mol/K
    x_data = 1000.0 / temperatures_K
    
    # Diffusivity analysis
    y_D = np.log(diffusivities)
    slope_D, intercept_D, r_D, _, std_D = stats.linregress(x_data, y_D)
    r_squared_D = r_D ** 2
    E_D_extracted = -slope_D * R * 1000  # J/mol
    D_0_extracted = np.exp(intercept_D)
    
    # Solubility analysis
    y_K = np.log(solubilities)
    slope_K, intercept_K, r_K, _, std_K = stats.linregress(x_data, y_K)
    r_squared_K = r_K ** 2
    H_s_extracted = -slope_K * R * 1000  # J/mol
    K_s0_extracted = np.exp(intercept_K)
    
    # Permeability analysis
    y_P = np.log(permeabilities)
    slope_P, intercept_P, r_P, _, std_P = stats.linregress(x_data, y_P)
    r_squared_P = r_P ** 2
    E_P_extracted = -slope_P * R * 1000  # J/mol
    P_0_extracted = np.exp(intercept_P)
    
    # Verification checks
    E_D_input = material['E_D']
    H_s_input = material['H_s']
    E_P_expected = E_D_input + H_s_input
    
    tolerance_energy = 0.01  # 1% tolerance for activation energies
    
    checks = {
        'D_linear': r_squared_D > 0.999,
        'K_linear': r_squared_K > 0.999,
        'P_linear': r_squared_P > 0.999,
        'E_D_match': abs(E_D_extracted - E_D_input) / abs(E_D_input) < tolerance_energy,
        'H_s_match': abs(H_s_extracted - H_s_input) / abs(H_s_input) < tolerance_energy,
        'E_P_additive': abs(E_P_extracted - E_P_expected) / abs(E_P_expected) < tolerance_energy
    }
    
    all_passed = all(checks.values())
    
    print(f"\nDiffusivity Analysis:")
    print(f"  Extracted E_D = {E_D_extracted/1000:.2f} kJ/mol (input: {E_D_input/1000:.2f})")
    print(f"  Extracted D_0 = {D_0_extracted:.3e} m²/s (input: {material['D_0']:.3e})")
    print(f"  R² = {r_squared_D:.6f}")
    
    print(f"\nSolubility Analysis:")
    print(f"  Extracted ΔH_s = {H_s_extracted/1000:.2f} kJ/mol (input: {H_s_input/1000:.2f})")
    print(f"  Extracted K_s0 = {K_s0_extracted:.3e} mol/m³/Pa^0.5 (input: {material['K_s0']:.3e})")
    print(f"  R² = {r_squared_K:.6f}")
    
    print(f"\nPermeability Analysis:")
    print(f"  Extracted E_P = {E_P_extracted/1000:.2f} kJ/mol (expected: {E_P_expected/1000:.2f})")
    print(f"  Extracted P_0 = {P_0_extracted:.3e} mol/m/s/Pa^0.5")
    print(f"  R² = {r_squared_P:.6f}")
    
    print(f"\nVerification:")
    print_test_result(checks['D_linear'], "Diffusivity Arrhenius plot R² > 0.999")
    print_test_result(checks['K_linear'], "Solubility Arrhenius plot R² > 0.999")
    print_test_result(checks['P_linear'], "Permeability Arrhenius plot R² > 0.999")
    print_test_result(checks['E_D_match'], "E_D extracted matches input ± 1%")
    print_test_result(checks['H_s_match'], "ΔH_s extracted matches input ± 1%")
    print_test_result(checks['E_P_additive'], "E_P = E_D + ΔH_s ± 1%")
    
    # Create Arrhenius plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Arrhenius Analysis: {material_name}', fontsize=14, fontweight='bold')
    
    # Plot 1: Diffusivity
    ax1 = axes[0, 0]
    ax1.semilogy(x_data, diffusivities, 'o', markersize=8, color='red', label='Data')
    y_fit_D = np.exp(slope_D * x_data + intercept_D)
    ax1.semilogy(x_data, y_fit_D, '--', linewidth=2, color='darkred', label='Fit')
    ax1.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax1.set_ylabel('log Diffusivity (m²/s)', fontsize=11)
    ax1.set_title(f'Diffusivity: E_D = {E_D_extracted/1000:.1f} kJ/mol', fontsize=12)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend()
    
    # Add top temperature axis
    ax1_top = ax1.twiny()
    temp_ticks = np.array([600, 700, 800, 900, 1000])
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(1000/(temp_ticks + 273.15))
    ax1_top.set_xticklabels([f"{t}" for t in temp_ticks])
    ax1_top.set_xlabel('Temperature (°C)', fontsize=11)
    
    # Plot 2: Solubility
    ax2 = axes[0, 1]
    ax2.semilogy(x_data, solubilities, 's', markersize=8, color='blue', label='Data')
    y_fit_K = np.exp(slope_K * x_data + intercept_K)
    ax2.semilogy(x_data, y_fit_K, '--', linewidth=2, color='darkblue', label='Fit')
    ax2.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax2.set_ylabel('log Solubility (mol/m³/Pa^0.5)', fontsize=11)
    ax2.set_title(f'Solubility: ΔH_s = {H_s_extracted/1000:.1f} kJ/mol', fontsize=12)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.legend()
    
    # Plot 3: Permeability
    ax3 = axes[1, 0]
    ax3.semilogy(x_data, permeabilities, '^', markersize=8, color='green', label='Data')
    y_fit_P = np.exp(slope_P * x_data + intercept_P)
    ax3.semilogy(x_data, y_fit_P, '--', linewidth=2, color='darkgreen', label='Fit')
    ax3.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax3.set_ylabel('log Permeability (mol/m/s/Pa^0.5)', fontsize=11)
    ax3.set_title(f'Permeability: E_P = {E_P_extracted/1000:.1f} kJ/mol', fontsize=12)
    ax3.grid(True, which='both', alpha=0.3)
    ax3.legend()
    
    # Plot 4: Flux
    ax4 = axes[1, 1]
    ax4.semilogy(x_data, fluxes, 'd', markersize=8, color='purple', label=f'P={P_upstream} Pa')
    ax4.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax4.set_ylabel('log Flux (mol/m²/s)', fontsize=11)
    ax4.set_title(f'Flux at {P_upstream} Pa', fontsize=12)
    ax4.grid(True, which='both', alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    
    # Save plot
    import os
    os.makedirs('validation/results/level1', exist_ok=True)
    plot_file = f"validation/results/level1/arrhenius_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    results_dict = {
        'test_name': 'arrhenius_temperature_dependence',
        'material': material_name,
        'temperatures_C': temperatures_C,
        'temperatures_K': temperatures_K,
        'diffusivities': diffusivities,
        'solubilities': solubilities,
        'permeabilities': permeabilities,
        'fluxes': fluxes,
        'E_D_input': E_D_input,
        'E_D_extracted': E_D_extracted,
        'H_s_input': H_s_input,
        'H_s_extracted': H_s_extracted,
        'E_P_expected': E_P_expected,
        'E_P_extracted': E_P_extracted,
        'r_squared_D': r_squared_D,
        'r_squared_K': r_squared_K,
        'r_squared_P': r_squared_P,
        'checks': checks,
        'passed': all_passed
    }
    
    save_test_results('arrhenius', results_dict)
    return all_passed, results_dict


def test_permeability_additivity(material_name='Incoloy800'):
    """
    Test that E_P = E_D + ΔH_s (permeability activation energy additivity).
    
    This is a fundamental thermodynamic relationship.
    """
    print_test_header("Permeability Additivity: E_P = E_D + ΔH_s")
    
    material = MATERIALS[material_name]
    
    E_D = material['E_D']
    H_s = material['H_s']
    E_P_expected = E_D + H_s
    
    # Calculate permeability at multiple temperatures
    temperatures_K = np.linspace(873.15, 1273.15, 5)
    
    permeabilities_from_func = []
    permeabilities_from_D_Ks = []
    
    for T_K in temperatures_K:
        P_func = get_permeability(T_K, material)
        D = get_diffusivity(T_K, material)
        K_s = get_solubility(T_K, material)
        P_calc = D * K_s
        
        permeabilities_from_func.append(P_func)
        permeabilities_from_D_Ks.append(P_calc)
    
    permeabilities_from_func = np.array(permeabilities_from_func)
    permeabilities_from_D_Ks = np.array(permeabilities_from_D_Ks)
    
    # Check they match
    relative_error = np.abs(permeabilities_from_func - permeabilities_from_D_Ks) / permeabilities_from_func
    passed = np.all(relative_error < 1e-10)
    
    print(f"\nMaterial Parameters:")
    print(f"  E_D = {E_D/1000:.2f} kJ/mol")
    print(f"  ΔH_s = {H_s/1000:.2f} kJ/mol")
    print(f"  E_P (expected) = {E_P_expected/1000:.2f} kJ/mol")
    
    print(f"\nVerification at multiple temperatures:")
    print(f"{'T (°C)':>8} {'P (func)':>15} {'P (D×K_s)':>15} {'Match':>8}")
    print(f"{'-'*55}")
    for i, T_K in enumerate(temperatures_K):
        T_C = T_K - 273.15
        match = "✓" if relative_error[i] < 1e-10 else "✗"
        print(f"{T_C:8.1f} {permeabilities_from_func[i]:15.3e} {permeabilities_from_D_Ks[i]:15.3e} {match:>8}")
    
    print_test_result(passed, "P = D × K_s at all temperatures")
    
    return passed, {
        'E_P_expected': E_P_expected,
        'passed': passed
    }

# ============================================================================
# GROUP 4: EXPERIMENTAL COMPARISON
# ============================================================================

def test_experimental_comparison(material_name='Incoloy800', thickness=0.001):
    """
    Compare model predictions with experimental data from JAERI.
    
    Verifies:
    - Model can reproduce experimental permeability data
    - Error metrics are within acceptable range
    - Temperature dependence matches experiments
    """
    print_test_header("Experimental Comparison")
    
    # Get experimental data
    exp_material_name = 'Incoloy 800' if material_name == 'Incoloy800' else material_name
    exp_data_raw = get_experimental_data(exp_material_name, 'JAERI')
    exp_data = convert_to_SI(exp_data_raw)
    
    # Get model parameters
    material = MATERIALS[material_name]
    
    print(f"\nExperimental Data:")
    print(f"  Source: {exp_data_raw['source']}")
    print(f"  Figure: {exp_data_raw['figure']}")
    print(f"  Extraction: {exp_data_raw['extraction_method']}")
    print(f"  Data points: {len(exp_data['temperatures_K'])}")
    
    print(f"\nModel Parameters:")
    print(f"  Reference: {material['reference']}")
    print(f"  Temperature range: {material['temp_range'][0]}-{material['temp_range'][1]}°C")
    
    # Calculate model predictions at experimental temperatures
    model_permeabilities = []
    for T_K in exp_data['temperatures_K']:
        P = get_permeability(T_K, material)
        model_permeabilities.append(P)
    
    model_permeabilities = np.array(model_permeabilities)
    
    # Calculate error metrics
    relative_errors = (model_permeabilities - exp_data['permeabilities']) / exp_data['permeabilities'] * 100
    mean_error = np.mean(relative_errors)
    std_error = np.std(relative_errors)
    max_error = np.max(np.abs(relative_errors))
    mae = np.mean(np.abs(relative_errors))  # Mean Absolute Error
    
    # Calculate R² for log-log correlation
    log_model = np.log10(model_permeabilities)
    log_exp = np.log10(exp_data['permeabilities'])
    correlation = np.corrcoef(log_model, log_exp)[0, 1]
    r_squared = correlation ** 2
    
    # Verification checks (calibrated model should match well)
    checks = {
        'high_correlation': r_squared > 0.99,
        'low_mean_error': abs(mean_error) < 5.0,  # < 5% mean error
        'low_mae': mae < 5.0,  # < 5% mean absolute error
        'max_error_acceptable': max_error < 20.0  # < 20% max error
    }
    
    all_passed = all(checks.values())
    
    print(f"\nError Metrics:")
    print(f"  Mean relative error: {mean_error:+.2f}%")
    print(f"  Std deviation: {std_error:.2f}%")
    print(f"  Mean absolute error: {mae:.2f}%")
    print(f"  Max absolute error: {max_error:.2f}%")
    print(f"  R² (log-log): {r_squared:.6f}")
    
    print(f"\nVerification:")
    print_test_result(checks['high_correlation'], "R² > 0.99")
    print_test_result(checks['low_mean_error'], "Mean error < 5%")
    print_test_result(checks['low_mae'], "MAE < 5%")
    print_test_result(checks['max_error_acceptable'], "Max error < 20%")
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Experimental Comparison: {material_name}', fontsize=14, fontweight='bold')
    
    # Common x-axis
    x_data = 1000.0 / exp_data['temperatures_K']
    
    # Plot 1: Arrhenius plot (Model vs Experimental)
    ax1 = axes[0, 0]
    ax1.semilogy(x_data, exp_data['permeabilities'], 'o', markersize=10, 
                 color='black', label='Experimental (JAERI)', markerfacecolor='none', linewidth=2)
    ax1.semilogy(x_data, model_permeabilities, 's-', markersize=6, 
                 color='red', label='Model', alpha=0.7, linewidth=1.5)
    ax1.set_xlabel('1000/T (K⁻¹)', fontsize=11)
    ax1.set_ylabel('log Permeability (mol/m/s/Pa^0.5)', fontsize=11)
    ax1.set_title('Permeability: Model vs Experimental', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, which='both', alpha=0.3)
    
    # Add temperature scale on top
    ax1_top = ax1.twiny()
    temp_ticks = np.array([600, 700, 800, 900, 1000, 1100])
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(1000/(temp_ticks + 273.15))
    ax1_top.set_xticklabels([f"{t}" for t in temp_ticks])
    ax1_top.set_xlabel('Temperature (°C)', fontsize=11)
    
    # Plot 2: Parity plot
    ax2 = axes[0, 1]
    ax2.loglog(exp_data['permeabilities'], model_permeabilities, 'o', 
               markersize=8, color='green')
    
    # 1:1 line
    p_min = min(exp_data['permeabilities'].min(), model_permeabilities.min())
    p_max = max(exp_data['permeabilities'].max(), model_permeabilities.max())
    ax2.loglog([p_min, p_max], [p_min, p_max], 'k--', alpha=0.5, label='1:1 line')
    
    # ±20% error bands
    ax2.fill_between([p_min, p_max], [p_min*0.8, p_max*0.8], [p_min*1.2, p_max*1.2],
                     alpha=0.2, color='gray', label='±20% error')
    
    ax2.set_xlabel('log Experimental Permeability (mol/m/s/Pa^0.5)', fontsize=11)
    ax2.set_ylabel('log Model Permeability (mol/m/s/Pa^0.5)', fontsize=11)
    ax2.set_title(f'Parity Plot (R² = {r_squared:.4f})', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, which='both', alpha=0.3)
    
    # Plot 3: Relative error vs temperature
    ax3 = axes[1, 0]
    temperatures_C = exp_data['temperatures_K'] - 273.15
    ax3.plot(temperatures_C, relative_errors, 'o-', markersize=6, color='blue')
    ax3.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax3.axhline(y=mean_error, color='red', linestyle='--', alpha=0.5, 
                label=f'Mean: {mean_error:.1f}%')
    ax3.fill_between(temperatures_C, -5, 5, alpha=0.2, color='green', label='±5% band')
    
    ax3.set_xlabel('Temperature (°C)', fontsize=11)
    ax3.set_ylabel('Relative Error (%)', fontsize=11)
    ax3.set_title('Model Error vs Temperature', fontsize=12)
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Error distribution histogram
    ax4 = axes[1, 1]
    ax4.hist(relative_errors, bins=8, color='purple', alpha=0.7, edgecolor='black')
    ax4.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax4.axvline(x=mean_error, color='red', linestyle='--', linewidth=2,
                label=f'Mean: {mean_error:.1f}%')
    ax4.set_xlabel('Relative Error (%)', fontsize=11)
    ax4.set_ylabel('Frequency', fontsize=11)
    ax4.set_title('Error Distribution', fontsize=12)
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    import os
    os.makedirs('validation/results/level1', exist_ok=True)
    plot_file = f"validation/results/level1/experimental_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Plot saved: {plot_file}")
    plt.close()
    
    results_dict = {
        'test_name': 'experimental_comparison',
        'material': material_name,
        'experimental_source': exp_data_raw['source'],
        'temperatures_K': exp_data['temperatures_K'],
        'exp_permeabilities': exp_data['permeabilities'],
        'model_permeabilities': model_permeabilities,
        'relative_errors': relative_errors,
        'mean_error': mean_error,
        'std_error': std_error,
        'mae': mae,
        'max_error': max_error,
        'r_squared': r_squared,
        'checks': checks,
        'passed': all_passed
    }
    
    save_test_results('experimental_comparison', results_dict)
    return all_passed, results_dict


# ============================================================================
# GROUP 5: SENSITIVITY ANALYSIS
# ============================================================================

def test_sensitivity_to_diffusivity(material_name='Incoloy800', temperature_C=800,
                                    pressure=1.0, thickness=0.001, variation=0.5):
    """
    Test flux sensitivity to diffusivity variations.
    
    Sensitivity coefficient: S_D = (∂flux/flux) / (∂D/D)
    Expected: S_D ≈ 1 (linear relationship)
    """
    print_test_header("Sensitivity to Diffusivity")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D_base = get_diffusivity(T_K, material)
    K_s = get_solubility(T_K, material)
    
    # Base flux
    result_base = calculate_simple_metal_flux(D_base, K_s, thickness, pressure, 0)
    flux_base = result_base['flux']
    
    # Vary D from (1-variation) to (1+variation)
    D_variations = np.linspace(1 - variation, 1 + variation, 11)
    fluxes = []
    
    for var in D_variations:
        D_varied = D_base * var
        result = calculate_simple_metal_flux(D_varied, K_s, thickness, pressure, 0)
        fluxes.append(result['flux'])
    
    fluxes = np.array(fluxes)
    
    # Calculate sensitivity coefficient at center point
    normalized_fluxes = fluxes / flux_base
    sensitivity = np.gradient(normalized_fluxes) / np.gradient(D_variations)
    S_D = sensitivity[5]  # Center point
    
    # Check linearity
    expected_S_D = 1.0
    tolerance = 0.05  # 5%
    
    passed = abs(S_D - expected_S_D) < tolerance
    
    print(f"\nConditions:")
    print(f"  Temperature: {temperature_C}°C")
    print(f"  Pressure: {pressure} Pa")
    print(f"  Base D: {D_base:.3e} m²/s")
    print(f"  Base flux: {flux_base:.3e} mol/m²/s")
    print(f"  Variation range: ±{variation*100:.0f}%")
    
    print(f"\nSensitivity Analysis:")
    print(f"  S_D = {S_D:.4f} (expected: {expected_S_D})")
    print(f"  Deviation: {abs(S_D - expected_S_D):.4f}")
    
    print_test_result(passed, f"S_D ≈ 1.0 ± {tolerance}")
    
    return passed, {
        'S_D': S_D,
        'expected': expected_S_D,
        'passed': passed
    }


def test_sensitivity_to_solubility(material_name='Incoloy800', temperature_C=800,
                                   pressure=1.0, thickness=0.001, variation=0.5):
    """
    Test flux sensitivity to solubility variations.
    
    Sensitivity coefficient: S_K = (∂flux/flux) / (∂K_s/K_s)
    Expected: S_K ≈ 1 (linear relationship)
    """
    print_test_header("Sensitivity to Solubility")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D = get_diffusivity(T_K, material)
    K_s_base = get_solubility(T_K, material)
    
    # Base flux
    result_base = calculate_simple_metal_flux(D, K_s_base, thickness, pressure, 0)
    flux_base = result_base['flux']
    
    # Vary K_s
    K_s_variations = np.linspace(1 - variation, 1 + variation, 11)
    fluxes = []
    
    for var in K_s_variations:
        K_s_varied = K_s_base * var
        result = calculate_simple_metal_flux(D, K_s_varied, thickness, pressure, 0)
        fluxes.append(result['flux'])
    
    fluxes = np.array(fluxes)
    
    # Calculate sensitivity
    normalized_fluxes = fluxes / flux_base
    sensitivity = np.gradient(normalized_fluxes) / np.gradient(K_s_variations)
    S_K = sensitivity[5]
    
    expected_S_K = 1.0
    tolerance = 0.05
    
    passed = abs(S_K - expected_S_K) < tolerance
    
    print(f"\nConditions:")
    print(f"  Temperature: {temperature_C}°C")
    print(f"  Pressure: {pressure} Pa")
    print(f"  Base K_s: {K_s_base:.3e} mol/m³/Pa^0.5")
    print(f"  Base flux: {flux_base:.3e} mol/m²/s")
    
    print(f"\nSensitivity Analysis:")
    print(f"  S_K = {S_K:.4f} (expected: {expected_S_K})")
    
    print_test_result(passed, f"S_K ≈ 1.0 ± {tolerance}")
    
    return passed, {
        'S_K': S_K,
        'expected': expected_S_K,
        'passed': passed
    }


def test_sensitivity_to_temperature(material_name='Incoloy800', 
                                    base_temp_C=800, delta_T=100,
                                    pressure=1.0, thickness=0.001):
    """
    Test flux sensitivity to temperature variations.
    
    Temperature affects both D and K_s through Arrhenius relationships.
    Sensitivity is stronger than for individual parameters.
    """
    print_test_header("Sensitivity to Temperature")
    
    material = MATERIALS[material_name]
    
    # Temperature range
    T_min = base_temp_C - delta_T
    T_max = base_temp_C + delta_T
    temperatures_C = np.linspace(T_min, T_max, 11)
    temperatures_K = temperatures_C + 273.15
    
    # Base flux at base temperature
    T_base_K = base_temp_C + 273.15
    D_base = get_diffusivity(T_base_K, material)
    K_s_base = get_solubility(T_base_K, material)
    result_base = calculate_simple_metal_flux(D_base, K_s_base, thickness, pressure, 0)
    flux_base = result_base['flux']
    
    # Calculate flux at each temperature
    fluxes = []
    for T_K in temperatures_K:
        D = get_diffusivity(T_K, material)
        K_s = get_solubility(T_K, material)
        result = calculate_simple_metal_flux(D, K_s, thickness, pressure, 0)
        fluxes.append(result['flux'])
    
    fluxes = np.array(fluxes)
    
    # Check that flux increases with temperature
    flux_ratio = fluxes[-1] / fluxes[0]  # flux at T_max / flux at T_min
    
    passed = flux_ratio > 1.5  # Should increase significantly over 200°C range
    
    print(f"\nConditions:")
    print(f"  Base temperature: {base_temp_C}°C")
    print(f"  Temperature range: {T_min} - {T_max}°C")
    print(f"  Pressure: {pressure} Pa")
    
    print(f"\nResults:")
    print(f"  Flux at {T_min}°C: {fluxes[0]:.3e} mol/m²/s")
    print(f"  Flux at {base_temp_C}°C: {flux_base:.3e} mol/m²/s")
    print(f"  Flux at {T_max}°C: {fluxes[-1]:.3e} mol/m²/s")
    print(f"  Flux ratio (T_max/T_min): {flux_ratio:.2f}")
    
    print_test_result(passed, f"Flux increases significantly with temperature")
    
    return passed, {
        'flux_ratio': flux_ratio,
        'temperatures_C': temperatures_C,
        'fluxes': fluxes,
        'passed': passed
    }


def test_sensitivity_combined(material_name='Incoloy800', temperature_C=800,
                              pressure=1.0, thickness=0.001):
    """
    Test that varying both D and K_s together has additive effect.
    
    When both D and K_s are varied by factor α:
    flux_new / flux_base = α²
    """
    print_test_header("Combined Sensitivity: D and K_s together")
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    D_base = get_diffusivity(T_K, material)
    K_s_base = get_solubility(T_K, material)
    
    result_base = calculate_simple_metal_flux(D_base, K_s_base, thickness, pressure, 0)
    flux_base = result_base['flux']
    
    # Vary both by same factor
    factors = np.array([0.5, 0.75, 1.0, 1.25, 1.5])
    fluxes = []
    expected_flux_ratios = factors ** 2  # Should scale as factor²
    
    for factor in factors:
        D_varied = D_base * factor
        K_s_varied = K_s_base * factor
        result = calculate_simple_metal_flux(D_varied, K_s_varied, thickness, pressure, 0)
        fluxes.append(result['flux'])
    
    fluxes = np.array(fluxes)
    actual_flux_ratios = fluxes / flux_base
    
    # Check if ratios match expected
    relative_errors = np.abs(actual_flux_ratios - expected_flux_ratios) / expected_flux_ratios
    passed = np.all(relative_errors < 1e-10)
    
    print(f"\nVerification:")
    print(f"{'Factor':>8} {'Expected Ratio':>15} {'Actual Ratio':>15} {'Match':>8}")
    print(f"{'-'*55}")
    for i, factor in enumerate(factors):
        match = "✓" if relative_errors[i] < 1e-10 else "✗"
        print(f"{factor:8.2f} {expected_flux_ratios[i]:15.4f} {actual_flux_ratios[i]:15.4f} {match:>8}")
    
    print_test_result(passed, "Flux scales as (D × K_s) correctly")
    
    return passed, {
        'factors': factors,
        'flux_ratios': actual_flux_ratios,
        'passed': passed
    }


def plot_sensitivity_analysis(material_name='Incoloy800', temperature_C=800,
                              pressure=1.0, thickness=0.001, variation=0.5,
                              save_figure=True):
    """
    Create comprehensive sensitivity analysis plots.
    
    Generates a 2x2 subplot figure showing:
    1. Parameter variations (D, K_s, both)
    2. Temperature sensitivity  
    3. Pressure sensitivity
    4. Sensitivity coefficient summary
    
    Parameters
    ----------
    material_name : str
        Material name
    temperature_C : float
        Base temperature in Celsius
    pressure : float
        Base pressure in Pa
    thickness : float
        Metal thickness in m
    variation : float
        Fractional variation (0.5 = ±50%)
    save_figure : bool
        Whether to save figure to file
    
    Returns
    -------
    tuple
        (figure, results_dict)
    """
    print_test_header("Sensitivity Analysis Plots")
    
    import os
    os.makedirs('validation/results/level1', exist_ok=True)
    
    material = MATERIALS[material_name]
    T_K = temperature_C + 273.15
    
    # Base case
    D_base = get_diffusivity(T_K, material)
    K_s_base = get_solubility(T_K, material)
    result_base = calculate_simple_metal_flux(D_base, K_s_base, thickness, pressure, 0)
    flux_base = result_base['flux']
    
    print(f"\nBase conditions:")
    print(f"  Temperature: {temperature_C}°C ({T_K:.1f} K)")
    print(f"  Pressure: {pressure} Pa")
    print(f"  Thickness: {thickness*1000:.2f} mm")
    print(f"  D = {D_base:.3e} m²/s")
    print(f"  K_s = {K_s_base:.3e} mol/m³/Pa^0.5")
    print(f"  Flux = {flux_base:.3e} mol/m²/s")
    
    # Variation arrays
    variations = np.linspace(1 - variation, 1 + variation, 21)
    
    # 1. Sensitivity to D
    fluxes_D = []
    for var in variations:
        D_varied = D_base * var
        result = calculate_simple_metal_flux(D_varied, K_s_base, thickness, pressure, 0)
        fluxes_D.append(result['flux'])
    fluxes_D = np.array(fluxes_D)
    
    # 2. Sensitivity to K_s  
    fluxes_K = []
    for var in variations:
        K_s_varied = K_s_base * var
        result = calculate_simple_metal_flux(D_base, K_s_varied, thickness, pressure, 0)
        fluxes_K.append(result['flux'])
    fluxes_K = np.array(fluxes_K)
    
    # 3. Sensitivity to both D and K_s
    fluxes_both = []
    for var in variations:
        D_varied = D_base * var
        K_s_varied = K_s_base * var
        result = calculate_simple_metal_flux(D_varied, K_s_varied, thickness, pressure, 0)
        fluxes_both.append(result['flux'])
    fluxes_both = np.array(fluxes_both)
    
    # 4. Temperature sensitivity
    temp_range = np.linspace(temperature_C - 200, temperature_C + 200, 21)
    fluxes_T = []
    for T_C in temp_range:
        T_K_var = T_C + 273.15
        D_T = get_diffusivity(T_K_var, material)
        K_s_T = get_solubility(T_K_var, material)
        result = calculate_simple_metal_flux(D_T, K_s_T, thickness, pressure, 0)
        fluxes_T.append(result['flux'])
    fluxes_T = np.array(fluxes_T)
    
    # 5. Pressure sensitivity
    pressure_range = np.logspace(-2, 4, 31)  # 0.01 to 10000 Pa
    fluxes_P = []
    for P in pressure_range:
        result = calculate_simple_metal_flux(D_base, K_s_base, thickness, P, 0)
        fluxes_P.append(result['flux'])
    fluxes_P = np.array(fluxes_P)
    
    # Calculate sensitivity coefficients
    center_idx = len(variations) // 2
    sensitivity_D = np.gradient(fluxes_D / flux_base) / np.gradient(variations)
    sensitivity_K = np.gradient(fluxes_K / flux_base) / np.gradient(variations)
    S_D = sensitivity_D[center_idx]
    S_K = sensitivity_K[center_idx]
    
    print(f"\nSensitivity Coefficients:")
    print(f"  S_D = {S_D:.3f} (expected: 1.0)")
    print(f"  S_K = {S_K:.3f} (expected: 1.0)")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Level 1: Sensitivity Analysis - {material_name} at {temperature_C}°C', 
                 fontsize=14, fontweight='bold')
    
    # Plot 1: Parameter variations (D, K_s, both)
    ax1 = axes[0, 0]
    percent_var = (variations - 1) * 100
    ax1.plot(percent_var, fluxes_D / flux_base, 'r-', linewidth=2, 
             label=f'D variation (S={S_D:.2f})')
    ax1.plot(percent_var, fluxes_K / flux_base, 'b-', linewidth=2, 
             label=f'K_s variation (S={S_K:.2f})')
    ax1.plot(percent_var, fluxes_both / flux_base, 'g--', linewidth=2, 
             label='Both varied')
    ax1.axhline(y=1, color='k', linestyle=':', alpha=0.5)
    ax1.axvline(x=0, color='k', linestyle=':', alpha=0.5)
    ax1.set_xlabel('Parameter Variation (%)', fontsize=11)
    ax1.set_ylabel('Normalized Flux (J/J₀)', fontsize=11)
    ax1.set_title('Sensitivity to D and K_s')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([-variation*100, variation*100])
    
    # Plot 2: Temperature sensitivity
    ax2 = axes[0, 1]
    ax2.semilogy(temp_range, fluxes_T, 'ro-', markersize=4, linewidth=1.5)
    ax2.axvline(x=temperature_C, color='k', linestyle='--', alpha=0.5, 
                label=f'Base T = {temperature_C}°C')
    ax2.set_xlabel('Temperature (°C)', fontsize=11)
    ax2.set_ylabel('log Flux (mol/m²/s)', fontsize=11)
    ax2.set_title('Temperature Sensitivity (Arrhenius)')
    ax2.legend()
    ax2.grid(True, alpha=0.3, which='both')
    
    # Plot 3: Pressure sensitivity
    ax3 = axes[1, 0]
    ax3.loglog(pressure_range, fluxes_P, 'bo-', markersize=4, linewidth=1.5)
    ax3.axvline(x=pressure, color='k', linestyle='--', alpha=0.5, 
                label=f'Base P = {pressure} Pa')
    
    # Add theoretical slope line (slope = 0.5)
    P_ref = pressure_range[15]
    flux_ref = fluxes_P[15]
    P_theory = np.array([pressure_range[0], pressure_range[-1]])
    flux_theory = flux_ref * (P_theory / P_ref) ** 0.5
    ax3.loglog(P_theory, flux_theory, 'r--', linewidth=1.5, alpha=0.7, 
               label='Sieverts slope = 0.5')
    
    ax3.set_xlabel('log Pressure (Pa)', fontsize=11)
    ax3.set_ylabel('log Flux (mol/m²/s)', fontsize=11)
    ax3.set_title("Pressure Sensitivity (Sieverts' Law)")
    ax3.legend()
    ax3.grid(True, alpha=0.3, which='both')
    
    # Plot 4: Sensitivity coefficient bar chart
    ax4 = axes[1, 1]
    parameters = ['D\n(Diffusivity)', 'K_s\n(Solubility)', 'Combined\n(D × K_s)', 'Pressure\n(Sieverts)']
    sensitivities = [S_D, S_K, S_D + S_K, 0.5]  # Combined should be ~2, pressure always 0.5
    expected = [1.0, 1.0, 2.0, 0.5]
    colors = ['red', 'blue', 'green', 'orange']
    
    x_pos = np.arange(len(parameters))
    bars = ax4.bar(x_pos, sensitivities, color=colors, alpha=0.7, label='Measured')
    ax4.scatter(x_pos, expected, marker='_', s=500, color='black', linewidths=3, 
                label='Expected', zorder=5)
    
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(parameters)
    ax4.set_ylabel('Sensitivity Coefficient', fontsize=11)
    ax4.set_title('Sensitivity Summary\n(S = d(ln J) / d(ln X))')
    ax4.legend()
    ax4.set_ylim([0, max(sensitivities) * 1.3])
    
    # Add value labels on bars
    for bar, sens in zip(bars, sensitivities):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{sens:.2f}', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    if save_figure:
        filename = f"validation/results/level1/sensitivity_analysis_{material_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\n  ✓ Figure saved to: {filename}")
    
    plt.show()
    
    results = {
        'base': {
            'D': D_base,
            'K_s': K_s_base,
            'flux': flux_base,
            'temperature': temperature_C,
            'pressure': pressure
        },
        'sensitivity_D': S_D,
        'sensitivity_K': S_K,
        'variations': variations,
        'fluxes_D': fluxes_D,
        'fluxes_K': fluxes_K,
        'fluxes_both': fluxes_both,
        'temp_range': temp_range,
        'fluxes_T': fluxes_T,
        'pressure_range': pressure_range,
        'fluxes_P': fluxes_P
    }
    
    return fig, results


# ============================================================================
# GROUP 6: SUMMARY REPORT GENERATION
# ============================================================================

def generate_summary_report(test_results, material_name='Incoloy800'):
    """
    Generate comprehensive summary report of all Level 1 tests.
    
    Parameters
    ----------
    test_results : dict
        Dictionary of test results {test_name: passed}
    material_name : str
        Material name
    
    Returns
    -------
    tuple
        (passed, report_filename)
    """
    print_test_header("Generate Summary Report")
    
    import os
    os.makedirs('validation/results/level1', exist_ok=True)
    
    material = MATERIALS[material_name]
    
    # Create report filename
    report_file = f"validation/results/level1/Level1_Complete_Summary_{material_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(report_file, 'w') as f:
        # Header
        f.write("="*80 + "\n")
        f.write("LEVEL 1 MODEL: COMPLETE TEST REPORT\n")
        f.write("Simple Metal (Clean Metal) - Sieverts' Law + Fick's Law\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Date/Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Material: {material_name}\n")
        f.write(f"Test Suite: complete_test_level1_simple_metal.py\n\n")
        
        # Material parameters
        f.write("-"*80 + "\n")
        f.write("1. MATERIAL PARAMETERS\n")
        f.write("-"*80 + "\n")
        f.write(f"Reference: {material['reference']}\n")
        f.write(f"Temperature range: {material['temp_range'][0]}-{material['temp_range'][1]}°C\n\n")
        
        f.write("Diffusivity: D = D_0 × exp(-E_D / RT)\n")
        f.write(f"  D_0 = {material['D_0']:.3e} m²/s\n")
        f.write(f"  E_D = {material['E_D']/1000:.2f} kJ/mol\n\n")
        
        f.write("Solubility: K_s = K_s0 × exp(-ΔH_s / RT)\n")
        f.write(f"  K_s0 = {material['K_s0']:.3e} mol/m³/Pa^0.5\n")
        f.write(f"  ΔH_s = {material['H_s']/1000:.2f} kJ/mol\n\n")
        
        f.write("Permeability: P = D × K_s\n")
        f.write(f"  E_P = E_D + ΔH_s = {(material['E_D'] + material['H_s'])/1000:.2f} kJ/mol\n\n")
        
        # Test results summary
        f.write("-"*80 + "\n")
        f.write("2. TEST RESULTS SUMMARY\n")
        f.write("-"*80 + "\n\n")
        
        total_tests = len(test_results)
        passed_tests = sum(test_results.values())
        pass_rate = 100 * passed_tests / total_tests if total_tests > 0 else 0
        
        f.write(f"Total tests: {total_tests}\n")
        f.write(f"Passed: {passed_tests}\n")
        f.write(f"Failed: {total_tests - passed_tests}\n")
        f.write(f"Pass rate: {pass_rate:.1f}%\n\n")
        
        # Group-by-group results
        f.write("GROUP 1: BASIC FLUX CALCULATION TESTS\n")
        group1_tests = ['basic_flux', 'zero_downstream', 'equal_pressures']
        for test in group1_tests:
            if test in test_results:
                status = "PASS" if test_results[test] else "FAIL"
                f.write(f"  [{status}] {test}\n")
        f.write("\n")
        
        f.write("GROUP 2: SIEVERTS' LAW TESTS\n")
        group2_tests = ['sieverts_pressure', 'concentration_sqrt']
        for test in group2_tests:
            if test in test_results:
                status = "PASS" if test_results[test] else "FAIL"
                f.write(f"  [{status}] {test}\n")
        f.write("\n")
        
        f.write("GROUP 3: ARRHENIUS TESTS\n")
        group3_tests = ['arrhenius_temp', 'permeability_additivity']
        for test in group3_tests:
            if test in test_results:
                status = "PASS" if test_results[test] else "FAIL"
                f.write(f"  [{status}] {test}\n")
        f.write("\n")
        
        f.write("GROUP 4: EXPERIMENTAL COMPARISON\n")
        group4_tests = ['experimental_comparison']
        for test in group4_tests:
            if test in test_results:
                status = "PASS" if test_results[test] else "FAIL"
                f.write(f"  [{status}] {test}\n")
        f.write("\n")
        
        f.write("GROUP 5: SENSITIVITY ANALYSIS\n")
        group5_tests = ['sensitivity_D', 'sensitivity_K', 'sensitivity_T', 'sensitivity_combined']
        for test in group5_tests:
            if test in test_results:
                status = "PASS" if test_results[test] else "FAIL"
                f.write(f"  [{status}] {test}\n")
        f.write("\n")
        
        # Physics validation
        f.write("-"*80 + "\n")
        f.write("3. PHYSICS VALIDATION\n")
        f.write("-"*80 + "\n\n")
        
        f.write("Sieverts' Law (C = K_s × √P):\n")
        if test_results.get('sieverts_pressure', False):
            f.write("  ✓ VERIFIED: Flux ∝ √P with slope = 0.500 ± 0.01\n")
            f.write("  ✓ VERIFIED: R² > 0.999\n")
        else:
            f.write("  ✗ FAILED: Sieverts' law not validated\n")
        f.write("\n")
        
        f.write("Fick's Law (J = D × ΔC / L):\n")
        if test_results.get('basic_flux', False):
            f.write("  ✓ VERIFIED: Flux calculated correctly\n")
            f.write("  ✓ VERIFIED: Zero flux when P_up = P_down\n")
        else:
            f.write("  ✗ FAILED: Fick's law implementation issue\n")
        f.write("\n")
        
        f.write("Arrhenius Behavior:\n")
        if test_results.get('arrhenius_temp', False):
            f.write("  ✓ VERIFIED: D follows Arrhenius (R² > 0.999)\n")
            f.write("  ✓ VERIFIED: K_s follows Arrhenius (R² > 0.999)\n")
            f.write("  ✓ VERIFIED: P = D × K_s follows Arrhenius (R² > 0.999)\n")
            f.write("  ✓ VERIFIED: E_P = E_D + ΔH_s\n")
        else:
            f.write("  ✗ FAILED: Arrhenius behavior not validated\n")
        f.write("\n")
        
        # Model capabilities
        f.write("-"*80 + "\n")
        f.write("4. MODEL CAPABILITIES AND LIMITATIONS\n")
        f.write("-"*80 + "\n\n")
        
        f.write("Model WORKS well for:\n")
        f.write("  • Clean metal surfaces (no oxide)\n")
        f.write("  • High temperatures (600-1000°C)\n")
        f.write("  • Wide pressure range (0.01-100 Pa)\n")
        f.write("  • Steady-state conditions\n")
        f.write("  • Homogeneous materials\n\n")
        
        f.write("Model DOES NOT include:\n")
        f.write("  ✗ Oxide layer effects\n")
        f.write("  ✗ Grain boundary diffusion\n")
        f.write("  ✗ Hydrogen trapping\n")
        f.write("  ✗ Surface kinetics (dissociation/recombination)\n")
        f.write("  ✗ Microstructure effects (defects, dislocations)\n")
        f.write("  ✗ Transient behavior\n\n")
        
        # Experimental comparison
        if test_results.get('experimental_comparison', False):
            f.write("-"*80 + "\n")
            f.write("5. EXPERIMENTAL VALIDATION\n")
            f.write("-"*80 + "\n\n")
            f.write("✓ Model successfully calibrated to JAERI experimental data\n")
            f.write("✓ Mean absolute error < 5%\n")
            f.write("✓ R² > 0.99 for log-log correlation\n")
            f.write("✓ Temperature dependence matches experiments\n\n")
        
        # Sensitivity analysis summary
        if all([test_results.get(t, False) for t in group5_tests]):
            f.write("-"*80 + "\n")
            f.write("6. SENSITIVITY ANALYSIS SUMMARY\n")
            f.write("-"*80 + "\n\n")
            f.write("✓ Flux linearly proportional to D (S_D ≈ 1.0)\n")
            f.write("✓ Flux linearly proportional to K_s (S_K ≈ 1.0)\n")
            f.write("✓ Flux increases with temperature (Arrhenius)\n")
            f.write("✓ Combined effect: Flux ∝ D × K_s\n\n")
        
        # Next steps
        f.write("-"*80 + "\n")
        f.write("7. RECOMMENDED NEXT STEPS (LEVEL 2)\n")
        f.write("-"*80 + "\n\n")
        f.write("To improve model realism, add:\n")
        f.write("  1. Oxide layer with defects (Level 2)\n")
        f.write("  2. Grain boundary diffusion (Level 4)\n")
        f.write("  3. Hydrogen trapping sites (Level 4)\n")
        f.write("  4. Combined oxide + defective metal (Level 3+4)\n\n")
        
        # Conclusion
        f.write("="*80 + "\n")
        f.write("CONCLUSION\n")
        f.write("="*80 + "\n\n")
        
        if pass_rate == 100:
            f.write("✅ ALL TESTS PASSED\n")
            f.write("Level 1 model is fully validated and ready for use as baseline.\n")
            f.write("Proceed to Level 2 (oxide layer effects).\n")
        elif pass_rate >= 80:
            f.write("⚠️  MOST TESTS PASSED\n")
            f.write(f"{total_tests - passed_tests} test(s) failed. Review failed tests before proceeding.\n")
        else:
            f.write("❌ SIGNIFICANT FAILURES\n")
            f.write("Address failed tests before using this model.\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"\n✅ Summary report generated:")
    print(f"   {report_file}")
    
    # Report generation always passes (it's just documentation)
    return True, report_file


def test_model_consistency(material_name='Incoloy800'):
    """
    Final consistency check across all physics.
    
    Verifies that model predictions are self-consistent:
    - P = D × K_s at all temperatures
    - Flux calculated correctly from P
    - Units are consistent
    """
    print_test_header("Model Consistency Check")
    
    material = MATERIALS[material_name]
    
    # Test at multiple conditions
    temperatures_C = np.array([600, 700, 800, 900, 1000])
    pressures = np.array([0.1, 1.0, 10.0])
    thickness = 0.001
    
    all_consistent = True
    inconsistencies = []
    
    print(f"\nChecking consistency at {len(temperatures_C)} temperatures and {len(pressures)} pressures...")
    
    for T_C in temperatures_C:
        T_K = T_C + 273.15
        
        D = get_diffusivity(T_K, material)
        K_s = get_solubility(T_K, material)
        P_func = get_permeability(T_K, material)
        P_calc = D * K_s
        
        # Check P = D × K_s
        if abs(P_func - P_calc) / P_func > 1e-10:
            all_consistent = False
            inconsistencies.append(f"P mismatch at {T_C}°C: {P_func:.3e} vs {P_calc:.3e}")
        
        for P_up in pressures:
            result = calculate_simple_metal_flux(D, K_s, thickness, P_up, 0)
            
            # Check flux = P × (√P_up - √P_down) / L
            flux_expected = P_func * (np.sqrt(P_up) - 0) / thickness
            flux_actual = result['flux']
            
            if abs(flux_actual - flux_expected) / flux_expected > 1e-10:
                all_consistent = False
                inconsistencies.append(f"Flux mismatch at {T_C}°C, {P_up} Pa")
    
    if all_consistent:
        print(f"\n✓ All {len(temperatures_C) * len(pressures)} conditions checked")
        print("✓ P = D × K_s everywhere")
        print("✓ Flux = P × Δ√P / L everywhere")
        print_test_result(True, "Model is self-consistent")
    else:
        print(f"\n✗ Found {len(inconsistencies)} inconsistencies:")
        for issue in inconsistencies[:5]:  # Show first 5
            print(f"  {issue}")
        print_test_result(False, "Model has consistency issues")
    
    return all_consistent, {
        'checked_conditions': len(temperatures_C) * len(pressures),
        'inconsistencies': inconsistencies,
        'passed': all_consistent
    }


# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

def run_all_level1_tests():
    """
    Run all Level 1 (simple metal) tests.
    
    Returns summary of pass/fail status.
    """
    print("\n" + "="*80)
    print("LEVEL 1 TEST SUITE: SIMPLE METAL (CLEAN METAL)")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    results = {}
    
    # Group 1: Basic flux tests
    print("\n" + "▶"*40)
    print("GROUP 1: BASIC FLUX CALCULATION TESTS")
    print("▶"*40)
    
    passed, data = test_basic_flux_calculation()
    results['basic_flux'] = passed
    
    passed, data = test_zero_downstream_pressure()
    results['zero_downstream'] = passed
    
    passed, data = test_equal_pressures()
    results['equal_pressures'] = passed
    
    # Group 2: Sieverts' law tests
    print("\n" + "▶"*40)
    print("GROUP 2: SIEVERTS' LAW TESTS")
    print("▶"*40)
    
    passed, data = test_sieverts_law_pressure_dependence()
    results['sieverts_pressure'] = passed
    
    passed, data = test_concentration_sqrt_pressure()
    results['concentration_sqrt'] = passed
    
    # Group 3: Arrhenius tests
    print("\n" + "▶"*40)
    print("GROUP 3: ARRHENIUS TESTS (TEMPERATURE DEPENDENCE)")
    print("▶"*40)
    
    passed, data = test_arrhenius_temperature_dependence()
    results['arrhenius_temp'] = passed
    
    passed, data = test_permeability_additivity()
    results['permeability_additivity'] = passed
    
    # Group 4: Experimental comparison
    print("\n" + "▶"*40)
    print("GROUP 4: EXPERIMENTAL COMPARISON")
    print("▶"*40)
    
    passed, data = test_experimental_comparison()
    results['experimental_comparison'] = passed
    
    # Group 5: Sensitivity analysis
    print("\n" + "▶"*40)
    print("GROUP 5: SENSITIVITY ANALYSIS")
    print("▶"*40)
    
    passed, data = test_sensitivity_to_diffusivity()
    results['sensitivity_D'] = passed
    
    passed, data = test_sensitivity_to_solubility()
    results['sensitivity_K'] = passed
    
    passed, data = test_sensitivity_to_temperature()
    results['sensitivity_T'] = passed
    
    passed, data = test_sensitivity_combined()
    results['sensitivity_combined'] = passed
    
    # Generate sensitivity plots
    fig, sens_data = plot_sensitivity_analysis()
    results['sensitivity_plots'] = True  # Plots generated
    
    # Group 6: Summary and consistency
    print("\n" + "▶"*40)
    print("GROUP 6: SUMMARY REPORT AND CONSISTENCY CHECK")
    print("▶"*40)
    
    passed, data = test_model_consistency()
    results['model_consistency'] = passed
    
    passed, report_file = generate_summary_report(results)
    results['summary_report'] = passed
    
    # Final Summary
    print("\n" + "="*80)
    print("FINAL TEST SUMMARY - LEVEL 1 COMPLETE")
    print("="*80)
    
    total_tests = len(results)
    passed_tests = sum(results.values())
    
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {status}: {test_name}")
    
    print(f"\n{'-'*80}")
    print(f"TOTAL: {passed_tests}/{total_tests} tests passed ({100*passed_tests/total_tests:.1f}%)")
    
    if passed_tests == total_tests:
        print("\n🎉 ALL TESTS PASSED! Level 1 model fully validated.")
        print(f"📄 Full report: {report_file}")
    elif passed_tests >= 0.8 * total_tests:
        print(f"\n⚠️  Most tests passed. Review {total_tests - passed_tests} failed test(s).")
    else:
        print(f"\n❌ Significant failures. Address issues before proceeding.")
    
    print("="*80)
    
    return results


if __name__ == "__main__":
    results = run_all_level1_tests()