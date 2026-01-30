"""
Complete Test Suite for Level 4: Defective Metal (Microstructure Effects)
=========================================================================
Covers:
- Basic flux (should match Level 1 if no traps/GBs)
- Trapping effects (trap_density, trap_energy)
- Grain boundary effects (grain_boundary_fraction, gb_enhancement)
- Mixed microstructure
- Pressure, temperature, thickness dependence
- Regime classification (lattice_limited, traps_defect_limited)
"""

import numpy as np
import os
import json
from datetime import datetime
from calculations.permeation_calc import calculate_defective_metal_flux, calculate_simple_metal_flux
import matplotlib.pyplot as plt
# Import microstructure parameters
from data.microstruture_parameters import PROCESSING_CONDITIONS, TRAP_PROPERTIES

# Results directory
RESULTS_DIR = 'validation/results/Azeezlevel4'
os.makedirs(RESULTS_DIR, exist_ok=True)

# Physical constants
R_GAS = 8.314  # J/mol/K

# ============================================================================
# DEFAULT PARAMETERS - Incoloy 800 at 600-1000°C
# ============================================================================

# Temperature range
DEFAULT_T_MIN = 873.15   # K (600°C)
DEFAULT_T_MAX = 1273.15  # K (1000°C)

# Arrhenius parameters for lattice diffusivity: D(T) = D_0 * exp(-E_D / RT)
DEFAULT_D0 = 1.5e-6      # Pre-exponential factor (m²/s) - Incoloy 800
DEFAULT_E_D = 55000      # Activation energy for diffusion (J/mol) ~0.57 eV

# Sieverts' constant
DEFAULT_K_S = 2e-4       # mol/m³/Pa^0.5

# Geometry
DEFAULT_THICKNESS = 1e-3  # 1 mm
DEFAULT_P_UP = 100        # 100 Pa (more realistic H2 pressure)
DEFAULT_P_DOWN = 0.0      # Vacuum downstream

# Pressure range for sweeps
DEFAULT_P_MIN = 1         # Pa
DEFAULT_P_MAX = 1000      # Pa

# Thickness range for sweeps
DEFAULT_L_MIN = 0.1e-3    # 0.1 mm
DEFAULT_L_MAX = 10e-3     # 10 mm

# Microstructure parameters
DEFAULT_GRAIN_SIZE = 50e-6    # 50 μm (solution annealed)
DEFAULT_TRAP_DENSITY = 1e24   # m⁻³ (moderate cold work)
DEFAULT_E_BINDING = 0.5 * 96485  # 0.5 eV binding energy in J/mol

def get_D_lattice(T, D_0=DEFAULT_D0, E_D=DEFAULT_E_D):
    """Calculate temperature-dependent lattice diffusivity using Arrhenius equation.
    
    D(T) = D_0 * exp(-E_D / RT)
    
    Args:
        T: Temperature (K)
        D_0: Pre-exponential factor (m²/s)
        E_D: Activation energy for diffusion (J/mol)
    
    Returns:
        D_lattice at temperature T (m²/s)
    """
    return D_0 * np.exp(-E_D / (R_GAS * T))

# Utility to build microstructure_params from PROCESSING_CONDITIONS
def get_microstructure_params(condition='solution_annealed', extra=None):
    if condition not in PROCESSING_CONDITIONS:
        raise ValueError(f"Unknown condition: {condition}")
    cond = PROCESSING_CONDITIONS[condition]
    gb_fraction = 0.5 if condition == 'nanocrystalline' else 0.1
    # Build a realistic trap_list from TRAP_PROPERTIES
    trap_list = []
    for name, prop in TRAP_PROPERTIES.items():
        # Determine default density for this trap type
        default_density = None
        if 'density_range' in prop:
            default_density = prop['density_range'].get('typical', None)
            if default_density is None:
                # Fallback to any numeric entry in density_range
                for v in prop['density_range'].values():
                    if isinstance(v, (int, float)):
                        default_density = v
                        break
        # Allow overriding all trap densities via extra['trap_density']
        override_density = None
        if extra and 'trap_density' in extra:
            override_density = extra['trap_density']

        trap = {
            'name': name,
            'type': name,
            'binding_energy': prop.get('binding_energy', None),
            'density': override_density if override_density is not None else (default_density if default_density is not None else 0.0),
        }
        trap_list.append(trap)
    params = {
        'trap_density': cond['dislocation_density'],
        'grain_boundary_fraction': gb_fraction,
        'grain_size': cond['grain_size'],
        'grain_shape': 'equiaxed',
        'gb_type': 'HAGB',
        'trap_list': trap_list,
    }
    if extra:
        params.update(extra)
    return params
# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def save_test_results(test_name, results, results_dir='validation/results/Azeezlevel4'):
    os.makedirs(results_dir, exist_ok=True)
    filename = f"{results_dir}/{test_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    def convert_to_serializable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        if isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        if isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert_to_serializable(x) for x in obj]
        return obj
    json_results = convert_to_serializable(results)
    with open(filename, 'w') as f:
        json.dump(json_results, f, indent=2)
    print(f"  ✓ Results saved to: {filename}")

def print_test_header(test_name):
    print(f"\n{'='*80}")
    print(f"TEST: {test_name}")
    print(f"{'='*80}")

def print_test_result(passed, message=""):
    status = "✅ PASSED" if passed else "❌ FAILED"
    print(f"{status}: {message}")

# ============================================================================
# GROUP 1: BASIC FLUX CALCULATION TESTS
# ============================================================================

def test_basic_flux_level4(D=1e-8, K_s=1e-4, thickness=1e-3, P_up=1.0, P_down=0.0, temperature=1073):
    print_test_header("Level 4: Basic Flux Calculation (No Traps/GBs)")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 0, 'grain_boundary_fraction': 0})
    result = calculate_defective_metal_flux(
        D_lattice=D, K_s=K_s, thickness=thickness,
        P_up=P_up, P_down=P_down, temperature=temperature,
        microstructure_params=params
    )
    expected_flux = D * K_s / thickness * (np.sqrt(P_up) - np.sqrt(P_down))
    checks = {
        'flux_positive': result['flux'] > 0,
        'matches_level1': np.isclose(result['flux'], expected_flux, rtol=1e-3),
        'lattice_limited': result['regime_classification']['regime_detail'] == 'lattice_limited'
    }
    for check, passed in checks.items():
        print_test_result(passed, check)
    save_test_results('basic_flux_level4', {'params': locals(), 'result': result, 'checks': checks})

def test_zero_downstream_pressure_level4():
    print_test_header("Level 4: Zero Downstream Pressure")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 0})
    result = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    passed = abs(result['C_down']) < 1e-15
    print_test_result(passed, "C_down should be exactly zero")
    save_test_results('zero_downstream_level4', {'result': result, 'passed': passed})

def test_equal_pressures_level4():
    print_test_header("Level 4: Equal Pressures (No Driving Force)")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 0})
    result = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=1.0, temperature=1073,
        microstructure_params=params
    )
    passed = abs(result['flux']) < 1e-20
    print_test_result(passed, "Flux should be zero when no pressure gradient")
    save_test_results('equal_pressures_level4', {'result': result, 'passed': passed})

# ============================================================================
# GROUP 2: MICROSTRUCTURE EFFECTS
# ============================================================================

def test_trapping_effect_level4():
    print_test_header("Level 4: Trapping Effect")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 1e25, 'trap_energy': 0.6})
    result = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    passed = result['modification_factor'] < 0.5 and result['regime_classification']['regime_detail'] == 'traps_defect_limited'
    print_test_result(passed, "Trapping reduces flux and regime is traps_defect_limited")
    save_test_results('trapping_effect_level4', {'result': result, 'passed': passed})

def test_grain_boundary_effect_level4():
    print_test_header("Level 4: Grain Boundary Effect")
    params = get_microstructure_params('nanocrystalline', extra={'grain_boundary_fraction': 0.5, 'gb_enhancement': 10})
    result = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    passed = result['modification_factor'] > 1.0
    print_test_result(passed, "Grain boundary enhancement increases flux")
    save_test_results('grain_boundary_effect_level4', {'result': result, 'passed': passed})

def test_mixed_microstructure_level4():
    print_test_header("Level 4: Mixed Microstructure (Traps + GBs)")
    params = get_microstructure_params('nanocrystalline', extra={'trap_density': 1e25, 'trap_energy': 0.6, 'grain_boundary_fraction': 0.3, 'gb_enhancement': 5})
    result = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    print("Regime:", result['regime_classification']['regime_hierarchy'])
    save_test_results('mixed_microstructure_level4', {'result': result})

# ============================================================================
# GROUP 3: PHYSICAL DEPENDENCE TESTS
# ============================================================================

def test_thickness_dependence_level4():
    print_test_header("Level 4: Thickness Dependence")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 0})
    result_thin = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-4,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    result_thick = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-2,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    passed = result_thin['flux'] > result_thick['flux']
    print_test_result(passed, "Flux decreases with increasing thickness")
    save_test_results('thickness_dependence_level4', {'result_thin': result_thin, 'result_thick': result_thick, 'passed': passed})

def test_pressure_dependence_level4():
    print_test_header("Level 4: Pressure Dependence")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 0})
    result_low = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    result_high = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1e3, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    passed = result_high['flux'] > result_low['flux']
    print_test_result(passed, "Flux increases with pressure")
    save_test_results('pressure_dependence_level4', {'result_low': result_low, 'result_high': result_high, 'passed': passed})

def test_temperature_dependence_level4():
    print_test_header("Level 4: Temperature Dependence")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 0})
    result_cold = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=300,
        microstructure_params=params
    )
    result_hot = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=1200,
        microstructure_params=params
    )
    passed = result_hot['flux'] > result_cold['flux']
    print_test_result(passed, "Flux increases with temperature")
    save_test_results('temperature_dependence_level4', {'result_cold': result_cold, 'result_hot': result_hot, 'passed': passed})

# ============================================================================
# GROUP 4: REGIME CLASSIFICATION
# ============================================================================

def test_regime_classification_fields_level4():
    print_test_header("Level 4: Regime Classification Fields")
    params = get_microstructure_params('solution_annealed', extra={'trap_density': 1e25, 'trap_energy': 0.6})
    result = calculate_defective_metal_flux(
        D_lattice=1e-8, K_s=1e-4, thickness=1e-3,
        P_up=1.0, P_down=0.0, temperature=1073,
        microstructure_params=params
    )
    rc = result['regime_classification']
    checks = {
        'model_level': rc['model_level'] == 'Level 1,4',
        'base_regime': rc['base_regime'] == 'metal_limited',
        'regime_hierarchy': rc['regime_hierarchy'].startswith('metal_limited')
    }
    for check, passed in checks.items():
        print_test_result(passed, check)
    save_test_results('regime_classification_fields_level4', {'result': result, 'checks': checks})


def plot_flux_vs_pressure_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                  temperature=1073, P_min=DEFAULT_P_MIN, 
                                  P_max=DEFAULT_P_MAX, n_points=20):
    """
    Plot 1: Flux vs Pressure - Level 1 vs Level 4 (all modes) continuous comparison.
    """
    print_test_header("Plot 1: Flux vs Pressure (Level 1 vs Level 4)")
    
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    D_T = get_D_lattice(temperature)  # D at this temperature
    
    # Level 1: Perfect lattice (no microstructure)
    fluxes_level1 = []
    for P in pressures:
        result = calculate_simple_metal_flux(D_T, K_s, thickness, P, 0.0)
        fluxes_level1.append(result['flux'])
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only', 'none']
    mode_labels = {
        'both': 'Level 4: Both (GB + Trapping)',
        'gb_only': 'Level 4: GB Enhancement Only',
        'trapping_only': 'Level 4: Trapping Only',
        'none': 'Level 4: No Effects (=Level 1)'
    }
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange', 'none': 'purple'}
    
    fluxes_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')  # Use solution annealed
    
    for P in pressures:
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode=mode
            )
            fluxes_level4[mode].append(result['flux'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    
    # Level 1 baseline (dashed black)
    plt.loglog(pressures, fluxes_level1, 'k--', linewidth=2, label='Level 1: Perfect Lattice')
    
    # Level 4 modes
    for mode in modes:
        plt.loglog(pressures, fluxes_level4[mode], 'o-', color=mode_colors[mode], 
                   label=mode_labels[mode], markersize=4)
    
    plt.xlabel('Upstream Pressure (Pa)', fontsize=12)
    plt.ylabel('Flux (mol/m²/s)', fontsize=12)
    plt.title('Flux vs Pressure: Level 1 vs Level 4 Comparison', fontsize=14)
    plt.grid(True, which='both', alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot1_flux_vs_pressure_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'pressures': pressures, 'fluxes_level1': fluxes_level1, 'fluxes_level4': fluxes_level4}

def plot_flux_vs_temperature_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                     P_up=DEFAULT_P_UP, T_min=DEFAULT_T_MIN, 
                                     T_max=DEFAULT_T_MAX, n_points=20):
    """
    Plot 2: Flux vs Temperature - Level 1 vs Level 4 (all modes) continuous comparison.
    Uses Arrhenius D(T) = D_0 * exp(-E_D/RT) for temperature dependence.
    """
    print_test_header("Plot 2: Flux vs Temperature (Level 1 vs Level 4)")
    
    temperatures = np.linspace(T_min, T_max, n_points)
    
    # Level 1: Perfect lattice with Arrhenius D(T)
    fluxes_level1 = []
    for T in temperatures:
        D_T = get_D_lattice(T)
        result = calculate_simple_metal_flux(D_T, K_s, thickness, P_up, 0.0)
        fluxes_level1.append(result['flux'])
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only', 'none']
    mode_labels = {
        'both': 'Level 4: Both (GB + Trapping)',
        'gb_only': 'Level 4: GB Enhancement Only',
        'trapping_only': 'Level 4: Trapping Only',
        'none': 'Level 4: No Effects (=Level 1)'
    }
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange', 'none': 'purple'}
    
    fluxes_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for T in temperatures:
        D_T = get_D_lattice(T)
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=0.0, temperature=T,
                microstructure_params=params, mode=mode
            )
            fluxes_level4[mode].append(result['flux'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    
    # Level 1 baseline (dashed black)
    plt.semilogy(temperatures, fluxes_level1, 'k--', linewidth=2, label='Level 1: Perfect Lattice')
    
    # Level 4 modes
    for mode in modes:
        plt.semilogy(temperatures, fluxes_level4[mode], 'o-', color=mode_colors[mode], 
                     label=mode_labels[mode], markersize=4)
    
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Flux (mol/m²/s)', fontsize=12)
    plt.title(f'Flux vs Temperature: Level 1 vs Level 4 (P={P_up} Pa, L={thickness*1e3:.1f} mm)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot2_flux_vs_temperature_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'temperatures': temperatures, 'fluxes_level1': fluxes_level1, 'fluxes_level4': fluxes_level4}

def plot_sensitivity_to_trap_density_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                             P_up=DEFAULT_P_UP, temperature=1073, n_points=20):
    """
    Plot 3: Flux vs Trap Density - Level 1 (constant) vs Level 4 continuous comparison.
    """
    print_test_header("Plot 3: Sensitivity to Trap Density (Level 1 vs Level 4)")
    
    trap_densities = np.logspace(20, 26, n_points)  # 1e20 to 1e26 m^-3
    D_T = get_D_lattice(temperature)
    
    # Level 1: Constant (no traps)
    result_level1 = calculate_simple_metal_flux(D_T, K_s, thickness, P_up, 0.0)
    fluxes_level1 = [result_level1['flux']] * n_points
    
    # Level 4 with varying trap density
    fluxes_level4 = []
    mod_factors = []
    
    for td in trap_densities:
        params = get_microstructure_params('solution_annealed', extra={'trap_density': td})
        result = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=temperature,
            microstructure_params=params, mode='trapping_only'
        )
        fluxes_level4.append(result['flux'])
        mod_factors.append(result['modification_factor'])
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Flux comparison
    ax1.semilogx(trap_densities, fluxes_level1, 'k--', linewidth=2, label='Level 1: Perfect Lattice')
    ax1.semilogx(trap_densities, fluxes_level4, 'o-', color='orange', label='Level 4: Trapping Only', markersize=4)
    ax1.set_xlabel('Trap Density (m⁻³)', fontsize=12)
    ax1.set_ylabel('Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux vs Trap Density', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Right: Modification factor
    ax2.semilogx(trap_densities, mod_factors, 's-', color='blue', markersize=4)
    ax2.axhline(y=1.0, color='k', linestyle='--', label='No modification')
    ax2.set_xlabel('Trap Density (m⁻³)', fontsize=12)
    ax2.set_ylabel('Modification Factor (D_eff/D_lattice)', fontsize=12)
    ax2.set_title('D_eff Reduction vs Trap Density', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    fname = f'{RESULTS_DIR}/plot3_sensitivity_trap_density_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'trap_densities': trap_densities, 'fluxes_level1': fluxes_level1, 
            'fluxes_level4': fluxes_level4, 'modification_factors': mod_factors}

def plot_sensitivity_to_grain_size_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                           P_up=DEFAULT_P_UP, temperature=1073, n_points=20):
    """
    Plot 4: Flux vs Grain Size - Level 1 (constant) vs Level 4 continuous comparison.
    Smaller grain size = more GB = more enhancement.
    """
    print_test_header("Plot 4: Sensitivity to Grain Size (Level 1 vs Level 4)")
    
    grain_sizes = np.logspace(-8, -4, n_points)  # 10 nm to 100 μm
    D_T = get_D_lattice(temperature)
    
    # Level 1: Constant (no GB effects)
    result_level1 = calculate_simple_metal_flux(D_T, K_s, thickness, P_up, 0.0)
    fluxes_level1 = [result_level1['flux']] * n_points
    
    # Level 4 with varying grain size
    fluxes_level4 = []
    mod_factors = []
    
    for gs in grain_sizes:
        params = get_microstructure_params('solution_annealed', extra={
            'grain_size': gs,
            'trap_density': 0  # No traps, only GB effect
        })
        result = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=temperature,
            microstructure_params=params, mode='gb_only'
        )
        fluxes_level4.append(result['flux'])
        mod_factors.append(result['modification_factor'])
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Flux comparison
    ax1.semilogx(grain_sizes * 1e6, fluxes_level1, 'k--', linewidth=2, label='Level 1: Perfect Lattice')
    ax1.semilogx(grain_sizes * 1e6, fluxes_level4, 'o-', color='green', label='Level 4: GB Only', markersize=4)
    ax1.set_xlabel('Grain Size (μm)', fontsize=12)
    ax1.set_ylabel('Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux vs Grain Size', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Right: Modification factor
    ax2.semilogx(grain_sizes * 1e6, mod_factors, 's-', color='blue', markersize=4)
    ax2.axhline(y=1.0, color='k', linestyle='--', label='No modification')
    ax2.set_xlabel('Grain Size (μm)', fontsize=12)
    ax2.set_ylabel('Modification Factor (D_eff/D_lattice)', fontsize=12)
    ax2.set_title('D_eff Enhancement vs Grain Size', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    fname = f'{RESULTS_DIR}/plot4_sensitivity_grain_size_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'grain_sizes': grain_sizes, 'fluxes_level1': fluxes_level1, 
            'fluxes_level4': fluxes_level4, 'modification_factors': mod_factors}


# ============================================================================
# GROUP 5: ADDITIONAL COMPARISON PLOTS (Level 1 vs Level 4)
# ============================================================================

def plot_flux_vs_thickness_level4(K_s=DEFAULT_K_S, P_up=DEFAULT_P_UP, temperature=1073,
                                   L_min=DEFAULT_L_MIN, L_max=DEFAULT_L_MAX, n_points=20):
    """
    Plot 5: Flux vs Thickness - Level 1 vs Level 4 (all modes) continuous comparison.
    """
    print_test_header("Plot 5: Flux vs Thickness (Level 1 vs Level 4)")
    
    thicknesses = np.logspace(np.log10(L_min), np.log10(L_max), n_points)
    D_T = get_D_lattice(temperature)
    
    # Level 1: Perfect lattice
    fluxes_level1 = []
    for L in thicknesses:
        result = calculate_simple_metal_flux(D_T, K_s, L, P_up, 0.0)
        fluxes_level1.append(result['flux'])
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only', 'none']
    mode_labels = {
        'both': 'Level 4: Both',
        'gb_only': 'Level 4: GB Only',
        'trapping_only': 'Level 4: Trapping Only',
        'none': 'Level 4: None'
    }
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange', 'none': 'purple'}
    
    fluxes_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for L in thicknesses:
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=L,
                P_up=P_up, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode=mode
            )
            fluxes_level4[mode].append(result['flux'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.loglog(thicknesses * 1e3, fluxes_level1, 'k--', linewidth=2, label='Level 1: Perfect Lattice')
    for mode in modes:
        plt.loglog(thicknesses * 1e3, fluxes_level4[mode], 'o-', color=mode_colors[mode], 
                   label=mode_labels[mode], markersize=4)
    
    plt.xlabel('Thickness (mm)', fontsize=12)
    plt.ylabel('Flux (mol/m²/s)', fontsize=12)
    plt.title('Flux vs Thickness: Level 1 vs Level 4 Comparison', fontsize=14)
    plt.grid(True, which='both', alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot5_flux_vs_thickness_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'thicknesses': thicknesses, 'fluxes_level1': fluxes_level1, 'fluxes_level4': fluxes_level4}


def plot_arrhenius_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, P_up=DEFAULT_P_UP,
                          T_min=DEFAULT_T_MIN, T_max=DEFAULT_T_MAX, n_points=20):
    """
    Plot 6: Arrhenius Plot - ln(Flux) vs 1000/T for Level 1 vs Level 4.
    Uses Arrhenius D(T) - should show linear behavior for Level 1.
    """
    print_test_header("Plot 6: Arrhenius Plot (Level 1 vs Level 4)")
    
    temperatures = np.linspace(T_min, T_max, n_points)
    inv_T = 1000.0 / temperatures  # 1000/T for better scale
    
    # Level 1 with Arrhenius D(T)
    fluxes_level1 = []
    for T in temperatures:
        D_T = get_D_lattice(T)
        result = calculate_simple_metal_flux(D_T, K_s, thickness, P_up, 0.0)
        fluxes_level1.append(result['flux'])
    ln_flux_level1 = np.log(fluxes_level1)
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only']
    mode_labels = {'both': 'Level 4: Both', 'gb_only': 'Level 4: GB Only', 'trapping_only': 'Level 4: Trapping Only'}
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange'}
    
    ln_flux_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for T in temperatures:
        D_T = get_D_lattice(T)
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=0.0, temperature=T,
                microstructure_params=params, mode=mode
            )
            ln_flux_level4[mode].append(np.log(result['flux']))
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.plot(inv_T, ln_flux_level1, 'k--', linewidth=2, label=f'Level 1: Perfect Lattice (E_D={DEFAULT_E_D/1000:.0f} kJ/mol)')
    for mode in modes:
        plt.plot(inv_T, ln_flux_level4[mode], 'o-', color=mode_colors[mode], 
                 label=mode_labels[mode], markersize=4)
    
    plt.xlabel('1000/T (K⁻¹)', fontsize=12)
    plt.ylabel('ln(Flux)', fontsize=12)
    plt.title('Arrhenius Plot: Level 1 vs Level 4 Comparison', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot6_arrhenius_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'inv_T': inv_T, 'temperatures': temperatures, 'ln_flux_level1': ln_flux_level1, 'ln_flux_level4': ln_flux_level4}


def plot_sqrt_P_scaling_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, temperature=1073,
                                P_min=DEFAULT_P_MIN, P_max=DEFAULT_P_MAX, n_points=20):
    """
    Plot 7: √P Scaling Verification - Flux vs √P should be linear for Sieverts' law.
    """
    print_test_header("Plot 7: √P Scaling Verification (Level 1 vs Level 4)")
    
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    sqrt_P = np.sqrt(pressures)
    D_T = get_D_lattice(temperature)
    
    # Level 1
    fluxes_level1 = []
    for P in pressures:
        result = calculate_simple_metal_flux(D_T, K_s, thickness, P, 0.0)
        fluxes_level1.append(result['flux'])
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only']
    mode_labels = {'both': 'Level 4: Both', 'gb_only': 'Level 4: GB Only', 'trapping_only': 'Level 4: Trapping Only'}
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange'}
    
    fluxes_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for P in pressures:
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode=mode
            )
            fluxes_level4[mode].append(result['flux'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.plot(sqrt_P, fluxes_level1, 'k--', linewidth=2, label='Level 1: Perfect Lattice')
    for mode in modes:
        plt.plot(sqrt_P, fluxes_level4[mode], 'o-', color=mode_colors[mode], 
                 label=mode_labels[mode], markersize=4)
    
    plt.xlabel('√P (Pa⁰·⁵)', fontsize=12)
    plt.ylabel('Flux (mol/m²/s)', fontsize=12)
    plt.title('√P Scaling Verification: Should be Linear', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot7_sqrt_P_scaling_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'sqrt_P': sqrt_P, 'pressures': pressures, 'fluxes_level1': fluxes_level1, 'fluxes_level4': fluxes_level4}


def plot_D_eff_vs_pressure_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, temperature=1073,
                                   P_min=DEFAULT_P_MIN, P_max=DEFAULT_P_MAX, n_points=20):
    """
    Plot 8: D_eff vs Pressure - Shows how effective diffusivity changes with pressure.
    Level 1: D_eff = D_lattice (constant)
    Level 4 gb_only: constant (no concentration dependence)
    Level 4 trapping_only: varies with concentration (pressure)
    """
    print_test_header("Plot 8: D_eff vs Pressure (Level 1 vs Level 4)")
    
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    D_T = get_D_lattice(temperature)
    
    # Level 1: D_eff = D_lattice(T) at this temperature
    D_eff_level1 = [D_T] * n_points
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only']
    mode_labels = {'both': 'Level 4: Both', 'gb_only': 'Level 4: GB Only', 'trapping_only': 'Level 4: Trapping Only'}
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange'}
    
    D_eff_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for P in pressures:
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode=mode
            )
            D_eff_level4[mode].append(result['D_eff'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.semilogx(pressures, D_eff_level1, 'k--', linewidth=2, label=f'Level 1: D_lattice(T={temperature}K)')
    for mode in modes:
        plt.semilogx(pressures, D_eff_level4[mode], 'o-', color=mode_colors[mode], 
                     label=mode_labels[mode], markersize=4)
    
    plt.xlabel('Upstream Pressure (Pa)', fontsize=12)
    plt.ylabel('D_eff (m²/s)', fontsize=12)
    plt.title('Effective Diffusivity vs Pressure: Level 1 vs Level 4', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot8_D_eff_vs_pressure_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'pressures': pressures, 'D_eff_level1': D_eff_level1, 'D_eff_level4': D_eff_level4}


def plot_D_eff_vs_temperature_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, P_up=DEFAULT_P_UP,
                                      T_min=DEFAULT_T_MIN, T_max=DEFAULT_T_MAX, n_points=20):
    """
    Plot 9: D_eff vs Temperature - Shows how effective diffusivity changes with T.
    Uses Arrhenius D(T) for temperature dependence.
    """
    print_test_header("Plot 9: D_eff vs Temperature (Level 1 vs Level 4)")
    
    temperatures = np.linspace(T_min, T_max, n_points)
    
    # Level 1: D_eff = D_lattice(T) with Arrhenius dependence
    D_eff_level1 = [get_D_lattice(T) for T in temperatures]
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only']
    mode_labels = {'both': 'Level 4: Both', 'gb_only': 'Level 4: GB Only', 'trapping_only': 'Level 4: Trapping Only'}
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange'}
    
    D_eff_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for T in temperatures:
        D_T = get_D_lattice(T)
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=0.0, temperature=T,
                microstructure_params=params, mode=mode
            )
            D_eff_level4[mode].append(result['D_eff'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.semilogy(temperatures, D_eff_level1, 'k--', linewidth=2, label='Level 1: D_lattice(T) Arrhenius')
    for mode in modes:
        plt.semilogy(temperatures, D_eff_level4[mode], 'o-', color=mode_colors[mode], 
                     label=mode_labels[mode], markersize=4)
    
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('D_eff (m²/s)', fontsize=12)
    plt.title(f'Effective Diffusivity vs Temperature (P={P_up} Pa)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot9_D_eff_vs_temperature_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'temperatures': temperatures, 'D_eff_level1': D_eff_level1, 'D_eff_level4': D_eff_level4}


def plot_modification_factor_vs_pressure_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                                 temperature=1073, P_min=DEFAULT_P_MIN, 
                                                 P_max=DEFAULT_P_MAX, n_points=20):
    """
    Plot 10: Modification Factor (D_eff/D_lattice) vs Pressure.
    Level 1: Always = 1 (baseline)
    """
    print_test_header("Plot 10: Modification Factor vs Pressure")
    
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    D_T = get_D_lattice(temperature)
    
    # Level 1: Modification factor = 1 (constant)
    mod_factor_level1 = [1.0] * n_points
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only']
    mode_labels = {'both': 'Level 4: Both', 'gb_only': 'Level 4: GB Only', 'trapping_only': 'Level 4: Trapping Only'}
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange'}
    
    mod_factor_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for P in pressures:
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode=mode
            )
            mod_factor_level4[mode].append(result['modification_factor'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.semilogx(pressures, mod_factor_level1, 'k--', linewidth=2, label='Level 1: No modification (=1)')
    for mode in modes:
        plt.semilogx(pressures, mod_factor_level4[mode], 'o-', color=mode_colors[mode], 
                     label=mode_labels[mode], markersize=4)
    
    plt.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
    plt.xlabel('Upstream Pressure (Pa)', fontsize=12)
    plt.ylabel('Modification Factor (D_eff/D_lattice)', fontsize=12)
    plt.title('Modification Factor vs Pressure', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot10_mod_factor_vs_pressure_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'pressures': pressures, 'mod_factor_level1': mod_factor_level1, 'mod_factor_level4': mod_factor_level4}


def plot_mode_comparison_continuous_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                            temperature=1073, P_min=DEFAULT_P_MIN, 
                                            P_max=DEFAULT_P_MAX, n_points=20):
    """
    Plot 11: Mode Comparison - Continuous flux ratio to Level 1 across pressure.
    Shows how each Level 4 mode deviates from Level 1 baseline.
    """
    print_test_header("Plot 11: Mode Comparison - Flux Ratio to Level 1 (Continuous)")
    
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    D_T = get_D_lattice(temperature)
    
    # Level 1 baseline
    fluxes_level1 = []
    for P in pressures:
        result = calculate_simple_metal_flux(D_T, K_s, thickness, P, 0.0)
        fluxes_level1.append(result['flux'])
    
    # Level 4 modes
    modes = ['none', 'gb_only', 'trapping_only', 'both']
    mode_labels = {
        'none': 'Level 4: None (=Level 1)',
        'gb_only': 'Level 4: GB Only',
        'trapping_only': 'Level 4: Trapping Only',
        'both': 'Level 4: Both'
    }
    mode_colors = {'none': 'purple', 'gb_only': 'green', 'trapping_only': 'orange', 'both': 'red'}
    
    flux_ratios = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for i, P in enumerate(pressures):
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode=mode
            )
            flux_ratios[mode].append(result['flux'] / fluxes_level1[i])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.axhline(y=1.0, color='k', linestyle='--', linewidth=2, label='Level 1 Baseline (=1)')
    for mode in modes:
        plt.semilogx(pressures, flux_ratios[mode], 'o-', color=mode_colors[mode], 
                     label=mode_labels[mode], markersize=4)
    
    plt.xlabel('Upstream Pressure (Pa)', fontsize=12)
    plt.ylabel('Flux Ratio (J_Level4 / J_Level1)', fontsize=12)
    plt.title('Mode Comparison: Flux Ratio to Level 1 Baseline', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot11_mode_comparison_continuous_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'pressures': pressures, 'flux_ratios': flux_ratios}


def plot_modification_factor_vs_temperature_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                                    P_up=DEFAULT_P_UP, T_min=DEFAULT_T_MIN, 
                                                    T_max=DEFAULT_T_MAX, n_points=20):
    """
    Plot 12: Modification Factor (D_eff/D_lattice) vs Temperature.
    Uses Arrhenius D(T) for temperature dependence.
    """
    print_test_header("Plot 12: Modification Factor vs Temperature")
    
    temperatures = np.linspace(T_min, T_max, n_points)
    
    # Level 1: Modification factor = 1 (constant - by definition)
    mod_factor_level1 = [1.0] * n_points
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only']
    mode_labels = {'both': 'Level 4: Both', 'gb_only': 'Level 4: GB Only', 'trapping_only': 'Level 4: Trapping Only'}
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange'}
    
    mod_factor_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for T in temperatures:
        D_T = get_D_lattice(T)
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=0.0, temperature=T,
                microstructure_params=params, mode=mode
            )
            mod_factor_level4[mode].append(result['modification_factor'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.plot(temperatures, mod_factor_level1, 'k--', linewidth=2, label='Level 1: No modification (=1)')
    for mode in modes:
        plt.plot(temperatures, mod_factor_level4[mode], 'o-', color=mode_colors[mode], 
                 label=mode_labels[mode], markersize=4)
    
    plt.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Modification Factor (D_eff/D_lattice)', fontsize=12)
    plt.title(f'Modification Factor vs Temperature (P={P_up} Pa)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot12_mod_factor_vs_temperature_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'temperatures': temperatures, 'mod_factor_level1': mod_factor_level1, 'mod_factor_level4': mod_factor_level4}


def plot_permeability_check_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, temperature=1073,
                                    P_min=DEFAULT_P_MIN, P_max=DEFAULT_P_MAX, n_points=20):
    """
    Plot 13: J×L vs P (Permeability Check).
    For Sieverts' law: J = D*K_s/L * sqrt(P), so J*L = D*K_s*sqrt(P)
    J*L should scale linearly with sqrt(P).
    """
    print_test_header("Plot 13: J×L vs P (Permeability Check)")
    
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    D_T = get_D_lattice(temperature)
    
    # Level 1
    JL_level1 = []
    for P in pressures:
        result = calculate_simple_metal_flux(D_T, K_s, thickness, P, 0.0)
        JL_level1.append(result['flux'] * thickness)
    
    # Level 4 modes
    modes = ['both', 'gb_only', 'trapping_only']
    mode_labels = {'both': 'Level 4: Both', 'gb_only': 'Level 4: GB Only', 'trapping_only': 'Level 4: Trapping Only'}
    mode_colors = {'both': 'red', 'gb_only': 'green', 'trapping_only': 'orange'}
    
    JL_level4 = {mode: [] for mode in modes}
    params = get_microstructure_params('solution_annealed')
    
    for P in pressures:
        for mode in modes:
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode=mode
            )
            JL_level4[mode].append(result['flux'] * thickness)
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.loglog(pressures, JL_level1, 'k--', linewidth=2, label='Level 1: Perfect Lattice')
    for mode in modes:
        plt.loglog(pressures, JL_level4[mode], 'o-', color=mode_colors[mode], 
                   label=mode_labels[mode], markersize=4)
    
    plt.xlabel('Upstream Pressure (Pa)', fontsize=12)
    plt.ylabel('J × L (mol/m/s)', fontsize=12)
    plt.title('Permeability Check: J×L vs P', fontsize=14)
    plt.grid(True, which='both', alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot13_permeability_check_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'pressures': pressures, 'JL_level1': JL_level1, 'JL_level4': JL_level4}


def plot_trap_occupancy_vs_pressure_level4(K_s=1e-4, thickness=1e-3, temperature=1073,
                                            P_min=0.01, P_max=100, n_points=20):
    """
    Plot 14: Trap Occupancy (θ) vs Pressure.
    θ = K*C_L / (1 + K*C_L) where K = exp(E_b/RT)
    Level 1: θ = 0 (no traps)
    Note: θ depends on K_s, P, and T, not on D(T) directly.
    """
    print_test_header("Plot 14: Trap Occupancy (θ) vs Pressure")
    
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    
    # Level 1: θ = 0 (no traps)
    theta_level1 = [0.0] * n_points
    
    # Level 4 with trapping
    theta_level4 = []
    
    # Physical constants
    K_trap = np.exp(DEFAULT_E_BINDING / (R_GAS * temperature))
    
    for P in pressures:
        # Lattice concentration from Sieverts' law
        C_L = K_s * np.sqrt(P)
        # Trap occupancy
        theta = K_trap * C_L / (1 + K_trap * C_L)
        theta_level4.append(theta)
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.semilogx(pressures, theta_level1, 'k--', linewidth=2, label='Level 1: θ = 0 (no traps)')
    plt.semilogx(pressures, theta_level4, 'o-', color='orange', label='Level 4: θ (Oriani equilibrium)', markersize=4)
    
    plt.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5, label='θ = 0.5 (half-filled)')
    plt.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
    plt.xlabel('Upstream Pressure (Pa)', fontsize=12)
    plt.ylabel('Trap Occupancy θ', fontsize=12)
    plt.title(f'Trap Occupancy vs Pressure (T = {temperature} K)', fontsize=14)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot14_trap_occupancy_vs_pressure_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'pressures': pressures, 'theta_level1': theta_level1, 'theta_level4': theta_level4}


def plot_trap_occupancy_vs_temperature_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                               P_up=DEFAULT_P_UP, T_min=DEFAULT_T_MIN, 
                                               T_max=DEFAULT_T_MAX, n_points=20):
    """
    Plot 15: Trap Occupancy (θ) vs Temperature.
    θ decreases at high T as thermal detrapping increases.
    Level 1: θ = 0 (no traps)
    Note: θ depends on K_s and T, not on D(T) directly.
    """
    print_test_header("Plot 15: Trap Occupancy (θ) vs Temperature")
    
    temperatures = np.linspace(T_min, T_max, n_points)
    
    # Level 1: θ = 0 (no traps)
    theta_level1 = [0.0] * n_points
    
    # Level 4 with trapping
    theta_level4 = []
    
    for T in temperatures:
        # Lattice concentration from Sieverts' law
        C_L = K_s * np.sqrt(P_up)
        # Temperature-dependent trap equilibrium constant
        K_trap = np.exp(DEFAULT_E_BINDING / (R_GAS * T))
        # Trap occupancy
        theta = K_trap * C_L / (1 + K_trap * C_L)
        theta_level4.append(theta)
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.plot(temperatures, theta_level1, 'k--', linewidth=2, label='Level 1: θ = 0 (no traps)')
    plt.plot(temperatures, theta_level4, 'o-', color='orange', label='Level 4: θ (Oriani equilibrium)', markersize=4)
    
    plt.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5, label='θ = 0.5 (half-filled)')
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Trap Occupancy θ', fontsize=12)
    plt.title(f'Trap Occupancy vs Temperature (P = {P_up} Pa)', fontsize=14)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot15_trap_occupancy_vs_temperature_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'temperatures': temperatures, 'theta_level1': theta_level1, 'theta_level4': theta_level4}


def plot_mode_dominance_vs_temperature_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                               P_up=DEFAULT_P_UP, T_min=DEFAULT_T_MIN, 
                                               T_max=DEFAULT_T_MAX, n_points=20):
    """
    Plot 16: Mode Effect Magnitude vs Temperature.
    Shows which mode (GB enhancement or trapping) dominates at different temperatures.
    Uses Arrhenius D(T) for temperature dependence.
    """
    print_test_header("Plot 16: Mode Dominance vs Temperature")
    
    temperatures = np.linspace(T_min, T_max, n_points)
    
    # Level 4 modes
    params = get_microstructure_params('solution_annealed')
    
    gb_effect = []  # (J_gb_only - J_level1) / J_level1
    trap_effect = []  # (J_trap_only - J_level1) / J_level1
    combined_effect = []  # (J_both - J_level1) / J_level1
    
    for T in temperatures:
        D_T = get_D_lattice(T)
        
        # Level 1 baseline at this temperature
        result_level1 = calculate_simple_metal_flux(D_T, K_s, thickness, P_up, 0.0)
        flux_level1 = result_level1['flux']
        
        result_gb = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=T,
            microstructure_params=params, mode='gb_only'
        )
        result_trap = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=T,
            microstructure_params=params, mode='trapping_only'
        )
        result_both = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=T,
            microstructure_params=params, mode='both'
        )
        
        gb_effect.append((result_gb['flux'] - flux_level1) / flux_level1 * 100)
        trap_effect.append((result_trap['flux'] - flux_level1) / flux_level1 * 100)
        combined_effect.append((result_both['flux'] - flux_level1) / flux_level1 * 100)
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.axhline(y=0, color='k', linestyle='--', linewidth=2, label='Level 1 Baseline (0%)')
    plt.plot(temperatures, gb_effect, 'o-', color='green', label='GB Only Effect', markersize=4)
    plt.plot(temperatures, trap_effect, 'o-', color='orange', label='Trapping Only Effect', markersize=4)
    plt.plot(temperatures, combined_effect, 'o-', color='red', label='Combined Effect', markersize=4)
    
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('% Change from Level 1 Baseline', fontsize=12)
    plt.title('Mode Effect Magnitude vs Temperature', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot16_mode_dominance_vs_temperature_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'temperatures': temperatures, 'gb_effect': gb_effect, 'trap_effect': trap_effect, 'combined_effect': combined_effect}


def plot_D_eff_profile_through_thickness_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                                  P_up=DEFAULT_P_UP, temperature=1073, n_points=50):
    """
    Plot 17: D_eff Profile Through Thickness.
    Shows position-dependent D_eff(x) due to concentration-dependent trapping.
    Level 1: D_eff = D_lattice(T) (constant through thickness)
    Uses Arrhenius D(T) for temperature dependence.
    """
    print_test_header("Plot 17: D_eff Profile Through Thickness")
    
    positions = np.linspace(0, thickness, n_points)
    x_normalized = positions / thickness  # 0 to 1
    
    # Get D_lattice at this temperature
    D = get_D_lattice(temperature)
    
    # Level 1: D_eff = D_lattice (constant through thickness)
    D_eff_level1 = [D] * n_points
    
    # Level 4: D_eff varies with position due to concentration gradient
    # C(x) = C_up - (C_up - C_down) * x/L
    C_up = K_s * np.sqrt(P_up)
    C_down = 0.0
    
    # Physical constants
    K_trap = np.exp(DEFAULT_E_BINDING / (R_GAS * temperature))
    N_t = DEFAULT_TRAP_DENSITY  # Trap density
    
    D_eff_level4_trap = []
    D_eff_level4_gb = []
    D_eff_level4_both = []
    
    # GB enhancement factor (constant, not position-dependent)
    f_gb = 0.1  # GB fraction
    D_gb_factor = 10  # GB enhancement
    D_eff_gb = D * (1 + f_gb * (D_gb_factor - 1))
    
    for x in positions:
        # Local concentration
        C_local = C_up - (C_up - C_down) * (x / thickness)
        
        # Trap occupancy at this position
        theta = K_trap * C_local / (1 + K_trap * C_local) if C_local > 0 else 0
        
        # Effective D with trapping: D_eff = D / (1 + N_t/N_L * dθ/dC_L)
        # Simplified: D_eff ≈ D / (1 + N_t*K/(N_L*(1+K*C_L)^2))
        N_L = 1e29  # Lattice site density
        denom = 1 + (N_t / N_L) * K_trap / (1 + K_trap * C_local)**2 if C_local > 0 else 1
        D_eff_trap = D / denom
        
        D_eff_level4_trap.append(D_eff_trap)
        D_eff_level4_gb.append(D_eff_gb)
        D_eff_level4_both.append(D_eff_gb / denom)
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.plot(x_normalized, D_eff_level1, 'k--', linewidth=2, label='Level 1: D_lattice (constant)')
    plt.plot(x_normalized, D_eff_level4_gb, 'o-', color='green', label='Level 4: GB Only', markersize=3)
    plt.plot(x_normalized, D_eff_level4_trap, 'o-', color='orange', label='Level 4: Trapping Only', markersize=3)
    plt.plot(x_normalized, D_eff_level4_both, 'o-', color='red', label='Level 4: Both', markersize=3)
    
    plt.xlabel('Normalized Position x/L (0=upstream, 1=downstream)', fontsize=12)
    plt.ylabel('D_eff (m²/s)', fontsize=12)
    plt.title('D_eff Profile Through Thickness', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot17_D_eff_profile_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'positions': positions, 'x_normalized': x_normalized, 'D_eff_level1': D_eff_level1,
            'D_eff_level4_trap': D_eff_level4_trap, 'D_eff_level4_gb': D_eff_level4_gb, 'D_eff_level4_both': D_eff_level4_both}


def plot_concentration_profile_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, 
                                       P_up=DEFAULT_P_UP, temperature=1073, n_points=50):
    """
    Plot 18: Concentration Profile C(x) Through Thickness.
    Level 1: Linear profile C(x) = C_up * (1 - x/L)
    Level 4: May deviate due to trapping effects
    Note: Concentration profile depends on K_s, not D(T) directly.
    """
    print_test_header("Plot 18: Concentration Profile Through Thickness")
    
    positions = np.linspace(0, thickness, n_points)
    x_normalized = positions / thickness  # 0 to 1
    
    # Boundary concentrations
    C_up = K_s * np.sqrt(P_up)
    C_down = 0.0
    
    # Level 1: Linear profile
    C_level1 = [C_up * (1 - x/thickness) for x in positions]
    
    # Level 4: With trapping, total hydrogen = C_L + C_T
    # For steady-state, C_L still follows linear profile
    # But total H includes trapped H: C_total = C_L + N_t * θ
    
    K_trap = np.exp(DEFAULT_E_BINDING / (R_GAS * temperature))
    N_t = DEFAULT_TRAP_DENSITY  # mol/m³
    
    C_lattice_level4 = []
    C_total_level4 = []
    
    for x in positions:
        C_L = C_up * (1 - x/thickness)  # Lattice H
        theta = K_trap * C_L / (1 + K_trap * C_L) if C_L > 0 else 0
        C_T = N_t * theta  # Trapped H (in mol/m³)
        
        C_lattice_level4.append(C_L)
        C_total_level4.append(C_L + C_T)
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Lattice concentration
    ax1.plot(x_normalized, C_level1, 'k--', linewidth=2, label='Level 1: Linear')
    ax1.plot(x_normalized, C_lattice_level4, 'o-', color='blue', label='Level 4: Lattice C_L', markersize=3)
    ax1.set_xlabel('Normalized Position x/L', fontsize=12)
    ax1.set_ylabel('Lattice Concentration C_L (mol/m³)', fontsize=12)
    ax1.set_title('Lattice Hydrogen Concentration', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Right: Total concentration (lattice + trapped)
    ax2.plot(x_normalized, C_level1, 'k--', linewidth=2, label='Level 1: C_total = C_L')
    ax2.plot(x_normalized, C_total_level4, 'o-', color='red', label='Level 4: C_total = C_L + C_T', markersize=3)
    ax2.plot(x_normalized, C_lattice_level4, 'o-', color='blue', label='Level 4: C_L only', markersize=3, alpha=0.5)
    ax2.set_xlabel('Normalized Position x/L', fontsize=12)
    ax2.set_ylabel('Total Concentration (mol/m³)', fontsize=12)
    ax2.set_title('Total Hydrogen Concentration', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    fname = f'{RESULTS_DIR}/plot18_concentration_profile_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'positions': positions, 'x_normalized': x_normalized, 'C_level1': C_level1,
            'C_lattice_level4': C_lattice_level4, 'C_total_level4': C_total_level4}


def plot_competing_effects_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, P_up=DEFAULT_P_UP,
                                   T_min=DEFAULT_T_MIN, T_max=DEFAULT_T_MAX, n_points=20):
    """
    Plot 19: Competing Effects - GB Enhancement vs Trapping Reduction.
    Shows modification factors separately across temperature.
    Uses Arrhenius D(T) for temperature dependence.
    """
    print_test_header("Plot 19: Competing Effects (GB vs Trapping)")
    
    temperatures = np.linspace(T_min, T_max, n_points)
    
    # Level 1: mod_factor = 1 (baseline - by definition)
    mod_factor_level1 = [1.0] * n_points
    
    # Level 4
    params = get_microstructure_params('solution_annealed')
    
    mod_gb = []
    mod_trap = []
    mod_combined = []
    
    for T in temperatures:
        D_T = get_D_lattice(T)
        result_gb = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=T,
            microstructure_params=params, mode='gb_only'
        )
        result_trap = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=T,
            microstructure_params=params, mode='trapping_only'
        )
        result_both = calculate_defective_metal_flux(
            D_lattice=D_T, K_s=K_s, thickness=thickness,
            P_up=P_up, P_down=0.0, temperature=T,
            microstructure_params=params, mode='both'
        )
        
        mod_gb.append(result_gb['modification_factor'])
        mod_trap.append(result_trap['modification_factor'])
        mod_combined.append(result_both['modification_factor'])
    
    # Plot
    plt.figure(figsize=(10, 7))
    plt.axhline(y=1.0, color='k', linestyle='--', linewidth=2, label='Level 1: No modification (=1)')
    plt.plot(temperatures, mod_gb, 'o-', color='green', label='GB Enhancement Factor', markersize=4)
    plt.plot(temperatures, mod_trap, 'o-', color='orange', label='Trapping Reduction Factor', markersize=4)
    plt.plot(temperatures, mod_combined, 'o-', color='red', label='Combined Factor', markersize=4)
    
    # Fill regions
    plt.fill_between(temperatures, 1, mod_gb, alpha=0.2, color='green', label='_GB region')
    plt.fill_between(temperatures, mod_trap, 1, alpha=0.2, color='orange', label='_Trap region')
    
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Modification Factor (D_eff/D_lattice)', fontsize=12)
    plt.title('Competing Effects: GB Enhancement vs Trapping Reduction', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    fname = f'{RESULTS_DIR}/plot19_competing_effects_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'temperatures': temperatures, 'mod_gb': mod_gb, 'mod_trap': mod_trap, 'mod_combined': mod_combined}


def plot_2D_regime_map_level4(K_s=DEFAULT_K_S, thickness=DEFAULT_THICKNESS, P_up=DEFAULT_P_UP, 
                               temperature=1073, n_grain=15, n_trap=15):
    """
    Plot 20: 2D Regime Map - Grain Size vs Trap Density.
    Shows modification factor as a heatmap.
    Level 1 annotation: Entire map would be mod_factor = 1 for Level 1.
    """
    print_test_header("Plot 20: 2D Regime Map (Grain Size vs Trap Density)")
    
    grain_sizes = np.logspace(-8, -4, n_grain)  # 10 nm to 100 μm
    trap_densities = np.logspace(20, 26, n_trap)  # 1e20 to 1e26 m^-3
    D_T = get_D_lattice(temperature)
    
    # Calculate modification factor for each combination
    mod_factor_map = np.zeros((n_trap, n_grain))
    
    for i, td in enumerate(trap_densities):
        for j, gs in enumerate(grain_sizes):
            params = get_microstructure_params('solution_annealed', extra={
                'grain_size': gs,
                'trap_density': td
            })
            result = calculate_defective_metal_flux(
                D_lattice=D_T, K_s=K_s, thickness=thickness,
                P_up=P_up, P_down=0.0, temperature=temperature,
                microstructure_params=params, mode='both'
            )
            mod_factor_map[i, j] = result['modification_factor']
    
    # Plot
    plt.figure(figsize=(12, 9))
    
    # Create meshgrid for pcolormesh
    GS, TD = np.meshgrid(grain_sizes * 1e6, trap_densities)  # Convert grain size to μm
    
    # Heatmap
    pcm = plt.pcolormesh(GS, TD, mod_factor_map, cmap='RdYlBu_r', shading='auto')
    plt.colorbar(pcm, label='Modification Factor (D_eff/D_lattice)')
    
    # Contour line at mod_factor = 1 (Level 1 equivalent)
    CS = plt.contour(GS, TD, mod_factor_map, levels=[1.0], colors='black', linewidths=2)
    plt.clabel(CS, inline=True, fontsize=10, fmt='Level 1 = 1.0')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Grain Size (μm)', fontsize=12)
    plt.ylabel('Trap Density (m⁻³)', fontsize=12)
    plt.title('2D Regime Map: Modification Factor\n(Black contour = Level 1 baseline)', fontsize=14)
    
    # Annotation for Level 1
    plt.annotate('Level 1: Entire map = 1.0\n(No microstructure effects)', 
                 xy=(0.02, 0.98), xycoords='axes fraction',
                 fontsize=10, ha='left', va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    fname = f'{RESULTS_DIR}/plot20_2D_regime_map_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png'
    plt.savefig(fname, dpi=300)
    print(f"  ✓ Plot saved: {fname}")
    plt.close()
    
    return {'grain_sizes': grain_sizes, 'trap_densities': trap_densities, 'mod_factor_map': mod_factor_map}


# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

def run_all_level4_tests():
    """Run all Level 4 validation tests."""
    print("\n" + "="*80)
    print("RUNNING ALL LEVEL 4 TESTS")
    print("="*80)
    
    # Group 1: Basic flux tests
    test_basic_flux_level4()
    test_zero_downstream_pressure_level4()
    test_equal_pressures_level4()
    
    # Group 2: Microstructure effects
    test_trapping_effect_level4()
    test_grain_boundary_effect_level4()
    test_mixed_microstructure_level4()
    
    # Group 3: Physical dependence
    test_thickness_dependence_level4()
    test_pressure_dependence_level4()
    test_temperature_dependence_level4()
    
    # Group 4: Regime classification
    test_regime_classification_fields_level4()
    
    print("\nAll Level 4 validation tests completed.")


def run_all_level4_plots():
    """Run all 20 Level 4 comparison plots (Level 1 vs Level 4)."""
    print("\n" + "="*80)
    print("GENERATING ALL 20 LEVEL 1 vs LEVEL 4 COMPARISON PLOTS")
    print("="*80)
    
    # Plot 1-2: Basic flux vs parameter sweeps
    plot_flux_vs_pressure_level4()
    plot_flux_vs_temperature_level4()
    
    # Plot 3-4: Sensitivity to microstructure parameters
    plot_sensitivity_to_trap_density_level4()
    plot_sensitivity_to_grain_size_level4()
    
    # Plot 5: Flux vs thickness
    plot_flux_vs_thickness_level4()
    
    # Plot 6-7: Scaling verification
    plot_arrhenius_level4()
    plot_sqrt_P_scaling_level4()
    
    # Plot 8-10: D_eff and modification factor vs pressure
    plot_D_eff_vs_pressure_level4()
    plot_D_eff_vs_temperature_level4()
    plot_modification_factor_vs_pressure_level4()
    
    # Plot 11: Mode comparison (continuous)
    plot_mode_comparison_continuous_level4()
    
    # Plot 12: Modification factor vs temperature
    plot_modification_factor_vs_temperature_level4()
    
    # Plot 13: Permeability check (J×L vs P)
    plot_permeability_check_level4()
    
    # Plot 14-15: Trap occupancy
    plot_trap_occupancy_vs_pressure_level4()
    plot_trap_occupancy_vs_temperature_level4()
    
    # Plot 16: Mode dominance vs temperature
    plot_mode_dominance_vs_temperature_level4()
    
    # Plot 17-18: Profiles through thickness
    plot_D_eff_profile_through_thickness_level4()
    plot_concentration_profile_level4()
    
    # Plot 19: Competing effects visualization
    plot_competing_effects_level4()
    
    # Plot 20: 2D regime map
    plot_2D_regime_map_level4()
    
    print("\nAll 20 Level 1 vs Level 4 comparison plots generated.")


def run_all():
    """Run all tests and generate all plots."""
    run_all_level4_tests()
    run_all_level4_plots()
    print("\n" + "="*80)
    print("ALL LEVEL 4 TESTS AND PLOTS COMPLETED")
    print("="*80)


if __name__ == "__main__":
    run_all()