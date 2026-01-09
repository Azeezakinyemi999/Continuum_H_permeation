"""
Phase 2: Pressure Dependence Study
Validates Sieverts' law by analyzing flux vs pressure relationship
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys
import os
from datetime import datetime

# Add parent directory to path to import your modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.utils import get_diffusivity, get_solubility
from data.material_data import MATERIALS

def pressure_sweep(material_name='Incoloy800', temperature_C=800, thickness=0.001, 
                   P_min=0.01, P_max=100, n_points=20, P_down=0):
    """
    Perform pressure sweep and calculate flux at each pressure.
    
    Parameters
    ----------
    material_name : str
        Material key from MATERIALS dictionary
    temperature_C : float
        Temperature in Celsius
    thickness : float
        Material thickness in meters
    P_min : float
        Minimum upstream pressure in Pa
    P_max : float
        Maximum upstream pressure in Pa
    n_points : int
        Number of pressure points (log-spaced)
    P_down : float
        Downstream pressure in Pa (default 0)
    
    Returns
    -------
    dict
        Dictionary containing pressures, fluxes, and calculation parameters
    """
    # Get material properties
    material = MATERIALS[material_name]
    T_kelvin = temperature_C + 273.15
    
    # Calculate material properties at temperature
    D = get_diffusivity(T_kelvin, material)
    K_s = get_solubility(T_kelvin, material)
    permeability = D * K_s
    
    print(f"\n{'='*60}")
    print(f"Pressure Sweep for {material_name}")
    print(f"{'='*60}")
    print(f"Temperature: {temperature_C}°C ({T_kelvin:.1f} K)")
    print(f"Thickness: {thickness*1000:.2f} mm")
    print(f"D = {D:.2e} m²/s")
    print(f"K_s = {K_s:.2e} mol/m³/Pa^0.5")
    print(f"Permeability = {permeability:.2e} mol/m/s/Pa^0.5")
    print(f"{'='*60}\n")
    
    # Generate pressure points (log-spaced)
    pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
    
    # Calculate flux at each pressure
    fluxes = []
    concentrations_up = []
    concentrations_down = []
    
    for P in pressures:
        result = calculate_simple_metal_flux(D, K_s, thickness, P, P_down)
        fluxes.append(result['flux'])
        concentrations_up.append(result['C_up'])
        concentrations_down.append(result['C_down'])
    
    # Store results
    results = {
        'pressures': pressures,
        'fluxes': np.array(fluxes),
        'C_up': np.array(concentrations_up),
        'C_down': np.array(concentrations_down),
        'material': material_name,
        'temperature_C': temperature_C,
        'thickness': thickness,
        'D': D,
        'K_s': K_s,
        'permeability': permeability
    }
    
    return results

def validate_sieverts_law(pressures, fluxes, tolerance = 0.01):
    """
    Calculate slope of log(flux) vs log(pressure) and verify Sieverts' law.
    
    Parameters
    ----------
    pressures : array-like
        Pressure values in Pa
    fluxes : array-like
        Flux values in mol/m²/s
    
    Returns
    -------
    dict
        Dictionary containing slope, R², and validation status
    """
    # Calculate logarithms
    log_p = np.log10(pressures)
    log_f = np.log10(fluxes)
    
    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_p, log_f)
    r_squared = r_value**2
    
    # Check if slope is close to 0.5 (Sieverts' law)
    expected_slope = 0.5
    is_valid = abs(slope - expected_slope) < tolerance
    
    validation_results = {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'std_error': std_err,
        'expected_slope': expected_slope,
        'is_valid': is_valid,
        'deviation': abs(slope - expected_slope)
    }
    
    print(f"\nSieverts' Law Validation:")
    print(f"{'='*40}")
    print(f"Calculated slope: {slope:.4f}")
    print(f"Expected slope: {expected_slope}")
    print(f"Deviation: {validation_results['deviation']:.4f}")
    print(f"R² value: {r_squared:.6f}")
    print(f"Validation: {'✅ PASSED' if is_valid else '❌ FAILED'}")
    print(f"{'='*40}\n")
    
    return validation_results

def plot_results(results, validation, save_figure=False):
    """
    Create comprehensive plots for pressure study.
    
    Parameters
    ----------
    results : dict
        Results from pressure_sweep
    validation : dict
        Results from validate_sieverts_law
    save_figure : bool
        Whether to save the figure to file
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Pressure Dependence Study - {results['material']} at {results['temperature_C']}°C", 
                 fontsize=14, fontweight='bold')
    
    # Plot 1: Flux vs Pressure (log-log)
    ax1 = axes[0, 0]
    ax1.loglog(results['pressures'], results['fluxes'], 'o-', 
               color='blue', markersize=6, linewidth=2, label='Calculated')
    
    # Add fitted line ::: flux is proportional to the sqrt of p 
    log_p = np.log10(results['pressures'])
    fitted_log_flux = validation['slope'] * log_p + validation['intercept']
    fitted_flux = 10**fitted_log_flux
    ax1.loglog(results['pressures'], fitted_flux, '--', 
               color='red', linewidth=1, alpha=0.7, label='Fitted')
    
    ax1.set_xlabel('Upstream Pressure (Pa)')
    ax1.set_ylabel('Flux (mol/m²/s)')
    ax1.set_title('Flux vs Pressure (log-log scale)')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.legend()
    
    # Add slope annotation
    ax1.text(0.05, 0.95, f'Slope = {validation["slope"]:.3f}\nExpected = 0.500\nR² = {validation["r_squared"]:.4f}',
             transform=ax1.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Plot 2: Linear scale
    ax2 = axes[0, 1]
    ax2.plot(results['pressures'], results['fluxes'], 'o-', 
             color='green', markersize=6, linewidth=2)
    ax2.set_xlabel('Upstream Pressure (Pa)')
    ax2.set_ylabel('Flux (mol/m²/s)')
    ax2.set_title('Flux vs Pressure (linear scale)')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Concentration profiles
    ax3 = axes[1, 0]
    ax3.plot(results['pressures'], results['C_up']*1000, 'o-', 
             color='red', markersize=6, linewidth=2, label='Upstream')
    ax3.plot(results['pressures'], results['C_down']*1000, 's-', 
             color='blue', markersize=6, linewidth=2, label='Downstream')
    ax3.set_xlabel('Upstream Pressure (Pa)')
    ax3.set_ylabel('Concentration (mmol/m³)')
    ax3.set_title('Surface Concentrations')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Plot 4: Residuals
    # ax4 = axes[1, 1]
    # log_flux_actual = np.log10(results['fluxes'])
    # log_flux_fitted = validation['slope'] * np.log10(results['pressures']) + validation['intercept']
    # residuals = log_flux_actual - log_flux_fitted
    
    # ax4.plot(results['pressures'], residuals, 'o', color='purple', markersize=6)
    # ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    # ax4.set_xlabel('Pressure (Pa)')
    # ax4.set_ylabel('Log(Flux) Residuals')
    # ax4.set_title('Fit Residuals')
    # ax4.grid(True, alpha=0.3)
    # ax4.set_xscale('log')
    
    # plt.tight_layout()
    
    if save_figure:
        filename = f"pressure_study_{results['material']}_{results['temperature_C']}C_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filename}")
    
    plt.show()

def save_results_to_file(results, validation, filename=None):
    """
    Save results to a text file for future reference.
    
    Parameters
    ----------
    results : dict
        Results from pressure_sweep
    validation : dict
        Results from validate_sieverts_law
    filename : str, optional
        Output filename
    """
    if filename is None:
        filename = f"pressure_study_{results['material']}_{results['temperature_C']}C_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(filename, 'w') as f:
        f.write(f"Pressure Dependence Study Results\n")
        f.write(f"{'='*60}\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Material: {results['material']}\n")
        f.write(f"Temperature: {results['temperature_C']}°C\n")
        f.write(f"Thickness: {results['thickness']*1000:.2f} mm\n")
        f.write(f"D: {results['D']:.2e} m²/s\n")
        f.write(f"K_s: {results['K_s']:.2e} mol/m³/Pa^0.5\n")
        f.write(f"Permeability: {results['permeability']:.2e} mol/m/s/Pa^0.5\n")
        f.write(f"\nSieverts' Law Validation:\n")
        f.write(f"Slope: {validation['slope']:.4f} (expected: 0.500)\n")
        f.write(f"R²: {validation['r_squared']:.6f}\n")
        f.write(f"Validation: {'PASSED' if validation['is_valid'] else 'FAILED'}\n")
        f.write(f"\nPressure (Pa)\tFlux (mol/m²/s)\tC_up (mol/m³)\n")
        f.write(f"{'-'*60}\n")
        
        for i, P in enumerate(results['pressures']):
            f.write(f"{P:.2e}\t{results['fluxes'][i]:.2e}\t{results['C_up'][i]:.2e}\n")
    
    print(f"Results saved to: {filename}")

def main():
    """
    Run the complete pressure dependence study.
    """
    print("\n" + "="*60)
    print("PHASE 2: PRESSURE DEPENDENCE STUDY")
    print("="*60)
    
    # Perform pressure sweep
    results = pressure_sweep(
        material_name='Incoloy800',
        temperature_C=800,
        thickness=0.001,  # 1 mm
        P_min=0.01,       # Pa
        P_max=100,        # Pa
        n_points=20
    )
    
    # Validate Sieverts' law
    validation = validate_sieverts_law(results['pressures'], results['fluxes'])
    
    # Plot results
    plot_results(results, validation, save_figure=True)
    
    # Save results to file
    save_results_to_file(results, validation)
    
    # Print summary
    print("\n" + "="*60)
    print("STUDY COMPLETE")
    print("="*60)
    print(f"✅ Pressure sweep completed: {len(results['pressures'])} points")
    print(f"✅ Sieverts' law validation: {'PASSED' if validation['is_valid'] else 'FAILED'}")
    print(f"✅ Results plotted and saved")
    print("="*60 + "\n")
    
    return results, validation

if __name__ == "__main__":
    results, validation = main()