"""
Phase 3: Temperature Dependence Study
Validates Arrhenius behavior and extracts activation energies
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys
import os
from datetime import datetime

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.utils import get_diffusivity, get_solubility, get_permeability
from data.material_data import MATERIALS

def temperature_sweep(material_name='Incoloy800', T_min_C=600, T_max_C=1000, 
                     n_points=9, P_upstream=1.0, P_downstream=0, thickness=0.001):
    """
    Perform temperature sweep and calculate transport properties.
    
    Parameters
    ----------
    material_name : str
        Material key from MATERIALS dictionary
    T_min_C : float
        Minimum temperature in Celsius
    T_max_C : float
        Maximum temperature in Celsius
    n_points : int
        Number of temperature points
    P_upstream : float
        Upstream pressure in Pa
    P_downstream : float
        Downstream pressure in Pa
    thickness : float
        Material thickness in meters
    
    Returns
    -------
    dict
        Dictionary containing temperatures, properties, and fluxes
    """
    # Get material
    material = MATERIALS[material_name]
    
    # Check temperature range
    if T_min_C < material['temp_range'][0] or T_max_C > material['temp_range'][1]:
        print(f"⚠️  Warning: Temperature range ({T_min_C}-{T_max_C}°C) exceeds")
        print(f"   material data range {material['temp_range']}°C")
    
    print(f"\n{'='*60}")
    print(f"Temperature Sweep for {material_name}")
    print(f"{'='*60}")
    print(f"Temperature range: {T_min_C} - {T_max_C}°C")
    print(f"Pressure: {P_upstream} Pa (upstream), {P_downstream} Pa (downstream)")
    print(f"Thickness: {thickness*1000:.2f} mm")
    print(f"Reference: {material['reference']}")
    print(f"{'='*60}\n")
    
    # Generate temperature points
    temperatures_C = np.linspace(T_min_C, T_max_C, n_points)
    temperatures_K = temperatures_C + 273.15
    
    # Calculate properties at each temperature
    diffusivities = []
    solubilities = []
    permeabilities = []
    fluxes = []
    
    print(f"{'T (°C)':>8} {'T (K)':>8} {'D (m²/s)':>12} {'K_s (mol/m³/Pa^0.5)':>18} {'P (mol/m/s/Pa^0.5)':>18}")
    print("-" * 80)
    
    for T_C, T_K in zip(temperatures_C, temperatures_K):
        D = get_diffusivity(T_K, material)
        K_s = get_solubility(T_K, material)
        P = get_permeability(T_K, material)
        
        # Calculate flux
        result = calculate_simple_metal_flux(D, K_s, thickness, P_upstream, P_downstream)
        flux = result['flux']
        
        diffusivities.append(D)
        solubilities.append(K_s)
        permeabilities.append(P)
        fluxes.append(flux)
        
        print(f"{T_C:8.1f} {T_K:8.1f} {D:12.3e} {K_s:18.3e} {P:18.3e}")
    
    results = {
        'temperatures_C': temperatures_C,
        'temperatures_K': temperatures_K,
        'diffusivities': np.array(diffusivities),
        'solubilities': np.array(solubilities),
        'permeabilities': np.array(permeabilities),
        'fluxes': np.array(fluxes),
        'material': material_name,
        'material_data': material,
        'P_upstream': P_upstream,
        'P_downstream': P_downstream,
        'thickness': thickness
    }
    
    return results

def arrhenius_analysis(temperatures_K, property_values, property_name="Property"):
    """
    Perform Arrhenius analysis to extract activation energy.
    
    Parameters
    ----------
    temperatures_K : array-like
        Temperatures in Kelvin
    property_values : array-like
        Property values (D, K_s, or P)
    property_name : str
        Name of property for labeling
    
    Returns
    -------
    dict
        Dictionary containing activation energy and pre-exponential factor
    """
    # Arrhenius plot: ln(property) vs 1000/T
    x_data = 1000.0 / temperatures_K  # 1000/T for better scaling
    y_data = np.log(property_values)
    
    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
    r_squared = r_value**2
    
    # Extract activation energy
    # slope = -E_a / (R * 1000) since we used 1000/T
    R = 8.314  # J/mol/K
    E_a = -slope * R * 1000  # J/mol
    E_a_kJ = E_a / 1000  # kJ/mol
    
    # Pre-exponential factor
    pre_exp = np.exp(intercept)
    
    results = {
        'property_name': property_name,
        'E_a': E_a,
        'E_a_kJ': E_a_kJ,
        'pre_exp': pre_exp,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'std_error': std_err
    }
    
    print(f"\n{property_name} Arrhenius Analysis:")
    print(f"  Activation energy: {E_a_kJ:.2f} kJ/mol")
    print(f"  Pre-exponential: {pre_exp:.3e}")
    print(f"  R²: {r_squared:.6f}")
    
    return results

def plot_arrhenius(results, analyses, save_figure=False):
    """
    Create comprehensive Arrhenius plots.
    
    Parameters
    ----------
    results : dict
        Results from temperature_sweep
    analyses : dict
        Dictionary of Arrhenius analyses
    save_figure : bool
        Whether to save the figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Temperature Dependence - {results['material']} (Arrhenius Analysis)", 
                 fontsize=14, fontweight='bold')
    
    # Common x-axis data
    x_data = 1000.0 / results['temperatures_K']
    
    # Plot 1: Diffusivity
    ax1 = axes[0, 0]
    y_D = results['diffusivities']
    ax1.semilogy(x_data, y_D, 'o', color='red', markersize=8, label='Data')
    
    # Add fitted line
    analysis_D = analyses['diffusivity']
    y_fit_D = analysis_D['pre_exp'] * np.exp(analysis_D['slope'] * x_data)
    ax1.semilogy(x_data, y_fit_D, '--', color='darkred', linewidth=2, label='Fit')
    
    ax1.set_xlabel('1000/T (K⁻¹)')
    ax1.set_ylabel('Diffusivity (m²/s)')
    ax1.set_title(f'Diffusivity: E_D = {analysis_D["E_a_kJ"]:.1f} kJ/mol')
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend()
    
    # Add temperature scale on top
    ax1_top = ax1.twiny()
    temp_ticks = np.array([600, 700, 800, 900, 1000])
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(1000/(temp_ticks+273.15))
    ax1_top.set_xticklabels([f"{t}" for t in temp_ticks])
    ax1_top.set_xlabel('Temperature (°C)')
    
    # Plot 2: Solubility
    ax2 = axes[0, 1]
    y_K = results['solubilities']
    ax2.semilogy(x_data, y_K, 's', color='blue', markersize=8, label='Data')
    
    analysis_K = analyses['solubility']
    y_fit_K = analysis_K['pre_exp'] * np.exp(analysis_K['slope'] * x_data)
    ax2.semilogy(x_data, y_fit_K, '--', color='darkblue', linewidth=2, label='Fit')
    
    ax2.set_xlabel('1000/T (K⁻¹)')
    ax2.set_ylabel('Solubility (mol/m³/Pa^0.5)')
    ax2.set_title(f'Solubility: ΔH_s = {analysis_K["E_a_kJ"]:.1f} kJ/mol')
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend()
    
    # Plot 3: Permeability
    ax3 = axes[1, 0]
    y_P = results['permeabilities']
    ax3.semilogy(x_data, y_P, '^', color='green', markersize=8, label='Data')
    
    analysis_P = analyses['permeability']
    y_fit_P = analysis_P['pre_exp'] * np.exp(analysis_P['slope'] * x_data)
    ax3.semilogy(x_data, y_fit_P, '--', color='darkgreen', linewidth=2, label='Fit')
    
    ax3.set_xlabel('1000/T (K⁻¹)')
    ax3.set_ylabel('Permeability (mol/m/s/Pa^0.5)')
    ax3.set_title(f'Permeability: E_P = {analysis_P["E_a_kJ"]:.1f} kJ/mol')
    ax3.grid(True, which="both", alpha=0.3)
    ax3.legend()
    
    # Plot 4: Flux at fixed pressure
    ax4 = axes[1, 1]
    y_flux = results['fluxes']
    ax4.semilogy(x_data, y_flux, 'd', color='purple', markersize=8, label=f'P={results["P_upstream"]} Pa')
    ax4.set_xlabel('1000/T (K⁻¹)')
    ax4.set_ylabel('Flux (mol/m²/s)')
    ax4.set_title(f'Flux at {results["P_upstream"]} Pa')
    ax4.grid(True, which="both", alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    
    if save_figure:
        filename = f"temperature_study_{results['material']}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved as: {filename}")
    
    plt.show()

def validate_against_input(results, analyses):
    """
    Compare extracted values with input material data.
    
    Parameters
    ----------
    results : dict
        Results from temperature_sweep
    analyses : dict
        Arrhenius analysis results
    """
    material = results['material_data']
    
    print("\n" + "="*60)
    print("VALIDATION: Extracted vs Input Parameters")
    print("="*60)
    
    # Diffusivity
    print("\nDiffusivity:")
    print(f"  Input E_D:      {material['E_D']/1000:.2f} kJ/mol")
    print(f"  Extracted E_D:  {analyses['diffusivity']['E_a_kJ']:.2f} kJ/mol")
    print(f"  Difference:     {abs(material['E_D']/1000 - analyses['diffusivity']['E_a_kJ']):.2f} kJ/mol")
    
    # Solubility
    print("\nSolubility:")
    print(f"  Input ΔH_s:     {material['H_s']/1000:.2f} kJ/mol")
    print(f"  Extracted ΔH_s: {analyses['solubility']['E_a_kJ']:.2f} kJ/mol")
    print(f"  Difference:     {abs(material['H_s']/1000 - analyses['solubility']['E_a_kJ']):.2f} kJ/mol")
    
    # Permeability (should equal E_D + H_s)
    expected_E_P = (material['E_D'] + material['H_s'])/1000
    print("\nPermeability:")
    print(f"  Expected E_P (E_D + ΔH_s): {expected_E_P:.2f} kJ/mol")
    print(f"  Extracted E_P:              {analyses['permeability']['E_a_kJ']:.2f} kJ/mol")
    print(f"  Difference:                 {abs(expected_E_P - analyses['permeability']['E_a_kJ']):.2f} kJ/mol")
    
    print("="*60)

def save_results_to_file(results, analyses, filename=None):
    """
    Save comprehensive results to a text file.
    
    Parameters
    ----------
    results : dict
        Results from temperature_sweep
    analyses : dict
        Arrhenius analysis results
    filename : str, optional
        Output filename
    """
    if filename is None:
        filename = f"temperature_study_{results['material']}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    material = results['material_data']
    
    with open(filename, 'w') as f:
        # Header
        f.write("="*80 + "\n")
        f.write("TEMPERATURE DEPENDENCE STUDY RESULTS\n")
        f.write("="*80 + "\n\n")
        
        # Study parameters
        f.write(f"Date/Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Material: {results['material']}\n")
        f.write(f"Reference: {material['reference']}\n")
        f.write(f"Temperature range: {results['temperatures_C'][0]:.1f} - {results['temperatures_C'][-1]:.1f}°C\n")
        f.write(f"Pressure: {results['P_upstream']} Pa (upstream), {results['P_downstream']} Pa (downstream)\n")
        f.write(f"Thickness: {results['thickness']*1000:.2f} mm\n\n")
        
        # Input material properties
        f.write("-"*80 + "\n")
        f.write("INPUT MATERIAL PROPERTIES\n")
        f.write("-"*80 + "\n")
        f.write(f"D_0 = {material['D_0']:.3e} m²/s\n")
        f.write(f"E_D = {material['E_D']/1000:.2f} kJ/mol\n")
        f.write(f"K_s0 = {material['K_s0']:.3e} mol/m³/Pa^0.5\n")
        f.write(f"ΔH_s = {material['H_s']/1000:.2f} kJ/mol\n\n")
        
        # Calculated values table
        f.write("-"*80 + "\n")
        f.write("CALCULATED VALUES\n")
        f.write("-"*80 + "\n")
        f.write(f"{'T(°C)':>8} {'T(K)':>8} {'1000/T':>10} {'D(m²/s)':>12} {'K_s(mol/m³/Pa^0.5)':>18} {'P(mol/m/s/Pa^0.5)':>18} {'Flux(mol/m²/s)':>15}\n")
        f.write("-"*80 + "\n")
        
        for i in range(len(results['temperatures_C'])):
            f.write(f"{results['temperatures_C'][i]:8.1f} ")
            f.write(f"{results['temperatures_K'][i]:8.1f} ")
            f.write(f"{1000/results['temperatures_K'][i]:10.4f} ")
            f.write(f"{results['diffusivities'][i]:12.3e} ")
            f.write(f"{results['solubilities'][i]:18.3e} ")
            f.write(f"{results['permeabilities'][i]:18.3e} ")
            f.write(f"{results['fluxes'][i]:15.3e}\n")
        
        # Arrhenius analysis results
        f.write("\n" + "-"*80 + "\n")
        f.write("ARRHENIUS ANALYSIS RESULTS\n")
        f.write("-"*80 + "\n\n")
        
        # Diffusivity analysis
        f.write("DIFFUSIVITY:\n")
        f.write(f"  Extracted E_D = {analyses['diffusivity']['E_a_kJ']:.2f} kJ/mol\n")
        f.write(f"  Extracted D_0 = {analyses['diffusivity']['pre_exp']:.3e} m²/s\n")
        f.write(f"  R² = {analyses['diffusivity']['r_squared']:.6f}\n")
        f.write(f"  Standard error = {analyses['diffusivity']['std_error']:.3e}\n\n")
        
        # Solubility analysis
        f.write("SOLUBILITY:\n")
        f.write(f"  Extracted ΔH_s = {analyses['solubility']['E_a_kJ']:.2f} kJ/mol\n")
        f.write(f"  Extracted K_s0 = {analyses['solubility']['pre_exp']:.3e} mol/m³/Pa^0.5\n")
        f.write(f"  R² = {analyses['solubility']['r_squared']:.6f}\n")
        f.write(f"  Standard error = {analyses['solubility']['std_error']:.3e}\n\n")
        
        # Permeability analysis
        f.write("PERMEABILITY:\n")
        f.write(f"  Extracted E_P = {analyses['permeability']['E_a_kJ']:.2f} kJ/mol\n")
        f.write(f"  Extracted P_0 = {analyses['permeability']['pre_exp']:.3e} mol/m/s/Pa^0.5\n")
        f.write(f"  R² = {analyses['permeability']['r_squared']:.6f}\n")
        f.write(f"  Standard error = {analyses['permeability']['std_error']:.3e}\n\n")
        
        # Validation
        f.write("-"*80 + "\n")
        f.write("VALIDATION: Input vs Extracted Parameters\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Parameter':<20} {'Input':<15} {'Extracted':<15} {'Difference':<15} {'Error (%)':<10}\n")
        f.write("-"*80 + "\n")
        
        # Diffusivity validation
        input_ED = material['E_D']/1000
        extracted_ED = analyses['diffusivity']['E_a_kJ']
        diff_ED = abs(input_ED - extracted_ED)
        error_ED = 100 * diff_ED / abs(input_ED) if input_ED != 0 else 0
        f.write(f"{'E_D (kJ/mol)':<20} {input_ED:<15.2f} {extracted_ED:<15.2f} {diff_ED:<15.2f} {error_ED:<10.2f}\n")
        
        # Solubility validation
        input_Hs = material['H_s']/1000
        extracted_Hs = analyses['solubility']['E_a_kJ']
        diff_Hs = abs(input_Hs - extracted_Hs)
        error_Hs = 100 * diff_Hs / abs(input_Hs) if input_Hs != 0 else 0
        f.write(f"{'ΔH_s (kJ/mol)':<20} {input_Hs:<15.2f} {extracted_Hs:<15.2f} {diff_Hs:<15.2f} {error_Hs:<10.2f}\n")
        
        # Permeability validation
        expected_EP = (material['E_D'] + material['H_s'])/1000
        extracted_EP = analyses['permeability']['E_a_kJ']
        diff_EP = abs(expected_EP - extracted_EP)
        error_EP = 100 * diff_EP / abs(expected_EP) if expected_EP != 0 else 0
        f.write(f"{'E_P (kJ/mol)':<20} {expected_EP:<15.2f} {extracted_EP:<15.2f} {diff_EP:<15.2f} {error_EP:<10.2f}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"\nResults saved to: {filename}")

def main():
    """
    Run the complete temperature dependence study.
    """
    print("\n" + "="*60)
    print("PHASE 3: TEMPERATURE DEPENDENCE STUDY")
    print("="*60)
    
    # Perform temperature sweep
    results = temperature_sweep(
        material_name='Incoloy800',
        T_min_C=600,
        T_max_C=1000,
        n_points=9,
        P_upstream=1.0,  # Pa
        thickness=0.001   # 1 mm
    )
    
    # Perform Arrhenius analysis for each property
    analyses = {
        'diffusivity': arrhenius_analysis(
            results['temperatures_K'], 
            results['diffusivities'], 
            "Diffusivity"
        ),
        'solubility': arrhenius_analysis(
            results['temperatures_K'], 
            results['solubilities'], 
            "Solubility"
        ),
        'permeability': arrhenius_analysis(
            results['temperatures_K'], 
            results['permeabilities'], 
            "Permeability"
        )
    }
    
    # Validate against input
    validate_against_input(results, analyses)

    # Save results to text file
    save_results_to_file(results, analyses)
    
    # Plot results
    plot_arrhenius(results, analyses, save_figure=True)
    
    print("\n" + "="*60)
    print("STUDY COMPLETE")
    print("="*60)
    print(f"✅ Temperature sweep completed: {len(results['temperatures_K'])} points")
    print(f"✅ Activation energies extracted")
    print(f"✅ Results validated against input")
    print("="*60 + "\n")
    
    return results, analyses

if __name__ == "__main__":
    results, analyses = main()