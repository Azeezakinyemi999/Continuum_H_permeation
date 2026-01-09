"""
Phase 4: Experimental Comparison
Compares model predictions with experimental data from literature
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
from data.experimental_data import get_experimental_data, convert_to_SI

def compare_with_experiment(material_name='Incoloy800', thickness=0.001):
    """
    Compare model predictions with experimental data.
    
    Parameters
    ----------
    material_name : str
        Material name from MATERIALS
    thickness : float
        Assumed thickness for flux calculations (m)
    
    Returns
    -------
    dict
        Comparison results
    """
    # Get experimental data
    # Handle material name formatting (Incoloy800 vs Incoloy 800)
    exp_material_name = 'Incoloy 800' if material_name == 'Incoloy800' else material_name
    exp_data_raw = get_experimental_data(exp_material_name, 'JAERI')
    exp_data = convert_to_SI(exp_data_raw)
    
    # Get model parameters
    material = MATERIALS[material_name]
    
    print(f"\n{'='*70}")
    print(f"EXPERIMENTAL COMPARISON: {material_name}")
    print(f"{'='*70}")
    print(f"Experimental data source: {exp_data_raw['source']}")
    print(f"Figure: {exp_data_raw['figure']}")
    print(f"Extraction method: {exp_data_raw['extraction_method']}")
    print(f"Model parameters from: {material['reference']}")
    print(f"{'='*70}\n")
    
    # Calculate model predictions at experimental temperatures
    model_permeabilities = []
    model_diffusivities = []
    model_solubilities = []
    
    for T_K in exp_data['temperatures_K']:
        D = get_diffusivity(T_K, material)
        K_s = get_solubility(T_K, material)
        P = get_permeability(T_K, material)
        
        model_permeabilities.append(P)
        model_diffusivities.append(D)
        model_solubilities.append(K_s)
    
    model_permeabilities = np.array(model_permeabilities)
    
    # Calculate error metrics
    relative_errors = (model_permeabilities - exp_data['permeabilities']) / exp_data['permeabilities'] * 100
    mean_error = np.mean(relative_errors)
    std_error = np.std(relative_errors)
    max_error = np.max(np.abs(relative_errors))
    
    # Calculate R² for log-log fit
    log_model = np.log10(model_permeabilities)
    log_exp = np.log10(exp_data['permeabilities'])
    correlation = np.corrcoef(log_model, log_exp)[0, 1]
    r_squared = correlation ** 2
    
    print(f"ERROR ANALYSIS:")
    print(f"  Mean relative error: {mean_error:+.1f}%")
    print(f"  Std deviation: {std_error:.1f}%")
    print(f"  Max absolute error: {max_error:.1f}%")
    print(f"  R² (log-log): {r_squared:.4f}")
    
    results = {
        'experimental': exp_data,
        'model': {
            'temperatures_K': exp_data['temperatures_K'],
            'permeabilities': model_permeabilities,
            'diffusivities': np.array(model_diffusivities),
            'solubilities': np.array(model_solubilities)
        },
        'errors': {
            'relative_errors': relative_errors,
            'mean_error': mean_error,
            'std_error': std_error,
            'max_error': max_error,
            'r_squared': r_squared
        },
        'material': material_name,
        'source': exp_data_raw
    }
    
    return results

def plot_comparison(results, save_figure=False):
    """
    Create comprehensive comparison plots.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Model vs Experimental Comparison - {results['material']}", 
                 fontsize=14, fontweight='bold')
    
    # Common x-axis data
    x_data = 1000.0 / results['experimental']['temperatures_K']
    
    # Plot 1: Permeability comparison (Arrhenius plot)
    ax1 = axes[0, 0]
    ax1.semilogy(x_data, results['experimental']['permeabilities'], 
                 'o', color='black', markersize=8, label='Experimental (JAERI)', markerfacecolor='none')
    ax1.semilogy(x_data, results['model']['permeabilities'], 
                 's-', color='red', markersize=6, linewidth=1, label='Model', alpha=0.7)
    
    ax1.set_xlabel('1000/T (K⁻¹)')
    ax1.set_ylabel('Permeability (mol/m/s/Pa^0.5)')
    ax1.set_title('Permeability: Model vs Experimental')
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)
    
    # Add temperature scale on top
    ax1_top = ax1.twiny()
    temp_ticks = np.array([600, 700, 800, 900, 1000, 1100])
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(1000/(temp_ticks+273.15))
    ax1_top.set_xticklabels([f"{t}" for t in temp_ticks])
    ax1_top.set_xlabel('Temperature (°C)')
    
    # Plot 2: Relative error vs temperature
    ax2 = axes[0, 1]
    temperatures_C = results['experimental']['temperatures_K'] - 273.15
    ax2.plot(temperatures_C, results['errors']['relative_errors'], 
             'o-', color='blue', markersize=6)
    ax2.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax2.axhline(y=results['errors']['mean_error'], color='red', 
                linestyle='--', alpha=0.5, label=f"Mean: {results['errors']['mean_error']:.1f}%")
    
    ax2.set_xlabel('Temperature (°C)')
    ax2.set_ylabel('Relative Error (%)')
    ax2.set_title('Model Error vs Temperature')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Parity plot (Model vs Experimental)
    ax3 = axes[1, 0]
    ax3.loglog(results['experimental']['permeabilities'], 
               results['model']['permeabilities'], 
               'o', color='green', markersize=8)
    
    # Add 1:1 line
    p_min = min(results['experimental']['permeabilities'].min(), 
                results['model']['permeabilities'].min())
    p_max = max(results['experimental']['permeabilities'].max(), 
                results['model']['permeabilities'].max())
    ax3.loglog([p_min, p_max], [p_min, p_max], 'k--', alpha=0.5, label='1:1 line')
    
    # Add ±50% error bands
    ax3.fill_between([p_min, p_max], [p_min*0.5, p_max*0.5], [p_min*1.5, p_max*1.5], 
                     alpha=0.2, color='gray', label='±50% error')
    
    ax3.set_xlabel('Experimental P (mol/m/s/Pa^0.5)')
    ax3.set_ylabel('Model P (mol/m/s/Pa^0.5)')
    ax3.set_title(f"Parity Plot (R² = {results['errors']['r_squared']:.4f})")
    ax3.legend()
    ax3.grid(True, which="both", alpha=0.3)
    
    # Plot 4: Error distribution histogram
    ax4 = axes[1, 1]
    ax4.hist(results['errors']['relative_errors'], bins=10, 
             color='purple', alpha=0.7, edgecolor='black')
    ax4.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax4.axvline(x=results['errors']['mean_error'], color='red', 
                linestyle='--', linewidth=2, label=f"Mean: {results['errors']['mean_error']:.1f}%")
    
    ax4.set_xlabel('Relative Error (%)')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Error Distribution')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_figure:
        filename = f"experimental_comparison_{results['material']}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved as: {filename}")
    
    plt.show()

def save_comparison_results(results, filename=None):
    """
    Save comparison results to text file with proper citation.
    """
    if filename is None:
        filename = f"experimental_comparison_{results['material']}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(filename, 'w') as f:
        f.write("="*80 + "\n")
        f.write("EXPERIMENTAL COMPARISON RESULTS\n")
        f.write("="*80 + "\n\n")
        
        # Important citation
        f.write("DATA CITATION:\n")
        f.write(f"Data for {results['source']['material']} [{results['source']['source']}] were extracted\n")
        f.write(f"from {results['source']['figure']} using {results['source']['extraction_method']}.\n\n")
        
        f.write(f"Date/Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Material: {results['material']}\n\n")
        
        # Error summary
        f.write("-"*80 + "\n")
        f.write("ERROR ANALYSIS SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write(f"Mean relative error: {results['errors']['mean_error']:+.2f}%\n")
        f.write(f"Standard deviation: {results['errors']['std_error']:.2f}%\n")
        f.write(f"Maximum absolute error: {results['errors']['max_error']:.2f}%\n")
        f.write(f"R² (log-log correlation): {results['errors']['r_squared']:.4f}\n\n")
        
        # Detailed comparison table
        f.write("-"*80 + "\n")
        f.write("DETAILED COMPARISON\n")
        f.write("-"*80 + "\n")
        f.write(f"{'T(°C)':>8} {'T(K)':>8} {'P_exp':>14} {'P_model':>14} {'Error(%)':>12}\n")
        f.write("-"*80 + "\n")
        
        for i in range(len(results['experimental']['temperatures_K'])):
            T_K = results['experimental']['temperatures_K'][i]
            T_C = T_K - 273.15
            P_exp = results['experimental']['permeabilities'][i]
            P_model = results['model']['permeabilities'][i]
            error = results['errors']['relative_errors'][i]
            
            f.write(f"{T_C:8.1f} {T_K:8.1f} {P_exp:14.3e} {P_model:14.3e} {error:12.1f}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"Results saved to: {filename}")

def main():
    """
    Run experimental comparison analysis.
    """
    print("\n" + "="*70)
    print("PHASE 4: EXPERIMENTAL COMPARISON")
    print("="*70)
    
    # Perform comparison
    results = compare_with_experiment('Incoloy800')
    
    # Create plots
    plot_comparison(results, save_figure=True)
    
    # Save results
    save_comparison_results(results)
    
    print("\n" + "="*70)
    print("COMPARISON COMPLETE")
    print("="*70)
    
    return results

if __name__ == "__main__":
    results = main()