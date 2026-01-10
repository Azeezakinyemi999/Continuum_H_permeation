"""
Level 3 Experimental Comparison

Compares Level 1, Level 2, and Level 3 predictions against experimental data.
Key goal: Determine if defect model can explain discrepancies between 
perfect oxide predictions and measured permeation rates.

Based on:
- Strehlow & Savage (1974): Parallel path model
- Zhang et al. (2018): PRF values 10-3828
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize
from datetime import datetime
import os
import sys

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import models
from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.interface_solver import calculate_oxide_metal_system
from calculations.parallel_oxide_defect_paths import (
    calculate_parallel_path_flux,
    calculate_PRF
)

# Import data
from data.material_data import MATERIALS
from data.oxide_properties import OXIDE_PROPERTIES
from data.experimental_data import get_experimental_data, convert_to_SI
from data.oxide_defect_parameters import DEFECT_CONFIGURATIONS

# Constants
R = 8.314  # J/mol/K

# Output directory
OUTPUT_DIR = "analysis/results/level3_comparison"


def ensure_output_dir():
    """Create output directory if needed."""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)


def get_properties_at_temperature(T_K, oxide_name='Cr2O3', metal_name='Incoloy800'):
    """Calculate temperature-dependent properties."""
    oxide_data = OXIDE_PROPERTIES[oxide_name]
    metal_data = MATERIALS[metal_name]
    
    # Oxide properties
    D_ox = oxide_data['D_ox_0'] * np.exp(-oxide_data['E_D_ox'] / (R * T_K))
    K_ox = oxide_data['K_ox_0'] * np.exp(-oxide_data['H_sol_ox'] / (R * T_K))
    
    oxide_props = {
        'D_ox': D_ox,
        'K_ox': K_ox,
        'thickness': oxide_data['thickness']
    }
    
    # Metal properties
    D_metal = metal_data['D_0'] * np.exp(-metal_data['E_D'] / (R * T_K))
    K_s_metal = metal_data['K_s0'] * np.exp(-metal_data['H_s'] / (R * T_K))
    
    metal_props = {
        'D_metal': D_metal,
        'K_s_metal': K_s_metal,
        'thickness': 1e-3  # 1 mm standard
    }
    
    return oxide_props, metal_props


def calculate_model_permeability(T_K, P_test, model_level, defect_params=None,
                                  oxide_name='Cr2O3', metal_name='Incoloy800'):
    """
    Calculate permeability at given temperature for specified model level.
    
    Permeability P = J * L / sqrt(P_up) for Sieverts-type behavior
    For oxide systems, we report effective permeability.
    """
    oxide_props, metal_props = get_properties_at_temperature(T_K, oxide_name, metal_name)
    L = metal_props['thickness']
    
    if model_level == 1:
        # Level 1: Bare metal
        result = calculate_simple_metal_flux(
            metal_props['D_metal'],
            metal_props['K_s_metal'],
            L, P_test, 0
        )
        flux = result['flux']
        # Permeability = D * K_s
        permeability = metal_props['D_metal'] * metal_props['K_s_metal']
        
    elif model_level == 2:
        # Level 2: Perfect oxide + metal
        result = calculate_oxide_metal_system(P_test, 0, oxide_props, metal_props)
        flux = result['flux']
        # Effective permeability from flux
        permeability = flux * L / np.sqrt(P_test) if P_test > 0 else 0
        
    elif model_level == 3:
        # Level 3: Defective oxide + metal
        if defect_params is None:
            defect_params = {'area_fraction': 0.01, 'type': 'pinhole'}
        result = calculate_parallel_path_flux(
            P_test, 0, oxide_props, metal_props, defect_params
        )
        flux = result['flux_total']
        permeability = flux * L / np.sqrt(P_test) if P_test > 0 else 0
        
    else:
        raise ValueError(f"Unknown model level: {model_level}")
    
    return permeability, flux


def fit_defect_fraction(exp_temperatures, exp_permeabilities, P_test=1000,
                        defect_type='pinhole', oxide_name='Cr2O3'):
    """
    Find the defect fraction that best fits experimental data.
    
    Parameters:
    -----------
    exp_temperatures : array
        Experimental temperatures (K)
    exp_permeabilities : array
        Experimental permeabilities (mol/m/s/Pa^0.5)
    P_test : float
        Test pressure for flux calculations (Pa)
    defect_type : str
        Type of defect model
    
    Returns:
    --------
    dict
        Optimal defect fraction and fit quality metrics
    """
    
    def objective(f_defect):
        """Calculate sum of squared log errors."""
        if f_defect <= 0 or f_defect >= 1:
            return 1e20
        
        defect_params = {'area_fraction': f_defect, 'type': defect_type}
        
        if defect_type == 'crack':
            defect_params['thickness_factor'] = 0.1
        elif defect_type == 'grain_boundary':
            defect_params['diffusivity_factor'] = 10
        
        errors = []
        for T_K, P_exp in zip(exp_temperatures, exp_permeabilities):
            try:
                P_model, _ = calculate_model_permeability(
                    T_K, P_test, model_level=3, 
                    defect_params=defect_params,
                    oxide_name=oxide_name
                )
                if P_model > 0 and P_exp > 0:
                    log_error = (np.log(P_model) - np.log(P_exp))**2
                    errors.append(log_error)
            except Exception:
                errors.append(100)  # Penalty for failed calculations
        
        return np.mean(errors) if errors else 1e20
    
    # Optimize defect fraction
    result = minimize_scalar(objective, bounds=(1e-6, 0.5), method='bounded')
    
    optimal_f = result.x
    
    # Calculate fit quality with optimal f
    defect_params = {'area_fraction': optimal_f, 'type': defect_type}
    if defect_type == 'crack':
        defect_params['thickness_factor'] = 0.1
    elif defect_type == 'grain_boundary':
        defect_params['diffusivity_factor'] = 10
    
    model_perms = []
    for T_K in exp_temperatures:
        P_model, _ = calculate_model_permeability(
            T_K, P_test, model_level=3, defect_params=defect_params
        )
        model_perms.append(P_model)
    
    model_perms = np.array(model_perms)
    
    # Calculate R²
    log_exp = np.log(exp_permeabilities)
    log_model = np.log(model_perms)
    ss_res = np.sum((log_exp - log_model)**2)
    ss_tot = np.sum((log_exp - np.mean(log_exp))**2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    
    # Mean absolute percentage error
    mape = np.mean(np.abs(model_perms - exp_permeabilities) / exp_permeabilities) * 100
    
    return {
        'optimal_f_defect': optimal_f,
        'defect_type': defect_type,
        'r_squared': r_squared,
        'mape': mape,
        'model_permeabilities': model_perms,
        'optimization_result': result
    }


def compare_all_levels(material='Incoloy800', P_test=1000, save_plots=True):
    """
    Compare Level 1, 2, and 3 predictions against experimental data.
    """
    print("\n" + "=" * 70)
    print("LEVEL 3 EXPERIMENTAL COMPARISON")
    print("=" * 70)
    
    # Get experimental data
    exp_data_raw = get_experimental_data('Incoloy 800', 'JAERI')
    exp_data = convert_to_SI(exp_data_raw)
    
    T_exp = exp_data['temperatures_K']
    P_exp = exp_data['permeabilities']
    
    print(f"\nExperimental data: {len(T_exp)} points")
    print(f"Temperature range: {T_exp.min()-273:.0f}°C to {T_exp.max()-273:.0f}°C")
    print(f"Test pressure: {P_test} Pa")
    
    # Calculate predictions for each level
    results = {
        'Level 1': {'permeabilities': [], 'label': 'Level 1 (Bare Metal)'},
        'Level 2': {'permeabilities': [], 'label': 'Level 2 (Perfect Oxide)'},
        'Level 3 (1% pinhole)': {'permeabilities': [], 'label': 'Level 3 (1% pinholes)'},
    }
    
    for T_K in T_exp:
        # Level 1
        P1, _ = calculate_model_permeability(T_K, P_test, model_level=1)
        results['Level 1']['permeabilities'].append(P1)
        
        # Level 2
        P2, _ = calculate_model_permeability(T_K, P_test, model_level=2)
        results['Level 2']['permeabilities'].append(P2)
        
        # Level 3 with 1% pinholes
        defect_params = {'area_fraction': 0.01, 'type': 'pinhole'}
        P3, _ = calculate_model_permeability(T_K, P_test, model_level=3, defect_params=defect_params)
        results['Level 3 (1% pinhole)']['permeabilities'].append(P3)
    
    for key in results:
        results[key]['permeabilities'] = np.array(results[key]['permeabilities'])
    
    # Fit optimal defect fraction
    print("\n" + "-" * 70)
    print("FITTING OPTIMAL DEFECT FRACTION")
    print("-" * 70)
    
    fit_results = {}
    for defect_type in ['pinhole', 'crack', 'grain_boundary']:
        fit = fit_defect_fraction(T_exp, P_exp, P_test, defect_type)
        fit_results[defect_type] = fit
        print(f"\n{defect_type.upper()}:")
        print(f"  Optimal f_defect = {fit['optimal_f_defect']*100:.3f}%")
        print(f"  R² = {fit['r_squared']:.4f}")
        print(f"  MAPE = {fit['mape']:.1f}%")
    
    # Find best fit
    best_type = max(fit_results.keys(), key=lambda x: fit_results[x]['r_squared'])
    best_fit = fit_results[best_type]
    
    print(f"\n→ Best fit: {best_type} with f = {best_fit['optimal_f_defect']*100:.3f}%")
    
    # Add best fit to results
    results['Level 3 (fitted)'] = {
        'permeabilities': best_fit['model_permeabilities'],
        'label': f"Level 3 (fitted: {best_fit['optimal_f_defect']*100:.2f}% {best_type})"
    }
    
    # Calculate error metrics for each level
    print("\n" + "-" * 70)
    print("ERROR ANALYSIS")
    print("-" * 70)
    print(f"\n{'Model':<30} | {'R²':>8} | {'MAPE (%)':>10} | {'Max Error (%)':>12}")
    print("-" * 70)
    
    for name, data in results.items():
        P_model = data['permeabilities']
        
        # R² on log scale
        log_exp = np.log(P_exp)
        log_model = np.log(P_model)
        ss_res = np.sum((log_exp - log_model)**2)
        ss_tot = np.sum((log_exp - np.mean(log_exp))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        
        # MAPE
        mape = np.mean(np.abs(P_model - P_exp) / P_exp) * 100
        
        # Max error
        max_err = np.max(np.abs(P_model - P_exp) / P_exp) * 100
        
        data['r_squared'] = r2
        data['mape'] = mape
        data['max_error'] = max_err
        
        print(f"{name:<30} | {r2:>8.4f} | {mape:>10.1f} | {max_err:>12.1f}")
    
    # Create comparison plot
    if save_plots:
        ensure_output_dir()
        create_comparison_plot(T_exp, P_exp, results, exp_data_raw, P_test)
        create_arrhenius_plot(T_exp, P_exp, results, exp_data_raw)
        create_prf_analysis_plot(T_exp, P_exp, best_fit, P_test)
    
    return results, fit_results


def create_comparison_plot(T_exp, P_exp, results, exp_data_raw, P_test):
    """Create main comparison plot."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Arrhenius plot
    ax1 = axes[0]
    x_exp = 1000 / T_exp
    
    # Experimental data
    ax1.semilogy(x_exp, P_exp, 'ko', markersize=10, markerfacecolor='none',
                 linewidth=2, label='Experimental (JAERI)', zorder=10)
    
    # Model predictions
    colors = ['blue', 'red', 'green', 'purple']
    linestyles = ['--', ':', '-', '-']
    
    for i, (name, data) in enumerate(results.items()):
        ax1.semilogy(x_exp, data['permeabilities'], 
                    color=colors[i % len(colors)], 
                    linestyle=linestyles[i % len(linestyles)],
                    linewidth=2, label=data['label'])
    
    ax1.set_xlabel('1000/T (K⁻¹)', fontsize=12)
    ax1.set_ylabel('log Permeability (mol/m/s/Pa⁰·⁵)', fontsize=12)
    ax1.set_title('Model Comparison: Permeability vs Temperature', fontsize=14)
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Add temperature axis on top
    ax1_top = ax1.twiny()
    temp_ticks = np.array([600, 700, 800, 900, 1000, 1100])
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(1000/(temp_ticks+273.15))
    ax1_top.set_xticklabels([f"{t}°C" for t in temp_ticks])
    
    # Plot 2: Parity plot for best model
    ax2 = axes[1]
    
    for name, data in results.items():
        if 'fitted' in name:
            ax2.loglog(P_exp, data['permeabilities'], 'go', markersize=8,
                      label=f"{name} (R²={data['r_squared']:.3f})")
        elif 'Level 1' in name:
            ax2.loglog(P_exp, data['permeabilities'], 'b^', markersize=6, alpha=0.5,
                      label=f"{name}")
    
    # 1:1 line
    p_range = [P_exp.min()*0.5, P_exp.max()*2]
    ax2.loglog(p_range, p_range, 'k--', linewidth=1, label='1:1 line')
    
    # ±50% bounds
    ax2.fill_between(p_range, [p*0.5 for p in p_range], [p*2 for p in p_range],
                     alpha=0.1, color='gray', label='±factor of 2')
    
    ax2.set_xlabel('log Experimental Permeability (mol/m/s/Pa⁰·⁵)', fontsize=12)
    ax2.set_ylabel('log Model Permeability (mol/m/s/Pa⁰·⁵)', fontsize=12)
    ax2.set_title('Parity Plot', fontsize=14)
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    
    filename = f"{OUTPUT_DIR}/level3_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"\n✓ Comparison plot saved: {filename}")
    plt.close()


def create_arrhenius_plot(T_exp, P_exp, results, exp_data_raw):
    """Create detailed Arrhenius plot."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    x_exp = 1000 / T_exp
    
    # Experimental
    ax.semilogy(x_exp, P_exp, 'ko', markersize=12, markerfacecolor='none',
               linewidth=2, label=f"Experimental ({exp_data_raw['source']})", zorder=10)
    
    # Models
    colors = {'Level 1': 'blue', 'Level 2': 'red', 
              'Level 3 (1% pinhole)': 'orange', 'Level 3 (fitted)': 'green'}
    
    for name, data in results.items():
        color = colors.get(name, 'gray')
        lw = 3 if 'fitted' in name else 2
        ax.semilogy(x_exp, data['permeabilities'], 
                   color=color, linewidth=lw, 
                   label=f"{data['label']} (R²={data['r_squared']:.3f})")
    
    ax.set_xlabel('1000/T (K⁻¹)', fontsize=14)
    ax.set_ylabel('log Permeability (mol/m/s/Pa⁰·⁵)', fontsize=14)
    ax.set_title('Level 1, 2, 3 Comparison with Experimental Data\n(Incoloy 800)', fontsize=16)
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Temperature axis
    ax_top = ax.twiny()
    temp_ticks = np.array([600, 700, 800, 900, 1000, 1100])
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(1000/(temp_ticks+273.15))
    ax_top.set_xticklabels([f"{t}°C" for t in temp_ticks])
    ax_top.set_xlabel('Temperature', fontsize=12)
    
    plt.tight_layout()
    
    filename = f"{OUTPUT_DIR}/arrhenius_all_levels_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"✓ Arrhenius plot saved: {filename}")
    plt.close()


def create_prf_analysis_plot(T_exp, P_exp, best_fit, P_test):
    """Create PRF analysis plot."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: PRF vs Temperature
    ax1 = axes[0]
    
    PRF_values = []
    for T_K in T_exp:
        oxide_props, metal_props = get_properties_at_temperature(T_K)
        
        defect_params = {
            'area_fraction': best_fit['optimal_f_defect'],
            'type': best_fit['defect_type']
        }
        
        result = calculate_PRF(P_test, oxide_props, metal_props, defect_params)
        PRF_values.append(result['PRF'])
    
    PRF_values = np.array(PRF_values)
    T_C = T_exp - 273.15
    
    ax1.semilogy(T_C, PRF_values, 'go-', linewidth=2, markersize=8,
                label=f"Level 3 ({best_fit['optimal_f_defect']*100:.2f}% {best_fit['defect_type']})")
    
    # Literature range
    ax1.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range (10-3828)')
    
    ax1.set_xlabel('Temperature (°C)', fontsize=12)
    ax1.set_ylabel('log PRF (Permeation Reduction Factor)', fontsize=12)
    ax1.set_title('PRF vs Temperature', fontsize=14)
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: PRF vs defect fraction at 800°C
    ax2 = axes[1]
    
    T_ref = 1073  # 800°C
    oxide_props, metal_props = get_properties_at_temperature(T_ref)
    
    f_defects = np.logspace(-4, -0.5, 30)
    PRF_pinhole = []
    PRF_crack = []
    
    for f in f_defects:
        # Pinhole
        result = calculate_PRF(P_test, oxide_props, metal_props, 
                              {'area_fraction': f, 'type': 'pinhole'})
        PRF_pinhole.append(result['PRF'])
        
        # Crack
        result = calculate_PRF(P_test, oxide_props, metal_props,
                              {'area_fraction': f, 'type': 'crack', 'thickness_factor': 0.1})
        PRF_crack.append(result['PRF'])
    
    ax2.loglog(f_defects * 100, PRF_pinhole, 'r-', linewidth=2, label='Pinholes')
    ax2.loglog(f_defects * 100, PRF_crack, 'b--', linewidth=2, label='Cracks (α=0.1)')
    
    # Mark fitted value
    ax2.axvline(best_fit['optimal_f_defect'] * 100, color='green', linestyle=':',
               linewidth=2, label=f"Fitted: {best_fit['optimal_f_defect']*100:.2f}%")
    
    # Literature range
    ax2.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature PRF range')
    
    ax2.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_ylabel('log PRF', fontsize=12)
    ax2.set_title(f'PRF vs Defect Fraction at 800°C', fontsize=14)
    ax2.legend(loc='best')
    ax2.grid(True, which='both', alpha=0.3)
    
    plt.tight_layout()
    
    filename = f"{OUTPUT_DIR}/PRF_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"✓ PRF analysis plot saved: {filename}")
    plt.close()


def save_results_report(results, fit_results, filename=None):
    """Save detailed results report."""
    ensure_output_dir()
    
    if filename is None:
        filename = f"{OUTPUT_DIR}/level3_comparison_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(filename, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("LEVEL 3 EXPERIMENTAL COMPARISON REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Material: Incoloy 800 with Cr2O3 oxide\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("MODEL COMPARISON SUMMARY\n")
        f.write("-" * 80 + "\n\n")
        
        f.write(f"{'Model':<35} | {'R²':>8} | {'MAPE (%)':>10}\n")
        f.write("-" * 60 + "\n")
        
        for name, data in results.items():
            f.write(f"{name:<35} | {data['r_squared']:>8.4f} | {data['mape']:>10.1f}\n")
        
        f.write("\n" + "-" * 80 + "\n")
        f.write("DEFECT FRACTION FITTING RESULTS\n")
        f.write("-" * 80 + "\n\n")
        
        for defect_type, fit in fit_results.items():
            f.write(f"\n{defect_type.upper()}:\n")
            f.write(f"  Optimal defect fraction: {fit['optimal_f_defect']*100:.4f}%\n")
            f.write(f"  R² (log scale): {fit['r_squared']:.4f}\n")
            f.write(f"  MAPE: {fit['mape']:.1f}%\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("CONCLUSIONS\n")
        f.write("=" * 80 + "\n\n")
        
        best_type = max(fit_results.keys(), key=lambda x: fit_results[x]['r_squared'])
        best = fit_results[best_type]
        
        f.write(f"Best fitting defect model: {best_type}\n")
        f.write(f"Optimal defect fraction: {best['optimal_f_defect']*100:.3f}%\n")
        f.write(f"This is {'within' if 0.001 <= best['optimal_f_defect'] <= 0.1 else 'outside'} ")
        f.write("the expected range (0.1% - 10%) from literature.\n\n")
        
        f.write("Physical interpretation:\n")
        if best['optimal_f_defect'] < 0.01:
            f.write("  - Low defect fraction suggests a well-formed oxide barrier\n")
            f.write("  - Permeation is primarily through intact oxide\n")
        elif best['optimal_f_defect'] < 0.05:
            f.write("  - Moderate defect fraction (typical for operating conditions)\n")
            f.write("  - Both oxide and defect paths contribute to permeation\n")
        else:
            f.write("  - High defect fraction indicates degraded oxide\n")
            f.write("  - Defects dominate permeation pathway\n")
    
    print(f"✓ Report saved: {filename}")


def main():
    """Run complete Level 3 experimental comparison."""
    print("\n" + "=" * 70)
    print("LEVEL 3: EXPERIMENTAL COMPARISON ANALYSIS")
    print("Comparing model predictions with JAERI experimental data")
    print("=" * 70)
    
    # Run comparison
    results, fit_results = compare_all_levels(
        material='Incoloy800',
        P_test=1000,  # 1000 Pa test pressure
        save_plots=True
    )
    
    # Save report
    save_results_report(results, fit_results)
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    
    return results, fit_results


if __name__ == "__main__":
    results, fit_results = main()