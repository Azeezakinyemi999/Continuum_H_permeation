import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from SALib.sample import morris as morris_sampler
from SALib.analyze import morris as morris_analyzer
from SALib.sample import saltelli as sobol_sampler
from SALib.analyze import sobol as sobol_analyzer

from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.utils import get_diffusivity, get_solubility, get_permeability
from data.material_data import MATERIALS

# =============================================================================
# DEFAULT PARAMETERS 
# =============================================================================

DEFAULT_PARAMS_LEVEL1 = {
    # Material: Incoloy800 at 800°C
    'D': 1.43e-9,              # m²/s (calculated at 800°C)
    'K_s': 5.92e-2,            # mol/m³/Pa^0.5 (calculated at 800°C)
    'thickness': 1e-3,         # m (1 mm)
    'P_up': 1e5,               # Pa (100 kPa - realistic value)
    'P_down': 0.0,             # Pa
    
    # Arrhenius parameters for Incoloy800
    'D_0': 6.40e-7,            # m²/s
    'E_D': 54000,              # J/mol
    'K_s0': 1.80e2,            # mol/m³/Pa^0.5
    'H_s': -21924,             # J/mol
    'temperature': 1073.15,    # K (800°C)
    
    # Material name
    'material_name': 'Incoloy800'
}

# Valid output metrics for sensitivity analysis
VALID_OUTPUT_METRICS = ['flux', 'permeability', 'C_up', 'C_down', 'D', 'K_s']

# =============================================================================
# MODEL WRAPPER FUNCTIONS
# =============================================================================
def level1_model_wrapper(params_dict):
    """
    Wrapper for LEVEL 1 simple metal permeation model.
    Works with direct parameters (D, K_s) OR Arrhenius parameters.
    
    Now properly handles temperature sensitivity:
    - If 'temperature' is varied, D and K_s are ALWAYS calculated from Arrhenius
    - Uses default Arrhenius parameters (D_0, E_D, K_s0, H_s) if not provided
    
    Parameters:
    -----------
    params_dict : dict
        Can contain:
        - Direct: 'D', 'K_s', 'thickness', 'P_up', 'P_down'
        - Arrhenius: 'D_0', 'E_D', 'K_s0', 'H_s', 'temperature'
        
    Returns:
    --------
    dict : Output metrics
        - 'flux': Hydrogen flux (mol/m²/s)
        - 'C_up': Upstream concentration (mol/m³)
        - 'C_down': Downstream concentration (mol/m³)
        - 'permeability': D × K_s (mol/m/s/Pa^0.5)
        - 'D': Diffusivity used (m²/s)
        - 'K_s': Solubility used (mol/m³/Pa^0.5)
        - 'temperature': Temperature (K)
    """
    # Merge with defaults
    full_params = DEFAULT_PARAMS_LEVEL1.copy()
    full_params.update(params_dict)
    
    try:
        R = 8.314  # J/mol/K
        T = full_params['temperature']
        
        # Determine if we should use Arrhenius equations
        # Use Arrhenius if:
        #   1. Temperature is being varied (in params_dict), OR
        #   2. Any Arrhenius parameter is being varied
        use_arrhenius_D = (
            'temperature' in params_dict or 
            'D_0' in params_dict or 
            'E_D' in params_dict
        )
        
        use_arrhenius_Ks = (
            'temperature' in params_dict or 
            'K_s0' in params_dict or 
            'H_s' in params_dict
        )
        
        # Calculate D
        if use_arrhenius_D:
            D_0 = full_params['D_0']
            E_D = full_params['E_D']
            D = D_0 * np.exp(-E_D / (R * T))
        else:
            D = full_params['D']
        
        # Calculate K_s
        if use_arrhenius_Ks:
            K_s0 = full_params['K_s0']
            H_s = full_params['H_s']
            K_s = K_s0 * np.exp(-H_s / (R * T))
        else:
            K_s = full_params['K_s']
        
        # Get other parameters
        thickness = full_params['thickness']
        P_up = full_params['P_up']
        P_down = full_params['P_down']
        
        # Validate inputs
        if D <= 0 or K_s <= 0 or thickness <= 0:
            raise ValueError(f"D ({D}), K_s ({K_s}), and thickness ({thickness}) must be positive")
        if P_up < 0 or P_down < 0:
            raise ValueError("Pressures must be non-negative")
        if P_up < P_down:
            raise ValueError("P_up must be >= P_down")
        
        # Calculate flux using your function
        result = calculate_simple_metal_flux(D, K_s, thickness, P_up, P_down)
        
        # Add D, K_s, and temperature to result
        result['D'] = D
        result['K_s'] = K_s
        result['temperature'] = T
        
        return result
        
    except Exception as e:
        print(f"Error in Level 1 model: {e}")
        print(f"  params_dict: {params_dict}")
        # Return safe defaults
        return {
            'flux': 1e-20,
            'C_up': 0.0,
            'C_down': 0.0,
            'permeability': 1e-20,
            'D': 1e-12,
            'K_s': 1e-6,
            'temperature': full_params.get('temperature', 1073.15)
        }


# =============================================================================
# MORRIS SENSITIVITY ANALYSIS
# =============================================================================

def morris_sensitivity_level1(
    param_ranges,
    N_trajectories=10,
    num_levels=4,
    output_metric='flux',
    use_arrhenius=False
):
    """
    Morris sensitivity analysis for LEVEL 1 clean metal model.
    
    Parameters:
    -----------
    param_ranges : dict
        Parameter ranges.
        
    N_trajectories : int
        Number of trajectories (default: 10)
    num_levels : int
        Number of grid levels (default: 4)
    output_metric : str
        Output to analyze: 'flux', 'permeability', 'C_up', 'D', 'K_s'
    use_arrhenius : bool
        Whether ranges include Arrhenius parameters (for documentation)
        
    Returns:
    --------
    Si : Morris sensitivity indices
    problem : Problem definition
    Y : Model outputs
    """
    # Validate output metric
    if output_metric not in VALID_OUTPUT_METRICS:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS}, got '{output_metric}'")
    
    # Check if output makes sense for the parameter ranges
    if output_metric in ['D', 'K_s']:
        has_temp_or_arrhenius = (
            'temperature' in param_ranges or 
            'D_0' in param_ranges or 
            'E_D' in param_ranges or
            'K_s0' in param_ranges or
            'H_s' in param_ranges
        )
        if not has_temp_or_arrhenius:
            print(f"\n⚠️  WARNING: Analyzing '{output_metric}' but no temperature or Arrhenius ")
            print(f"    parameters in param_ranges. '{output_metric}' will be constant!")
            print(f"    Include 'temperature' or Arrhenius parameters to see variation.\n")
    
    # Define problem
    problem = {
        'num_vars': len(param_ranges),
        'names': list(param_ranges.keys()),
        'bounds': list(param_ranges.values())
    }
    
    # Generate Morris samples
    param_values = morris_sampler.sample(
        problem, 
        N=N_trajectories, 
        num_levels=num_levels
    )
    
    # Run model for each sample
    Y = np.zeros(param_values.shape[0])
    
    print(f"\n{'='*70}")
    print(f"MORRIS SENSITIVITY ANALYSIS - LEVEL 1 (Clean Metal)")
    print(f"{'='*70}")
    print(f"Running {param_values.shape[0]} Morris samples...")
    print(f"Varying parameters: {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print(f"{'='*70}\n")
    
    for i, X in enumerate(param_values):
        sample_params = dict(zip(problem['names'], X))
        result = level1_model_wrapper(sample_params)
        Y[i] = result[output_metric]
        
        if (i + 1) % 10 == 0:
            print(f"  Completed {i + 1}/{param_values.shape[0]} samples")
    
    # Handle invalid values
    valid_idx = np.isfinite(Y) & (Y > 0)
    n_invalid = np.sum(~valid_idx)
    if n_invalid > 0:
        print(f"\nWarning: {n_invalid} invalid outputs detected")
        if np.any(valid_idx):
            median_val = np.median(Y[valid_idx])
            Y[~valid_idx] = median_val
            print(f"  Replaced with median: {median_val:.2e}")
        else:
            Y[:] = 1e-15
    
    # Analyze
    Si = morris_analyzer.analyze(
        problem, 
        param_values, 
        Y, 
        conf_level=0.95,
        print_to_console=False
    )
    
    print(f"\n✓ Morris analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    
    return Si, problem, Y


# =============================================================================
# SOBOL SENSITIVITY ANALYSIS
# =============================================================================

def sobol_sensitivity_level1(
    param_ranges,
    N_samples=1024,
    output_metric='flux',
    calc_second_order=True
):
    """
    Sobol variance-based sensitivity analysis for LEVEL 1.
    
    Parameters:
    -----------
    param_ranges : dict
        Parameter ranges (same format as Morris)
    N_samples : int
        Number of samples (must be power of 2: 256, 512, 1024, etc.)
    output_metric : str
        Output to analyze: 'flux', 'permeability', 'C_up', 'D', 'K_s'
    calc_second_order : bool
        Calculate second-order interactions (slower but more informative)
        
    Returns:
    --------
    Si : Sobol indices (S1, ST, S2 if calc_second_order=True)
    problem : Problem definition
    Y : Model outputs
    """
    # Validate output metric
    if output_metric not in VALID_OUTPUT_METRICS:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS}, got '{output_metric}'")
    
    # Check if output makes sense for the parameter ranges
    if output_metric in ['D', 'K_s']:
        has_temp_or_arrhenius = (
            'temperature' in param_ranges or 
            'D_0' in param_ranges or 
            'E_D' in param_ranges or
            'K_s0' in param_ranges or
            'H_s' in param_ranges
        )
        if not has_temp_or_arrhenius:
            print(f"\n⚠️  WARNING: Analyzing '{output_metric}' but no temperature or Arrhenius ")
            print(f"    parameters in param_ranges. '{output_metric}' will be constant!")
            print(f"    Include 'temperature' or Arrhenius parameters to see variation.\n")
    
    problem = {
        'num_vars': len(param_ranges),
        'names': list(param_ranges.keys()),
        'bounds': list(param_ranges.values())
    }
    
    # Generate Sobol samples (Saltelli)
    param_values = sobol_sampler.sample(
        problem, 
        N_samples,
        calc_second_order=calc_second_order
    )
    
    Y = np.zeros(param_values.shape[0])
    
    print(f"\n{'='*70}")
    print(f"SOBOL SENSITIVITY ANALYSIS - LEVEL 1 (Clean Metal)")
    print(f"{'='*70}")
    print(f"Running {param_values.shape[0]} Sobol samples...")
    print(f"Varying parameters: {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print("This may take a while...")
    print(f"{'='*70}\n")
    
    for i, X in enumerate(param_values):
        sample_params = dict(zip(problem['names'], X))
        result = level1_model_wrapper(sample_params)
        Y[i] = result[output_metric]
        
        if (i + 1) % 500 == 0:
            print(f"  Completed {i + 1}/{param_values.shape[0]} samples")
    
    # Handle invalid values
    valid_idx = np.isfinite(Y) & (Y > 0)
    if not np.all(valid_idx):
        print(f"\nWarning: {np.sum(~valid_idx)} invalid outputs")
        if np.any(valid_idx):
            Y[~valid_idx] = np.median(Y[valid_idx])
    
    # Analyze
    Si = sobol_analyzer.analyze(
        problem, 
        Y,
        calc_second_order=calc_second_order,
        conf_level=0.95,
        print_to_console=False
    )
    
    print(f"\n✓ Sobol analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    
    return Si, problem, Y
# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================
def plot_morris_results(Si, problem, output_metric='Model Output'):
    """
    Visualize Morris sensitivity results.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    param_names = problem['names']
    mu_star = np.nan_to_num(Si['mu_star'], nan=0.0)
    sigma = np.nan_to_num(Si['sigma'], nan=0.0)
    
    # Sort by importance
    sorted_idx = np.argsort(mu_star)[::-1]
    sorted_names = [param_names[i] for i in sorted_idx]
    sorted_mu_star = mu_star[sorted_idx]
    
    # Bar chart
    ax1.barh(sorted_names, sorted_mu_star, color='steelblue', alpha=0.7)
    ax1.set_xlabel('μ* (Mean Absolute Elementary Effect)', fontsize=12)
    ax1.set_title(f'Morris Sensitivity Analysis\n{output_metric}', fontsize=14, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Scatter plot
    ax2.scatter(mu_star, sigma, s=150, alpha=0.6, c='steelblue', edgecolors='black')
    
    for i, name in enumerate(param_names):
        ax2.annotate(name, (mu_star[i], sigma[i]), 
                    fontsize=10, ha='right', va='bottom',
                    xytext=(-5, 5), textcoords='offset points')
    
    ax2.set_xlabel('μ* (Importance)', fontsize=12)
    ax2.set_ylabel('σ (Nonlinearity/Interactions)', fontsize=12)
    ax2.set_title('Elementary Effects', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add interpretation zones
    if np.any(mu_star > 0):
        max_mu = np.max(mu_star)
        max_sigma = np.max(sigma)
        if max_mu > 0:
            ax2.axvline(max_mu * 0.3, color='red', linestyle='--', alpha=0.3, 
                       label='High importance')
        if max_sigma > 0:
            ax2.axhline(max_sigma * 0.3, color='blue', linestyle='--', alpha=0.3, 
                       label='High nonlinearity')
        ax2.legend(fontsize=10)
    
    plt.tight_layout()
    plt.show()
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"MORRIS SENSITIVITY RESULTS - {output_metric}")
    print(f"{'='*70}")
    df = pd.DataFrame({
        'Parameter': sorted_names,
        'μ*': sorted_mu_star,
        'σ': [sigma[i] for i in sorted_idx],
        'μ*_conf': [Si['mu_star_conf'][i] for i in sorted_idx]
    })
    print(df.to_string(index=False))
    print(f"{'='*70}")
    print("\nInterpretation:")
    print("  μ* : Overall importance (higher = more influential)")
    print("  σ  : Nonlinearity/Interactions (higher = more complex)")
    print(f"{'='*70}\n")


def plot_sobol_results(Si, problem, output_metric='Model Output'):
    """
    Visualize Sobol sensitivity results.
    """
    param_names = problem['names']
    S1 = Si['S1']
    ST = Si['ST']
    
    has_S2 = 'S2' in Si and Si['S2'] is not None
    n_plots = 2 if has_S2 else 1
    
    fig, axes = plt.subplots(1, n_plots, figsize=(7*n_plots, 6))
    if n_plots == 1:
        axes = [axes]
    
    # First-order and Total-order indices
    ax1 = axes[0]
    x = np.arange(len(param_names))
    width = 0.35
    
    ax1.bar(x - width/2, S1, width, label='S1 (First-order)', alpha=0.8, color='steelblue')
    ax1.bar(x + width/2, ST, width, label='ST (Total-order)', alpha=0.8, color='coral')
    
    ax1.set_xlabel('Parameters', fontsize=12)
    ax1.set_ylabel('Sobol Index', fontsize=12)
    ax1.set_title(f'Sobol Sensitivity\n{output_metric}', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(param_names, rotation=45, ha='right')
    ax1.legend(fontsize=11)
    ax1.grid(axis='y', alpha=0.3)
    
    # Second-order interactions
    if has_S2:
        ax2 = axes[1]
        S2 = Si['S2']
        
        im = ax2.imshow(S2, cmap='RdYlGn', aspect='auto', vmin=0, vmax=np.max(S2))
        
        ax2.set_xticks(np.arange(len(param_names)))
        ax2.set_yticks(np.arange(len(param_names)))
        ax2.set_xticklabels(param_names, rotation=45, ha='right')
        ax2.set_yticklabels(param_names)
        
        cbar = plt.colorbar(im, ax=ax2)
        cbar.set_label('S2', fontsize=11)
        
        # Annotate significant interactions
        for i in range(len(param_names)):
            for j in range(len(param_names)):
                if i != j and S2[i, j] > 0.01:
                    ax2.text(j, i, f'{S2[i, j]:.3f}',
                           ha="center", va="center", color="black", fontsize=9)
        
        ax2.set_title('Second-Order Interactions', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"SOBOL SENSITIVITY RESULTS - {output_metric}")
    print(f"{'='*70}")
    df = pd.DataFrame({
        'Parameter': param_names,
        'S1': S1,
        'S1_conf': Si['S1_conf'],
        'ST': ST,
        'ST_conf': Si['ST_conf']
    })
    df['Interactions'] = df['ST'] - df['S1']
    print(df.to_string(index=False))
    print(f"{'='*70}")
    print("\nInterpretation:")
    print("  S1 : First-order effect (direct influence)")
    print("  ST : Total-order effect (direct + interactions)")
    print("  ST-S1 : Pure interaction effects")
    print(f"{'='*70}\n")
    
    # Top parameters
    important = df.nlargest(3, 'ST')
    print("Top 3 Most Influential Parameters:")
    for idx, row in important.iterrows():
        print(f"  {row['Parameter']}: ST = {row['ST']:.3f} (S1 = {row['S1']:.3f})")
    print(f"{'='*70}\n")


###########################################################################
##################### LEVEL 2 SENSITIVITY ANALYSIS ########################
###########################################################################
###########################################################################











###########################################################################
##################### LEVEL 3 SENSITIVITY ANALYSIS ########################
###########################################################################
###########################################################################

###########################################################################
##################### LEVEL 4 SENSITIVITY ANALYSIS ########################
###########################################################################
###########################################################################



###########################################################################
##################### LEVEL 5 SENSITIVITY ANALYSIS ########################
###########################################################################
###########################################################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from SALib.sample import morris as morris_sampler
from SALib.analyze import morris as morris_analyzer
from SALib.sample import saltelli as sobol_sampler
from SALib.analyze import sobol as sobol_analyzer

# =============================================================================
# DEFAULT PARAMETERS FOR LEVEL 5 (Complete Hierarchical Model)
# =============================================================================

DEFAULT_PARAMS_LEVEL5 = {
    # =========================================================================
    # LEVEL 1: Base Metal Properties (Incoloy 800)
    # =========================================================================
    # Arrhenius diffusivity: D = D_0 * exp(-E_D / RT)
    'D_0': 5.00e-7,            # m²/s (pre-exponential)
    'E_D': 52000,              # J/mol (activation energy)
    
    # Arrhenius solubility: K_s = K_s0 * exp(-H_s / RT)
    'K_s0': 1.0e-9,            # mol/m³/Pa^0.5 (pre-exponential)
    'H_s': -20000,             # J/mol (heat of solution, negative = exothermic)
    
    # Geometry
    'metal_thickness': 1e-3,   # m (1 mm)
    
    # Operating conditions
    'P_upstream': 1e5,         # Pa (100 kPa)
    'P_downstream': 0.0,       # Pa
    'temperature': 1073.15,    # K (800°C)
    
    # =========================================================================
    # LEVEL 2: Oxide Layer Properties (Cr2O3)
    # =========================================================================
    # Arrhenius diffusivity: D_ox = D_ox_0 * exp(-E_D_ox / RT)
    'D_ox_0': 1e-6,            # m²/s (pre-exponential)
    'E_D_ox': 1.55e5,          # J/mol (activation energy)
    
    # Arrhenius solubility: K_ox = K_ox_0 * exp(-H_sol_ox / RT)
    'K_ox_0': 1e-4,            # mol/m³/Pa (pre-exponential)
    'H_sol_ox': 1.85e5,        # J/mol (solution enthalpy)
    
    # Oxide geometry
    'oxide_thickness': 1e-6,   # m (1 μm)
    
    # =========================================================================
    # LEVEL 3: Oxide Defect Properties
    # =========================================================================
    'defect_fraction': 0.01,   # Area fraction with defects (1%)
    'defect_type': 'pinhole',  # 'pinhole', 'crack', 'grain_boundary', 'mixed'
    'crack_thickness_factor': 0.1,    # For cracks: oxide thickness = factor × nominal
    'gb_diffusivity_factor': 10.0,    # For GB: D_defect = factor × D_oxide
    
    # =========================================================================
    # LEVEL 4: Metal Microstructure Properties
    # =========================================================================
    # Grain structure
    'grain_size': 50e-6,       # m (50 μm) - large enough for trapping to dominate
    'grain_shape': 'equiaxed', # 'equiaxed', 'columnar', or 'planar'
    'gb_type': 'HAGB',         # 'HAGB', 'LAGB', 'twin', 'special'
    'gb_thickness': 0.5e-9,    # m (0.5 nm)
    
    # GB diffusion enhancement: D_gb = gb_enhancement_factor × D_lattice
    'gb_enhancement_factor': 100,  # At 800°C (temperature-dependent in full model)
    
    # Lattice site density
    'lattice_density': 1.06e29,  # m⁻³ (FCC interstitial sites)
    
    # =========================================================================
    # LEVEL 4: Trap Properties
    # =========================================================================
    # Dislocation traps
    'trap_dislocation_E_b': 15e3,     # J/mol (15 kJ/mol binding energy)
    'trap_dislocation_N_T': 1e15,     # m⁻³ (trap density)
    
    # Grain boundary traps (separate from GB fast diffusion)
    'trap_gb_E_b': 20e3,              # J/mol (20 kJ/mol)
    'sites_per_area': 1e19,           # trap sites/m² of GB area
    # N_T for GB calculated from grain_size internally
    
    # Vacancy traps
    'trap_vacancy_E_b': 45e3,         # J/mol (45 kJ/mol)
    'trap_vacancy_N_T': 1e23,         # m⁻³ (T-dependent in full model)
    
    # Carbide interface traps
    'trap_carbide_E_b': 72e3,         # J/mol (72 kJ/mol - strong trap)
    'trap_carbide_N_T': 1e21,         # m⁻³ (aged condition)
    
    # =========================================================================
    # MODEL OPTIONS
    # =========================================================================
    'include_gb_enhancement': True,   # Include GB fast diffusion path
    'include_trapping': True,         # Include Oriani trapping
    'D_eff_method': 'average',        # 'average', 'harmonic', 'inlet', 'outlet'
}

# =============================================================================
# VALID OUTPUT METRICS
# =============================================================================
VALID_OUTPUT_METRICS_L5 = [
    'flux',              # Total permeation flux [mol/m²/s]
    'PRF',               # Permeation Reduction Factor [-]
    'D_eff',             # Effective metal diffusivity [m²/s]
    'D_modification',    # D_eff / D_lattice [-]
    'permeability',      # Effective permeability [mol/m/s/Pa^0.5]
    'P_interface',       # Oxide-metal interface pressure [Pa]
    'flux_intact',       # Flux through intact oxide [mol/m²/s]
    'flux_defect',       # Flux through defect paths [mol/m²/s]
]

# =============================================================================
# PARAMETER GROUPS (for organized sensitivity analysis)
# =============================================================================
PARAM_GROUPS = {
    'level1_metal': ['D_0', 'E_D', 'K_s0', 'H_s', 'metal_thickness'],
    'level2_oxide': ['D_ox_0', 'E_D_ox', 'K_ox_0', 'H_sol_ox', 'oxide_thickness'],
    'level3_defects': ['defect_fraction', 'crack_thickness_factor', 'gb_diffusivity_factor'],
    'level4_microstructure': ['grain_size', 'gb_enhancement_factor'],
    'level4_traps': ['trap_dislocation_E_b', 'trap_dislocation_N_T',
                     'trap_vacancy_E_b', 'trap_vacancy_N_T',
                     'trap_carbide_E_b', 'trap_carbide_N_T'],
    'operating': ['temperature', 'P_upstream'],
}

# =============================================================================
# SUGGESTED PARAMETER RANGES FOR SENSITIVITY ANALYSIS (FULL 30 PARAMS)
# =============================================================================
SUGGESTED_RANGES_LEVEL5 = {
    # =========================================================================
    # LEVEL 1: Base Metal Properties (6 params)
    # =========================================================================
    'D_0': [1e-8, 1e-5],               # m²/s (pre-exponential diffusivity)
    'E_D': [40000, 70000],             # J/mol (activation energy)
    'K_s0': [1e-10, 1e-7],             # mol/m³/Pa^0.5 (pre-exponential solubility)
    'H_s': [-40000, 0],                # J/mol (heat of solution)
    'metal_thickness': [0.5e-3, 5e-3], # m (0.5-5 mm)
    'P_upstream': [1e3, 1e6],          # Pa (1 kPa to 1 MPa)
    
    # =========================================================================
    # LEVEL 2: Oxide Layer Properties (5 params)
    # =========================================================================
    'D_ox_0': [1e-8, 1e-4],            # m²/s (oxide diffusivity pre-exp)
    'E_D_ox': [1.2e5, 2.0e5],          # J/mol (oxide activation energy)
    'K_ox_0': [1e-6, 1e-2],            # mol/m³/Pa (oxide solubility pre-exp)
    'H_sol_ox': [1.5e5, 2.2e5],        # J/mol (oxide solution enthalpy)
    'oxide_thickness': [1e-8, 1e-5],   # m (10 nm to 10 μm)
    
    # =========================================================================
    # LEVEL 3: Oxide Defect Properties (3 params)
    # =========================================================================
    'defect_fraction': [0.0001, 0.1],      # 0.01% to 10%
    'crack_thickness_factor': [0.01, 0.5], # Crack thinning factor
    'gb_diffusivity_factor': [1.0, 100],   # GB enhancement in oxide
    
    # =========================================================================
    # LEVEL 4: Grain Structure (5 params)
    # =========================================================================
    'grain_size': [1e-7, 1e-4],            # m (0.1-100 μm)
    'gb_thickness': [0.3e-9, 1e-9],        # m (0.3-1 nm)
    'gb_enhancement_factor': [10, 1000],   # D_gb/D_lattice ratio
    'lattice_density': [8e28, 1.2e29],     # m⁻³ (lattice sites)
    'sites_per_area': [1e18, 1e20],        # trap sites/m² of GB area
    
    # =========================================================================
    # LEVEL 4: Dislocation Traps (2 params)
    # =========================================================================
    'trap_dislocation_E_b': [10e3, 35e3],  # J/mol (10-35 kJ/mol)
    'trap_dislocation_N_T': [1e13, 1e17],  # m⁻³ (annealed to cold-worked)
    
    # =========================================================================
    # LEVEL 4: Grain Boundary Traps (1 param)
    # =========================================================================
    'trap_gb_E_b': [15e3, 40e3],           # J/mol (15-40 kJ/mol)
    
    # =========================================================================
    # LEVEL 4: Vacancy Traps (2 params)
    # =========================================================================
    'trap_vacancy_E_b': [30e3, 60e3],      # J/mol (30-60 kJ/mol)
    'trap_vacancy_N_T': [1e20, 1e25],      # m⁻³ (T-dependent)
    
    # =========================================================================
    # LEVEL 4: Carbide/Precipitate Traps (2 params)
    # =========================================================================
    'trap_carbide_E_b': [50e3, 100e3],     # J/mol (50-100 kJ/mol, deep traps)
    'trap_carbide_N_T': [1e19, 1e23],      # m⁻³ (aging dependent)
    
    # =========================================================================
    # Operating Conditions (1 param)
    # =========================================================================
    'temperature': [573, 1273],            # K (300-1000°C)
}
# Total: 30 parameters


# =============================================================================
# MODEL WRAPPER FUNCTIONS
# =============================================================================

def level5_model_wrapper(params_dict):
    """
    Wrapper for LEVEL 5 complete hierarchical model.
    
    Combines all levels:
    - Level 1: Base metal (Arrhenius D, K_s)
    - Level 2: Oxide layer (barrier)
    - Level 3: Oxide defects (parallel paths)
    - Level 4: Metal microstructure (GB enhancement + trapping)
    
    Parameters
    ----------
    params_dict : dict
        Any subset of parameters from DEFAULT_PARAMS_LEVEL5.
        Missing parameters are filled from defaults.
        
    Returns
    -------
    dict
        Output metrics:
        - 'flux': Total permeation flux [mol/m²/s]
        - 'PRF': Permeation Reduction Factor [-]
        - 'D_eff': Effective metal diffusivity [m²/s]
        - 'D_modification': D_eff / D_lattice [-]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
        - 'P_interface': Oxide-metal interface pressure [Pa]
        - 'flux_intact': Flux through intact oxide path [mol/m²/s]
        - 'flux_defect': Flux through defect paths [mol/m²/s]
        - 'flux_bare_metal': Flux without oxide (Level 4 only) [mol/m²/s]
        - 'regime': Operating regime classification
    """
    from calculations.parallel_oxide_defect_paths import (
        calculate_parallel_path_flux_defective_metal,
        calculate_PRF_defective_metal
    )
    from calculations.permeation_calc import calculate_defective_metal_flux
    
    # Merge with defaults
    full_params = DEFAULT_PARAMS_LEVEL5.copy()
    full_params.update(params_dict)
    
    try:
        R = 8.314  # J/mol/K
        T = full_params['temperature']
        
        # =====================================================================
        # Calculate temperature-dependent properties
        # =====================================================================
        
        # Metal diffusivity: D = D_0 * exp(-E_D / RT)
        D_0 = full_params['D_0']
        E_D = full_params['E_D']
        D_metal = D_0 * np.exp(-E_D / (R * T))
        
        # Metal solubility: K_s = K_s0 * exp(-H_s / RT)
        K_s0 = full_params['K_s0']
        H_s = full_params['H_s']
        K_s_metal = K_s0 * np.exp(-H_s / (R * T))
        
        # Oxide diffusivity: D_ox = D_ox_0 * exp(-E_D_ox / RT)
        D_ox_0 = full_params['D_ox_0']
        E_D_ox = full_params['E_D_ox']
        D_ox = D_ox_0 * np.exp(-E_D_ox / (R * T))
        
        # Oxide solubility: K_ox = K_ox_0 * exp(-H_sol_ox / RT)
        K_ox_0 = full_params['K_ox_0']
        H_sol_ox = full_params['H_sol_ox']
        K_ox = K_ox_0 * np.exp(-H_sol_ox / (R * T))
        
        # =====================================================================
        # Build property dictionaries for calculation functions
        # =====================================================================
        
        oxide_props = {
            'D_ox': D_ox,
            'K_ox': K_ox,
            'thickness': full_params['oxide_thickness']
        }
        
        metal_props = {
            'D_metal': D_metal,  # This is D_lattice for Level 4
            'K_s_metal': K_s_metal,
            'thickness': full_params['metal_thickness']
        }
        
        # Build trap list from individual trap parameters
        trap_list = []
        
        # Dislocation traps
        if full_params.get('trap_dislocation_N_T', 0) > 0:
            trap_list.append({
                'name': 'dislocations',
                'binding_energy': full_params['trap_dislocation_E_b'],
                'density': full_params['trap_dislocation_N_T']
            })
        
        # Vacancy traps
        if full_params.get('trap_vacancy_N_T', 0) > 0:
            trap_list.append({
                'name': 'vacancies',
                'binding_energy': full_params['trap_vacancy_E_b'],
                'density': full_params['trap_vacancy_N_T']
            })
        
        # Carbide traps
        if full_params.get('trap_carbide_N_T', 0) > 0:
            trap_list.append({
                'name': 'carbides',
                'binding_energy': full_params['trap_carbide_E_b'],
                'density': full_params['trap_carbide_N_T']
            })
        
        # Grain boundary traps (density calculated from grain size)
        grain_size = full_params['grain_size']
        gb_thickness = full_params['gb_thickness']
        sites_per_area = full_params['sites_per_area']  # trap sites/m² of GB area
        # GB area per unit volume ≈ 3/grain_size for equiaxed grains
        gb_area_per_volume = 3.0 / grain_size
        N_T_gb = sites_per_area * gb_area_per_volume * gb_thickness
        
        if full_params.get('trap_gb_E_b', 0) > 0 and N_T_gb > 0:
            trap_list.append({
                'name': 'grain_boundaries',
                'binding_energy': full_params['trap_gb_E_b'],
                'density': N_T_gb
            })
        
        microstructure_params = {
            'grain_size': grain_size,
            'grain_shape': full_params['grain_shape'],
            'gb_type': full_params['gb_type'],
            'gb_thickness': gb_thickness,
            'trap_list': trap_list,
            'gb_enhancement_factor': full_params.get('gb_enhancement_factor', 100)
        }
        
        # Defect parameters for Level 3
        defect_params = {
            'area_fraction': full_params['defect_fraction'],
            'type': full_params['defect_type'],
            'thickness_factor': full_params.get('crack_thickness_factor', 0.1),
            'diffusivity_factor': full_params.get('gb_diffusivity_factor', 10.0)
        }
        
        # Operating conditions
        P_upstream = full_params['P_upstream']
        P_downstream = full_params['P_downstream']
        lattice_density = full_params['lattice_density']
        method = full_params.get('D_eff_method', 'average')
        
        # Determine mode based on flags
        include_gb = full_params.get('include_gb_enhancement', True)
        include_trap = full_params.get('include_trapping', True)
        if include_gb and include_trap:
            mode = 'both'
        elif include_gb:
            mode = 'gb_only'
        elif include_trap:
            mode = 'trapping_only'
        else:
            mode = 'none'
        
        # =====================================================================
        # Calculate Level 5: Full hierarchical model
        # =====================================================================
        
        # Full Level 3+4 calculation
        result_l5 = calculate_parallel_path_flux_defective_metal(
            P_upstream=P_upstream,
            P_downstream=P_downstream,
            oxide_props=oxide_props,
            metal_props=metal_props,
            defect_params=defect_params,
            temperature=T,
            microstructure_params=microstructure_params,
            lattice_density=lattice_density,
            method=method,
            n_points=10,
            mode=mode
        )
        
        # Calculate bare metal flux for PRF (Level 4 only, no oxide)
        result_bare = calculate_defective_metal_flux(
            D_lattice=D_metal,
            K_s=K_s_metal,
            thickness=metal_props['thickness'],
            P_up=P_upstream,
            P_down=P_downstream,
            temperature=T,
            microstructure_params=microstructure_params,
            lattice_density=lattice_density,
            method=method,
            n_points=10,
            mode=mode
        )
        
        flux_total = result_l5['flux_total']
        flux_bare = result_bare['flux']
        
        # PRF = flux_bare / flux_with_oxide
        PRF = flux_bare / flux_total if flux_total > 0 else float('inf')
        
        # Effective permeability
        D_eff = result_l5.get('D_eff_metal', D_metal)
        permeability = D_eff * K_s_metal
        
        # D modification factor
        D_modification = D_eff / D_metal if D_metal > 0 else 1.0
        
        # =====================================================================
        # Return comprehensive results
        # =====================================================================
        
        return {
            # Primary outputs
            'flux': flux_total,
            'PRF': PRF,
            'D_eff': D_eff,
            'D_modification': D_modification,
            'permeability': permeability,
            
            # Interface and path details
            'P_interface': result_l5.get('P_interface_intact', 0),
            'flux_intact': result_l5['flux_intact_contribution'],
            'flux_defect': result_l5['flux_defect_contribution'],
            'flux_bare_metal': flux_bare,
            
            # Regime classification
            'regime': result_l5.get('regime', 'unknown'),
            'dominant_path': result_l5.get('dominant_path', 'unknown'),
            
            # Temperature-dependent calculated values
            'D_metal': D_metal,
            'K_s_metal': K_s_metal,
            'D_ox': D_ox,
            'K_ox': K_ox,
            'temperature': T,
            
            # Microstructure details
            'modification_factor': result_l5.get('modification_factor', 1.0),
            'defect_enhancement': result_l5.get('defect_enhancement_factor', 1.0),
        }
        
    except Exception as e:
        print(f"Error in Level 5 model: {e}")
        print(f"  params_dict keys: {list(params_dict.keys())}")
        import traceback
        traceback.print_exc()
        
        # Return safe defaults
        return {
            'flux': 1e-20,
            'PRF': 1.0,
            'D_eff': 1e-12,
            'D_modification': 1.0,
            'permeability': 1e-20,
            'P_interface': 0,
            'flux_intact': 1e-20,
            'flux_defect': 0,
            'flux_bare_metal': 1e-20,
            'regime': 'error',
            'dominant_path': 'error',
            'D_metal': 1e-12,
            'K_s_metal': 1e-6,
            'D_ox': 1e-15,
            'K_ox': 1e-10,
            'temperature': full_params.get('temperature', 1073.15),
            'modification_factor': 1.0,
            'defect_enhancement': 1.0,
        }


# (Add after level5_model_wrapper function)

# =============================================================================
# MORRIS SENSITIVITY ANALYSIS - LEVEL 5
# =============================================================================

def morris_sensitivity_level5(
    param_ranges,
    N_trajectories=10,
    num_levels=4,
    seed=42,
    output_metric='flux'
):
    """
    Morris sensitivity analysis for LEVEL 5 complete hierarchical model.
    
    Parameters
    ----------
    param_ranges : dict
        Parameter ranges as {param_name: [min, max]}.
        Use SUGGESTED_RANGES_LEVEL5 as starting point.
        
    N_trajectories : int
        Number of Morris trajectories (default: 10).
        Total samples = N_trajectories × (num_params + 1)
        
    num_levels : int
        Number of grid levels (default: 4)
        
    output_metric : str
        Output to analyze. Options:
        - 'flux': Total permeation flux [mol/m²/s]
        - 'PRF': Permeation Reduction Factor [-]
        - 'D_eff': Effective metal diffusivity [m²/s]
        - 'D_modification': D_eff / D_lattice [-]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
        - 'flux_intact': Flux through intact oxide [mol/m²/s]
        - 'flux_defect': Flux through defect paths [mol/m²/s]
        
    Returns
    -------
    Si : dict
        Morris sensitivity indices:
        - 'mu_star': Mean absolute elementary effect (importance)
        - 'sigma': Standard deviation (nonlinearity/interactions)
        - 'mu_star_conf': Confidence interval for mu_star
    problem : dict
        SALib problem definition
    Y : ndarray
        Model outputs for all samples
        
    Examples
    --------
    >>> # Quick screening of key parameters
    >>> param_ranges = {
    ...     'temperature': [873, 1273],
    ...     'defect_fraction': [0.001, 0.1],
    ...     'grain_size': [1e-6, 100e-6],
    ...     'trap_vacancy_N_T': [1e21, 1e24],
    ... }
    >>> Si, problem, Y = morris_sensitivity_level5(param_ranges, N_trajectories=10)
    >>> plot_morris_results(Si, problem, 'Flux (Level 5)')
    """
    # Validate output metric
    if output_metric not in VALID_OUTPUT_METRICS_L5:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS_L5}, got '{output_metric}'")
    
    # Define problem for SALib
    problem = {
        'num_vars': len(param_ranges),
        'names': list(param_ranges.keys()),
        'bounds': list(param_ranges.values())
    }
    
    # Generate Morris samples
    param_values = morris_sampler.sample(
        problem, 
        N=N_trajectories, 
        num_levels=num_levels,
        seed=seed
    )
    
    n_samples = param_values.shape[0]
    Y = np.zeros(n_samples)
    
    print(f"\n{'='*70}")
    print(f"MORRIS SENSITIVITY ANALYSIS - LEVEL 5 (Complete Hierarchical Model)")
    print(f"{'='*70}")
    print(f"Running {n_samples} Morris samples...")
    print(f"Parameters ({len(problem['names'])}): {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print(f"{'='*70}\n")
    
    # Run model for each sample
    for i, X in enumerate(param_values):
        sample_params = dict(zip(problem['names'], X))
        result = level5_model_wrapper(sample_params)
        Y[i] = result[output_metric]
        
        if (i + 1) % 10 == 0 or (i + 1) == n_samples:
            print(f"  Completed {i + 1}/{n_samples} samples")
    
    # Handle invalid values
    valid_idx = np.isfinite(Y) & (Y > 0)
    n_invalid = np.sum(~valid_idx)
    
    if n_invalid > 0:
        print(f"\n⚠️  Warning: {n_invalid} invalid outputs detected")
        if np.any(valid_idx):
            median_val = np.median(Y[valid_idx])
            Y[~valid_idx] = median_val
            print(f"   Replaced with median: {median_val:.2e}")
        else:
            print("   All outputs invalid! Check parameter ranges.")
            Y[:] = 1e-15
    
    # Analyze with Morris method
    Si = morris_analyzer.analyze(
        problem, 
        param_values, 
        Y, 
        conf_level=0.95,
        print_to_console=False
    )
    
    print(f"\n✓ Morris analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    print(f"  Output span: {np.max(Y)/np.min(Y):.1f}× variation")
    
    # Print quick summary
    mu_star = Si['mu_star']
    sorted_idx = np.argsort(mu_star)[::-1]
    print(f"\n  Top 10 most influential parameters:")
    for rank, idx in enumerate(sorted_idx[:10]):
        print(f"    {rank+1}. {problem['names'][idx]}: μ* = {mu_star[idx]:.2e}")
    
    return Si, problem, Y   



# =============================================================================
# SOBOL SENSITIVITY ANALYSIS - LEVEL 5
# =============================================================================

def sobol_sensitivity_level5(
    param_ranges,
    N_samples=1024,
    output_metric='flux',
    seed=42,
    calc_second_order=True
):
    """
    Sobol variance-based sensitivity analysis for LEVEL 5 complete hierarchical model.
    
    Provides quantitative decomposition of output variance into contributions
    from each parameter and their interactions.
    
    Parameters
    ----------
    param_ranges : dict
        Parameter ranges as {param_name: [min, max]}.
        Use top parameters from Morris screening to reduce computation.
        
    N_samples : int
        Base number of samples (must be power of 2: 256, 512, 1024, 2048).
        Total model evaluations = N_samples × (2×num_params + 2) for second-order
        or N_samples × (num_params + 2) for first-order only.
        
    output_metric : str
        Output to analyze. Options:
        - 'flux': Total permeation flux [mol/m²/s]
        - 'PRF': Permeation Reduction Factor [-]
        - 'D_eff': Effective metal diffusivity [m²/s]
        - 'D_modification': D_eff / D_lattice [-]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
        - 'flux_intact': Flux through intact oxide [mol/m²/s]
        - 'flux_defect': Flux through defect paths [mol/m²/s]
        
    calc_second_order : bool
        If True, calculate pairwise interaction indices (S2).
        More expensive but reveals parameter interactions.
        
    Returns
    -------
    Si : dict
        Sobol sensitivity indices:
        - 'S1': First-order indices (direct effect)
        - 'S1_conf': Confidence intervals for S1
        - 'ST': Total-order indices (direct + all interactions)
        - 'ST_conf': Confidence intervals for ST
        - 'S2': Second-order indices (pairwise interactions), if calc_second_order=True
    problem : dict
        SALib problem definition
    Y : ndarray
        Model outputs for all samples
        
    Notes
    -----
    Interpretation:
    - S1 ≈ ST: Parameter acts independently (no interactions)
    - ST >> S1: Strong interactions with other parameters
    - Sum(S1) ≈ 1: Model is additive
    - Sum(ST) > 1: Significant interactions exist
    
    Computational cost:
    - 5 params, N=1024, second_order=True: ~12,000 evaluations
    - 10 params, N=1024, second_order=True: ~22,000 evaluations
    - Consider using Morris first to screen down to 5-8 key parameters
        
    Examples
    --------
    >>> # After Morris screening, analyze top parameters
    >>> param_ranges = {
    ...     'temperature': [873, 1273],
    ...     'defect_fraction': [0.001, 0.1],
    ...     'oxide_thickness': [1e-7, 1e-5],
    ...     'trap_vacancy_N_T': [1e21, 1e24],
    ... }
    >>> Si, problem, Y = sobol_sensitivity_level5(param_ranges, N_samples=1024)
    >>> plot_sobol_results(Si, problem, 'Flux (Level 5)')
    """
    # Validate output metric
    if output_metric not in VALID_OUTPUT_METRICS_L5:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS_L5}, got '{output_metric}'")
    
    # Define problem for SALib
    problem = {
        'num_vars': len(param_ranges),
        'names': list(param_ranges.keys()),
        'bounds': list(param_ranges.values())
    }
    
    # Set random seed for reproducibility (older SALib versions don't support seed in sampler)
    if seed is not None:
        np.random.seed(seed)
    
    # Generate Sobol samples using Saltelli's extension
    param_values = sobol_sampler.sample(
        problem, 
        N_samples,
        #seed=seed,
        calc_second_order=calc_second_order
        
    )
    
    n_samples = param_values.shape[0]
    Y = np.zeros(n_samples)
    
    # Estimate time
    n_params = len(param_ranges)
    if calc_second_order:
        expected_samples = N_samples * (2 * n_params + 2)
    else:
        expected_samples = N_samples * (n_params + 2)
    
    print(f"\n{'='*70}")
    print(f"SOBOL SENSITIVITY ANALYSIS - LEVEL 5 (Complete Hierarchical Model)")
    print(f"{'='*70}")
    print(f"Running {n_samples} Sobol samples (N_base={N_samples}, {n_params} params)")
    print(f"Parameters: {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print(f"Second-order interactions: {'Yes' if calc_second_order else 'No'}")
    print(f"{'='*70}")
    print(f"⏳ This may take a while...\n")
    
    # Run model for each sample with progress updates
    for i, X in enumerate(param_values):
        sample_params = dict(zip(problem['names'], X))
        result = level5_model_wrapper(sample_params)
        Y[i] = result[output_metric]
        
        # Progress updates at 10%, 25%, 50%, 75%, 100%
        progress = (i + 1) / n_samples
        if (i + 1) % max(1, n_samples // 10) == 0 or (i + 1) == n_samples:
            print(f"  Completed {i + 1}/{n_samples} samples ({progress*100:.0f}%)")
    
    # Handle invalid values
    valid_idx = np.isfinite(Y) & (Y > 0)
    n_invalid = np.sum(~valid_idx)
    
    if n_invalid > 0:
        print(f"\n⚠️  Warning: {n_invalid} invalid outputs detected ({n_invalid/n_samples*100:.1f}%)")
        if np.any(valid_idx):
            median_val = np.median(Y[valid_idx])
            Y[~valid_idx] = median_val
            print(f"   Replaced with median: {median_val:.2e}")
        else:
            print("   All outputs invalid! Check parameter ranges.")
            Y[:] = 1e-15
    
    # Analyze with Sobol method
    Si = sobol_analyzer.analyze(
        problem, 
        Y,
        calc_second_order=calc_second_order,
        conf_level=0.95,
        print_to_console=False
    )
    
    print(f"\n✓ Sobol analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    print(f"  Output span: {np.max(Y)/np.min(Y):.1f}× variation")
    
    # Print summary table
    print(f"\n  {'Parameter':<25} {'S1':>8} {'ST':>8} {'ST-S1':>8}")
    print(f"  {'-'*25} {'-'*8} {'-'*8} {'-'*8}")
    
    S1 = Si['S1']
    ST = Si['ST']
    sorted_idx = np.argsort(ST)[::-1]
    
    for idx in sorted_idx:
        name = problem['names'][idx]
        s1 = S1[idx]
        st = ST[idx]
        interaction = st - s1
        print(f"  {name:<25} {s1:>8.3f} {st:>8.3f} {interaction:>8.3f}")
    
    # Interpretation
    print(f"\n  Sum(S1) = {np.sum(S1):.3f} (should be ~1 if additive)")
    print(f"  Sum(ST) = {np.sum(ST):.3f} (>1 indicates interactions)")
    
    # Top interactions if S2 available
    if calc_second_order and 'S2' in Si and Si['S2'] is not None:
        S2 = Si['S2']
        print(f"\n  Top parameter interactions (S2 > 0.01):")
        interactions_found = False
        for i in range(n_params):
            for j in range(i+1, n_params):
                if S2[i, j] > 0.01:
                    interactions_found = True
                    print(f"    {problem['names'][i]} × {problem['names'][j]}: S2 = {S2[i,j]:.3f}")
        if not interactions_found:
            print(f"    No significant pairwise interactions detected")
    
    return Si, problem, Y