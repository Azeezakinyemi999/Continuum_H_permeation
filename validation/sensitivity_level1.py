import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from SALib.sample import morris as morris_sampler
from SALib.analyze import morris as morris_analyzer
from SALib.sample import saltelli as sobol_sampler
from SALib.analyze import sobol as sobol_analyzer

# =============================================================================
# LEVEL 5: Complete Hierarchical Model (L3 + L4)
# Defective oxide (L3) + defective metal microstructure (L4)
#
# Reference-point Arrhenius: X(T) = X_ref * exp(-E/R * (1/T - 1/T_ref))
# Reference values from model_config.py:
#   Incoloy 802 (X40 NiCrAlTi) — T_ref_metal  = 1223 K  (Schmidt 1985)
#   Cr2O3_sample4               — T_ref_oxide  =  673 K  (Nemanic 2023 / Stover 1986)
# =============================================================================

DEFAULT_PARAMS_LEVEL5 = {
    # -------------------------------------------------------------------------
    # Metal transport  (Incoloy 802 at T_ref_metal = 1223 K, Schmidt 1985)
    # -------------------------------------------------------------------------
    'D_ref':       4.285e-9,    # m²/s           D_metal at T_ref_metal
    'E_D':         72200,       # J/mol           activation energy for diffusion
    'K_s_ref':     0.05093,     # mol/m³/Pa^0.5   K_s at T_ref_metal
    'H_s':         7200,        # J/mol           heat of solution
    'T_ref_metal': 1223,        # K               Incoloy 802 measurement reference

    # -------------------------------------------------------------------------
    # Oxide transport  (Cr2O3_sample4 at T_ref_oxide = 673 K)
    # -------------------------------------------------------------------------
    'D_ox_ref':    7.800e-19,   # m²/s            D_ox at T_ref_oxide
    'E_D_ox':      70434,       # J/mol
    'K_ox_ref':    0.35417,     # mol/m³/Pa^0.5
    'H_sol_ox':    163566,      # J/mol
    'T_ref_oxide': 673,         # K               Cr2O3 measurement reference

    # -------------------------------------------------------------------------
    # Geometry
    # -------------------------------------------------------------------------
    'metal_thickness': 1e-3,    # m (1 mm)
    'oxide_thickness': 4.8e-8,  # m (48 nm, Cr2O3_sample4)

    # -------------------------------------------------------------------------
    # Operating conditions
    # -------------------------------------------------------------------------
    'P_upstream':   1e5,        # Pa (100 kPa)
    'P_downstream': 0.0,        # Pa
    'temperature':  1073,       # K (800°C)

    # -------------------------------------------------------------------------
    # Level 3: Oxide defect properties  (structure mirrors model_config OXIDE_DEFECTS)
    # -------------------------------------------------------------------------
    'f_pinhole':              0.001,    # area fraction — pinholes (config: 0.001)
    'f_crack':                0.001,    # area fraction — cracks   (config: 0.001)
    'f_gb_defect':            0.001,    # area fraction — GB paths (config: 0.001)
    'crack_thickness_factor': 0.1,      # L_crack = factor × L_oxide
    'gb_diffusivity_factor':  10.0,     # D_gb_defect = factor × D_oxide
    'use_sieverts_pinhole':   False,    # False = use metal surface kinetics at pinholes

    # -------------------------------------------------------------------------
    # Level 4: Grain structure  (Zhu 2021 for Hastelloy N baseline)
    # -------------------------------------------------------------------------
    'grain_size':            100e-6,    # m (100 μm, Zhu 2021 midpoint)
    'grain_shape':           'equiaxed',
    'gb_type':               'LAGB',
    'gb_thickness':          0.5e-9,    # m (0.5 nm)
    'gb_enhancement_factor': 100,
    'lattice_density':       8.774e28,  # m⁻³ (Fe BCC estimate — update for alloy)

    # -------------------------------------------------------------------------
    # Level 4: Trap properties  (Lu 2022 + Young 1997 via model_config MICROSTRUCTURE)
    # -------------------------------------------------------------------------
    'trap_dislocation_E_b': 19297,   # J/mol  0.20 eV (Lu 2022: 0.186–0.215 eV)
    'trap_dislocation_N_T': 8.16e12, # m⁻³   (Zhu 2021 GND density)
    'trap_gb_E_b':          26051,   # J/mol  0.27 eV (Lu 2022 Peak 2: 0.258–0.281 eV)
    'trap_gb_N_T':          6e14,    # m⁻³   (geometric, 100 μm grain)
    'trap_vacancy_E_b':     41489,   # J/mol  0.43 eV (Ni estimate)
    'trap_vacancy_N_T':     1e26,    # m⁻³
    'trap_carbide_E_b':     26051,   # J/mol  0.27 eV (M6C, Lu 2022)
    'trap_carbide_N_T':     2e25,    # m⁻³   (Young 1997 upper bound)

    # -------------------------------------------------------------------------
    # Model options
    # -------------------------------------------------------------------------
    'include_gb_enhancement': True,
    'include_trapping':       True,
    'D_eff_method':           'average',
}

# =============================================================================
# VALID OUTPUT METRICS
# =============================================================================
VALID_OUTPUT_METRICS_L5 = [
    'flux',        # Total permeation flux [mol/m²/s]
    'permeability', # Effective permeability [mol/m/s/Pa^0.5]
]

# =============================================================================
# PARAMETER GROUPS (for organised sensitivity analysis)
# =============================================================================
PARAM_GROUPS = {
    'metal_transport': ['D_ref', 'E_D', 'K_s_ref', 'H_s', 'T_ref_metal', 'metal_thickness'],
    'oxide_transport': ['D_ox_ref', 'E_D_ox', 'K_ox_ref', 'H_sol_ox', 'T_ref_oxide', 'oxide_thickness'],
    'oxide_defects':   ['f_pinhole', 'f_crack', 'f_gb_defect',
                        'crack_thickness_factor', 'gb_diffusivity_factor'],
    'microstructure':  ['grain_size', 'gb_thickness', 'lattice_density'],
    'traps':           ['trap_dislocation_E_b', 'trap_dislocation_N_T',
                        'trap_gb_E_b', 'trap_gb_N_T', 'trap_vacancy_E_b', 'trap_vacancy_N_T',
                        'trap_carbide_E_b', 'trap_carbide_N_T'],
    'operating':       ['temperature', 'P_upstream'],
}

# =============================================================================
# SUGGESTED PARAMETER RANGES FOR SENSITIVITY ANALYSIS  (28 parameters)
# Ref-point values span ±2 orders of magnitude around the reference value.
# Activation energies and T_ref span physically meaningful bounds.
# =============================================================================
SUGGESTED_RANGES_LEVEL5 = {
    # Metal transport (Incoloy 802, Schmidt 1985)
    'D_ref':        [4.285e-11, 4.285e-7],
    'E_D':          [60000, 80000],
    'K_s_ref':      [1e-4, 0.1],
    'H_s':          [1000, 50000],
    'T_ref_metal':  [600, 1400],

    # Oxide transport (Cr2O3_sample4)
    'D_ox_ref':     [7.800e-21, 7.800e-17],
    'E_D_ox':       [69000, 71000],
    'K_ox_ref':     [0.08, 0.4],
    'H_sol_ox':     [20000, 170000],
    'T_ref_oxide':  [600, 1400],

    # Geometry
    'metal_thickness': [5e-4, 5e-3],
    'oxide_thickness': [1e-10, 1e-5],
    'P_upstream':      [1e-7, 1e7],

    # Oxide defects
    'f_pinhole':              [1e-5, 0.1],
    'f_crack':                [1e-5, 0.1],
    'f_gb_defect':            [1e-5, 0.1],
    'crack_thickness_factor': [0.01, 0.5],
    'gb_diffusivity_factor':  [1.0, 100.0],

    # Grain structure
    'grain_size':      [1e-5, 5e-4],
    'gb_thickness':    [3e-10, 1e-9],
    'lattice_density': [8e28, 1.2e29],

    # Traps (Lu 2022 + Young 1997; ranges narrowed to physically motivated bounds)
    'trap_dislocation_E_b': [15000, 25000],
    'trap_dislocation_N_T': [8.16e10, 8.16e14],
    'trap_gb_E_b':          [25000, 30000],
    'trap_gb_N_T':          [6e12, 6e16],
    'trap_vacancy_E_b':     [35000, 50000],
    'trap_vacancy_N_T':     [1e24, 1e28],
    'trap_carbide_E_b':     [20000, 40000],
    'trap_carbide_N_T':     [2e23, 2e27],

    # Operating conditions
    'temperature': [573, 1273],
}


# =============================================================================
# MODEL WRAPPER — LEVEL 5
# =============================================================================

def level5_model_wrapper(params_dict):
    """
    Wrapper for LEVEL 5 complete hierarchical model (L3 + L4).

    Combines defective oxide (L3) and defective metal microstructure (L4).
    Uses reference-point Arrhenius: X(T) = X_ref * exp(-E/R * (1/T - 1/T_ref)).

    Parameters
    ----------
    params_dict : dict
        Any subset of parameters from DEFAULT_PARAMS_LEVEL5.
        Missing parameters are filled from defaults.

    Returns
    -------
    dict
        - 'flux': Total permeation flux [mol/m²/s]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
        - 'PRF': Permeation Reduction Factor [-]
        - 'D_eff', 'D_modification', 'P_interface'
        - 'flux_intact', 'flux_defect', 'flux_bare_metal'
    """
    from calculations.parallel_oxide_defect_paths import (
        calculate_parallel_path_flux_defective_metal,
    )
    from calculations.permeation_calc import calculate_defective_metal_flux

    full_params = DEFAULT_PARAMS_LEVEL5.copy()
    full_params.update(params_dict)

    try:
        R = 8.314  # J/mol/K
        T = full_params['temperature']

        # Reference-point Arrhenius
        T_ref_m  = full_params['T_ref_metal']
        T_ref_ox = full_params['T_ref_oxide']
        inv_m    = 1.0 / T - 1.0 / T_ref_m
        inv_ox   = 1.0 / T - 1.0 / T_ref_ox

        D_metal   = full_params['D_ref']    * np.exp((-full_params['E_D']      / R) * inv_m)
        K_s_metal = full_params['K_s_ref']  * np.exp((-full_params['H_s']      / R) * inv_m)
        D_ox      = full_params['D_ox_ref'] * np.exp((-full_params['E_D_ox']   / R) * inv_ox)
        K_ox      = full_params['K_ox_ref'] * np.exp((-full_params['H_sol_ox'] / R) * inv_ox)

        oxide_props = {
            'D_ox':      D_ox,
            'K_ox':      K_ox,
            'thickness': full_params['oxide_thickness'],
        }
        metal_props = {
            'D_metal':   D_metal,
            'K_s_metal': K_s_metal,
            'thickness': full_params['metal_thickness'],
        }

        # Build trap list
        trap_list = []
        if full_params.get('trap_dislocation_N_T', 0) > 0:
            trap_list.append({
                'name': 'dislocations',
                'binding_energy': full_params['trap_dislocation_E_b'],
                'density': full_params['trap_dislocation_N_T'],
            })
        if full_params.get('trap_vacancy_N_T', 0) > 0:
            trap_list.append({
                'name': 'vacancies',
                'binding_energy': full_params['trap_vacancy_E_b'],
                'density': full_params['trap_vacancy_N_T'],
            })

        grain_size   = full_params['grain_size']
        gb_thickness = full_params['gb_thickness']
        N_T_gb       = full_params.get('trap_gb_N_T', 0)
        if full_params.get('trap_gb_E_b', 0) > 0 and N_T_gb > 0:
            trap_list.append({
                'name': 'grain_boundaries',
                'binding_energy': full_params['trap_gb_E_b'],
                'density': N_T_gb,
            })
        N_T_carbide = full_params.get('trap_carbide_N_T', 0)
        if full_params.get('trap_carbide_E_b', 0) > 0 and N_T_carbide > 0:
            trap_list.append({
                'name': 'carbides',
                'binding_energy': full_params['trap_carbide_E_b'],
                'density': N_T_carbide,
            })

        microstructure_params = {
            'grain_size':            grain_size,
            'grain_shape':           full_params['grain_shape'],
            'gb_type':               full_params['gb_type'],
            'gb_thickness':          gb_thickness,
            'trap_list':             trap_list,
            'gb_enhancement_factor': full_params.get('gb_enhancement_factor', 100),
        }

        f_pin = full_params.get('f_pinhole', 0.0)
        f_cra = full_params.get('f_crack',   0.0)
        f_gb  = full_params.get('f_gb_defect', 0.0)
        defect_params = {
            'area_fraction': f_pin + f_cra + f_gb,
            'type': 'mixed',
            'components': {
                'pinholes':         f_pin,
                'cracks':           f_cra,
                'grain_boundaries': f_gb,
            },
            'thickness_factor':   full_params.get('crack_thickness_factor', 0.1),
            'diffusivity_factor': full_params.get('gb_diffusivity_factor', 10.0),
        }

        P_upstream      = full_params['P_upstream']
        P_downstream    = full_params['P_downstream']
        lattice_density = full_params['lattice_density']
        method          = full_params.get('D_eff_method', 'average')

        include_gb   = full_params.get('include_gb_enhancement', True)
        include_trap = full_params.get('include_trapping', True)
        if include_gb and include_trap:
            mode = 'both'
        elif include_gb:
            mode = 'gb_only'
        elif include_trap:
            mode = 'trapping_only'
        else:
            mode = 'none'

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
            mode=mode,
        )

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
            mode=mode,
        )

        flux_total = result_l5['flux_total']
        flux_bare  = result_bare['flux']
        PRF        = flux_bare / flux_total if flux_total > 0 else float('inf')
        D_eff      = result_l5.get('D_eff_metal', D_metal)

        # Effective permeability via harmonic mean (series resistances)
        # 1/Φ_eff = 1/Φ_oxide + 1/Φ_metal (no length dependence)
        Phi_oxide = D_ox * K_ox
        Phi_metal = D_eff * K_s_metal

        if Phi_oxide > 0 and Phi_metal > 0:
            permeability = 1.0 / (1.0/Phi_oxide + 1.0/Phi_metal)
        else:
            permeability = np.nan

        return {
            'flux':         flux_total,
            'PRF':          PRF,
            'D_eff':        D_eff,
            'D_modification': D_eff / D_metal if D_metal > 0 else 1.0,
            'permeability': permeability,
            'P_interface':  result_l5.get('P_interface_intact', 0),
            'flux_intact':  result_l5['flux_intact_contribution'],
            'flux_defect':  result_l5['flux_defect_contribution'],
            'flux_bare_metal': flux_bare,
            'regime':       result_l5.get('regime', 'unknown'),
            'dominant_path': result_l5.get('dominant_path', 'unknown'),
            'D_metal':      D_metal,
            'K_s_metal':    K_s_metal,
            'D_ox':         D_ox,
            'K_ox':         K_ox,
            'temperature':  T,
            'modification_factor': result_l5.get('modification_factor', 1.0),
            'defect_enhancement': result_l5.get('defect_enhancement_factor', 1.0),
        }

    except Exception as e:
        print(f"Error in Level 5 model: {e}")
        import traceback; traceback.print_exc()
        return {
            'flux': 1e-20, 'PRF': 1.0, 'D_eff': 1e-12, 'D_modification': 1.0,
            'permeability': 1e-20, 'P_interface': 0,
            'flux_intact': 1e-20, 'flux_defect': 0, 'flux_bare_metal': 1e-20,
            'regime': 'error', 'dominant_path': 'error',
            'D_metal': 1e-12, 'K_s_metal': 1e-6, 'D_ox': 1e-15, 'K_ox': 1e-10,
            'temperature': full_params.get('temperature', 1073.15),
            'modification_factor': 1.0, 'defect_enhancement': 1.0,
        }


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_morris_results(Si, problem, output_metric='Model Output'):
    """Visualize Morris sensitivity results."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    param_names = problem['names']
    mu    = np.nan_to_num(Si['mu'], nan=0.0)
    sigma = np.nan_to_num(Si['sigma'], nan=0.0)

    # Sort by |μ|
    sorted_idx   = np.argsort(np.abs(mu))[::-1]
    sorted_names = [param_names[i] for i in sorted_idx]
    sorted_mu    = mu[sorted_idx]

    # Bar chart — signed μ, colour shows direction
    colors = ['steelblue' if v >= 0 else 'tomato' for v in sorted_mu]
    ax1.barh(sorted_names, sorted_mu, color=colors, alpha=0.7)
    ax1.axvline(0, color='black', linewidth=0.8)
    ax1.set_xlabel('μ (Mean Elementary Effect)', fontsize=12)
    ax1.set_title(f'Morris Sensitivity Analysis\n{output_metric}', fontsize=14, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)

    # μ vs σ scatter (using |μ| for normalized range starting at 0)
    mu_abs = np.abs(mu)
    ax2.scatter(mu_abs, sigma, s=150, alpha=0.6, c='steelblue', edgecolors='black')
    for i, name in enumerate(param_names):
        ax2.annotate(name, (mu_abs[i], sigma[i]),
                     fontsize=10, ha='right', va='bottom',
                     xytext=(-5, 5), textcoords='offset points')
    ax2.set_xlabel('|μ| (Mean Elementary Effect — magnitude)', fontsize=12)
    ax2.set_ylabel('σ (Nonlinearity/Interactions)', fontsize=12)
    ax2.set_title('Elementary Effects', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    print(f"\n{'='*70}")
    print(f"MORRIS SENSITIVITY RESULTS - {output_metric}")
    print(f"{'='*70}")
    df = pd.DataFrame({
        'Parameter': sorted_names,
        'μ':      sorted_mu,
        'σ':      [sigma[i] for i in sorted_idx],
        'μ_conf': [Si['mu_star_conf'][i] for i in sorted_idx],
    })
    print(df.to_string(index=False))
    print(f"{'='*70}")
    print("\nInterpretation:")
    print("  μ > 0 : Parameter increases output on average")
    print("  μ < 0 : Parameter decreases output on average")
    print("  σ     : Nonlinearity/Interactions (higher = more complex)")
    print(f"{'='*70}\n")


def plot_sobol_results(Si, problem, output_metric='Model Output'):
    """Visualize Sobol sensitivity results."""
    param_names = problem['names']
    S1 = Si['S1']
    ST = Si['ST']

    has_S2 = 'S2' in Si and Si['S2'] is not None
    n_plots = 2 if has_S2 else 1

    fig, axes = plt.subplots(1, n_plots, figsize=(7 * n_plots, 6))
    if n_plots == 1:
        axes = [axes]

    ax1 = axes[0]
    x     = np.arange(len(param_names))
    width = 0.35
    ax1.bar(x - width / 2, S1, width, label='S1 (First-order)',  alpha=0.8, color='steelblue')
    ax1.bar(x + width / 2, ST, width, label='ST (Total-order)',  alpha=0.8, color='coral')
    ax1.set_xlabel('Parameters', fontsize=12)
    ax1.set_ylabel('Sobol Index', fontsize=12)
    ax1.set_title(f'Sobol Sensitivity\n{output_metric}', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(param_names, rotation=45, ha='right')
    ax1.legend(fontsize=11)
    ax1.grid(axis='y', alpha=0.3)

    if has_S2:
        ax2 = axes[1]
        S2  = Si['S2']
        im  = ax2.imshow(S2, cmap='RdYlGn', aspect='auto', vmin=0, vmax=np.max(S2))
        ax2.set_xticks(np.arange(len(param_names)))
        ax2.set_yticks(np.arange(len(param_names)))
        ax2.set_xticklabels(param_names, rotation=45, ha='right')
        ax2.set_yticklabels(param_names)
        cbar = plt.colorbar(im, ax=ax2)
        cbar.set_label('S2', fontsize=11)
        for i in range(len(param_names)):
            for j in range(len(param_names)):
                if i != j and S2[i, j] > 0.01:
                    ax2.text(j, i, f'{S2[i, j]:.3f}', ha='center', va='center',
                             color='black', fontsize=9)
        ax2.set_title('Second-Order Interactions', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.show()

    print(f"\n{'='*70}")
    print(f"SOBOL SENSITIVITY RESULTS - {output_metric}")
    print(f"{'='*70}")
    df = pd.DataFrame({
        'Parameter': param_names,
        'S1': S1,
        'S1_conf': Si['S1_conf'],
        'ST': ST,
        'ST_conf': Si['ST_conf'],
    })
    df['Interactions'] = df['ST'] - df['S1']
    print(df.to_string(index=False))
    print(f"{'='*70}")
    print("\nInterpretation:")
    print("  S1   : First-order effect (direct influence)")
    print("  ST   : Total-order effect (direct + interactions)")
    print("  ST-S1: Pure interaction effects")
    print(f"{'='*70}\n")
    important = df.nlargest(3, 'ST')
    print("Top 3 Most Influential Parameters:")
    for _, row in important.iterrows():
        print(f"  {row['Parameter']}: ST = {row['ST']:.3f} (S1 = {row['S1']:.3f})")
    print(f"{'='*70}\n")


# =============================================================================
# MORRIS SENSITIVITY ANALYSIS — LEVEL 5
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
    num_levels : int
        Number of grid levels (default: 4).
    output_metric : str
        One of VALID_OUTPUT_METRICS_L5.

    Returns
    -------
    Si : dict
        Morris sensitivity indices (mu, mu_star, sigma, mu_star_conf).
    problem : dict
        SALib problem definition.
    Y : ndarray
        Model outputs for all samples.
    """
    if output_metric not in VALID_OUTPUT_METRICS_L5:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS_L5}, got '{output_metric}'")

    problem = {
        'num_vars': len(param_ranges),
        'names':    list(param_ranges.keys()),
        'bounds':   list(param_ranges.values()),
    }

    param_values = morris_sampler.sample(problem, N=N_trajectories, num_levels=num_levels, seed=seed)
    n_samples = param_values.shape[0]
    Y = np.zeros(n_samples)

    print(f"\n{'='*70}")
    print(f"MORRIS SENSITIVITY ANALYSIS - LEVEL 5 (Complete Hierarchical Model)")
    print(f"{'='*70}")
    print(f"Running {n_samples} Morris samples...")
    print(f"Parameters ({len(problem['names'])}): {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print(f"{'='*70}\n")

    for i, X in enumerate(param_values):
        result = level5_model_wrapper(dict(zip(problem['names'], X)))
        Y[i] = result[output_metric]
        if (i + 1) % 10 == 0 or (i + 1) == n_samples:
            print(f"  Completed {i + 1}/{n_samples} samples")

    valid_idx = np.isfinite(Y) & (Y > 0)
    n_invalid = np.sum(~valid_idx)
    if n_invalid > 0:
        print(f"\n⚠️  Warning: {n_invalid} invalid outputs detected")
        Y[~valid_idx] = np.median(Y[valid_idx]) if np.any(valid_idx) else 1e-15

    Si = morris_analyzer.analyze(problem, param_values, Y, conf_level=0.95, print_to_console=False)

    print(f"\n✓ Morris analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    print(f"  Output span: {np.max(Y)/np.min(Y):.1f}× variation")

    mu = Si['mu']
    sorted_idx = np.argsort(np.abs(mu))[::-1]
    print(f"\n  Top 10 most influential parameters:")
    for rank, idx in enumerate(sorted_idx[:10]):
        print(f"    {rank+1}. {problem['names'][idx]}: μ = {mu[idx]:.2e}")

    return Si, problem, Y


# =============================================================================
# SOBOL SENSITIVITY ANALYSIS — LEVEL 5
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

    Parameters
    ----------
    param_ranges : dict
        Parameter ranges as {param_name: [min, max]}.
    N_samples : int
        Base number of samples (power of 2: 256, 512, 1024, 2048).
    output_metric : str
        One of VALID_OUTPUT_METRICS_L5.
    calc_second_order : bool
        Whether to compute pairwise S2 indices.

    Returns
    -------
    Si, problem, Y
    """
    if output_metric not in VALID_OUTPUT_METRICS_L5:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS_L5}, got '{output_metric}'")

    problem = {
        'num_vars': len(param_ranges),
        'names':    list(param_ranges.keys()),
        'bounds':   list(param_ranges.values()),
    }

    if seed is not None:
        np.random.seed(seed)

    param_values = sobol_sampler.sample(problem, N_samples, calc_second_order=calc_second_order)
    n_samples = param_values.shape[0]
    n_params  = len(param_ranges)
    Y = np.zeros(n_samples)

    print(f"\n{'='*70}")
    print(f"SOBOL SENSITIVITY ANALYSIS - LEVEL 5 (Complete Hierarchical Model)")
    print(f"{'='*70}")
    print(f"Running {n_samples} Sobol samples (N_base={N_samples}, {n_params} params)")
    print(f"Parameters: {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print(f"Second-order interactions: {'Yes' if calc_second_order else 'No'}")
    print(f"{'='*70}")
    print(f"⏳ This may take a while...\n")

    for i, X in enumerate(param_values):
        result = level5_model_wrapper(dict(zip(problem['names'], X)))
        Y[i] = result[output_metric]
        if (i + 1) % max(1, n_samples // 10) == 0 or (i + 1) == n_samples:
            print(f"  Completed {i + 1}/{n_samples} samples ({(i+1)/n_samples*100:.0f}%)")

    valid_idx = np.isfinite(Y) & (Y > 0)
    n_invalid = np.sum(~valid_idx)
    if n_invalid > 0:
        print(f"\n⚠️  Warning: {n_invalid} invalid outputs ({n_invalid/n_samples*100:.1f}%)")
        Y[~valid_idx] = np.median(Y[valid_idx]) if np.any(valid_idx) else 1e-15

    Si = sobol_analyzer.analyze(problem, Y, calc_second_order=calc_second_order,
                                conf_level=0.95, print_to_console=False)

    print(f"\n✓ Sobol analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    print(f"  Output span: {np.max(Y)/np.min(Y):.1f}× variation")

    S1, ST = Si['S1'], Si['ST']
    sorted_idx = np.argsort(ST)[::-1]
    print(f"\n  {'Parameter':<25} {'S1':>8} {'ST':>8} {'ST-S1':>8}")
    print(f"  {'-'*25} {'-'*8} {'-'*8} {'-'*8}")
    for idx in sorted_idx:
        name       = problem['names'][idx]
        interaction = ST[idx] - S1[idx]
        print(f"  {name:<25} {S1[idx]:>8.3f} {ST[idx]:>8.3f} {interaction:>8.3f}")
    print(f"\n  Sum(S1) = {np.sum(S1):.3f} (should be ~1 if additive)")
    print(f"  Sum(ST) = {np.sum(ST):.3f} (>1 indicates interactions)")

    if calc_second_order and 'S2' in Si and Si['S2'] is not None:
        S2 = Si['S2']
        print(f"\n  Top parameter interactions (S2 > 0.01):")
        found = False
        for i in range(n_params):
            for j in range(i + 1, n_params):
                if S2[i, j] > 0.01:
                    found = True
                    print(f"    {problem['names'][i]} × {problem['names'][j]}: S2 = {S2[i,j]:.3f}")
        if not found:
            print(f"    No significant pairwise interactions detected")

    return Si, problem, Y


# =============================================================================
# LEVEL 5L6: SURFACE KINETICS EXTENSION
# Adds oxide + metal surface dissociation/recombination kinetics.
# Maps to Surface_proposal.ipynb → calculate_full_model_flux_L346_v2.
# Cleanly inherits all Level 5 parameters via **DEFAULT_PARAMS_LEVEL5.
# =============================================================================

DEFAULT_PARAMS_LEVEL5L6 = {
    **DEFAULT_PARAMS_LEVEL5,
    # -------------------------------------------------------------------------
    # Oxide surface kinetics  (Cr2O3_sample4, T_ref_surface = 1623 K)
    # -------------------------------------------------------------------------
    'k_diss_ref':    9.487e-8,  # mol/m²/s/Pa   Grant 1988 oxidized surface condition
    'E_diss':        57950,     # J/mol          Cr2O3 surface activation energy
    'K_eq_ref':      1e-10,     # Pa⁻¹          PLACEHOLDER
    'H_eq':          20000,     # J/mol          PLACEHOLDER
    'T_ref_surface': 1623,      # K              oxide surface kinetics reference
    # -------------------------------------------------------------------------
    # Metal surface kinetics (Grant 1988 clean surface, T_ref = 965 K)
    # -------------------------------------------------------------------------
    'k_diss_metal_ref':  1.346e-6,  # mol/m²/s/Pa   Grant 1988 clean surface
    'E_diss_metal':      81560,     # J/mol
    'K_eq_metal_ref':    1e-8,      # Pa⁻¹          PLACEHOLDER
    'H_eq_metal':        15000,     # J/mol          PLACEHOLDER
    'T_ref_surface_metal': 965,     # K              Grant 1988 reference temperature
}

VALID_OUTPUT_METRICS_L5L6 = [m for m in VALID_OUTPUT_METRICS_L5 if m != 'PRF'] + [
    'frac_surface',  # Flux-weighted surface resistance fraction
    'frac_oxide',    # Flux-weighted oxide resistance fraction
    'frac_metal',    # Flux-weighted metal resistance fraction
    'theta',         # Surface coverage on intact path at steady state
]

SUGGESTED_RANGES_LEVEL5L6 = {
    **SUGGESTED_RANGES_LEVEL5,
    # Oxide surface kinetics (Cr2O3)
    'k_diss_ref':          [9.487e-10, 9.487e-6],
    'E_diss':              [40000, 70000],
    'K_eq_ref':            [1e-14, 1e-6],
    'H_eq':                [5000, 50000],
    'T_ref_surface':       [600, 1800],
    # Metal surface kinetics (Grant 1988)
    'k_diss_metal_ref':    [1.346e-8, 1.346e-4],
    'E_diss_metal':        [70000, 100000],
    'K_eq_metal_ref':      [1e-10, 1e-3],
    'H_eq_metal':          [5000, 40000],
    'T_ref_surface_metal': [600, 1100],
}


# =============================================================================
# MODEL WRAPPER — LEVEL 5L6
# =============================================================================

def level5L6_model_wrapper(params_dict):
    """
    Wrapper for LEVEL 5+L6 model: defective oxide + defective metal + surface kinetics.

    Extends Level 5 by adding Langmuir-Hinshelwood surface kinetics at both the
    gas-oxide interface (k_diss, K_eq) and the metal surface at pinholes
    (k_diss_metal, K_eq_metal). Calls calculate_full_model_flux_L346_v2.

    Parameters
    ----------
    params_dict : dict
        Any subset of parameters from DEFAULT_PARAMS_LEVEL5L6.

    Returns
    -------
    dict
        All Level 5 outputs plus:
        - 'frac_surface': Flux-weighted surface resistance fraction
        - 'frac_oxide':   Flux-weighted oxide resistance fraction
        - 'frac_metal':   Flux-weighted metal resistance fraction
        - 'theta':        Surface coverage on intact path
    """
    from calculations.surface_kinetics import calculate_full_model_flux_L346_v2

    full_params = DEFAULT_PARAMS_LEVEL5L6.copy()
    full_params.update(params_dict)

    try:
        R_gas = 8.314
        T = full_params['temperature']

        # Separate reference temperatures per property group
        T_ref_m  = full_params['T_ref_metal']          # metal transport
        T_ref_ox = full_params['T_ref_oxide']          # oxide transport
        T_ref_s  = full_params['T_ref_surface']        # oxide surface kinetics
        T_ref_sm = full_params['T_ref_surface_metal']  # metal surface kinetics (Grant 1988)

        inv_m  = 1.0 / T - 1.0 / T_ref_m
        inv_ox = 1.0 / T - 1.0 / T_ref_ox
        inv_s  = 1.0 / T - 1.0 / T_ref_s
        inv_sm = 1.0 / T - 1.0 / T_ref_sm

        # Metal transport
        D_metal = full_params['D_ref']   * np.exp((-full_params['E_D']      / R_gas) * inv_m)
        K_s_met = full_params['K_s_ref'] * np.exp((-full_params['H_s']      / R_gas) * inv_m)

        # Oxide transport
        D_ox = full_params['D_ox_ref'] * np.exp((-full_params['E_D_ox']   / R_gas) * inv_ox)
        K_ox = full_params['K_ox_ref'] * np.exp((-full_params['H_sol_ox'] / R_gas) * inv_ox)

        # Oxide surface kinetics
        k_diss = full_params['k_diss_ref'] * np.exp((-full_params['E_diss'] / R_gas) * inv_s)
        K_eq   = full_params['K_eq_ref']   * np.exp((-full_params['H_eq']   / R_gas) * inv_s)

        # Metal surface kinetics (pinholes) — own T_ref from Grant 1988 (965 K)
        k_diss_metal = full_params['k_diss_metal_ref'] * np.exp((-full_params['E_diss_metal'] / R_gas) * inv_sm)
        K_eq_metal   = full_params['K_eq_metal_ref']   * np.exp((-full_params['H_eq_metal']   / R_gas) * inv_sm)

        L_ox = full_params['oxide_thickness']
        L_m  = full_params['metal_thickness']
        P_up = full_params['P_upstream']
        lattice_density = full_params['lattice_density']
        method = full_params.get('D_eff_method', 'average')

        grain_size   = full_params['grain_size']
        gb_thickness = full_params['gb_thickness']
        N_T_gb       = full_params.get('trap_gb_N_T', 0)

        trap_list = []
        if full_params.get('trap_dislocation_N_T', 0) > 0:
            trap_list.append({'name': 'dislocations',
                              'binding_energy': full_params['trap_dislocation_E_b'],
                              'density': full_params['trap_dislocation_N_T']})
        if full_params.get('trap_vacancy_N_T', 0) > 0:
            trap_list.append({'name': 'vacancies',
                              'binding_energy': full_params['trap_vacancy_E_b'],
                              'density': full_params['trap_vacancy_N_T']})
        if full_params.get('trap_gb_E_b', 0) > 0 and N_T_gb > 0:
            trap_list.append({'name': 'grain_boundaries',
                              'binding_energy': full_params['trap_gb_E_b'],
                              'density': N_T_gb})
        N_T_carbide = full_params.get('trap_carbide_N_T', 0)
        if full_params.get('trap_carbide_E_b', 0) > 0 and N_T_carbide > 0:
            trap_list.append({'name': 'carbides',
                              'binding_energy': full_params['trap_carbide_E_b'],
                              'density': N_T_carbide})

        microstructure_params = {
            'grain_size':            grain_size,
            'grain_shape':           full_params.get('grain_shape', 'equiaxed'),
            'gb_type':               full_params.get('gb_type', 'LAGB'),
            'gb_thickness':          gb_thickness,
            'trap_list':             trap_list,
            'gb_enhancement_factor': full_params.get('gb_enhancement_factor', 100),
        }

        include_gb   = full_params.get('include_gb_enhancement', True)
        include_trap = full_params.get('include_trapping', True)
        if include_gb and include_trap:
            mode = 'both'
        elif include_gb:
            mode = 'gb_only'
        elif include_trap:
            mode = 'trapping_only'
        else:
            mode = 'none'

        f_pin = full_params.get('f_pinhole',   0.0)
        f_cra = full_params.get('f_crack',     0.0)
        f_gb  = full_params.get('f_gb_defect', 0.0)
        defect_config = {}
        if f_pin > 0:
            defect_config['pinhole'] = {'area_fraction': f_pin}
        if f_cra > 0:
            defect_config['crack'] = {
                'area_fraction':  f_cra,
                'thickness_factor': full_params.get('crack_thickness_factor', 0.1),
            }
        if f_gb > 0:
            defect_config['grain_boundary'] = {
                'area_fraction':     f_gb,
                'diffusivity_factor': full_params.get('gb_diffusivity_factor', 10.0),
            }
        if not defect_config:
            defect_config['crack'] = {'area_fraction': 1e-6}

        r = calculate_full_model_flux_L346_v2(
            P_up=P_up, P_down=0.0, L_m=L_m, temperature=T,
            k_diss=k_diss, K_eq=K_eq,
            D_ox=D_ox, K_ox=K_ox, L_ox=L_ox,
            D_lattice=D_metal, K_s_m=K_s_met,
            microstructure_params=microstructure_params,
            defect_config=defect_config,
            lattice_density=lattice_density,
            method=method,
            n_points=10,
            mode=mode,
            k_diss_metal=k_diss_metal,
            K_eq_metal=K_eq_metal,
        )

        fw           = r.get('flux_weighted_resistances', {})
        intact_theta = r.get('intact_path', {}).get('theta', np.nan)

        # Effective permeability via harmonic mean (series resistances)
        # 1/Φ_eff = 1/Φ_oxide + 1/Φ_metal (no length dependence)
        Phi_oxide_6 = D_ox * K_ox
        Phi_metal_6 = r.get('D_eff_avg', D_metal) * K_s_met

        if Phi_oxide_6 > 0 and Phi_metal_6 > 0 and not np.isnan(Phi_oxide_6) and not np.isnan(Phi_metal_6):
            permeability_6 = 1.0 / (1.0/Phi_oxide_6 + 1.0/Phi_metal_6)
        else:
            permeability_6 = np.nan

        return {
            'flux':           r['J_total'],
            'PRF':            np.nan,
            'D_eff':          r.get('D_eff_avg', np.nan),
            'D_modification': r.get('overall_modification_factor', np.nan),
            'permeability':   permeability_6,
            'P_interface':    r.get('intact_path', {}).get('P_int', np.nan),
            'flux_intact':    r['flux_breakdown'].get('intact', {}).get('contribution', np.nan),
            'flux_defect':    (r['J_total']
                               - r['flux_breakdown'].get('intact', {}).get('contribution', 0)),
            'frac_surface':   fw.get('fraction_surface', np.nan),
            'frac_oxide':     fw.get('fraction_oxide',   np.nan),
            'frac_metal':     fw.get('fraction_metal',   np.nan),
            'theta':          intact_theta,
            'D_metal':        D_metal,
            'K_s_metal':      K_s_met,
            'D_ox':           D_ox,
            'K_ox':           K_ox,
            'k_diss':         k_diss,
            'K_eq':           K_eq,
            'k_diss_metal':   k_diss_metal,
            'K_eq_metal':     K_eq_metal,
            'temperature':    T,
        }

    except Exception as e:
        print(f"Error in Level 5L6 model: {e}")
        import traceback; traceback.print_exc()
        return {
            'flux': 1e-20, 'PRF': np.nan, 'D_eff': 1e-12, 'D_modification': np.nan,
            'permeability': 1e-20, 'P_interface': np.nan,
            'flux_intact': 1e-20, 'flux_defect': 0.0,
            'frac_surface': np.nan, 'frac_oxide': np.nan, 'frac_metal': np.nan,
            'theta': np.nan, 'D_metal': 1e-12, 'K_s_metal': 1e-6,
            'D_ox': 1e-15, 'K_ox': 1e-10, 'k_diss': 1e-15, 'K_eq': 1e-10,
            'k_diss_metal': 1e-12, 'K_eq_metal': 1e-8,
            'temperature': full_params.get('temperature', 1073.15),
        }


# =============================================================================
# MORRIS SENSITIVITY ANALYSIS — LEVEL 5L6
# =============================================================================

def morris_sensitivity_level5L6(
    param_ranges,
    N_trajectories=10,
    num_levels=4,
    seed=42,
    output_metric='flux'
):
    """
    Morris sensitivity analysis for LEVEL 5+L6 (surface kinetics) model.

    Same methodology as morris_sensitivity_level5 but uses level5L6_model_wrapper.
    Additional output metrics: frac_surface, frac_oxide, frac_metal, theta.

    Parameters
    ----------
    param_ranges : dict
        Use SUGGESTED_RANGES_LEVEL5L6 as starting point.
    N_trajectories : int
        Number of Morris trajectories (default: 10).
    num_levels : int
        Number of grid levels (default: 4).
    output_metric : str
        One of VALID_OUTPUT_METRICS_L5L6.

    Returns
    -------
    Si, problem, Y : same structure as morris_sensitivity_level5
    """
    if output_metric not in VALID_OUTPUT_METRICS_L5L6:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS_L5L6}, got '{output_metric}'")

    problem = {
        'num_vars': len(param_ranges),
        'names':    list(param_ranges.keys()),
        'bounds':   list(param_ranges.values()),
    }

    param_values = morris_sampler.sample(problem, N=N_trajectories, num_levels=num_levels, seed=seed)
    n_samples = param_values.shape[0]
    Y = np.zeros(n_samples)

    print(f"\n{'='*70}")
    print(f"MORRIS SENSITIVITY ANALYSIS - LEVEL 5L6 (Surface Kinetics Model)")
    print(f"{'='*70}")
    print(f"Running {n_samples} Morris samples...")
    print(f"Parameters ({len(problem['names'])}): {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print(f"{'='*70}\n")

    for i, X in enumerate(param_values):
        result = level5L6_model_wrapper(dict(zip(problem['names'], X)))
        Y[i] = result.get(output_metric, np.nan)
        if (i + 1) % 10 == 0 or (i + 1) == n_samples:
            print(f"  Completed {i + 1}/{n_samples} samples")

    valid_idx = np.isfinite(Y) & (Y != 0)
    n_invalid = np.sum(~valid_idx)
    if n_invalid > 0:
        print(f"\n  Warning: {n_invalid} invalid outputs replaced with median")
        Y[~valid_idx] = np.median(Y[valid_idx]) if np.any(valid_idx) else 1e-15

    Si = morris_analyzer.analyze(problem, param_values, Y, conf_level=0.95, print_to_console=False)

    print(f"\n  Morris analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    mu = Si['mu']
    for rank, idx in enumerate(np.argsort(np.abs(mu))[::-1][:10]):
        print(f"    {rank+1}. {problem['names'][idx]}: μ = {mu[idx]:.2e}")

    return Si, problem, Y


# =============================================================================
# SOBOL SENSITIVITY ANALYSIS — LEVEL 5L6
# =============================================================================

def sobol_sensitivity_level5L6(
    param_ranges,
    N_samples=1024,
    output_metric='flux',
    seed=42,
    calc_second_order=True
):
    """
    Sobol variance-based sensitivity analysis for LEVEL 5+L6 (surface kinetics) model.

    Same methodology as sobol_sensitivity_level5 but uses level5L6_model_wrapper.

    Parameters
    ----------
    param_ranges : dict
        Use top parameters from Morris screening (SUGGESTED_RANGES_LEVEL5L6 subset).
    N_samples : int
        Base number of samples (power of 2: 256, 512, 1024).
    output_metric : str
        One of VALID_OUTPUT_METRICS_L5L6.
    calc_second_order : bool
        Whether to compute pairwise S2 indices.

    Returns
    -------
    Si, problem, Y : same structure as sobol_sensitivity_level5
    """
    if output_metric not in VALID_OUTPUT_METRICS_L5L6:
        raise ValueError(f"output_metric must be one of {VALID_OUTPUT_METRICS_L5L6}, got '{output_metric}'")

    problem = {
        'num_vars': len(param_ranges),
        'names':    list(param_ranges.keys()),
        'bounds':   list(param_ranges.values()),
    }

    if seed is not None:
        np.random.seed(seed)

    param_values = sobol_sampler.sample(problem, N_samples, calc_second_order=calc_second_order)
    n_samples = param_values.shape[0]
    Y = np.zeros(n_samples)

    print(f"\n{'='*70}")
    print(f"SOBOL SENSITIVITY ANALYSIS - LEVEL 5L6 (Surface Kinetics Model)")
    print(f"{'='*70}")
    print(f"Running {n_samples} Sobol samples (N_base={N_samples}, {len(param_ranges)} params)")
    print(f"Parameters: {', '.join(problem['names'])}")
    print(f"Output metric: {output_metric}")
    print(f"{'='*70}\n")

    for i, X in enumerate(param_values):
        result = level5L6_model_wrapper(dict(zip(problem['names'], X)))
        Y[i] = result.get(output_metric, np.nan)
        if (i + 1) % max(1, n_samples // 10) == 0 or (i + 1) == n_samples:
            print(f"  Completed {i + 1}/{n_samples} samples ({(i+1)/n_samples*100:.0f}%)")

    valid_idx = np.isfinite(Y) & (Y != 0)
    n_invalid = np.sum(~valid_idx)
    if n_invalid > 0:
        print(f"\n  Warning: {n_invalid} invalid outputs replaced with median")
        Y[~valid_idx] = np.median(Y[valid_idx]) if np.any(valid_idx) else 1e-15

    Si = sobol_analyzer.analyze(problem, Y, calc_second_order=calc_second_order,
                                conf_level=0.95, print_to_console=False)

    print(f"\n  Sobol analysis complete")
    print(f"  Output range: [{np.min(Y):.2e}, {np.max(Y):.2e}]")
    S1, ST = Si['S1'], Si['ST']
    print(f"\n  {'Parameter':<25} {'S1':>8} {'ST':>8}")
    print(f"  {'-'*25} {'-'*8} {'-'*8}")
    for idx in np.argsort(ST)[::-1]:
        print(f"  {problem['names'][idx]:<25} {S1[idx]:>8.3f} {ST[idx]:>8.3f}")
    print(f"\n  Sum(S1)={np.sum(S1):.3f}  Sum(ST)={np.sum(ST):.3f}")

    return Si, problem, Y
