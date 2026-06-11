import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from SALib.sample import morris as morris_sampler
from SALib.analyze import morris as morris_analyzer
from SALib.sample import saltelli as sobol_sampler
from SALib.analyze import sobol as sobol_analyzer
from SALib.sample import latin as latin_sampler
from SALib.analyze import pawn as pawn_analyzer
from SALib.analyze import delta as delta_analyzer

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
    #'T_ref_metal':  [600, 1400],

    # Oxide transport (Cr2O3_sample4)
    'D_ox_ref':     [7.800e-21, 7.800e-17],
    'E_D_ox':       [69000, 71000],
    'K_ox_ref':     [0.08, 0.4],
    'H_sol_ox':     [20000, 170000],
    #'T_ref_oxide':  [600, 1400],

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

VALID_OUTPUT_METRICS_L5L6  = VALID_OUTPUT_METRICS_L5
# VALID_OUTPUT_METRICS_L5L6 = [m for m in VALID_OUTPUT_METRICS_L5 if m != 'PRF'] + [
#     'frac_surface',  # Flux-weighted surface resistance fraction
#     'frac_oxide',    # Flux-weighted oxide resistance fraction
#     'frac_metal',    # Flux-weighted metal resistance fraction
#     'theta',         # Surface coverage on intact path at steady state
# ]

SUGGESTED_RANGES_LEVEL5L6 = {
    **SUGGESTED_RANGES_LEVEL5,
    # Oxide surface kinetics (Cr2O3)
    'k_diss_ref':          [9.487e-10, 9.487e-6],
    'E_diss':              [40000, 70000],
    'K_eq_ref':            [1e-14, 1e-6],
    'H_eq':                [5000, 50000],
    #'T_ref_surface':       [600, 1800],
    # Metal surface kinetics (Grant 1988)
    'k_diss_metal_ref':    [1.346e-8, 1.346e-4],
    'E_diss_metal':        [70000, 100000],
    'K_eq_metal_ref':      [1e-10, 1e-3],
    'H_eq_metal':          [5000, 40000],
    #'T_ref_surface_metal': [600, 1100],
}


# =============================================================================
# REGIME LABELLING — single source of truth
# =============================================================================

REGIME_LABELS = ('surface', 'oxide', 'metal')


def assign_regime(frac_surface, frac_oxide, frac_metal, rule='argmax', threshold=0.5):
    """
    Canonical rate-limiting-regime label from the flux-weighted resistance fractions.

    This is the single source of truth used by the regime-stratified sensitivity
    analysis (and anywhere else a regime label is needed), so the definition stays
    consistent across the codebase.

    Parameters
    ----------
    frac_surface, frac_oxide, frac_metal : float
        Flux-weighted resistance fractions (from
        calculate_full_model_flux_L346_v2 -> 'flux_weighted_resistances').
        They should sum to ~1 on a valid solve.
    rule : {'argmax', 'threshold'}
        'argmax'     : regime = whichever fraction is largest (default; always a
                       single label, no 'mixed' bucket).
        'threshold'  : regime = the largest fraction only if it exceeds `threshold`,
                       otherwise 'mixed' (matches the notebook's dominant-step rule).
    threshold : float
        Dominance cutoff used only when rule='threshold'.

    Returns
    -------
    str
        One of 'surface', 'oxide', 'metal' (or 'mixed' for rule='threshold'),
        or 'undefined' if any fraction is NaN/None (e.g. a failed solve).
    """
    fracs = {'surface': frac_surface, 'oxide': frac_oxide, 'metal': frac_metal}

    # Failed / incomplete solve -> undefined; excluded from regime clusters.
    for v in fracs.values():
        if v is None or (isinstance(v, (float, np.floating)) and np.isnan(v)):
            return 'undefined'

    label, val = max(fracs.items(), key=lambda kv: kv[1])

    if rule == 'argmax':
        return label
    elif rule == 'threshold':
        return label if val > threshold else 'mixed'
    else:
        raise ValueError(f"rule must be 'argmax' or 'threshold', got '{rule}'")


# =============================================================================
# MODEL WRAPPER — LEVEL 5L6
# =============================================================================

def level5L6_model_wrapper(params_dict, return_full_record=False):
    """
    Wrapper for LEVEL 5+L6 model: defective oxide + defective metal + surface kinetics.

    Extends Level 5 by adding Langmuir-Hinshelwood surface kinetics at both the
    gas-oxide interface (k_diss, K_eq) and the metal surface at pinholes
    (k_diss_metal, K_eq_metal). Calls calculate_full_model_flux_L346_v2.

    Parameters
    ----------
    params_dict : dict
        Any subset of parameters from DEFAULT_PARAMS_LEVEL5L6.
    return_full_record : bool
        When True, also attach the string diagnostics 'system_rate_limiting' and
        'dominant_path' to the returned dict. Used by the regime-stratified global
        scan (run_global_lhs_scan). Default False leaves the existing scalar-metric
        callers (Morris/Sobol loops) unaffected.

    Returns
    -------
    dict
        All Level 5 outputs plus:
        - 'frac_surface': Flux-weighted surface resistance fraction
        - 'frac_oxide':   Flux-weighted oxide resistance fraction
        - 'frac_metal':   Flux-weighted metal resistance fraction
        - 'theta':        Surface coverage on intact path
        - 'regime':       Rate-limiting regime label via assign_regime() (argmax of
                          the flux-weighted fractions; 'undefined' on a failed solve)
        and, when return_full_record=True, 'system_rate_limiting' and 'dominant_path'.
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

        # Canonical regime label (argmax of the flux-weighted fractions).
        regime = assign_regime(
            fw.get('fraction_surface', np.nan),
            fw.get('fraction_oxide',   np.nan),
            fw.get('fraction_metal',   np.nan),
        )

        # Effective permeability via harmonic mean (series resistances)
        # 1/Φ_eff = 1/Φ_oxide + 1/Φ_metal (no length dependence)
        Phi_oxide_6 = D_ox * K_ox
        Phi_metal_6 = r.get('D_eff_avg', D_metal) * K_s_met

        if Phi_oxide_6 > 0 and Phi_metal_6 > 0 and not np.isnan(Phi_oxide_6) and not np.isnan(Phi_metal_6):
            permeability_6 = 1.0 / (1.0/Phi_oxide_6 + 1.0/Phi_metal_6)
        else:
            permeability_6 = np.nan

        record = {
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
            'regime':         regime,
        }
        if return_full_record:
            record['system_rate_limiting'] = r.get('system_rate_limiting', 'unknown')
            record['dominant_path']        = r.get('dominant_path', 'unknown')
        return record

    except Exception as e:
        print(f"Error in Level 5L6 model: {e}")
        import traceback; traceback.print_exc()
        record = {
            'flux': 1e-20, 'PRF': np.nan, 'D_eff': 1e-12, 'D_modification': np.nan,
            'permeability': 1e-20, 'P_interface': np.nan,
            'flux_intact': 1e-20, 'flux_defect': 0.0,
            'frac_surface': np.nan, 'frac_oxide': np.nan, 'frac_metal': np.nan,
            'theta': np.nan, 'D_metal': 1e-12, 'K_s_metal': 1e-6,
            'D_ox': 1e-15, 'K_ox': 1e-10, 'k_diss': 1e-15, 'K_eq': 1e-10,
            'k_diss_metal': 1e-12, 'K_eq_metal': 1e-8,
            'temperature': full_params.get('temperature', 1073.15),
            'regime': 'undefined',
        }
        if return_full_record:
            record['system_rate_limiting'] = 'error'
            record['dominant_path']        = 'error'
        return record


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


# =============================================================================
# REGIME-STRATIFIED SENSITIVITY ANALYSIS — LEVEL 5L6
# Unified Global LHS Foundation + global Morris screen + per-regime Route A/B.
# See plan: argmax regime labelling; Route B (given-data on true clusters) is
# primary, Route A (accept-and-report Saltelli over a reduced k-box) is secondary.
# =============================================================================

# Scalar model outputs stored as columns in the global-scan DataFrame.
SCAN_OUTPUT_FIELDS = [
    'flux', 'permeability', 'theta',
    'frac_surface', 'frac_oxide', 'frac_metal',
    'PRF', 'D_eff', 'P_interface', 'flux_intact', 'flux_defect',
]
# String diagnostics stored alongside.
SCAN_LABEL_FIELDS = ['regime', 'system_rate_limiting', 'dominant_path']
# Default per-cluster output metrics for the regime SA. 'flux' is primary (the
# only metric that responds to surface kinetics); 'permeability' is bulk-only
# (surface-insensitive by construction); 'theta' is the surface coverage.
# The regime-defining frac_* are intentionally excluded as per-cluster metrics
# (degenerate within their own cluster).
REGIME_SA_METRICS = ['flux', 'permeability', 'theta']


def _make_problem(param_ranges):
    """SALib problem dict from an ordered {name: [lo, hi]} mapping."""
    return {
        'num_vars': len(param_ranges),
        'names':    list(param_ranges.keys()),
        'bounds':   [list(v) for v in param_ranges.values()],
    }


# Metrics that span orders of magnitude — analysed on a log10 scale so the
# variance/density estimators behind Morris/Sobol/δ aren't dominated by a few
# huge values (raw flux ranges ~10 decades). theta (∈[0,1]) stays linear.
LOG_METRICS_DEFAULT = ('flux', 'permeability')


def _prep_Y_full(Y, metric, log_metrics):
    """
    Full-length Y for structured estimators (Morris/Sobol) — cannot drop rows.
    Applies log10 for heavy-tailed metrics (non-positive -> NaN) and imputes any
    NaN with the median so the sampling design stays intact.
    """
    Y = np.asarray(Y, dtype=float).copy()
    if metric in log_metrics:
        good = np.isfinite(Y) & (Y > 0)
        out = np.full(Y.shape, np.nan)
        out[good] = np.log10(Y[good])
    else:
        out = np.where(np.isfinite(Y), Y, np.nan)
    bad = ~np.isfinite(out)
    if bad.any():
        out[bad] = np.median(out[~bad]) if (~bad).any() else 0.0
    return out


def _save_scan(df, path):
    """Persist the scan DataFrame; parquet if requested+available, else CSV."""
    path = str(path)
    if path.endswith('.parquet'):
        try:
            df.to_parquet(path)
            print(f"  saved scan -> {path}")
            return
        except Exception as e:
            print(f"  parquet save failed ({e}); falling back to CSV")
            path = path.rsplit('.', 1)[0] + '.csv'
    df.to_csv(path, index=False)
    print(f"  saved scan -> {path}")


# -----------------------------------------------------------------------------
# Phase 1 — Step 3: Global LHS foundation (one evaluation pass, full records)
# -----------------------------------------------------------------------------

def run_global_lhs_scan(param_ranges, N, seed=42, wrapper=None,
                        save_path=None, verbose=True):
    """
    Space-filling LHS over all parameters, evaluated once with the full record.

    This single pass feeds BOTH the descriptive partitioning and Route B
    (given-data SA on the true regime clusters) — zero further model evaluations
    are needed for Route B.

    Parameters
    ----------
    param_ranges : dict
        {name: [lo, hi]}; e.g. SUGGESTED_RANGES_LEVEL5L6.
    N : int
        LHS sample size. Choose so the rarest regime cluster has >= a few hundred
        points (set by the rarest regime's frequency).
    seed : int
        Reproducibility seed for latin.sample.
    wrapper : callable
        Model wrapper (default level5L6_model_wrapper); called with
        return_full_record=True.
    save_path : str or None
        If given, persist the tidy DataFrame (CSV, or parquet by extension).

    Returns
    -------
    (df, problem)
        df : tidy DataFrame — one row per LHS draw; columns = sampled parameter
             values (the SA inputs X) + SCAN_OUTPUT_FIELDS + SCAN_LABEL_FIELDS.
        problem : SALib problem dict (names/bounds) for reuse downstream.
    """
    if wrapper is None:
        wrapper = level5L6_model_wrapper

    problem = _make_problem(param_ranges)
    names   = problem['names']
    X       = latin_sampler.sample(problem, N, seed=seed)
    n       = X.shape[0]

    if verbose:
        print(f"\n{'='*70}")
        print(f"GLOBAL LHS SCAN — {n} samples × {len(names)} params (seed={seed})")
        print(f"{'='*70}")

    rows = []
    for i in range(n):
        params = dict(zip(names, X[i]))
        rec    = wrapper(params, return_full_record=True)
        row    = dict(params)  # sampled values ARE the SA inputs
        for f in SCAN_OUTPUT_FIELDS:
            row[f] = rec.get(f, np.nan)
        for f in SCAN_LABEL_FIELDS:
            row[f] = rec.get(f, 'undefined')
        rows.append(row)
        if verbose and ((i + 1) % max(1, n // 10) == 0 or (i + 1) == n):
            print(f"  scan {i + 1}/{n} ({(i + 1) / n * 100:.0f}%)")

    df = pd.DataFrame(rows)

    if save_path:
        _save_scan(df, save_path)

    if verbose:
        print("\n  Regime counts:")
        for reg, c in df['regime'].value_counts().items():
            print(f"    {reg:<10s}: {c:6d}  ({c / n * 100:5.1f}%)")

    return df, problem


# -----------------------------------------------------------------------------
# Phase 1 — Step 4: Partition into regime clusters
# -----------------------------------------------------------------------------

def partition_by_regime(df, regime_col='regime', min_cluster=300,
                        drop_undefined=True, verbose=True):
    """
    Split the scan DataFrame into {regime: sub_df}.

    'undefined' (failed-solve) rows are dropped from the clusters by default and
    reported separately. Warns when a cluster is too small for stable given-data
    indices (δ/PAWN).
    """
    counts = df[regime_col].value_counts()
    n = len(df)
    if verbose:
        print("\n  Regime cluster sizes:")
        for reg, c in counts.items():
            print(f"    {reg:<10s}: {c:6d}  ({c / n * 100:5.1f}%)")

    work = df
    n_undef = int((df[regime_col] == 'undefined').sum())
    if drop_undefined and n_undef:
        if verbose:
            print(f"  dropping {n_undef} 'undefined' (failed-solve) rows from clusters")
        work = df[df[regime_col] != 'undefined']

    partition = {reg: sub.reset_index(drop=True)
                 for reg, sub in work.groupby(regime_col)}

    if verbose:
        for reg, sub in partition.items():
            if len(sub) < min_cluster:
                print(f"  ⚠ cluster '{reg}' has {len(sub)} < {min_cluster} pts — "
                      f"given-data δ/PAWN may be noisy; consider a targeted top-up draw.")

    return partition


# -----------------------------------------------------------------------------
# Phase 1 — Step 5: Descriptive plots
# -----------------------------------------------------------------------------

def plot_regime_exploration(partition, output_metrics=('flux', 'permeability', 'theta'),
                            scan_df=None, show=True):
    """
    Per-regime descriptive view: cluster sizes + output distributions by regime,
    and (if scan_df given) a regime × system_rate_limiting agreement cross-tab.
    """
    regimes = list(partition.keys())
    sizes   = [len(partition[r]) for r in regimes]

    n_panels = 1 + len(output_metrics)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    axes[0].bar(regimes, sizes, color='steelblue', alpha=0.85)
    axes[0].set_title('Regime cluster sizes')
    axes[0].set_ylabel('count')
    axes[0].tick_params(axis='x', rotation=30)

    for ax, m in zip(axes[1:], output_metrics):
        data = []
        for r in regimes:
            vals = partition[r][m].replace([np.inf, -np.inf], np.nan).dropna()
            data.append(vals.values)
        ax.boxplot(data, showfliers=False)
        ax.set_xticks(range(1, len(regimes) + 1))
        ax.set_xticklabels(regimes, rotation=30)
        all_pos = all(len(d) and np.all(d > 0) for d in data)
        ax.set_yscale('log' if all_pos else 'linear')
        ax.set_title(f'{m} by regime')

    fig.tight_layout()
    if show:
        plt.show()

    # Agreement check: argmax-regime vs the solver's own system_rate_limiting.
    if scan_df is not None and 'system_rate_limiting' in scan_df.columns:
        valid = scan_df[scan_df['regime'] != 'undefined']
        ct = pd.crosstab(valid['regime'], valid['system_rate_limiting'])
        print("\n  regime (argmax)  ×  system_rate_limiting (solver):")
        print(ct.to_string())

    return fig


# -----------------------------------------------------------------------------
# Phase 1 — Step 6: Global Morris screen (multi-metric, single eval pass) + rank
# -----------------------------------------------------------------------------

def morris_screen_global(param_ranges, output_metrics=REGIME_SA_METRICS,
                         N_trajectories=500, num_levels=4, seed=42,
                         wrapper=None, log_metrics=LOG_METRICS_DEFAULT,
                         cache_csv=None, verbose=True):
    """
    Standard 40-D Morris screen over the FULL space (no regime filtering — clean,
    valid, structured). Evaluates the Morris design once and analyses every output
    metric on the shared sample (the multi-metric efficiency win).

    cache_csv : if given and the file exists, the (X, Y) model evaluations are
        loaded from it instead of re-running the model (the analysis is always
        recomputed — it is cheap). If it does not exist, the evaluations are saved
        there after running. Lets reruns skip the expensive model calls.

    Returns
    -------
    dict : {metric: {'Si': Si, 'problem': problem, 'Y': Y, 'log': bool}}
    """
    if wrapper is None:
        wrapper = level5L6_model_wrapper

    problem = _make_problem(param_ranges)
    names   = problem['names']
    X       = morris_sampler.sample(problem, N=N_trajectories,
                                    num_levels=num_levels, seed=seed)
    n = X.shape[0]

    if cache_csv and os.path.exists(cache_csv):
        cdf = pd.read_csv(cache_csv)
        if all(c in cdf.columns for c in names) and len(cdf) == n:
            X = cdf[names].to_numpy(dtype=float)
        Ys = {m: cdf[f'Y_{m}'].to_numpy(dtype=float) for m in output_metrics
              if f'Y_{m}' in cdf.columns}
        if verbose:
            print(f"  [cache] loaded Morris evals from {cache_csv} ({len(cdf)} rows)")
    else:
        if verbose:
            print(f"\n{'='*70}")
            print(f"GLOBAL MORRIS SCREEN — {n} samples × {len(names)} params, "
                  f"metrics={list(output_metrics)}")
            print(f"{'='*70}")
        Ys = {m: np.full(n, np.nan) for m in output_metrics}
        for i in range(n):
            rec = wrapper(dict(zip(names, X[i])), return_full_record=True)
            for m in output_metrics:
                Ys[m][i] = rec.get(m, np.nan)
            if verbose and ((i + 1) % max(1, n // 10) == 0 or (i + 1) == n):
                print(f"  morris-screen {i + 1}/{n}")
        if cache_csv:
            d = pd.DataFrame(X, columns=names)
            for m in output_metrics:
                d[f'Y_{m}'] = Ys[m]
            d.to_csv(cache_csv, index=False)
            if verbose:
                print(f"  [cache] saved Morris evals -> {cache_csv}")

    results = {}
    for m in output_metrics:
        Y = _prep_Y_full(Ys[m], m, log_metrics)
        Si = morris_analyzer.analyze(problem, X, Y, conf_level=0.95,
                                     print_to_console=False)
        results[m] = {'Si': Si, 'problem': problem, 'Y': Ys[m],
                      'log': m in log_metrics}

    return results


def rank_and_select_top_k(screen_results, output_metrics=REGIME_SA_METRICS,
                          k=10, verbose=True):
    """
    Rank parameters by μ* normalised per metric then averaged across metrics, and
    select the top-k. Mirrors the existing notebook ranking logic (cell 6).

    Returns
    -------
    (df_rank, top_k_names)
    """
    names = screen_results[output_metrics[0]]['problem']['names']
    df = pd.DataFrame({'Parameter': names})

    for m in output_metrics:
        df[f'mu_star_{m}'] = screen_results[m]['Si']['mu_star']

    norm_cols = []
    for m in output_metrics:
        col = f'mu_star_{m}'
        mx  = df[col].abs().max()
        nc  = f'{col}_norm'
        df[nc] = df[col] / mx if mx > 0 else 0.0
        norm_cols.append(nc)

    df['Avg_Importance'] = df[norm_cols].mean(axis=1)
    df = df.sort_values('Avg_Importance', ascending=False).reset_index(drop=True)
    top_k = df.head(k)['Parameter'].tolist()

    if verbose:
        show_cols = ['Parameter'] + [f'mu_star_{m}' for m in output_metrics] + ['Avg_Importance']
        print("\n  Global Morris ranking (μ* normalised & averaged):")
        print(df[show_cols].to_string(index=False))
        print(f"\n  Top {k}: {top_k}")

    return df, top_k


# -----------------------------------------------------------------------------
# Shared helper: top-k percentile sub-box (used by the review gate and Route A)
# -----------------------------------------------------------------------------

def extract_topk_subbox(cluster_df, top_k, low=5, high=95, fallback_ranges=None):
    """
    Axis-aligned percentile-band box on the top-k params for one regime cluster.
    Degenerate (zero-width) bounds fall back to the global range when provided.
    """
    box = {}
    for p in top_k:
        lo = float(np.percentile(cluster_df[p], low))
        hi = float(np.percentile(cluster_df[p], high))
        if hi <= lo:
            if fallback_ranges and p in fallback_ranges:
                lo, hi = float(fallback_ranges[p][0]), float(fallback_ranges[p][1])
            else:
                span = abs(lo) if lo != 0 else 1.0
                lo, hi = lo - 1e-6 * span, hi + 1e-6 * span
        box[p] = [lo, hi]
    return box


def estimate_subbox_purity(scan_df, partition, top_k, regime_col='regime',
                           low=5, high=95, fallback_ranges=None, verbose=True):
    """
    REVIEW-GATE estimator: for each regime, build its top-k percentile sub-box and
    measure purity p = fraction of ALL scan points inside that box actually labelled
    as the target regime. Uses existing scan data — zero new model evaluations.

    Low purity => a rectangular sub-box leaks heavily (curved boundary) => Route A
    is 'neighbourhood' sensitivity, and per-regime Morris-with-discard would starve
    (survival ~ p^(k+1)).

    Returns
    -------
    dict : {regime: {'n_in_box': int, 'purity': float, 'box': {param: [lo, hi]}}}
    """
    out = {}
    for reg, sub in partition.items():
        box = extract_topk_subbox(sub, top_k, low, high, fallback_ranges)
        mask = np.ones(len(scan_df), dtype=bool)
        for p, (lo, hi) in box.items():
            mask &= (scan_df[p] >= lo) & (scan_df[p] <= hi)
        inside = scan_df[mask]
        n_in   = int(len(inside))
        purity = float((inside[regime_col] == reg).mean()) if n_in else float('nan')
        out[reg] = {'n_in_box': n_in, 'purity': purity, 'box': box}

    if verbose:
        print("\n  Sub-box purity (top-k percentile box; estimated from scan):")
        print(f"    {'regime':<10s} {'n_in_box':>9s} {'purity':>8s}  "
              f"~Morris-traj survival p^(k+1)")
        kk = len(top_k)
        for reg, d in out.items():
            surv = d['purity'] ** (kk + 1) if np.isfinite(d['purity']) else float('nan')
            print(f"    {reg:<10s} {d['n_in_box']:>9d} {d['purity']:>8.3f}  {surv:>10.2e}")

    return out


# -----------------------------------------------------------------------------
# Regime-targeted sampling presets (Option A)
# Each preset = the default ranges with ONLY the regime-controlling parameters
# overridden, left as broad as possible so the per-regime SA reflects genuine
# sensitivity rather than the targeting itself. Verified yields (uniform LHS):
#   metal   ~97%   (default ranges)
#   surface ~32%   (low pressure + slow surface kinetics)
#   oxide   ~11%   (suppressed defects + slow/thick oxide + fast surface & metal)
# -----------------------------------------------------------------------------

def _ranges_with(overrides):
    R = {k: list(v) for k, v in SUGGESTED_RANGES_LEVEL5L6.items()}
    R.update({k: list(v) for k, v in overrides.items()})
    return R


REGIME_PRESETS = {
    'metal': _ranges_with({}),  # default ranges already ~97% metal-limited

    'surface': _ranges_with({
        'P_upstream':       [1e-7, 1e1],     # low pressure
        'k_diss_ref':       [9.5e-12, 9.5e-9],   # slow oxide-surface dissociation
        'k_diss_metal_ref': [1.3e-10, 1.3e-7],   # slow metal-surface dissociation
    }),

    'oxide': _ranges_with({
        'f_pinhole':   [1e-7, 1e-5],   # suppress defect bypass so flux uses the intact oxide
        'f_crack':     [1e-7, 1e-5],
        'f_gb_defect': [1e-7, 1e-5],
        'oxide_thickness': [1e-7, 1e-5],     # thick, slow oxide -> oxide is the bottleneck
        'D_ox_ref':        [7.8e-21, 7.8e-20],
        'P_upstream':      [1e4, 1e7],       # fast surface (high P + high k_diss)
        'k_diss_ref':      [9.5e-8, 9.5e-6],
        'metal_thickness': [5e-5, 5e-4],     # thin, fast metal -> metal not limiting
        'D_ref':           [4.3e-8, 4.3e-7],
    }),
}

# Default LHS sizes per regime, chosen from the verified yields so each cluster
# clears the ~300-500 floor for stable given-data indices.
DEFAULT_N_PER_REGIME = {'metal': 1200, 'surface': 3000, 'oxide': 6000}


def run_targeted_regime_scans(N_per_regime=None, regimes=('metal', 'surface', 'oxide'),
                              seed=42, wrapper=None, save_dir=None,
                              min_cluster=300, verbose=True):
    """
    Phase 1 (Option A): one LHS scan per regime over its targeted preset, keeping
    the in-regime rows (argmax label). Feeds Route B (given-data) and Route A.

    Returns
    -------
    clusters : {regime: cluster_df}   in-regime rows only (regime == target preset)
    scans    : {regime: full_scan_df} every draw from that preset (for yield diagnostics)
    master   : DataFrame              all clusters concatenated (+ 'target_preset' col)
    yields   : {regime: in_regime_fraction}
    """
    if N_per_regime is None:
        N_per_regime = dict(DEFAULT_N_PER_REGIME)
    if isinstance(N_per_regime, int):
        N_per_regime = {r: N_per_regime for r in regimes}

    clusters, scans, yields = {}, {}, {}
    for r in regimes:
        if verbose:
            print(f"\n--- targeted scan: {r} (N={N_per_regime[r]}) ---")
        save_path = (f"{save_dir.rstrip('/')}/scan_{r}.csv" if save_dir else None)
        df, _ = run_global_lhs_scan(REGIME_PRESETS[r], N_per_regime[r], seed=seed,
                                    wrapper=wrapper, save_path=save_path, verbose=False)
        df['target_preset'] = r
        cluster = df[df['regime'] == r].reset_index(drop=True)
        clusters[r], scans[r] = cluster, df
        yields[r] = len(cluster) / len(df) if len(df) else 0.0
        if verbose:
            print(f"  yield: {len(cluster)}/{len(df)} in-regime ({yields[r]*100:.1f}%)")
            if len(cluster) < min_cluster:
                print(f"  ⚠ '{r}' cluster {len(cluster)} < {min_cluster} — raise N_per_regime['{r}'].")

    master = pd.concat([clusters[r] for r in regimes], ignore_index=True)
    return clusters, scans, master, yields


def load_regime_scans(save_dir, regimes=('metal', 'surface', 'oxide')):
    """
    Reload scans saved by run_targeted_regime_scans(save_dir=...) and rebuild the
    clusters/master structures — so an analysis can reuse cached model evaluations
    instead of re-running the (expensive) scans.

    Returns (clusters, scans, master) matching run_targeted_regime_scans' first three.
    """
    scans, clusters = {}, {}
    for r in regimes:
        path = os.path.join(save_dir, f'scan_{r}.csv')
        df = pd.read_csv(path)
        scans[r] = df
        clusters[r] = df[df['regime'] == r].reset_index(drop=True)
    master = pd.concat([clusters[r] for r in regimes], ignore_index=True)
    return clusters, scans, master


# -----------------------------------------------------------------------------
# Phase 2 Route B — given-data sensitivity on the true clusters (PAWN + delta)
# -----------------------------------------------------------------------------

def _givendata_problem(X, names):
    """Problem dict for given-data analyzers; bounds from the data (epsilon-guarded)."""
    lo, hi = X.min(axis=0), X.max(axis=0)
    bounds = []
    for j in range(len(names)):
        a, b = float(lo[j]), float(hi[j])
        if b <= a:
            b = a + (abs(a) if a != 0 else 1.0) * 1e-9
        bounds.append([a, b])
    return {'num_vars': len(names), 'names': list(names), 'bounds': bounds}


def givendata_sensitivity_by_regime(clusters, param_names, output_metrics=REGIME_SA_METRICS,
                                    methods=('pawn', 'delta'), pawn_S=10,
                                    delta_resamples=100, seed=42,
                                    log_metrics=LOG_METRICS_DEFAULT, verbose=True):
    """
    Route B (primary): for each (regime, metric), run PAWN and Borgonovo-delta on the
    TRUE in-regime cluster — valid on arbitrary scattered points, zero new model evals.

    Returns
    -------
    {regime: {metric: {'n': int, 'pawn': Si_pawn, 'delta': Si_delta}}}
    (an entry carries 'skipped' instead of indices if the cluster is too small/degenerate)
    """
    param_names = list(param_names)
    results = {}
    for reg, cdf in clusters.items():
        X = cdf[param_names].to_numpy(dtype=float)
        problem = _givendata_problem(X, param_names)
        results[reg] = {}
        for m in output_metrics:
            Y = cdf[m].to_numpy(dtype=float)
            use_log = m in log_metrics
            mask = (np.isfinite(Y) & (Y > 0)) if use_log else np.isfinite(Y)
            Xv = X[mask]
            Yv = np.log10(Y[mask]) if use_log else Y[mask]
            entry = {'n': int(mask.sum()), 'log': use_log}
            if entry['n'] < 10 or np.ptp(Yv) == 0:
                entry['skipped'] = 'insufficient or degenerate output'
                results[reg][m] = entry
                if verbose:
                    print(f"  [{reg}/{m}] skipped ({entry['n']} pts)")
                continue
            if 'pawn' in methods:
                entry['pawn'] = pawn_analyzer.analyze(problem, Xv, Yv, S=pawn_S, seed=seed)
            if 'delta' in methods:
                entry['delta'] = delta_analyzer.analyze(problem, Xv, Yv,
                                                        num_resamples=delta_resamples,
                                                        conf_level=0.95, seed=seed,
                                                        print_to_console=False)
            results[reg][m] = entry
            if verbose:
                top = ''
                if 'delta' in entry:
                    d = np.asarray(entry['delta']['delta'])
                    top = param_names[int(np.argmax(d))]
                print(f"  [{reg}/{m}] n={entry['n']:5d}  top-δ param: {top}")
    return results


def summarize_givendata(results, param_names):
    """Flatten Route-B results into a tidy long DataFrame for tables/heatmaps."""
    param_names = list(param_names)
    rows = []
    for reg, by_m in results.items():
        for m, entry in by_m.items():
            if 'skipped' in entry:
                continue
            pawn_med = (np.asarray(entry['pawn']['median']) if 'pawn' in entry
                        else np.full(len(param_names), np.nan))
            delta = (np.asarray(entry['delta']['delta']) if 'delta' in entry
                     else np.full(len(param_names), np.nan))
            s1 = (np.asarray(entry['delta']['S1']) if 'delta' in entry
                  else np.full(len(param_names), np.nan))
            for j, p in enumerate(param_names):
                rows.append({'regime': reg, 'metric': m, 'parameter': p,
                             'pawn_median': pawn_med[j], 'delta': delta[j],
                             'S1_givendata': s1[j], 'n': entry['n']})
    return pd.DataFrame(rows)


def plot_givendata_results(results, regime, output_metric, param_names, top_n=15, show=True):
    """Bar plots of PAWN median and Borgonovo-delta (with conf) for one regime+metric."""
    entry = results[regime][output_metric]
    if 'skipped' in entry:
        print(f"[{regime}/{output_metric}] skipped: {entry['skipped']}")
        return None
    names = list(param_names)
    pawn_med = np.asarray(entry['pawn']['median']) if 'pawn' in entry else None
    delta = np.asarray(entry['delta']['delta']) if 'delta' in entry else None
    delta_conf = np.asarray(entry['delta']['delta_conf']) if 'delta' in entry else None

    order = np.argsort(delta if delta is not None else pawn_med)[::-1][:top_n]
    n_panels = sum(x is not None for x in (pawn_med, delta))
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, max(4, 0.35 * len(order))))
    if n_panels == 1:
        axes = [axes]
    ax_i = 0
    if delta is not None:
        ax = axes[ax_i]; ax_i += 1
        ax.barh([names[i] for i in order][::-1], delta[order][::-1],
                xerr=delta_conf[order][::-1], color='coral', alpha=0.8)
        ax.set_title(f'Borgonovo δ — {regime} / {output_metric}  (n={entry["n"]})')
        ax.set_xlabel('δ (moment-independent)')
    if pawn_med is not None:
        ax = axes[ax_i]
        ax.barh([names[i] for i in order][::-1], pawn_med[order][::-1],
                color='steelblue', alpha=0.8)
        ax.set_title(f'PAWN median — {regime} / {output_metric}')
        ax.set_xlabel('PAWN KS (median)')
    fig.tight_layout()
    if show:
        plt.show()
    return fig


# -----------------------------------------------------------------------------
# Phase 2 Route A — accept-and-report Saltelli over the regime's top-k sub-box
# -----------------------------------------------------------------------------

def sobol_regime_subbox(cluster_df, top_k, target_regime, output_metrics=REGIME_SA_METRICS,
                        N=1024, seed=42, low=5, high=95, fallback_ranges=None,
                        calc_second_order=True, wrapper=None,
                        log_metrics=LOG_METRICS_DEFAULT, cache_csv=None, verbose=True):
    """
    Route A (secondary): Saltelli S1/ST over the regime's top-k percentile sub-box.
    Non-top-k parameters are held at the cluster median (focused-Sobol style). Every
    sample is re-labelled and the off-regime contamination fraction is reported
    (accept-and-report — NO row deletion, which would invalidate the estimator).

    cache_csv : if given and present, the (X, Y, regime_hit) model evaluations are
        loaded instead of re-running the model; otherwise they are saved there after
        running. The Sobol analysis is always recomputed (cheap).

    Returns
    -------
    {'sobol': {metric: {'Si': Si, 'problem': problem, 'log': bool}},
     'contamination': float, 'regime_counts': {regime: n}, 'subbox': {...}, 'n': int}
    """
    if wrapper is None:
        wrapper = level5L6_model_wrapper

    subbox  = extract_topk_subbox(cluster_df, top_k, low, high, fallback_ranges)
    problem = _make_problem(subbox)
    names   = problem['names']

    all_params = list(SUGGESTED_RANGES_LEVEL5L6.keys())
    fixed = {p: float(cluster_df[p].median()) for p in all_params if p not in top_k}

    if seed is not None:
        np.random.seed(seed)
    X = sobol_sampler.sample(problem, N, calc_second_order=calc_second_order)
    n = X.shape[0]

    if cache_csv and os.path.exists(cache_csv):
        cdf = pd.read_csv(cache_csv)
        if all(c in cdf.columns for c in names) and len(cdf) == n:
            X = cdf[names].to_numpy(dtype=float)
        Ys = {m: cdf[f'Y_{m}'].to_numpy(dtype=float) for m in output_metrics
              if f'Y_{m}' in cdf.columns}
        regime_hits = (cdf['regime_hit'].to_numpy() if 'regime_hit' in cdf.columns
                       else np.array([target_regime] * n))
        if verbose:
            print(f"  [cache] loaded Route-A evals ('{target_regime}') from {cache_csv} ({len(cdf)} rows)")
    else:
        if verbose:
            print(f"\n  Route A Sobol — regime '{target_regime}', {len(names)} vars, {n} samples")
        Ys = {m: np.full(n, np.nan) for m in output_metrics}
        regime_hits = []
        for i in range(n):
            params = dict(fixed)
            params.update(dict(zip(names, X[i])))
            rec = wrapper(params, return_full_record=True)
            for m in output_metrics:
                Ys[m][i] = rec.get(m, np.nan)
            regime_hits.append(rec.get('regime', 'undefined'))
            if verbose and ((i + 1) % max(1, n // 5) == 0 or (i + 1) == n):
                print(f"    {i + 1}/{n}")
        regime_hits = np.array(regime_hits)
        if cache_csv:
            d = pd.DataFrame(X, columns=names)
            for m in output_metrics:
                d[f'Y_{m}'] = Ys[m]
            d['regime_hit'] = regime_hits
            d.to_csv(cache_csv, index=False)
            if verbose:
                print(f"  [cache] saved Route-A evals -> {cache_csv}")

    regime_hits = np.asarray(regime_hits)
    regime_counts = {r: int((regime_hits == r).sum()) for r in np.unique(regime_hits)}
    contamination = float(1.0 - (regime_hits == target_regime).mean())

    sobol = {}
    for m in output_metrics:
        Y = _prep_Y_full(Ys[m], m, log_metrics)
        Si = sobol_analyzer.analyze(problem, Y, calc_second_order=calc_second_order,
                                    conf_level=0.95, print_to_console=False)
        sobol[m] = {'Si': Si, 'problem': problem, 'log': m in log_metrics}

    if verbose:
        print(f"  off-regime contamination: {contamination*100:.1f}%  "
              f"(in-regime {regime_counts.get(target_regime, 0)}/{n})")

    return {'sobol': sobol, 'contamination': contamination,
            'regime_counts': regime_counts, 'subbox': subbox, 'n': n}


# -----------------------------------------------------------------------------
# Phase 3 — cross-regime comparison reporting
# -----------------------------------------------------------------------------

def regime_comparison_matrix(givendata_results, output_metric, param_names,
                             index='delta', top_union=12):
    """
    Build a (parameter × regime) matrix of a chosen Route-B index for one metric,
    restricted to the union of each regime's top-`top_union` parameters.

    index : 'delta' | 'pawn_median' | 'S1_givendata'
    Returns the DataFrame (params as rows, regimes as columns).
    """
    df = summarize_givendata(givendata_results, param_names)
    df = df[df['metric'] == output_metric]
    if df.empty:
        return df

    union = set()
    for reg, sub in df.groupby('regime'):
        union |= set(sub.nlargest(top_union, index)['parameter'])

    mat = (df[df['parameter'].isin(union)]
           .pivot(index='parameter', columns='regime', values=index))
    # order rows by overall importance (row mean)
    mat = mat.reindex(mat.mean(axis=1).sort_values(ascending=False).index)
    return mat


def plot_regime_comparison_heatmap(givendata_results, output_metric, param_names,
                                   index='delta', top_union=12, show=True):
    """Heatmap of a Route-B index across regimes (rows=params, cols=regimes)."""
    mat = regime_comparison_matrix(givendata_results, output_metric, param_names,
                                   index=index, top_union=top_union)
    if mat.empty:
        print(f"No data for metric '{output_metric}'.")
        return None
    fig, ax = plt.subplots(figsize=(1.6 * mat.shape[1] + 3, 0.45 * mat.shape[0] + 2))
    im = ax.imshow(mat.values, cmap='YlOrRd', aspect='auto')
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=0)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(f'{index}')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=8,
                        color='black' if v < np.nanmax(mat.values) * 0.6 else 'white')
    ax.set_title(f'Cross-regime sensitivity ({index}) — {output_metric}')
    fig.tight_layout()
    if show:
        plt.show()
    return fig


def contamination_summary(route_a_results):
    """Tidy table of Route-A off-regime contamination + in-regime counts per regime."""
    rows = []
    for reg, res in route_a_results.items():
        rows.append({'regime': reg, 'n': res['n'],
                     'in_regime': res['regime_counts'].get(reg, 0),
                     'contamination_%': round(res['contamination'] * 100, 1)})
    return pd.DataFrame(rows)


def summarize_route_a(route_a_results, param_names=None):
    """
    Flatten Route-A Sobol results into a tidy long DataFrame (one row per
    regime × metric × parameter) for CSV export and plotting from disk.
    Columns: regime, metric, parameter, S1, S1_conf, ST, ST_conf, log, contamination_%.
    """
    rows = []
    for reg, res in route_a_results.items():
        cont = round(res['contamination'] * 100, 1)
        for m, entry in res['sobol'].items():
            Si    = entry['Si']
            names = entry['problem']['names']
            for j, p in enumerate(names):
                rows.append({'regime': reg, 'metric': m, 'parameter': p,
                             'S1': Si['S1'][j], 'S1_conf': Si['S1_conf'][j],
                             'ST': Si['ST'][j], 'ST_conf': Si['ST_conf'][j],
                             'log': entry.get('log', False),
                             'contamination_%': cont})
    return pd.DataFrame(rows)


# =============================================================================
# PARALLEL-COORDINATES PLOTS (Plotly) — interactive views of the regime data
# Plotly is imported lazily so the module never hard-requires it.
# =============================================================================

REGIME_COLOR_CODE = {'surface': 0, 'oxide': 1, 'metal': 2}
# Discrete 3-band colorscale (used with cmin=-0.5, cmax=2.5 so codes 0/1/2 centre in each band)
_REGIME_COLORSCALE = [[0.0, '#1f77b4'], [0.333, '#1f77b4'],
                      [0.334, '#ff7f0e'], [0.666, '#ff7f0e'],
                      [0.667, '#2ca02c'], [1.0, '#2ca02c']]


def top_drivers(routeB_df, regime, metric='flux', k=6, by='delta'):
    """Pick a regime's top-k driver parameters from the Route-B tidy table."""
    sub = routeB_df[(routeB_df['regime'] == regime) & (routeB_df['metric'] == metric)]
    return sub.nlargest(k, by)['parameter'].tolist()


def _pcp_axis(df, col, log_dims):
    """Return (values, label) for one PCP axis; log10-transform if col in log_dims."""
    v = df[col].to_numpy(dtype=float)
    if col in log_dims:
        pos = v > 0
        lv = np.full_like(v, np.nan)
        lv[pos] = np.log10(v[pos])
        if not pos.all():
            lv[~pos] = np.nanmin(lv) if pos.any() else 0.0
        return lv, f'log10({col})'
    return v, col


def parallel_coordinates_samples(df, dimensions, color_by='regime', log_dims=(),
                                 title='', save_html=None,
                                 width=None, height=520, labelangle=30):
    """
    Per-sample parallel-coordinates plot (covers views A/B/D).

    dimensions : list of column names → one axis each (heavy-tailed ones listed in
                 `log_dims` are shown as log10).
    color_by   : 'regime' → discrete colouring (surface/oxide/metal); otherwise a
                 numeric column ('flux'/'permeability' are coloured on log10).
    save_html  : optional path to write a standalone interactive HTML.
    width      : figure width in px; default scales with the number of axes so the
                 axis titles don't overlap. height/labelangle tune layout.
    """
    import plotly.graph_objects as go

    log_dims = set(log_dims)
    dims = [dict(label=lbl, values=vals)
            for vals, lbl in (_pcp_axis(df, d, log_dims) for d in dimensions)]

    if color_by == 'regime':
        codes = df['regime'].map(REGIME_COLOR_CODE).to_numpy(dtype=float)
        line = dict(color=codes, colorscale=_REGIME_COLORSCALE, cmin=-0.5, cmax=2.5,
                    showscale=True,
                    colorbar=dict(title='regime', tickvals=[0, 1, 2],
                                  ticktext=['surface', 'oxide', 'metal']))
    else:
        clog = {color_by} if color_by in ('flux', 'permeability') else set()
        cvals, clbl = _pcp_axis(df, color_by, clog)
        line = dict(color=cvals, colorscale='Viridis', showscale=True,
                    colorbar=dict(title=clbl))

    # Only pin an explicit (wide) width when there are MANY axes (else the titles
    # collide). With few axes, leave width unset so the figure stays responsive and
    # fills the cell/page like before.
    if width is None and len(dims) > 7:
        width = 135 * len(dims)
    fig = go.Figure(data=go.Parcoords(
        line=line, dimensions=dims,
        labelangle=labelangle, labelside='top',
        labelfont=dict(size=11), tickfont=dict(size=9), rangefont=dict(size=9)))
    layout = dict(title=title, font=dict(size=12), margin=dict(l=100, r=120, t=110, b=40))
    if width is not None:
        layout['width'] = width
        layout['height'] = height
    fig.update_layout(**layout)
    if save_html:
        fig.write_html(save_html)
    return fig


def parallel_coordinates_sensitivity(matrix_df, regimes=('metal', 'oxide', 'surface'),
                                     param_col='parameter', title='', save_html=None,
                                     top_n=None):
    """
    Sensitivity-shift PCP (view C): one line per PARAMETER, one axis per regime,
    axis value = the index in `matrix_df` (e.g. δ from compare_delta_flux.csv).

    A leftmost categorical axis labels each line with its parameter name (sorted by
    mean importance, highest at top) so individual lines are identifiable —
    go.Parcoords has no per-line legend/hover. Line colour = mean importance.
    `top_n` keeps only the most important N parameters to reduce clutter.
    """
    import plotly.graph_objects as go

    regs = [r for r in regimes if r in matrix_df.columns]
    df = matrix_df.copy()
    df['_imp'] = df[regs].mean(axis=1)
    df = df.sort_values('_imp', ascending=False)
    if top_n:
        df = df.head(top_n)

    params = df[param_col].astype(str).tolist()
    n = len(params)
    M = df[regs].to_numpy(dtype=float)
    lo, hi = float(np.nanmin(M)), float(np.nanmax(M))

    # Leftmost axis = parameter name (rank position; highest importance at the top),
    # so every line can be traced back to which parameter it is.
    rank = np.arange(n)[::-1]
    param_axis = dict(label=param_col, values=rank,
                      tickvals=list(rank), ticktext=params, range=[-0.5, n - 0.5])
    dims = [param_axis] + [dict(label=r, values=df[r].to_numpy(dtype=float),
                                range=[lo, hi]) for r in regs]
    line = dict(color=df['_imp'].to_numpy(dtype=float), colorscale='YlOrRd',
                showscale=True, colorbar=dict(title='mean<br>importance'))

    fig = go.Figure(data=go.Parcoords(
        line=line, dimensions=dims,
        labelangle=0, labelside='top',
        labelfont=dict(size=12), tickfont=dict(size=9), rangefont=dict(size=9)))
    fig.update_layout(title=title, font=dict(size=12),
                      height=max(500, 26 * n + 150),
                      margin=dict(l=175, r=120, t=80, b=40))
    if save_html:
        fig.write_html(save_html)
    return fig


"""regime_viz.py  —  one function, three modes.

    plot_route_a(clusters, axes=...)                              -> Phase 1 view
    plot_route_a(clusters, axes=..., focus='metal')               -> Route A view
    plot_route_a(clusters, axes=..., focus='metal', animate=True) -> frame sequence
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def plot_route_a(
    clusters,
    *,
    axes,
    focus=None,
    route_a=None,
    results_dir="sa_results",
    evals_df=None,
    regime_col="regime_hit",
    animate=False,
    step=200,
    log_x=False,
    log_y=False,
    figsize=(9, 6.5),
    save_path=None,
    title=None,
):
    """Visualize Route A regime separation in parameter space.
 
    Parameters
    ----------
    clusters : dict[str, pd.DataFrame]
        {regime: cluster_df} from Phase 1.
    axes : (str, str)
        Two parameter names to project on.
    focus : str, optional
        Regime whose Saltelli sample to overlay. None -> Phase 1 view only.
    route_a : dict, optional
        The full route_a result dict. If given, the library's recorded sub-box
        (route_a[r]['subbox']) is used for each regime instead of the cluster
        percentile band. Falls back to percentile bands for regimes / axes the
        sub-box doesn't cover.
    results_dir : str
        Directory holding routeA_evals_<focus>.csv (default 'sa_results').
    evals_df : pd.DataFrame, optional
        Pre-loaded Saltelli evaluations DataFrame; overrides the CSV read.
    regime_col : str
        Column in the eval CSV that holds the post-wrapper regime label
        (default 'regime_hit').
    animate : bool
        Write a frame sequence to <results_dir>/anim/<focus>/frame_NNNNN.png.
    step : int
        Animation frame cadence (rows per frame).
    log_x, log_y : bool
        Log-scale the corresponding axis.
    figsize, save_path, title : usual matplotlib knobs.
 
    Examples
    --------
    >>> plot_route_a(clusters, axes=('k_diss_ref', 'P_upstream'))
    >>> plot_route_a(clusters, axes=regime_topk('metal')[:2],
    ...              focus='metal', route_a=route_a)
    >>> plot_route_a(clusters, axes=regime_topk('metal')[:2],
    ...              focus='metal', animate=True, step=200)
    """
    REGIME_COLORS = {"surface": "#1D9E75", "oxide": "#7F77DD", "metal": "#D85A30"}
    CONTAM_EDGE = "#A32D2D"
    P1, P2 = axes
 
    def percentile_box(cluster, params, lo_q=0.10, hi_q=0.90):
        return {p: (cluster[p].quantile(lo_q), cluster[p].quantile(hi_q))
                for p in params}
 
    def box_bounds_for(regime):
        """Pick the best available (lo, hi) for both axes for one regime.
        Prefer the library's recorded sub-box; fall back to percentile band."""
        out = {}
        sb_lib = None
        if route_a is not None and regime in route_a:
            sb_lib = route_a[regime].get("subbox")
        pb = percentile_box(clusters[regime], [P1, P2])
        for ax_name in (P1, P2):
            if sb_lib is not None and ax_name in sb_lib:
                lo, hi = sb_lib[ax_name]
                out[ax_name] = (float(lo), float(hi))
            else:
                out[ax_name] = pb[ax_name]
        return out
 
    def draw(ax, evals_subset=None, n_total=None):
        bg_alpha = 0.22 if (focus is not None) else 0.45
        bg_size  = 12 if (focus is not None) else 14
        for r, c in REGIME_COLORS.items():
            if r not in clusters:
                continue
            sub = clusters[r]
            ax.scatter(sub[P1], sub[P2], c=c, s=bg_size, alpha=bg_alpha,
                       linewidths=0, label=f"{r} cluster  n={len(sub)}")
            b = box_bounds_for(r)
            (lo1, hi1), (lo2, hi2) = b[P1], b[P2]
            ax.add_patch(Rectangle(
                (lo1, lo2), hi1 - lo1, hi2 - lo2,
                fill=False, edgecolor=c, linestyle="--", linewidth=1.3,
            ))
 
        if evals_subset is not None:
            in_ = evals_subset[regime_col] == focus
            ax.scatter(
                evals_subset.loc[in_, P1], evals_subset.loc[in_, P2],
                facecolors="none", edgecolors="#444",
                s=20, linewidths=0.7,
                label=f"Saltelli -> in-regime  n={int(in_.sum())}",
            )
            for r, c in REGIME_COLORS.items():
                if r == focus:
                    continue
                m = evals_subset[regime_col] == r
                if not m.any():
                    continue
                ax.scatter(
                    evals_subset.loc[m, P1], evals_subset.loc[m, P2],
                    facecolors=c, edgecolors=CONTAM_EDGE,
                    s=30, linewidths=1.3,
                    label=f"Saltelli -> {r} (contam)  n={int(m.sum())}",
                )
            contam_pct = 100 * (~in_).mean()
            n_str = f"n={len(evals_subset)}" + (f"/{n_total}" if n_total else "")
            ax.set_title(title or (
                f"Route A geometry  -  focus = {focus}  "
                f"({P1}, {P2})  .  {n_str}  .  contam {contam_pct:.1f}%"
            ))
        else:
            ax.set_title(title or
                         f"Phase 1 geometry  -  projection on ({P1}, {P2})")
 
        if log_x:
            ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel(P1)
        ax.set_ylabel(P2)
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, alpha=0.15)
 
    # --- Phase 1 view -------------------------------------------------------
    if focus is None:
        fig, ax = plt.subplots(figsize=figsize)
        draw(ax)
        fig.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=120, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()
        return None
 
    # --- Load Saltelli evaluations -----------------------------------------
    if evals_df is None:
        path = os.path.join(results_dir, f"routeA_evals_{focus}.csv")
        evals = pd.read_csv(path)
    else:
        evals = evals_df.copy()
    if regime_col not in evals.columns:
        raise KeyError(
            f"column {regime_col!r} not in eval CSV. Found: "
            f"{list(evals.columns)}. Pass regime_col=<name> if your library "
            f"uses a different column for the post-wrapper regime label."
        )
 
    # --- Animation mode ----------------------------------------------------
    if animate:
        outdir = os.path.join(results_dir, "anim", focus)
        os.makedirs(outdir, exist_ok=True)
        n = len(evals)
        ends = list(range(step, n + 1, step))
        if not ends or ends[-1] != n:
            ends.append(n)
        for end in ends:
            fig, ax = plt.subplots(figsize=figsize)
            draw(ax, evals_subset=evals.iloc[:end], n_total=n)
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, f"frame_{end:05d}.png"), dpi=110)
            plt.close(fig)
        print(f"wrote {len(ends)} frames to {outdir}")
        print(f"stitch:  ffmpeg -framerate 8 -i "
              f"{outdir}/frame_%05d.png anim_{focus}.gif")
        return outdir
 
    # --- Static Route A view -----------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    draw(ax, evals_subset=evals)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
    return None
 