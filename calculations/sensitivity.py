import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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

def level5_model_wrapper(params_dict, return_full_record=False):
    """
    Wrapper for LEVEL 5 complete hierarchical model (L3 + L4).

    Combines defective oxide (L3) and defective metal microstructure (L4).
    Uses reference-point Arrhenius: X(T) = X_ref * exp(-E/R * (1/T - 1/T_ref)).

    Parameters
    ----------
    params_dict : dict
        Any subset of parameters from DEFAULT_PARAMS_LEVEL5.
        Missing parameters are filled from defaults.
    return_full_record : bool
        When True, also attach the string diagnostics 'regime_hierarchy' and
        'dominant_path'. Used by the regime-stratified scan (run_global_lhs_scan).

    Returns
    -------
    dict
        - 'flux': Total permeation flux [mol/m²/s]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
        - 'PRF': Permeation Reduction Factor [-]
        - 'D_eff', 'D_modification', 'P_interface'
        - 'flux_intact', 'flux_defect', 'flux_bare_metal'
        - 'frac_oxide', 'frac_metal', 'frac_defect': regime fractions (sum ~1)
        - 'regime': argmax label — 'oxide' | 'metal' | 'defect' | 'undefined'

    Note
    ----
    'regime' holds the argmax label used to build the regime clusters. The older
    hierarchical L3/L4 string (e.g. 'metal_limited/traps_defect_limited') moved to
    'regime_hierarchy' so that partition_by_regime() can key on 'regime' unchanged.
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

        # --- Regime fractions -------------------------------------------------
        # Parallel axis: what share of the flux bypasses the intact oxide.
        # Series axis: how the driving pressure splits across oxide and metal.
        # Taken straight from the solved interface pressure, so it is exact — no
        # linearisation of the sqrt-P metal law (calculate_metal_resistance's
        # Taylor form is only valid for small dP, and these sweeps run P_down=0).
        # The series split is then scaled by the share of flux that actually
        # travels the series path, so all three fractions sum to 1 and the
        # parallel axis competes directly against the two series terms.
        f_int = result_l5['flux_intact_contribution']
        f_def = result_l5['flux_defect_contribution']
        f_tot = f_int + f_def
        frac_defect = f_def / f_tot if f_tot > 0 else np.nan

        P_int = result_l5.get('P_interface_intact', np.nan)
        if P_upstream > 0 and np.isfinite(P_int) and np.isfinite(frac_defect):
            # clip guards the rare solver overshoot P_int > P_up
            drop_ox    = float(np.clip((P_upstream - P_int) / P_upstream, 0.0, 1.0))
            series     = 1.0 - frac_defect
            frac_oxide = drop_ox * series
            frac_metal = (1.0 - drop_ox) * series
        else:
            frac_oxide = frac_metal = np.nan

        record = {
            'flux':         flux_total,
            'PRF':          PRF,
            'D_eff':        D_eff,
            'D_modification': D_eff / D_metal if D_metal > 0 else 1.0,
            'permeability': permeability,
            'P_interface':  result_l5.get('P_interface_intact', 0),
            'flux_intact':  f_int,
            'flux_defect':  f_def,
            'flux_bare_metal': flux_bare,
            'frac_oxide':   frac_oxide,
            'frac_metal':   frac_metal,
            'frac_defect':  frac_defect,
            'regime':       assign_regime_L5(frac_oxide, frac_metal, frac_defect),
            'D_metal':      D_metal,
            'K_s_metal':    K_s_metal,
            'D_ox':         D_ox,
            'K_ox':         K_ox,
            'temperature':  T,
            'modification_factor': result_l5.get('modification_factor', 1.0),
            'defect_enhancement': result_l5.get('defect_enhancement_factor', 1.0),
        }
        if return_full_record:
            record['regime_hierarchy'] = result_l5.get('regime', 'unknown')
            record['dominant_path']    = result_l5.get('dominant_path', 'unknown')
        return record

    except Exception as e:
        print(f"Error in Level 5 model: {e}")
        import traceback; traceback.print_exc()
        record = {
            'flux': 1e-20, 'PRF': 1.0, 'D_eff': 1e-12, 'D_modification': 1.0,
            'permeability': 1e-20, 'P_interface': 0,
            'flux_intact': 1e-20, 'flux_defect': 0, 'flux_bare_metal': 1e-20,
            'frac_oxide': np.nan, 'frac_metal': np.nan, 'frac_defect': np.nan,
            'regime': 'undefined',
            'D_metal': 1e-12, 'K_s_metal': 1e-6, 'D_ox': 1e-15, 'K_ox': 1e-10,
            'temperature': full_params.get('temperature', 1073.15),
            'modification_factor': 1.0, 'defect_enhancement': 1.0,
        }
        if return_full_record:
            record['regime_hierarchy'] = 'error'
            record['dominant_path']    = 'error'
        return record


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

REGIME_LABELS    = ('surface', 'oxide', 'metal')   # L5L6, surface kinetics
REGIME_LABELS_L5 = ('oxide', 'metal', 'defect')    # L5, no surface kinetics


def _argmax_label(fracs, rule='argmax', threshold=0.5):
    """Shared labelling rule behind assign_regime / assign_regime_L5.

    Kept in one place so the two model levels cannot drift apart on how a set of
    fractions becomes a regime string.

    Parameters
    ----------
    fracs : dict[str, float]
        {label: fraction}; should sum to ~1 on a valid solve.
    rule : {'argmax', 'threshold'}
    threshold : float
        Dominance cutoff, used only when rule='threshold'.

    Returns
    -------
    str
        The winning label, 'mixed' (threshold rule, no clear winner), or
        'undefined' when any fraction is NaN/None — i.e. a failed solve, which
        partition_by_regime() then drops from the clusters.
    """
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


def assign_regime_L5(frac_oxide, frac_metal, frac_defect, rule='argmax', threshold=0.5):
    """
    Regime label for the LEVEL 5 model (no surface kinetics).

    Without a surface step there are only two series resistances, so the third
    competitor is the *parallel* one: flux that bypasses the intact oxide through
    pinholes/cracks/grain-boundary paths. The three fractions come from
    level5_model_wrapper and sum to 1 by construction.

    Parameters
    ----------
    frac_oxide, frac_metal : float
        Share of the driving pressure dropped across the oxide / the metal,
        scaled by the fraction of flux taking the intact (series) path.
    frac_defect : float
        J_defect / J_total — the parallel bypass share.

    Returns
    -------
    str
        One of 'oxide', 'metal', 'defect' (or 'mixed' for rule='threshold'),
        or 'undefined' if any fraction is NaN/None.

    See Also
    --------
    assign_regime : the L5L6 counterpart (surface / oxide / metal).
    """
    return _argmax_label(
        {'oxide': frac_oxide, 'metal': frac_metal, 'defect': frac_defect},
        rule=rule, threshold=threshold)


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
    return _argmax_label(
        {'surface': frac_surface, 'oxide': frac_oxide, 'metal': frac_metal},
        rule=rule, threshold=threshold)


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
        scan (run_global_lhs_scan). Default False returns the scalar metrics only.

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
# REGIME-STRATIFIED SENSITIVITY ANALYSIS — LEVEL 5L6
# Global LHS foundation + per-regime given-data SA (Route B).
# Argmax regime labelling; PAWN + Borgonovo-delta on the true in-regime clusters,
# which need no structured design and cost no extra model evaluations.
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
# density estimators behind PAWN/δ aren't dominated by a few huge values (raw
# flux ranges ~10 decades). theta (∈[0,1]) stays linear.
LOG_METRICS_DEFAULT = ('flux', 'permeability')


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
                        output_fields=None, label_fields=None,
                        fixed_params=None, save_path=None, verbose=True):
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
    if output_fields is None:
        output_fields = SCAN_OUTPUT_FIELDS
    if label_fields is None:
        label_fields = SCAN_LABEL_FIELDS
    fixed_params = dict(fixed_params or {})

    clash = set(fixed_params) & set(param_ranges)
    if clash:
        raise ValueError(f"{sorted(clash)} is both sampled and fixed — drop it from "
                         f"param_ranges (see presets_without) so it is not an SA input")

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
        rec    = wrapper({**params, **fixed_params}, return_full_record=True)
        row    = dict(params)  # ONLY the sampled values are SA inputs
        # fixed values are recorded for provenance but must stay out of the SA
        # parameter list — a constant column has zero variance, which PAWN and
        # delta cannot condition on.
        row.update(fixed_params)
        for f in output_fields:
            row[f] = rec.get(f, np.nan)
        for f in label_fields:
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
                            scan_df=None, cross_col=None, show=True):
    """
    Per-regime descriptive view: cluster sizes + output distributions by regime,
    and (if scan_df given) a cross-tab of the argmax regime against the solver's
    own categorical verdict.

    cross_col : str, optional
        Column holding that verdict. Defaults to whichever of
        'system_rate_limiting' (L5L6) / 'regime_hierarchy' (L5) is present.
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

    # Agreement check: our argmax label vs the solver's own categorical verdict.
    # L5L6 calls that 'system_rate_limiting'; L5 calls it 'regime_hierarchy'.
    if scan_df is not None:
        col = next((c for c in (cross_col, 'system_rate_limiting', 'regime_hierarchy')
                    if c and c in scan_df.columns), None)
        if col:
            valid = scan_df[scan_df['regime'] != 'undefined']
            print(f"\n  regime (argmax)  ×  {col} (solver):")
            print(pd.crosstab(valid['regime'], valid[col]).to_string())

    return fig




# -----------------------------------------------------------------------------
# Regime-targeted sampling presets
# Each preset = the default ranges with ONLY the regime-controlling parameters
# overridden, left as broad as possible so the per-regime SA reflects genuine
# sensitivity rather than the targeting itself. Targeting is necessary because an
# untargeted sweep over the full ranges yields essentially no surface-limited rows.
#
# Yields measured on the 36-param production scans (Application/sa_results/scans):
#   metal    813/1000  = 81.3%   (default ranges)
#   surface  714/2500  = 28.6%   (low pressure + slow surface kinetics)
#   oxide   3725/5000  = 74.5%   (suppressed defects + slow/thick oxide + fast surface & metal)
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

# Default LHS draws per regime = target cluster size / that preset's measured yield,
# so all three clusters land near 1500 points. Balance matters as much as absolute
# size here: delta/PAWN estimator variance depends on n, so unequal clusters make the
# cross-regime comparison in regime_comparison_matrix() compare noise levels as well
# as sensitivities. Column-normalisation in the heatmap rescales magnitude but cannot
# undo that. 1500 also leaves headroom over the ~300-500 stability floor.
DEFAULT_N_PER_REGIME = {'metal': 1850, 'surface': 5250, 'oxide': 2010}


# -----------------------------------------------------------------------------
# LEVEL 5 (no surface kinetics) — scan fields, metrics and regime-targeted presets
#
# The L5 regimes are oxide / metal / defect: two series resistances plus the
# PARALLEL bypass through pinholes, cracks and grain-boundary paths. There is no
# 'theta' (no surface coverage) and no surface-kinetics parameters to target with,
# so the oxide preset has to lean entirely on oxide/metal transport and geometry.
#
# Measured on 300-draw LHS probes (seed 42):
#   default ranges   -> 98% metal, max frac_defect 0.415  (i.e. ONE regime;
#                       targeting is not optional here, it is what makes the study)
#   metal preset     -> 97.0% in-regime
#   oxide preset     -> 77.7%
#   defect preset    -> 87.3%,  max frac_defect 0.999
# -----------------------------------------------------------------------------

SCAN_OUTPUT_FIELDS_L5 = [
    'flux', 'permeability', 'PRF',
    'frac_oxide', 'frac_metal', 'frac_defect',
    'D_eff', 'D_modification', 'P_interface',
    'flux_intact', 'flux_defect', 'flux_bare_metal',
]
SCAN_LABEL_FIELDS_L5 = ['regime', 'regime_hierarchy', 'dominant_path']

# 'flux' is primary. 'PRF' (bare-metal flux / coated flux) is the coating-
# effectiveness metric and replaces L5L6's 'theta'. NOTE 'permeability' here is a
# harmonic mean of the oxide and metal permeabilities only — it ignores the defect
# paths entirely, so it is bulk-only by construction and will look insensitive to
# the very parameters that define the 'defect' regime. Kept for continuity with
# L5L6, but read 'flux' and 'PRF' for the defect cluster.
REGIME_SA_METRICS_L5 = ['flux', 'permeability', 'PRF']
LOG_METRICS_L5       = ('flux', 'permeability', 'PRF')   # PRF spans decades too


def _ranges_with_L5(overrides):
    R = {k: list(v) for k, v in SUGGESTED_RANGES_LEVEL5.items()}
    R.update({k: list(v) for k, v in overrides.items()})
    return R


REGIME_PRESETS_L5 = {
    # default ranges are already ~97% metal-limited
    'metal': _ranges_with_L5({}),

    'oxide': _ranges_with_L5({
        'f_pinhole':   [1e-7, 1e-5],    # suppress the bypass so flux uses the oxide
        'f_crack':     [1e-7, 1e-5],
        'f_gb_defect': [1e-7, 1e-5],
        'oxide_thickness': [1e-7, 1e-5],      # thick, slow oxide -> the bottleneck
        'D_ox_ref':        [7.8e-21, 7.8e-20],
        'metal_thickness': [5e-5, 5e-4],      # thin, fast metal -> not limiting
        'D_ref':           [4.3e-8, 4.3e-7],
    }),

    'defect': _ranges_with_L5({
        # Large defect area AND a slow intact oxide, so the bypass actually wins.
        # NOTE this is a heavily degraded coating (total defect area up to ~0.9) —
        # a wider claim than the <=0.3 the default ranges imply. State it in any
        # write-up rather than leaving it buried here.
        'f_pinhole':   [0.05, 0.30],
        'f_crack':     [0.05, 0.30],
        'f_gb_defect': [0.05, 0.30],
        'oxide_thickness': [1e-7, 1e-5],
        'D_ox_ref':        [7.8e-21, 7.8e-20],
    }),
}

# draws = target cluster size / measured yield, balanced at ~1500 per cluster so
# delta/PAWN estimator noise is comparable across the regime columns (see the
# note on DEFAULT_N_PER_REGIME above).
DEFAULT_N_PER_REGIME_L5 = {'metal': 1550, 'oxide': 1930, 'defect': 1720}


def presets_without(presets, *drop):
    """Copy of `presets` with the named parameters removed from every range dict.

    Used to take a parameter OUT of the sampled set so it can be pinned via
    run_global_lhs_scan(fixed_params=...) instead. The motivating case is
    `temperature`: it enters D, K_s, D_ox and K_ox exponentially over 573-1273 K
    and so monopolises the variance of log10(flux) — a dummy-parameter test on the
    varying-T clusters put every other parameter at or below the noise floor
    (delta ~0.09), while temperature scored ~0.5. Pinning it lets the remaining
    parameters compete; in a T-banded probe `H_sol_ox` then reached 5x the floor in
    the oxide regime and `f_crack` surfaced in the defect regime.

    >>> presets_without(REGIME_PRESETS_L5, 'temperature')
    """
    return {r: {k: v for k, v in rng.items() if k not in drop}
            for r, rng in presets.items()}


# Temperatures for the isothermal L5 sweep. Spans the operating window; the probe
# showed the driver ranking genuinely shifts between the low and high ends.
DEFAULT_SWEEP_TEMPERATURES = (773.0, 1073.0, 1273.0)

# In-regime yield of each preset AT FIXED TEMPERATURE (250-draw LHS probes, seed 42).
# These differ from the varying-T yields, so draw counts must be sized per temperature —
# a single dict would leave the oxide cluster ~40% short at 1273 K and reintroduce the
# cross-regime imbalance the isothermal design exists to remove.
#
# The oxide yield falls steadily with temperature (0.92 -> 0.59) because the oxide's
# activation energies (E_D_ox ~70 kJ/mol, H_sol_ox up to ~164 kJ/mol) exceed the metal's,
# so it speeds up faster than the metal and stops being the bottleneck.
MEASURED_YIELDS_L5_BY_T = {
    773.0:  {'oxide': 0.916, 'metal': 0.992, 'defect': 0.916},
    1073.0: {'oxide': 0.700, 'metal': 0.988, 'defect': 0.820},
    1273.0: {'oxide': 0.592, 'metal': 0.972, 'defect': 0.800},
}


def run_targeted_regime_scans(N_per_regime=None, regimes=('metal', 'surface', 'oxide'),
                              seed=42, wrapper=None, save_dir=None,
                              presets=None, min_cluster=300, verbose=True,
                              **scan_kw):
    """
    Phase 1: one LHS scan per regime over its targeted preset, keeping
    the in-regime rows (argmax label). Feeds Route B (given-data).

    Parameters
    ----------
    presets : {regime: {param: [lo, hi]}}, optional
        Targeted ranges per regime. Defaults to REGIME_PRESETS (L5L6); pass
        REGIME_PRESETS_L5 together with regimes=REGIME_LABELS_L5 and
        wrapper=level5_model_wrapper to run the no-surface model instead.
    **scan_kw
        Forwarded to run_global_lhs_scan — e.g. output_fields / label_fields,
        which differ between the two model levels.

    Returns
    -------
    clusters : {regime: cluster_df}   in-regime rows only (regime == target preset)
    scans    : {regime: full_scan_df} every draw from that preset (for yield diagnostics)
    master   : DataFrame              all clusters concatenated (+ 'target_preset' col)
    yields   : {regime: in_regime_fraction}
    """
    if presets is None:
        presets = REGIME_PRESETS
    if N_per_regime is None:
        N_per_regime = dict(DEFAULT_N_PER_REGIME)
    if isinstance(N_per_regime, int):
        N_per_regime = {r: N_per_regime for r in regimes}

    missing = [r for r in regimes if r not in presets]
    if missing:
        raise KeyError(f"no preset for regime(s) {missing}; presets has {list(presets)}")

    clusters, scans, yields = {}, {}, {}
    for r in regimes:
        if verbose:
            print(f"\n--- targeted scan: {r} (N={N_per_regime[r]}) ---")
        save_path = (f"{save_dir.rstrip('/')}/scan_{r}.csv" if save_dir else None)
        df, _ = run_global_lhs_scan(presets[r], N_per_regime[r], seed=seed,
                                    wrapper=wrapper, save_path=save_path,
                                    verbose=False, **scan_kw)
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
                                    log_metrics=LOG_METRICS_DEFAULT, n_dummy=3,
                                    verbose=True):
    """
    Route B (primary): for each (regime, metric), run PAWN and Borgonovo-delta on the
    TRUE in-regime cluster — valid on arbitrary scattered points, zero new model evals.

    n_dummy : int
        Number of DUMMY inputs (uniform random columns the model never saw) appended
        to X purely to measure the estimator's noise floor. delta is biased positive
        for an irrelevant input, and that bias does NOT shrink with n — measured flat
        at ~0.09 from n=400 to n=1500 — so without this an entire table of noise reads
        as a ranking. Each entry gains 'floor' = max dummy delta; treat any parameter
        not clearly above it as unresolved. Set 0 to disable.

    Returns
    -------
    {regime: {metric: {'n': int, 'floor': float, 'pawn': Si_pawn, 'delta': Si_delta}}}
    (an entry carries 'skipped' instead of indices if the cluster is too small/degenerate)
    """
    param_names = list(param_names)
    dummy_names = [f'_dummy_{i}' for i in range(n_dummy)]
    results = {}
    for reg, cdf in clusters.items():
        X = cdf[param_names].to_numpy(dtype=float)
        if n_dummy:
            rng = np.random.default_rng(seed)
            X = np.column_stack([X, rng.uniform(0.0, 1.0, (len(X), n_dummy))])
        problem = _givendata_problem(X, param_names + dummy_names)
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
            # Noise floor = the best score any dummy achieved. Indices for the real
            # parameters are trimmed back to param_names so every downstream
            # consumer still sees exactly the parameters it asked for.
            if n_dummy and 'delta' in entry:
                d_all = np.asarray(entry['delta']['delta'])
                entry['floor'] = float(np.max(d_all[len(param_names):]))
            for meth, keys in (('pawn', ('median', 'mean', 'minimum', 'CV')),
                               ('delta', ('delta', 'delta_conf', 'S1', 'S1_conf'))):
                if n_dummy and meth in entry:
                    for k in keys:
                        if k in entry[meth]:
                            entry[meth][k] = np.asarray(entry[meth][k])[:len(param_names)]
                    entry[meth]['names'] = list(param_names)

            results[reg][m] = entry
            if verbose:
                top, floor = '', entry.get('floor', np.nan)
                if 'delta' in entry:
                    d = np.asarray(entry['delta']['delta'])
                    top = f"{param_names[int(np.argmax(d))]} ({np.max(d):.3f})"
                print(f"  [{reg}/{m}] n={entry['n']:5d}  floor={floor:.3f}  top-δ: {top}")
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
            floor = entry.get('floor', np.nan)
            for j, p in enumerate(param_names):
                rows.append({'regime': reg, 'metric': m, 'parameter': p,
                             'pawn_median': pawn_med[j], 'delta': delta[j],
                             'S1_givendata': s1[j], 'n': entry['n'],
                             # dummy-parameter noise floor and this parameter's
                             # margin over it — <=1 means "not resolved"
                             'floor': floor,
                             'delta_over_floor': (delta[j] / floor
                                                  if floor and np.isfinite(floor)
                                                  else np.nan)})
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
# Phase 3 — cross-regime comparison reporting
# -----------------------------------------------------------------------------

def _column_normalize(mat):
    """
    Divide every column of `mat` by its own maximum, so each column peaks at 1.0.

    Columns whose max is non-finite or <= 0 are left untouched — max-division is
    meaningless there (relevant for 'S1_givendata', which can come out slightly
    negative on scattered data).
    """
    out = mat.astype(float).copy()
    for c in out.columns:
        col = out[c].to_numpy(dtype=float)
        mx = np.nanmax(col) if np.isfinite(col).any() else np.nan
        if np.isfinite(mx) and mx > 0:
            out[c] = out[c] / mx
    return out


def regime_comparison_matrix(givendata_results, output_metric, param_names,
                             index='delta', top_union=12, normalize=None):
    """
    Build a (parameter × regime) matrix of a chosen Route-B index for one metric,
    restricted to the union of each regime's top-`top_union` parameters.

    index : 'delta' | 'pawn_median' | 'S1_givendata'
    normalize : None (raw values) | 'column' (divide each column by its own max,
        so every regime peaks at 1.0 and rows read as *relative* importance
        within that regime). The top-`top_union` selection always uses the raw
        values, so the same parameters appear either way; only the row ordering
        (by row mean) shifts, since it is computed on whatever is returned.

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
    if normalize == 'column':
        mat = _column_normalize(mat)
    elif normalize not in (None, False):
        raise ValueError(f"normalize must be None or 'column', got {normalize!r}")
    # order rows by overall importance (row mean)
    mat = mat.reindex(mat.mean(axis=1).sort_values(ascending=False).index)
    return mat


def plot_regime_comparison_heatmap(givendata_results, output_metric, param_names,
                                   index='delta', top_union=12, normalize='column',
                                   show_raw=True, show=True):
    """
    Heatmap of a Route-B index across regimes (rows=params, cols=regimes).

    normalize : 'column' (default) colours each cell by its column-normalized
        value, so every regime spans a full 0->1 scale and within-regime contrast
        is visible. Without this, the dominant parameter (usually `temperature`)
        saturates the single global colour scale and everything else reads white.
        Pass None for the raw, un-normalized heatmap.
    show_raw : when normalizing, also print the raw value in a small box in each
        cell's upper-right corner, so absolute magnitude is not lost. The large
        centred number is the normalized value.
    """
    raw = regime_comparison_matrix(givendata_results, output_metric, param_names,
                                   index=index, top_union=top_union)
    if raw.empty:
        print(f"No data for metric '{output_metric}'.")
        return None

    normalized = (normalize == 'column')
    if normalized:
        colour_mat = _column_normalize(raw)
        # re-order rows on the normalized values, then keep raw in lockstep
        colour_mat = colour_mat.reindex(
            colour_mat.mean(axis=1).sort_values(ascending=False).index)
        raw = raw.reindex(colour_mat.index)
        vmin, vmax = 0.0, 1.0
    else:
        colour_mat = raw
        vmin = vmax = None
        show_raw = False

    # dual-annotated cells need a little more room for the corner box
    col_w, row_h = (1.9, 0.55) if show_raw else (1.6, 0.45)
    fig, ax = plt.subplots(figsize=(col_w * colour_mat.shape[1] + 3,
                                    row_h * colour_mat.shape[0] + 2))
    im = ax.imshow(colour_mat.values, cmap='YlOrRd', aspect='auto',
                   vmin=vmin, vmax=vmax)
    ax.set_xticks(range(colour_mat.shape[1]))
    ax.set_xticklabels(colour_mat.columns, rotation=0)
    ax.set_yticks(range(colour_mat.shape[0]))
    ax.set_yticklabels(colour_mat.index)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(f'{index} (fraction of column max)' if normalized else f'{index}')

    flip = (vmax if normalized else np.nanmax(colour_mat.values)) * 0.6
    for i in range(colour_mat.shape[0]):
        for j in range(colour_mat.shape[1]):
            v = colour_mat.values[i, j]
            if not np.isfinite(v):
                continue
            ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=9,
                    color='black' if v < flip else 'white')
            if show_raw:
                r = raw.values[i, j]
                if np.isfinite(r):
                    # translucent white box keeps black text legible on any fill
                    ax.text(j + 0.45, i - 0.39, f'{r:.3f}', ha='right', va='top',
                            fontsize=6.5, color='black',
                            bbox=dict(boxstyle='round,pad=0.15', fc='white',
                                      ec='none', alpha=0.65))

    title = f'Cross-regime sensitivity ({index}) — {output_metric}'
    if normalized:
        title += ('\ncolour + large number = normalized to column max'
                  + ('; corner = raw value' if show_raw else ''))
    ax.set_title(title)
    fig.tight_layout()
    if show:
        plt.show()
    return fig


# =============================================================================
# PARALLEL-COORDINATES PLOTS (Plotly) — interactive views of the regime data
# Plotly is imported lazily so the module never hard-requires it.
# =============================================================================

# Canonical ordering/colour per regime across BOTH model levels: 'surface' only
# occurs at L5L6, 'defect' only at L5.
REGIME_COLOR_CODE   = {'surface': 0, 'oxide': 1, 'metal': 2, 'defect': 3}
_REGIME_BAND_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']


def _regime_colorbar(labels):
    """Discrete colorscale + colorbar ticks for exactly the regimes present.

    Codes are compacted to 0..n-1 rather than reusing REGIME_COLOR_CODE's absolute
    codes: an L5 run ('oxide','metal','defect' -> 1,2,3) would otherwise leave an
    empty band at 0 and shift every colour off its label.

    Returns (code_map, colorscale, tickvals, ticktext); pair with
    cmin=-0.5, cmax=len(tickvals)-0.5 so each code centres in its band.
    """
    present = [r for r in REGIME_COLOR_CODE if r in set(labels)]
    n = max(len(present), 1)
    code, scale = {}, []
    for i, r in enumerate(present):
        code[r] = i
        c = _REGIME_BAND_COLORS[REGIME_COLOR_CODE[r] % len(_REGIME_BAND_COLORS)]
        scale += [[i / n, c], [(i + 1) / n - 1e-9, c]]
    if scale:
        scale[-1][0] = 1.0
    return code, scale, list(range(n)), present


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
        cmap, cscale, tickvals, ticktext = _regime_colorbar(df['regime'])
        codes = df['regime'].map(cmap).to_numpy(dtype=float)
        line = dict(color=codes, colorscale=cscale,
                    cmin=-0.5, cmax=len(tickvals) - 0.5, showscale=True,
                    colorbar=dict(title='regime', tickvals=tickvals,
                                  ticktext=ticktext))
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


# =============================================================================
# CONDITIONAL SWEEP — mechanistic complement to the given-data indices
#
# delta/PAWN answer "how much of the output spread does this input account for
# while everything else also varies". They cannot answer "if I increase this one
# parameter, does the output move, and by how much" — and Puy, Lo Piano & Saltelli
# (2020) show PAWN in particular cannot separate an input whose effect is purely
# interactive from a genuinely inert one. The conditional sweep answers the second
# question directly, so a below-floor index can be checked rather than assumed.
#
# Baselines are drawn FROM the regime cluster, not from one nominal point: a single
# one-at-a-time sweep explores a vanishing fraction of the space (Saltelli & Annoni
# 2010). The SPREAD of the response across baselines is the informative part — a
# wide spread IS the interaction, made visible.
#
# This is conditional, not a global index. Never report it as a sensitivity index.
# =============================================================================

def _sweep_grid(lo, hi, n):
    """Log grid when the range spans a decade or more, linear otherwise."""
    lo, hi = float(lo), float(hi)
    if lo > 0 and hi / lo >= 10.0:
        return np.logspace(np.log10(lo), np.log10(hi), n)
    return np.linspace(lo, hi, n)


def conditional_sweep(cluster_df, param_names, ranges, wrapper=None, fixed_params=None,
                      n_baselines=25, n_grid=9, metric='flux', regime_col='regime',
                      seed=42, verbose=True):
    """
    Sweep each parameter across `ranges` from baselines drawn out of `cluster_df`.

    Parameters
    ----------
    cluster_df : DataFrame
        One regime's cluster; supplies the baseline parameter vectors.
    param_names : list[str]
        Parameters to sweep (one sweep each, all others held at the baseline).
    ranges : {param: [lo, hi]}
        Bounds to sweep over. Pass the global ranges to ask "what if this took any
        plausible value", or the regime's own preset to stay inside the regime —
        the two answer different questions, so run both if it matters.
    fixed_params : dict, optional
        Merged into every model call (e.g. a pinned temperature) — kept out of the
        swept set.
    n_baselines, n_grid : int
        Baselines per parameter, and grid points per sweep. Cost is
        len(param_names) * n_baselines * n_grid model calls.

    Returns
    -------
    DataFrame, one row per (parameter, baseline):
        swing        log10(max metric / min metric) across the grid — the response size
        n_regimes    distinct regime labels seen along the sweep (>1 = crossed a boundary)
        regime_start / regime_end
        n_valid      grid points that produced a usable value
    """
    if wrapper is None:
        wrapper = level5_model_wrapper
    fixed_params = dict(fixed_params or {})
    param_names = [p for p in param_names if p in ranges]

    rng = np.random.default_rng(seed)
    idx = rng.choice(len(cluster_df), size=min(n_baselines, len(cluster_df)),
                     replace=False)
    baselines = cluster_df.iloc[idx].reset_index(drop=True)
    base_cols = [c for c in ranges if c in baselines.columns]

    rows = []
    for k, p in enumerate(param_names):
        grid = _sweep_grid(*ranges[p], n_grid)
        for b in range(len(baselines)):
            base = {c: float(baselines.loc[b, c]) for c in base_cols}
            vals, regs = [], []
            for g in grid:
                rec = wrapper({**base, **{p: float(g)}, **fixed_params},
                              return_full_record=True)
                y = rec.get(metric, np.nan)
                if np.isfinite(y) and y > 0:
                    vals.append(y)
                    regs.append(rec.get(regime_col, 'undefined'))
            if len(vals) >= 2:
                swing = float(np.log10(max(vals) / min(vals)))
            else:
                swing = np.nan
            rows.append({'parameter': p, 'baseline': b, 'swing': swing,
                         'n_regimes': len(set(regs)) if regs else 0,
                         'regime_start': regs[0] if regs else 'undefined',
                         'regime_end': regs[-1] if regs else 'undefined',
                         'n_valid': len(vals)})
        if verbose:
            print(f"    swept {k+1}/{len(param_names)}: {p}")
    return pd.DataFrame(rows)


def summarize_conditional_sweep(sweep_df):
    """Median/IQR swing per parameter, plus how often the sweep left its regime.

    Sort by `swing_med`. A large `swing_iqr` relative to `swing_med` means the
    response depends strongly on the rest of the parameters — i.e. the effect is
    interaction-driven, which is exactly what a global index can miss.
    """
    g = sweep_df.groupby('parameter')['swing']
    out = pd.DataFrame({
        'swing_med': g.median(),
        'swing_iqr': g.quantile(0.75) - g.quantile(0.25),
        'swing_max': g.max(),
        'frac_left_regime': sweep_df.groupby('parameter')['n_regimes']
                                    .apply(lambda s: float((s > 1).mean())),
    })
    return out.sort_values('swing_med', ascending=False)


# =============================================================================
# REGIME GEOMETRY PLOT — cluster projection in parameter space
# =============================================================================

def plot_regime_geometry(
    clusters,
    *,
    axes,
    log_x=False,
    log_y=False,
    figsize=(9, 6.5),
    save_path=None,
    title=None,
):
    """Project the Phase 1 regime clusters onto two parameters.

    Each regime's points are scattered in its own colour and outlined with the
    10-90 percentile band of the two plotted parameters, so cluster overlap (and
    the curvature of the true regime boundary) is visible.

    Parameters
    ----------
    clusters : dict[str, pd.DataFrame]
        {regime: cluster_df} from run_targeted_regime_scans / load_regime_scans.
    axes : (str, str)
        Two parameter names to project on.
    log_x, log_y : bool
        Log-scale the corresponding axis.
    figsize, save_path, title : usual matplotlib knobs.

    Examples
    --------
    >>> plot_regime_geometry(clusters, axes=('k_diss_ref', 'P_upstream'), log_x=True)
    >>> plot_regime_geometry(clusters, axes=regime_topk('metal')[:2])
    """
    from matplotlib.patches import Rectangle

    # 'surface' only occurs at L5L6, 'defect' only at L5; regimes absent from
    # `clusters` are skipped, so one map serves both levels.
    REGIME_COLORS = {"surface": "#1D9E75", "oxide": "#7F77DD",
                     "metal": "#D85A30", "defect": "#B5179E"}
    P1, P2 = axes

    fig, ax = plt.subplots(figsize=figsize)
    for r, c in REGIME_COLORS.items():
        if r not in clusters:
            continue
        sub = clusters[r]
        ax.scatter(sub[P1], sub[P2], c=c, s=14, alpha=0.45, linewidths=0,
                   label=f"{r} cluster  n={len(sub)}")
        lo1, hi1 = sub[P1].quantile(0.10), sub[P1].quantile(0.90)
        lo2, hi2 = sub[P2].quantile(0.10), sub[P2].quantile(0.90)
        ax.add_patch(Rectangle((lo1, lo2), hi1 - lo1, hi2 - lo2,
                               fill=False, edgecolor=c, linestyle="--",
                               linewidth=1.3))

    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(P1)
    ax.set_ylabel(P2)
    ax.set_title(title or f"Regime geometry  -  projection on ({P1}, {P2})")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.15)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
    return fig
 