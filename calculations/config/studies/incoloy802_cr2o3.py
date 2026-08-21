"""
================================================================================
UNIFIED CONFIGURATION FOR HIERARCHICAL HYDROGEN PERMEATION MODEL
================================================================================
All material properties, operating conditions, and model parameters in one place.
"""

import numpy as np

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================
R = 8.314      # J/mol/K — Universal gas constant
F = 96485      # C/mol  — Faraday constant (for eV conversions)

# =============================================================================
# METAL PROPERTIES (Level 1, 4)
# =============================================================================
METALS = {
        'metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985': {
        'T_ref':   1223,           # K
        'D_ref':   4.285e-09,     # m²/s at T_ref
        'E_D':     72200,         # J/mol
        'D_0':     5.196e-06,      # m²/s
        'K_s_ref': 0.05093,      # mol/m³/Pa^0.5 at T_ref
        'K_s0':    0.1034,      # mol/m³/Pa^0.5
        'H_s':     7200,         # J/mol
        'Phi_ref': 2.182e-10,     # mol/m/s/Pa^0.5 at T_ref
        'Q_p':     79400,         # J/mol
        'Phi_0':   5.373e-07,      # mol/m/s/Pa^0.5

        # ---------------------------------------------------------------
        # SURFACE KINETICS — metal surface (for pinhole paths, L1+L6)
        # ---------------------------------------------------------------
        'surface_kinetics': {
            'T_ref':            965,        # K
            'k_diss_metal_ref': 1.346e-06,      # mol/m²/s/Pa at T_ref
            'E_diss_metal':     81560,      # J/mol
            'K_eq_metal_ref':   1e-03,       # Pa⁻¹ at T_ref PLACEHOLDER
            'H_eq_metal':       15000,      # J/mol PLACEHOLDER
            'reference':        'clean metal surface condition(metal_316_steel_Grant1988): Grant 1988: Grant DM, Cummings DL, Blackburn DA	J. Nuclear Materials 152	1988	DOI: 10.1016/0022-3115(88)90128-7	316 steel H₂ transport. ONLY source with full surface kinetics k1 (3 conditions).	Metals + Surface Kinetics',
        },

        'pressure':       [],
        'reference':      'Schmidt 1985	Schmidt H, Batfalsky P	J. Nuclear Materials 131	1985		H₂ and T₂ permeation in 30+ heat-resistant alloys. Full D+Ks+Phi for many Ni-Cr alloys.	Metals: X40 NiCrAlTi (Incoloy 802) — use as primary reference for Incoloy 800.',
        'temp_range':     [1023, 1223],
        'pressure_range': [1.33e2, 1e5],
        'metal_thickness': [5e-3],
        'gas':            'H₂',
        'Notes':          'Phi_ref = D_ref * Ks_ref.',
    },

    # 'metal_316L_Heat_treated_ref_cast': {
    #     'T_ref':   873,           # K
    #     'D_ref':   2.194e-10,     # m²/s at T_ref
    #     'E_D':     42500,         # J/mol
    #     'D_0':     7.66e-08,      # m²/s
    #     'K_s_ref': 2.723e-2,      # mol/m³/Pa^0.5 at T_ref
    #     'K_s0':    4.640e-1,      # mol/m³/Pa^0.5
    #     'H_s':     20585,         # J/mol
    #     'Phi_ref': 1.881e-11,     # mol/m/s/Pa^0.5 at T_ref
    #     'Q_p':     63085,         # J/mol
    #     'Phi_0':   1.12e-07,      # mol/m/s/Pa^0.5

    #     # ---------------------------------------------------------------
    #     # SURFACE KINETICS — metal surface (for pinhole paths, L1+L6)
    #     # ---------------------------------------------------------------
    #     'surface_kinetics': {
    #         'T_ref':            873,        # K
    #         'k_diss_metal_ref': 1e-12,      # mol/m²/s/Pa at T_ref
    #         'E_diss_metal':     40000,      # J/mol
    #         'K_eq_metal_ref':   1e-8,       # Pa⁻¹ at T_ref
    #         'H_eq_metal':       15000,      # J/mol
    #         'reference':        'Placeholder — replace with literature values',
    #     },

    #     'pressure':       [],
    #     'reference':      'Forcey et al. 1988 — Heat treated reference cast 316L',
    #     'temp_range':     [523, 873],
    #     'pressure_range': [1.33e2, 1e5],
    #     'metal_thickness': [5e-3],
    #     'gas':            'H₂',
    #     'Notes':          'Phi_ref = D_ref * Ks_ref.',
    # },
}


# =============================================================================
# OXIDE PROPERTIES (Level 2a, 2b, 3)
# =============================================================================
OXIDES = {

        'Cr2O3_sample4': {
        # ---- Reference conditions ----
        'T_ref':       673,        # K (400 °C)
        'D_ox_ref':    7.800e-19,  # m²/s 
        'K_ox_ref':    0.35417,      # mol/m³/Pa^0.5 
        # ---- Activation energies ----
        'E_D_ox':      70434,      # J/mol (1.24 eV) 
        'H_sol_ox':    163566,      # J/mol — Q_p − E_D 

        # ---- Pre-exponential factors (Sieverts, Pa^0.5) ----
        'D_ox_0':      2.813e-13,   # m²/s
        'K_ox_0':      1.757e+12,   # mol/m³/Pa^0.5

        # ---- Permeability at T_ref (Sieverts, Pa^0.5) ----
        'Phi_ox_ref':  3.4e-19,  # mol/m/s/Pa^0.5
        'Q_p_ox_J_per_mol': 234000,


        # ---- Default geometry ----
        'thickness':       4.8e-08,
        'thickness_range': [1e-9, 1e-5],

        # ---------------------------------------------------------------
        # SURFACE KINETICS — oxide surface (gas–oxide interface, L6)
        # Used by: get_all_properties → solve_steady_state_flux,
        #          solve_steady_state_flux_L2aL6, calculate_path_flux_L6,
        #          calculate_path_flux_L346_v2
        # Arrhenius form:
        #   k_diss(T) = k_diss_ref × exp(-E_diss/R  × (1/T - 1/T_ref))
        #   K_eq(T)   = K_eq_ref   × exp(-H_eq/R    × (1/T - 1/T_ref))
        #   k_recomb  = k_diss / K_eq  (derived)
        # ---------------------------------------------------------------
        'surface_kinetics': {
            'T_ref':       1623,        # K — same as transport T_ref
            'k_diss_ref':  9.487e-08,       # mol/m²/s/Pa at T_ref
            'E_diss':      57950,       # J/mol
            'K_eq_ref':    1e-4,       # Pa⁻¹ at T_ref PLACEHOLDER. Changed from 1e-10 to 1e-4
            'H_eq':        20000,       # J/mol PLACEHOLDER
            'reference':   'Grant et al. 1988: Grant DM, Cummings DL, Blackburn DA	J. Nuclear Materials 152	1988	DOI: 10.1016/0022-3115(88)90128-7	316 steel H₂ transport. ONLY source with full surface kinetics k1 (3 conditions).	Oxide Surface Kinetics',
        },

        # ---- Metadata ----
        'reference': (
            "Nemanic et al. 2023 (Phi, D), Stover 1986 (Q_p) and Chen 2011 (E_D): Undamaged chromia limit from Stover 1986 Fig.7 (sample Cr2)"
        ),
        'temp_range_K':       [473, 773],
        'temperature_range':  [473, 773],
    },
    

    # 'Al2O3': {
    #     # ---- Reference conditions ----
    #     'T_ref':       1623,        # K (1350°C)
    #     'D_ox_ref':    3.0647e-11,  # m²/s — Belonoshko 2004
    #     'K_ox_ref':    4.9077,      # mol/m³/Pa^0.5 — Sieverts re-assumed at 1 atm

    #     # ---- Activation energies ----
    #     'E_D_ox':      119641,      # J/mol (1.24 eV) — Belonoshko 2004
    #     'H_sol_ox':    198559,      # J/mol — Q_p − E_D = 318200 − 119641

    #     # ---- Pre-exponential factors (Sieverts, Pa^0.5) ----
    #     'D_ox_0':      2.1730e-7,   # m²/s
    #     'K_ox_0':      1.2065e+7,   # mol/m³/Pa^0.5

    #     # ---- Permeability at T_ref (Sieverts, Pa^0.5) ----
    #     'Phi_ox_ref':  1.5040e-10,  # mol/m/s/Pa^0.5
    #     'Phi_ox_0':    2.6216,      # mol/m/s/Pa^0.5

    #     # ---- Raw Roberts values (Pa^0.43) — retained for traceability ----
    #     'Phi_ox_0_raw_Pa043':       5.8744,
    #     'K_ox_ref_raw_Pa043':       1.0997e+1,
    #     'K_ox_0_raw_Pa043':         2.7034e+7,
    #     'Phi_ox_ref_raw_Pa043':     3.3701e-10,
    #     'pressure_exponent_Roberts': 0.43,
    #     'P_reassumption_Pa':        101325,

    #     # ---- Default geometry ----
    #     'thickness':       1e-6,
    #     'thickness_range': [1e-7, 1e-5],

    #     # ---------------------------------------------------------------
    #     # SURFACE KINETICS — oxide surface (gas–oxide interface, L6)
    #     # Used by: get_all_properties → solve_steady_state_flux,
    #     #          solve_steady_state_flux_L2aL6, calculate_path_flux_L6,
    #     #          calculate_path_flux_L346_v2
    #     # Arrhenius form:
    #     #   k_diss(T) = k_diss_ref × exp(-E_diss/R  × (1/T - 1/T_ref))
    #     #   K_eq(T)   = K_eq_ref   × exp(-H_eq/R    × (1/T - 1/T_ref))
    #     #   k_recomb  = k_diss / K_eq  (derived)
    #     # ---------------------------------------------------------------
    #     'surface_kinetics': {
    #         'T_ref':       1623,        # K — same as transport T_ref
    #         'k_diss_ref':  1e-15,       # mol/m²/s/Pa at T_ref
    #         'E_diss':      50000,       # J/mol
    #         'K_eq_ref':    1e-10,       # Pa⁻¹ at T_ref
    #         'H_eq':        20000,       # J/mol (van't Hoff)
    #         'reference':   'Placeholder — replace with literature values',
    #     },

    #     # ---- Metadata ----
    #     'reference': (
    #         'Roberts RM et al. J. Am. Ceram. Soc. 62(9–10), 495–499 (1979). '
    #         'E_D_ox: Belonoshko et al. 2004, Phys. Rev. B 69, 024302. '
    #         'Sieverts re-assumption at P_ref=101325 Pa.'
    #     ),
    #     'temp_range_K':       [1473, 1723],
    #     'temperature_range':  [1473, 1723],
    #     'uncertainty_factor': 10,
    # },
}


# =============================================================================
# MICROSTRUCTURE PARAMETERS (Level 4)
# =============================================================================
# MICROSTRUCTURE = {
#     'grain_size':          50e-6,       # m (50 μm default)
#     'gb_thickness':        0.5e-9,      # m (0.5 nm)
#     'grain_shape':         'equiaxed',  # 'equiaxed', 'columnar'
#     'gb_type':             'LAGB',      # 'HAGB', 'LAGB'
#     'include_gb_trapping': False,

#     'trap_list': [
#         {
#             'name':           'vacancies',
#             'binding_energy': 0.5 * F,  # J/mol (0.5 eV)
#             'density':        1e26,     # m⁻³
#         },
#         {
#             'name':           'dislocations',
#             'binding_energy': 0.7 * F,  # J/mol (0.7 eV)
#             'density':        1e22,     # m⁻³
#         },
#         {
#             'name':           'grain_boundaries',
#             'binding_energy': 0.9 * F,  # J/mol (0.9 eV)
#             'density':        1e24,     # m⁻³
#         },

#         {
#             'name':           'Carbides',
#             'binding_energy': 0 * F,  # J/mol (0 eV)
#             'density':        0,     # m⁻³
#         },
#     ],

MICROSTRUCTURE = {
    'grain_size':          100e-6,      # m  (midpoint 50-150 um, Zhu 2021)
    'gb_thickness':        0.5e-9,      # m  (unchanged -- no Has-N measurement)
    'grain_shape':         'equiaxed',  # confirmed Zhu 2021
    'gb_type':             'LAGB',     # ~50% Sigma3 (0 eV) + ~50% HAGB (0.27 eV)
    'include_gb_trapping': False,       # keep False for baseline; True for sensitivity

    'trap_list': [
        {
            'name':           'vacancies',
            'binding_energy': 0.43 * F,  # J/mol  (Ni est. vs placeholder 0.50)
            'density':        1e26,       # m^-3   (no Has-N data -- keep)
        },
        {
            'name':           'dislocations',
            'binding_energy': 0.20 * F,  # J/mol  (Lu 2022: 0.186-0.215 eV, was 0.70)
            'density':        8.16e12,   # m^-3   (Zhu 2021 GND, was 1e22 -- 4 orders too high)
        },
        {
            'name':           'grain_boundaries',
            'binding_energy': 0.27 * F,  # J/mol  (Lu 2022 Peak 2: 0.258-0.281, was 0.90)
            'density':        6e14,      # m^-3   (geometric from 100 um grain, was 1e24)
        },
        {
            'name':           'Carbides',
            'binding_energy': 0.27 * F,  # J/mol  (M6C Lu 2022: 0.258-0.281, was 0.0)
            'density':        2e25,      # m^-3   (Young 1997 upper bound, was 0)
        },
    ],


    # -------------------------------------------------------------------
    # LATTICE SITE DENSITY
    # -------------------------------------------------------------------
    'lattice_density':     8.774e28,     # m⁻³ — Fe BCC; update for your alloy
    'N_L':                 8.774e28,     # m⁻³ — alias retained for backward compatibility
}


# =============================================================================
# OXIDE DEFECT PARAMETERS (Level 3, 5)
# =============================================================================
OXIDE_DEFECTS = {
    'area_fraction':    0.02,
    'type':             'mixed',
    'components': {
        'pinholes':         0.01,
        'cracks':           0.005,
        'grain_boundaries': 0.005,
    },
    'thickness_factor':   0.1,    # crack: L_crack = 0.1 × L_oxide
    'diffusivity_factor': 10,     # GB: D_gb = 10 × D_oxide
    'use_sieverts_pinhole': False,   # ← ADD THIS
    # False = use metal surface kinetics from METALS['surface_kinetics']
    # True  = assume fast metal kinetics (Sieverts limit, no k_diss_metal needed)
}


# =============================================================================
# OPERATING CONDITIONS
# =============================================================================
CONDITIONS = {
    'T_operating':    873,              # K
    'T_range':        (623, 1200),      # K
    'n_T_points':     20,

    'P_upstream':     1e5,              # Pa
    'P_downstream':   0,                # Pa
    'P_range':        (1e-4, 1e12),      # Pa
    'n_P_points':     40,

    'L_metal':        1e-3,             # m
    'L_oxide':        1e-6,             # m
    'L_metal_range':  (0.1e-3, 5e-3),  # m
    'L_oxide_range':  (1e-7, 1e-5),    # m
    'n_L_points':     20,
}


# # =============================================================================
# # SURFACE KINETICS — flat defaults (Level 6)
# # Used when no per-material Arrhenius data is available.
# # K_eq is derived: K_eq = k_diss / k_recomb
# # =============================================================================
# SURFACE_KINETICS = {
#     'k_diss':         1e-15,                        # mol/m²/s/Pa
#     'k_recomb':       1e-3,                         # m⁴/mol/s
#     'K_eq':           1e-15 / 1e-3,                 # Pa⁻¹ — derived: k_diss / k_recomb
#     'coverage_mode':  'steady_state',               # 'steady_state', 'langmuir'
# }



# =============================================================================
# =============================================================================
# SENSITIVITY ANALYSIS PARAMETERS
#
# Everything the regime-stratified sensitivity analysis lets you tune lives here,
# so there is one file to open when you want to change what the study explores.
# Consumed by calculations/sensitivity.py; nothing here imports it back.
#
# Two kinds of thing, and the distinction matters:
#
#   DEFAULT_PARAMS_*   are DERIVED from METALS / OXIDES / MICROSTRUCTURE /
#                      OXIDE_DEFECTS / CONDITIONS above. Do not hardcode a number
#                      into them — change the source dict instead. They were once
#                      hand-copied literals and had already drifted: the defect
#                      area fractions, the operating temperature and both surface
#                      equilibrium constants disagreed with the values above by up
#                      to six orders of magnitude.
#
#   SUGGESTED_RANGES_*, REGIME_PRESETS_*, draw counts and sweep temperatures are
#                      STUDY DESIGN — deliberately literal. A preset that
#                      suppresses f_pinhole to force an oxide-limited cluster
#                      expresses intent about the analysis, not a material
#                      property, so it is written out explicitly.
#
# Run calculations.sensitivity.check_against_config() after editing.
# =============================================================================
# =============================================================================

# Which METALS / OXIDES entry the sensitivity analysis characterises: THE FIRST
# ENTRY of each dict. Taken positionally rather than spelled out, so adding or
# renaming a material needs no edit here and the name cannot fall out of sync with
# the dict it indexes.
#
# The convention this encodes: the first entry is the active one. Dicts preserve
# insertion order (Python 3.7+), so that means the first literal in the file. If you
# add a second material, put it BELOW the active one — inserting above silently
# switches every derived default.


def _first_key(d, what):
    """First key of `d`, with a clear error instead of a bare StopIteration."""
    try:
        return next(iter(d))
    except StopIteration:
        raise ValueError(f"{what} is empty — cannot determine the active material") from None


ACTIVE_METAL = _first_key(METALS, 'METALS')
ACTIVE_OXIDE = _first_key(OXIDES, 'OXIDES')

_SA_M  = METALS[ACTIVE_METAL]
_SA_O  = OXIDES[ACTIVE_OXIDE]
_SA_MS = _SA_M['surface_kinetics']          # metal surface kinetics
_SA_OS = _SA_O['surface_kinetics']          # oxide surface kinetics
_SA_DC = OXIDE_DEFECTS['components']


def _sa_trap(name):
    """MICROSTRUCTURE trap entry by name (case-insensitive — config uses 'Carbides')."""
    for t in MICROSTRUCTURE['trap_list']:
        if t['name'].lower() == name.lower():
            return t
    raise KeyError(f"no trap named {name!r} in MICROSTRUCTURE['trap_list']; "
                   f"have {[t['name'] for t in MICROSTRUCTURE['trap_list']]}")


# -----------------------------------------------------------------------------
# LEVEL 5 defaults — defective oxide (L3) + defective metal microstructure (L4)
#
# Reference-point Arrhenius: X(T) = X_ref * exp(-E/R * (1/T - 1/T_ref))
#   Incoloy 802 (X40 NiCrAlTi) — T_ref_metal = 1223 K (Schmidt 1985)
#   Cr2O3_sample4              — T_ref_oxide =  673 K (Nemanic 2023 / Stover 1986)
# -----------------------------------------------------------------------------

DEFAULT_PARAMS_LEVEL5 = {
    # Metal transport  (METALS[ACTIVE_METAL])
    'D_ref':       _SA_M['D_ref'],      # m²/s           D_metal at T_ref_metal
    'E_D':         _SA_M['E_D'],        # J/mol          activation energy for diffusion
    'K_s_ref':     _SA_M['K_s_ref'],    # mol/m³/Pa^0.5  K_s at T_ref_metal
    'H_s':         _SA_M['H_s'],        # J/mol          heat of solution
    'T_ref_metal': _SA_M['T_ref'],      # K              measurement reference

    # Oxide transport  (OXIDES[ACTIVE_OXIDE])
    'D_ox_ref':    _SA_O['D_ox_ref'],   # m²/s
    'E_D_ox':      _SA_O['E_D_ox'],     # J/mol
    'K_ox_ref':    _SA_O['K_ox_ref'],   # mol/m³/Pa^0.5
    'H_sol_ox':    _SA_O['H_sol_ox'],   # J/mol
    'T_ref_oxide': _SA_O['T_ref'],      # K

    # Geometry — metal thickness comes from CONDITIONS, NOT METALS['metal_thickness'],
    # which is a LIST of sweep values and would silently poison every model call.
    'metal_thickness': CONDITIONS['L_metal'],   # m
    'oxide_thickness': _SA_O['thickness'],      # m

    # Operating conditions  (CONDITIONS)
    'P_upstream':   CONDITIONS['P_upstream'],   # Pa
    'P_downstream': CONDITIONS['P_downstream'], # Pa
    'temperature':  CONDITIONS['T_operating'],  # K

    # Level 3: oxide defects  (OXIDE_DEFECTS)
    'f_pinhole':              _SA_DC['pinholes'],          # area fraction
    'f_crack':                _SA_DC['cracks'],            # area fraction
    'f_gb_defect':            _SA_DC['grain_boundaries'],  # area fraction
    'crack_thickness_factor': OXIDE_DEFECTS['thickness_factor'],    # L_crack = f × L_oxide
    'gb_diffusivity_factor':  OXIDE_DEFECTS['diffusivity_factor'],  # D_gb = f × D_oxide
    'use_sieverts_pinhole':   OXIDE_DEFECTS['use_sieverts_pinhole'],

    # Level 4: grain structure  (MICROSTRUCTURE; Zhu 2021 baseline)
    'grain_size':      MICROSTRUCTURE['grain_size'],       # m
    'grain_shape':     MICROSTRUCTURE['grain_shape'],
    'gb_type':         MICROSTRUCTURE['gb_type'],
    'gb_thickness':    MICROSTRUCTURE['gb_thickness'],     # m
    'lattice_density': MICROSTRUCTURE['lattice_density'],  # m⁻³

    # No config source — an SA modelling choice, not a material property.
    'gb_enhancement_factor': 100,

    # Level 4: traps  (MICROSTRUCTURE['trap_list'])
    # Lu 2022 binding energies + Young 1997 carbide density
    'trap_dislocation_E_b': _sa_trap('dislocations')['binding_energy'],      # J/mol
    'trap_dislocation_N_T': _sa_trap('dislocations')['density'],             # m⁻³
    'trap_gb_E_b':          _sa_trap('grain_boundaries')['binding_energy'],  # J/mol
    'trap_gb_N_T':          _sa_trap('grain_boundaries')['density'],         # m⁻³
    'trap_vacancy_E_b':     _sa_trap('vacancies')['binding_energy'],         # J/mol
    'trap_vacancy_N_T':     _sa_trap('vacancies')['density'],                # m⁻³
    'trap_carbide_E_b':     _sa_trap('carbides')['binding_energy'],          # J/mol
    'trap_carbide_N_T':     _sa_trap('carbides')['density'],                 # m⁻³

    # Model options — how to run the model, not what it is; no config source
    'include_gb_enhancement': True,
    'include_trapping':       True,
    'D_eff_method':           'average',
}

# -----------------------------------------------------------------------------
# LEVEL 5L6 defaults — adds oxide + metal surface dissociation/recombination.
# Inherits every Level 5 parameter via **DEFAULT_PARAMS_LEVEL5, which is what
# keeps the two levels sharing one set of values for everything they have in common.
# -----------------------------------------------------------------------------

DEFAULT_PARAMS_LEVEL5L6 = {
    **DEFAULT_PARAMS_LEVEL5,
    # Oxide surface kinetics — OXIDES[ACTIVE_OXIDE]['surface_kinetics'] (Grant 1988)
    'k_diss_ref':    _SA_OS['k_diss_ref'],   # mol/m²/s/Pa
    'E_diss':        _SA_OS['E_diss'],       # J/mol
    'K_eq_ref':      _SA_OS['K_eq_ref'],     # Pa⁻¹
    'H_eq':          _SA_OS['H_eq'],         # J/mol
    'T_ref_surface': _SA_OS['T_ref'],        # K

    # Metal surface kinetics — METALS[ACTIVE_METAL]['surface_kinetics'] (Grant 1988)
    'k_diss_metal_ref':    _SA_MS['k_diss_metal_ref'],  # mol/m²/s/Pa
    'E_diss_metal':        _SA_MS['E_diss_metal'],      # J/mol
    'K_eq_metal_ref':      _SA_MS['K_eq_metal_ref'],    # Pa⁻¹
    'H_eq_metal':          _SA_MS['H_eq_metal'],        # J/mol
    'T_ref_surface_metal': _SA_MS['T_ref'],             # K
}


# -----------------------------------------------------------------------------
# PARAMETER GROUPS (for organised reporting)
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# SUGGESTED PARAMETER RANGES — Level 5 (28 parameters)   ** TUNE ME **
#
# Reference-point values span roughly ±2 decades around the reference; activation
# energies span physically meaningful bounds. The four T_ref_* are commented out:
# they are measurement metadata, redundant with the *_ref prefactors, and varying
# them injects an artificial 10^3-10^4 prefactor sweep.
# -----------------------------------------------------------------------------
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

    # Traps (Lu 2022 + Young 1997; narrowed to physically motivated bounds)
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

# -----------------------------------------------------------------------------
# SUGGESTED PARAMETER RANGES — Level 5L6 (36 parameters)   ** TUNE ME **
# Inherits all 28 Level 5 ranges and adds the surface-kinetics parameters.
# -----------------------------------------------------------------------------
SUGGESTED_RANGES_LEVEL5L6 = {
    **SUGGESTED_RANGES_LEVEL5,
    # Oxide surface kinetics (Cr2O3)
    'k_diss_ref':          [9.487e-10, 9.487e-6],
    'E_diss':              [40000, 70000],
    # K_eq ranges are centred on the config values (Grant 1988) at ±2 decades, the
    # same convention as every other *_ref parameter. They previously spanned
    # [1e-14, 1e-6] and [1e-10, 1e-3], written around placeholder defaults of 1e-10
    # and 1e-8: the config value 1e-4 fell 100x ABOVE its own sweep range, and
    # K_eq_metal_ref's 1e-3 sat exactly ON its upper bound, so that parameter was
    # only ever sampled below its nominal value.
    'K_eq_ref':            [1e-6, 1e-2],
    'H_eq':                [5000, 50000],
    #'T_ref_surface':       [600, 1800],
    # Metal surface kinetics (Grant 1988)
    'k_diss_metal_ref':    [1.346e-8, 1.346e-4],
    'E_diss_metal':        [70000, 100000],
    'K_eq_metal_ref':      [1e-5, 1e-1],
    'H_eq_metal':          [5000, 40000],
    #'T_ref_surface_metal': [600, 1100],
}


# -----------------------------------------------------------------------------
# REGIME-TARGETED SAMPLING PRESETS   ** TUNE ME **
#
# Targeting is not a convenience: over the default ranges L5L6 yields essentially
# no surface-limited rows, and L5 is ~98% metal-limited with max frac_defect 0.415
# — a single regime, nothing to stratify. The presets are what create the study.
#
# Each preset = the base ranges with ONLY the regime-controlling parameters
# overridden, left as broad as possible so the per-regime sensitivity reflects
# genuine physics rather than the targeting itself.
# -----------------------------------------------------------------------------

def _sa_ranges_with(base, overrides):
    """`base` ranges with `overrides` applied — the single preset constructor.

    Takes the base explicitly so L5 and L5L6 share one implementation instead of
    two identical helpers differing only in which ranges dict they close over.
    """
    R = {k: list(v) for k, v in base.items()}
    R.update({k: list(v) for k, v in overrides.items()})
    return R


# Preset building blocks — SHARED BETWEEN L5 AND L5L6. Each names one physical
# intent, defined exactly once. Both levels compose their presets from these, so
# the shared parts cannot drift: the slow-oxide block alone was previously written
# out three times (L5L6 oxide, L5 oxide, L5 defect), and the suppressed-defect and
# fast-metal blocks twice each. Editing one copy and missing the others would
# silently make the two levels' oxide clusters non-comparable while every test
# still passed. Surface blocks are L5L6-only — there is no surface step in L5.

# thick, slow oxide -> the oxide is the bottleneck
_SA_SLOW_OXIDE    = {'oxide_thickness': [1e-7, 1e-5],
                     'D_ox_ref':        [7.8e-21, 7.8e-20]}

# thin, fast metal -> the metal is NOT limiting
_SA_FAST_METAL    = {'metal_thickness': [5e-5, 5e-4],
                     'D_ref':           [4.3e-8, 4.3e-7]}

# suppress the defect bypass so flux must cross the intact oxide
_SA_DEFECTS_OFF   = {'f_pinhole':   [1e-7, 1e-5],
                     'f_crack':     [1e-7, 1e-5],
                     'f_gb_defect': [1e-7, 1e-5]}

# large defect area -> the bypass can carry the majority of the flux.
# NOTE total defect area reaches ~0.9 here: a heavily degraded coating, a wider
# claim than the <=0.3 the default ranges imply. State it in any write-up.
_SA_DEFECTS_LARGE = {'f_pinhole':   [0.05, 0.30],
                     'f_crack':     [0.05, 0.30],
                     'f_gb_defect': [0.05, 0.30]}

# L5L6 only — fast surface (high pressure + fast dissociation)
_SA_FAST_SURFACE  = {'P_upstream':  [1e4, 1e7],
                     'k_diss_ref':  [9.5e-8, 9.5e-6]}

# L5L6 only — slow surface (low pressure + slow dissociation on oxide and metal)
_SA_SLOW_SURFACE  = {'P_upstream':       [1e-7, 1e1],
                     'k_diss_ref':       [9.5e-12, 9.5e-9],
                     'k_diss_metal_ref': [1.3e-10, 1.3e-7]}

# Level 5L6 presets. Yields measured on the 36-param production scans:
#   metal    813/1000  = 81.3%   (default ranges)
#   surface  714/2500  = 28.6%
#   oxide   3725/5000  = 74.5%
REGIME_PRESETS = {
    'metal':   _sa_ranges_with(SUGGESTED_RANGES_LEVEL5L6, {}),
    'surface': _sa_ranges_with(SUGGESTED_RANGES_LEVEL5L6, _SA_SLOW_SURFACE),
    'oxide':   _sa_ranges_with(SUGGESTED_RANGES_LEVEL5L6,
                               {**_SA_DEFECTS_OFF, **_SA_SLOW_OXIDE,
                                **_SA_FAST_METAL, **_SA_FAST_SURFACE}),
}

# Level 5 presets (no surface kinetics). The L5 regimes are oxide / metal / defect,
# the third being the PARALLEL bypass rather than a series resistance.
# Yields measured on 300-draw LHS probes (seed 42):
#   metal 97.0%,  oxide 77.7%,  defect 87.3% (max frac_defect 0.999)
REGIME_PRESETS_L5 = {
    'metal':  _sa_ranges_with(SUGGESTED_RANGES_LEVEL5, {}),
    # same oxide-limiting intent as the L5L6 oxide preset, minus the surface block
    'oxide':  _sa_ranges_with(SUGGESTED_RANGES_LEVEL5,
                              {**_SA_DEFECTS_OFF, **_SA_SLOW_OXIDE, **_SA_FAST_METAL}),
    # large defect area AND a slow intact oxide, so the bypass actually wins
    'defect': _sa_ranges_with(SUGGESTED_RANGES_LEVEL5,
                              {**_SA_DEFECTS_LARGE, **_SA_SLOW_OXIDE}),
}


# -----------------------------------------------------------------------------
# SAMPLE SIZES AND SWEEP POINTS   ** TUNE ME **
#
# draws = target cluster size / that preset's measured yield, balanced at ~1500 per
# cluster. Balance matters as much as absolute size: delta/PAWN estimator variance
# depends on n, so unequal clusters make the cross-regime comparison compare noise
# levels as well as sensitivities. Column-normalisation in the heatmap rescales
# magnitude but cannot undo that. 1500 also leaves headroom over the ~300-500
# stability floor.
# -----------------------------------------------------------------------------
DEFAULT_N_PER_REGIME    = {'metal': 1850, 'surface': 5250, 'oxide': 2010}
DEFAULT_N_PER_REGIME_L5 = {'metal': 1550, 'oxide': 1930, 'defect': 1720}

# Temperatures for the isothermal L5 sweep. With temperature varying it monopolises
# the variance of log10(flux) and leaves every other parameter at or below the noise
# floor; sampling harder does not help, because the floor is flat in n (0.085 at
# n=400, 0.092 at n=1500). The driver ranking genuinely shifts across this span.
DEFAULT_SWEEP_TEMPERATURES = (773.0, 1073.0, 1273.0)

# In-regime yield of each L5 preset AT FIXED TEMPERATURE (250-draw probes, seed 42).
# These differ from the varying-T yields, so draw counts must be sized per
# temperature — a single dict would leave the oxide cluster ~40% short at 1273 K and
# reintroduce the cross-regime imbalance the isothermal design exists to remove.
#
# The oxide yield falls steadily with temperature (0.92 -> 0.59) because the oxide's
# activation energies (E_D_ox ~70 kJ/mol, H_sol_ox up to ~164 kJ/mol) exceed the
# metal's, so it speeds up faster than the metal and stops being the bottleneck.
MEASURED_YIELDS_L5_BY_T = {
    773.0:  {'oxide': 0.916, 'metal': 0.992, 'defect': 0.916},
    1073.0: {'oxide': 0.700, 'metal': 0.988, 'defect': 0.820},
    1273.0: {'oxide': 0.592, 'metal': 0.972, 'defect': 0.800},
}

# Same, for the L5L6 presets (250-draw probes, seed 42). L5L6 is run isothermally for
# the same reason as L5, and the case is if anything stronger: a dummy-parameter test on
# varying-T L5L6 clusters left temperature at 3.4-4.3x the noise floor and, in the metal
# cluster, NOTHING else above 1.0x — the runner-up scored 0.84x, i.e. below a random
# input. L5L6 spreads 36 parameters over the variance left after temperature takes its
# share, versus L5's 28, so it has less signal per parameter, not more.
#
# The surface yield RISES steeply with temperature (0.23 -> 0.57) — the opposite trend to
# L5's oxide yield — so per-temperature sizing matters even more here. At 773 K the
# surface preset needs ~6.5k draws for a 1500-point cluster against ~2.6k at 1273 K.
#
# NOTE for 'theta': temperature does NOT dominate that metric (1.54x oxide, 0.88x metal),
# because theta is bounded [0,1] and analysed linearly, so it cannot absorb an
# exponential response. L5L6 theta results from a varying-T run are far less compromised
# than the flux ones.
MEASURED_YIELDS_L5L6_BY_T = {
    773.0:  {'surface': 0.232, 'oxide': 0.780, 'metal': 0.772},
    1073.0: {'surface': 0.476, 'oxide': 0.836, 'metal': 0.852},
    1273.0: {'surface': 0.568, 'oxide': 0.820, 'metal': 0.864},
}


# =============================================================================
# PLOTTING STYLE
# =============================================================================
PLOT_STYLE = {
    'figsize':              (14, 12),
    'fontsize_title':       10,
    'fontsize_suptitle':    18,
    'fontsize_axis':        14,
    'fontsize_tick':        12,
    'fontsize_legend':      12,
    'fontsize_annotation':  11,
    'linewidth':            2.5,
    'markersize':           8,
    'grid_alpha':           0.3,
}

COLORS = {
    'L1':       'black',
    'L2a':      'blue',
    'L2b':      'purple',
    'L3':       'cyan',
    'L4_gb':    'green',
    'L4_trap':  'orange',
    'L4_both':  'crimson',
    'L5':       'red',
}

# =============================================================================
# CURVE STYLE SYSTEM
# Each level/mode gets a unique combination of color + linestyle + marker
# Distinguishable by color alone, linestyle alone, or marker alone
# =============================================================================
CURVE_STYLES = {
    # Model levels
    'L1':       {'color': 'black',   'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L1_L6':    {'color': 'black',   'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L2a':      {'color': '#1f77b4', 'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},  # blue
    'L2a_L6':   {'color': '#1f77b4', 'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L2b':      {'color': 'purple',  'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L2_L6':    {'color': 'purple',  'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L3':       {'color': 'cyan',    'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L3_L6':    {'color': 'teal',    'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},
    'L4_gb':    {'color': 'green',   'ls': '-',    'marker': '^',  'lw': 2.5, 'ms': 6},
    'L4_trap':  {'color': 'orange',  'ls': '-',    'marker': 'v',  'lw': 2.5, 'ms': 6},
    'L4_both':  {'color': 'crimson', 'ls': '-',    'marker': 'D',  'lw': 2.5, 'ms': 6},
    'L5':       {'color': 'red',     'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 6},
    'L3_L4_L6': {'color': 'red',     'ls': '--',   'marker': 's',  'lw': 2.5, 'ms': 6},

    # Defect fractions (used in L3 sweep plots)
    'f_0':      {'color': 'black',   'ls': '-',    'marker': None, 'lw': 2.5, 'ms': 0},
    'f_001':    {'color': '#2196F3', 'ls': '--',   'marker': None, 'lw': 2.0, 'ms': 0},  # 0.1%
    'f_01':     {'color': '#FF9800', 'ls': '-.',   'marker': None, 'lw': 2.0, 'ms': 0},  # 1%
    'f_10':     {'color': '#F44336', 'ls': ':',    'marker': None, 'lw': 2.0, 'ms': 0},  # 10%

    # Defect types
    'intact':   {'color': 'black',   'ls': '-',    'marker': 'o',  'lw': 2.5, 'ms': 5},
    'pinhole':  {'color': '#F44336', 'ls': '--',   'marker': '^',  'lw': 2.0, 'ms': 5},  # red
    'crack':    {'color': '#FF9800', 'ls': '-.',   'marker': 's',  'lw': 2.0, 'ms': 5},  # orange
    'gb':       {'color': '#4CAF50', 'ls': ':',    'marker': 'D',  'lw': 2.0, 'ms': 5},  # green

    # Rate-limiting regions (shading)
    'surface_region': {'color': '#F44336', 'alpha': 0.9},  # red
    'oxide_region':   {'color': '#FF9800', 'alpha': 0.9},  # orange
    'metal_region':   {'color': '#2196F3', 'alpha': 0.9},  # blue
    'mixed_region':   {'color': '#4CAF50', 'alpha': 0.9},  # green

    # Reference lines (slopes, theory)
    'slope_1':  {'color': '#F44336', 'ls': '--',   'lw': 1.5, 'alpha': 0.6},
    'slope_05': {'color': '#4CAF50', 'ls': '--',   'lw': 1.5, 'alpha': 0.6},
    'parity':   {'color': 'red',     'ls': '--',   'lw': 2.0, 'alpha': 0.8},
    'threshold':{'color': 'gray',    'ls': '--',   'lw': 1.5, 'alpha': 0.5},
}


# =============================================================================
# PLOT HELPER — apply a curve style in one call
# =============================================================================
def apply_style(ax, x, y, style_key, label=None, log=False):
    """
    Plot a curve using a named style from CURVE_STYLES.

    Parameters
    ----------
    ax        : matplotlib Axes
    x, y      : array-like
    style_key : str — key in CURVE_STYLES
    label     : str — legend label (optional)
    log       : bool — if True use loglog, else plot

    Returns
    -------
    matplotlib Line2D
    """
    s = CURVE_STYLES[style_key]
    kwargs = dict(
        color=s['color'],
        linestyle=s['ls'],
        linewidth=s['lw'],
        alpha=s.get('alpha', 1.0),
    )
    if s.get('marker') and s.get('ms', 0) > 0:
        kwargs['marker']     = s['marker']
        kwargs['markersize'] = s['ms']
        kwargs['markevery']  = 5   # don't crowd markers
    if label:
        kwargs['label'] = label

    if log:
        return ax.loglog(x, y, **kwargs)[0]
    else:
        return ax.plot(x, y, **kwargs)[0]

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_surface_kinetics_from_config(material_key, temperature_K, material_dict):
    """
    Compute temperature-dependent surface kinetics from a METALS or OXIDES entry.

    Uses the 'surface_kinetics' sub-dict within the material entry.
    Arrhenius form:
        k_diss(T) = k_diss_ref × exp(-E_diss/R × (1/T - 1/T_ref))
        K_eq(T)   = K_eq_ref   × exp(-H_eq/R   × (1/T - 1/T_ref))
        k_recomb  = k_diss / K_eq

    Parameters
    ----------
    material_key  : str   — key in METALS or OXIDES
    temperature_K : float — temperature [K]
    material_dict : dict  — either METALS or OXIDES

    Returns
    -------
    dict : k_diss, k_recomb, K_eq
    """
    entry = material_dict[material_key]
    sk    = entry['surface_kinetics']
    T_ref = sk['T_ref']

    # Oxide surface kinetics keys
    if 'k_diss_ref' in sk:
        k_diss_ref = sk['k_diss_ref']
        E_diss     = sk['E_diss']
        K_eq_ref   = sk['K_eq_ref']
        H_eq       = sk['H_eq']
    # Metal surface kinetics keys
    else:
        k_diss_ref = sk['k_diss_metal_ref']
        E_diss     = sk['E_diss_metal']
        K_eq_ref   = sk['K_eq_metal_ref']
        H_eq       = sk['H_eq_metal']

    inv_T_diff = 1.0 / temperature_K - 1.0 / T_ref
    k_diss  = k_diss_ref * np.exp((-E_diss / R) * inv_T_diff)
    K_eq    = K_eq_ref   * np.exp((-H_eq   / R) * inv_T_diff)
    k_recomb = k_diss / K_eq

    return {'k_diss': k_diss, 'k_recomb': k_recomb, 'K_eq': K_eq}


def build_simulation_config(
    metal=ACTIVE_METAL,   # first entry of METALS — see ACTIVE_METAL above
    oxide=ACTIVE_OXIDE,   # first entry of OXIDES
    T_operating=None,
    P_upstream=None,
    L_metal=None,
    L_oxide=None,
    microstructure_overrides=None,
    defect_overrides=None,
):
    od = OXIDE_DEFECTS.copy()
    if defect_overrides:
        od.update(defect_overrides)

    # Build solver-ready defect config from OXIDE_DEFECTS structure
    defect_config = {}
    f_pin = od['components'].get('pinholes', 0.0)
    f_cra = od['components'].get('cracks', 0.0)
    f_gb  = od['components'].get('grain_boundaries', 0.0)
    if f_pin > 0:
        defect_config['pinhole'] = {'area_fraction': f_pin}
    if f_cra > 0:
        defect_config['crack'] = {
            'area_fraction':    f_cra,
            'thickness_factor': od.get('thickness_factor', 0.1),
        }
    if f_gb > 0:
        defect_config['grain_boundary'] = {
            'area_fraction':      f_gb,
            'diffusivity_factor': od.get('diffusivity_factor', 100.0),
        }

    config = {
        'metal_name':  metal,
        'oxide_name':  oxide,
        'metal_props': METALS[metal].copy(),
        'oxide_props': OXIDES[oxide].copy(),

        'T_operating':  T_operating or CONDITIONS['T_operating'],
        'P_upstream':   P_upstream  or CONDITIONS['P_upstream'],
        'P_downstream': CONDITIONS['P_downstream'],

        'L_metal': L_metal or CONDITIONS['L_metal'],
        'L_oxide': L_oxide or OXIDES[oxide].get('thickness', CONDITIONS['L_oxide']),

        'microstructure': MICROSTRUCTURE.copy(),
        'oxide_defects':  od,              # raw OXIDE_DEFECTS (with any overrides)
        'defect_config':  defect_config,   # ← solver-ready format, built automatically

        'T_range':       CONDITIONS['T_range'],
        'P_range':       CONDITIONS['P_range'],
        'L_metal_range': CONDITIONS['L_metal_range'],
        'L_oxide_range': CONDITIONS['L_oxide_range'],
        'n_T_points':    CONDITIONS['n_T_points'],
        'n_P_points':    CONDITIONS['n_P_points'],
        'n_L_points':    CONDITIONS['n_L_points'],
    }

    if microstructure_overrides:
        config['microstructure'].update(microstructure_overrides)

    return config


# =============================================================================
# VALIDATION / TEST PARAMETERS
# =============================================================================
VALIDATION = {
    'L1': {
        'P_ref':         1e5,
        'n_test':        50,
        'test_pressures': np.logspace(-2, 2, 50),
        'test_temps_C':  np.array([700, 775, 850, 925]),
        'test_temps_K':  np.array([700, 775, 850, 925]) + 273.15,
    },

    'L2a': {
        'P_ref':         1e5,
        'n_test':        50,
        'test_pressures': np.logspace(-2, 5, 50),
        'test_temps_C':  np.array([700, 800, 900, 1000]),
        'test_temps_K':  np.array([700, 800, 900, 1000]) + 273.15,
    },

    'L2b': {
        'P_fixed':              1e5,
        'P_fixed_metal_sweep':  1e5,
        'oxide_thickness_thin': 1e-8,
        'L_metal_sweep':        np.logspace(-5, -1, 30),
        'P_fixed_oxide_sweep':  1e5,
        'delta_ox_sweep':       np.logspace(-4, -12, 40),
    },

    'L3': {
        'P_fixed':                1e5,
        'P_down':                 0,
        'defect_fractions_compare': [0.0, 0.001, 0.01, 0.1],
        'f_defect_min':           1e-5,
        'f_defect_max':           0.5,
        'n_defect_points':        40,
        'f_defect_limit_check':   1e-100,
        # New for L3+L6
        'defect_types':             ['pinhole', 'crack', 'grain_boundary'],
        'gamma_sweep':              np.logspace(-2, 0, 20),   # crack thickness factors
        'delta_sweep':              np.logspace(0, 4, 20),    # GB diffusivity factors
        'n_enhancement_P_points':   30,        # pressure points for enhancement map
        'enhancement_P_range':      (1e0, 1e10),  # Pa — for enhancement heatmap
    },

    'L4': {
        'P_fixed':            1e5,
        'P_down':             0,
        'T_mode_comparison':  873,
        'T_min':              473,
        'T_max':              1173,
        'n_T_points':         40,
        'T_grain_size':       873,
        'grain_size':         np.logspace(-8, -3, 30),
        'f_gb_max':           0.5,
    },

    'L5': {
        'P_down':         0,
        'T_comparison':   873,
        'T_min':          573,
        'T_max':          1173,
        'n_T_points':     20,
        'P_fixed':        1e5,
        'T_sensitivity':  673,
        'P_sensitivity':  1e5,
        'f_defect_sweep': np.array([0.001, 0.01, 0.05, 0.1]),
        'T_limit':        873,
    },

    'L6': {
        # ========================================================================
        # SURFACE KINETICS (L6) — Limit Checks and Model Recovery
        # ========================================================================
        # Fast-kinetics limit: verify recovery of limiting cases when k_diss → ∞
        #   - L1+L6 → Sieverts (pure bulk diffusion, no surface limit)
        #   - L2a+L6 → L2a (pure oxide diffusion, no surface limit)
        #   - L2+L6 → series resistance (oxide + metal, no surface limit)
        #   - L3+L6 → L3 (defective oxide, no surface limit)
        #   - L4+L6 → L4 (defective metal, no surface limit)
        #   - L3+L4+L6 → L5 (full defective system, no surface limit)
        # Parity plots show J_L*+L6(k_diss=1e-3) vs J_L* or J_analytical
        # ========================================================================
        'k_diss_sweep_exp_range': (-15, -3),    # exponent range for log sweep
        'k_diss_sweep_n_points': 7,              # number of points in sweep
        # Fast kinetics limit value (used in parity plots)
        'k_diss_fast_limit': 1e-3,               # mol/m²/s/Pa — "fast" regime
    },
}
