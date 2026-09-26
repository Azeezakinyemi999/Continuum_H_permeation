"""
================================================================================
316L STAINLESS / Cr2O3 — Guo et al. 2025
================================================================================
THE INPUT FILE for the notebooks in this folder. Every material property,
operating condition and sensitivity-study parameter they run on is here.

calculations.config.model_config picks this file up when a notebook runs with this
folder as its working directory — nothing else selects it, and editing it changes
this study alone. To read it from elsewhere:

    PERMEATION_STUDY=Guo_etal_2025/study_config.py python your_script.py
"""



# ==============================================================================
# 1. IDENTITY AND CONSTANTS
#
# Who this study is, and the two physical constants the rest of the file uses.
# ==============================================================================


STUDY_NAME = 'Guo_etal_2025_316L'

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================
R = 8.314      # J/mol/K — Universal gas constant

F = 96485      # C/mol  — Faraday constant (for eV conversions)


# ==============================================================================
# 2. WHAT IS BEING STUDIED
#   OPERATING CONDITIONS
# ==============================================================================

CONDITIONS = {
    # ---- the operating point ----
    'T_operating':    873,              # K
    'P_upstream':     1e5,              # Pa
    'P_downstream':   0,                # Pa
    'L_metal':        1e-3,             # m
    'L_oxide':        1e-6,             # m

    # ---- the range of conditions being studied ----
    'T_range':        (573, 1273),      # K
    'P_range':        (1e-7, 1e7),      # Pa   1e-7 = UHV, 1e7 = 100 bar
    'L_metal_range':  (1e-4, 5e-3),     # m    0.1 - 5 mm
    'L_oxide_range':  (1e-9, 1e-5),     # m    1 nm - 10 um

    # ---- what a defect study varies ----
    # Defect area fraction: the fraction of the oxide area carrying a bypass path.
    'defect_fractions_compare': [0.0, 0.001, 0.01, 0.1],   # the curves a comparison draws
    'defect_fraction_range':    (1e-5, 0.5),               # swept axis, dimensionless
    'defect_limit_check':       1e-100,                    # f -> 0: L3 must collapse onto L2b
    'f_gb_max':                 0.5,                       # cap on GB volume fraction 3*d_gb/d

    # ---- sweep resolution: a plotting choice, not a condition ----
    'n_points':       40,
}



# -----------------------------------------------------------------------------
# HOW THE MODEL IS RUN — not data, and not a property of anything.
# -----------------------------------------------------------------------------
MODEL_OPTIONS = {
    # Level 4: trap grain boundaries as well as enhancing diffusion along them.
    'include_gb_trapping':    False,
    # Level 3: treat a pinhole as fast metal kinetics (Sieverts) rather than using the
    # metal's own surface kinetics.
    'use_sieverts_pinhole':   False,
    # Level 4, read by the model wrappers.
    'include_gb_enhancement': True,
    'include_trapping':       True,
    'D_eff_method':           'average',
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
    'use_sieverts_pinhole': MODEL_OPTIONS['use_sieverts_pinhole'],
    # False = use metal surface kinetics from METALS['surface_kinetics']
    # True  = assume fast metal kinetics (Sieverts limit, no k_diss_metal needed)
}

# -----------------------------------------------------------------------------
# SUGGESTED PARAMETER RANGES — Level 5 (28 parameters)   ** TUNE ME **
# -----------------------------------------------------------------------------
SUGGESTED_RANGES_LEVEL5 = {
    # Metal transport (Incoloy 802, Schmidt 1985)
    'D_ref':        [4.285e-11, 4.285e-7],
    'E_D':          [50000, 80000],
    'K_s_ref':      [1e-4, 0.1],
    'H_s':          [1000, 50000],
  

    # Oxide transport (Cr2O3_sample4)
    'D_ox_ref':     [7.800e-21, 7.800e-17],
    'E_D_ox':       [69000, 71000],
    'K_ox_ref':     [0.08, 0.4],
    'H_sol_ox':     [20000, 170000],

    # Geometry and operating point — from CONDITIONS
    'metal_thickness': list(CONDITIONS['L_metal_range']),
    'oxide_thickness': list(CONDITIONS['L_oxide_range']),
    'P_upstream':      list(CONDITIONS['P_range']),

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
    'temperature': list(CONDITIONS['T_range']),
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
    'K_eq_ref':            [1e-6, 1e-2],
    'H_eq':                [5000, 50000],
    # Metal surface kinetics (Grant 1988)
    'k_diss_metal_ref':    [1.346e-8, 1.346e-4],
    'E_diss_metal':        [70000, 100000],
    'K_eq_metal_ref':      [1e-5, 1e-1],
    'H_eq_metal':          [5000, 40000],
}
# The floor is the sensitivity range's floor times 1e-5.
CONDITIONS['k_diss_sweep_range'] = (SUGGESTED_RANGES_LEVEL5L6['k_diss_ref'][0] * 1e-5,
                                    SUGGESTED_RANGES_LEVEL5L6['k_diss_ref'][1])


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
    """``base`` ranges with ``overrides`` applied — the single preset constructor.

    Takes the base explicitly so L5 and L5L6 share one implementation instead of
    two identical helpers differing only in which ranges dict they close over.
    """
    R = {k: list(v) for k, v in base.items()}
    R.update({k: list(v) for k, v in overrides.items()})
    return R

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

REGIME_PRESETS = {
    'metal':   _sa_ranges_with(SUGGESTED_RANGES_LEVEL5L6, {}),
    'surface': _sa_ranges_with(SUGGESTED_RANGES_LEVEL5L6, _SA_SLOW_SURFACE),
    'oxide':   _sa_ranges_with(SUGGESTED_RANGES_LEVEL5L6,
                               {**_SA_DEFECTS_OFF, **_SA_SLOW_OXIDE,
                                **_SA_FAST_METAL, **_SA_FAST_SURFACE}),
}

# Level 5 presets (no surface kinetics). The L5 regimes are oxide / metal / defect
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
# SWEEP POINTS
# What remains here is study design, which no probe can infer.
# -----------------------------------------------------------------------------
# Temperatures for the isothermal L5 sweep. 
DEFAULT_SWEEP_TEMPERATURES = (773.0, 1073.0, 1273.0)

# ==============================================================================
# 3. WHAT THE MATERIALS ARE
#
# Measured properties, with the provenance of each measurement. Nothing here describes the
# wall being simulated — no thickness, no operating point. A material entry states what the
# substance does at a temperature; CONDITIONS above states what is being built from it.
# ==============================================================================

# =============================================================================
# METAL PROPERTIES (Level 1, 4)
# =============================================================================
METALS = {
        '316L_Guo_2025': {
        'T_ref':   673,           # K
        'D_ref':   2.8e-11,     # m²/s at T_ref
        'E_D':     52102,         # J/mol
        'D_0':     3.1e-07,      # m²/s
        'K_s_ref': 0.039,      # mol/m³/Pa^0.5 at T_ref
        'K_s0':    0.219,      # mol/m³/Pa^0.5
        'H_s':     9648,         # J/mol
        'Phi_ref': 1.0e-12,     # mol/m/s/Pa^0.5 at T_ref
        'Q_p':     61750,         # J/mol
        'Phi_0':   6.79e-08,      # mol/m/s/Pa^0.5
        # Guo et al. 2025
        # ---------------------------------------------------------------
        # SURFACE KINETICS — metal surface (for pinhole paths, L1+L6)
        # ---------------------------------------------------------------
        'surface_kinetics': {
            'T_ref':            673,        # K
            'k_diss_metal_ref': 2.6287e-11,      # mol/m²/s/Pa at T_ref
            'E_diss_metal':     54996,      # J/mol
            'K_eq_metal_ref':   1.52e-03,       # Pa⁻¹ at T_ref 
            'H_eq_metal':       19297,      # J/mol 
            'reference':        'Guo et al. 2025',
        },

        'pressure':       [],
        'reference':      'Guo et al. 2025',
        'temperature_range': [468, 823],
        'pressure_range': [1.33e2, 1e5],
        'sample_thickness': [5e-3],
        'gas':            'H₂',
        'Notes':          'Phi_ref = D_ref * Ks_ref.',
    }
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


        # ---- Provenance, NOT geometry ----
        'sample_thickness': 4.8e-08,

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
        # Validity range of the Nemanic/Stover data
        'temperature_range':  [473, 773],
    },
    
}

MICROSTRUCTURE = {
    'grain_size':          100e-6,      # m  (midpoint 50-150 um, Zhu 2021)
    'gb_thickness':        0.5e-9,      # m  (unchanged -- no Has-N measurement)
    'grain_shape':         'equiaxed',  # confirmed Zhu 2021
    'gb_type':             'LAGB',     # ~50% Sigma3 (0 eV) + ~50% HAGB (0.27 eV)
    'include_gb_trapping': MODEL_OPTIONS['include_gb_trapping'],

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
}


MICROSTRUCTURE['N_L'] = MICROSTRUCTURE['lattice_density']
