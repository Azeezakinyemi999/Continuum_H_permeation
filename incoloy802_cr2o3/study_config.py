"""
================================================================================
INCOLOY 802 (X40 NiCrAlTi 31/19) / Cr2O3 — Schmidt 1985
================================================================================
THE INPUT FILE for the notebooks in this folder. Every material property,
operating condition and sensitivity-study parameter they run on is here.

calculations.config.model_config picks this file up when a notebook runs with this
folder as its working directory — nothing else selects it, and editing it changes
this study alone. To read it from elsewhere:

    MHI_STUDY=incoloy802_cr2o3/study_config.py python your_script.py
"""

import numpy as np


# ==============================================================================
# 1. IDENTITY AND CONSTANTS
#
# Who this study is, and the two physical constants the rest of the file uses.
# ==============================================================================

# Identity of this study: labels figures and CSVs as ACTIVE_STUDY_NAME, and selects
# the reference set in container/model_smoke_test.py. Independent of the folder name,
# so renaming the folder does not silently relabel the results.
STUDY_NAME = 'incoloy802_cr2o3'

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================
R = 8.314      # J/mol/K — Universal gas constant

F = 96485      # C/mol  — Faraday constant (for eV conversions)


# ==============================================================================
# 2. WHAT IS BEING STUDIED
#
# The scenario: the operating point, the state of the wall, how the model is run, and the
# ranges the sensitivity analysis explores around it. These are the entries you edit to
# ask a different question of the same materials.
#
# CONDITIONS is the authority for the operating point. Nothing below may restate it.
# ==============================================================================

# =============================================================================
# OPERATING CONDITIONS
# =============================================================================
CONDITIONS = {
    # ---- the operating point ----
    'T_operating':    873,              # K
    'P_upstream':     1e5,              # Pa
    'P_downstream':   0,                # Pa
    'L_metal':        1e-3,             # m
    'L_oxide':        1e-6,             # m

    # ---- the range of conditions being studied ----
    # These ARE operating conditions, so they are stated here and nowhere else: the
    # sensitivity analysis reads them rather than restating them. It used to keep its
    # own set, and each contained a bound that is nonsense read as a scenario — this
    # dict reached 1e12 Pa (10 million bar) and the SA floor was 1e-10 m, sub-atomic.
    'T_range':        (573, 1273),      # K
    'P_range':        (1e-7, 1e7),      # Pa   1e-7 = UHV, 1e7 = 100 bar
    'L_metal_range':  (1e-4, 5e-3),     # m    0.1 - 5 mm
    # absorbed OXIDES[...]['thickness_range'], which stated the same quantity
    'L_oxide_range':  (1e-9, 1e-5),     # m    1 nm - 10 um

    # ---- sweep resolution: a plotting choice, not a condition ----
    'n_T_points':     20,
    'n_P_points':     40,
    'n_L_points':     20,
}

# -----------------------------------------------------------------------------
# HOW THE MODEL IS RUN — not data, and not a property of anything.
#
# These were scattered: include_gb_trapping inside MICROSTRUCTURE, use_sieverts_pinhole
# inside OXIDE_DEFECTS (with a "# <- ADD THIS" comment still attached), and the last
# three inside DEFAULT_PARAMS_LEVEL5. A reader looking for "what does this study switch
# on" had to know all three places.
#
# The dicts below still carry the flags the solvers read out of them — the consumers take
# a microstructure dict or a defect dict and look inside — but they carry a REFERENCE to
# this one, so there is a single authored value.
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
  

    # Oxide transport (Cr2O3_sample4)
    'D_ox_ref':     [7.800e-21, 7.800e-17],
    'E_D_ox':       [69000, 71000],
    'K_ox_ref':     [0.08, 0.4],
    'H_sol_ox':     [20000, 170000],
   

    # Geometry and operating point — from CONDITIONS, which owns the range of
    # conditions being studied. Material-parameter ranges below stay here:
    # they describe materials, not conditions.
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
    # K_eq ranges are centred on the config values (Grant 1988) at ±2 decades, the
    # same convention as every other *_ref parameter. They previously spanned
    # [1e-14, 1e-6] and [1e-10, 1e-3], written around placeholder defaults of 1e-10
    # and 1e-8: the config value 1e-4 fell 100x ABOVE its own sweep range, and
    # K_eq_metal_ref's 1e-3 sat exactly ON its upper bound, so that parameter was
    # only ever sampled below its nominal value.
    'K_eq_ref':            [1e-6, 1e-2],
    'H_eq':                [5000, 50000],
    # Metal surface kinetics (Grant 1988)
    'k_diss_metal_ref':    [1.346e-8, 1.346e-4],
    'E_diss_metal':        [70000, 100000],
    'K_eq_metal_ref':      [1e-5, 1e-1],
    'H_eq_metal':          [5000, 40000],
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
    """``base`` ranges with ``overrides`` applied — the single preset constructor.

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
# SWEEP POINTS
#
# Draw counts and per-preset yields used to live here as measured constants. They
# are gone: `calculations.sensitivity.size_draws_for_target()` now probes each
# preset and sizes the run from the yield it actually observes. A cached yield
# cannot notice that a range or a preset block was edited, or that the study
# switched material, and silently unbalanced clusters make a cross-regime
# comparison partly a comparison of estimator noise. TARGET_CLUSTER_SIZE lives in
# sensitivity.py beside the code that consumes it.
#
# What remains here is study design, which no probe can infer.
# -----------------------------------------------------------------------------

# Temperatures for the isothermal L5 sweep. With temperature varying it monopolises
# the variance of log10(flux) and leaves every other parameter at or below the noise
# floor; sampling harder does not help, because the floor is flat in n (0.085 at
# n=400, 0.092 at n=1500). The driver ranking genuinely shifts across this span.
DEFAULT_SWEEP_TEMPERATURES = (773.0, 1073.0, 1273.0)

# =============================================================================
# VALIDATION / TEST PARAMETERS
# =============================================================================
VALIDATION = {
    # -------------------------------------------------------------------------
    # OPERATING POINT — shared by every level, read from CONDITIONS.
    #
    # These used to be restated per level under three different names for the same
    # quantity: P_ref (L1, L2a), P_fixed (L2b, L3, L4, L5) and P_sensitivity (L5)
    # were all 1e5, nine literals for one pressure. T_mode_comparison, T_grain_size,
    # T_comparison and T_limit were all 873. A check that runs at a different pressure
    # than the model it is validating is checking something else.
    # -------------------------------------------------------------------------
    'P_up':   CONDITIONS['P_upstream'],
    'P_down': CONDITIONS['P_downstream'],
    'T':      CONDITIONS['T_operating'],
    'T_min':  CONDITIONS['T_range'][0],
    'T_max':  CONDITIONS['T_range'][1],

    # -------------------------------------------------------------------------
    # SWEEP GRIDS — what each level varies, prefixed by the level that owns it.
    #
    # Several are deliberately far wider than the nominal operating window, because
    # that span is the point of the plot; each says so. Point counts stay literal
    # too: resolution is a plotting choice, not an operating condition.
    # -------------------------------------------------------------------------

    # ---- L1: perfect metal ----
    'L1_n_test':         50,
    # 0.01-100 Pa, far below the operating point on purpose: this is where the
    # sqrt(P) Sieverts law is checked against its analytic limit.
    'L1_test_pressures': np.logspace(-2, 2, 50),
    'L1_test_temps_C':   np.array([700, 775, 850, 925]),
    'L1_test_temps_K':   np.array([700, 775, 850, 925]) + 273.15,

    # ---- L2a: perfect oxide ----
    'L2a_n_test':         50,
    # seven decades, to give the oxide's pressure exponent a long lever.
    'L2a_test_pressures': np.logspace(-2, 5, 50),
    'L2a_test_temps_C':   np.array([700, 800, 900, 1000]),
    'L2a_test_temps_K':   np.array([700, 800, 900, 1000]) + 273.15,

    # ---- L2b: oxide + metal in series ----
    'L2b_oxide_thickness_thin': 1e-8,
    # 10 um - 10 cm: four decades to show the 1/L scaling, deliberately wider than
    # CONDITIONS['L_metal_range'].
    'L2b_L_metal_sweep':        np.logspace(-5, -1, 30),
    # down to 1e-12 m: the thin-film limit, where the oxide stops mattering at all.
    # Narrowing this to CONDITIONS['L_oxide_range'] would delete the effect.
    'L2b_delta_ox_sweep':       np.logspace(-4, -12, 40),

    # ---- L3: defective oxide ----
    'L3_defect_fractions_compare': [0.0, 0.001, 0.01, 0.1],
    'L3_f_defect_min':             1e-5,
    'L3_f_defect_max':             0.5,
    'L3_n_defect_points':          40,
    'L3_f_defect_limit_check':     1e-100,
    'L3_defect_types':             ['pinhole', 'crack', 'grain_boundary'],
    'L3_gamma_sweep':              np.logspace(-2, 0, 20),   # crack thickness factors
    'L3_delta_sweep':              np.logspace(0, 4, 20),    # GB diffusivity factors
    'L3_n_enhancement_P_points':   30,
    'L3_enhancement_P_range':      (1e0, 1e10),  # Pa — for the enhancement heatmap

    # ---- L4: defective metal ----
    'L4_n_T_points': 40,                          # finer than L5: a smoother curve
    'L4_grain_size': np.logspace(-8, -3, 30),
    'L4_f_gb_max':   0.5,

    # ---- L5: full system ----
    'L5_n_T_points':    20,
    'L5_f_defect_sweep': np.array([0.001, 0.01, 0.05, 0.1]),
    # 673 K, not T: inside the oxide's validated [473, 773] K window, so the
    # sensitivity check is not reading an extrapolation.
    'L5_T_sensitivity': 673,

    # ---- L6: surface kinetics ----
    # Fast-kinetics limit: each Level 6 model must collapse onto its non-surface
    # counterpart as k_diss → ∞ (L1+L6 → Sieverts, L2a+L6 → L2a, L3+L4+L6 → L5, ...).
    # Parity plots show J_L*+L6(k_diss=k_diss_fast_limit) against J_L*.
    'L6_k_diss_sweep_exp_range': (-15, -3),
    'L6_k_diss_sweep_n_points':  7,
    'L6_k_diss_fast_limit':      1e-3,      # mol/m²/s/Pa — the "fast" regime
}


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
        # The scale thickness of the sample these data were measured on. The model does
        # not read it: the wall being simulated is described by CONDITIONS['L_oxide'],
        # exactly as the metal is. This used to be the value every model call actually
        # used, 20x thinner than the CONDITIONS entry that claimed to set it.
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
        # Validity range of the Nemanic/Stover data; read by
        # oxide_permeation.get_oxide_properties_at_T to warn on extrapolation.
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

# 'N_L' is the same quantity under the symbol the derivations use, and
# Application/Proposal.ipynb reads it 14 times. Aliased rather than written out
# twice: two literals of the same number drift the moment one is edited.
MICROSTRUCTURE['N_L'] = MICROSTRUCTURE['lattice_density']


# ==============================================================================
# 4. DERIVED
#
# Assembled from sections 2 and 3 by reference, never by retyping. This is why the section
# comes last: Python executes top to bottom, so everything these read must already exist.
#
# A literal here is a bug. It is a second copy of a number that lives above, and the two
# will diverge — they already did once, by up to six orders of magnitude.
# ==============================================================================

# The flat parameter sets the model wrappers and the sensitivity analysis consume,
# assembled from sections 2 and 3. Run calculations.sensitivity.check_against_config()
# after editing anything they read.

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
    """First key of ``d``, with a clear error instead of a bare StopIteration."""
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
    'oxide_thickness': CONDITIONS['L_oxide'],   # m

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
    # 'use_sieverts_pinhole' is NOT here: neither model wrapper reads it from the
    # parameter dict. Surface_proposal.ipynb takes it from SIM['oxide_defects'],
    # i.e. straight from OXIDE_DEFECTS, so that is where it lives.

    # Level 4: grain structure  (MICROSTRUCTURE; Zhu 2021 baseline)
    'grain_size':      MICROSTRUCTURE['grain_size'],       # m
    'grain_shape':     MICROSTRUCTURE['grain_shape'],
    'gb_type':         MICROSTRUCTURE['gb_type'],
    'gb_thickness':    MICROSTRUCTURE['gb_thickness'],     # m
    'lattice_density': MICROSTRUCTURE['lattice_density'],  # m⁻³

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
    'include_gb_enhancement': MODEL_OPTIONS['include_gb_enhancement'],
    'include_trapping':       MODEL_OPTIONS['include_trapping'],
    'D_eff_method':           MODEL_OPTIONS['D_eff_method'],
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

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_surface_kinetics_from_config(material_key, temperature_K, material_dict):
    r"""
    Compute temperature-dependent surface kinetics from a METALS or OXIDES entry.

    Uses the 'surface_kinetics' sub-dict within the material entry.
    Arrhenius form:

    .. math::

        k_{diss}(T) &= k_{diss,ref}\exp\!\left[\frac{-E_{diss}}{R}
            \left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right] \\
        K_{eq}(T) &= K_{eq,ref}\exp\!\left[\frac{-H_{eq}}{R}
            \left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right] \\
        k_{recomb} &= \frac{k_{diss}}{K_{eq}}

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
    """Assemble a solver-ready configuration dict for one metal/oxide pair.

    Pulls the chosen entries out of METALS and OXIDES and merges them with the
    operating conditions in CONDITIONS, so a caller hands the model a single dict
    instead of threading a dozen parameters through by hand.

    The oxide defect fractions in OXIDE_DEFECTS are also translated here into the
    per-path ``defect_config`` format the Level 3 solvers expect. A pinhole, crack
    or grain-boundary path is included only when its area fraction is non-zero, so
    setting a fraction to zero removes that transport path entirely rather than
    leaving a zero-weighted one in the parallel sum.

    Parameters
    ----------
    metal, oxide : str
        Keys into METALS and OXIDES. Default to ACTIVE_METAL / ACTIVE_OXIDE, which
        are the first entry of each dict.
    T_operating, P_upstream, L_metal, L_oxide : float, optional
        Override the corresponding CONDITIONS value. L_oxide falls back to the
        oxide entry's own 'thickness' before CONDITIONS. Note these use ``or``, so
        a falsy 0 also takes the default — pass a small positive number instead.
    microstructure_overrides, defect_overrides : dict, optional
        Shallow-merged onto MICROSTRUCTURE and OXIDE_DEFECTS respectively. The
        defect overrides are applied BEFORE ``defect_config`` is derived, so they
        do change which paths exist.

    Returns
    -------
    dict
        Material properties, operating conditions and sweep ranges, plus both the
        raw 'oxide_defects' block and the derived 'defect_config'.
    """
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
        'L_oxide': L_oxide or CONDITIONS['L_oxide'],

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


# ==============================================================================
# 5. PRESENTATION
#
# Figure styling. Affects no number.
# ==============================================================================

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
