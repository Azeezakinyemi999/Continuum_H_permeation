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
    metal='metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985',
    oxide='Cr2O3_sample4',
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
