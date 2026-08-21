# """
# ================================================================================
# UNIFIED CONFIGURATION FOR HIERARCHICAL HYDROGEN PERMEATION MODEL
# ================================================================================
# All material properties, operating conditions, and model parameters in one place.
# """

# import numpy as np

# # =============================================================================
# # PHYSICAL CONSTANTS
# # =============================================================================
# R = 8.314          # J/mol/K - Universal gas constant
# F = 96485          # C/mol - Faraday constant (for eV conversions)

# # =============================================================================
# # METAL PROPERTIES (Level 1, 4)
# # =============================================================================
# METALS = {
#     'Fe_alpha': {
#         # Reference conditions
#         'T_ref': 625,                # K
#         'D_ref': 1e-11,              # m²/s at T_ref
#         'K_s_ref': 181.615,          # mol/m³/Pa^0.5 at T_ref
        
#         # Activation energies
#         'E_D': 54000,                # J/mol (diffusion)
#         'H_s': 5900,                 # J/mol (solution enthalpy)
        
#         # Metadata
#         'reference': 'San Marchi 2007',
#         'temp_range_C': [149.87, 626.87],
#     },
    
#     'Incoloy800': {
#         'T_ref': 1073.15,            # K (800°C)
#         'D_ref': 1.43e-9,            # m²/s
#         'K_s_ref': 5.92e-2,          # mol/m³/Pa^0.5
#         'E_D': 52000,                # J/mol
#         'H_s': -20000,               # J/mol
#         'reference': 'JAERI-Tech 2002-090, Forcey et al. (1988)',
#         'temp_range_C': [600, 1000],
#     },
# }

# # =============================================================================
# # OXIDE PROPERTIES (Level 2a, 2b, 3)
# # =============================================================================
# OXIDES = {
#     'Cr2O3': {
#         # Reference conditions
#         'T_ref': 1073.15,            # K (800°C)
#         'D_ox_ref': 1e-9,            # m²/s
#         'K_ox_ref': 1e-6,            # mol/m³/Pa (Henry's law)
        
#         # Activation energies
#         'E_D_ox': 1.55e5,            # J/mol
#         'H_sol_ox': 1.85e5,          # J/mol
        
#         # Default geometry
#         'thickness': 1e-6,   # m (1 μm)
#         'thickness_range': [1e-7, 1e-5],
        
#         # Metadata
#         'reference': 'Strehlow & Savage (1974), Serra (1998)',
#         'temp_range_K': [873, 1273],
#         'temperature_range': [873, 1273],
#         'uncertainty_factor': 10,
#     },
    
#     'Cr2O3_thin': {
#         'T_ref': 1073.15,
#         'D_ox_ref': 1e-18,
#         'K_ox_ref': 1e-12,
#         'E_D_ox': 50000,
#         'H_sol_ox': 30000,
#         'thickness': 6e-10,  # 6 Å
#         'reference': 'Zarchy & Axtmann (1979)',
#     },
# }

# # =============================================================================
# # MICROSTRUCTURE PARAMETERS (Level 4)
# # =============================================================================
# MICROSTRUCTURE = {
#     # Grain boundary parameters
#     'grain_size': 50e-6,             # m (50 μm default, use 100e-9 for nanocrystalline)
#     'gb_thickness': 0.5e-9,          # m (0.5 nm)
#     'grain_shape': 'equiaxed',       # 'equiaxed', 'columnar'
#     'gb_type': 'LAGB',               # 'HAGB', 'LAGB'
#     'include_gb_trapping': False,
    
#     # Trap parameters
#     'trap_list': [
#         {
#             'name': 'vacancies',
#             'binding_energy': 0.5 * F,   # J/mol (0.5 eV)
#             'density': 1e26,              # m⁻³
#         },
#     ],
    
#     # Lattice site density (for Oriani model)
#     'N_L': 8.46e28,                  # m⁻³ (Fe BCC lattice sites)
# }

# # =============================================================================
# # OXIDE DEFECT PARAMETERS (Level 3, 5)
# # =============================================================================
# OXIDE_DEFECTS = {
#     'area_fraction': 0.01,           # 1% pinholes
#     'type': 'pinhole',               # 'pinhole', 'crack'
# }

# # =============================================================================
# # OPERATING CONDITIONS
# # =============================================================================
# CONDITIONS = {
#     # Temperature
#     'T_operating': 625,                    # K - Reference temperature
#     'T_range': (423, 700),           # K - Sweep range
#     'n_T_points': 20,
    
#     # Pressure
#     'P_upstream': 1e4,               # Pa (10 kPa)
#     'P_downstream': 0,               # Pa
#     'P_range': (1e-7, 1e26),           # Pa - Sweep range
#     'n_P_points': 30,
    
#     # Geometry
#     'L_metal': 1e-3,                 # m (1 mm)
#     'L_oxide': 1e-6,                 # m (1 μm)
#     'L_metal_range': (0.1e-3, 5e-3), # m
#     'L_oxide_range': (1e-7, 1e-5),   # m
#     'n_L_points': 20,
# }

# # =============================================================================
# # SURFACE KINETICS PARAMETERS (Level 6 - if used)
# # =============================================================================
# SURFACE_KINETICS = {
#     'k_diss': 1e-15,                 # mol/m²/s/Pa
#     'k_recomb': 1e-3,                # m⁴/mol/s
#     'coverage_mode': 'steady_state', # 'steady_state', 'langmuir'
# }


# # =============================================================================
# # PLOTTING STYLE
# # =============================================================================
# PLOT_STYLE = {
#     'figsize': (14, 12),
#     'fontsize_title': 16,
#     'fontsize_suptitle': 18,
#     'fontsize_axis': 14,
#     'fontsize_tick': 12,
#     'fontsize_legend': 12,
#     'fontsize_annotation': 11,
#     'linewidth': 2.5,
#     'markersize': 8,
#     'grid_alpha': 0.3,
# }

# # Color scheme by level
# COLORS = {
#     'L1': 'black',           # Perfect Metal (baseline)
#     'L2a': 'blue',           # Perfect Oxide Only
#     'L2b': 'purple',         # Oxide + Metal
#     'L3': 'cyan',            # Defective Oxide + Metal
#     'L4_gb': 'green',        # Defective Metal - GB mode
#     'L4_trap': 'orange',     # Defective Metal - Trapping mode
#     'L4_both': 'crimson',    # Defective Metal - Combined
#     'L5': 'red',             # Full System
# }


# # =============================================================================
# # HELPER FUNCTION: Build complete config dict for a simulation
# # =============================================================================
# def build_simulation_config(
#     metal='Fe_alpha',
#     oxide='Cr2O3',
#     T_operating=None,
#     P_upstream=None,
#     L_metal=None,
#     L_oxide=None,
#     microstructure_overrides=None,
#     defect_overrides=None,
# ):
#     """
#     Build a complete configuration dictionary for a simulation.
    
#     Parameters
#     ----------
#     metal : str
#         Key from METALS dict
#     oxide : str
#         Key from OXIDES dict
#     T_operating : float, optional
#         Operating temperature [K] (NOT the reference temperature for properties)
#     P_upstream : float, optional
#         Override upstream pressure (Pa)
#     L_metal : float, optional
#         Override metal thickness (m)
#     L_oxide : float, optional
#         Override oxide thickness (m)
#     microstructure_overrides : dict, optional
#         Override specific microstructure parameters
#     defect_overrides : dict, optional
#         Override specific defect parameters
        
#     Returns
#     -------
#     dict
#         Complete configuration for simulation

#     Notes
#     -----
#     T_ref for Arrhenius calculations is stored WITHIN each material/oxide dict,
#     NOT in the operating conditions. Each material may have a different T_ref.
#     """
#     config = {
#         # Material selections
#         'metal_name': metal,
#         'oxide_name': oxide,
#         'metal_props': METALS[metal].copy(),
#         'oxide_props': OXIDES[oxide].copy(),
        
#         # Operating conditions
#         'T_operating': T_operating or CONDITIONS['T_operating'],
#         'P_upstream': P_upstream or CONDITIONS['P_upstream'],
#         'P_downstream': CONDITIONS['P_downstream'],
        
#         # Geometry
#         'L_metal': L_metal or CONDITIONS['L_metal'],
#         'L_oxide': L_oxide or OXIDES[oxide].get('thickness_default', CONDITIONS['L_oxide']),
        
#         # Microstructure
#         'microstructure': MICROSTRUCTURE.copy(),
        
#         # Oxide defects
#         'oxide_defects': OXIDE_DEFECTS.copy(),
        
#         # Sweep ranges
#         'T_range': CONDITIONS['T_range'],
#         'P_range': CONDITIONS['P_range'],
#         'L_metal_range': CONDITIONS['L_metal_range'],
#         'L_oxide_range': CONDITIONS['L_oxide_range'],
#         'n_T_points': CONDITIONS['n_T_points'],
#         'n_P_points': CONDITIONS['n_P_points'],
#         'n_L_points': CONDITIONS['n_L_points'],
#     }
    
#     # Apply overrides
#     if microstructure_overrides:
#         config['microstructure'].update(microstructure_overrides)
#     if defect_overrides:
#         config['oxide_defects'].update(defect_overrides)
    
#     return config


# # =============================================================================
# # VALIDATION / TEST PARAMETERS
# # =============================================================================
# """
# Standardized parameters for model validation and testing.
# Organized by hierarchical level to ensure consistency across all simulations.
# """

# VALIDATION = {
#     # -------------------------------------------------------------------------
#     # LEVEL 1: Perfect Metal Validation
#     # -------------------------------------------------------------------------
#     'L1': {
#         'P_ref': 1.0,                                    # Pa - reference pressure for thickness tests
#         'n_test': 50,                                    # number of test points
#         'test_pressures': np.logspace(-2, 2, 50),        # Pa - pressure sweep range
#         'test_temps_C': np.array([700, 775, 850, 925]),  # °C - 4 test temperatures
#         'test_temps_K': np.array([700, 775, 850, 925]) + 273.15,  # K - should match test_temps_C
#     },
    
#     # -------------------------------------------------------------------------
#     # LEVEL 2a: Perfect Oxide Only Validation
#     # -------------------------------------------------------------------------
#     'L2a': {
#         'P_ref': 1.0,                                    # Pa - reference pressure
#         'n_test': 50,                                    # number of test points
#         'test_pressures': np.logspace(-2, 5, 50),        # Pa - wider range for oxide
#         'test_temps_C': np.array([700, 800, 900, 1000]), # °C
#         'test_temps_K': np.array([700, 800, 900, 1000]) + 273.15,  # K
#     },
    
#     # -------------------------------------------------------------------------
#     # LEVEL 2b: Perfect Oxide + Perfect Metal Validation
#     # -------------------------------------------------------------------------
#     'L2b': {
#         # Temperature sweep parameters
#         'P_fixed': 1e6,                                  # Pa (1 MPa) - for T sweep
        
#         # Metal thickness sweep parameters
#         'P_fixed_metal_sweep': 1e8,                      # Pa (100 MPa) - high P for metal regime
#         'oxide_thickness_thin': 1e-8,                    # m (10 nm) - thin oxide for transition
#         'L_metal_sweep': np.logspace(-5, -1, 30),        # m (0.01 mm to 100 mm)
        
#         # Oxide thickness sweep (limit check)
#         'P_fixed_oxide_sweep': 1e6,                      # Pa (1 MPa)
#         'delta_ox_sweep': np.logspace(-4, -12, 40),      # m (100 μm to 1 pm)
#     },
#     # -------------------------------------------------------------------------
#     # LEVEL 3: Defective Oxide + Perfect Metal Validation
#     # -------------------------------------------------------------------------
#     'L3': {
#         # Fixed pressure for defect analysis
#         'P_fixed': 1e6,                                  # Pa (1 MPa)
#         'P_down': 0,                                     # Pa (vacuum)
        
#         # Defect fraction comparisons (Panel A)
#         'defect_fractions_compare': [0.0, 0.001, 0.01, 0.1],  # 0%, 0.1%, 1%, 10%
        
#         # Defect fraction sweep (Panel B, C)
#         'f_defect_min': 1e-5,                           # 0.001% minimum
#         'f_defect_max': 0.5,                            # 50% maximum (~10^-0.3)
#         'n_defect_points': 40,                          # number of points in sweep
        
#         # Limit check (Panel D)
#         'f_defect_limit_check': 1e-10,                  # Tiny defect fraction for validation
#     },
    
#     # -------------------------------------------------------------------------
#     # LEVEL 4: Defective Metal Validation
#     # -------------------------------------------------------------------------
#     'L4': {
#         # Pressure 
#         'P_fixed': 1e3,                                  # Pa (1 MPa) - fixed pressure for L4 validation
#         'P_down': 0,                                    # Pa (vacuum) - downstream pressure

#         # Mode comparison (Panel A)
#         'T_mode_comparison': 873,                        # K (600°C) - shows trapping effect
        
#         # Temperature crossover (Panel B)
#         'T_min': 473,                                    # K (200°C) - strong trapping
#         'T_max': 1173,                                   # K (900°C) - GB dominates
#         'n_T_points': 40,                                # temperature sweep points
        
#         # Grain size sensitivity (Panel C)
#         'T_grain_size': 873,                             # K (600°C)
#         'grain_size': np.logspace(-8, -3, 30),          # m (10 nm to 1 mm)
#         'f_gb_max': 0.5,                                 # maximum GB volume fraction
#     },
    
#     # -------------------------------------------------------------------------
#     # LEVEL 5: Full System Validation
#     # -------------------------------------------------------------------------
#     'L5': {
#         # All levels comparison (Panel A)
#         'P_down': 0,                                    # Pa (vacuum) - downstream pressure
#         'T_comparison': 873,                             # K (600°C)
        
#         # Arrhenius plot (Panel B)
#         'T_min': 573,                                    # K (300°C)
#         'T_max': 1173,                                   # K (900°C)
#         'n_T_points': 20,                                # temperature points
#         'P_fixed': 1e5,                                  # Pa (100 kPa)
        
#         # Sensitivity map (Panel C)
#         'T_sensitivity': 673,                            # K (400°C)
#         'P_sensitivity': 1e5,                            # Pa (100 kPa)
#         'f_defect_sweep': np.array([0.001, 0.01, 0.05, 0.1]),  # 0.1%, 1%, 5%, 10%
        
#         # Limit checks (Panel D)
#         'T_limit': 873,                            # K (600°C)
#     },
# }

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
R = 8.314          # J/mol/K - Universal gas constant
F = 96485          # C/mol - Faraday constant (for eV conversions)

# =============================================================================
# METAL PROPERTIES (Level 1, 4)
# =============================================================================
METALS = {
    'Fe_alpha': {
        # Reference conditions
        'T_ref': 625,                # K (351.85°C)
        'D_ref': 1.000e-11,          # m²/s at T_ref
        'K_s_ref': 1.070e-05,        # mol/m³/Pa^0.5 at T_ref
        
        # Activation energies (from literature)
        'E_D': 54000,                # J/mol (diffusion)
        'H_s': 5900,                 # J/mol (solution enthalpy)
        
        # Pre-exponential factors (calculated to match target Φ)
        'D_0': 1.292e-08,            # m²/s
        'K_s0': 1.822e-02,           # mol/m³/Pa^0.5
        
        # Target permeability verification
        'Phi_ref': 1.070e-08,        # mol/m/s/Pa^0.5 at T_ref (= D_ref × K_s_ref)
        
        # Metadata
        'reference': 'Calculated to match target permeability (activation energies from San Marchi 2007)',
        'temp_range_C': [149.87, 626.87],
    },

    'metal_316L_Heat_treated_ref_cast': {
        #'Diffusion parameters': {
            'T_ref': 873,           # K
            'D_ref': 2.194e-10,     # m²/s at T_ref
            'E_D':   42500,         # J/mol
            'D_0':   7.66e-08,      # m²/s
        #},
        #'Solubility parameters': {
            'K_s_ref': 2.723e-2,     # mol/m³/Pa^0.5 at T_ref
            'K_s0':   4.640e-1,     # mol/m³/Pa^0.5
            'H_s':    20585,        # J/mol
        #},
        #'Permeation parameters': {
            'Phi_ref': 1.881e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
            'Q_p':     63085,       # J/mol (= E_D + H_s)
            'Phi_0':   1.12e-07,    # mol/m/s/Pa^0.5 (= D_0 * K_s0)
        #},
        #'Other common parameters': {
            'pressure': [],
            'reference': 'Forcey et al. 1988 — Heat treated reference cast 316L austenitic steel',
            'temp_range': [523, 873],           # K
            'pressure_range': [1.33e2, 1e5],    # Pa
            'metal_thickness': [5e-3],          # m
            'gas': 'H₂',
            'Notes': 'Data in SI units. Phi_ref = D_ref * Ks_ref. Ks_0 = Phi_0/D_0 = 1.462.',
        #},
    }
    
}

# =============================================================================
# OXIDE PROPERTIES (Level 2a, 2b, 3)
# =============================================================================
OXIDES = {
    # 'Cr2O3': {
    #     # Reference conditions
    #     'T_ref': 625,                # K (351.85°C) - same as metal for consistency
    #     'D_ox_ref': 1.000e-18,       # m²/s at T_ref
    #     'K_ox_ref': 1.070e-07,       # mol/m³/Pa at T_ref (Henry's law)
        
    #     # Activation energies (from literature)
    #     'E_D_ox': 1.55e5,            # J/mol (155 kJ/mol - diffusion)
    #     'H_sol_ox': 1.85e5,          # J/mol (185 kJ/mol - solution)
        
    #     # Pre-exponential factors (calculated to match target Φ_ox)
    #     'D_ox_0': 8.414e+10,         # m²/s
    #     'K_ox_0': 1.514e+16,         # mol/m³/Pa
        
    #     # Target permeability verification
    #     'Phi_ox_ref': 1.070e-18,     # mol/m/s/Pa at T_ref (= D_ox_ref × K_ox_ref)
        
    #     # Default geometry
    #     'thickness': 1e-12,           # m (1 μm)
    #     'thickness_range': [1e-7, 1e-5],
        
    #     # Metadata
    #     'reference': 'Calculated to match target permeability (activation energies from Strehlow & Savage 1974, Serra 1998)',
    #     'temp_range_K': [625, 1273],  # Extended to lower T
    #     'temperature_range': [625, 1273],
    #     'uncertainty_factor': 10,
    # },
    
        'Al2O3': {
        # ---- Reference conditions ----
        'T_ref':       1623,        # K (1350°C) — centre of Roberts 1979 range [R79]
        'D_ox_ref':    3.0647e-11,  # m²/s — Belonoshko 2004 [B04], Step 3 above
                                    # Eq: D_ox_0 × exp(−E_D / R T_ref)
        'K_ox_ref':    4.9077,      # mol/m³/Pa^0.5 — Sieverts re-assumed at 1 atm, Step 7
                                    # Eq: (Phi_R_ref / D_ox_ref) × 101325^(−0.07)

        # ---- Activation energies ----
        'E_D_ox':      119641,      # J/mol (1.24 eV) — Belonoshko 2004 DFT [B04 via C11]
                                    # Cross-validated vs Chen 2011 Table 3: within 8% ✓
        'H_sol_ox':    198559,      # J/mol — Step 2: Q_p − E_D = 318200 − 119641
                                    # Q_p from Roberts 1979 Eq.(4) [R79]; E_D from [B04]

        # ---- Pre-exponential factors (Sieverts, Pa^0.5) ----
        'D_ox_0':      2.1730e-7,   # m²/s — Belonoshko 2004 [B04 via C11]
                                    # Eq: D(T) = D_ox_0 × exp(−E_D_ox / RT)
        'K_ox_0':      1.2065e+7,   # mol/m³/Pa^0.5 — Sieverts re-assumed at 1 atm, Step 7
                                    # Eq: K_R_0 × 101325^(−0.07)
                                    # where K_R_0 = K_R_ref × exp(+H_sol / R T_ref)

        # ---- Permeability at T_ref (Sieverts, Pa^0.5) ----
        'Phi_ox_ref':  1.5040e-10,  # mol/m/s/Pa^0.5 — Sieverts re-assumed at 1 atm, Step 7
                                    # Eq: Phi_R_ref × 101325^(−0.07)
                                    # VERIFY: D_ox_ref × K_ox_ref = 3.065e-11 × 4.908 = 1.504e-10 ✓
        'Phi_ox_0':    2.6216,      # mol/m/s/Pa^0.5 — Sieverts re-assumed at 1 atm, Step 7
                                    # Eq: Phi_0_SI × 101325^(−0.07) = 5.8744 × 0.4463

        # ---- Raw Roberts values (Pa^0.43) — retained for traceability ----
        'Phi_ox_0_raw_Pa043':  5.8744,      # mol/m/s/Pa^0.43 — Roberts 1979 Eq.(4) converted to SI
        'K_ox_ref_raw_Pa043':  1.0997e+1,   # mol/m³/Pa^0.43  — before pressure re-assumption
        'K_ox_0_raw_Pa043':    2.7034e+7,   # mol/m³/Pa^0.43  — before pressure re-assumption
        'Phi_ox_ref_raw_Pa043':3.3701e-10,  # mol/m/s/Pa^0.43 — before pressure re-assumption
        'pressure_exponent_Roberts': 0.43,  # measured in Roberts 1979 Fig 3 at 1350°C, 2–50 kPa
        'P_reassumption_Pa':   101325,      # Pa — reference pressure used for Sieverts re-assumption

        # ---- Default geometry ----
        'thickness':       1e-6,    # m (1 µm) — typical CVD/PVD Al2O3 barrier coating [H95]
        'thickness_range': [1e-7, 1e-5],  # m — 0.1 µm to 10 µm [H95]

        # ---- Metadata ----
        'reference': (
            'Phi_ox_0, Q_p: Roberts RM, Elleman TS, Palmour H III, Verghese K. '
            '"Hydrogen Permeability of Sintered Aluminum Oxide." '
            'J. Am. Ceram. Soc. 62(9–10), 495–499 (1979). Eq.(4)/(5). '
            'Directly from original paper. '
            'Specimen: sintered Al2O3 tubes, 99.8% purity, SG=3.85, gastight. '
            'Measurement: T=1200–1450°C, P=2–50 kPa, H2+T2 tracer. '
            'E_D_ox (1.24 eV): Belonoshko et al. 2004, Phys. Rev. B 69, 024302. '
            'Cited as Ref 13 in Chen CF et al. 2011. DOI: 10.1007/s11431-010-4112-3. '
            'H_sol_ox: derived Q_p − E_D = 318200 − 119641 = 198559 J/mol. '
            'K_ox_0, K_ox_ref, Phi_ox_ref: Sieverts re-assumed at P_ref=101325 Pa '
            'using scale factor = 101325^(0.43−0.50) = 0.4463. '
            'Full equation chain documented in entry header above.'
        ),
        'temp_range_K':       [1473, 1723],  # K — Roberts 1979 (1200–1450°C); extrapolation beyond is unvalidated
        'temperature_range':  [1473, 1723],
        'uncertainty_factor': 10,
        # Justification for ×10:
        # Phi_ox_0 and Q_p from Roberts 1979 Eq.(4) directly — no estimation.
        # E_D_ox from Belonoshko 2004 DFT, cross-validated vs Roberts D(T) via Chen 2011.
        # H_sol_ox exact (arithmetic from two primary sources).
        # Remaining uncertainties:
        # (1) Sieverts re-assumption introduces ~10% error at pressures far from 101325 Pa —
        #     exact error = |P^0.43 / P^0.5 - (P/P_ref)^(-0.07)| grows with |P - P_ref|.
        # (2) Roberts data at 1200–1450°C only; lower-T extrapolation unvalidated.
        # (3) Roberts warns (Section IV): microstructural changes expected above 1450°C.
        # (4) Hollenberg 1995: grain-boundary diffusion at lower T may raise Phi by
        #     "several orders of magnitude" relative to bulk extrapolation.

        'notes': (
            '=== EQUATION CHAIN (summary) ===\n'
            'All steps numbered to match entry header.\n'
            '\n'
            'Step 0 — Roberts 1979 Eq.(4), paper units:\n'
            '  φ = exp(48.95±0.61) × exp(−318200±18800 / RT)\n'
            '  units: H-atom·cm·cm⁻²·s⁻¹·kPa⁻⁰·⁴³\n'
            '\n'
            'Step 1 — Unit conversion Roberts → SI (Pa^0.43):\n'
            '  Phi_0_SI = exp(48.95) / N_A × (cm→m) × (cm⁻²→m⁻²) × (kPa→Pa)^0.43\n'
            '           = 1.814e21 / 6.022e23 × 1e-2 × 1e4 × (1e3)^0.43\n'
            '           = 5.8744 mol/m/s/Pa^0.43\n'
            '\n'
            'Step 2 — H_sol:\n'
            '  H_sol = Q_p − E_D = 318200 − 119641 = 198559 J/mol\n'
            '  Q_p from Roberts 1979 Eq.(4); E_D from Belonoshko 2004 (1.24 eV)\n'
            '\n'
            'Step 3 — D_ox_ref at T_ref = 1623 K:\n'
            '  D_ox_ref = 2.173e-7 × exp(−119641 / (8.314 × 1623)) = 3.065e-11 m²/s\n'
            '\n'
            'Step 4 — Phi_R_ref at T_ref (Pa^0.43):\n'
            '  Phi_R_ref = 5.8744 × exp(−318200 / (8.314 × 1623)) = 3.370e-10 mol/m/s/Pa^0.43\n'
            '\n'
            'Step 5 — K_R_ref (Pa^0.43):\n'
            '  K_R_ref = Phi_R_ref / D_ox_ref = 3.370e-10 / 3.065e-11 = 10.997 mol/m³/Pa^0.43\n'
            '\n'
            'Step 6 — K_R_0 (Pa^0.43):\n'
            '  K_R_0 = K_R_ref × exp(+H_sol / R T_ref)\n'
            '        = 10.997 × exp(198559 / (8.314 × 1623)) = 2.703e7 mol/m³/Pa^0.43\n'
            '\n'
            'Step 7 — Sieverts re-assumption at P_ref = 101325 Pa:\n'
            '  Derivation: equate fluxes J = Phi_R × P^0.43 / d = Phi_S × P^0.5 / d\n'
            '  At P = P_ref: Phi_S = Phi_R × P_ref^(0.43−0.50) = Phi_R × P_ref^(−0.07)\n'
            '  scale = 101325^(−0.07) = exp(−0.07 × 11.526) = exp(−0.8068) = 0.4463\n'
            '  Applied uniformly to Phi_0, K_ref, K_0 — same T-independent scale factor.\n'
            '  E_D_ox and H_sol_ox are NOT scaled (they live in the exponent).\n'
            '  Phi_S_0   = 5.8744  × 0.4463 = 2.6216  mol/m/s/Pa^0.5\n'
            '  K_S_ref   = 10.997  × 0.4463 = 4.9077  mol/m³/Pa^0.5\n'
            '  K_S_0     = 2.703e7 × 0.4463 = 1.207e7 mol/m³/Pa^0.5\n'
            '  Phi_S_ref = 3.370e-10 × 0.4463 = 1.504e-10 mol/m/s/Pa^0.5\n'
            '\n'
            'Step 8 — Identity checks (all verified at 0.0000% error):\n'
            '  D_ox_ref × K_S_ref = 3.065e-11 × 4.908 = 1.504e-10 = Phi_S_ref ✓\n'
            '  E_D + H_sol = 119641 + 198559 = 318200 = Q_p ✓\n'
            '  K_S_0 × exp(−H_sol/RT_ref) = K_S_ref ✓\n'
            '  Spot check 1473K: D×K = Phi_Roberts ✓\n'
            '  Spot check 1723K: D×K = Phi_Roberts ✓\n'
            '  Spot check  973K: D×K = Phi_Roberts ✓\n'
            '\n'
            '=== PRESSURE RE-ASSUMPTION APPLICABILITY ===\n'
            'Valid when operating pressure P is near P_ref = 101325 Pa (1 atm).\n'
            'Error from re-assumption grows as |P − P_ref| increases:\n'
            '  At P = 10 kPa  (Roberts lower range): scale = 0.5248, error ~18%\n'
            '  At P = 50 kPa:  scale = 0.4796, error ~7%\n'
            '  At P = 101 kPa: scale = 0.4463, error = 0% (by definition)\n'
            '  At P = 500 kPa: scale = 0.4073, error ~9%\n'
            'For high-accuracy work at pressures far from 1 atm, use raw Pa^0.43\n'
            'values (stored as Phi_ox_0_raw_Pa043, K_ox_ref_raw_Pa043, etc.).\n'
            '\n'
            '=== TEMPERATURE EXTRAPOLATION WARNING ===\n'
            'Roberts 1979 data range: 1200–1450°C (1473–1723 K).\n'
            'Roberts (Section IV) explicitly warns against extrapolation above 1450°C\n'
            'due to microstructural changes (bloating, grain growth).\n'
            'Hollenberg 1995: extrapolation to lower T may underestimate Phi by\n'
            '"several orders of magnitude" due to grain-boundary diffusion crossover.\n'
            'Roberts conclusion 3 confirms: at 1200–1450°C, grain boundaries do NOT\n'
            'contribute — making this a clean defect-free anchor in that range.\n'
            '\n'
            '=== D(T) CROSS-VALIDATION (Belonoshko 2004 vs Chen 2011 Table 3) ===\n'
            '  273K: D_Bel=2.78e-30, D_Chen=3.04e-30 m²/s  (8% diff)\n'
            '  773K: D_Bel=1.79e-15, D_Chen=1.84e-15 m²/s  (3% diff)\n'
            ' 1273K: D_Bel=2.68e-12, D_Chen=2.73e-12 m²/s  (2% diff)\n'
            '\n'
            '=== DEFECT-FREE MODEL APPLICABILITY ===\n'
            'Roberts (conclusion 3): grain boundaries and closed pores (<1%)\n'
            'do not serve as rapid diffusion paths — consistent with defect-free assumption.\n'
            'Hollenberg 1995 area-defect caveat (Q_p ≈ 65 kJ/mol in aluminized SS316)\n'
            'does NOT apply here — that result reflects coating defects in real systems,\n'
            'not the intrinsic Al2O3 bulk property used in this entry.'
        ),
    }

}

# =============================================================================
# MICROSTRUCTURE PARAMETERS (Level 4)
# =============================================================================
# MICROSTRUCTURE = {
#     # Grain boundary parameters
#     'grain_size': 50e-6,             # m (50 μm default, use 100e-9 for nanocrystalline)
#     'gb_thickness': 0.5e-9,          # m (0.5 nm)
#     'grain_shape': 'equiaxed',       # 'equiaxed', 'columnar'
#     'gb_type': 'LAGB',               # 'HAGB', 'LAGB'
#     'include_gb_trapping': False,
    
#     # Trap parameters
#     'trap_list': [
#         {
#             'name': 'vacancies',
#             'binding_energy': 0.3 * F,   # J/mol (0.5 eV)
#             'density': 1e26,              # m⁻³
#         },
#     ],
    
#     # Lattice site density (for Oriani model)
#     'N_L': 8.46e28,                  # m⁻³ (Fe BCC lattice sites)
# }

# Extended version — multiple traps simultaneously
MICROSTRUCTURE = {
    'grain_size': 50e-6,
    'gb_thickness': 0.5e-9,
    'grain_shape': 'equiaxed',
    'gb_type': 'LAGB',
    'include_gb_trapping': False,

    'trap_list': [
        {
            'name': 'vacancies',
            'binding_energy': 0.5 * F,   # J/mol (0.5 eV) — weak trap
            'density': 1e26,              # m⁻³
        },
        {
            'name': 'dislocations',
            'binding_energy': 0.7 * F,   # J/mol (0.7 eV) — medium trap
            'density': 1e22,              # m⁻³ (line density × length)
        },
        {
            'name': 'grain_boundaries',
            'binding_energy': 0.9 * F,   # J/mol (0.9 eV) — strong trap
            'density': 1e24,              # m⁻³
        },
    ],

    'N_L': 8.46e28,                  # m⁻³ (Fe BCC lattice sites)
}

# =============================================================================
# OXIDE DEFECT PARAMETERS (Level 3, 5)
# =============================================================================
# OXIDE_DEFECTS = {
#     'area_fraction': 0.01,           # 1% pinholes
#     'type': 'pinhole',               # 'pinhole', 'crack'
# }
OXIDE_DEFECTS = {
    'area_fraction': 0.015,   # Total defect fraction (sum of all components)
    'type': 'mixed',
    'components': {
        'pinholes': 0.005,        # 0.5% — no oxide barrier
        'cracks': 0.005,          # 0.5% — thin oxide path
        'grain_boundaries': 0.005 # 0.5% — enhanced diffusion
    },
    'thickness_factor': 0.1,   # For crack component: L_crack = 0.1 × L_oxide
    'diffusivity_factor': 10   # For GB component: D_gb = 10 × D_oxide
}
# =============================================================================
# OPERATING CONDITIONS
# =============================================================================
CONDITIONS = {
    # Temperature
    'T_operating': 625,              # K - Reference temperature (NOT T_ref for Arrhenius!)
    'T_range': (423, 1200),           # K - Sweep range
    'n_T_points': 20,
    
    # Pressure
    'P_upstream': 1e4,               # Pa (10 kPa)
    'P_downstream': 0,               # Pa
    'P_range': (1e-8, 1e8),           # Pa - Sweep range
    'n_P_points': 100,
    
    # Geometry
    'L_metal': 1e-3,                 # m (1 mm)
    'L_oxide': 1e-6,                 # m (1 μm)
    'L_metal_range': (0.1e-3, 5e-3), # m
    'L_oxide_range': (1e-7, 1e-5),   # m
    'n_L_points': 20,
}

# =============================================================================
# SURFACE KINETICS PARAMETERS (Level 6 - if used)
# =============================================================================
SURFACE_KINETICS = {
    'k_diss': 1e-15,                 # mol/m²/s/Pa
    'k_recomb': 1e-3,                # m⁴/mol/s
    'coverage_mode': 'steady_state', # 'steady_state', 'langmuir'
}


# =============================================================================
# PLOTTING STYLE
# =============================================================================
PLOT_STYLE = {
    'figsize': (14, 12),
    'fontsize_title': 16,
    'fontsize_suptitle': 18,
    'fontsize_axis': 14,
    'fontsize_tick': 12,
    'fontsize_legend': 12,
    'fontsize_annotation': 11,
    'linewidth': 2.5,
    'markersize': 8,
    'grid_alpha': 0.3,
}

# Color scheme by level
COLORS = {
    'L1': 'black',           # Perfect Metal (baseline)
    'L2a': 'blue',           # Perfect Oxide Only
    'L2b': 'purple',         # Oxide + Metal
    'L3': 'cyan',            # Defective Oxide + Metal
    'L4_gb': 'green',        # Defective Metal - GB mode
    'L4_trap': 'orange',     # Defective Metal - Trapping mode
    'L4_both': 'crimson',    # Defective Metal - Combined
    'L5': 'red',             # Full System
}


# =============================================================================
# HELPER FUNCTION: Build complete config dict for a simulation
# =============================================================================
def build_simulation_config(
    metal='None',
    oxide='None',
    T_operating=None,
    P_upstream=None,
    L_metal=None,
    L_oxide=None,
    microstructure_overrides=None,
    defect_overrides=None,
):
    """
    Build a complete configuration dictionary for a simulation.
    
    Parameters
    ----------
    metal : str
        Key from METALS dict
    oxide : str
        Key from OXIDES dict
    T_operating : float, optional
        Operating temperature [K] (NOT the reference temperature for properties)
    P_upstream : float, optional
        Override upstream pressure (Pa)
    L_metal : float, optional
        Override metal thickness (m)
    L_oxide : float, optional
        Override oxide thickness (m)
    microstructure_overrides : dict, optional
        Override specific microstructure parameters
    defect_overrides : dict, optional
        Override specific defect parameters
        
    Returns
    -------
    dict
        Complete configuration for simulation

    Notes
    -----
    T_ref for Arrhenius calculations is stored WITHIN each material/oxide dict,
    NOT in the operating conditions. Each material may have a different T_ref.
    """
    config = {
        # Material selections
        'metal_name': metal,
        'oxide_name': oxide,
        'metal_props': METALS[metal].copy(),
        'oxide_props': OXIDES[oxide].copy(),
        
        # Operating conditions
        'T_operating': T_operating or CONDITIONS['T_operating'],
        'P_upstream': P_upstream or CONDITIONS['P_upstream'],
        'P_downstream': CONDITIONS['P_downstream'],
        
        # Geometry
        'L_metal': L_metal or CONDITIONS['L_metal'],
        'L_oxide': L_oxide or OXIDES[oxide].get('thickness_default', CONDITIONS['L_oxide']),
        
        # Microstructure
        'microstructure': MICROSTRUCTURE.copy(),
        
        # Oxide defects
        'oxide_defects': OXIDE_DEFECTS.copy(),
        
        # Sweep ranges
        'T_range': CONDITIONS['T_range'],
        'P_range': CONDITIONS['P_range'],
        'L_metal_range': CONDITIONS['L_metal_range'],
        'L_oxide_range': CONDITIONS['L_oxide_range'],
        'n_T_points': CONDITIONS['n_T_points'],
        'n_P_points': CONDITIONS['n_P_points'],
        'n_L_points': CONDITIONS['n_L_points'],
    }
    
    # Apply overrides
    if microstructure_overrides:
        config['microstructure'].update(microstructure_overrides)
    if defect_overrides:
        config['oxide_defects'].update(defect_overrides)
    
    return config


# =============================================================================
# VALIDATION / TEST PARAMETERS
# =============================================================================
"""
Standardized parameters for model validation and testing.
Organized by hierarchical level to ensure consistency across all simulations.
"""

VALIDATION = {
    # -------------------------------------------------------------------------
    # LEVEL 1: Perfect Metal Validation
    # -------------------------------------------------------------------------
    'L1': {
        'P_ref': 1.0,                                    # Pa - reference pressure for thickness tests
        'n_test': 50,                                    # number of test points
        'test_pressures': np.logspace(-2, 2, 50),        # Pa - pressure sweep range
        'test_temps_C': np.array([700, 775, 850, 925]),  # °C - 4 test temperatures
        'test_temps_K': np.array([700, 775, 850, 925]) + 273.15,  # K - should match test_temps_C
    },
    
    # -------------------------------------------------------------------------
    # LEVEL 2a: Perfect Oxide Only Validation
    # -------------------------------------------------------------------------
    'L2a': {
        'P_ref': 1.0,                                    # Pa - reference pressure
        'n_test': 50,                                    # number of test points
        'test_pressures': np.logspace(-2, 5, 50),        # Pa - wider range for oxide
        'test_temps_C': np.array([700, 800, 900, 1000]), # °C
        'test_temps_K': np.array([700, 800, 900, 1000]) + 273.15,  # K
    },
    
    # -------------------------------------------------------------------------
    # LEVEL 2b: Perfect Oxide + Perfect Metal Validation
    # -------------------------------------------------------------------------
    'L2b': {
        # Temperature sweep parameters
        'P_fixed': 1e3,                                  # Pa (1 MPa) - for T sweep
        
        # Metal thickness sweep parameters
        'P_fixed_metal_sweep': 1e3,                      # Pa (100 MPa) - high P for metal regime
        'oxide_thickness_thin': 1e-8,                    # m (10 nm) - thin oxide for transition
        'L_metal_sweep': np.logspace(-5, -1, 30),        # m (0.01 mm to 100 mm)
        
        # Oxide thickness sweep (limit check)
        'P_fixed_oxide_sweep': 1e3,                      # Pa (1 MPa)
        'delta_ox_sweep': np.logspace(-4, -12, 40),      # m (100 μm to 1 pm)
    },

    # -------------------------------------------------------------------------
    # LEVEL 3: Defective Oxide + Perfect Metal Validation
    # -------------------------------------------------------------------------
    'L3': {
        # Fixed pressure for defect analysis
        'P_fixed': 1e3,                                  # Pa (1 kPa)
        'P_down': 0,                                     # Pa (vacuum)
        
        # Defect fraction comparisons (Panel A)
        'defect_fractions_compare': [0.0, 0.001, 0.01, 0.1],  # 0%, 0.1%, 1%, 10%
        
        # Defect fraction sweep (Panel B, C)
        'f_defect_min': 1e-5,                           # 0.001% minimum
        'f_defect_max': 0.5,                            # 50% maximum (~10^-0.3)
        'n_defect_points': 40,                          # number of points in sweep
        
        # Limit check (Panel D)
        'f_defect_limit_check': 1e-100,                  # Tiny defect fraction for validation
    },
    
    # -------------------------------------------------------------------------
    # LEVEL 4: Defective Metal Validation
    # -------------------------------------------------------------------------
    'L4': {
        # Pressure 
        'P_fixed': 1e3,                                  # Pa (1 MPa) - fixed pressure for L4 validation
        'P_down': 0,                                     # Pa (vacuum) - downstream pressure

        # Mode comparison (Panel A)
        'T_mode_comparison': 873,                        # K (600°C) - shows trapping effect
        
        # Temperature crossover (Panel B)
        'T_min': 473,                                    # K (200°C) - strong trapping
        'T_max': 1173,                                   # K (900°C) - GB dominates
        'n_T_points': 40,                                # temperature sweep points
        
        # Grain size sensitivity (Panel C)
        'T_grain_size': 873,                             # K (600°C)
        'grain_size': np.logspace(-8, -3, 30),           # m (10 nm to 1 mm)
        'f_gb_max': 0.5,                                 # maximum GB volume fraction
    },
    
    # -------------------------------------------------------------------------
    # LEVEL 5: Full System Validation
    # -------------------------------------------------------------------------
    'L5': {
        # All levels comparison (Panel A)
        'P_down': 0,                                     # Pa (vacuum) - downstream pressure
        'T_comparison': 873,                             # K (600°C)
        
        # Arrhenius plot (Panel B)
        'T_min': 573,                                    # K (300°C)
        'T_max': 1173,                                   # K (900°C)
        'n_T_points': 20,                                # temperature points
        'P_fixed': 1e5,                                  # Pa (100 kPa)
        
        # Sensitivity map (Panel C)
        'T_sensitivity': 673,                            # K (400°C)
        'P_sensitivity': 1e5,                            # Pa (100 kPa)
        'f_defect_sweep': np.array([0.001, 0.01, 0.05, 0.1]),  # 0.1%, 1%, 5%, 10%
        
        # Limit checks (Panel D)
        'T_limit': 873,                                  # K (600°C)
    },
}


# =============================================================================
# TEMPLATE
# =============================================================================
"""
'metal_NAME': {
    'Diffusion parameters': {
        'T_ref': ,        # K
        'D_ref': ,        # m²/s at T_ref
        'E_D':   ,        # J/mol (diffusion activation energy)
        'D_0':   ,        # m²/s (pre-exponential factor)
    },
    'Solubility parameters': {
        'Ks_ref': ,       # mol/m³/Pa^0.5 at T_ref
        'K_s0':   ,       # mol/m³/Pa^0.5 (pre-exponential factor)
        'H_s':    ,       # J/mol (heat of solution)
    },
    'Permeation parameters': {
        'Phi_ref': ,      # mol/m/s/Pa^0.5 at T_ref
        'Q_p':     ,      # J/mol (permeation activation energy)
        'Phi_0':   ,      # mol/m/s/Pa^0.5 (pre-exponential factor)
    },
    'Other common parameters': {
        'pressure': [],           # Pa
        'reference': '',
        'temp_range': [],         # K validity range
        'pressure_range': [],     # Pa
        'metal_thickness': [],    # m
        'gas': '',                # e.g. H₂, D₂, T₂
        'Notes': '',
    },
},

'oxide_Cr2O3_NAME': {
    'Diffusion parameters': {
        'T_ref': ,        # K
        'D_ref': ,        # m²/s at T_ref
        'E_D':   ,        # J/mol (diffusion activation energy)
        'D_0':   ,        # m²/s (pre-exponential factor)
    },
    'Solubility parameters': {
        'Ks_ref': ,       # mol/m³/Pa^0.5 at T_ref
        'K_s0':   ,       # mol/m³/Pa^0.5 (pre-exponential factor)
        'H_s':    ,       # J/mol (heat of solution)
    },
    'Permeation parameters': {
        'Phi_ref': ,      # mol/m/s/Pa^0.5 at T_ref
        'Q_p':     ,      # J/mol (permeation activation energy)
        'Phi_0':   ,      # mol/m/s/Pa^0.5 (pre-exponential factor)
    },
    'Permeation Reduction Factor': {
        'PRF':     ,      # dimensionless (Phi_clean_metal / Phi_oxide_coated)
        'PRF_ref': ,      # K (temperature at which PRF is reported)
    },
    'Other common parameters': {
        'oxide_thickness': [],    # m
        'pressure': [],           # Pa
        'reference': '',
        'temp_range': [],         # K validity range
        'pressure_range': [],     # Pa
        'substrate': '',          # e.g. 316L, Incoloy800
        'gas': '',                # e.g. H₂, D₂, T₂
        'Notes': '',
    },
},

'surface_kinetics_NAME': {
    'Dissociation (upstream surface)': {
        'K_d':     ,      # m⁴/mol/s or appropriate units — dissociative adsorption rate constant at T_ref
        'K_d0':    ,      # pre-exponential factor
        'E_d':     ,      # J/mol (activation energy for dissociation)
        'T_ref':   ,      # K
    },
    'Recombination (downstream surface)': {
        'K_r':     ,      # m⁴/mol/s or appropriate units — recombinative desorption rate constant at T_ref
        'K_r0':    ,      # pre-exponential factor
        'E_r':     ,      # J/mol (activation energy for recombination)
        'T_ref':   ,      # K
    },
    'Sticking coefficient': {
        'S_0':     ,      # dimensionless (0–1)
        'E_s':     ,      # J/mol (activation energy, if temperature dependent)
    },
    'Other common parameters': {
        'pressure': [],           # Pa
        'reference': '',
        'temp_range': [],         # K validity range
        'pressure_range': [],     # Pa
        'surface_condition': '',  # e.g. clean, oxidized, coated
        'substrate': '',          # e.g. 316L, Ni, Incoloy800
        'gas': '',                # e.g. H₂, D₂, T₂
        'Notes': '',
    },
},





#### OXIDE LAYER PROPERTIES - LITERATURE DATA COLLECTION ####
# 
# DATA QUALITY TIERS:
# -------------------
# Tier 1 — Full back-calculation possible:
#           Phi_composite, Phi_bare_metal, oxide_thickness, metal_thickness all known
#           → standalone Phi_oxide can be computed via series resistance formula
#
# Tier 2 — Partial back-calculation possible:
#           PRF and oxide_thickness known, but Phi_bare_metal taken from METALS dict
#           → back-calculation possible but depends on choice of metal baseline
#
# Tier 3 — PRF only:
#           Only PRF reported, oxide_thickness unknown
#           → cannot back-calculate standalone Phi_oxide without assumptions
#
#
# BACK-CALCULATION NOTE (Series Resistance Model):
# -------------------------------------------------
# For a composite membrane (oxide + metal in series), permeation resistance adds as:
#
#   R_total   = R_oxide   + R_metal
#   d_total   = d_oxide   + d_metal          (thicknesses, m)
#
# In terms of permeability (Phi, mol/m/s/Pa^0.5):
#
#   d_total / Phi_total = d_oxide / Phi_oxide + d_metal / Phi_metal
#
# Therefore, to back-calculate standalone Phi_oxide:
#
#   Phi_oxide = d_oxide / (d_total/Phi_total - d_metal/Phi_metal)
#
# Where:
#   Phi_total  = Phi_composite_ref   (measured, oxide+metal together)
#   Phi_metal  = from METALS dict, key = metal_baseline_key, field = Phi_ref
#                evaluated at same T_ref using Arrhenius:
#                Phi_metal(T) = Phi_metal_ref * exp((-Q_p/R) * (1/T - 1/T_ref_metal))
#   d_oxide    = oxide_thickness     (m)
#   d_metal    = metal_thickness     (m)
#   d_total    = d_oxide + d_metal   (m)
#
# For Arrhenius parameters of standalone oxide (if composite Arrhenius is known):
#   This requires back-calculation at multiple temperatures to fit:
#   Phi_oxide(T) = Phi_oxide_0 * exp(-Q_p_oxide / RT)
#   → fit Phi_oxide(T) vs 1/T to extract Phi_oxide_0 and Q_p_oxide
#
# IMPORTANT ASSUMPTIONS:
#   - Oxide and metal are in series (1D diffusion, no lateral transport)
#   - Phi_metal baseline must be chosen carefully — use same gas, similar T range
#   - If oxide_thickness is not reported, Tier 3 only; no back-calculation possible
#   - PRF = Phi_bare / Phi_composite (dimensionless, always >= 1 for a barrier oxide)
#
# =============================================================================

OXIDES = {

    # =========================================================================
    # TEMPLATE — copy this block for each new oxide entry
    # =========================================================================

    'oxide_TYPE_on_SUBSTRATE_SOURCE': {
        'Oxide identity': {
            'oxide_type': '',               # e.g. 'Cr₂O₃', 'Al₂O₃', 'Fe₃O₄', 'SiO₂'
            'substrate': '',                # e.g. '316L', 'Incoloy800' — matches METALS dict key
            'formation_condition': '',      # e.g. 'air oxidized 800°C 1h', 'steam 500°C',
                                            #       'electrochemical', 'native/as-received'
            'oxide_thickness': None,        # m (None if not reported)
            'thickness_known': False,       # True / False
            'data_tier': None,              # 1, 2, or 3 — set based on what is available
        },
        'As-reported parameters (composite)': {
            'T_ref': None,                  # K — temperature at which Phi_composite_ref is reported
            'PRF': None,                    # dimensionless — Phi_bare / Phi_composite (>= 1)
            'PRF_T_ref': None,              # K — temperature at which PRF was evaluated
            'Phi_composite_ref': None,      # mol/m/s/Pa^0.5 at T_ref (oxide + metal together)
            'Q_p_composite': None,          # J/mol — permeation activation energy (composite)
            'Phi_composite_0': None,        # mol/m/s/Pa^0.5 — pre-exponential (composite)
            'Phi_bare_ref': None,           # mol/m/s/Pa^0.5 — bare metal baseline used for PRF
            'metal_baseline_key': None,     # str — key in METALS dict used as baseline
                                            #        e.g. 'metal_316L_Heat_treated_ref_cast'
        },
        'Back-calculated standalone oxide parameters': {
            # All None by default — see back-calculation note at top of file
            # Tier 1 or 2: compute using series resistance formula above
            # Tier 3: cannot be computed without oxide_thickness
            'Phi_oxide_ref': None,          # mol/m/s/Pa^0.5 at T_ref — back-calculated
            'Q_p_oxide': None,              # J/mol — requires multi-T fit (see note)
            'Phi_oxide_0': None,            # mol/m/s/Pa^0.5 — requires multi-T fit
            'D_ref': None,                  # m²/s — None if D and Ks not separable
            'E_D': None,                    # J/mol
            'D_0': None,                    # m²/s
            'Ks_ref': None,                 # mol/m³/Pa^0.5
            'K_s0': None,                   # mol/m³/Pa^0.5
            'H_s': None,                    # J/mol
            'back_calc_method': None,       # e.g. 'series resistance — Tier 1',
                                            #       'series resistance — Tier 2 (METALS baseline)',
                                            #       'not possible — Tier 3'
        },
        'Other common parameters': {
            'pressure': [],                 # Pa
            'reference': '',
            'temp_range': [],               # K validity range
            'pressure_range': [],           # Pa
            'metal_thickness': [],          # m — thickness of substrate
            'gas': '',                      # e.g. 'H₂', 'D₂', 'T₂'
            'Notes': '',
        },
    },

}
"""
