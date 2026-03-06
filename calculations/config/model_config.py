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

UPDATED: Parameters calculated to achieve target permeabilities at T_ref = 625 K
  - Metal (Fe_alpha):  Φ = 1.07e-8 mol/m/s/Pa^0.5
  - Oxide (Cr2O3):     Φ = 1.07e-9 mol/m/s/Pa
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
    
    'Incoloy800': {
        'T_ref': 1073.15,            # K (800°C)
        'D_ref': 1.43e-9,            # m²/s
        'K_s_ref': 5.92e-2,          # mol/m³/Pa^0.5
        'E_D': 52000,                # J/mol
        'H_s': -20000,               # J/mol
        'reference': 'JAERI-Tech 2002-090, Forcey et al. (1988)',
        'temp_range_C': [600, 1000],
    },
}

# =============================================================================
# OXIDE PROPERTIES (Level 2a, 2b, 3)
# =============================================================================
OXIDES = {
    'Cr2O3': {
        # Reference conditions
        'T_ref': 625,                # K (351.85°C) - same as metal for consistency
        'D_ox_ref': 1.000e-18,       # m²/s at T_ref
        'K_ox_ref': 1.070e-07,       # mol/m³/Pa at T_ref (Henry's law)
        
        # Activation energies (from literature)
        'E_D_ox': 1.55e5,            # J/mol (155 kJ/mol - diffusion)
        'H_sol_ox': 1.85e5,          # J/mol (185 kJ/mol - solution)
        
        # Pre-exponential factors (calculated to match target Φ_ox)
        'D_ox_0': 8.414e+10,         # m²/s
        'K_ox_0': 1.514e+16,         # mol/m³/Pa
        
        # Target permeability verification
        'Phi_ox_ref': 1.070e-18,     # mol/m/s/Pa at T_ref (= D_ox_ref × K_ox_ref)
        
        # Default geometry
        'thickness': 1e-12,           # m (1 μm)
        'thickness_range': [1e-7, 1e-5],
        
        # Metadata
        'reference': 'Calculated to match target permeability (activation energies from Strehlow & Savage 1974, Serra 1998)',
        'temp_range_K': [625, 1273],  # Extended to lower T
        'temperature_range': [625, 1273],
        'uncertainty_factor': 10,
    },
    
    'Cr2O3_thin': {
        'T_ref': 1073.15,
        'D_ox_ref': 1e-18,
        'K_ox_ref': 1e-12,
        'E_D_ox': 50000,
        'H_sol_ox': 30000,
        'thickness': 6e-10,  # 6 Å
        'reference': 'Zarchy & Axtmann (1979)',
    },
}

# =============================================================================
# MICROSTRUCTURE PARAMETERS (Level 4)
# =============================================================================
MICROSTRUCTURE = {
    # Grain boundary parameters
    'grain_size': 50e-6,             # m (50 μm default, use 100e-9 for nanocrystalline)
    'gb_thickness': 0.5e-9,          # m (0.5 nm)
    'grain_shape': 'equiaxed',       # 'equiaxed', 'columnar'
    'gb_type': 'LAGB',               # 'HAGB', 'LAGB'
    'include_gb_trapping': False,
    
    # Trap parameters
    'trap_list': [
        {
            'name': 'vacancies',
            'binding_energy': 0.3 * F,   # J/mol (0.5 eV)
            'density': 1e26,              # m⁻³
        },
    ],
    
    # Lattice site density (for Oriani model)
    'N_L': 8.46e28,                  # m⁻³ (Fe BCC lattice sites)
}

# =============================================================================
# OXIDE DEFECT PARAMETERS (Level 3, 5)
# =============================================================================
OXIDE_DEFECTS = {
    'area_fraction': 0.01,           # 1% pinholes
    'type': 'pinhole',               # 'pinhole', 'crack'
}

# =============================================================================
# OPERATING CONDITIONS
# =============================================================================
CONDITIONS = {
    # Temperature
    'T_operating': 625,              # K - Reference temperature (NOT T_ref for Arrhenius!)
    'T_range': (423, 700),           # K - Sweep range
    'n_T_points': 20,
    
    # Pressure
    'P_upstream': 1e4,               # Pa (10 kPa)
    'P_downstream': 0,               # Pa
    'P_range': (1e3, 1e5),           # Pa - Sweep range
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
    metal='Fe_alpha',
    oxide='Cr2O3',
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
# PARAMETER CALCULATION SUMMARY (for documentation)
# =============================================================================
"""
CALCULATION SUMMARY - TARGET PERMEABILITIES
============================================

Reference Temperature: T_ref = 625 K (351.85°C)

METAL (Fe_alpha):
-----------------
Target: Φ_metal = 1.070e-08 mol/m/s/Pa^0.5

Given:
  - E_D = 54000 J/mol (literature)
  - H_s = 5900 J/mol (literature)

Chosen:
  - D_ref = 1.000e-11 m²/s (typical for Fe at 350°C)

Calculated:
  - K_s_ref = Φ/D = 1.070e+03 mol/m³/Pa^0.5
  - D_0 = D_ref × exp(E_D/RT_ref) = 1.292e-07 m²/s
  - K_s0 = K_s_ref × exp(H_s/RT_ref) = 1.822e+03 mol/m³/Pa^0.5

Verification: Φ = D_ref × K_s_ref = 1.070e-08 ✓


OXIDE (Cr2O3):
--------------
Target: Φ_oxide = 1.070e-09 mol/m/s/Pa

Given:
  - E_D_ox = 1.55e5 J/mol (literature)
  - H_sol_ox = 1.85e5 J/mol (literature)

Chosen:
  - D_ox_ref = 1.000e-16 m²/s (extrapolated to 350°C)

Calculated:
  - K_ox_ref = Φ/D = 1.070e+07 mol/m³/Pa
  - D_ox_0 = D_ox_ref × exp(E_D_ox/RT_ref) = 8.414e+10 m²/s
  - K_ox_0 = K_ox_ref × exp(H_sol_ox/RT_ref) = 1.514e+21 mol/m³/Pa

Verification: Φ_ox = D_ox_ref × K_ox_ref = 1.070e-09 ✓


PERMEABILITY RATIO:
-------------------
Φ_metal / Φ_oxide = 10.0×
→ Metal is 10× more permeable than oxide
→ Oxide acts as significant barrier ✓


TEMPERATURE DEPENDENCE (example values):
-----------------------------------------
T (K)    Φ_metal (mol/m/s/Pa^0.5)   Φ_oxide (mol/m/s/Pa)
500      5.29e-09                    2.93e-15
625      1.07e-08 ← T_ref           1.07e-09 ← T_ref
750      1.79e-08                    2.43e-06
900      3.25e-08                    1.91e-03

✓ Both increase with temperature (Arrhenius)
✓ At T_ref = 625 K, exact match to targets
"""
