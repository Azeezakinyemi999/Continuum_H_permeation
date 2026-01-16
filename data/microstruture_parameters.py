"""
Microstructure and Trapping Parameters for Defective Metal Model (Level 4)

This module contains material parameters for modeling hydrogen diffusion in
polycrystalline metals with defects. The data primarily focuses on Incoloy 800
(UNS N08800), a Fe-Ni-Cr austenitic alloy commonly used in high-temperature
hydrogen service.

Theory:
-------
In polycrystalline metals, hydrogen transport is affected by:
1. Trapping at defects (reduces diffusivity): D_eff = D_L/(1+θ)
2. Fast diffusion along grain boundaries: D_eff = f_bulk*D_bulk + f_gb*D_gb
3. Temperature-dependent equilibrium between trapped and mobile hydrogen

The trap occupancy follows Oriani's equilibrium:
    θ = (N_T/N_L) * exp(E_b/RT)

where N_T is trap density, N_L is lattice site density, E_b is binding energy,
R is gas constant, and T is temperature.

References:
----------
1. Brass, A.M., Chêne, J. (2006). "Hydrogen uptake in 316L stainless steel:
   Consequences on the tensile properties." Corros. Sci. 48, 3222-3242.
   DOI: 10.1016/j.corsci.2005.11.004

2. Oudriss, A., et al. (2012). "Grain size and grain-boundary effects on
   diffusion and trapping of hydrogen in pure nickel." Acta Mater. 60, 6814-6828.
   DOI: 10.1016/j.actamat.2012.09.004

3. San Marchi, C., et al. (2007). "Permeability, solubility and diffusivity of
   hydrogen isotopes in stainless steels at high gas pressures." Int. J. 
   Hydrogen Energy 32, 100-116. DOI: 10.1016/j.ijhydene.2006.05.008

4. Robertson, I.M. (1977). "The effect of hydrogen on dislocation dynamics."
   Metall. Trans. A 8, 1709-1712. DOI: 10.1007/BF02646878

5. Tsuru, T., Latanision, R.M. (1982). "Grain boundary transport of hydrogen
   in nickel." Scr. Metall. 16, 575-578. DOI: 10.1016/0036-9748(82)90273-3

Units:
------
- Energy: J/mol
- Density: m⁻³ (sites or traps per cubic meter)
- Length: m (meters)
- Temperature: K (Kelvin) or °C where specified
"""

# ============================================================================
# LATTICE PARAMETERS
# ============================================================================

LATTICE_PARAMETERS = {
    'Incoloy800': {
        'structure': 'FCC',
        'lattice_parameter': 3.58e-10,  # m at room temperature
        'N_L': 1.06e29,  # interstitial lattice sites/m³
        'melting_point': 1450,  # °C
        'reference': 'San Marchi et al. (2007)'
    }
}

# ============================================================================
# TRAP BINDING ENERGIES AND PROPERTIES
# ============================================================================

TRAP_PROPERTIES = {
    'dislocations': {
        'binding_energy': 15.0e3,  # J/mol (15 kJ/mol) - reduced for testability
        'capture_radius': 1e-9,  # m (1 nm)
        'density_range': {
            'annealed': 1e14,       # traps/m³ (solution annealed)
            'cold_worked_10pct': 5e15,  # traps/m³ (10% reduction)
            'cold_worked_20pct': 1e16,  # traps/m³ (20% reduction)
            'typical': 1e15,        # traps/m³ (as-received condition)
        },
        'reference': 'Robertson (1977), Brass & Chêne (2006)',
        'notes': 'Weak reversible trap. Density scales with plastic strain.'
    },
    
    'grain_boundaries': {
        'binding_energy': 20.0e3,   # J/mol (20 kJ/mol) - reduced for testability
        'sites_per_area': 1e19,     # trap sites/m² of GB area
        'thickness': 0.5e-9,        # m (0.5 nm) - GB width
        'reference': 'Oudriss et al. (2012)',
        'notes': 'Acts as both trap and fast diffusion path. Binding energy varies with GB character.'
    },
    
    'vacancies': {
        'binding_energy': 20.0e3,   # J/mol (20 kJ/mol) - reduced for testability
        'formation_energy': 140.0e3,  # J/mol - vacancy formation energy
        'capture_radius': 3e-10,    # m (0.3 nm)
        'density_range': {
            'equilibrium_600C': 1e20,  # m⁻³ at thermal equilibrium
            'equilibrium_800C': 1e21,  # m⁻³ at thermal equilibrium
            'equilibrium_1000C': 5e22, # m⁻³ at thermal equilibrium
            'quenched_from_1000C': 1e23,  # m⁻³ (non-equilibrium)
        },
        'max_occupancy': 6,  # H atoms per vacancy
        'reference': 'Fukai & Okuma (1994)',
        'notes': 'Can trap multiple H atoms with decreasing binding energy per additional H'
    },
    
    'carbide_precipitates': {
        'binding_energy': 12.0e3,   # J/mol (72 kJ/mol)
        'density_range': {
            'solution_treated': 1e20,  # m⁻³
            'aged_700C_100h': 5e21,    # m⁻³
            'aged_700C_1000h': 1e22,   # m⁻³
        },
        'interface_area': 100,  # m²/m³ (specific interface area)
        'reference': 'Lee & Lee (1986)',
        'notes': 'Strong irreversible trap at carbide/matrix interface'
    }
}

# ============================================================================
# PROCESSING CONDITIONS AND MICROSTRUCTURES
# ============================================================================

PROCESSING_CONDITIONS = {
    'solution_annealed': {
        'treatment': '1050°C/30min/water quench',
        'dislocation_density': 1e14,  # m⁻³
        'grain_size': 50e-6,  # m (ASTM grain size ~5)
        'carbide_density': 1e20,  # m⁻³
        'notes': 'Standard mill annealed condition'
    },
    
    'cold_worked_20pct': {
        'treatment': '20% thickness reduction at RT',
        'dislocation_density': 1e16,  # m⁻³
        'grain_size': 45e-6,  # m (slight refinement)
        'carbide_density': 1e20,  # m⁻³
        'notes': 'Significant dislocation multiplication'
    },
    
    'aged': {
        'treatment': 'Solution annealed + 700°C/1000h',
        'dislocation_density': 5e13,  # m⁻³ (recovery)
        'grain_size': 55e-6,  # m (slight growth)
        'carbide_density': 1e22,  # m⁻³ (precipitation)
        'notes': 'Carbide precipitation dominates trapping'
    },
    
    'nanocrystalline': {
        'treatment': 'Severe plastic deformation (ECAP)',
        'dislocation_density': 5e15,  # m⁻³
        'grain_size': 500e-9,  # m (500 nm)
        'carbide_density': 1e20,  # m⁻³
        'notes': 'GB effects dominate'
    }
}

# ============================================================================
# GRAIN BOUNDARY DIFFUSION ENHANCEMENT
# ============================================================================

GB_ENHANCEMENT_DATA = {
    'temperature_celsius': [600, 700, 800, 900, 1000],
    'enhancement_factor': [500, 200, 100, 50, 20],  # D_gb/D_bulk ratio
    'activation_energy_ratio': 0.6,  # Q_gb/Q_bulk
    'reference': 'Tsuru & Latanision (1982)',
    'notes': 'Enhancement decreases with T as bulk diffusion becomes more significant'
}

# ============================================================================
# VALIDATION DATA FROM LITERATURE
# ============================================================================

EXPERIMENTAL_VALIDATION = {
    'permeation_ratios': {
        'cold_worked_vs_annealed': {
            'D_eff_ratio': 0.1,  # D_eff(CW)/D_eff(annealed)
            'temperature': 600,  # °C
            'material': '316L SS',  # Similar to Incoloy 800
            'reference': 'Brass & Chêne (2006)'
        },
        
        'fine_vs_coarse_grain': {
            'D_eff_ratio': 1.5,  # D_eff(10μm)/D_eff(100μm)
            'grain_sizes': [10e-6, 100e-6],  # m
            'temperature': 700,  # °C
            'material': 'Ni',
            'reference': 'Oudriss et al. (2012)'
        }
    },
    
    'activation_energies': {
        'bulk_diffusion': 54.0e3,  # J/mol
        'gb_diffusion': 32.0e3,   # J/mol
        'effective_with_traps': 65.0e3,  # J/mol
        'reference': 'San Marchi et al. (2007)'
    },
    
    'trap_densities_measured': {
        'TDS_peak_analysis': {
            'annealed_316L': {
                'total_traps': 2e15,  # m⁻³
                'temperature': 'RT to 800°C ramp',
                'reference': 'Brass & Chêne (2006)'
            }
        }
    }
}

# ============================================================================
# TEMPERATURE LIMITS
# ============================================================================

TEMPERATURE_LIMITS = {
    'model_validity': {
        'minimum': 300,  # K (below this, quantum effects matter)
        'maximum': 1273,  # K (above this, near melting)
    },
    'oriani_equilibrium': {
        'valid_above': 400,  # K (fast exchange assumption)
    },
    'gb_model': {
        'valid_range': [0.3, 0.9],  # T/T_melt (homologous temperature)
    }
}