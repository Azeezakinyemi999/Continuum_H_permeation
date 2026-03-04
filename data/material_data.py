from calculations.config.model_config import METALS

MATERIALS = METALS

# # Start with ONE material (Incoloy 800):
# # •	Find D₀ and E_D from literature
# # •	Find K_s₀ and ΔH_s from literature
# # •	Document your source
# # •	Include temperature range of validity


# """
# Material properties for hydrogen permeation calculations.

# References:
# - Forcey et al. (1988): J. Nucl. Mater. 160, 117-124
# - JAERI-Tech 2002-090: Hydrogen permeability data compilation
# - San Marchi et al. (2007): Technical Reference for Hydrogen Compatibility

# IMPORTANT: Parameters have been calibrated to match JAERI experimental data.
# The original Forcey1988 values gave permeability ~10^6 too low.
# """

# # import numpy as np

# # # Gas constant
# # R = 8.314  # J/mol/K

# # MATERIALS = {
# #     # =========================================================================
# #     # CALIBRATED to JAERI-Tech 2002-090 experimental data
# #     # Arrhenius fit: P = 1.1521e-04 × exp(-32076 / RT)
# #     # R² = 0.9993, MAPE = 1.0%
# #     # =========================================================================
# #     'Incoloy800': {
# #         # Diffusivity: D = D_0 * exp(-E_D / RT)
# #         'D_0': 5.00e-07,    # m²/s
# #         'E_D': 52000,       # J/mol (literature value)
        
# #         # Solubility: K_s = K_s0 * exp(-H_s / RT)
# #         # H_s < 0 means exothermic (K_s decreases with T)
# #         'K_s0': 1.0e-09,   # mol/m³/Pa^0.5
# #         'H_s': -20000,      # J/mol
        
# #         # Derived: P = D * K_s = P_0 * exp(-E_p / RT)
# #         # P_0 = D_0 * K_s0 = 1.15e-04 mol/m/s/Pa^0.5
# #         # E_p = E_D + H_s = 52000 + (-20000) = 32000 J/mol
        
# #         # Metadata
# #         'reference': 'Calibrated to JAERI-Tech 2002-090 Fig 2.2',
# #         'temp_range': [600, 1000],  # °C
# #         'notes': 'Arrhenius fit R²=0.9993, MAPE=1.0%'
# #     },
    
# #     # =========================================================================
# #     # Pure Iron (for validation/comparison)
# #     # =========================================================================
# #     'Fe_alpha': {
# #         'D_0': 4.0e-8,      # m²/s
# #         'E_D': 4200,        # J/mol
# #         'K_s0': 1.9e-1,     # mol/m³/Pa^0.5
# #         'H_s': 28600,       # J/mol (positive - endothermic)
# #         'reference': 'San Marchi 2007',
# #         'temp_range': [25, 900],
# #         'notes': 'Alpha (BCC) iron'
# #     }
# # }

# """
# Material properties for hydrogen permeation calculations.

# Uses REFERENCE-TEMPERATURE Arrhenius format:
#     k(T) = k_ref × exp((-E/R) × (1/T - 1/T_ref))

# where k_ref is the MEASURED property value at T_ref.
# """

# import numpy as np

# R = 8.314  # J/mol/K

MATERIALSsss = {
    'Incoloy800': {
        # Reference temperature
        'T_ref': 1073.15,       # K (800°C)
        
        # Reference values at T_ref (from literature/experiment)
        'D_ref': 1.43e-9,       # m²/s - Diffusivity at 800°C
        'K_s_ref': 5.92e-2,     # mol/m³/Pa^0.5 - Solubility at 800°C
        
        # Activation energies
        'E_D': 52000,           # J/mol (diffusion activation energy)
        'H_s': -20000,          # J/mol (heat of solution)
        
        # Metadata
        'reference': 'JAERI-Tech 2002-090, Forcey et al. (1988)',
        'temp_range': [600, 1000],  # °C validity range
        'notes': 'Calibrated to experimental permeability data'
    },
    
    'Fe_alpha': {
        'T_ref': 1073.15,
        'D_ref': 2.5e-8,        # m²/s at 800°C
        'K_s_ref': 0.12,        # mol/m³/Pa^0.5 at 800°C
        'E_D': 4200,
        'H_s': 28600,
        'reference': 'San Marchi 2007',
        'temp_range': [25, 900],
    }
}