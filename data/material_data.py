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


#### LITERATURE DATA COLLECTIONS ####


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
    },

        'Austenitic_SS_1': {
        'T_ref': 700.0,             # K (427°C)
        'D_ref': 6.16e-11,           # m²/s at 427°C
        'K_s_ref': 6.495e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 4.004e-12,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 54000,              # J/mol (diffusion activation energy)
        'H_s': 5900,               # J/mol (heat of solution)
        'reference': 'San Marchi 2007',
        'temp_range': [423, 700], # K validity range
        'pressure_range': [1e5, 3e5],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from Table 1 of San Marchi 2007 data. doi:10.1016/j.ijhydene.2006.05.008'
    },


        'Austenitic_SS_2': {
        'T_ref': 700.0,             # K (427°C)
        'D_ref': 5.734e-11,           # m²/s at 427°C
        'K_s_ref': 1.104e-1,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 6.339e-12,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 53620,              # J/mol (diffusion activation energy)
        'H_s': 8650,               # J/mol (heat of solution)
        'reference': 'San Marchi 2007',
        'temp_range': [423, 703], # K validity range
        'pressure_range': [1e5],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from Table 1 of San Marchi 2007 data. doi:10.1016/j.ijhydene.2006.05.008'
    },

        'Austenitic_SS_3': {
        'T_ref': 623.0,             # K (427°C)
        'D_ref': 1.478e-11,           # m²/s at 427°C
        'K_s_ref': 7.074e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 1.045e-12,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 49300,              # J/mol (diffusion activation energy)
        'H_s': 6860,               # J/mol (heat of solution)
        'reference': 'San Marchi 2007',
        'temp_range': [373, 623], # K validity range
        'pressure_range': [1e2, 3e4],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from Table 1 of San Marchi 2007 data. doi:10.1016/j.ijhydene.2006.05.008'
    },
    'Austenitic_SS_300series': {
        'T_ref': 623.0,             # K (427°C)
        'D_ref': 2.69e-11,           # m²/s at 427°C
        'K_s_ref': 4.322e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 1.163e-12,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 53900,              # J/mol (diffusion activation energy)
        'H_s': 5900,               # J/mol (heat of solution)
        'reference': 'San Marchi 2007',
        'temp_range': [373, 623], # K validity range
        'pressure_range': [1e2, 3e4 ],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from Table 3 of San Marchi 2007 data. doi:10.1016/j.ijhydene.2006.05.008'
    },
    
    'Austenitic_SS_21Cr–6Ni–9Mn, 22Cr–13Ni–5Mn ': {
        'T_ref': 623.0,             # K (427°C)
        'D_ref': 1.633e-11,           # m²/s at 427°C
        'K_s_ref': 7.107e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 1.161e-12,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 53900,              # J/mol (diffusion activation energy)
        'H_s': 5900,               # J/mol (heat of solution)
        'reference': 'San Marchi 2007',
        'temp_range': [373, 623], # K validity range
        'pressure_range': [1e2, 3e4 ],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from Table 3 of San Marchi 2007 data. doi:10.1016/j.ijhydene.2006.05.008'
    },

    'Austenitic_SS_A-286/JBK-75 aged ': {
        'T_ref': 623.0,             # K (427°C)
        'D_ref': 3.629e-11,           # m²/s at 427°C
        'K_s_ref': 3.297e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 1.197e-12,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 53900,              # J/mol (diffusion activation energy)
        'H_s': 5900,               # J/mol (heat of solution)
        'reference': 'San Marchi 2007',
        'temp_range': [373, 623], # K validity range
        'pressure_range': [1e2, 3e4 ],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from Table 3 of San Marchi 2007 data. doi:10.1016/j.ijhydene.2006.05.008'
    },

    'Heat_treated_316L steel': {
        'T_ref': 873.0,             # K (427°C)
        'D_ref': 2.194e-10,           # m²/s at 427°C
        'K_s_ref': 2.723e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 5.975e-12,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 42500,              # J/mol (diffusion activation energy)
        'H_s': 20580,               # J/mol (heat of solution)
        'reference': 'Forcey 1988',
        'temp_range': [523, 873], # K validity range
        'pressure_range': [1.33e2, 1e5],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from 5. Experimental results of Forcey 1988 data. doi:'
    },

    'Commercial_316L steel': {
        'T_ref': 873.0,             # K (427°C)
        'D_ref': 7.237e-10,           # m²/s at 427°C
        'K_s_ref': 1.171e-1,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 8.474e-11,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 45500,              # J/mol (diffusion activation energy)
        'H_s': 18510,               # J/mol (heat of solution)
        'reference': 'Forcey 1988',
        'temp_range': [523, 873], # K validity range
        'pressure_range': [1.33e2, 1e5],  # Pa
        'Notes': 'Austenitic stainless steel, calibrated from 5. Experimental results of Forcey 1988 data. doi:'
    },
    'Heat_treated_1.4914 steel': {
        'T_ref': 873.0,             # K (427°C)
        'D_ref': 1.118e-8,           # m²/s at 427°C
        'K_s_ref': 2.179e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 2.436e-10,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 13490,              # J/mol (diffusion activation energy)
        'H_s': 29620,               # J/mol (heat of solution)
        'reference': 'Forcey 1988',
        'temp_range': [523, 873], # K validity range
        'pressure_range': [1.33e2, 1e5],  # Pa
        'Notes': 'Martensitic stainless steel, calibrated from 5. Experimental results of Forcey 1988 data. doi:'
    },
    'Commercial_1.4914 steel': {
        'T_ref': 873.0,             # K (427°C)
        'D_ref': 1.047e-8,           # m²/s at 427°C
        'K_s_ref': 2.903e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 3.039e-10,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 15470,              # J/mol (diffusion activation energy)
        'H_s': 26890,               # J/mol (heat of solution)
        'reference': 'Forcey 1988',
        'temp_range': [523, 873], # K validity range
        'pressure_range': [1.33e2, 1e5],  # Pa
        'Notes': 'Martensitic stainless steel, calibrated from 5. Experimental results of Forcey 1988 data. doi:'
    },

    'Commercial_1.4914 steel': {
        'T_ref': 873.0,             # K (427°C)
        'D_ref': 1.047e-8,           # m²/s at 427°C
        'K_s_ref': 2.903e-2,        # mol/m³/Pa^0.5 at 427°C
        'Phi_ref': 3.039e-10,        # mol/m/s/Pa^0.5 at 427°C
        'E_D': 15470,              # J/mol (diffusion activation energy)
        'H_s': 26890,               # J/mol (heat of solution)
        'reference': 'Forcey 1988',
        'temp_range': [523, 873], # K validity range
        'pressure_range': [1.33e2, 1e5],  # Pa
        'Notes': 'Martensitic stainless steel, calibrated from 5. Experimental results of Forcey 1988 data. doi:'
    },

    'Name': {
        'T_ref':             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':        # mol/m/s/Pa^0.5 at Tref
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'reference': '',
        'temp_range': [], # K validity range
        'pressure_range': [],  # Pa
        'Notes': ''
    },

    'Nickel_unoxidized': {
        'T_ref': 896,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':1.872e-9 ,  # mol/m/s/Pa^0.5 at Tref
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 53550         #J/mol usually E_D + H_s
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [873, 973], # K validity range
        'pressure_range': [6.67e-2,1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282'
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Permeation rate = 1.5e-3 to 1.7e-3 cm3 (STP)/(h cm2) =  1.56e-10 to 1.77e-10  mol/m/s/Pa^0.5
        At pressure = 1 Torr = 133.322pa
        temperature range = [873, 973] kelvin
        Activation energy ranges= [12.3, 12.6 kcal/mol], 1kcal = 4.184kj, [51.463, 52.718]KJ/mol

        '''
    },
    'Nickel_unoxidized': {
        'T_ref': 896,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':1.66e-9,  # mol/m/s/Pa^0.5 at Tref
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 53550         #J/mol usually E_D + H_s
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [873, 973], # K validity range
        'pressure_range': [1.2e-1, 1e5],  # Pa
        'metal_thickness':[1e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II'
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.51,0.55] #n in P^n
        Permeation rate = 1.58e-2 to 1.62e-2 cm3 (STP)/(h cm2) =  1.64e-9 to 1.68e-9  mol/m/s/Pa^0.5
        At pressure = 1 Torr = 133.322pa
        '''
    },
    'Hastelloy_N_unoxidized': {
        'T_ref': 896,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':        # mol/m/s/Pa^0.5 at Tref
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [5.33e-1,4e3],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282'
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Permeation rate = 1e-3 to 1e-2 cm3 (STP)/(h cm2) =  1.04e-10 to 1.04e-9  mol/m/s/Pa^0.5
        At pressure = 1 Torr = 133.322pa
        temperature range = [873, 973] kelvin
        Activation enerfy ranges= [17, 18 kcal/mol], 1kcal = 4.184kj, [71.128,75.312]KJ/mol
        '''
    },

    'Hastelloy_N_unoxidized': {
        'T_ref': 842,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':2.808e-10 ,        # mol/m/s/Pa^0.5 at Tref
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [5.33e-1,4.186e3],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II'
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.5,0.54] #n in P^n
        Permeation rate = 1e-3 to 1e-2 cm3 (STP)/(h cm2) =  1.04e-10 to 1.04e-9  mol/m/s/Pa^0.5
        At pressure = 1 Torr = 133.322pa
        temperature range = [873, 973] kelvin
        Activation enerfy ranges= [17, 18 kcal/mol], 1kcal = 4.184kj, [71.128,75.312]KJ/mol
        '''
    },
    'Type304L_SS_unoxidized': {
        'T_ref': 970,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':2.392e-10 ,        # mol/m/s/Pa^0.5 at Tref #uncertainty = not available
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [1.33e2, 4.186e3],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974

    },
    'Type304L_SS_unoxidized': {
        'T_ref': 1058,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':2.392e-10 ,        # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.4, +0.4)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [4.0, 4.186e3],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.65,0.75] #n in P^n
        '''

    },

    'Type304L_SS_oxidized': {
        'T_ref': 970,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':2.704e-11 ,        # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.5, +0.5)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [6.67, 4.186e3],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.49,0.59] #n in P^n
        '''

    },


    'Type406L_SS_oxidized': {
        'T_ref': 915,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':1.664e-12, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.3, +0.3)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [133.322, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.4,0.6] #n in P^n
        '''

    },

    'Incoloy800_unoxidized': {
        'T_ref': 922,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':1.664e-10, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.2, +0.2)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [133.322, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.56,0.65] #n in P^n
        '''

    },

    'Incoloy800_oxidized': {
        'T_ref': 922,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':4.888e-13, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.4, +0.4)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [1333.22, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.91,1.10] #n in P^n
        '''

    },
    'Incoloy800_oxidized': {
        'T_ref': 804,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':6.552e-11, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.4, +0.4)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [173.319, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.49,0.53] #n in P^n
        '''

    },

    'CroloyT9_oxidized': {
        'T_ref': 846,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':6.864e-12, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.5, +0.5)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [933.254, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.63,0.73] #n in P^n
        '''

    },

    'CroloyT9_oxidized': {
        'T_ref': 843,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':5.928e-13, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.4, +0.4)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [1333.22, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.83,0.88] #n in P^n
        '''

    },
    'CroloyT22_unoxidized': {
        'T_ref': 558,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':1.872e-11, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.1, +0.1)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [133.322, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.66,0.68] #n in P^n
        '''

    },

    'CroloyT22_unoxidized': {
        'T_ref': 853,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':8.32e-10, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.1, +0.1)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [133.322, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.55,0.57] #n in P^n
        '''

    },

        'CroloyT22_oxidized': {
        'T_ref': 853,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':6.24e-10, # mol/m/s/Pa^0.5 at Tref #uncertainty = (-0.1, +0.1)times1.04e-7
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 
        'pressure': [133.322],  # Pa
        'reference': 'Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282',
        'temp_range': [], # K validity range
        'pressure_range': [4e2, 1e5],  # Pa
        'metal_thickness':[9e-4,1.65e-3] # m
        'Notes': 'This result is for Deuterium by Strehlow and Savage 1974 DOI: 10.13182/NT74-A16282. TABLE II '
        #Expected Permeation rate according to Strehlow and Savage 1974
        '''
        Pressure dependence = [0.57,0.61] #n in P^n
        '''

    },

    'Incoloy800_unoxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':2.397e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 74056 #J/mol
        'pressure': [],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [], # K validity range
        'pressure_range': [1e2, 1e5],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 4.4415e-7 # mol/m/s/Pa^0.5
        '''

    },

    'Incoloy800_oxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':2.397e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 74056 #J/mol
        'pressure': [],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [], # K validity range
        'pressure_range': [1e2, 1e5],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 4.4415e-7 # mol/m/s/Pa^0.5
        '''

    },

    'Incoloy800_oxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':5.499e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 74056 #J/mol
        'pressure': [],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [873, 1223], # K validity range
        'pressure_range': [5e1,5e4],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 1.049e-6 # mol/m/s/Pa^0.5
        '''

    },

    'Incoloy600_oxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':6.627e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 64015.2 #J/mol
        'pressure': [],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [973, 1223], # K validity range
        'pressure_range': [1e5,5e5],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 4.5402e-7 # mol/m/s/Pa^0.5
        '''

    },

    '304L_SS': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':1.1844e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 73220 #J/mol
        'pressure': [4e3],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [973, 1173], # K validity range
        'pressure_range': [],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 2.0445e-7 # mol/m/s/Pa^0.5
        '''

    },
        'X 15 Cr-Ni-Si 2520_oxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':3.807e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 64852 #J/mol
        'pressure': [2e5],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [773, 1223], # K validity range
        'pressure_range': [],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 2.82e-10 # mol/m/s/Pa^0.5
        '''

    },

    'Hastelloy_N_oxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':4.23e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 77947.92 #J/mol
        'pressure': [],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [773, 1223], # K validity range
        'pressure_range': [5e-1, 4e3],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 1.173e-6 # mol/m/s/Pa^0.5
        '''
        

    },

    'Hastelloy_N_oxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':4.23e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 77947.92 #J/mol
        'pressure': [],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [773, 1223], # K validity range
        'pressure_range': [5e-1, 4e3],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 1.173e-6 # mol/m/s/Pa^0.5
        '''    
    },


    'Nickel_unoxidized': {
        'T_ref': 1173,             # K 
        'D_ref':           # m²/s at Tref
        'K_s_ref':        # mol/m³/Pa^0.5 at Tref
        'Phi_ref':9.306e-10, # mol/m/s/Pa^0.5 at Tref 
        'E_D':              # J/mol (diffusion activation energy)
        'H_s':               # J/mol (heat of solution)
        'Activation energy of permeation': 53555.2#J/mol
        'pressure': [],  # Pa
        'reference': ' Rohrig et al 1975 ',
        'temp_range': [773, 1223], # K validity range
        'pressure_range': [1e-1, 1e5],  # Pa
        'metal_thickness':[3.3e-3] # m
        'Notes': 'This result is for Deuterium by Rohrig et al 1975 '
        '''
        phi_0 = 2.1714e-7 # mol/m/s/Pa^0.5
        '''    
    },

    #### Grant et al 1988: HYDROGEN IN 316 STEEL - DIFFUSION, PERMEATION AND SURFACE REACTION DOI:
    '316_steel': {
        'Diffusion parameters':{
            'T_ref': 965 # K
            'D_ref': 1.069e-9          # m²/s at Tref
            'E_D': 53292.74              # J/mol (diffusion activation energy)
            'D_o':8.2e-7             # m²/s pre-exponential factor for diffusion
        },
        'Solubility parameters':{
            'K_s_ref': 0.1295,        # mol/m³/Pa^0.5 at Tref
            'K_s0': 1.19,        # mol/m³/Pa^0.5
            'H_s': 17791.96,       # J/mol (heat of solution)
        },
        'Permeation parameters':{
            'Phi_ref': 1.384e-10, # mol/m/s/Pa^0.5 at Tref 
            'Activation energy of permeation': 68756.78 #J/mol
            'Phi_0': 8.9e-7 # mol/m/s/Pa^0.5 pre-exponential factor for permeation
        },
        'surface parameters':{
            'Clean surface':{
                'k1_ref': 1.346e-06, #  mol/s/m^-2/Pa^-1 at Tref
                'E_k1':  81560.34, # J/mol activation energy for surface reaction
                'k1_0': 3.5e-2 #  at Tref

            },
            'Surface subjected to partial ion beamed cleaning':{
                'k1_ref': 6.726e-08, #  mol/s/m^-2/Pa^-1 at Tref
                'E_k1':  61856.16, # J/mol activation energy for surface reaction
                'k1_0': 1.5e-4 #  at Tref
            },
            'Foil oxidized on both surface':{
                'k1_ref': 9.2e-08, #  mol/s/m^-2/Pa^-1 at Tref
                'E_k1':  59860.8, # J/mol activation energy for surface reaction
                'k1_0': 1.6e-4 #  at Tref
            },
    },
    'Other common parameters':{
        'pressure': [],  # Pa
        'reference': 'Grant et al 1988: HYDROGEN IN 316 STEEL - DIFFUSION, PERMEATION AND SURFACE REACTION DOI: 10.1016/0022-3115(88)90128-7',
        'temp_range': [502, 965], # K validity range
        'pressure_range': [1e-2, 1e4],  # Pa
        '_thickness_range':[5e-5,1e-4] # m
        '''
        'Note': 'This result is for Hydrogen by Grant et al 1988: HYDROGEN IN 316 STEEL - DIFFUSION, PERMEATION AND SURFACE REACTION DOI: 10.1016/0022-3115(88)90128-7'
        All these results have an uncertainty attached to it but I added the highest bound of the assigned uncertainty e.g  for E_D, the uncertainty is 9(-0.11 to +0.11), 
        so, i added the +0.11 to the literature value 

        '''
    }
    }
}