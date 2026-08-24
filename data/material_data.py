from calculations.config.model_config import METALS
MATERIALS = METALS


# #### METAL PROPERTIES - LITERATURE DATA COLLECTION ####
# # All entries filled into clean template format
# # Sources: San Marchi 2007, Forcey 1988, Strehlow & Savage 1974,
# #          Rohrig et al. 1975, Schmidt et al. 1985, Perng & Altstetter,
# #          Brass & Chene, Louthan et al.

# METALSssssssss = {

#     # =========================================================================
#     # SOURCE: San Marchi 2007
#     # doi:10.1016/j.ijhydene.2006.05.008
#     # =========================================================================

#     'Austenitic_SS_1': {
#         'Diffusion parameters': {
#             'T_ref': 700.0,         # K (427°C)
#             'D_ref': 6.16e-11,      # m²/s at T_ref
#             'E_D':   54000,         # J/mol (diffusion activation energy)
#             'D_0':   None,          # m²/s (pre-exponential factor) — not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 6.495e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # mol/m³/Pa^0.5 (pre-exponential factor) — not provided
#             'H_s':    5900,         # J/mol (heat of solution)
#         },
#         'Permeation parameters': {
#             'Phi_ref': 4.004e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # J/mol — not provided (Q_p = E_D + H_s = 59900)
#             'Phi_0':   None,        # mol/m/s/Pa^0.5 — not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'San Marchi 2007, Table 1. doi:10.1016/j.ijhydene.2006.05.008',
#             'temp_range': [423, 700],       # K validity range
#             'pressure_range': [1e5, 3e5],   # Pa
#             'metal_thickness': [],          # m — not provided
#             'gas': 'H₂',
#             'Notes': 'Austenitic stainless steel, calibrated from Table 1.',
#         },
#     },

#     'Austenitic_SS_2': {
#         'Diffusion parameters': {
#             'T_ref': 700.0,         # K (427°C)
#             'D_ref': 5.734e-11,     # m²/s at T_ref
#             'E_D':   53620,         # J/mol
#             'D_0':   None,          # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.104e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # not provided
#             'H_s':    8650,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 6.339e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided (Q_p = E_D + H_s = 62270)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'San Marchi 2007, Table 1. doi:10.1016/j.ijhydene.2006.05.008',
#             'temp_range': [423, 703],   # K
#             'pressure_range': [1e5],    # Pa
#             'metal_thickness': [],
#             'gas': 'H₂',
#             'Notes': 'Austenitic stainless steel, calibrated from Table 1.',
#         },
#     },

#     'Austenitic_SS_3': {
#         'Diffusion parameters': {
#             'T_ref': 623.0,         # K
#             'D_ref': 1.478e-11,     # m²/s at T_ref
#             'E_D':   49300,         # J/mol
#             'D_0':   None,          # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 7.074e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # not provided
#             'H_s':    6860,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.045e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided (Q_p = E_D + H_s = 56160)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'San Marchi 2007, Table 1. doi:10.1016/j.ijhydene.2006.05.008',
#             'temp_range': [373, 623],       # K
#             'pressure_range': [1e2, 3e4],   # Pa
#             'metal_thickness': [],
#             'gas': 'H₂',
#             'Notes': 'Austenitic stainless steel, calibrated from Table 1.',
#         },
#     },

#     'Austenitic_SS_300series': {
#         'Diffusion parameters': {
#             'T_ref': 623.0,         # K
#             'D_ref': 2.69e-11,      # m²/s at T_ref
#             'E_D':   53900,         # J/mol
#             'D_0':   None,          # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 4.322e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # not provided
#             'H_s':    5900,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.163e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided (Q_p = E_D + H_s = 59800)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'San Marchi 2007, Table 3. doi:10.1016/j.ijhydene.2006.05.008',
#             'temp_range': [373, 623],       # K
#             'pressure_range': [1e2, 3e4],   # Pa
#             'metal_thickness': [],
#             'gas': 'H₂',
#             'Notes': 'Austenitic stainless steel 300 series, calibrated from Table 3.',
#         },
#     },

#     'Austenitic_SS_21Cr_6Ni_9Mn_22Cr_13Ni_5Mn': {
#         'Diffusion parameters': {
#             'T_ref': 623.0,         # K
#             'D_ref': 1.633e-11,     # m²/s at T_ref
#             'E_D':   53900,         # J/mol
#             'D_0':   None,          # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 7.107e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # not provided
#             'H_s':    5900,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.161e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided (Q_p = E_D + H_s = 59800)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'San Marchi 2007, Table 3. doi:10.1016/j.ijhydene.2006.05.008',
#             'temp_range': [373, 623],       # K
#             'pressure_range': [1e2, 3e4],   # Pa
#             'metal_thickness': [],
#             'gas': 'H₂',
#             'Notes': 'Austenitic SS grades 21Cr–6Ni–9Mn and 22Cr–13Ni–5Mn, calibrated from Table 3.',
#         },
#     },

#     'Austenitic_SS_A286_JBK75_aged': {
#         'Diffusion parameters': {
#             'T_ref': 623.0,         # K
#             'D_ref': 3.629e-11,     # m²/s at T_ref
#             'E_D':   53900,         # J/mol
#             'D_0':   None,          # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 3.297e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # not provided
#             'H_s':    5900,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.197e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided (Q_p = E_D + H_s = 59800)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'San Marchi 2007, Table 3. doi:10.1016/j.ijhydene.2006.05.008',
#             'temp_range': [373, 623],       # K
#             'pressure_range': [1e2, 3e4],   # Pa
#             'metal_thickness': [],
#             'gas': 'H₂',
#             'Notes': 'Austenitic SS A-286/JBK-75 aged, calibrated from Table 3.',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Forcey et al. 1988
#     # HYDROGEN TRANSPORT AND SOLUBILITY IN 316L AND 1.4914 STEELS
#     # =========================================================================

#     'metal_316L_Heat_treated_ref_cast': {
#         'Diffusion parameters': {
#             'T_ref': 873,           # K
#             'D_ref': 2.194e-10,     # m²/s at T_ref
#             'E_D':   42500,         # J/mol
#             'D_0':   7.66e-08,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 2.723e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   4.640e-1,     # mol/m³/Pa^0.5
#             'H_s':    20585,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.881e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     63085,       # J/mol (= E_D + H_s)
#             'Phi_0':   1.12e-07,    # mol/m/s/Pa^0.5 (= D_0 * K_s0)
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Forcey et al. 1988 — Heat treated reference cast 316L austenitic steel',
#             'temp_range': [523, 873],           # K
#             'pressure_range': [1.33e2, 1e5],    # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Data in SI units. Phi_ref = D_ref * Ks_ref. Ks_0 = Phi_0/D_0 = 1.462.',
#         },
#     },

#     'metal_316L_Commercial': {
#         'Diffusion parameters': {
#             'T_ref': 873,           # K
#             'D_ref': 7.237e-10,     # m²/s at T_ref
#             'E_D':   45500,         # J/mol
#             'D_0':   3.82e-07,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.171e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.500e+0,     # mol/m³/Pa^0.5
#             'H_s':    18510,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.655e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     64030,       # J/mol (= E_D + H_s, note: 18530 used vs 18510 reported — rounding)
#             'Phi_0':   1.80e-07,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Forcey et al. 1988 — Commercial 316L austenitic steel',
#             'temp_range': [523, 873],           # K
#             'pressure_range': [1.33e2, 1e5],    # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Data in SI units. Phi_ref = D_ref * Ks_ref. Ks_0 = Phi_0/D_0 = 0.4712.',
#         },
#     },

#     'metal_1.4914_Heat_treated_ref_cast': {
#         'Diffusion parameters': {
#             'T_ref': 873,           # K
#             'D_ref': 1.118e-8,      # m²/s at T_ref
#             'E_D':   13490,         # J/mol
#             'D_0':   7.17e-08,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 2.179e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.290e+0,     # mol/m³/Pa^0.5
#             'H_s':    29620,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 7.700e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     43100,       # J/mol (= E_D + H_s)
#             'Phi_0':   2.92e-08,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Forcey et al. 1988 — Heat-treated reference cast 1.4914 martensitic steel',
#             'temp_range': [523, 873],           # K
#             'pressure_range': [1.33e2, 1e5],    # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Martensitic steel. Phi_ref = D_ref * Ks_ref. Ks_0 = Phi_0/D_0 = 0.4073.',
#         },
#     },

#     'metal_1.4914_Commercial': {
#         'Diffusion parameters': {
#             'T_ref': 873,           # K
#             'D_ref': 1.047e-8,      # m²/s at T_ref
#             'E_D':   15470,         # J/mol
#             'D_0':   8.82e-08,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 2.903e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.180e+0,     # mol/m³/Pa^0.5
#             'H_s':    26890,        # J/mol (reported); Qs = 27800 — ~910 J/mol rounding in source
#         },
#         'Permeation parameters': {
#             'Phi_ref': 9.737e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     43270,       # J/mol
#             'Phi_0':   3.78e-08,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Forcey et al. 1988 — Commercial 1.4914 martensitic steel',
#             'temp_range': [523, 873],           # K
#             'pressure_range': [1.33e2, 1e5],    # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Martensitic steel. Phi_ref = D_ref * Ks_ref. Ks_0 = Phi_0/D_0 = 0.4286. Qs vs H_s diff = 910 J/mol (rounding in source).',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Strehlow & Savage 1974
#     # DOI: 10.13182/NT74-A16282
#     # Gas: Deuterium (D₂) unless stated
#     # =========================================================================

#     'Nickel_unoxidized_Strehlow1974_tableI': {
#         'Diffusion parameters': {
#             'T_ref': 896,       # K
#             'D_ref': None,      # m²/s — not provided
#             'E_D':   None,      # J/mol — not provided separately
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.872e-9,    # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     53550,       # J/mol (activation energy of permeation)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE I',
#             'temp_range': [873, 973],               # K
#             'pressure_range': [6.67e-2, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Permeation rate at 1 Torr = 1.5e-3 to 1.7e-3 cm³(STP)/(h·cm²) '
#                 '= 1.56e-10 to 1.77e-10 mol/m/s/Pa^0.5. '
#                 'Activation energy range: [12.3, 12.6] kcal/mol = [51.463, 52.718] kJ/mol.'
#             ),
#         },
#     },

#     'Nickel_unoxidized_Strehlow1974_tableII': {
#         'Diffusion parameters': {
#             'T_ref': 896,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.66e-9,     # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     53550,       # J/mol
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [873, 973],               # K
#             'pressure_range': [1.2e-1, 1e5],        # Pa
#             'metal_thickness': [1e-3],              # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.51, 0.55]. '
#                 'Permeation rate at 1 Torr = 1.58e-2 to 1.62e-2 cm³(STP)/(h·cm²) '
#                 '= 1.64e-9 to 1.68e-9 mol/m/s/Pa^0.5.'
#             ),
#         },
#     },

#     'Hastelloy_N_unoxidized_Strehlow1974_tableI': {
#         'Diffusion parameters': {
#             'T_ref': 896,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': None,    # not provided — range only: 1.04e-10 to 1.04e-9 mol/m/s/Pa^0.5
#             'Q_p':     None,    # range: [71128, 75312] J/mol
#             'Phi_0':   None,    # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE I',
#             'temp_range': [873, 973],               # K
#             'pressure_range': [5.33e-1, 4e3],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Permeation rate at 1 Torr = 1e-3 to 1e-2 cm³(STP)/(h·cm²) '
#                 '= 1.04e-10 to 1.04e-9 mol/m/s/Pa^0.5. '
#                 'Activation energy range: [17, 18] kcal/mol = [71.128, 75.312] kJ/mol.'
#             ),
#         },
#     },

#     'Hastelloy_N_unoxidized_Strehlow1974_tableII': {
#         'Diffusion parameters': {
#             'T_ref': 842,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.808e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # range: [71128, 75312] J/mol — not single value
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [873, 973],               # K
#             'pressure_range': [5.33e-1, 4.186e3],   # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.50, 0.54]. '
#                 'Activation energy range: [17, 18] kcal/mol = [71.128, 75.312] kJ/mol.'
#             ),
#         },
#     },

#     'Type304L_SS_unoxidized_Strehlow1974_tableII_970K': {
#         'Diffusion parameters': {
#             'T_ref': 970,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.392e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [1.33e2, 4.186e3],    # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': 'Uncertainty on Phi_ref not available.',
#         },
#     },

#     'Type304L_SS_unoxidized_Strehlow1974_tableII_1058K': {
#         'Diffusion parameters': {
#             'T_ref': 1058,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.392e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [4.0, 4.186e3],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.65, 0.75]. '
#                 'Uncertainty = (-0.4, +0.4) × 1.04e-7.'
#             ),
#         },
#     },

#     'Type304L_SS_oxidized_Strehlow1974_tableII_970K': {
#         'Diffusion parameters': {
#             'T_ref': 970,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.704e-11,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [6.67, 4.186e3],      # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.49, 0.59]. '
#                 'Uncertainty = (-0.5, +0.5) × 1.04e-7.'
#             ),
#         },
#     },

#     'Type406L_SS_oxidized_Strehlow1974_tableII_915K': {
#         'Diffusion parameters': {
#             'T_ref': 915,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.664e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [133.322, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.4, 0.6]. '
#                 'Uncertainty = (-0.3, +0.3) × 1.04e-7.'
#             ),
#         },
#     },

#     'Incoloy800_unoxidized_Strehlow1974_tableII_922K': {
#         'Diffusion parameters': {
#             'T_ref': 922,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.664e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [133.322, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.56, 0.65]. '
#                 'Uncertainty = (-0.2, +0.2) × 1.04e-7.'
#             ),
#         },
#     },

#     'Incoloy800_oxidized_Strehlow1974_tableII_922K': {
#         'Diffusion parameters': {
#             'T_ref': 922,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 4.888e-13,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [1333.22, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.91, 1.10]. '
#                 'Uncertainty = (-0.4, +0.4) × 1.04e-7.'
#             ),
#         },
#     },

#     'Incoloy800_oxidized_Strehlow1974_tableII_804K': {
#         'Diffusion parameters': {
#             'T_ref': 804,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 6.552e-11,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [173.319, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.49, 0.53]. '
#                 'Uncertainty = (-0.4, +0.4) × 1.04e-7.'
#             ),
#         },
#     },

#     'CroloyT9_oxidized_Strehlow1974_tableII_846K': {
#         'Diffusion parameters': {
#             'T_ref': 846,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 6.864e-12,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [933.254, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.63, 0.73]. '
#                 'Uncertainty = (-0.5, +0.5) × 1.04e-7.'
#             ),
#         },
#     },

#     'CroloyT9_oxidized_Strehlow1974_tableII_843K': {
#         'Diffusion parameters': {
#             'T_ref': 843,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 5.928e-13,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [1333.22, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.83, 0.88]. '
#                 'Uncertainty = (-0.4, +0.4) × 1.04e-7.'
#             ),
#         },
#     },

#     'CroloyT22_unoxidized_Strehlow1974_tableII_558K': {
#         'Diffusion parameters': {
#             'T_ref': 558,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.872e-11,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [133.322, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.66, 0.68]. '
#                 'Uncertainty = (-0.1, +0.1) × 1.04e-7.'
#             ),
#         },
#     },

#     'CroloyT22_unoxidized_Strehlow1974_tableII_853K': {
#         'Diffusion parameters': {
#             'T_ref': 853,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 8.32e-10,    # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [133.322, 1e5],       # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.55, 0.57]. '
#                 'Uncertainty = (-0.1, +0.1) × 1.04e-7.'
#             ),
#         },
#     },

#     'CroloyT22_oxidized_Strehlow1974_tableII_853K': {
#         'Diffusion parameters': {
#             'T_ref': 853,       # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 6.24e-10,    # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [133.322],                  # Pa (1 Torr)
#             'reference': 'Strehlow & Savage 1974. DOI: 10.13182/NT74-A16282. TABLE II',
#             'temp_range': [],                       # K — not provided
#             'pressure_range': [4e2, 1e5],           # Pa
#             'metal_thickness': [9e-4, 1.65e-3],     # m
#             'gas': 'D₂',
#             'Notes': (
#                 'Pressure dependence exponent n = [0.57, 0.61]. '
#                 'Uncertainty = (-0.1, +0.1) × 1.04e-7.'
#             ),
#         },
#     },

#     # =========================================================================
#     # SOURCE: Rohrig et al. 1975
#     # Gas: Deuterium (D₂)
#     # =========================================================================

#     'Incoloy800_unoxidized_Rohrig1975_1173K': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.397e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     74056,       # J/mol
#             'Phi_0':   4.4415e-7,   # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [],                   # K — not provided
#             'pressure_range': [1e2, 1e5],       # Pa
#             'metal_thickness': [3.3e-3],        # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     'Incoloy800_oxidized_Rohrig1975_1173K_v1': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.397e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     74056,       # J/mol
#             'Phi_0':   4.4415e-7,   # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [],                   # K — not provided
#             'pressure_range': [1e2, 1e5],       # Pa
#             'metal_thickness': [3.3e-3],        # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     'Incoloy800_oxidized_Rohrig1975_1173K_v2': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 5.499e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     74056,       # J/mol
#             'Phi_0':   1.049e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [873, 1223],          # K
#             'pressure_range': [5e1, 5e4],       # Pa
#             'metal_thickness': [3.3e-3],        # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     'Incoloy600_oxidized_Rohrig1975_1173K': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 6.627e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     64015.2,     # J/mol
#             'Phi_0':   4.5402e-7,   # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [973, 1223],          # K
#             'pressure_range': [1e5, 5e5],       # Pa
#             'metal_thickness': [3.3e-3],        # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     '304L_SS_Rohrig1975_1173K': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.1844e-10,  # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     73220,       # J/mol
#             'Phi_0':   2.0445e-7,   # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [4e3],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [973, 1173],          # K
#             'pressure_range': [],               # not provided
#             'metal_thickness': [3.3e-3],        # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     'X15_CrNiSi_2520_oxidized_Rohrig1975_1173K': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.807e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     64852,       # J/mol
#             'Phi_0':   2.82e-10,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [2e5],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [773, 1223],          # K
#             'pressure_range': [],               # not provided
#             'metal_thickness': [3.3e-3],        # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     'Hastelloy_N_oxidized_Rohrig1975_1173K_v1': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 4.23e-10,    # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     77947.92,    # J/mol
#             'Phi_0':   1.173e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [773, 1223],              # K
#             'pressure_range': [5e-1, 4e3],          # Pa
#             'metal_thickness': [3.3e-3],            # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     'Hastelloy_N_oxidized_Rohrig1975_1173K_v2': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 4.23e-10,    # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     77947.92,    # J/mol
#             'Phi_0':   1.173e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [773, 1223],              # K
#             'pressure_range': [5e-1, 4e3],          # Pa
#             'metal_thickness': [3.3e-3],            # m
#             'gas': 'D₂',
#             'Notes': 'Duplicate entry in source data — identical parameters.',
#         },
#     },

#     'Nickel_unoxidized_Rohrig1975_1173K': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 9.306e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     53555.2,     # J/mol
#             'Phi_0':   2.1714e-7,   # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Rohrig et al. 1975',
#             'temp_range': [773, 1223],          # K
#             'pressure_range': [1e-1, 1e5],      # Pa
#             'metal_thickness': [3.3e-3],        # m
#             'gas': 'D₂',
#             'Notes': '',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Grant et al. 1988
#     # HYDROGEN IN 316 STEEL - DIFFUSION, PERMEATION AND SURFACE REACTION
#     # DOI: 10.1016/0022-3115(88)90128-7
#     # Gas: H₂
#     # Note: All values include the +upper bound of stated uncertainty
#     #       e.g. E_D uncertainty ±0.11 → value = literature + 0.11
#     # =========================================================================

#     'metal_316_steel_Grant1988': {
#         'Diffusion parameters': {
#             'T_ref': 965,           # K
#             'D_ref': 1.069e-9,      # m²/s at T_ref
#             'E_D':   53292.74,      # J/mol (diffusion activation energy, +upper uncertainty bound)
#             'D_0':   8.2e-7,        # m²/s (pre-exponential factor)
#         },
#         'Solubility parameters': {
#             'Ks_ref': 0.1295,       # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.19,         # mol/m³/Pa^0.5 (pre-exponential factor)
#             'H_s':    17791.96,     # J/mol (heat of solution)
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.384e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     68756.78,    # J/mol (permeation activation energy = E_D + H_s)
#             'Phi_0':   8.9e-7,      # mol/m/s/Pa^0.5 (pre-exponential factor)
#         },
#         'Surface kinetics parameters': {
#             'Clean surface': {
#                 'k1_ref': 1.346e-6,     # mol/s/m²/Pa at T_ref
#                 'E_k1':   81560.34,     # J/mol (activation energy for surface reaction)
#                 'k1_0':   3.5e-2,       # mol/s/m²/Pa (pre-exponential factor)
#             },
#             'Surface subjected to partial ion beam cleaning': {
#                 'k1_ref': 6.726e-8,     # mol/s/m²/Pa at T_ref
#                 'E_k1':   61856.16,     # J/mol
#                 'k1_0':   1.5e-4,       # mol/s/m²/Pa
#             },
#             'Foil oxidized on both surfaces': {
#                 'k1_ref': 9.2e-8,       # mol/s/m²/Pa at T_ref
#                 'E_k1':   59860.8,      # J/mol
#                 'k1_0':   1.6e-4,       # mol/s/m²/Pa
#             },
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Grant et al. 1988 — Hydrogen in 316 steel: diffusion, permeation and surface reaction. DOI: 10.1016/0022-3115(88)90128-7',
#             'temp_range': [502, 965],           # K
#             'pressure_range': [1e-2, 1e4],      # Pa
#             'metal_thickness': [5e-5, 1e-4],    # m
#             'gas': 'H₂',
#             'Notes': (
#                 'All parameter values include upper bound of stated uncertainty. '
#                 'k1 = dissociative adsorption / surface recombination rate constant. '
#                 'Three surface conditions measured: clean, partially ion-beam cleaned, and oxidized.'
#             ),
#         },
#     },

#     # =========================================================================
#     # SOURCE: Schmidt et al. 1985
#     # Studies on the permeation of hydrogen and tritium through heat resistant alloys
#     # Gas: T₂ unless stated H₂
#     # =========================================================================

#     'metal_1.4841_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.046e-9,      # m²/s at T_ref
#             'E_D':   35400,         # J/mol
#             'D_0':   3.400e-8,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 9.010e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   7.465e-1,     # mol/m³/Pa^0.5
#             'H_s':    21500,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 9.423e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     56900,       # J/mol
#             'Phi_0':   2.538e-8,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985 — Studies on the permeation of hydrogen and tritium through heat resistant alloys',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4876_Incoloy800_Schmidt1985_v1': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 3.783e-9,      # m²/s at T_ref
#             'E_D':   56400,         # J/mol
#             'D_0':   9.700e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 2.138e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   2.180e-2,     # mol/m³/Pa^0.5
#             'H_s':    200,          # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 8.088e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     56600,       # J/mol
#             'Phi_0':   2.115e-8,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X2_NiCr_30_25_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 3.333e-9,      # m²/s at T_ref
#             'E_D':   69500,         # J/mol
#             'D_0':   3.100e-6,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 6.505e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.774e-1,     # mol/m³/Pa^0.5
#             'H_s':    10200,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.168e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     79700,       # J/mol
#             'Phi_0':   5.499e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X2_NiCrSi_30_25_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 4.067e-9,      # m²/s at T_ref
#             'E_D':   102900,        # J/mol
#             'D_0':   1.010e-4,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.843e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.569e+2,     # mol/m³/Pa^0.5
#             'H_s':    68600,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 7.497e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     171500,      # J/mol
#             'Phi_0':   1.585e-2,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X2_NiCrTi_30_25_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 3.056e-9,      # m²/s at T_ref
#             'E_D':   62300,         # J/mol
#             'D_0':   1.400e-6,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 8.059e-3,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.410e-1,     # mol/m³/Pa^0.5
#             'H_s':    29100,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.463e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     91400,       # J/mol
#             'Phi_0':   1.974e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X10_NiCrTi_30_25_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.152e-9,      # m²/s at T_ref
#             'E_D':   42600,         # J/mol
#             'D_0':   7.600e-8,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 7.315e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   6.679e+0,     # mol/m³/Pa^0.5
#             'H_s':    45900,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 8.424e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     88500,       # J/mol
#             'Phi_0':   5.076e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X40_NiCrTi_30_25_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 6.909e-10,     # m²/s at T_ref
#             'E_D':   23400,         # J/mol
#             'D_0':   6.900e-9,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 9.948e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   2.657e+0,     # mol/m³/Pa^0.5
#             'H_s':    33400,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 6.873e-11,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     56800,       # J/mol
#             'Phi_0':   1.833e-8,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X10_NiCr_30_20_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 2.117e-9,      # m²/s at T_ref
#             'E_D':   48100,         # J/mol
#             'D_0':   2.400e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 7.104e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   3.525e+0,     # mol/m³/Pa^0.5
#             'H_s':    39700,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.504e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     87800,       # J/mol
#             'Phi_0':   8.460e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X10_NiCr_30_22_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.764e-9,      # m²/s at T_ref
#             'E_D':   48600,         # J/mol
#             'D_0':   2.100e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 6.537e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   5.304e+0,     # mol/m³/Pa^0.5
#             'H_s':    44700,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.153e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     93300,       # J/mol
#             'Phi_0':   1.114e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X10_NiCr_30_26_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 2.555e-9,      # m²/s at T_ref
#             'E_D':   112800,        # J/mol
#             'D_0':   1.680e-4,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.502e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   7.008e+2,     # mol/m³/Pa^0.5
#             'H_s':    85900,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.838e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     198700,      # J/mol
#             'Phi_0':   1.177e-1,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X10_NiCrAlTi_32_12_Incoloy800H_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 8.050e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     75200,       # J/mol
#             'Phi_0':   1.311e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_NiCr_22_CoMoAl_Inconel617_Schmidt1985_1023K': {
#         'Diffusion parameters': {
#             'T_ref': 1023,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.349e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     68100,       # J/mol
#             'Phi_0':   7.050e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K (note: source had reversed range — corrected)
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_NiCr_22_CoMoAl_Inconel617_Schmidt1985_1173K': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.973e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     83000,       # J/mol
#             'Phi_0':   1.974e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1173],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Trapping effects observed. Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_NiCr_22_MoFeNb_Inconel625_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1173,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.183e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     6.0e4,       # J/mol
#             'Phi_0':   2.397e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1173],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Trapping effects observed. Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 2.474e-9,      # m²/s at T_ref
#             'E_D':   72200,         # J/mol
#             'D_0':   3.000e-6,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 5.093e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.034e-1,     # mol/m³/Pa^0.5
#             'H_s':    7200,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.260e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     79400,       # J/mol
#             'Phi_0':   3.102e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4762_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1123,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': None,    # not provided
#             'Q_p':     None,    # not provided
#             'Phi_0':   None,    # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [973, 1123],          # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phase transformation in temperature range. No permeation data provided.',
#         },
#     },

#     'metal_1.4763_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1123,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.731e-11,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     141300,      # J/mol
#             'Phi_0':   1.021e-4,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [973, 1123],          # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phase transformation in temperature range. Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_1.4301_Schmidt1985_v1': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.620e-9,      # m²/s at T_ref
#             'E_D':   46700,         # J/mol
#             'D_0':   1.600e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 2.091e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   2.882e+2,     # mol/m³/Pa^0.5
#             'H_s':    73500,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.387e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     120200,      # J/mol
#             'Phi_0':   4.611e-5,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4301_Schmidt1985_v2': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 2.857e-9,      # m²/s at T_ref
#             'E_D':   54200,         # J/mol
#             'D_0':   5.900e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.061e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   4.780e-1,     # mol/m³/Pa^0.5
#             'H_s':    15300,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.032e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     69500,       # J/mol
#             'Phi_0':   2.820e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4401_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 3.229e-9,      # m²/s at T_ref
#             'E_D':   50200,         # J/mol
#             'D_0':   4.500e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 3.966e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   5.327e-2,     # mol/m³/Pa^0.5
#             'H_s':    3000,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.281e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     53200,       # J/mol
#             'Phi_0':   2.397e-8,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4541_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 2.100e-9,      # m²/s at T_ref
#             'E_D':   56800,         # J/mol
#             'D_0':   5.600e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 8.383e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.158e+0,     # mol/m³/Pa^0.5
#             'H_s':    26700,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.760e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     83500,       # J/mol
#             'Phi_0':   6.486e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4571_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.830e-9,      # m²/s at T_ref
#             'E_D':   50000,         # J/mol
#             'D_0':   2.500e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.027e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   9.024e-1,     # mol/m³/Pa^0.5
#             'H_s':    22100,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.879e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     72100,       # J/mol
#             'Phi_0':   2.256e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4550_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.442e-9,      # m²/s at T_ref
#             'E_D':   53200,         # J/mol
#             'D_0':   2.700e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.024e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   3.394e+0,     # mol/m³/Pa^0.5
#             'H_s':    35600,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.477e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     88800,       # J/mol
#             'Phi_0':   9.165e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4580_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.703e-9,      # m²/s at T_ref
#             'E_D':   62300,         # J/mol
#             'D_0':   7.800e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.036e+0,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.645e+0,     # mol/m³/Pa^0.5
#             'H_s':    4700,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.764e-9,    # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     67000,       # J/mol
#             'Phi_0':   1.283e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_X10_NiCrNbAl_30_25_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 9.468e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     63900,       # J/mol
#             'Phi_0':   5.076e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_X10_NiCrNbAlTi_30_25_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.156e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     63900,       # J/mol
#             'Phi_0':   1.692e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_NiCr_22_FeW_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 3.854e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     46700,       # J/mol
#             'Phi_0':   3.807e-8,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_X25_NiCrNb_30_27_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.425e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     70700,       # J/mol
#             'Phi_0':   2.538e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_X25_NiCrNbCe_30_30_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 7.573e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     135400,      # J/mol
#             'Phi_0':   4.597e-4,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Oxidation observed. Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_X25_NiCrCoNb_30_27_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.787e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     78300,       # J/mol
#             'Phi_0':   3.948e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_X25_NiCrSiNbCe_35_27_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 6.278e-10,   # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     62100,       # J/mol
#             'Phi_0':   2.820e-7,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     'metal_Nimonic_PE_13_Schmidt1985': {
#         'Diffusion parameters': {
#             'T_ref': 1223,          # K
#             'D_ref': 1.908e-9,      # m²/s at T_ref
#             'E_D':   33800,         # J/mol
#             'D_0':   5.300e-8,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 6.469e-2,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   3.991e-1,     # mol/m³/Pa^0.5
#             'H_s':    18500,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 1.234e-10,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     52300,       # J/mol
#             'Phi_0':   2.115e-8,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Schmidt et al. 1985',
#             'temp_range': [1023, 1223],         # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     'metal_1.4876_Incoloy800_Schmidt1985_v2': {
#         'Diffusion parameters': {
#             'T_ref': 1223,      # K
#             'D_ref': None,      # not provided
#             'E_D':   None,      # not provided
#             'D_0':   None,      # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,     # not provided
#             'K_s0':   None,     # not provided
#             'H_s':    None,     # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': 2.884e-9,    # mol/m/s/Pa^0.5 at T_ref
#             'Q_p':     59900,       # J/mol
#             'Phi_0':   1.043e-6,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Forcey et al. — H₂ 600–950°C',
#             'temp_range': [873, 1223],          # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref computed from Phi_0 and Q_p via Arrhenius (D_0/Ks_0 unavailable).',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Various authors — D₂, 100–450°C
#     # =========================================================================

#     'metal_various_D2_100_450C': {
#         'Diffusion parameters': {
#             'T_ref': 723,           # K
#             'D_ref': 5.897e-11,     # m²/s at T_ref
#             'E_D':   54000,         # J/mol
#             'D_0':   4.700e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 1.334e-1,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.800e-1,     # mol/m³/Pa^0.5
#             'H_s':    1800,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 7.868e-12,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     55800,       # J/mol
#             'Phi_0':   8.460e-8,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Various authors — D₂, 100–450°C',
#             'temp_range': [373, 723],           # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'D₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Perng & Altstetter — AISI 304, H₂ 100–600°C
#     # =========================================================================

#     'metal_AISI_304_Perng_Altstetter': {
#         'Diffusion parameters': {
#             'T_ref': 873,           # K
#             'D_ref': 1.501e-9,      # m²/s at T_ref
#             'E_D':   54400,         # J/mol
#             'D_0':   2.700e-6,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': 3.200e-4,     # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   1.201e-3,     # mol/m³/Pa^0.5
#             'H_s':    9600,         # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': 4.802e-13,   # mol/m/s/Pa^0.5 at T_ref (= D_ref * Ks_ref)
#             'Q_p':     64000,       # J/mol
#             'Phi_0':   3.243e-9,    # mol/m/s/Pa^0.5
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Perng & Altstetter — AISI 304, H₂ 100–600°C',
#             'temp_range': [373, 873],           # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'H₂',
#             'Notes': 'Phi_ref = D_ref * Ks_ref.',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Brass & Chene — AISI 304/316, grain boundary, T₂ -78–185°C
#     # =========================================================================

#     'metal_AISI_304_316_Brass_Chene_grain_boundary': {
#         'Diffusion parameters': {
#             'T_ref': 458,           # K
#             'D_ref': 1.284e-8,      # m²/s at T_ref
#             'E_D':   42400,         # J/mol
#             'D_0':   8.800e-4,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,         # not provided
#             'K_s0':   None,         # not provided
#             'H_s':    None,         # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': None,        # not provided (Ks_ref unavailable)
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Brass & Chene — AISI 304/316, grain boundary, T₂ -78–185°C',
#             'temp_range': [195, 458],           # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Grain boundary diffusion. Solubility/permeation data not provided.',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Louthan & Derrick — AISI 304/316, T₂ 25–222°C
#     # =========================================================================

#     'metal_AISI_304_316_Louthan_Derrick': {
#         'Diffusion parameters': {
#             'T_ref': 495,           # K
#             'D_ref': 1.178e-12,     # m²/s at T_ref
#             'E_D':   58600,         # J/mol
#             'D_0':   1.800e-6,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,         # not provided
#             'K_s0':   None,         # not provided
#             'H_s':    None,         # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': None,        # not provided (Ks_ref unavailable)
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Louthan & Derrick — AISI 304/316, T₂ 25–222°C',
#             'temp_range': [298, 495],           # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Solubility/permeation data not provided.',
#         },
#     },

#     # =========================================================================
#     # SOURCE: Louthan et al. — AISI 304L, T₂ 22–157°C
#     # =========================================================================

#     'metal_AISI_304L_Louthan': {
#         'Diffusion parameters': {
#             'T_ref': 430,           # K
#             'D_ref': 7.438e-14,     # m²/s at T_ref
#             'E_D':   54000,         # J/mol
#             'D_0':   2.700e-7,      # m²/s
#         },
#         'Solubility parameters': {
#             'Ks_ref': None,         # not provided
#             'K_s0':   None,         # not provided
#             'H_s':    None,         # not provided
#         },
#         'Permeation parameters': {
#             'Phi_ref': None,        # not provided (Ks_ref unavailable)
#             'Q_p':     None,        # not provided
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'Louthan et al. — AISI 304L, T₂ 22–157°C',
#             'temp_range': [295, 430],           # K
#             'pressure_range': [1e5, 4e6],       # Pa
#             'metal_thickness': [5e-3],          # m
#             'gas': 'T₂',
#             'Notes': 'Solubility/permeation data not provided.',
#         },
#     },

#     # =========================================================================
#     # SOURCE: JAERI-Tech 2002-090 / Forcey et al. 1988
#     # =========================================================================

#     'Incoloy800_JAERI_Forcey': {
#         'Diffusion parameters': {
#             'T_ref': 1073.15,       # K (800°C)
#             'D_ref': 1.43e-9,       # m²/s at T_ref
#             'E_D':   52000,         # J/mol
#             'D_0':   None,          # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 5.92e-2,      # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # not provided
#             'H_s':    -20000,       # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': None,        # not computed — Phi_0 not provided
#             'Q_p':     None,        # not provided (Q_p = E_D + H_s = 32000)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'JAERI-Tech 2002-090; Forcey et al. 1988',
#             'temp_range': [873, 1273],          # K (600–1000°C)
#             'pressure_range': [],               # not provided
#             'metal_thickness': [],              # not provided
#             'gas': 'H₂',
#             'Notes': 'Calibrated to experimental permeability data.',
#         },
#     },

#     # =========================================================================
#     # SOURCE: San Marchi 2007
#     # =========================================================================

#     'Fe_alpha_SanMarchi2007': {
#         'Diffusion parameters': {
#             'T_ref': 1073.15,       # K (800°C)
#             'D_ref': 2.5e-8,        # m²/s at T_ref
#             'E_D':   4200,          # J/mol
#             'D_0':   None,          # not provided
#         },
#         'Solubility parameters': {
#             'Ks_ref': 0.12,         # mol/m³/Pa^0.5 at T_ref
#             'K_s0':   None,         # not provided
#             'H_s':    28600,        # J/mol
#         },
#         'Permeation parameters': {
#             'Phi_ref': None,        # not provided
#             'Q_p':     None,        # not provided (Q_p = E_D + H_s = 32800)
#             'Phi_0':   None,        # not provided
#         },
#         'Other common parameters': {
#             'pressure': [],
#             'reference': 'San Marchi 2007',
#             'temp_range': [298, 1173],          # K (25–900°C)
#             'pressure_range': [],               # not provided
#             'metal_thickness': [],              # not provided
#             'gas': 'H₂',
#             'Notes': '',
#         },
#     },

# }

# #### OXIDE LAYER PROPERTIES - LITERATURE DATA COLLECTION ####
# #
# # DATA QUALITY TIERS:
# # -------------------
# # Tier 1: Phi_composite, Phi_bare, oxide_thickness, metal_thickness all known
# #         → standalone Phi_oxide computable via series resistance
# # Tier 2: oxide_thickness known, Phi_bare from METALS dict
# #         → back-calculation possible but baseline-dependent
# # Tier 3: PRF only, thickness unknown → cannot back-calculate
# #
# # BACK-CALCULATION (Series Resistance Model):
# # -------------------------------------------
# #   Phi_oxide = d_oxide / (d_total/Phi_total - d_metal/Phi_metal)
# #
# #   where d_total = d_oxide + d_metal
# #   Phi_metal from METALS dict at same T:
# #     Phi_metal(T) = Phi_ref * exp((-Q_p/R)*(1/T - 1/T_ref))
# #
# #   For Arrhenius parameters of standalone oxide:
# #     Back-calculate Phi_oxide at multiple T, then fit:
# #     Phi_oxide(T) = Phi_oxide_0 * exp(-Q_p_oxide / RT)
# #
# #   PRF = Phi_bare / Phi_composite  (dimensionless, >= 1)
# # =============================================================================

