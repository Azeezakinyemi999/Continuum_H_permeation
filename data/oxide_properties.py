from calculations.config.model_config import OXIDES

OXIDE_PROPERTIES = OXIDES
# """
# Oxide properties for hydrogen permeation calculations.

# References:
# - Strehlow & Savage (1974): Nuclear Technology, 22:127-137
# - Zarchy & Axtmann (1979): J. Nuclear Materials, 79:110-117
# - Serra et al. (1998): J. Nuclear Materials, 258-263:1028-1032
# IMPORTANT: The oxide properties have been adjusted to yield realistic
# Permeation Reduction Factors (PRF) in the range of 10 to 3800, as suggested
# by literature for dense oxide barriers.
# The previous parameters led to unrealistically high PRF (~10^22).
# """

# OXIDE_PROPERTIES = {
#     'Cr2O3': {
#         # ===================================================================
#         # CORRECTED VALUES for realistic PRF (10-3800 range)
#         # ===================================================================
#         # The previous values gave PRF ~ 10^22 (unrealistic)
#         # Literature suggests PRF = 10-3800 for oxide barriers
#         # 
#         # Key insight: Dense oxides are still permeable to H2
#         # The barrier effect comes from ~100-1000x reduction, not 10^22x
#         # ===================================================================
        
#         # Temperature-dependent parameters
#         # D = D_0 * exp(-E_D / RT)
#         'D_ox_0': 1e-6,  # m²/s (pre-exponential) - INCREASED from 1e-15
#         'E_D_ox': 1.55e5,  # J/mol (activation energy) - DECREASED from 80000
        
#         # K = K_0 * exp(-H_sol / RT)  
#         # Note: For molecular H2 dissolution, H_sol is typically positive
#         'K_ox_0': 1e-4,   # mol/m³/Pa (pre-exponential) - INCREASED from 1e-20
#         'H_sol_ox': 1.85e5,  # J/mol (solution enthalpy) - INCREASED from 20000
        
#         # Geometric properties
#         'thickness': 1e-6,  # m (1 μm) - INCREASED from 6Å for more realistic barrier
#         'thickness_range': [1e-7, 1e-5],  # m (0.1 μm to 10 μm typical range)
        
#         # Metadata
#         'reference': 'Strehlow & Savage (1974), Zarchy & Axtmann (1979), Serra (1998)',
#         'temperature_range': [873, 1273],  # K (600-1000°C)
#         'uncertainty_factor': 10,  # Properties uncertain within 10x
#         'notes': 'Values adjusted to give PRF in 10-3800 range per literature'
#     },
    
#     # Alternative: Ultra-thin native oxide (6 Å as in original)
#     'Cr2O3_thin': {
#         'D_ox': 1e-15,
#         'K_ox': 1e-10,
#         'D_ox_0': 1e-10,
#         'E_D_ox': 50000,
#         'K_ox_0': 1e-8,
#         'H_sol_ox': 30000,
#         'thickness': 6e-10,  # 6 Å (native oxide)
#         'thickness_range': [3e-10, 2e-9],
#         'reference': 'Zarchy & Axtmann (1979) - 6Å oxide',
#         'temperature_range': [873, 1273],
#         'uncertainty_factor': 100,
#         'notes': 'Ultra-thin native oxide - small PRF expected'
#     }
# }


"""
Oxide properties for hydrogen permeation calculations.

Uses REFERENCE-TEMPERATURE Arrhenius format.
"""

# OXIDE_PROPERTIES = {
#     'Cr2O3': {
#         # Reference temperature
#         'T_ref': 1073.15,       # K (800°C)
        
#         # Reference values at T_ref
#         'D_ox_ref': 1e-15,      # m²/s - Oxide diffusivity at 800°C
#         'K_ox_ref': 1e-10,      # mol/m³/Pa - Oxide solubility at 800°C
        
#         # Activation energies
#         'E_D_ox': 1.55e5,       # J/mol
#         'H_sol_ox': 1.85e5,     # J/mol
        
#         # Geometry
#         'thickness': 1e-6,      # m (1 μm)
#         'thickness_range': [1e-7, 1e-5],
        
#         # Metadata
#         'reference': 'Strehlow & Savage (1974), Serra (1998)',
#         'temperature_range': [873, 1273],  # K
#         'uncertainty_factor': 10,
#     },
    
#     'Cr2O3_thin': {
#         'T_ref': 1073.15,
#         'D_ox_ref': 1e-18,
#         'K_ox_ref': 1e-12,
#         'E_D_ox': 50000,
#         'H_sol_ox': 30000,
#         'thickness': 6e-10,     # 6 Å
#         'reference': 'Zarchy & Axtmann (1979)',
#     }
# }


"""
OXIDES — Combined Hydrogen Transport Properties Dictionary
===========================================================
All units SI throughout:
    D          : m²/s
    K (Ks)     : mol/m³/Pa^0.5    [Sieverts convention — all entries use Pa^0.5]
    Phi        : mol/m/s/Pa^0.5
    E_D, H_sol : J/mol
    thickness  : m
    T          : K
    PRF        : dimensionless

SOLUBILITY CONVENTION:
    All oxide solubility (K_ox) values use the SIEVERTS convention (Pa^0.5),
    consistent with the Nemanic 2023 experimental baseline.
    Phi_ox_ref = D_ox_ref × K_ox_ref  [identity check at T_ref]
    K_ox(T) = K_ox_0 × exp(−H_sol_ox / RT)
    D_ox(T) = D_ox_0  × exp(−E_D_ox   / RT)
    Phi_ox(T) = D_ox(T) × K_ox(T) = Phi_ox_0 × exp(−Q_p_ox / RT)
    where Q_p_ox = E_D_ox + H_sol_ox

DATA QUALITY FLAGS (uncertainty_factor):
    ×5   : Two or more independent measurements agree at same T, or
           direct Phi_ox from series-resistance back-calc (Tier 1, thick oxide).
    ×10  : Single measurement point + DFT activation energy cross-check (this study).
    ×100 : Phi estimated from digitized figure with unknown specimen thickness.
    ×1000: PRF-only data, no full back-calculation possible.

PRIMARY SOURCES
---------------
[N23]  Nemanic V, Zajec B, McGuiness P, Kovac J, Zuzek B, Spiler J.
       "Impact of surface oxide on hydrogen permeability of chromium membranes."
       Int. J. Hydrogen Energy 48 (2023) 9723–9733.
       DOI: 10.1016/j.ijhydene.2022.11.267

[C11]  Chen CF, Yu HB, Zheng SQ.
       "First-principles study of hydrogen diffusion mechanism in Cr2O3."
       Sci. China Tech. Sci. 54(1), 88–94 (2011).
       DOI: 10.1007/s11431-010-4112-3

[H95]  Hollenberg GW, Simonen EP, Kalinin G, Terlain A.
       "Tritium/hydrogen barrier development."
       Fusion Eng. Design 28 (1995) 190–208.
       DOI: 10.1016/0920-3796(95)90039-X

[R79]  Roberts RM, Elleman TS, Palour H III, Verghese K.
       "Hydrogen permeability of sintered aluminum oxide."
       J. Am. Ceram. Soc. 62 (1979) 495–499.
       [Cited as Ref 22 in H95; Q_p = 318 kJ/mol taken from H95 p.200]

[B04]  Belonoshko AB, Rosengren A, Dong Q, Hultquist G, Leygraf C.
       "First-principles study of H diffusion in corundum Al2O3."
       Phys. Rev. B 69 (2004) 024302.
       [Cited as Ref 13 in C11; D(T) formula: D = 2.173e-7 exp(−1.24 eV/kT)]

[S79]  Swansiger WA, Bastasz R.
       "Tritium and deuterium permeation in stainless steels:
        Influence of thin oxide films."
       J. Nucl. Mater. 85–86 (1979) 335–339.
       DOI: 10.1016/0022-3115(79)90512-9

[M19]  Medasani BK, Sushko ML, Rosso KM, Schreiber DK, Bruemmer SM.
       "Temperature Dependence of Self-Diffusion in Cr2O3 from First Principles."
       J. Phys. Chem. C 123(36) (2019) 22139–22150.
       DOI: 10.1021/acs.jpcc.9b03218
       [NOTE: reports Cr/O ION self-diffusion ONLY — not H atom diffusion.
        Used here for oxide microstructure context, not for H transport parameters.]

CALCULATION NOTES (reproduced for traceability)
------------------------------------------------
Cr2O3 back-calculation:
    D_ox_0 = D_ox_ref × exp(+E_D_ox / R T_ref)
           = 9.6e-19 × exp(70434 / (8.314 × 673))
           = 2.813e-13 m²/s

    K_ox_ref = Phi_ox_ref / D_ox_ref = 3.4e-19 / 9.6e-19 = 0.3542 mol/m³/Pa^0.5
    K_ox_0   = K_ox_ref (H_sol = 0 assumed)

Al2O3 back-calculation:
    D_ox_0 = 2.173e-7 m²/s  [Belonoshko 2004, D(T) = D_0 exp(−1.24 eV/kT)]
    D_ox_ref(973K) = 2.173e-7 × exp(−119641 / (8.314 × 973)) = 8.203e-14 m²/s
    H_sol_ox = Q_p − E_D_ox = 318000 − 119641 = 198359 J/mol
    K_ox_0 = K_ox_ref × exp(+H_sol / R T_ref) = 3.830e-10 × exp(198359/(8.314×973))
           = 17.07 mol/m³/Pa^0.5
    Phi_ox_ref estimated: J_flux(973K,101325Pa) ~ 1e-10 mol/m²/s [Hollenberg Fig 2],
        d_specimen ~ 1e-10 m [Roberts "~10^-8 cm"], P = 101325 Pa
        Phi = J×d/sqrt(P) = 1e-10 × 1e-10 / sqrt(101325) = 3.14e-23 mol/m/s/Pa^0.5
        *** THIS IS AN ESTIMATE — Roberts 1979 original paper needed for exact value ***
"""

import numpy as np
R_GAS = 8.314  # J/mol/K


OXIDES = {

    # =========================================================================
    # Cr2O3 — PRIMARY ENTRY
    # D_ox_ref, Phi_ox_ref : Nemanic 2023 [N23], Table 3, Sample 4 (48 nm oxide on 1 mm Cr)
    # E_D_ox               : Chen 2011 [C11] DFT GGA/PW91, 0.73 eV (rate-limiting step)
    # D_ox_0               : back-calculated (see module docstring)
    # H_sol_ox             : assumed 0 (H solubility in Cr2O3 not measured; isothermal basis)
    # K_ox_0               : = K_ox_ref (H_sol=0 → temperature-independent solubility)
    # =========================================================================
    'Cr2O3': {
        # ---- Reference conditions ----
        'T_ref':       673,         # K (400°C) — Nemanic 2023 experimental temperature [N23]
        'D_ox_ref':    9.6e-19,     # m²/s     — Nemanic 2023, Table 3, Sample 4 [N23]
        'K_ox_ref':    3.5417e-1,   # mol/m³/Pa^0.5 — back-calc: Phi_ox_ref / D_ox_ref

        # ---- Activation energies ----
        'E_D_ox':      70434,       # J/mol  (0.73 eV) — Chen 2011 DFT, rate-limiting hop [C11]
        'H_sol_ox':    0.0,         # J/mol  — NOT MEASURED; isothermal assumption;
                                    # H solubility in Cr2O3 not available in literature

        # ---- Pre-exponential factors ----
        # D_ox_0 = D_ox_ref × exp(+E_D_ox / R T_ref)  →  verified: D_ox_0×exp(-E_D_ox/RT_ref) = D_ox_ref ✓
        'D_ox_0':      2.8131e-13,  # m²/s   — back-calc from D_ox_ref and E_D_ox [N23+C11]
        'K_ox_0':      3.5417e-1,   # mol/m³/Pa^0.5 — = K_ox_ref (H_sol=0)

        # ---- Permeability at T_ref (identity check: = D_ox_ref × K_ox_ref) ----
        'Phi_ox_ref':  3.4e-19,     # mol/m/s/Pa^0.5 — Nemanic 2023, Table 3, Sample 4 [N23]
                                    # VERIFY: 9.6e-19 × 3.5417e-1 = 3.400e-19 ✓

        # ---- Default geometry ----
        'thickness':        48e-9,  # m (48 nm) — Nemanic 2023 Sample 4, XPS+FIB confirmed [N23]
        'thickness_range':  [2.7e-9, 48e-9],  # m — range across Nemanic 2023 samples (S6 to S4)

        # ---- Metadata ----
        'reference': (
            'D_ox_ref, Phi_ox_ref: Nemanic et al. 2023, Table 3, Sample 4 (48 nm Cr2O3 on 1 mm Cr). '
            'DOI: 10.1016/j.ijhydene.2022.11.267. '
            'E_D_ox (0.73 eV): Chen et al. 2011 DFT GGA/PW91, Fig 5(a), rate-limiting hop A→B. '
            'DOI: 10.1007/s11431-010-4112-3. '
            'D_ox_0, K_ox_0: back-calculated (see module docstring). '
            'H_sol_ox = 0: assumed; not measured in either source. '
            'Oxide identity (Cr2O3 alpha-corundum): XPS Cr 2p confirmed in [N23].'
        ),
        'temp_range_K':       [473, 773],   # K — Nemanic 2023 experimental range (200–500°C)
        'temperature_range':  [473, 773],
        'uncertainty_factor': 10,
        # Justification for ×10: D_ox_ref and Phi_ox_ref from single T measurement (Nemanic Table 3).
        # E_D_ox from DFT on perfect bulk crystal (GGA underestimates barriers).
        # Real defective oxide coatings will have higher effective E_D.
        # PRF_Sample4 not separately tabulated; single temperature only.
        # Bounds: expect factor 3–30× on Phi_ox at any given T.

        'notes': (
            'IMPORTANT: E_D_ox = 70.4 kJ/mol from DFT perfect bulk crystal. '
            'Real polycrystalline oxide with grain boundaries and defects will have '
            'higher effective activation energy — treat as lower bound. '
            'H_sol_ox = 0 (isothermal K_s assumption) because H solubility in Cr2O3 '
            'has not been experimentally determined. '
            'Chen 2011 Table 3 H-diffusion values (converted to m²/s): '
            '  273K: 2.34e-22,  473K: 6.55e-17,  773K: 5.03e-14,  1273K: 3.08e-12. '
            'Medasani 2019 [M19] provides Cr/O ION self-diffusion — do NOT substitute '
            'for H diffusion coefficient. '
            'Swansiger 1979 [S79] finds SS-native Cr-rich oxide acts as SURFACE-KINETICS '
            'barrier (not bulk diffusion barrier) — series-resistance model invalid '
            'for thin native SS oxides; use this entry for engineered Cr2O3 coatings only.'
        ),
    },


    # =========================================================================
    # Cr2O3 — THIN COATING (28 nm, Nemanic 2023 Sample 3)
    # Phi_ox_ref : Nemanic 2023 [N23], Table 3, Sample 3 (28 nm)
    # D_ox_ref   : NOT directly measured for Sample 3; estimated using Chen 2011 D_0 [C11]
    # PRF        : 3900 (best barrier in study) at 673 K [N23]
    # =========================================================================
    'Cr2O3_28nm': {
        # ---- Reference conditions ----
        'T_ref':       673,         # K — Nemanic 2023 experimental temperature [N23]
        'D_ox_ref':    9.6e-19,     # m²/s — ESTIMATED using Chen 2011 D_0 (same T as S4) [C11]
                                    # NOT directly measured for this sample — use with caution
        'K_ox_ref':    7.2917e-2,   # mol/m³/Pa^0.5 — back-calc: Phi_ox_ref / D_ox_ref(estimated)

        # ---- Activation energies ----
        'E_D_ox':      70434,       # J/mol (0.73 eV) — Chen 2011 DFT [C11], same as Cr2O3 entry
        'H_sol_ox':    0.0,         # J/mol — assumed isothermal (same limitation as above)

        # ---- Pre-exponential factors ----
        'D_ox_0':      2.8131e-13,  # m²/s — same Chen 2011 D_0 (shared between S3 and S4)
        'K_ox_0':      7.2917e-2,   # mol/m³/Pa^0.5 — = K_ox_ref (H_sol=0)

        # ---- Permeability at T_ref ----
        'Phi_ox_ref':  7.0e-20,     # mol/m/s/Pa^0.5 — Nemanic 2023, Table 3, Sample 3 [N23]
                                    # NOTE: lower Phi_ox than 48 nm sample (better barrier)
                                    # PRF_Sample3 = 3900 (highest in study)

        # ---- Default geometry ----
        'thickness':       28e-9,   # m (28 nm) — Nemanic 2023 Sample 3, XPS+FIB, ±4 nm [N23]
        'thickness_range': [20e-9, 35e-9],  # m — ±4 nm uncertainty from XPS

        # ---- Metadata ----
        'reference': (
            'Phi_ox_ref: Nemanic et al. 2023, Table 3, Sample 3 (28 nm Cr2O3). PRF = 3900. '
            'DOI: 10.1016/j.ijhydene.2022.11.267. '
            'D_ox_ref: ESTIMATED — not directly measured; Chen 2011 D_0 used [C11]. '
            'DOI: 10.1007/s11431-010-4112-3. '
            'K_ox_ref = Phi_ox_ref / D_ox_ref (ESTIMATED D). '
            'E_D_ox: same source as Cr2O3 entry. '
            'Single-T data; Arrhenius parameters should not be extrapolated far from 673 K.'
        ),
        'temp_range_K':       [673, 673],   # K — single measurement temperature only
        'temperature_range':  [673, 673],
        'uncertainty_factor': 30,
        # Justification for ×30: Phi_ox_ref measured [N23]; D_ox_ref estimated (not measured).
        # K_ox_ref therefore doubly uncertain. Single T only — no Arrhenius fit possible.

        'notes': (
            'Best-performing oxide in Nemanic 2023 study (PRF = 3900 vs PRF = baseline). '
            'Pressure dependence shifted from p^0.5 (bare Cr, Sieverts) to p^1 (surface-limited) '
            'after oxidation — only one high-pressure data point measurable. '
            'D_ox_ref at 673 K ESTIMATED from Chen 2011 D_0; not independently validated. '
            'Use Cr2O3 (48 nm) entry for modelling unless specific 28 nm thickness required.'
        ),
    },


    # =========================================================================
    # Al2O3 — PRIMARY ENTRY
    # =========================================================================
    # SOURCE BREAKDOWN:
    #   Phi_ox_0, Q_p        : Roberts 1979 [R79], Eq.(4)/(5) — DIRECTLY FROM ORIGINAL PAPER
    #   E_D_ox, D_ox_0       : Belonoshko 2004 [B04] via Chen 2011 [C11] — DFT, cross-validated
    #   H_sol_ox             : derived Q_p − E_D_ox (both primary sources)
    #   K_ox_0, K_ox_ref     : derived from Phi/D split (see equation chain below)
    #
    # PRESSURE EXPONENT: Roberts 1979 measured n = 0.43, not 0.5 (Sieverts).
    # Raw Roberts values therefore have units Pa^0.43.
    # All K_ox and Phi_ox values stored here are RE-ASSUMED to Sieverts (Pa^0.5)
    # at P_ref = 101325 Pa (1 atm) using the equation:
    #   Phi_S(T, P_ref) = Phi_R(T) × P_ref^(n_R - n_S)
    #                   = Phi_R(T) × 101325^(0.43 - 0.5)
    #                   = Phi_R(T) × 101325^(-0.07)
    #                   = Phi_R(T) × 0.4463
    # This scaling is T-INDEPENDENT — it applies equally to Phi_0, K_0, K_ref.
    # E_D_ox and H_sol_ox are NOT affected (they are in the exponent, not pre-factor).
    #
    # FULL EQUATION CHAIN (all steps, reproduced for traceability):
    # ---------------------------------------------------------------
    # Step 0  — Roberts 1979 Eq.(4), paper units:
    #   φ = exp(48.95) × exp(−318200/RT)  [H-atom·cm·cm⁻²·s⁻¹·kPa⁻⁰·⁴³]
    #
    # Step 1  — Unit conversion to SI (Pa^0.43):
    #   Phi_0_SI = exp(48.95) / N_A × (1e-2) × (1e4) × (1e3)^0.43
    #            = 1.8143e21 / 6.022e23 × 1e-2 × 1e4 × 19.498
    #            = 5.8744  mol/m/s/Pa^0.43
    #   factors: /N_A = H-atoms→mol; ×1e-2 = cm→m (thickness);
    #            ×1e4 = cm⁻²→m⁻² (area); ×(1e3)^0.43 = kPa→Pa
    #
    # Step 2  — H_sol from Roberts Q_p and Belonoshko E_D:
    #   H_sol = Q_p − E_D = 318200 − 119641 = 198559  J/mol
    #
    # Step 3  — D_ox_ref at T_ref (Belonoshko 2004 formula):
    #   D_ox_ref = D_ox_0 × exp(−E_D / R T_ref)
    #            = 2.173e-7 × exp(−119641 / (8.314 × 1623))
    #            = 3.065e-11  m²/s
    #
    # Step 4  — Phi_ox_ref at T_ref (Roberts formula, Pa^0.43):
    #   Phi_R_ref = Phi_0_SI × exp(−Q_p / R T_ref)
    #             = 5.8744 × exp(−318200 / (8.314 × 1623))
    #             = 3.370e-10  mol/m/s/Pa^0.43
    #
    # Step 5  — K_ox_ref (Pa^0.43):
    #   K_R_ref = Phi_R_ref / D_ox_ref = 3.370e-10 / 3.065e-11 = 10.997  mol/m³/Pa^0.43
    #
    # Step 6  — K_ox_0 (Pa^0.43):
    #   K_R_0 = K_R_ref × exp(+H_sol / R T_ref)
    #         = 10.997 × exp(198559 / (8.314 × 1623))
    #         = 10.997 × 2.458e6 = 2.703e7  mol/m³/Pa^0.43
    #
    # Step 7  — Sieverts re-assumption at P_ref = 101325 Pa:
    #   scale = P_ref^(n_R − n_S) = 101325^(0.43 − 0.50) = 101325^(−0.07)
    #         = exp(−0.07 × ln(101325)) = exp(−0.07 × 11.526) = exp(−0.8068) = 0.4463
    #   Phi_S_0   = Phi_0_SI   × scale = 5.8744  × 0.4463 = 2.6216  mol/m/s/Pa^0.5
    #   K_S_ref   = K_R_ref    × scale = 10.997  × 0.4463 = 4.9077  mol/m³/Pa^0.5
    #   K_S_0     = K_R_0      × scale = 2.703e7 × 0.4463 = 1.207e7 mol/m³/Pa^0.5
    #   Phi_S_ref = Phi_R_ref  × scale = 3.370e-10 × 0.4463 = 1.504e-10  mol/m/s/Pa^0.5
    #
    # Step 8  — Identity checks (all verified ✓):
    #   D_ox_ref × K_S_ref = 3.065e-11 × 4.908 = 1.504e-10 = Phi_S_ref ✓
    #   E_D + H_sol = 119641 + 198559 = 318200 = Q_p ✓
    #   K_S_0 × exp(−H_sol/RT_ref) = K_S_ref ✓
    #   All three routes to Phi_S_ref agree to 0.0000% ✓
    # =========================================================================
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
    },


    # =========================================================================
    # Cr2O3 on SS (surface-kinetics limited) — Swansiger 1979
    # This entry represents the EFFECTIVE composite behaviour, NOT standalone oxide D/K.
    # Series-resistance model INVALID: Swansiger confirmed surface-kinetics controls.
    # Q_p_composite = 87 kJ/mol (elevated vs bare ~62–66 kJ/mol) [S79]
    # =========================================================================
    'Cr2O3_SS_surface_kinetics': {
        # ---- Reference conditions ----
        # No standalone D_ox or K_ox — surface-kinetics limited, not bulk diffusion limited.
        # Providing composite effective parameters and PRF only.
        'T_ref':       500,         # K — approximate mid-range of Swansiger 1979 data [S79]
        'D_ox_ref':    None,        # NOT APPLICABLE — surface-limited mechanism
        'K_ox_ref':    None,        # NOT APPLICABLE

        # ---- Activation energies ----
        'E_D_ox':      None,        # NOT APPLICABLE — bulk diffusion not rate-limiting
        'H_sol_ox':    None,        # NOT APPLICABLE

        # ---- Pre-exponential factors ----
        'D_ox_0':      None,
        'K_ox_0':      None,

        # ---- Permeability at T_ref ----
        'Phi_ox_ref':  None,        # NOT EXTRACTABLE — series resistance invalid

        # ---- Composite effective parameters (what CAN be used) ----
        'Q_p_composite':    87000,  # J/mol — Swansiger 1979, 21-6-9 SS, Fig 4 [S79]
                                    # cf. bare SS Q_p ~ 62–66 kJ/mol; oxide raises apparent Ea
        'PRF_range':        [100, 1000],  # dimensionless — 2–3 orders of magnitude [S79]
        'PRF_T_ref':        500,    # K
        'oxide_thickness':  15e-9,  # m — 21-6-9 SS, ~15 nm Cr-rich oxide, AES confirmed [S79]
        'oxide_thickness_A286': 4e-9,  # m — modified A-286, ~4 nm Ti-enriched oxide [S79]

        # ---- Default geometry ----
        'thickness':       15e-9,   # m — 21-6-9 reference case
        'thickness_range': [4e-9, 15e-9],

        # ---- Metadata ----
        'reference': (
            'Swansiger WA, Bastasz R. '
            '"Tritium and deuterium permeation in stainless steels: '
            'Influence of thin oxide films." '
            'J. Nucl. Mater. 85–86 (1979) 335–339. '
            'DOI: 10.1016/0022-3115(79)90512-9. '
            'Substrates: 21-6-9 SS (15 nm Cr-rich oxide) and modified A-286 (4 nm Ti-enriched). '
            'Oxide formed by Nitradd-nitric acid etch. '
            'Q_p_composite from Fig 4 text. PRF from Figs 4–5.'
        ),
        'temp_range_K':       [325, 700],
        'temperature_range':  [325, 700],
        'uncertainty_factor': 1000,
        # Justification for ×1000: surface-kinetics mechanism — series resistance invalid.
        # Phi_oxide not extractable. PRF only with order-of-magnitude bounds.

        'notes': (
            'MECHANISM: surface-kinetics limited (NOT bulk diffusion limited). '
            'Evidence from Swansiger 1979: '
            '(1) Pressure exponent → p^1 (linear) after oxidation vs p^0.5 (Sieverts) bare. '
            '    Linear p-dependence indicates upstream surface recombination rate-limiting. '
            '(2) Placing Pd over oxide on UPSTREAM face eliminates barrier — confirms '
            '    upstream surface (not bulk or downstream) controls transport. '
            '(3) Modified A-286 (4 nm oxide) and 21-6-9 (15 nm) show similar PRF despite '
            '    3.75× thickness difference — thickness does NOT control barrier. '
            '(4) Q_p_composite = 87 kJ/mol > bare SS ~62–66 kJ/mol: elevated activation '
            '    energy consistent with surface H2 dissociation/recombination kinetics. '
            'CONSEQUENCE FOR MODELLING: '
            'Series-resistance model (Phi_ox from oxide layer thickness) is INVALID here. '
            'Use a surface-kinetics boundary condition (surface rate constant k_s) instead. '
            'Oxide AES composition [S79]: Cr ~19%, O ~43%, Fe ~29%, C ~20% (contamination).'
        ),
    },


    # =========================================================================
    # Cr2O3 — DEFECT CONTEXT (Medasani 2019)
    # NOT an H permeation entry. Provides Cr/O vacancy data for oxide microstructure modelling.
    # =========================================================================
    'Cr2O3_defect_context': {
        # ---- Reference conditions ----
        'T_ref':       None,
        'D_ox_ref':    None,        # Cr/O ion self-diffusion — NOT H atom diffusion
        'K_ox_ref':    None,

        # ---- NOT APPLICABLE for H transport ----
        'E_D_ox':      None,
        'H_sol_ox':    None,
        'D_ox_0':      None,
        'K_ox_0':      None,
        'Phi_ox_ref':  None,

        # ---- Defect formation energies (HSE, O-rich) from Medasani 2019 Table 3 ----
        'V_Cr_formation_eV':    0.90,   # eV — Cr vacancy, lowest EF; makes Cr2O3 p-type [M19]
        'V_O_formation_eV':     5.55,   # eV — O vacancy, much higher EF under O-rich [M19]
        'Cr_i_formation_eV':   11.65,   # eV — Cr interstitial, very high; negligible [M19]

        # ---- Ion self-diffusion at ~1000 K, 1 atm O2 (approximate, Fig 9a of M19) ----
        'D_Cr_ion_1000K':   1e-26,  # m²/s — Cr ion in Cr2O3 bulk at 1 atm [M19 Fig 9a]
        'D_O_ion_1000K':    1e-34,  # m²/s — O ion in Cr2O3 bulk at 1 atm [M19 Fig 9a]

        # ---- Default geometry ----
        'thickness':       None,
        'thickness_range': None,

        # ---- Metadata ----
        'reference': (
            'Medasani BK, Sushko ML, Rosso KM, Schreiber DK, Bruemmer SM. '
            '"Temperature Dependence of Self-Diffusion in Cr2O3 from First Principles." '
            'J. Phys. Chem. C 123(36) (2019) 22139–22150. '
            'DOI: 10.1021/acs.jpcc.9b03218. '
            'HSE hybrid DFT (VASP, PAW, AFM). Defect formation energies: Table 3. '
            'Ion self-diffusion coefficients: Fig 9.'
        ),
        'temp_range_K':       [300, 1600],
        'temperature_range':  [300, 1600],
        'uncertainty_factor': None,

        'notes': (
            'WARNING: This entry contains Cr and O ION self-diffusion data. '
            'DO NOT use D_Cr_ion or D_O_ion as H atom diffusion coefficients. '
            'Ion diffusion (10^-26 to 10^-34 m²/s) is entirely different physics '
            'from H interstitial diffusion (10^-19 m²/s from Nemanic at 673 K). '
            'USE FOR: '
            '(a) Passive film growth rate estimation on Hastelloy X surface. '
            '(b) Oxide layer thickness evolution with time (parabolic oxidation kinetics). '
            '(c) Defect chemistry context: V_Cr dominates at O-rich conditions → '
            '    Cr2O3 is p-type; V_O dominates at Cr/Cr2O3 interface and low pO2. '
            '    This determines grain boundary character and effective H trapping depth. '
            'Dominant defect transitions (Fig 9 of M19): '
            '    300–1100 K at 1 atm: V_Cr dominant (Cr diffuses faster than O). '
            '    >1100 K or low pO2:  V_O dominant. '
            '    Cr/Cr2O3 interface:  V_O dominant at all T.'
        ),
    },

}


# =============================================================================
# HELPER: evaluate any entry at arbitrary T (returns D, K, Phi)
# =============================================================================

def eval_oxide(entry_dict, T_K):
    """
    Evaluate D_ox(T), K_ox(T), Phi_ox(T) for an OXIDES entry at temperature T_K.
    Returns None for entries with surface-kinetics or defect-context data.

    Parameters
    ----------
    entry_dict : dict   — one entry from OXIDES
    T_K        : float  — temperature in Kelvin

    Returns
    -------
    dict with keys: T_K, D_ox, K_ox, Phi_ox  (all None if not applicable)
    """
    D0  = entry_dict.get('D_ox_0')
    K0  = entry_dict.get('K_ox_0')
    ED  = entry_dict.get('E_D_ox')
    Hs  = entry_dict.get('H_sol_ox')

    if any(v is None for v in [D0, K0, ED, Hs]):
        return {'T_K': T_K, 'D_ox': None, 'K_ox': None, 'Phi_ox': None}

    D_ox   = D0  * np.exp(-ED / (R_GAS * T_K))
    K_ox   = K0  * np.exp(-Hs / (R_GAS * T_K))
    Phi_ox = D_ox * K_ox
    return {'T_K': T_K, 'D_ox': D_ox, 'K_ox': K_ox, 'Phi_ox': Phi_ox}


# =============================================================================
# QUICK SELF-VERIFICATION (run module directly: python OXIDES_combined.py)
# =============================================================================

if __name__ == '__main__':
    print("=" * 60)
    print("OXIDES dict — self-verification")
    print("=" * 60)

    for name, ox in OXIDES.items():
        T = ox.get('T_ref')
        D = ox.get('D_ox_ref')
        K = ox.get('K_ox_ref')
        P = ox.get('Phi_ox_ref')
        D0 = ox.get('D_ox_0')
        K0 = ox.get('K_ox_0')
        ED = ox.get('E_D_ox')
        Hs = ox.get('H_sol_ox')

        print(f"\n[{name}]")
        if None in (T, D, K, P, D0, K0, ED, Hs):
            print("  Skipping (surface-kinetics or defect-context entry)")
            continue

        # Identity check 1: D*K = Phi at T_ref
        DK = D * K
        err1 = abs(DK / P - 1.0)
        flag1 = "✓" if err1 < 1e-4 else f"FAIL (err={err1:.2e})"
        print(f"  D_ox_ref × K_ox_ref = Phi_ox_ref?  {flag1}  [{DK:.4e} vs {P:.4e}]")

        # Identity check 2: D_ox_0 × exp(−E_D/RT_ref) = D_ox_ref
        D_recalc = D0 * np.exp(-ED / (R_GAS * T))
        err2 = abs(D_recalc / D - 1.0)
        flag2 = "✓" if err2 < 1e-4 else f"FAIL (err={err2:.2e})"
        print(f"  D_ox_0 × exp(−E_D/RT_ref) = D_ox_ref?  {flag2}  [{D_recalc:.4e} vs {D:.4e}]")

        # eval_oxide round-trip
        out = eval_oxide(ox, T)
        err3 = abs(out['Phi_ox'] / P - 1.0)
        flag3 = "✓" if err3 < 1e-4 else f"FAIL (err={err3:.2e})"
        print(f"  eval_oxide(T_ref) = Phi_ox_ref?  {flag3}  [{out['Phi_ox']:.4e}]")





OXIDES = {

    # =========================================================================
    # SOURCE: Nemanic et al. 2023
    # "Impact of surface oxide on hydrogen permeability of chromium membranes"
    # Int. J. Hydrogen Energy 48 (2023) 9723-9733
    # DOI: 10.1016/j.ijhydene.2022.11.267
    # Substrate: Pure bulk Cr (99.95%), 1 mm thick
    # Oxide: Cr2O3 (XPS-confirmed)
    # Gas: H2 (Pd-coated reference) / D2 (dynamic QMS mode)
    # =========================================================================

    'Cr_bare_Pd_coated_Nemanic2023': {
        'Oxide identity': {
            'oxide_type': None,
            'substrate': 'Pure Cr (99.95%)',
            'formation_condition': 'Pd-coated both sides (~100 nm Pd, magnetron sputtering) — oxide-free reference',
            'oxide_thickness': None,
            'thickness_known': False,
            'data_tier': None,
        },
        'As-reported parameters (bare metal)': {
            'T_ref': 673,                   # K (400 C)
            'Phi_bare_ref': 2.6058e-12,    # mol/m/s/Pa^0.5 — computed from Arrhenius at 673 K
            'Q_p': 65620,                   # J/mol (0.68 eV * 96485)
            'Phi_0': 3.23e-7,              # mol/m/s/Pa^0.5
            'D_ref': 3.4352e-9,            # m^2/s — computed from Arrhenius at 673 K
            'Q_D': 56924,                   # J/mol (0.59 eV * 96485)
            'D_0': 9.0e-5,                 # m2/s
            'Ks_ref': 7.5858e-4,           # mol/m^3/Pa^0.5 — back-calc: Phi_ref / D_ref at 673 K
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Nemanic et al. 2023. DOI: 10.1016/j.ijhydene.2022.11.267',
            'temp_range': [473, 773],       # K (200-500 C)
            'pressure_range': [2.5e3, 1e5], # Pa (25-1000 mbar)
            'metal_thickness': [1e-3],      # m
            'gas': 'H2',
            'Notes': (
                'Pd-coated reference (Sample 5). Sieverts law (p^0.5) confirmed. '
                'Data at 500C corrected x2 for Pd migration into Cr. '
                'These parameters are for PURE Cr bulk, not the oxide layer.'
            ),
        },
    },

    'CrOx_28nm_on_Cr_Nemanic2023_Sample3': {
        'Oxide identity': {
            'oxide_type': 'Cr2O3',
            'substrate': 'Pure Cr (99.95%)',
            'formation_condition': 'Controlled oxidation: pure O2 900 mbar at 400 C, 2h + 20h',
            'oxide_thickness': 28e-9,       # m (28 +/- 4 nm, XPS depth profiling + FIB cross-section)
            'thickness_known': True,
            'data_tier': 1,
        },
        'As-reported parameters (composite)': {
            'T_ref': 673,                   # K (400 C)
            'PRF': 3900,                    # dimensionless — highest PRF in study
            'PRF_T_ref': 673,               # K
            'Phi_composite_ref': None,      # only one high-pressure point recorded
            'Q_p_composite': None,          # single point — not extractable
            'Phi_composite_0': None,
            'Phi_bare_ref': None,           # from Pd-coated Sample 5
            'metal_baseline_key': 'Cr_bare_Pd_coated_Nemanic2023',
        },
        'Back-calculated standalone oxide parameters': {
            # Phi_oxide taken directly from Nemanic 2023 Table 3.
            # Independent series-resistance check using Arrhenius Phi_bare gives 1.87e-20,
            # a factor ~3.7x lower. Discrepancy because paper uses MEASURED j_Cr(400C,1bar)
            # from Fig 4, not the Arrhenius extrapolation. Paper Table 3 value is preferred.
            'Phi_oxide_ref': 7.0e-20,      # mol/m/s/Pa^0.5 at 673 K — Table 3, Nemanic 2023
            'Q_p_oxide': None,             # single T only — not extractable
            'Phi_oxide_0': None,
            'D_ref': None,                  # D_CrOx not in Table 3 for Sample 3
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,                # cannot compute without D_CrOx
            'K_s0': None,
            'H_s': None,
            'back_calc_method': (
                'Tier 1 — Phi_oxide from Nemanic 2023 Table 3 (paper Eqs. 3-5 applied '
                'to measured j_Cr and j_CrOx at 400C, 1 bar). '
                'D_CrOx not tabulated for Sample 3 — Ks_oxide not derivable. '
                'Single T only — Arrhenius parameters not extractable. '
                'Cross-check: series resistance with Arrhenius Phi_bare gives 1.87e-20 '
                '(factor ~3.7x lower — due to Arrhenius vs measured Phi_bare).'
            ),
        },
        'Other common parameters': {
            'pressure': [1e5],              # Pa — single highest-pressure point
            'reference': 'Nemanic et al. 2023. DOI: 10.1016/j.ijhydene.2022.11.267. Tables 2 & 3.',
            'temp_range': [673],
            'pressure_range': [1e5],
            'metal_thickness': [1e-3],
            'gas': 'H2',
            'Notes': (
                'Best barrier result. PRF ~3900 at 28 nm. '
                'Pressure dependence shifts from p^0.5 to p^1 after oxidation: surface-limited regime. '
                'Oxide: Cr2O3-type confirmed by XPS Cr 2p spectra. '
                'Only one data point measurable after oxidation.'
            ),
        },
    },

    'CrOx_48nm_on_Cr_Nemanic2023_Sample4': {
        'Oxide identity': {
            'oxide_type': 'Cr2O3',
            'substrate': 'Pure Cr (99.95%)',
            'formation_condition': 'Two successive oxidation cycles: pure O2 900 mbar at 400 C',
            'oxide_thickness': 48e-9,       # m (48 +/- 6 nm, XPS + FIB)
            'thickness_known': True,
            'data_tier': 1,
        },
        'As-reported parameters (composite)': {
            'T_ref': 673,
            'PRF': None,                    # not separately tabulated for Sample 4
            'PRF_T_ref': 673,
            'Phi_composite_ref': None,
            'Q_p_composite': None,
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': 'Cr_bare_Pd_coated_Nemanic2023',
        },
        'Back-calculated standalone oxide parameters': {
            # Phi_oxide and D_CrOx from Nemanic 2023 Table 3 (paper's own back-calculation).
            # Ks_oxide back-calculated here as Phi_oxide / D_oxide.
            # PRF cannot be computed without measured j_Cr(400C,1bar) from Fig 4.
            'Phi_oxide_ref': 3.4e-19,      # mol/m/s/Pa^0.5 at 673 K — Table 3, Nemanic 2023
            'Q_p_oxide': None,             # single T only — not extractable
            'Phi_oxide_0': None,
            'D_ref': 9.6e-19,              # m^2/s at 673 K — Table 3, Nemanic 2023
            'E_D': None,
            'D_0': None,
            'Ks_ref': 3.54e-1,             # mol/m^3/Pa^0.5 at 673 K — back-calc: Phi_oxide/D_oxide
            'K_s0': None,
            'H_s': None,
            'back_calc_method': (
                'Tier 1 — Phi_oxide and D_CrOx from Nemanic 2023 Table 3 '
                '(paper Eqs. 3-5 applied to measured j_Cr and j_CrOx). '
                'Ks_oxide = Phi_oxide / D_oxide = 3.54e-1 mol/m^3/Pa^0.5 (back-calc). '
                'PRF not computed — requires measured j_Cr(400C,1bar) from Fig 4. '
                'Single T only — Arrhenius parameters not extractable.'
            ),
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Nemanic et al. 2023. DOI: 10.1016/j.ijhydene.2022.11.267. Table 3.',
            'temp_range': [673],
            'pressure_range': [],
            'metal_thickness': [1e-3],
            'gas': 'H2',
            'Notes': (
                'Thicker oxide (48 nm) than Sample 3 (28 nm) but higher Phi_oxide — '
                'consistent with non-uniform oxide quality / scatter between samples. '
                'Reference DiHea value from literature: Phi_CrOx = 1.5e-17 mol/m/s/Pa^0.5 (Table 3).'
            ),
        },
    },

    'CrOx_2p7nm_native_on_Cr_Nemanic2023_Sample6': {
        'Oxide identity': {
            'oxide_type': 'Cr2O3',
            'substrate': 'Pure Cr (99.95%)',
            'formation_condition': 'Native oxide: polishing + vacuum thermal treatment 18h at 400 C',
            'oxide_thickness': 2.7e-9,      # m (2.7 +/- 0.2 nm, XPS depth profiling)
            'thickness_known': True,
            'data_tier': 2,
        },
        'As-reported parameters (composite)': {
            'T_ref': 673,
            'PRF': 5.2,                     # dimensionless — average iPRF_DL (range 3-9, ~73% scatter)
            'PRF_T_ref': 673,
            'Phi_composite_ref': None,
            'Q_p_composite': None,
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': 'Cr_bare_Pd_coated_Nemanic2023',
        },
        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,          # unreliable due to PRF scatter
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': (
                'Tier 2 — thickness known (2.7 nm) but back-calculation unreliable. '
                'iPRF_DL scatter ~73% across identically-prepared samples. '
                'Attributed to uncontrolled water adsorption during heating period. '
                'Phi_oxide not determinable from this data.'
            ),
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Nemanic et al. 2023. DOI: 10.1016/j.ijhydene.2022.11.267. Table 2.',
            'temp_range': [673],
            'pressure_range': [],
            'metal_thickness': [1e-3],
            'gas': 'H2',
            'Notes': (
                'Native oxide at XPS detection limit. Cr2O3-type confirmed. '
                'iPRF scatter (range 3-9, avg 5.2) prevents reliable back-calculation.'
            ),
        },
    },

    # =========================================================================
    # SOURCE: Swansiger & Bastasz 1979
    # "Tritium and Deuterium Permeation in Stainless Steels: Influence of Thin Oxide Films"
    # J. Nuclear Materials 85-86 (1979) 335-339
    # DOI: 10.1016/0022-3115(79)90512-9
    # Substrates: 309S, 21-6-9, modified A-286 SS
    # Oxide: Thin mixed Cr/Fe oxide from Nitradd-nitric acid etch (~4-15 nm)
    # Gas: T2/D2 (T2 data normalized to D2 by multiply by 1/sqrt(5))
    # T range: 325-700 K
    # =========================================================================

    'SS_bare_PdCoated_Swansiger1979': {
        'Oxide identity': {
            'oxide_type': None,
            'substrate': '309S / 21-6-9 / modified A-286 SS (combined)',
            'formation_condition': 'Sputter-cleaned + Pd-coated both sides — oxide-free baseline',
            'oxide_thickness': None,
            'thickness_known': False,
            'data_tier': None,
        },
        'As-reported parameters (bare metal)': {
            'T_ref': None,                  # Arrhenius fit over full range 325-700 K
            # Individual grade Arrhenius (units: mol D2/m/s/MPa^0.5 in paper):
            # 309S:     Ea=62300 J/mol, Phi_0=0.12e-3 mol D2/m/s/MPa^0.5
            # 21-6-9:   Ea=66000 J/mol, Phi_0=0.36e-3
            # mod A-286:Ea=62100 J/mol, Phi_0=1.14e-3
            # Combined: Ea=59900 J/mol, Phi_0=8.0e-5 mol D2/m/s/MPa^0.5
            # Conversion: 1 mol D2/m/s/MPa^0.5 = 1e-3 mol/m/s/Pa^0.5
            #             (divide by sqrt(1e6) = 1e3 for Pa -> MPa)
            #             Actually: Phi [mol/m/s/Pa^0.5] = Phi [mol D2/m/s/MPa^0.5] * 1e-3
            'Q_p_309S': 62300,              # J/mol
            'Phi_0_309S': 1.2e-7,          # mol/m/s/Pa^0.5 (converted from 0.12e-3 mol D2/m/s/MPa^0.5)
            'Q_p_21_6_9': 66000,            # J/mol
            'Phi_0_21_6_9': 3.6e-7,        # mol/m/s/Pa^0.5
            'Q_p_A286': 62100,              # J/mol
            'Phi_0_A286': 1.14e-6,         # mol/m/s/Pa^0.5
            'Q_p_combined': 59900,          # J/mol
            'Phi_0_combined': 8.0e-8,      # mol/m/s/Pa^0.5 (combined all three)
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Swansiger & Bastasz 1979. DOI: 10.1016/0022-3115(79)90512-9. Fig 1.',
            'temp_range': [325, 700],       # K
            'pressure_range': [1.3e3, 1e5], # Pa (1.3-100 kPa)
            'metal_thickness': [9e-5, 2.8e-4], # m (0.09-0.28 mm)
            'gas': 'D2',
            'Notes': (
                'All three SS bulk permeabilities within 30% of each other. '
                'T2 data (below 473 K) normalized to D2 by x 1/sqrt(5). '
                'Combined Arrhenius: Phi = 8.0e-5*exp(-59900/RT) mol D2/m/s/MPa^0.5.'
            ),
        },
    },

    'SS_oxide_21_6_9_NitradEtch_Swansiger1979': {
        'Oxide identity': {
            'oxide_type': 'Cr-rich mixed oxide (Fe/Cr/Ni)',  # AES: Cr~19%, O~43%, Fe~29%
            'substrate': '21-6-9 austenitic SS',
            'formation_condition': 'Nitradd-nitric acid chemical etch (20% Nitradd + 30% HNO3 in DI water)',
            'oxide_thickness': 15e-9,       # m (~15 nm from AES sputter profiling)
            'thickness_known': True,
            'data_tier': 2,
        },
        'As-reported parameters (composite)': {
            'T_ref': None,                  # Arrhenius fit over 325-700 K
            'PRF': None,                    # 2-3 orders of magnitude — exact value from Fig 4 (not tabulated)
            'PRF_T_ref': None,
            'Phi_composite_ref': None,      # from Fig 4, not tabulated
            'Q_p_composite': 87000,         # J/mol — up to 87 kJ/mol (significantly above bare ~62-66 kJ/mol)
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': 'SS_bare_PdCoated_Swansiger1979',
        },
        'Back-calculated standalone oxide parameters': {
            # Series resistance back-calculation NOT physically meaningful here.
            # Paper explicitly confirms upstream surface effects are rate-limiting —
            # oxide acts as a surface kinetics barrier, not a bulk diffusion barrier.
            # Numerical check: at T=500K with PRF=100-1000, series resistance gives
            # Phi_oxide ~ 4.4e-20 to 4.4e-21 mol/m/s/Pa^0.5, but these values are
            # unreliable because the transport mechanism is surface-limited not diffusion-limited.
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': (
                'Tier 2 — series resistance back-calculation NOT applicable. '
                'Paper confirms upstream surface effects are rate-limiting (not bulk oxide diffusion). '
                'Oxide functions as surface kinetics barrier. Series resistance assumes '
                'diffusion-limited transport — invalid here. '
                'Phi_composite requires digitizing Fig 4 of Swansiger & Bastasz 1979. '
                'Even if extracted, Phi_oxide would represent effective surface conductance, '
                'not true diffusion permeability.'
            ),
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Swansiger & Bastasz 1979. DOI: 10.1016/0022-3115(79)90512-9. Figs 1 & 4.',
            'temp_range': [325, 700],
            'pressure_range': [1.3e3, 1e5],
            'metal_thickness': [9e-5, 2.8e-4],
            'gas': 'D2',
            'Notes': (
                'Permeability reduced 2-3 orders of magnitude vs bare SS. '
                'Q_p_composite up to 87 kJ/mol (vs ~62-66 kJ/mol bare). '
                'Upstream surface effects rate-limiting (confirmed by Pd-over-oxide experiments). '
                'Downstream oxide has negligible effect on permeation. '
                'Oxide surface composition from AES: Cr~19%, O~43%, Fe~29%, C~20% (carbon contamination).'
            ),
        },
    },

    'SS_oxide_modA286_NitradEtch_Swansiger1979': {
        'Oxide identity': {
            'oxide_type': 'Ti-enriched mixed oxide (Cr/Fe/Ti)',  # Ti enrichment from heat treatment
            'substrate': 'Modified A-286 (precipitation hardened SS)',
            'formation_condition': 'Nitradd-nitric acid etch. Ti surface enrichment from prior heat treatment (1172 K, 7.2 ks).',
            'oxide_thickness': 4e-9,        # m (~4 nm from AES sputter profiling)
            'thickness_known': True,
            'data_tier': 2,
        },
        'As-reported parameters (composite)': {
            'T_ref': None,
            'PRF': None,                    # 2-3 orders of magnitude — from Fig 5 (not tabulated)
            'PRF_T_ref': None,
            'Phi_composite_ref': None,
            'Q_p_composite': None,          # not quoted for A-286 separately
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': 'SS_bare_PdCoated_Swansiger1979',
        },
        'Back-calculated standalone oxide parameters': {
            # Series resistance back-calculation NOT physically meaningful.
            # Same conclusion as 21-6-9 entry above: surface-limited mechanism confirmed.
            # Key evidence: A-286 (~4 nm oxide) shows similar PRF to 21-6-9 (~15 nm oxide),
            # proving oxide THICKNESS does not control the barrier — surface kinetics does.
            # Numerical check at T=500K: series resistance gives Phi_oxide ~ 1.2e-20 to 1.2e-21
            # mol/m/s/Pa^0.5 (PRF=100-1000) but these are physically invalid estimates.
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': (
                'Tier 2 — series resistance back-calculation NOT applicable. '
                'Surface-limited transport confirmed: A-286 (4 nm oxide) and 21-6-9 (15 nm) '
                'show similar PRF despite 3.75x thickness difference — '
                'oxide thickness does not control the barrier, surface kinetics does. '
                'Series resistance model (diffusion-limited) is invalid for this system.'
            ),
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Swansiger & Bastasz 1979. DOI: 10.1016/0022-3115(79)90512-9. Fig 5.',
            'temp_range': [325, 700],
            'pressure_range': [1.3e3, 1e5],
            'metal_thickness': [9e-5, 2.8e-4],
            'gas': 'D2',
            'Notes': (
                'Apparent scatter between individual A-286 samples despite identical preparation. '
                'Ti surface enrichment from prior heat treatment observed in AES (Fig 3). '
                'Upstream surface effects confirmed rate-limiting. '
                'Paper general finding: thick oxides (>100 nm) show T-invariant PRF — '
                'suggests fraction of oxide cracked/defective, allowing partial H access to substrate.'
            ),
        },
    },

    # =========================================================================
    # SOURCE: Yamawaki et al. 1989
    # "Effect of surface impurities on the hydrogen recombination coefficient"
    # J. Nuclear Materials 162-164 (1989) 1071-1076
    # NOTE: This paper reports SURFACE KINETICS (recombination coefficient kR),
    #       not oxide permeability. Included as surface oxide effect data.
    # kR units in paper: cm^4/s  (to convert: 1 cm^4/s = 1e-8 m^4/s)
    # kR relates to recombination: release rate ~ kR * C_s^2
    # =========================================================================

    'SS304_surface_oxide_O_impurity_Yamawaki1989': {
        'Oxide identity': {
            'oxide_type': 'O-contaminated surface (adsorbed O -> oxide on heating)',
            'substrate': 'Type 304 SS',
            'formation_condition': 'O adsorbed from atmosphere at low T, forms oxide layer on heating to 773 K',
            'oxide_thickness': None,
            'thickness_known': False,
            'data_tier': 3,
        },
        'As-reported parameters (composite)': {
            'T_ref': 773,
            'PRF': None,                    # not reported as PRF
            'PRF_T_ref': 773,
            'Phi_composite_ref': None,
            'Q_p_composite': None,
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': None,
        },
        'Surface kinetics parameters': {
            # kR functional form: kR(T, Co) = k0(T) * exp(-b * Co)  [Eq. 6]
            'T_ref': 773,                   # K
            'dominant_impurity': 'Oxygen (O)',
            'kR_units_paper': 'cm^4/s',
            'kR_at_Co_zero_approx': 1e-16,  # cm^4/s — read from Fig 1 extrapolated to Co=0
            'kR_at_Co_35at_approx': 1e-20,  # cm^4/s — read from Fig 1 at Co~35 at% O
            'kR_range': [1e-20, 1e-16],     # cm^4/s
            'b_coefficient': None,          # slope of ln(kR) vs Co — read from Fig 1
            'functional_form': 'kR(T,Co) = k0(T) * exp(-b * Co)',
            'comparison_kR_Baskes': None,   # (kR)B from Baskes model (Eq. 7)
            'comparison_kR_PickSonnenberg': None,  # (kR)PS from Pick & Sonnenberg (Eq. 8)
        },
        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': 'Tier 3 — kR data only, no Phi_composite or oxide_thickness.',
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Yamawaki et al. 1989. J. Nucl. Mater. 162-164, 1071-1076.',
            'temp_range': [773],
            'pressure_range': [],
            'metal_thickness': [],
            'gas': 'D2',
            'Notes': (
                'kR decreases from ~1e-16 to ~1e-20 cm^4/s as Co increases 0 to 35 at% O. '
                'O is most correlated impurity (then Si weakly). '
                'Mechanism: O increases surface electronegativity -> suppresses H2 dissociation. '
                'Result closest to Pick & Sonnenberg model prediction.'
            ),
        },
    },

    'Ni_surface_S_impurity_Yamawaki1989': {
        'Oxide identity': {
            'oxide_type': 'S-contaminated surface (S segregation from bulk)',
            'substrate': 'Nickel',
            'formation_condition': 'Sulfur segregation from bulk on heating — up to 40 at% S on surface',
            'oxide_thickness': None,
            'thickness_known': False,
            'data_tier': 3,
        },
        'As-reported parameters (composite)': {
            'T_ref': None,
            'PRF': None,
            'PRF_T_ref': None,
            'Phi_composite_ref': None,
            'Q_p_composite': None,
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': None,
        },
        'Surface kinetics parameters': {
            'dominant_impurity': 'Sulfur (S)',
            'kR_units_paper': 'cm^4/s',
            'kR_range': [1e-19, 1e-16],     # cm^4/s — approximate from context
            '2Ec_max_eV': 0.7,              # eV — max chemisorption activation barrier at high Cs
            '2Ec_max_Jmol': 67540,          # J/mol (0.7 eV * 96485)
            'functional_form': 'kR(T,Ci) = k0(T) * exp(-b * Ci)',
        },
        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': 'Tier 3 — surface kinetics study only.',
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Yamawaki et al. 1989. J. Nucl. Mater. 162-164, 1071-1076.',
            'temp_range': [],
            'pressure_range': [],
            'metal_thickness': [],
            'gas': 'D2',
            'Notes': (
                'kR decreases as Cs (sulfur) increases 0 to 40 at% S. '
                '2Ec increases up to 0.7 eV at high Cs. '
                'Mechanism: S increases surface electronegativity -> suppresses H2 dissociative adsorption. '
                'kR values orders of magnitude larger than other literature for Ni.'
            ),
        },
    },

    'V_surface_S_impurity_Yamawaki1989': {
        'Oxide identity': {
            'oxide_type': 'S-contaminated surface (S segregation from bulk)',
            'substrate': 'Vanadium',
            'formation_condition': 'Sulfur segregation from bulk on heating — up to 40 at% S',
            'oxide_thickness': None,
            'thickness_known': False,
            'data_tier': 3,
        },
        'As-reported parameters (composite)': {
            'T_ref': None,
            'PRF': None,
            'PRF_T_ref': None,
            'Phi_composite_ref': None,
            'Q_p_composite': None,
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': None,
        },
        'Surface kinetics parameters': {
            'dominant_impurity': 'Sulfur (S)',
            'kR_units_paper': 'cm^4/s',
            'kR_approx': 1e-25,             # cm^4/s — ~10^-25 cm^4/s (close to Holland & Anderl ~10^-25)
            'Notes': (
                'kR for V much smaller than SS and Ni. '
                'kR-Cs functional dependence differs from Ni (not simple exp form). '
                'Permeation reduced by order of magnitude when S present.'
            ),
        },
        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': 'Tier 3 — surface kinetics study only.',
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Yamawaki et al. 1989. J. Nucl. Mater. 162-164, 1071-1076.',
            'temp_range': [],
            'pressure_range': [],
            'metal_thickness': [],
            'gas': 'D2',
            'Notes': 'Surface impurity effect on recombination coefficient. Not an oxide permeability study.',
        },
    },

    # =========================================================================
    # SOURCE: Ronnebro et al. 2022 (Review)
    # Molecules 2022, 27, 6528
    # DOI: 10.3390/molecules27196528
    # NOTE: Review paper — no primary experimental data.
    # =========================================================================

    'Review_oxide_barriers_Ronnebro2022': {
        'Oxide identity': {
            'oxide_type': 'Various — Al2O3, Cr2O3, Y2O3, Er2O3 and nitrides/carbides',
            'substrate': 'Various',
            'formation_condition': 'Various — PVD, CVD, cold spray, sol-gel, pack cementation',
            'oxide_thickness': None,
            'thickness_known': False,
            'data_tier': 3,
        },
        'As-reported parameters (composite)': {
            'T_ref': None,
            'PRF': None,
            'PRF_T_ref': None,
            'Phi_composite_ref': None,
            'Q_p_composite': None,
            'Phi_composite_0': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': None,
        },
        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': 'Tier 3 — review paper, no primary data.',
        },
        'Other common parameters': {
            'pressure': [],
            'reference': 'Ronnebro et al. 2022. Molecules 27, 6528. DOI: 10.3390/molecules27196528.',
            'temp_range': [],
            'pressure_range': [],
            'metal_thickness': [],
            'gas': 'H2/D2/T2',
            'Notes': (
                'Review only. Key classifications: oxide barriers Al2O3, Cr2O3, Y2O3, Er2O3. '
                'Permeation Pm = D0*exp(-E/RT) * S0*exp(-deltaH/RT). '
                'Barrier requirements: PRF, mechanical strength, radiation resistance, '
                'corrosion resistance, fabrication feasibility. '
                'Irradiation causes amorphous-to-crystalline transformation in Al2O3 -> grain growth. '
                'Use as background reference only — no numerical data.'
            ),
        },
    },

}

"""
OXIDE TRANSPORT PROPERTIES DICTIONARY
======================================
Entries extracted from three source papers:

[A] Chen CF, Yu HB, Zheng SQ.
    "First-principles study of hydrogen diffusion mechanism in Cr2O3"
    Science China Technological Sciences, 54(1), 88-94 (2011).
    DOI: 10.1007/s11431-010-4112-3
    Content: DFT/GGA-PW91 + molecular dynamics calculation of H atom
    diffusion in bulk alpha-Cr2O3 crystal. Provides D_0 and E_D for H.

[B] Hollenberg GW, Simonen EP, Kalinin G, Terlain A.
    "Tritium/hydrogen barrier development"
    Fusion Engineering and Design, 28, 190-208 (1995).
    DOI: 10.1016/0920-3796(95)90039-X
    Content: Review of hydrogen permeation barriers for fusion reactors.
    Provides PRF values for Cr2O3, Al2O3, TiC/TiN on structural steels,
    and activation energy for H permeation through Al2O3.
    Source of Al2O3 data: Roberts et al. [Ref 22 in paper]:
      R.M. Roberts, T.S. Elleman, H. Palour III, K. Verghese,
      "Hydrogen permeability of sintered aluminum oxide",
      J. Am. Ceram. Soc., 62 (1979) 495-499.

[C] Medasani BK, Sushko ML, Rosso KM, Schreiber DK, Bruemmer SM.
    "Temperature Dependence of Self-Diffusion in Cr2O3 from First Principles"
    J. Physical Chemistry C, 123(36), 22139-22150 (2019).
    DOI: 10.1021/acs.jpcc.9b03218
    Content: HSE hybrid DFT self-diffusion of Cr and O IONS in Cr2O3
    (NOT H atom diffusion). Relevant to oxide growth kinetics only.
    NOT used for H permeation Arrhenius parameters.

DATA QUALITY TIERS (inherited from OXIDES dict convention):
    Tier 1: Phi_composite, Phi_bare, oxide_thickness all known
    Tier 2: oxide_thickness known, Phi_bare from metals dict
    Tier 3: PRF only or activation energy only — no full back-calculation
    Tier DFT: Standalone DFT bulk crystal calculation (no composite data)

UNIT CONVENTIONS (SI throughout):
    D      : m²/s
    Phi    : mol/m/s/Pa^0.5
    E_D    : J/mol
    T      : K
    d_ox   : m
"""

OXIDES_FROM_PAPERS = {

    # =========================================================================
    # ENTRY 1: H diffusion in bulk alpha-Cr2O3 (DFT)
    # Source: Paper [A] — Chen et al. 2011
    # =========================================================================

    'alpha_Cr2O3_H_diffusion_bulk_Chen2011': {
        'Oxide identity': {
            'oxide_type': 'alpha-Cr2O3 (corundum, R-3C space group)',
            'substrate': 'None — standalone bulk crystal DFT calculation',
            'crystal_structure': {
                'space_group': 'R-3C trigonal',
                'a': 4.9598e-10,           # m (calculated, matches expt 4.953 Å, error 0.14%)
                'c': 13.5894e-10,           # m (calculated, matches expt 13.578 Å, error 0.08%)
                'Cr_per_cell': 12,
                'O_per_cell': 18,
            },
            'formation_condition': (
                'Perfect stoichiometric bulk crystal — DFT geometry optimised. '
                'No grain boundaries, vacancies, or defects. '
                'H stable position: bilateral positions of center of unoccupied '
                'O octahedral interstice (NOT at center, unlike alpha-Al2O3).'
            ),
            'oxide_thickness': None,        # standalone bulk — no coating
            'thickness_known': False,
            'data_tier': 'DFT',
        },

        'Diffusion parameters (H atom)': {
            # From Section 3.5 and 3.6 of Chen et al. 2011
            # Two-step diffusion A->B->C through nearest-neighbour octahedral interstices:
            #   Step A->B: Ea = 0.73 eV (rate-limiting, H crosses O atomic plane)
            #   Step B->C: Ea = 0.68 eV (H moves between stable positions within interstice)
            #   Total transition: Ea = 0.73 eV (rate-limiting step dominates)
            'E_D': 70434,                   # J/mol  (0.73 eV × 96485 J/mol/eV)
            'E_D_eV': 0.73,                 # eV — directly from paper Fig 5(a)
            'E_D_step2_eV': 0.68,           # eV — step B->C, paper Fig 5(b) (not rate-limiting)
            'D_0': 1.78e-9,                 # m²/s  (converted from 1.78×10⁻⁵ cm²/s)
            'attempt_frequency_nu': 3.74e10, # s⁻¹ — from first-principles MD at 1400 K
            'jump_distance_l': 3.63e-10,    # m  (3.63 Å — A to C transition distance)
            'method': 'DFT GGA/PW91 (CASTEP) + LST/QST transition state search + ab initio MD',
            'diffusion_path': (
                'H hops between nearest-neighbour unoccupied O octahedral interstices. '
                'Transition state is the O atomic plane between interstices. '
                'Two stable positions per interstice (bilateral, symmetric about center). '
                'Preferred path: along [001] direction through O layer.'
            ),
        },

        'Computed D(T) values — Table 3 of Chen et al. 2011': {
            # Units converted from cm²/s to m²/s (×1e-4)
            273:  2.34e-22,   # m²/s at   0°C
            373:  6.65e-19,   # m²/s at 100°C
            473:  6.55e-17,   # m²/s at 200°C
            773:  5.03e-14,   # m²/s at 500°C
            1273: 3.08e-12,   # m²/s at 1000°C
        },

        'Comparison values from Table 3 of Chen et al. 2011': {
            # D at 0°C (273 K) for context — all converted to m²/s
            'alpha_Al2O3_273K': 3.04e-30,   # m²/s — Belonoshko et al. (ref [13] in paper)
            'steel_16MnR_273K': 2.11e-10,   # m²/s
            'steel_20G_273K':   6.49e-10,   # m²/s
            'steel_Q235_273K':  3.75e-11,   # m²/s
            'steel_HR1_SS_273K': 6.76e-17,  # m²/s — HR-1 stainless
            'note': (
                'Cr2O3 D is ~8 orders of magnitude lower than structural steels '
                'but ~8 orders of magnitude higher than alpha-Al2O3. '
                'Both oxides are effective H barriers relative to base metals.'
            ),
        },

        'Permeation parameters': {
            'Phi_0': None,      # not calculated — no solubility data in this paper
            'Q_p': None,        # not calculable — only D_0 and E_D available
            'Phi_ref': None,
            'T_ref': None,
            'Ks_0': None,       # H solubility in Cr2O3 not determined
            'H_s': None,
        },

        'Other common parameters': {
            'reference': (
                'Chen CF, Yu HB, Zheng SQ. '
                '"First-principles study of hydrogen diffusion mechanism in Cr2O3." '
                'Sci. China Tech. Sci. 54(1), 88-94 (2011). '
                'DOI: 10.1007/s11431-010-4112-3'
            ),
            'temp_range_validity': [273, 1273],  # K — range of Table 3 values
            'gas': 'H (atomic)',
            'pressure_range': None,
            'Notes': (
                '[IMPORTANT LIMITATIONS] '
                '(1) DFT calculation on PERFECT bulk crystal — no grain boundaries, '
                'vacancies, or defects. Real oxide coatings will have higher effective '
                'E_D due to trapping at grain boundaries and structural defects. '
                '(2) Only D_0 and E_D computed — Ks and Phi NOT available from this paper. '
                'Cannot construct full Arrhenius permeability without separate solubility data. '
                '(3) GGA/PW91 functional — known to underestimate activation barriers in '
                'transition metal oxides. HSE-corrected values would likely be higher. '
                '(4) H most stable position in Cr2O3 is bilateral (off-center in octahedral '
                'interstice) — distinct from Al2O3 where H sits at interstice center. '
                '(5) This E_D (70.4 kJ/mol) is much lower than Al2O3 (318 kJ/g-atom '
                'per Hollenberg 1995), consistent with Cr2O3 being a weaker H barrier.'
            ),
        },
    },


    # =========================================================================
    # ENTRY 2: Cr2O3 coating on SS316 — PRF from irradiation test (TREXMAN)
    # Source: Paper [B] — Hollenberg et al. 1995, Table 5
    # =========================================================================

    'Cr2O3_on_SS316_TREXMAN_Hollenberg1995': {
        'Oxide identity': {
            'oxide_type': 'Cr2O3 (chemically densified coating via SiO2 slurry + CrO3 + firing at 450°C)',
            'substrate': 'SS316 austenitic stainless steel',
            'formation_condition': (
                'Terai et al. method [Ref 14 in Hollenberg 1995]: '
                'SiO2 slurry deposition followed by chemical densification with chromia, '
                'fired at 450°C. Reference: T. Terai et al., J. Nucl. Mater. (in press at 1995).'
            ),
            'oxide_thickness': None,        # not stated in Hollenberg Table 5
            'thickness_known': False,
            'data_tier': 3,
        },

        'TREXMAN in-reactor test results — Table 5 of Hollenberg 1995': {
            # Two configurations tested at 600°C in YAYOI reactor, Pb-17Li source
            'config_A': {
                'barrier_system': 'Cr2O3/SS316',
                'barrier_orientation': 'Cr2O3 upstream (facing Pb-17Li)',
                'T_test': 873,               # K (600°C)
                'PRF': 10,
                'PRF_note': (
                    'POOR result — Cr2O3 in contact with molten Pb-17Li. '
                    'Chemical reduction of chromia by liquid metal anticipated '
                    '(Tanabe et al. noted Cr2O3 unstable in reducing Pb-17Li environment). '
                    'Hollenberg explicitly states chromia next to molten alloy '
                    'leads to ineffective barrier.'
                ),
            },
            'config_B': {
                'barrier_system': 'SS316/Cr2O3',
                'barrier_orientation': 'Cr2O3 downstream (exterior, no liquid metal contact)',
                'T_test': 873,               # K (600°C)
                'PRF': 100,
                'PRF_note': (
                    'IMPROVED result — Cr2O3 on exterior with no Pb-17Li contact. '
                    'Effective PRF ~100 achieved. Confirms Cr2O3 is chemically '
                    'unstable in Pb-17Li but functions as barrier when protected.'
                ),
            },
        },

        'As-reported parameters': {
            'T_ref': 873,                   # K (600°C test temperature)
            'PRF_upstream': 10,
            'PRF_downstream': 100,
            'Phi_composite_ref': None,      # not tabulated
            'Q_p_composite': None,
            'Phi_bare_ref': None,
            'metal_baseline_key': None,
        },

        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'Phi_oxide_0': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': (
                'Tier 3 — PRF only, no Phi_composite or oxide thickness given. '
                'Back-calculation not possible from available data.'
            ),
        },

        'Other common parameters': {
            'reference': (
                'Hollenberg GW, Simonen EP, Kalinin G, Terlain A. '
                '"Tritium/hydrogen barrier development." '
                'Fusion Engineering and Design 28, 190-208 (1995). '
                'DOI: 10.1016/0920-3796(95)90039-X. '
                'Table 5 (TREXMAN irradiation test). '
                'Primary Cr2O3 coating data from: T. Terai et al., '
                'J. Nucl. Mater. (to be published at time of Hollenberg 1995).'
            ),
            'temp_range': [873],            # K — single test temperature
            'pressure_range': None,
            'tritium_source': 'Pb-17Li eutectic',
            'reactor': 'YAYOI',
            'gas': 'T2 (tritium)',
            'Notes': (
                'In-reactor test result. Hollenberg Table 1 also lists '
                'Cr (or Cr2O3) on SS316 with PRF ~10 in lab experiments [Ref 14]. '
                'Cr2O3 is explicitly noted as LESS THERMODYNAMICALLY STABLE than '
                'Al2O3 and TiO2 because its oxide formation free energy is lower '
                '(see Fig 1 of Hollenberg 1995 — Cr2O3 + 3H2 = 2Cr + 3H2O '
                'line is above FeO+H2 line, meaning Cr2O3 is reduced at lower '
                'H2O/H2 ratios than iron oxides). '
                'Hollenberg concludes Cr2O3 is NOT recommended for Pb-17Li '
                'upstream contact; Al2O3 and TiC/TiN are preferred.'
            ),
        },
    },


    # =========================================================================
    # ENTRY 3: Al2O3 (sintered) — activation energy and permeability trend
    # Source: Paper [B] — Hollenberg et al. 1995 (citing Roberts et al. 1979)
    # =========================================================================

    'Al2O3_sintered_Roberts1979_via_Hollenberg1995': {
        'Oxide identity': {
            'oxide_type': 'Al2O3 (sintered polycrystalline alumina)',
            'substrate': 'None — standalone sintered alumina membrane',
            'formation_condition': (
                'Sintered Al2O3 ceramic specimens. '
                'Primary source: Roberts RM, Elleman TS, Palour H III, Verghese K. '
                '"Hydrogen permeability of sintered aluminum oxide." '
                'J. Am. Ceram. Soc. 62 (1979) 495-499. [Ref 22 in Hollenberg 1995]'
            ),
            'oxide_thickness': None,
            'thickness_known': False,
            'data_tier': 3,
        },

        'As-reported parameters — from Hollenberg 1995 text and Fig 2': {
            'T_ref': None,
            'Q_p': 318000,                  # J/mol (318 kJ/g-atom, stated p.200 of Hollenberg)
            'Q_p_source': (
                'Explicitly stated in Hollenberg 1995 p.200: '
                '"The activation energy for permeation through alumina at higher '
                'temperatures is 318 kJ g-atom-1" — attributed to Roberts et al. 1979.'
            ),
            'Phi_0': None,                  # not extracted from Fig 2 — no tabulated value
            'Phi_composite_ref': None,
            'PRF': None,                    # standalone membrane — no composite system
            'metal_baseline_key': None,
        },

        'Permeability trend from Fig 2 of Hollenberg 1995': {
            # Read from Fig 2 (Arrhenius plot, 700-300°C range)
            # Units in figure: moles/cm²/s at 101 kPa
            # Approximate values read from graph:
            'Phi_at_700C_approx': 1e-20,    # mol/cm²/s — from Fig 2 (very approximate)
            'Phi_at_500C_approx': 1e-25,    # mol/cm²/s — from Fig 2 (very approximate)
            'Phi_at_300C_approx': 1e-30,    # mol/cm²/s — extrapolated from Fig 2
            'read_method': (
                'Approximate visual read from Fig 2 of Hollenberg 1995. '
                'These are AREA-NORMALIZED fluxes (mol/cm²/s), NOT permeability '
                '(mol/m/s/Pa^0.5). Thickness of the Roberts et al. specimen not '
                'stated in Hollenberg — conversion to permeability not possible '
                'without original Roberts et al. 1979 paper.'
            ),
            'Hollenberg_note': (
                'Hollenberg explicitly states: "It is accepted that the extrapolation '
                'of the Al2O3 data to lower temperatures may lead to a significant '
                'underestimation of the permeability, because of the change in mechanism '
                'from bulk to either grain boundary or defect diffusion, which could '
                'lead to an increase of several orders of magnitude." '
                'Bulk Al2O3 permeability is so low that even with this uncertainty, '
                'a thin coating provides sufficient barrier for PRF > 10,000.'
            ),
        },

        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': 318000,            # J/mol — only value available
            'Phi_oxide_0': None,            # not extractable without thickness
            'D_ref': None,
            'E_D': None,                    # D_0 and E_D not separated in Hollenberg
            'D_0': None,
            'Ks_ref': None,
            'K_s0': None,
            'H_s': None,
            'back_calc_method': (
                'Tier 3 — Q_p = 318 kJ/mol extracted from Hollenberg 1995 text (p.200). '
                'Phi_0 not available without Roberts et al. 1979 original paper. '
                'Fig 2 of Hollenberg gives approximate flux values but specimen thickness '
                'not given, preventing conversion to permeability. '
                'To get full Arrhenius parameters, obtain: '
                'Roberts RM et al., J. Am. Ceram. Soc. 62 (1979) 495-499.'
            ),
        },

        'Other common parameters': {
            'reference': (
                'Hollenberg GW et al. Fusion Eng. Design 28, 190-208 (1995). '
                'DOI: 10.1016/0920-3796(95)90039-X. '
                'Al2O3 data originally from: '
                'Roberts RM, Elleman TS, Palour H III, Verghese K. '
                '"Hydrogen permeability of sintered aluminum oxide." '
                'J. Am. Ceram. Soc. 62, 495-499 (1979). [Ref 22 in Hollenberg]'
            ),
            'temp_range': [573, 973],       # K — approximate range of Fig 2 data
            'gas': 'H2',
            'Notes': (
                'Al2O3 is identified in Hollenberg as thermodynamically superior to Cr2O3 '
                'because Al has higher free energy of oxide formation (Fig 1). '
                'Al2O3 and TiO2 can re-form even at low O2 partial pressures. '
                'Practical barrier coatings using aluminide/Al2O3 achieve PRF = '
                '1000 to >100,000 in laboratory (Table 1), far exceeding Cr2O3 (~10). '
                'In-reactor PRF drops to 3-150 due to radiation-induced diffusion (REID) '
                'and radiation-electric-field-induced diffusion (REID) effects (Table 5). '
                'Hollenberg: activation energy for permeation through SS316 ~65 kJ/g-atom; '
                'through Al2O3 ~318 kJ/g-atom — factor ~5x higher explains barrier mechanism.'
            ),
        },
    },


    # =========================================================================
    # ENTRY 4: Al2O3 aluminide coatings on SS316 — composite system
    # Source: Paper [B] — Hollenberg et al. 1995, Fig 7, Table 4
    # Substrate data from: Forcey et al. [Refs 2,3]; Gilbert et al. [Ref 5];
    #                      McGuire [Ref 1]; Van Deventer & Maroni [Ref 4]
    # =========================================================================

    'Al2O3_aluminide_on_SS316_composite_Hollenberg1995': {
        'Oxide identity': {
            'oxide_type': 'Aluminide intermetallic / Al2O3 surface oxide',
            'substrate': 'SS316 austenitic stainless steel',
            'formation_condition': (
                'Packed-bed cementation (pack aluminizing): substrate immersed in '
                'Al/Al2O3/NH4Cl powder at 850-950°C. Forms FeAl, FeAl2, Fe2Al5, FeAl3 '
                'intermetallics + outer Al2O3 surface oxide on exposure to air/moisture. '
                'Coating thicknesses: 17-330 µm (see Table 4 of Hollenberg 1995). '
                'Al/M (aluminum to metal) ratio at surface: 0.42 to >4 depending on process.'
            ),
            'oxide_thickness': None,        # varies — see Table 4
            'thickness_known': False,
            'data_tier': 3,
        },

        'As-reported parameters (composite system)': {
            'T_ref': None,
            'PRF_lab_range': [1000, 1000000],   # dimensionless — from Hollenberg Table 1 and text
            'PRF_inreactor_range': [3, 150],     # from Table 5 (LIBRETTO, Loop-1, WC-1 tests)
            'Q_p_composite': None,               # see note below on activation energy
            'Q_p_note': (
                'Forcey et al. [Refs 2,3] and Gilbert et al. [Ref 5] observe that '
                'the activation energy for permeation through aluminized SS316 '
                'remains approximately the SAME as bare SS316 (~65 kJ/g-atom), '
                'NOT the Al2O3 bulk value (318 kJ/g-atom). '
                'Hollenberg interprets this as evidence that permeation is controlled '
                'by diffusion through the BASE METAL at defect sites (area-defect model), '
                'NOT by diffusion through the Al2O3 layer itself. '
                'This means the observed PRF reflects DEFECT DENSITY in the coating, '
                'not the intrinsic Al2O3 permeability.'
            ),
            'pressure_exponent': {
                'Forcey_high_P': 0.46,      # Sievert-like (p^0.5) at 1000-10000 Pa
                'McGuire_low_P': 1.1,       # linear (p^1) at 200-1000 Pa
                'Perujo_high_P': 0.5,       # Sievert at >20000 Pa
                'Perujo_low_P': 1.0,        # linear at <20000 Pa
                'interpretation': (
                    'p^0.5 at high pressure = area-defect / bulk diffusion control. '
                    'p^1 at low pressure = surface-limited / recombination control. '
                    'Transition consistent with surface recombination limiting at '
                    'low flux conditions (few H surface sites available).'
                ),
            },
        },

        'Coating parameter summary from Table 4 of Hollenberg 1995': {
            # Reference: Table 4
            'McGuire_1980': {
                'aluminide_thickness_um': 40,
                'metal_thickness_cm': 0.04,
                'substrate': '304SS',
                'Al_M_ratio': '>2',
            },
            'Forcey_tubes_1991': {
                'aluminide_thickness_um': '275-330',
                'metal_thickness_cm': 0.2,
                'substrate': '316L',
                'Al_M_ratio': 0.5,
                'coating_sides': 'both',
            },
            'Forcey_disks_1991': {
                'aluminide_thickness_um': 30,
                'metal_thickness_cm': 0.16,
                'substrate': '316L',
                'Al_M_ratio': 0.9,
                'coating_sides': 'both',
            },
            'Gilbert_1992': {
                'aluminide_thickness_um': '17-120',
                'metal_thickness_cm': 0.05,
                'substrate': 'SS316',
                'Al_M_ratio': '0.42-4',
                'coating_sides': 'one',
            },
            'VanDeventer_1980': {
                'aluminide_thickness_um': 4,
                'metal_thickness_cm': 0.13,
                'substrate': 'Pure Al sputtered on SS',
                'Al_M_ratio': None,
            },
        },

        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'back_calc_method': (
                'Tier 3 — composite PRF data only. '
                'Permeation controlled by base metal at defect sites (area-defect model), '
                'NOT by Al2O3 layer. Extracting standalone Al2O3 Phi_oxide '
                'from these data would be physically incorrect.'
            ),
        },

        'Other common parameters': {
            'reference': (
                'Hollenberg GW et al. Fusion Eng. Design 28, 190-208 (1995). '
                'DOI: 10.1016/0920-3796(95)90039-X. '
                'Fig 7, Tables 1, 4, 5. '
                'Composite permeation data from: '
                'McGuire JG (1980); '
                'Forcey KS, Ross DK, Wu CH, J. Nucl. Mater. 182 (1991) 36; '
                'Forcey KS et al., J. Nucl. Mater. 161 (1989) 108-116; '
                'Gilbert ER et al., Fusion Technol. 21 (1992) 739-744; '
                'Van Deventer EH et al., J. Nucl. Mater. 88 (1980) 168-173.'
            ),
            'temp_range': [573, 923],       # K — Fig 7 data range 300-650°C
            'pressure_range': [200, 101325], # Pa
            'gas': 'H2 / D2 / T2',
            'Notes': (
                'KEY FINDING from Hollenberg: Activation energy of aluminized coating '
                'matches BASE METAL (~65 kJ/mol), not Al2O3 (~318 kJ/mol). '
                'This supports area-defect mechanism: H permeates through coating '
                'defects (pinholes, cracks) and rate is controlled by diffusion through '
                'the underlying metal at those defect sites, not by bulk oxide diffusion. '
                'PRF = 1000-1,000,000 in lab. PRF = 3-150 in-reactor (Table 5). '
                'Gap attributed to: REID, corrosion by Pb-17Li, coating scale-up defects. '
                'Eddy current NDE (Gilbert et al.) essential for quality control.'
            ),
        },
    },


    # =========================================================================
    # ENTRY 5: TiC / TiN coatings on SS316 — composite system
    # Source: Paper [B] — Hollenberg et al. 1995, Fig 13, Table 1
    # Primary data: Forcey et al. [Ref 11]; Shan et al. [Ref 10]
    # =========================================================================

    'TiC_TiN_on_SS316_composite_Hollenberg1995': {
        'Oxide identity': {
            'oxide_type': 'TiC / TiN ceramic coatings (carbide/nitride, not oxide)',
            'substrate': 'SS316 austenitic stainless steel',
            'formation_condition': (
                'CVD (chemical vapor deposition) or PVD (plasma vapor deposition). '
                'Forcey et al. [Ref 11]: CVD TiN(1.5µm)/TiC(1.5µm) bilayer, '
                'and Al2O3(4µm) over TiC/TiN. '
                'Shan et al. [Ref 10]: PVD Ti(0.01µm)/TiC(2.5µm) or '
                'Ti(0.01µm)/TiN(1µm)/TiC/TiN composite (2.5µm total).'
            ),
            'oxide_thickness': {
                'Forcey_TiN_TiC': [1.5e-6, 1.5e-6],   # m each layer
                'Shan_TiC': 2.5e-6,                     # m
                'Shan_TiN_TiC_composite': 2.5e-6,       # m total
            },
            'thickness_known': True,
            'data_tier': 3,
        },

        'As-reported parameters (composite system) — Fig 13 of Hollenberg 1995': {
            # Approximate permeability values read from Fig 13
            # Units in Fig 13: mol/cm²/s (area-normalized flux, NOT permeability)
            # Temperature range: ~550-830 K (280-560°C)
            'Forcey_TiN_TiC': {
                'PRF_vs_bare_SS316': 10,            # ~1 order of magnitude reduction
                'effective_permeability_at_550K_approx': 1e-14,  # mol/cm²/s from Fig 13
            },
            'Shan_TiC': {
                'PRF_vs_bare_SS316': 1e5,           # ~5-6 orders of magnitude
                'effective_permeability_at_625K_approx': 1e-17,  # mol/cm²/s from Fig 13
            },
            'Shan_TiN_TiC_composite': {
                'PRF_vs_bare_SS316': 1e6,           # ~6 orders of magnitude
                'effective_permeability_at_625K_approx': 1e-18,  # mol/cm²/s from Fig 13
            },
            'Hollenberg_interpretation': (
                'Vastly different PRFs (10 vs 1,000,000) between Forcey and Shan '
                'attributed by Hollenberg to DEFECT DENSITY differences, not '
                'intrinsic material property differences. Shan et al. fabricated '
                'a coating with reduced defect size and population. '
                'Hollenberg notes Shan results show evidence of chemical reactivity '
                'and hydrocarbon formation during testing (not discussed by Shan).'
            ),
        },

        'Back-calculated standalone oxide parameters': {
            'Phi_oxide_ref': None,
            'Q_p_oxide': None,
            'D_ref': None,
            'E_D': None,
            'D_0': None,
            'Ks_ref': None,
            'back_calc_method': (
                'Tier 3 — PRF/flux data only from Fig 13. '
                'Specimen thicknesses not consistently reported. '
                'Back-calculation not possible.'
            ),
        },

        'Other common parameters': {
            'reference': (
                'Hollenberg GW et al. Fusion Eng. Design 28, 190-208 (1995). '
                'DOI: 10.1016/0920-3796(95)90039-X. '
                'Fig 13, Table 1. '
                'Forcey data: Forcey KS, Perujo A, Reiter F, Lolli Ceroni PL, '
                '"The formation of tritium permeation barriers by CVD," '
                'J. Nucl. Mater. 200 (1993) 417-420. '
                'Shan data: Shan C et al., '
                '"The behaviour of diffusion and permeation of tritium through '
                '316L stainless steel with coating of TiC and TiN/TiC," '
                'J. Nucl. Mater. 191-194A (1992) 221-225.'
            ),
            'temp_range': [553, 833],       # K — approximate from Fig 13 x-axis
            'gas': 'D2 (Forcey) / T2 (Shan)',
            'Notes': (
                'Hollenberg Table 1 entry: TiC, TiN, TiO2 on SS316/MANET/TZM/Ti '
                'achieves PRF = less than 10 to above 10,000. '
                'TiC/TiN are CVD-deposited non-oxide ceramics — included here '
                'because they appear in Hollenberg alongside oxide barriers and '
                'are directly compared in Fig 13.'
            ),
        },
    },


    # =========================================================================
    # ENTRY 6: Cr2O3 self-diffusion (Cr and O ions) — NOT H diffusion
    # Source: Paper [C] — Medasani et al. 2019
    # =========================================================================

    'Cr2O3_self_diffusion_Cr_O_ions_Medasani2019': {
        'Oxide identity': {
            'oxide_type': 'alpha-Cr2O3 (R-3c corundum) — bulk crystal',
            'substrate': 'None — standalone bulk DFT calculation',
            'crystal_parameters': {
                'a': 4.925e-10,             # m (HSE optimised, Table 2)
                'c': 13.54e-10,             # m (HSE optimised, Table 2)
                'magnetic_moment_Cr': 2.81, # µ_B per Cr (HSE, Table 2)
                'bandgap': 3.38,            # eV (HSE with screening=0.6, Table 1+2)
                'expt_bandgap': 3.4,        # eV
            },
            'formation_condition': (
                'Perfect stoichiometric bulk crystal. '
                'DFT method: HSE hybrid functional (25% HF exchange, screening=0.6), '
                'VASP code, PAW pseudopotentials, AFM spin ordering. '
                'Supercell: 2×2×1 for defect calculations.'
            ),
            'thickness_known': False,
            'data_tier': 'DFT',
        },

        'CRITICAL NOTE — SCOPE OF THIS PAPER': {
            'subject': 'Cr and O ION self-diffusion (oxide growth kinetics)',
            'NOT_subject': 'H atom diffusion through Cr2O3 for permeation modelling',
            'relevance_to_H_permeation': (
                'INDIRECT ONLY. This paper cannot provide D_0 or E_D for H transport. '
                'It is useful for: '
                '(a) understanding Cr2O3 passive film growth rates on alloy surfaces, '
                '(b) estimating oxide layer thickness evolution over time, '
                '(c) understanding point defect chemistry (V_Cr, V_O, Cr_i, O_i) '
                'which controls oxide microstructure and thus affects real '
                'H permeation barriers (grain boundaries, defect density).'
            ),
        },

        'Self-diffusion coefficients from Fig 9 of Medasani 2019': {
            # All at vacuum/oxide interface (representative condition for passive films)
            # Read from Fig 9 at approximately 1000 K and 1 atm O2
            # Units: cm²/s in original figure
            'conditions_for_H_permeation_relevance': {
                'pO2': '1 atm (most relevant for passive film in air)',
                'T_approx': 1000,           # K
                'D_Cr_approx': 1e-22,       # cm²/s ~ 1e-26 m²/s (from Fig 9a at 1000K)
                'D_O_approx': 1e-30,        # cm²/s ~ 1e-34 m²/s (from Fig 9a at 1000K)
                'dominant_species': 'Cr (higher mobility at 1 atm, T < 1100 K)',
                'note': (
                    'Read from Fig 9a (log D vs 10^4/T). Values are approximate. '
                    'At 1 atm O2 and T < 1100K: Cr diffuses faster than O. '
                    'At low pO2 or high T: O diffusion dominates. '
                    'Experimental values from Sabioni et al. (1992): '
                    'D_Cr = D_O ~ 2-8 × 10^-18 cm²/s at 1300-1100°C (high vacuum). '
                    'Paper reproduces these at high T within factor ~2-5.'
                ),
            },
            'dominant_defects_by_condition': {
                '300K_1atm': 'V_Cr (Cr vacancies dominant)',
                '800K_1atm': 'V_Cr (Cr vacancies dominant)',
                '1200K_low_pO2': 'V_O (O vacancies become dominant)',
                'Cr_Cr2O3_interface': 'V_O dominant at all T (Fig 9d)',
            },
        },

        'Defect formation energies at O-rich conditions (Table 3)': {
            # EF at middle of bandgap, O-rich conditions
            # HSE functional values
            'V_O_0': 5.55,      # eV
            'V_O_1+': 5.62,     # eV
            'V_O_2+': 6.16,     # eV
            'V_Cr_0': 0.90,     # eV — lowest formation energy
            'V_Cr_1-': 0.93,    # eV
            'V_Cr_2-': 1.30,    # eV
            'V_Cr_3-': 2.00,    # eV
            'O_i_0': 2.94,      # eV
            'Cr_i_0': 11.65,    # eV
            'note': (
                'V_Cr (Cr vacancies) have lowest formation energy at O-rich conditions. '
                'This makes Cr2O3 intrinsically p-type at standard conditions '
                '(Fermi level below midgap). V_O forms less readily but dominates '
                'at reducing/low-pO2 conditions. '
                'These defect densities control grain boundary character and '
                'therefore affect real H permeation barriers.'
            ),
        },

        'Other common parameters': {
            'reference': (
                'Medasani BK, Sushko ML, Rosso KM, Schreiber DK, Bruemmer SM. '
                '"Temperature Dependence of Self-Diffusion in Cr2O3 from First Principles." '
                'J. Phys. Chem. C 123(36), 22139-22150 (2019). '
                'DOI: 10.1021/acs.jpcc.9b03218. '
                'Data/scripts: https://github.com/mbkumar/Cr2O3_diffusion'
            ),
            'temp_range': [300, 1600],      # K — range of Fig 9 plots
            'gas': None,
            'Notes': (
                'DO NOT use D_Cr or D_O from this paper as H diffusion coefficients. '
                'Cr and O ion diffusion (10^-22 to 10^-34 m²/s range) is completely '
                'different physics from H atom interstitial diffusion '
                '(2×10^-22 to 3×10^-12 m²/s range from Chen et al. 2011). '
                'Use this paper for: passive film growth rate estimates, '
                'oxide microstructure evolution modelling, '
                'understanding defect chemistry that governs H trapping in oxide. '
                'Downloaded via Northeastern University on March 16, 2026 '
                '(confirmed in PDF header).'
            ),
        },
    },

}


# =============================================================================
# SUMMARY TABLE of extractable Arrhenius parameters from all three papers
# =============================================================================

OXIDE_ARRHENIUS_SUMMARY = {
    'alpha_Cr2O3_H_atom_bulk': {
        'D_0': 1.78e-9,     # m²/s
        'E_D': 70434,       # J/mol (0.73 eV)
        'Phi_0': None,      # not available
        'Q_p': None,        # not available
        'source': 'Chen et al. 2011 [A]',
        'notes': 'DFT bulk perfect crystal. Real coatings will have higher effective E_D.',
    },
    'Al2O3_H_permeation_Q_p_only': {
        'D_0': None,        # not available in Hollenberg
        'E_D': None,        # not separated from Q_p in Hollenberg
        'Phi_0': None,      # not tabulated
        'Q_p': 318000,      # J/mol — from Hollenberg 1995 p.200 citing Roberts 1979
        'source': 'Hollenberg et al. 1995 [B] citing Roberts et al. 1979',
        'notes': (
            'Activation energy for bulk permeation through sintered Al2O3 at high T. '
            'Full Arrhenius Phi_0 requires Roberts et al. 1979 original paper. '
            'At lower T, grain boundary/defect diffusion may increase effective permeability.'
        ),
    },
    'Cr2O3_SS316_composite_PRF': {
        'D_0': None,
        'E_D': None,
        'Phi_0': None,
        'Q_p': None,
        'PRF_lab': 10,      # from Hollenberg Table 1
        'PRF_TREXMAN_downstream': 100,  # from Hollenberg Table 5
        'source': 'Hollenberg et al. 1995 [B]',
        'notes': 'PRF only. Composite/in-reactor data. No standalone Phi_oxide extractable.',
    },
    'Cr2O3_self_diffusion_ions': {
        'D_Cr_at_1000K_approx': 1e-26,  # m²/s at 1 atm (from Fig 9a Medasani 2019)
        'D_O_at_1000K_approx': 1e-34,   # m²/s at 1 atm
        'source': 'Medasani et al. 2019 [C]',
        'notes': 'Cr and O ION diffusion — NOT H atom diffusion. Use for oxide growth only.',
    },
}
