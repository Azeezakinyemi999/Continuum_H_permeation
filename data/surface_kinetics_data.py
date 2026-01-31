"""
Surface Kinetics Parameters for Level 6: Surface Chemistry

This module contains material parameters for modeling hydrogen surface kinetics
including dissociative adsorption and recombinative desorption. The data focuses
on metals relevant to high-temperature hydrogen service.

Theory:
-------
At metal surfaces, hydrogen undergoes:
1. Dissociative adsorption: H₂(g) + 2S → 2H(ads)
   Rate: J_diss = k_diss × P × (1 - θ)²
   
2. Recombinative desorption: 2H(ads) → H₂(g) + 2S
   Rate: J_recomb = k_recomb × θ²

where:
- θ = surface coverage (fraction of occupied sites, 0 ≤ θ ≤ 1)
- S = empty surface site
- k_diss, k_recomb follow Arrhenius form: k = k₀ × exp(-E_act/RT)

At equilibrium (Langmuir isotherm):
    θ/(1-θ) = √(K_eq × P)
    where K_eq = k_diss/k_recomb

Relationship to Sieverts' constant:
    K_s = √(k_diss/k_recomb) × N_surf/N_L
    (connects surface equilibrium to bulk solubility)

The surface site density N_surf determines maximum coverage:
    N_surf ~ 10¹⁵ atoms/cm² ~ 10¹⁹ m⁻² (typical close-packed metal)

Units:
------
- k_diss_0: m⁴/(mol·s) - dissociation pre-exponential
- k_recomb_0: m²/s - recombination pre-exponential  
- E_diss: J/mol - dissociation activation energy
- E_recomb: J/mol - recombination activation energy
- N_surf: m⁻² - surface site density (sites per unit area)

References:
-----------
1. Pick, M.A. & Sonnenberg, K. (1985). "A model for atomic hydrogen-metal 
   interactions - Application to recycling, recombination and permeation."
   J. Nucl. Mater. 131, 208-220. DOI:doi.org/10.1016/0022-3115(85)90459-3

2. Baskes, M.I. (1980). "A calculation of the surface recombination rate 
   constant for hydrogen isotopes on metals." J. Nucl. Mater. 92, 318-324.
   DOI: doi.org/10.1016/0022-3115(80)90117-8

3. Andrew, P.L. & Haasz, A.A. (1992). "Models for hydrogen permeation in 
   metals." J. Appl. Phys. 72, 2749-2757. DOI: 10.1063/1.351526

4. Causey, R.A. (2002). "Hydrogen isotope retention and recycling in fusion 
   reactor plasma-facing components." J. Nucl. Mater. 300, 91-117.
   DOI: 10.1016/S0022-3115(01)00732-2

5. Wampler, W.R. (1986). "Surface recombination of hydrogen on clean nickel."
   Appl. Phys. Lett. 48, 405-407. DOI: 10.1063/1.96521
"""

import numpy as np

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

R_GAS = 8.314  # J/(mol·K) - Universal gas constant

# =============================================================================
# SURFACE SITE DENSITIES
# =============================================================================

SURFACE_SITE_DENSITY = {
    'FCC_111': {
        'N_surf': 1.5e19,  # m⁻² (sites per unit area)
        'description': 'Close-packed FCC (111) surface',
        'reference': 'Crystallographic calculation'
    },
    'FCC_100': {
        'N_surf': 1.2e19,  # m⁻²
        'description': 'FCC (100) surface',
        'reference': 'Crystallographic calculation'
    },
    'BCC_110': {
        'N_surf': 1.4e19,  # m⁻²
        'description': 'Close-packed BCC (110) surface',
        'reference': 'Crystallographic calculation'
    },
    'polycrystalline': {
        'N_surf': 1.0e19,  # m⁻² (effective average)
        'description': 'Polycrystalline average',
        'reference': 'Typical value for mixed orientations'
    }
}

# =============================================================================
# SURFACE KINETICS PARAMETERS
# =============================================================================

SURFACE_KINETICS = {
    # =========================================================================
    # Incoloy 800 (Fe-Ni-Cr austenitic alloy)
    # Primary material for this project
    # =========================================================================
    'Incoloy800': {
        # Dissociative adsorption: H₂ + 2S → 2H(ads)
        # k_diss = k_diss_0 × exp(-E_diss/RT)
        'k_diss_0': 1.0e-2,     # m⁴/(mol·s)
        'E_diss': 15000.0,       # J/mol (~15 kJ/mol)
        
        # Recombinative desorption: 2H(ads) → H₂ + 2S
        # k_recomb = k_recomb_0 × exp(-E_recomb/RT)
        'k_recomb_0': 1.0e-7,   # m²/s
        'E_recomb': 80000.0,     # J/mol (~80 kJ/mol)
        
        # Surface properties
        'N_surf': 1.0e19,        # m⁻² (polycrystalline)
        'surface_type': 'polycrystalline',
        
        # Metadata
        'reference': 'Estimated from Ni data (Wampler 1986), adjusted for alloy',
        'temp_range': [600, 1100],  # K
        'notes': 'Ni-Fe-Cr alloy; Ni dominates surface behavior. Values estimated.',
        'uncertainty_factor': 5  # Parameters uncertain within 5×
    },
    
    # =========================================================================
    # Pure Iron (BCC α-Fe) - Well characterized reference
    # =========================================================================
    'Fe_alpha': {
        'k_diss_0': 5.0e-2,     # m⁴/(mol·s)
        'E_diss': 10000.0,       # J/mol (~10 kJ/mol) - small barrier on clean Fe
        
        'k_recomb_0': 5.0e-8,   # m²/s
        'E_recomb': 70000.0,     # J/mol (~70 kJ/mol)
        
        'N_surf': 1.4e19,        # m⁻² (BCC 110)
        'surface_type': 'BCC_110',
        
        'reference': 'Pick & Sonnenberg (1985), Baskes (1980)',
        'temp_range': [300, 900],  # K
        'notes': 'Clean iron surface; oxide-free. Well-characterized system.',
        'uncertainty_factor': 2
    },
    
    # =========================================================================
    # Pure Nickel (FCC) - Reference data
    # =========================================================================
    'Ni': {
        'k_diss_0': 2.0e-3,     # m⁴/(mol·s)
        'E_diss': 20000.0,       # J/mol (~20 kJ/mol)
        
        'k_recomb_0': 2.0e-8,   # m²/s
        'E_recomb': 85000.0,     # J/mol (~85 kJ/mol)
        
        'N_surf': 1.5e19,        # m⁻² (FCC 111)
        'surface_type': 'FCC_111',
        
        'reference': 'Wampler (1986), Baskes (1980)',
        'temp_range': [300, 1200],  # K
        'notes': 'Clean Ni surface. Higher barrier than Fe.',
        'uncertainty_factor': 2
    },
    
    # =========================================================================
    # SS316L (Austenitic stainless steel)
    # =========================================================================
    'SS316L': {
        'k_diss_0': 5.0e-3,     # m⁴/(mol·s)
        'E_diss': 18000.0,       # J/mol (~18 kJ/mol)
        
        'k_recomb_0': 5.0e-8,   # m²/s
        'E_recomb': 82000.0,     # J/mol (~82 kJ/mol)
        
        'N_surf': 1.0e19,        # m⁻² (polycrystalline)
        'surface_type': 'polycrystalline',
        
        'reference': 'Causey (2002), estimated from Fe/Ni data',
        'temp_range': [300, 1000],  # K
        'notes': 'Cr-Ni-Mo alloy; surface composition differs from bulk.',
        'uncertainty_factor': 3
    }
}

# =============================================================================
# STICKING COEFFICIENTS (for reference/future use)
# =============================================================================

STICKING_COEFFICIENTS = {
    # Sticking coefficient s = probability that incident H₂ adsorbs
    # s₀ is the zero-coverage limit; decreases as (1-θ)² with coverage
    
    'Fe_alpha': {
        's_0': 0.1,  # Clean surface, low T
        'reference': 'Pick & Sonnenberg (1985)'
    },
    'Ni': {
        's_0': 0.05,  # Clean surface
        'reference': 'Wampler (1986)'
    },
    'Incoloy800': {
        's_0': 0.03,  # Estimated
        'reference': 'Estimated from Ni'
    }
}


# =============================================================================
# GETTER FUNCTIONS
# =============================================================================

def get_surface_kinetics(material_name, temperature):
    """
    Get surface kinetics rate constants for a material at given temperature.
    
    Parameters
    ----------
    material_name : str
        Material name (e.g., 'Incoloy800', 'Ni', 'Fe_alpha', 'SS316L')
    temperature : float
        Temperature in K
        
    Returns
    -------
    dict with keys:
        - k_diss: float, dissociation rate constant at T
        - k_recomb: float, recombination rate constant at T
        - N_surf: float, surface site density (m⁻²)
        - K_eq: float, equilibrium constant k_diss/k_recomb
        - reference: str, data source
        - notes: str, additional information
        
    Raises
    ------
    ValueError
        If material not found in database
    """
    if material_name not in SURFACE_KINETICS:
        available = list(SURFACE_KINETICS.keys())
        raise ValueError(f"Material '{material_name}' not found. Available: {available}")
    
    params = SURFACE_KINETICS[material_name]
    
    # Calculate rate constants at temperature using Arrhenius form
    k_diss = params['k_diss_0'] * np.exp(-params['E_diss'] / (R_GAS * temperature))
    k_recomb = params['k_recomb_0'] * np.exp(-params['E_recomb'] / (R_GAS * temperature))
    
    # Equilibrium constant
    K_eq = k_diss / k_recomb if k_recomb > 0 else np.inf
    
    return {
        'k_diss': k_diss,
        'k_recomb': k_recomb,
        'K_eq': K_eq,
        'N_surf': params['N_surf'],
        'E_diss': params['E_diss'],
        'E_recomb': params['E_recomb'],
        'reference': params.get('reference', 'Unknown'),
        'notes': params.get('notes', ''),
        'uncertainty_factor': params.get('uncertainty_factor', 1)
    }


def list_available_materials():
    """Return list of materials with surface kinetics data."""
    return list(SURFACE_KINETICS.keys())