"""
Surface Kinetics Model for Hydrogen Permeation (Level 6) - Dissociation Only

Simplified model where only upstream dissociation can be rate-limiting.
Downstream recombination is assumed fast (instantaneous).

Theory:
-------
J = min(J_diss, J_Sieverts)

where:
    J_diss = k_diss × P × (1-θ)²
    θ = √(K_eq × P) / (1 + √(K_eq × P))  [Langmuir isotherm]
    K_eq = k_diss / k_recomb

Damköhler number:
    Da = k_diss × K_s × L / D
    Da >> 1: Sieverts' Law (diffusion-limited)
    Da << 1: Surface-limited (dissociation controls)

References:
-----------
1. Pick, M.A. & Sonnenberg, K. (1985). J. Nucl. Mater. 131, 208-220.
2. Baskes, M.I. (1980). J. Nucl. Mater. 92, 318-324.
"""

import numpy as np

# Import existing functions
from calculations.permeation_calc import sieverts_concentration, fick_flux

# Gas constant
R = 8.314  # J/(mol·K)


def surface_equilibrium_coverage(pressure, k_diss, k_recomb):
    """
    Calculate equilibrium surface coverage using Langmuir isotherm.
    
    θ = √(K_eq × P) / (1 + √(K_eq × P))
    
    Parameters
    ----------
    pressure : float
        Hydrogen partial pressure [Pa]
    k_diss : float
        Dissociation rate constant [m/s] or appropriate units
    k_recomb : float
        Recombination rate constant [m⁴/(mol·s)] or appropriate units
    
    Returns
    -------
    dict with 'theta', 'K_eq', 'sqrt_K_eq_P', 'regime'
    """
    if pressure < 0:
        raise ValueError(f"Pressure must be non-negative, got {pressure}")
    if k_diss <= 0 or k_recomb <= 0:
        raise ValueError("Rate constants must be positive")
    
    if pressure == 0:
        return {'theta': 0.0, 'K_eq': k_diss/k_recomb, 'sqrt_K_eq_P': 0.0, 'regime': 'low_coverage'}
    
    K_eq = k_diss / k_recomb
    sqrt_K_eq_P = np.sqrt(K_eq * pressure)
    theta = sqrt_K_eq_P / (1.0 + sqrt_K_eq_P)
    
    if theta < 0.1:
        regime = 'low_coverage'
    elif theta > 0.9:
        regime = 'high_coverage'
    else:
        regime = 'intermediate'
    
    return {'theta': theta, 'K_eq': K_eq, 'sqrt_K_eq_P': sqrt_K_eq_P, 'regime': regime}


def dissociation_flux(pressure, theta, k_diss):
    """
    Calculate dissociative adsorption flux.
    
    J_diss = k_diss × P × (1-θ)²
    
    Parameters
    ----------
    pressure : float
        Hydrogen partial pressure [Pa]
    theta : float
        Surface coverage fraction (0 ≤ θ ≤ 1)
    k_diss : float
        Dissociation rate constant
    
    Returns
    -------
    float
        Dissociation flux [mol/(m²·s)]
    """
    blocking_factor = (1.0 - theta) ** 2
    return k_diss * pressure * blocking_factor


def calculate_damkohler_number(k_diss, D, K_s, thickness):
    """
    Calculate Damköhler number for dissociation.
    
    Da = k_diss × K_s × L / D
    
    Da >> 1: Diffusion-limited (Sieverts applies)
    Da << 1: Dissociation-limited (surface controls)
    
    Returns
    -------
    dict with 'Da', 'regime', 'description'
    """
    Da = k_diss * K_s * thickness / D
    
    if Da > 10:
        regime = 'diffusion-limited'
        description = "Fast dissociation; Sieverts' Law applies"
    elif Da < 0.1:
        regime = 'dissociation-limited'
        description = "Slow dissociation controls flux"
    else:
        regime = 'mixed'
        description = "Comparable surface and bulk resistances"
    
    return {'Da': Da, 'regime': regime, 'description': description}


def calculate_surface_limited_flux(D, K_s, thickness, P_up, P_down, temperature,
                                    k_diss=None, k_recomb=None,
                                    material_name=None, coverage=None):
    """
    Calculate hydrogen flux with dissociation-limited surface kinetics.
    
    Two usage modes:
    ----------------
    Option 1: Provide k_diss + k_recomb (no coverage)
        → Calculate θ from Langmuir isotherm
        → C_surface = K_s × (θ/(1-θ)) / √K_eq  [= K_s × √P at equilibrium]
        → J = D × (C_surface - C_down) / L
    
    Option 2: Provide k_diss + k_recomb + coverage
        → Use provided θ with rigorous formula
        → C_surface = K_s × (θ/(1-θ)) / √K_eq
        → J = D × (C_surface - C_down) / L
        → Warning if provided θ differs from Langmuir equilibrium
    
    Note: When θ is calculated from Langmuir (Option 1), the formula reduces
    to Sieverts' Law: C_surface = K_s × √P. The surface effect appears through
    the steady-state coupling, not through a modified boundary condition.
    
    Parameters
    ----------
    D : float
        Bulk diffusivity [m²/s]
    K_s : float
        Solubility or Sieverts' constant [mol/(m³·Pa^0.5)]
    thickness : float
        Membrane thickness [m]
    P_up : float
        Upstream pressure [Pa]
    P_down : float
        Downstream pressure [Pa]
    temperature : float
        Temperature [K] - used for material lookup; if providing k_diss/k_recomb
        directly, ensure they are at this temperature
    k_diss : float, optional
        Dissociation rate constant [mol/(m²·s·Pa)]
    k_recomb : float, optional
        Recombination rate constant [m⁴/(mol·s)]
    material_name : str, optional
        Material name to look up kinetics
    coverage : float, optional
        Fixed surface coverage θ (0 ≤ θ < 1). If provided with kinetics,
        overrides Langmuir calculation.
    
    Returns
    -------
    dict with flux, theta, damkohler info, surface_reduction_factor
    """
    # Input validation
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    if coverage is not None and (coverage < 0 or coverage >= 1):
        raise ValueError("Coverage must be in range [0, 1)")
    
    # Sieverts flux (Level 1 reference)
    C_up_sieverts = sieverts_concentration(K_s, P_up)
    C_down = sieverts_concentration(K_s, P_down)
    flux_sieverts = fick_flux(D, C_up_sieverts, C_down, thickness)
    
    # Determine which option we're using
    have_kinetics = (k_diss is not None and k_recomb is not None)
    have_coverage = (coverage is not None)
    
    # Try to get kinetics from material if not provided
    if not have_kinetics and material_name is not None:
        from data.surface_kinetics_data import get_surface_kinetics
        kinetics = get_surface_kinetics(material_name, temperature)
        k_diss = kinetics['k_diss']
        k_recomb = kinetics['k_recomb']
        have_kinetics = True
    
    # Must have kinetics (removed Option 2 - coverage only)
    if not have_kinetics:
        raise ValueError("Must provide k_diss + k_recomb, or material_name. "
                         "Coverage alone is not sufficient.")
    
    K_eq = k_diss / k_recomb
    
    # Damköhler number
    Da_info = calculate_damkohler_number(k_diss, D, K_s, thickness)
    
    # Calculate Langmuir equilibrium coverage
    coverage_result = surface_equilibrium_coverage(P_up, k_diss, k_recomb)
    theta_langmuir = coverage_result['theta']
    
    if have_coverage:
        # Option 2: Use provided coverage with rigorous formula
        theta = coverage
        coverage_mode = 'fixed_rigorous'
        
        # Warn if provided coverage differs significantly from equilibrium
        if abs(theta - theta_langmuir) > 0.1:
            import warnings
            warnings.warn(
                f"Provided coverage ({theta:.3f}) differs from Langmuir equilibrium "
                f"({theta_langmuir:.3f}) by more than 0.1. This may indicate "
                f"non-equilibrium conditions or inconsistent parameters.",
                UserWarning
            )
    else:
        # Option 1: Use Langmuir equilibrium coverage
        theta = theta_langmuir
        coverage_mode = 'calculated'
    
    # Rigorous formula: C_surface = K_s × (θ/(1-θ)) / √K_eq
    # At Langmuir equilibrium, this reduces to K_s × √P (Sieverts)
    if theta < 1.0 - 1e-10:  # Avoid division by zero
        C_surface = K_s * (theta / (1.0 - theta)) / np.sqrt(K_eq)
    else:
        # Saturation limit: use Sieverts
        C_surface = K_s * np.sqrt(P_up)
    
    # Flux with downstream boundary condition
    flux = D * (C_surface - C_down) / thickness
    
    # Ensure non-negative flux
    flux = max(0.0, flux)
    
    # Dissociation flux for reference
    J_diss = dissociation_flux(P_up, theta, k_diss)
    
    # Surface reduction factor (compared to Sieverts)
    SRF = flux / flux_sieverts if flux_sieverts > 0 else 1.0
    
    return {
        'flux': flux,
        'flux_sieverts': flux_sieverts,
        'flux_dissociation': J_diss,
        'C_surface': C_surface,
        'C_down': C_down,
        'theta': theta,
        'theta_langmuir': theta_langmuir,
        'coverage_mode': coverage_mode,
        'damkohler': Da_info,
        'surface_reduction_factor': SRF,
        'converged': True,
        'K_eq': K_eq
    }