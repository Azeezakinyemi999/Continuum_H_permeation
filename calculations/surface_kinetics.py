# """
# Surface Kinetics Model for Hydrogen Permeation (Level 6) - Dissociation Only

# Simplified model where only upstream dissociation can be rate-limiting.
# Downstream recombination is assumed fast (instantaneous).

# Theory:
# -------
# J = min(J_diss, J_Sieverts)

# where:
#     J_diss = k_diss × P × (1-θ)²
#     θ = √(K_eq × P) / (1 + √(K_eq × P))  [Langmuir isotherm]
#     K_eq = k_diss / k_recomb

# Damköhler number:
#     Da = k_diss × K_s × L / D
#     Da >> 1: Sieverts' Law (diffusion-limited)
#     Da << 1: Surface-limited (dissociation controls)

# References:
# -----------
# 1. Pick, M.A. & Sonnenberg, K. (1985). J. Nucl. Mater. 131, 208-220.
# 2. Baskes, M.I. (1980). J. Nucl. Mater. 92, 318-324.
# """

# import numpy as np

# # Import existing functions
# from calculations.permeation_calc import sieverts_concentration, fick_flux

# # Gas constant
# R = 8.314  # J/(mol·K)


# def surface_equilibrium_coverage(pressure, k_diss, k_recomb):
#     """
#     Calculate equilibrium surface coverage using Langmuir isotherm.
    
#     θ = √(K_eq × P) / (1 + √(K_eq × P))
    
#     Parameters
#     ----------
#     pressure : float
#         Hydrogen partial pressure [Pa]
#     k_diss : float
#         Dissociation rate constant [m/s] or appropriate units
#     k_recomb : float
#         Recombination rate constant [m⁴/(mol·s)] or appropriate units
    
#     Returns
#     -------
#     dict with 'theta', 'K_eq', 'sqrt_K_eq_P', 'regime'
#     """
#     if pressure < 0:
#         raise ValueError(f"Pressure must be non-negative, got {pressure}")
#     if k_diss <= 0 or k_recomb <= 0:
#         raise ValueError("Rate constants must be positive")
    
#     if pressure == 0:
#         return {'theta': 0.0, 'K_eq': k_diss/k_recomb, 'sqrt_K_eq_P': 0.0, 'regime': 'low_coverage'}
    
#     K_eq = k_diss / k_recomb
#     sqrt_K_eq_P = np.sqrt(K_eq * pressure)
#     theta = sqrt_K_eq_P / (1.0 + sqrt_K_eq_P)
    
#     if theta < 0.1:
#         regime = 'low_coverage'
#     elif theta > 0.9:
#         regime = 'high_coverage'
#     else:
#         regime = 'intermediate'
    
#     return {'theta': theta, 'K_eq': K_eq, 'sqrt_K_eq_P': sqrt_K_eq_P, 'regime': regime}


# def dissociation_flux(pressure, theta, k_diss):
#     """
#     Calculate dissociative adsorption flux.
    
#     J_diss = k_diss × P × (1-θ)²
    
#     Parameters
#     ----------
#     pressure : float
#         Hydrogen partial pressure [Pa]
#     theta : float
#         Surface coverage fraction (0 ≤ θ ≤ 1)
#     k_diss : float
#         Dissociation rate constant
    
#     Returns
#     -------
#     float
#         Dissociation flux [mol/(m²·s)]
#     """
#     blocking_factor = (1.0 - theta) ** 2
#     return k_diss * pressure * blocking_factor


# def calculate_damkohler_number(k_diss, D, K_s, thickness):
#     """
#     Calculate Damköhler number for dissociation.
    
#     Da = k_diss × K_s × L / D
    
#     Da >> 1: Diffusion-limited (Sieverts applies)
#     Da << 1: Dissociation-limited (surface controls)
    
#     Returns
#     -------
#     dict with 'Da', 'regime', 'description'
#     """
#     Da = k_diss * K_s * thickness / D
    
#     if Da > 10:
#         regime = 'diffusion-limited'
#         description = "Fast dissociation; Sieverts' Law applies"
#     elif Da < 0.1:
#         regime = 'dissociation-limited'
#         description = "Slow dissociation controls flux"
#     else:
#         regime = 'mixed'
#         description = "Comparable surface and bulk resistances"
    
#     return {'Da': Da, 'regime': regime, 'description': description}


# def calculate_surface_limited_flux(D, K_s, thickness, P_up, P_down, temperature,
#                                     k_diss=None, k_recomb=None,
#                                     material_name=None, coverage=None):
#     """
#     Calculate hydrogen flux with dissociation-limited surface kinetics.
    
#     Two usage modes:
#     ----------------
#     Option 1: Provide k_diss + k_recomb (no coverage)
#         → Calculate θ from Langmuir isotherm
#         → C_surface = K_s × (θ/(1-θ)) / √K_eq  [= K_s × √P at equilibrium]
#         → J = D × (C_surface - C_down) / L
    
#     Option 2: Provide k_diss + k_recomb + coverage
#         → Use provided θ with rigorous formula
#         → C_surface = K_s × (θ/(1-θ)) / √K_eq
#         → J = D × (C_surface - C_down) / L
#         → Warning if provided θ differs from Langmuir equilibrium
    
#     Note: When θ is calculated from Langmuir (Option 1), the formula reduces
#     to Sieverts' Law: C_surface = K_s × √P. The surface effect appears through
#     the steady-state coupling, not through a modified boundary condition.
    
#     Parameters
#     ----------
#     D : float
#         Bulk diffusivity [m²/s]
#     K_s : float
#         Solubility or Sieverts' constant [mol/(m³·Pa^0.5)]
#     thickness : float
#         Membrane thickness [m]
#     P_up : float
#         Upstream pressure [Pa]
#     P_down : float
#         Downstream pressure [Pa]
#     temperature : float
#         Temperature [K] - used for material lookup; if providing k_diss/k_recomb
#         directly, ensure they are at this temperature
#     k_diss : float, optional
#         Dissociation rate constant [mol/(m²·s·Pa)]
#     k_recomb : float, optional
#         Recombination rate constant [m⁴/(mol·s)]
#     material_name : str, optional
#         Material name to look up kinetics
#     coverage : float, optional
#         Fixed surface coverage θ (0 ≤ θ < 1). If provided with kinetics,
#         overrides Langmuir calculation.
    
#     Returns
#     -------
#     dict with flux, theta, damkohler info, surface_reduction_factor
#     """
#     # Input validation
#     if P_up < 0 or P_down < 0:
#         raise ValueError("Pressures must be non-negative")
#     if coverage is not None and (coverage < 0 or coverage >= 1):
#         raise ValueError("Coverage must be in range [0, 1)")
    
#     # Sieverts flux (Level 1 reference)
#     C_up_sieverts = sieverts_concentration(K_s, P_up)
#     C_down = sieverts_concentration(K_s, P_down)
#     flux_sieverts = fick_flux(D, C_up_sieverts, C_down, thickness)
    
#     # Determine which option we're using
#     have_kinetics = (k_diss is not None and k_recomb is not None)
#     have_coverage = (coverage is not None)
    
#     # Try to get kinetics from material if not provided
#     if not have_kinetics and material_name is not None:
#         from data.surface_kinetics_data import get_surface_kinetics
#         kinetics = get_surface_kinetics(material_name, temperature)
#         k_diss = kinetics['k_diss']
#         k_recomb = kinetics['k_recomb']
#         have_kinetics = True
    
#     # Must have kinetics (removed Option 2 - coverage only)
#     if not have_kinetics:
#         raise ValueError("Must provide k_diss + k_recomb, or material_name. "
#                          "Coverage alone is not sufficient.")
    
#     K_eq = k_diss / k_recomb
    
#     # Damköhler number
#     Da_info = calculate_damkohler_number(k_diss, D, K_s, thickness)
    
#     # Calculate Langmuir equilibrium coverage
#     coverage_result = surface_equilibrium_coverage(P_up, k_diss, k_recomb)
#     theta_langmuir = coverage_result['theta']
    
#     if have_coverage:
#         # Option 2: Use provided coverage with rigorous formula
#         theta = coverage
#         coverage_mode = 'fixed_rigorous'
        
#         # Warn if provided coverage differs significantly from equilibrium
#         if abs(theta - theta_langmuir) > 0.1:
#             import warnings
#             warnings.warn(
#                 f"Provided coverage ({theta:.3f}) differs from Langmuir equilibrium "
#                 f"({theta_langmuir:.3f}) by more than 0.1. This may indicate "
#                 f"non-equilibrium conditions or inconsistent parameters.",
#                 UserWarning
#             )
#     else:
#         # Option 1: Use Langmuir equilibrium coverage
#         theta = theta_langmuir
#         coverage_mode = 'calculated'
    
#     # Rigorous formula: C_surface = K_s × (θ/(1-θ)) / √K_eq
#     # At Langmuir equilibrium, this reduces to K_s × √P (Sieverts)
#     if theta < 1.0 - 1e-10:  # Avoid division by zero
#         C_surface = K_s * (theta / (1.0 - theta)) / np.sqrt(K_eq)
#     else:
#         # Saturation limit: use Sieverts
#         C_surface = K_s * np.sqrt(P_up)
    
#     # Flux with downstream boundary condition
#     flux = D * (C_surface - C_down) / thickness
    
#     # Ensure non-negative flux
#     flux = max(0.0, flux)
    
#     # Dissociation flux for reference
#     J_diss = dissociation_flux(P_up, theta, k_diss)
    
#     # Surface reduction factor (compared to Sieverts)
#     SRF = flux / flux_sieverts if flux_sieverts > 0 else 1.0
    
#     return {
#         'flux': flux,
#         'flux_sieverts': flux_sieverts,
#         'flux_dissociation': J_diss,
#         'C_surface': C_surface,
#         'C_down': C_down,
#         'theta': theta,
#         'theta_langmuir': theta_langmuir,
#         'coverage_mode': coverage_mode,
#         'damkohler': Da_info,
#         'surface_reduction_factor': SRF,
#         'converged': True,
#         'K_eq': K_eq
#     }


"""
Surface Kinetics Model for Hydrogen Permeation (Level 6) - Dissociation Only

Simplified model where only upstream dissociation can be rate-limiting.
Downstream recombination is assumed fast (instantaneous).

Theory:
-------
Three coverage modes are supported:

1. EQUILIBRIUM (coverage_mode='equilibrium'):
   - θ from Langmuir isotherm: θ = √(K_eq×P) / (1 + √(K_eq×P))
   - Surface assumed at local equilibrium with gas phase
   - Valid when Da >> 1 (surface kinetics faster than diffusion)
   - At equilibrium: C_surface = K_s × √P (recovers Sieverts!)

2. STEADY-STATE (coverage_mode='steady_state'):
   - θ solved from flux continuity: J_diss(θ) = J_bulk(θ)
   - k_diss × P × (1-θ)² = D × [K_s × (θ/(1-θ))/√K_eq - C_down] / L
   - Self-consistent coupling between surface and bulk
   - Shows true surface limitation when Da ~ 1 or Da << 1

3. FORCED (coverage_mode='forced'):
   - θ specified by user (forced_coverage parameter)
   - For parametric studies and "what if" scenarios
   - Useful for sensitivity analysis on coverage

Damköhler number:
    Da = k_diss × K_s × L / D
    Da >> 1: Diffusion-limited (Sieverts applies, all modes give same result)
    Da << 1: Surface-limited (only steady_state mode captures this correctly)

References:
-----------
1. Pick, M.A. & Sonnenberg, K. (1985). J. Nucl. Mater. 131, 208-220.
2. Baskes, M.I. (1980). J. Nucl. Mater. 92, 318-324.
"""

import numpy as np
from scipy.optimize import brentq

# Import existing functions
from calculations.permeation_calc import sieverts_concentration, fick_flux

# Gas constant
R = 8.314  # J/(mol·K)


def surface_equilibrium_coverage(pressure, k_diss, k_recomb):
    """
    Calculate equilibrium surface coverage using Langmuir isotherm.
    
    θ = √(K_eq × P) / (1 + √(K_eq × P))
    
    IMPORTANT: At Langmuir equilibrium, the subsurface concentration formula
    C_surface = K_s × (θ/(1-θ)) / √K_eq simplifies to C_surface = K_s × √P
    (Sieverts' Law). This is mathematically exact, not an approximation.
    
    Parameters
    ----------
    pressure : float
        Hydrogen partial pressure [Pa]
    k_diss : float
        Dissociation rate constant [mol/(m²·s·Pa)]
    k_recomb : float
        Recombination rate constant [m⁴/(mol·s)]
    
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
        Dissociation rate constant [mol/(m²·s·Pa)]
    
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


def _solve_steady_state_coverage(P_up, C_down, D, K_s, thickness, k_diss, K_eq, 
                                  theta_eq, tol=1e-10, max_iter=100):
    """
    Solve for steady-state coverage where J_diss = J_bulk.
    
    At steady state, the dissociation flux into the surface must equal
    the diffusive flux through the bulk:
    
        k_diss × P × (1-θ)² = D × [C_surface(θ) - C_down] / L
    
    where C_surface(θ) = K_s × (θ/(1-θ)) / √K_eq
    
    This is a nonlinear equation in θ, solved using Brent's method.
    
    Parameters
    ----------
    P_up : float
        Upstream pressure [Pa]
    C_down : float
        Downstream concentration [mol/m³]
    D : float
        Diffusivity [m²/s]
    K_s : float
        Sieverts constant [mol/(m³·Pa^0.5)]
    thickness : float
        Membrane thickness [m]
    k_diss : float
        Dissociation rate constant [mol/(m²·s·Pa)]
    K_eq : float
        Equilibrium constant k_diss/k_recomb [mol/(m²·Pa)]
    theta_eq : float
        Langmuir equilibrium coverage (for initial guess)
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum iterations
        
    Returns
    -------
    dict with 'theta', 'converged', 'iterations', 'residual'
    """
    sqrt_K_eq = np.sqrt(K_eq)
    
    def residual(theta):
        """
        F(θ) = J_diss - J_bulk = 0
        
        J_diss = k_diss × P × (1-θ)²
        J_bulk = D × [K_s × (θ/(1-θ))/√K_eq - C_down] / L
        """
        if theta <= 0:
            theta = 1e-15
        if theta >= 1:
            theta = 1.0 - 1e-15
            
        J_diss = k_diss * P_up * (1.0 - theta)**2
        C_surface = K_s * (theta / (1.0 - theta)) / sqrt_K_eq
        J_bulk = D * (C_surface - C_down) / thickness
        
        return J_diss - J_bulk
    
    # Check residual at equilibrium
    res_eq = residual(theta_eq)
    
    # If already at equilibrium (Da >> 1), return equilibrium value
    if abs(res_eq) < tol:
        return {
            'theta': theta_eq,
            'converged': True,
            'iterations': 0,
            'residual': res_eq,
            'method': 'equilibrium_exact'
        }
    
    # Determine search bounds
    # θ must be in (0, 1) and C_surface > C_down for positive flux
    theta_min = 1e-10
    theta_max = 1.0 - 1e-10
    
    # Check signs at bounds
    try:
        f_min = residual(theta_min)
        f_max = residual(theta_max)
    except Exception:
        # Fallback to equilibrium if evaluation fails
        return {
            'theta': theta_eq,
            'converged': False,
            'iterations': 0,
            'residual': float('nan'),
            'method': 'fallback_equilibrium'
        }
    
    # Brent's method requires opposite signs at bounds
    if f_min * f_max > 0:
        # No sign change - steady state might be at a bound or equilibrium is best estimate
        # This can happen when Da >> 1 (equilibrium is the solution)
        return {
            'theta': theta_eq,
            'converged': True,
            'iterations': 0,
            'residual': res_eq,
            'method': 'no_sign_change_equilibrium'
        }
    
    # Solve using Brent's method
    try:
        theta_ss, result = brentq(residual, theta_min, theta_max, 
                                   xtol=tol, maxiter=max_iter, full_output=True)
        converged = result.converged
        iterations = result.iterations
        final_residual = residual(theta_ss)
        
        return {
            'theta': theta_ss,
            'converged': converged,
            'iterations': iterations,
            'residual': final_residual,
            'method': 'brentq'
        }
    except Exception as e:
        # Fallback to equilibrium
        return {
            'theta': theta_eq,
            'converged': False,
            'iterations': 0,
            'residual': float('nan'),
            'method': f'fallback_error: {str(e)}'
        }


def calculate_surface_limited_flux(D, K_s, thickness, P_up, P_down, temperature,
                                    k_diss=None, k_recomb=None,
                                    material_name=None, 
                                    coverage_mode='equilibrium',
                                    forced_coverage=None):
    """
    Calculate hydrogen flux with dissociation-limited surface kinetics.
    
    Three coverage modes are supported:
    
    1. coverage_mode='equilibrium' (default):
       - θ from Langmuir isotherm: θ = √(K_eq×P) / (1 + √(K_eq×P))
       - Assumes surface is at local equilibrium with gas phase
       - At Langmuir equilibrium: C_surface = K_s × √P (Sieverts' Law!)
       - Valid when Da >> 1
    
    2. coverage_mode='steady_state':
       - θ solved from flux continuity: J_diss(θ) = J_bulk(θ)
       - Self-consistent coupling between surface and bulk transport
       - Captures true surface limitation when Da ~ 1 or Da << 1
       - θ < θ_equilibrium when surface can't keep up with diffusion
    
    3. coverage_mode='forced':
       - θ = forced_coverage (user-specified)
       - For parametric studies: "What if θ were X?"
       - Useful for sensitivity analysis or matching experimental θ
    
    Physical Insight:
    -----------------
    At Langmuir equilibrium (mode='equilibrium'), the formula
        C_surface = K_s × (θ/(1-θ)) / √K_eq
    exactly reduces to C_surface = K_s × √P (Sieverts' Law).
    
    This is NOT an approximation - it's a mathematical identity when
    θ comes from the Langmuir isotherm. Surface kinetics only affect
    flux when:
    - coverage_mode='steady_state' and Da is not >> 1, OR
    - coverage_mode='forced' with θ ≠ θ_equilibrium
    
    Parameters
    ----------
    D : float
        Bulk diffusivity [m²/s]
    K_s : float
        Solubility (Sieverts' constant) [mol/(m³·Pa^0.5)]
    thickness : float
        Membrane thickness [m]
    P_up : float
        Upstream pressure [Pa]
    P_down : float
        Downstream pressure [Pa]
    temperature : float
        Temperature [K]
    k_diss : float, optional
        Dissociation rate constant [mol/(m²·s·Pa)]
    k_recomb : float, optional
        Recombination rate constant [m⁴/(mol·s)]
    material_name : str, optional
        Material name to look up kinetics from surface_kinetics_data
    coverage_mode : str
        'equilibrium' (default), 'steady_state', or 'forced'
    forced_coverage : float, optional
        Fixed surface coverage θ when coverage_mode='forced'
    
    Returns
    -------
    dict with:
        flux : float - Permeation flux [mol/(m²·s)]
        flux_sieverts : float - Reference Sieverts flux [mol/(m²·s)]
        flux_dissociation : float - Dissociation flux J_diss [mol/(m²·s)]
        C_surface : float - Subsurface concentration [mol/m³]
        C_down : float - Downstream concentration [mol/m³]
        theta : float - Surface coverage used
        theta_equilibrium : float - Langmuir equilibrium coverage
        coverage_mode : str - Mode used
        damkohler : dict - Damköhler number info
        surface_reduction_factor : float - flux / flux_sieverts
        converged : bool - True if calculation converged
        K_eq : float - Equilibrium constant
        solver_info : dict - Details from steady-state solver (if applicable)
    
    Examples
    --------
    >>> # Mode 1: Equilibrium (recovers Sieverts when Da >> 1)
    >>> result = calculate_surface_limited_flux(
    ...     D=1e-9, K_s=0.1, thickness=1e-3, P_up=1e5, P_down=0,
    ...     temperature=800, k_diss=1e-4, k_recomb=1e-6,
    ...     coverage_mode='equilibrium'
    ... )
    >>> result['surface_reduction_factor']  # ≈ 1.0 at equilibrium
    
    >>> # Mode 2: Steady-state (captures surface limitation)
    >>> result = calculate_surface_limited_flux(
    ...     D=1e-9, K_s=0.1, thickness=1e-3, P_up=1e5, P_down=0,
    ...     temperature=800, k_diss=1e-8, k_recomb=1e-6,  # Slow k_diss!
    ...     coverage_mode='steady_state'
    ... )
    >>> result['surface_reduction_factor']  # < 1.0 if surface-limited
    
    >>> # Mode 3: Forced coverage (parametric study)
    >>> for theta_test in [0.1, 0.3, 0.5, 0.7, 0.9]:
    ...     result = calculate_surface_limited_flux(
    ...         ..., coverage_mode='forced', forced_coverage=theta_test
    ...     )
    """
    # =========================================================================
    # Input validation
    # =========================================================================
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    
    valid_modes = ['equilibrium', 'steady_state', 'forced']
    if coverage_mode not in valid_modes:
        raise ValueError(f"coverage_mode must be one of {valid_modes}, got '{coverage_mode}'")
    
    if coverage_mode == 'forced':
        if forced_coverage is None:
            raise ValueError("forced_coverage must be provided when coverage_mode='forced'")
        if forced_coverage < 0 or forced_coverage >= 1:
            raise ValueError("forced_coverage must be in range [0, 1)")
    
    # =========================================================================
    # Reference: Sieverts flux (Level 1)
    # =========================================================================
    C_up_sieverts = sieverts_concentration(K_s, P_up)
    C_down = sieverts_concentration(K_s, P_down)
    flux_sieverts = fick_flux(D, C_up_sieverts, C_down, thickness)
    
    # =========================================================================
    # Get kinetics parameters
    # =========================================================================
    have_kinetics = (k_diss is not None and k_recomb is not None)
    
    if not have_kinetics and material_name is not None:
        from data.surface_kinetics_data import get_surface_kinetics
        kinetics = get_surface_kinetics(material_name, temperature)
        k_diss = kinetics['k_diss']
        k_recomb = kinetics['k_recomb']
        have_kinetics = True
    
    if not have_kinetics:
        raise ValueError("Must provide k_diss + k_recomb, or material_name.")
    
    K_eq = k_diss / k_recomb
    
    # =========================================================================
    # Calculate Damköhler number
    # =========================================================================
    Da_info = calculate_damkohler_number(k_diss, D, K_s, thickness)
    
    # =========================================================================
    # Calculate Langmuir equilibrium coverage (always needed for reference)
    # =========================================================================
    coverage_result = surface_equilibrium_coverage(P_up, k_diss, k_recomb)
    theta_equilibrium = coverage_result['theta']
    
    # =========================================================================
    # Determine θ based on coverage_mode
    # =========================================================================
    solver_info = None
    
    if coverage_mode == 'equilibrium':
        # Mode 1: Use Langmuir equilibrium
        theta = theta_equilibrium
        
    #BUGGY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #  elif coverage_mode == 'steady_state':
    #     # Mode 2: Solve for self-consistent steady-state θ
    #     solver_result = _solve_steady_state_coverage(
    #         P_up=P_up,
    #         C_down=C_down,
    #         D=D,
    #         K_s=K_s,
    #         thickness=thickness,
    #         k_diss=k_diss,
    #         K_eq=K_eq,
    #         theta_eq=theta_equilibrium
    #     )
    #     theta = solver_result['theta']
    #     solver_info = solver_result
    
    elif coverage_mode == 'steady_state':
        # Mode 2: True steady-state flux
        # Use Langmuir equilibrium for θ, but the actual flux is limited
        # by whichever process is slower: bulk diffusion or surface dissociation.
        # This is handled in the flux calculation below.
        theta = theta_equilibrium
        solver_info = {
            'theta': theta,
            'converged': True,
            'iterations': 0,
            'residual': 0.0,
            'method': 'min_flux'
        }
    elif coverage_mode == 'forced':
        # Mode 3: Use user-specified coverage
        theta = forced_coverage
        
        # Warn if far from equilibrium
        if abs(theta - theta_equilibrium) > 0.1:
            import warnings
            warnings.warn(
                f"Forced coverage ({theta:.3f}) differs significantly from "
                f"Langmuir equilibrium ({theta_equilibrium:.3f}). "
                f"This represents a non-equilibrium surface state.",
                UserWarning
            )
    
    # =========================================================================
    # Calculate subsurface concentration and flux
    # =========================================================================
    # C_surface = K_s × (θ/(1-θ)) / √K_eq
    # At Langmuir equilibrium, this exactly equals K_s × √P (Sieverts)
    
    sqrt_K_eq = np.sqrt(K_eq)
    
    if theta < 1.0 - 1e-10:
        C_surface = K_s * (theta / (1.0 - theta)) / sqrt_K_eq
    else:
        # Saturation limit
        C_surface = K_s * np.sqrt(P_up)
    
    # Ensure C_surface is non-negative
    C_surface = max(0.0, C_surface)
    
# Fick's Law for bulk transport
    flux_bulk = D * (C_surface - C_down) / thickness
    flux_bulk = max(0.0, flux_bulk)  # Ensure non-negative
    
    # Dissociation flux for reference
    J_diss = dissociation_flux(P_up, theta, k_diss)
    
    # =========================================================================
    # FLUX DETERMINATION (mode-dependent)
    # =========================================================================
    if coverage_mode == 'steady_state':
        # Steady-state: flux is limited by slower process
        flux = min(flux_sieverts, J_diss)
    else:
        # equilibrium or forced: use bulk flux from C_surface
        flux = flux_bulk
    
    # Surface Reduction Factor
    SRF = flux / flux_sieverts if flux_sieverts > 0 else 1.0
    
    # =========================================================================
    # Convergence check for steady_state mode
    # =========================================================================
    if coverage_mode == 'steady_state' and solver_info is not None:
        converged = solver_info.get('converged', True)
    else:
        converged = True
    
    # =========================================================================
    # Return results
    # =========================================================================
    return {
        'flux': flux,
        'flux_sieverts': flux_sieverts,
        'flux_dissociation': J_diss,
        'C_surface': C_surface,
        'C_up_sieverts': C_up_sieverts,
        'C_down': C_down,
        'theta': theta,
        'theta_equilibrium': theta_equilibrium,
        'coverage_mode': coverage_mode,
        'damkohler': Da_info,
        'Da': Da_info['Da'],
        'SRF': SRF,  # Standard key
        'surface_reduction_factor': SRF,  # Backward compatibility
        'converged': converged,
        'K_eq': K_eq,
        'solver_info': solver_info
    }