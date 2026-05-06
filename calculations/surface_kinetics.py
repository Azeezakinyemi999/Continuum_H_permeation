"""
surface_kinetics.py
===================
Importable module for hydrogen permeation surface kinetics models.

Hierarchy
---------
L1+L6  : Perfect metal + surface kinetics (no oxide)
L2a+L6 : Perfect oxide + surface kinetics (no metal block)
L2+L6  : Perfect oxide + perfect metal + surface kinetics
L2+L4+L6 : Perfect oxide + defective metal + surface kinetics
L3+L6  : Defective oxide + perfect metal + surface kinetics
L3+L4+L6 : Defective oxide + defective metal + surface kinetics (full model)

Usage
-----
    from surface_kinetics import solve_steady_state_flux_L1L6
    from surface_kinetics import calculate_full_model_flux_L346_v2
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
from scipy.optimize import brentq

from data.surface_kinetics_data import SURFACE_KINETICS, get_surface_kinetics
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS
from calculations.defective_metal import combined_microstructure_model

R = 8.314  # J/(mol·K)


# =============================================================================
# SECTION 1: Shared Helpers
# Used by all model levels
# =============================================================================

def get_all_properties(oxide_name, metal_name, temperature_K):
    """
    Look up and compute all transport properties from data modules.

    Returns surface kinetics, oxide transport, and metal transport
    at the given temperature via Arrhenius correction.
    """
    surfraction_kin = get_surface_kinetics(oxide_name, temperature_K)

    oxide_data = OXIDE_PROPERTIES[oxide_name]
    T_ref_ox = oxide_data['T_ref']
    D_ox = oxide_data['D_ox_ref'] * np.exp((-oxide_data['E_D_ox'] / R) * (1/temperature_K - 1/T_ref_ox))
    K_ox = oxide_data['K_ox_ref'] * np.exp((-oxide_data['H_sol_ox'] / R) * (1/temperature_K - 1/T_ref_ox))
    L_ox = oxide_data['thickness']

    metal_data = MATERIALS[metal_name]
    T_ref_m = metal_data['T_ref']
    D_m   = metal_data['D_ref']  * np.exp((-metal_data['E_D'] / R) * (1/temperature_K - 1/T_ref_m))
    K_s_m = metal_data['K_s_ref'] * np.exp((-metal_data['H_s']  / R) * (1/temperature_K - 1/T_ref_m))

    return {
        'k_diss':   surfraction_kin['k_diss'],
        'k_recomb': surfraction_kin['k_recomb'],
        'K_eq':     surfraction_kin['K_eq'],
        'D_ox': D_ox, 'K_ox': K_ox, 'L_ox': L_ox,
        'D_m':  D_m,  'K_s_m': K_s_m,
    }


def compute_permeances(props, L_m):
    """Return oxide permeance α and metal permeance β from a props dict."""
    alpha = props['D_ox'] * props['K_ox'] / props['L_ox']
    beta  = props['D_m']  * props['K_s_m'] / L_m
    return alpha, beta


def g_theta(theta, K_eq):
    """
    Surface isotherm: converts coverage θ to effective √P at the surface.

        g(θ) = θ / ((1 - θ) × √K_eq)
    """
    if theta >= 1.0:
        return np.inf
    return theta / ((1.0 - theta) * np.sqrt(K_eq))


def sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down):
    """
    Analytical oxide-metal interface pressure from flux balance (L2+L6).

        √P_int = (α × g(θ) + β × √P_down) / (α + β)
    """
    g = g_theta(theta, K_eq)
    return (alpha * g + beta * np.sqrt(P_down)) / (alpha + beta)


def smooth_surface_resistance(k_diss_eff, P_up, theta, eps_theta=1e-12, eps_rate=1e-60):
    """
    Smooth surface resistance with no hard cutoff.

        R_surface ~ 1 / (k_diss × P_up × (1 - θ)²)
    """
    theta_clamped = min(max(theta, 0.0), 1.0 - eps_theta)
    vacant = max(1.0 - theta_clamped, eps_theta)
    rate = k_diss_eff * max(P_up, eps_theta) * (vacant ** 2)
    return 1.0 / max(rate, eps_rate)


def surface_flux(theta, P_up, k_diss, K_eq):
    """Langmuir-Hinshelwood net dissociation flux at the gas-surface interface."""
    k_recomb = k_diss / K_eq
    return k_diss * P_up * (1 - theta)**2 - k_recomb * theta**2


def oxide_flux(theta, alpha, beta, K_eq, P_down):
    """Oxide diffusion flux using the oxide-metal interface pressure (L2+L6 form)."""
    g        = g_theta(theta, K_eq)
    sqrt_P_int = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
    return alpha * (g - sqrt_P_int)


def metal_flux(theta, alpha, beta, K_eq, P_down):
    """Metal diffusion flux using the oxide-metal interface pressure (L2+L6 form)."""
    sqrt_P_int  = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
    sqrt_P_down = np.sqrt(P_down)
    return beta * (sqrt_P_int - sqrt_P_down)


def surface_flux_residual(theta, P_up, alpha, beta, P_down, k_diss, K_eq):
    """Root-finding residual for L2+L6: J_surface - J_oxide = 0."""
    return surface_flux(theta, P_up, k_diss, K_eq) - oxide_flux(theta, alpha, beta, K_eq, P_down)


# =============================================================================
# SECTION 2: L1+L6 — Perfect Metal + Surface Kinetics (No Oxide)
# =============================================================================

def metal_flux_L1(theta, beta, K_eq, P_down):
    """Metal flux for L1+L6: √P_int = g(θ) directly, no oxide layer."""
    return beta * (g_theta(theta, K_eq) - np.sqrt(P_down))


def surface_flux_residual_L1(theta, P_up, beta, P_down, k_diss, K_eq):
    """Residual for L1+L6: f(θ) = J_surface - J_metal = 0."""
    return surface_flux(theta, P_up, k_diss, K_eq) - metal_flux_L1(theta, beta, K_eq, P_down)


def solve_steady_state_flux_L1L6(P_up, P_down, L_m, k_diss, K_eq, D_m, K_s_m):
    """
    Solve L1+L6 steady-state flux: perfect metal + surface kinetics, no oxide.

    Parameters
    ----------
    P_up, P_down : float — upstream/downstream H2 pressure [Pa]
    L_m          : float — metal thickness [m]
    k_diss       : float — dissociation rate constant [mol/m²/s/Pa]
    K_eq         : float — surface equilibrium constant [Pa⁻¹]
    D_m          : float — metal diffusivity [m²/s]
    K_s_m        : float — metal Sieverts constant [mol/m³/Pa^0.5]

    Returns
    -------
    dict : theta, P_int, J_ss, beta, rate_limiting, resistances
    """
    beta = D_m * K_s_m / L_m

    theta_ss = brentq(
        surface_flux_residual_L1,
        1e-10, 1.0 - 1e-10,
        args=(P_up, beta, P_down, k_diss, K_eq)
    )

    sqrt_P_int = g_theta(theta_ss, K_eq)
    J_ss       = metal_flux_L1(theta_ss, beta, K_eq, P_down)

    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
    R_metal   = 1.0 / beta
    R_total   = R_surface + R_metal

    fraction_surface = R_surface / R_total
    fraction_metal   = R_metal   / R_total

    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'

    return {
        'theta':         theta_ss,
        'P_int':         sqrt_P_int**2,
        'J_ss':          J_ss,
        'beta':          beta,
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface':        R_surface,
            'R_metal':          R_metal,
            'R_total':          R_total,
            'fraction_surface': fraction_surface,
            'fraction_metal':   fraction_metal,
        },
    }


# =============================================================================
# SECTION 3: L2a+L6 — Perfect Oxide + Surface Kinetics (No Metal Block)
# =============================================================================

def oxide_flux_L2a(theta, alpha, K_eq, P_down):
    """Oxide flux for L2a+L6: √P_surf = g(θ) directly, no metal block."""
    return alpha * (g_theta(theta, K_eq) - np.sqrt(P_down))


def surface_flux_residual_L2a(theta, P_up, alpha, P_down, k_diss, K_eq):
    """Residual for L2a+L6: f(θ) = J_surface - J_oxide = 0."""
    return surface_flux(theta, P_up, k_diss, K_eq) - oxide_flux_L2a(theta, alpha, K_eq, P_down)


def solve_steady_state_flux_L2aL6(P_up, P_down, k_diss, K_eq, D_ox, K_ox, L_ox):
    """
    Solve L2a+L6 steady-state flux: perfect oxide + surface kinetics, no metal block.

    Parameters
    ----------
    P_up, P_down : float — upstream/downstream H2 pressure [Pa]
    k_diss       : float — dissociation rate constant [mol/m²/s/Pa]
    K_eq         : float — surface equilibrium constant [Pa⁻¹]
    D_ox         : float — oxide diffusivity [m²/s]
    K_ox         : float — oxide Sieverts constant [mol/m³/Pa^0.5]
    L_ox         : float — oxide thickness [m]

    Returns
    -------
    dict : theta, P_surf, J_ss, alpha, rate_limiting, resistances
    """
    alpha = D_ox * K_ox / L_ox

    theta_ss = brentq(
        surface_flux_residual_L2a,
        1e-10, 1.0 - 1e-10,
        args=(P_up, alpha, P_down, k_diss, K_eq)
    )

    sqrt_P_surf = g_theta(theta_ss, K_eq)
    J_ss        = oxide_flux_L2a(theta_ss, alpha, K_eq, P_down)

    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
    R_oxide   = 1.0 / alpha
    R_total   = R_surface + R_oxide

    fraction_surface = R_surface / R_total
    fraction_oxide   = R_oxide   / R_total

    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    else:
        rate_limiting = 'mixed'

    return {
        'theta':         theta_ss,
        'P_surf':        sqrt_P_surf**2,
        'J_ss':          J_ss,
        'alpha':         alpha,
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface':        R_surface,
            'R_oxide':          R_oxide,
            'R_total':          R_total,
            'fraction_surface': fraction_surface,
            'fraction_oxide':   fraction_oxide,
        },
    }


# =============================================================================
# SECTION 4: L2+L6 — Perfect Oxide + Perfect Metal + Surface Kinetics
# =============================================================================

def solve_steady_state_flux(P_up, P_down, L_m, oxide_name, metal_name, temperature_K):
    """
    Solve L2+L6 using data module lookup (Arrhenius-corrected properties).

    Parameters
    ----------
    P_up, P_down  : float — upstream/downstream H2 pressure [Pa]
    L_m           : float — metal thickness [m]
    oxide_name    : str   — key in OXIDE_PROPERTIES
    metal_name    : str   — key in MATERIALS
    temperature_K : float — temperature [K]

    Returns
    -------
    dict : theta, P_int, J_ss, J_surface, J_oxide, J_metal,
           alpha, beta, rate_limiting, resistances
    """
    props = get_all_properties(oxide_name, metal_name, temperature_K)
    alpha   = props['D_ox'] * props['K_ox'] / props['L_ox']
    beta    = props['D_m']  * props['K_s_m'] / L_m
    k_diss  = props['k_diss']
    K_eq    = props['K_eq']

    theta_ss = brentq(
        surface_flux_residual,
        1e-10, 1.0 - 1e-10,
        args=(P_up, alpha, beta, P_down, k_diss, K_eq)
    )

    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
    P_int      = sqrt_P_int**2
    J_ss       = metal_flux(theta_ss, alpha, beta, K_eq, P_down)
    J_surf     = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_ox       = oxide_flux(theta_ss, alpha, beta, K_eq, P_down)

    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
    R_oxide   = 1.0 / alpha
    R_metal   = 1.0 / beta
    R_total   = R_surface + R_oxide + R_metal

    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide   = R_oxide   / R_total if R_total > 0 else 0
    fraction_metal   = R_metal   / R_total if R_total > 0 else 0

    if fraction_surface > 0.5:
        rate_limiting = 'surface (dissociation)'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide (diffusion)'
    else:
        rate_limiting = 'metal (diffusion)'

    return {
        'theta':     theta_ss,
        'P_int':     P_int,
        'J_ss':      J_ss,
        'J_surface': J_surf,
        'J_oxide':   J_ox,
        'J_metal':   J_ss,
        'alpha':     alpha,
        'beta':      beta,
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface':        R_surface,
            'R_oxide':          R_oxide,
            'R_metal':          R_metal,
            'R_total':          R_total,
            'fraction_surface': fraction_surface,
            'fraction_oxide':   fraction_oxide,
            'fraction_metal':   fraction_metal,
        },
    }


def solve_steady_state_flux_direct(P_up, P_down, L_m, k_diss, K_eq,
                                    D_ox, K_ox, L_ox, D_m, K_s_m):
    """
    Solve L2+L6 with direct parameter input (no data module lookup).

    Parameters
    ----------
    P_up, P_down : float — upstream/downstream H2 pressure [Pa]
    L_m          : float — metal thickness [m]
    k_diss       : float — dissociation rate constant [mol/m²/s/Pa]
    K_eq         : float — surface equilibrium constant [Pa⁻¹]
    D_ox         : float — oxide diffusivity [m²/s]
    K_ox         : float — oxide Sieverts constant [mol/m³/Pa^0.5]
    L_ox         : float — oxide thickness [m]
    D_m          : float — metal diffusivity [m²/s]
    K_s_m        : float — metal Sieverts constant [mol/m³/Pa^0.5]

    Returns
    -------
    dict : theta, P_int, J_ss, J_surface, J_oxide, J_metal,
           alpha, beta, rate_limiting, resistances
    """
    alpha = D_ox * K_ox / L_ox
    beta  = D_m  * K_s_m / L_m

    theta_ss = brentq(
        surface_flux_residual,
        1e-10, 1.0 - 1e-10,
        args=(P_up, alpha, beta, P_down, k_diss, K_eq)
    )

    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
    P_int      = sqrt_P_int**2
    J_ss       = metal_flux(theta_ss, alpha, beta, K_eq, P_down)
    J_surf     = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_ox       = oxide_flux(theta_ss, alpha, beta, K_eq, P_down)

    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
    R_oxide   = 1.0 / alpha
    R_metal   = 1.0 / beta
    R_total   = R_surface + R_oxide + R_metal

    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide   = R_oxide   / R_total if R_total > 0 else 0
    fraction_metal   = R_metal   / R_total if R_total > 0 else 0

    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'

    return {
        'theta':     theta_ss,
        'P_int':     P_int,
        'J_ss':      J_ss,
        'J_surface': J_surf,
        'J_oxide':   J_ox,
        'J_metal':   J_ss,
        'alpha':     alpha,
        'beta':      beta,
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface':        R_surface,
            'R_oxide':          R_oxide,
            'R_metal':          R_metal,
            'R_total':          R_total,
            'fraction_surface': fraction_surface,
            'fraction_oxide':   fraction_oxide,
            'fraction_metal':   fraction_metal,
        },
    }


# =============================================================================
# SECTION 5: L2+L4+L6 — Perfect Oxide + Defective Metal + Surface Kinetics
# =============================================================================

def calculate_defective_metal_flux_L6(
    P_up, P_down, thickness, temperature,
    k_diss, K_eq,
    D_ox, K_ox, L_ox,
    D_lattice, K_s_m,
    microstructure_params,
    lattice_density=1.06e29,
    method='average', n_points=50, mode='both',
    max_iterations=15, tolerance=1e-6,
):
    """
    Solve L2+L4+L6: perfect oxide + defective metal + surface kinetics.

    Iterates to find self-consistent D_eff and θ_ss, since D_eff depends
    on the concentration profile which depends on θ_ss.

    Parameters
    ----------
    P_up, P_down      : float — upstream/downstream H2 pressure [Pa]
    thickness         : float — metal thickness [m]
    temperature       : float — temperature [K]
    k_diss, K_eq      : float — surface kinetics parameters
    D_ox, K_ox, L_ox  : float — oxide transport properties
    D_lattice, K_s_m  : float — metal lattice properties
    microstructure_params : dict — grain size, traps, etc.
    lattice_density   : float — lattice site density [sites/m³]
    method            : str  — D_eff averaging: 'average', 'harmonic', 'inlet', 'outlet'
    n_points          : int  — concentration profile discretisation points
    mode              : str  — microstructure mode: 'both', 'gb_only', 'trapping_only', 'none'
    max_iterations    : int  — convergence iteration limit
    tolerance         : float — relative tolerance for D_eff convergence

    Returns
    -------
    dict : flux, theta_surface, P_int, D_eff, profiles, rate_limiting, resistances, convergence
    """
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    if thickness <= 0:
        raise ValueError(f"Thickness must be positive: {thickness} m")
    if D_lattice <= 0:
        raise ValueError("Diffusion coefficient must be positive")
    if k_diss <= 0 or K_eq <= 0:
        raise ValueError("Surface kinetics parameters must be positive")

    alpha      = D_ox * K_ox / L_ox
    D_eff      = D_lattice
    sqrt_P_down = np.sqrt(max(P_down, 0))
    convergence_history = []

    for iteration in range(max_iterations):
        beta = D_eff * K_s_m / thickness

        try:
            theta_ss = brentq(
                surface_flux_residual,
                1e-10, 1.0 - 1e-10,
                args=(P_up, alpha, beta, P_down, k_diss, K_eq)
            )
        except ValueError as e:
            return {
                'flux': np.nan,
                'error': f'Failed to solve for theta at iteration {iteration}: {str(e)}',
                'convergence_history': convergence_history,
            }

        sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)

        x_array    = np.linspace(0, thickness, n_points)
        sqrt_P_arr = sqrt_P_int - (sqrt_P_int - sqrt_P_down) * x_array / thickness
        C_array    = K_s_m * sqrt_P_arr

        D_array          = np.zeros(n_points)
        theta_trap_array = np.zeros(n_points)
        gb_factor_array  = np.zeros(n_points)

        for i, C_local in enumerate(C_array):
            C_local = max(C_local, 1e-20)
            result_i = combined_microstructure_model(
                D_lattice=D_lattice,
                temperature=temperature,
                microstructure_params=microstructure_params,
                lattice_concentration=C_local,
                lattice_density=lattice_density,
                mode=mode,
            )
            D_array[i] = result_i['D_eff']
            if 'trapping' in result_i and result_i['trapping'] is not None:
                theta_trap_array[i] = result_i['trapping'].get('theta_total', 0.0)
            elif 'theta_total' in result_i:
                theta_trap_array[i] = result_i['theta_total']
            if 'gb_enhancement' in result_i and result_i['gb_enhancement'] is not None:
                gb_factor_array[i] = result_i['gb_enhancement'].get('factor', 1.0)
            elif 'gb_enhancement_factor' in result_i:
                gb_factor_array[i] = result_i['gb_enhancement_factor']
            else:
                gb_factor_array[i] = 1.0

        if method == 'harmonic':
            D_eff_new = len(D_array) / np.sum(1.0 / D_array)
        elif method == 'inlet':
            D_eff_new = D_array[0]
        elif method == 'outlet':
            D_eff_new = D_array[-1]
        else:
            D_eff_new = np.mean(D_array)

        rel_change = abs(D_eff_new - D_eff) / D_eff if D_eff > 0 else np.inf
        convergence_history.append({
            'iteration': iteration,
            'D_eff': D_eff_new,
            'theta': theta_ss,
            'relative_change': rel_change,
        })

        if rel_change < tolerance:
            D_eff = D_eff_new
            break
        D_eff = D_eff_new

    beta_eff   = D_eff * K_s_m / thickness
    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta_eff, K_eq, P_down)
    P_int      = sqrt_P_int**2
    C_up       = K_s_m * sqrt_P_int
    C_down     = K_s_m * sqrt_P_down if P_down > 0 else 0.0
    J_metal    = metal_flux(theta_ss, alpha, beta_eff, K_eq, P_down)
    J_surface  = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_oxide    = oxide_flux(theta_ss, alpha, beta_eff, K_eq, P_down)

    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
    R_oxide   = 1.0 / alpha
    R_metal   = 1.0 / beta_eff
    R_total   = R_surface + R_oxide + R_metal

    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide   = R_oxide   / R_total if R_total > 0 else 0
    fraction_metal   = R_metal   / R_total if R_total > 0 else 0

    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'

    modification_factor = D_eff / D_lattice
    avg_gb_factor       = np.mean(gb_factor_array)
    avg_theta_trap      = np.mean(theta_trap_array)
    trap_reduction      = 1.0 / (1.0 + avg_theta_trap) if avg_theta_trap > 0 else 1.0

    if avg_gb_factor > 1.5 and trap_reduction > 0.5:
        dominant_effect = 'gb_enhancement'
    elif avg_gb_factor < 1.5 and trap_reduction < 0.5:
        dominant_effect = 'trapping'
    else:
        dominant_effect = 'balanced'

    return {
        'flux':      J_metal,
        'J_surface': J_surface,
        'J_oxide':   J_oxide,
        'J_metal':   J_metal,
        'theta_surface': theta_ss,
        'P_int':     P_int,
        'C_up':      C_up,
        'C_down':    C_down,
        'alpha':     alpha,
        'beta_lattice': D_lattice * K_s_m / thickness,
        'beta_eff':  beta_eff,
        'D_eff':     D_eff,
        'D_lattice': D_lattice,
        'modification_factor': modification_factor,
        'microstructure_details': {
            'average_theta_trap':    avg_theta_trap,
            'average_gb_enhancement': avg_gb_factor,
            'trap_reduction_factor': trap_reduction,
            'dominant_effect':       dominant_effect,
            'method_used':           method,
        },
        'flux_balance': {
            'J_surface': J_surface,
            'J_oxide':   J_oxide,
            'J_metal':   J_metal,
            'balanced':  np.allclose(J_surface, J_oxide, rtol=1e-4) and
                         np.allclose(J_oxide, J_metal, rtol=1e-4),
        },
        'profiles': {
            'x': x_array, 'D': D_array, 'C': C_array,
            'theta_trap': theta_trap_array, 'gb_factor': gb_factor_array,
        },
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface':        R_surface,
            'R_oxide':          R_oxide,
            'R_metal':          R_metal,
            'R_total':          R_total,
            'fraction_surface': fraction_surface,
            'fraction_oxide':   fraction_oxide,
            'fraction_metal':   fraction_metal,
        },
        'convergence': {
            'iterations': len(convergence_history),
            'converged':  len(convergence_history) < max_iterations,
            'history':    convergence_history,
        },
    }


# =============================================================================
# SECTION 6: L3+L6 — Defective Oxide + Perfect Metal + Surface Kinetics
# =============================================================================

def calculate_path_flux_L6(
    P_up, P_down, L_m,
    k_diss, K_eq,
    alpha,
    D_m, K_s_m,
    path_type='intact',
    k_diss_metal=None,
    K_eq_metal=None,
):
    """
    Calculate steady-state flux through a single oxide path (L3+L6).

    Supports intact, pinhole, crack, and grain boundary paths.
    For pinhole: uses metal surface kinetics if k_diss_metal/K_eq_metal
    are provided, otherwise falls back to Sieverts' law limit.

    Parameters
    ----------
    P_up, P_down     : float      — upstream/downstream H2 pressure [Pa]
    L_m              : float      — metal thickness [m]
    k_diss, K_eq     : float      — oxide surface kinetics parameters
    alpha            : float      — oxide permeance (np.inf for pinhole)
    D_m, K_s_m       : float      — metal transport properties
    path_type        : str        — 'intact', 'pinhole', 'crack', 'grain_boundary'
    k_diss_metal     : float|None — metal surface dissociation rate [mol/m²/s/Pa]
    K_eq_metal       : float|None — metal surface equilibrium constant [Pa⁻¹]

    Returns
    -------
    dict : flux, theta, P_int, path_type, alpha, beta,
           kinetics_used, flux_balance, resistances, rate_limiting
    """
    beta        = D_m * K_s_m / L_m
    sqrt_P_down = np.sqrt(max(P_down, 0))
    sqrt_P_up   = np.sqrt(max(P_up,   0))
    is_pinhole  = (path_type == 'pinhole' or alpha == np.inf or alpha > 1e10)

    # -------------------------------------------------------------------------
    # Pinhole path
    # -------------------------------------------------------------------------
    if is_pinhole:
        if k_diss_metal is not None and K_eq_metal is not None:
            k_diss_eff      = k_diss_metal
            K_eq_eff        = K_eq_metal
            use_sieverts_limit = False
        else:
            use_sieverts_limit = True

        if use_sieverts_limit:
            sqrt_P_int = sqrt_P_up
            P_int      = P_up
            J_ss       = beta * (sqrt_P_int - sqrt_P_down)
            theta_ss   = 0.0
            J_surf     = J_ss
            J_ox       = np.nan
            R_surface  = 0.0
            R_oxide    = 0.0
            R_metal    = 1.0 / beta
            R_total    = R_metal
            fraction_surface = 0.0
            fraction_oxide   = 0.0
            fraction_metal   = 1.0
            rate_limiting    = 'metal'
        else:
            def residual(theta):
                if theta <= 0 or theta >= 1:
                    return np.inf
                J_surf_loc = surface_flux(theta, P_up, k_diss_eff, K_eq_eff)
                sqrt_P_int_loc = g_theta(theta, K_eq_eff)
                J_metal_loc    = beta * (sqrt_P_int_loc - sqrt_P_down)
                return J_surf_loc - J_metal_loc

            try:
                theta_ss = brentq(residual, 1e-10, 1.0 - 1e-10)
            except ValueError as e:
                return {
                    'flux': np.nan, 'theta': np.nan, 'P_int': np.nan,
                    'path_type': path_type, 'alpha': np.inf, 'beta': beta,
                    'error': str(e),
                }

            sqrt_P_int = g_theta(theta_ss, K_eq_eff)
            P_int      = sqrt_P_int**2
            J_ss       = beta * (sqrt_P_int - sqrt_P_down)
            J_surf     = surface_flux(theta_ss, P_up, k_diss_eff, K_eq_eff)
            J_ox       = np.nan
            R_surface  = smooth_surface_resistance(k_diss_eff, P_up, theta_ss)
            R_oxide    = 0.0
            R_metal    = 1.0 / beta
            R_total    = R_surface + R_metal
            fraction_surface = R_surface / R_total if R_total > 0 else 0
            fraction_oxide   = 0.0
            fraction_metal   = R_metal / R_total if R_total > 0 else 0
            if fraction_surface > 0.5:
                rate_limiting = 'surface'
            elif fraction_metal > 0.5:
                rate_limiting = 'metal'
            else:
                rate_limiting = 'mixed'

        return {
            'flux': J_ss, 'theta': theta_ss, 'P_int': P_int,
            'path_type': path_type, 'alpha': np.inf, 'beta': beta,
            'kinetics_used': 'sieverts' if use_sieverts_limit else 'metal',
            'flux_balance': {'J_surface': J_surf, 'J_oxide': J_ox, 'J_metal': J_ss, 'balanced': True},
            'resistances': {
                'R_surface': R_surface, 'R_oxide': R_oxide, 'R_metal': R_metal,
                'R_total': R_total, 'fraction_surface': fraction_surface,
                'fraction_oxide': fraction_oxide, 'fraction_metal': fraction_metal,
            },
            'rate_limiting': rate_limiting,
        }

    # -------------------------------------------------------------------------
    # Non-pinhole paths (intact, crack, GB)
    # -------------------------------------------------------------------------
    def residual(theta):
        if theta <= 0 or theta >= 1:
            return np.inf
        J_surf_loc     = surface_flux(theta, P_up, k_diss, K_eq)
        sqrt_P_int_loc = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
        J_metal_loc    = beta * (sqrt_P_int_loc - sqrt_P_down)
        return J_surf_loc - J_metal_loc

    try:
        theta_ss = brentq(residual, 1e-10, 1.0 - 1e-10)
    except ValueError as e:
        return {
            'flux': np.nan, 'theta': np.nan, 'P_int': np.nan,
            'path_type': path_type, 'alpha': alpha, 'beta': beta,
            'error': str(e),
        }

    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
    P_int      = sqrt_P_int**2
    J_ss       = beta * (sqrt_P_int - sqrt_P_down)
    J_surf     = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_ox       = oxide_flux(theta_ss, alpha, beta, K_eq, P_down)

    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
    R_oxide   = 1.0 / alpha
    R_metal   = 1.0 / beta
    R_total   = R_surface + R_oxide + R_metal

    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide   = R_oxide   / R_total if R_total > 0 else 0
    fraction_metal   = R_metal   / R_total if R_total > 0 else 0

    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'

    return {
        'flux': J_ss, 'theta': theta_ss, 'P_int': P_int,
        'path_type': path_type, 'alpha': alpha, 'beta': beta,
        'kinetics_used': 'oxide',
        'flux_balance': {
            'J_surface': J_surf, 'J_oxide': J_ox, 'J_metal': J_ss,
            'balanced': np.allclose(J_surf, J_ox, rtol=1e-6) and
                        np.allclose(J_ox, J_ss, rtol=1e-6),
        },
        'resistances': {
            'R_surface': R_surface, 'R_oxide': R_oxide, 'R_metal': R_metal,
            'R_total': R_total, 'fraction_surface': fraction_surface,
            'fraction_oxide': fraction_oxide, 'fraction_metal': fraction_metal,
        },
        'rate_limiting': rate_limiting,
    }


def calculate_parallel_path_flux_L6(
    P_up, P_down, L_m,
    k_diss, K_eq,
    D_ox, K_ox, L_ox,
    D_m, K_s_m,
    defect_area_fraction,
    defect_type='crack',
    thickness_factor=0.1,
    diffusivity_factor=100.0,
    k_diss_metal=None,
    K_eq_metal=None,
):
    """
    Area-weighted flux for one intact + one defect path (L3+L6).

        J_total = (1 - fraction_defect) × J_intact + fraction_defect × J_defect

    Parameters
    ----------
    defect_area_fraction : float — area fraction with defects (0 to 1)
    defect_type          : str  — 'pinhole', 'crack', or 'grain_boundary'
    thickness_factor     : float — γ for crack (L_crack = γ × L_ox)
    diffusivity_factor   : float — δ for GB (D_gb = δ × D_ox)
    k_diss_metal, K_eq_metal : float|None — metal kinetics for pinhole

    Returns
    -------
    dict : J_total, enhancement_factor, dominant_path,
           intact_path, defect_path, alpha_intact, alpha_defect
    """
    if not 0 <= defect_area_fraction <= 1:
        raise ValueError(f"defect_area_fraction must be in [0, 1], got {defect_area_fraction}")

    alpha_intact = D_ox * K_ox / L_ox

    if defect_type == 'pinhole':
        alpha_defect = np.inf
    elif defect_type == 'crack':
        alpha_defect = alpha_intact / thickness_factor
    elif defect_type == 'grain_boundary':
        alpha_defect = diffusivity_factor * alpha_intact
    else:
        raise ValueError(f"Unknown defect_type: {defect_type}")

    intact_result = calculate_path_flux_L6(
        P_up, P_down, L_m, k_diss, K_eq, alpha_intact, D_m, K_s_m, path_type='intact'
    )
    defect_result = calculate_path_flux_L6(
        P_up, P_down, L_m, k_diss, K_eq, alpha_defect, D_m, K_s_m,
        path_type=defect_type, k_diss_metal=k_diss_metal, K_eq_metal=K_eq_metal,
    )

    fraction_intact = 1.0 - defect_area_fraction
    J_intact        = intact_result['flux']
    J_defect        = defect_result['flux']
    J_total         = fraction_intact * J_intact + defect_area_fraction * J_defect

    enhancement_factor = J_total / J_intact if J_intact > 0 else np.inf
    flux_from_intact   = fraction_intact * J_intact
    flux_from_defect   = defect_area_fraction * J_defect

    if flux_from_defect > 0.5 * J_total:
        dominant_path = 'defect'
    elif flux_from_intact > 0.5 * J_total:
        dominant_path = 'intact'
    else:
        dominant_path = 'mixed'

    return {
        'J_total':            J_total,
        'enhancement_factor': enhancement_factor,
        'dominant_path':      dominant_path,
        'fraction_intact':    fraction_intact,
        'fraction_defect':    defect_area_fraction,
        'flux_from_intact':   flux_from_intact,
        'flux_from_defect':   flux_from_defect,
        'fraction_from_defect': flux_from_defect / J_total if J_total > 0 else 0,
        'intact_path':  intact_result,
        'defect_path':  defect_result,
        'alpha_intact': alpha_intact,
        'alpha_defect': alpha_defect,
        'alpha_ratio':  alpha_defect / alpha_intact if alpha_defect != np.inf else np.inf,
        'defect_type':        defect_type,
        'thickness_factor':   thickness_factor if defect_type == 'crack' else None,
        'diffusivity_factor': diffusivity_factor if defect_type == 'grain_boundary' else None,
    }


def calculate_mixed_defect_flux_L6(
    P_up, P_down, L_m,
    k_diss, K_eq,
    D_ox, K_ox, L_ox,
    D_m, K_s_m,
    defect_config,
    k_diss_metal=None,
    K_eq_metal=None,
):
    """
    Flux through oxide with multiple simultaneous defect types (L3+L6).

        J_total = fraction_intact × J_intact + Σᵢ (fraction_i × J_i)

    Parameters
    ----------
    defect_config : dict — per-defect configuration, e.g.:
        {
            'pinhole':        {'area_fraction': 0.001},
            'crack':          {'area_fraction': 0.005, 'thickness_factor': 0.1},
            'grain_boundary': {'area_fraction': 0.02,  'diffusivity_factor': 100},
        }
    k_diss_metal, K_eq_metal : float|None — metal kinetics for pinhole path

    Returns
    -------
    dict : J_total, enhancement_factor, dominant_path, flux_breakdown,
           intact_path, defect_paths, system_rate_limiting
    """
    valid_defect_types  = ['pinhole', 'crack', 'grain_boundary']
    total_defect_fraction = 0.0

    for defect_type, config in defect_config.items():
        if defect_type not in valid_defect_types:
            raise ValueError(f"Unknown defect type: {defect_type}")
        total_defect_fraction += config.get('area_fraction', 0.0)

    if total_defect_fraction > 1.0:
        raise ValueError(f"Total defect fraction ({total_defect_fraction:.3f}) exceeds 1.0")

    fraction_intact = 1.0 - total_defect_fraction
    alpha_intact    = D_ox * K_ox / L_ox

    try:
        intact_result = calculate_path_flux_L6(
            P_up, P_down, L_m, k_diss, K_eq, alpha_intact, D_m, K_s_m, path_type='intact'
        )
        if 'error' in intact_result:
            return {'J_total': np.nan, 'error': f"Intact path failed: {intact_result['error']}"}
    except Exception as e:
        return {'J_total': np.nan, 'error': f"Intact path calculation failed: {e}"}

    defect_results = {}

    for defect_type, config in defect_config.items():
        fraction_defect = config.get('area_fraction', 0.0)
        if fraction_defect <= 0:
            continue

        if defect_type == 'pinhole':
            alpha_defect = np.inf
        elif defect_type == 'crack':
            alpha_defect = alpha_intact / config.get('thickness_factor', 0.1)
        elif defect_type == 'grain_boundary':
            alpha_defect = config.get('diffusivity_factor', 100.0) * alpha_intact

        try:
            path_result = calculate_path_flux_L6(
                P_up, P_down, L_m, k_diss, K_eq, alpha_defect, D_m, K_s_m,
                path_type=defect_type,
                k_diss_metal=k_diss_metal if defect_type == 'pinhole' else None,
                K_eq_metal=K_eq_metal   if defect_type == 'pinhole' else None,
            )
            if 'error' in path_result:
                continue
        except Exception:
            continue

        defect_results[defect_type] = {
            'area_fraction':     fraction_defect,
            'alpha':             alpha_defect,
            'alpha_ratio':       alpha_defect / alpha_intact if alpha_defect != np.inf else np.inf,
            'path_result':       path_result,
            'flux_contribution': fraction_defect * path_result['flux'],
        }

    J_intact_contribution = fraction_intact * intact_result['flux']
    J_total = J_intact_contribution + sum(d['flux_contribution'] for d in defect_results.values())

    flux_breakdown = {
        'intact': {
            'area_fraction':    fraction_intact,
            'flux':             intact_result['flux'],
            'contribution':     J_intact_contribution,
            'fraction_of_total': J_intact_contribution / J_total if J_total > 0 else 0,
        }
    }
    for defect_type, data in defect_results.items():
        flux_breakdown[defect_type] = {
            'area_fraction':    data['area_fraction'],
            'flux':             data['path_result']['flux'],
            'contribution':     data['flux_contribution'],
            'fraction_of_total': data['flux_contribution'] / J_total if J_total > 0 else 0,
        }

    max_contribution = J_intact_contribution
    dominant_path    = 'intact'
    for defect_type, data in defect_results.items():
        if data['flux_contribution'] > max_contribution:
            max_contribution = data['flux_contribution']
            dominant_path    = defect_type

    dominant_fraction = max_contribution / J_total if J_total > 0 else 0
    if dominant_fraction < 0.7:
        dominant_path = 'mixed'

    if dominant_path == 'intact':
        system_rate_limiting = intact_result['rate_limiting']
        dominant_resistances = intact_result['resistances']
    elif dominant_path in defect_results:
        system_rate_limiting = defect_results[dominant_path]['path_result']['rate_limiting']
        dominant_resistances = defect_results[dominant_path]['path_result']['resistances']
    else:
        system_rate_limiting = 'mixed (multiple paths)'
        dominant_resistances = None

    enhancement_factor = J_total / intact_result['flux'] if intact_result['flux'] > 0 else np.inf

    return {
        'J_total':              J_total,
        'enhancement_factor':   enhancement_factor,
        'dominant_path':        dominant_path,
        'dominant_fraction':    dominant_fraction,
        'flux_breakdown':       flux_breakdown,
        'fraction_intact':      fraction_intact,
        'total_defect_fraction': total_defect_fraction,
        'intact_path':          intact_result,
        'defect_paths':         defect_results,
        'system_rate_limiting': system_rate_limiting,
        'dominant_resistances': dominant_resistances,
        'alpha_intact':         alpha_intact,
        'defect_config':        defect_config,
    }


# =============================================================================
# SECTION 7: L3+L4+L6 — Defective Oxide + Defective Metal + Surface Kinetics
# =============================================================================

def calculate_path_flux_L346_v2(
    P_up, P_down, L_m, temperature,
    k_diss, K_eq,
    alpha,
    D_lattice, K_s_m,
    microstructure_params,
    path_type='intact',
    lattice_density=1.06e29,
    method='average',
    n_points=20,
    mode='both',
    max_iterations=10,
    tolerance=1e-5,
    k_diss_metal=None,
    K_eq_metal=None,
):
    """
    Flux through a single oxide path with defective metal (L3+L4+L6 building block).

    Combines L6 surface kinetics, variable α for oxide path type, and L4
    microstructure effects in the metal via iterative D_eff convergence.

    Parameters
    ----------
    P_up, P_down         : float      — upstream/downstream H2 pressure [Pa]
    L_m                  : float      — metal thickness [m]
    temperature          : float      — temperature [K]
    k_diss, K_eq         : float      — oxide surface kinetics
    alpha                : float      — oxide permeance (np.inf for pinhole)
    D_lattice, K_s_m     : float      — metal lattice properties
    microstructure_params : dict      — grain size, traps, etc.
    path_type            : str        — 'intact', 'pinhole', 'crack', 'grain_boundary'
    lattice_density      : float      — lattice site density [sites/m³]
    method               : str        — D_eff averaging method
    n_points             : int        — profile discretisation points
    mode                 : str        — microstructure mode
    max_iterations       : int        — D_eff convergence limit
    tolerance            : float      — relative tolerance for D_eff convergence
    k_diss_metal, K_eq_metal : float|None — metal kinetics for pinhole path

    Returns
    -------
    dict : flux, theta, P_int, D_eff, modification_factor, profiles,
           resistances, rate_limiting, convergence
    """
    is_pinhole  = (path_type == 'pinhole' or alpha == np.inf or alpha > 1e10)
    sqrt_P_down = np.sqrt(max(P_down, 0))
    sqrt_P_up   = np.sqrt(max(P_up,   0))

    if is_pinhole:
        if k_diss_metal is not None and K_eq_metal is not None:
            k_diss_eff      = k_diss_metal
            K_eq_eff        = K_eq_metal
            use_sieverts_limit = False
        else:
            use_sieverts_limit = True
    else:
        k_diss_eff      = k_diss
        K_eq_eff        = K_eq
        use_sieverts_limit = False

    D_eff               = D_lattice
    convergence_history = []

    for iteration in range(max_iterations):
        beta = D_eff * K_s_m / L_m

        if is_pinhole and use_sieverts_limit:
            sqrt_P_int = sqrt_P_up
            theta_ss   = 0.0
        else:
            def residual(theta):
                if theta <= 0 or theta >= 1:
                    return np.inf
                J_surf_loc = surface_flux(theta, P_up, k_diss_eff, K_eq_eff)
                if is_pinhole:
                    sqrt_P_int_loc = g_theta(theta, K_eq_eff)
                else:
                    sqrt_P_int_loc = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
                J_metal_loc = beta * (sqrt_P_int_loc - sqrt_P_down)
                return J_surf_loc - J_metal_loc

            try:
                theta_ss = brentq(residual, 1e-12, 1.0 - 1e-12)
            except ValueError as e:
                return {
                    'flux': np.nan, 'theta': np.nan, 'P_int': np.nan, 'D_eff': np.nan,
                    'path_type': path_type,
                    'error': f'Failed to solve for theta: {str(e)}',
                    'iteration': iteration,
                }

            if is_pinhole:
                sqrt_P_int = g_theta(theta_ss, K_eq_eff)
            else:
                sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)

        x_array    = np.linspace(0, L_m, n_points)
        sqrt_P_arr = sqrt_P_int - (sqrt_P_int - sqrt_P_down) * x_array / L_m
        C_array    = K_s_m * sqrt_P_arr

        D_array          = np.zeros(n_points)
        theta_trap_array = np.zeros(n_points)
        gb_factor_array  = np.zeros(n_points)

        for i, C_local in enumerate(C_array):
            C_local = max(C_local, 1e-20)
            result_i = combined_microstructure_model(
                D_lattice=D_lattice,
                temperature=temperature,
                microstructure_params=microstructure_params,
                lattice_concentration=C_local,
                lattice_density=lattice_density,
                mode=mode,
            )
            D_array[i] = result_i['D_eff']
            if 'trapping' in result_i and result_i['trapping'] is not None:
                theta_trap_array[i] = result_i['trapping'].get('theta_total', 0.0)
            elif 'theta_total' in result_i:
                theta_trap_array[i] = result_i['theta_total']
            if 'gb_enhancement' in result_i and result_i['gb_enhancement'] is not None:
                gb_factor_array[i] = result_i['gb_enhancement'].get('factor', 1.0)
            elif 'gb_enhancement_factor' in result_i:
                gb_factor_array[i] = result_i['gb_enhancement_factor']
            else:
                gb_factor_array[i] = 1.0

        if method == 'harmonic':
            D_eff_new = len(D_array) / np.sum(1.0 / D_array)
        elif method == 'inlet':
            D_eff_new = D_array[0]
        elif method == 'outlet':
            D_eff_new = D_array[-1]
        else:
            D_eff_new = np.mean(D_array)

        rel_change = abs(D_eff_new - D_eff) / D_eff if D_eff > 0 else np.inf
        convergence_history.append({
            'iteration': iteration, 'D_eff': D_eff_new,
            'theta': theta_ss, 'relative_change': rel_change,
        })

        if rel_change < tolerance:
            D_eff = D_eff_new
            break
        D_eff = D_eff_new

    beta_eff = D_eff * K_s_m / L_m
    P_int    = sqrt_P_int**2
    J_metal  = beta_eff * (sqrt_P_int - sqrt_P_down)

    if is_pinhole and use_sieverts_limit:
        J_surface    = J_metal
        J_oxide      = np.nan
        R_oxide      = 0.0
        R_surface    = 0.0
        alpha_out    = np.inf
        kinetics_used = 'sieverts'
    elif is_pinhole:
        J_surface    = surface_flux(theta_ss, P_up, k_diss_eff, K_eq_eff)
        J_oxide      = np.nan
        R_oxide      = 0.0
        R_surface    = smooth_surface_resistance(k_diss_eff, P_up, theta_ss)
        alpha_out    = np.inf
        kinetics_used = 'metal'
    else:
        J_surface    = surface_flux(theta_ss, P_up, k_diss, K_eq)
        J_oxide      = oxide_flux(theta_ss, alpha, beta_eff, K_eq, P_down)
        R_oxide      = 1.0 / alpha
        R_surface    = smooth_surface_resistance(k_diss, P_up, theta_ss)
        alpha_out    = alpha
        kinetics_used = 'oxide'

    R_metal   = 1.0 / beta_eff
    R_total   = R_surface + R_oxide + R_metal

    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide   = R_oxide   / R_total if R_total > 0 else 0
    fraction_metal   = R_metal   / R_total if R_total > 0 else 0

    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'

    modification_factor = D_eff / D_lattice

    return {
        'flux':      J_metal,
        'theta':     theta_ss,
        'P_int':     P_int,
        'path_type': path_type,
        'kinetics_used': kinetics_used,
        'alpha':     alpha_out,
        'beta_lattice': D_lattice * K_s_m / L_m,
        'beta_eff':  beta_eff,
        'D_eff':     D_eff,
        'D_lattice': D_lattice,
        'modification_factor': modification_factor,
        'microstructure': {
            'avg_gb_factor':   np.mean(gb_factor_array),
            'avg_theta_trap':  np.mean(theta_trap_array),
        },
        'flux_balance': {
            'J_surface': J_surface, 'J_oxide': J_oxide, 'J_metal': J_metal,
        },
        'resistances': {
            'R_surface': R_surface, 'R_oxide': R_oxide, 'R_metal': R_metal,
            'R_total': R_total, 'fraction_surface': fraction_surface,
            'fraction_oxide': fraction_oxide, 'fraction_metal': fraction_metal,
        },
        'rate_limiting': rate_limiting,
        'convergence': {
            'iterations': len(convergence_history),
            'converged':  len(convergence_history) < max_iterations,
            'history':    convergence_history,
        },
        'profiles': {
            'x': x_array, 'D': D_array, 'C': C_array,
            'theta_trap': theta_trap_array, 'gb_factor': gb_factor_array,
        },
    }


def calculate_full_model_flux_L346_v2(
    P_up, P_down, L_m, temperature,
    k_diss, K_eq,
    D_ox, K_ox, L_ox,
    D_lattice, K_s_m,
    microstructure_params,
    defect_config,
    lattice_density=1.06e29,
    method='average',
    n_points=20,
    mode='both',
    k_diss_metal=None,
    K_eq_metal=None,
):
    """
    Full L3+L4+L6 model: defective oxide + defective metal + surface kinetics.

        J_total = fraction_intact × J_intact + Σᵢ (fraction_i × J_i)

    Parameters
    ----------
    P_up, P_down         : float — upstream/downstream H2 pressure [Pa]
    L_m                  : float — metal thickness [m]
    temperature          : float — temperature [K]
    k_diss, K_eq         : float — oxide surface kinetics
    D_ox, K_ox, L_ox     : float — intact oxide properties
    D_lattice, K_s_m     : float — metal lattice properties
    microstructure_params : dict — grain size, traps, etc.
    defect_config        : dict — oxide defect configuration (see below)
    k_diss_metal, K_eq_metal : float|None — metal surface kinetics for pinhole

    defect_config example
    ---------------------
    {
        'pinhole':        {'area_fraction': 0.001},
        'crack':          {'area_fraction': 0.005, 'thickness_factor': 0.1},
        'grain_boundary': {'area_fraction': 0.02,  'diffusivity_factor': 100},
    }

    Returns
    -------
    dict : J_total, enhancement_vs_intact, dominant_path, flux_breakdown,
           D_eff_avg, system_rate_limiting, alpha_intact
    """
    valid_defect_types    = ['pinhole', 'crack', 'grain_boundary']
    total_defect_fraction = 0.0

    for defect_type, config in defect_config.items():
        if defect_type not in valid_defect_types:
            raise ValueError(f"Unknown defect type: {defect_type}")
        total_defect_fraction += config.get('area_fraction', 0.0)

    if total_defect_fraction > 1.0:
        raise ValueError(f"Total defect fraction ({total_defect_fraction:.3f}) exceeds 1.0")

    fraction_intact = 1.0 - total_defect_fraction
    alpha_intact    = D_ox * K_ox / L_ox

    intact_result = calculate_path_flux_L346_v2(
        P_up, P_down, L_m, temperature,
        k_diss, K_eq, alpha_intact,
        D_lattice, K_s_m, microstructure_params,
        path_type='intact',
        lattice_density=lattice_density, method=method,
        n_points=n_points, mode=mode,
    )

    defect_results = {}

    for defect_type, config in defect_config.items():
        fraction_defect = config.get('area_fraction', 0.0)
        if fraction_defect <= 0:
            continue

        if defect_type == 'pinhole':
            alpha_defect = np.inf
        elif defect_type == 'crack':
            alpha_defect = alpha_intact / config.get('thickness_factor', 0.1)
        elif defect_type == 'grain_boundary':
            alpha_defect = config.get('diffusivity_factor', 100.0) * alpha_intact

        path_result = calculate_path_flux_L346_v2(
            P_up, P_down, L_m, temperature,
            k_diss, K_eq, alpha_defect,
            D_lattice, K_s_m, microstructure_params,
            path_type=defect_type,
            lattice_density=lattice_density, method=method,
            n_points=n_points, mode=mode,
            k_diss_metal=k_diss_metal if defect_type == 'pinhole' else None,
            K_eq_metal=K_eq_metal   if defect_type == 'pinhole' else None,
        )

        defect_results[defect_type] = {
            'area_fraction':     fraction_defect,
            'alpha':             alpha_defect,
            'alpha_ratio':       alpha_defect / alpha_intact if alpha_defect != np.inf else np.inf,
            'path_result':       path_result,
            'flux_contribution': fraction_defect * path_result['flux'],
        }

    J_intact_contribution = fraction_intact * intact_result['flux']
    J_total = J_intact_contribution + sum(d['flux_contribution'] for d in defect_results.values())

    flux_breakdown = {
        'intact': {
            'area_fraction':    fraction_intact,
            'flux':             intact_result['flux'],
            'contribution':     J_intact_contribution,
            'fraction_of_total': J_intact_contribution / J_total if J_total > 0 else 0,
            'D_eff':            intact_result['D_eff'],
            'modification_factor': intact_result['modification_factor'],
            'theta':            intact_result['theta'],
            'P_int':            intact_result['P_int'],
        }
    }
    for defect_type, data in defect_results.items():
        pr = data['path_result']
        flux_breakdown[defect_type] = {
            'area_fraction':    data['area_fraction'],
            'flux':             pr['flux'],
            'contribution':     data['flux_contribution'],
            'fraction_of_total': data['flux_contribution'] / J_total if J_total > 0 else 0,
            'D_eff':            pr['D_eff'],
            'modification_factor': pr['modification_factor'],
            'theta':            pr['theta'],
            'P_int':            pr['P_int'],
            'kinetics_used':    pr.get('kinetics_used', 'oxide'),
        }

    max_contribution = J_intact_contribution
    dominant_path    = 'intact'
    for defect_type, data in defect_results.items():
        if data['flux_contribution'] > max_contribution:
            max_contribution = data['flux_contribution']
            dominant_path    = defect_type

    dominant_fraction = max_contribution / J_total if J_total > 0 else 0
    if dominant_fraction < 0.7:
        dominant_path = 'mixed'

    if dominant_path == 'intact':
        system_rate_limiting = intact_result['rate_limiting']
        system_resistances   = intact_result['resistances']
    elif dominant_path in defect_results:
        system_rate_limiting = defect_results[dominant_path]['path_result']['rate_limiting']
        system_resistances   = defect_results[dominant_path]['path_result']['resistances']
    else:
        system_rate_limiting = 'mixed'
        system_resistances   = None

    D_eff_avg = sum(
        data['fraction_of_total'] * data['D_eff']
        for data in flux_breakdown.values()
    )

    # ── Flux-weighted layer fractions (true system-level view) ────────────────
    all_path_data = [
        (J_intact_contribution / J_total if J_total > 0 else 0,
        intact_result['resistances'])
    ]
    for data in defect_results.values():
        all_path_data.append((
            data['flux_contribution'] / J_total if J_total > 0 else 0,
            data['path_result']['resistances']
        ))

    fw_resistances = {
        'fraction_surface': sum(w * rs.get('fraction_surface', 0) for w, rs in all_path_data),
        'fraction_oxide':   sum(w * rs.get('fraction_oxide',   0) for w, rs in all_path_data),
        'fraction_metal':   sum(w * rs.get('fraction_metal',   0) for w, rs in all_path_data),
    }

    return {
        'J_total':              J_total,
        'enhancement_vs_intact': J_total / intact_result['flux'] if intact_result['flux'] > 0 else np.inf,
        'dominant_path':        dominant_path,
        'dominant_fraction':    dominant_fraction,
        'flux_breakdown':       flux_breakdown,
        'fraction_intact':      fraction_intact,
        'total_defect_fraction': total_defect_fraction,
        'intact_path':          intact_result,
        'defect_paths':         defect_results,
        'D_eff_avg':            D_eff_avg,
        'D_lattice':            D_lattice,
        'overall_modification_factor': D_eff_avg / D_lattice,
        'system_rate_limiting': system_rate_limiting,
        'system_resistances':   system_resistances,
        'alpha_intact':         alpha_intact,
        'defect_config':        defect_config,
        'microstructure_params': microstructure_params,
        'flux_weighted_resistances': fw_resistances
    }
