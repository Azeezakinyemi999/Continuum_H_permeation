# # calculations/interface_solver.py

# import numpy as np
# from scipy.optimize import brentq, root_scalar
# from calculations.oxide_permeation import (
#     molecular_diffusion_flux,
#     calculate_oxide_resistance,
#     calculate_metal_resistance
# )

# def calculate_metal_flux_sieverts(D_metal, K_s_metal, thickness, P_interface, P_downstream):
#     """
#     Calculate flux through metal using Sieverts' law.
    
#     Parameters:
#     -----------
#     D_metal : float
#         Diffusion coefficient in metal (m²/s)
#     K_s_metal : float
#         Sieverts' constant (mol/m³/Pa^0.5)
#     thickness : float
#         Metal thickness (m)
#     P_interface : float
#         Pressure at oxide/metal interface (Pa)
#     P_downstream : float
#         Downstream pressure (Pa)
    
#     Returns:
#     --------
#     float
#         Flux through metal (mol/m²/s)
#     """
#     if P_interface < 0 or P_downstream < 0:
#         raise ValueError("Pressures must be non-negative")
    
#     C_interface = K_s_metal * np.sqrt(P_interface)
#     C_downstream = K_s_metal * np.sqrt(P_downstream)
    
#     flux = D_metal * (C_interface - C_downstream) / thickness
#     return flux


# def flux_balance_equation(P_interface, P_upstream, P_downstream, oxide_props, metal_props):
#     """
#     Flux balance equation that equals zero when fluxes match.
    
#     This is the key equation: flux_oxide - flux_metal = 0
    
#     Parameters:
#     -----------
#     P_interface : float
#         Interface pressure to solve for (Pa)
#     P_upstream : float
#         Upstream pressure (Pa)
#     P_downstream : float
#         Downstream pressure (Pa)
#     oxide_props : dict
#         Contains D_ox, K_ox, thickness
#     metal_props : dict
#         Contains D_metal, K_s_metal, thickness
    
#     Returns:
#     --------
#     float
#         Flux difference (should be zero at solution)
#     """
#     # Oxide flux (molecular diffusion)
#     flux_oxide = molecular_diffusion_flux(
#         oxide_props['D_ox'],
#         oxide_props['K_ox'],
#         oxide_props['thickness'],
#         P_upstream,
#         P_interface
#     )
    
#     # Metal flux (Sieverts' law)
#     flux_metal = calculate_metal_flux_sieverts(
#         metal_props['D_metal'],
#         metal_props['K_s_metal'],
#         metal_props['thickness'],
#         P_interface,
#         P_downstream
#     )
    
#     return flux_oxide - flux_metal


# def solve_interface_pressure(P_upstream, P_downstream, oxide_props, metal_props, method='brentq'):
#     """
#     Solve for interface pressure where oxide and metal fluxes match.
    
#     Parameters:
#     -----------
#     P_upstream : float
#         Upstream pressure (Pa)
#     P_downstream : float
#         Downstream pressure (Pa)
#     oxide_props : dict
#         Oxide layer properties
#     metal_props : dict
#         Metal layer properties
#     method : str
#         Solver method ('brentq' or 'root_scalar')
    
#     Returns:
#     --------
#     dict
#         Contains P_interface, flux, convergence info
#     """
#     # Physical bounds: P_interface must be between downstream and upstream
#     P_min = P_downstream + 1e-10 * P_upstream  # Slightly above downstream
#     P_max = P_upstream * (1 - 1e-10)  # Slightly below upstream
    
#     if P_min >= P_max:
#         raise ValueError(f"Invalid pressure bounds: P_min={P_min} >= P_max={P_max}")
    
#     try:
#         if method == 'brentq':
#             # Brent's method - robust and fast
#             P_interface = brentq(
#                 flux_balance_equation,
#                 P_min, P_max,
#                 args=(P_upstream, P_downstream, oxide_props, metal_props),
#                 xtol=1e-12,
#                 rtol=1e-12
#             )
#             converged = True
            
#         elif method == 'root_scalar':
#             # Alternative solver with more diagnostics
#             sol = root_scalar(
#                 flux_balance_equation,
#                 args=(P_upstream, P_downstream, oxide_props, metal_props),
#                 bracket=[P_min, P_max],
#                 method='brentq'
#             )
#             P_interface = sol.root
#             converged = sol.converged
#         else:
#             raise ValueError(f"Unknown method: {method}")
            
#     except ValueError as e:
#         print(f"Solver failed: {e}")
#         print(f"P_upstream={P_upstream:.2e}, P_downstream={P_downstream:.2e}")
#         # Return a reasonable guess if solver fails
#         P_interface = np.sqrt(P_upstream * P_downstream)  # Geometric mean
#         converged = False
    
#     # Calculate actual flux at solution
#     flux = molecular_diffusion_flux(
#         oxide_props['D_ox'],
#         oxide_props['K_ox'],
#         oxide_props['thickness'],
#         P_upstream,
#         P_interface
#     )
    
#     # Verify flux continuity
#     flux_metal_check = calculate_metal_flux_sieverts(
#         metal_props['D_metal'],
#         metal_props['K_s_metal'],
#         metal_props['thickness'],
#         P_interface,
#         P_downstream
#     )
    
#     flux_error = abs(flux - flux_metal_check) / flux if flux > 0 else 0
    
#     return {
#         'P_interface': P_interface,
#         'P_upstream': P_upstream,
#         'P_downstream': P_downstream,
#         'flux': flux,
#         'flux_error': flux_error,
#         'converged': converged,
#         'P_interface_normalized': (P_interface - P_downstream) / (P_upstream - P_downstream)
#     }


# def calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props, T_K=None):
#     """
#     Main function to calculate flux through oxide+metal system.
    
#     Parameters:
#     -----------
#     P_upstream : float
#         Upstream pressure (Pa)
#     P_downstream : float
#         Downstream pressure (Pa)
#     oxide_props : dict
#         Oxide properties (can include T-dependent parameters)
#     metal_props : dict
#         Metal properties (can include T-dependent parameters)
#     T_K : float, optional
#         Temperature in Kelvin (for T-dependent properties)
    
#     Returns:
#     --------
#     dict
#         Complete system solution including flux, pressures, regime
#     """
#     # Solve for interface pressure
#     solution = solve_interface_pressure(P_upstream, P_downstream, oxide_props, metal_props)
    
#     # Calculate resistances for regime identification
#     R_oxide = calculate_oxide_resistance(
#         oxide_props['D_ox'],
#         oxide_props['K_ox'],
#         oxide_props['thickness']
#     )
    
#     R_metal = calculate_metal_resistance(
#         metal_props['D_metal'],
#         metal_props['K_s_metal'],
#         metal_props['thickness'],
#         solution['P_interface']
#     )
    
#     # Identify limiting mechanism
#     ratio = R_oxide / R_metal
#     if ratio > 10:
#         regime = "oxide_limited"
#     elif ratio < 0.1:
#         regime = "metal_limited"
#     else:
#         regime = "transition"
    
#     # Add additional information
#     solution.update({
#         'R_oxide': R_oxide,
#         'R_metal': R_metal,
#         'resistance_ratio': ratio,
#         'regime': regime,
#         'temperature': T_K
#     })
    
#     return solution


# def calculate_concentration_profile(P_upstream, P_downstream, oxide_props, metal_props, n_points=100):
#     """
#     Calculate concentration profile through oxide and metal layers.
    
#     Parameters:
#     -----------
#     P_upstream, P_downstream : float
#         Boundary pressures (Pa)
#     oxide_props, metal_props : dict
#         Layer properties
#     n_points : int
#         Number of points in each layer
    
#     Returns:
#     --------
#     dict
#         Contains position and concentration arrays
#     """
#     # First solve for interface pressure
#     solution = solve_interface_pressure(P_upstream, P_downstream, oxide_props, metal_props)
#     P_interface = solution['P_interface']
    
#     # Oxide layer positions and concentrations
#     x_oxide = np.linspace(0, oxide_props['thickness'], n_points)
#     # Linear profile in oxide
#     C_oxide_up = oxide_props['K_ox'] * P_upstream
#     C_oxide_interface = oxide_props['K_ox'] * P_interface
#     C_oxide = C_oxide_up - (C_oxide_up - C_oxide_interface) * x_oxide / oxide_props['thickness']
    
#     # Metal layer positions and concentrations
#     x_metal_start = oxide_props['thickness']
#     x_metal = np.linspace(x_metal_start, 
#                           x_metal_start + metal_props['thickness'], 
#                           n_points)
#     # Square root profile in metal
#     C_metal_interface = metal_props['K_s_metal'] * np.sqrt(P_interface)
#     C_metal_down = metal_props['K_s_metal'] * np.sqrt(P_downstream)
    
#     # Linear interpolation of sqrt(P), then calculate C
#     x_normalized = (x_metal - x_metal_start) / metal_props['thickness']
#     sqrt_P = np.sqrt(P_interface) - (np.sqrt(P_interface) - np.sqrt(P_downstream)) * x_normalized
#     C_metal = metal_props['K_s_metal'] * sqrt_P
    
#     return {
#         'x_oxide': x_oxide,
#         'C_oxide': C_oxide,
#         'x_metal': x_metal,
#         'C_metal': C_metal,
#         'x_all': np.concatenate([x_oxide, x_metal]),
#         'C_all': np.concatenate([C_oxide, C_metal]),
#         'P_interface': P_interface,
#         'C_discontinuity': C_oxide_interface - C_metal_interface
#     }

# calculations/interface_solver.py

import numpy as np
from scipy.optimize import brentq, root_scalar
from calculations.oxide_permeation import (
    molecular_diffusion_flux,
    calculate_oxide_resistance,
    calculate_metal_resistance
)

def calculate_metal_flux_sieverts(D_metal, K_s_metal, thickness, P_interface, P_downstream):
    """
    Calculate flux through metal using Sieverts' law.
    
    Parameters:
    -----------
    D_metal : float
        Diffusion coefficient in metal (m²/s)
    K_s_metal : float
        Sieverts' constant (mol/m³/Pa^0.5)
    thickness : float
        Metal thickness (m)
    P_interface : float
        Pressure at oxide/metal interface (Pa)
    P_downstream : float
        Downstream pressure (Pa)
    
    Returns:
    --------
    float
        Flux through metal (mol/m²/s)
    """
    if P_interface < 0 or P_downstream < 0:
        raise ValueError("Pressures must be non-negative")
    
    C_interface = K_s_metal * np.sqrt(P_interface)
    C_downstream = K_s_metal * np.sqrt(P_downstream)
    
    flux = D_metal * (C_interface - C_downstream) / thickness
    return flux


def flux_balance_equation(P_interface, P_upstream, P_downstream, oxide_props, metal_props):
    """
    Flux balance equation that equals zero when fluxes match.
    
    This is the key equation: flux_oxide - flux_metal = 0
    
    Parameters:
    -----------
    P_interface : float
        Interface pressure to solve for (Pa)
    P_upstream : float
        Upstream pressure (Pa)
    P_downstream : float
        Downstream pressure (Pa)
    oxide_props : dict
        Contains D_ox, K_ox, thickness
    metal_props : dict
        Contains D_metal, K_s_metal, thickness
    
    Returns:
    --------
    float
        Flux difference (should be zero at solution)
    """
    # Oxide flux (molecular diffusion)
    flux_oxide = molecular_diffusion_flux(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness'],
        P_upstream,
        P_interface
    )
    
    # Metal flux (Sieverts' law)
    flux_metal = calculate_metal_flux_sieverts(
        metal_props['D_metal'],
        metal_props['K_s_metal'],
        metal_props['thickness'],
        P_interface,
        P_downstream
    )
    
    return flux_oxide - flux_metal


def solve_interface_pressure(P_upstream, P_downstream, oxide_props, metal_props, method='brentq'):
    """
    Solve for interface pressure where oxide and metal fluxes match.
    
    Parameters:
    -----------
    P_upstream : float
        Upstream pressure (Pa)
    P_downstream : float
        Downstream pressure (Pa)
    oxide_props : dict
        Oxide layer properties
    metal_props : dict
        Metal layer properties
    method : str
        Solver method ('brentq' or 'root_scalar')
    
    Returns:
    --------
    dict
        Contains P_interface, flux, convergence info
    """
    # Handle edge case where upstream pressure is too low
    min_pressure = 1e-20  # Minimum meaningful pressure
    
    if P_upstream <= min_pressure:
        # For extremely low pressures, interface pressure ≈ downstream pressure
        P_interface = P_downstream + min_pressure
        flux = molecular_diffusion_flux(
            oxide_props['D_ox'],
            oxide_props['K_ox'],
            oxide_props['thickness'],
            P_upstream,
            P_interface
        )
        
        return {
            'P_interface': P_interface,
            'P_upstream': P_upstream,
            'P_downstream': P_downstream,
            'flux': flux,
            'flux_error': 0,
            'converged': False,  # Not actually solved, just approximated
            'P_interface_normalized': 0
        }
    
    # Physical bounds: P_interface must be between downstream and upstream
    P_min = max(P_downstream + P_upstream * 1e-10, min_pressure)
    P_max = P_upstream * (1 - 1e-10)
    
    if P_min >= P_max:
        # If bounds are invalid, use geometric mean as guess
        P_interface = np.sqrt(max(P_upstream * P_downstream, min_pressure))
        flux = molecular_diffusion_flux(
            oxide_props['D_ox'],
            oxide_props['K_ox'],
            oxide_props['thickness'],
            P_upstream,
            P_interface
        )
        
        return {
            'P_interface': P_interface,
            'P_upstream': P_upstream,
            'P_downstream': P_downstream,
            'flux': flux,
            'flux_error': 0,
            'converged': False,
            'P_interface_normalized': (P_interface - P_downstream) / max(P_upstream - P_downstream, min_pressure)
        }
    
    try:
        # Check if the function has different signs at boundaries
        f_min = flux_balance_equation(P_min, P_upstream, P_downstream, oxide_props, metal_props)
        f_max = flux_balance_equation(P_max, P_upstream, P_downstream, oxide_props, metal_props)
        
        if f_min * f_max > 0:
            # Same sign at boundaries - no solution in this interval
            # This typically means oxide is completely dominant
            P_interface = P_min  # Interface pressure drops to minimum
            converged = False
        else:
            if method == 'brentq':
                # Brent's method - robust and fast
                P_interface = brentq(
                    flux_balance_equation,
                    P_min, P_max,
                    args=(P_upstream, P_downstream, oxide_props, metal_props),
                    xtol=1e-12,
                    rtol=1e-12
                )
                converged = True
                
            elif method == 'root_scalar':
                # Alternative solver with more diagnostics
                sol = root_scalar(
                    flux_balance_equation,
                    args=(P_upstream, P_downstream, oxide_props, metal_props),
                    bracket=[P_min, P_max],
                    method='brentq'
                )
                P_interface = sol.root
                converged = sol.converged
            else:
                raise ValueError(f"Unknown method: {method}")
            
    except (ValueError, RuntimeError) as e:
        # If solver fails, use approximation
        # For oxide-dominated case, interface pressure is very low
        P_interface = P_min
        converged = False
    
    # Calculate actual flux at solution
    flux = molecular_diffusion_flux(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness'],
        P_upstream,
        P_interface
    )
    
    # Verify flux continuity (only if we have a valid solution)
    if P_interface > 0:
        flux_metal_check = calculate_metal_flux_sieverts(
            metal_props['D_metal'],
            metal_props['K_s_metal'],
            metal_props['thickness'],
            P_interface,
            P_downstream
        )
        flux_error = abs(flux - flux_metal_check) / flux if flux > 0 else 0
    else:
        flux_error = 1.0  # Maximum error
    
    # Calculate normalized position safely
    if P_upstream > P_downstream:
        P_interface_normalized = (P_interface - P_downstream) / (P_upstream - P_downstream)
    else:
        P_interface_normalized = 0
    
    return {
        'P_interface': P_interface,
        'P_upstream': P_upstream,
        'P_downstream': P_downstream,
        'flux': flux,
        'flux_error': flux_error,
        'converged': converged,
        'P_interface_normalized': P_interface_normalized
    }


def calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props, T_K=None):
    """
    Main function to calculate flux through oxide+metal system.
    
    Parameters:
    -----------
    P_upstream : float
        Upstream pressure (Pa)
    P_downstream : float
        Downstream pressure (Pa)
    oxide_props : dict
        Oxide properties (can include T-dependent parameters)
    metal_props : dict
        Metal properties (can include T-dependent parameters)
    T_K : float, optional
        Temperature in Kelvin (for T-dependent properties)
    
    Returns:
    --------
    dict
        Complete system solution including flux, pressures, regime
    """
    # Solve for interface pressure
    solution = solve_interface_pressure(P_upstream, P_downstream, oxide_props, metal_props)
    
    # Calculate resistances for regime identification
    R_oxide = calculate_oxide_resistance(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness']
    )
    
    # Use a minimum pressure for resistance calculation to avoid errors
    P_interface_for_resistance = max(solution['P_interface'], min_pressure)
    
    R_metal = calculate_metal_resistance(
        metal_props['D_metal'],
        metal_props['K_s_metal'],
        metal_props['thickness'],
        P_interface_for_resistance
    )
    
    # Identify limiting mechanism
    ratio = R_oxide / R_metal
    if ratio > 10:
        regime = "oxide_limited"
    elif ratio < 0.5:
        regime = "metal_limited"
    else:
        regime = "transition"
    
    # Add additional information
    solution.update({
        'R_oxide': R_oxide,
        'R_metal': R_metal,
        'resistance_ratio': ratio,
        'regime': regime,
        'temperature': T_K
    })
    
    return solution


def calculate_concentration_profile(P_upstream, P_downstream, oxide_props, metal_props, n_points=100):
    """
    Calculate concentration profile through oxide and metal layers.
    
    Parameters:
    -----------
    P_upstream, P_downstream : float
        Boundary pressures (Pa)
    oxide_props, metal_props : dict
        Layer properties
    n_points : int
        Number of points in each layer
    
    Returns:
    --------
    dict
        Contains position and concentration arrays
    """
    # First solve for interface pressure
    solution = solve_interface_pressure(P_upstream, P_downstream, oxide_props, metal_props)
    P_interface = max(solution['P_interface'], min_pressure)  # Ensure positive for sqrt
    
    # Oxide layer positions and concentrations
    x_oxide = np.linspace(0, oxide_props['thickness'], n_points)
    # Linear profile in oxide
    C_oxide_up = oxide_props['K_ox'] * P_upstream
    C_oxide_interface = oxide_props['K_ox'] * P_interface
    C_oxide = C_oxide_up - (C_oxide_up - C_oxide_interface) * x_oxide / oxide_props['thickness']
    
    # Metal layer positions and concentrations
    x_metal_start = oxide_props['thickness']
    x_metal = np.linspace(x_metal_start, 
                          x_metal_start + metal_props['thickness'], 
                          n_points)
    # Square root profile in metal
    C_metal_interface = metal_props['K_s_metal'] * np.sqrt(P_interface)
    C_metal_down = metal_props['K_s_metal'] * np.sqrt(P_downstream)
    
    # Linear interpolation of sqrt(P), then calculate C
    x_normalized = (x_metal - x_metal_start) / metal_props['thickness']
    sqrt_P = np.sqrt(P_interface) - (np.sqrt(P_interface) - np.sqrt(P_downstream)) * x_normalized
    C_metal = metal_props['K_s_metal'] * sqrt_P
    
    return {
        'x_oxide': x_oxide,
        'C_oxide': C_oxide,
        'x_metal': x_metal,
        'C_metal': C_metal,
        'x_all': np.concatenate([x_oxide, x_metal]),
        'C_all': np.concatenate([C_oxide, C_metal]),
        'P_interface': P_interface,
        'C_discontinuity': C_oxide_interface - C_metal_interface
    }