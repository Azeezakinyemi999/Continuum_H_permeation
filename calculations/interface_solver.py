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

# Handle edge case where upstream pressure is too low
min_pressure = 1e-20  # Minimum meaningful pressure
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



# =============================================================================
# LEVEL 4: Interface Solvers for Defective Metal
# =============================================================================

def calculate_defective_metal_flux_sieverts(D_lattice, K_s_metal, thickness, 
                                            P_interface, P_downstream,
                                            temperature, microstructure_params,
                                            lattice_density=1.06e29,
                                            method='average', n_points=10, mode='both'):
    """
    Calculate flux through defective metal using Sieverts' law boundary conditions.
    
    This is the Level 4 equivalent of calculate_metal_flux_sieverts().
    Designed for use in iterative interface pressure solving.
    
    Parameters
    ----------
    D_lattice : float
        Intrinsic lattice diffusion coefficient [m²/s]
    K_s_metal : float
        Sieverts' constant [mol/m³/Pa^0.5]
    thickness : float
        Metal thickness [m]
    P_interface : float
        Pressure at oxide/metal interface [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification (grain_size, grain_shape, gb_type, trap_list)
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method ('average', 'harmonic', 'inlet', 'outlet')
    n_points : int, optional
        Points for integration
    
    Returns
    -------
    float
        Flux through defective metal [mol/m²/s]
    """
    from calculations.permeation_calc import calculate_defective_metal_flux
    
    result = calculate_defective_metal_flux(
        D_lattice=D_lattice,
        K_s=K_s_metal,
        thickness=thickness,
        P_up=P_interface,
        P_down=P_downstream,
        temperature=temperature,
        microstructure_params=microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode
    )
    
    return result['flux']


def flux_balance_equation_defective_metal(P_interface, P_upstream, P_downstream, 
                                          oxide_props, metal_props, temperature,
                                          microstructure_params, lattice_density=1.06e29,
                                          method='average', n_points=10, mode='both'):
    """
    Flux balance equation for oxide + defective metal system.
    
    Returns zero when oxide flux equals defective metal flux.
    Used by brentq solver.
    
    Parameters
    ----------
    P_interface : float
        Interface pressure to solve for [Pa]
    P_upstream : float
        Upstream pressure [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal/D_lattice, K_s_metal, thickness)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification
    lattice_density : float
        Lattice site density [m⁻³]
    method : str
        D_eff averaging method
    n_points : int
        Points for integration
    
    Returns
    -------
    float
        Flux difference (oxide - metal), should be zero at solution
    """
    from calculations.oxide_permeation import molecular_diffusion_flux
    
    # Oxide flux (molecular diffusion - same as Level 2)
    flux_oxide = molecular_diffusion_flux(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness'],
        P_upstream,
        P_interface
    )
    
    # Defective metal flux (Level 4)
    flux_metal = calculate_defective_metal_flux_sieverts(
        D_lattice=metal_props['D_metal'],  # This is D_lattice for Level 4
        K_s_metal=metal_props['K_s_metal'],
        thickness=metal_props['thickness'],
        P_interface=P_interface,
        P_downstream=P_downstream,
        temperature=temperature,
        microstructure_params=microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode
    )
    
    return flux_oxide - flux_metal


def solve_interface_pressure_defective_metal(P_upstream, P_downstream, oxide_props, 
                                             metal_props, temperature, microstructure_params,
                                             lattice_density=1.06e29, method='average',
                                             n_points=10, solver_method='brentq',
                                             max_iterations=10, tolerance=1e-6, mode='both'):
    """
    Solve for interface pressure with Level 4 defective metal (iterative).
    
    This is the Level 4 equivalent of solve_interface_pressure().
    Uses iterative approach since D_eff depends on concentration which
    depends on P_interface.
    
    Parameters
    ----------
    P_upstream : float
        Upstream pressure [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    oxide_props : dict
        Oxide layer properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification:
        - 'grain_size': Average grain diameter [m]
        - 'grain_shape': 'equiaxed', 'columnar', or 'planar'
        - 'gb_type': 'HAGB', 'LAGB', 'twin', or 'special'
        - 'trap_list': List of trap dictionaries
    lattice_density : float, optional
        Lattice site density [m⁻³] (default 1.06e29 for FCC)
    method : str, optional
        D_eff averaging method ('average', 'harmonic', 'inlet', 'outlet')
    n_points : int, optional
        Points for D_eff integration (default 10)
    solver_method : str, optional
        Root finding method ('brentq' or 'root_scalar')
    max_iterations : int, optional
        Maximum outer iterations for D_eff convergence (default 10)
    tolerance : float, optional
        Relative convergence tolerance for P_interface (default 1e-6)
    
    Returns
    -------
    dict
        Solution dictionary containing:
        - 'P_interface': Solved interface pressure [Pa]
        - 'P_upstream': Input upstream pressure [Pa]
        - 'P_downstream': Input downstream pressure [Pa]
        - 'flux': Steady-state flux [mol/m²/s]
        - 'flux_error': Relative flux mismatch [-]
        - 'converged': True if solver converged
        - 'P_interface_normalized': (P_int - P_down)/(P_up - P_down)
        
        Level 4 specific:
        - 'D_eff': Final effective metal diffusivity [m²/s]
        - 'D_lattice': Input lattice diffusivity [m²/s]
        - 'modification_factor': D_eff/D_lattice [-]
        - 'level4_iterations': Number of outer iterations
        - 'level4_converged': True if D_eff iteration converged
        - 'microstructure_details': Dict with GB and trapping info
    
    Notes
    -----
    The iterative approach works as follows:
    1. Initial guess: Use D_lattice to estimate P_interface
    2. Calculate D_eff at estimated concentration profile
    3. Re-solve for P_interface with updated D_eff
    4. Repeat until P_interface converges
    
    This accounts for the coupling between interface pressure and
    concentration-dependent effective diffusivity.
    """
    from calculations.oxide_permeation import molecular_diffusion_flux
    from calculations.permeation_calc import calculate_defective_metal_flux
    
    # Handle edge cases
    min_pressure = 1e-20
    
    if P_upstream <= min_pressure:
        return {
            'P_interface': P_downstream + min_pressure,
            'P_upstream': P_upstream,
            'P_downstream': P_downstream,
            'flux': 0.0,
            'flux_error': 0,
            'converged': False,
            'P_interface_normalized': 0,
            'D_eff': metal_props['D_metal'],
            'D_lattice': metal_props['D_metal'],
            'modification_factor': 1.0,
            'level4_iterations': 0,
            'level4_converged': False,
            'microstructure_details': {}
        }
    
    # Initial guess for P_interface (from Level 2 behavior)
    P_interface_old = P_upstream / 2
    D_eff_current = metal_props['D_metal']  # Start with lattice value
    
    # Outer iteration loop for Level 4 coupling
    for iteration in range(max_iterations):
        
        # Create effective metal props with current D_eff estimate
        effective_metal_props = {
            'D_metal': D_eff_current,
            'K_s_metal': metal_props['K_s_metal'],
            'thickness': metal_props['thickness']
        }
        
        # Solve for interface pressure using brentq
        P_min = max(P_downstream + P_upstream * 1e-10, min_pressure)
        P_max = P_upstream * (1 - 1e-10)
        
        if P_min >= P_max:
            P_interface_new = np.sqrt(max(P_upstream * P_downstream, min_pressure))
            converged = False
        else:
            try:
                # Check if function has different signs at boundaries
                f_min = flux_balance_equation_defective_metal(
                    P_min, P_upstream, P_downstream, oxide_props, metal_props,
                    temperature, microstructure_params, lattice_density, method, n_points, mode
                )
                f_max = flux_balance_equation_defective_metal(
                    P_max, P_upstream, P_downstream, oxide_props, metal_props,
                    temperature, microstructure_params, lattice_density, method, n_points, mode
                )
                
                if f_min * f_max > 0:
                    # Same sign - no solution in interval
                    P_interface_new = P_min
                    converged = False
                else:
                    P_interface_new = brentq(
                        flux_balance_equation_defective_metal,
                        P_min, P_max,
                        args=(P_upstream, P_downstream, oxide_props, metal_props,
                              temperature, microstructure_params, lattice_density, 
                              method, n_points, mode),
                        xtol=1e-12,
                        rtol=1e-12
                    )
                    converged = True
                    
            except (ValueError, RuntimeError):
                P_interface_new = P_min
                converged = False
        
        # Update D_eff based on new P_interface
        metal_result = calculate_defective_metal_flux(
            D_lattice=metal_props['D_metal'],
            K_s=metal_props['K_s_metal'],
            thickness=metal_props['thickness'],
            P_up=P_interface_new,
            P_down=P_downstream,
            temperature=temperature,
            microstructure_params=microstructure_params,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode
        )
        
        D_eff_new = metal_result['D_eff']
        
        # Check convergence of P_interface
        if P_interface_old > 0:
            rel_change = abs(P_interface_new - P_interface_old) / P_interface_old
        else:
            rel_change = abs(P_interface_new - P_interface_old)
        
        if rel_change < tolerance and iteration > 0:
            # Converged
            break
        
        # Update for next iteration
        P_interface_old = P_interface_new
        D_eff_current = D_eff_new
    
    level4_converged = (rel_change < tolerance) if iteration > 0 else True
    
    # Calculate final flux
    flux = molecular_diffusion_flux(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness'],
        P_upstream,
        P_interface_new
    )
    
    # Verify flux continuity
    flux_metal_check = metal_result['flux']
    flux_error = abs(flux - flux_metal_check) / flux if flux > 0 else 0
    
    # Calculate normalized position
    if P_upstream > P_downstream:
        P_interface_normalized = (P_interface_new - P_downstream) / (P_upstream - P_downstream)
    else:
        P_interface_normalized = 0
    
    return {
        # Standard interface solver outputs
        'P_interface': P_interface_new,
        'P_upstream': P_upstream,
        'P_downstream': P_downstream,
        'flux': flux,
        'flux_error': flux_error,
        'converged': converged,
        'P_interface_normalized': P_interface_normalized,
        
        # Level 4 specific outputs
        'D_eff': D_eff_new,
        'D_lattice': metal_props['D_metal'],
        'modification_factor': D_eff_new / metal_props['D_metal'],
        'level4_iterations': iteration + 1,
        'level4_converged': level4_converged,
        'microstructure_details': metal_result['microstructure_details'],
        
        # Temperature for reference
        'temperature': temperature
    }


def calculate_oxide_defective_metal_system(P_upstream, P_downstream, oxide_props, 
                                           metal_props, temperature, microstructure_params,
                                           lattice_density=1.06e29, method='average',
                                           n_points=10, max_iterations=10, tolerance=1e-6,mode='both'):
    """
    Main function to calculate flux through oxide + defective metal system.
    
    This is the Level 4 equivalent of calculate_oxide_metal_system().
    Includes regime identification and resistance calculations.
    
    Parameters
    ----------
    P_upstream : float
        Upstream pressure [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    max_iterations : int, optional
        Max iterations for D_eff convergence
    tolerance : float, optional
        Convergence tolerance
    
    Returns
    -------
    dict
        Complete system solution including:
        - All outputs from solve_interface_pressure_defective_metal()
        - 'R_oxide': Oxide resistance
        - 'R_metal': Effective metal resistance
        - 'resistance_ratio': R_oxide/R_metal
        - 'regime': Operating regime classification
    """
    from calculations.oxide_permeation import calculate_oxide_resistance, calculate_metal_resistance
    
    # Solve for interface pressure with Level 4
    solution = solve_interface_pressure_defective_metal(
        P_upstream, P_downstream, oxide_props, metal_props,
        temperature, microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        max_iterations=max_iterations,
        tolerance=tolerance,
        mode=mode
    )
    
    # Calculate resistances for regime identification
    R_oxide = calculate_oxide_resistance(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness']
    )
    
    # Use effective D_metal for resistance calculation
    P_interface_for_resistance = max(solution['P_interface'], 1e-20)
    
    R_metal = calculate_metal_resistance(
        solution['D_eff'],  # Use effective diffusivity!
        metal_props['K_s_metal'],
        metal_props['thickness'],
        P_interface_for_resistance
    )
    
    # Identify limiting mechanism
    ratio = R_oxide / R_metal if R_metal > 0 else float('inf')
    
    if ratio > 10:
        regime = "oxide_limited"
    elif ratio < 0.1:
        regime = "metal_limited"
    else:
        regime = "transition"
    
    # Add resistance info to solution
    solution.update({
        'R_oxide': R_oxide,
        'R_metal': R_metal,
        'resistance_ratio': ratio,
        'regime': regime
    })
    
    return solution
