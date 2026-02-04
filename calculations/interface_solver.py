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


##############################################################################
##############################################################################
##############################################################################

# =============================================================================
# LEVEL 6: Interface Solvers with Surface Kinetics
# =============================================================================

def calculate_metal_flux_with_surface(D_metal, K_s_metal, thickness,
                                                P_interface, P_downstream,
                                                temperature, k_diss=None,
                                                k_recomb=None, material_name=None,
                                                coverage_mode='equilibrium',
                                                forced_coverage=None):
    """
    Calculate flux through metal with surface kinetics at oxide-metal interface.
    
    Uses calculate_surface_limited_flux() from surface_kinetics.py.
    Dissociation occurs at P_interface (oxide-metal interface).
    
    Parameters
    ----------
    D_metal : float
        Diffusion coefficient in metal [m²/s]
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
    k_diss, k_recomb : float, optional
        Surface kinetics parameters
    material_name : str, optional
        Material name for kinetics lookup
    coverage_mode : str, optional
        Mode for surface coverage calculation:
        - 'equilibrium': Langmuir isotherm (adsorption-desorption equilibrium)
        - 'steady_state': Solve J_diss = J_bulk iteratively
        - 'forced': Use forced_coverage value
    forced_coverage : float, optional
        Required when coverage_mode='forced', surface coverage (0-1)
    
    Returns
    -------
    dict
        Contains flux, flux_sieverts, theta, damkohler, SRF, coverage_mode
    """
    from calculations.surface_kinetics import calculate_surface_limited_flux
    
    # Call L6 function with P_interface as the upstream pressure for metal
    result = calculate_surface_limited_flux(
        D=D_metal,
        K_s=K_s_metal,
        thickness=thickness,
        P_up=P_interface,
        P_down=P_downstream,
        temperature=temperature,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    
    return result


def flux_balance_equation_with_surface(P_interface, P_upstream, P_downstream,
                                        oxide_props, metal_props, temperature,
                                        k_diss, k_recomb, material_name,
                                        coverage_mode='equilibrium',
                                        forced_coverage=None):
    """
    Flux balance for oxide + metal with surface kinetics (L2+L6).
    
    Returns zero when oxide flux equals surface-limited metal flux.
    
    Parameters
    ----------
    coverage_mode : str
        'equilibrium', 'steady_state', or 'forced'
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    """
    # Oxide flux (molecular diffusion - same as L2)
    flux_oxide = molecular_diffusion_flux(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness'],
        P_upstream,
        P_interface
    )
    
    # Metal flux with surface kinetics (L6)
    metal_result = calculate_metal_flux_with_surface(
        D_metal=metal_props['D_metal'],
        K_s_metal=metal_props['K_s_metal'],
        thickness=metal_props['thickness'],
        P_interface=P_interface,
        P_downstream=P_downstream,
        temperature=temperature,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    
    return flux_oxide - metal_result['flux']


def solve_interface_pressure_with_surface(P_upstream, P_downstream, oxide_props,
                                           metal_props, temperature,
                                           k_diss=None, k_recomb=None,
                                           material_name=None, method='brentq',
                                           coverage_mode='equilibrium',
                                           forced_coverage=None):
    """
    Solve for interface pressure with Level 6 surface kinetics.
    
    This is the L2+L6 equivalent of solve_interface_pressure().
    
    Parameters
    ----------
    P_upstream : float
        Upstream pressure [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    oxide_props : dict
        Oxide layer properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    temperature : float
        Temperature [K]
    k_diss, k_recomb : float, optional
        Surface kinetics parameters
    material_name : str, optional
        Material name for kinetics lookup
    method : str
        Solver method ('brentq')
    coverage_mode : str
        Mode for surface coverage:
        - 'equilibrium': Langmuir isotherm
        - 'steady_state': Solve J_diss = J_bulk
        - 'forced': Use forced_coverage
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    dict
        Solution with P_interface, flux, theta, Da, SRF, coverage_mode
    """
    from calculations.surface_kinetics import calculate_damkohler_number, surface_equilibrium_coverage
    
    min_pressure = 1e-20
    
    if P_upstream <= min_pressure:
        return {
            'P_interface': P_downstream + min_pressure,
            'P_upstream': P_upstream,
            'P_downstream': P_downstream,
            'flux': 0,
            'converged': False,
            'SRF': 1.0,
            'coverage_mode': coverage_mode
        }
    
    # Get kinetics if not provided
    if k_diss is None or k_recomb is None:
        if material_name is None:
            raise ValueError("Must provide k_diss/k_recomb or material_name")
        from data.surface_kinetics_data import get_surface_kinetics
        kinetics = get_surface_kinetics(material_name, temperature)
        k_diss = kinetics['k_diss']
        k_recomb = kinetics['k_recomb']
    
    # Physical bounds
    P_min = max(P_downstream + P_upstream * 1e-10, min_pressure)
    P_max = P_upstream * (1 - 1e-10)
    
    if P_min >= P_max:
        P_interface = np.sqrt(max(P_upstream * P_downstream, min_pressure))
        converged = False
    else:
        try:
            f_min = flux_balance_equation_with_surface(
                P_min, P_upstream, P_downstream, oxide_props, metal_props,
                temperature, k_diss, k_recomb, material_name,
                coverage_mode, forced_coverage
            )
            f_max = flux_balance_equation_with_surface(
                P_max, P_upstream, P_downstream, oxide_props, metal_props,
                temperature, k_diss, k_recomb, material_name,
                coverage_mode, forced_coverage
            )
            
            if f_min * f_max > 0:
                P_interface = P_min
                converged = False
            else:
                P_interface = brentq(
                    flux_balance_equation_with_surface,
                    P_min, P_max,
                    args=(P_upstream, P_downstream, oxide_props, metal_props,
                          temperature, k_diss, k_recomb, material_name,
                          coverage_mode, forced_coverage),
                    xtol=1e-12, rtol=1e-12
                )
                converged = True
                
        except (ValueError, RuntimeError):
            P_interface = P_min
            converged = False
    
    # Calculate final flux at solution
    metal_result = calculate_metal_flux_with_surface(
        D_metal=metal_props['D_metal'],
        K_s_metal=metal_props['K_s_metal'],
        thickness=metal_props['thickness'],
        P_interface=P_interface,
        P_downstream=P_downstream,
        temperature=temperature,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    
    # Damköhler at interface
    Da_info = calculate_damkohler_number(
        k_diss, metal_props['D_metal'],
        metal_props['K_s_metal'], metal_props['thickness']
    )
    
    return {
        'P_interface': P_interface,
        'P_upstream': P_upstream,
        'P_downstream': P_downstream,
        'flux': metal_result['flux'],
        'flux_sieverts': metal_result['flux_sieverts'],
        'SRF': metal_result['SRF'],
        'theta': metal_result['theta'],
        'theta_equilibrium': metal_result.get('theta_equilibrium'),
        'Da': Da_info['Da'],
        'damkohler': Da_info,
        'coverage_mode': coverage_mode,
        'converged': converged,
        'P_interface_normalized': (P_interface - P_downstream) / (P_upstream - P_downstream) if P_upstream > P_downstream else 0
    }


def calculate_oxide_metal_system_with_surface(P_upstream, P_downstream, oxide_props,
                                               metal_props, temperature,
                                               k_diss=None, k_recomb=None,
                                               material_name=None,
                                               coverage_mode='equilibrium',
                                               forced_coverage=None):
    """
    Main function for oxide + metal + surface kinetics (L2+L6).
    
    Analogous to calculate_oxide_metal_system() but with L6.
    
    Parameters
    ----------
    P_upstream : float
        Upstream pressure [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    oxide_props : dict
        Oxide layer properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    temperature : float
        Temperature [K]
    k_diss, k_recomb : float, optional
        Surface kinetics parameters
    material_name : str, optional
        Material name for kinetics lookup
    coverage_mode : str
        Mode for surface coverage:
        - 'equilibrium': Langmuir isotherm
        - 'steady_state': Solve J_diss = J_bulk
        - 'forced': Use forced_coverage
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    dict
        Complete L2+L6 solution with regime classification
    """
    solution = solve_interface_pressure_with_surface(
        P_upstream, P_downstream, oxide_props, metal_props,
        temperature, k_diss, k_recomb, material_name,
        coverage_mode=coverage_mode, forced_coverage=forced_coverage
    )
    
    # Calculate resistances (same as L2)
    R_oxide = calculate_oxide_resistance(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness']
    )
    
    P_interface_for_resistance = max(solution['P_interface'], 1e-20)
    R_metal = calculate_metal_resistance(
        metal_props['D_metal'],
        metal_props['K_s_metal'],
        metal_props['thickness'],
        P_interface_for_resistance
    )
    
    # Regime with surface effects
    ratio = R_oxide / R_metal
    Da = solution['Da']
    SRF = solution['SRF']
    
    # Determine dominant limitation
    if Da < 0.1 and SRF < 0.5:
        regime = "surface-limited"
    elif ratio > 10:
        regime = "oxide-limited"
    elif ratio < 0.1:
        regime = "metal-limited"
    else:
        regime = "transition"
    
    solution.update({
        'R_oxide': R_oxide,
        'R_metal': R_metal,
        'resistance_ratio': ratio,
        'regime': regime,
        'temperature': temperature
    })
    
    return solution


# =============================================================================
# LEVEL 5+6: Oxide + Defective Metal + Surface Kinetics (COUPLED)
# =============================================================================

def calculate_defective_metal_flux_sieverts_with_surface(D_lattice, K_s_metal, thickness,
                                                          P_interface, P_downstream,
                                                          temperature, microstructure_params,
                                                          k_diss=None, k_recomb=None,
                                                          material_name=None,
                                                          lattice_density=1.06e29,
                                                          method='average', n_points=10,
                                                          mode='both',
                                                          coverage_mode='equilibrium',
                                                          forced_coverage=None):
    """
    Calculate flux through defective metal with surface kinetics (L4+L6 for interface use).
    
    This is a wrapper that calls calculate_defective_metal_flux_with_surface()
    from permeation_calc.py, designed for use in interface pressure solving.
    
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
        Microstructure specification
    k_diss, k_recomb : float, optional
        Surface kinetics parameters
    material_name : str, optional
        Material name for kinetics lookup
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    mode : str, optional
        'both', 'gb_only', 'trapping_only'
    coverage_mode : str, optional
        Mode for surface coverage:
        - 'equilibrium': Langmuir isotherm
        - 'steady_state': Solve J_diss = J_bulk
        - 'forced': Use forced_coverage
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    dict
        Contains flux, flux_sieverts, theta, Da, SRF, D_eff, coverage_mode
    """
    from calculations.permeation_calc import calculate_defective_metal_flux_with_surface
    
    result = calculate_defective_metal_flux_with_surface(
        D_lattice=D_lattice,
        K_s=K_s_metal,
        thickness=thickness,
        P_up=P_interface,
        P_down=P_downstream,
        temperature=temperature,
        microstructure_params=microstructure_params,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    
    return result


def flux_balance_equation_defective_metal_with_surface(P_interface, P_upstream, P_downstream,
                                                        oxide_props, metal_props, temperature,
                                                        microstructure_params, k_diss, k_recomb,
                                                        material_name, lattice_density=1.06e29,
                                                        method='average', n_points=10, mode='both',
                                                        coverage_mode='equilibrium',
                                                        forced_coverage=None):
    """
    Flux balance equation for oxide + defective metal + surface kinetics (L5+L6).
    
    This is the COUPLED flux balance where surface kinetics is included in the
    metal flux calculation, not applied as a post-processing step.
    
    At steady state: J_oxide = J_metal_with_surface
    
    This function returns the residual (J_oxide - J_metal_with_surface) which
    equals zero at the correct interface pressure.
    
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
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification
    k_diss, k_recomb : float
        Surface kinetics parameters
    material_name : str
        Material name (can be None if k_diss/k_recomb provided)
    lattice_density : float
        Lattice site density [m⁻³]
    method : str
        D_eff averaging method
    n_points : int
        Points for integration
    mode : str
        'both', 'gb_only', 'trapping_only'
    coverage_mode : str
        'equilibrium', 'steady_state', or 'forced'
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    float
        Residual (J_oxide - J_metal_with_surface), zero at solution
    """
    # Oxide flux (molecular diffusion through oxide)
    flux_oxide = molecular_diffusion_flux(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness'],
        P_upstream,
        P_interface
    )
    
    # Defective metal flux WITH surface kinetics (L4+L6)
    # Surface dissociation occurs at P_interface (oxide-metal interface)
    metal_result = calculate_defective_metal_flux_sieverts_with_surface(
        D_lattice=metal_props['D_metal'],
        K_s_metal=metal_props['K_s_metal'],
        thickness=metal_props['thickness'],
        P_interface=P_interface,
        P_downstream=P_downstream,
        temperature=temperature,
        microstructure_params=microstructure_params,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    
    return flux_oxide - metal_result['flux']


def solve_interface_pressure_defective_metal_with_surface(P_upstream, P_downstream, oxide_props,
                                                           metal_props, temperature, microstructure_params,
                                                           k_diss=None, k_recomb=None, material_name=None,
                                                           lattice_density=1.06e29, method='average',
                                                           n_points=10, mode='both',
                                                           max_iterations=10, tolerance=1e-6,
                                                           coverage_mode='equilibrium',
                                                           forced_coverage=None):
    """
    Solve for interface pressure with L5+L6: defective metal + surface kinetics (COUPLED).
    
    This is the properly COUPLED solver where surface kinetics is included in the
    flux balance equation, not applied as a post-processing step.
    
    Physics:
    --------
    At steady state, flux continuity requires:
        J_oxide(P_up, P_int) = J_metal_with_surface(P_int, P_down, θ)
    
    The key difference from the decoupled approach:
    - Decoupled (WRONG): Solve P_interface ignoring surface, then apply surface separately
    - Coupled (CORRECT): Solve P_interface where J_oxide = J_metal_with_surface
    
    This matters because:
    1. Surface kinetics affects the effective "resistance" of the metal
    2. This changes where P_interface settles
    3. Which in turn affects surface coverage θ
    4. The decoupled approach calculates θ at the wrong pressure
    
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
    k_diss, k_recomb : float, optional
        Surface kinetics parameters [mol/m²/s/Pa] and [m⁴/mol/s]
    material_name : str, optional
        Material name for kinetics lookup (alternative to k_diss/k_recomb)
    lattice_density : float, optional
        Lattice site density [m⁻³] (default 1.06e29 for FCC)
    method : str, optional
        D_eff averaging method ('average', 'harmonic', 'inlet', 'outlet')
    n_points : int, optional
        Points for D_eff integration
    mode : str, optional
        Microstructure mode: 'both', 'gb_only', 'trapping_only'
    max_iterations : int, optional
        Max iterations (kept for API compatibility, not used in brentq)
    tolerance : float, optional
        Convergence tolerance (kept for API compatibility)
    coverage_mode : str
        Mode for surface coverage:
        - 'equilibrium': Langmuir isotherm (fast adsorption-desorption)
        - 'steady_state': Solve J_diss = J_bulk iteratively
        - 'forced': Use forced_coverage value
    forced_coverage : float, optional
        Required when coverage_mode='forced', surface coverage (0-1)
    
    Returns
    -------
    dict
        Complete solution dictionary:
        
        Standard interface solver outputs:
        - 'P_interface': Solved interface pressure [Pa]
        - 'P_upstream': Input upstream pressure [Pa]
        - 'P_downstream': Input downstream pressure [Pa]
        - 'flux': Steady-state flux [mol/m²/s]
        - 'flux_error': Relative flux mismatch [-]
        - 'converged': True if solver converged
        - 'P_interface_normalized': (P_int - P_down)/(P_up - P_down)
        
        Level 4 outputs (defective metal):
        - 'D_eff': Effective metal diffusivity [m²/s]
        - 'D_lattice': Input lattice diffusivity [m²/s]
        - 'modification_factor': D_eff/D_lattice [-]
        - 'microstructure_details': Dict with GB and trapping info
        
        Level 6 outputs (surface kinetics):
        - 'flux_sieverts': Flux without surface limitation [mol/m²/s]
        - 'theta': Surface coverage at P_interface [-]
        - 'theta_equilibrium': Equilibrium coverage (if steady_state mode)
        - 'Da': Damköhler number [-]
        - 'damkohler': Full Damköhler info dict
        - 'SRF': Surface Reduction Factor [-]
        - 'coverage_mode': Mode used for calculation
        
        Kinetics reference:
        - 'k_diss': Dissociation rate constant used
        - 'k_recomb': Recombination rate constant used
        - 'temperature': Temperature [K]
    """
    from calculations.surface_kinetics import calculate_damkohler_number
    from calculations.permeation_calc import calculate_defective_metal_flux_with_surface
    
    min_pressure = 1e-20
    
    # Handle edge case: very low upstream pressure
    if P_upstream <= min_pressure:
        return {
            'P_interface': P_downstream + min_pressure,
            'P_upstream': P_upstream,
            'P_downstream': P_downstream,
            'flux': 0,
            'flux_sieverts': 0,
            'flux_error': 0,
            'converged': False,
            'P_interface_normalized': 0,
            'D_eff': metal_props['D_metal'],
            'D_lattice': metal_props['D_metal'],
            'modification_factor': 1.0,
            'theta': 0,
            'theta_equilibrium': 0,
            'Da': float('inf'),
            'damkohler': {'Da': float('inf')},
            'SRF': 1.0,
            'coverage_mode': coverage_mode,
            'microstructure_details': {},
            'k_diss': k_diss,
            'k_recomb': k_recomb,
            'temperature': temperature
        }
    
    # Get kinetics parameters if not provided
    if k_diss is None or k_recomb is None:
        if material_name is None:
            raise ValueError("Must provide k_diss/k_recomb or material_name")
        from data.surface_kinetics_data import get_surface_kinetics
        kinetics = get_surface_kinetics(material_name, temperature)
        k_diss = kinetics['k_diss']
        k_recomb = kinetics['k_recomb']
    
    # Physical bounds for P_interface
    # Must be between downstream and upstream pressures
    P_min = max(P_downstream + P_upstream * 1e-10, min_pressure)
    P_max = P_upstream * (1 - 1e-10)
    
    if P_min >= P_max:
        # Invalid bounds - use geometric mean as fallback
        P_interface = np.sqrt(max(P_upstream * P_downstream, min_pressure))
        converged = False
    else:
        try:
            # Check if function has different signs at boundaries
            # This ensures a root exists in the interval
            f_min = flux_balance_equation_defective_metal_with_surface(
                P_min, P_upstream, P_downstream, oxide_props, metal_props,
                temperature, microstructure_params, k_diss, k_recomb, material_name,
                lattice_density, method, n_points, mode, coverage_mode, forced_coverage
            )
            f_max = flux_balance_equation_defective_metal_with_surface(
                P_max, P_upstream, P_downstream, oxide_props, metal_props,
                temperature, microstructure_params, k_diss, k_recomb, material_name,
                lattice_density, method, n_points, mode, coverage_mode, forced_coverage
            )
            
            if f_min * f_max > 0:
                # Same sign at boundaries - no root in interval
                # This typically means oxide is completely dominant
                P_interface = P_min
                converged = False
            else:
                # Use Brent's method to find root
                P_interface = brentq(
                    flux_balance_equation_defective_metal_with_surface,
                    P_min, P_max,
                    args=(P_upstream, P_downstream, oxide_props, metal_props,
                          temperature, microstructure_params, k_diss, k_recomb,
                          material_name, lattice_density, method, n_points, mode,
                          coverage_mode, forced_coverage),
                    xtol=1e-12, rtol=1e-12
                )
                converged = True
                
        except (ValueError, RuntimeError):
            # Solver failed - use minimum as fallback
            P_interface = P_min
            converged = False
    
    # Calculate all outputs at the solved interface pressure
    metal_result = calculate_defective_metal_flux_with_surface(
        D_lattice=metal_props['D_metal'],
        K_s=metal_props['K_s_metal'],
        thickness=metal_props['thickness'],
        P_up=P_interface,
        P_down=P_downstream,
        temperature=temperature,
        microstructure_params=microstructure_params,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    
    # The system flux equals the oxide flux (which equals metal flux at solution)
    flux = molecular_diffusion_flux(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness'],
        P_upstream,
        P_interface
    )
    
    # Verify flux continuity
    flux_metal = metal_result['flux']
    flux_error = abs(flux - flux_metal) / flux if flux > 0 else 0
    
    # Calculate Damköhler number at interface
    D_eff = metal_result['D_eff']
    Da_info = calculate_damkohler_number(
        k_diss, D_eff, metal_props['K_s_metal'], metal_props['thickness']
    )
    
    # Normalized interface position
    if P_upstream > P_downstream:
        P_interface_normalized = (P_interface - P_downstream) / (P_upstream - P_downstream)
    else:
        P_interface_normalized = 0
    
    return {
        # Standard interface solver outputs
        'P_interface': P_interface,
        'P_upstream': P_upstream,
        'P_downstream': P_downstream,
        'flux': flux,
        'flux_error': flux_error,
        'converged': converged,
        'P_interface_normalized': P_interface_normalized,
        
        # Level 4 outputs (defective metal)
        'D_eff': D_eff,
        'D_lattice': metal_props['D_metal'],
        'modification_factor': D_eff / metal_props['D_metal'],
        'microstructure_details': metal_result.get('microstructure_details', {}),
        
        # Level 6 outputs (surface kinetics)
        'flux_sieverts': metal_result.get('flux_sieverts', flux),
        'theta': metal_result.get('theta', 0),
        'theta_equilibrium': metal_result.get('theta_equilibrium'),
        'Da': Da_info['Da'],
        'damkohler': Da_info,
        'SRF': metal_result.get('SRF', 1.0),
        'coverage_mode': coverage_mode,
        
        # Kinetics for reference
        'k_diss': k_diss,
        'k_recomb': k_recomb,
        'temperature': temperature
    }


def calculate_oxide_defective_metal_system_with_surface(P_upstream, P_downstream,
                                                         oxide_props, metal_props,
                                                         temperature, microstructure_params,
                                                         k_diss=None, k_recomb=None,
                                                         material_name=None,
                                                         lattice_density=1.06e29,
                                                         method='average', n_points=10,
                                                         max_iterations=10, tolerance=1e-6,
                                                         mode='both',
                                                         coverage_mode='equilibrium',
                                                         forced_coverage=None):
    """
    Complete oxide + defective metal + surface kinetics system (L5+L6) - COUPLED.
    
    This is the main function for the full L5+L6 model, analogous to
    calculate_oxide_defective_metal_system() but with properly coupled
    Level 6 surface kinetics.
    
    The coupling ensures that surface kinetics is included in the interface
    pressure solver, not applied as a post-processing step. This is physically
    correct because surface kinetics affects the effective metal resistance,
    which determines where P_interface settles.
    
    Parameters
    ----------
    P_upstream : float
        Upstream pressure [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    oxide_props : dict
        Oxide layer properties:
        - 'D_ox': Oxide diffusivity [m²/s]
        - 'K_ox': Henry's law constant [mol/m³/Pa]
        - 'thickness': Oxide thickness [m]
    metal_props : dict
        Metal properties:
        - 'D_metal': Lattice diffusivity (D_lattice) [m²/s]
        - 'K_s_metal': Sieverts' constant [mol/m³/Pa^0.5]
        - 'thickness': Metal thickness [m]
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification:
        - 'grain_size': Average grain diameter [m]
        - 'grain_shape': 'equiaxed', 'columnar', or 'planar'
        - 'gb_type': 'HAGB', 'LAGB', 'twin', or 'special'
        - 'trap_list': List of trap dictionaries (optional)
    k_diss, k_recomb : float, optional
        Surface kinetics parameters
    material_name : str, optional
        Material name for kinetics lookup
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    max_iterations : int, optional
        Max iterations for interface solver
    tolerance : float, optional
        Convergence tolerance
    mode : str, optional
        Microstructure mode: 'both', 'gb_only', 'trapping_only'
    coverage_mode : str, optional
        Mode for surface coverage:
        - 'equilibrium': Langmuir isotherm
        - 'steady_state': Solve J_diss = J_bulk
        - 'forced': Use forced_coverage
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    dict
        Complete L5+L6 solution including:
        
        All outputs from solve_interface_pressure_defective_metal_with_surface()
        plus:
        - 'R_oxide': Oxide resistance [m²·s/mol]
        - 'R_metal': Effective metal resistance [m²·s/mol]
        - 'resistance_ratio': R_oxide/R_metal [-]
        - 'regime': Operating regime classification
        - 'flux_L5': Flux without surface kinetics (for comparison)
    """
    # Use the COUPLED solver
    solution = solve_interface_pressure_defective_metal_with_surface(
        P_upstream, P_downstream, oxide_props, metal_props,
        temperature, microstructure_params,
        k_diss=k_diss, k_recomb=k_recomb, material_name=material_name,
        lattice_density=lattice_density, method=method, n_points=n_points,
        mode=mode, max_iterations=max_iterations, tolerance=tolerance,
        coverage_mode=coverage_mode, forced_coverage=forced_coverage
    )
    
    # Calculate resistances for regime identification
    R_oxide = calculate_oxide_resistance(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness']
    )
    
    P_interface_for_resistance = max(solution['P_interface'], 1e-20)
    D_eff = solution['D_eff']
    
    R_metal = calculate_metal_resistance(
        D_eff,  # Use effective diffusivity
        metal_props['K_s_metal'],
        metal_props['thickness'],
        P_interface_for_resistance
    )
    
    # Calculate L5 flux for comparison (without surface kinetics)
    # Call the full L5 system solver to get true system-level comparison
    result_L5_system = calculate_oxide_defective_metal_system(
        P_upstream=P_upstream,
        P_downstream=P_downstream,
        oxide_props=oxide_props,
        metal_props=metal_props,
        temperature=temperature,
        microstructure_params=microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode
    )
    flux_L5 = result_L5_system['flux']

    # Calculate SYSTEM-LEVEL SRF (correct definition)
    # SRF = flux(with surface) / flux(without surface)
    # This compares complete systems, not fluxes at same P_interface
    SRF_system = solution['flux'] / flux_L5 if flux_L5 > 0 else 1.0
      
    # Regime classification considering all three effects:
    # 1. Oxide vs metal resistance
    # 2. Surface kinetics (Da, SRF)
    # 3. Microstructure (modification_factor)
    ratio = R_oxide / R_metal if R_metal > 0 else float('inf')
    Da = solution['Da']
    SRF = SRF_system  # Use system-level SRF for regime classification
    
    
    # Determine dominant limitation
    if Da < 0.1 and SRF < 0.5:
        # Surface kinetics is rate-limiting
        regime = "surface-limited"
    elif ratio > 10:
        # Oxide is rate-limiting
        regime = "oxide-limited"
    elif ratio < 0.1:
        # Metal diffusion is rate-limiting
        regime = "metal-limited"
    else:
        # Mixed/transition regime
        regime = "transition"
    
    # Add additional info to solution
    solution.update({
        'R_oxide': R_oxide,
        'R_metal': R_metal,
        'resistance_ratio': ratio,
        'regime': regime,
        'flux_L5': flux_L5,
        'SRF': SRF_system 
    })
    
    return solution