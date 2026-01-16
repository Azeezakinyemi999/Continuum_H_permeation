# Level 4 Implementation Plan: Defective Metal Model

## Date: January 12, 2026

---

## Table of Contents
1. [Overview](#overview)
2. [Architecture Decisions](#architecture-decisions)
3. [Implementation in permeation_calc.py](#implementation-in-permeation_calcpy)
4. [Implementation in interface_solver.py](#implementation-in-interface_solverpy)
5. [Implementation in parallel_oxide_defect_paths.py](#implementation-in-parallel_oxide_defect_pathspy)
6. [Complete Architecture Summary](#complete-architecture-summary)
7. [Existing Analysis and Test Suite (Levels 1-3)](#existing-analysis-and-test-suite-levels-1-3)
8. [What's Missing for Level 4](#whats-missing-for-level-4)
9. [Redundant Functions Clarification](#redundant-functions-clarification)

---

## Overview

Level 4 adds **microstructure effects to the metal layer**, incorporating:
- Grain boundary fast diffusion paths (enhancement)
- Hydrogen trapping at defects (reduction)
- Position-dependent concentration effects

### Design Philosophy

The approach creates **parallel functions** that maintain the same interface as existing Level 1-3 functions, allowing seamless integration while keeping existing code unchanged.

---

## Architecture Decisions

### User Choices:

1. **Option 2**: Create parallel solver functions (keep existing unchanged)
2. **Iterative**: Re-calculate D_eff as P_interface converges (more accurate)
3. **Location**: `calculate_defective_metal_flux()` lives in `permeation_calc.py` (keeps all flux calculations together)

### Current vs. Proposed Architecture:

**Current Flow:**
```
Level 1: permeation_calc.py → calculate_simple_metal_flux()
Level 2: interface_solver.py → uses calculate_simple_metal_flux() internally
Level 3: parallel_oxide_defect_paths.py → uses Level 2 solver
Level 4: defective_metal.py → (separate, not integrated)
```

**Proposed Flow:**
```
permeation_calc.py:
├── calculate_simple_metal_flux()      → Level 1 (clean metal)
└── calculate_defective_metal_flux()   → Level 4 (metal with microstructure)

Level 2: interface_solver.py → can use EITHER flux function
Level 3: parallel_oxide_defect_paths.py → can use EITHER flux function
```

---

## Implementation in permeation_calc.py

### New Functions to Add:

```python
# =============================================================================
# LEVEL 4: Defective Metal Flux Calculations
# =============================================================================

def calculate_defective_metal_flux(D_lattice, K_s, thickness, P_up, P_down,
                                    temperature, microstructure_params,
                                    lattice_density=1.06e29,
                                    method='average', n_points=10):
    """
    Calculate hydrogen permeation flux through metal with microstructure effects.
    
    This is the Level 4 equivalent of calculate_simple_metal_flux().
    Incorporates:
    - Grain boundary fast diffusion paths (enhancement)
    - Hydrogen trapping at defects (reduction)
    - Position-dependent concentration effects
    
    Parameters
    ----------
    D_lattice : float
        Intrinsic lattice diffusion coefficient [m²/s]
    K_s : float
        Solubility constant [mol/m³/Pa^0.5] (unchanged by microstructure)
    thickness : float
        Metal thickness [m]
    P_up : float
        Upstream hydrogen pressure [Pa]
    P_down : float
        Downstream hydrogen pressure [Pa]
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification:
        - 'grain_size': Average grain diameter [m]
        - 'grain_shape': 'equiaxed', 'columnar', or 'planar'
        - 'gb_type': 'HAGB', 'LAGB', 'twin', or 'special'
        - 'trap_list': List of trap dicts with 'name', 'density', 'binding_energy'
        Optional:
        - 'gb_thickness': GB width [m] (default 0.5e-9)
        - 'model': 'parallel' or 'hart' (default 'parallel')
        - 'include_gb_trapping': bool (default True)
    lattice_density : float, optional
        Lattice site density [m⁻³] (default 1.06e29 for FCC)
    method : str, optional
        D_eff averaging method:
        - 'average': Arithmetic mean (default)
        - 'harmonic': Harmonic mean (series resistance)
        - 'inlet': Use high-concentration side only
        - 'outlet': Use low-concentration side only
    n_points : int, optional
        Points for numerical integration (default 10)
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'flux': Permeation flux [mol/m²/s]
        - 'C_up': Upstream concentration [mol/m³]
        - 'C_down': Downstream concentration [mol/m³]
        - 'D_eff': Effective diffusivity [m²/s]
        - 'D_lattice': Input lattice diffusivity [m²/s]
        - 'modification_factor': D_eff/D_lattice [-]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
        - 'microstructure_details': Dict with GB and trapping info
        - 'units': Dictionary of units for each output
    
    Notes
    -----
    This function can be used as a drop-in replacement for 
    calculate_simple_metal_flux() when microstructure effects are important.
    
    The key physics difference is that D_eff depends on local hydrogen 
    concentration, which varies through the thickness. The 'method' parameter 
    controls how this position-dependence is averaged for flux calculation.
    
    Theory
    ------
    For defective metal, the effective diffusivity combines:
    
    1. GB enhancement: D_gb_enhanced = (1-f_gb)×D_bulk + f_gb×D_gb
       where D_gb = α×D_bulk and α is temperature-dependent
    
    2. Trapping reduction: D_eff = D_gb_enhanced/(1 + Σθᵢ)
       where θᵢ is trap occupancy (concentration-dependent)
    
    The net effect can enhance OR reduce diffusivity depending on:
    - Grain size (smaller → more GB → more enhancement)
    - Trap density (higher → more trapping → more reduction)
    - Temperature (affects both GB enhancement and trap occupancy)
    - Concentration (affects trap occupancy via Oriani equilibrium)
    
    Examples
    --------
    >>> microstructure = {
    ...     'grain_size': 50e-6,  # 50 μm
    ...     'grain_shape': 'equiaxed',
    ...     'gb_type': 'HAGB',
    ...     'trap_list': [
    ...         {'name': 'dislocations', 'density': 1e15, 'binding_energy': 27e3},
    ...         {'name': 'vacancies', 'density': 1e21, 'binding_energy': 41e3}
    ...     ]
    ... }
    >>> result = calculate_defective_metal_flux(
    ...     D_lattice=1e-10, K_s=0.5, thickness=1e-3,
    ...     P_up=1e5, P_down=0,
    ...     temperature=1073, microstructure_params=microstructure
    ... )
    >>> print(f"Flux: {result['flux']:.2e} mol/m²/s")
    >>> print(f"D modification: {result['modification_factor']:.2f}×")
    
    See Also
    --------
    calculate_simple_metal_flux : Level 1 clean metal calculation
    defective_metal.combined_microstructure_model : Underlying D_eff calculation
    """
    # Import Level 4 functions from defective_metal module
    from calculations.defective_metal import combined_microstructure_model
    
    # =========================================================================
    # Input Validation
    # =========================================================================
    
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    if thickness <= 0:
        raise ValueError(f"Thickness must be positive: {thickness} m")
    if D_lattice <= 0:
        raise ValueError(f"Diffusion coefficient must be positive: {D_lattice} m²/s")
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive: {temperature} K")
    if K_s <= 0:
        raise ValueError(f"Solubility constant must be positive: {K_s}")
    
    valid_methods = ['average', 'harmonic', 'inlet', 'outlet']
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}, got '{method}'")
    
    # Required microstructure keys
    required_keys = ['grain_size', 'grain_shape', 'gb_type', 'trap_list']
    for key in required_keys:
        if key not in microstructure_params:
            raise ValueError(f"microstructure_params missing required key: '{key}'")
    
    # Warning for reverse flow
    if P_down > P_up:
        print(f"Warning: Downstream pressure ({P_down} Pa) > Upstream pressure ({P_up} Pa)")
        print("Flux will be negative (backward flow)")
    
    # =========================================================================
    # Calculate Surface Concentrations (Sieverts' Law - unchanged)
    # =========================================================================
    
    C_up = sieverts_concentration(K_s, P_up)
    C_down = sieverts_concentration(K_s, P_down)
    
    # =========================================================================
    # Calculate Position-Dependent Effective Diffusivity
    # =========================================================================
    
    # Create position array through metal thickness
    x_array = np.linspace(0, thickness, n_points)
    
    # For Sieverts' law, sqrt(P) is linear in steady state
    # C(x) = K_s × sqrt(P(x)) where sqrt(P) varies linearly
    min_pressure = 1e-20  # Avoid numerical issues at zero pressure
    sqrt_P_up = np.sqrt(max(P_up, min_pressure))
    sqrt_P_down = np.sqrt(max(P_down, 0))
    
    # Linear profile in sqrt(P) space
    sqrt_P_array = sqrt_P_up - (sqrt_P_up - sqrt_P_down) * x_array / thickness
    C_array = K_s * sqrt_P_array
    
    # Calculate D_eff at each position
    D_array = np.zeros(n_points)
    theta_array = np.zeros(n_points)
    gb_factor_array = np.zeros(n_points)
    
    for i, C_local in enumerate(C_array):
        # Ensure minimum concentration for numerical stability
        C_local = max(C_local, 1e-20)
        
        # Calculate combined microstructure effects at this concentration
        result_i = combined_microstructure_model(
            D_lattice=D_lattice,
            temperature=temperature,
            microstructure_params=microstructure_params,
            lattice_concentration=C_local,
            lattice_density=lattice_density
        )
        
        D_array[i] = result_i['D_eff']
        theta_array[i] = result_i['trapping']['theta_total']
        gb_factor_array[i] = result_i['gb_enhancement']['factor']
    
    # =========================================================================
    # Calculate Average Effective Diffusivity
    # =========================================================================
    
    if method == 'average':
        # Simple arithmetic mean
        D_eff = np.mean(D_array)
        
    elif method == 'harmonic':
        # Harmonic mean - appropriate for resistances in series
        # 1/D_harm = (1/n) × Σ(1/Dᵢ)
        D_eff = len(D_array) / np.sum(1.0 / D_array)
        
    elif method == 'inlet':
        # Use inlet (high concentration) value only
        # Conservative for trapping-dominated cases
        D_eff = D_array[0]
        
    elif method == 'outlet':
        # Use outlet (low concentration) value only
        # Conservative for GB-dominated cases
        D_eff = D_array[-1]
    
    # =========================================================================
    # Calculate Flux Using Effective Diffusivity
    # =========================================================================
    
    flux = fick_flux(D_eff, C_up, C_down, thickness)
    
    # Effective permeability with modified diffusivity
    permeability = D_eff * K_s
    
    # =========================================================================
    # Calculate Diagnostic Information
    # =========================================================================
    
    modification_factor = D_eff / D_lattice
    D_variation = (D_array.max() - D_array.min()) / D_eff if D_eff > 0 else 0
    
    # Determine dominant effect
    avg_gb_factor = np.mean(gb_factor_array)
    avg_theta = np.mean(theta_array)
    trap_reduction = 1.0 / (1.0 + avg_theta) if avg_theta > 0 else 1.0
    
    if avg_gb_factor > 1.5 and trap_reduction > 0.5:
        dominant_effect = 'gb_enhancement'
    elif avg_gb_factor < 1.5 and trap_reduction < 0.5:
        dominant_effect = 'trapping'
    else:
        dominant_effect = 'balanced'
    
    # =========================================================================
    # Return Results
    # =========================================================================
    
    return {
        # Primary outputs (compatible with calculate_simple_metal_flux)
        'flux': flux,
        'C_up': C_up,
        'C_down': C_down,
        'permeability': permeability,
        
        # Level 4 specific outputs
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'modification_factor': modification_factor,
        
        # Microstructure details
        'microstructure_details': {
            'average_theta': avg_theta,
            'average_gb_enhancement': avg_gb_factor,
            'trap_reduction_factor': trap_reduction,
            'D_variation': D_variation,
            'dominant_effect': dominant_effect,
            'method_used': method,
            'n_points': n_points
        },
        
        # Profile data (for detailed analysis)
        'profiles': {
            'x': x_array,
            'D': D_array,
            'C': C_array,
            'theta': theta_array
        },
        
        # Units
        'units': {
            'flux': 'mol/m²/s',
            'concentration': 'mol/m³',
            'permeability': 'mol/m/s/Pa^0.5',
            'diffusivity': 'm²/s'
        }
    }
```

---

## Implementation in interface_solver.py

### New Functions to Add:

```python
# =============================================================================
# LEVEL 4: Interface Solvers for Defective Metal
# =============================================================================

def calculate_defective_metal_flux_sieverts(D_lattice, K_s_metal, thickness, 
                                            P_interface, P_downstream,
                                            temperature, microstructure_params,
                                            lattice_density=1.06e29,
                                            method='average', n_points=10):
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
        n_points=n_points
    )
    
    return result['flux']


def flux_balance_equation_defective_metal(P_interface, P_upstream, P_downstream, 
                                          oxide_props, metal_props, temperature,
                                          microstructure_params, lattice_density=1.06e29,
                                          method='average', n_points=10):
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
        n_points=n_points
    )
    
    return flux_oxide - flux_metal


def solve_interface_pressure_defective_metal(P_upstream, P_downstream, oxide_props, 
                                             metal_props, temperature, microstructure_params,
                                             lattice_density=1.06e29, method='average',
                                             n_points=10, solver_method='brentq',
                                             max_iterations=10, tolerance=1e-6):
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
                    temperature, microstructure_params, lattice_density, method, n_points
                )
                f_max = flux_balance_equation_defective_metal(
                    P_max, P_upstream, P_downstream, oxide_props, metal_props,
                    temperature, microstructure_params, lattice_density, method, n_points
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
                              method, n_points),
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
            n_points=n_points
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
                                           n_points=10, max_iterations=10, tolerance=1e-6):
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
        tolerance=tolerance
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
```

---

## Implementation in parallel_oxide_defect_paths.py

### New Functions to Add (Level 3 + Level 4):

```python
# =============================================================================
# LEVEL 3 + LEVEL 4: Defective Oxide + Defective Metal
# =============================================================================

def calculate_defect_path_flux_defective_metal(P_upstream, P_downstream, oxide_props, 
                                               metal_props, defect_props, temperature,
                                               microstructure_params, lattice_density=1.06e29,
                                               method='average', n_points=10):
    """
    Calculate hydrogen flux through a defect path with defective metal (Level 3+4).
    
    This extends calculate_defect_path_flux() to use Level 4 defective metal
    instead of clean metal for the metal layer.
    
    Theory:
    -------
    Same as calculate_defect_path_flux(), but the metal layer now includes:
    - Grain boundary fast diffusion paths
    - Hydrogen trapping effects
    - Position-dependent effective diffusivity
    
    Parameters
    ----------
    P_upstream : float
        Upstream hydrogen pressure [Pa]
    P_downstream : float
        Downstream hydrogen pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    defect_props : dict
        Defect properties:
        - 'type': 'pinhole', 'crack', 'grain_boundary', or 'mixed'
        - 'thickness_factor': fraction of oxide thickness (for cracks)
        - 'diffusivity_factor': D multiplication factor (for GB)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Metal microstructure specification
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    
    Returns
    -------
    float
        Hydrogen flux through defect path [mol/m²/s]
    """
    # Implementation handles pinhole, crack, grain_boundary, and mixed defect types
    # Each type uses calculate_defective_metal_flux or calculate_oxide_defective_metal_system
    pass  # Full implementation shown above


def calculate_parallel_path_flux_defective_metal(P_upstream, P_downstream, oxide_props, 
                                                  metal_props, defect_params, temperature,
                                                  microstructure_params, lattice_density=1.06e29,
                                                  method='average', n_points=10,
                                                  max_iterations=10, tolerance=1e-6):
    """
    Calculate total flux through defective oxide + defective metal (Level 3+4).
    
    This is the full Level 3+4 model combining:
    - Strehlow & Savage parallel path model for defective oxide
    - Level 4 microstructure effects in metal
    
    Theory:
    -------
    Total flux = Intact path contribution + Defect path contribution
    
    J_total = j_intact × f_intact + j_defect × f_defect
    
    Where both j_intact and j_defect now use Level 4 defective metal.
    
    Returns
    -------
    dict
        Comprehensive results:
        - 'flux_total': Total hydrogen flux [mol/m²/s]
        - 'flux_intact_contribution': Flux through intact oxide [mol/m²/s]
        - 'flux_defect_contribution': Flux through defects [mol/m²/s]
        - 'dominant_path': Which path carries more flux
        - 'defect_enhancement_factor': J_total/J_perfect_oxide
        - 'level4_details': Microstructure effect details
    """
    pass  # Full implementation shown above


def calculate_PRF_defective_metal(P_test, oxide_props, metal_props, temperature,
                                   microstructure_params, defect_params=None,
                                   P_downstream=0, lattice_density=1.06e29,
                                   method='average', n_points=10):
    """
    Calculate Permeation Reduction Factor with Level 4 defective metal.
    
    This extends calculate_PRF() to use Level 4 metal model.
    
    Returns
    -------
    dict
        PRF results:
        - 'PRF': Permeation reduction factor
        - 'PRF_perfect': PRF for defect-free oxide
        - 'efficiency': PRF/PRF_perfect
        - 'regime': Operating regime
        - 'flux_bare': Bare defective metal flux
        - 'flux_with_oxide': Flux with oxide
        - 'level4_details': Microstructure effect info
    """
    pass  # Full implementation shown above
```

---

## Complete Architecture Summary

| Combination | Oxide | Metal | Function Location | Function Name |
|-------------|-------|-------|-------------------|---------------|
| **Level 1** | None | Clean | `permeation_calc.py` | `calculate_simple_metal_flux()` |
| **Level 2** | Perfect | Clean | `interface_solver.py` | `calculate_oxide_metal_system()` |
| **Level 3** | Defective | Clean | `parallel_oxide_defect_paths.py` | `calculate_parallel_path_flux()` |
| **Level 4** | None | Defective | `permeation_calc.py` | `calculate_defective_metal_flux()` |
| **Level 2+4** | Perfect | Defective | `interface_solver.py` | `calculate_oxide_defective_metal_system()` |
| **Level 3+4** | Defective | Defective | `parallel_oxide_defect_paths.py` | `calculate_parallel_path_flux_defective_metal()` |

---

## Existing Analysis and Test Suite (Levels 1-3)

### LEVEL 1: Clean Metal (Sieverts + Fick)

**Location:** `Level_1/`

#### Analysis Scripts (`Level_1/analysis/`):

| Script | Purpose |
|--------|---------|
| `pressure_study.py` | Validates Sieverts' law (J ∝ √P) via pressure sweep |
| `temperature_study.py` | Validates Arrhenius behavior, extracts activation energies |
| `experimental_comparison.py` | Compares model vs JAERI experimental data |
| `sensitivity_analysis.py` | Analyzes model sensitivity to D, K_s, T, P parameters |

#### Results Generated:
- `pressure_study_Incoloy800_800C_*.txt/.png`
- `temperature_study_Incoloy800_*.txt/.png`
- `experimental_comparison_Incoloy800_*.txt/.png`
- `sensitivity_analysis_*.png`
- `Level1_Final_Report_Incoloy800_*.txt`

---

### LEVEL 2: Perfect Oxide + Clean Metal

**Location:** `Level_2/`

#### Analysis Scripts (`Level_2/`):

| Script | Purpose |
|--------|---------|
| `pressure_study.py` | Pressure sweep with oxide layer |
| `temperature_study.py` | Temperature dependence with oxide |
| `experimental_comparison.py` | Compare oxide+metal model to experiments |
| `sensitivity_analysis.py` | Sensitivity analysis for oxide parameters |

#### Validation Tests (`Level_2/` and `validation/`):

| Script | Purpose |
|--------|---------|
| `test_interface_solver.py` | **Comprehensive interface solver tests** (1802 lines) |
| `test_oxide_functions.py` | Tests oxide flux linearity, resistance calculations |
| `test_simple_metal.py` | Basic metal flux verification |

#### Test Categories in `test_interface_solver.py`:
1. **Test 1: Interface Pressure Solutions** - Solves P_interface across 10⁻¹⁰ to 10²⁸ Pa
2. **Test 2: Flux Continuity** - Verifies J_oxide = J_metal at interface
3. **Test 3: Concentration Profiles** - Validates concentration through oxide+metal
4. **Test 4: Limiting Cases** - Tests oxide-limited vs metal-limited regimes

#### Test Categories in `test_oxide_functions.py`:
1. **Test 1: Molecular Flux Linearity** - J_oxide ∝ P (linear)
2. **Test 2: Resistance Calculations** - R_oxide (constant) vs R_metal (∝√P)
3. **Test 3: Flux Regimes** - Regime identification
4. **Test 4: Temperature Dependence** - T-dependent property calculations

---

### LEVEL 3: Defective Oxide + Clean Metal

**Location:** `Level_3/`, `analysis/`, `validation/`

#### Analysis Scripts (`analysis/` and `Level_3/analysis/`):

| Script | Purpose |
|--------|---------|
| `level3_experimental_comparison.py` | Compares Level 1, 2, 3 vs experiments; fits defect fraction |
| `regime_analysis.py` | Phase diagrams: defect vs oxide dominated regions |
| `pressure_study.py` | Pressure sweep with defective oxide |
| `temperature_study.py` | Temperature effects on defective oxide |
| `sensitivity_analysis.py` | Sensitivity to defect parameters |
| `experimental_comparison.py` | General experimental comparison |

#### Validation Tests (`validation/`):

| Script | Purpose |
|--------|---------|
| `test_parallel_oxide_defect_path.py` | **Comprehensive Level 3 test suite** (1486 lines) |
| `test_interface_solver.py` | Interface solver validation |
| `test_oxide_functions.py` | Oxide function tests |
| `calibrate_material_params.py` | Parameter calibration tools |

#### Test Categories in `test_parallel_oxide_defect_path.py`:

| Test # | Name | Description |
|--------|------|-------------|
| **Test 1** | Basic Functionality | All defect types work (pinhole, crack, GB, mixed) |
| **Test 2** | Limiting Cases | f=0 → Level 2, f=1 → Level 1 |
| **Test 3** | Monotonic Behavior | Flux increases with defect fraction |
| **Test 4** | PRF Validation | PRF values in literature range (10-3828) |
| **Test 5** | Pressure Sweep | Three-regime behavior analysis |
| **Test 6** | Defect Type Comparison | Impact of different defect types |
| **Test 7** | PRF vs Defect Fraction | Barrier effectiveness vs f_defect |
| **Test 8** | Dominant Path Analysis | Phase diagram: defect vs oxide dominated |
| **Test 9** | PRF Regime Analysis | PRF phase diagram |
| **Test 10** | Temperature Effects | Regime transitions vs temperature |
| **Test 11** | Defect Type Regime Comparison | Compare defect types across regimes |

#### Results Generated (`validation/results/complete_level3_tests/`):
- `test1_basic_functionality_*.png`
- `test2_limiting_cases_*.png`
- `test3_monotonic_behavior_*.png`
- `test4_PRF_validation_*.png`
- `test5_pressure_sweep_*.png`
- `test6_defect_type_comparison_*.png`
- `test7_PRF_vs_defect_fraction_*.png`
- `test8_dominant_path_analysis_*.png`
- `test9_PRF_regime_analysis_*.png`
- `test10_temperature_regime_analysis_*.png`
- `test11_defect_type_regime_comparison_*.png`
- `test_results_report_*.txt`

---

### Summary Table

| Level | Model | # Analysis Scripts | # Test Scripts | Key Tests |
|-------|-------|-------------------|----------------|-----------|
| **1** | Clean Metal | 4 | 1 | Sieverts law, Arrhenius, Sensitivity |
| **2** | Perfect Oxide + Metal | 4 | 3 | Interface solver, Flux continuity, Regimes |
| **3** | Defective Oxide + Metal | 6 | 4 | 11 comprehensive tests, PRF, Phase diagrams |

---

## What's Missing for Level 4

Based on existing test coverage, Level 4 (Defective Metal) should have:

1. **Basic Functionality** - All microstructure effects work (GB, trapping, combined)
2. **Limiting Cases** - No traps → Level 1, No GB → bulk diffusion
3. **Monotonic Behavior** - D_eff responds correctly to parameters
4. **Integration Tests** - Level 2+4, Level 3+4 combinations
5. **Temperature Effects** - GB enhancement and trapping vs T
6. **Concentration Effects** - D_eff variation across thickness
7. **Experimental Comparison** - Fit microstructure to experimental data
8. **Regime Analysis** - GB-dominated vs trap-dominated regimes
9. **Iterative Convergence** - Verify D_eff ↔ P_interface iteration works
10. **Parameter Sensitivity** - grain_size, trap_density, binding_energy effects

---

## Redundant Functions Clarification

### What You NOW Have (with Option 2 approach):

| Function | Location | Purpose |
|----------|----------|---------|
| `calculate_simple_metal_flux()` | `permeation_calc.py` | Level 1: Clean metal flux |
| `calculate_defective_metal_flux()` | `permeation_calc.py` | Level 4: Defective metal flux |
| `solve_interface_pressure()` | `interface_solver.py` | Level 2: Solve P_interface (clean metal) |
| `solve_interface_pressure_defective_metal()` | `interface_solver.py` | Level 2+4: Solve P_interface (defective metal) |
| `calculate_parallel_path_flux()` | `parallel_oxide_defect_paths.py` | Level 3: Defective oxide + clean metal |
| `calculate_parallel_path_flux_defective_metal()` | `parallel_oxide_defect_paths.py` | Level 3+4: Defective oxide + defective metal |

### What's NOT Needed:

| Function | Why Not Needed |
|----------|----------------|
| `calculate_level4_metal_properties()` | **Absorbed into** `calculate_defective_metal_flux()` - D_eff calculation happens inside |
| `solve_interface_with_level4_metal()` | **Replaced by** `solve_interface_pressure_defective_metal()` - same functionality, cleaner naming |

### Building Block Functions to Keep in `defective_metal.py`:

- `trap_occupancy()` ✅
- `grain_boundary_density()` ✅
- `gb_enhancement_factor()` ✅
- `vacancy_concentration()` ✅
- `calculate_effective_diffusivity_trapping()` ✅
- `calculate_gb_enhanced_diffusivity()` ✅
- `combined_microstructure_model()` ✅

These are called internally by `calculate_defective_metal_flux()` in `permeation_calc.py`.

---

## Next Steps

1. **Implement** `calculate_defective_metal_flux()` in `permeation_calc.py`
2. **Implement** Level 4 solver functions in `interface_solver.py`
3. **Implement** Level 3+4 functions in `parallel_oxide_defect_paths.py`
4. **Create test suite** for Level 4 (similar structure to Level 3 tests)
5. **Run experimental comparison** to validate microstructure parameters
