import numpy as np
from calculations.classify_regime import classify_regime_level14

def sieverts_concentration(K_s, pressure):
    """
    Calculate hydrogen concentration at metal surface using Sieverts' law.
    
    Parameters
    ----------
    K_s : float
        Solubility constant [mol/m³/Pa^0.5]
    pressure : float
        Hydrogen partial pressure [Pa]
    
    Returns
    -------
    float
        Hydrogen concentration [mol/m³]
    
    Notes
    -----
    Sieverts' law: C = K_s * sqrt(P)
    Valid for molecular hydrogen dissociating at surface
    """
    if pressure < 0:
        raise ValueError(f"Pressure cannot be negative: {pressure} Pa")
    
    concentration = K_s * np.sqrt(pressure)
    return concentration


def fick_flux(D, C_up, C_down, thickness):
    """
    Calculate diffusive flux using Fick's first law..
    
    Parameters
    ----------
    D : float
        Diffusion coefficient [m²/s]
    C_up : float
        Upstream concentration [mol/m³]
    C_down : float
        Downstream concentration [mol/m³]
    thickness : float
        Material thickness [m]
    
    Returns
    -------
    float
        Diffusive flux [mol/m²/s]
    
    Notes
    -----
    Fick's first law: J = -D * dC/dx
    For steady-state through flat plate: J = D * (C_up - C_down) / thickness
    Positive flux is from upstream to downstream
    """
    if thickness <= 0:
        raise ValueError(f"Thickness must be positive: {thickness} m")
    if D < 0:
        raise ValueError(f"Diffusion coefficient cannot be negative: {D} m²/s")
    
    flux = D * (C_up - C_down) / thickness
    return flux


def calculate_simple_metal_flux(D, K_s, thickness, P_up, P_down):
    """
    Calculate hydrogen permeation flux through clean metal..
    
    Combines Sieverts' law for surface concentration with Fick's law
    for bulk diffusion.
    
    Parameters
    ----------
    D : float
        Diffusion coefficient [m²/s]
    K_s : float
        Solubility constant [mol/m³/Pa^0.5]
    thickness : float
        Metal thickness [m]
    P_up : float
        Upstream hydrogen pressure [Pa]
    P_down : float
        Downstream hydrogen pressure [Pa]
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'flux': Permeation flux [mol/m²/s]
        - 'C_up': Upstream concentration [mol/m³]
        - 'C_down': Downstream concentration [mol/m³]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
    
    Notes
    -----
    This assumes:
    - Clean metal surfaces (no oxide)
    - Equilibrium at surfaces (fast dissociation/recombination)
    - Steady-state diffusion
    - No trapping or other complications
    """
    # Input validation
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    if P_down > P_up:
        print(f"Warning: Downstream pressure ({P_down} Pa) > Upstream pressure ({P_up} Pa)")
        print("Flux will be negative (backward flow)")
    
    # Calculate surface concentrations using Sieverts' law
    C_up = sieverts_concentration(K_s, P_up)
    C_down = sieverts_concentration(K_s, P_down)
    
    # Calculate flux using Fick's law
    flux = fick_flux(D, C_up, C_down, thickness)
    
    # Calculate effective permeability
    # P = D * K_s for metal membranes
    permeability = D * K_s
    
    # Return comprehensive results
    return {
        'flux': flux,
        'C_up': C_up,
        'C_down': C_down,
        'permeability': permeability,
        'Diffusivity': D,
        'solubility': K_s,
        'units': {
            'flux': 'mol/m²/s',
            'concentration': 'mol/m³',
            'permeability': 'mol/m/s/Pa^0.5',
            'Diffusivity': 'm²/s',
            'solubility': 'mol/m³/Pa^0.5'
        }
    }


# =============================================================================
# LEVEL 4: Defective Metal Flux Calculations
# =============================================================================

def calculate_defective_metal_flux(D_lattice, K_s, thickness, P_up, P_down,
                                    temperature, microstructure_params,
                                    lattice_density=1.06e29,
                                    method='average', n_points=10, mode='both'):
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
    # Validate method
    valid_methods = ['average', 'harmonic', 'inlet', 'outlet']
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}, got '{method}'")
    # Validate mode
    valid_modes = ['both', 'gb_only', 'trapping_only', 'none']
    if mode not in valid_modes:
        raise ValueError(f"mode must be one of {valid_modes}, got '{mode}'")
    
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
            lattice_density=lattice_density,
            mode=mode
        )
        
        D_array[i] = result_i['D_eff']
        # Handle None cases for trapping/gb_enhancement
        # theta_array[i] = result_i['trapping']['theta_total'] if result_i['trapping'] else 0.0
        # gb_factor_array[i] = result_i['gb_enhancement']['factor'] if result_i['gb_enhancement'] else 1.0
            # Handle different return structures based on mode
        # The keys depend on what combined_microstructure_model returns
        if 'trapping' in result_i and result_i['trapping'] is not None:
            theta_array[i] = result_i['trapping'].get('theta_total', 0.0)
        elif 'theta_total' in result_i:
            theta_array[i] = result_i['theta_total']
        else:
            theta_array[i] = 0.0
            
        if 'gb_enhancement' in result_i and result_i['gb_enhancement'] is not None:
            gb_factor_array[i] = result_i['gb_enhancement'].get('factor', 1.0)
        elif 'gb_enhancement_factor' in result_i:
            gb_factor_array[i] = result_i['gb_enhancement_factor']
        else:
            gb_factor_array[i] = 1.0
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

    # Regime classification for Level 1,4 (bare metal / defective metal)
    regime_classification = classify_regime_level14(modification_factor)
    regime = regime_classification.get('regime_hierarchy')
    regime_base = regime_classification.get('base_regime')
    regime_detail = regime_classification.get('regime_detail')
    
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
        # Regime classification (Level 1,4)
        'regime_classification': regime_classification,
        'regime': regime,
        'regime_base': regime_base,
        'regime_detail': regime_detail,
        
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


def calculate_defective_metal_flux_sieverts(D_lattice, K_s, thickness, P_interface, P_downstream,
                                            temperature, microstructure_params,
                                            lattice_density=1.06e29,
                                            method='average', n_points=10, mode='both'):
    """
    Calculate flux through defective metal using Sieverts' law boundary conditions.
    
    This is the Level 4 equivalent of calculate_metal_flux_sieverts() from
    interface_solver.py. Designed for use in iterative interface pressure solving.
    
    Parameters
    ----------
    D_lattice : float
        Intrinsic lattice diffusion coefficient [m²/s]
    K_s : float
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
        Microstructure specification (see calculate_defective_metal_flux)
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    
    Returns
    -------
    float
        Flux through defective metal [mol/m²/s]
    
    Notes
    -----
    This function returns just the flux value (float) for compatibility with
    the brentq solver in interface_solver.py. For full details, use
    calculate_defective_metal_flux().
    """
    result = calculate_defective_metal_flux(
        D_lattice=D_lattice,
        K_s=K_s,
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