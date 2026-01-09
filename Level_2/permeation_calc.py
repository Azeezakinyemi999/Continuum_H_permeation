import numpy as np

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
        'units': {
            'flux': 'mol/m²/s',
            'concentration': 'mol/m³',
            'permeability': 'mol/m/s/Pa^0.5'
        }
    }


