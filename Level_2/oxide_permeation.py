import numpy as np
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS

def molecular_diffusion_flux(D_ox, K_ox, thickness, P_up, P_down):
    """
    Calculate flux through oxide layer via molecular diffusion.
    
    This is fundamentally different from metal permeation:
    - H2 molecules diffuse intact (don't dissociate)
    - Concentration is linear in pressure (Henry's law)
    - Results in flux linear in pressure difference
    
    Parameters:
    -----------
    D_ox : float
        Molecular diffusion coefficient in oxide (m²/s)
    K_ox : float
        Henry's law constant for H2 in oxide (mol/m³/Pa)
    thickness : float
        Oxide layer thickness (m)
    P_up : float
        Upstream pressure (Pa)
    P_down : float
        Downstream pressure (Pa)
    
    Returns:
    --------
    float
        Flux through oxide (mol/m²/s)
    
    Physics Note:
    ------------
    C = K_ox * P (Henry's law, NOT Sieverts' law)
    J = -D * dC/dx = D * (C_up - C_down) / thickness
    """

    if thickness <= 0:
        raise ValueError("Oxide thickness must be positive")
    if P_up < P_down:
        raise ValueError("Upstream pressure must be >= downstream pressure")


    # Henry's law: concentration linear in pressure
    C_up = K_ox * P_up      # mol/m³
    C_down = K_ox * P_down  # mol/m³
    
    # Fick's first law with linear concentration gradient
    flux = D_ox * (C_up - C_down) / thickness  # mol/m²/s
    
    return flux



# def calculate_oxide_flux_from_material(material_name, temperature_K, P_up, P_down):
#     """
#     Calculate flux using material properties from database.
    
#     Parameters:
#     -----------
#     material_name : str
#         Name of oxide (e.g., 'Cr2O3')
#     temperature_K : float
#         Temperature in Kelvin
#     P_up, P_down : float
#         Upstream and downstream pressures (Pa)
    
#     Returns:
#     --------
#     float
#         Flux through oxide (mol/m²/s)
#     """
#     # Get material properties
#     if material_name not in OXIDE_PROPERTIES:
#         raise ValueError(f"Unknown oxide material: {material_name}")
    
#     oxide_props = OXIDE_PROPERTIES[material_name]
    
#     # Calculate temperature-dependent properties
#     R = 8.314  # J/mol/K
#     D_ox = oxide_props['D_ox_0'] * np.exp(-oxide_props['E_D_ox'] / (R * temperature_K))
#     K_ox = oxide_props['K_ox_0'] * np.exp(-oxide_props['H_sol_ox'] / (R * temperature_K))
    
#     # Use the base function
#     flux = molecular_diffusion_flux(
#         D_ox, K_ox, 
#         oxide_props['thickness'],
#         P_up, P_down
#     )
    
#     return flux


def calculate_oxide_resistance(D_ox, K_ox, thickness):
    """
    Calculate permeation resistance of oxide layer.
    
    Resistance is defined such that:
    Flux = ΔP / Resistance
    
    For molecular diffusion: R = thickness / (D_ox * K_ox)
    
    Parameters:
    -----------
    D_ox : float
        Molecular diffusion coefficient in oxide (m²/s)
    K_ox : float
        Henry's law constant (mol/m³/Pa)
    thickness : float
        Oxide thickness (m)
    
    Returns:
    --------
    float
        Oxide resistance (Pa·s·m²/mol)
    
    Note:
    -----
    This resistance is pressure-independent (linear transport)
    """
    if D_ox <= 0 or K_ox <= 0:
        raise ValueError("D_ox and K_ox must be positive")
    if thickness <= 0:
        raise ValueError("Thickness must be positive")
    resistance = thickness / (D_ox * K_ox)  # Pa·s·m²/mol
    
    return resistance


def calculate_metal_resistance(D_metal, K_s_metal, thickness, P_interface):
    """
    Calculate permeation resistance of metal layer.
    
    Metal resistance depends on interface pressure because of √P solubility!
    This is a KEY DIFFERENCE from oxide resistance.
    
    Parameters:
    -----------
    D_metal : float
        Atomic diffusion coefficient in metal (m²/s)
    K_s_metal : float
        Sieverts' constant (mol/m³/Pa^0.5)
    thickness : float
        Metal thickness (m)
    P_interface : float
        Pressure at oxide/metal interface (Pa)
    
    Returns:
    --------
    float
        Metal resistance at given interface pressure (Pa·s·m²/mol)
    
    Physics Note:
    ------------
    For metal with Sieverts' law:
    J = (D * K_s / thickness) * (√P_up - √P_down)
    
    Linearizing around P_interface:
    R_metal ≈ (thickness * 2 * √P_interface) / (D * K_s )
    
    This is an approximation valid for small ΔP across metal.
    """
    if P_interface <= 0:
        raise ValueError("Interface pressure must be positive")
    if D_metal <= 0 or K_s_metal <= 0:
        raise ValueError("D_metal and K_s_metal must be positive")
    if thickness <= 0:
        raise ValueError("Thickness must be positive")

    # Effective permeability at interface pressure
    # Factor of 2 comes from derivative of √P: i used taylor expansion for this calculation
    effective_permeability = (D_metal * K_s_metal) / (2* np.sqrt(P_interface))
    
    resistance = thickness / effective_permeability  # Pa·s·m²/mol
    
    return resistance



def get_oxide_properties_at_T(oxide_name, temperature_K):
    """
    Calculate temperature-dependent oxide properties.
    
    Parameters:
    -----------
    oxide_name : str
        Name of oxide material (e.g., 'Cr2O3')
    temperature_K : float
        Temperature in Kelvin
    
    Returns:
    --------
    dict
        Contains D_ox, K_ox, thickness at specified temperature
    """
    if oxide_name not in OXIDE_PROPERTIES:
        raise ValueError(f"Unknown oxide material: {oxide_name}")
    
    oxide_data = OXIDE_PROPERTIES[oxide_name]
    R = 8.314  # J/mol/K
    
    # Check temperature range
    T_min, T_max = oxide_data['temperature_range']
    if not (T_min <= temperature_K <= T_max):
        print(f"Warning: Temperature {temperature_K}K outside validated range [{T_min}, {T_max}]K")
    
    # Calculate temperature-dependent properties
    D_ox = oxide_data['D_ox_0'] * np.exp(-oxide_data['E_D_ox'] / (R * temperature_K))
    K_ox = oxide_data['K_ox_0'] * np.exp(-oxide_data['H_sol_ox'] / (R * temperature_K))
    
    return {
        'D_ox': D_ox,
        'K_ox': K_ox,
        'thickness': oxide_data['thickness']
    }


def get_metal_properties_at_T(metal_name, temperature_K):
    """
    Calculate temperature-dependent metal properties.
    Using your Level 1 material data.
    
    Parameters:
    -----------
    metal_name : str
        Name of metal (e.g., 'Incoloy800')
    temperature_K : float
        Temperature in Kelvin
    
    Returns:
    --------
    dict
        Contains D_metal, K_s_metal at specified temperature
    """
    if metal_name not in MATERIALS:
        raise ValueError(f"Unknown metal material: {metal_name}")
    
    metal_data = MATERIALS[metal_name]
    R = 8.314  # J/mol/K
    
    # From your Level 1 implementation
    D_metal = metal_data['D_0'] * np.exp(-metal_data['E_D'] / (R * temperature_K))
    K_s_metal = metal_data['K_s0'] * np.exp(-metal_data['H_s'] / (R * temperature_K))
    
    return {
        'D_metal': D_metal,
        'K_s_metal': K_s_metal
    }

def compare_resistances(oxide_props, metal_props, P_interface):
    """
    Compare oxide and metal resistances to identify limiting mechanism.
    
    Parameters:
    -----------
    oxide_props : dict
        Contains D_ox, K_ox, thickness
    metal_props : dict
        Contains D_metal, K_s_metal, thickness  
    P_interface : float
        Interface pressure (Pa)
    
    Returns:
    --------
    dict
        Contains R_oxide, R_metal, ratio, limiting_mechanism
    """
    R_oxide = calculate_oxide_resistance(
        oxide_props['D_ox'],
        oxide_props['K_ox'],
        oxide_props['thickness']
    )
    
    R_metal = calculate_metal_resistance(
        metal_props['D_metal'],
        metal_props['K_s_metal'],
        metal_props['thickness'],
        P_interface
    )
    
    ratio = R_oxide / R_metal
    
    # Classify regime based on resistance ratio
    if ratio > 100:
        limiting = "oxide_limited"
    elif ratio < 0.1:
        limiting = "metal_limited"
    else:
        limiting = "transition_regime"
    
    return {
        'R_oxide': R_oxide,
        'R_metal': R_metal,
        'ratio': ratio,
        'limiting_mechanism': limiting,
        'P_interface': P_interface
    }

def calculate_transition_pressure(oxide_props, metal_props):
    """
    Estimate the pressure where oxide and metal resistances are equal. This is approximately where the transition occurs.
    
    Parameters:
    -----------
    oxide_props : dict
        Contains D_ox, K_ox, thickness
    metal_props : dict
        Contains D_metal, K_s_metal, thickness
    
    Returns:
    --------
    float
        Approximate transition pressure (Pa)
    """
    # At transition: R_oxide = R_metal
    # X_ox/(D_ox*K_ox) = (X_metal*2*√P_trans)/(D_metal*K_s_metal)
    # Solving for P_trans:
    
    numerator = (oxide_props['thickness'] * metal_props['D_metal'] * metal_props['K_s_metal'])**2
    denominator = (2 * oxide_props['D_ox'] * oxide_props['K_ox'] * metal_props['thickness'])**2
    
    P_transition = numerator / denominator
    
    return P_transition

def pressure_dependence_analysis(P_range, oxide_props, metal_props, T_K):
    """
    Analyze how system behavior changes with pressure.
    Useful for understanding regime transitions.
    
    Parameters:
    -----------
    P_range : array-like
        Range of upstream pressures to analyze (Pa)
    oxide_props : dict or str
        Oxide properties dict or material name
    metal_props : dict or str
        Metal properties dict or material name
    T_K : float
        Temperature in Kelvin
    
    Returns:
    --------
    dict
        Arrays of fluxes, resistances, and regime classifications
    """
    # Handle material names or property dicts
    if isinstance(oxide_props, str):
        oxide_props = get_oxide_properties_at_T(oxide_props, T_K)
    if isinstance(metal_props, str):
        metal_data = get_metal_properties_at_T(metal_props, T_K)
        metal_props = {
            'D_metal': metal_data['D_metal'],
            'K_s_metal': metal_data['K_s_metal'],
            'thickness': 1e-3  # Default 1mm metal thickness
        }
    
    results = {
        'pressures': P_range,
        'oxide_fluxes': [],
        'regimes': [],
        'R_oxide_values': [],
        'R_metal_values': []
    }
    
    for P in P_range:
        # Calculate oxide-only flux for comparison
        flux_oxide = molecular_diffusion_flux(
            oxide_props['D_ox'],
            oxide_props['K_ox'],
            oxide_props['thickness'],
            P, P_down
        )
        results['oxide_fluxes'].append(flux_oxide)
        
        # Compare resistances
        comparison = compare_resistances(oxide_props, metal_props, P)
        results['regimes'].append(comparison['limiting_mechanism'])
        results['R_oxide_values'].append(comparison['R_oxide'])
        results['R_metal_values'].append(comparison['R_metal'])
    
    return results