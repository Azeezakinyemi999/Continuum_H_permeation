import numpy as np

def arrhenius(pre_exp, activation_energy, temperature):
    """
    Calculate temperature-dependent property using Arrhenius equation.
    
    Parameters
    ----------
    pre_exp : float
        Pre-exponential factor (same units as property)
    activation_energy : float
        Activation energy [J/mol]
    temperature : float
        Absolute temperature [K]
    
    Returns
    -------
    float
        Temperature-dependent property value
    
    Notes
    -----
    Arrhenius equation: property = pre_exp * exp(-E_a / R*T)
    R = 8.314 J/mol/K (universal gas constant)
    
    Common uses:
    - Diffusivity: D = D_0 * exp(-E_D / RT)
    - Solubility: K_s = K_s0 * exp(-H_s / RT)
    - Permeability: P = P_0 * exp(-E_p / RT)
    """
    R = 8.314  # J/mol/K
    
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive: {temperature} K")
    
    # Handle both endothermic (positive E) and exothermic (negative E) processes
    property_value = pre_exp * np.exp(-activation_energy / (R * temperature))
    
    return property_value

def get_diffusivity(temperature, material_dict):
    """
    Calculate diffusivity at given temperature.
    
    Parameters
    ----------
    temperature : float
        Temperature [K]
    material_dict : dict
        Dictionary containing 'D_0' and 'E_D'
    
    Returns
    -------
    float
        Diffusivity [m²/s]
    """
    return arrhenius(material_dict['D_0'], material_dict['E_D'], temperature)


def get_solubility(temperature, material_dict):
    """
    Calculate Sieverts' solubility at given temperature.
    
    Parameters
    ----------
    temperature : float
        Temperature [K]
    material_dict : dict
        Dictionary containing 'K_s0' and 'H_s'
    
    Returns
    -------
    float
        Solubility constant [mol/m³/Pa^0.5]
    """
    return arrhenius(material_dict['K_s0'], material_dict['H_s'], temperature)


def get_permeability(temperature, material_dict):
    """
    Calculate permeability at given temperature.
    
    Parameters
    ----------
    temperature : float
        Temperature [K]
    material_dict : dict
        Dictionary with D_0, E_D, K_s0, H_s
    
    Returns
    -------
    float
        Permeability [mol/m/s/Pa^0.5]
    """
    D = get_diffusivity(temperature, material_dict)
    K_s = get_solubility(temperature, material_dict)
    return D * K_s