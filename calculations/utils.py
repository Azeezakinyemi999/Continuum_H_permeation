"""
Utility functions for hydrogen permeation calculations.

Uses REFERENCE-TEMPERATURE Arrhenius format:
    k(T) = k_ref × exp((-E/R) × (1/T - 1/T_ref))

Reference values are read from material/oxide data dictionaries.
"""

import numpy as np

# Gas constant
R = 8.314  # J/mol/K


def arrhenius(value_ref, activation_energy, temperature, T_ref):
    """
    Calculate temperature-dependent property using reference-temperature Arrhenius equation.
    
    Parameters
    ----------
    value_ref : float
        Property value at reference temperature T_ref
    activation_energy : float
        Activation energy [J/mol]
    temperature : float
        Target temperature [K] - the operating condition
    T_ref : float
        Reference temperature [K] - where value_ref was measured
    
    Returns
    -------
    float
        Property value at the target temperature
    
    Notes
    -----
    Reference-temperature Arrhenius: k(T) = k_ref × exp((-E/R) × (1/T - 1/T_ref))
    """
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive: {temperature} K")
    if T_ref <= 0:
        raise ValueError(f"Reference temperature must be positive: {T_ref} K")
    
    return value_ref * np.exp((-activation_energy / R) * (1/temperature - 1/T_ref))


def get_diffusivity(temperature, material_dict):
    """
    Calculate diffusivity at given temperature.
    
    Parameters
    ----------
    temperature : float
        Operating temperature [K]
    material_dict : dict
        Dictionary containing 'D_ref', 'E_D', and 'T_ref'
    
    Returns
    -------
    float
        Diffusivity [m²/s]
    """
    T_ref = material_dict['T_ref']
    D_ref = material_dict['D_ref']
    E_D = material_dict['E_D']
    return arrhenius(D_ref, E_D, temperature, T_ref)


def get_solubility(temperature, material_dict):
    """
    Calculate Sieverts' solubility at given temperature.
    
    Parameters
    ----------
    temperature : float
        Operating temperature [K]
    material_dict : dict
        Dictionary containing 'K_s_ref', 'H_s', and 'T_ref'
    
    Returns
    -------
    float
        Solubility constant [mol/m³/Pa^0.5]
    """
    T_ref = material_dict['T_ref']
    K_s_ref = material_dict['K_s_ref']
    H_s = material_dict['H_s']
    return arrhenius(K_s_ref, H_s, temperature, T_ref)


def get_permeability(temperature, material_dict):
    """
    Calculate permeability at given temperature.
    
    Parameters
    ----------
    temperature : float
        Operating temperature [K]
    material_dict : dict
        Dictionary with D_ref, E_D, K_s_ref, H_s, T_ref
    
    Returns
    -------
    float
        Permeability [mol/m/s/Pa^0.5]
    """
    D = get_diffusivity(temperature, material_dict)
    K_s = get_solubility(temperature, material_dict)
    return D * K_s