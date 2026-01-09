"""
Experimental permeation data extracted from literature
Proper citation is crucial for academic integrity
"""

import numpy as np

# Data for Incoloy 800 from JAERI-Tech 2002-090
# Data extracted from Fig. 2.2 using WebPlotDigitizer
# Original units: Permeability in cm³(NTP)·cm⁻¹·s⁻¹·Pa⁻⁰·⁵
# Citation: "Data for Incoloy 800 [JAERI-Tech 2002-090] were extracted from Fig. 2.2 using WebPlotDigitizer."

INCOLOY_800_JAERI = {
    'material': 'Incoloy 800',
    'source': 'JAERI-Tech 2002-090',
    'figure': 'Fig. 2.2',
    'extraction_method': 'WebPlotDigitizer',
    'original_units': {
        'temperature': '1000/T (K^-1)',
        'permeability': 'ln(P) where P in cm³(NTP)·cm⁻¹·s⁻¹·Pa⁻⁰·⁵'
    },
    'si_units': {
        'temperature': '1000/T (K^-1)',
        'permeability': 'ln(P) where P in mol·m⁻¹·s⁻¹·Pa⁻⁰·⁵'
    },
    'data': {
        '1000/T': np.array([
            0.81875, 0.826388889, 0.826388889, 0.839583333, 0.858333333,
            0.875694444, 0.901388889, 0.923611111, 0.949305556, 0.976388889,
            1.0, 1.025694444, 1.049305556, 1.070833333, 1.088194444,
            1.108333333, 1.129166667, 1.14375
        ]),
        'ln_permeability_NTP': np.array([  # Original extracted values in ln(P)
            -6.826361584, -6.860371831, -6.860371831, -6.903277893, -6.959525597,
            -7.024811789, -7.130454546, -7.210858862, -7.305800117, -7.406610942,
            -7.513274282, -7.593483584, -7.698099596, -7.789663426, -7.857900144,
            -7.936631741, -8.030998435, -8.081294884
        ])
    }
}

def convert_NTP_to_SI(log_P_NTP):
    """
    Convert permeability from NTP volumetric units to molar SI units.
    
    From: cm³(NTP)·cm⁻¹·s⁻¹·Pa⁻⁰·⁵
    To: mol·m⁻¹·s⁻¹·Pa⁻⁰·⁵
    
    Parameters
    ----------
    log_P_NTP : float or array
        Natural log of permeability in cm³(NTP)·cm⁻¹·s⁻¹·Pa⁻⁰·⁵
    
    Returns
    -------
    float or array
        Natural log of permeability in mol·m⁻¹·s⁻¹·Pa⁻⁰·⁵
    
    Notes
    -----
    Conversion factors:
    - 1 cm³(NTP) = 1/22414 mol (at NTP: 0°C, 1 atm)
    - 1 cm⁻¹ = 100 m⁻¹
    - Overall factor: (1/22414) × 100 = 0.00446
    - ln(P_SI) = ln(P_NTP × 0.00446) = ln(P_NTP) + ln(0.00446)
    """
    conversion_factor = (1.0 / 22414) * 100  # = 0.00446
    log_conversion = np.log(conversion_factor)  # = -5.41
    
    return log_P_NTP + log_conversion


# Convert to SI units and add to dictionary
INCOLOY_800_JAERI['data']['ln_permeability_SI'] = convert_NTP_to_SI(
    INCOLOY_800_JAERI['data']['ln_permeability_NTP']
)

def get_experimental_data(material='Incoloy 800', source='JAERI'):
    """
    Retrieve experimental data for specified material and source.
    
    Parameters
    ----------
    material : str
        Material name
    source : str
        Data source identifier
    
    Returns
    -------
    dict
        Experimental data dictionary with SI units
    """
    if material == 'Incoloy 800' and 'JAERI' in source:
        return INCOLOY_800_JAERI
    else:
        raise ValueError(f"No experimental data available for {material} from {source}")

def convert_to_SI(data_dict):
    """
    Convert experimental data to SI units for comparison.
    
    Parameters
    ----------
    data_dict : dict
        Experimental data dictionary
    
    Returns
    -------
    dict
        Data in SI units with temperature in K and permeability in mol/m/s/Pa^0.5
    """
    # Convert 1000/T to T
    temperatures_K = 1000.0 / data_dict['data']['1000/T']
    
    # Convert ln(P) to P (already in SI units)
    permeabilities = np.exp(data_dict['data']['ln_permeability_SI'])
    
    return {
        'temperatures_K': temperatures_K,
        'permeabilities': permeabilities,
        'material': data_dict['material'],
        'source': data_dict['source']
    }

# Test the conversion
if __name__ == "__main__":
    data = get_experimental_data('Incoloy 800', 'JAERI')
    si_data = convert_to_SI(data)
    
    print("Unit Conversion Check:")
    print("-" * 50)
    print(f"Original (NTP): ln(P) = {data['data']['ln_permeability_NTP'][5]:.3f}")
    print(f"Converted (SI): ln(P) = {data['data']['ln_permeability_SI'][5]:.3f}")
    print(f"Final P (SI): {si_data['permeabilities'][5]:.3e} mol/m/s/Pa^0.5")
    print(f"Temperature: {si_data['temperatures_K'][5]:.1f} K ({si_data['temperatures_K'][5]-273.15:.1f}°C)")

# Add more experimental datasets as needed
# For example:
# INCOLOY_800_FORCEY = { ... }
# SS316_SERRA = { ... }