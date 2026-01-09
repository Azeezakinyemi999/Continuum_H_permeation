OXIDE_PROPERTIES = {
    'Cr2O3': {
        # Fixed temperature values (for initial testing)
        'D_ox': 5e-21,  # m²/s at 800°C
        'K_ox': 1e-22,  # mol/m³/Pa at 800°C
        
        # Temperature-dependent parameters
        'D_ox_0': 1e-15,  # m²/s (pre-exponential)
        'E_D_ox': 80000,  # J/mol (activation energy)
        'K_ox_0': 1e-20,  # mol/m³/Pa (pre-exponential)
        'H_sol_ox': 20000,  # J/mol (solution enthalpy)
        
        # Geometric properties
        'thickness': 6e-10,  # m (6 Angstroms)
        'thickness_range': [3e-10, 2e-9],  # m (uncertainty)
        
        # Metadata
        'reference': 'Strehlow & Savage (1974), Zarchy & Axtmann (1979)',
        'temperature_range': [873, 1273],  # K (600-1000°C)
        'uncertainty_factor': 100
    }
}