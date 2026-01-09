# Start with ONE material (Incoloy 800):
# •	Find D₀ and E_D from literature
# •	Find K_s₀ and ΔH_s from literature
# •	Document your source
# •	Include temperature range of validity


# Start with simple dictionaries
## Made up data by Azeez: I believe i use be able to get these from DFT or MACE/LAMMPS
MATERIALS = {
    'Incoloy800': {
        'D_0': 3.4e-7,  # m2/s
        'E_D': 54000,  # J/mol
        'K_s0': 2.1e-3,  # mol/m3/Pa^0.5
        'H_s': -15000,  # J/mol
        'reference': 'Forcey1988',
        'temp_range':[600, 1000] # celsius
    }

    #     'Incoloy800': {
    #     'D_0': 3.4e-8,  # m2/s
    #     'E_D': 56400,  # J/mol
    #     'K_s0': 2.1e-8,  # mol/m3/Pa^0.5
    #     'H_s': -56600,  # J/mol
    #     'reference': 'Forcey1988',
    #     'temp_range':[750, 950] # celsius
    # }
}