from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.utils import get_diffusivity, get_solubility
from data.material_data import MATERIALS

# Your analysis code
# from calculations.permeation_calc import calculate_simple_metal_flux
# from calculations.utils import get_diffusivity, get_solubility
# from data.material_data import MATERIALS

# Test at 800°C
T_celsius = 800
T_kelvin = T_celsius + 273.15

# Get material properties at temperature
incoloy = MATERIALS['Incoloy800']
D = get_diffusivity(T_kelvin, incoloy)
K_s = get_solubility(T_kelvin, incoloy)

# Calculate flux
thickness = 0.001  # 1 mm
P_up = 1  # Pa
P_down = 0  # Pa

result = calculate_simple_metal_flux(D, K_s, thickness, P_up, P_down)
print(f"\nFor P_up = {P_up} Pa:")
print(f"C_up = {result['C_up']:.2e} mol/m³")
print(f"C_down = {result['C_down']:.2e} mol/m³")
print(f"Flux = {result['flux']:.2e} mol/m²/s")
print(f"permeability at {T_celsius}°C: {result['permeability']:.2e} mol/m/s/Pa^0.5")
