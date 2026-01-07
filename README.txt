Start level 1:
Phase1: Setup and First calculation
- step 1.1: Create project structure 
        |   permeation_project/
        ├── calculations/
        │   ├── permeation_calc.py
        │   └── utils.py
        ├── data/
        │   ├── materials_data.py
        │   └── experimental_data.py
        ├── validation/
        │   └── test_simple_metal.py
        └── run_analysis.py

- step 1.2 Define Your Material Data
       - Start with ONE material (Incoloy 800):
       - Find D₀ and E_D from literature
      - Find K_s₀ and ΔH_s from literature
        Document your source
        Include temperature range of validity

- Step 1.3: Write Core Calculation Function
        File: calculations/permeation_calc.py
        Functions to write:
        sieverts_concentration(K_s, pressure)
        fick_flux(D, C_up, C_down, thickness)
        calculate_simple_metal_flux(D, K_s, thickness, P_up, P_down)
        arrhenius(pre_exp, activation_energy, temperature)

- Step 1.4: First Validation
    Calculate flux at T=800°C, P=1 Pa, thickness=1mm
    Check if your result has reasonable magnitude (typically 10⁻¹⁰ to 10⁻⁶ mol/m²/s)
    Verify units throughout calculation


Phase 2: Pressure Dependence Study (Days 4-5)
Step 2.1: Pressure Sweep

Step 2.2: Validation Plot
    Plot log(flux) vs log(pressure)
    Calculate slope using linear regression
    Success criterion: Slope should be 0.5 ± 0.01

Step 2.3: Document Findings
    Create a results dictionary:


Phase 3: Temperature Dependence (Days 6-7)
Step 3.1: Temperature Functions
        Add functions:

        get_diffusivity(T, material_dict)
        get_solubility(T, material_dict)
        get_permeability(T, material_dict)

Step 3.2: Temperature Sweep

Step 3.3: Arrhenius Validation
        Plot ln(permeability) vs 1000/T
        Should be linear
        Extract activation energy from slope
        Compare to literature value

Phase 4: Experimental Comparison (Days 8-10)
Step 4.1: Get Experimental Data
    Find ONE experimental curve from the papers (e.g., Incoloy 800 from JAERI-Tech 2002-090)
    Extract data points (use WebPlotDigitizer if needed)
    Save in experimental_data.py

Phase 5: Analysis and Documentation (Days 11-12)
    Step 5.1: Sensitivity Analysis
    Test sensitivity to parameters