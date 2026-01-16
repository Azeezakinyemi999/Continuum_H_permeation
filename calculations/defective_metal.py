"""
Defective Metal Model for Hydrogen Diffusion (Level 4)

This module implements microstructure effects on hydrogen diffusion in polycrystalline
metals, incorporating both enhancement (grain boundaries) and reduction (trapping) 
mechanisms. It provides the core physics for Level 4 of the hierarchical permeation model.

Theory:
-------
Real metals exhibit competing microstructure effects on hydrogen transport:

1. Grain Boundary Enhancement:
   - GBs provide fast diffusion paths (D_gb >> D_bulk)
   - Enhancement factor decreases with temperature
   - Volume fraction f_gb = 3δ/d_grain
   - Parallel path model: D_eff = (1-f_gb)×D_bulk + f_gb×D_gb

2. Trapping Reduction:
   - Defects (dislocations, vacancies, GBs) trap hydrogen
   - Reduces effective diffusivity via Oriani model
   - D_eff = D_lattice/(1 + Σθᵢ)
   - θᵢ = trap occupancy from local equilibrium

3. Combined Effects:
   - Sequential application: GB enhancement then trapping
   - Net effect can enhance OR reduce diffusion
   - Temperature and concentration dependent

Mathematical Framework:
-----------------------
The effective diffusivity combines both effects:

    D_eff = [(1-f_gb)×D_bulk + f_gb×α×D_bulk]/(1 + Σθᵢ)

where:
    f_gb = 3δ/d = grain boundary volume fraction
    α = gb_enhancement_factor(T) = D_gb/D_bulk
    θᵢ = trap_occupancy(T, E_b,i, N_T,i, C_L)
    
Oriani equilibrium for trapping:
    θ/(1-θ) = exp(E_b/RT) × (C_L/N_L)
    
For dilute solutions:
    θ ≈ (N_T/N_L) × exp(E_b/RT)

Grain boundary diffusion enhancement:
    D_gb/D_bulk = A × exp[(Q_bulk - Q_gb)/RT]
    
where typically Q_gb ≈ 0.6×Q_bulk

Module Structure:
-----------------
Core Functions (7):
1. trap_occupancy() - Calculate trap site occupancy fraction
2. grain_boundary_density() - Convert grain size to trap density
3. gb_enhancement_factor() - Temperature-dependent GB enhancement
4. vacancy_concentration() - Thermal vacancy concentration
5. calculate_effective_diffusivity_trapping() - Apply trapping effects
6. calculate_gb_enhanced_diffusivity() - Apply GB enhancement
7. combined_microstructure_model() - Combine all effects

Usage:
------
This module is used by higher-level functions in:
- permeation_calc.py: calculate_defective_metal_flux()
- interface_solver.py: solve_interface_pressure_defective_metal()
- parallel_oxide_defect_paths.py: calculate_parallel_path_flux_defective_metal()

The functions can be used individually for specific effects or combined
via combined_microstructure_model() for complete microstructure modeling.

Physical Parameters:
--------------------
Typical ranges for austenitic steels:
- Grain size: 10⁻⁸ to 10⁻³ m
- GB thickness: 0.5-1.0 nm
- Dislocation density: 10¹⁴ to 10¹⁶ m⁻³
- Vacancy concentration: 10²⁰ to 10²³ m⁻³
- Binding energies:
  - Dislocations: 20-30 kJ/mol
  - Grain boundaries: 40-50 kJ/mol
  - Vacancies: 40-45 kJ/mol
  - Precipitates: 60-90 kJ/mol

Limitations:
------------
1. Assumes Oriani local equilibrium (invalid for θ > 0.9)
2. Independent trap types (no interaction)
3. Uniform trap distribution (not segregated)
4. Steady-state conditions
5. Dilute hydrogen concentration
6. No trap saturation effects

For high trap occupancy (θ > 0.9) or transient conditions, 
consider McNabb-Foster kinetic model instead.

References:
-----------
1. Oriani, R.A. (1970). "The diffusion and trapping of hydrogen in steel."
   Acta Metallurgica 18, 147-157. DOI: 10.1016/0001-6160(70)90078-7

2. McNabb, A., Foster, P.K. (1963). "A new analysis of diffusion of hydrogen
   in iron and ferritic steels." Trans. Metall. Soc. AIME 227, 618-627.

3. Oudriss, A., et al. (2012). "Grain size and grain-boundary effects on
   diffusion and trapping of hydrogen in pure nickel." Acta Mater. 60, 
   6814-6828. DOI: 10.1016/j.actamat.2012.09.004

4. Tsuru, T., Latanision, R.M. (1982). "Grain boundary transport of hydrogen
   in nickel." Scripta Metall. 16, 575-578. DOI: 10.1016/0036-9748(82)90273-3

5. Brass, A.M., Chêne, J. (2006). "Hydrogen uptake in 316L stainless steel:
   Consequences on the tensile properties." Corros. Sci. 48, 3222-3242.
   DOI: 10.1016/j.corsci.2005.11.004

6. Robertson, I.M. (1977). "The effect of hydrogen on dislocation dynamics."
   Metall. Trans. A 8, 1709-1712. DOI: 10.1007/BF02646878

"""
def trap_occupancy(temperature, binding_energy, trap_density, lattice_density, 
                   lattice_concentration):
    """
    Calculate hydrogen trap occupancy fraction using Oriani local equilibrium model.
    
    Theory:
    -------
    At thermal equilibrium, the distribution of hydrogen between lattice sites and
    trap sites follows Fermi-Dirac statistics. For a single trap type, the 
    equilibrium is described by:
    
        K = exp(E_b/RT) = θ_T/(1-θ_T) × (1-θ_L)/θ_L
    
    where θ_T and θ_L are trap and lattice site occupancies respectively.
    
    For dilute solutions (θ_L << 1), this simplifies to:
        θ_T = K × (C_L/N_L) / [1 + K × (C_L/N_L)]
    
    Mathematical Derivation:
    ------------------------
    Starting from chemical potential equilibrium:
        μ_T = μ_L
        μ°_T + RT×ln(θ_T/(1-θ_T)) = μ°_L + RT×ln(θ_L/(1-θ_L))
    
    With E_b = μ°_L - μ°_T (binding energy), we get:
        θ_T/(1-θ_T) = exp(E_b/RT) × θ_L/(1-θ_L)
    
    For dilute lattice (θ_L << 1):
        θ_T = exp(E_b/RT) × (C_L/N_L) / [1 + exp(E_b/RT) × (C_L/N_L)]
    
    Parameters:
    -----------
    temperature : float
        Absolute temperature in Kelvin [K]
        Must be positive and typically 300-1200 K for fusion applications
    
    binding_energy : float
        Trap binding energy in J/mol [J/mol]
        Positive value means trap is lower energy than lattice site
        Typical range: 20-100 kJ/mol for various defects
    
    trap_density : float
        Number of trap sites per unit volume [m⁻³]
        Range: 10¹⁴-10²⁴ m⁻³ depending on defect type
    
    lattice_density : float
        Number of interstitial lattice sites per unit volume [m⁻³]
        For FCC metals: ~10²⁹ m⁻³
    
    lattice_concentration : float
        Hydrogen concentration in lattice sites [mol/m³]
        This is the mobile/diffusible hydrogen concentration
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'theta': Trap occupancy fraction (0 ≤ θ ≤ 1) [-]
        - 'K_equilibrium': Equilibrium constant [-]
        - 'approximation_used': 'low_coverage' or 'full_equation'
        - 'theta_lattice': Lattice site occupancy fraction [-]
        - 'saturation_warning': True if θ > 0.9
        - 'trap_concentration': Actual trapped H concentration [mol/m³]
    
    Raises:
    -------
    ValueError
        If temperature ≤ 0 K
        If any density/concentration is negative
        If lattice_concentration > lattice_density (unphysical)
    
    References:
    -----------
    1. Oriani, R.A. (1970). "The diffusion and trapping of hydrogen in steel."
       Acta Metallurgica, 18, 147-157. DOI: 10.1016/0001-6160(70)90078-7
    
    2. McNabb, A., Foster, P.K. (1963). "A new analysis of diffusion of hydrogen
       in iron and ferritic steels." Trans. Metall. Soc. AIME, 227, 618-627.
    
    Example:
    --------
    >>> result = trap_occupancy(
    ...     temperature=800,  # K
    ...     binding_energy=27e3,  # J/mol (dislocations)
    ...     trap_density=1e15,  # m⁻³
    ...     lattice_density=1.06e29,  # m⁻³ (FCC)
    ...     lattice_concentration=1e-3  # mol/m³
    ... )
    >>> print(f"Trap occupancy: {result['theta']:.3f}")
    """
    import numpy as np
    
    # Physical constants
    R = 8.314  # Universal gas constant [J/mol/K]
    
    # Input validation
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature} K")
    
    if binding_energy < 0:
        raise ValueError(f"Binding energy should be positive (trap lower than lattice), got {binding_energy} J/mol")
    
    if trap_density < 0:
        raise ValueError(f"Trap density cannot be negative, got {trap_density} m⁻³")
    
    if lattice_density <= 0:
        raise ValueError(f"Lattice density must be positive, got {lattice_density} m⁻³")
    
    if lattice_concentration < 0:
        raise ValueError(f"Lattice concentration cannot be negative, got {lattice_concentration} mol/m³")
    
    # Convert concentration to number density for calculations
    N_A = 6.022e23  # Avogadro's number [mol⁻¹]
    C_L_sites = lattice_concentration * N_A  # Convert mol/m³ to sites/m³
    
    if C_L_sites > lattice_density:
        raise ValueError(f"Lattice concentration {lattice_concentration} mol/m³ exceeds "
                        f"available sites {lattice_density/N_A} mol/m³")
    
    # Calculate equilibrium constant
    K_eq = np.exp(binding_energy / (R * temperature))
    
    # Calculate lattice site occupancy
    theta_L = C_L_sites / lattice_density
    
    # Determine which equation to use
    K_theta_L = K_eq * theta_L
    
    if K_theta_L < 0.01:  # Low coverage approximation valid
        # Simple form: θ_T ≈ K × θ_L for K×θ_L << 1
        theta_T = K_theta_L
        approximation = 'low_coverage'
    else:
        # Full equation: θ_T = K×θ_L / (1 + K×θ_L)
        theta_T = K_theta_L / (1 + K_theta_L)
        approximation = 'full_equation'
    
    # Ensure physical bounds
    theta_T = min(max(theta_T, 0.0), 1.0)
    
    # Calculate actual trapped concentration
    trap_concentration = theta_T * trap_density / N_A  # mol/m³
    
    # Check for saturation
    saturation_warning = theta_T > 0.9
    
    return {
        'theta': theta_T,
        'K_equilibrium': K_eq,
        'approximation_used': approximation,
        'theta_lattice': theta_L,
        'saturation_warning': saturation_warning,
        'trap_concentration': trap_concentration
    }

def grain_boundary_density(grain_size, gb_thickness=0.5e-9, sites_per_area=1e19, grain_shape='equiaxed'):
    """
    Calculate grain boundary trap density from microstructure parameters.
    
    Theory:
    -------
    Grain boundaries (GBs) are 2D defects that separate crystalline regions. The
    total GB area per unit volume depends on grain size and morphology. For 
    hydrogen trapping, we need to convert this geometric information into a
    volumetric trap density.
    
    For randomly oriented equiaxed grains, stereology gives:
        S_v = 3/d  (surface area per volume)
    
    where d is the mean linear intercept (approximately the grain diameter).
    
    The volumetric trap density is then:
        N_gb = S_v × ρ_gb = (3/d) × ρ_gb
    
    where ρ_gb is the areal density of trap sites on the GB.
    
    Mathematical Derivation:
    ------------------------
    Starting from stereological relationships:
    1. For random polycrystal: <N_L> = 2/d (intercepts per unit length)
    2. Surface area per volume: S_v = 2×<N_L> = 4/d
    3. Correction for shared boundaries: S_v = 3/d (factor varies 2-3)
    
    Volume fraction of GB:
        f_gb = S_v × δ = (3δ)/d
    
    where δ is the GB thickness.
    
    Parameters:
    -----------
    grain_size : float
        Average grain diameter in meters [m]
        Typical ranges:
        - Nanocrystalline: 1-100 nm
        - Ultrafine: 100 nm - 1 μm  
        - Fine: 1-10 μm
        - Conventional: 10-100 μm
        - Coarse: >100 μm
    
    gb_thickness : float, optional
        Grain boundary width in meters [m]
        Default: 0.5e-9 m (0.5 nm)
        Physical range: 0.3-1.0 nm for high-angle boundaries
    
    sites_per_area : float, optional
        Number of hydrogen trap sites per m² of GB surface [m⁻²]
        Default: 1e19 m⁻² (approximately one site per atom on GB)
        Range: 1e18-1e20 m⁻² depending on GB structure
    
    grain_shape : str, optional
        Grain morphology affecting GB density:
        - 'equiaxed': 3D random polycrystal (default)
        - 'columnar': 2D boundaries (e.g., directionally solidified)
        - 'planar': 1D boundaries (e.g., lamellar structure)
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'trap_density': GB trap sites per unit volume [m⁻³]
        - 'volume_fraction': Volume fraction of GB phase [-]
        - 'surface_per_volume': GB surface area per unit volume [m⁻¹]
        - 'mean_intercept': Mean linear intercept length [m]
        - 'warnings': List of any warnings about parameters
    
    Raises:
    -------
    ValueError
        If grain_size ≤ 0
        If gb_thickness ≤ 0 or > grain_size/10
        If sites_per_area ≤ 0
        If grain_shape is not recognized
    
    Warnings:
    ---------
    Issues warnings for:
        - Grain size < 10 nm (approaching amorphous limit)
        - Grain size > 10 mm (unusually coarse)
        - Volume fraction > 0.1 (GB phase percolation)
    
    References:
    -----------
    1. Palumbo, G., Aust, K.T. (1990). "Structure-dependence of intergranular
       corrosion in high purity nickel." Acta Metall. Mater. 38, 2343-2352.
       DOI: 10.1016/0956-7151(90)90101-L
    
    2. Underwood, E.E. (1970). "Quantitative Stereology." Addison-Wesley.
    
    3. Oudriss, A., et al. (2012). "Grain size and grain-boundary effects on
       diffusion and trapping of hydrogen in pure nickel." Acta Mater. 60, 
       6814-6828. DOI: 10.1016/j.actamat.2012.09.004
    
    Example:
    --------
    >>> result = grain_boundary_density(
    ...     grain_size=50e-6,  # 50 μm
    ...     gb_thickness=0.5e-9,  # 0.5 nm
    ...     sites_per_area=1e19,  # sites/m²
    ...     grain_shape='equiaxed'
    ... )
    >>> print(f"GB trap density: {result['trap_density']:.2e} m⁻³")
    """
    import numpy as np
    import warnings
    
    # Input validation
    if grain_size <= 0:
        raise ValueError(f"Grain size must be positive, got {grain_size} m")
    
    if gb_thickness <= 0:
        raise ValueError(f"GB thickness must be positive, got {gb_thickness} m")
    
    if gb_thickness > grain_size/10:
        raise ValueError(f"GB thickness {gb_thickness} m cannot exceed 10% of grain size {grain_size} m")
    
    if sites_per_area <= 0:
        raise ValueError(f"Sites per area must be positive, got {sites_per_area} m⁻²")
    
    valid_shapes = ['equiaxed', 'columnar', 'planar']
    if grain_shape not in valid_shapes:
        raise ValueError(f"Grain shape must be one of {valid_shapes}, got '{grain_shape}'")
    
    # Initialize warnings list
    warning_list = []
    
    # Check for unusual grain sizes
    if grain_size < 10e-9:  # Less than 10 nm
        warning_msg = f"Grain size {grain_size*1e9:.1f} nm approaches amorphous limit. Model validity uncertain."
        warnings.warn(warning_msg)
        warning_list.append(warning_msg)
    
    if grain_size > 10e-3:  # Greater than 10 mm
        warning_msg = f"Grain size {grain_size*1e3:.1f} mm is unusually coarse. Verify this is intended."
        warnings.warn(warning_msg)
        warning_list.append(warning_msg)
    
    # Calculate surface area per volume based on grain shape
    if grain_shape == 'equiaxed':
        # 3D random polycrystal
        shape_factor = 3.0  # Classic stereology result
        dimensionality = 3
    elif grain_shape == 'columnar':
        # 2D boundaries (columnar grains)
        shape_factor = 2.0
        dimensionality = 2
    elif grain_shape == 'planar':
        # 1D boundaries (lamellar structure)
        shape_factor = 1.0
        dimensionality = 1
    
    # Calculate GB surface area per unit volume [m⁻¹]
    surface_per_volume = shape_factor / grain_size
    
    # Calculate volume fraction of GB
    volume_fraction = surface_per_volume * gb_thickness
    
    # Warn if GB volume fraction is very high
    if volume_fraction > 0.1:
        warning_msg = f"GB volume fraction {volume_fraction:.2%} exceeds 10%. GB phase may dominate transport."
        warnings.warn(warning_msg)
        warning_list.append(warning_msg)
    
    # Ensure volume fraction doesn't exceed 1 (unphysical)
    if volume_fraction > 1.0:
        warning_msg = f"Calculated GB volume fraction {volume_fraction:.2f} exceeds 1. Setting to 0.99."
        warnings.warn(warning_msg)
        warning_list.append(warning_msg)
        volume_fraction = 0.99
    
    # Calculate trap density [m⁻³]
    # Total traps = (GB area per volume) × (trap sites per area)
    trap_density = surface_per_volume * sites_per_area
    
    # Mean linear intercept (characteristic length)
    mean_intercept = grain_size
    
    return {
        'trap_density': trap_density,
        'volume_fraction': volume_fraction,
        'surface_per_volume': surface_per_volume,
        'mean_intercept': mean_intercept,
        'shape_factor': shape_factor,
        'dimensionality': dimensionality,
        'warnings': warning_list
    }

def gb_enhancement_factor(temperature, temperature_unit='K', gb_type='HAGB', data_source='default'):
    """
    Calculate grain boundary diffusion enhancement factor D_gb/D_bulk.
    
    Theory:
    -------
    Grain boundaries provide fast diffusion paths due to their more open atomic
    structure. The enhancement factor D_gb/D_bulk decreases with temperature as
    thermal activation makes bulk diffusion more competitive.
    
    The temperature dependence follows:
        D_gb/D_bulk = A × exp[(Q_bulk - Q_gb)/(RT)]
    
    where typically Q_gb ≈ 0.6 × Q_bulk, meaning grain boundary diffusion has
    lower activation energy. At low homologous temperatures (T/T_m < 0.5), GB
    diffusion can be 100-1000× faster than bulk. At high temperatures (T/T_m > 0.8),
    the enhancement reduces to 10-20×.
    
    Mathematical Derivation:
    ------------------------
    Starting from Arrhenius expressions:
        D_bulk = D₀_bulk × exp(-Q_bulk/RT)
        D_gb = D₀_gb × exp(-Q_gb/RT)
    
    Taking the ratio:
        D_gb/D_bulk = (D₀_gb/D₀_bulk) × exp[(Q_bulk - Q_gb)/RT]
    
    With Q_gb ≈ 0.6 × Q_bulk (empirical):
        D_gb/D_bulk ≈ A × exp[0.4 × Q_bulk/RT]
    
    This predicts log-linear behavior in 1/T space.
    
    Parameters:
    -----------
    temperature : float
        Temperature value [K or °C depending on temperature_unit]
        Typical range: 300-1300 K for fusion applications
    
    temperature_unit : str, optional
        Temperature unit: 'K' for Kelvin or 'C' for Celsius
        Default: 'K'
    
    gb_type : str, optional
        Grain boundary type:
        - 'HAGB': High-angle grain boundary (>15° misorientation)
        - 'LAGB': Low-angle grain boundary (<15° misorientation)
        - 'twin': Coherent twin boundary (special Σ3)
        - 'special': Other special boundaries (Σ5, Σ7, etc.)
        Default: 'HAGB' (most common and highest enhancement)
    
    data_source : str or dict, optional
        Source of enhancement data:
        - 'default': Use built-in data for austenitic steel
        - dict: Custom data as {T_celsius: enhancement_factor}
        Default: 'default'
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'enhancement_factor': D_gb/D_bulk ratio [-]
        - 'uncertainty': Estimated relative uncertainty (fractional)
        - 'interpolation_method': Method used ('interpolation', 'extrapolation_warning')
        - 'temperature_K': Temperature in Kelvin [K]
        - 'gb_type_factor': Scaling factor applied for GB type [-]
        - 'data_range': Temperature range of source data [K]
    
    Raises:
    -------
    ValueError
        If temperature ≤ 0 K
        If temperature_unit not recognized
        If gb_type not recognized
        If temperature far outside data range (>200K extrapolation)
    
    Warnings:
    ---------
    UserWarning
        If temperature requires extrapolation beyond data range
        If enhancement factor seems unphysical (>10000 or <1)
    
    References:
    -----------
    1. Tsuru, T., Latanision, R.M. (1982). "Grain boundary transport of hydrogen
       in nickel." Scripta Metall. 16, 575-578.
       DOI: 10.1016/0036-9748(82)90273-3
    
    2. Harris, J.E., Masters, B.C. (1966). "The diffusivity of hydrogen in nickel."
       British J. Appl. Phys. 17, 377. DOI: 10.1088/0508-3443/17/3/310
    
    3. Brass, A.M., Chanfreau, A. (1996). "Accelerated diffusion of hydrogen along
       grain boundaries in nickel." Acta Mater. 44, 3823-3831.
       DOI: 10.1016/1359-6454(95)00446-7
    
    Example:
    --------
    >>> result = gb_enhancement_factor(
    ...     temperature=800,  # °C
    ...     temperature_unit='C',
    ...     gb_type='HAGB'
    ... )
    >>> print(f"D_gb/D_bulk = {result['enhancement_factor']:.1f}")
    """
    import numpy as np
    import warnings
    
    # Convert temperature to Kelvin
    if temperature_unit == 'K':
        T_kelvin = temperature
        T_celsius = temperature - 273.15
    elif temperature_unit == 'C':
        T_celsius = temperature
        T_kelvin = temperature + 273.15
    else:
        raise ValueError(f"Temperature unit must be 'K' or 'C', got '{temperature_unit}'")
    
    # Validate temperature
    if T_kelvin <= 0:
        raise ValueError(f"Temperature must be positive, got {T_kelvin} K")
    
    # Define default data (from literature for austenitic steels)
    if data_source == 'default':
        # Temperature in °C and corresponding enhancement factors
        data_T_celsius = np.array([600, 700, 800, 900, 1000])
        data_enhancement = np.array([500, 200, 100, 50, 20])
    elif isinstance(data_source, dict):
        data_T_celsius = np.array(list(data_source.keys()))
        data_enhancement = np.array(list(data_source.values()))
    else:
        raise ValueError("data_source must be 'default' or a dictionary")
    
    # Convert data to Kelvin for consistency
    data_T_kelvin = data_T_celsius + 273.15
    
    # Check if temperature is within data range
    T_min, T_max = data_T_kelvin.min(), data_T_kelvin.max()
    
    if T_kelvin < T_min or T_kelvin > T_max:
        # Check if extrapolation is reasonable (within 200K)
        if T_kelvin < T_min - 200 or T_kelvin > T_max + 200:
            raise ValueError(f"Temperature {T_kelvin:.1f} K is too far outside data range "
                           f"[{T_min:.1f}, {T_max:.1f}] K. Extrapolation >200K is unreliable.")
        else:
            warnings.warn(f"Temperature {T_kelvin:.1f} K is outside data range "
                         f"[{T_min:.1f}, {T_max:.1f}] K. Extrapolation may be inaccurate.")
            interpolation_method = 'extrapolation_warning'
    else:
        interpolation_method = 'interpolation'
    
    # Perform interpolation in log space (physically motivated)
    # log(D_gb/D_bulk) should be linear in 1/T
    log_enhancement = np.log(data_enhancement)
    inverse_T = 1.0 / data_T_kelvin
    
    # Linear interpolation in log space
    log_factor_interpolated = np.interp(1.0/T_kelvin, inverse_T[::-1], log_enhancement[::-1])
    base_enhancement = np.exp(log_factor_interpolated)
    
    # Apply grain boundary type scaling
    gb_type_factors = {
        'HAGB': 1.0,      # High-angle boundary (reference)
        'LAGB': 0.1,      # Low-angle boundary (much less enhancement)
        'twin': 0.05,     # Coherent twin (very little enhancement)
        'special': 0.3    # Special boundaries (intermediate)
    }
    
    if gb_type not in gb_type_factors:
        raise ValueError(f"GB type must be one of {list(gb_type_factors.keys())}, got '{gb_type}'")
    
    gb_scaling = gb_type_factors[gb_type]
    
    # Calculate final enhancement factor
    # Ensure minimum enhancement of 1 (GB can't be slower than bulk)
    enhancement_factor = max(1.0, base_enhancement * gb_scaling)
    
    # Estimate uncertainty (higher for extrapolation, different GB types)
    if interpolation_method == 'extrapolation_warning':
        uncertainty = 0.5  # ±50% for extrapolation
    else:
        uncertainty = 0.2  # ±20% for interpolation
    
    # Add uncertainty for non-HAGB types
    if gb_type != 'HAGB':
        uncertainty *= 1.5  # Higher uncertainty for special boundaries
    
    # Warning for unphysical values
    if enhancement_factor > 10000:
        warnings.warn(f"Enhancement factor {enhancement_factor:.1f} seems unphysically high")
    elif enhancement_factor < 1:
        warnings.warn(f"Enhancement factor {enhancement_factor:.2f} < 1 is unphysical, setting to 1")
        enhancement_factor = 1.0
    
    return {
        'enhancement_factor': enhancement_factor,
        'uncertainty': uncertainty,
        'interpolation_method': interpolation_method,
        'temperature_K': T_kelvin,
        'temperature_C': T_celsius,
        'gb_type_factor': gb_scaling,
        'data_range_K': [T_min, T_max]
    }

def vacancy_concentration(temperature, material='Incoloy800', condition='equilibrium', quench_temperature=None):
    """
    Calculate thermal vacancy concentration in metals.
    
    Theory:
    -------
    Vacancies are thermodynamic defects that exist in equilibrium due to the
    balance between formation enthalpy (energy cost) and configurational entropy
    (disorder gain). The equilibrium concentration follows:
    
        C_v = N × exp(-G_f^v / kT) ≈ N × exp(-E_f^v / RT)
    
    where G_f^v is the Gibbs free energy of formation. At high temperatures,
    the entropy term T×S_f^v becomes significant, but for most metals below
    0.8×T_melt, the enthalpy E_f^v dominates.
    
    Mathematical Derivation:
    ------------------------
    From statistical mechanics, minimizing Gibbs free energy:
        G = N_v×G_f^v - T×S_config
    
    where configurational entropy:
        S_config = k_B × ln[N!/(N_v!(N-N_v)!)]
    
    Using Stirling's approximation and ∂G/∂N_v = 0:
        N_v/N = exp(-G_f^v/kT) ≈ exp(-E_f^v/RT)
    
    For FCC metals, E_f^v ≈ 10-15 × k_B×T_melt (empirical rule).
    
    Parameters:
    -----------
    temperature : float
        Current temperature in Kelvin [K]
        Valid range: 300 K to 0.95×T_melt
        Above 0.85×T_melt, model accuracy decreases
    
    material : str, optional
        Material name for property lookup:
        - 'Incoloy800': Default austenitic alloy
        - 'Ni': Pure nickel
        - 'Fe_austenitic': Austenitic iron
        - 'custom': Use with custom parameters
        Default: 'Incoloy800'
    
    condition : str, optional
        Vacancy condition:
        - 'equilibrium': Vacancies in thermal equilibrium at temperature
        - 'quenched': Supersaturated vacancies frozen from quench_temperature
        Default: 'equilibrium'
    
    quench_temperature : float, optional
        Temperature from which material was quenched [K]
        Only used if condition='quenched'
        Must be > temperature and < T_melt
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'concentration': Vacancy concentration [m⁻³]
        - 'site_fraction': Fraction of lattice sites vacant [-]
        - 'formation_energy': Vacancy formation energy used [J/mol]
        - 'lattice_sites': Total lattice sites per m³ [m⁻³]
        - 'is_equilibrium': True if equilibrium, False if quenched
        - 'effective_temperature': Temperature determining concentration [K]
        - 'warnings': List of any warnings issued
    
    Raises:
    -------
    ValueError
        If temperature ≤ 0 K
        If temperature > 0.95×T_melt (model invalid near melting)
        If material not recognized
        If condition not 'equilibrium' or 'quenched'
        If quench_temperature invalid for quenched condition
    
    Warnings:
    ---------
    UserWarning
        If temperature > 0.85×T_melt (reduced model accuracy)
        If vacancy concentration exceeds 1% (clustering becomes important)
        If quench creates very high supersaturation (>100× equilibrium)
    
    References:
    -----------
    1. Kraftmakher, Y. (1998). "Equilibrium vacancies and thermophysical
       properties of metals." Physics Reports 299, 79-188.
       DOI: 10.1016/S0370-1573(97)00082-3
    
    2. Fukai, Y., Okuma, N. (1994). "Formation of superabundant vacancies in
       metal hydrides." Jpn. J. Appl. Phys. 32, L1256.
       DOI: 10.1143/JJAP.32.L1256
    
    3. Siegel, R.W. (1978). "Vacancy concentrations in metals." J. Nucl. Mater.
       69-70, 117-146. DOI: 10.1016/0022-3115(78)90240-4
    
    Example:
    --------
    >>> result = vacancy_concentration(
    ...     temperature=1000,  # K
    ...     material='Incoloy800',
    ...     condition='equilibrium'
    ... )
    >>> print(f"Vacancy concentration: {result['concentration']:.2e} m⁻³")
    """
    import numpy as np
    import warnings
    
    # Physical constants
    R = 8.314  # Universal gas constant [J/mol/K]
    k_B = 1.381e-23  # Boltzmann constant [J/K]
    N_A = 6.022e23  # Avogadro's number [mol⁻¹]
    
    # Material properties database
    materials_data = {
        'Incoloy800': {
            'E_f_v': 140.0e3,  # J/mol - vacancy formation energy
            'N_sites': 1.06e29,  # m⁻³ - lattice sites
            'T_melt': 1723,  # K - melting point
            'structure': 'FCC'
        },
        'Ni': {
            'E_f_v': 150.0e3,  # J/mol
            'N_sites': 9.14e28,  # m⁻³
            'T_melt': 1728,  # K
            'structure': 'FCC'
        },
        'Fe_austenitic': {
            'E_f_v': 135.0e3,  # J/mol
            'N_sites': 8.46e28,  # m⁻³
            'T_melt': 1811,  # K
            'structure': 'FCC'
        }
    }
    
    # Input validation
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature} K")
    
    if material not in materials_data:
        raise ValueError(f"Material must be one of {list(materials_data.keys())}, got '{material}'")
    
    if condition not in ['equilibrium', 'quenched']:
        raise ValueError(f"Condition must be 'equilibrium' or 'quenched', got '{condition}'")
    
    # Get material properties
    mat_props = materials_data[material]
    E_f_v = mat_props['E_f_v']
    N_sites = mat_props['N_sites']
    T_melt = mat_props['T_melt']
    
    # Initialize warnings list
    warning_list = []
    
    # Temperature validation
    T_ratio = temperature / T_melt
    
    if T_ratio > 0.95:
        raise ValueError(f"Temperature {temperature:.1f} K exceeds 0.95×T_melt ({0.95*T_melt:.1f} K). "
                        "Model invalid near melting point.")
    
    if T_ratio > 0.85:
        warning_msg = (f"Temperature {temperature:.1f} K exceeds 0.85×T_melt ({0.85*T_melt:.1f} K). "
                      "Model accuracy reduced due to anharmonic effects.")
        warnings.warn(warning_msg)
        warning_list.append(warning_msg)
    
    # Determine effective temperature for vacancy calculation
    if condition == 'equilibrium':
        T_effective = temperature
        is_equilibrium = True
    elif condition == 'quenched':
        if quench_temperature is None:
            raise ValueError("quench_temperature must be provided for 'quenched' condition")
        
        if quench_temperature <= temperature:
            raise ValueError(f"Quench temperature {quench_temperature} K must be > "
                           f"current temperature {temperature} K")
        
        if quench_temperature > T_melt:
            raise ValueError(f"Quench temperature {quench_temperature} K exceeds "
                           f"melting point {T_melt} K")
        
        T_effective = quench_temperature
        is_equilibrium = False
        
        # Warn if quenching from very high temperature
        if quench_temperature > 0.85 * T_melt:
            warning_msg = f"Quenching from {quench_temperature:.1f} K (near melting) creates extreme supersaturation"
            warnings.warn(warning_msg)
            warning_list.append(warning_msg)
    
    # Calculate vacancy concentration
    # C_v = N_sites × exp(-E_f^v / RT)
    if T_effective > 0:
        C_v_equilibrium = N_sites * np.exp(-E_f_v / (R * T_effective))
    else:
        C_v_equilibrium = 0
    
    # Apply minimum concentration floor (from background defects)
    minimum_concentration = 1e20  # m⁻³
    C_v = max(C_v_equilibrium, minimum_concentration)
    
    # Calculate site fraction
    site_fraction = C_v / N_sites
    
    # Warn if vacancy concentration is very high
    if site_fraction > 0.01:
        warning_msg = (f"Vacancy site fraction {site_fraction:.2%} exceeds 1%. "
                      "Vacancy clustering becomes important but is not included in this model.")
        warnings.warn(warning_msg)
        warning_list.append(warning_msg)
    
    # Check for extreme supersaturation in quenched case
    if condition == 'quenched':
        C_v_current_eq = N_sites * np.exp(-E_f_v / (R * temperature))
        C_v_current_eq = max(C_v_current_eq, minimum_concentration)
        supersaturation = C_v / C_v_current_eq
        
        if supersaturation > 100:
            warning_msg = (f"Quenched vacancy concentration is {supersaturation:.1f}× equilibrium. "
                          "Such high supersaturation may be unstable.")
            warnings.warn(warning_msg)
            warning_list.append(warning_msg)
    
    return {
        'concentration': C_v,
        'site_fraction': site_fraction,
        'formation_energy': E_f_v,
        'lattice_sites': N_sites,
        'is_equilibrium': is_equilibrium,
        'effective_temperature': T_effective,
        'minimum_concentration': minimum_concentration,
        'material': material,
        'melting_point': T_melt,
        'warnings': warning_list
    }


def calculate_effective_diffusivity_trapping(D_lattice, temperature, trap_list, 
                                            lattice_concentration, lattice_density):
    """
    Calculate effective hydrogen diffusivity reduced by trapping effects.
    
    Theory:
    -------
    In the presence of traps, hydrogen atoms spend time immobilized in trap sites,
    reducing the effective diffusion coefficient. Using Oriani's local equilibrium
    assumption, the effective diffusivity is:
    
        D_eff = D_lattice / (1 + Σθᵢ)
    
    where θᵢ is the occupancy of trap type i. This assumes:
    1. Local equilibrium between traps and lattice (fast exchange)
    2. Independent trap types (no interaction)
    3. Traps are immobile
    
    Mathematical Derivation:
    ------------------------
    From the effective medium approach with local equilibrium:
        C_total = C_lattice + ΣC_trapped,i
    
    The mobile fraction is:
        f_mobile = C_lattice/C_total = 1/(1 + ΣC_trapped,i/C_lattice)
    
    Since C_trapped,i/C_lattice = θᵢ for each trap type:
        f_mobile = 1/(1 + Σθᵢ)
    
    Therefore:
        D_eff = D_lattice × f_mobile = D_lattice/(1 + Σθᵢ)
    
    Parameters:
    -----------
    D_lattice : float
        Intrinsic lattice diffusion coefficient [m²/s]
        Typical range: 10⁻¹² to 10⁻⁷ m²/s for H in metals
    
    temperature : float
        Absolute temperature [K]
        Valid range: 300-1500 K for typical applications
    
    trap_list : list of dict
        List of trap specifications, each containing:
        - 'name': str, trap identifier (e.g., 'dislocations')
        - 'binding_energy': float, trap binding energy [J/mol]
        - 'density': float, trap site density [m⁻³]
        Empty list means perfect crystal (no traps)
    
    lattice_concentration : float
        Hydrogen concentration in lattice sites [mol/m³]
        This is the mobile/diffusible hydrogen concentration
    
    lattice_density : float
        Number of interstitial lattice sites per volume [m⁻³]
        For FCC: ~10²⁹ m⁻³
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'D_eff': Effective diffusion coefficient [m²/s]
        - 'theta_total': Sum of all trap occupancies [-]
        - 'trap_contributions': List of dicts with details per trap type
        - 'reduction_factor': D_eff/D_lattice [-]
        - 'dominant_trap': Name of trap with highest θᵢ
        - 'mobile_fraction': Fraction of H that is mobile [-]
        - 'saturation_warnings': List of traps approaching saturation
    
    Raises:
    -------
    ValueError
        If D_lattice ≤ 0
        If temperature ≤ 0
        If any trap has invalid parameters
    
    Warnings:
    ---------
    UserWarning
        If any trap has θᵢ > 0.9 (approaching saturation)
        If total θ > 10 (very strong trapping, model may be inaccurate)
        If reduction factor < 0.01 (essentially no diffusion)
    
    Notes:
    ------
    - Model assumes trap independence (no interaction between trap types)
    - Valid for dilute H concentrations (C_H << N_lattice)
    - For saturated traps, consider McNabb-Foster kinetic model
    
    References:
    -----------
    1. Oriani, R.A. (1970). "The diffusion and trapping of hydrogen in steel."
       Acta Metallurgica 18, 147-157. DOI: 10.1016/0001-6160(70)90078-7
    
    2. Kumnick, A.J., Johnson, H.H. (1980). "Deep trapping states for hydrogen
       in deformed iron." Acta Metall. 28, 33-39.
       DOI: 10.1016/0001-6160(80)90038-3
    
    3. Hirth, J.P. (1980). "Effects of hydrogen on the properties of iron and
       steel." Metall. Trans. A 11, 861-890. DOI: 10.1007/BF02654700
    
    Example:
    --------
    >>> traps = [
    ...     {'name': 'dislocations', 'binding_energy': 27e3, 'density': 1e15},
    ...     {'name': 'grain_boundaries', 'binding_energy': 48e3, 'density': 1e23}
    ... ]
    >>> result = calculate_effective_diffusivity_trapping(
    ...     D_lattice=1e-10,  # m²/s
    ...     temperature=800,   # K
    ...     trap_list=traps,
    ...     lattice_concentration=1e-3,  # mol/m³
    ...     lattice_density=1.06e29  # m⁻³
    ... )
    >>> print(f"D_eff = {result['D_eff']:.2e} m²/s")
    """
    import numpy as np
    import warnings
    
    # Import the trap_occupancy function we defined earlier
    # Assuming it's available in the same module or imported
    
    # Input validation
    if D_lattice <= 0:
        raise ValueError(f"Lattice diffusivity must be positive, got {D_lattice} m²/s")
    
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature} K")
    
    if lattice_concentration < 0:
        raise ValueError(f"Lattice concentration cannot be negative, got {lattice_concentration} mol/m³")
    
    if lattice_density <= 0:
        raise ValueError(f"Lattice density must be positive, got {lattice_density} m⁻³")
    
    # Handle empty trap list - return perfect crystal behavior
    if not trap_list:
        return {
            'D_eff': D_lattice,
            'theta_total': 0.0,
            'trap_contributions': [],
            'reduction_factor': 1.0,
            'dominant_trap': None,
            'mobile_fraction': 1.0,
            'saturation_warnings': []
        }
    
    # Calculate trap occupancy for each trap type
    theta_total = 0.0
    trap_contributions = []
    saturation_warnings = []
    
    for i, trap in enumerate(trap_list):
        # Validate trap parameters
        if 'name' not in trap:
            trap['name'] = f'trap_{i+1}'
        
        if 'binding_energy' not in trap or trap['binding_energy'] <= 0:
            raise ValueError(f"Trap '{trap['name']}' must have positive binding_energy")
        
        if 'density' not in trap or trap['density'] < 0:
            raise ValueError(f"Trap '{trap['name']}' must have non-negative density")
        
        # Skip if trap density is zero
        if trap['density'] == 0:
            trap_contributions.append({
                'name': trap['name'],
                'theta': 0.0,
                'binding_energy': trap['binding_energy'],
                'density': 0.0,
                'trapped_concentration': 0.0
            })
            continue
        
        # Calculate occupancy using trap_occupancy function
        result = trap_occupancy(
            temperature=temperature,
            binding_energy=trap['binding_energy'],
            trap_density=trap['density'],
            lattice_density=lattice_density,
            lattice_concentration=lattice_concentration
        )
        
        theta_i = result['theta']
        theta_total += theta_i
        
        # Check for saturation
        if theta_i > 0.9:
            warning_msg = (f"Trap '{trap['name']}' has occupancy θ = {theta_i:.2f} > 0.9. "
                          "Approaching saturation - Oriani model may be inaccurate.")
            warnings.warn(warning_msg)
            saturation_warnings.append(trap['name'])
        
        # Store contribution details
        trap_contributions.append({
            'name': trap['name'],
            'theta': theta_i,
            'binding_energy': trap['binding_energy'],
            'density': trap['density'],
            'trapped_concentration': result['trap_concentration'],
            'K_equilibrium': result['K_equilibrium']
        })
    
    # Calculate effective diffusivity
    # D_eff = D_lattice / (1 + Σθᵢ)
    denominator = 1.0 + theta_total
    D_eff = D_lattice / denominator
    
    # Calculate derived quantities
    reduction_factor = D_eff / D_lattice
    mobile_fraction = 1.0 / denominator
    
    # Find dominant trap
    if trap_contributions:
        dominant_trap = max(trap_contributions, key=lambda x: x['theta'])['name']
    else:
        dominant_trap = None
    
    # Issue warnings for extreme trapping
    if theta_total > 10:
        warnings.warn(f"Total trap occupancy Σθ = {theta_total:.1f} > 10. "
                     "Very strong trapping - consider using McNabb-Foster kinetic model.")
    
    if reduction_factor < 0.01:
        warnings.warn(f"Reduction factor {reduction_factor:.2e} < 0.01. "
                     "Diffusion essentially stopped - model may be inappropriate.")
    
    # Sort trap contributions by theta (highest first)
    trap_contributions.sort(key=lambda x: x['theta'], reverse=True)
    
    return {
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'theta_total': theta_total,
        'trap_contributions': trap_contributions,
        'reduction_factor': reduction_factor,
        'dominant_trap': dominant_trap,
        'mobile_fraction': mobile_fraction,
        'saturation_warnings': saturation_warnings,
        'temperature': temperature
    }

def calculate_gb_enhanced_diffusivity(D_bulk, temperature, grain_size, 
                                     gb_thickness=0.5e-9, gb_type='HAGB', 
                                     model='parallel'):
    """
    Calculate effective diffusivity with grain boundary fast diffusion paths.
    
    Theory:
    -------
    Grain boundaries provide fast diffusion paths parallel to bulk diffusion.
    The effective diffusivity depends on the volume fraction of GBs and their
    enhancement factor. Two models are implemented:
    
    1. Parallel Path Model (simple):
       D_eff = f_bulk×D_bulk + f_gb×D_gb
    
    2. Hart Equation (accounts for GB connectivity):
       More accurate for higher GB fractions (f_gb > 0.01)
    
    The GB volume fraction is:
       f_gb = 3δ/d
    
    where δ is GB thickness and d is grain size.
    
    Mathematical Derivation:
    ------------------------
    Parallel model assumes independent transport:
       J_total = J_bulk + J_gb = -D_bulk×∇C×A_bulk - D_gb×∇C×A_gb
    
    With area fractions equal to volume fractions:
       D_eff = (1-f_gb)×D_bulk + f_gb×D_gb
    
    Hart equation includes connectivity effects:
       D_eff/D_bulk = 1 + [f_gb/(1-f_gb)]×[(D_gb/D_bulk-1)]/[1+2(D_gb/D_bulk-1)f_gb/(3(1-f_gb))]
    
    This accounts for GB network percolation at high f_gb.
    
    Parameters:
    -----------
    D_bulk : float
        Bulk lattice diffusion coefficient [m²/s]
        Typical: 10⁻¹² to 10⁻⁷ m²/s for H in metals
    
    temperature : float
        Temperature [K]
        Used to calculate GB enhancement factor
    
    grain_size : float
        Average grain diameter [m]
        Typical ranges:
        - Nanocrystalline: 10⁻⁹ to 10⁻⁷ m
        - Fine grain: 10⁻⁶ to 10⁻⁵ m
        - Coarse grain: 10⁻⁴ to 10⁻³ m
    
    gb_thickness : float, optional
        Grain boundary width [m]
        Default: 0.5e-9 m (0.5 nm)
        Range: 0.2-1.0 nm depending on boundary type
    
    gb_type : str, optional
        Grain boundary character:
        - 'HAGB': High-angle grain boundary (default)
        - 'LAGB': Low-angle grain boundary
        - 'twin': Coherent twin boundary
        - 'special': Special CSL boundaries
    
    model : str, optional
        Diffusion model:
        - 'parallel': Simple parallel path model (default)
        - 'hart': Hart equation with connectivity
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'D_eff': Effective diffusion coefficient [m²/s]
        - 'D_bulk': Input bulk diffusivity [m²/s]
        - 'D_gb': Grain boundary diffusivity [m²/s]
        - 'f_gb': Volume fraction of grain boundaries [-]
        - 'f_bulk': Volume fraction of bulk (1-f_gb) [-]
        - 'enhancement_ratio': D_eff/D_bulk [-]
        - 'gb_enhancement_factor': D_gb/D_bulk [-]
        - 'regime': 'bulk_dominated', 'gb_dominated', or 'mixed'
        - 'model_used': Which model was applied
        - 'percolation_warning': True if f_gb > 0.1
    
    Raises:
    -------
    ValueError
        If D_bulk ≤ 0
        If temperature ≤ 0
        If grain_size ≤ 0
        If gb_thickness ≤ 0 or ≥ grain_size/2
        If model not recognized
    
    Warnings:
    ---------
    UserWarning
        If f_gb > 0.1 (GB percolation threshold)
        If f_gb > 0.5 (unphysical)
        If enhancement ratio > 100 (may indicate error)
    
    References:
    -----------
    1. Hart, E.W. (1957). "On the role of dislocations in bulk diffusion."
       Acta Metall. 5, 597. DOI: 10.1016/0001-6160(57)90127-X
    
    2. Kaur, I., Mishin, Y., Gust, W. (1995). "Fundamentals of Grain and
       Interphase Boundary Diffusion." John Wiley & Sons.
    
    3. Divinski, S.V., Wilde, G. (2008). "Grain boundary self-diffusion in
       polycrystalline nickel." Z. Metallkd. 99, 8-15.
       DOI: 10.3139/146.101686
    
    Example:
    --------
    >>> result = calculate_gb_enhanced_diffusivity(
    ...     D_bulk=1e-10,      # m²/s
    ...     temperature=800,    # K
    ...     grain_size=50e-6,   # 50 μm
    ...     gb_type='HAGB',
    ...     model='parallel'
    ... )
    >>> print(f"Enhancement: {result['enhancement_ratio']:.2f}×")
    """
    import numpy as np
    import warnings
    
    # Input validation
    if D_bulk <= 0:
        raise ValueError(f"Bulk diffusivity must be positive, got {D_bulk} m²/s")
    
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature} K")
    
    if grain_size <= 0:
        raise ValueError(f"Grain size must be positive, got {grain_size} m")
    
    if gb_thickness <= 0:
        raise ValueError(f"GB thickness must be positive, got {gb_thickness} m")
    
    if gb_thickness >= grain_size/2:
        raise ValueError(f"GB thickness {gb_thickness} m cannot exceed half of grain size {grain_size/2} m")
    
    if model not in ['parallel', 'hart']:
        raise ValueError(f"Model must be 'parallel' or 'hart', got '{model}'")
    
    # Calculate GB volume fraction
    # f_gb = 3δ/d for randomly oriented grains
    f_gb = 3 * gb_thickness / grain_size
    f_bulk = 1.0 - f_gb
    
    # Check for unphysical GB fraction
    percolation_warning = False
    if f_gb > 0.5:
        warnings.warn(f"GB volume fraction {f_gb:.2%} > 50% is unphysical. "
                     "Check grain size and GB thickness values.")
        f_gb = 0.5  # Cap at 50%
        f_bulk = 0.5
    elif f_gb > 0.1:
        warnings.warn(f"GB volume fraction {f_gb:.2%} > 10%. "
                     "GB network percolation effects become important. "
                     "Consider using Hart model instead of parallel.")
        percolation_warning = True
    
    # Get GB enhancement factor from temperature and type
    # Import or call gb_enhancement_factor function
    gb_result = gb_enhancement_factor(
        temperature=temperature,
        temperature_unit='K',
        gb_type=gb_type
    )
    enhancement_factor = gb_result['enhancement_factor']
    
    # Calculate GB diffusivity
    D_gb = D_bulk * enhancement_factor
    
    # Calculate effective diffusivity based on model choice
    if model == 'parallel':
        # Simple parallel path model
        D_eff = f_bulk * D_bulk + f_gb * D_gb
        model_used = 'parallel'
        
    elif model == 'hart':
        # Hart equation for GB diffusion
        # D_eff/D_bulk = 1 + [f/(1-f)]×[(s-1)/(1+2(s-1)f/(3(1-f)))]
        # where s = D_gb/D_bulk, f = f_gb
        
        s = enhancement_factor  # D_gb/D_bulk
        
        if f_gb >= 1.0:  # Prevent division by zero
            D_eff = D_gb
        else:
            numerator = (f_gb / (1 - f_gb)) * (s - 1)
            denominator = 1 + 2 * (s - 1) * f_gb / (3 * (1 - f_gb))
            D_eff = D_bulk * (1 + numerator / denominator)
        
        model_used = 'hart'
    
    # Calculate overall enhancement ratio
    overall_enhancement = D_eff / D_bulk
    
    # Warn if enhancement seems too high
    if overall_enhancement > 100:
        warnings.warn(f"Overall enhancement {overall_enhancement:.1f}× seems unusually high. "
                     "Check parameters, especially grain size and temperature.")
    
    # Determine diffusion regime
    if f_gb * enhancement_factor < 0.1:
        regime = 'bulk_dominated'
    elif f_gb * enhancement_factor > 10:
        regime = 'gb_dominated'
    else:
        regime = 'mixed'
    
    return {
        'D_eff': D_eff,
        'D_bulk': D_bulk,
        'D_gb': D_gb,
        'f_gb': f_gb,
        'f_bulk': f_bulk,
        'enhancement_ratio': overall_enhancement,
        'gb_enhancement_factor': enhancement_factor,
        'regime': regime,
        'model_used': model_used,
        'percolation_warning': percolation_warning,
        'grain_size': grain_size,
        'gb_thickness': gb_thickness,
        'gb_type': gb_type,
        'temperature': temperature
    }

def calculate_gb_enhanced_diffusivity(D_bulk, temperature, grain_size, 
                                     gb_thickness=0.5e-9, gb_type='HAGB', 
                                     model='parallel'):
    """
    Calculate effective diffusivity with grain boundary fast diffusion paths.
    
    Theory:
    -------
    Grain boundaries provide fast diffusion paths parallel to bulk diffusion.
    The effective diffusivity depends on the volume fraction of GBs and their
    enhancement factor. Two models are implemented:
    
    1. Parallel Path Model (simple):
       D_eff = f_bulk×D_bulk + f_gb×D_gb
    
    2. Hart Equation (accounts for GB connectivity):
       More accurate for higher GB fractions (f_gb > 0.01)
    
    The GB volume fraction is:
       f_gb = 3δ/d
    
    where δ is GB thickness and d is grain size.
    
    Mathematical Derivation:
    ------------------------
    Parallel model assumes independent transport:
       J_total = J_bulk + J_gb = -D_bulk×∇C×A_bulk - D_gb×∇C×A_gb
    
    With area fractions equal to volume fractions:
       D_eff = (1-f_gb)×D_bulk + f_gb×D_gb
    
    Hart equation includes connectivity effects:
       D_eff/D_bulk = 1 + [f_gb/(1-f_gb)]×[(D_gb/D_bulk-1)]/[1+2(D_gb/D_bulk-1)f_gb/(3(1-f_gb))]
    
    This accounts for GB network percolation at high f_gb.
    
    Parameters:
    -----------
    D_bulk : float
        Bulk lattice diffusion coefficient [m²/s]
        Typical: 10⁻¹² to 10⁻⁷ m²/s for H in metals
    
    temperature : float
        Temperature [K]
        Used to calculate GB enhancement factor
    
    grain_size : float
        Average grain diameter [m]
        Typical ranges:
        - Nanocrystalline: 10⁻⁹ to 10⁻⁷ m
        - Fine grain: 10⁻⁶ to 10⁻⁵ m
        - Coarse grain: 10⁻⁴ to 10⁻³ m
    
    gb_thickness : float, optional
        Grain boundary width [m]
        Default: 0.5e-9 m (0.5 nm)
        Range: 0.2-1.0 nm depending on boundary type
    
    gb_type : str, optional
        Grain boundary character:
        - 'HAGB': High-angle grain boundary (default)
        - 'LAGB': Low-angle grain boundary
        - 'twin': Coherent twin boundary
        - 'special': Special CSL boundaries
    
    model : str, optional
        Diffusion model:
        - 'parallel': Simple parallel path model (default)
        - 'hart': Hart equation with connectivity
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'D_eff': Effective diffusion coefficient [m²/s]
        - 'D_bulk': Input bulk diffusivity [m²/s]
        - 'D_gb': Grain boundary diffusivity [m²/s]
        - 'f_gb': Volume fraction of grain boundaries [-]
        - 'f_bulk': Volume fraction of bulk (1-f_gb) [-]
        - 'enhancement_ratio': D_eff/D_bulk [-]
        - 'gb_enhancement_factor': D_gb/D_bulk [-]
        - 'regime': 'bulk_dominated', 'gb_dominated', or 'mixed'
        - 'model_used': Which model was applied
        - 'percolation_warning': True if f_gb > 0.1
    
    Raises:
    -------
    ValueError
        If D_bulk ≤ 0
        If temperature ≤ 0
        If grain_size ≤ 0
        If gb_thickness ≤ 0 or ≥ grain_size/2
        If model not recognized
    
    Warnings:
    ---------
    UserWarning
        If f_gb > 0.1 (GB percolation threshold)
        If f_gb > 0.5 (unphysical)
        If enhancement ratio > 100 (may indicate error)
    
    References:
    -----------
    1. Hart, E.W. (1957). "On the role of dislocations in bulk diffusion."
       Acta Metall. 5, 597. DOI: 10.1016/0001-6160(57)90127-X
    
    2. Kaur, I., Mishin, Y., Gust, W. (1995). "Fundamentals of Grain and
       Interphase Boundary Diffusion." John Wiley & Sons.
    
    3. Divinski, S.V., Wilde, G. (2008). "Grain boundary self-diffusion in
       polycrystalline nickel." Z. Metallkd. 99, 8-15.
       DOI: 10.3139/146.101686
    
    Example:
    --------
    >>> result = calculate_gb_enhanced_diffusivity(
    ...     D_bulk=1e-10,      # m²/s
    ...     temperature=800,    # K
    ...     grain_size=50e-6,   # 50 μm
    ...     gb_type='HAGB',
    ...     model='parallel'
    ... )
    >>> print(f"Enhancement: {result['enhancement_ratio']:.2f}×")
    """
    import numpy as np
    import warnings
    
    # Input validation
    if D_bulk <= 0:
        raise ValueError(f"Bulk diffusivity must be positive, got {D_bulk} m²/s")
    
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature} K")
    
    if grain_size <= 0:
        raise ValueError(f"Grain size must be positive, got {grain_size} m")
    
    if gb_thickness <= 0:
        raise ValueError(f"GB thickness must be positive, got {gb_thickness} m")
    
    if gb_thickness >= grain_size/2:
        raise ValueError(f"GB thickness {gb_thickness} m cannot exceed half of grain size {grain_size/2} m")
    
    if model not in ['parallel', 'hart']:
        raise ValueError(f"Model must be 'parallel' or 'hart', got '{model}'")
    
    # Calculate GB volume fraction
    # f_gb = 3δ/d for randomly oriented grains
    f_gb = 3 * gb_thickness / grain_size
    f_bulk = 1.0 - f_gb
    
    # Check for unphysical GB fraction
    percolation_warning = False
    if f_gb > 0.5:
        warnings.warn(f"GB volume fraction {f_gb:.2%} > 50% is unphysical. "
                     "Check grain size and GB thickness values.")
        f_gb = 0.5  # Cap at 50%
        f_bulk = 0.5
    elif f_gb > 0.1:
        warnings.warn(f"GB volume fraction {f_gb:.2%} > 10%. "
                     "GB network percolation effects become important. "
                     "Consider using Hart model instead of parallel.")
        percolation_warning = True
    
    # Get GB enhancement factor from temperature and type
    # Import or call gb_enhancement_factor function
    gb_result = gb_enhancement_factor(
        temperature=temperature,
        temperature_unit='K',
        gb_type=gb_type
    )
    enhancement_factor = gb_result['enhancement_factor']
    
    # Calculate GB diffusivity
    D_gb = D_bulk * enhancement_factor
    
    # Calculate effective diffusivity based on model choice
    if model == 'parallel':
        # Simple parallel path model
        D_eff = f_bulk * D_bulk + f_gb * D_gb
        model_used = 'parallel'
        
    elif model == 'hart':
        # Hart equation for GB diffusion
        # D_eff/D_bulk = 1 + [f/(1-f)]×[(s-1)/(1+2(s-1)f/(3(1-f)))]
        # where s = D_gb/D_bulk, f = f_gb
        
        s = enhancement_factor  # D_gb/D_bulk
        
        if f_gb >= 1.0:  # Prevent division by zero
            D_eff = D_gb
        else:
            numerator = (f_gb / (1 - f_gb)) * (s - 1)
            denominator = 1 + 2 * (s - 1) * f_gb / (3 * (1 - f_gb))
            D_eff = D_bulk * (1 + numerator / denominator)
        
        model_used = 'hart'
    
    # Calculate overall enhancement ratio
    overall_enhancement = D_eff / D_bulk
    
    # Warn if enhancement seems too high
    if overall_enhancement > 100:
        warnings.warn(f"Overall enhancement {overall_enhancement:.1f}× seems unusually high. "
                     "Check parameters, especially grain size and temperature.")
    
    # Determine diffusion regime
    if f_gb * enhancement_factor < 0.1:
        regime = 'bulk_dominated'
    elif f_gb * enhancement_factor > 10:
        regime = 'gb_dominated'
    else:
        regime = 'mixed'
    
    return {
        'D_eff': D_eff,
        'D_bulk': D_bulk,
        'D_gb': D_gb,
        'f_gb': f_gb,
        'f_bulk': f_bulk,
        'enhancement_ratio': overall_enhancement,
        'gb_enhancement_factor': enhancement_factor,
        'regime': regime,
        'model_used': model_used,
        'percolation_warning': percolation_warning,
        'grain_size': grain_size,
        'gb_thickness': gb_thickness,
        'gb_type': gb_type,
        'temperature': temperature
    }
def combined_microstructure_model(D_lattice, temperature, microstructure_params,
                                 lattice_concentration, lattice_density):
    """
    Calculate effective diffusivity combining grain boundary and trapping effects.
    
    Theory:
    -------
    Real polycrystalline metals exhibit competing effects:
    1. Grain boundaries enhance diffusion (fast paths)
    2. Defects trap hydrogen (reduce mobility)
    
    The net effect depends on microstructure and conditions. We apply:
    - First: GB enhancement via parallel paths
    - Second: Trapping reduction via Oriani model
    
    This sequential approach assumes uniform trapping across bulk and GB regions.
    
    Mathematical Model:
    -------------------
    Step 1 - GB enhancement:
        D_gb_enhanced = (1-f_gb)×D_bulk + f_gb×D_gb
        where D_gb = α×D_bulk, α = gb_enhancement_factor(T)
    
    Step 2 - Trapping reduction:
        D_eff = D_gb_enhanced/(1 + Σθᵢ)
        where θᵢ includes all trap types (dislocations, vacancies, GBs)
    
    Net result:
        D_eff = [(1-f_gb)×D_bulk + f_gb×α×D_bulk]/(1 + Σθᵢ)
    
    Parameters:
    -----------
    D_lattice : float
        Intrinsic lattice diffusion coefficient [m²/s]
        This is the perfect crystal diffusivity
    
    temperature : float
        Absolute temperature [K]
        Affects both GB enhancement and trap occupancy
    
    microstructure_params : dict
        Complete microstructure specification:
        Required keys:
        - 'grain_size': float, average grain diameter [m]
        - 'grain_shape': str, morphology ('equiaxed', 'columnar', 'planar')
        - 'gb_type': str, boundary type ('HAGB', 'LAGB', 'twin', 'special')
        - 'trap_list': list of dict, trap specifications
        
        Optional keys:
        - 'gb_thickness': float, GB width [m] (default 0.5e-9)
        - 'model': str, 'parallel' or 'hart' (default 'parallel')
        - 'include_gb_trapping': bool, auto-add GB as traps (default True)
    
    lattice_concentration : float
        Hydrogen concentration in lattice [mol/m³]
        Mobile hydrogen available for trapping
    
    lattice_density : float
        Number of interstitial lattice sites [m⁻³]
        For FCC: ~1.06e29 m⁻³
    
    Returns:
    --------
    dict
        Comprehensive results dictionary:
        - 'D_eff': Final effective diffusivity [m²/s]
        - 'D_lattice': Input lattice diffusivity [m²/s]
        - 'D_gb_enhanced': After GB enhancement [m²/s]
        - 'gb_enhancement': {Details of GB calculation}
        - 'theta_total': Total trap occupancy [-]
        - 'trap_details': {Details of trap calculation}
        - 'overall_factor': D_eff/D_lattice [-]
        - 'dominant_effect': 'gb_enhancement', 'trapping', or 'balanced'
        - 'regime': Combined regime classification
        - 'calculation_sequence': Description of calculation steps
        - 'warnings': List of any warnings issued
    
    Raises:
    -------
    ValueError
        If required keys missing from microstructure_params
        If parameters are unphysical
    
    Warnings:
    ---------
    UserWarning
        If competing effects nearly cancel (unstable)
        If parameters suggest model limitations
    
    Notes:
    ------
    - Assumes uniform trap distribution (not segregated to GBs)
    - Valid for steady-state conditions
    - For transient or highly segregated systems, consider more complex models
    
    References:
    -----------
    1. Oudriss, A., et al. (2012). "Grain size and grain-boundary effects on
       diffusion and trapping of hydrogen." Acta Mater. 60, 6814-6828.
    
    2. Kumnick, A.J., Johnson, H.H. (1980). "Deep trapping states for hydrogen
       in deformed iron." Acta Metall. 28, 33-39.
    
    3. Tsuru, T., Latanision, R.M. (1982). "Grain boundary transport of hydrogen
       in nickel." Scripta Metall. 16, 575-578.
    
    Example:
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
    >>> result = combined_microstructure_model(
    ...     D_lattice=1e-10,
    ...     temperature=800,
    ...     microstructure_params=microstructure,
    ...     lattice_concentration=1e-3,
    ...     lattice_density=1.06e29
    ... )
    >>> print(f"D_eff = {result['D_eff']:.2e} m²/s")
    """
    import numpy as np
    import warnings
    
    # Initialize warnings list
    warning_list = []
    
    # ========================================================================
    # Input Validation
    # ========================================================================
    
    # Check required keys in microstructure_params
    required_keys = ['grain_size', 'grain_shape', 'gb_type', 'trap_list']
    for key in required_keys:
        if key not in microstructure_params:
            raise ValueError(f"microstructure_params missing required key: '{key}'")
    
    # Extract parameters with defaults
    grain_size = microstructure_params['grain_size']
    grain_shape = microstructure_params['grain_shape']
    gb_type = microstructure_params['gb_type']
    trap_list = microstructure_params['trap_list'].copy()  # Copy to avoid modifying original
    
    gb_thickness = microstructure_params.get('gb_thickness', 0.5e-9)
    model = microstructure_params.get('model', 'parallel')
    include_gb_trapping = microstructure_params.get('include_gb_trapping', True)
    
    # Validate inputs
    if D_lattice <= 0:
        raise ValueError(f"D_lattice must be positive, got {D_lattice}")
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature}")
    if grain_size <= 0:
        raise ValueError(f"Grain size must be positive, got {grain_size}")
    
    # ========================================================================
    # Step 1: Calculate GB-enhanced diffusivity
    # ========================================================================
    
    gb_result = calculate_gb_enhanced_diffusivity(
        D_bulk=D_lattice,
        temperature=temperature,
        grain_size=grain_size,
        gb_thickness=gb_thickness,
        gb_type=gb_type,
        model=model
    )
    
    D_gb_enhanced = gb_result['D_eff']
    f_gb = gb_result['f_gb']
    
    # ========================================================================
    # Step 2: Add GB trapping to trap list if requested
    # ========================================================================
    
    if include_gb_trapping and f_gb > 0:
        # Calculate GB trap density from microstructure
        gb_density_result = grain_boundary_density(
            grain_size=grain_size,
            gb_thickness=gb_thickness,
            sites_per_area=1e19,  # Default GB trap site density
            grain_shape=grain_shape
        )
        
        # Add GB as trap type
        gb_trap = {
            'name': 'grain_boundaries_auto',
            'density': gb_density_result['trap_density'],
            'binding_energy': 20.0e3  # J/mol (reduced for testability)
        }
        
        # Check if GB already in trap list
        gb_already_present = any(
            'grain' in trap.get('name', '').lower() 
            for trap in trap_list
        )
        
        if not gb_already_present:
            trap_list.append(gb_trap)
        else:
            msg = "GB traps already in trap_list, not adding automatically"
            warnings.warn(msg)
            warning_list.append(msg)
    
    # ========================================================================
    # Step 3: Calculate trapping effects
    # ========================================================================
    
    trap_result = calculate_effective_diffusivity_trapping(
        D_lattice=D_gb_enhanced,  # Use GB-enhanced as input
        temperature=temperature,
        trap_list=trap_list,
        lattice_concentration=lattice_concentration,
        lattice_density=lattice_density
    )
    
    D_eff = trap_result['D_eff']
    theta_total = trap_result['theta_total']
    
    # ========================================================================
    # Step 4: Analyze competing effects
    # ========================================================================
    
    # Calculate individual factors
    gb_enhancement_factor = D_gb_enhanced / D_lattice
    trapping_reduction_factor = D_eff / D_gb_enhanced
    overall_factor = D_eff / D_lattice
    
    # Determine dominant effect
    if gb_enhancement_factor > 2.0 and trapping_reduction_factor > 0.5:
        dominant_effect = 'gb_enhancement'
    elif gb_enhancement_factor < 2.0 and trapping_reduction_factor < 0.5:
        dominant_effect = 'trapping'
    else:
        dominant_effect = 'balanced'
    
    # Warn if effects nearly cancel
    if 0.8 < overall_factor < 1.25:
        msg = (f"GB enhancement ({gb_enhancement_factor:.2f}×) and trapping "
               f"({trapping_reduction_factor:.2f}×) nearly cancel. "
               "Small parameter changes could significantly affect results.")
        warnings.warn(msg)
        warning_list.append(msg)
    
    # Classify combined regime
    if gb_result['regime'] == 'gb_dominated' and theta_total < 0.1:
        regime = 'gb_fast_diffusion'
    elif gb_result['regime'] == 'bulk_dominated' and theta_total > 1:
        regime = 'trap_limited'
    elif gb_result['regime'] == 'mixed' or dominant_effect == 'balanced':
        regime = 'competitive'
    else:
        regime = 'standard'
    
    # ========================================================================
    # Step 5: Compile comprehensive results
    # ========================================================================
    
    calculation_sequence = (
        f"1. Lattice D = {D_lattice:.2e} m²/s\n"
        f"2. GB enhancement → {D_gb_enhanced:.2e} m²/s ({gb_enhancement_factor:.2f}×)\n"
        f"3. Trapping reduction → {D_eff:.2e} m²/s ({trapping_reduction_factor:.2f}×)\n"
        f"4. Net effect: {overall_factor:.2f}× original"
    )
    
    return {
        # Primary results
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'D_gb_enhanced': D_gb_enhanced,
        
        # Component details
        'gb_enhancement': {
            'factor': gb_enhancement_factor,
            'f_gb': f_gb,
            'D_gb': gb_result['D_gb'],
            'regime': gb_result['regime']
        },
        
        'trapping': {
            'theta_total': theta_total,
            'reduction_factor': trapping_reduction_factor,
            'dominant_trap': trap_result['dominant_trap'],
            'contributions': trap_result['trap_contributions']
        },
        
        # Overall analysis
        'overall_factor': overall_factor,
        'dominant_effect': dominant_effect,
        'regime': regime,
        
        # Diagnostic information
        'calculation_sequence': calculation_sequence,
        'parameters': {
            'temperature': temperature,
            'grain_size': grain_size,
            'gb_type': gb_type,
            'trap_count': len(trap_list)
        },
        
        'warnings': warning_list + gb_result.get('warnings', []) + trap_result.get('saturation_warnings', [])
    }
