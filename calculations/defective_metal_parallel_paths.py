"""
Defective Metal Model - Parallel Paths Implementation (Level 4 Alternative)

This module implements an ALTERNATIVE approach to microstructure effects where
ALL defect types (grain boundaries, dislocations, etc.) are modeled as parallel
diffusion paths, similar to the oxide defect model.

Key Difference from defective_metal.py:
---------------------------------------
- **Original (defective_metal.py):** GB parallel + all traps homogenized
  D_eff = [(1-f_gb)×D_L + f_gb×D_gb] / (1 + Σθᵢ)

- **This file:** ALL paths parallel, each with its own trapping
  D_eff = Σ[f_i × D_i/(1+θᵢ)]

Usage:
------
To switch from homogenized to parallel path model, simply change:
    from calculations.defective_metal import combined_microstructure_model
to:
    from calculations.defective_metal_parallel_paths import combined_microstructure_model

All function signatures remain identical!

When to Use This Model:
-----------------------
1. Nanocrystalline materials (grain size < 1 μm)
2. Very high dislocation densities (>10¹⁶ m⁻³) forming networks
3. Severe plastic deformation microstructures
4. When consistency with oxide parallel path model is desired

When to Use Original Model:
----------------------------
1. Conventional grain sizes (>10 μm)
2. Moderate dislocation densities
3. Standard metallurgical conditions
4. Validated against most literature data

Theory:
-------
Each microstructural feature provides an independent diffusion pathway:

1. Lattice path: Perfect crystal diffusion (slowest, largest volume)
2. Grain boundary path: Fast diffusion along GBs (2D network)
3. Dislocation path: Pipe diffusion along dislocation cores (1D)
4. Carbide interface path: Diffusion along precipitate/matrix interfaces

Volume conservation: f_lattice + f_gb + f_disloc + ... = 1

Each path experiences its own trapping environment, reducing its contribution.

Mathematical Framework:
-----------------------
For N parallel paths:
    D_eff = Σ(i=1 to N) [f_i × D_i / (1 + θ_i)]

where:
    f_i = volume fraction of path i
    D_i = intrinsic diffusivity of path i
    θ_i = trap occupancy specific to path i

References:
-----------
1. Hart, E.W. (1957). "On the role of dislocations in bulk diffusion."
   Acta Metallurgica 5, 597-606.

2. Divinski, S.V., et al. (2011). "Grain boundary self-diffusion in 
   polycrystalline nickel of different purity levels." Acta Mater. 59, 
   1974-1985.

3. Oudriss, A., et al. (2012). "Grain size and grain-boundary effects on
   diffusion and trapping of hydrogen in pure nickel." Acta Mater. 60, 
   6814-6828.

Limitations:
------------
1. Requires geometric models for each path type
2. Need diffusivity data for each path (scarce for dislocations)
3. More parameters than homogenized model
4. Less validated against experimental data
5. Assumes non-interacting parallel paths (no cross-path trapping)

Author: Auto-generated parallel path implementation
Date: January 16, 2026
"""

import numpy as np
import warnings


# ============================================================================
# CORE FUNCTIONS (Same signatures as defective_metal.py)
# ============================================================================

def trap_occupancy(temperature, binding_energy, trap_density, lattice_density, 
                   lattice_concentration):
    """
    Calculate hydrogen trap occupancy fraction using Oriani local equilibrium model.
    
    [Same implementation as original - trapping physics unchanged]
    """
    import numpy as np
    
    # Physical constants
    R = 8.314  # Universal gas constant [J/mol/K]
    
    # Input validation
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature} K")
    
    if binding_energy < 0:
        raise ValueError(f"Binding energy cannot be negative, got {binding_energy} J/mol")
    
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
        raise ValueError(f"Lattice concentration ({lattice_concentration:.2e} mol/m³) exceeds "
                        f"lattice site density ({lattice_density:.2e} sites/m³)")
    
    # Calculate equilibrium constant
    K_eq = np.exp(binding_energy / (R * temperature))
    
    # Calculate lattice site occupancy
    theta_L = C_L_sites / lattice_density
    
    # Determine which equation to use
    K_theta_L = K_eq * theta_L
    
    if K_theta_L < 0.01:  # Low coverage - use simplified form
        theta_T = K_theta_L
        approximation = 'low_coverage'
    else:
        # Full Oriani equation
        theta_T = K_theta_L / (1.0 + K_theta_L)
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
    
    [Same implementation as original]
    """
    import numpy as np
    
    # Input validation
    if grain_size <= 0:
        raise ValueError(f"Grain size must be positive, got {grain_size} m")
    
    if gb_thickness <= 0:
        raise ValueError(f"GB thickness must be positive, got {gb_thickness} m")
    
    if gb_thickness > grain_size/10:
        raise ValueError(f"GB thickness ({gb_thickness:.2e} m) too large for grain size ({grain_size:.2e} m)")
    
    if sites_per_area <= 0:
        raise ValueError(f"Sites per area must be positive, got {sites_per_area} m⁻²")
    
    valid_shapes = ['equiaxed', 'columnar', 'planar']
    if grain_shape not in valid_shapes:
        raise ValueError(f"grain_shape must be one of {valid_shapes}, got '{grain_shape}'")
    
    # Initialize warnings list
    warning_list = []
    
    # Check for unusual grain sizes
    if grain_size < 10e-9:  # Less than 10 nm
        warning_list.append(f"Very small grain size ({grain_size*1e9:.1f} nm) - approaching amorphous limit")
    
    if grain_size > 10e-3:  # Greater than 10 mm
        warning_list.append(f"Very large grain size ({grain_size*1e3:.1f} mm) - unusually coarse")
    
    # Calculate surface area per volume based on grain shape
    if grain_shape == 'equiaxed':
        shape_factor = 3.0  # 3D random polycrystal
        dimensionality = 3
    elif grain_shape == 'columnar':
        shape_factor = 2.0  # 2D (no boundaries perpendicular to columnar direction)
        dimensionality = 2
    elif grain_shape == 'planar':
        shape_factor = 1.0  # 1D (lamellar structure)
        dimensionality = 1
    
    # Calculate GB surface area per unit volume [m⁻¹]
    surface_per_volume = shape_factor / grain_size
    
    # Calculate volume fraction of GB
    volume_fraction = surface_per_volume * gb_thickness
    
    # Warn if GB volume fraction is very high
    if volume_fraction > 0.1:
        warning_list.append(f"GB volume fraction ({volume_fraction:.1%}) > 10% - percolation effects may be important")
    
    # Ensure volume fraction doesn't exceed 1 (unphysical)
    if volume_fraction > 1.0:
        warning_list.append(f"GB volume fraction ({volume_fraction:.2f}) > 1 - unphysical! Check grain size and GB thickness.")
        volume_fraction = min(volume_fraction, 0.99)  # Cap at 99%
    
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
    
    [Same implementation as original]
    """
    import numpy as np
    
    # Convert temperature to Kelvin
    if temperature_unit == 'K':
        T_kelvin = temperature
        T_celsius = temperature - 273.15
    elif temperature_unit == 'C':
        T_celsius = temperature
        T_kelvin = temperature + 273.15
    else:
        raise ValueError(f"temperature_unit must be 'K' or 'C', got '{temperature_unit}'")
    
    # Validate temperature
    if T_kelvin <= 0:
        raise ValueError(f"Temperature must be positive, got {T_kelvin} K")
    
    # Define default data (from literature for austenitic steels)
    if data_source == 'default':
        data_T_celsius = np.array([600, 700, 800, 900, 1000])
        data_enhancement = np.array([500, 200, 100, 50, 20])
    elif isinstance(data_source, dict):
        data_T_celsius = np.array(list(data_source.keys()))
        data_enhancement = np.array(list(data_source.values()))
    else:
        raise ValueError("data_source must be 'default' or dict")
    
    # Convert data to Kelvin for consistency
    data_T_kelvin = data_T_celsius + 273.15
    
    # Check if temperature is within data range
    T_min, T_max = data_T_kelvin.min(), data_T_kelvin.max()
    
    if T_kelvin < T_min or T_kelvin > T_max:
        warnings.warn(f"Temperature {T_kelvin} K outside data range [{T_min}, {T_max}] K. "
                     "Extrapolating - results may be unreliable.")
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
        'HAGB': 1.0,      # High-angle grain boundary (reference)
        'LAGB': 0.1,      # Low-angle grain boundary (much lower)
        'twin': 0.05,     # Coherent twin boundary (very low)
        'special': 0.3    # Special boundaries (Σ5, Σ7, etc.)
    }
    
    if gb_type not in gb_type_factors:
        raise ValueError(f"gb_type must be one of {list(gb_type_factors.keys())}, got '{gb_type}'")
    
    gb_scaling = gb_type_factors[gb_type]
    
    # Calculate final enhancement factor
    # Ensure minimum enhancement of 1 (GB can't be slower than bulk)
    enhancement_factor = max(1.0, base_enhancement * gb_scaling)
    
    # Estimate uncertainty (higher for extrapolation, different GB types)
    if interpolation_method == 'extrapolation_warning':
        uncertainty = 0.5  # 50% uncertainty for extrapolation
    else:
        uncertainty = 0.2  # 20% uncertainty for interpolation
    
    # Add uncertainty for non-HAGB types
    if gb_type != 'HAGB':
        uncertainty += 0.3  # Additional uncertainty for scaling
    
    # Warning for unphysical values
    if enhancement_factor > 10000:
        warnings.warn(f"Enhancement factor {enhancement_factor:.0f} > 10000 - seems unphysically high")
    elif enhancement_factor < 1:
        warnings.warn(f"Enhancement factor {enhancement_factor:.2f} < 1 was capped at 1 "
                     "(GB can't be slower than bulk)")
    
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
    
    [Same implementation as original]
    """
    import numpy as np
    
    # Physical constants
    R = 8.314  # Universal gas constant [J/mol/K]
    k_B = 1.381e-23  # Boltzmann constant [J/K]
    N_A = 6.022e23  # Avogadro's number [mol⁻¹]
    
    # Material properties database
    materials_data = {
        'Incoloy800': {
            'E_f_v': 140.0e3,  # J/mol - vacancy formation energy
            'N_sites': 1.06e29,  # m⁻³ - lattice site density
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
        raise ValueError(f"condition must be 'equilibrium' or 'quenched', got '{condition}'")
    
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
        raise ValueError(f"Temperature {temperature} K > 0.95×T_melt ({0.95*T_melt:.0f} K) - model invalid near melting")
    
    if T_ratio > 0.85:
        warning_list.append(f"Temperature {temperature} K > 0.85×T_melt ({0.85*T_melt:.0f} K) - "
                           "reduced model accuracy, entropy effects become significant")
    
    # Determine effective temperature for vacancy calculation
    if condition == 'equilibrium':
        T_effective = temperature
        is_equilibrium = True
    elif condition == 'quenched':
        if quench_temperature is None:
            raise ValueError("quench_temperature must be provided when condition='quenched'")
        if quench_temperature <= temperature:
            raise ValueError(f"quench_temperature ({quench_temperature} K) must be > temperature ({temperature} K)")
        if quench_temperature > T_melt:
            raise ValueError(f"quench_temperature ({quench_temperature} K) must be < T_melt ({T_melt} K)")
        
        T_effective = quench_temperature
        is_equilibrium = False
        warning_list.append(f"Using non-equilibrium vacancy concentration frozen from {quench_temperature} K")
    
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
        warning_list.append(f"Vacancy site fraction {site_fraction:.2%} > 1% - "
                           "clustering and multi-vacancy complexes may be important")
    
    # Check for extreme supersaturation in quenched case
    if condition == 'quenched':
        C_v_eq_at_T = N_sites * np.exp(-E_f_v / (R * temperature))
        supersaturation = C_v / max(C_v_eq_at_T, minimum_concentration)
        if supersaturation > 100:
            warning_list.append(f"Supersaturation ratio {supersaturation:.0f}× - "
                               "vacancies may anneal out quickly")
    
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


# ============================================================================
# PARALLEL PATH SPECIFIC FUNCTIONS
# ============================================================================

def calculate_dislocation_path_properties(dislocation_density, D_lattice, temperature):
    """
    Calculate dislocation pipe diffusion contribution.
    
    Theory:
    -------
    Dislocations provide 1D fast diffusion paths ("pipe diffusion").
    The dislocation core has enhanced diffusivity due to disordered structure.
    
    Volume fraction: f_disloc = π × r_pipe² × ρ_disloc
    
    Parameters:
    -----------
    dislocation_density : float
        Dislocation line density [m⁻²]
    D_lattice : float
        Lattice diffusivity [m²/s]
    temperature : float
        Temperature [K]
    
    Returns:
    --------
    dict with 'fraction', 'D_intrinsic', 'name'
    """
    # Dislocation core radius (typical value)
    r_pipe = 0.5e-9  # m (0.5 nm)
    
    # Volume fraction occupied by dislocation cores
    f_disloc = np.pi * r_pipe**2 * dislocation_density
    
    # Dislocation core diffusivity enhancement
    # Typically 10-100× faster than lattice, but slower than GB
    # Use temperature-dependent model similar to GB
    Q_ratio = 0.7  # Q_disloc/Q_lattice (between lattice and GB)
    enhancement_base = 50  # at reference temperature
    
    # Simple temperature scaling (could be made more sophisticated)
    D_disloc = D_lattice * enhancement_base
    
    return {
        'name': 'dislocation_pipe',
        'fraction': f_disloc,
        'D_intrinsic': D_disloc,
        'core_radius': r_pipe
    }


def calculate_effective_diffusivity_trapping(D_lattice, temperature, trap_list, 
                                            lattice_concentration, lattice_density):
    """
    Calculate effective diffusivity for a SINGLE PATH with trapping.
    
    This function is path-specific in the parallel model.
    
    Uses the CORRECT Oriani equilibrium model:
        D_eff = D_lattice / (1 + Σ(N_T,i × K_i / N_L))
    
    where:
        - N_T,i = trap density for trap type i [m⁻³]
        - K_i = exp(E_b,i / RT) = equilibrium constant for trap i
        - N_L = lattice site density [m⁻³]
        - The term (N_T × K / N_L) represents the ratio of trapped to mobile H
    
    Parameters
    ----------
    D_lattice : float
        Lattice diffusivity [m²/s]
    temperature : float
        Temperature [K]
    trap_list : list of dict
        Each trap: {'name': str, 'binding_energy': float [J/mol], 'density': float [m⁻³]}
    lattice_concentration : float
        Mobile H concentration in lattice [mol/m³]
    lattice_density : float
        Lattice interstitial site density [m⁻³], typically ~1e29 for FCC metals
    
    Returns
    -------
    dict with:
        - 'D_eff': Effective diffusivity [m²/s]
        - 'theta_total': Sum of trap occupancy fractions (for info)
        - 'trapping_term': Sum of (N_T × K / N_L) - the actual denominator factor
        - 'trap_contributions': List of per-trap details including trapping_contribution
        - 'reduction_factor': D_eff / D_lattice
        - 'dominant_trap': Name of trap with largest trapping_contribution
        - 'mobile_fraction': Fraction of H that is mobile
        - 'saturation_warnings': List of traps near saturation
    """
    import numpy as np
    import warnings
    
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
            'trapping_term': 0.0,
            'trap_contributions': [],
            'reduction_factor': 1.0,
            'dominant_trap': None,
            'mobile_fraction': 1.0,
            'saturation_warnings': []
        }
    
    # Calculate trap occupancy for each trap type
    theta_total = 0.0
    trapping_term = 0.0  # NEW: Sum of (N_T × K / N_L) for correct Oriani formula
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
                'trapped_concentration': 0.0,
                'trapping_contribution': 0.0
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
        
        # Calculate the TRAPPING CONTRIBUTION using correct Oriani formula
        # trapping_contribution = (N_T × K) / N_L
        K_eq = result['K_equilibrium']
        trapping_contribution = (trap['density'] * K_eq) / lattice_density
        trapping_term += trapping_contribution
        
        # Check for saturation (still useful for warnings)
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
            'K_equilibrium': K_eq,
            'trapping_contribution': trapping_contribution  # NEW: N_T × K / N_L
        })
    
    # Calculate effective diffusivity using CORRECT Oriani formula
    # D_eff = D_lattice / (1 + Σ(N_T,i × K_i / N_L))
    # NOT the old wrong formula: D_eff = D_lattice / (1 + Σθᵢ)
    denominator = 1.0 + trapping_term
    D_eff = D_lattice / denominator
    
    # Calculate derived quantities
    reduction_factor = D_eff / D_lattice
    mobile_fraction = 1.0 / denominator
    
    # Find dominant trap (now based on trapping_contribution, not theta)
    if trap_contributions:
        dominant_trap = max(trap_contributions, key=lambda x: x.get('trapping_contribution', 0))['name']
    else:
        dominant_trap = None
    
    # Issue warnings for extreme trapping
    if trapping_term > 100:
        warnings.warn(f"Total trapping term Σ(N_T×K/N_L) = {trapping_term:.1f} > 100. "
                     "Very strong trapping - consider using McNabb-Foster kinetic model.")
    
    if reduction_factor < 0.01:
        warnings.warn(f"Reduction factor {reduction_factor:.2e} < 0.01. "
                     "Diffusion essentially stopped - model may be inappropriate.")
    
    # Sort trap contributions by trapping_contribution (highest first)
    trap_contributions.sort(key=lambda x: x.get('trapping_contribution', 0), reverse=True)
    
    return {
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'theta_total': theta_total,
        'trapping_term': trapping_term,  # NEW: The actual factor in denominator
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
    Calculate GB path properties for parallel path model.
    
    In parallel path model, this returns GB path characteristics, not total D_eff.
    
    Returns:
    --------
    dict with:
        - 'fraction': GB volume fraction
        - 'D_intrinsic': GB diffusivity
        - 'name': 'grain_boundary'
    """
    import numpy as np
    
    # Input validation
    if D_bulk <= 0:
        raise ValueError(f"Bulk diffusivity must be positive, got {D_bulk} m²/s")
    
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature} K")
    
    if grain_size <= 0:
        raise ValueError(f"Grain size must be positive, got {grain_size} m")
    
    # Calculate GB volume fraction
    # f_gb = (surface area per volume) × thickness
    # For equiaxed grains: S_v = 3/d
    f_gb = (3.0 * gb_thickness) / grain_size
    
    # Cap at reasonable maximum
    if f_gb > 0.5:
        warnings.warn(f"GB volume fraction {f_gb:.2%} > 50% - very fine grain size, "
                     "model assumptions may break down")
        f_gb = min(f_gb, 0.5)
    
    # Get GB enhancement factor
    gb_result = gb_enhancement_factor(
        temperature=temperature,
        temperature_unit='K',
        gb_type=gb_type
    )
    
    alpha = gb_result['enhancement_factor']
    D_gb = alpha * D_bulk
    
    # Classify regime based on GB importance
    if f_gb * alpha > 10:
        regime = 'gb_dominated'
    elif f_gb * alpha < 0.1:
        regime = 'bulk_dominated'
    else:
        regime = 'mixed'
    
    return {
        'name': 'grain_boundary',
        'fraction': f_gb,
        'D_intrinsic': D_gb,
        'D_gb': D_gb,
        'D_bulk': D_bulk,
        'f_gb': f_gb,
        'enhancement_alpha': alpha,
        'regime': regime,
        'D_eff': D_gb,  # For compatibility with original interface
        'warnings': gb_result.get('warnings', [])
    }


def combined_microstructure_model(D_lattice, temperature, microstructure_params,
                                 lattice_concentration, lattice_density):
    """
    PARALLEL PATH implementation: All defects as independent parallel paths.
    
    Key Difference from Original:
    ------------------------------
    Original: D_eff = [(1-f_gb)×D_L + f_gb×D_gb] / (1 + Σθᵢ)
    This:     D_eff = Σ[f_i × D_i / (1 + θ_i)]
    
    Each path has:
    - Volume fraction f_i
    - Intrinsic diffusivity D_i
    - Path-specific trapping θ_i
    
    Paths:
    ------
    1. Lattice (intact crystal)
    2. Grain boundaries
    3. Dislocations (if density > threshold)
    4. Other paths (carbides, etc.) if specified
    
    [Same function signature as original for drop-in replacement]
    """
    import numpy as np
    import warnings
    
    # Initialize warnings list
    warning_list = []
    
    # ========================================================================
    # Input Validation
    # ========================================================================
    
    required_keys = ['grain_size', 'grain_shape', 'gb_type', 'trap_list']
    for key in required_keys:
        if key not in microstructure_params:
            raise ValueError(f"microstructure_params missing required key: '{key}'")
    
    # Extract parameters
    grain_size = microstructure_params['grain_size']
    grain_shape = microstructure_params['grain_shape']
    gb_type = microstructure_params['gb_type']
    trap_list = microstructure_params['trap_list'].copy()
    
    gb_thickness = microstructure_params.get('gb_thickness', 0.5e-9)
    include_gb_trapping = microstructure_params.get('include_gb_trapping', True)
    dislocation_density = microstructure_params.get('dislocation_density', 0)
    include_dislocation_path = microstructure_params.get('include_dislocation_path', True)
    
    # Validate inputs
    if D_lattice <= 0:
        raise ValueError(f"D_lattice must be positive, got {D_lattice}")
    if temperature <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature}")
    if grain_size <= 0:
        raise ValueError(f"Grain size must be positive, got {grain_size}")
    
    # ========================================================================
    # Step 1: Define all parallel paths
    # ========================================================================
    
    paths = []
    
    # Path 1: Grain Boundary
    gb_result = calculate_gb_enhanced_diffusivity(
        D_bulk=D_lattice,
        temperature=temperature,
        grain_size=grain_size,
        gb_thickness=gb_thickness,
        gb_type=gb_type,
        model='parallel'
    )
    
    paths.append({
        'name': 'grain_boundary',
        'fraction': gb_result['f_gb'],
        'D_intrinsic': gb_result['D_gb'],
        'trap_list': []  # Will populate with GB-specific traps
    })
    
    # Path 2: Dislocation pipe (if significant)
    if include_dislocation_path and dislocation_density > 1e15:
        disloc_result = calculate_dislocation_path_properties(
            dislocation_density=dislocation_density,
            D_lattice=D_lattice,
            temperature=temperature
        )
        
        paths.append({
            'name': 'dislocation_pipe',
            'fraction': disloc_result['fraction'],
            'D_intrinsic': disloc_result['D_intrinsic'],
            'trap_list': []  # Will populate with dislocation-specific traps
        })
    
    # Path 3: Lattice (remaining volume)
    f_lattice = 1.0 - sum(p['fraction'] for p in paths)
    
    if f_lattice < 0:
        warning_msg = f"Lattice fraction became negative ({f_lattice:.2e}). Renormalizing path fractions."
        warnings.warn(warning_msg)
        warning_list.append(warning_msg)
        
        # Renormalize
        total_f = sum(p['fraction'] for p in paths)
        for p in paths:
            p['fraction'] /= total_f
        f_lattice = 0.01  # Keep small lattice contribution
    
    paths.append({
        'name': 'lattice',
        'fraction': f_lattice,
        'D_intrinsic': D_lattice,
        'trap_list': []  # Will populate with bulk traps
    })
    
    # ========================================================================
    # Step 2: Distribute traps to paths
    # ========================================================================
    
    # Add GB traps automatically if requested
    if include_gb_trapping:
        gb_density_result = grain_boundary_density(
            grain_size=grain_size,
            gb_thickness=gb_thickness,
            grain_shape=grain_shape
        )
        
        gb_trap = {
            'name': 'grain_boundaries_auto',
            'density': gb_density_result['trap_density'],
            'binding_energy': 20.0e3  # J/mol
        }
        
        # Check if GB already in trap list
        gb_already_present = any('grain' in trap.get('name', '').lower() for trap in trap_list)
        
        if not gb_already_present:
            trap_list.append(gb_trap)
    
    # Distribute traps to appropriate paths
    # Simple rule: GB traps go to GB path, others to lattice path
    for trap in trap_list:
        trap_name_lower = trap.get('name', '').lower()
        
        if 'grain' in trap_name_lower or 'gb' in trap_name_lower:
            # Assign to GB path
            for p in paths:
                if p['name'] == 'grain_boundary':
                    p['trap_list'].append(trap)
        elif 'disloc' in trap_name_lower:
            # Assign to dislocation path if it exists
            for p in paths:
                if p['name'] == 'dislocation_pipe':
                    p['trap_list'].append(trap)
                    break
            else:
                # If no dislocation path, assign to lattice
                for p in paths:
                    if p['name'] == 'lattice':
                        p['trap_list'].append(trap)
        else:
            # Default: assign to lattice path
            for p in paths:
                if p['name'] == 'lattice':
                    p['trap_list'].append(trap)
    
    # ========================================================================
    # Step 3: Calculate trapping for each path
    # ========================================================================
    
    for path in paths:
        if path['trap_list']:
            trap_result = calculate_effective_diffusivity_trapping(
                D_lattice=path['D_intrinsic'],
                temperature=temperature,
                trap_list=path['trap_list'],
                lattice_concentration=lattice_concentration,
                lattice_density=lattice_density
            )
            
            path['D_eff'] = trap_result['D_eff']
            path['theta_total'] = trap_result['theta_total']
            path['trap_contributions'] = trap_result['trap_contributions']
            path['reduction_factor'] = trap_result['reduction_factor']
        else:
            # No traps for this path
            path['D_eff'] = path['D_intrinsic']
            path['theta_total'] = 0.0
            path['trap_contributions'] = []
            path['reduction_factor'] = 1.0
    
    # ========================================================================
    # Step 4: Combine paths via rule of mixtures
    # ========================================================================
    
    D_eff = sum(p['fraction'] * p['D_eff'] for p in paths)
    
    # Overall enhancement/reduction factor
    overall_factor = D_eff / D_lattice
    
    # Find dominant path
    path_contributions = [(p['fraction'] * p['D_eff'], p['name']) for p in paths]
    dominant_path = max(path_contributions, key=lambda x: x[0])[1]
    
    # Calculate weighted average theta
    theta_total_weighted = sum(p['fraction'] * p['theta_total'] for p in paths)
    
    # Classify regime
    gb_contribution = next((p['fraction'] * p['D_eff'] for p in paths if p['name'] == 'grain_boundary'), 0)
    gb_fraction_of_flux = gb_contribution / D_eff if D_eff > 0 else 0
    
    if gb_fraction_of_flux > 0.8:
        regime = 'gb_dominated'
    elif gb_fraction_of_flux < 0.2:
        regime = 'lattice_dominated'
    else:
        regime = 'mixed_paths'
    
    # ========================================================================
    # Step 5: Compile results
    # ========================================================================
    
    # Create calculation sequence description
    calc_lines = ["PARALLEL PATH MODEL:", ""]
    for i, p in enumerate(paths, 1):
        calc_lines.append(f"{i}. {p['name']}:")
        calc_lines.append(f"   f = {p['fraction']:.4f}, D₀ = {p['D_intrinsic']:.2e} m²/s")
        calc_lines.append(f"   θ = {p['theta_total']:.3f}, D_eff = {p['D_eff']:.2e} m²/s")
        calc_lines.append(f"   Contribution: {p['fraction']*p['D_eff']:.2e} m²/s")
    calc_lines.append("")
    calc_lines.append(f"Total: D_eff = {D_eff:.2e} m²/s ({overall_factor:.2f}× lattice)")
    
    calculation_sequence = "\n".join(calc_lines)
    
    # Determine dominant effect
    if overall_factor > 2.0:
        dominant_effect = 'fast_path_enhancement'
    elif overall_factor < 0.5:
        dominant_effect = 'trapping_reduction'
    else:
        dominant_effect = 'balanced'
    
    return {
        # Primary results
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'D_gb_enhanced': D_eff,  # For compatibility
        
        # Parallel path details
        'paths': paths,
        'dominant_path': dominant_path,
        
        # Component details (for compatibility with original interface)
        'gb_enhancement': {
            'factor': gb_result['D_gb'] / D_lattice,
            'f_gb': gb_result['f_gb'],
            'D_gb': gb_result['D_gb'],
            'regime': gb_result['regime']
        },
        
        'trapping': {
            'theta_total': theta_total_weighted,
            'reduction_factor': D_eff / D_lattice,
            'dominant_trap': None,  # Multiple paths make this ambiguous
            'contributions': []  # Would need to aggregate
        },
        
        # Overall analysis
        'overall_factor': overall_factor,
        'dominant_effect': dominant_effect,
        'regime': regime,
        'model_type': 'parallel_paths',
        
        # Diagnostic information
        'calculation_sequence': calculation_sequence,
        'parameters': {
            'temperature': temperature,
            'grain_size': grain_size,
            'gb_type': gb_type,
            'num_paths': len(paths),
            'dislocation_density': dislocation_density
        },
        
        'warnings': warning_list
    }
