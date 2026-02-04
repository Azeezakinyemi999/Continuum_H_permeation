# calculations/parallel_oxide_defect_paths.py
"""
Level 3: Parallel Path Model for Oxide with Defects

This module implements the Strehlow & Savage (1974) parallel path model
for hydrogen permeation through metals with defective oxide films.

References:
- Strehlow & Savage (1974), Nuclear Technology, 22:127-137
- Zarchy & Axtmann (1979), J. Nuclear Materials, 79:110-117
- Zhang et al. (2018), Int. J. Hydrogen Energy, 43:3353-3365
"""

import numpy as np
from scipy.optimize import brentq

# Import from Level 1
from calculations.permeation_calc import calculate_simple_metal_flux

# Import from Level 2
from calculations.interface_solver import calculate_oxide_metal_system
from calculations.oxide_permeation import molecular_diffusion_flux

# Centralized regime classification utilities
from calculations.classify_regime import (
    classify_regime_level2,
    classify_regime_level3,
    classify_regime_level4_metal,
    classify_regime_level14,
    classify_regime_level24,
    classify_regime_level34,
)

def calculate_defect_path_flux(P_upstream, P_downstream, oxide_props, metal_props, defect_props):
    """
    Calculate hydrogen flux through a defect in the oxide layer..
    
    Theory:
    -------
    Based on Strehlow & Savage (1974) parallel path model for defective oxides.
    Defects create alternative permeation paths with reduced resistance compared
    to intact oxide. The defect can be:
    1. Pinhole: Direct metal exposure (no oxide barrier)
    2. Crack: Partial oxide with reduced thickness
    3. Grain boundary: Modified transport properties
    
    Reference:
    ----------
    Strehlow, R.A. and Savage, H.C., "The permeation of hydrogen isotopes through 
    structural metals at low pressures and through metals with oxide film barriers," 
    Nuclear Technology, 22:127-137 (1974). DOI: 10.13182/NT74-A31383
    
    Mathematical Derivation:
    ------------------------
    For a pinhole (complete oxide absence):
        J_defect = (D_metal * K_s_metal / L_metal) * (sqrt(P_up) - sqrt(P_down))
        This reduces to Level 1 metal-only permeation (Sieverts' law)
    
    For a crack with thin oxide:
        The oxide thickness in crack: t_crack = α * t_oxide, where α < 1
        Uses Level 2 model with modified oxide thickness
        Interface pressure solved automatically by Level 2 solver
    
    For grain boundaries:
        D_gb = β * D_oxide, where β > 1 (enhanced diffusion)
        Otherwise similar to oxide calculation with modified properties
    
    Parameters:
    -----------
    P_upstream : float
        Upstream hydrogen pressure (Pa)
    P_downstream : float
        Downstream hydrogen pressure (Pa)
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    defect_props : dict
        Defect properties:
        - 'type': 'pinhole', 'crack', or 'grain_boundary'
        - 'thickness_factor': fraction of oxide thickness (for cracks)
        - 'diffusivity_factor': D multiplication factor (for GB)
    
    Returns:
    --------
    float
        Hydrogen flux through defect path (mol/m²/s)
    """
    
    defect_type = defect_props.get('type', 'pinhole')
    
    if defect_type == 'pinhole':
        # Direct metal exposure - use Level 1 model
        # No oxide barrier at all
        result = calculate_simple_metal_flux(
            metal_props['D_metal'],
            metal_props['K_s_metal'],
            metal_props['thickness'],
            P_upstream,
            P_downstream
        )
        flux_defect = result['flux']
        
    elif defect_type == 'crack':
        # Crack has thin oxide layer
        # Modify oxide thickness: t_crack = α * t_oxide
        alpha = defect_props.get('thickness_factor', 0.1)  # Default 10% thickness
        
        # Create modified oxide properties for crack
        crack_oxide_props = oxide_props.copy()
        crack_oxide_props['thickness'] *= alpha
        
        # Use Level 2 oxide+metal model with thin oxide
        result = calculate_oxide_metal_system(
            P_upstream, 
            P_downstream, 
            crack_oxide_props, 
            metal_props
        )
        flux_defect = result['flux']
        
    elif defect_type == 'grain_boundary':
        # Enhanced diffusion through oxide grain boundaries
        # D_gb = β * D_oxide where β > 1
        beta = defect_props.get('diffusivity_factor', 10)  # Default 10x faster
        
        gb_oxide_props = oxide_props.copy()
        gb_oxide_props['D_ox'] *= beta
        
        # Use Level 2 model with enhanced oxide diffusivity
        result = calculate_oxide_metal_system(
            P_upstream,
            P_downstream,
            gb_oxide_props,
            metal_props
        )
        flux_defect = result['flux']
    
    elif defect_type == 'mixed':
        # Mixed defects: combination of pinholes, cracks, and grain boundaries
        # Calculate weighted average flux based on component area fractions
        components = defect_props.get('components', {})
        total_component_fraction = sum(components.values())
        
        if total_component_fraction == 0:
            # No components defined, default to pinhole behavior
            result = calculate_simple_metal_flux(
                metal_props['D_metal'],
                metal_props['K_s_metal'],
                metal_props['thickness'],
                P_upstream,
                P_downstream
            )
            flux_defect = result['flux']
        else:
            # Calculate flux for each component type and weight by relative fraction
            flux_defect = 0.0
            
            # Pinhole component
            if 'pinholes' in components and components['pinholes'] > 0:
                pinhole_result = calculate_simple_metal_flux(
                    metal_props['D_metal'],
                    metal_props['K_s_metal'],
                    metal_props['thickness'],
                    P_upstream,
                    P_downstream
                )
                weight = components['pinholes'] / total_component_fraction
                flux_defect += pinhole_result['flux'] * weight
            
            # Crack component
            if 'cracks' in components and components['cracks'] > 0:
                alpha = defect_props.get('thickness_factor', 0.1)
                crack_oxide_props = oxide_props.copy()
                crack_oxide_props['thickness'] *= alpha
                crack_result = calculate_oxide_metal_system(
                    P_upstream, 
                    P_downstream, 
                    crack_oxide_props, 
                    metal_props
                )
                weight = components['cracks'] / total_component_fraction
                flux_defect += crack_result['flux'] * weight
            
            # Grain boundary component
            if 'grain_boundaries' in components and components['grain_boundaries'] > 0:
                beta = defect_props.get('diffusivity_factor', 10)
                gb_oxide_props = oxide_props.copy()
                gb_oxide_props['D_ox'] *= beta
                gb_result = calculate_oxide_metal_system(
                    P_upstream,
                    P_downstream,
                    gb_oxide_props,
                    metal_props
                )
                weight = components['grain_boundaries'] / total_component_fraction
                flux_defect += gb_result['flux'] * weight
    
    else:
        raise ValueError(f"Unknown defect type: {defect_type}")
    
    return flux_defect


def calculate_parallel_path_flux(P_upstream, P_downstream, oxide_props, metal_props, 
                                 defect_params):
    """
    Calculate total hydrogen flux through oxide with defects using parallel path model..
    
    Theory:
    -------
    The Strehlow & Savage (1974) parallel path model treats defective oxide as
    parallel permeation paths with different resistances. Unlike the simplified
    area-defect model that assumes intact oxide is impermeable, this model 
    recognizes that both intact oxide and defects contribute to permeation.
    
    The electrical resistance analogy:
    For parallel resistors: 1/R_total = 1/R_1 + 1/R_2 + ... + 1/R_n
    For permeation: J_total = J_intact * A_intact + J_defect * A_defect
    where A represents area fractions.
    
    Reference:
    ----------
    1. Strehlow & Savage (1974), Nuclear Technology, 22:127-137
    2. Zarchy & Axtmann (1979), J. Nuclear Materials, 79:110-117
       - Showed even 6Å oxide affects permeation with ~1% defects
    
    Mathematical Derivation:
    ------------------------
    Total surface area: A_total = A_intact + A_defect
    Area fractions: f_intact = A_intact/A_total, f_defect = A_defect/A_total
    
    Total flux: J_total = ∫∫ j(x,y) dA over total area
    
    Assuming uniform flux in each region:
    J_total = (1/A_total) * [∫∫_intact j_intact dA + ∫∫_defect j_defect dA]
    J_total = j_intact * (A_intact/A_total) + j_defect * (A_defect/A_total)
    J_total = j_intact * f_intact + j_defect * f_defect
    
    Where:
    - j_intact: flux density through intact oxide (from Level 2)
    - j_defect: flux density through defects
    - f_intact + f_defect = 1
    
    Parameters:
    -----------
    P_upstream : float
        Upstream hydrogen pressure (Pa)
    P_downstream : float  
        Downstream hydrogen pressure (Pa)
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    defect_params : dict
        - 'area_fraction': fraction of surface with defects (0-1)
        - 'type': defect type ('pinhole', 'crack', 'grain_boundary')
        - Additional parameters for specific defect types
    
    Returns:
    --------
    dict
        'flux_total': Total hydrogen flux (mol/m²/s)
        'flux_intact_contribution': Flux through intact oxide (mol/m²/s)
        'flux_defect_contribution': Flux through defects (mol/m²/s)
        'dominant_path': Which path carries more flux
        'defect_enhancement_factor': J_total/J_perfect_oxide
    """
    
    # Extract area fractions
    f_defect = defect_params.get('area_fraction', 0.01)  # Default 1% defects
    f_intact = 1.0 - f_defect
    
    # Validate area fraction
    if not 0 <= f_defect <= 1:
        raise ValueError(f"Defect area fraction must be 0-1, got {f_defect}")
    
    # Path 1: Through intact oxide+metal (Level 2 model)
    # This path is always present and permeable
    intact_result = calculate_oxide_metal_system(
        P_upstream, 
        P_downstream, 
        oxide_props, 
        metal_props
    )
    j_intact = intact_result['flux']
    
    # Path 2: Through defects
    j_defect = calculate_defect_path_flux(
        P_upstream, 
        P_downstream,
        oxide_props, 
        metal_props,
        defect_params
    )
    
    # Calculate contributions
    flux_intact_contribution = j_intact * f_intact
    flux_defect_contribution = j_defect * f_defect
    
    # Total flux (area-weighted parallel paths)
    flux_total = flux_intact_contribution + flux_defect_contribution
    
    # Determine dominant path
    if flux_defect_contribution > flux_intact_contribution:
        dominant = 'defects'
    else:
        dominant = 'intact_oxide'
    
    # # Calculate enhancement factor vs perfect oxide
    # enhancement = flux_total / j_intact if j_intact > 0 else float('inf')
    
    # return {
    #     'flux_total': flux_total,
    #     'flux_intact_contribution': flux_intact_contribution,
    #     'flux_defect_contribution': flux_defect_contribution,
    #     'flux_intact_per_area': j_intact,
    #     'flux_defect_per_area': j_defect,
    #     'dominant_path': dominant,
    #     'defect_enhancement_factor': enhancement,
    #     'area_fraction_defect': f_defect,
    #     'P_interface_intact': intact_result.get('P_interface', None),
    #     'regime_intact': intact_result.get('regime', None)
    # }
    # Enhancement factor vs perfect oxide with defective metal
    enhancement = flux_total / j_intact if j_intact > 0 else float('inf')
    
    # Level 3,4 hierarchical regime classification
    regime_classification = classify_regime_level34(
        base_regime=intact_result.get('regime', 'unknown'),
        flux_intact_contribution=flux_intact_contribution,
        flux_defect_contribution=flux_defect_contribution,
        modification_factor=intact_result.get('modification_factor', 1.0)
    )
    
    return {
        'flux_total': flux_total,
        'flux_intact_contribution': flux_intact_contribution,
        'flux_defect_contribution': flux_defect_contribution,
        'flux_intact_per_area': j_intact,
        'flux_defect_per_area': j_defect,
        'dominant_path': dominant,
        'defect_enhancement_factor': enhancement,
        'area_fraction_defect': f_defect,
        
        # Level 4 details from intact path
        'D_eff_metal': intact_result.get('D_eff'),
        'modification_factor': intact_result.get('modification_factor'),
        'level4_converged': intact_result.get('level4_converged'),
        
        # Interface and regime info
        'P_interface_intact': intact_result.get('P_interface'),
        'regime_intact': intact_result.get('regime'),
        
        # Level 3,4 hierarchical regime classification
        'regime_classification': regime_classification,
        'regime': regime_classification['regime_hierarchy'],
        'regime_base': regime_classification['base_regime'],
        'regime_detail': regime_classification['regime_detail'],
        
        # Microstructure details
        'microstructure_details': intact_result.get('microstructure_details', {})
    }

def calculate_PRF(P_test, oxide_props, metal_props, defect_params=None, P_downstream=0):
    """
    Calculate Permeation Reduction Factor (PRF) for defective oxide barrier..
    
    Theory:
    -------
    PRF quantifies the effectiveness of an oxide barrier in reducing hydrogen
    permeation compared to bare metal. It's defined as the ratio of flux through
    bare metal to flux through oxide-covered metal.
    
    PRF > 1: Oxide reduces permeation (desired)
    PRF >> 1: Very effective barrier
    PRF ~ 1: Oxide has little effect
    
    Reference:
    ----------
    Zhang, Q. et al., "Effects of surface oxide films on hydrogen permeation 
    and susceptibility to embrittlement of X80 steel under hydrogen atmosphere,"
    Int. J. Hydrogen Energy, 43(7):3353-3365 (2018). DOI: 10.1016/j.ijhydene.2017.12.170
    - Reported PRF up to 3828 for high-temperature oxidized steel
    
    Mathematical Derivation:
    ------------------------
    PRF = J_bare_metal / J_oxide_covered
    
    Where:
    J_bare_metal = Result from Level 1 model (Sieverts' law)
    J_oxide_covered = Result from Level 2 (perfect) or Level 3 (defective)
    
    For perfect oxide (no defects):
    PRF_max represents theoretical maximum barrier effectiveness
    
    For defective oxide:
    PRF_actual < PRF_max due to defect short-circuits
    
    The ratio PRF_actual/PRF_max indicates oxide quality
    
    Parameters:
    -----------
    P_test : float
        Test pressure for PRF calculation (Pa)
        Note: PRF is pressure-dependent due to regime transitions
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    defect_params : dict or None
        Defect parameters. If None, calculates PRF for perfect oxide
    P_downstream : float
        Downstream pressure (Pa), default 0 for standard PRF test
    
    Returns:
    --------
    dict
        'PRF': Permeation reduction factor
        'PRF_perfect': PRF for defect-free oxide (theoretical max)
        'efficiency': PRF/PRF_perfect (0-1, oxide quality metric)
        'regime': Operating regime at test pressure
        'flux_bare': Bare metal flux for reference
        'flux_with_oxide': Flux with oxide (perfect or defective)
    """
    
    # Calculate bare metal flux (Level 1)
    bare_result = calculate_simple_metal_flux(
        metal_props['D_metal'],
        metal_props['K_s_metal'],
        metal_props['thickness'],
        P_test,
        P_downstream
    )
    flux_bare_metal = bare_result['flux']
    
    # Calculate flux with perfect oxide (Level 2)
    perfect_result = calculate_oxide_metal_system(
        P_test, 
        P_downstream,
        oxide_props, 
        metal_props
    )
    flux_perfect_oxide = perfect_result['flux']
    
    PRF_perfect = flux_bare_metal / flux_perfect_oxide if flux_perfect_oxide > 0 else float('inf')
    
    # Calculate flux with defective oxide if parameters provided
    if defect_params is not None:
        # Level 3 model with defects
        defect_result = calculate_parallel_path_flux(
            P_test, 
            P_downstream,
            oxide_props, 
            metal_props,
            defect_params
        )
        flux_defective = defect_result['flux_total']
        
        PRF_actual = flux_bare_metal / flux_defective if flux_defective > 0 else float('inf')
        efficiency = PRF_actual / PRF_perfect if PRF_perfect != float('inf') else 0
        
        # Determine regime based on dominant path
        if defect_result['dominant_path'] == 'defects':
            regime = 'defect_limited'
        else:
            regime = perfect_result['regime']
            
        flux_with_oxide = flux_defective
    else:
        # No defects - perfect oxide
        PRF_actual = PRF_perfect
        efficiency = 1.0
        regime = perfect_result['regime']
        flux_with_oxide = flux_perfect_oxide
    
    return {
        'PRF': PRF_actual,
        'PRF_perfect': PRF_perfect,
        'efficiency': efficiency,
        'regime': regime,
        'test_pressure': P_test,
        'flux_bare_metal': flux_bare_metal,
        'flux_with_oxide': flux_with_oxide,
        'flux_reduction_factor': flux_bare_metal / flux_with_oxide if flux_with_oxide > 0 else float('inf')
    }


# =============================================================================
# LEVEL 3 + LEVEL 4: Defective Oxide + Defective Metal
# =============================================================================

def calculate_defect_path_flux_defective_metal(P_upstream, P_downstream, oxide_props, 
                                               metal_props, defect_props, temperature,
                                               microstructure_params, lattice_density=1.06e29,
                                               method='average', n_points=10, mode='both'):
    """
    Calculate hydrogen flux through a defect path with defective metal (Level 3+4).
    
    This extends calculate_defect_path_flux() to use Level 4 defective metal
    instead of clean metal for the metal layer.
    
    Theory:
    -------
    Same as calculate_defect_path_flux(), but the metal layer now includes:
    - Grain boundary fast diffusion paths
    - Hydrogen trapping effects
    - Position-dependent effective diffusivity
    
    Parameters
    ----------
    P_upstream : float
        Upstream hydrogen pressure [Pa]
    P_downstream : float
        Downstream hydrogen pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    defect_props : dict
        Defect properties:
        - 'type': 'pinhole', 'crack', 'grain_boundary', or 'mixed'
        - 'thickness_factor': fraction of oxide thickness (for cracks)
        - 'diffusivity_factor': D multiplication factor (for GB)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Metal microstructure specification
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    
    Returns
    -------
    float
        Hydrogen flux through defect path [mol/m²/s]
    """
    from calculations.permeation_calc import calculate_defective_metal_flux
    from calculations.interface_solver import calculate_oxide_defective_metal_system
    
    defect_type = defect_props.get('type', 'pinhole')
    
    if defect_type == 'pinhole':
        # Direct metal exposure - use Level 4 metal-only model
        # No oxide barrier at all
        result = calculate_defective_metal_flux(
            D_lattice=metal_props['D_metal'],
            K_s=metal_props['K_s_metal'],
            thickness=metal_props['thickness'],
            P_up=P_upstream,
            P_down=P_downstream,
            temperature=temperature,
            microstructure_params=microstructure_params,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode
        )
        flux_defect = result['flux']
        
    elif defect_type == 'crack':
        # Crack has thin oxide layer
        alpha = defect_props.get('thickness_factor', 0.1)
        
        crack_oxide_props = oxide_props.copy()
        crack_oxide_props['thickness'] *= alpha
        
        # Use Level 2+4: thin oxide + defective metal
        result = calculate_oxide_defective_metal_system(
            P_upstream, P_downstream,
            crack_oxide_props, metal_props,
            temperature, microstructure_params,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode
        )
        flux_defect = result['flux']
        
    elif defect_type == 'grain_boundary':
        # Enhanced diffusion through oxide grain boundaries
        beta = defect_props.get('diffusivity_factor', 10)
        
        gb_oxide_props = oxide_props.copy()
        gb_oxide_props['D_ox'] *= beta
        
        # Use Level 2+4: GB-enhanced oxide + defective metal
        result = calculate_oxide_defective_metal_system(
            P_upstream, P_downstream,
            gb_oxide_props, metal_props,
            temperature, microstructure_params,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode
        )
        flux_defect = result['flux']
        
    elif defect_type == 'mixed':
        # Mixed defects: combination weighted by component fractions
        components = defect_props.get('components', {})
        total_component_fraction = sum(components.values())
        
        if total_component_fraction == 0:
            result = calculate_defective_metal_flux(
                D_lattice=metal_props['D_metal'],
                K_s=metal_props['K_s_metal'],
                thickness=metal_props['thickness'],
                P_up=P_upstream,
                P_down=P_downstream,
                temperature=temperature,
                microstructure_params=microstructure_params,
                lattice_density=lattice_density,
                method=method,
                n_points=n_points,
                mode=mode
            )
            flux_defect = result['flux']
        else:
            flux_defect = 0.0
            
            # Pinhole component
            if 'pinholes' in components and components['pinholes'] > 0:
                pinhole_result = calculate_defective_metal_flux(
                    D_lattice=metal_props['D_metal'],
                    K_s=metal_props['K_s_metal'],
                    thickness=metal_props['thickness'],
                    P_up=P_upstream,
                    P_down=P_downstream,
                    temperature=temperature,
                    microstructure_params=microstructure_params,
                    lattice_density=lattice_density,
                    method=method,
                    n_points=n_points,
                    mode=mode
                )
                weight = components['pinholes'] / total_component_fraction
                flux_defect += pinhole_result['flux'] * weight
            
            # Crack component
            if 'cracks' in components and components['cracks'] > 0:
                alpha = defect_props.get('thickness_factor', 0.1)
                crack_oxide_props = oxide_props.copy()
                crack_oxide_props['thickness'] *= alpha
                crack_result = calculate_oxide_defective_metal_system(
                    P_upstream, P_downstream,
                    crack_oxide_props, metal_props,
                    temperature, microstructure_params,
                    lattice_density=lattice_density,
                    method=method,
                    n_points=n_points,
                    mode=mode
                )
                weight = components['cracks'] / total_component_fraction
                flux_defect += crack_result['flux'] * weight
            
            # Grain boundary component
            if 'grain_boundaries' in components and components['grain_boundaries'] > 0:
                beta = defect_props.get('diffusivity_factor', 10)
                gb_oxide_props = oxide_props.copy()
                gb_oxide_props['D_ox'] *= beta
                gb_result = calculate_oxide_defective_metal_system(
                    P_upstream, P_downstream,
                    gb_oxide_props, metal_props,
                    temperature, microstructure_params,
                    lattice_density=lattice_density,
                    method=method,
                    n_points=n_points,
                    mode=mode
                )
                weight = components['grain_boundaries'] / total_component_fraction
                flux_defect += gb_result['flux'] * weight
    else:
        raise ValueError(f"Unknown defect type: {defect_type}")
    
    return flux_defect


def calculate_parallel_path_flux_defective_metal(P_upstream, P_downstream, oxide_props, 
                                                  metal_props, defect_params, temperature,
                                                  microstructure_params, lattice_density=1.06e29,
                                                  method='average', n_points=10,
                                                  max_iterations=10, tolerance=1e-6, mode='both'):
    """
    Calculate total flux through defective oxide + defective metal (Level 3+4).
    
    This is the full Level 3+4 model combining:
    - Strehlow & Savage parallel path model for defective oxide
    - Level 4 microstructure effects in metal
    
    Theory:
    -------
    Total flux = Intact path contribution + Defect path contribution
    
    J_total = j_intact × f_intact + j_defect × f_defect
    
    Where both j_intact and j_defect now use Level 4 defective metal.
    
    Parameters
    ----------
    P_upstream : float
        Upstream hydrogen pressure [Pa]
    P_downstream : float
        Downstream hydrogen pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    defect_params : dict
        Defect parameters:
        - 'area_fraction': fraction of surface with defects (0-1)
        - 'type': defect type
        - Additional parameters for specific defect types
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Metal microstructure specification
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    max_iterations : int, optional
        Max iterations for interface solver
    tolerance : float, optional
        Convergence tolerance
    
    Returns
    -------
    dict
        Comprehensive results:
        - 'flux_total': Total hydrogen flux [mol/m²/s]
        - 'flux_intact_contribution': Flux through intact oxide [mol/m²/s]
        - 'flux_defect_contribution': Flux through defects [mol/m²/s]
        - 'dominant_path': Which path carries more flux
        - 'defect_enhancement_factor': J_total/J_perfect_oxide
        - 'level4_details': Microstructure effect details
    """
    from calculations.interface_solver import calculate_oxide_defective_metal_system
    
    # Extract area fractions
    f_defect = defect_params.get('area_fraction', 0.01)
    f_intact = 1.0 - f_defect
    
    if not 0 <= f_defect <= 1:
        raise ValueError(f"Defect area fraction must be 0-1, got {f_defect}")
    
    # Path 1: Through intact oxide + defective metal (Level 2+4)
    intact_result = calculate_oxide_defective_metal_system(
        P_upstream, P_downstream,
        oxide_props, metal_props,
        temperature, microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        max_iterations=max_iterations,
        tolerance=tolerance,
        mode=mode
    )
    j_intact = intact_result['flux']
    
    # Path 2: Through defects + defective metal (Level 3+4)
    j_defect = calculate_defect_path_flux_defective_metal(
        P_upstream, P_downstream,
        oxide_props, metal_props,
        defect_params, temperature,
        microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode
    )
    
    # Calculate contributions (area-weighted)
    flux_intact_contribution = j_intact * f_intact
    flux_defect_contribution = j_defect * f_defect
    
    # Total flux
    flux_total = flux_intact_contribution + flux_defect_contribution
    
    # Determine dominant path
    if flux_defect_contribution > flux_intact_contribution:
        dominant = 'defects'
    else:
        dominant = 'intact_oxide'
    
    # # Enhancement factor vs perfect oxide with defective metal
    # enhancement = flux_total / j_intact if j_intact > 0 else float('inf')
    
    # return {
    #     'flux_total': flux_total,
    #     'flux_intact_contribution': flux_intact_contribution,
    #     'flux_defect_contribution': flux_defect_contribution,
    #     'flux_intact_per_area': j_intact,
    #     'flux_defect_per_area': j_defect,
    #     'dominant_path': dominant,
    #     'defect_enhancement_factor': enhancement,
    #     'area_fraction_defect': f_defect,
        
    #     # Level 4 details from intact path
    #     'D_eff_metal': intact_result.get('D_eff'),
    #     'modification_factor': intact_result.get('modification_factor'),
    #     'level4_converged': intact_result.get('level4_converged'),
        
    #     # Interface and regime info
    #     'P_interface_intact': intact_result.get('P_interface'),
    #     'regime_intact': intact_result.get('regime'),
        
    #     # Microstructure details
    #     'microstructure_details': intact_result.get('microstructure_details', {})
    # }
    # Enhancement factor vs perfect oxide with defective metal
    enhancement = flux_total / j_intact if j_intact > 0 else float('inf')
    
    # Level 3,4 hierarchical regime classification
    regime_classification = classify_regime_level34(
        base_regime=intact_result.get('regime', 'unknown'),
        flux_intact_contribution=flux_intact_contribution,
        flux_defect_contribution=flux_defect_contribution,
        modification_factor=intact_result.get('modification_factor', 1.0)
    )
    
    return {
        'flux_total': flux_total,
        'flux_intact_contribution': flux_intact_contribution,
        'flux_defect_contribution': flux_defect_contribution,
        'flux_intact_per_area': j_intact,
        'flux_defect_per_area': j_defect,
        'dominant_path': dominant,
        'defect_enhancement_factor': enhancement,
        'area_fraction_defect': f_defect,
        
        # Level 4 details from intact path
        'D_eff_metal': intact_result.get('D_eff'),
        'modification_factor': intact_result.get('modification_factor'),
        'level4_converged': intact_result.get('level4_converged'),
        
        # Interface and regime info
        'P_interface_intact': intact_result.get('P_interface'),
        'regime_intact': intact_result.get('regime'),
        
        # Level 3,4 hierarchical regime classification
        'regime_classification': regime_classification,
        'regime': regime_classification['regime_hierarchy'],
        'regime_base': regime_classification['base_regime'],
        'regime_detail': regime_classification['regime_detail'],
        
        # Microstructure details
        'microstructure_details': intact_result.get('microstructure_details', {})
    }

def calculate_PRF_defective_metal(P_test, oxide_props, metal_props, temperature,
                                   microstructure_params, defect_params=None,
                                   P_downstream=0, lattice_density=1.06e29,
                                   method='average', n_points=10, mode='both'):
    """
    Calculate Permeation Reduction Factor with Level 4 defective metal.
    
    This extends calculate_PRF() to use Level 4 metal model.
    
    Parameters
    ----------
    P_test : float
        Test pressure [Pa]
    oxide_props : dict
        Oxide properties
    metal_props : dict
        Metal properties (D_metal as D_lattice)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Metal microstructure specification
    defect_params : dict or None
        Oxide defect parameters. If None, calculates PRF for perfect oxide
    P_downstream : float
        Downstream pressure [Pa]
    lattice_density : float
        Lattice site density [m⁻³]
    method : str
        D_eff averaging method
    n_points : int
        Points for integration
    
    Returns
    -------
    dict
        PRF results:
        - 'PRF': Permeation reduction factor
        - 'PRF_perfect': PRF for defect-free oxide
        - 'efficiency': PRF/PRF_perfect
        - 'regime': Operating regime
        - 'flux_bare': Bare defective metal flux
        - 'flux_with_oxide': Flux with oxide
        - 'level4_details': Microstructure effect info
    """
    from calculations.permeation_calc import calculate_defective_metal_flux
    from calculations.interface_solver import calculate_oxide_defective_metal_system
    
    # Calculate bare defective metal flux (Level 4 only)
    bare_result = calculate_defective_metal_flux(
        D_lattice=metal_props['D_metal'],
        K_s=metal_props['K_s_metal'],
        thickness=metal_props['thickness'],
        P_up=P_test,
        P_down=P_downstream,
        temperature=temperature,
        microstructure_params=microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode
    )
    flux_bare_metal = bare_result['flux']
    
    # Calculate flux with perfect oxide + defective metal (Level 2+4)
    perfect_result = calculate_oxide_defective_metal_system(
        P_test, P_downstream,
        oxide_props, metal_props,
        temperature, microstructure_params,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode
    )
    flux_perfect_oxide = perfect_result['flux']
    
    PRF_perfect = flux_bare_metal / flux_perfect_oxide if flux_perfect_oxide > 0 else float('inf')
    
    # Calculate flux with defective oxide if parameters provided
    if defect_params is not None:
        # Level 3+4: defective oxide + defective metal
        defect_result = calculate_parallel_path_flux_defective_metal(
            P_test, P_downstream,
            oxide_props, metal_props,
            defect_params, temperature,
            microstructure_params,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode
        )
        flux_defective = defect_result['flux_total']
        
        PRF_actual = flux_bare_metal / flux_defective if flux_defective > 0 else float('inf')
        efficiency = PRF_actual / PRF_perfect if PRF_perfect != float('inf') else 0
        
        if defect_result['dominant_path'] == 'defects':
            regime = 'defect_limited'
        else:
            regime = perfect_result['regime']
            
        flux_with_oxide = flux_defective
    else:
        PRF_actual = PRF_perfect
        efficiency = 1.0
        regime = perfect_result['regime']
        flux_with_oxide = flux_perfect_oxide
    
    return {
        'PRF': PRF_actual,
        'PRF_perfect': PRF_perfect,
        'efficiency': efficiency,
        'regime': regime,
        'test_pressure': P_test,
        'flux_bare_metal': flux_bare_metal,
        'flux_with_oxide': flux_with_oxide,
        'flux_reduction_factor': flux_bare_metal / flux_with_oxide if flux_with_oxide > 0 else float('inf'),
        
        # Level 4 details
        'D_eff_bare_metal': bare_result['D_eff'],
        'modification_factor_bare': bare_result['modification_factor'],
        'D_eff_with_oxide': perfect_result.get('D_eff'),
        'modification_factor_oxide': perfect_result.get('modification_factor'),
        
        # Microstructure info
        'microstructure_details': bare_result['microstructure_details']
    }



####################################################################################
####################################################################################
####################################################################################

# =============================================================================
# LEVEL 3+6: Defective Oxide + Surface Kinetics
# =============================================================================

def calculate_defect_path_flux_with_surface(P_upstream, P_downstream, oxide_props, metal_props,
                                             defect_props, temperature,
                                             k_diss=None, k_recomb=None, material_name=None,
                                             coverage_mode='equilibrium', forced_coverage=None):
    """
    Calculate hydrogen flux through a defect in the oxide layer WITH surface kinetics (L3+L6).
    
    This extends calculate_defect_path_flux() to include Level 6 surface kinetics.
    
    Theory:
    -------
    Same defect types as Level 3, but now surface dissociation limits flux:
    
    1. Pinhole: Direct metal exposure → L1+L6 (surface-limited bare metal)
    2. Crack: Thin oxide → L2+L6 (surface-limited at oxide-metal interface)
    3. Grain boundary: Enhanced oxide diffusion → L2+L6 with modified oxide
    
    The surface kinetics apply at ANY gas-metal interface:
    - For pinholes: H₂ dissociates at exposed metal surface
    - For cracks/GBs: H₂ arrives as molecular species in oxide, then dissociates
      at the oxide-metal interface
    
    Parameters
    ----------
    P_upstream : float
        Upstream hydrogen pressure [Pa]
    P_downstream : float
        Downstream hydrogen pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    defect_props : dict
        Defect properties:
        - 'type': 'pinhole', 'crack', 'grain_boundary', or 'mixed'
        - 'thickness_factor': fraction of oxide thickness (for cracks)
        - 'diffusivity_factor': D multiplication factor (for GB)
    temperature : float
        Temperature [K]
    k_diss : float, optional
        Dissociation rate constant [mol/m²/s/Pa]
    k_recomb : float, optional
        Recombination rate constant [m⁴/mol/s]
    material_name : str, optional
        Material name for kinetics lookup
    coverage_mode : str
        Mode for surface coverage:
        - 'equilibrium': Langmuir isotherm
        - 'steady_state': Solve J_diss = J_bulk
        - 'forced': Use forced_coverage
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    dict
        Contains flux and surface kinetics info (theta, Da, SRF)
    """
    from calculations.permeation_calc import calculate_simple_metal_flux_with_surface
    from calculations.interface_solver import calculate_oxide_metal_system_with_surface
    
    defect_type = defect_props.get('type', 'pinhole')
    
    if defect_type == 'pinhole':
        # Direct metal exposure - use L1+L6 (bare metal with surface kinetics)
        result = calculate_simple_metal_flux_with_surface(
            D=metal_props['D_metal'],
            K_s=metal_props['K_s_metal'],
            thickness=metal_props['thickness'],
            P_up=P_upstream,
            P_down=P_downstream,
            temperature=temperature,
            k_diss=k_diss,
            k_recomb=k_recomb,
            material_name=material_name,
            coverage_mode=coverage_mode,
            forced_coverage=forced_coverage
        )
        return {
            'flux': result['flux'],
            'flux_sieverts': result.get('flux_sieverts', result['flux']),
            'theta': result.get('theta', 0),
            'Da': result.get('Da', float('inf')),
            'SRF': result.get('SRF', 1.0),
            'coverage_mode': coverage_mode,
            'defect_type': 'pinhole'
        }
        
    elif defect_type == 'crack':
        # Crack has thin oxide layer - use L2+L6 with modified thickness
        alpha = defect_props.get('thickness_factor', 0.1)
        
        crack_oxide_props = oxide_props.copy()
        crack_oxide_props['thickness'] *= alpha
        
        result = calculate_oxide_metal_system_with_surface(
            P_upstream=P_upstream,
            P_downstream=P_downstream,
            oxide_props=crack_oxide_props,
            metal_props=metal_props,
            temperature=temperature,
            k_diss=k_diss,
            k_recomb=k_recomb,
            material_name=material_name,
            coverage_mode=coverage_mode,
            forced_coverage=forced_coverage
        )
        return {
            'flux': result['flux'],
            'flux_sieverts': result.get('flux_sieverts', result['flux']),
            'theta': result.get('theta', 0),
            'Da': result.get('Da', float('inf')),
            'SRF': result.get('SRF', 1.0),
            'coverage_mode': coverage_mode,
            'defect_type': 'crack',
            'P_interface': result.get('P_interface')
        }
        
    elif defect_type == 'grain_boundary':
        # Enhanced diffusion through oxide grain boundaries - L2+L6 with modified D_ox
        beta = defect_props.get('diffusivity_factor', 10)
        
        gb_oxide_props = oxide_props.copy()
        gb_oxide_props['D_ox'] *= beta
        
        result = calculate_oxide_metal_system_with_surface(
            P_upstream=P_upstream,
            P_downstream=P_downstream,
            oxide_props=gb_oxide_props,
            metal_props=metal_props,
            temperature=temperature,
            k_diss=k_diss,
            k_recomb=k_recomb,
            material_name=material_name,
            coverage_mode=coverage_mode,
            forced_coverage=forced_coverage
        )
        return {
            'flux': result['flux'],
            'flux_sieverts': result.get('flux_sieverts', result['flux']),
            'theta': result.get('theta', 0),
            'Da': result.get('Da', float('inf')),
            'SRF': result.get('SRF', 1.0),
            'coverage_mode': coverage_mode,
            'defect_type': 'grain_boundary',
            'P_interface': result.get('P_interface')
        }
        
    elif defect_type == 'mixed':
        # Mixed defects: weighted average by component fractions
        components = defect_props.get('components', {})
        total_component_fraction = sum(components.values())
        
        if total_component_fraction == 0:
            # Default to pinhole behavior
            result = calculate_simple_metal_flux_with_surface(
                D=metal_props['D_metal'],
                K_s=metal_props['K_s_metal'],
                thickness=metal_props['thickness'],
                P_up=P_upstream,
                P_down=P_downstream,
                temperature=temperature,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            return {
                'flux': result['flux'],
                'flux_sieverts': result.get('flux_sieverts', result['flux']),
                'theta': result.get('theta', 0),
                'Da': result.get('Da', float('inf')),
                'SRF': result.get('SRF', 1.0),
                'coverage_mode': coverage_mode,
                'defect_type': 'mixed'
            }
        
        # Calculate flux-weighted averages
        flux_defect = 0.0
        theta_avg = 0.0
        SRF_avg = 0.0
        Da_min = float('inf')  # Track minimum Da (most surface-limited)
        
        # Pinhole component
        if 'pinholes' in components and components['pinholes'] > 0:
            pinhole_result = calculate_simple_metal_flux_with_surface(
                D=metal_props['D_metal'],
                K_s=metal_props['K_s_metal'],
                thickness=metal_props['thickness'],
                P_up=P_upstream,
                P_down=P_downstream,
                temperature=temperature,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            weight = components['pinholes'] / total_component_fraction
            flux_defect += pinhole_result['flux'] * weight
            theta_avg += pinhole_result.get('theta', 0) * weight
            SRF_avg += pinhole_result.get('SRF', 1.0) * weight
            Da_min = min(Da_min, pinhole_result.get('Da', float('inf')))
        
        # Crack component
        if 'cracks' in components and components['cracks'] > 0:
            alpha = defect_props.get('thickness_factor', 0.1)
            crack_oxide_props = oxide_props.copy()
            crack_oxide_props['thickness'] *= alpha
            
            crack_result = calculate_oxide_metal_system_with_surface(
                P_upstream=P_upstream,
                P_downstream=P_downstream,
                oxide_props=crack_oxide_props,
                metal_props=metal_props,
                temperature=temperature,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            weight = components['cracks'] / total_component_fraction
            flux_defect += crack_result['flux'] * weight
            theta_avg += crack_result.get('theta', 0) * weight
            SRF_avg += crack_result.get('SRF', 1.0) * weight
            Da_min = min(Da_min, crack_result.get('Da', float('inf')))
        
        # Grain boundary component
        if 'grain_boundaries' in components and components['grain_boundaries'] > 0:
            beta = defect_props.get('diffusivity_factor', 10)
            gb_oxide_props = oxide_props.copy()
            gb_oxide_props['D_ox'] *= beta
            
            gb_result = calculate_oxide_metal_system_with_surface(
                P_upstream=P_upstream,
                P_downstream=P_downstream,
                oxide_props=gb_oxide_props,
                metal_props=metal_props,
                temperature=temperature,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            weight = components['grain_boundaries'] / total_component_fraction
            flux_defect += gb_result['flux'] * weight
            theta_avg += gb_result.get('theta', 0) * weight
            SRF_avg += gb_result.get('SRF', 1.0) * weight
            Da_min = min(Da_min, gb_result.get('Da', float('inf')))
        
        return {
            'flux': flux_defect,
            'theta': theta_avg,
            'Da': Da_min,
            'SRF': SRF_avg,
            'coverage_mode': coverage_mode,
            'defect_type': 'mixed'
        }
    
    else:
        raise ValueError(f"Unknown defect type: {defect_type}")


def calculate_parallel_path_flux_with_surface_L3L6(P_upstream, P_downstream, oxide_props, 
                                                    metal_props, defect_params, temperature,
                                                    k_diss=None, k_recomb=None, 
                                                    material_name=None,
                                                    coverage_mode='equilibrium',
                                                    forced_coverage=None):
    """
    Calculate total flux through defective oxide WITH surface kinetics (L3+L6).
    
    This is the parallel path model (Strehlow & Savage 1974) combined with
    Level 6 surface kinetics at all gas-metal interfaces.
    
    Theory:
    -------
    Total flux = Intact path contribution + Defect path contribution
    
    J_total = j_intact × f_intact + j_defect × f_defect
    
    Where both paths now include surface kinetics:
    - j_intact: L2+L6 flux through intact oxide (surface kinetics at oxide-metal interface)
    - j_defect: L1+L6 or L2+L6 flux through defects (surface kinetics at exposed/interface)
    
    The surface coverage (θ) and Damköhler number (Da) are calculated at
    each interface. The effective SRF represents the overall surface limitation.
    
    Parameters
    ----------
    P_upstream : float
        Upstream hydrogen pressure [Pa]
    P_downstream : float
        Downstream hydrogen pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    defect_params : dict
        Defect parameters:
        - 'area_fraction': fraction of surface with defects (0-1)
        - 'type': defect type ('pinhole', 'crack', 'grain_boundary', 'mixed')
        - Additional parameters for specific defect types
    temperature : float
        Temperature [K]
    k_diss : float, optional
        Dissociation rate constant [mol/m²/s/Pa]
    k_recomb : float, optional
        Recombination rate constant [m⁴/mol/s]
    material_name : str, optional
        Material name for kinetics lookup
    coverage_mode : str
        Mode for surface coverage:
        - 'equilibrium': Langmuir isotherm (default)
        - 'steady_state': Solve J_diss = J_bulk
        - 'forced': Use forced_coverage
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    dict
        Comprehensive L3+L6 results:
        - 'flux_total': Total hydrogen flux [mol/m²/s]
        - 'flux_intact_contribution': Flux through intact oxide [mol/m²/s]
        - 'flux_defect_contribution': Flux through defects [mol/m²/s]
        - 'dominant_path': Which path carries more flux
        - 'defect_enhancement_factor': J_total/J_perfect_oxide
        - Surface kinetics outputs (theta, Da, SRF for each path)
        - 'regime': Combined regime classification
    """
    from calculations.interface_solver import calculate_oxide_metal_system_with_surface
    from calculations.surface_kinetics import calculate_damkohler_number
    
    # Extract area fractions
    f_defect = defect_params.get('area_fraction', 0.01)
    f_intact = 1.0 - f_defect
    
    if not 0 <= f_defect <= 1:
        raise ValueError(f"Defect area fraction must be 0-1, got {f_defect}")
    
    # Path 1: Through intact oxide+metal with surface kinetics (L2+L6)
    intact_result = calculate_oxide_metal_system_with_surface(
        P_upstream=P_upstream,
        P_downstream=P_downstream,
        oxide_props=oxide_props,
        metal_props=metal_props,
        temperature=temperature,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    j_intact = intact_result['flux']
    
    # Path 2: Through defects with surface kinetics
    defect_result = calculate_defect_path_flux_with_surface(
        P_upstream=P_upstream,
        P_downstream=P_downstream,
        oxide_props=oxide_props,
        metal_props=metal_props,
        defect_props=defect_params,
        temperature=temperature,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    j_defect = defect_result['flux']
    
    # Calculate contributions (area-weighted)
    flux_intact_contribution = j_intact * f_intact
    flux_defect_contribution = j_defect * f_defect
    
    # Total flux
    flux_total = flux_intact_contribution + flux_defect_contribution
    
    # Determine dominant path
    if flux_defect_contribution > flux_intact_contribution:
        dominant = 'defects'
    else:
        dominant = 'intact_oxide'
    
    # Enhancement factor vs perfect oxide with surface kinetics
    enhancement = flux_total / j_intact if j_intact > 0 else float('inf')
    
    # Calculate effective surface parameters (flux-weighted)
    theta_intact = intact_result.get('theta', 0)
    theta_defect = defect_result.get('theta', 0)
    
    Da_intact = intact_result.get('Da', float('inf'))
    Da_defect = defect_result.get('Da', float('inf'))
    
    SRF_intact = intact_result.get('SRF', 1.0)
    SRF_defect = defect_result.get('SRF', 1.0)
    
    # Effective values (flux-weighted average)
    if flux_total > 0:
        theta_eff = (theta_intact * flux_intact_contribution + 
                     theta_defect * flux_defect_contribution) / flux_total
        SRF_eff = (SRF_intact * flux_intact_contribution + 
                   SRF_defect * flux_defect_contribution) / flux_total
    else:
        theta_eff = theta_intact
        SRF_eff = SRF_intact
    
    # Use minimum Da (most surface-limited path)
    Da_eff = min(Da_intact, Da_defect)
    
    # Regime classification with surface effects
    base_regime = intact_result.get('regime', 'transition')
    
    if Da_eff < 0.1 and SRF_eff < 0.5:
        regime = "surface-limited"
    elif dominant == 'defects':
        regime = f"defect-limited ({defect_params.get('type', 'pinhole')})"
    else:
        regime = base_regime
    
    return {
        # Primary flux outputs
        'flux_total': flux_total,
        'flux_intact_contribution': flux_intact_contribution,
        'flux_defect_contribution': flux_defect_contribution,
        'flux_intact_per_area': j_intact,
        'flux_defect_per_area': j_defect,
        
        # Path analysis
        'dominant_path': dominant,
        'defect_enhancement_factor': enhancement,
        'area_fraction_defect': f_defect,
        
        # Surface kinetics - intact path (L2+L6)
        'theta_intact': theta_intact,
        'Da_intact': Da_intact,
        'SRF_intact': SRF_intact,
        
        # Surface kinetics - defect path
        'theta_defect': theta_defect,
        'Da_defect': Da_defect,
        'SRF_defect': SRF_defect,
        
        # Effective (flux-weighted) surface parameters
        'theta': theta_eff,
        'Da': Da_eff,
        'SRF': SRF_eff,
        
        # Interface info
        'P_interface_intact': intact_result.get('P_interface'),
        'P_interface_defect': defect_result.get('P_interface'),
        
        # Regime classification
        'regime': regime,
        'regime_intact': base_regime,
        'coverage_mode': coverage_mode,
        
        # Temperature
        'temperature': temperature
    }


def calculate_PRF_with_surface(P_test, oxide_props, metal_props, temperature,
                                k_diss=None, k_recomb=None, material_name=None,
                                defect_params=None, P_downstream=0,
                                coverage_mode='equilibrium', forced_coverage=None):
    """
    Calculate Permeation Reduction Factor with surface kinetics (L3+L6 or L2+L6).
    
    This extends calculate_PRF() to include Level 6 surface kinetics.
    
    Parameters
    ----------
    P_test : float
        Test pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal, K_s_metal, thickness)
    temperature : float
        Temperature [K]
    k_diss : float, optional
        Dissociation rate constant
    k_recomb : float, optional
        Recombination rate constant
    material_name : str, optional
        Material name for kinetics lookup
    defect_params : dict or None
        Defect parameters. If None, calculates PRF for perfect oxide
    P_downstream : float
        Downstream pressure [Pa]
    coverage_mode : str
        Mode for surface coverage
    forced_coverage : float, optional
        Fixed coverage when coverage_mode='forced'
    
    Returns
    -------
    dict
        PRF results with surface kinetics info
    """
    from calculations.permeation_calc import calculate_simple_metal_flux_with_surface
    from calculations.interface_solver import calculate_oxide_metal_system_with_surface
    
    # Calculate bare metal flux with surface kinetics (L1+L6)
    bare_result = calculate_simple_metal_flux_with_surface(
        D=metal_props['D_metal'],
        K_s=metal_props['K_s_metal'],
        thickness=metal_props['thickness'],
        P_up=P_test,
        P_down=P_downstream,
        temperature=temperature,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    flux_bare_metal = bare_result['flux']
    
    # Calculate flux with perfect oxide + surface kinetics (L2+L6)
    perfect_result = calculate_oxide_metal_system_with_surface(
        P_upstream=P_test,
        P_downstream=P_downstream,
        oxide_props=oxide_props,
        metal_props=metal_props,
        temperature=temperature,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    flux_perfect_oxide = perfect_result['flux']
    
    PRF_perfect = flux_bare_metal / flux_perfect_oxide if flux_perfect_oxide > 0 else float('inf')
    
    # Calculate flux with defective oxide if parameters provided (L3+L6)
    if defect_params is not None:
        defect_result = calculate_parallel_path_flux_with_surface_L3L6(
            P_upstream=P_test,
            P_downstream=P_downstream,
            oxide_props=oxide_props,
            metal_props=metal_props,
            defect_params=defect_params,
            temperature=temperature,
            k_diss=k_diss,
            k_recomb=k_recomb,
            material_name=material_name,
            coverage_mode=coverage_mode,
            forced_coverage=forced_coverage
        )
        flux_defective = defect_result['flux_total']
        
        PRF_actual = flux_bare_metal / flux_defective if flux_defective > 0 else float('inf')
        efficiency = PRF_actual / PRF_perfect if PRF_perfect != float('inf') else 0
        
        regime = defect_result['regime']
        flux_with_oxide = flux_defective
        
        # Surface parameters from defect result
        theta = defect_result.get('theta', 0)
        Da = defect_result.get('Da', float('inf'))
        SRF = defect_result.get('SRF', 1.0)
    else:
        PRF_actual = PRF_perfect
        efficiency = 1.0
        regime = perfect_result.get('regime', 'transition')
        flux_with_oxide = flux_perfect_oxide
        
        theta = perfect_result.get('theta', 0)
        Da = perfect_result.get('Da', float('inf'))
        SRF = perfect_result.get('SRF', 1.0)
    
    return {
        'PRF': PRF_actual,
        'PRF_perfect': PRF_perfect,
        'efficiency': efficiency,
        'regime': regime,
        'test_pressure': P_test,
        'flux_bare_metal': flux_bare_metal,
        'flux_with_oxide': flux_with_oxide,
        'flux_reduction_factor': flux_bare_metal / flux_with_oxide if flux_with_oxide > 0 else float('inf'),
        
        # Surface kinetics outputs
        'theta': theta,
        'Da': Da,
        'SRF': SRF,
        'coverage_mode': coverage_mode,
        
        # Bare metal surface info
        'theta_bare': bare_result.get('theta', 0),
        'Da_bare': bare_result.get('Da', float('inf')),
        'SRF_bare': bare_result.get('SRF', 1.0)
    }

# =============================================================================
# LEVEL 3+4+6: Full System with Surface Kinetics (COUPLED)
# =============================================================================

def calculate_defect_path_flux_defective_metal_with_surface(P_upstream, P_downstream, 
                                                            oxide_props, metal_props,
                                                            defect_props, temperature,
                                                            microstructure_params,
                                                            k_diss=None, k_recomb=None,
                                                            material_name=None,
                                                            lattice_density=1.06e29,
                                                            method='average', n_points=10,
                                                            mode='both',
                                                            coverage_mode='equilibrium',
                                                            forced_coverage=None):
    """
    Calculate flux through a defect path with defective metal AND surface kinetics (L3+L4+L6).
    
    This extends calculate_defect_path_flux_defective_metal() to include properly
    coupled Level 6 surface kinetics.
    
    Theory:
    -------
    For each defect type, uses the appropriate coupled model:
    1. Pinhole: L4+L6 (defective metal with surface kinetics)
    2. Crack: L5+L6 (thin oxide + defective metal + surface kinetics)
    3. Grain boundary: L5+L6 (enhanced oxide + defective metal + surface kinetics)
    
    Parameters
    ----------
    P_upstream : float
        Upstream hydrogen pressure [Pa]
    P_downstream : float
        Downstream hydrogen pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    defect_props : dict
        Defect properties (type, thickness_factor, diffusivity_factor)
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Metal microstructure specification
    k_diss, k_recomb : float, optional
        Surface kinetics parameters
    material_name : str, optional
        Material name for kinetics lookup
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    mode : str, optional
        Microstructure mode
    coverage_mode : str
        Surface coverage mode ('equilibrium', 'steady_state', 'forced')
    forced_coverage : float, optional
        Required when coverage_mode='forced'
    
    Returns
    -------
    dict
        Contains flux, D_eff, theta, Da, SRF, and microstructure details
    """
    from calculations.permeation_calc import calculate_defective_metal_flux_with_surface
    from calculations.interface_solver import calculate_oxide_defective_metal_system_with_surface
    
    defect_type = defect_props.get('type', 'pinhole')
    
    if defect_type == 'pinhole':
        # Direct metal exposure - use L4+L6 (defective metal with surface kinetics)
        result = calculate_defective_metal_flux_with_surface(
            D_lattice=metal_props['D_metal'],
            K_s=metal_props['K_s_metal'],
            thickness=metal_props['thickness'],
            P_up=P_upstream,
            P_down=P_downstream,
            temperature=temperature,
            microstructure_params=microstructure_params,
            k_diss=k_diss,
            k_recomb=k_recomb,
            material_name=material_name,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode,
            coverage_mode=coverage_mode,
            forced_coverage=forced_coverage
        )
        return {
            'flux': result['flux'],
            'flux_sieverts': result.get('flux_sieverts', result['flux']),
            'D_eff': result.get('D_eff'),
            'modification_factor': result.get('modification_factor', 1.0),
            'theta': result.get('theta', 0),
            'Da': result.get('Da', float('inf')),
            'SRF': result.get('SRF', 1.0),
            'coverage_mode': coverage_mode,
            'defect_type': 'pinhole',
            'microstructure_details': result.get('microstructure_details', {})
        }
        
    elif defect_type == 'crack':
        # Crack with thin oxide - use L5+L6 COUPLED
        alpha = defect_props.get('thickness_factor', 0.1)
        
        crack_oxide_props = oxide_props.copy()
        crack_oxide_props['thickness'] *= alpha
        
        result = calculate_oxide_defective_metal_system_with_surface(
            P_upstream=P_upstream,
            P_downstream=P_downstream,
            oxide_props=crack_oxide_props,
            metal_props=metal_props,
            temperature=temperature,
            microstructure_params=microstructure_params,
            k_diss=k_diss,
            k_recomb=k_recomb,
            material_name=material_name,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode,
            coverage_mode=coverage_mode,
            forced_coverage=forced_coverage
        )
        return {
            'flux': result['flux'],
            'flux_sieverts': result.get('flux_sieverts', result['flux']),
            'D_eff': result.get('D_eff'),
            'modification_factor': result.get('modification_factor', 1.0),
            'theta': result.get('theta', 0),
            'Da': result.get('Da', float('inf')),
            'SRF': result.get('SRF', 1.0),
            'coverage_mode': coverage_mode,
            'defect_type': 'crack',
            'P_interface': result.get('P_interface'),
            'microstructure_details': result.get('microstructure_details', {})
        }
        
    elif defect_type == 'grain_boundary':
        # Enhanced oxide diffusion - use L5+L6 COUPLED
        beta = defect_props.get('diffusivity_factor', 10)
        
        gb_oxide_props = oxide_props.copy()
        gb_oxide_props['D_ox'] *= beta
        
        result = calculate_oxide_defective_metal_system_with_surface(
            P_upstream=P_upstream,
            P_downstream=P_downstream,
            oxide_props=gb_oxide_props,
            metal_props=metal_props,
            temperature=temperature,
            microstructure_params=microstructure_params,
            k_diss=k_diss,
            k_recomb=k_recomb,
            material_name=material_name,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode,
            coverage_mode=coverage_mode,
            forced_coverage=forced_coverage
        )
        return {
            'flux': result['flux'],
            'flux_sieverts': result.get('flux_sieverts', result['flux']),
            'D_eff': result.get('D_eff'),
            'modification_factor': result.get('modification_factor', 1.0),
            'theta': result.get('theta', 0),
            'Da': result.get('Da', float('inf')),
            'SRF': result.get('SRF', 1.0),
            'coverage_mode': coverage_mode,
            'defect_type': 'grain_boundary',
            'P_interface': result.get('P_interface'),
            'microstructure_details': result.get('microstructure_details', {})
        }
        
    elif defect_type == 'mixed':
        # Mixed defects: weighted average by component fractions
        components = defect_props.get('components', {})
        total_component_fraction = sum(components.values())
        
        if total_component_fraction == 0:
            # Default to pinhole behavior
            result = calculate_defective_metal_flux_with_surface(
                D_lattice=metal_props['D_metal'],
                K_s=metal_props['K_s_metal'],
                thickness=metal_props['thickness'],
                P_up=P_upstream,
                P_down=P_downstream,
                temperature=temperature,
                microstructure_params=microstructure_params,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                lattice_density=lattice_density,
                method=method,
                n_points=n_points,
                mode=mode,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            return {
                'flux': result['flux'],
                'flux_sieverts': result.get('flux_sieverts', result['flux']),
                'D_eff': result.get('D_eff'),
                'modification_factor': result.get('modification_factor', 1.0),
                'theta': result.get('theta', 0),
                'Da': result.get('Da', float('inf')),
                'SRF': result.get('SRF', 1.0),
                'coverage_mode': coverage_mode,
                'defect_type': 'mixed',
                'microstructure_details': result.get('microstructure_details', {})
            }
        
        # Calculate flux-weighted averages
        flux_defect = 0.0
        flux_sieverts_defect = 0.0
        theta_avg = 0.0
        SRF_avg = 0.0
        Da_min = float('inf')
        D_eff_avg = 0.0
        
        # Pinhole component
        if 'pinholes' in components and components['pinholes'] > 0:
            pinhole_result = calculate_defective_metal_flux_with_surface(
                D_lattice=metal_props['D_metal'],
                K_s=metal_props['K_s_metal'],
                thickness=metal_props['thickness'],
                P_up=P_upstream,
                P_down=P_downstream,
                temperature=temperature,
                microstructure_params=microstructure_params,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                lattice_density=lattice_density,
                method=method,
                n_points=n_points,
                mode=mode,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            weight = components['pinholes'] / total_component_fraction
            flux_defect += pinhole_result['flux'] * weight
            flux_sieverts_defect += pinhole_result.get('flux_sieverts', pinhole_result['flux']) * weight
            theta_avg += pinhole_result.get('theta', 0) * weight
            SRF_avg += pinhole_result.get('SRF', 1.0) * weight
            D_eff_avg += pinhole_result.get('D_eff', metal_props['D_metal']) * weight
            Da_min = min(Da_min, pinhole_result.get('Da', float('inf')))
        
        # Crack component
        if 'cracks' in components and components['cracks'] > 0:
            alpha = defect_props.get('thickness_factor', 0.1)
            crack_oxide_props = oxide_props.copy()
            crack_oxide_props['thickness'] *= alpha
            
            crack_result = calculate_oxide_defective_metal_system_with_surface(
                P_upstream=P_upstream,
                P_downstream=P_downstream,
                oxide_props=crack_oxide_props,
                metal_props=metal_props,
                temperature=temperature,
                microstructure_params=microstructure_params,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                lattice_density=lattice_density,
                method=method,
                n_points=n_points,
                mode=mode,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            weight = components['cracks'] / total_component_fraction
            flux_defect += crack_result['flux'] * weight
            flux_sieverts_defect += crack_result.get('flux_sieverts', crack_result['flux']) * weight
            theta_avg += crack_result.get('theta', 0) * weight
            SRF_avg += crack_result.get('SRF', 1.0) * weight
            D_eff_avg += crack_result.get('D_eff', metal_props['D_metal']) * weight
            Da_min = min(Da_min, crack_result.get('Da', float('inf')))
        
        # Grain boundary component
        if 'grain_boundaries' in components and components['grain_boundaries'] > 0:
            beta = defect_props.get('diffusivity_factor', 10)
            gb_oxide_props = oxide_props.copy()
            gb_oxide_props['D_ox'] *= beta
            
            gb_result = calculate_oxide_defective_metal_system_with_surface(
                P_upstream=P_upstream,
                P_downstream=P_downstream,
                oxide_props=gb_oxide_props,
                metal_props=metal_props,
                temperature=temperature,
                microstructure_params=microstructure_params,
                k_diss=k_diss,
                k_recomb=k_recomb,
                material_name=material_name,
                lattice_density=lattice_density,
                method=method,
                n_points=n_points,
                mode=mode,
                coverage_mode=coverage_mode,
                forced_coverage=forced_coverage
            )
            weight = components['grain_boundaries'] / total_component_fraction
            flux_defect += gb_result['flux'] * weight
            flux_sieverts_defect += gb_result.get('flux_sieverts', gb_result['flux']) * weight
            theta_avg += gb_result.get('theta', 0) * weight
            SRF_avg += gb_result.get('SRF', 1.0) * weight
            D_eff_avg += gb_result.get('D_eff', metal_props['D_metal']) * weight
            Da_min = min(Da_min, gb_result.get('Da', float('inf')))
        
        return {
            'flux': flux_defect,
            'flux_sieverts': flux_sieverts_defect,
            'D_eff': D_eff_avg,
            'modification_factor': D_eff_avg / metal_props['D_metal'] if metal_props['D_metal'] > 0 else 1.0,
            'theta': theta_avg,
            'Da': Da_min,
            'SRF': SRF_avg,
            'coverage_mode': coverage_mode,
            'defect_type': 'mixed',
            'microstructure_details': {}
        }
    
    else:
        raise ValueError(f"Unknown defect type: {defect_type}")


def calculate_parallel_path_flux_with_surface(P_upstream, P_downstream, oxide_props, metal_props,
                                               defect_params, temperature, microstructure_params,
                                               k_diss=None, k_recomb=None, material_name=None,
                                               lattice_density=1.06e29,
                                               method='average', n_points=10,
                                               mode='both',
                                               coverage_mode='equilibrium',
                                               forced_coverage=None):
    """
    Calculate flux through defective oxide + defective metal WITH surface kinetics (COUPLED).
    
    This is Level 3+4+6: The complete hierarchical model combining:
    - Surface dissociation limitation (Level 6) - PROPERLY COUPLED
    - Defective oxide with parallel paths (Level 3)
    - Defective metal with GB + trapping (Level 4)
    
    Physics:
    --------
    Total flux = Intact path contribution + Defect path contribution
    
    J_total = j_intact × f_intact + j_defect × f_defect
    
    Where:
    - j_intact: L5+L6 flux through intact oxide (COUPLED surface kinetics)
    - j_defect: L4+L6 or L5+L6 flux through defects (COUPLED surface kinetics)
    
    The coupling ensures surface kinetics is included in the interface pressure
    solver, not applied as a post-processing step.
    
    Parameters
    ----------
    P_upstream : float
        Upstream pressure [Pa]
    P_downstream : float
        Downstream pressure [Pa]
    oxide_props : dict
        Oxide properties (D_ox, K_ox, thickness)
    metal_props : dict
        Metal properties (D_metal as D_lattice, K_s_metal, thickness)
    defect_params : dict
        Defect parameters:
        - 'area_fraction': fraction of surface with defects (0-1)
        - 'type': defect type
        - Additional parameters for specific defect types
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Microstructure specification
    k_diss : float, optional
        Dissociation rate constant [mol/m²/s/Pa]
    k_recomb : float, optional
        Recombination rate constant [m⁴/mol/s]
    material_name : str, optional
        Material name for kinetics lookup
    lattice_density : float, optional
        Lattice site density [m⁻³]
    method : str, optional
        D_eff averaging method
    n_points : int, optional
        Points for integration
    mode : str, optional
        Microstructure mode ('both', 'gb_only', 'trapping_only', 'none')
    coverage_mode : str
        Surface coverage mode ('equilibrium', 'steady_state', 'forced')
    forced_coverage : float, optional
        Required when coverage_mode='forced'
        
    Returns
    -------
    dict
        Complete Level 3+4+6 results with properly coupled surface kinetics
    """
    from calculations.interface_solver import calculate_oxide_defective_metal_system_with_surface
    from calculations.classify_regime import classify_regime_level56
    
    # Extract area fractions
    f_defect = defect_params.get('area_fraction', 0.01)
    f_intact = 1.0 - f_defect
    
    if not 0 <= f_defect <= 1:
        raise ValueError(f"Defect area fraction must be 0-1, got {f_defect}")
    
    # Path 1: Through intact oxide + defective metal + surface kinetics (L5+L6 COUPLED)
    intact_result = calculate_oxide_defective_metal_system_with_surface(
        P_upstream=P_upstream,
        P_downstream=P_downstream,
        oxide_props=oxide_props,
        metal_props=metal_props,
        temperature=temperature,
        microstructure_params=microstructure_params,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    j_intact = intact_result['flux']
    
    # Path 2: Through defects + defective metal + surface kinetics (L4+L6 or L5+L6 COUPLED)
    defect_result = calculate_defect_path_flux_defective_metal_with_surface(
        P_upstream=P_upstream,
        P_downstream=P_downstream,
        oxide_props=oxide_props,
        metal_props=metal_props,
        defect_props=defect_params,
        temperature=temperature,
        microstructure_params=microstructure_params,
        k_diss=k_diss,
        k_recomb=k_recomb,
        material_name=material_name,
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode,
        coverage_mode=coverage_mode,
        forced_coverage=forced_coverage
    )
    j_defect = defect_result['flux']
    
    # Calculate contributions (area-weighted)
    flux_intact_contribution = j_intact * f_intact
    flux_defect_contribution = j_defect * f_defect
    
    # Total flux
    flux_total = flux_intact_contribution + flux_defect_contribution
    
    # Determine dominant path
    if flux_defect_contribution > flux_intact_contribution:
        dominant = 'defects'
    else:
        dominant = 'intact_oxide'
    
    # Enhancement factor vs perfect oxide with defective metal and surface kinetics
    enhancement = flux_total / j_intact if j_intact > 0 else float('inf')
    
    # Surface kinetics from each path
    theta_intact = intact_result.get('theta', 0)
    theta_defect = defect_result.get('theta', 0)
    
    Da_intact = intact_result.get('Da', float('inf'))
    Da_defect = defect_result.get('Da', float('inf'))
    
    SRF_intact = intact_result.get('SRF', 1.0)
    SRF_defect = defect_result.get('SRF', 1.0)
    
    # Effective values (flux-weighted average)
    if flux_total > 0:
        theta_eff = (theta_intact * flux_intact_contribution + 
                     theta_defect * flux_defect_contribution) / flux_total
        SRF_eff = (SRF_intact * flux_intact_contribution + 
                   SRF_defect * flux_defect_contribution) / flux_total
    else:
        theta_eff = theta_intact
        SRF_eff = SRF_intact
    
    # Use minimum Da (most surface-limited path)
    Da_eff = min(Da_intact, Da_defect)
    
    # Level 4 outputs
    D_eff_intact = intact_result.get('D_eff', metal_props['D_metal'])
    D_eff_defect = defect_result.get('D_eff', metal_props['D_metal'])
    modification_factor_intact = intact_result.get('modification_factor', 1.0)
    modification_factor_defect = defect_result.get('modification_factor', 1.0)
    
    # Base regime from intact path
    base_regime = intact_result.get('regime', 'transition')
    
    # Full regime classification
    regime_class = classify_regime_level56(
        base_regime=base_regime,
        Da=Da_eff, theta=theta_eff, SRF=SRF_eff,
        modification_factor=modification_factor_intact,
        flux_intact_contribution=flux_intact_contribution,
        flux_defect_contribution=flux_defect_contribution
    )
    
    return {
        # Primary outputs
        'flux_total': flux_total,
        'flux_intact_contribution': flux_intact_contribution,
        'flux_defect_contribution': flux_defect_contribution,
        'flux_intact_per_area': j_intact,
        'flux_defect_per_area': j_defect,
        
        # Path analysis
        'dominant_path': dominant,
        'defect_enhancement_factor': enhancement,
        'area_fraction_defect': f_defect,
        
        # Surface kinetics - intact path (L5+L6)
        'theta_intact': theta_intact,
        'Da_intact': Da_intact,
        'SRF_intact': SRF_intact,
        
        # Surface kinetics - defect path (L4+L6 or L5+L6)
        'theta_defect': theta_defect,
        'Da_defect': Da_defect,
        'SRF_defect': SRF_defect,
        
        # Effective (flux-weighted) surface parameters
        'theta': theta_eff,
        'Da': Da_eff,
        'SRF': SRF_eff,
        
        # Level 4 outputs
        'D_eff_intact': D_eff_intact,
        'D_eff_defect': D_eff_defect,
        'D_eff': D_eff_intact,  # Primary D_eff from intact path
        'modification_factor': modification_factor_intact,
        
        # Interface info
        'P_interface_intact': intact_result.get('P_interface'),
        'P_interface_defect': defect_result.get('P_interface'),
        
        # Regime classification
        'regime_classification': regime_class,
        'regime': regime_class['regime_hierarchy'],
        'regime_base': base_regime,
        'surface_significant': regime_class['surface_significant'],
        'dominant_limitation': regime_class['dominant_limitation'],
        
        # Microstructure details
        'microstructure_details': intact_result.get('microstructure_details', {}),
        
        # Coverage mode
        'coverage_mode': coverage_mode,
        'temperature': temperature,
        
        'units': {
            'flux': 'mol/m²/s',
            'diffusivity': 'm²/s'
        }
    }