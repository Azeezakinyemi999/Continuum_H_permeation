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
    
    # Calculate enhancement factor vs perfect oxide
    enhancement = flux_total / j_intact if j_intact > 0 else float('inf')
    
    return {
        'flux_total': flux_total,
        'flux_intact_contribution': flux_intact_contribution,
        'flux_defect_contribution': flux_defect_contribution,
        'flux_intact_per_area': j_intact,
        'flux_defect_per_area': j_defect,
        'dominant_path': dominant,
        'defect_enhancement_factor': enhancement,
        'area_fraction_defect': f_defect,
        'P_interface_intact': intact_result.get('P_interface', None),
        'regime_intact': intact_result.get('regime', None)
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