"""
calculations/classify_regime.py

Hierarchical regime classification utilities used by Level 1-4 models.

This module centralizes the various `classify_regime_*` helper functions so
that other calculation modules can import and reuse a single implementation.
"""

def classify_regime_level2(base_regime):
    return {
        'model_level': 'Level 2',
        'base_regime': base_regime,
        'regime_hierarchy': base_regime,
        'regime_detail': None,
        'regime_subdetail': None,
        'classification_depth': 1
    }


def classify_regime_level3(base_regime, flux_intact_contribution, flux_defect_contribution):
    if base_regime == 'oxide_limited':
        if flux_defect_contribution > flux_intact_contribution:
            regime_detail = 'defect_limited'
        else:
            regime_detail = 'regime_intact_oxide'
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
    else:
        regime_detail = None
        regime_hierarchy = base_regime
        classification_depth = 1

    return {
        'model_level': 'Level 3',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': classification_depth,
        'flux_ratio_defect_to_intact': flux_defect_contribution / flux_intact_contribution if flux_intact_contribution > 0 else float('inf')
    }


def classify_regime_level4_metal(modification_factor, threshold_traps=0.5):
    if modification_factor < threshold_traps:
        metal_regime_detail = 'traps_defect_limited'
        trapping_significant = True
    else:
        metal_regime_detail = 'lattice_limited'
        trapping_significant = False

    return {
        'metal_regime_detail': metal_regime_detail,
        'modification_factor': modification_factor,
        'trapping_significant': trapping_significant,
        'trapping_reduction_percent': (1 - modification_factor) * 100
    }


def classify_regime_level14(modification_factor, threshold_traps=0.5):
    metal_detail = classify_regime_level4_metal(modification_factor, threshold_traps)

    base_regime = 'metal_limited'
    regime_detail = metal_detail['metal_regime_detail']
    regime_hierarchy = f"{base_regime}/{regime_detail}"

    return {
        'model_level': 'Level 1,4',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': 2,
        'modification_factor': modification_factor,
        'trapping_significant': metal_detail['trapping_significant'],
        'trapping_reduction_percent': metal_detail['trapping_reduction_percent']
    }


def classify_regime_level24(base_regime, modification_factor, threshold_traps=0.5):
    if base_regime == 'metal_limited':
        metal_detail = classify_regime_level4_metal(modification_factor, threshold_traps)
        regime_detail = metal_detail['metal_regime_detail']
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
        trapping_significant = metal_detail['trapping_significant']
        trapping_reduction = metal_detail['trapping_reduction_percent']
    else:
        regime_detail = None
        regime_hierarchy = base_regime
        classification_depth = 1
        trapping_significant = False
        trapping_reduction = 0.0

    return {
        'model_level': 'Level 2,4',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': classification_depth,
        'modification_factor': modification_factor,
        'trapping_significant': trapping_significant,
        'trapping_reduction_percent': trapping_reduction
    }


def classify_regime_level34(base_regime, flux_intact_contribution, flux_defect_contribution,
                            modification_factor, threshold_traps=0.5):
    if base_regime == 'oxide_limited':
        if flux_defect_contribution > flux_intact_contribution:
            regime_detail = 'defect_limited'
        else:
            regime_detail = 'regime_intact_oxide'
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
        trapping_significant = False
        trapping_reduction = 0.0
        flux_ratio = flux_defect_contribution / flux_intact_contribution if flux_intact_contribution > 0 else float('inf')

    elif base_regime == 'metal_limited':
        metal_detail = classify_regime_level4_metal(modification_factor, threshold_traps)
        regime_detail = metal_detail['metal_regime_detail']
        regime_hierarchy = f"{base_regime}/{regime_detail}"
        classification_depth = 2
        trapping_significant = metal_detail['trapping_significant']
        trapping_reduction = metal_detail['trapping_reduction_percent']
        flux_ratio = None

    else:
        regime_detail = None
        regime_hierarchy = base_regime
        classification_depth = 1
        trapping_significant = False
        trapping_reduction = 0.0
        flux_ratio = None

    return {
        'model_level': 'Level 3,4',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': None,
        'classification_depth': classification_depth,
        'modification_factor': modification_factor,
        'trapping_significant': trapping_significant,
        'trapping_reduction_percent': trapping_reduction,
        'flux_ratio_defect_to_intact': flux_ratio
    }


# =============================================================================
# LEVEL 6: Surface Kinetics Regime Classification
# =============================================================================

def classify_regime_level6(Da, theta=None, coverage_mode='equilibrium', theta_equilibrium=None):
    """
    Classify surface kinetics regime based on Damköhler number and coverage mode.
    
    Parameters
    ----------
    Da : float
        Damköhler number for dissociation = k_diss × K_s × L / D
    theta : float, optional
        Surface coverage (0-1)
    coverage_mode : str, optional
        Mode used to calculate coverage:
        - 'equilibrium': Langmuir isotherm (adsorption-desorption equilibrium)
        - 'steady_state': Solved from J_diss = J_bulk balance
        - 'forced': User-specified coverage
    theta_equilibrium : float, optional
        Equilibrium coverage from Langmuir (for comparison with steady_state)
        
    Returns
    -------
    dict
        Surface kinetics regime classification
    """
    # Da-based regime classification
    if Da > 10:
        surface_regime = 'diffusion-limited'
        surface_significant = False
    elif Da < 0.1:
        surface_regime = 'dissociation-limited'
        surface_significant = True
    else:
        surface_regime = 'mixed-surface-control'
        surface_significant = True
    
    # Coverage regime
    if theta is not None:
        if theta > 0.9:
            coverage_regime = 'high_coverage'
        elif theta < 0.1:
            coverage_regime = 'low_coverage'
        else:
            coverage_regime = 'intermediate_coverage'
    else:
        coverage_regime = 'unknown'
    
    # Coverage mode interpretation
    coverage_mode_info = {
        'equilibrium': 'Surface at adsorption-desorption equilibrium (Langmuir)',
        'steady_state': 'Surface depleted by bulk flux (J_diss = J_bulk)',
        'forced': 'User-specified coverage (non-equilibrium)'
    }.get(coverage_mode, 'Unknown coverage mode')
    
    # Check for significant depletion (steady_state vs equilibrium)
    coverage_depletion = None
    if coverage_mode == 'steady_state' and theta_equilibrium is not None and theta is not None:
        if theta_equilibrium > 0:
            coverage_depletion = (theta_equilibrium - theta) / theta_equilibrium
    
    return {
        'model_level': 'Level 6',
        'surface_regime': surface_regime,
        'coverage_regime': coverage_regime,
        'coverage_mode': coverage_mode,
        'coverage_mode_info': coverage_mode_info,
        'coverage_depletion': coverage_depletion,
        'Da': Da,
        'theta': theta,
        'theta_equilibrium': theta_equilibrium,
        'surface_significant': surface_significant,
        'regime_hierarchy': f"L6:{surface_regime}({coverage_mode})"
    }


def classify_regime_level16(Da, theta=None, SRF=1.0, coverage_mode='equilibrium', 
                           theta_equilibrium=None):
    """
    Classify regime for Level 1 + Level 6 (Perfect Metal with Surface Kinetics).
    
    Parameters
    ----------
    Da : float
        Damköhler number
    theta : float, optional
        Surface coverage used in calculation
    SRF : float
        Surface Reduction Factor (flux/flux_sieverts)
    coverage_mode : str
        'equilibrium', 'steady_state', or 'forced'
    theta_equilibrium : float, optional
        Equilibrium coverage (Langmuir) for comparison
    """
    surface_class = classify_regime_level6(Da, theta, coverage_mode, theta_equilibrium)
    
    base_regime = 'metal_limited'
    if surface_class['surface_significant']:
        regime_detail = surface_class['surface_regime']
        regime_hierarchy = f"{base_regime}/{regime_detail}"
    else:
        regime_detail = 'sieverts_equilibrium'
        regime_hierarchy = f"{base_regime}/{regime_detail}"
    
    # Add coverage mode context to regime
    if coverage_mode != 'equilibrium':
        regime_hierarchy = f"{regime_hierarchy}({coverage_mode})"
    
    return {
        'model_level': 'Level 1+6',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'surface_regime': surface_class['surface_regime'],
        'coverage_regime': surface_class['coverage_regime'],
        'coverage_mode': coverage_mode,
        'coverage_mode_info': surface_class['coverage_mode_info'],
        'coverage_depletion': surface_class['coverage_depletion'],
        'Da': Da,
        'theta': theta,
        'theta_equilibrium': theta_equilibrium,
        'SRF': SRF,
        'surface_significant': surface_class['surface_significant'],
        'classification_depth': 2
    }


def classify_regime_level46(Da, theta, SRF, modification_factor, threshold_traps=0.5,
                           coverage_mode='equilibrium', theta_equilibrium=None):
    """
    Classify regime for Level 4 + Level 6 (Defective Metal with Surface Kinetics).
    
    Parameters
    ----------
    Da : float
        Damköhler number (using D_eff from L4)
    theta : float
        Surface coverage used in calculation
    SRF : float
        Surface Reduction Factor
    modification_factor : float
        D_eff/D_lattice from trapping effects
    threshold_traps : float
        Threshold for significant trapping
    coverage_mode : str
        'equilibrium', 'steady_state', or 'forced'
    theta_equilibrium : float, optional
        Equilibrium coverage for comparison
    """
    surface_class = classify_regime_level6(Da, theta, coverage_mode, theta_equilibrium)
    metal_class = classify_regime_level4_metal(modification_factor, threshold_traps)
    
    base_regime = 'metal_limited'
    
    if surface_class['surface_significant'] and SRF < 0.5:
        regime_detail = surface_class['surface_regime']
        regime_subdetail = metal_class['metal_regime_detail']
        dominant = 'surface'
    elif metal_class['trapping_significant']:
        regime_detail = metal_class['metal_regime_detail']
        regime_subdetail = surface_class['surface_regime']
        dominant = 'trapping'
    else:
        regime_detail = 'combined_effects'
        regime_subdetail = f"surface:{surface_class['surface_regime']}/metal:{metal_class['metal_regime_detail']}"
        dominant = 'balanced'
    
    regime_hierarchy = f"{base_regime}/{regime_detail}"
    if coverage_mode != 'equilibrium':
        regime_hierarchy = f"{regime_hierarchy}({coverage_mode})"
    
    return {
        'model_level': 'Level 4+6',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': regime_subdetail,
        'dominant_limitation': dominant,
        'surface_regime': surface_class['surface_regime'],
        'metal_regime': metal_class['metal_regime_detail'],
        'coverage_mode': coverage_mode,
        'coverage_mode_info': surface_class['coverage_mode_info'],
        'coverage_depletion': surface_class['coverage_depletion'],
        'Da': Da,
        'theta': theta,
        'theta_equilibrium': theta_equilibrium,
        'SRF': SRF,
        'modification_factor': modification_factor,
        'surface_significant': surface_class['surface_significant'],
        'trapping_significant': metal_class['trapping_significant'],
        'classification_depth': 3
    }


def classify_regime_level56(base_regime, Da, theta, SRF, modification_factor,
                            flux_intact_contribution=None, flux_defect_contribution=None,
                            threshold_traps=0.5, coverage_mode='equilibrium', 
                            theta_equilibrium=None):
    """
    Classify regime for Level 5 + Level 6 (Full System with Surface Kinetics).
    
    Parameters
    ----------
    base_regime : str
        'oxide_limited' or 'metal_limited'
    Da : float
        Damköhler number
    theta : float
        Surface coverage
    SRF : float
        Surface Reduction Factor
    modification_factor : float
        D_eff/D_lattice from trapping
    flux_intact_contribution : float, optional
        Flux through intact oxide
    flux_defect_contribution : float, optional
        Flux through oxide defects
    threshold_traps : float
        Threshold for significant trapping
    coverage_mode : str
        'equilibrium', 'steady_state', or 'forced'
    theta_equilibrium : float, optional
        Equilibrium coverage for comparison
    """
    surface_class = classify_regime_level6(Da, theta, coverage_mode, theta_equilibrium)
    surface_is_primary = surface_class['surface_significant'] and SRF < 0.5
    
    if surface_is_primary:
        regime_detail = surface_class['surface_regime']
        regime_subdetail = f"base:{base_regime}"
        dominant = 'surface'
        
    elif base_regime == 'oxide_limited':
        if flux_defect_contribution is not None and flux_intact_contribution is not None:
            if flux_defect_contribution > flux_intact_contribution:
                regime_detail = 'defect_limited'
            else:
                regime_detail = 'intact_oxide'
        else:
            regime_detail = 'oxide_barrier'
        regime_subdetail = surface_class['surface_regime']
        dominant = 'oxide'
        
    else:  # metal_limited
        metal_class = classify_regime_level4_metal(modification_factor, threshold_traps)
        regime_detail = metal_class['metal_regime_detail']
        regime_subdetail = surface_class['surface_regime']
        dominant = 'metal' if not metal_class['trapping_significant'] else 'trapping'
    
    regime_hierarchy = f"{base_regime}/{regime_detail}"
    if surface_class['surface_significant']:
        regime_hierarchy = f"L6:{surface_class['surface_regime']}({coverage_mode})|{regime_hierarchy}"
    elif coverage_mode != 'equilibrium':
        regime_hierarchy = f"{regime_hierarchy}({coverage_mode})"
    
    return {
        'model_level': 'Level 5+6',
        'base_regime': base_regime,
        'regime_hierarchy': regime_hierarchy,
        'regime_detail': regime_detail,
        'regime_subdetail': regime_subdetail,
        'dominant_limitation': dominant,
        'surface_regime': surface_class['surface_regime'],
        'coverage_mode': coverage_mode,
        'coverage_mode_info': surface_class['coverage_mode_info'],
        'coverage_depletion': surface_class['coverage_depletion'],
        'Da': Da,
        'theta': theta,
        'theta_equilibrium': theta_equilibrium,
        'SRF': SRF,
        'modification_factor': modification_factor,
        'surface_significant': surface_class['surface_significant'],
        'classification_depth': 3
    }