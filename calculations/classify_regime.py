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
