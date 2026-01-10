# data/oxide_defect_parameters.py
"""
Defect parameters for Level 3 parallel path model.

Based on literature values from:
- Strehlow & Savage (1974): Original parallel path model
- Zarchy & Axtmann (1979): 6Å oxide with 1% defects
- Zhang et al. (2018): PRF measurements
"""

# Standard defect configurations based on literature
DEFECT_CONFIGURATIONS = {
    'perfect': {
        'description': 'No defects - baseline for comparison',
        'area_fraction': 0.0,
        'type': 'pinhole',
        'reference': 'Level 2 baseline'
    },
    
    'minimal_pinholes': {
        'description': 'Fresh oxide with minimal pinholes',
        'area_fraction': 0.001,  # 0.1%
        'type': 'pinhole',
        'reference': 'Typical for well-prepared surfaces'
    },
    
    'typical_pinholes': {
        'description': 'Normal operating conditions with pinholes',
        'area_fraction': 0.01,   # 1%
        'type': 'pinhole',
        'reference': 'Zarchy & Axtmann (1979) - 1% defects observed'
    },
    
    'degraded_pinholes': {
        'description': 'Degraded oxide with significant pinholes',
        'area_fraction': 0.1,    # 10%
        'type': 'pinhole',
        'reference': 'Upper limit for damaged oxides'
    },
    
    'thin_cracks': {
        'description': 'Cracks with thin oxide (10% normal thickness)',
        'area_fraction': 0.01,   # 1%
        'type': 'crack',
        'thickness_factor': 0.1,  # Oxide in crack is 10% of normal
        'reference': 'Thermal cycling damage'
    },
    
    'thick_cracks': {
        'description': 'Cracks with partial oxide (50% normal thickness)',
        'area_fraction': 0.01,   # 1%
        'type': 'crack',
        'thickness_factor': 0.5,  # Oxide in crack is 50% of normal
        'reference': 'Mechanical stress cracks'
    },
    
    'grain_boundaries': {
        'description': 'Enhanced diffusion through grain boundaries',
        'area_fraction': 0.02,   # 2%
        'type': 'grain_boundary',
        'diffusivity_factor': 10,  # 10x faster diffusion
        'reference': 'Typical polycrystalline oxide'
    },
    
    'mixed_defects': {
        'description': 'Realistic mixed defect population',
        'area_fraction': 0.015,  # 1.5% total
        'type': 'mixed',
        'components': {
            'pinholes': 0.005,    # 0.5%
            'cracks': 0.005,      # 0.5%
            'grain_boundaries': 0.005  # 0.5%
        },
        'reference': 'Real-world oxide'
    }
}

# Parameter ranges for sensitivity studies
PARAMETER_RANGES = {
    'area_fraction': {
        'min': 0.0,
        'max': 0.2,
        'typical': 0.01,
        'units': 'fraction',
        'description': 'Fraction of surface area with defects'
    },
    
    'thickness_factor': {
        'min': 0.0,
        'max': 1.0,
        'typical': 0.1,
        'units': 'fraction',
        'description': 'Oxide thickness in crack relative to intact'
    },
    
    'diffusivity_factor': {
        'min': 1,
        'max': 100,
        'typical': 10,
        'units': 'multiplier',
        'description': 'Enhancement of diffusivity in grain boundaries'
    }
}

# Expected PRF ranges from literature
PRF_RANGES = {
    'perfect_oxide': {
        'min': 100,
        'max': 10000,
        'reference': 'Theoretical maximum'
    },
    
    'minimal_defects': {
        'min': 50,
        'max': 1000,
        'reference': 'Well-maintained barriers'
    },
    
    'typical_defects': {
        'min': 10,
        'max': 100,
        'reference': 'Normal operating conditions'
    },
    
    'severe_defects': {
        'min': 2,
        'max': 10,
        'reference': 'Degraded barriers'
    }
}