r"""
Defective metal model for hydrogen diffusion (Level 4): grain-boundary enhancement.

Grain boundaries are fast diffusion paths. In a polycrystal they raise the metal's
effective diffusivity above the lattice value, and because they carry mobile
hydrogen they raise the steady-state permeation flux.

Theory
------
Parallel-path model, with GB volume fraction :math:`f_{gb} = 3\delta/d`:

.. math::

    D_{metal} = (1-f_{gb})\,D_{bulk} + f_{gb}\,D_{gb},
    \qquad
    \frac{D_{gb}}{D_{bulk}} = A\exp\!\left[\frac{Q_{bulk}-Q_{gb}}{RT}\right]

with :math:`Q_{gb} \approx 0.6\,Q_{bulk}`. The enhancement depends on temperature,
grain size and boundary type, but not on hydrogen concentration, so a wall has one
metal diffusivity and the Level 1 solvers apply to it unchanged.

Why there is no trapping here
-----------------------------
Hydrogen trapping (Oriani, McNabb-Foster) was removed from this model. At steady
state the traps neither fill nor empty, so only lattice hydrogen carries flux and
:math:`J = D_L\,\Delta C_L / L` whatever the trap density, binding energy or
saturation. Traps set the transient (time lag) and the hydrogen inventory, which
are outside the scope of a steady-state permeation model. See ``STEADY_STATE.md``.

Module structure
----------------
1. :func:`gb_enhancement_factor` -- :math:`D_{gb}/D_{bulk}` from the active metal's E_D
2. :func:`calculate_gb_enhanced_diffusivity` -- parallel-path (or Hart) D_metal
3. :func:`metal_diffusivity` -- the diffusivity every Level 4/5 solver uses

References
----------
1. Hart, E.W. (1957). "On the role of dislocations in bulk diffusion."
   Acta Metall. 5, 597. DOI: 10.1016/0001-6160(57)90127-X

2. Kaur, I., Mishin, Y., Gust, W. (1995). "Fundamentals of Grain and
   Interphase Boundary Diffusion." John Wiley & Sons.

3. Oudriss, A., et al. (2012). "Grain size and grain-boundary effects on
   diffusion and trapping of hydrogen in pure nickel." Acta Mater. 60,
   6814-6828. DOI: 10.1016/j.actamat.2012.09.004

4. Tsuru, T., Latanision, R.M. (1982). "Grain boundary transport of hydrogen
   in nickel." Scripta Metall. 16, 575-578. DOI: 10.1016/0036-9748(82)90273-3
"""




def gb_enhancement_factor(temperature, temperature_unit='K', gb_type='HAGB', data_source='default'):
    r"""
    Calculate grain boundary diffusion enhancement factor D_gb/D_bulk.
    
    Theory
    ------
    Grain boundaries provide fast diffusion paths due to their more open atomic
    structure. The enhancement factor D_gb/D_bulk decreases with temperature as
    thermal activation makes bulk diffusion more competitive.
    
    The temperature dependence follows:
        .. math::

        \frac{D_{gb}}{D_{bulk}} = A\exp\!\left[\frac{Q_{bulk}-Q_{gb}}{RT}\right]
    
    where typically Q_gb ≈ 0.6 × Q_bulk, meaning grain boundary diffusion has
    lower activation energy. At low homologous temperatures (T/T_m < 0.5), GB
    diffusion can be 100-1000× faster than bulk. At high temperatures (T/T_m > 0.8),
    the enhancement reduces to 10-20×.

    **This is what the function computes.** Q_bulk is the active metal's ``E_D`` from
    the study config; the two dimensionless ratios (Q_gb/Q_bulk and the prefactor
    ratio A) are fixed in the body. Until 2026-09-28 the derivation below was prose
    only: the code interpolated a hardcoded austenitic-steel table instead, whatever
    metal was loaded.
    
    Mathematical Derivation
    -----------------------
    Starting from Arrhenius expressions:
        D_bulk = D₀_bulk × exp(-Q_bulk/RT)
        D_gb = D₀_gb × exp(-Q_gb/RT)
    
    Taking the ratio:
        D_gb/D_bulk = (D₀_gb/D₀_bulk) × exp[(Q_bulk - Q_gb)/RT]
    
    With Q_gb ≈ 0.6 × Q_bulk (empirical):
        .. math::

        \frac{D_{gb}}{D_{bulk}} \approx A\exp\!\left[\frac{0.4\,Q_{bulk}}{RT}\right]
    
    This predicts log-linear behavior in 1/T space.
    
    Parameters
    ----------
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
        Where the enhancement comes from. Default: ``'default'``.

        - ``'default'``: compute alpha from the ACTIVE METAL's own ``E_D``, read from
          the METALS entry of the study config. No enhancement data is stored here.
        - a dict: an explicit table, ``{T_celsius: enhancement_factor}``, interpolated
          in log space. Note ``np.interp`` clamps outside the supplied range rather
          than extrapolating, so the nearest endpoint is returned there.
    
    Returns
    -------
    dict
        Dictionary containing:

        - 'enhancement_factor': D_gb/D_bulk ratio [-]
        - 'uncertainty': Estimated relative uncertainty (fractional); 0.2 inside the
          validated window, 0.5 outside it, x1.5 for a GB type other than HAGB
        - 'interpolation_method': one of 'arrhenius_from_config',
          'arrhenius_extrapolated' (default source), or 'interpolation',
          'extrapolation_warning' (explicit table)
        - 'temperature_K': Temperature in Kelvin [K]
        - 'temperature_C': Temperature in Celsius [C]
        - 'gb_type_factor': Scaling factor applied for GB type [-]
        - 'data_range_K': validity window — the metal's own ``temperature_range`` for
          the default source, or the supplied table's span for a dict
    
    Raises
    ------
    ValueError
        If temperature ≤ 0 K
        If temperature_unit not recognized
        If gb_type not recognized
        If data_source is neither 'default' nor a dict

    .. note::

       Being outside the validated temperature range **warns**, it does not raise.
       The docstring previously promised a ValueError beyond 200 K of extrapolation;
       no such check ever existed.
    
    Warnings
    --------
    UserWarning
        If temperature requires extrapolation beyond data range
        If enhancement factor seems unphysical (>10000 or <1)
    
    References
    ----------
    1. Tsuru, T., Latanision, R.M. (1982). "Grain boundary transport of hydrogen
       in nickel." Scripta Metall. 16, 575-578.
       DOI: 10.1016/0036-9748(82)90273-3
    
    2. Harris, J.E., Masters, B.C. (1966). "The diffusivity of hydrogen in nickel."
       British J. Appl. Phys. 17, 377. DOI: 10.1088/0508-3443/17/3/310
    
    3. Brass, A.M., Chanfreau, A. (1996). "Accelerated diffusion of hydrogen along
       grain boundaries in nickel." Acta Mater. 44, 3823-3831.
       DOI: 10.1016/1359-6454(95)00446-7
    
    Example
    -------
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
    
    if data_source == 'default':
        from calculations.config.model_config import ACTIVE_METAL, METALS, R
        
        Q_RATIO         = 0.6   # Q_gb / Q_bulk, the standard empirical value
        PREFACTOR_RATIO = 1.0   # A = D0_gb / D0_bulk

        metal  = METALS[ACTIVE_METAL]
        Q_bulk = metal['E_D']

        # D_gb/D_bulk = A exp[(Q_bulk - Q_gb)/RT], with Q_gb = Q_RATIO * Q_bulk
        base_enhancement = PREFACTOR_RATIO * np.exp(
            (1.0 - Q_RATIO) * Q_bulk / (R * T_kelvin))

        # Validity is now the window the metal's OWN diffusion data were measured over,
        # from the same METALS entry as E_D, instead of an unrelated table's range.
        window = metal.get('temperature_range')
        if window:
            T_min, T_max = float(window[0]), float(window[1])
            if T_min <= T_kelvin <= T_max:
                interpolation_method = 'arrhenius_from_config'
            else:
                interpolation_method = 'arrhenius_extrapolated'
                warnings.warn(
                    f"Temperature {T_kelvin:.1f} K is outside the validated range "
                    f"[{T_min:.1f}, {T_max:.1f}] K of {ACTIVE_METAL}'s diffusion data. "
                    "The GB enhancement is extrapolated.")
        else:
            T_min = T_max = float('nan')
            interpolation_method = 'arrhenius_from_config'

    elif isinstance(data_source, dict):
        # Explicit table override: {T_celsius: enhancement}. Interpolated in log space,
        # because log(D_gb/D_bulk) is linear in 1/T.
        data_T_celsius = np.array(list(data_source.keys()))
        data_enhancement = np.array(list(data_source.values()))
        data_T_kelvin = data_T_celsius + 273.15
        T_min, T_max = data_T_kelvin.min(), data_T_kelvin.max()

        if T_min <= T_kelvin <= T_max:
            interpolation_method = 'interpolation'
        else:
            interpolation_method = 'extrapolation_warning'
            warnings.warn(f"Temperature {T_kelvin:.1f} K is outside the supplied data "
                          f"range [{T_min:.1f}, {T_max:.1f}] K. np.interp clamps rather "
                          "than extrapolates, so the nearest endpoint value is returned.")

        log_enhancement = np.log(data_enhancement)
        inverse_T = 1.0 / data_T_kelvin
        log_factor_interpolated = np.interp(
            1.0/T_kelvin, inverse_T[::-1], log_enhancement[::-1])
        base_enhancement = np.exp(log_factor_interpolated)

    else:
        raise ValueError("data_source must be 'default' or a dictionary")
    
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
    
    # Estimate uncertainty (higher outside the validated window, and for GB types
    # whose scaling is itself a rule of thumb).
    #
    # One tier for "outside", deliberately. There used to be two, and only
    # 'extrapolation_warning' was tested here, so the FAR case fell through to the
    # interpolation value: going further from the data lowered the reported
    # uncertainty from 0.5 to 0.2.
    if interpolation_method in ('extrapolation_warning', 'arrhenius_extrapolated'):
        uncertainty = 0.5  # ±50% outside the validated range
    else:
        uncertainty = 0.2  # ±20% inside it
    
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


def calculate_gb_enhanced_diffusivity(D_bulk, temperature, grain_size, 
                                     gb_thickness=0.5e-9, gb_type='HAGB', 
                                     model='parallel'):
    r"""
    Calculate effective diffusivity with grain boundary fast diffusion paths.
    
    Theory
    ------
    Grain boundaries provide fast diffusion paths parallel to bulk diffusion.
    The effective diffusivity depends on the volume fraction of GBs and their
    enhancement factor. Two models are implemented:
    
    1. Parallel Path Model (simple):

       .. math::

           D_{eff} = f_{bulk}D_{bulk} + f_{gb}D_{gb}
    
    2. Hart Equation (accounts for GB connectivity):
       More accurate for higher GB fractions (f_gb > 0.01)
    
    The GB volume fraction is:

    .. math::

        f_{gb} = \frac{3\delta}{d}
    
    where δ is GB thickness and d is grain size.
    
    Mathematical Derivation
    -----------------------
    Parallel model assumes independent transport:
       .. math::

       J_{total} = J_{bulk} + J_{gb}
           = -D_{bulk}\nabla C\,A_{bulk} - D_{gb}\nabla C\,A_{gb}
    
    With area fractions equal to volume fractions:
       .. math::

       D_{eff} = (1-f_{gb})D_{bulk} + f_{gb}D_{gb}
    
    Hart equation includes connectivity effects:
       D_eff/D_bulk = 1 + [f_gb/(1-f_gb)]×[(D_gb/D_bulk-1)]/[1+2(D_gb/D_bulk-1)f_gb/(3(1-f_gb))]
    
    This accounts for GB network percolation at high f_gb.
    
    Parameters
    ----------
    D_bulk : float
        Bulk lattice diffusion coefficient [m²/s]
        Typical: 10⁻¹² to 10⁻⁷ m²/s for H in metals
    
    temperature : float
        Operating temperature [K] (all properties should be evaluated at this T)
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
    
    Returns
    -------
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
    
    Raises
    ------
    ValueError
        If D_bulk ≤ 0
        If temperature ≤ 0
        If grain_size ≤ 0
        If gb_thickness ≤ 0 or ≥ grain_size/2
        If model not recognized
    
    Warnings
    --------
    UserWarning
        If f_gb > 0.1 (GB percolation threshold)
        If f_gb > 0.5 (unphysical)
        If enhancement ratio > 100 (may indicate error)
    
    References
    ----------
    1. Hart, E.W. (1957). "On the role of dislocations in bulk diffusion."
       Acta Metall. 5, 597. DOI: 10.1016/0001-6160(57)90127-X
    
    2. Kaur, I., Mishin, Y., Gust, W. (1995). "Fundamentals of Grain and
       Interphase Boundary Diffusion." John Wiley & Sons.
    
    3. Divinski, S.V., Wilde, G. (2008). "Grain boundary self-diffusion in
       polycrystalline nickel." Z. Metallkd. 99, 8-15.
       DOI: 10.3139/146.101686
    
    Example
    -------
    >>> result = calculate_gb_enhanced_diffusivity(
    ...     D_bulk=1e-10,      # m²/s
    ...     temperature=800,    # K
    ...     grain_size=50e-6,   # 50 μm
    ...     gb_type='HAGB',
    ...     model='parallel'
    ... )
    >>> print(f"Enhancement: {result['enhancement_ratio']:.2f}×")
    """
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
    # f_gb = 3δ/d for randomly oriented grains.
    #
    # This is the METAL's grain-boundary VOLUME fraction, of order 1e-5 at a 100 um
    # grain. It is not the oxide's grain-boundary DEFECT area fraction, which is of
    # order 1e-3 and is called f_gb_defect in the config and the sensitivity layer.
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


def metal_diffusivity(D_lattice, temperature, microstructure_params,
                      enable_gb_enhancement=True):
    """
    The metal diffusivity every Level 4/5 solver uses.

    The lattice value, raised by grain-boundary fast paths when
    ``enable_gb_enhancement`` is true. It does not depend on hydrogen
    concentration, so it is computed once and handed to the Level 1-3 solvers in
    place of ``D_lattice``.

    Parameters
    ----------
    D_lattice : float
        Lattice diffusion coefficient [m²/s]
    temperature : float
        Temperature [K]
    microstructure_params : dict
        Needs ``grain_size`` and ``gb_type``; ``gb_thickness`` (default 0.5 nm)
        and ``model`` (``'parallel'`` or ``'hart'``, default ``'parallel'``) are
        optional.
    enable_gb_enhancement : bool, optional
        False returns ``D_lattice`` unchanged. Default True.

    Returns
    -------
    float
        Metal diffusivity [m²/s], never below ``D_lattice``.
    """
    if not enable_gb_enhancement:
        return D_lattice
    return calculate_gb_enhanced_diffusivity(
        D_bulk=D_lattice,
        temperature=temperature,
        grain_size=microstructure_params['grain_size'],
        gb_thickness=microstructure_params.get('gb_thickness', 0.5e-9),
        gb_type=microstructure_params['gb_type'],
        model=microstructure_params.get('model', 'parallel'),
    )['D_eff']
