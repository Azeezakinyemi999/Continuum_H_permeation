import numpy as np

def sieverts_concentration(K_s, pressure):
    r"""
    Calculate hydrogen concentration at metal surface using Sieverts' law.
    
    Parameters
    ----------
    K_s : float
        Solubility constant [mol/m³/Pa^0.5]
    pressure : float
        Hydrogen partial pressure [Pa]
    
    Returns
    -------
    float
        Hydrogen concentration [mol/m³]
    
    Notes
    -----
    Sieverts' law:

    .. math::

        C = K_s \sqrt{P}

    Valid for molecular hydrogen dissociating at surface
    """
    if pressure < 0:
        raise ValueError(f"Pressure cannot be negative: {pressure} Pa")
    
    concentration = K_s * np.sqrt(pressure)
    return concentration


def fick_flux(D, C_up, C_down, layer_thickness):
    r"""
    Calculate diffusive flux using Fick's first law..
    
    Parameters
    ----------
    D : float
        Diffusion coefficient [m²/s]
    C_up : float
        Upstream concentration [mol/m³]
    C_down : float
        Downstream concentration [mol/m³]
    layer_thickness : float
        Material layer_thickness [m]
    
    Returns
    -------
    float
        Diffusive flux [mol/m²/s]
    
    Notes
    -----
    Fick's first law, and its steady-state form through a flat plate where the
    concentration gradient is constant:

    .. math::

        J = -D\frac{dC}{dx}
        \qquad\Longrightarrow\qquad
        J = \frac{D\,(C_{up} - C_{down})}{L}

    Positive flux is from upstream to downstream
    """
    if layer_thickness <= 0:
        raise ValueError(f"Thickness must be positive: {layer_thickness} m")
    if D < 0:
        raise ValueError(f"Diffusion coefficient cannot be negative: {D} m²/s")
    
    flux = D * (C_up - C_down) / layer_thickness
    return flux


def calculate_simple_metal_flux(D, K_s, metal_thickness, P_up, P_down):
    """
    Calculate hydrogen permeation flux through clean metal..
    
    Combines Sieverts' law for surface concentration with Fick's law
    for bulk diffusion.
    
    Parameters
    ----------
    D : float
        Diffusion coefficient [m²/s]
    K_s : float
        Solubility constant [mol/m³/Pa^0.5]
    metal_thickness : float
        Metal metal_thickness [m]
    P_up : float
        Upstream hydrogen pressure [Pa]
    P_down : float
        Downstream hydrogen pressure [Pa]
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'flux': Permeation flux [mol/m²/s]
        - 'C_up': Upstream concentration [mol/m³]
        - 'C_down': Downstream concentration [mol/m³]
        - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
    
    Notes
    -----
    This assumes:
    - Clean metal surfaces (no oxide)
    - Equilibrium at surfaces (fast dissociation/recombination)
    - Steady-state diffusion
    - No trapping or other complications
    """
    # Input validation
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    if P_down > P_up:
        print(f"Warning: Downstream pressure ({P_down} Pa) > Upstream pressure ({P_up} Pa)")
        print("Flux will be negative (backward flow)")
    
    # Calculate surface concentrations using Sieverts' law
    C_up = sieverts_concentration(K_s, P_up)
    C_down = sieverts_concentration(K_s, P_down)
    
    # Calculate flux using Fick's law
    flux = fick_flux(D, C_up, C_down, metal_thickness)
    
    # Calculate effective permeability
    # P = D * K_s for metal membranes
    permeability = D * K_s
    
    # Return comprehensive results
    return {
        'flux': flux,
        'C_up': C_up,
        'C_down': C_down,
        'permeability': permeability,
        'Diffusivity': D,
        'solubility': K_s,
        'units': {
            'flux': 'mol/m²/s',
            'concentration': 'mol/m³',
            'permeability': 'mol/m/s/Pa^0.5',
            'Diffusivity': 'm²/s',
            'solubility': 'mol/m³/Pa^0.5'
        }
    }


# =============================================================================
# LEVEL 4: Defective Metal Flux Calculations
# =============================================================================


def calculate_defective_metal_flux(D_lattice, K_s, metal_thickness, P_up, P_down,
                                   temperature, microstructure_params,
                                   enable_gb_enhancement=True):
    r"""
    Steady-state hydrogen flux through a polycrystalline metal (Level 4).

    Level 1 with the lattice diffusivity replaced by the grain-boundary-enhanced
    one from :func:`calculations.defective_metal.metal_diffusivity`:

    .. math::

        J = \frac{D_{metal}\,K_s\left(\sqrt{P_{up}}-\sqrt{P_{down}}\right)}{L}

    Trapping is deliberately absent. At steady state the traps are in equilibrium
    with the lattice and exchange no net hydrogen, so only lattice hydrogen carries
    flux and :math:`J = D\,\Delta C_L/L` regardless of trap density or binding
    energy. See ``STEADY_STATE.md``.

    Parameters
    ----------
    D_lattice : float
        Lattice diffusion coefficient [m²/s]
    K_s : float
        Sieverts' constant [mol/m³/Pa^0.5]
    metal_thickness : float
        Metal thickness [m]
    P_up, P_down : float
        Upstream and downstream hydrogen pressure [Pa]
    temperature : float
        Temperature [K]
    microstructure_params : dict
        ``grain_size`` and ``gb_type``; optionally ``gb_thickness`` and ``model``.
    enable_gb_enhancement : bool, optional
        False reduces this to Level 1. Default True.

    Returns
    -------
    dict
        ``flux`` [mol/m²/s], ``C_up``, ``C_down`` [mol/m³], ``permeability``
        (D_eff·K_s) [mol/m/s/Pa^0.5], ``D_eff`` (the metal diffusivity used)
        and ``D_lattice`` [m²/s], ``modification_factor`` (D_eff/D_lattice, ≥ 1)
        and ``regime`` (always ``'metal_limited'``: a bare wall has nothing else
        to limit it).
    """
    from calculations.defective_metal import metal_diffusivity

    if temperature <= 0:
        raise ValueError(f"Temperature must be positive: {temperature} K")

    if D_lattice <= 0:
        raise ValueError(f"Diffusion coefficient must be positive: {D_lattice} m²/s")
    if K_s <= 0:
        raise ValueError(f"Solubility constant must be positive: {K_s}")
    D_eff = metal_diffusivity(D_lattice, temperature, microstructure_params,
                              enable_gb_enhancement)
    result = calculate_simple_metal_flux(D_eff, K_s, metal_thickness, P_up, P_down)

    return {
        'flux': result['flux'],
        'C_up': result['C_up'],
        'C_down': result['C_down'],
        'permeability': result['permeability'],
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'modification_factor': D_eff / D_lattice,
        'regime': 'metal_limited',
        'units': {
            'flux': 'mol/m²/s',
            'concentration': 'mol/m³',
            'permeability': 'mol/m/s/Pa^0.5',
            'diffusivity': 'm²/s'
        }
    }
