"""
Permeability reporting, tiered by what each quantity actually depends on.

Hydrogen changes chemical identity crossing a coated wall, and *where* it
dissociates decides whether a permeability is a material constant at all. This
codebase holds two deliberate hypotheses:

**Model 1** (L1-L5) puts dissociation at the oxide/metal interface. Molecular H2
crosses the oxide, so the oxide obeys Henry's law, :math:`C = K_{ox}P`, while the
metal obeys Sieverts, :math:`C = K_s\\sqrt{P}`. The two layers therefore have
*incommensurable* permeabilities -- mol/m/s/Pa against mol/m/s/Pa^0.5 -- and no
weighting combines them.

**Model 2** (L6 family) puts dissociation on the gas/oxide surface. Atomic H
crosses the oxide, so both layers obey Sieverts and the transport stack collapses
to a single permeance.

Four tiers follow, distinguished by what each value needs to be reproduced:

==========================  ==============  ==========================================
tier                        depends on      quantity
==========================  ==============  ==========================================
1 intrinsic                 nothing         ``D_m*K_s``, ``D_ox*K_ox``
2 effective layer           nothing         ``D_eff*K_s``, ``Phi_ox_eff`` (cracks/GB)
3 stack                     ``L``           ``Phi_transport`` -- Model 2 only
4 apparent                  ``L`` and ``P`` ``Phi_app`` from the solved flux
==========================  ==============  ==========================================

Tiers 1-3 are reproducible from parameters alone. Tier 4 is a property of the
wall *and* the experiment, and is always returned with its operating point.

Nothing here changes any flux. Every solved flux is taken from the existing
solvers unmodified; this module only decides how to report it.

See also
--------
docs/theory/permeability.md : the physics, with worked numbers
"""

import numpy as np

# -----------------------------------------------------------------------------
# Model and unit tags
# -----------------------------------------------------------------------------

MODEL_1 = 'model_1'   # Henry oxide + Sieverts metal, dissociation at the interface
MODEL_2 = 'model_2'   # surface kinetics + Sieverts oxide + Sieverts metal

#: Permeability units. These are not interchangeable and the module refuses to
#: mix them -- that refusal is the point, not an inconvenience.
UNITS_HENRY    = 'mol/m/s/Pa'
UNITS_SIEVERTS = 'mol/m/s/Pa^0.5'

_FLUX_KEYS = ('flux_total', 'J_total', 'J_ss', 'flux')


class UnitMismatchError(ValueError):
    """Raised when two permeabilities with different sorption laws are combined."""


class TierError(ValueError):
    """Raised when a tier is requested that does not exist for the given model."""


# -----------------------------------------------------------------------------
# Primitives
# -----------------------------------------------------------------------------

def _d_sqrt_P(P_up, P_down):
    return np.sqrt(max(P_up, 0.0)) - np.sqrt(max(P_down, 0.0))


def permeance(flux, P_up, P_down):
    r"""Permeance :math:`\Pi = J/\Delta\sqrt{P}` [mol/m^2/s/Pa^0.5].

    Needs no thickness convention, which makes it the safest quantity to compare
    across levels. Relates exactly to the permeation reduction factor:
    ``PRF = permeance(bare) / permeance(coated)``.
    """
    d = _d_sqrt_P(P_up, P_down)
    return flux / d if d > 0 else np.nan


def apparent_permeability(flux, P_up, P_down, L_total):
    r"""Tier-4 apparent permeability :math:`\Phi_{app} = J L_{tot}/\Delta\sqrt{P}`.

    Derived from the flux the level actually solved, so it inherits every defect
    path, every microstructural correction and the surface kinetics. Only
    meaningful when quoted with ``(T, P_up, P_down)``.
    """
    if L_total <= 0:
        raise ValueError(f"L_total must be positive, got {L_total}")
    return permeance(flux, P_up, P_down) * L_total


def pressure_exponent(flux_fn, P_up, rel_step=1e-2):
    r""":math:`n = \mathrm{d}\ln J/\mathrm{d}\ln P` by a two-point log derivative.

    Diagnoses which step controls the wall without needing any model internals,
    which makes it an independent check on ``classify_regime``:

    - :math:`n \to 1/2` -- a Sieverts transport step controls
    - :math:`n \to 1` -- a Henry layer (Model 1) or the surface rate (Model 2)
    - :math:`n < 1/2` -- Model 2 at high pressure, where coverage saturation
      makes :math:`g(\theta)` grow more slowly than :math:`\sqrt{P}`

    ``flux_fn`` takes an upstream pressure and returns a flux.
    """
    J1, J2 = flux_fn(P_up), flux_fn(P_up * (1.0 + rel_step))
    if not (J1 > 0 and J2 > 0):
        return np.nan
    return float(np.log(J2 / J1) / np.log(1.0 + rel_step))


def extract_flux(result, flux_key=None):
    """Pull the whole-wall flux out of any solver result.

    Resolution order is ``flux_total``, ``J_total``, ``J_ss``, ``flux`` -- most
    specific first, so a parallel-path result never returns a single-path flux.
    """
    if flux_key is not None:
        return result[flux_key]
    if not isinstance(result, dict):
        return float(result)          # molecular_diffusion_flux returns a bare float
    for k in _FLUX_KEYS:
        if k in result:
            return result[k]
    raise KeyError(f"no flux key in result; looked for {_FLUX_KEYS}")


def permeability_from_result(result, P_up, P_down, L_total, flux_key=None):
    """``apparent_permeability`` applied to a solver result dict."""
    return apparent_permeability(extract_flux(result, flux_key), P_up, P_down, L_total)


# -----------------------------------------------------------------------------
# Tier 1 -- intrinsic material constants
# -----------------------------------------------------------------------------

def intrinsic_metal_permeability(D_m, K_s):
    r"""Tier 1: :math:`\Phi_m = D_m K_s` [mol/m/s/Pa^0.5]."""
    return {'permeability': D_m * K_s, 'units': UNITS_SIEVERTS, 'tier': 1}


def intrinsic_oxide_permeability(D_ox, K_ox, model):
    r"""Tier 1: :math:`\Phi_{ox} = D_{ox}K_{ox}`.

    Units depend on the model, and that is physical rather than bookkeeping:
    Model 1 dissolves molecular H2 (Henry, per Pa), Model 2 dissolves atomic H
    (Sieverts, per Pa^0.5).
    """
    if model not in (MODEL_1, MODEL_2):
        raise ValueError(f"model must be {MODEL_1!r} or {MODEL_2!r}, got {model!r}")
    return {
        'permeability': D_ox * K_ox,
        'units': UNITS_HENRY if model == MODEL_1 else UNITS_SIEVERTS,
        'tier': 1,
    }


# -----------------------------------------------------------------------------
# Tier 2 -- effective layer properties
# -----------------------------------------------------------------------------

def effective_metal_permeability(D_eff, K_s):
    r"""Tier 2: :math:`D_{eff}K_s`, metal GB enhancement and trapping folded in.

    The sorption law is untouched, so this stays thickness- and
    pressure-independent.
    """
    return {'permeability': D_eff * K_s, 'units': UNITS_SIEVERTS, 'tier': 2}


def oxide_only_permeability(Phi_ox, f_crack=0.0, gamma=0.1,
                            f_gb=0.0, beta=10.0, f_pinhole=0.0, model=MODEL_1):
    r"""Tier 2: free-standing defective oxide, cracks and grain boundaries only.

    Parallel paths under one sorption law, so conductances add by area:

    .. math::

        \frac{\Phi_{ox,eff}}{\Phi_{ox}}
        = (1 - f_{cr} - f_{gb}) + \frac{f_{cr}}{\gamma} + f_{gb}\beta

    Exact -- it is the per-path :math:`\alpha` the solvers already use
    (:math:`\alpha_{cr} = \alpha_{int}/\gamma`, :math:`\alpha_{gb} = \beta\alpha_{int}`)
    summed by area rather than solved path by path. Depends only on the *ratio*
    :math:`\gamma`, never on :math:`L_{ox}`.

    Raises
    ------
    ValueError
        If ``f_pinhole > 0``. A free-standing oxide with a pinhole has infinite
        permeance, and what crosses it is aperture flow -- viscous or Knudsen --
        not permeation. A pinhole's resistance is supplied entirely by whatever
        sits behind it, so it belongs to the assembly and never to the oxide.
        Use :func:`lumped_oxide_permeability` for the coupled wall.
    """
    if f_pinhole > 0:
        raise ValueError(
            f"f_pinhole={f_pinhole} is not defined for a free-standing oxide: with no "
            "metal behind it a pinhole has zero oxide resistance, so the permeance is "
            "infinite and the transport is aperture flow, not permeation. Use "
            "lumped_oxide_permeability() for the coupled wall."
        )
    if not 0.0 <= f_crack + f_gb <= 1.0:
        raise ValueError(f"f_crack + f_gb must lie in [0, 1], got {f_crack + f_gb}")
    if not 0.0 < gamma <= 1.0:
        raise ValueError(f"gamma = L_crack/L_ox must lie in (0, 1], got {gamma}")
    if beta < 1.0:
        raise ValueError(f"beta = D_gb/D_ox must be >= 1, got {beta}")

    ratio = (1.0 - f_crack - f_gb) + f_crack / gamma + f_gb * beta
    return {
        'permeability': Phi_ox * ratio,
        'enhancement_ratio': ratio,
        'units': UNITS_HENRY if model == MODEL_1 else UNITS_SIEVERTS,
        'tier': 2,
        'contrast': {'crack': 1.0 / gamma, 'grain_boundary': beta},
    }


def lumped_oxide_permeability(Phi_ox, L_ox, Phi_m, L_m, f_pinhole=0.0,
                              f_crack=0.0, gamma=0.1, f_gb=0.0, beta=10.0,
                              model=MODEL_1, metal_dominance_threshold=10.0):
    r"""Coupled wall: may the defect branches be replaced by one lumped oxide?

    In general no. Each parallel branch carries its own interface pressure -- the
    pinhole branch sits at :math:`P_{int}\approx P_{up}`, the intact branch far
    below -- whereas a lumped oxide forces a single interface pressure on the
    whole wall. Parallel-of-series is not series-of-parallel.

    The validity condition is physical and independent of the sorption law:
    lumping works when the **metal** dominates, because then every branch shares
    nearly the same interface pressure. Measured error against the true
    parallel-path flux at ``f_pinhole = 1e-6``:

    ========================  ==============
    :math:`R_{ox}/R_m`        lumped / true
    ========================  ==============
    0.11                      1.0002
    1.14                      1.0018
    3.42                      1.0054
    11.4                      1.0182
    34.2                      1.0555
    114                       1.196
    342                       1.678
    3.4e6                     up to 1.07e4
    ========================  ==============

    The error grows roughly as 0.16% per unit of :math:`R_{ox}/R_m`, hence the
    default threshold of 10 -- about 2% worst case. Raise it only with the table
    above in view.

    Returns
    -------
    dict
        ``permeability`` (``inf`` when ``f_pinhole > 0``), ``valid``,
        ``R_ox_over_R_m`` and ``reason``.
    """
    R_ox = L_ox / Phi_ox if Phi_ox > 0 else np.inf
    R_m  = L_m / Phi_m if Phi_m > 0 else np.inf
    ratio = R_ox / R_m if R_m > 0 else np.inf
    valid = ratio <= metal_dominance_threshold

    if f_pinhole > 0:
        value, reason = np.inf, (
            'pinhole present: oxide resistance is zero over that area, so a lumped '
            'oxide makes the whole wall metal-limited instead of f_pinhole of it'
        )
    else:
        value = oxide_only_permeability(
            Phi_ox, f_crack=f_crack, gamma=gamma, f_gb=f_gb, beta=beta, model=model
        )['permeability']
        reason = 'cracks and grain boundaries only: homogenises exactly'

    return {
        'permeability': value,
        'valid': bool(valid),
        'R_ox_over_R_m': float(ratio),
        'reason': reason,
        'units': UNITS_HENRY if model == MODEL_1 else UNITS_SIEVERTS,
        'tier': 2,
    }


# -----------------------------------------------------------------------------
# Tier 3 -- stack permeability (Model 2 only)
# -----------------------------------------------------------------------------

def transport_permeability(Phi_ox, L_ox, Phi_m, L_m, model=MODEL_2):
    r"""Tier 3: thickness-weighted stack permeability. **Model 2 only.**

    .. math::

        \Phi_{transport} = \frac{L_{ox}+L_m}{L_{ox}/\Phi_{ox} + L_m/\Phi_m}

    Exact and pressure-independent, because both layers share :math:`n = 1/2`.
    It is the harmonic mean of the two *permeances* -- note that the harmonic
    mean of the two *permeabilities* is a different and wrong quantity, off by
    the thickness ratio.

    Raises
    ------
    TierError
        For Model 1. There the oxide permeability is per Pa and the metal's per
        Pa^0.5, so no weighting combines them and no stack permeability exists.
        This is not a gap to be filled -- it is what a molecular-transport
        barrier on an atomic-transport substrate costs.
    """
    if model != MODEL_2:
        raise TierError(
            "tier 3 does not exist for Model 1: Phi_ox is mol/m/s/Pa and Phi_m is "
            "mol/m/s/Pa^0.5, so they cannot be combined by any weighting. Report "
            "them separately (tier 1) and use apparent_permeability() for the wall."
        )
    R_total = L_ox / Phi_ox + L_m / Phi_m
    return {
        'permeability': (L_ox + L_m) / R_total,
        'permeance': 1.0 / R_total,
        'units': UNITS_SIEVERTS,
        'tier': 3,
        'L_total': L_ox + L_m,
    }


def surface_efficiency(theta, K_eq, P_up, P_down=0.0):
    r"""Model 2 surface factor :math:`\eta_{surf}\in(0,1]`.

    .. math::

        \eta_{surf} = \frac{g(\theta) - \sqrt{P_{down}}}{\sqrt{P_{up}} - \sqrt{P_{down}}}

    The surface adds **no resistance** to the transport stack; it acts by
    replacing :math:`\sqrt{P_{up}}` with the virtual pressure :math:`g(\theta)`.
    Hence :math:`\Phi_{app} = \Phi_{transport}\,\eta_{surf}` exactly.

    This is the same quantity the Level 6 chapter reports as
    :math:`P_{up}/P_{virtual}`, where its square root reproduces the measured
    flux ratio -- written here as a factor on permeability.
    """
    from calculations.surface_kinetics import g_theta
    d = _d_sqrt_P(P_up, P_down)
    if d <= 0:
        return np.nan
    return float((g_theta(theta, K_eq) - np.sqrt(max(P_down, 0.0))) / d)


# -----------------------------------------------------------------------------
# Guard against combining incommensurable quantities
# -----------------------------------------------------------------------------

def combine_series(a, b, L_a, L_b):
    """Thickness-weighted series combination, refusing to mix sorption laws.

    Raises
    ------
    UnitMismatchError
        If the two operands carry different units. Model 1 always hits this, by
        construction.
    """
    if a['units'] != b['units']:
        raise UnitMismatchError(
            f"cannot combine {a['units']} with {b['units']}: the two layers obey "
            "different sorption laws, so their permeabilities are not commensurable"
        )
    R = L_a / a['permeability'] + L_b / b['permeability']
    return {'permeability': (L_a + L_b) / R, 'units': a['units'],
            'tier': 3, 'L_total': L_a + L_b}


# -----------------------------------------------------------------------------
# Model 1 closed form for the pristine bilayer
# -----------------------------------------------------------------------------

def model1_bilayer_closed_form(Phi_ox, L_ox, Phi_m, L_m, P_up, P_down=0.0):
    r"""Exact Model 1 bilayer flux -- a Henry layer matched to a Sieverts layer.

    Equating the two fluxes with :math:`\alpha = \Phi_{ox}/L_{ox}` and
    :math:`\beta = \Phi_m/L_m`,

    .. math::

        \alpha(P_{up} - P_{int}) = \beta(\sqrt{P_{int}} - \sqrt{P_{down}})

    is a quadratic in :math:`\sqrt{P_{int}}`:

    .. math::

        \sqrt{P_{int}} = \frac{-\beta + \sqrt{\beta^2
                         + 4\alpha(\alpha P_{up} + \beta\sqrt{P_{down}})}}{2\alpha}

    Replaces a ``brentq`` root-find with an expression, and gives
    ``interface_solver`` a free regression test.
    """
    alpha, beta = Phi_ox / L_ox, Phi_m / L_m
    sd = np.sqrt(max(P_down, 0.0))
    sqrt_P_int = (-beta + np.sqrt(beta**2 + 4.0 * alpha * (alpha * P_up + beta * sd))) / (2.0 * alpha)
    return {'flux': beta * (sqrt_P_int - sd),
            'P_interface': sqrt_P_int**2,
            'alpha': alpha, 'beta': beta}


# -----------------------------------------------------------------------------
# Per-level reporting
# -----------------------------------------------------------------------------

def _row(level, model, tier, flux, P_up, P_down, L_total, units,
         temperature=None, Phi_transport=None, eta_surf=None,
         pressure_exponent_value=None, closed_form=None, valid=True, **extra):
    """Assemble one level's permeability record.

    Tier 1-2 rows carry no operating point and no exponent; tier 3-4 rows always
    do. That asymmetry is enforced here so a caller cannot mistake an apparent
    value for a material constant.
    """
    row = {
        'level': level, 'model': model, 'tier': tier,
        'flux': flux,
        'permeability': apparent_permeability(flux, P_up, P_down, L_total),
        'permeance': permeance(flux, P_up, P_down),
        'units': units, 'L_total': L_total,
        'Phi_transport': Phi_transport, 'eta_surf': eta_surf,
        'pressure_exponent': pressure_exponent_value,
        'P_up': P_up, 'P_down': P_down, 'temperature': temperature,
        'valid': valid, 'closed_form': closed_form,
    }
    row.update(extra)
    return row


# --- Model 1 ------------------------------------------------------------------

def level1(D_m, K_s, L_m, P_up, P_down=0.0, temperature=None):
    """L1 -- pristine metal. Tier 1: the apparent value *is* ``D_m*K_s``."""
    from calculations.permeation_calc import calculate_simple_metal_flux
    flux = calculate_simple_metal_flux(D_m, K_s, L_m, P_up, P_down)['flux']
    return _row('L1', MODEL_1, 1, flux, P_up, P_down, L_m, UNITS_SIEVERTS,
                temperature, closed_form=D_m * K_s)


def level2a(D_ox, K_ox, L_ox, P_up, P_down=0.0, temperature=None):
    """L2a -- pristine oxide alone. Tier 1, and in Henry units (per Pa).

    ``permeability`` here is normalised by ``d(sqrt P)`` for cross-level
    comparability, so it is *not* ``D_ox*K_ox``; the intrinsic value is returned
    separately as ``closed_form``. The two differ by ``sqrt(P_up)`` -- which is
    exactly the Henry-vs-Sieverts mismatch, made visible rather than hidden.
    """
    from calculations.oxide_permeation import molecular_diffusion_flux
    flux = molecular_diffusion_flux(D_ox, K_ox, L_ox, P_up, P_down)
    return _row('L2a', MODEL_1, 1, flux, P_up, P_down, L_ox, UNITS_HENRY,
                temperature, closed_form=D_ox * K_ox,
                intrinsic_permeability=D_ox * K_ox)


def level4(D_lattice, K_s, L_m, P_up, P_down, temperature,
           microstructure_params, lattice_density=1.06e29,
           method='average', mode='both'):
    """L4 -- bare defective metal. Tier 2: apparent value equals ``D_eff*K_s``."""
    from calculations.permeation_calc import calculate_defective_metal_flux
    r = calculate_defective_metal_flux(
        D_lattice=D_lattice, K_s=K_s, thickness=L_m, P_up=P_up, P_down=P_down,
        temperature=temperature, microstructure_params=microstructure_params,
        lattice_density=lattice_density, method=method, mode=mode)
    return _row('L4', MODEL_1, 2, r['flux'], P_up, P_down, L_m, UNITS_SIEVERTS,
                temperature, closed_form=r['D_eff'] * K_s, D_eff=r['D_eff'])


def level2b(oxide_props, metal_props, P_up, P_down=0.0, temperature=None):
    """L2b -- pristine bilayer. Tier 4: mixed exponents, so pressure-dependent."""
    from calculations.interface_solver import calculate_oxide_metal_system
    r = calculate_oxide_metal_system(P_up, P_down, oxide_props, metal_props)
    L_tot = oxide_props['thickness'] + metal_props['thickness']
    cf = model1_bilayer_closed_form(
        oxide_props['D_ox'] * oxide_props['K_ox'], oxide_props['thickness'],
        metal_props['D_metal'] * metal_props['K_s_metal'], metal_props['thickness'],
        P_up, P_down)['flux']
    n = pressure_exponent(
        lambda p: calculate_oxide_metal_system(p, P_down, oxide_props, metal_props)['flux'], P_up)
    return _row('L2b', MODEL_1, 4, r['flux'], P_up, P_down, L_tot, UNITS_SIEVERTS,
                temperature, pressure_exponent_value=n, closed_form=cf,
                P_interface=r.get('P_interface'))


def level3(oxide_props, metal_props, defect_params, P_up, P_down=0.0, temperature=None):
    """L3 -- defective oxide + pristine metal. Tier 4: branch topology."""
    from calculations.parallel_oxide_defect_paths import calculate_parallel_path_flux
    r = calculate_parallel_path_flux(P_up, P_down, oxide_props, metal_props, defect_params)
    L_tot = oxide_props['thickness'] + metal_props['thickness']
    n = pressure_exponent(
        lambda p: calculate_parallel_path_flux(
            p, P_down, oxide_props, metal_props, defect_params)['flux_total'], P_up)
    return _row('L3', MODEL_1, 4, r['flux_total'], P_up, P_down, L_tot,
                UNITS_SIEVERTS, temperature, pressure_exponent_value=n,
                flux_intact=r.get('flux_intact_contribution'),
                flux_defect=r.get('flux_defect_contribution'))


def level5(oxide_props, metal_props, defect_params, P_up, P_down, temperature,
           microstructure_params, lattice_density=1.06e29,
           method='average', mode='both'):
    """L5 -- defective oxide + defective metal. Tier 4."""
    from calculations.parallel_oxide_defect_paths import (
        calculate_parallel_path_flux_defective_metal as solve)
    kw = dict(oxide_props=oxide_props, metal_props=metal_props,
              defect_params=defect_params, temperature=temperature,
              microstructure_params=microstructure_params,
              lattice_density=lattice_density, method=method, n_points=10, mode=mode)
    r = solve(P_upstream=P_up, P_downstream=P_down, **kw)
    L_tot = oxide_props['thickness'] + metal_props['thickness']
    n = pressure_exponent(
        lambda p: solve(P_upstream=p, P_downstream=P_down, **kw)['flux_total'], P_up)
    return _row('L5', MODEL_1, 4, r['flux_total'], P_up, P_down, L_tot,
                UNITS_SIEVERTS, temperature, pressure_exponent_value=n,
                D_eff=r.get('D_eff_metal'),
                flux_intact=r.get('flux_intact_contribution'),
                flux_defect=r.get('flux_defect_contribution'))


# --- Model 2 ------------------------------------------------------------------

def _m2_row(level, r, P_up, P_down, L_total, temperature, K_eq,
            Phi_ox=None, L_ox=None, Phi_m=None, L_m=None, flux_key=None, **extra):
    """Model 2 row, with the tier-3 stack value and the surface factor attached."""
    flux = extract_flux(r, flux_key)
    theta = r.get('theta', r.get('theta_surface'))
    Phi_t = None
    if None not in (Phi_ox, L_ox, Phi_m, L_m):
        Phi_t = transport_permeability(Phi_ox, L_ox, Phi_m, L_m, MODEL_2)['permeability']
    eta = surface_efficiency(theta, K_eq, P_up, P_down) if theta is not None else None
    return _row(level, MODEL_2, 4, flux, P_up, P_down, L_total, UNITS_SIEVERTS,
                temperature, Phi_transport=Phi_t, eta_surf=eta,
                theta=theta, **extra)


def level1_L6(P_up, P_down, L_m, k_diss, K_eq, D_m, K_s_m, temperature=None):
    """L1+L6 -- pristine metal + surface. Tier 4: a rate constant demotes it.

    L1 alone is tier 1. Switching the surface on costs it that status with no
    change to the transport physics at all.
    """
    from calculations.surface_kinetics import solve_steady_state_flux_L1L6
    r = solve_steady_state_flux_L1L6(P_up, P_down, L_m, k_diss, K_eq, D_m, K_s_m)
    n = pressure_exponent(
        lambda p: solve_steady_state_flux_L1L6(p, P_down, L_m, k_diss, K_eq, D_m, K_s_m)['J_ss'], P_up)
    row = _m2_row('L1L6', r, P_up, P_down, L_m, temperature, K_eq)
    row['pressure_exponent'] = n
    row['Phi_transport'] = D_m * K_s_m          # single layer: the stack is the metal
    return row


def level2a_L6(P_up, P_down, k_diss, K_eq, D_ox, K_ox, L_ox, temperature=None):
    """L2a+L6 -- pristine oxide + surface, no metal block. Tier 4."""
    from calculations.surface_kinetics import solve_steady_state_flux_L2aL6
    r = solve_steady_state_flux_L2aL6(P_up, P_down, k_diss, K_eq, D_ox, K_ox, L_ox)
    n = pressure_exponent(
        lambda p: solve_steady_state_flux_L2aL6(p, P_down, k_diss, K_eq, D_ox, K_ox, L_ox)['J_ss'], P_up)
    row = _m2_row('L2aL6', r, P_up, P_down, L_ox, temperature, K_eq)
    row['pressure_exponent'] = n
    row['Phi_transport'] = D_ox * K_ox
    return row


def level2b_L6(P_up, P_down, L_m, k_diss, K_eq, D_ox, K_ox, L_ox, D_m, K_s_m,
               temperature=None):
    """L2b+L6 -- pristine bilayer + surface. Tier 3 stack times eta_surf."""
    from calculations.surface_kinetics import solve_steady_state_flux_direct
    args = (L_m, k_diss, K_eq, D_ox, K_ox, L_ox, D_m, K_s_m)
    r = solve_steady_state_flux_direct(P_up, P_down, *args)
    n = pressure_exponent(lambda p: solve_steady_state_flux_direct(p, P_down, *args)['J_ss'], P_up)
    row = _m2_row('L2bL6', r, P_up, P_down, L_ox + L_m, temperature, K_eq,
                  Phi_ox=D_ox * K_ox, L_ox=L_ox, Phi_m=D_m * K_s_m, L_m=L_m)
    row['pressure_exponent'] = n
    return row


def level3_L6(P_up, P_down, L_m, k_diss, K_eq, D_ox, K_ox, L_ox, D_m, K_s_m,
              defect_config, k_diss_metal=None, K_eq_metal=None, temperature=None):
    """L3+L6 -- defective oxide + surface. Tier 4: branch topology.

    ``eta_surf`` is reported from the intact path only. The pinhole branch runs
    on the *metal's* surface kinetics (``k_diss_metal``, ``K_eq_metal``), so a
    wall-level eta would blend two different surface chemistries; ``eta_blended``
    flags when that is the case.
    """
    from calculations.surface_kinetics import calculate_mixed_defect_flux_L6
    r = calculate_mixed_defect_flux_L6(
        P_up, P_down, L_m, k_diss, K_eq, D_ox, K_ox, L_ox, D_m, K_s_m,
        defect_config, k_diss_metal=k_diss_metal, K_eq_metal=K_eq_metal)
    theta_intact = r.get('intact_path', {}).get('theta')
    eta = surface_efficiency(theta_intact, K_eq, P_up, P_down) if theta_intact is not None else None
    has_pinhole = defect_config.get('pinhole', {}).get('area_fraction', 0) > 0
    return _row('L3L6', MODEL_2, 4, r['J_total'], P_up, P_down, L_ox + L_m,
                UNITS_SIEVERTS, temperature, eta_surf=eta,
                theta=theta_intact, eta_blended=bool(has_pinhole),
                dominant_path=r.get('dominant_path'))


def level4_L6(P_up, P_down, L_m, temperature, k_diss, K_eq, D_ox, K_ox, L_ox,
              D_lattice, K_s_m, microstructure_params,
              lattice_density=1.06e29, method='average', mode='both'):
    """L4+L6 -- pristine oxide + defective metal + surface. Tier 3 x eta_surf."""
    from calculations.surface_kinetics import calculate_defective_metal_flux_L6
    r = calculate_defective_metal_flux_L6(
        P_up=P_up, P_down=P_down, thickness=L_m, temperature=temperature,
        k_diss=k_diss, K_eq=K_eq, D_ox=D_ox, K_ox=K_ox, L_ox=L_ox,
        D_lattice=D_lattice, K_s_m=K_s_m,
        microstructure_params=microstructure_params,
        lattice_density=lattice_density, method=method, mode=mode)
    row = _m2_row('L4L6', r, P_up, P_down, L_ox + L_m, temperature, K_eq,
                  Phi_ox=D_ox * K_ox, L_ox=L_ox,
                  Phi_m=r['D_eff'] * K_s_m, L_m=L_m, flux_key='flux',
                  D_eff=r['D_eff'])
    return row


def level5_L6(P_up, P_down, L_m, temperature, k_diss, K_eq, D_ox, K_ox, L_ox,
              D_lattice, K_s_m, microstructure_params, defect_config,
              lattice_density=1.06e29, method='average', mode='both',
              k_diss_metal=None, K_eq_metal=None):
    """L5+L6 -- the full model. Tier 4."""
    from calculations.surface_kinetics import calculate_full_model_flux_L346_v2
    r = calculate_full_model_flux_L346_v2(
        P_up=P_up, P_down=P_down, L_m=L_m, temperature=temperature,
        k_diss=k_diss, K_eq=K_eq, D_ox=D_ox, K_ox=K_ox, L_ox=L_ox,
        D_lattice=D_lattice, K_s_m=K_s_m,
        microstructure_params=microstructure_params, defect_config=defect_config,
        lattice_density=lattice_density, method=method, mode=mode,
        k_diss_metal=k_diss_metal, K_eq_metal=K_eq_metal)
    theta_intact = r.get('intact_path', {}).get('theta')
    eta = surface_efficiency(theta_intact, K_eq, P_up, P_down) if theta_intact is not None else None
    has_pinhole = defect_config.get('pinhole', {}).get('area_fraction', 0) > 0
    return _row('L5L6', MODEL_2, 4, r['J_total'], P_up, P_down, L_ox + L_m,
                UNITS_SIEVERTS, temperature, eta_surf=eta, theta=theta_intact,
                eta_blended=bool(has_pinhole), D_eff=r.get('D_eff_avg'),
                dominant_path=r.get('dominant_path'))
