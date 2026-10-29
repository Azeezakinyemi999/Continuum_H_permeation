r"""The Example Wall: a small fictitious system that every figure and check runs on.

Nothing here is a real material. The numbers are round powers of ten chosen so that

* every intermediate quantity can be checked by hand,
* the oxide and the metal carry *comparable* resistance, so the series solve in
  :mod:`calculations.interface_solver` is non-trivial rather than collapsing to one
  layer, and
* the Level 2b interface pressure has an **exact closed form** (see
  :func:`exact_interface_pressure`), which gives :mod:`check_equations` an analytical
  target instead of comparing the code against itself.

Why a fictitious wall at all. The documentation's figures must not report results from
a real run, and must not change meaning when the active study changes. A real study
would do both. See ``docs/CHECKPOINT.md`` decision D5.

Run it directly to print the whole worked example::

    python docs/figures/example_wall.py
"""

from __future__ import annotations

import math

# -----------------------------------------------------------------------------
# The wall
# -----------------------------------------------------------------------------

#: Operating point. T is used only by the temperature-dependent figures; the
#: transport numbers below are stated AT this temperature rather than extrapolated
#: to it, so the example needs no Arrhenius reference data.
T_K = 800.0
P_UP = 1.0e5      # Pa
P_DOWN = 0.0      # Pa

#: Metal: Sieverts sorption, atomic H. Round numbers, not a real alloy.
D_M = 1.0e-10     # m^2/s
K_S_M = 1.0       # mol/m^3/Pa^0.5
L_M = 1.0e-3      # m   (1 mm wall)

#: Oxide: Henry sorption, molecular H2 (Model 1). Also not a real oxide.
D_OX = 1.0e-16    # m^2/s
K_OX = 1.0        # mol/m^3/Pa
L_OX = 1.0e-6     # m   (1 um coating)

#: Activation energies, for the temperature figures only. Round numbers again, and
#: positive means "increases with T" under the reference-temperature Arrhenius form in
#: calculations.utils. The permeability activation energy is their sum,
#: Q_perm = E_D + H_s = 60 kJ/mol, which is what an Arrhenius plot of flux recovers.
E_D_M = 50.0e3     # J/mol, metal diffusivity
H_S_M = 10.0e3     # J/mol, metal solubility
E_D_OX = 100.0e3   # J/mol, oxide diffusivity
H_OX = 20.0e3      # J/mol, oxide solubility

#: Microstructure, for the Level 4 figures.
N_L = 1.0e29      # m^-3, interstitial lattice sites (FCC order of magnitude)
TRAP = {
    "name": "example_trap",
    "binding_energy": 80.0e3,   # J/mol -- strong enough that the effect is visible
    "density": 1.0e24,          # m^-3
}
GRAIN_SIZE = 20.0e-6          # m
GB_THICKNESS = 0.5e-9         # m
GB_DIFFUSIVITY_RATIO = 100.0  # D_gb / D_lattice

#: Defects in the oxide, for the Level 3 figures.
F_DEFECT = 1.0e-3             # area fraction
DEFECT_TYPE = "pinhole"

#: Surface kinetics, for the Level 6 figures.
K_DISS = 1.0e-8               # mol/m^2/s/Pa
K_EQ = 1.0e-3                 # Pa^-1

R_GAS = 8.314                 # J/mol/K -- matches calculations.utils.R


# -----------------------------------------------------------------------------
# Derived quantities
# -----------------------------------------------------------------------------

def permeances():
    r"""Oxide and metal permeances, :math:`\alpha` and :math:`\beta`.

    .. math::

        \alpha = \frac{D_{ox}K_{ox}}{L_{ox}}
        \qquad
        \beta = \frac{D_m K_{s,m}}{L_m}

    Note the units do **not** match: :math:`\alpha` is per Pa (Henry) and
    :math:`\beta` per :math:`\sqrt{\mathrm{Pa}}` (Sieverts). That mismatch is
    physical, not a bookkeeping slip -- it is the whole subject of the two-models
    chapter -- and it is why they may appear in one flux balance but never be
    summed.

    Returns
    -------
    (float, float)
        ``(alpha, beta)``; for the Example Wall, exactly ``(1e-10, 1e-7)``.
    """
    alpha = D_OX * K_OX / L_OX
    beta = D_M * K_S_M / L_M
    return alpha, beta


def exact_interface_pressure(P_up=P_UP, P_down=P_DOWN):
    r"""The Level 2b interface pressure in closed form, for :math:`P_{down}=0`.

    Flux balance across the two layers is

    .. math::

        \underbrace{\alpha\,(P_{up} - P_{int})}_{\text{oxide, Henry}}
        =
        \underbrace{\beta\,(\sqrt{P_{int}} - \sqrt{P_{down}})}_{\text{metal, Sieverts}}

    With :math:`P_{down}=0` and :math:`x=\sqrt{P_{int}}` this is an ordinary
    quadratic, :math:`\alpha x^2 + \beta x - \alpha P_{up} = 0`, whose positive root is

    .. math::

        \sqrt{P_{int}} = \frac{-\beta + \sqrt{\beta^2 + 4\alpha^2 P_{up}}}{2\alpha}

    This is the analytical value :mod:`check_equations` tests ``brentq`` against. It
    exists only because the oxide is linear in :math:`P` and the metal in
    :math:`\sqrt{P}`; add a defect path or surface kinetics and the closed form is gone.

    Raises
    ------
    ValueError
        If ``P_down`` is non-zero, where no such closed form applies.
    """
    if P_down != 0.0:
        raise ValueError(
            "the closed form holds only for P_down = 0; "
            f"got {P_down}. Use the solver for the general case."
        )
    alpha, beta = permeances()
    sqrt_P_int = (-beta + math.sqrt(beta**2 + 4.0 * alpha**2 * P_up)) / (2.0 * alpha)
    return sqrt_P_int**2


def exact_level1_flux(P_up=P_UP, P_down=P_DOWN):
    r""":math:`J = \frac{D_m K_{s,m}}{L_m}\left(\sqrt{P_{up}}-\sqrt{P_{down}}\right)`."""
    _, beta = permeances()
    return beta * (math.sqrt(P_up) - math.sqrt(P_down))


def exact_level2a_flux(P_up=P_UP, P_down=P_DOWN):
    r""":math:`J = \frac{D_{ox}K_{ox}}{L_{ox}}\left(P_{up}-P_{down}\right)`."""
    alpha, _ = permeances()
    return alpha * (P_up - P_down)


def exact_level2b_flux(P_up=P_UP, P_down=P_DOWN):
    """Series flux, evaluated on the metal side of the closed-form interface pressure."""
    _, beta = permeances()
    P_int = exact_interface_pressure(P_up, P_down)
    return beta * (math.sqrt(P_int) - math.sqrt(P_down))


def exact_trapping_term(temperature=T_K, trap=None, lattice_density=N_L):
    r"""Oriani trapping term :math:`\sum_i N_{T,i}K_i/N_L`, with :math:`K=e^{E_b/RT}`.

    The effective diffusivity is :math:`D_{eff} = D_L/(1 + \text{this})`.
    """
    trap = trap or TRAP
    K = math.exp(trap["binding_energy"] / (R_GAS * temperature))
    return trap["density"] * K / lattice_density


def exact_D_eff_trapping(D_lattice=D_M, temperature=T_K):
    r""":math:`D_{eff} = D_L\big/\left(1 + \sum_i N_{T,i}K_i/N_L\right)`."""
    return D_lattice / (1.0 + exact_trapping_term(temperature))


def exact_g_theta(theta, K_eq=K_EQ):
    r""":math:`g(\theta)=\theta/\left[(1-\theta)\sqrt{K_{eq}}\right] = \sqrt{P_{virtual}}`."""
    if theta >= 1.0:
        return math.inf
    return theta / ((1.0 - theta) * math.sqrt(K_eq))


# -----------------------------------------------------------------------------
# Dictionaries in the shape the solvers expect
# -----------------------------------------------------------------------------

def oxide_props():
    """``oxide_props`` as :mod:`calculations.interface_solver` expects it."""
    return {"D_ox": D_OX, "K_ox": K_OX, "thickness": L_OX, "name": "ExampleOxide"}


def metal_props():
    """``metal_props`` as :mod:`calculations.interface_solver` expects it."""
    return {"D_metal": D_M, "K_s_metal": K_S_M, "thickness": L_M, "name": "ExampleMetal"}


def defect_params(area_fraction=F_DEFECT, defect_type=DEFECT_TYPE):
    """``defect_params`` for the Level 3 parallel-path solvers."""
    return {"area_fraction": area_fraction, "type": defect_type}


def trap_list():
    """``trap_list`` for :func:`calculations.defective_metal.calculate_effective_diffusivity_trapping`."""
    return [dict(TRAP)]


# -----------------------------------------------------------------------------
# Palette
# -----------------------------------------------------------------------------
#
# The production palette lives in calculations/config/model_config.py, chosen against
# a colour-vision validator rather than eyeballed. Importing it requires a resolved
# study -- model_config raises when none is in scope -- and the documentation figures
# must build with no study selected. So it is vendored here and check_equations.py
# asserts the two agree whenever a study IS resolvable, which fails the check rather
# than letting the docs drift to a different set of hues.
#
# Keep in sync with model_config.COLORS.

COLORS = {
    "L1":      "black",      # baseline reference ink
    "L2a":     "#eb6834",    # orange
    "L2b":     "#4a3aa7",    # violet
    "L3":      "#1baf7a",    # aqua
    "L4_gb":   "#008300",    # green
    "L4_trap": "#eda100",    # yellow
    "L4_both": "#2a78d6",    # blue
    "L5":      "#e34948",    # red
}

#: Each level carries a linestyle and marker as well as a hue, so a figure stays
#: readable in greyscale and under colour-vision deficiency. Mirrors
#: model_config.CURVE_STYLES for the keys the documentation actually draws.
CURVE_STYLES = {
    "L1":      {"color": COLORS["L1"],      "ls": "-",  "marker": "o", "lw": 2.5, "ms": 6},
    "L1_L6":   {"color": COLORS["L1"],      "ls": "--", "marker": "s", "lw": 2.5, "ms": 6},
    "L2a":     {"color": COLORS["L2a"],     "ls": "-",  "marker": "o", "lw": 2.5, "ms": 6},
    "L2b":     {"color": COLORS["L2b"],     "ls": "-",  "marker": "o", "lw": 2.5, "ms": 6},
    "L2_L6":   {"color": COLORS["L2b"],     "ls": "--", "marker": "s", "lw": 2.5, "ms": 6},
    "L3":      {"color": COLORS["L3"],      "ls": "-",  "marker": "o", "lw": 2.5, "ms": 6},
    "L4_gb":   {"color": COLORS["L4_gb"],   "ls": "-",  "marker": "^", "lw": 2.5, "ms": 6},
    "L4_trap": {"color": COLORS["L4_trap"], "ls": "-",  "marker": "v", "lw": 2.5, "ms": 6},
    "L4_both": {"color": COLORS["L4_both"], "ls": "-",  "marker": "D", "lw": 2.5, "ms": 6},
    "L5":      {"color": COLORS["L5"],      "ls": "-",  "marker": "o", "lw": 2.5, "ms": 6},
    "L5_L6":   {"color": COLORS["L5"],      "ls": "--", "marker": "s", "lw": 2.5, "ms": 6},
}

#: Text and frame settings. Lighter than the notebooks' ``TEXT_STYLE`` because a
#: documentation figure is read inline at a fraction of a notebook panel's size.
DOC_STYLE = {
    "figsize": (7.0, 4.4),
    "fontsize_title": 11,
    "fontsize_axis": 10,
    "fontsize_tick": 9,
    "fontsize_legend": 9,
    "fontsize_annotation": 8.5,
    "linewidth": 2.0,
    "markersize": 5,
    "grid_alpha": 0.25,
    "dpi": 200,
}


def summary():
    """Every derived quantity of the Example Wall, as a dict. Used by the checks."""
    alpha, beta = permeances()
    P_int = exact_interface_pressure()
    return {
        "alpha": alpha,
        "beta": beta,
        "J_L1": exact_level1_flux(),
        "J_L2a": exact_level2a_flux(),
        "P_int": P_int,
        "J_L2b": exact_level2b_flux(),
        "trapping_term": exact_trapping_term(),
        "D_eff_trap": exact_D_eff_trapping(),
    }


if __name__ == "__main__":
    print(__doc__.split("\n")[0])
    print()
    print(f"  T            {T_K:>14.4g} K")
    print(f"  P_up/P_down  {P_UP:>14.4g} / {P_DOWN:.4g} Pa")
    print(f"  metal        D={D_M:.4g} m2/s  K_s={K_S_M:.4g}  L={L_M:.4g} m")
    print(f"  oxide        D={D_OX:.4g} m2/s  K_ox={K_OX:.4g}  L={L_OX:.4g} m")
    print()
    for key, value in summary().items():
        print(f"  {key:<14} {value:>14.6g}")
