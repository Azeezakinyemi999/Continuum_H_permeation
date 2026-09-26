"""Physical invariants of the two model wrappers.

``level5_model_wrapper`` and ``level5L6_model_wrapper`` are what the entire sensitivity
analysis calls, and until now they were pinned only by the thirteen values in the
container gate. Pinned values catch a change; they do not say which changes would be
*wrong*. These tests state the physics instead, so they stay meaningful when the study,
the materials or the operating point change.

Each invariant here is one the model must satisfy for any study, not a property of the
316L numbers that happen to be active. The exception is the Sieverts-exponent test,
which asserts a half-power only because the active study is metal-limited — it checks
the regime first and says so.

Deliberately NOT tested: that L5L6 collapses onto L5 as the surface becomes fast. The
two are different models rather than refinements of one another — L5 puts intact H2
through the oxide under Henry's law and L5L6 puts atomic H through it under Sieverts —
so they do not converge, and asserting they should would pin a bug in place. See
docs/theory/two-models.md.
"""
import math

import pytest

from calculations.config.model_config import (
    DEFAULT_PARAMS_LEVEL5, DEFAULT_PARAMS_LEVEL5L6,
)
from calculations.sensitivity import level5_model_wrapper, level5L6_model_wrapper

# (label, wrapper, its default parameters). Parametrising on this pair keeps every
# invariant that applies to both models written exactly once.
MODELS = [
    ("L5", level5_model_wrapper, DEFAULT_PARAMS_LEVEL5),
    ("L5L6", level5L6_model_wrapper, DEFAULT_PARAMS_LEVEL5L6),
]
IDS = [m[0] for m in MODELS]


def run(wrapper, params, **overrides):
    """The wrapper's result with ``overrides`` applied to a copy of ``params``."""
    p = dict(params)
    p.update(overrides)
    return wrapper(p)


def flux(wrapper, params, **overrides):
    return run(wrapper, params, **overrides)["flux"]


@pytest.mark.parametrize("label,wrapper,params", MODELS, ids=IDS)
def test_nominal_flux_is_positive_and_finite(label, wrapper, params):
    """The baseline every other test perturbs from."""
    J = flux(wrapper, params)
    assert math.isfinite(J), f"{label} returned {J} at its own defaults"
    assert J > 0


@pytest.mark.parametrize("label,wrapper,params", MODELS, ids=IDS)
def test_flux_rises_with_upstream_pressure(label, wrapper, params):
    """Raising the driving force cannot reduce the flux."""
    fluxes = [flux(wrapper, params, P_upstream=P) for P in (1e3, 1e4, 1e5, 1e6)]
    assert all(a < b for a, b in zip(fluxes, fluxes[1:])), f"{label}: {fluxes}"


@pytest.mark.parametrize("label,wrapper,params", MODELS, ids=IDS)
def test_flux_falls_with_metal_thickness(label, wrapper, params):
    """More wall is more resistance, whatever else the wall is doing."""
    fluxes = [flux(wrapper, params, metal_thickness=L) for L in (5e-4, 1e-3, 2e-3)]
    assert all(a > b for a, b in zip(fluxes, fluxes[1:])), f"{label}: {fluxes}"


@pytest.mark.parametrize("label,wrapper,params", MODELS, ids=IDS)
def test_flux_falls_with_oxide_thickness(label, wrapper, params):
    """A thicker barrier is a better barrier.

    Worth stating even though it looks obvious: the oxide thickness reaches the
    solvers from ``CONDITIONS``, and a regression that disconnected it would leave
    this flat rather than raise.
    """
    fluxes = [flux(wrapper, params, oxide_thickness=L) for L in (1e-7, 1e-6, 1e-5)]
    assert all(a > b for a, b in zip(fluxes, fluxes[1:])), f"{label}: {fluxes}"


@pytest.mark.parametrize("label,wrapper,params", MODELS, ids=IDS)
def test_flux_rises_with_temperature(label, wrapper, params):
    """Every activated process in the stack has a positive activation energy."""
    fluxes = [flux(wrapper, params, temperature=T) for T in (673, 773, 873, 973)]
    assert all(a < b for a, b in zip(fluxes, fluxes[1:])), f"{label}: {fluxes}"


@pytest.mark.parametrize("label,wrapper,params", MODELS, ids=IDS)
def test_resistance_fractions_partition_the_whole(label, wrapper, params):
    """The frac_* keys are a decomposition, so they must sum to one.

    This is what makes `regime` meaningful: it names the largest share, and a share of
    something that does not add to 1 is not a share.
    """
    result = run(wrapper, params)
    fracs = {k: v for k, v in result.items() if k.startswith("frac_")}
    assert fracs, f"{label} reports no frac_* keys"
    assert sum(fracs.values()) == pytest.approx(1.0, abs=1e-9), fracs
    assert all(0.0 <= v <= 1.0 for v in fracs.values()), fracs


@pytest.mark.parametrize("label,wrapper,params", MODELS, ids=IDS)
def test_regime_names_the_largest_resistance(label, wrapper, params):
    """`regime` is read straight off these plots, so it must agree with the numbers."""
    result = run(wrapper, params)
    fracs = {k[len("frac_"):]: v for k, v in result.items() if k.startswith("frac_")}
    assert result["regime"] == max(fracs, key=fracs.get), (result["regime"], fracs)


def test_zero_driving_force_gives_zero_flux():
    """No gradient, no transport — the cleanest check there is."""
    assert flux(level5_model_wrapper, DEFAULT_PARAMS_LEVEL5,
                P_upstream=DEFAULT_PARAMS_LEVEL5["P_downstream"]) == 0.0


@pytest.mark.xfail(strict=True, reason=(
    "L5L6 returns nan rather than 0 when P_upstream == P_downstream. The surface "
    "solve cannot bracket a root at zero driving force. Recorded in CHECKLIST.md; "
    "flip this to a plain assertion when it is fixed."))
def test_zero_driving_force_gives_zero_flux_L5L6():
    assert flux(level5L6_model_wrapper, DEFAULT_PARAMS_LEVEL5L6,
                P_upstream=DEFAULT_PARAMS_LEVEL5L6["P_downstream"]) == 0.0


# A wall with a negligible oxide and no defects is metal-limited by construction, for
# any study. Forcing the limit rather than detecting it is deliberate: an earlier
# version of the test below skipped itself unless the nominal point happened to be
# metal-limited, and every mutation that broke Sieverts' law also moved the regime — so
# the guard silently skipped the test instead of failing it.
BARE_METAL = dict(oxide_thickness=1e-9, f_pinhole=0.0, f_crack=0.0, f_gb_defect=0.0)


def test_the_forced_bare_metal_limit_really_is_metal_limited():
    """The premise of the Sieverts test below, asserted rather than assumed."""
    result = run(level5_model_wrapper, DEFAULT_PARAMS_LEVEL5, **BARE_METAL)
    assert result["regime"] == "metal"
    assert result["frac_metal"] > 0.99


def test_metal_limited_flux_follows_a_square_root_in_pressure():
    """Sieverts' law, visible end to end through the wrapper.

    Hydrogen arrives as H2 but dissolves as atoms, so a metal-limited wall carries
    J proportional to sqrt(P): the log-log slope against pressure is 1/2. A Henry law
    would give 1, which is what a dropped square root anywhere in the metal path looks
    like from here.
    """
    low = flux(level5_model_wrapper, DEFAULT_PARAMS_LEVEL5, P_upstream=1e4, **BARE_METAL)
    high = flux(level5_model_wrapper, DEFAULT_PARAMS_LEVEL5, P_upstream=4e4, **BARE_METAL)
    exponent = math.log(high / low) / math.log(4.0)
    assert exponent == pytest.approx(0.5, abs=0.01)


def test_prf_is_at_least_one():
    """A coating cannot make the wall more permeable than the bare metal.

    Level 5 only: `PRF` is nan from the Level 5L6 wrapper at default parameters, which
    is a known defect rather than a value to assert on.
    """
    result = run(level5_model_wrapper, DEFAULT_PARAMS_LEVEL5)
    assert result["PRF"] >= 1.0


def test_removing_the_defects_cannot_increase_the_flux():
    """Pinholes, cracks and grain boundaries are short circuits through the barrier.

    Closing them is the one change guaranteed to help, so it is the sharpest test of
    whether the Level 3 paths are wired in at all: if the defect fractions were being
    dropped somewhere, both fluxes would come back equal.
    """
    leaky = flux(level5_model_wrapper, DEFAULT_PARAMS_LEVEL5)
    sealed = flux(level5_model_wrapper, DEFAULT_PARAMS_LEVEL5,
                  f_pinhole=0.0, f_crack=0.0, f_gb_defect=0.0)
    assert sealed < leaky
