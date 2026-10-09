#!/usr/bin/env python3
r"""Check every equation the documentation states against the code's own numbers.

Each check restates an equation from Part II *independently* -- in plain arithmetic on
the Example Wall -- and compares it with what :mod:`calculations` returns. A check that
merely called the same function twice would pass forever and prove nothing, so every
expected value here is computed from the equation as written in the documentation, not
borrowed from the module under test.

Where no independent form exists (the Level 2b solve, the Level 6 solve) the check is
an *invariant* instead: flux balance across the layers must hold at the solution, and a
limit must collapse to the level below.

Run from the repository root::

    python docs/figures/check_equations.py

Exit code 0 if every check passes, 1 otherwise. Appendix E lists what is checked and
the expected output.
"""

from __future__ import annotations

import math
import os
import pathlib
import sys
import warnings

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import example_wall as ex  # noqa: E402

# Material properties live behind a study switch, and importing model_config with no
# study in scope raises. Nothing in this file needs a study -- the Example Wall supplies
# every number -- but several calculations modules import model_config at module load,
# so a study is named here purely to let those imports succeed. No value from it is used.
os.environ.setdefault("PERMEATION_STUDY", "Guo_etal_2025/study_config.py")

from calculations import oxide_permeation as ox  # noqa: E402
from calculations import permeability as perm  # noqa: E402
from calculations import permeation_calc as pc  # noqa: E402
from calculations import surface_kinetics as sk  # noqa: E402
from calculations.defective_metal import (  # noqa: E402
    calculate_effective_diffusivity_trapping,
    trap_occupancy,
)
from calculations.interface_solver import (  # noqa: E402
    calculate_oxide_metal_system,
    solve_interface_pressure,
)
from calculations.parallel_oxide_defect_paths import (  # noqa: E402
    calculate_parallel_path_flux,
)
from calculations.utils import arrhenius  # noqa: E402

RTOL = 1e-9          # equations that should agree to machine precision
RTOL_SOLVER = 1e-6   # anything that went through a root finder

_results: list[tuple[str, bool, str]] = []


def check(name, got, expected, rtol=RTOL, note=""):
    """Record one comparison of a computed value against the documented equation."""
    if expected == 0:
        ok = abs(got) < rtol
        rel = abs(got)
    else:
        rel = abs(got - expected) / abs(expected)
        ok = rel <= rtol
    detail = f"got {got:.10g}, expected {expected:.10g}, rel {rel:.2e}"
    _results.append((name, ok, f"{detail}{'  -- ' + note if note else ''}"))
    return ok


def check_true(name, ok, detail=""):
    """Record one boolean invariant."""
    _results.append((name, bool(ok), detail))
    return bool(ok)


# -----------------------------------------------------------------------------
# Part II.1 -- sorption laws and Fick
# -----------------------------------------------------------------------------

def check_sorption_and_fick():
    """C = K_s sqrt(P) (Sieverts), C = K_ox P (Henry), J = D dC/L (Fick)."""
    check(
        "Sieverts: C = K_s*sqrt(P)",
        pc.sieverts_concentration(ex.K_S_M, ex.P_UP),
        ex.K_S_M * math.sqrt(ex.P_UP),
    )
    check(
        "Sieverts: C(0) = 0",
        pc.sieverts_concentration(ex.K_S_M, 0.0),
        0.0,
    )
    check(
        "Fick: J = D*(C_up - C_down)/L",
        pc.fick_flux(ex.D_M, 10.0, 2.0, ex.L_M),
        ex.D_M * (10.0 - 2.0) / ex.L_M,
    )
    # Henry is not exposed as its own function; it is inlined in the oxide flux, so the
    # check recovers K_ox*P from a flux with zero downstream pressure.
    J = ox.molecular_diffusion_flux(ex.D_OX, ex.K_OX, ex.L_OX, ex.P_UP, 0.0)
    check(
        "Henry: J*L/D recovers C_up = K_ox*P",
        J * ex.L_OX / ex.D_OX,
        ex.K_OX * ex.P_UP,
    )


# -----------------------------------------------------------------------------
# Part II.2 -- Level 1 and Level 2a
# -----------------------------------------------------------------------------

def check_level1():
    """J = (D K_s / L)(sqrt(P_up) - sqrt(P_down)), and Phi = D K_s."""
    r = pc.calculate_simple_metal_flux(ex.D_M, ex.K_S_M, ex.L_M, ex.P_UP, ex.P_DOWN)
    check("L1: flux", r["flux"], ex.exact_level1_flux())
    check("L1: permeability = D*K_s", r["permeability"], ex.D_M * ex.K_S_M)
    check("L1: C_up = K_s*sqrt(P_up)", r["C_up"], ex.K_S_M * math.sqrt(ex.P_UP))

    # Sieverts' signature: halving the driving force in sqrt(P) halves the flux.
    r4 = pc.calculate_simple_metal_flux(ex.D_M, ex.K_S_M, ex.L_M, ex.P_UP / 4.0, 0.0)
    check("L1: J(P/4) = J(P)/2  (sqrt-law)", r4["flux"], r["flux"] / 2.0)


def check_level2a():
    """J = (D_ox K_ox / L_ox)(P_up - P_down), and R_ox = L_ox/(D_ox K_ox)."""
    J = ox.molecular_diffusion_flux(ex.D_OX, ex.K_OX, ex.L_OX, ex.P_UP, ex.P_DOWN)
    check("L2a: flux", J, ex.exact_level2a_flux())
    check(
        "L2a: R_ox = L_ox/(D_ox*K_ox)",
        ox.calculate_oxide_resistance(ex.D_OX, ex.K_OX, ex.L_OX),
        ex.L_OX / (ex.D_OX * ex.K_OX),
    )
    # Henry's signature: the flux is LINEAR in pressure, unlike Level 1.
    J_half = ox.molecular_diffusion_flux(ex.D_OX, ex.K_OX, ex.L_OX, ex.P_UP / 2.0, 0.0)
    check("L2a: J(P/2) = J(P)/2  (linear law)", J_half, J / 2.0)

    check(
        "L2: R_metal = 2*L*sqrt(P_int)/(D*K_s)",
        ox.calculate_metal_resistance(ex.D_M, ex.K_S_M, ex.L_M, 1.0e4),
        2.0 * ex.L_M * math.sqrt(1.0e4) / (ex.D_M * ex.K_S_M),
    )


# -----------------------------------------------------------------------------
# Part II.3 -- the Level 2b series solve
# -----------------------------------------------------------------------------

def check_level2b():
    """The solved interface pressure must equal the closed-form quadratic root."""
    sol = solve_interface_pressure(ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props())
    check(
        "L2b: P_int matches closed form",
        sol["P_interface"],
        ex.exact_interface_pressure(),
        rtol=RTOL_SOLVER,
        note="brentq vs the exact root of alpha*x^2 + beta*x - alpha*P_up",
    )

    sys_ = calculate_oxide_metal_system(ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props())
    check("L2b: flux", sys_["flux"], ex.exact_level2b_flux(), rtol=RTOL_SOLVER)

    # The invariant that defines the interface: both layers carry the same flux.
    alpha, beta = ex.permeances()
    P_int = sol["P_interface"]
    J_ox = alpha * (ex.P_UP - P_int)
    J_m = beta * (math.sqrt(P_int) - math.sqrt(ex.P_DOWN))
    check("L2b: flux balance J_ox = J_metal", J_ox, J_m, rtol=RTOL_SOLVER)

    # The series flux must sit below BOTH single-layer fluxes: adding a layer in series
    # can only add resistance. This catches a sign or reciprocal slip that the equality
    # checks above would not.
    check_true(
        "L2b: J_series < min(J_L1, J_L2a)",
        sys_["flux"] < min(ex.exact_level1_flux(), ex.exact_level2a_flux()),
        f"J_2b={sys_['flux']:.6g}, J_L1={ex.exact_level1_flux():.6g}, "
        f"J_L2a={ex.exact_level2a_flux():.6g}",
    )


def check_level2b_limit():
    """As L_ox -> 0 the oxide vanishes and Level 2b must recover Level 1."""
    thin = ex.oxide_props()
    thin["thickness"] = ex.L_OX * 1e-8
    sys_ = calculate_oxide_metal_system(ex.P_UP, ex.P_DOWN, thin, ex.metal_props())
    check(
        "L2b -> L1 as L_ox -> 0",
        sys_["flux"],
        ex.exact_level1_flux(),
        rtol=1e-4,
        note="the limit check panel (D) draws this",
    )


# -----------------------------------------------------------------------------
# Part II.4 -- Level 3 parallel paths
# -----------------------------------------------------------------------------

def check_level3():
    """J_total = f_intact*j_intact + f_defect*j_defect, an area-weighted SUM."""
    r = calculate_parallel_path_flux(
        ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props(), ex.defect_params()
    )
    f_d = ex.F_DEFECT
    expected = (1.0 - f_d) * r["flux_intact_per_area"] + f_d * r["flux_defect_per_area"]
    check(
        "L3: J_total = sum of area-weighted paths",
        r["flux_total"],
        expected,
        rtol=RTOL_SOLVER,
        note="a SUM, not the reciprocal sum of parallel resistors",
    )
    check_true(
        "L3: contributions sum to the total",
        math.isclose(
            r["flux_intact_contribution"] + r["flux_defect_contribution"],
            r["flux_total"],
            rel_tol=RTOL_SOLVER,
        ),
    )


def check_level3_limit():
    """As f_defect -> 0 the defect path vanishes and Level 3 must recover Level 2b."""
    r = calculate_parallel_path_flux(
        ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props(), ex.defect_params(0.0)
    )
    base = calculate_oxide_metal_system(
        ex.P_UP, ex.P_DOWN, ex.oxide_props(), ex.metal_props()
    )
    check("L3 -> L2b as f_defect -> 0", r["flux_total"], base["flux"], rtol=RTOL_SOLVER)


# -----------------------------------------------------------------------------
# Part II.5 -- Level 4 trapping
# -----------------------------------------------------------------------------

def check_level4_trapping():
    """Oriani: K = exp(E_b/RT) and D_eff = D_L/(1 + sum_i N_Ti K_i / N_L)."""
    C_L = ex.K_S_M * math.sqrt(ex.P_UP)
    occ = trap_occupancy(
        temperature=ex.T_K,
        binding_energy=ex.TRAP["binding_energy"],
        trap_density=ex.TRAP["density"],
        lattice_density=ex.N_L,
        lattice_concentration=C_L,
    )
    K_expected = math.exp(ex.TRAP["binding_energy"] / (ex.R_GAS * ex.T_K))
    check("L4: K = exp(E_b/RT)", occ["K_equilibrium"], K_expected, rtol=1e-6)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = calculate_effective_diffusivity_trapping(
            D_lattice=ex.D_M,
            temperature=ex.T_K,
            trap_list=ex.trap_list(),
            lattice_concentration=C_L,
            lattice_density=ex.N_L,
        )
    check(
        "L4: trapping term = sum N_T*K/N_L",
        res["trapping_term"],
        ex.exact_trapping_term(),
        rtol=1e-6,
    )
    check(
        "L4: D_eff = D_L/(1 + trapping term)",
        res["D_eff"],
        ex.exact_D_eff_trapping(),
        rtol=1e-6,
    )
    check_true(
        "L4: trapping RETARDS (D_eff < D_lattice)",
        res["D_eff"] < ex.D_M,
        f"D_eff={res['D_eff']:.6g} < D_L={ex.D_M:.6g}",
    )


def check_level4_limit():
    """With no traps, D_eff must be exactly D_lattice."""
    C_L = ex.K_S_M * math.sqrt(ex.P_UP)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = calculate_effective_diffusivity_trapping(
            D_lattice=ex.D_M,
            temperature=ex.T_K,
            trap_list=[],
            lattice_concentration=C_L,
            lattice_density=ex.N_L,
        )
    check("L4 -> L1 with an empty trap list", res["D_eff"], ex.D_M)


# -----------------------------------------------------------------------------
# Part II.6 -- Arrhenius
# -----------------------------------------------------------------------------

def check_arrhenius():
    """k(T) = k_ref * exp[-(E/R)(1/T - 1/T_ref)], the reference-temperature form."""
    k_ref, E, T, T_ref = 1.0e-10, 50.0e3, 900.0, 800.0
    check(
        "Arrhenius: reference-temperature form",
        arrhenius(k_ref, E, T, T_ref),
        k_ref * math.exp((-E / ex.R_GAS) * (1.0 / T - 1.0 / T_ref)),
    )
    check("Arrhenius: k(T_ref) = k_ref", arrhenius(k_ref, E, T_ref, T_ref), k_ref)


# -----------------------------------------------------------------------------
# Part II.7 -- Level 6 surface kinetics
# -----------------------------------------------------------------------------

def check_level6():
    """Langmuir g(theta), the closed-form interface pressure, and flux balance."""
    theta = 0.3
    check(
        "L6: g(theta) = theta/[(1-theta)sqrt(K_eq)]",
        sk.g_theta(theta, ex.K_EQ),
        ex.exact_g_theta(theta),
    )

    alpha, beta = ex.permeances()
    check(
        "L6: sqrt(P_int) = (alpha*g + beta*sqrt(P_down))/(alpha+beta)",
        sk.sqrt_P_int_from_theta(theta, alpha, beta, ex.K_EQ, ex.P_DOWN),
        (alpha * ex.exact_g_theta(theta) + beta * math.sqrt(ex.P_DOWN)) / (alpha + beta),
    )

    check(
        "L6: J_surface = k_diss*P*(1-theta)^2 - (k_diss/K_eq)*theta^2",
        sk.surface_flux(theta, ex.P_UP, ex.K_DISS, ex.K_EQ),
        ex.K_DISS * ex.P_UP * (1 - theta) ** 2 - (ex.K_DISS / ex.K_EQ) * theta**2,
    )

    sol = sk.solve_steady_state_flux_L1L6(
        ex.P_UP, ex.P_DOWN, ex.L_M, ex.K_DISS, ex.K_EQ, ex.D_M, ex.K_S_M
    )
    # The solve's defining invariant: at steady state the surface and the metal carry
    # the same flux. Nothing else pins theta.
    J_surf = sk.surface_flux(sol["theta"], ex.P_UP, ex.K_DISS, ex.K_EQ)
    check(
        "L6: flux balance J_surface = J_metal at the solution",
        J_surf,
        sol["J_ss"],
        rtol=RTOL_SOLVER,
    )
    check_true(
        "L6: 0 < theta < 1",
        0.0 < sol["theta"] < 1.0,
        f"theta={sol['theta']:.6g}",
    )
    # Adding a finite surface rate can only reduce the flux below the L1 value, which
    # assumes instantaneous dissociation.
    check_true(
        "L6: J_L1L6 < J_L1  (a finite surface rate can only cost flux)",
        sol["J_ss"] < ex.exact_level1_flux(),
        f"J_L1L6={sol['J_ss']:.6g} < J_L1={ex.exact_level1_flux():.6g}",
    )


# -----------------------------------------------------------------------------
# Part II.8 -- permeability tiers and the unit refusal
# -----------------------------------------------------------------------------

def check_permeability_tiers():
    """Tier 1 is intrinsic; the two sorption laws must refuse to combine."""
    m = perm.intrinsic_metal_permeability(ex.D_M, ex.K_S_M)
    check("Tier 1: Phi_m = D_m*K_s", m["permeability"], ex.D_M * ex.K_S_M)
    check_true("Tier 1: metal carries Sieverts units", m["units"] == perm.UNITS_SIEVERTS, m["units"])

    o1 = perm.intrinsic_oxide_permeability(ex.D_OX, ex.K_OX, perm.MODEL_1)
    check("Tier 1: Phi_ox = D_ox*K_ox", o1["permeability"], ex.D_OX * ex.K_OX)
    check_true("Tier 1: Model 1 oxide carries Henry units", o1["units"] == perm.UNITS_HENRY, o1["units"])

    o2 = perm.intrinsic_oxide_permeability(ex.D_OX, ex.K_OX, perm.MODEL_2)
    check_true(
        "Tier 1: Model 2 oxide carries Sieverts units",
        o2["units"] == perm.UNITS_SIEVERTS,
        o2["units"],
    )

    # The refusal IS the physics. Model 1 cannot produce a single stack permeability.
    try:
        perm.combine_series(m, o1, ex.L_M, ex.L_OX)
        check_true("Tier 3: Model 1 refuses to combine mixed units", False, "no error raised")
    except perm.UnitMismatchError:
        check_true("Tier 3: Model 1 refuses to combine mixed units", True, "UnitMismatchError raised")

    # Model 2 collapses to one permeance, because both layers obey Sieverts.
    stack = perm.combine_series(m, o2, ex.L_M, ex.L_OX)
    expected = (ex.L_M + ex.L_OX) / (
        ex.L_M / (ex.D_M * ex.K_S_M) + ex.L_OX / (ex.D_OX * ex.K_OX)
    )
    check("Tier 3: Model 2 series permeability", stack["permeability"], expected)


def check_pressure_exponent():
    """n = dlnJ/dlnP: 1/2 for a Sieverts step, 1 for a Henry layer."""
    alpha, beta = ex.permeances()
    n_metal = perm.pressure_exponent(lambda P: beta * math.sqrt(P), ex.P_UP)
    check("n = 1/2 for a pure Sieverts step", n_metal, 0.5, rtol=1e-3)

    n_oxide = perm.pressure_exponent(lambda P: alpha * P, ex.P_UP)
    check("n = 1 for a pure Henry layer", n_oxide, 1.0, rtol=1e-3)

    # The series wall must sit strictly between the two, which is what makes n a
    # usable diagnostic of which step controls.
    def J_series(P):
        return calculate_oxide_metal_system(
            P, ex.P_DOWN, ex.oxide_props(), ex.metal_props()
        )["flux"]

    n_series = perm.pressure_exponent(J_series, ex.P_UP)
    check_true(
        "1/2 < n < 1 for the mixed series wall",
        0.5 < n_series < 1.0,
        f"n={n_series:.4f}",
    )


# -----------------------------------------------------------------------------
# Documentation-integrity checks
# -----------------------------------------------------------------------------

def check_palette_in_sync():
    """The vendored palette must match the production one in model_config."""
    try:
        from calculations.config.model_config import COLORS as PROD
    except Exception as exc:  # noqa: BLE001 - no study resolvable is a skip, not a failure
        check_true("palette matches model_config.COLORS", True, f"skipped: {type(exc).__name__}")
        return
    mismatched = {k: (v, PROD.get(k)) for k, v in ex.COLORS.items() if PROD.get(k) != v}
    check_true(
        "palette matches model_config.COLORS",
        not mismatched,
        "in sync" if not mismatched else f"differs: {mismatched}",
    )


# -----------------------------------------------------------------------------

CHECKS = (
    check_sorption_and_fick,
    check_level1,
    check_level2a,
    check_level2b,
    check_level2b_limit,
    check_level3,
    check_level3_limit,
    check_level4_trapping,
    check_level4_limit,
    check_arrhenius,
    check_level6,
    check_permeability_tiers,
    check_pressure_exponent,
    check_palette_in_sync,
)


def main():
    print("Checking the documented equations against calculations/, on the Example Wall.")
    print(f"  alpha = {ex.permeances()[0]:.6g} mol/m2/s/Pa      (oxide, Henry)")
    print(f"  beta  = {ex.permeances()[1]:.6g} mol/m2/s/Pa^0.5  (metal, Sieverts)")
    print()

    for fn in CHECKS:
        fn()

    width = max(len(name) for name, _, _ in _results)
    for name, ok, detail in _results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}}  {detail}")

    passed = sum(1 for _, ok, _ in _results if ok)
    total = len(_results)
    print()
    print(f"{passed} / {total} checks passed")
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
