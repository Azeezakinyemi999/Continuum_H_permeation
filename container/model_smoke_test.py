#!/usr/bin/env python3
"""Gate 2 — does the scientific application still compute the right numbers?

Environment reproducibility (gate 1) does not imply scientific reproducibility.
This asserts pinned reference values through the real solvers.

**Study-aware.** The model computes whichever study `ACTIVE_STUDY` names, so a
single set of pinned numbers could only ever test one of them. REFERENCE below
holds a measured set per study and the test selects on ACTIVE_STUDY; a study with
no entry fails with a message saying to measure one rather than silently passing.

Model INPUTS are derived from the active study too — its own metal/oxide
properties at T, its own surface kinetics — so the test exercises the study it
claims to. That matters: the earlier fixed-input version paired the metal's
k_diss with the OXIDE's K_eq (1e-4) rather than the metal's own K_eq_metal_ref
(1e-3), which is why the L1L6 and microstructure numbers here differ from the
values pinned before. The model is unchanged — re-running it with the old
hardcoded inputs still returns the old numbers exactly.

Reference values were measured on Python 3.9.23 / macOS / x86_64 (mace_env) and
confirmed bit-identical on Python 3.12.14 / linux-amd64 from
container/requirements.lock.txt. See PACKAGING_1.0.0.md.

Comparison uses rtol=1e-9, not equality: macOS libm and glibc round np.log10
differently in the last bit (observed 2.1e-16 relative on a plot colour channel),
and exact comparison would make this flaky for no scientific reason. A
disagreement larger than rtol=1e-9 is a real finding and blocks the release.

Run from the workspace root so `calculations` is importable:
    docker run --rm -v "$PWD/workspace:/workspace" -w /workspace \
        hydrogen-model:1.0.0 python /opt/model_smoke_test.py

Exit 0 = pass, 1 = fail.
"""
import math
import os
import sys
import warnings

# Running `python /opt/model_smoke_test.py` puts /opt on sys.path[0], not the
# cwd, so `import calculations` would fail even with -w /workspace. The image
# sets PYTHONPATH=/workspace; this makes the script work regardless.
for _cand in (os.getcwd(), "/workspace"):
    if _cand not in sys.path and os.path.isdir(os.path.join(_cand, "calculations")):
        sys.path.insert(0, _cand)

RTOL = 1e-9

# Fixed operating point. T=700 K is deliberate: it sits inside the Cr2O3
# validated range [473, 773] K, so no extrapolation warning perturbs the run.
T, P_UP, P_DOWN = 700.0, 1.0e5, 1.0e2
L_METAL = 1.0e-3
# L_OXIDE is NOT fixed here: it comes from the active study's oxide, so the test
# follows the study rather than asserting one study's scale onto another.

REFERENCE = {
    'incoloy802_cr2o3': {
        'metal_name':             'metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985',
        'oxide_name':             'Cr2O3_sample4',
        'L1_flux':                1.955109516956673e-07,
        'L2b_flux':               1.888003770552079e-07,
        'L2b_P_interface':        93462.901558976,
        'L2b_resistance_ratio':   0.036154225782795736,
        'L1L6_J_ss':              1.9549321314132636e-07,
        'L1L6_P_int':             99982.42875170171,
        'L1L6_theta':             0.9090836473455054,
        'L1L6_rate_limiting':     'metal',
        'micro_D_eff':            8.720200643474982e-12,
        'micro_overall_factor':   0.4098345562673666,
        'micro_theta_total':      0.09282336531113745,
    },
    'fuerst_etal_2024_model_config': {
        'metal_name':             'Hastelloy_N_fuerst_2024',
        'oxide_name':             'Cr2O3_sample4',
        'L1_flux':                2.074929655645745e-06,
        'L2b_flux':               1.4462213335730932e-06,
        'L2b_P_interface':        49925.464278538966,
        'L2b_resistance_ratio':   0.5249886661138037,
        'L1L6_J_ss':              2.0729356704106247e-06,
        'L1L6_P_int':             99813.96656812982,
        'L1L6_theta':             0.9090139348659355,
        'L1L6_rate_limiting':     'metal',
        'micro_D_eff':            8.008430918590877e-11,
        'micro_overall_factor':   0.4098345562673666,
        'micro_theta_total':      0.09282336531113745,
    },
    'Guo_etal_2025_316L': {
        'metal_name':             '316L_Guo_2025',
        'oxide_name':             'Cr2O3_sample4',
        'L1_flux':                5.118412424120975e-07,
        'L2b_flux':               4.672050860445382e-07,
        'L2b_P_interface':        83823.30752058215,
        'L2b_resistance_ratio':   0.09994486510402441,
        'L1L6_J_ss':              1.276041708643705e-07,
        'L1L6_P_int':             7455.263517774021,
        'L1L6_theta':             0.7709733288080891,
        'L1L6_rate_limiting':     'metal',
        'micro_D_eff':            1.6434267493955927e-11,
        'micro_overall_factor':   0.40983455626736665,
        'micro_theta_total':      0.09282336531113745,
    },
}

failures: list[str] = []


def close(label, got, want, rtol=RTOL):
    ok = isinstance(got, (int, float)) and math.isclose(got, want, rel_tol=rtol)
    if ok:
        print(f"OK    {label:<34} {got!r}")
    else:
        rel = "n/a"
        if isinstance(got, (int, float)) and want:
            rel = f"{abs(got - want) / abs(want):.3e}"
        print(f"FAIL  {label:<34} got {got!r}  want {want!r}  rel={rel}")
        failures.append(label)
    return ok


def equal(label, got, want):
    ok = got == want
    print(f"{'OK  ' if ok else 'FAIL'}  {label:<34} {got!r}"
          f"{'' if ok else f'  want {want!r}'}")
    if not ok:
        failures.append(label)
    return ok


print("=" * 78)
print("HYDROGEN MODEL — SCIENTIFIC SMOKE TEST (gate 2 of 2)")
print("=" * 78)

try:
    from calculations.classify_regime import classify_regime_level14
    from calculations.config.model_config import (
        ACTIVE_STUDY, METALS, MICROSTRUCTURE, OXIDES, build_simulation_config,
    )
    from calculations.defective_metal import combined_microstructure_model
    from calculations.interface_solver import calculate_oxide_metal_system
    from calculations.oxide_permeation import (
        get_metal_properties_at_T, get_oxide_properties_at_T,
    )
    from calculations.permeation_calc import calculate_simple_metal_flux
    from calculations.surface_kinetics import solve_steady_state_flux_L1L6
    from calculations.utils import arrhenius
except ImportError as exc:
    print(f"FAIL  cannot import calculations: {exc}")
    print("      Run with the workspace mounted and cwd=/workspace, e.g.")
    print('      docker run -v "$PWD/workspace:/workspace" -w /workspace ...')
    sys.exit(1)

warnings.simplefilter("ignore")

print(f"\nactive study: {ACTIVE_STUDY}")
if ACTIVE_STUDY not in REFERENCE:
    print(f"FAIL  no reference values for study {ACTIVE_STUDY!r}.")
    print(f"      known: {sorted(REFERENCE)}")
    print( "      Measure a set for it and add it to REFERENCE — passing without")
    print( "      one would mean the gate checks nothing about this study.")
    sys.exit(1)
REF = REFERENCE[ACTIVE_STUDY]
print(f"OK    reference set found     {len(REF)} pinned values")

SIM = build_simulation_config()
equal("metal resolves in METALS", SIM["metal_name"], REF["metal_name"])
equal("oxide resolves in OXIDES", SIM["oxide_name"], REF["oxide_name"])

# ---------------------------------------------------------------- analytic
print("\nanalytic:")
close("arrhenius", arrhenius(1.0e-11, 50_000.0, 800.0, 700.0),
      2.926830409475082e-11)

# ------------------------------------------------------- L1: perfect metal
print("\nLevel 1 (perfect metal):")
# NB: the T-evaluated getters return D_metal / K_s_metal, not D / K_s.
L_OXIDE = OXIDES[SIM["oxide_name"]]["thickness"]
mp = dict(get_metal_properties_at_T(SIM["metal_name"], T), thickness=L_METAL)
op = dict(get_oxide_properties_at_T(SIM["oxide_name"], T), thickness=L_OXIDE)
r1 = calculate_simple_metal_flux(mp["D_metal"], mp["K_s_metal"],
                                 L_METAL, P_UP, P_DOWN)
close("L1 flux", r1["flux"], REF["L1_flux"])

# --------------------------------------------- L2b: oxide + metal (brentq)
print("\nLevel 2b (oxide+metal, brentq):")
# calculate_oxide_metal_system needs `thickness` in BOTH props dicts; the
# T-evaluated getter does not supply it for the metal.
r2 = calculate_oxide_metal_system(P_UP, P_DOWN, op, mp, T_K=T)
close("L2b flux", r2["flux"], REF["L2b_flux"])
close("L2b P_interface", r2["P_interface"], REF["L2b_P_interface"])
close("L2b resistance_ratio", r2["resistance_ratio"], REF["L2b_resistance_ratio"])
if not (isinstance(r2.get("flux_error"), float) and abs(r2["flux_error"]) < 1e-12):
    print(f"FAIL  L2b solver converged           flux_error={r2.get('flux_error')!r}")
    failures.append("L2b solver converged")
else:
    print(f"OK    {'L2b solver converged':<34} flux_error={r2['flux_error']:.3e}")

# ------------------------------------- L1L6: surface kinetics (brentq #2)
print("\nLevel 1+6 (surface kinetics, brentq):")
# Returns J_ss / P_int / theta / beta / rate_limiting — there is no 'flux' key.
# Inputs from the study's OWN metal surface kinetics, not literals -- pairing one
# study's k_diss with another's K_eq is exactly the mismatch this replaced.
_sk = METALS[SIM["metal_name"]]["surface_kinetics"]
r3 = solve_steady_state_flux_L1L6(
    P_up=P_UP, P_down=P_DOWN, L_m=L_METAL,
    k_diss=_sk["k_diss_metal_ref"], K_eq=_sk["K_eq_metal_ref"],
    D_m=mp["D_metal"], K_s_m=mp["K_s_metal"])
close("L1L6 J_ss", r3["J_ss"], REF["L1L6_J_ss"])
close("L1L6 P_int", r3["P_int"], REF["L1L6_P_int"])
close("L1L6 theta", r3["theta"], REF["L1L6_theta"])
equal("L1L6 rate_limiting", r3["rate_limiting"], REF["L1L6_rate_limiting"])

# ------------------------------------------- microstructure (4 trap types)
print("\nLevel 4 (microstructure, 4 traps):")
m = combined_microstructure_model(
    mp["D_metal"], T, dict(MICROSTRUCTURE),
    lattice_concentration=10.0, lattice_density=MICROSTRUCTURE["lattice_density"])
close("micro D_eff", m["D_eff"], REF["micro_D_eff"])
close("micro overall_factor", m["overall_factor"], REF["micro_overall_factor"])
close("micro theta_total", m["trapping"]["theta_total"], REF["micro_theta_total"])

# ------------------------------------------------- regime classification
print("\nregime classification:")
equal("classify(0.2)", classify_regime_level14(0.2)["regime_hierarchy"],
      "metal_limited/traps_defect_limited")
equal("classify(0.9)", classify_regime_level14(0.9)["regime_hierarchy"],
      "metal_limited/lattice_limited")

# -------------------------------------------------- SALib (seeded, exact)
print("\nSALib (seeded Latin hypercube -> PAWN + Borgonovo delta):")
import numpy as np
from SALib.analyze import delta as delta_an, pawn
from SALib.sample import latin

problem = {"num_vars": 3, "names": ["a", "b", "c"],
           "bounds": [[1.0, 10.0], [0.1, 1.0], [100.0, 1000.0]]}
X = latin.sample(problem, 256, seed=12345)
Y = np.log10(X[:, 0]) * X[:, 1] + np.sqrt(X[:, 2])
close("SALib X sum", float(X.sum()), 142364.52517734436)
close("SALib Y mean", float(Y.mean()), 23.05501556524375)
for i, want in enumerate([0.17518028846153846, 0.14708533653846154, 0.6875]):
    close(f"pawn_median[{i}]",
          float(pawn.analyze(problem, X, Y, S=10, seed=12345,
                             print_to_console=False)["median"][i]), want)
_dl = delta_an.analyze(problem, X, Y, num_resamples=10, seed=12345,
                       print_to_console=False)
for i, want in enumerate([0.10706087473576083, 0.06905348258240432,
                          0.6947500046510651]):
    close(f"delta[{i}]", float(_dl["delta"][i]), want)

# ---------------------------------- real SA artefacts (needs the workspace)
print("\nshipped SA artefacts:")
R = os.path.join("Application", "sa_results")
needed = ["routeB_givendata.csv", "master_clusters.csv", "compare_delta_flux.csv"]
if not all(os.path.exists(os.path.join(R, f)) for f in needed):
    # Expected for the shipped workspace: the SA result CSVs are deliberately
    # not distributed (see PACKAGING_1.0.0.md), so this section only runs in the
    # development tree. It is not a failure.
    print(f"SKIP  {R}/ CSVs absent — not shipped; run the sensitivity notebooks first")
else:
    import pandas as pd
    from calculations.sensitivity import (
        parallel_coordinates_samples, parallel_coordinates_sensitivity,
        top_drivers,
    )
    rb = pd.read_csv(os.path.join(R, "routeB_givendata.csv"))
    mc = pd.read_csv(os.path.join(R, "master_clusters.csv"))
    cf = pd.read_csv(os.path.join(R, "compare_delta_flux.csv"))
    equal("routeB shape", rb.shape, (324, 7))
    equal("master_clusters shape", mc.shape, (4508, 51))
    equal("top_drivers[metal]", list(top_drivers(rb, "metal", "flux", k=6)),
          ["temperature", "P_upstream", "f_pinhole", "H_sol_ox",
           "lattice_density", "trap_vacancy_N_T"])
    equal("top_drivers[surface]", list(top_drivers(rb, "surface", "flux", k=6)),
          ["temperature", "P_upstream", "k_diss_ref", "f_pinhole",
           "K_s_ref", "D_ref"])
    # Figure HTML length is NOT asserted: it varies with the embedded plotly
    # version string. Only that a figure builds at all.
    n1 = len(parallel_coordinates_sensitivity(cf, title="C").to_html(
        include_plotlyjs="cdn"))
    dims = list(top_drivers(rb, "metal", "flux", k=5)) + ["flux"]
    n2 = len(parallel_coordinates_samples(mc[mc["regime"] == "metal"], dims,
                                          color_by="flux").to_html(
        include_plotlyjs="cdn"))
    print(f"OK    {'parcoords figures build':<34} {n1} / {n2} chars")

print("=" * 78)
if failures:
    print(f"SCIENTIFIC SMOKE TEST FAILED — {len(failures)} check(s):")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("SCIENTIFIC SMOKE TEST PASSED")
