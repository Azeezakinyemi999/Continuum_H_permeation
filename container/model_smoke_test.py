#!/usr/bin/env python3
"""Gate 2 — does the scientific application still compute the right numbers?

Environment reproducibility (gate 1) does not imply scientific reproducibility.
This asserts pinned reference values through the real solvers.

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
L_METAL, L_OXIDE = 1.0e-3, 4.8e-08

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
        ACTIVE_STUDY, MICROSTRUCTURE, build_simulation_config,
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
equal("ACTIVE_STUDY", ACTIVE_STUDY, "incoloy802_cr2o3")

SIM = build_simulation_config()
equal("metal resolves in METALS", SIM["metal_name"],
      "metal_X40_NiCrAlTi_31_19_Incoloy802_Schmidt1985")
equal("oxide resolves in OXIDES", SIM["oxide_name"], "Cr2O3_sample4")

# ---------------------------------------------------------------- analytic
print("\nanalytic:")
close("arrhenius", arrhenius(1.0e-11, 50_000.0, 800.0, 700.0),
      2.926830409475082e-11)

# ------------------------------------------------------- L1: perfect metal
print("\nLevel 1 (perfect metal):")
# NB: the T-evaluated getters return D_metal / K_s_metal, not D / K_s.
mp = dict(get_metal_properties_at_T(SIM["metal_name"], T), thickness=L_METAL)
op = dict(get_oxide_properties_at_T(SIM["oxide_name"], T), thickness=L_OXIDE)
r1 = calculate_simple_metal_flux(mp["D_metal"], mp["K_s_metal"],
                                 L_METAL, P_UP, P_DOWN)
close("L1 flux", r1["flux"], 1.955109516956673e-07)

# --------------------------------------------- L2b: oxide + metal (brentq)
print("\nLevel 2b (oxide+metal, brentq):")
# calculate_oxide_metal_system needs `thickness` in BOTH props dicts; the
# T-evaluated getter does not supply it for the metal.
r2 = calculate_oxide_metal_system(P_UP, P_DOWN, op, mp, T_K=T)
close("L2b flux", r2["flux"], 1.888003770552079e-07)
close("L2b P_interface", r2["P_interface"], 93462.901558976)
close("L2b resistance_ratio", r2["resistance_ratio"], 0.036154225782795736)
if not (isinstance(r2.get("flux_error"), float) and abs(r2["flux_error"]) < 1e-12):
    print(f"FAIL  L2b solver converged           flux_error={r2.get('flux_error')!r}")
    failures.append("L2b solver converged")
else:
    print(f"OK    {'L2b solver converged':<34} flux_error={r2['flux_error']:.3e}")

# ------------------------------------- L1L6: surface kinetics (brentq #2)
print("\nLevel 1+6 (surface kinetics, brentq):")
# Returns J_ss / P_int / theta / beta / rate_limiting — there is no 'flux' key.
r3 = solve_steady_state_flux_L1L6(
    P_up=P_UP, P_down=P_DOWN, L_m=L_METAL, k_diss=1.346e-06, K_eq=1.0e-4,
    D_m=6.881385163348364e-11, K_s_m=0.03373208680855995)
close("L1L6 J_ss", r3["J_ss"], 7.10792976386876e-07)
close("L1L6 P_int", r3["P_int"], 99990.85191948501)
close("L1L6 theta", r3["theta"], 0.7597385771009892)
equal("L1L6 rate_limiting", r3["rate_limiting"], "metal")

# ------------------------------------------- microstructure (4 trap types)
print("\nLevel 4 (microstructure, 4 traps):")
m = combined_microstructure_model(
    6.881385163348364e-11, 873.0, dict(MICROSTRUCTURE),
    lattice_concentration=10.0, lattice_density=MICROSTRUCTURE["lattice_density"])
close("micro D_eff", m["D_eff"], 5.084388686889583e-11)
close("micro overall_factor", m["overall_factor"], 0.7388612272380938)
close("micro theta_total", m["trapping"]["theta_total"], 0.026370841516212183)

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
    print(f"SKIP  {R}/ not present — mount the workspace to exercise this")
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
