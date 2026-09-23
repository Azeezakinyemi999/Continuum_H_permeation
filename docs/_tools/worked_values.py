#!/usr/bin/env python3
"""Generate the worked values used in the theory chapters.

Every number quoted in docs/theory/ comes from here, so a model change shows up
as a diff in this output rather than as a quietly wrong figure in the prose.

    python docs/_tools/worked_values.py            # active study
    python docs/_tools/worked_values.py --check     # diff against the stored copy

The stored copy is docs/_tools/worked_values.txt. Regenerate it with:
    python docs/_tools/worked_values.py > docs/_tools/worked_values.txt
"""
from __future__ import annotations
import argparse
import contextlib
import io
import pathlib
import subprocess
import sys
import warnings

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
STORED = pathlib.Path(__file__).with_suffix(".txt")


def fmt(v, w=28):
    if isinstance(v, float):
        return f"  {'':<0}{v:.6e}"
    return f"  {v}"


def report() -> str:
    # Some model paths report range violations with print(), not warnings.warn,
    # so simplefilter alone does not silence them and they would otherwise be
    # interleaved into this report. Captured and discarded; the caveat is
    # documented in the theory chapters instead.
    warnings.simplefilter("ignore")
    out: list[str] = []
    def w(s=""): out.append(s)
    _noise = io.StringIO()

    from calculations.config.model_config import (
        ACTIVE_STUDY, METALS, OXIDES, MICROSTRUCTURE, CONDITIONS,
        build_simulation_config, DEFAULT_PARAMS_LEVEL5, DEFAULT_PARAMS_LEVEL5L6)
    from calculations.oxide_permeation import (
        get_metal_properties_at_T, get_oxide_properties_at_T,
        molecular_diffusion_flux, calculate_oxide_resistance)
    from calculations.permeation_calc import calculate_simple_metal_flux
    from calculations.interface_solver import calculate_oxide_metal_system
    from calculations.parallel_oxide_defect_paths import (
        calculate_parallel_path_flux, calculate_PRF)
    from calculations.defective_metal import combined_microstructure_model
    from calculations.surface_kinetics import solve_steady_state_flux_L1L6
    from calculations.sensitivity import level5_model_wrapper, level5L6_model_wrapper

    SIM = build_simulation_config()
    T   = CONDITIONS["T_operating"]
    PU  = CONDITIONS["P_upstream"]
    PD  = CONDITIONS["P_downstream"]
    LM  = CONDITIONS["L_metal"]
    LOX = OXIDES[SIM["oxide_name"]]["thickness"]

    w(f"study            {ACTIVE_STUDY}")
    w(f"metal / oxide    {SIM['metal_name']} / {SIM['oxide_name']}")
    w(f"T                {T} K")
    w(f"P_up / P_down    {PU:g} / {PD:g} Pa")
    w(f"L_metal/L_oxide  {LM:g} / {LOX:g} m")

    mp = dict(get_metal_properties_at_T(SIM["metal_name"], T), thickness=LM)
    op = dict(get_oxide_properties_at_T(SIM["oxide_name"], T), thickness=LOX)
    w()
    w("[properties at T]")
    for k in ("D_metal", "K_s_metal"):
        w(f"  {k:<26} {mp[k]:.6e}")
    for k in ("D_ox", "K_ox"):
        w(f"  {k:<26} {op[k]:.6e}")

    w()
    w("[L1 perfect metal] permeation_calc.calculate_simple_metal_flux")
    r1 = calculate_simple_metal_flux(mp["D_metal"], mp["K_s_metal"], LM, PU, PD)
    for k in ("flux", "C_up", "C_down", "permeability"):
        if k in r1:
            w(f"  {k:<26} {r1[k]:.6e}")

    w()
    w("[L2a perfect oxide] oxide_permeation.molecular_diffusion_flux")
    w(f"  {'flux':<26} {molecular_diffusion_flux(op['D_ox'], op['K_ox'], LOX, PU, PD):.6e}")
    w(f"  {'R_oxide':<26} {calculate_oxide_resistance(op['D_ox'], op['K_ox'], LOX):.6e}")

    w()
    w("[L2b series] interface_solver.calculate_oxide_metal_system")
    r2 = calculate_oxide_metal_system(PU, PD, op, mp, T_K=T)
    for k in ("flux", "P_interface", "resistance_ratio", "regime"):
        v = r2.get(k)
        w(f"  {k:<26} {v:.6e}" if isinstance(v, float) else f"  {k:<26} {v}")

    w()
    w("[L3 defective oxide] parallel_oxide_defect_paths.calculate_parallel_path_flux")
    r3 = calculate_parallel_path_flux(PU, PD, op, mp, SIM["oxide_defects"])
    for k in ("flux_total", "flux_intact_contribution", "flux_defect_contribution"):
        if k in r3:
            w(f"  {k:<26} {r3[k]:.6e}")
    prf = calculate_PRF(PU, op, mp, defect_params=SIM["oxide_defects"], P_downstream=PD)
    for k in ("PRF", "PRF_perfect"):
        if isinstance(prf, dict) and k in prf:
            w(f"  {k:<26} {prf[k]:.6e}")

    w()
    w("[L4 microstructure] defective_metal.combined_microstructure_model")
    w("  (lattice_concentration=10.0, as in container/model_smoke_test.py)")
    m = combined_microstructure_model(mp["D_metal"], T, dict(MICROSTRUCTURE),
                                      lattice_concentration=10.0,
                                      lattice_density=MICROSTRUCTURE["lattice_density"])
    for k in ("D_lattice", "D_eff", "overall_factor"):
        w(f"  {k:<26} {m[k]:.6e}")
    w(f"  {'trapping.theta_total':<26} {m['trapping']['theta_total']:.6e}")
    w(f"  {'dominant_effect':<26} {m['dominant_effect']}")
    w(f"  {'regime':<26} {m['regime']}")

    w()
    w("[L1+L6 surface kinetics] surface_kinetics.solve_steady_state_flux_L1L6")
    sk = METALS[SIM["metal_name"]]["surface_kinetics"]
    r6 = solve_steady_state_flux_L1L6(
        P_up=PU, P_down=PD, L_m=LM,
        k_diss=sk["k_diss_metal_ref"], K_eq=sk["K_eq_metal_ref"],
        D_m=mp["D_metal"], K_s_m=mp["K_s_metal"])
    w("  (k_diss/K_eq are the study's *_ref values, as pinned by the smoke test)")
    for k in ("J_ss", "P_int", "theta", "beta", "rate_limiting"):
        v = r6.get(k)
        w(f"  {k:<26} {v:.6e}" if isinstance(v, float) else f"  {k:<26} {v}")

    w()
    w("[permeability tiers] calculations.permeability")
    from calculations import permeability as pmod
    Phi_ox_1 = pmod.intrinsic_oxide_permeability(op["D_ox"], op["K_ox"], pmod.MODEL_1)
    Phi_m_1  = pmod.intrinsic_metal_permeability(mp["D_metal"], mp["K_s_metal"])
    w(f"  tier 1 Phi_oxide           {Phi_ox_1['permeability']:.6e}  {Phi_ox_1['units']}")
    w(f"  tier 1 Phi_metal           {Phi_m_1['permeability']:.6e}  {Phi_m_1['units']}")
    w("  Model 1 has no tier 3: the two units above are not commensurable.")

    # tier 2 -- free-standing defective oxide, cracks and oxide GB only
    for fc, fg in ((1e-3, 0.0), (0.0, 1e-3), (1e-2, 1e-2)):
        ox_eff = pmod.oxide_only_permeability(
            Phi_ox_1["permeability"], f_crack=fc, gamma=0.1, f_gb=fg, beta=10.0)
        w(f"  tier 2 Phi_ox_eff/Phi_ox   {ox_eff['enhancement_ratio']:.6f}"
          f"   (f_crack={fc:g}, f_gb={fg:g})")

    # Model 1 bilayer: closed form against the brentq solve
    cf = pmod.model1_bilayer_closed_form(
        Phi_ox_1["permeability"], LOX, Phi_m_1["permeability"], LM, PU, PD)
    l2b = calculate_oxide_metal_system(PU, PD, op, mp)
    w(f"  Model 1 bilayer closed form {cf['flux']:.6e}  mol/m2/s")
    w(f"  Model 1 bilayer brentq      {l2b['flux']:.6e}  mol/m2/s")
    w(f"  relative difference         {abs(cf['flux']/l2b['flux'] - 1):.3e}")

    # Model 2 decomposition on the pristine bilayer
    sk_ox = OXIDES[SIM["oxide_name"]]["surface_kinetics"]
    from calculations.surface_kinetics import solve_steady_state_flux_direct
    r2 = solve_steady_state_flux_direct(
        PU, PD, LM, sk_ox["k_diss_ref"], sk_ox["K_eq_ref"],
        op["D_ox"], op["K_ox"], LOX, mp["D_metal"], mp["K_s_metal"])
    t3 = pmod.transport_permeability(
        op["D_ox"] * op["K_ox"], LOX,
        mp["D_metal"] * mp["K_s_metal"], LM, pmod.MODEL_2)
    eta = pmod.surface_efficiency(r2["theta"], sk_ox["K_eq_ref"], PU, PD)
    app = pmod.apparent_permeability(r2["J_ss"], PU, PD, LOX + LM)
    w(f"  tier 3 Phi_transport       {t3['permeability']:.6e}  {t3['units']}")
    w(f"  Model 2 eta_surf           {eta:.6f}")
    w(f"  tier 4 Phi_app             {app:.6e}")
    w(f"  Phi_transport * eta_surf   {t3['permeability'] * eta:.6e}  (must equal Phi_app)")
    Phi_ox_v = op["D_ox"] * op["K_ox"]
    Phi_m_v  = mp["D_metal"] * mp["K_s_metal"]
    old_val  = 1.0 / (1.0 / Phi_ox_v + 1.0 / Phi_m_v)
    w(f"  old unweighted harmonic    {old_val:.6e}")
    w(f"  understated by             {t3['permeability'] / old_val:.0f}x")

    for label, fn, params in (("L5", level5_model_wrapper, DEFAULT_PARAMS_LEVEL5),
                              ("L5+L6", level5L6_model_wrapper, DEFAULT_PARAMS_LEVEL5L6)):
        w()
        w(f"[{label} canonical wrapper] sensitivity.{fn.__name__} on its DEFAULT_PARAMS")
        r = fn(params)
        for k in sorted(r):
            v = r[k]
            if isinstance(v, bool) or isinstance(v, str):
                w(f"  {k:<26} {v}")
            elif isinstance(v, (int, float)):
                w(f"  {k:<26} {v:.6e}")
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="diff against the stored copy; exit 1 if it drifted")
    args = ap.parse_args()
    # calculations/oxide_permeation.py:181 reports an out-of-range temperature
    # with print(), not warnings.warn, so it cannot be silenced by simplefilter
    # and would otherwise be captured into the stored report. The active study
    # runs at 873 K, outside Cr2O3's validated [473, 773] K, so this fires on
    # every default call — see the caveat in the theory chapters.
    _noise = io.StringIO()
    with contextlib.redirect_stdout(_noise):
        text = report()
    if not args.check:
        sys.stdout.write(text)
        return 0
    if not STORED.exists():
        print(f"no stored copy at {STORED}; generate one first")
        return 1
    if STORED.read_text() == text:
        print(f"worked values unchanged ({STORED.name})")
        return 0
    print(f"worked values DRIFTED from {STORED.name}:\n")
    import difflib
    for line in difflib.unified_diff(STORED.read_text().splitlines(),
                                     text.splitlines(),
                                     "stored", "computed", lineterm="", n=1):
        print(f"  {line}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
