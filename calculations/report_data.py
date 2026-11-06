"""The numbers of the study report that no notebook saves: the model at the study's
operating point, level by level, and the facts of the study itself.

The notebooks save figures, and the sweep and sensitivity notebooks save tables; the
Level 1-5 and Level 6 walkthroughs save figures only. This fills that gap with the
same calls :mod:`container.model_smoke_test` and the sensitivity wrappers make, at
``CONDITIONS`` of the active study, so the report can quote each level's flux,
regime and resistance split. Everything is closed-form; a run takes well under a
minute, most of it imports.

Run it in the study's own process, because :mod:`calculations.config.model_config`
binds one study at import::

    PERMEATION_STUDY=Guo_etal_2025/study_config.py python -m calculations.report_data Guo_etal_2025

It writes ``<study>/report/data/operating_point.csv`` (level, quantity, value, units,
source) and ``<study>/report/data/study_facts.json``. :mod:`calculations.report` runs
it before every build, in a separate process, so the table always matches the
current ``study_config.py``.
"""
import contextlib
import csv
import datetime
import hashlib
import io
import json
import math
import os
import sys
import warnings

FLUX = 'mol m-2 s-1'
PERM = 'mol m-1 s-1 Pa-0.5'


def _plain(v):
    """A value JSON and CSV can hold."""
    if isinstance(v, (str, bool)) or v is None:
        return v
    if isinstance(v, (int, float)):
        return None if isinstance(v, float) and not math.isfinite(v) else v
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items()}
    try:
        return float(v)
    except (TypeError, ValueError):
        return str(v)


def operating_point():
    """[(level, quantity, value, units, source)] at CONDITIONS, and the study facts."""
    from calculations.config.model_config import (
        ACTIVE_STUDY_NAME, ACTIVE_STUDY_PATH, ACTIVE_STUDY_SOURCE, CONDITIONS, METALS,
        build_simulation_config)
    from calculations.interface_solver import calculate_oxide_metal_system
    from calculations.oxide_permeation import get_metal_properties_at_T, get_oxide_properties_at_T
    from calculations.permeation_calc import calculate_simple_metal_flux
    from calculations.sensitivity import (
        DEFAULT_PARAMS_LEVEL5, DEFAULT_PARAMS_LEVEL5L6, check_against_config,
        level5_model_wrapper, level5L6_model_wrapper)
    from calculations.surface_kinetics import solve_steady_state_flux_L1L6

    T = CONDITIONS['T_operating']
    P_up, P_down = CONDITIONS['P_upstream'], CONDITIONS['P_downstream']
    L_m, L_ox = CONDITIONS['L_metal'], CONDITIONS['L_oxide']
    sim = build_simulation_config()
    rows, errors = [], []

    def add(level, quantity, value, units, source):
        rows.append((level, quantity, _plain(value), units, source))

    def attempt(level, fn):
        try:
            fn()
        except Exception as exc:                  # noqa: BLE001 -- one level failing must not hide the others
            errors.append(f'{level}: {type(exc).__name__}: {exc}')
            add(level, 'error', f'{type(exc).__name__}: {exc}', '', '')

    mp = dict(get_metal_properties_at_T(sim['metal_name'], T), thickness=L_m)
    op = dict(get_oxide_properties_at_T(sim['oxide_name'], T), thickness=L_ox)

    def l1():
        r = calculate_simple_metal_flux(mp['D_metal'], mp['K_s_metal'], L_m, P_up, P_down)
        add('L1', 'flux', r['flux'], FLUX, 'permeation_calc.calculate_simple_metal_flux')

    def l2b():
        r = calculate_oxide_metal_system(P_up, P_down, op, mp, T_K=T)
        src = 'interface_solver.calculate_oxide_metal_system'
        add('L2b', 'flux', r['flux'], FLUX, src)
        add('L2b', 'P_interface', r['P_interface'], 'Pa', src)
        add('L2b', 'resistance_ratio', r['resistance_ratio'], '-', src)
        if 'flux_error' in r:
            add('L2b', 'flux_error', r['flux_error'], FLUX, src)

    def l1l6():
        sk = METALS[sim['metal_name']]['surface_kinetics']
        r = solve_steady_state_flux_L1L6(P_up=P_up, P_down=P_down, L_m=L_m, k_diss=sk['k_diss_metal_ref'],
                                         K_eq=sk['K_eq_metal_ref'], D_m=mp['D_metal'], K_s_m=mp['K_s_metal'])
        src = 'surface_kinetics.solve_steady_state_flux_L1L6'
        add('L1L6', 'flux', r['J_ss'], FLUX, src)
        add('L1L6', 'P_interface', r['P_int'], 'Pa', src)
        add('L1L6', 'theta', r['theta'], '-', src)
        add('L1L6', 'rate_limiting', r['rate_limiting'], '', src)

    record5 = {}

    def l5():
        r = level5_model_wrapper(dict(DEFAULT_PARAMS_LEVEL5), return_full_record=True)
        record5.update(r)
        src = 'sensitivity.level5_model_wrapper'
        for k, u in (('flux', FLUX), ('permeability', PERM), ('PRF', '-'), ('P_interface', 'Pa'),
                     ('flux_intact', FLUX), ('flux_defect', FLUX), ('flux_bare_metal', FLUX),
                     ('D_eff', 'm2 s-1'), ('D_modification', '-'), ('frac_oxide', '-'), ('frac_metal', '-'),
                     ('frac_defect', '-'), ('regime', ''), ('regime_hierarchy', ''), ('dominant_path', '')):
            if k in r:
                add('L5', k, r[k], u, src)

    def l5l6():
        r = level5L6_model_wrapper(dict(DEFAULT_PARAMS_LEVEL5L6), return_full_record=True)
        src = 'sensitivity.level5L6_model_wrapper'
        for k, u in (('flux', FLUX), ('permeability', PERM), ('PRF', '-'), ('P_interface', 'Pa'),
                     ('flux_intact', FLUX), ('flux_defect', FLUX), ('theta', '-'), ('eta_surf', '-'),
                     ('frac_surface', '-'), ('frac_oxide', '-'), ('frac_metal', '-'), ('regime', ''),
                     ('system_rate_limiting', ''), ('dominant_path', '')):
            if k in r:
                add('L5L6', k, r[k], u, src)

    for level, fn in (('L1', l1), ('L2b', l2b), ('L1L6', l1l6), ('L5', l5), ('L5L6', l5l6)):
        attempt(level, fn)

    with open(ACTIVE_STUDY_PATH, 'rb') as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    try:
        problems = check_against_config(raise_on_fail=False, verbose=False)
    except Exception as exc:                      # noqa: BLE001
        problems = [f'check_against_config itself failed: {type(exc).__name__}: {exc}']
    try:
        from importlib.metadata import version
        pkg = version('continuum-h-permeation')
    except Exception:                             # noqa: BLE001 -- not installed: run from the tree
        pkg = None
    properties = {k: record5.get(k) for k in ('D_metal', 'K_s_metal', 'D_ox', 'K_ox')}
    facts = {
        'study_name': ACTIVE_STUDY_NAME, 'study_path': str(ACTIVE_STUDY_PATH),
        'study_source': str(ACTIVE_STUDY_SOURCE), 'study_sha256': sha,
        'metal_name': sim.get('metal_name'), 'oxide_name': sim.get('oxide_name'),
        'conditions': _plain(dict(CONDITIONS)), 'properties_at_T': _plain(properties),
        'config_problems': list(problems or []), 'errors': errors,
        'computed_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'python': sys.version.split()[0], 'package_version': pkg,
    }
    return rows, facts


def write(study_dir, out_dir=None):
    """Compute and write the two files; returns their folder."""
    from calculations.config.model_config import ACTIVE_STUDY_PATH
    want = os.path.realpath(os.path.join(study_dir, 'study_config.py'))
    if os.path.realpath(str(ACTIVE_STUDY_PATH)) != want:
        raise SystemExit(f'report_data: the active study is {ACTIVE_STUDY_PATH}, not {want}; '
                         f'run with PERMEATION_STUDY={want}')
    out_dir = out_dir or os.path.join(study_dir, 'report', 'data')
    os.makedirs(out_dir, exist_ok=True)
    printed = io.StringIO()
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(printed):
        warnings.simplefilter('always')
        rows, facts = operating_point()
    # what the model said while computing: extrapolation notices and the like, once each
    said = [line.strip() for line in printed.getvalue().splitlines() if 'warning' in line.lower()]
    said += [f'{w.category.__name__}: {w.message}' for w in caught]
    facts['model_warnings'] = list(dict.fromkeys(s for s in said if s))
    with open(os.path.join(out_dir, 'operating_point.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['level', 'quantity', 'value', 'units', 'source'])
        w.writerows(rows)
    with open(os.path.join(out_dir, 'study_facts.json'), 'w') as fh:
        json.dump(facts, fh, indent=2, ensure_ascii=False)
    return out_dir


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__)
        return 2
    out = write(os.path.abspath(argv[0]))
    print(f'report_data: {out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
