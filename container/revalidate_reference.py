"""Cross-interpreter revalidation for the 3.9 -> 3.12 migration.

Prints a fixed set of solver / SALib / pandas results. Run under both
interpreters and diff the output; anything other than an empty diff blocks the
release. Confirmed identical on Python 3.9.23/macOS-x86_64 and
Python 3.12.14/linux-amd64 with numpy 1.26.4, scipy 1.13.1, pandas 2.3.3,
SALib 1.5.0 (2026-08-24).

Phase D will fold these into container/model_smoke_test.py with assertions.
"""
import io, warnings, numpy as np, pandas as pd
warnings.simplefilter('ignore')

out = {}

# 1. surface kinetics solver (brentq on a different equation than L2b)
from calculations.surface_kinetics import solve_steady_state_flux_L1L6
r = solve_steady_state_flux_L1L6(P_up=1.0e5, P_down=1.0e2, L_m=1.0e-3,
                                 k_diss=1.346e-06, K_eq=1.0e-4,
                                 D_m=6.881385163348364e-11, K_s_m=0.03373208680855995)
for k in ('J_ss', 'P_int', 'theta', 'beta', 'rate_limiting'):
    out[f'L1L6_{k}'] = r[k]

# 2. defective metal / microstructure (numpy-heavy, many traps)
from calculations.config.model_config import MICROSTRUCTURE
from calculations.defective_metal import combined_microstructure_model
m = combined_microstructure_model(6.881385163348364e-11, 873.0, dict(MICROSTRUCTURE),
                                  lattice_concentration=10.0,
                                  lattice_density=MICROSTRUCTURE['lattice_density'])
out['micro_D_eff']          = m['D_eff']
out['micro_overall_factor'] = m['overall_factor']
out['micro_theta_total']    = m['trapping']['theta_total']

# 3. SALib: seeded Latin hypercube -> PAWN + delta (the SA path)
from SALib.sample import latin
from SALib.analyze import pawn, delta as delta_an
problem = {'num_vars': 3, 'names': ['a','b','c'],
           'bounds': [[1.0, 10.0], [0.1, 1.0], [100.0, 1000.0]]}
X = latin.sample(problem, 256, seed=12345)
Y = np.log10(X[:,0]) * X[:,1] + np.sqrt(X[:,2])
out['salib_X_sum']    = float(X.sum())
out['salib_Y_mean']   = float(Y.mean())
pw = pawn.analyze(problem, X, Y, S=10, seed=12345, print_to_console=False)
dl = delta_an.analyze(problem, X, Y, num_resamples=10, seed=12345, print_to_console=False)
out['pawn_median']    = [float(v) for v in pw['median']]
out['delta']          = [float(v) for v in dl['delta']]
out['S1']             = [float(v) for v in dl['S1']]

# 4. pandas CSV round-trip (SA scans are written and re-read as CSV)
df = pd.DataFrame({'regime': ['metal','oxide','surface'],
                   'flux': [1.888003770552079e-07, 3.1e-12, 7.77e-09]})
buf = io.StringIO(); df.to_csv(buf, index=False); buf.seek(0)
out['csv_roundtrip_exact'] = bool((pd.read_csv(buf)['flux'] == df['flux']).all())

import sys
print('PYTHON', sys.version.split()[0], '| numpy', np.__version__, '| pandas', pd.__version__)
for k in sorted(out):
    print(f'{k:24s} {out[k]!r}')
