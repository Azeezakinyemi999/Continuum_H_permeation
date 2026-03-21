---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.7
  kernelspec:
    display_name: rmg_rms_env
    language: python
    name: python3
---

# L1+L6: Perfect Metal + Surface Kinetics — Analysis Notebook

**Physics:**
- Surface dissociation/recombination at the gas–metal interface (L6)
- Bulk lattice diffusion of atomic H through a perfect metal (L1)
- No oxide layer

**Two competing resistances:**
- $R_{surface} \sim 1 / (k_{diss} \cdot P_{up} \cdot (1-\theta)^2)$
- $R_{metal} = L_m / (D_m \cdot K_{s,m})$

**Key physics question:** At what pressure and temperature does the rate-limiting step shift from surface dissociation to bulk diffusion?


```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import sys

# Add parent directory to path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath('__file__')))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)
```

```python
# =============================================================================
# CELL 1 — SETUP
# Run once. All imports and helper functions.
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import brentq
import pandas as pd
from itertools import groupby
from operator import itemgetter
from IPython.display import display

from calculations.surface_kinetics import solve_steady_state_flux_L1L6
from calculations.config.model_config import (
    R, F,
    METALS, CONDITIONS, VALIDATION,
    PLOT_STYLE as STYLE,
    CURVE_STYLES,
    apply_style,
    build_simulation_config,
)

# ── Simulation config ────────────────────────────────────────────────────────
SIM       = build_simulation_config(
    metal='metal_316L_Heat_treated_ref_cast',
    oxide='Al2O3',
    T_operating=873,
    P_upstream=1e5,
    L_metal=1e-3,
)

METAL_KEY = SIM['metal_name']
METAL     = METALS[METAL_KEY]
SK        = METAL['surface_kinetics']     # surface kinetics sub-dict

# ── Pressure sweep (matches widget range) ────────────────────────────────────
P_RANGE_ARR = np.logspace(
    np.log10(CONDITIONS['P_range'][0]),
    np.log10(CONDITIONS['P_range'][1]),
    CONDITIONS['n_P_points']
)   #np.logspace(-3, 12, 100)
P_DOWN      = SIM['P_downstream']         # Pa
T_RANGE     = CONDITIONS['T_range']       # K
N_T         = CONDITIONS['n_T_points']

# ── Arrhenius helper functions ───────────────────────────────────────────────
def get_k_diss_at_T(T_K):
    """k_diss at temperature T_K via Arrhenius from config."""
    return SK['k_diss_metal_ref'] * np.exp(
        (-SK['E_diss_metal'] / R) * (1/T_K - 1/SK['T_ref'])
    )

def get_K_eq_at_T(T_K):
    """K_eq at temperature T_K via van't Hoff from config."""
    return SK['K_eq_metal_ref'] * np.exp(
        (-SK['H_eq_metal'] / R) * (1/T_K - 1/SK['T_ref'])
    )

def get_D_m_at_T(T_K):
    """D_m at temperature T_K via Arrhenius from config."""
    return METAL['D_ref'] * np.exp(
        (-METAL['E_D'] / R) * (1/T_K - 1/METAL['T_ref'])
    )

def get_K_s_m_at_T(T_K):
    """K_s_m at temperature T_K via Arrhenius from config."""
    return METAL['K_s_ref'] * np.exp(
        (-METAL['H_s'] / R) * (1/T_K - 1/METAL['T_ref'])
    )

# ── Known reference activation energies (for validation) ────────────────────
E_REF_SURFACE = SK['E_diss_metal']              # J/mol — surface dissociation
E_REF_METAL   = METAL['E_D'] + METAL['H_s']    # J/mol — bulk permeation

print("=" * 60)
print("SETUP COMPLETE")
print("=" * 60)
print(f"  Metal:           {METAL_KEY}")
print(f"  E_ref surface:   {E_REF_SURFACE/1000:.1f} kJ/mol  (E_diss)")
print(f"  E_ref metal:     {E_REF_METAL/1000:.1f} kJ/mol  (E_D + H_s)")
print(f"  P range:         {P_RANGE_ARR[0]:.0e} – {P_RANGE_ARR[-1]:.0e} Pa")
print(f"  T range:         {T_RANGE[0]} – {T_RANGE[1]} K")
print("=" * 60)
```

```python
P_RANGE_ARR
```

```python
# =============================================================================
# CELL 2 — PARAMETERS
# All physics values at operating temperature, pulled from config.
# Change placeholder values in config and rerun from this cell downward.
# =============================================================================

# ── Operating conditions ─────────────────────────────────────────────────────
T_op  = SIM['T_operating']    # K
L_m   = SIM['L_metal']        # m
P_down = SIM['P_downstream']  # Pa

# ── Surface kinetics at T_op ─────────────────────────────────────────────────
k_diss = get_k_diss_at_T(T_op)
K_eq   = get_K_eq_at_T(T_op)
k_recomb = k_diss / K_eq      # derived — not stored separately

# ── Metal transport at T_op ──────────────────────────────────────────────────
D_m   = get_D_m_at_T(T_op)
K_s_m = get_K_s_m_at_T(T_op)

# ── Derived ──────────────────────────────────────────────────────────────────
beta = D_m * K_s_m / L_m      # metal permeance [mol/m²/s/Pa^0.5]

# ── Print verification ───────────────────────────────────────────────────────
print("=" * 60)
print(f"L1+L6 PARAMETERS AT T = {T_op-273.15:.0f}°C  ({T_op} K)")
print("=" * 60)

print(f"\n  Geometry")
print(f"    L_m      = {L_m:.3e}  m")
print(f"    P_down   = {P_down:.1e}  Pa")

print(f"\n  Surface kinetics  [config: surface_kinetics sub-dict]")
print(f"    k_diss   = {k_diss:.3e}  mol/m²/s/Pa")
print(f"    K_eq     = {K_eq:.3e}  Pa⁻¹")
print(f"    k_recomb = {k_recomb:.3e}  mol/m²/s  (= k_diss / K_eq)")

print(f"\n  Metal transport  [config: Arrhenius]")
print(f"    D_m      = {D_m:.3e}  m²/s")
print(f"    K_s_m    = {K_s_m:.3e}  mol/m³/Pa^0.5")

print(f"\n  Derived")
print(f"    β        = {beta:.3e}  mol/m²/s/Pa^0.5  (= D_m × K_s_m / L_m)")

print(f"\n  Reference activation energies  [config]")
print(f"    E_diss   = {E_REF_SURFACE/1000:.1f} kJ/mol   (surface dissociation)")
print(f"    E_D+Hs   = {E_REF_METAL/1000:.1f} kJ/mol   (bulk permeation)")

# ── Placeholder warning ──────────────────────────────────────────────────────
if SK['k_diss_metal_ref'] <= 1e-12:
    print(f"\n  ⚠  WARNING: k_diss_metal_ref = {SK['k_diss_metal_ref']:.0e} looks like")
    print(f"     a placeholder. Update config surface_kinetics sub-dict")
    print(f"     with literature values before interpreting results.")

print("\n" + "=" * 60)
```

```python
# =============================================================================
# CELL 3 — COMPUTE
# Single-pass loop over full pressure range.
# Builds all arrays, regime summary, and Arrhenius pressure selection.
# =============================================================================

# ── Single-pass loop ─────────────────────────────────────────────────────────
rows             = []
J_arr            = []
theta_arr        = []
P_int_arr        = []
frac_surface_arr = []
frac_metal_arr   = []
R_surface_arr    = []
R_metal_arr      = []
R_total_arr      = []
rate_lim_arr     = []

for P_up in P_RANGE_ARR:
    try:
        r = solve_steady_state_flux_L1L6(
            P_up, P_down, L_m,
            k_diss, K_eq, D_m, K_s_m
        )
        J_arr.append(r['J_ss'])
        theta_arr.append(r['theta'])
        P_int_arr.append(r['P_int'])
        frac_surface_arr.append(r['resistances']['fraction_surface'])
        frac_metal_arr.append(r['resistances']['fraction_metal'])
        R_surface_arr.append(r['resistances']['R_surface'])
        R_metal_arr.append(r['resistances']['R_metal'])
        R_total_arr.append(r['resistances']['R_total'])
        rate_lim_arr.append(r['rate_limiting'])

        rows.append({
            'P_up (Pa)':            P_up,
            'P_int (Pa)':           r['P_int'],
            'J_ss (mol/m²/s)':      r['J_ss'],
            'θ_surface':            r['theta'],
            'β_metal':              r['beta'],
            'R_surface':            r['resistances']['R_surface'],
            'R_metal':              r['resistances']['R_metal'],
            'R_total':              r['resistances']['R_total'],
            'fraction_surface (%)': r['resistances']['fraction_surface'] * 100,
            'fraction_metal (%)':   r['resistances']['fraction_metal']   * 100,
            'Rate-Limiting':        r['rate_limiting'].upper(),
        })

    except Exception as e:
        J_arr.append(np.nan);            theta_arr.append(np.nan)
        P_int_arr.append(np.nan);        frac_surface_arr.append(np.nan)
        frac_metal_arr.append(np.nan);   R_surface_arr.append(np.nan)
        R_metal_arr.append(np.nan);      R_total_arr.append(np.nan)
        rate_lim_arr.append('error')
        rows.append({'P_up (Pa)': P_up, 'Rate-Limiting': 'ERROR', 'Error': str(e)})

# ── Convert to arrays ────────────────────────────────────────────────────────
J_arr            = np.array(J_arr)
theta_arr        = np.array(theta_arr)
P_int_arr        = np.array(P_int_arr)
frac_surface_arr = np.array(frac_surface_arr)
frac_metal_arr   = np.array(frac_metal_arr)
R_surface_arr    = np.array(R_surface_arr)
R_metal_arr      = np.array(R_metal_arr)
R_total_arr      = np.array(R_total_arr)
rate_lim_arr     = np.array(rate_lim_arr)
valid            = ~np.isnan(J_arr)

# ── Rate-limiting classification array (for plot shading) ────────────────────
rl_arr = np.where(frac_surface_arr > 0.5, 'surface',
         np.where(frac_metal_arr   > 0.5, 'metal', 'mixed'))

# ── Regime summary ───────────────────────────────────────────────────────────
# Count points in each regime
n_surface = np.sum((frac_surface_arr > 0.5)  & valid)
n_metal   = np.sum((frac_metal_arr   > 0.5)  & valid)
n_mixed   = np.sum(
    (frac_surface_arr <= 0.5) & (frac_metal_arr <= 0.5) & valid
)

# Pressure ranges for each dominant regime (>90% — clean, not transitional)
surf_dom_mask  = (frac_surface_arr > 0.90) & valid
metal_dom_mask = (frac_metal_arr   > 0.90) & valid

P_surf_dom  = P_RANGE_ARR[surf_dom_mask]
P_metal_dom = P_RANGE_ARR[metal_dom_mask]

# ── Arrhenius pressure selection ─────────────────────────────────────────────
# Use median of each dominant region — no dependence on crossover
P_for_surface_Ea = np.median(P_surf_dom)  if len(P_surf_dom)  > 3 else None
P_for_metal_Ea   = np.median(P_metal_dom) if len(P_metal_dom) > 3 else None

# ── Crossover pressure (purely diagnostic, not used for Ea selection) ────────
P_crossover = None
valid_cross = valid & ~np.isnan(frac_surface_arr)
if np.sum(valid_cross) > 2:
    fs_v   = frac_surface_arr[valid_cross]
    P_v    = P_RANGE_ARR[valid_cross]
    sc     = np.where(np.diff(np.sign(fs_v - 0.5)))[0]
    if len(sc) > 0:
        idx = sc[0]
        # Direction check: surface fraction must drop (not rise) at crossover
        if idx + 1 < len(fs_v) and fs_v[idx] > 0.5 and fs_v[idx + 1] < 0.5:
            P_crossover = P_v[idx]

# ── Print regime summary ─────────────────────────────────────────────────────
print("=" * 60)
print(f"L1+L6 REGIME SUMMARY  (T = {T_op-273.15:.0f}°C, {np.sum(valid)} valid points)")
print("=" * 60)

print(f"\n  Regime breakdown (fraction > 50%)")
print(f"    Surface-limited:  {n_surface:3d} points")
print(f"    Metal-limited:    {n_metal:3d} points")
print(f"    Mixed:            {n_mixed:3d} points")

print(f"\n  Dominant regions (fraction > 90%) — used for Arrhenius")
if len(P_surf_dom) > 3:
    print(f"    Surface-dominant: {len(P_surf_dom):3d} points  "
          f"P = {P_surf_dom.min():.1e} – {P_surf_dom.max():.1e} Pa")
    print(f"    → P_for_surface_Ea = {P_for_surface_Ea:.2e} Pa")
else:
    print(f"    Surface-dominant: {len(P_surf_dom):3d} points  "
          f"(insufficient — surface Ea not extractable)")

if len(P_metal_dom) > 3:
    print(f"    Metal-dominant:   {len(P_metal_dom):3d} points  "
          f"P = {P_metal_dom.min():.1e} – {P_metal_dom.max():.1e} Pa")
    print(f"    → P_for_metal_Ea  = {P_for_metal_Ea:.2e} Pa")
else:
    print(f"    Metal-dominant:   {len(P_metal_dom):3d} points  "
          f"(insufficient — metal Ea not extractable)")

print(f"\n  Crossover pressure (diagnostic only)")
if P_crossover is not None:
    print(f"    P_crossover = {P_crossover:.2e} Pa")
else:
    print(f"    No crossover detected in physical pressure range")
    dominant = 'surface-limited' if n_surface > n_metal else 'metal-limited'
    print(f"    System is predominantly {dominant} across full range")

print(f"\n  θ range:  {np.nanmin(theta_arr):.4f} → {np.nanmax(theta_arr):.4f}")
print(f"  J range:  {np.nanmin(J_arr):.2e} → {np.nanmax(J_arr):.2e} mol/m²/s")
print("=" * 60)
```

```python
# =============================================================================
# CELL 4 — FIGURE 1: Core Validation (2×2)
# (A) Flux vs Pressure
# (B) Surface Coverage θ vs Pressure
# (C) Resistance Fractions vs Pressure
# (D) Limit Check: Sieverts Recovery
# =============================================================================


"""
Figure 1 validates the L1+L6 model from four angles at the operating temperature.

Panel (A) — Flux vs Pressure
    Question: Does the model produce the correct pressure-dependent behaviour?
    The hero plot. Shows J_ss vs P_up on a log-log scale. The curve should
    transition from slope=1 (surface-limited, low P) to slope=0.5
    (metal-limited, high P). Style follows the widget — no fill_between.
    Instead the curve itself is recolored by regime using thick colored
    segments, with a slope annotation box sitting directly on each segment.
    The black model curve sits underneath, the colored overlay sits on top.
    Reference lines for slope=1 and slope=0.5 are drawn as dashed lines
    anchored to the start and end of the curve respectively.
    Net slope annotation box in the bottom right corner.

Panel (B) — Surface Coverage θ vs Pressure
    Question: How does surface occupancy evolve with pressure, and does it
    match the Langmuir equilibrium isotherm?
    Plots the steady-state θ from the compute loop alongside the analytical
    Langmuir isotherm θ_eq = sqrt(K_eq * P) / (1 + sqrt(K_eq * P)) as a
    gray dash-dot line. The gap between θ_ss and θ_eq is physically
    meaningful — when they coincide the surface is in local equilibrium
    (metal-limited). When θ_ss < θ_eq the surface is being depleted by fast
    diffusion into the metal (surface-limited). Crossover pressure marked
    if it exists.

Panel (C) — Resistance Fractions vs Pressure
    Question: Which step controls the rate, and how does the balance shift
    with pressure?
    Plots R_surface/R_total and R_metal/R_total as percentages vs pressure
    on a semilog scale. The 50% threshold line marks the crossover. If no
    crossover exists in the physical range a text annotation states which
    regime dominates throughout. Crossover pressure annotated if present.

Panel (D) — Limit Check: Sieverts Recovery
    Question: As k_diss → ∞, does L1+L6 recover the pure Sieverts result?
    Runs the solver at k_diss = 1e-3 (fast kinetics limit) and compares
    J_L1+L6 against the analytical Sieverts flux J = beta*(sqrt(P_up) -
    sqrt(P_down)). A parity plot shows J_L1+L6 vs J_Sieverts — perfect
    agreement is a diagonal line. Max and mean relative errors are shown
    in an annotation box. Green box if max error < 1%, yellow otherwise.
    This validates the model correctly recovers L1 in the fast-kinetics limit.
"""

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle(
    f'L1+L6: Perfect Metal + Surface Kinetics — Core Validation\n'
    f'{METAL_KEY}  |  T = {T_op-273.15:.0f}°C  |  L = {L_m*1e3:.1f} mm',
    fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98
)

# ─────────────────────────────────────────────────────────────────────────────
# (A) FLUX vs PRESSURE — widget style: colored segments + slope annotations
# ─────────────────────────────────────────────────────────────────────────────
ax1 = axes[0, 0]

# Base black curve
ax1.loglog(
    P_RANGE_ARR[valid], J_arr[valid],
    color='black', lw=STYLE['linewidth'],
    label='L1+L6 Model'
)

# Slope reference lines
if np.any(valid):
    P_ref1  = P_RANGE_ARR[valid][0]
    J_ref1  = J_arr[valid][0]
    ax1.loglog(
        P_RANGE_ARR[valid],
        J_ref1 * (P_RANGE_ARR[valid] / P_ref1) ** 1.0,
        color=CURVE_STYLES['slope_1']['color'],
        ls=CURVE_STYLES['slope_1']['ls'],
        lw=CURVE_STYLES['slope_1']['lw'],
        alpha=CURVE_STYLES['slope_1']['alpha'],
        label='Slope = 1 (surface)'
    )
    P_ref05 = P_RANGE_ARR[valid][-1]
    J_ref05 = J_arr[valid][-1]
    ax1.loglog(
        P_RANGE_ARR[valid],
        J_ref05 * (P_RANGE_ARR[valid] / P_ref05) ** 0.5,
        color=CURVE_STYLES['slope_05']['color'],
        ls=CURVE_STYLES['slope_05']['ls'],
        lw=CURVE_STYLES['slope_05']['lw'],
        alpha=CURVE_STYLES['slope_05']['alpha'],
        label='Slope = 0.5 (metal)'
    )

# Colored regime segments + slope annotation boxes (widget style)
region_styles = [
    ('surface', CURVE_STYLES['surface_region']['color'], 'Surface-limited'),
    ('metal',   CURVE_STYLES['metal_region']['color'],   'Metal-limited'),
    ('mixed',   CURVE_STYLES['mixed_region']['color'],   'Mixed'),
]

for regime, color, label in region_styles:
    mask = (rl_arr == regime) & valid
    if not np.any(mask):
        continue
    idxs = np.where(mask)[0]
    for _, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
        group = list(map(itemgetter(1), g))
        if len(group) < 3:
            continue
        P_seg = P_RANGE_ARR[group]
        J_seg = J_arr[group]
        ax1.loglog(P_seg, J_seg, color=color, lw=4, alpha=0.7)
        slope_seg, *_ = stats.linregress(
            np.log10(P_seg), np.log10(np.abs(J_seg))
        )
        mid = len(group) // 2
        ax1.text(
            P_seg[mid], J_seg[mid],
            f'{label}\nSlope={slope_seg:.2f}',
            color=color,
            fontsize=STYLE['fontsize_annotation'],
            fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

# Net slope annotation
log_P_v = np.log10(P_RANGE_ARR[valid])
log_J_v = np.log10(np.abs(J_arr[valid]))
slope_net, *_ = stats.linregress(log_P_v, log_J_v)
ax1.text(
    0.98, 0.02, f'Net slope = {slope_net:.2f}',
    transform=ax1.transAxes, ha='right', va='bottom',
    fontsize=STYLE['fontsize_annotation'], fontweight='bold',
    bbox=dict(boxstyle='square', fc='wheat', ec='gray', alpha=1)
)
ax1.text(
    0.05, 0.95,
    r'$J = \beta\left(g(\theta) - \sqrt{P_{down}}\right)$',
    transform=ax1.transAxes, fontsize=10, va='top',
    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
)

ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux $J_{ss}$ (mol/m²/s)',         fontsize=STYLE['fontsize_axis'])
ax1.set_title('(A) Flux vs Pressure',               fontsize=STYLE['fontsize_title'])
ax1.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper left')
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (B) SURFACE COVERAGE θ vs PRESSURE
# ─────────────────────────────────────────────────────────────────────────────
ax2 = axes[0, 1]

# Steady-state θ from compute loop
ax2.semilogx(
    P_RANGE_ARR[valid], theta_arr[valid],
    color=CURVE_STYLES['L1_L6']['color'],
    ls=CURVE_STYLES['L1_L6']['ls'],
    lw=STYLE['linewidth'],
    marker=CURVE_STYLES['L1_L6']['marker'],
    ms=CURVE_STYLES['L1_L6']['ms'],
    markevery=8,
    label='θ steady-state'
)

# Equilibrium Langmuir isotherm overlay
theta_eq = (np.sqrt(K_eq * P_RANGE_ARR) /
            (1 + np.sqrt(K_eq * P_RANGE_ARR)))
ax2.semilogx(
    P_RANGE_ARR, theta_eq,
    color='gray', ls='-.', lw=1.5, alpha=0.7,
    label=r'θ$_{eq}$ (Langmuir isotherm)'
)

# Reference lines
ax2.axhline(
    1.0, color='red', ls='--', lw=1.5, alpha=0.5,
    label='θ = 1 (saturation)'
)
ax2.axhline(
    0.5, color='orange', ls='--', lw=1.5, alpha=0.5,
    label='θ = 0.5'
)

# Crossover line
if P_crossover is not None:
    ax2.axvline(
        P_crossover, color='gray', ls=':', lw=1.5, alpha=0.7,
        label=f'Crossover  P={P_crossover:.1e} Pa'
    )

ax2.text(
    0.05, 0.95,
    r'$\theta_{eq} = \frac{\sqrt{K_{eq} P}}{1 + \sqrt{K_{eq} P}}$',
    transform=ax2.transAxes, fontsize=10, va='top',
    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
)
ax2.text(
    0.95, 0.30,
    'θ_ss < θ_eq → surface depleted\n   by fast metal diffusion\n'
    'θ_ss ≈ θ_eq → local equilibrium\n   (metal-limited)',
    transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='top', ha='right',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
)

ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)',  fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('Surface Coverage $\\theta$',        fontsize=STYLE['fontsize_axis'])
ax2.set_title('(B) Surface Coverage vs Pressure',   fontsize=STYLE['fontsize_title'])
ax2.set_ylim(-0.05, 1.1)
ax2.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (C) RESISTANCE FRACTIONS vs PRESSURE
# ─────────────────────────────────────────────────────────────────────────────
ax3 = axes[1, 0]

ax3.semilogx(
    P_RANGE_ARR[valid], frac_surface_arr[valid] * 100,
    color=CURVE_STYLES['surface_region']['color'],
    lw=STYLE['linewidth'], label='Surface (dissociation)'
)
ax3.semilogx(
    P_RANGE_ARR[valid], frac_metal_arr[valid] * 100,
    color=CURVE_STYLES['metal_region']['color'],
    lw=STYLE['linewidth'], label='Metal (diffusion)'
)

# 50% threshold
ax3.axhline(
    50,
    color=CURVE_STYLES['threshold']['color'],
    ls=CURVE_STYLES['threshold']['ls'],
    lw=CURVE_STYLES['threshold']['lw'],
    alpha=CURVE_STYLES['threshold']['alpha'],
    label='50% threshold'
)

# Crossover annotation
if P_crossover is not None:
    ax3.axvline(P_crossover, color='gray', ls=':', lw=1.5, alpha=0.7)
    ax3.text(
        P_crossover * 1.5, 53,
        f'P = {P_crossover:.1e} Pa',
        fontsize=STYLE['fontsize_annotation']-1, color='gray'
    )
else:
    dominant = 'Surface-limited' if n_surface > n_metal else 'Metal-limited'
    ax3.text(
        0.05, 0.50,
        f'No crossover in range\n{dominant} throughout',
        transform=ax3.transAxes,
        fontsize=STYLE['fontsize_annotation'],
        va='center',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )

ax3.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Resistance Fraction (%)',          fontsize=STYLE['fontsize_axis'])
ax3.set_title('(C) Rate-Limiting Step Analysis',   fontsize=STYLE['fontsize_title'])
ax3.set_ylim(0, 105)
ax3.legend(fontsize=STYLE['fontsize_legend'], loc='center right')
ax3.grid(True, alpha=STYLE['grid_alpha'])
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (D) LIMIT CHECK — Sieverts Recovery as k_diss → ∞
# ─────────────────────────────────────────────────────────────────────────────
ax4 = axes[1, 1]

# Sieverts reference (L1, no surface kinetics)
J_sieverts = beta * (np.sqrt(P_RANGE_ARR) - np.sqrt(P_down))

# Sweep k_diss over many orders of magnitude
k_diss_sweep = np.logspace(-15, -3, 7)

# Compute L1+L6 at the highest k_diss only — for parity plot
J_high_k = []
for P_up in P_RANGE_ARR:
    try:
        r = solve_steady_state_flux_L1L6(
            P_up, P_down, L_m,
            k_diss_sweep[-1], K_eq, D_m, K_s_m
        )
        J_high_k.append(r['J_ss'])
    except Exception:
        J_high_k.append(np.nan)
J_high_k = np.array(J_high_k)

valid_parity = ~np.isnan(J_high_k) & (J_sieverts > 0) & (J_high_k > 0)

# Parity plot
ax4.loglog(
    J_sieverts[valid_parity], J_high_k[valid_parity],
    'o',
    color=CURVE_STYLES['L1_L6']['color'],
    markersize=6, alpha=0.8,
    label=f'L1+L6  (k_diss = {k_diss_sweep[-1]:.0e})'
)

# Perfect agreement line
J_min = min(J_sieverts[valid_parity].min(), J_high_k[valid_parity].min()) * 0.5
J_max = max(J_sieverts[valid_parity].max(), J_high_k[valid_parity].max()) * 2.0
ax4.loglog(
    [J_min, J_max], [J_min, J_max],
    color=CURVE_STYLES['parity']['color'],
    ls=CURVE_STYLES['parity']['ls'],
    lw=CURVE_STYLES['parity']['lw'],
    alpha=CURVE_STYLES['parity']['alpha'],
    label='Perfect agreement'
)

# Max relative error
rel_errors = (np.abs(J_high_k[valid_parity] - J_sieverts[valid_parity])
              / J_sieverts[valid_parity] * 100)
max_err    = rel_errors.max()
mean_err   = rel_errors.mean()
box_color  = 'lightgreen' if max_err < 1.0 else 'lightyellow'

ax4.text(
    0.05, 0.95,
    f'k_diss = {k_diss_sweep[-1]:.0e}\n'
    f'Max error  = {max_err:.2e}%\n'
    f'Mean error = {mean_err:.2e}%\n'
    f'{"✓ VALIDATED" if max_err < 1.0 else "⚠ CHECK k_diss"}',
    transform=ax4.transAxes,
    fontsize=STYLE['fontsize_annotation'],
    va='top',
    bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.9)
)

ax4.set_xlabel('Sieverts Flux $J_{L1}$ (mol/m²/s)',         fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('L1+L6 Flux $J_{L1+L6}$ (mol/m²/s)',         fontsize=STYLE['fontsize_axis'])
ax4.set_title(r'(D) Limit Check: $k_{diss} \to \infty$ → Sieverts', fontsize=STYLE['fontsize_title'])
ax4.set_aspect('equal', adjustable='box')
ax4.legend(fontsize=STYLE['fontsize_legend'], loc='lower right')
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# ── Validation summary ───────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("L1+L6 FIGURE 1 — VALIDATION SUMMARY")
print("=" * 60)
print(f"  Net log-log slope:       {slope_net:.4f}")
print(f"  θ range:                 "
      f"{np.nanmin(theta_arr):.4f} → {np.nanmax(theta_arr):.4f}")
print(f"  Sieverts recovery error: {max_err:.2e}%  "
      f"({'✓' if max_err < 1.0 else '⚠'})")
if P_crossover is not None:
    print(f"  Crossover pressure:      {P_crossover:.2e} Pa")
else:
    dominant = 'surface-limited' if n_surface > n_metal else 'metal-limited'
    print(f"  No crossover — system is {dominant} across full range")
print("=" * 60)
```

```python
# =============================================================================
# CELL 5 — FIGURE 2: Extended Analysis (1×2)
# (E) Arrhenius: apparent Ea per regime vs known reference
# (F) Rate-limiting map: surface fraction vs pressure at multiple temperatures
# =============================================================================

# =============================================================================
"""
Figure 2 examines how temperature changes the balance between surface and
metal resistance. While Figure 1 is at a single operating temperature,
Figure 2 asks what happens across the full temperature range.

Panel (E) — Arrhenius: apparent Ea per regime
    Question: Does the model extract the correct activation energy in each
    regime, and does it match the known reference values from config?
    Sweeps temperature from T_RANGE[0] to T_RANGE[1]. At each temperature
    all four properties (k_diss, K_eq, D_m, K_s_m) are recomputed via the
    Arrhenius helper functions from Cell 1. The solver is called at two fixed
    pressures — P_for_surface_Ea and P_for_metal_Ea — selected in Cell 3 as
    the median of each dominant region (fraction > 90%). A linear fit to
    ln(J) vs 1000/T gives the apparent activation energy for each regime.
    Extracted Ea is compared to the known reference from config:
      - Surface regime: extracted Ea vs E_diss (config surface_kinetics)
      - Metal regime:   extracted Ea vs E_D + H_s (config metal transport)
    Dual x-axis: 1000/T on the bottom (K-1), degrees C on top.
    If a regime has insufficient dominant points (P set to None in Cell 3),
    that fit is skipped cleanly with a printed message — no crash.

Panel (F) — Rate-limiting map: surface fraction vs pressure at multiple T
    Question: How does the balance between surface and metal resistance shift
    as temperature increases, and in which direction does the crossover move?
    Runs the full pressure sweep at 5 temperatures spread across T_RANGE.
    Plots surface resistance fraction (%) vs pressure for each temperature
    using a plasma colormap — cool for low T, warm for high T. The direction
    of the crossover shift reveals whether surface kinetics or metal diffusion
    activates faster with temperature, determined by comparing E_diss vs
    E_D + H_s. An annotation box states the observed direction and the
    physical reason derived from config activation energies.
"""




fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
fig2.suptitle(
    f'L1+L6: Extended Analysis — Temperature Effects\n'
    f'{METAL_KEY}  |  L = {L_m*1e3:.1f} mm',
    fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=1.01
)

# # ─────────────────────────────────────────────────────────────────────────────
# # (E) ARRHENIUS — apparent Ea per regime
# # ─────────────────────────────────────────────────────────────────────────────
# ax_E = axes2[0]

# temperatures_K = np.linspace(T_RANGE[0], T_RANGE[1], N_T)
# inv_T          = 1000.0 / temperatures_K

# # ── Compute J at each regime pressure across temperature sweep ────────────────
# results = {}

# for label, P_fixed, color in [
#     ('surface', P_for_surface_Ea, CURVE_STYLES['surface_region']['color']),
#     ('metal',   P_for_metal_Ea,   CURVE_STYLES['metal_region']['color']),
# ]:
#     if P_fixed is None:
#         print(f"  ⚠  {label}-regime: P_for_{label}_Ea is None — fit skipped")
#         results[label] = None
#         continue

#     J_row = []
#     for T_K in temperatures_K:
#         try:
#             r = solve_steady_state_flux_L1L6(
#                 P_fixed, P_down, L_m,
#                 get_k_diss_at_T(T_K),
#                 get_K_eq_at_T(T_K),
#                 get_D_m_at_T(T_K),
#                 get_K_s_m_at_T(T_K),
#             )
#             J_row.append(r['J_ss'])
#         except Exception:
#             J_row.append(np.nan)

#     J_arr_T = np.array(J_row)
#     valid_T  = ~np.isnan(J_arr_T) & (J_arr_T > 0)

#     if np.sum(valid_T) < 3:
#         print(f"  ⚠  {label}-regime: only {np.sum(valid_T)} valid points — fit skipped")
#         results[label] = None
#         continue

#     # Arrhenius fit
#     slope, intercept, r_val, *_ = stats.linregress(
#         inv_T[valid_T], np.log(J_arr_T[valid_T])
#     )
#     E_a_extracted = -slope * R * 1000   # J/mol

#     results[label] = {
#         'J':           J_arr_T,
#         'valid':       valid_T,
#         'slope':       slope,
#         'intercept':   intercept,
#         'r_sq':        r_val**2,
#         'E_a':         E_a_extracted,
#         'P_fixed':     P_fixed,
#         'color':       color,
#     }

#     # Plot data points
#     ax_E.semilogy(
#         inv_T[valid_T], J_arr_T[valid_T],
#         color=color, ls='-', lw=STYLE['linewidth'],
#         marker='o' if label == 'surface' else 's',
#         ms=5, markevery=2,
#         label=f'{label.capitalize()}-limited  (P = {P_fixed:.1e} Pa)'
#     )

#     # Plot fitted line
#     ax_E.semilogy(
#         inv_T[valid_T],
#         np.exp(slope * inv_T[valid_T] + intercept),
#         color=color, ls='--', lw=1.5, alpha=0.7
#     )

# # ── Dual x-axis: °C on top ────────────────────────────────────────────────────
# ax_E_top = ax_E.twiny()
# ax_E_top.set_xlim(ax_E.get_xlim())
# T_ticks_C  = np.array([200, 300, 400, 500, 600, 700, 800, 900])
# T_ticks_K  = T_ticks_C + 273.15
# in_range   = (T_ticks_K >= T_RANGE[0]) & (T_ticks_K <= T_RANGE[1])
# if np.any(in_range):
#     ax_E_top.set_xticks(1000.0 / T_ticks_K[in_range])
#     ax_E_top.set_xticklabels([f'{t}' for t in T_ticks_C[in_range]])
# ax_E_top.set_xlabel('Temperature (°C)',      fontsize=STYLE['fontsize_axis'])
# ax_E_top.tick_params(labelsize=STYLE['fontsize_tick'])

# # ── Annotation box ────────────────────────────────────────────────────────────
# annot_lines = []

# for label, E_ref, ref_name in [
#     ('surface', E_REF_SURFACE, 'E_diss'),
#     ('metal',   E_REF_METAL,   'E_D + H_s'),
# ]:
#     res = results.get(label)
#     if res is None:
#         annot_lines.append(f'{label.capitalize()}: not extractable')
#         annot_lines.append(f'  Expected ({ref_name}): {E_ref/1000:.1f} kJ/mol')
#     else:
#         E_a = res['E_a']
#         match = abs(E_a - E_ref) < 3000    # within 3 kJ/mol
#         annot_lines.append(
#             f'{label.capitalize()}-limited:'
#         )
#         annot_lines.append(
#             f'  Extracted:  {E_a/1000:.1f} kJ/mol'
#         )
#         annot_lines.append(
#             f'  Expected ({ref_name}): {E_ref/1000:.1f} kJ/mol'
#         )
#         annot_lines.append(
#             f'  R² = {res["r_sq"]:.6f}  '
#             f'{"✓" if match else "⚠ check regime"}'
#         )
#     annot_lines.append('')   # blank line between blocks

# ax_E.text(
#     0.97, 0.97, '\n'.join(annot_lines).rstrip(),
#     transform=ax_E.transAxes,
#     fontsize=STYLE['fontsize_annotation'],
#     va='top', ha='right',
#     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9)
# )

# ax_E.set_xlabel('1000/T (K⁻¹)',              fontsize=STYLE['fontsize_axis'])
# ax_E.set_ylabel('Flux $J_{ss}$ (mol/m²/s)',  fontsize=STYLE['fontsize_axis'])
# ax_E.set_title('(E) Arrhenius: Apparent $E_a$ per Regime',
#                fontsize=STYLE['fontsize_title'])
# ax_E.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower left')
# ax_E.grid(True, which='both', alpha=STYLE['grid_alpha'])
# ax_E.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (E) ARRHENIUS — apparent Ea per regime, classified by regime purity
# ─────────────────────────────────────────────────────────────────────────────
ax_E = axes2[0]

temperatures_K = np.linspace(T_RANGE[0], T_RANGE[1], N_T)
inv_T          = 1000.0 / temperatures_K

# ── Regime classification — two resistances only ──────────────────────────────
def classify_regime_L1(fs, fm, pure_thresh=0.90, mixed_thresh=0.50):
    if fs >= pure_thresh:    return 'pure_surface'
    elif fm >= pure_thresh:  return 'pure_metal'
    elif fs >= mixed_thresh: return 'mixed_surface'
    elif fm >= mixed_thresh: return 'mixed_metal'
    else:                    return 'mixed'

regime_class_styles = {
    'pure_surface':  {'color': CURVE_STYLES['surface_region']['color'], 'marker': 'o', 'ms': 8,  'alpha': 1.0},
    'mixed_surface': {'color': CURVE_STYLES['surface_region']['color'], 'marker': 'o', 'ms': 5,  'alpha': 0.35},
    'pure_metal':    {'color': CURVE_STYLES['metal_region']['color'],   'marker': '^', 'ms': 8,  'alpha': 1.0},
    'mixed_metal':   {'color': CURVE_STYLES['metal_region']['color'],   'marker': '^', 'ms': 5,  'alpha': 0.35},
    'mixed':         {'color': CURVE_STYLES['mixed_region']['color'],   'marker': 'x', 'ms': 6,  'alpha': 0.5},
}

regime_Ea_cases = [
    ('surface', P_for_surface_Ea, CURVE_STYLES['surface_region']['color']),
    ('metal',   P_for_metal_Ea,   CURVE_STYLES['metal_region']['color']),
]

results            = {}
plotted_cls_labels = set()

for label, P_fixed, color in regime_Ea_cases:
    if P_fixed is None:
        print(f"  ⚠  {label}: P_for_{label}_Ea is None — fit skipped")
        results[label] = None
        continue

    J_row     = []
    class_row = []

    for T_K in temperatures_K:
        try:
            r = solve_steady_state_flux_L1L6(
                P_fixed, P_down, L_m,
                get_k_diss_at_T(T_K),
                get_K_eq_at_T(T_K),
                get_D_m_at_T(T_K),
                get_K_s_m_at_T(T_K),
            )
            J_row.append(r['J_ss'])
            fs = r['resistances']['fraction_surface']
            fm = r['resistances']['fraction_metal']
            class_row.append(classify_regime_L1(fs, fm))
        except Exception:
            J_row.append(np.nan)
            class_row.append('mixed')

    J_arr_T   = np.array(J_row)
    class_arr = np.array(class_row)
    pure      = (class_arr == f'pure_{label}') & ~np.isnan(J_arr_T) & (J_arr_T > 0)

    # Plot all classified points
    for cls, sty in regime_class_styles.items():
        mask = (class_arr == cls) & ~np.isnan(J_arr_T) & (J_arr_T > 0)
        if not np.any(mask):
            continue
        cls_label = cls.replace('_', ' ').capitalize()
        ax_E.semilogy(
            inv_T[mask], J_arr_T[mask],
            ls='none',
            color=sty['color'], marker=sty['marker'],
            ms=sty['ms'], alpha=sty['alpha'],
            label=cls_label if cls_label not in plotted_cls_labels
                  else '_nolegend_'
        )
        plotted_cls_labels.add(cls_label)

    # Fit on pure points only
    if np.sum(pure) < 3:
        print(f"  ⚠  {label}: only {np.sum(pure)} pure T points — fit skipped")
        results[label] = None
        continue

    slope, intercept, r_val, *_ = stats.linregress(
        inv_T[pure], np.log(J_arr_T[pure])
    )
    E_a = -slope * R * 1000

    results[label] = {
        'J': J_arr_T, 'valid': pure,
        'slope': slope, 'intercept': intercept,
        'r_sq': r_val**2, 'E_a': E_a,
        'P_fixed': P_fixed, 'color': color,
        'n_pure': int(np.sum(pure)),
        'T_pure_range': (temperatures_K[pure].min() - 273.15,
                         temperatures_K[pure].max() - 273.15),
    }

    ax_E.semilogy(
        inv_T[pure],
        np.exp(slope * inv_T[pure] + intercept),
        color=color, ls='--', lw=2.0, alpha=0.9,
        label=f'Fit {label}  (pure pts)'
    )

# Dual x-axis
ax_E_top = ax_E.twiny()
ax_E_top.set_xlim(ax_E.get_xlim())
T_ticks_C = np.array([200, 300, 400, 500, 600, 700, 800, 900])
T_ticks_K = T_ticks_C + 273.15
in_range  = (T_ticks_K >= T_RANGE[0]) & (T_ticks_K <= T_RANGE[1])
if np.any(in_range):
    ax_E_top.set_xticks(1000.0 / T_ticks_K[in_range])
    ax_E_top.set_xticklabels([f'{t}' for t in T_ticks_C[in_range]])
ax_E_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax_E_top.tick_params(labelsize=STYLE['fontsize_tick'])

# Annotation box
annot = []
for lbl, E_ref, ref_name in [
    ('surface', E_REF_SURFACE, 'E_diss'),
    ('metal',   E_REF_METAL,   'E_D + H_s'),
]:
    res = results.get(lbl)
    if res is None:
        annot.append(f'{lbl.capitalize()}: not extractable')
        annot.append(f'  Expected ({ref_name}): {E_ref/1000:.1f} kJ/mol')
    else:
        match = abs(res['E_a'] - E_ref) < 3000
        annot.append(f'{lbl.capitalize()}-limited:')
        annot.append(f'  Extracted: {res["E_a"]/1000:.1f} kJ/mol')
        annot.append(f'  Expected ({ref_name}): {E_ref/1000:.1f} kJ/mol')
        annot.append(f'  R²={res["r_sq"]:.4f}  {"✓" if match else "⚠"}')
        annot.append(f'  T range: {res["T_pure_range"][0]:.0f}–'
                     f'{res["T_pure_range"][1]:.0f}°C  '
                     f'({res["n_pure"]} pure pts)')
    annot.append('')

ax_E.text(
    0.97, 0.97, '\n'.join(annot).rstrip(),
    transform=ax_E.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='top', ha='right',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9)
)

ax_E.set_xlabel('1000/T (K⁻¹)',             fontsize=STYLE['fontsize_axis'])
ax_E.set_ylabel('Flux $J_{ss}$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax_E.set_title('(E) Arrhenius: Apparent $E_a$ per Regime',
               fontsize=STYLE['fontsize_title'])
ax_E.legend(fontsize=STYLE['fontsize_legend']-2, loc='lower left',
            ncol=2, title='Symbol size = purity')
ax_E.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax_E.tick_params(labelsize=STYLE['fontsize_tick'])
# ─────────────────────────────────────────────────────────────────────────────
# (F) RATE-LIMITING MAP — surface fraction vs pressure at 5 temperatures
# ─────────────────────────────────────────────────────────────────────────────
ax_F = axes2[1]

T_map_values = np.linspace(T_RANGE[0], T_RANGE[1], 5)
T_colors     = plt.cm.plasma(np.linspace(0.1, 0.9, len(T_map_values)))

crossover_pressures = []   # collect crossover per T for annotation

for T_K, col in zip(T_map_values, T_colors):
    fs_row = []
    for P_up in P_RANGE_ARR:
        try:
            r = solve_steady_state_flux_L1L6(
                P_up, P_down, L_m,
                get_k_diss_at_T(T_K),
                get_K_eq_at_T(T_K),
                get_D_m_at_T(T_K),
                get_K_s_m_at_T(T_K),
            )
            fs_row.append(r['resistances']['fraction_surface'] * 100)
        except Exception:
            fs_row.append(np.nan)

    fs_row    = np.array(fs_row)
    valid_row = ~np.isnan(fs_row)

    if not np.any(valid_row):
        continue

    ax_F.semilogx(
        P_RANGE_ARR[valid_row], fs_row[valid_row],
        color=col, lw=STYLE['linewidth'],
        label=f'T = {T_K-273.15:.0f}°C'
    )

    # Find crossover for this temperature
    sc = np.where(np.diff(np.sign(fs_row[valid_row] - 50)))[0]
    if len(sc) > 0:
        idx = sc[0]
        fs_v = fs_row[valid_row]
        if idx + 1 < len(fs_v) and fs_v[idx] > 50 and fs_v[idx+1] < 50:
            crossover_pressures.append((T_K, P_RANGE_ARR[valid_row][idx]))

# 50% threshold
ax_F.axhline(
    50,
    color=CURVE_STYLES['threshold']['color'],
    ls=CURVE_STYLES['threshold']['ls'],
    lw=CURVE_STYLES['threshold']['lw'],
    alpha=CURVE_STYLES['threshold']['alpha'],
    label='50% threshold'
)

# ── Crossover shift annotation ────────────────────────────────────────────────
if len(crossover_pressures) >= 2:
    T_low,  P_co_low  = crossover_pressures[0]
    T_high, P_co_high = crossover_pressures[-1]

    if P_co_high < P_co_low:
        direction = 'lower P as T ↑'
        physical  = f'E_diss ({E_REF_SURFACE/1000:.0f}) > E_D+Hs ({E_REF_METAL/1000:.0f})'
        meaning   = 'Surface activates faster → less limiting at high T'
    else:
        direction = 'higher P as T ↑'
        physical  = f'E_D+Hs ({E_REF_METAL/1000:.0f}) > E_diss ({E_REF_SURFACE/1000:.0f})'
        meaning   = 'Metal activates faster → more limiting at high T'

    ax_F.text(
        0.05, 0.05,
        f'Crossover shifts to {direction}\n{physical}\n{meaning}',
        transform=ax_F.transAxes,
        fontsize=STYLE['fontsize_annotation']-1,
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )
elif len(crossover_pressures) == 0:
    dominant = 'surface-limited' if n_surface > n_metal else 'metal-limited'
    ax_F.text(
        0.05, 0.05,
        f'No crossover at any T\nSystem is {dominant}\nacross full pressure range',
        transform=ax_F.transAxes,
        fontsize=STYLE['fontsize_annotation']-1,
        va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )

ax_F.set_xlabel('Upstream Pressure $P_{up}$ (Pa)',  fontsize=STYLE['fontsize_axis'])
ax_F.set_ylabel('Surface Resistance Fraction (%)',  fontsize=STYLE['fontsize_axis'])
ax_F.set_title('(F) Rate-Limiting Map vs Temperature',
               fontsize=STYLE['fontsize_title'])
ax_F.set_ylim(0, 105)
ax_F.legend(
    fontsize=STYLE['fontsize_legend']-1,
    loc='upper right',
    title='Temperature'
)
ax_F.grid(True, alpha=STYLE['grid_alpha'])
ax_F.tick_params(labelsize=STYLE['fontsize_tick'])

plt.tight_layout()
plt.show()

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("L1+L6 FIGURE 2 — ARRHENIUS SUMMARY")
print("=" * 60)

for label, E_ref, ref_name in [
    ('surface', E_REF_SURFACE, 'E_diss'),
    ('metal',   E_REF_METAL,   'E_D + H_s'),
]:
    res = results.get(label)
    if res is None:
        print(f"  {label.capitalize()}-limited:  not extractable")
    else:
        match = abs(res['E_a'] - E_ref) < 3000
        print(f"  {label.capitalize()}-limited:")
        print(f"    Extracted  Ea = {res['E_a']/1000:.1f} kJ/mol")
        print(f"    Expected ({ref_name}) = {E_ref/1000:.1f} kJ/mol")
        print(f"    R²         = {res['r_sq']:.6f}")
        print(f"    Match:       {'✓' if match else '⚠ — check regime purity'}")

if crossover_pressures:
    print(f"\n  Crossover pressures across T sweep:")
    for T_K, P_co in crossover_pressures:
        print(f"    T = {T_K-273.15:.0f}°C  →  P_crossover = {P_co:.2e} Pa")

print("=" * 60)
```

```python
# =============================================================================
# CELL 6 — DATAFRAME OUTPUT
# =============================================================================
df = pd.DataFrame(rows)

# ── Formatting ────────────────────────────────────────────────────────────────
df_display = df.copy()

for col in df_display.columns:
    if col == 'Rate-Limiting' or col == 'Error':
        continue
    if df_display[col].dtype == float:
        if 'fraction' in col.lower():
            df_display[col] = df_display[col].round(2)
        elif 'theta' in col.lower():
            df_display[col] = df_display[col].round(6)
        else:
            df_display[col] = df_display[col].apply(
                lambda x: f'{x:.4e}' if pd.notna(x) else 'NaN'
            )

# ── Summary line ──────────────────────────────────────────────────────────────
n_surface_df = (df['Rate-Limiting'] == 'SURFACE').sum()
n_metal_df   = (df['Rate-Limiting'] == 'METAL').sum()
n_mixed_df   = (df['Rate-Limiting'] == 'MIXED').sum()
n_error_df   = (df['Rate-Limiting'] == 'ERROR').sum()

print("=" * 60)
print("L1+L6 RESULTS DATAFRAME")
print("=" * 60)
print(f"  Metal:       {METAL_KEY}")
print(f"  T:           {T_op-273.15:.0f}°C  ({T_op} K)")
print(f"  L_m:         {L_m*1e3:.1f} mm")
print(f"  k_diss:      {k_diss:.3e}  mol/m²/s/Pa")
print(f"  K_eq:        {K_eq:.3e}  Pa⁻¹")
print(f"  D_m:         {D_m:.3e}  m²/s")
print(f"  K_s_m:       {K_s_m:.3e}  mol/m³/Pa^0.5")
print(f"\n  Total rows:        {len(df)}")
print(f"  Surface-limited:   {n_surface_df}")
print(f"  Metal-limited:     {n_metal_df}")
print(f"  Mixed:             {n_mixed_df}")
if n_error_df > 0:
    print(f"  Errors:            {n_error_df}  ⚠")
print("=" * 60)

display(df_display)
```

### L2a+L6: Perfect Oxide + Surface Kinetics 

```python
# =============================================================================
# CELL 1 — SETUP
# Run once. All imports and helper functions.
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import brentq
import pandas as pd
from itertools import groupby
from operator import itemgetter
from IPython.display import display

from calculations.surface_kinetics import solve_steady_state_flux_L2aL6
from calculations.config.model_config import (
    R, F,
    METALS, OXIDES, CONDITIONS, VALIDATION,
    PLOT_STYLE as STYLE,
    CURVE_STYLES,
    apply_style,
    build_simulation_config,
)

# ── Simulation config ─────────────────────────────────────────────────────────
SIM      = build_simulation_config(
    metal='metal_316L_Heat_treated_ref_cast',
    oxide='Al2O3',
    T_operating=873,
    P_upstream=1e5,
    L_metal=1e-3,
)

OXIDE_KEY = SIM['oxide_name']
OXIDE     = OXIDES[OXIDE_KEY]
SK        = OXIDE['surface_kinetics']     # oxide surface kinetics sub-dict

# ── Pressure sweep (matches widget range) ─────────────────────────────────────
P_RANGE_ARR = np.logspace(np.log10(CONDITIONS['P_range'][0]), np.log10(CONDITIONS['P_range'][1]), CONDITIONS['n_P_points'])
P_DOWN      = SIM['P_downstream']
T_RANGE     = CONDITIONS['T_range']
N_T         = CONDITIONS['n_T_points']

# ── Arrhenius helper functions ────────────────────────────────────────────────
def get_k_diss_at_T(T_K):
    """k_diss at temperature T_K via Arrhenius from config."""
    return SK['k_diss_ref'] * np.exp(
        (-SK['E_diss'] / R) * (1/T_K - 1/SK['T_ref'])
    )

def get_K_eq_at_T(T_K):
    """K_eq at temperature T_K via van't Hoff from config."""
    return SK['K_eq_ref'] * np.exp(
        (-SK['H_eq'] / R) * (1/T_K - 1/SK['T_ref'])
    )

def get_D_ox_at_T(T_K):
    """D_ox at temperature T_K via Arrhenius from config."""
    return OXIDE['D_ox_ref'] * np.exp(
        (-OXIDE['E_D_ox'] / R) * (1/T_K - 1/OXIDE['T_ref'])
    )

def get_K_ox_at_T(T_K):
    """K_ox at temperature T_K via Arrhenius from config."""
    return OXIDE['K_ox_ref'] * np.exp(
        (-OXIDE['H_sol_ox'] / R) * (1/T_K - 1/OXIDE['T_ref'])
    )

# ── Known reference activation energies ──────────────────────────────────────
E_REF_SURFACE = SK['E_diss']                          # J/mol — surface dissociation
E_REF_OXIDE   = OXIDE['E_D_ox'] + OXIDE['H_sol_ox']  # J/mol — oxide permeation

print("=" * 60)
print("SETUP COMPLETE")
print("=" * 60)
print(f"  Oxide:           {OXIDE_KEY}")
print(f"  E_ref surface:   {E_REF_SURFACE/1000:.1f} kJ/mol  (E_diss)")
print(f"  E_ref oxide:     {E_REF_OXIDE/1000:.1f} kJ/mol  (E_D_ox + H_sol_ox)")
print(f"  P range:         {P_RANGE_ARR[0]:.0e} – {P_RANGE_ARR[-1]:.0e} Pa")
print(f"  T range:         {T_RANGE[0]} – {T_RANGE[1]} K")
print("=" * 60)
```

```python
# =============================================================================
# CELL 2 — PARAMETERS
# All physics values at operating temperature, pulled from config.
# Change placeholder values in config and rerun from this cell downward.
# =============================================================================

# ── Operating conditions ──────────────────────────────────────────────────────
T_op   = SIM['T_operating']
L_ox   = SIM['L_oxide']
P_down = SIM['P_downstream']

# ── Surface kinetics at T_op ──────────────────────────────────────────────────
k_diss   = get_k_diss_at_T(T_op)
K_eq     = get_K_eq_at_T(T_op)
k_recomb = k_diss / K_eq

# ── Oxide transport at T_op ───────────────────────────────────────────────────
D_ox = get_D_ox_at_T(T_op)
K_ox = get_K_ox_at_T(T_op)

# ── Derived ───────────────────────────────────────────────────────────────────
alpha = D_ox * K_ox / L_ox    # oxide permeance [mol/m²/s/Pa^0.5]

# ── Print verification ────────────────────────────────────────────────────────
print("=" * 60)
print(f"L2a+L6 PARAMETERS AT T = {T_op-273.15:.0f}°C  ({T_op} K)")
print("=" * 60)

print(f"\n  Geometry")
print(f"    L_ox     = {L_ox:.3e}  m")
print(f"    P_down   = {P_down:.1e}  Pa")

print(f"\n  Surface kinetics  [config: oxide surface_kinetics sub-dict]")
print(f"    k_diss   = {k_diss:.3e}  mol/m²/s/Pa")
print(f"    K_eq     = {K_eq:.3e}  Pa⁻¹")
print(f"    k_recomb = {k_recomb:.3e}  mol/m²/s  (= k_diss / K_eq)")

print(f"\n  Oxide transport  [config: Arrhenius]")
print(f"    D_ox     = {D_ox:.3e}  m²/s")
print(f"    K_ox     = {K_ox:.3e}  mol/m³/Pa^0.5")

print(f"\n  Derived")
print(f"    α        = {alpha:.3e}  mol/m²/s/Pa^0.5  (= D_ox × K_ox / L_ox)")

print(f"\n  Reference activation energies  [config]")
print(f"    E_diss        = {E_REF_SURFACE/1000:.1f} kJ/mol   (surface dissociation)")
print(f"    E_D_ox+H_sol  = {E_REF_OXIDE/1000:.1f} kJ/mol   (oxide permeation)")

if SK['k_diss_ref'] <= 1e-12:
    print(f"\n  ⚠  WARNING: k_diss_ref = {SK['k_diss_ref']:.0e} looks like a placeholder.")
    print(f"     Update config oxide surface_kinetics sub-dict with literature values.")

print("\n" + "=" * 60)
```

```python
# =============================================================================
# CELL 3 — COMPUTE
# Single-pass loop over full pressure range.
# Builds all arrays, regime summary, and Arrhenius pressure selection.
# =============================================================================

# ── Single-pass loop ──────────────────────────────────────────────────────────
rows             = []
J_arr            = []
theta_arr        = []
P_surf_arr       = []
frac_surface_arr = []
frac_oxide_arr   = []
R_surface_arr    = []
R_oxide_arr      = []
R_total_arr      = []
rate_lim_arr     = []

for P_up in P_RANGE_ARR:
    try:
        r = solve_steady_state_flux_L2aL6(
            P_up, P_down,
            k_diss, K_eq, D_ox, K_ox, L_ox
        )
        J_arr.append(r['J_ss'])
        theta_arr.append(r['theta'])
        P_surf_arr.append(r['P_surf'])
        frac_surface_arr.append(r['resistances']['fraction_surface'])
        frac_oxide_arr.append(r['resistances']['fraction_oxide'])
        R_surface_arr.append(r['resistances']['R_surface'])
        R_oxide_arr.append(r['resistances']['R_oxide'])
        R_total_arr.append(r['resistances']['R_total'])
        rate_lim_arr.append(r['rate_limiting'])

        rows.append({
            'P_up (Pa)':            P_up,
            'P_surf (Pa)':          r['P_surf'],
            'J_ss (mol/m²/s)':      r['J_ss'],
            'θ_surface':            r['theta'],
            'α_oxide':              r['alpha'],
            'R_surface':            r['resistances']['R_surface'],
            'R_oxide':              r['resistances']['R_oxide'],
            'R_total':              r['resistances']['R_total'],
            'fraction_surface (%)': r['resistances']['fraction_surface'] * 100,
            'fraction_oxide (%)':   r['resistances']['fraction_oxide']   * 100,
            'Rate-Limiting':        r['rate_limiting'].upper(),
        })

    except Exception as e:
        J_arr.append(np.nan);            theta_arr.append(np.nan)
        P_surf_arr.append(np.nan);       frac_surface_arr.append(np.nan)
        frac_oxide_arr.append(np.nan);   R_surface_arr.append(np.nan)
        R_oxide_arr.append(np.nan);      R_total_arr.append(np.nan)
        rate_lim_arr.append('error')
        rows.append({'P_up (Pa)': P_up, 'Rate-Limiting': 'ERROR', 'Error': str(e)})

# ── Convert to arrays ─────────────────────────────────────────────────────────
J_arr            = np.array(J_arr)
theta_arr        = np.array(theta_arr)
P_surf_arr       = np.array(P_surf_arr)
frac_surface_arr = np.array(frac_surface_arr)
frac_oxide_arr   = np.array(frac_oxide_arr)
R_surface_arr    = np.array(R_surface_arr)
R_oxide_arr      = np.array(R_oxide_arr)
R_total_arr      = np.array(R_total_arr)
rate_lim_arr     = np.array(rate_lim_arr)
valid            = ~np.isnan(J_arr)

# ── Rate-limiting classification array ───────────────────────────────────────
rl_arr = np.where(frac_surface_arr > 0.5, 'surface',
         np.where(frac_oxide_arr   > 0.5, 'oxide', 'mixed'))

# ── Regime summary ────────────────────────────────────────────────────────────
n_surface = np.sum((frac_surface_arr > 0.5) & valid)
n_oxide   = np.sum((frac_oxide_arr   > 0.5) & valid)
n_mixed   = np.sum(
    (frac_surface_arr <= 0.5) & (frac_oxide_arr <= 0.5) & valid
)

surf_dom_mask  = (frac_surface_arr > 0.90) & valid
oxide_dom_mask = (frac_oxide_arr   > 0.90) & valid
P_surf_dom     = P_RANGE_ARR[surf_dom_mask]
P_oxide_dom    = P_RANGE_ARR[oxide_dom_mask]

# ── Arrhenius pressure selection (no crossover dependence) ───────────────────
P_for_surface_Ea = np.median(P_surf_dom)  if len(P_surf_dom)  > 3 else None
P_for_oxide_Ea   = np.median(P_oxide_dom) if len(P_oxide_dom) > 3 else None

# ── Crossover pressure (diagnostic only) ─────────────────────────────────────
P_crossover = None
valid_cross = valid & ~np.isnan(frac_surface_arr)
if np.sum(valid_cross) > 2:
    fs_v = frac_surface_arr[valid_cross]
    P_v  = P_RANGE_ARR[valid_cross]
    sc   = np.where(np.diff(np.sign(fs_v - 0.5)))[0]
    if len(sc) > 0:
        idx = sc[0]
        if idx + 1 < len(fs_v) and fs_v[idx] > 0.5 and fs_v[idx+1] < 0.5:
            P_crossover = P_v[idx]

# ── Print regime summary ──────────────────────────────────────────────────────
print("=" * 60)
print(f"L2a+L6 REGIME SUMMARY  (T = {T_op-273.15:.0f}°C, {np.sum(valid)} valid points)")
print("=" * 60)

print(f"\n  Regime breakdown (fraction > 50%)")
print(f"    Surface-limited:  {n_surface:3d} points")
print(f"    Oxide-limited:    {n_oxide:3d} points")
print(f"    Mixed:            {n_mixed:3d} points")

print(f"\n  Dominant regions (fraction > 90%) — used for Arrhenius")
if len(P_surf_dom) > 3:
    print(f"    Surface-dominant: {len(P_surf_dom):3d} points  "
          f"P = {P_surf_dom.min():.1e} – {P_surf_dom.max():.1e} Pa")
    print(f"    → P_for_surface_Ea = {P_for_surface_Ea:.2e} Pa")
else:
    print(f"    Surface-dominant: {len(P_surf_dom):3d} points  "
          f"(insufficient — surface Ea not extractable)")

if len(P_oxide_dom) > 3:
    print(f"    Oxide-dominant:   {len(P_oxide_dom):3d} points  "
          f"P = {P_oxide_dom.min():.1e} – {P_oxide_dom.max():.1e} Pa")
    print(f"    → P_for_oxide_Ea  = {P_for_oxide_Ea:.2e} Pa")
else:
    print(f"    Oxide-dominant:   {len(P_oxide_dom):3d} points  "
          f"(insufficient — oxide Ea not extractable)")

print(f"\n  Crossover pressure (diagnostic only)")
if P_crossover is not None:
    print(f"    P_crossover = {P_crossover:.2e} Pa")
else:
    print(f"    No crossover detected in physical pressure range")
    dominant = 'surface-limited' if n_surface > n_oxide else 'oxide-limited'
    print(f"    System is predominantly {dominant} across full range")

print(f"\n  θ range:  {np.nanmin(theta_arr):.4f} → {np.nanmax(theta_arr):.4f}")
print(f"  J range:  {np.nanmin(J_arr):.2e} → {np.nanmax(J_arr):.2e} mol/m²/s")
print("=" * 60)
```

```python
# =============================================================================
# CELL 4 — FIGURE 1: Core Validation (2×2)
# =============================================================================
"""
Figure 1 validates the L2a+L6 model from four angles at operating temperature.

Panel (A) — Flux vs Pressure
    Question: Does the model produce the correct pressure-dependent behaviour?
    The hero plot. Shows J_ss vs P_up on a log-log scale. The curve should
    transition from slope=1 (surface-limited, low P) to slope=0.5
    (oxide-limited, high P). Note: in the pure L2a limit (no surface
    kinetics), the oxide follows Sieverts-like solubility so the slope is
    always 0.5. The addition of L6 introduces the slope=1 surface-limited
    regime at low P. Style follows the widget — no fill_between. The curve
    is recolored by regime using thick colored segments with slope annotation
    boxes sitting directly on each segment. Reference lines for slope=1 and
    slope=0.5 are drawn as dashed lines. Net slope annotation box in the
    bottom right corner.

Panel (B) — Surface Coverage θ vs Pressure
    Question: How does surface occupancy evolve with pressure, and does it
    match the Langmuir equilibrium isotherm?
    Plots steady-state θ from the compute loop alongside the equilibrium
    Langmuir isotherm θ_eq = sqrt(K_eq*P)/(1+sqrt(K_eq*P)) as a gray
    dash-dot line. In the oxide-limited regime the surface is in local
    equilibrium (θ_ss ≈ θ_eq). In the surface-limited regime the surface
    is being depleted faster than it can be replenished (θ_ss < θ_eq).
    Crossover pressure marked if it exists.

Panel (C) — Resistance Fractions vs Pressure
    Question: Which step controls the rate, and how does the balance shift
    with pressure?
    Plots R_surface/R_total and R_oxide/R_total as percentages vs pressure
    on a semilog scale. Note: there is no metal resistance term here — the
    two fractions always sum to 100%. The 50% threshold line marks the
    crossover. If no crossover exists a text annotation states which regime
    dominates throughout.

Panel (D) — Limit Check: Fast-Kinetics Recovery (L2a)
    Question: As k_diss → ∞, does L2a+L6 recover the pure L2a result?
    Runs the solver at k_diss = 1e-3 (fast kinetics limit) and compares
    J_L2a+L6 against the analytical L2a flux J = alpha*(sqrt(P_up) -
    sqrt(P_down)). A parity plot shows J_L2a+L6 vs J_L2a — perfect
    agreement is a diagonal line. This validates the model correctly
    recovers oxide-only Sieverts diffusion when surface kinetics is fast.
    Green box if max error < 1%, yellow otherwise.
"""

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle(
    f'L2a+L6: Perfect Oxide + Surface Kinetics — Core Validation\n'
    f'{OXIDE_KEY}  |  T = {T_op-273.15:.0f}°C  |  L_ox = {L_ox*1e6:.1f} μm',
    fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98
)

# ─────────────────────────────────────────────────────────────────────────────
# (A) FLUX vs PRESSURE
# ─────────────────────────────────────────────────────────────────────────────
ax1 = axes[0, 0]

ax1.loglog(
    P_RANGE_ARR[valid], J_arr[valid],
    color='black', lw=STYLE['linewidth'], label='L2a+L6 Model'
)

if np.any(valid):
    P_ref1 = P_RANGE_ARR[valid][0]
    J_ref1 = J_arr[valid][0]
    ax1.loglog(
        P_RANGE_ARR[valid],
        J_ref1 * (P_RANGE_ARR[valid] / P_ref1) ** 1.0,
        color=CURVE_STYLES['slope_1']['color'],
        ls=CURVE_STYLES['slope_1']['ls'],
        lw=CURVE_STYLES['slope_1']['lw'],
        alpha=CURVE_STYLES['slope_1']['alpha'],
        label='Slope = 1 (surface)'
    )
    P_ref05 = P_RANGE_ARR[valid][-1]
    J_ref05 = J_arr[valid][-1]
    ax1.loglog(
        P_RANGE_ARR[valid],
        J_ref05 * (P_RANGE_ARR[valid] / P_ref05) ** 0.5,
        color=CURVE_STYLES['slope_05']['color'],
        ls=CURVE_STYLES['slope_05']['ls'],
        lw=CURVE_STYLES['slope_05']['lw'],
        alpha=CURVE_STYLES['slope_05']['alpha'],
        label='Slope = 0.5 (oxide)'
    )

# Colored regime segments + slope annotation boxes
region_styles = [
    ('surface', CURVE_STYLES['surface_region']['color'], 'Surface-limited'),
    ('oxide',   CURVE_STYLES['oxide_region']['color'],   'Oxide-limited'),
    ('mixed',   CURVE_STYLES['mixed_region']['color'],   'Mixed'),
]

for regime, color, label in region_styles:
    mask = (rl_arr == regime) & valid
    if not np.any(mask):
        continue
    idxs = np.where(mask)[0]
    for _, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
        group = list(map(itemgetter(1), g))
        if len(group) < 3:
            continue
        P_seg = P_RANGE_ARR[group]
        J_seg = J_arr[group]
        ax1.loglog(P_seg, J_seg, color=color, lw=4, alpha=0.7)
        slope_seg, *_ = stats.linregress(
            np.log10(P_seg), np.log10(np.abs(J_seg))
        )
        mid = len(group) // 2
        ax1.text(
            P_seg[mid], J_seg[mid],
            f'{label}\nSlope={slope_seg:.2f}',
            color=color,
            fontsize=STYLE['fontsize_annotation'],
            fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

log_P_v   = np.log10(P_RANGE_ARR[valid])
log_J_v   = np.log10(np.abs(J_arr[valid]))
slope_net, *_ = stats.linregress(log_P_v, log_J_v)

ax1.text(
    0.98, 0.02, f'Net slope = {slope_net:.2f}',
    transform=ax1.transAxes, ha='right', va='bottom',
    fontsize=STYLE['fontsize_annotation'], fontweight='bold',
    bbox=dict(boxstyle='square', fc='wheat', ec='gray', alpha=1)
)
ax1.text(
    0.05, 0.95,
    r'$J = \alpha\left(g(\theta) - \sqrt{P_{down}}\right)$',
    transform=ax1.transAxes, fontsize=10, va='top',
    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
)

ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux $J_{ss}$ (mol/m²/s)',         fontsize=STYLE['fontsize_axis'])
ax1.set_title('(A) Flux vs Pressure',               fontsize=STYLE['fontsize_title'])
ax1.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper left')
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (B) SURFACE COVERAGE θ vs PRESSURE
# ─────────────────────────────────────────────────────────────────────────────
ax2 = axes[0, 1]

ax2.semilogx(
    P_RANGE_ARR[valid], theta_arr[valid],
    color=CURVE_STYLES['L2a_L6']['color'],
    ls=CURVE_STYLES['L2a_L6']['ls'],
    lw=STYLE['linewidth'],
    marker=CURVE_STYLES['L2a_L6']['marker'],
    ms=CURVE_STYLES['L2a_L6']['ms'],
    markevery=8,
    label='θ steady-state'
)

theta_eq = (np.sqrt(K_eq * P_RANGE_ARR) /
            (1 + np.sqrt(K_eq * P_RANGE_ARR)))
ax2.semilogx(
    P_RANGE_ARR, theta_eq,
    color='gray', ls='-.', lw=1.5, alpha=0.7,
    label=r'θ$_{eq}$ (Langmuir isotherm)'
)
ax2.axhline(1.0, color='red',    ls='--', lw=1.5, alpha=0.5, label='θ = 1 (saturation)')
ax2.axhline(0.5, color='orange', ls='--', lw=1.5, alpha=0.5, label='θ = 0.5')

if P_crossover is not None:
    ax2.axvline(P_crossover, color='gray', ls=':', lw=1.5, alpha=0.7,
                label=f'Crossover  P={P_crossover:.1e} Pa')

ax2.text(
    0.05, 0.95,
    r'$\theta_{eq} = \frac{\sqrt{K_{eq} P}}{1 + \sqrt{K_{eq} P}}$',
    transform=ax2.transAxes, fontsize=10, va='top',
    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
)
ax2.text(
    0.95, 0.30,
    'θ_ss < θ_eq → surface depleted\n   by fast oxide diffusion\n'
    'θ_ss ≈ θ_eq → local equilibrium\n   (oxide-limited)',
    transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='top', ha='right',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
)

ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)',  fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('Surface Coverage $\\theta$',        fontsize=STYLE['fontsize_axis'])
ax2.set_title('(B) Surface Coverage vs Pressure',   fontsize=STYLE['fontsize_title'])
ax2.set_ylim(-0.05, 1.1)
ax2.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (C) RESISTANCE FRACTIONS vs PRESSURE
# ─────────────────────────────────────────────────────────────────────────────
ax3 = axes[1, 0]

ax3.semilogx(
    P_RANGE_ARR[valid], frac_surface_arr[valid] * 100,
    color=CURVE_STYLES['surface_region']['color'],
    lw=STYLE['linewidth'], label='Surface (dissociation)'
)
ax3.semilogx(
    P_RANGE_ARR[valid], frac_oxide_arr[valid] * 100,
    color=CURVE_STYLES['oxide_region']['color'],
    lw=STYLE['linewidth'], label='Oxide (diffusion)'
)
ax3.axhline(
    50,
    color=CURVE_STYLES['threshold']['color'],
    ls=CURVE_STYLES['threshold']['ls'],
    lw=CURVE_STYLES['threshold']['lw'],
    alpha=CURVE_STYLES['threshold']['alpha'],
    label='50% threshold'
)

if P_crossover is not None:
    ax3.axvline(P_crossover, color='gray', ls=':', lw=1.5, alpha=0.7)
    ax3.text(
        P_crossover * 1.5, 53,
        f'P = {P_crossover:.1e} Pa',
        fontsize=STYLE['fontsize_annotation']-1, color='gray'
    )
else:
    dominant = 'Surface-limited' if n_surface > n_oxide else 'Oxide-limited'
    ax3.text(
        0.05, 0.50,
        f'No crossover in range\n{dominant} throughout',
        transform=ax3.transAxes,
        fontsize=STYLE['fontsize_annotation'],
        va='center',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )

ax3.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Resistance Fraction (%)',          fontsize=STYLE['fontsize_axis'])
ax3.set_title('(C) Rate-Limiting Step Analysis',   fontsize=STYLE['fontsize_title'])
ax3.set_ylim(0, 105)
ax3.legend(fontsize=STYLE['fontsize_legend'], loc='center right')
ax3.grid(True, alpha=STYLE['grid_alpha'])
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (D) LIMIT CHECK — Fast-Kinetics Recovery (L2a)
# ─────────────────────────────────────────────────────────────────────────────
ax4 = axes[1, 1]

# L2a reference flux (no surface kinetics)
J_L2a = alpha * (np.sqrt(P_RANGE_ARR) - np.sqrt(P_down))

# L2a+L6 at high k_diss (fast kinetics limit)
J_high_k = []
for P_up in P_RANGE_ARR:
    try:
        r = solve_steady_state_flux_L2aL6(
            P_up, P_down,
            1e-3, K_eq, D_ox, K_ox, L_ox    # k_diss → large
        )
        J_high_k.append(r['J_ss'])
    except Exception:
        J_high_k.append(np.nan)
J_high_k = np.array(J_high_k)

valid_parity = ~np.isnan(J_high_k) & (J_L2a > 0) & (J_high_k > 0)

ax4.loglog(
    J_L2a[valid_parity], J_high_k[valid_parity],
    'o', color=CURVE_STYLES['L2a_L6']['color'],
    markersize=6, alpha=0.8,
    label='L2a+L6  (k_diss = 1e-3)'
)

J_min = min(J_L2a[valid_parity].min(), J_high_k[valid_parity].min()) * 0.5
J_max = max(J_L2a[valid_parity].max(), J_high_k[valid_parity].max()) * 2.0
ax4.loglog(
    [J_min, J_max], [J_min, J_max],
    color=CURVE_STYLES['parity']['color'],
    ls=CURVE_STYLES['parity']['ls'],
    lw=CURVE_STYLES['parity']['lw'],
    alpha=CURVE_STYLES['parity']['alpha'],
    label='Perfect agreement'
)

rel_errors = (np.abs(J_high_k[valid_parity] - J_L2a[valid_parity])
              / J_L2a[valid_parity] * 100)
max_err   = rel_errors.max()
mean_err  = rel_errors.mean()
box_color = 'lightgreen' if max_err < 1.0 else 'lightyellow'

ax4.text(
    0.05, 0.95,
    f'k_diss = 1e-3\n'
    f'Max error  = {max_err:.2e}%\n'
    f'Mean error = {mean_err:.2e}%\n'
    f'{"✓ VALIDATED" if max_err < 1.0 else "⚠ CHECK k_diss"}',
    transform=ax4.transAxes, fontsize=STYLE['fontsize_annotation'],
    va='top',
    bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.9)
)

ax4.set_xlabel('L2a Flux $J_{L2a}$ (mol/m²/s)',              fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('L2a+L6 Flux $J_{L2a+L6}$ (mol/m²/s)',        fontsize=STYLE['fontsize_axis'])
ax4.set_title(r'(D) Limit Check: $k_{diss} \to \infty$ → L2a', fontsize=STYLE['fontsize_title'])
ax4.set_aspect('equal', adjustable='box')
ax4.legend(fontsize=STYLE['fontsize_legend'], loc='lower right')
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.tick_params(labelsize=STYLE['fontsize_tick'])

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

print("\n" + "=" * 60)
print("L2a+L6 FIGURE 1 — VALIDATION SUMMARY")
print("=" * 60)
print(f"  Net log-log slope:         {slope_net:.4f}")
print(f"  θ range:                   "
      f"{np.nanmin(theta_arr):.4f} → {np.nanmax(theta_arr):.4f}")
print(f"  L2a recovery error:        {max_err:.2e}%  "
      f"({'✓' if max_err < 1.0 else '⚠'})")
if P_crossover is not None:
    print(f"  Crossover pressure:        {P_crossover:.2e} Pa")
else:
    dominant = 'surface-limited' if n_surface > n_oxide else 'oxide-limited'
    print(f"  No crossover — system is {dominant} across full range")
print("=" * 60)
```

```python
# =============================================================================
# CELL 5 — FIGURE 2: Extended Analysis (1×3)
# =============================================================================
"""
Figure 2 examines three additional dimensions: oxide thickness sensitivity,
Arrhenius activation energy per regime, and temperature effects on the
rate-limiting balance.

Panel (E) — Flux vs Oxide Thickness
    Question: How does flux respond to oxide thickness, and does the
    response change between surface-limited and oxide-limited regimes?
    Sweeps L_ox from 0.1 μm to 10 μm at fixed T and P. In the oxide-limited
    regime flux scales as J ∝ 1/L_ox (slope = -1 on log-log). In the
    surface-limited regime the surface kinetics is the bottleneck so oxide
    thickness has no effect — the curve goes flat. Two pressures are used:
    one in each dominant regime (from Cell 3). Both curves are plotted on
    the same axes. A theoretical J ∝ 1/L_ox reference line is overlaid to
    validate the oxide-limited slope.

Panel (F) — Arrhenius: apparent Ea per regime, classified by regime purity
    Question: Does the model extract the correct activation energy in each
    regime, and does it match the known reference values from config?
    Sweeps temperature from T_RANGE[0] to T_RANGE[1]. At each temperature
    all properties are recomputed via helper functions. The solver is called
    at two fixed pressures — one per dominant regime from Cell 3. Each T
    point is classified into one of five categories based on actual regime
    fractions: pure_surface, mixed_surface, pure_oxide, mixed_oxide, mixed.
    Large opaque markers = pure (used for fit). Small faded markers = mixed
    (shown for context, excluded from fit). Dashed line = Arrhenius fit
    through pure points only. Extracted Ea compared to config references.
    Dual x-axis: 1000/T bottom, °C top.

Panel (G) — Rate-Limiting Map: surface fraction vs pressure at multiple T
    Question: How does the balance between surface and oxide resistance shift
    as temperature increases, and in which direction does the crossover move?
    Runs the full pressure sweep at 5 temperatures spread across T_RANGE.
    Plots surface resistance fraction (%) vs pressure for each temperature
    using a plasma colormap — cool for low T, warm for high T. The direction
    of the crossover shift reveals whether surface kinetics or oxide diffusion
    activates faster with temperature, determined by comparing E_diss vs
    E_D_ox + H_sol_ox from config.
"""

fig2, axes2 = plt.subplots(1, 3, figsize=(20, 6))
fig2.suptitle(
    f'L2a+L6: Extended Analysis — Thickness, Arrhenius & Temperature\n'
    f'{OXIDE_KEY}  |  T = {T_op-273.15:.0f}°C',
    fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=1.01
)

# ─────────────────────────────────────────────────────────────────────────────
# (E) FLUX vs OXIDE THICKNESS
# ─────────────────────────────────────────────────────────────────────────────
ax_E = axes2[0]

L_ox_sweep = np.logspace(
    np.log10(CONDITIONS['L_oxide_range'][0]),
    np.log10(CONDITIONS['L_oxide_range'][1]),
    CONDITIONS['n_L_points']
)

pressure_cases = []
if P_for_surface_Ea is not None:
    pressure_cases.append((
        P_for_surface_Ea,
        CURVE_STYLES['surface_region']['color'],
        f'P = {P_for_surface_Ea:.1e} Pa  (surface-dominated)'
    ))
if P_for_oxide_Ea is not None:
    pressure_cases.append((
        P_for_oxide_Ea,
        CURVE_STYLES['oxide_region']['color'],
        f'P = {P_for_oxide_Ea:.1e} Pa  (oxide-dominated)'
    ))
if not pressure_cases:
    pressure_cases = [
        (1e2, CURVE_STYLES['surface_region']['color'], 'P = 1e2 Pa'),
        (1e8, CURVE_STYLES['oxide_region']['color'],   'P = 1e8 Pa'),
    ]

for P_fixed, color, lbl in pressure_cases:
    J_thickness = []
    for L_test in L_ox_sweep:
        try:
            r = solve_steady_state_flux_L2aL6(
                P_fixed, P_down,
                k_diss, K_eq, D_ox, K_ox, L_test
            )
            J_thickness.append(r['J_ss'])
        except Exception:
            J_thickness.append(np.nan)

    J_thickness = np.array(J_thickness)
    valid_th    = ~np.isnan(J_thickness) & (J_thickness > 0)
    if not np.any(valid_th):
        continue

    ax_E.loglog(
        L_ox_sweep[valid_th] * 1e6, J_thickness[valid_th],
        color=color, lw=STYLE['linewidth'],
        marker='o', ms=4, markevery=3,
        label=lbl
    )

    slope_th, *_ = stats.linregress(
        np.log10(L_ox_sweep[valid_th]),
        np.log10(np.abs(J_thickness[valid_th]))
    )
    mid_idx = np.sum(valid_th) // 2
    ax_E.text(
        L_ox_sweep[valid_th][mid_idx] * 1e6,
        J_thickness[valid_th][mid_idx],
        f'Slope={slope_th:.2f}',
        color=color,
        fontsize=STYLE['fontsize_annotation'],
        fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )

if P_for_oxide_Ea is not None:
    J_theory_ref = solve_steady_state_flux_L2aL6(
        P_for_oxide_Ea, P_down,
        k_diss, K_eq, D_ox, K_ox, L_ox_sweep[len(L_ox_sweep)//2]
    )['J_ss']
    L_ref    = L_ox_sweep[len(L_ox_sweep)//2]
    J_theory = J_theory_ref * (L_ref / L_ox_sweep)
    ax_E.loglog(
        L_ox_sweep * 1e6, J_theory,
        color='gray', ls='--', lw=1.5, alpha=0.6,
        label='Theory: $J \\propto 1/L_{ox}$'
    )

ax_E.text(
    0.05, 0.05,
    'Oxide-limited: slope = -1\nSurface-limited: slope = 0 (flat)',
    transform=ax_E.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='bottom',
    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
)
ax_E.set_xlabel('Oxide Thickness $L_{ox}$ (μm)',  fontsize=STYLE['fontsize_axis'])
ax_E.set_ylabel('Flux $J_{ss}$ (mol/m²/s)',        fontsize=STYLE['fontsize_axis'])
ax_E.set_title('(E) Flux vs Oxide Thickness',      fontsize=STYLE['fontsize_title'])
ax_E.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper right')
ax_E.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax_E.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (F) ARRHENIUS — apparent Ea per regime, classified by regime purity
# ─────────────────────────────────────────────────────────────────────────────
ax_F = axes2[1]

temperatures_K = np.linspace(T_RANGE[0], 800, N_T)
inv_T          = 1000.0 / temperatures_K

def classify_regime_L2a(fs, fo, pure_thresh=0.90, mixed_thresh=0.50):
    if fs >= pure_thresh:    return 'pure_surface'
    elif fo >= pure_thresh:  return 'pure_oxide'
    elif fs >= mixed_thresh: return 'mixed_surface'
    elif fo >= mixed_thresh: return 'mixed_oxide'
    else:                    return 'mixed'

regime_class_styles = {
    'pure_surface':  {'color': CURVE_STYLES['surface_region']['color'], 'marker': 'o', 'ms': 8,  'alpha': 1.0},
    'mixed_surface': {'color': CURVE_STYLES['surface_region']['color'], 'marker': 'o', 'ms': 5,  'alpha': 0.35},
    'pure_oxide':    {'color': CURVE_STYLES['oxide_region']['color'],   'marker': 's', 'ms': 8,  'alpha': 1.0},
    'mixed_oxide':   {'color': CURVE_STYLES['oxide_region']['color'],   'marker': 's', 'ms': 5,  'alpha': 0.35},
    'mixed':         {'color': CURVE_STYLES['mixed_region']['color'],   'marker': 'x', 'ms': 6,  'alpha': 0.5},
}

regime_Ea_cases    = [
    ('surface', P_for_surface_Ea, CURVE_STYLES['surface_region']['color']),
    ('oxide',   P_for_oxide_Ea,   CURVE_STYLES['oxide_region']['color']),
]
results            = {}
plotted_cls_labels = set()

for label, P_fixed, color in regime_Ea_cases:
    if P_fixed is None:
        print(f"  ⚠  {label}: P_for_{label}_Ea is None — fit skipped")
        results[label] = None
        continue

    J_row     = []
    class_row = []

    for T_K in temperatures_K:
        try:
            r = solve_steady_state_flux_L2aL6(
                P_fixed, P_down,
                get_k_diss_at_T(T_K),
                get_K_eq_at_T(T_K),
                get_D_ox_at_T(T_K),
                get_K_ox_at_T(T_K),
                L_ox,
            )
            J_row.append(r['J_ss'])
            fs = r['resistances']['fraction_surface']
            fo = r['resistances']['fraction_oxide']
            class_row.append(classify_regime_L2a(fs, fo))
        except Exception:
            J_row.append(np.nan)
            class_row.append('mixed')

    J_arr_T   = np.array(J_row)
    class_arr = np.array(class_row)
    pure      = (class_arr == f'pure_{label}') & ~np.isnan(J_arr_T) & (J_arr_T > 0)

    # Plot all classified points onto ax_F
    for cls, sty in regime_class_styles.items():
        mask = (class_arr == cls) & ~np.isnan(J_arr_T) & (J_arr_T > 0)
        if not np.any(mask):
            continue
        cls_label = cls.replace('_', ' ').capitalize()
        ax_F.semilogy(
            inv_T[mask], J_arr_T[mask],
            ls='none',
            color=sty['color'], marker=sty['marker'],
            ms=sty['ms'], alpha=sty['alpha'],
            label=cls_label if cls_label not in plotted_cls_labels
                  else '_nolegend_'
        )
        plotted_cls_labels.add(cls_label)

    if np.sum(pure) < 3:
        print(f"  ⚠  {label}: only {np.sum(pure)} pure T points — fit skipped")
        results[label] = None
        continue

    slope, intercept, r_val, *_ = stats.linregress(
        inv_T[pure], np.log(J_arr_T[pure])
    )
    E_a = -slope * R * 1000

    results[label] = {
        'J': J_arr_T, 'valid': pure,
        'slope': slope, 'intercept': intercept,
        'r_sq': r_val**2, 'E_a': E_a,
        'P_fixed': P_fixed, 'color': color,
        'n_pure': int(np.sum(pure)),
        'T_pure_range': (temperatures_K[pure].min() - 273.15,
                         temperatures_K[pure].max() - 273.15),
    }

    # Fitted line onto ax_F
    ax_F.semilogy(
        inv_T[pure],
        np.exp(slope * inv_T[pure] + intercept),
        color=color, ls='--', lw=2.0, alpha=0.9,
        label=f'Fit {label}  (pure pts)'
    )

# Dual x-axis on ax_F
ax_F_top = ax_F.twiny()
ax_F_top.set_xlim(ax_F.get_xlim())
T_ticks_C = np.array([200, 300, 400, 500, 600, 700, 800, 900])
T_ticks_K = T_ticks_C + 273.15
in_range  = (T_ticks_K >= T_RANGE[0]) & (T_ticks_K <= T_RANGE[1])
if np.any(in_range):
    ax_F_top.set_xticks(1000.0 / T_ticks_K[in_range])
    ax_F_top.set_xticklabels([f'{t}' for t in T_ticks_C[in_range]])
ax_F_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax_F_top.tick_params(labelsize=STYLE['fontsize_tick'])

# Annotation box on ax_F
annot = []
for lbl, E_ref, ref_name in [
    ('surface', E_REF_SURFACE, 'E_diss'),
    ('oxide',   E_REF_OXIDE,   'E_D_ox + H_sol'),
]:
    res = results.get(lbl)
    if res is None:
        annot.append(f'{lbl.capitalize()}: not extractable')
        annot.append(f'  Expected ({ref_name}): {E_ref/1000:.1f} kJ/mol')
    else:
        match = abs(res['E_a'] - E_ref) < 3000
        annot.append(f'{lbl.capitalize()}-limited:')
        annot.append(f'  Extracted: {res["E_a"]/1000:.1f} kJ/mol')
        annot.append(f'  Expected ({ref_name}): {E_ref/1000:.1f} kJ/mol')
        annot.append(f'  R²={res["r_sq"]:.4f}  {"✓" if match else "⚠"}')
        annot.append(f'  T range: {res["T_pure_range"][0]:.0f}–'
                     f'{res["T_pure_range"][1]:.0f}°C  '
                     f'({res["n_pure"]} pure pts)')
    annot.append('')

ax_F.text(
    0.97, 0.97, '\n'.join(annot).rstrip(),
    transform=ax_F.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='top', ha='right',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9)
)
ax_F.set_xlabel('1000/T (K⁻¹)',             fontsize=STYLE['fontsize_axis'])
ax_F.set_ylabel('Flux $J_{ss}$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax_F.set_title('(F) Arrhenius: Apparent $E_a$ per Regime',
               fontsize=STYLE['fontsize_title'])
ax_F.legend(fontsize=STYLE['fontsize_legend']-2, loc='lower left',
            ncol=2, title='Symbol size = purity')
ax_F.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax_F.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (G) RATE-LIMITING MAP
# ─────────────────────────────────────────────────────────────────────────────
ax_G = axes2[2]

T_map_values        = np.linspace(T_RANGE[0], T_RANGE[1], 5)
T_colors            = plt.cm.plasma(np.linspace(0.1, 0.9, len(T_map_values)))
crossover_pressures = []

for T_K, col in zip(T_map_values, T_colors):
    fs_row = []
    for P_up in P_RANGE_ARR:
        try:
            r = solve_steady_state_flux_L2aL6(
                P_up, P_down,
                get_k_diss_at_T(T_K),
                get_K_eq_at_T(T_K),
                get_D_ox_at_T(T_K),
                get_K_ox_at_T(T_K),
                L_ox,
            )
            fs_row.append(r['resistances']['fraction_surface'] * 100)
        except Exception:
            fs_row.append(np.nan)

    fs_row    = np.array(fs_row)
    valid_row = ~np.isnan(fs_row)
    if not np.any(valid_row):
        continue

    ax_G.semilogx(
        P_RANGE_ARR[valid_row], fs_row[valid_row],
        color=col, lw=STYLE['linewidth'],
        label=f'T = {T_K-273.15:.0f}°C'
    )

    sc = np.where(np.diff(np.sign(fs_row[valid_row] - 50)))[0]
    if len(sc) > 0:
        idx  = sc[0]
        fs_v = fs_row[valid_row]
        if idx + 1 < len(fs_v) and fs_v[idx] > 50 and fs_v[idx+1] < 50:
            crossover_pressures.append((T_K, P_RANGE_ARR[valid_row][idx]))

ax_G.axhline(
    50,
    color=CURVE_STYLES['threshold']['color'],
    ls=CURVE_STYLES['threshold']['ls'],
    lw=CURVE_STYLES['threshold']['lw'],
    alpha=CURVE_STYLES['threshold']['alpha'],
    label='50% threshold'
)

if len(crossover_pressures) >= 2:
    T_low,  P_co_low  = crossover_pressures[0]
    T_high, P_co_high = crossover_pressures[-1]
    if P_co_high < P_co_low:
        direction = 'lower P as T ↑'
        physical  = (f'E_diss ({E_REF_SURFACE/1000:.0f}) > '
                     f'E_D_ox+H_sol ({E_REF_OXIDE/1000:.0f}) kJ/mol')
        meaning   = 'Surface activates faster → less limiting at high T'
    else:
        direction = 'higher P as T ↑'
        physical  = (f'E_D_ox+H_sol ({E_REF_OXIDE/1000:.0f}) > '
                     f'E_diss ({E_REF_SURFACE/1000:.0f}) kJ/mol')
        meaning   = 'Oxide activates faster → more limiting at high T'
    ax_G.text(
        0.05, 0.05,
        f'Crossover shifts to {direction}\n{physical}\n{meaning}',
        transform=ax_G.transAxes,
        fontsize=STYLE['fontsize_annotation']-1, va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )
elif len(crossover_pressures) == 0:
    dominant = 'surface-limited' if n_surface > n_oxide else 'oxide-limited'
    ax_G.text(
        0.05, 0.05,
        f'No crossover at any T\nSystem is {dominant}\nacross full pressure range',
        transform=ax_G.transAxes,
        fontsize=STYLE['fontsize_annotation']-1, va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )

ax_G.set_xlabel('Upstream Pressure $P_{up}$ (Pa)',  fontsize=STYLE['fontsize_axis'])
ax_G.set_ylabel('Surface Resistance Fraction (%)',  fontsize=STYLE['fontsize_axis'])
ax_G.set_title('(G) Rate-Limiting Map vs Temperature',
               fontsize=STYLE['fontsize_title'])
ax_G.set_ylim(0, 105)
ax_G.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper right',
            title='Temperature')
ax_G.grid(True, alpha=STYLE['grid_alpha'])
ax_G.tick_params(labelsize=STYLE['fontsize_tick'])

plt.tight_layout()
plt.show()

print("\n" + "=" * 60)
print("L2a+L6 FIGURE 2 — SUMMARY")
print("=" * 60)
for lbl in ('surface', 'oxide'):
    res = results.get(lbl)
    E_ref = E_REF_SURFACE if lbl == 'surface' else E_REF_OXIDE
    if res is None:
        print(f"  {lbl.capitalize()}: not extractable")
    else:
        match = abs(res['E_a'] - E_ref) < 3000
        print(f"  {lbl.capitalize()}: extracted={res['E_a']/1000:.1f}  "
              f"expected={E_ref/1000:.1f} kJ/mol  "
              f"({'✓' if match else '⚠'})  "
              f"{res['n_pure']} pure pts")
if crossover_pressures:
    print(f"\n  Crossover pressures across T sweep:")
    for T_K, P_co in crossover_pressures:
        print(f"    T = {T_K-273.15:.0f}°C  →  P_crossover = {P_co:.2e} Pa")
else:
    dominant = 'surface-limited' if n_surface > n_oxide else 'oxide-limited'
    print(f"  No crossover found — system predominantly {dominant}")
print("=" * 60)
```

```python
# =============================================================================
# CELL 6 — DATAFRAME OUTPUT
# =============================================================================
df = pd.DataFrame(rows)

df_display = df.copy()
for col in df_display.columns:
    if col in ('Rate-Limiting', 'Error'):
        continue
    if df_display[col].dtype == float:
        if 'fraction' in col.lower():
            df_display[col] = df_display[col].round(2)
        elif 'theta' in col.lower():
            df_display[col] = df_display[col].round(6)
        else:
            df_display[col] = df_display[col].apply(
                lambda x: f'{x:.4e}' if pd.notna(x) else 'NaN'
            )

n_surface_df = (df['Rate-Limiting'] == 'SURFACE').sum()
n_oxide_df   = (df['Rate-Limiting'] == 'OXIDE').sum()
n_mixed_df   = (df['Rate-Limiting'] == 'MIXED').sum()
n_error_df   = (df['Rate-Limiting'] == 'ERROR').sum()

print("=" * 60)
print("L2a+L6 RESULTS DATAFRAME")
print("=" * 60)
print(f"  Oxide:       {OXIDE_KEY}")
print(f"  T:           {T_op-273.15:.0f}°C  ({T_op} K)")
print(f"  L_ox:        {L_ox*1e6:.1f} μm")
print(f"  k_diss:      {k_diss:.3e}  mol/m²/s/Pa")
print(f"  K_eq:        {K_eq:.3e}  Pa⁻¹")
print(f"  D_ox:        {D_ox:.3e}  m²/s")
print(f"  K_ox:        {K_ox:.3e}  mol/m³/Pa^0.5")
print(f"\n  Total rows:        {len(df)}")
print(f"  Surface-limited:   {n_surface_df}")
print(f"  Oxide-limited:     {n_oxide_df}")
print(f"  Mixed:             {n_mixed_df}")
if n_error_df > 0:
    print(f"  Errors:            {n_error_df}  ⚠")
print("=" * 60)

display(df_display)
```

### L2a +L6

```python
# =============================================================================
# CELL 1 — SETUP
# Run once. All imports and helper functions.
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import brentq
import pandas as pd
from itertools import groupby
from operator import itemgetter
from IPython.display import display

from calculations.surface_kinetics import solve_steady_state_flux_direct
from calculations.config.model_config import (
    R, F,
    METALS, OXIDES, CONDITIONS,
    PLOT_STYLE as STYLE,
    CURVE_STYLES,
    apply_style,
    build_simulation_config,
)

# ── Simulation config ─────────────────────────────────────────────────────────
SIM       = build_simulation_config(
    metal='metal_316L_Heat_treated_ref_cast',
    oxide='Al2O3',
    T_operating=873,
    P_upstream=1e5,
    L_metal=1e-3,
)

METAL_KEY = SIM['metal_name']
OXIDE_KEY = SIM['oxide_name']
METAL     = METALS[METAL_KEY]
OXIDE     = OXIDES[OXIDE_KEY]
SK        = OXIDE['surface_kinetics']    # oxide surface kinetics sub-dict

# ── Pressure sweep ────────────────────────────────────────────────────────────
P_RANGE_ARR = np.logspace(0, 12, 100)
P_DOWN      = SIM['P_downstream']
T_RANGE     = CONDITIONS['T_range']
N_T         = CONDITIONS['n_T_points']
```

```python
# =============================================================================
# CELL 2 — PARAMETERS
# All physics values at operating temperature, pulled from config.
# Change placeholder values in config and rerun from this cell downward.
# =============================================================================

# ── Operating conditions ──────────────────────────────────────────────────────
T_op   = SIM['T_operating']
L_ox   = SIM['L_oxide']
L_m    = SIM['L_metal']
P_down = SIM['P_downstream']

# ── Surface kinetics at T_op ──────────────────────────────────────────────────
k_diss   = get_k_diss_at_T(T_op)
K_eq     = get_K_eq_at_T(T_op)
k_recomb = k_diss / K_eq

# ── Oxide transport at T_op ───────────────────────────────────────────────────
D_ox = get_D_ox_at_T(T_op)
K_ox = get_K_ox_at_T(T_op)

# ── Metal transport at T_op ───────────────────────────────────────────────────
D_m   = get_D_m_at_T(T_op)
K_s_m = get_K_s_m_at_T(T_op)

# ── Derived ───────────────────────────────────────────────────────────────────
alpha = D_ox * K_ox / L_ox
beta  = D_m  * K_s_m / L_m

# ── Print verification ────────────────────────────────────────────────────────
print("=" * 60)
print(f"L2+L6 PARAMETERS AT T = {T_op-273.15:.0f}°C  ({T_op} K)")
print("=" * 60)

print(f"\n  Geometry")
print(f"    L_ox     = {L_ox:.3e}  m")
print(f"    L_m      = {L_m:.3e}  m")
print(f"    P_down   = {P_down:.1e}  Pa")

print(f"\n  Surface kinetics  [config: oxide surface_kinetics sub-dict]")
print(f"    k_diss   = {k_diss:.3e}  mol/m²/s/Pa")
print(f"    K_eq     = {K_eq:.3e}  Pa⁻¹")
print(f"    k_recomb = {k_recomb:.3e}  mol/m²/s  (= k_diss / K_eq)")

print(f"\n  Oxide transport  [config: Arrhenius]")
print(f"    D_ox     = {D_ox:.3e}  m²/s")
print(f"    K_ox     = {K_ox:.3e}  mol/m³/Pa^0.5")

print(f"\n  Metal transport  [config: Arrhenius]")
print(f"    D_m      = {D_m:.3e}  m²/s")
print(f"    K_s_m    = {K_s_m:.3e}  mol/m³/Pa^0.5")

print(f"\n  Derived")
print(f"    α        = {alpha:.3e}  mol/m²/s/Pa^0.5  (= D_ox × K_ox / L_ox)")
print(f"    β        = {beta:.3e}  mol/m²/s/Pa^0.5  (= D_m × K_s_m / L_m)")
print(f"    α/β      = {alpha/beta:.3e}  (oxide/metal permeance ratio)")

print(f"\n  Reference activation energies  [config]")
print(f"    E_diss        = {E_REF_SURFACE/1000:.1f} kJ/mol   (surface)")
print(f"    E_D_ox+H_sol  = {E_REF_OXIDE/1000:.1f} kJ/mol   (oxide permeation)")
print(f"    E_D+H_s       = {E_REF_METAL/1000:.1f} kJ/mol   (metal permeation)")

if SK['k_diss_ref'] <= 1e-12:
    print(f"\n  ⚠  WARNING: k_diss_ref = {SK['k_diss_ref']:.0e} looks like a placeholder.")
    print(f"     Update config oxide surface_kinetics sub-dict with literature values.")

print("\n" + "=" * 60)
```

```python
# =============================================================================
# CELL 3 — COMPUTE
# Single-pass loop over full pressure range.
# Builds all arrays, regime summary, and Arrhenius pressure selection.
# =============================================================================

# ── Single-pass loop ──────────────────────────────────────────────────────────
rows             = []
J_arr            = []
theta_arr        = []
P_int_arr        = []
frac_surface_arr = []
frac_oxide_arr   = []
frac_metal_arr   = []
R_surface_arr    = []
R_oxide_arr      = []
R_metal_arr      = []
R_total_arr      = []
rate_lim_arr     = []

for P_up in P_RANGE_ARR:
    try:
        r = solve_steady_state_flux_direct(
            P_up, P_down, L_m,
            k_diss, K_eq,
            D_ox, K_ox, L_ox,
            D_m, K_s_m
        )
        J_arr.append(r['J_ss'])
        theta_arr.append(r['theta'])
        P_int_arr.append(r['P_int'])
        frac_surface_arr.append(r['resistances']['fraction_surface'])
        frac_oxide_arr.append(r['resistances']['fraction_oxide'])
        frac_metal_arr.append(r['resistances']['fraction_metal'])
        R_surface_arr.append(r['resistances']['R_surface'])
        R_oxide_arr.append(r['resistances']['R_oxide'])
        R_metal_arr.append(r['resistances']['R_metal'])
        R_total_arr.append(r['resistances']['R_total'])
        rate_lim_arr.append(r['rate_limiting'])

        rows.append({
            'P_up (Pa)':            P_up,
            'P_int (Pa)':           r['P_int'],
            'J_ss (mol/m²/s)':      r['J_ss'],
            'θ_surface':            r['theta'],
            'α_oxide':              r['alpha'],
            'β_metal':              r['beta'],
            'R_surface':            r['resistances']['R_surface'],
            'R_oxide':              r['resistances']['R_oxide'],
            'R_metal':              r['resistances']['R_metal'],
            'R_total':              r['resistances']['R_total'],
            'fraction_surface (%)': r['resistances']['fraction_surface'] * 100,
            'fraction_oxide (%)':   r['resistances']['fraction_oxide']   * 100,
            'fraction_metal (%)':   r['resistances']['fraction_metal']   * 100,
            'Rate-Limiting':        r['rate_limiting'].upper(),
        })

    except Exception as e:
        for lst in [J_arr, theta_arr, P_int_arr, frac_surface_arr,
                    frac_oxide_arr, frac_metal_arr, R_surface_arr,
                    R_oxide_arr, R_metal_arr, R_total_arr]:
            lst.append(np.nan)
        rate_lim_arr.append('error')
        rows.append({'P_up (Pa)': P_up, 'Rate-Limiting': 'ERROR', 'Error': str(e)})

# ── Convert to arrays ─────────────────────────────────────────────────────────
J_arr            = np.array(J_arr)
theta_arr        = np.array(theta_arr)
P_int_arr        = np.array(P_int_arr)
frac_surface_arr = np.array(frac_surface_arr)
frac_oxide_arr   = np.array(frac_oxide_arr)
frac_metal_arr   = np.array(frac_metal_arr)
R_surface_arr    = np.array(R_surface_arr)
R_oxide_arr      = np.array(R_oxide_arr)
R_metal_arr      = np.array(R_metal_arr)
R_total_arr      = np.array(R_total_arr)
rate_lim_arr     = np.array(rate_lim_arr)
valid            = ~np.isnan(J_arr)

# ── Rate-limiting classification array ───────────────────────────────────────
rl_arr = np.where(frac_surface_arr > 0.5, 'surface',
         np.where(frac_oxide_arr   > 0.5, 'oxide',
         np.where(frac_metal_arr   > 0.5, 'metal', 'mixed')))

# ── Regime summary ────────────────────────────────────────────────────────────
n_surface = np.sum((frac_surface_arr > 0.5) & valid)
n_oxide   = np.sum((frac_oxide_arr   > 0.5) & valid)
n_metal   = np.sum((frac_metal_arr   > 0.5) & valid)
n_mixed   = np.sum(
    (frac_surface_arr <= 0.5) & (frac_oxide_arr <= 0.5) &
    (frac_metal_arr   <= 0.5) & valid
)

surf_dom_mask  = (frac_surface_arr > 0.90) & valid
oxide_dom_mask = (frac_oxide_arr   > 0.90) & valid
metal_dom_mask = (frac_metal_arr   > 0.90) & valid
P_surf_dom     = P_RANGE_ARR[surf_dom_mask]
P_oxide_dom    = P_RANGE_ARR[oxide_dom_mask]
P_metal_dom    = P_RANGE_ARR[metal_dom_mask]

# ── Arrhenius pressure selection (no crossover dependence) ───────────────────
P_for_surface_Ea = np.median(P_surf_dom)  if len(P_surf_dom)  > 3 else None
P_for_oxide_Ea   = np.median(P_oxide_dom) if len(P_oxide_dom) > 3 else None
P_for_metal_Ea   = np.median(P_metal_dom) if len(P_metal_dom) > 3 else None

# ── Crossover pressures (diagnostic only) ────────────────────────────────────
def find_crossover(frac_a, frac_b, valid_mask):
    """Find pressure where frac_a crosses frac_b (from above)."""
    vc = valid_mask & ~np.isnan(frac_a) & ~np.isnan(frac_b)
    if np.sum(vc) < 2:
        return None
    diff = frac_a[vc] - frac_b[vc]
    Pv   = P_RANGE_ARR[vc]
    sc   = np.where(np.diff(np.sign(diff)))[0]
    if len(sc) == 0:
        return None
    idx = sc[0]
    if diff[idx] > 0 and diff[idx+1] < 0:
        return Pv[idx]
    return None

P_cross_surf_oxide = find_crossover(frac_surface_arr, frac_oxide_arr, valid)
P_cross_surf_metal = find_crossover(frac_surface_arr, frac_metal_arr, valid)
P_cross_oxide_metal= find_crossover(frac_oxide_arr,   frac_metal_arr, valid)

# ── Print regime summary ──────────────────────────────────────────────────────
print("=" * 60)
print(f"L2+L6 REGIME SUMMARY  (T = {T_op-273.15:.0f}°C, {np.sum(valid)} valid points)")
print("=" * 60)

print(f"\n  Regime breakdown (dominant fraction > 50%)")
print(f"    Surface-limited:  {n_surface:3d} points")
print(f"    Oxide-limited:    {n_oxide:3d} points")
print(f"    Metal-limited:    {n_metal:3d} points")
print(f"    Mixed:            {n_mixed:3d} points")

print(f"\n  Dominant regions (fraction > 90%) — used for Arrhenius")
for label, dom, P_ea in [
    ('Surface', P_surf_dom,  P_for_surface_Ea),
    ('Oxide',   P_oxide_dom, P_for_oxide_Ea),
    ('Metal',   P_metal_dom, P_for_metal_Ea),
]:
    if len(dom) > 3:
        print(f"    {label}-dominant: {len(dom):3d} points  "
              f"P = {dom.min():.1e} – {dom.max():.1e} Pa  "
              f"→ P_for_{label.lower()}_Ea = {P_ea:.2e} Pa")
    else:
        print(f"    {label}-dominant: {len(dom):3d} points  "
              f"(insufficient — {label.lower()} Ea not extractable)")

print(f"\n  Crossover pressures (diagnostic only)")
for desc, P_co in [
    ('Surface → Oxide',  P_cross_surf_oxide),
    ('Surface → Metal',  P_cross_surf_metal),
    ('Oxide   → Metal',  P_cross_oxide_metal),
]:
    val = f'{P_co:.2e} Pa' if P_co is not None else 'not found'
    print(f"    {desc}:  {val}")

print(f"\n  α/β ratio:  {alpha/beta:.3e}")
print(f"  θ range:    {np.nanmin(theta_arr):.4f} → {np.nanmax(theta_arr):.4f}")
print(f"  J range:    {np.nanmin(J_arr):.2e} → {np.nanmax(J_arr):.2e} mol/m²/s")
print("=" * 60)
```

```python
# =============================================================================
# CELL 4 — FIGURE 1: Core Validation (2×2)
# =============================================================================
"""
Figure 1 validates the L2+L6 model from four angles at operating temperature.
This is the first level with three competing resistances, making it the
canonical figure for the full series-resistance framework.

Panel (A) — Flux vs Pressure
    Question: Does the model produce the correct three-regime pressure
    behaviour?
    The hero plot. Shows J_ss vs P_up on a log-log scale. The curve should
    show three distinct regimes: slope≈1 (surface-limited, low P), slope≈0.5
    (oxide-limited or metal-limited, intermediate/high P), with possible
    transitions between them. Style follows the widget — no fill_between.
    The curve is recolored by regime using thick colored segments with slope
    annotation boxes sitting directly on each segment. The black model curve
    sits underneath, the colored overlay sits on top. Reference lines for
    slope=1 and slope=0.5 are drawn as dashed lines. Net slope annotation
    box in the bottom right corner.

Panel (B) — Arrhenius: apparent Ea per regime
    Question: Does the model extract the correct activation energy in each
    regime, and does it match the known reference values from config?
    Sweeps temperature from T_RANGE[0] to T_RANGE[1]. At each temperature
    all properties are recomputed via helper functions. The solver is called
    at up to three fixed pressures — one per dominant regime — selected in
    Cell 3 as the median of each dominant region (fraction > 90%). A linear
    fit to ln(J) vs 1000/T gives the apparent Ea for each regime. Extracted
    Ea is compared to config references: E_diss (surface), E_D_ox+H_sol_ox
    (oxide), E_D+H_s (metal). Dual x-axis: 1000/T bottom, °C top.

Panel (C) — Resistance Fractions vs Pressure
    Question: Which step controls the rate, and how does the balance shift
    with pressure?
    Plots R_surface/R_total, R_oxide/R_total, and R_metal/R_total as
    percentages vs pressure on a semilog scale — consistent with L1+L6
    and L2a+L6. With three resistances there can be two crossover
    pressures: surface→oxide (S→O) and oxide→metal (O→M). Both are
    marked with vertical lines if present. The 50% threshold line marks
    where each fraction is dominant. If no crossovers exist a text
    annotation states which regime dominates throughout.

Panel (D) — Limit Check: k_diss → ∞ recovers non-L6 result
    Question: As k_diss → ∞, does L2+L6 recover the oxide+metal result
    without surface kinetics?
    Runs the solver at k_diss = 1e-3 (fast kinetics limit) and compares
    J_L2+L6 against the analytical series-resistance flux:
    J = (sqrt(P_up) - sqrt(P_down)) / (1/alpha + 1/beta).
    A parity plot shows agreement — perfect match is a diagonal line.
    This validates the model correctly recovers the L2+L5 (no L6) limit
    when surface kinetics is fast. Green box if max error < 1%.
"""

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle(
    f'L2+L6: Perfect Oxide + Perfect Metal + Surface Kinetics — Core Validation\n'
    f'{OXIDE_KEY} + {METAL_KEY}  |  T = {T_op-273.15:.0f}°C  |  '
    f'L_ox = {L_ox*1e6:.1f} μm  |  L_m = {L_m*1e3:.1f} mm',
    fontsize=STYLE['fontsize_suptitle']-1, fontweight='bold', y=0.98
)

# ─────────────────────────────────────────────────────────────────────────────
# (A) FLUX vs PRESSURE — three regimes
# ─────────────────────────────────────────────────────────────────────────────
ax1 = axes[0, 0]

ax1.loglog(
    P_RANGE_ARR[valid], J_arr[valid],
    color='black', lw=STYLE['linewidth'], label='L2+L6 Model'
)

if np.any(valid):
    P_ref1  = P_RANGE_ARR[valid][0]
    J_ref1  = J_arr[valid][0]
    ax1.loglog(
        P_RANGE_ARR[valid],
        J_ref1 * (P_RANGE_ARR[valid] / P_ref1) ** 1.0,
        color=CURVE_STYLES['slope_1']['color'],
        ls=CURVE_STYLES['slope_1']['ls'],
        lw=CURVE_STYLES['slope_1']['lw'],
        alpha=CURVE_STYLES['slope_1']['alpha'],
        label='Slope = 1 (surface)'
    )
    P_ref05 = P_RANGE_ARR[valid][-1]
    J_ref05 = J_arr[valid][-1]
    ax1.loglog(
        P_RANGE_ARR[valid],
        J_ref05 * (P_RANGE_ARR[valid] / P_ref05) ** 0.5,
        color=CURVE_STYLES['slope_05']['color'],
        ls=CURVE_STYLES['slope_05']['ls'],
        lw=CURVE_STYLES['slope_05']['lw'],
        alpha=CURVE_STYLES['slope_05']['alpha'],
        label='Slope = 0.5 (diffusion)'
    )

# Colored regime segments + slope annotation boxes
region_styles = [
    ('surface', CURVE_STYLES['surface_region']['color'], 'Surface-limited'),
    ('oxide',   CURVE_STYLES['oxide_region']['color'],   'Oxide-limited'),
    ('metal',   CURVE_STYLES['metal_region']['color'],   'Metal-limited'),
    ('mixed',   CURVE_STYLES['mixed_region']['color'],   'Mixed'),
]

for regime, color, label in region_styles:
    mask = (rl_arr == regime) & valid
    if not np.any(mask):
        continue
    idxs = np.where(mask)[0]
    for _, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
        group = list(map(itemgetter(1), g))
        if len(group) < 3:
            continue
        P_seg = P_RANGE_ARR[group]
        J_seg = J_arr[group]
        ax1.loglog(P_seg, J_seg, color=color, lw=4, alpha=0.7)
        slope_seg, *_ = stats.linregress(
            np.log10(P_seg), np.log10(np.abs(J_seg))
        )
        mid = len(group) // 2
        ax1.text(
            P_seg[mid], J_seg[mid],
            f'{label}\nSlope={slope_seg:.2f}',
            color=color,
            fontsize=STYLE['fontsize_annotation'],
            fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

log_P_v   = np.log10(P_RANGE_ARR[valid])
log_J_v   = np.log10(np.abs(J_arr[valid]))
slope_net, *_ = stats.linregress(log_P_v, log_J_v)

ax1.text(
    0.98, 0.02, f'Net slope = {slope_net:.2f}',
    transform=ax1.transAxes, ha='right', va='bottom',
    fontsize=STYLE['fontsize_annotation'], fontweight='bold',
    bbox=dict(boxstyle='square', fc='wheat', ec='gray', alpha=1)
)

ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux $J_{ss}$ (mol/m²/s)',         fontsize=STYLE['fontsize_axis'])
ax1.set_title('(A) Flux vs Pressure',               fontsize=STYLE['fontsize_title'])
ax1.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper left')
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (B) ARRHENIUS — apparent Ea per regime, classified by regime purity
# ─────────────────────────────────────────────────────────────────────────────
ax2 = axes[0, 1]

temperatures_K = np.linspace(T_RANGE[0], T_RANGE[1], N_T)
inv_T          = 1000.0 / temperatures_K
arrhenius_results = {}

# ── Regime classification helper ─────────────────────────────────────────────
def classify_regime(fs, fo, fm, pure_thresh=0.90, mixed_thresh=0.50):
    """
    Classify a single T point into one of seven regime categories
    based on resistance fractions at that temperature and pressure.
    pure_thresh  : fraction above which a step is called 'pure'
    mixed_thresh : fraction above which a step is called 'dominant'
    """
    if fs >= pure_thresh:    return 'pure_surface'
    elif fo >= pure_thresh:  return 'pure_oxide'
    elif fm >= pure_thresh:  return 'pure_metal'
    elif fs >= mixed_thresh: return 'mixed_surface'
    elif fo >= mixed_thresh: return 'mixed_oxide'
    elif fm >= mixed_thresh: return 'mixed_metal'
    else:                    return 'mixed'

# ── Per-class plot styles ─────────────────────────────────────────────────────
regime_class_styles = {
    'pure_surface':  {'color': CURVE_STYLES['surface_region']['color'],
                      'marker': 'o', 'ms': 8,  'alpha': 1.0},
    'mixed_surface': {'color': CURVE_STYLES['surface_region']['color'],
                      'marker': 'o', 'ms': 5,  'alpha': 0.35},
    'pure_oxide':    {'color': CURVE_STYLES['oxide_region']['color'],
                      'marker': 's', 'ms': 8,  'alpha': 1.0},
    'mixed_oxide':   {'color': CURVE_STYLES['oxide_region']['color'],
                      'marker': 's', 'ms': 5,  'alpha': 0.35},
    'pure_metal':    {'color': CURVE_STYLES['metal_region']['color'],
                      'marker': '^', 'ms': 8,  'alpha': 1.0},
    'mixed_metal':   {'color': CURVE_STYLES['metal_region']['color'],
                      'marker': '^', 'ms': 5,  'alpha': 0.35},
    'mixed':         {'color': CURVE_STYLES['mixed_region']['color'],
                      'marker': 'x', 'ms': 6,  'alpha': 0.5},
}

regime_Ea_cases = [
    ('surface', P_for_surface_Ea, CURVE_STYLES['surface_region']['color'], 'o'),
    ('oxide',   P_for_oxide_Ea,   CURVE_STYLES['oxide_region']['color'],   's'),
    ('metal',   P_for_metal_Ea,   CURVE_STYLES['metal_region']['color'],   '^'),
]

plotted_class_labels = set()

for label, P_fixed, color, marker in regime_Ea_cases:
    if P_fixed is None:
        arrhenius_results[label] = None
        continue

    J_row     = []
    fs_row    = []
    fo_row    = []
    fm_row    = []
    class_row = []

    for T_K in temperatures_K:
        try:
            r = solve_steady_state_flux_direct(
                P_fixed, P_down, L_m,
                get_k_diss_at_T(T_K), get_K_eq_at_T(T_K),
                get_D_ox_at_T(T_K),   get_K_ox_at_T(T_K), L_ox,
                get_D_m_at_T(T_K),    get_K_s_m_at_T(T_K),
            )
            J_row.append(r['J_ss'])
            fs = r['resistances']['fraction_surface']
            fo = r['resistances']['fraction_oxide']
            fm = r['resistances']['fraction_metal']
            fs_row.append(fs)
            fo_row.append(fo)
            fm_row.append(fm)
            class_row.append(classify_regime(fs, fo, fm))
        except Exception:
            J_row.append(np.nan)
            fs_row.append(np.nan)
            fo_row.append(np.nan)
            fm_row.append(np.nan)
            class_row.append('mixed')

    J_arr_T   = np.array(J_row)
    class_arr = np.array(class_row)

    # Pure mask — only points classified as pure target regime
    pure = (class_arr == f'pure_{label}') & ~np.isnan(J_arr_T) & (J_arr_T > 0)

    # Plot all classified points with distinct size/alpha/marker
    for cls, sty in regime_class_styles.items():
        mask = (class_arr == cls) & ~np.isnan(J_arr_T) & (J_arr_T > 0)
        if not np.any(mask):
            continue
        cls_label = cls.replace('_', ' ').capitalize()
        ax2.semilogy(
            inv_T[mask], J_arr_T[mask],
            ls='none',
            color=sty['color'],
            marker=sty['marker'],
            ms=sty['ms'],
            alpha=sty['alpha'],
            label=cls_label if cls_label not in plotted_class_labels
                  else '_nolegend_'
        )
        plotted_class_labels.add(cls_label)

    # Arrhenius fit on pure points only
    if np.sum(pure) < 3:
        print(f"  ⚠  {label}: only {np.sum(pure)} pure T points "
              f"(fraction > 90%) — fit skipped")
        arrhenius_results[label] = None
        continue

    slope, intercept, r_val, *_ = stats.linregress(
        inv_T[pure], np.log(J_arr_T[pure])
    )
    E_a = -slope * R * 1000

    arrhenius_results[label] = {
        'J':            J_arr_T,
        'valid':        pure,
        'slope':        slope,
        'intercept':    intercept,
        'r_sq':         r_val**2,
        'E_a':          E_a,
        'P_fixed':      P_fixed,
        'color':        color,
        'n_pure':       int(np.sum(pure)),
        'T_pure_range': (temperatures_K[pure].min() - 273.15,
                         temperatures_K[pure].max() - 273.15),
        'class_arr':    class_arr,
    }

    # Fitted line through pure points only
    ax2.semilogy(
        inv_T[pure],
        np.exp(slope * inv_T[pure] + intercept),
        color=color, ls='--', lw=2.0, alpha=0.9,
        label=f'Fit {label}  (pure pts)'
    )

# ── Dual x-axis ───────────────────────────────────────────────────────────────
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
T_ticks_C = np.array([200, 300, 400, 500, 600, 700, 800, 900])
T_ticks_K = T_ticks_C + 273.15
in_range  = (T_ticks_K >= T_RANGE[0]) & (T_ticks_K <= T_RANGE[1])
if np.any(in_range):
    ax2_top.set_xticks(1000.0 / T_ticks_K[in_range])
    ax2_top.set_xticklabels([f'{t}' for t in T_ticks_C[in_range]])
ax2_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax2_top.tick_params(labelsize=STYLE['fontsize_tick'])

# ── Annotation box ────────────────────────────────────────────────────────────
E_refs    = {'surface': E_REF_SURFACE, 'oxide': E_REF_OXIDE, 'metal': E_REF_METAL}
ref_names = {'surface': 'E_diss', 'oxide': 'E_D_ox+H_sol', 'metal': 'E_D+H_s'}
annot     = []

for lbl in ('surface', 'oxide', 'metal'):
    res = arrhenius_results.get(lbl)
    E_r = E_refs[lbl]
    rn  = ref_names[lbl]
    if res is None:
        annot.append(f'{lbl.capitalize()}: not extractable')
        annot.append(f'  Expected ({rn}): {E_r/1000:.1f} kJ/mol')
    else:
        match = abs(res['E_a'] - E_r) < 3000
        annot.append(f'{lbl.capitalize()}-limited:')
        annot.append(f'  Extracted: {res["E_a"]/1000:.1f} kJ/mol')
        annot.append(f'  Expected ({rn}): {E_r/1000:.1f} kJ/mol')
        annot.append(f'  R²={res["r_sq"]:.4f}  {"✓" if match else "⚠"}')
        annot.append(f'  T range: {res["T_pure_range"][0]:.0f}–'
                     f'{res["T_pure_range"][1]:.0f}°C  '
                     f'({res["n_pure"]} pure pts)')
    annot.append('')

ax2.text(
    0.97, 0.97, '\n'.join(annot).rstrip(),
    transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='top', ha='right',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9)
)

ax2.set_xlabel('1000/T (K⁻¹)',             fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('Flux $J_{ss}$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax2.set_title('(B) Arrhenius: Apparent $E_a$ per Regime',
              fontsize=STYLE['fontsize_title'])
ax2.legend(fontsize=STYLE['fontsize_legend']-2, loc='lower left',
           ncol=2, title='Symbol size = purity')
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (C) RESISTANCE FRACTIONS vs PRESSURE
# ─────────────────────────────────────────────────────────────────────────────
ax3 = axes[1, 0]

ax3.semilogx(
    P_RANGE_ARR[valid], frac_surface_arr[valid] * 100,
    color=CURVE_STYLES['surface_region']['color'],
    lw=STYLE['linewidth'], label='Surface (dissociation)'
)
ax3.semilogx(
    P_RANGE_ARR[valid], frac_oxide_arr[valid] * 100,
    color=CURVE_STYLES['oxide_region']['color'],
    lw=STYLE['linewidth'], label='Oxide (diffusion)'
)
ax3.semilogx(
    P_RANGE_ARR[valid], frac_metal_arr[valid] * 100,
    color=CURVE_STYLES['metal_region']['color'],
    lw=STYLE['linewidth'], label='Metal (diffusion)'
)
ax3.axhline(
    50,
    color=CURVE_STYLES['threshold']['color'],
    ls=CURVE_STYLES['threshold']['ls'],
    lw=CURVE_STYLES['threshold']['lw'],
    alpha=CURVE_STYLES['threshold']['alpha'],
    label='50% threshold'
)

# Crossover annotations — up to two for three-resistance system
if P_cross_surf_oxide is not None:
    ax3.axvline(P_cross_surf_oxide, color='gray', ls=':', lw=1.5, alpha=0.7)
    ax3.text(
        P_cross_surf_oxide * 1.5, 53,
        f'S→O\n{P_cross_surf_oxide:.1e} Pa',
        fontsize=STYLE['fontsize_annotation']-2, color='gray'
    )
if P_cross_oxide_metal is not None:
    ax3.axvline(P_cross_oxide_metal, color='gray', ls=':', lw=1.5, alpha=0.7)
    ax3.text(
        P_cross_oxide_metal * 1.5, 53,
        f'O→M\n{P_cross_oxide_metal:.1e} Pa',
        fontsize=STYLE['fontsize_annotation']-2, color='gray'
    )

if P_cross_surf_oxide is None and P_cross_oxide_metal is None:
    dominant = ('Surface' if n_surface >= n_oxide and n_surface >= n_metal
                else 'Oxide' if n_oxide >= n_metal else 'Metal')
    ax3.text(
        0.05, 0.50,
        f'No crossover in range\n{dominant}-limited throughout',
        transform=ax3.transAxes,
        fontsize=STYLE['fontsize_annotation'],
        va='center',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )

ax3.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Resistance Fraction (%)',          fontsize=STYLE['fontsize_axis'])
ax3.set_title('(C) Rate-Limiting Step Analysis',   fontsize=STYLE['fontsize_title'])
ax3.set_ylim(0, 105)
ax3.legend(fontsize=STYLE['fontsize_legend'], loc='center right')
ax3.grid(True, alpha=STYLE['grid_alpha'])
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (D) LIMIT CHECK — k_diss → ∞ recovers series resistance result
# ─────────────────────────────────────────────────────────────────────────────
ax4 = axes[1, 1]

# Series resistance reference (no surface kinetics)
# J = (sqrt(P_up) - sqrt(P_down)) / (1/alpha + 1/beta)
J_series = ((np.sqrt(P_RANGE_ARR) - np.sqrt(P_down)) /
            (1.0/alpha + 1.0/beta))

# L2+L6 at high k_diss
J_high_k = []
for P_up in P_RANGE_ARR:
    try:
        r = solve_steady_state_flux_direct(
            P_up, P_down, L_m,
            1e-3, K_eq,
            D_ox, K_ox, L_ox,
            D_m, K_s_m
        )
        J_high_k.append(r['J_ss'])
    except Exception:
        J_high_k.append(np.nan)
J_high_k = np.array(J_high_k)

valid_parity = ~np.isnan(J_high_k) & (J_series > 0) & (J_high_k > 0)

ax4.loglog(
    J_series[valid_parity], J_high_k[valid_parity],
    'o', color=CURVE_STYLES['L2_L6']['color'],
    markersize=6, alpha=0.8,
    label='L2+L6  (k_diss = 1e-3)'
)

J_min = min(J_series[valid_parity].min(), J_high_k[valid_parity].min()) * 0.5
J_max = max(J_series[valid_parity].max(), J_high_k[valid_parity].max()) * 2.0
ax4.loglog(
    [J_min, J_max], [J_min, J_max],
    color=CURVE_STYLES['parity']['color'],
    ls=CURVE_STYLES['parity']['ls'],
    lw=CURVE_STYLES['parity']['lw'],
    alpha=CURVE_STYLES['parity']['alpha'],
    label='Perfect agreement'
)

rel_errors = (np.abs(J_high_k[valid_parity] - J_series[valid_parity])
              / J_series[valid_parity] * 100)
max_err   = rel_errors.max()
mean_err  = rel_errors.mean()
box_color = 'lightgreen' if max_err < 1.0 else 'lightyellow'

ax4.text(
    0.05, 0.95,
    f'k_diss = 1e-3\n'
    f'Max error  = {max_err:.2e}%\n'
    f'Mean error = {mean_err:.2e}%\n'
    f'{"✓ VALIDATED" if max_err < 1.0 else "⚠ CHECK k_diss"}',
    transform=ax4.transAxes, fontsize=STYLE['fontsize_annotation'],
    va='top',
    bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.9)
)

ax4.set_xlabel('Series-Resistance Flux $J_{ref}$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('L2+L6 Flux $J_{L2+L6}$ (mol/m²/s)',           fontsize=STYLE['fontsize_axis'])
ax4.set_title(r'(D) Limit Check: $k_{diss} \to \infty$ → Series Resistance',
              fontsize=STYLE['fontsize_title'])
ax4.set_aspect('equal', adjustable='box')
ax4.legend(fontsize=STYLE['fontsize_legend'], loc='lower right')
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.tick_params(labelsize=STYLE['fontsize_tick'])

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

print("\n" + "=" * 60)
print("L2+L6 FIGURE 1 — VALIDATION SUMMARY")
print("=" * 60)
print(f"  Net log-log slope:       {slope_net:.4f}")
print(f"  θ range:                 "
      f"{np.nanmin(theta_arr):.4f} → {np.nanmax(theta_arr):.4f}")
print(f"  Series limit error:      {max_err:.2e}%  "
      f"({'✓' if max_err < 1.0 else '⚠'})")
print(f"\n  Arrhenius results:")
for lbl in ('surface', 'oxide', 'metal'):
    res = arrhenius_results.get(lbl)
    E_r = E_refs[lbl]
    if res is None:
        print(f"    {lbl.capitalize()}: not extractable")
    else:
        match = abs(res['E_a'] - E_r) < 3000
        print(f"    {lbl.capitalize()}: extracted={res['E_a']/1000:.1f} "
              f"expected={E_r/1000:.1f} kJ/mol  "
              f"{'✓' if match else '⚠'}")
print("=" * 60)
```

```python
# =============================================================================
# CELL 5 — FIGURE 2: Extended Analysis (2×2)
# =============================================================================
"""
Figure 2 explores how system geometry and temperature affect the three-way
resistance balance.

Panel (E) — α/β Ratio Sweep
    Question: How does the oxide-to-metal permeance ratio determine which
    diffusion step is rate-limiting?
    Fixes surface kinetics and sweeps alpha/beta from 0.01 to 100 by varying
    alpha (oxide permeance) at fixed beta (metal permeance). At alpha/beta << 1
    the oxide is the bottleneck. At alpha/beta >> 1 the metal is the
    bottleneck. Plots J_ss vs P_up for several alpha/beta values on one axes.
    Each curve is labeled with its alpha/beta value. Shows that the
    surface-to-diffusion crossover pressure is independent of alpha/beta
    but the crossover between oxide-limited and metal-limited depends
    entirely on their ratio.

Panel (F) — Thickness Trade-off
    Question: If total membrane thickness is fixed, what is the optimal
    split between oxide and metal thickness?
    Co-varies L_ox and L_m at fixed total thickness L_total = L_ox + L_m.
    Plots J_ss vs the oxide fraction f_ox = L_ox / L_total. Shows there is
    a minimum flux (maximum resistance) at a specific oxide fraction, and
    that the minimum shifts with pressure and temperature.

Panel (G) — Surface Coverage θ vs Pressure
    Question: How does the addition of a metal block (vs L2a+L6) change
    the surface coverage behaviour?
    Plots steady-state θ alongside the Langmuir equilibrium isotherm.
    Compares θ_ss from L2+L6 against θ_ss from L2a+L6 at the same surface
    kinetics to isolate the effect of the metal block on surface coverage.
    The metal block provides an additional pathway for H to leave the
    surface, depleting θ further below equilibrium in the surface-limited
    regime.

Panel (H) — Rate-Limiting Map vs Temperature
    Question: How does the three-way balance shift with temperature?
    Runs the full pressure sweep at 5 temperatures spread across T_RANGE.
    Plots surface resistance fraction (%) vs pressure for each temperature
    using a plasma colormap. The crossover pressure shift direction is
    annotated with the physical reason from config activation energies.

Figure 3 (standalone) — Stacked Resistance Area
    Question: How do all three resistances share the total resistance burden
    visually across the full pressure range?
    Plots R_surface, R_oxide, and R_metal as a stacked area filling 100%
    at every pressure point. Unlike the line plot in panel (C), this shows
    all three fractions simultaneously as filled bands — any growth in one
    comes visually at the expense of the others. The dominant regime label
    is annotated at low, mid, and high pressure within the bands.
"""

from calculations.surface_kinetics import solve_steady_state_flux_L2aL6

fig2, axes2 = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig2.suptitle(
    f'L2+L6: Extended Analysis — Geometry & Temperature Effects\n'
    f'{OXIDE_KEY} + {METAL_KEY}',
    fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98
)

# ─────────────────────────────────────────────────────────────────────────────
# (E) α/β RATIO SWEEP
# ─────────────────────────────────────────────────────────────────────────────
ax_E = axes2[0, 0]

# Sweep alpha/beta by scaling alpha at fixed beta
ab_ratios  = [0.01, 0.1, 1.0, 10.0, 100.0]
ab_colors  = plt.cm.viridis(np.linspace(0.1, 0.9, len(ab_ratios)))

for ab, col in zip(ab_ratios, ab_colors):
    alpha_test = ab * beta    # scale alpha to hit target ratio
    J_ab = []
    for P_up in P_RANGE_ARR:
        try:
            # Equivalent D_ox for this alpha: D_ox_test = alpha_test * L_ox / K_ox
            D_ox_test = alpha_test * L_ox / K_ox
            r = solve_steady_state_flux_direct(
                P_up, P_down, L_m,
                k_diss, K_eq,
                D_ox_test, K_ox, L_ox,
                D_m, K_s_m
            )
            J_ab.append(r['J_ss'])
        except Exception:
            J_ab.append(np.nan)

    J_ab   = np.array(J_ab)
    v_ab   = ~np.isnan(J_ab) & (J_ab > 0)
    if not np.any(v_ab):
        continue

    ax_E.loglog(
        P_RANGE_ARR[v_ab], J_ab[v_ab],
        color=col, lw=STYLE['linewidth'],
        label=f'α/β = {ab:.2g}'
    )

ax_E.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax_E.set_ylabel('Flux $J_{ss}$ (mol/m²/s)',         fontsize=STYLE['fontsize_axis'])
ax_E.set_title('(E) α/β Ratio Sweep',               fontsize=STYLE['fontsize_title'])
ax_E.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper left',
            title='Oxide/Metal permeance')
ax_E.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax_E.tick_params(labelsize=STYLE['fontsize_tick'])
ax_E.text(
    0.05, 0.05,
    'α/β << 1 → oxide bottleneck\nα/β >> 1 → metal bottleneck',
    transform=ax_E.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='bottom',
    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
)

# ─────────────────────────────────────────────────────────────────────────────
# (F) THICKNESS TRADE-OFF
# ─────────────────────────────────────────────────────────────────────────────
ax_F = axes2[0, 1]

L_total     = L_ox + L_m
f_ox_sweep  = np.linspace(0.01, 0.99, 50)   # oxide fraction of total
P_cases     = [1e2, 1e5, 1e8]
P_colors_tf = [CURVE_STYLES['surface_region']['color'],
               CURVE_STYLES['oxide_region']['color'],
               CURVE_STYLES['metal_region']['color']]

for P_fixed, col in zip(P_cases, P_colors_tf):
    J_tf = []
    for f in f_ox_sweep:
        L_ox_t = f * L_total
        L_m_t  = (1 - f) * L_total
        try:
            r = solve_steady_state_flux_direct(
                P_fixed, P_down, L_m_t,
                k_diss, K_eq,
                D_ox, K_ox, L_ox_t,
                D_m, K_s_m
            )
            J_tf.append(r['J_ss'])
        except Exception:
            J_tf.append(np.nan)

    J_tf  = np.array(J_tf)
    v_tf  = ~np.isnan(J_tf) & (J_tf > 0)
    if not np.any(v_tf):
        continue

    ax_F.semilogy(
        f_ox_sweep[v_tf] * 100, J_tf[v_tf],
        color=col, lw=STYLE['linewidth'],
        label=f'P = {P_fixed:.0e} Pa'
    )

    # Mark minimum flux point
    min_idx = np.nanargmin(J_tf)
    if not np.isnan(J_tf[min_idx]):
        ax_F.plot(
            f_ox_sweep[min_idx] * 100, J_tf[min_idx],
            'v', color=col, ms=9
        )

ax_F.set_xlabel('Oxide Fraction $f_{ox}$ = $L_{ox}/L_{total}$ (%)',
                fontsize=STYLE['fontsize_axis'])
ax_F.set_ylabel('Flux $J_{ss}$ (mol/m²/s)',  fontsize=STYLE['fontsize_axis'])
ax_F.set_title(f'(F) Thickness Trade-off  ($L_{{total}}$ = {L_total*1e3:.1f} mm)',
               fontsize=STYLE['fontsize_title'])
ax_F.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper center')
ax_F.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax_F.tick_params(labelsize=STYLE['fontsize_tick'])
ax_F.text(
    0.05, 0.05,
    '▼ = minimum flux (maximum resistance)',
    transform=ax_F.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='bottom',
    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
)

# ─────────────────────────────────────────────────────────────────────────────
# (G) SURFACE COVERAGE θ — L2+L6 vs L2a+L6
# ─────────────────────────────────────────────────────────────────────────────
ax_G = axes2[1, 0]

# L2a+L6 coverage for comparison (no metal block)
theta_L2a = []
for P_up in P_RANGE_ARR:
    try:
        r2a = solve_steady_state_flux_L2aL6(
            P_up, P_down,
            k_diss, K_eq, D_ox, K_ox, L_ox
        )
        theta_L2a.append(r2a['theta'])
    except Exception:
        theta_L2a.append(np.nan)
theta_L2a = np.array(theta_L2a)

# Equilibrium Langmuir isotherm
theta_eq = (np.sqrt(K_eq * P_RANGE_ARR) /
            (1 + np.sqrt(K_eq * P_RANGE_ARR)))

ax_G.semilogx(
    P_RANGE_ARR[valid], theta_arr[valid],
    color=CURVE_STYLES['L2_L6']['color'],
    ls=CURVE_STYLES['L2_L6']['ls'],
    lw=STYLE['linewidth'],
    marker=CURVE_STYLES['L2_L6']['marker'],
    ms=CURVE_STYLES['L2_L6']['ms'],
    markevery=8,
    label='θ  L2+L6 (oxide + metal)'
)

v2a = ~np.isnan(theta_L2a)
ax_G.semilogx(
    P_RANGE_ARR[v2a], theta_L2a[v2a],
    color=CURVE_STYLES['L2a_L6']['color'],
    ls=CURVE_STYLES['L2a_L6']['ls'],
    lw=STYLE['linewidth'],
    marker=CURVE_STYLES['L2a_L6']['marker'],
    ms=CURVE_STYLES['L2a_L6']['ms'],
    markevery=8,
    label='θ  L2a+L6 (oxide only)'
)

ax_G.semilogx(
    P_RANGE_ARR, theta_eq,
    color='gray', ls='-.', lw=1.5, alpha=0.7,
    label=r'θ$_{eq}$ (Langmuir isotherm)'
)
ax_G.axhline(1.0, color='red',    ls='--', lw=1.5, alpha=0.4, label='θ = 1')
ax_G.axhline(0.5, color='orange', ls='--', lw=1.5, alpha=0.4, label='θ = 0.5')

ax_G.text(
    0.05, 0.95,
    'L2+L6 θ < L2a+L6 θ:\nmetal block provides extra\n'
    'pathway → more surface depletion',
    transform=ax_G.transAxes, fontsize=STYLE['fontsize_annotation']-1,
    va='top',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
)

ax_G.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax_G.set_ylabel('Surface Coverage $\\theta$',       fontsize=STYLE['fontsize_axis'])
ax_G.set_title('(G) Surface Coverage: L2+L6 vs L2a+L6',
               fontsize=STYLE['fontsize_title'])
ax_G.set_ylim(-0.05, 1.1)
ax_G.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax_G.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax_G.tick_params(labelsize=STYLE['fontsize_tick'])

# ─────────────────────────────────────────────────────────────────────────────
# (H) RATE-LIMITING MAP vs TEMPERATURE
# ─────────────────────────────────────────────────────────────────────────────
ax_H = axes2[1, 1]

T_map_values        = np.linspace(T_RANGE[0], T_RANGE[1], 5)
T_colors_map        = plt.cm.plasma(np.linspace(0.1, 0.9, len(T_map_values)))
crossover_pressures = []

for T_K, col in zip(T_map_values, T_colors_map):
    fs_row = []
    for P_up in P_RANGE_ARR:
        try:
            r = solve_steady_state_flux_direct(
                P_up, P_down, L_m,
                get_k_diss_at_T(T_K), get_K_eq_at_T(T_K),
                get_D_ox_at_T(T_K),   get_K_ox_at_T(T_K), L_ox,
                get_D_m_at_T(T_K),    get_K_s_m_at_T(T_K),
            )
            fs_row.append(r['resistances']['fraction_surface'] * 100)
        except Exception:
            fs_row.append(np.nan)

    fs_row    = np.array(fs_row)
    valid_row = ~np.isnan(fs_row)
    if not np.any(valid_row):
        continue

    ax_H.semilogx(
        P_RANGE_ARR[valid_row], fs_row[valid_row],
        color=col, lw=STYLE['linewidth'],
        label=f'T = {T_K-273.15:.0f}°C'
    )

    sc = np.where(np.diff(np.sign(fs_row[valid_row] - 50)))[0]
    if len(sc) > 0:
        idx  = sc[0]
        fs_v = fs_row[valid_row]
        if idx + 1 < len(fs_v) and fs_v[idx] > 50 and fs_v[idx+1] < 50:
            crossover_pressures.append((T_K, P_RANGE_ARR[valid_row][idx]))

ax_H.axhline(
    50,
    color=CURVE_STYLES['threshold']['color'],
    ls=CURVE_STYLES['threshold']['ls'],
    lw=CURVE_STYLES['threshold']['lw'],
    alpha=CURVE_STYLES['threshold']['alpha'],
    label='50% threshold'
)

if len(crossover_pressures) >= 2:
    T_low,  P_co_low  = crossover_pressures[0]
    T_high, P_co_high = crossover_pressures[-1]
    direction = 'lower P as T ↑' if P_co_high < P_co_low else 'higher P as T ↑'
    ax_H.text(
        0.05, 0.05,
        f'Crossover shifts to {direction}\n'
        f'E_diss={E_REF_SURFACE/1000:.0f} vs '
        f'E_D+Hs={E_REF_METAL/1000:.0f} kJ/mol',
        transform=ax_H.transAxes,
        fontsize=STYLE['fontsize_annotation']-1, va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )
elif len(crossover_pressures) == 0:
    ax_H.text(
        0.05, 0.05,
        'No surface/diffusion crossover\nfound across T range',
        transform=ax_H.transAxes,
        fontsize=STYLE['fontsize_annotation']-1, va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9)
    )

ax_H.set_xlabel('Upstream Pressure $P_{up}$ (Pa)',  fontsize=STYLE['fontsize_axis'])
ax_H.set_ylabel('Surface Resistance Fraction (%)',  fontsize=STYLE['fontsize_axis'])
ax_H.set_title('(H) Rate-Limiting Map vs Temperature',
               fontsize=STYLE['fontsize_title'])
ax_H.set_ylim(0, 105)
ax_H.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper right',
            title='Temperature')
ax_H.grid(True, alpha=STYLE['grid_alpha'])
ax_H.tick_params(labelsize=STYLE['fontsize_tick'])

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

print("\n" + "=" * 60)
print("L2+L6 FIGURE 2 — SUMMARY")
print("=" * 60)
print(f"  α/β at operating conditions: {alpha/beta:.3e}")
if crossover_pressures:
    print(f"  Surface/diffusion crossover across T sweep:")
    for T_K, P_co in crossover_pressures:
        print(f"    T = {T_K-273.15:.0f}°C  →  P_crossover = {P_co:.2e} Pa")
print("=" * 60)

# =============================================================================
# FIGURE 3 — Stacked Resistance Area (standalone)
# =============================================================================
fig3, ax_stack = plt.subplots(1, 1, figsize=(8, 5))
fig3.suptitle(
    f'L2+L6: Stacked Resistance Area\n'
    f'{OXIDE_KEY} + {METAL_KEY}  |  T = {T_op-273.15:.0f}°C',
    fontsize=STYLE['fontsize_title'], fontweight='bold'
)

fs = frac_surface_arr[valid] * 100
fo = frac_oxide_arr[valid]   * 100
fm = frac_metal_arr[valid]   * 100
Pv = P_RANGE_ARR[valid]

ax_stack.fill_between(
    Pv, 0, fs,
    color=CURVE_STYLES['surface_region']['color'],
    alpha=0.8, label='Surface'
)
ax_stack.fill_between(
    Pv, fs, fs + fo,
    color=CURVE_STYLES['oxide_region']['color'],
    alpha=0.8, label='Oxide'
)
ax_stack.fill_between(
    Pv, fs + fo, fs + fo + fm,
    color=CURVE_STYLES['metal_region']['color'],
    alpha=0.8, label='Metal'
)

ax_stack.set_xscale('log')
ax_stack.set_ylim(0, 100)
ax_stack.axhline(50, color='white', ls='--', lw=1.0, alpha=0.6)

# Dominant regime label at three representative pressures
for P_annot in [Pv[len(Pv)//8], Pv[len(Pv)//2], Pv[7*len(Pv)//8]]:
    idx = np.argmin(np.abs(Pv - P_annot))
    dom = rl_arr[valid][idx]
    ax_stack.text(
        P_annot, 50, dom.upper(),
        ha='center', va='center',
        fontsize=STYLE['fontsize_annotation'],
        color='white', fontweight='bold'
    )

ax_stack.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax_stack.set_ylabel('Resistance Fraction (%)',          fontsize=STYLE['fontsize_axis'])
ax_stack.set_title('Stacked Resistance Fractions',      fontsize=STYLE['fontsize_title'])
ax_stack.legend(fontsize=STYLE['fontsize_legend'], loc='upper right')
ax_stack.grid(True, alpha=STYLE['grid_alpha'])
ax_stack.tick_params(labelsize=STYLE['fontsize_tick'])

plt.tight_layout()
plt.show()
```

```python
# =============================================================================
# CELL 6 — DATAFRAME OUTPUT
# =============================================================================
df = pd.DataFrame(rows)

df_display = df.copy()
for col in df_display.columns:
    if col in ('Rate-Limiting', 'Error'):
        continue
    if df_display[col].dtype == float:
        if 'fraction' in col.lower():
            df_display[col] = df_display[col].round(2)
        elif 'theta' in col.lower():
            df_display[col] = df_display[col].round(6)
        else:
            df_display[col] = df_display[col].apply(
                lambda x: f'{x:.4e}' if pd.notna(x) else 'NaN'
            )

n_surface_df = (df['Rate-Limiting'] == 'SURFACE').sum()
n_oxide_df   = (df['Rate-Limiting'] == 'OXIDE').sum()
n_metal_df   = (df['Rate-Limiting'] == 'METAL').sum()
n_mixed_df   = (df['Rate-Limiting'] == 'MIXED').sum()
n_error_df   = (df['Rate-Limiting'] == 'ERROR').sum()

print("=" * 60)
print("L2+L6 RESULTS DATAFRAME")
print("=" * 60)
print(f"  Oxide:       {OXIDE_KEY}")
print(f"  Metal:       {METAL_KEY}")
print(f"  T:           {T_op-273.15:.0f}°C  ({T_op} K)")
print(f"  L_ox:        {L_ox*1e6:.1f} μm")
print(f"  L_m:         {L_m*1e3:.1f} mm")
print(f"  k_diss:      {k_diss:.3e}  mol/m²/s/Pa")
print(f"  K_eq:        {K_eq:.3e}  Pa⁻¹")
print(f"  α:           {alpha:.3e}  mol/m²/s/Pa^0.5")
print(f"  β:           {beta:.3e}  mol/m²/s/Pa^0.5")
print(f"  α/β:         {alpha/beta:.3e}")
print(f"\n  Total rows:        {len(df)}")
print(f"  Surface-limited:   {n_surface_df}")
print(f"  Oxide-limited:     {n_oxide_df}")
print(f"  Metal-limited:     {n_metal_df}")
print(f"  Mixed:             {n_mixed_df}")
if n_error_df > 0:
    print(f"  Errors:            {n_error_df}  ⚠")
print("=" * 60)

display(df_display)
```

```python

```

```python

```
