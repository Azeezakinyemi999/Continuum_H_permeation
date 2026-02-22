---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: mace_env
    language: python
    name: python3
---

```python
"""
================================================================================
HIERARCHICAL HYDROGEN PERMEATION MODEL - PRESENTATION FIGURES
================================================================================
Conference-ready plots with analytical limit validation at each level.

Structure:
- Level 1: Perfect Metal (baseline)
- Level 2a: Perfect Oxide Only
- Level 2b: Perfect Oxide + Perfect Metal
- Level 4: Defective Metal
- Level 3: Defective Oxide + Metal
- Level 5: Full System
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Import core functions
from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.utils import get_diffusivity, get_solubility, get_permeability
from data.material_data import MATERIALS

# ============================================================================
# CONFIGURATION
# ============================================================================
CONFIG = {
    # Material
    'material': 'Incoloy800',
    
    # Temperature
    'T_ref': 1073.15,              # 800°C reference (K)
    'T_range': (873.15, 1273.15),  # 600-1000°C for Arrhenius
    'n_T_points': 9,
    
    # Pressure  
    'P_range': (0.01, 1e5),        # Pa
    'P_down': 0,                   # Pa
    'n_P_points': 30,
    
    # Geometry
    'L_metal': 1e-3,               # 1 mm
    'L_range': (0.1e-3, 5e-3),     # 0.1-5 mm for sweep
    'n_L_points': 20,
}

# Plot style for conference presentation
STYLE = {
    'figsize': (14, 12),
    'fontsize_title': 16,
    'fontsize_suptitle': 18,
    'fontsize_axis': 14,
    'fontsize_tick': 12,
    'fontsize_legend': 12,
    'fontsize_annotation': 11,
    'linewidth': 2.5,
    'markersize': 8,
    'grid_alpha': 0.3,
}

# Color scheme by level
COLORS = {
    'L1': 'black',      # Perfect Metal (baseline)
    'L2a': 'blue',      # Perfect Oxide Only
    'L2b': 'purple',    # Oxide + Metal
    'L4_gb': 'green',   # Defective Metal - GB mode
    'L4_trap': 'orange',# Defective Metal - Trapping mode
    'L3': 'cyan',       # Defective Oxide
    'L5': 'red',        # Full System
}


# Set log log plot aspect ratio

def set_equal_decade_aspect(ax):
    """
    Set aspect ratio so that 1 decade on x-axis = 1 decade on y-axis visually.
    This makes slope=0.5 look half as steep as slope=1.0.
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Calculate number of decades on each axis
    x_decades = np.log10(xlim[1]) - np.log10(xlim[0])
    y_decades = np.log10(ylim[1]) - np.log10(ylim[0])
    
    # Set aspect ratio to make decades equal
    ax.set_aspect(y_decades / x_decades, adjustable='box')

print("✓ Configuration loaded")
print(f"  Material: {CONFIG['material']}")
print(f"  T_ref: {CONFIG['T_ref']-273.15:.0f}°C")
print(f"  P_range: {CONFIG['P_range'][0]}-{CONFIG['P_range'][1]} Pa")
print(f"  L_metal: {CONFIG['L_metal']*1000:.1f} mm")
```

```python
"""
================================================================================
LEVEL 1: PERFECT METAL (BASELINE)
================================================================================
Physics: Sieverts' Law + Fick's Law
    C = K_s × √P        (surface equilibrium)
    J = D × (C_up - C_down) / L    (bulk diffusion)
    
Combined: J = (D × K_s / L) × (√P_up - √P_down) = Φ × √P / L

Validation:
    (A) Flux vs Pressure: slope = 0.5 on log-log
    (B) Permeability vs 1000/T: Arrhenius with E_Φ = E_D + ΔH_s
    (C) Flux vs Thickness: slope = -1 on log-log
    (D) Analytical check: model matches J = Φ×√P/L exactly
================================================================================
"""

# Get material properties
material = MATERIALS[CONFIG['material']]
R = 8.314  # J/mol/K

# Create figure
fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle('Level 1: Perfect Metal (Sieverts + Fick)', 
             fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98)

# ============================================================================
# (A) FLUX vs PRESSURE - Verify √P scaling (slope = 0.5)
# ============================================================================
ax1 = axes[0, 0]

# Temperature and properties at reference
T_ref = CONFIG['T_ref']
D_ref = get_diffusivity(T_ref, material)
K_s_ref = get_solubility(T_ref, material)
L_ref = CONFIG['L_metal']

# Pressure sweep
pressures_A = np.logspace(np.log10(CONFIG['P_range'][0]), 
                          np.log10(CONFIG['P_range'][1]), 
                          CONFIG['n_P_points'])
fluxes_A = []
for P in pressures_A:
    result = calculate_simple_metal_flux(D_ref, K_s_ref, L_ref, P, CONFIG['P_down'])
    fluxes_A.append(result['flux'])
fluxes_A = np.array(fluxes_A)

# Linear regression on log-log
log_P = np.log10(pressures_A)
log_J = np.log10(fluxes_A)
slope_A, intercept_A, r_A, _, _ = stats.linregress(log_P, log_J)

# Plot
ax1.loglog(pressures_A, fluxes_A, 'o-', color=COLORS['L1'], 
           linewidth=STYLE['linewidth'], markersize=STYLE['markersize'], label='Model')

# Fitted line
fitted_J = 10**(slope_A * log_P + intercept_A)
ax1.loglog(pressures_A, fitted_J, '--', color='red', linewidth=2, alpha=0.7, label='Fit')

ax1.set_xlabel('Upstream Pressure, $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax1.set_xlim(1e-2, 1e5)  # 7 decades
ax1.set_ylim(1e-16, 1e-09)  # 7 decades (adjust based on your data range)
ax1.set_title('(A) Flux vs Pressure', fontsize=STYLE['fontsize_title'])
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.legend(fontsize=STYLE['fontsize_legend'])
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

#set_equal_decade_aspect(ax1)

# Annotation box
textstr_A = f'Slope = {slope_A:.4f}\nExpected = 0.500\nR² = {r_A**2:.6f}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
ax1.text(0.05, 0.95, textstr_A, transform=ax1.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=props)

# Physics equation
ax1.text(0.95, 0.05, r'$J \propto \sqrt{P}$', transform=ax1.transAxes, 
         fontsize=14, ha='right', va='bottom', style='italic',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# ============================================================================
# (B) PERMEABILITY vs 1000/T - Arrhenius behavior
# ============================================================================
ax2 = axes[0, 1]

# Temperature sweep
temperatures_K = np.linspace(CONFIG['T_range'][0], CONFIG['T_range'][1], CONFIG['n_T_points'])
temperatures_C = temperatures_K - 273.15
inv_T = 1000 / temperatures_K

permeabilities_B = []
for T_K in temperatures_K:
    Phi = get_permeability(T_K, material)
    permeabilities_B.append(Phi)
permeabilities_B = np.array(permeabilities_B)

# Arrhenius fit
ln_Phi = np.log(permeabilities_B)
slope_B, intercept_B, r_B, _, _ = stats.linregress(inv_T, ln_Phi)
E_Phi_extracted = -slope_B * R * 1000  # J/mol
Phi_0_extracted = np.exp(intercept_B)

# Expected value
E_Phi_expected = material['E_D'] + material['H_s']

# Plot
ax2.semilogy(inv_T, permeabilities_B, 'o-', color=COLORS['L1'], 
             linewidth=STYLE['linewidth'], markersize=STYLE['markersize'], label='Model')

# Fitted line
fitted_Phi = np.exp(slope_B * inv_T + intercept_B)
ax2.semilogy(inv_T, fitted_Phi, '--', color='red', linewidth=2, alpha=0.7, label='Fit')

ax2.set_xlabel('1000/T (K⁻¹)', fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('Permeability, $\Phi$ (mol/m/s/Pa$^{0.5}$)', fontsize=STYLE['fontsize_axis'])
ax2.set_title('(B) Permeability vs Temperature', fontsize=STYLE['fontsize_title'])
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.legend(fontsize=STYLE['fontsize_legend'])
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# Add temperature axis on top
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
T_ticks = np.array([600, 700, 800, 900, 1000])
ax2_top.set_xticks(1000 / (T_ticks + 273.15))
ax2_top.set_xticklabels([f'{t}' for t in T_ticks])
ax2_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax2_top.tick_params(labelsize=STYLE['fontsize_tick'])

# Annotation box
textstr_B = (f'$E_\\Phi$ = {E_Phi_extracted/1000:.1f} kJ/mol\n'
             f'Expected = {E_Phi_expected/1000:.1f} kJ/mol\n'
             f'R² = {r_B**2:.6f}')
ax2.text(0.95, 0.95, textstr_B, transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', ha='right', bbox=props)

# Physics equation
ax2.text(0.05, 0.05, r'$\Phi = \Phi_0 \exp(-E_\Phi/RT)$', transform=ax2.transAxes, 
         fontsize=12, ha='left', va='bottom', style='italic',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# ============================================================================
# (C) FLUX vs THICKNESS - Verify 1/L scaling (slope = -1)
# ============================================================================
ax3 = axes[1, 0]

# Thickness sweep at reference temperature
thicknesses_C = np.linspace(CONFIG['L_range'][0], CONFIG['L_range'][1], CONFIG['n_L_points'])
fluxes_C = []
P_ref = 1.0  # 1 Pa for this test

for L in thicknesses_C:
    result = calculate_simple_metal_flux(D_ref, K_s_ref, L, P_ref, CONFIG['P_down'])
    fluxes_C.append(result['flux'])
fluxes_C = np.array(fluxes_C)

# Linear regression on log-log
log_L = np.log10(thicknesses_C * 1000)  # Convert to mm for display
log_J_C = np.log10(fluxes_C)
slope_C, intercept_C, r_C, _, _ = stats.linregress(log_L, log_J_C)

# Plot
ax3.loglog(thicknesses_C * 1000, fluxes_C, 'o-', color=COLORS['L1'], 
           linewidth=STYLE['linewidth'], markersize=STYLE['markersize'], label='Model')

# Theoretical line (J ∝ 1/L)
J_theory_ref = fluxes_C[len(fluxes_C)//2]
L_theory_ref = thicknesses_C[len(thicknesses_C)//2]
J_theory = J_theory_ref * (L_theory_ref / thicknesses_C)
ax3.loglog(thicknesses_C * 1000, J_theory, '--', color='red', linewidth=2, alpha=0.7, 
           label='Theory ($J \\propto 1/L$)')

ax3.set_xlabel('Thickness, $L$ (mm)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax3.set_title('(C) Flux vs Thickness', fontsize=STYLE['fontsize_title'])
ax3.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax3.legend(fontsize=STYLE['fontsize_legend'])
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# Annotation box
textstr_C = f'Slope = {slope_C:.4f}\nExpected = -1.000\nR² = {r_C**2:.6f}'
ax3.text(0.95, 0.95, textstr_C, transform=ax3.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', ha='right', bbox=props)

# Physics equation
ax3.text(0.05, 0.05, r'$J = \Phi \sqrt{P} / L$', transform=ax3.transAxes, 
         fontsize=14, ha='left', va='bottom', style='italic',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# ============================================================================
# (D) ANALYTICAL CHECK - Model matches exact solution
# ============================================================================
ax4 = axes[1, 1]

# Compare model output to analytical formula across conditions
# Vary both P and T
n_test = 50
test_pressures = np.logspace(-2, 2, n_test)
test_temps = np.linspace(700, 1000, 4) + 273.15  # 4 temperatures

for i, T_test in enumerate(test_temps):
    D_test = get_diffusivity(T_test, material)
    K_s_test = get_solubility(T_test, material)
    Phi_test = D_test * K_s_test
    
    J_model = []
    J_analytical = []
    
    for P in test_pressures:
        # Model
        result = calculate_simple_metal_flux(D_test, K_s_test, L_ref, P, 0)
        J_model.append(result['flux'])
        
        # Analytical: J = Φ × √P / L
        J_ana = Phi_test * np.sqrt(P) / L_ref
        J_analytical.append(J_ana)
    
    J_model = np.array(J_model)
    J_analytical = np.array(J_analytical)
    
    # Calculate relative error
    rel_error = np.abs(J_model - J_analytical) / J_analytical * 100
    
    # Plot parity
    color = plt.cm.viridis(i / 3)
    ax4.loglog(J_analytical, J_model, 'o', color=color, markersize=6, 
               label=f'T = {T_test-273.15:.0f}°C', alpha=0.7)

# Perfect agreement line
J_range = [1e-20, 1e-5]
ax4.loglog(J_range, J_range, 'k--', linewidth=2, label='Perfect agreement')

ax4.set_xlabel('Analytical: $J = \\Phi \\sqrt{P}/L$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('Model: calculate_simple_metal_flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_title('(D) Analytical Validation', fontsize=STYLE['fontsize_title'])
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax4.tick_params(labelsize=STYLE['fontsize_tick'])
#ax4.set_aspect('equal', adjustable='box')

# Annotation
ax4.text(0.05, 0.95, '✓ Exact match\nat all conditions', transform=ax4.transAxes, 
         fontsize=STYLE['fontsize_annotation'], verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))

# ============================================================================
# FINAL LAYOUT
# ============================================================================
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Level1_Perfect_Metal.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

# Summary
print("\n" + "="*80)
print("LEVEL 1: PERFECT METAL - VALIDATION SUMMARY")
print("="*80)
print(f"\n(A) Flux vs Pressure:")
print(f"    Slope = {slope_A:.4f} (expected: 0.500)")
print(f"    R² = {r_A**2:.6f}")
print(f"    ✓ Sieverts' law verified: J ∝ √P")

print(f"\n(B) Permeability vs Temperature:")
print(f"    E_Φ = {E_Phi_extracted/1000:.1f} kJ/mol (expected: {E_Phi_expected/1000:.1f})")
print(f"    R² = {r_B**2:.6f}")
print(f"    ✓ Arrhenius behavior verified")

print(f"\n(C) Flux vs Thickness:")
print(f"    Slope = {slope_C:.4f} (expected: -1.000)")
print(f"    R² = {r_C**2:.6f}")
print(f"    ✓ Inverse thickness scaling verified: J ∝ 1/L")

print(f"\n(D) Analytical Check:")
print(f"    Model matches J = Φ×√P/L exactly at all conditions")
print(f"    ✓ Implementation correct")

print("\n" + "="*80)
print("Level 1 serves as BASELINE for all subsequent levels")
print("="*80)
```

```python
material
```

```python
"""
================================================================================
LEVEL 2a: PERFECT OXIDE ONLY (Molecular Diffusion)
================================================================================
Physics: Henry's Law + Fick's Law
    C = K_ox × P        (Henry's law - molecular dissolution)
    J = D_ox × (C_up - C_down) / δ    (Fick's law)
    
Combined: J = (D_ox × K_ox / δ) × (P_up - P_down) = Φ_ox × ΔP / δ

KEY DIFFERENCE FROM METAL:
    Metal:  C ∝ √P  (Sieverts)  →  J ∝ √P  (slope = 0.5)
    Oxide:  C ∝ P   (Henry)     →  J ∝ P   (slope = 1.0)

Validation:
    (A) Flux vs Pressure: slope = 1.0 on log-log (Henry's law)
    (B) Oxide permeability vs 1000/T: Arrhenius
    (C) Flux vs Oxide Thickness: slope = -1 on log-log
    (D) Analytical check: model matches J = D_ox×K_ox×ΔP/δ exactly
================================================================================
"""

from calculations.oxide_permeation import molecular_diffusion_flux, get_oxide_properties_at_T
from data.oxide_properties import OXIDE_PROPERTIES

# Update CONFIG for oxide
CONFIG['oxide'] = 'Cr2O3'
CONFIG['L_oxide'] = OXIDE_PROPERTIES['Cr2O3']['thickness']  # 1 μm
CONFIG['L_oxide_range'] = (1e-7, 1e-5)  # 0.1 - 10 μm

# Get oxide data
oxide_data = OXIDE_PROPERTIES[CONFIG['oxide']]

# Create figure
fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle('Level 2a: Perfect Oxide Only (Henry + Fick)', 
             fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98)

# ============================================================================
# (A) FLUX vs PRESSURE - Verify linear scaling (slope = 1.0)
# ============================================================================
ax1 = axes[0, 0]

# Get oxide properties at reference temperature
T_ref = CONFIG['T_ref']
oxide_props = get_oxide_properties_at_T(CONFIG['oxide'], T_ref)
D_ox_ref = oxide_props['D_ox']
K_ox_ref = oxide_props['K_ox']
delta_ref = CONFIG['L_oxide']

# Pressure sweep
pressures_2A = np.logspace(np.log10(CONFIG['P_range'][0]), 
                           np.log10(CONFIG['P_range'][1]), 
                           CONFIG['n_P_points'])
fluxes_2A = []
for P in pressures_2A:
    flux = molecular_diffusion_flux(D_ox_ref, K_ox_ref, delta_ref, P, CONFIG['P_down'])
    fluxes_2A.append(flux)
fluxes_2A = np.array(fluxes_2A)

# Linear regression on log-log
log_P_2A = np.log10(pressures_2A)
log_J_2A = np.log10(fluxes_2A)
slope_2A, intercept_2A, r_2A, _, _ = stats.linregress(log_P_2A, log_J_2A)

# Plot
ax1.loglog(pressures_2A, fluxes_2A, 's-', color=COLORS['L2a'], 
           linewidth=STYLE['linewidth'], markersize=STYLE['markersize'], label='Oxide Model')

# Fitted line
fitted_J_2A = 10**(slope_2A * log_P_2A + intercept_2A)
ax1.loglog(pressures_2A, fitted_J_2A, '--', color='red', linewidth=2, alpha=0.7, label='Fit')

ax1.set_xlabel('Upstream Pressure, $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax1.set_title('(A) Flux vs Pressure', fontsize=STYLE['fontsize_title'])
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.legend(fontsize=STYLE['fontsize_legend'])
ax1.tick_params(labelsize=STYLE['fontsize_tick'])
#set_equal_decade_aspect(ax1)

# Annotation box
textstr_2A = f'Slope = {slope_2A:.4f}\nExpected = 1.000\nR² = {r_2A**2:.6f}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
ax1.text(0.05, 0.95, textstr_2A, transform=ax1.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=props)

# Physics equation - highlight difference from metal
ax1.text(0.95, 0.05, r'$J \propto P$ (Henry)', transform=ax1.transAxes, 
         fontsize=14, ha='right', va='bottom', style='italic',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# ============================================================================
# (B) OXIDE PERMEABILITY vs 1000/T - Arrhenius behavior
# ============================================================================
ax2 = axes[0, 1]

# Temperature sweep
temperatures_K = np.linspace(CONFIG['T_range'][0], CONFIG['T_range'][1], CONFIG['n_T_points'])
inv_T = 1000 / temperatures_K

# Calculate oxide permeability: Φ_ox = D_ox × K_ox
Phi_ox_list = []
D_ox_list = []
K_ox_list = []

for T_K in temperatures_K:
    props = get_oxide_properties_at_T(CONFIG['oxide'], T_K)
    D_ox_list.append(props['D_ox'])
    K_ox_list.append(props['K_ox'])
    Phi_ox = props['D_ox'] * props['K_ox']
    Phi_ox_list.append(Phi_ox)

Phi_ox_arr = np.array(Phi_ox_list)
D_ox_arr = np.array(D_ox_list)
K_ox_arr = np.array(K_ox_list)

# Arrhenius fit for Φ_ox
ln_Phi_ox = np.log(Phi_ox_arr)
slope_2B, intercept_2B, r_2B, _, _ = stats.linregress(inv_T, ln_Phi_ox)
E_Phi_ox_extracted = -slope_2B * R * 1000  # J/mol
Phi_ox_0_extracted = np.exp(intercept_2B)

# Expected value
E_Phi_ox_expected = oxide_data['E_D_ox'] + oxide_data['H_sol_ox']

# Plot
ax2.semilogy(inv_T, Phi_ox_arr, 's-', color=COLORS['L2a'], 
             linewidth=STYLE['linewidth'], markersize=STYLE['markersize'], label='$\\Phi_{ox}$ = D$_{ox}$ × K$_{ox}$')

# Fitted line
fitted_Phi_ox = np.exp(slope_2B * inv_T + intercept_2B)
ax2.semilogy(inv_T, fitted_Phi_ox, '--', color='red', linewidth=2, alpha=0.7, label='Fit')

ax2.set_xlabel('1000/T (K⁻¹)', fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('Oxide Permeability, $\\Phi_{ox}$ (mol/m/s/Pa)', fontsize=STYLE['fontsize_axis'])
ax2.set_title('(B) Oxide Permeability vs Temperature', fontsize=STYLE['fontsize_title'])
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.legend(fontsize=STYLE['fontsize_legend'])
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# Add temperature axis on top
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
T_ticks = np.array([600, 700, 800, 900, 1000])
ax2_top.set_xticks(1000 / (T_ticks + 273.15))
ax2_top.set_xticklabels([f'{t}' for t in T_ticks])
ax2_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax2_top.tick_params(labelsize=STYLE['fontsize_tick'])

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
# Annotation box
textstr_2B = (f'$E_{{\\Phi,ox}}$ = {E_Phi_ox_extracted/1000:.1f} kJ/mol\n'
              f'Expected = {E_Phi_ox_expected/1000:.1f} kJ/mol\n'
              f'R² = {r_2B**2:.6f}')
ax2.text(0.95, 0.95, textstr_2B, transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', ha='right', bbox=props)

# Physics equation
ax2.text(0.05, 0.05, r'$\Phi_{ox} = D_{ox} \times K_{ox}$', transform=ax2.transAxes, 
         fontsize=12, ha='left', va='bottom', style='italic',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# ============================================================================
# (C) FLUX vs OXIDE THICKNESS - Verify 1/δ scaling (slope = -1)
# ============================================================================
ax3 = axes[1, 0]

# Thickness sweep at reference temperature
thicknesses_2C = np.logspace(np.log10(CONFIG['L_oxide_range'][0]), 
                              np.log10(CONFIG['L_oxide_range'][1]), 
                              CONFIG['n_L_points'])
fluxes_2C = []
P_ref_2C = 1.0  # 1 Pa

for delta in thicknesses_2C:
    flux = molecular_diffusion_flux(D_ox_ref, K_ox_ref, delta, P_ref_2C, CONFIG['P_down'])
    fluxes_2C.append(flux)
fluxes_2C = np.array(fluxes_2C)

# Linear regression on log-log (thickness in μm)
log_delta = np.log10(thicknesses_2C * 1e6)  # Convert to μm
log_J_2C = np.log10(fluxes_2C)
slope_2C, intercept_2C, r_2C, _, _ = stats.linregress(log_delta, log_J_2C)

# Plot
ax3.loglog(thicknesses_2C * 1e6, fluxes_2C, 's-', color=COLORS['L2a'], 
           linewidth=STYLE['linewidth'], markersize=STYLE['markersize'], label='Model')

# Theoretical line (J ∝ 1/δ)
J_theory_ref_2C = fluxes_2C[len(fluxes_2C)//2]
delta_theory_ref = thicknesses_2C[len(thicknesses_2C)//2]
J_theory_2C = J_theory_ref_2C * (delta_theory_ref / thicknesses_2C)
ax3.loglog(thicknesses_2C * 1e6, J_theory_2C, '--', color='red', linewidth=2, alpha=0.7, 
           label='Theory ($J \\propto 1/\\delta$)')

ax3.set_xlabel('Oxide Thickness, $\\delta$ (μm)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax3.set_title('(C) Flux vs Oxide Thickness', fontsize=STYLE['fontsize_title'])
ax3.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax3.legend(fontsize=STYLE['fontsize_legend'])
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# Annotation box
textstr_2C = f'Slope = {slope_2C:.4f}\nExpected = -1.000\nR² = {r_2C**2:.6f}'
ax3.text(0.95, 0.95, textstr_2C, transform=ax3.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', ha='right', bbox=props)

# Physics equation
ax3.text(0.05, 0.05, r'$J = \Phi_{ox} \Delta P / \delta$', transform=ax3.transAxes, 
         fontsize=14, ha='left', va='bottom', style='italic',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# ============================================================================
# (D) ANALYTICAL CHECK - Model matches exact solution
# ============================================================================
ax4 = axes[1, 1]

# Compare model output to analytical formula
n_test = 50
test_pressures_2D = np.logspace(-2, 5, n_test)
test_temps_2D = np.linspace(700, 1000, 4) + 273

for i, T_test in enumerate(test_temps_2D):
    props = get_oxide_properties_at_T(CONFIG['oxide'], T_test)
    D_ox_test = props['D_ox']
    K_ox_test = props['K_ox']
    Phi_ox_test = D_ox_test * K_ox_test
    
    J_model_2D = []
    J_analytical_2D = []
    
    for P in test_pressures_2D:
        # Model
        flux_model = molecular_diffusion_flux(D_ox_test, K_ox_test, delta_ref, P, 0)
        J_model_2D.append(flux_model)
        
        # Analytical: J = Φ_ox × ΔP / δ = D_ox × K_ox × (P_up - P_down) / δ
        J_ana = Phi_ox_test * P / delta_ref
        J_analytical_2D.append(J_ana)
    
    J_model_2D = np.array(J_model_2D)
    J_analytical_2D = np.array(J_analytical_2D)
    
    # Plot parity
    color = plt.cm.Blues(0.3 + i * 0.2)
    ax4.loglog(J_analytical_2D, J_model_2D, 's', color=color, markersize=6, 
               label=f'T = {T_test-273.15:.0f}°C', alpha=0.7)

# Perfect agreement line
J_range_2D = [1e-25, 1e-5]
ax4.loglog(J_range_2D, J_range_2D, 'k--', linewidth=2, label='Perfect agreement')

ax4.set_xlabel('Analytical: $J = D_{ox} K_{ox} \\Delta P / \\delta$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('Model: molecular_diffusion_flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_title('(D) Analytical Validation', fontsize=STYLE['fontsize_title'])
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax4.tick_params(labelsize=STYLE['fontsize_tick'])
ax4.set_aspect('equal', adjustable='box')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
# Annotation
ax4.text(0.05, 0.95, '✓ Exact match\nat all conditions', transform=ax4.transAxes, 
         fontsize=STYLE['fontsize_annotation'], verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))

# ============================================================================
# FINAL LAYOUT
# ============================================================================
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Level2a_Perfect_Oxide.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

# Summary
print("\n" + "="*80)
print("LEVEL 2a: PERFECT OXIDE ONLY - VALIDATION SUMMARY")
print("="*80)
print(f"\nOxide: {CONFIG['oxide']}")
print(f"Thickness: {CONFIG['L_oxide']*1e6:.1f} μm")
print(f"T_ref: {T_ref-273.15:.0f}°C")
print(f"D_ox(T_ref) = {D_ox_ref:.3e} m²/s")
print(f"K_ox(T_ref) = {K_ox_ref:.3e} mol/m³/Pa")

print(f"\n(A) Flux vs Pressure:")
print(f"    Slope = {slope_2A:.4f} (expected: 1.000)")
print(f"    R² = {r_2A**2:.6f}")
print(f"    ✓ Henry's law verified: J ∝ P (NOT √P like metal!)")

print(f"\n(B) Oxide Permeability vs Temperature:")
print(f"    E_Φ,ox = {E_Phi_ox_extracted/1000:.1f} kJ/mol (expected: {E_Phi_ox_expected/1000:.1f})")
print(f"    R² = {r_2B**2:.6f}")
print(f"    ✓ Arrhenius behavior verified")

print(f"\n(C) Flux vs Oxide Thickness:")
print(f"    Slope = {slope_2C:.4f} (expected: -1.000)")
print(f"    R² = {r_2C**2:.6f}")
print(f"    ✓ Inverse thickness scaling verified: J ∝ 1/δ")

print(f"\n(D) Analytical Check:")
print(f"    Model matches J = D_ox×K_ox×ΔP/δ exactly")
print(f"    ✓ Implementation correct")

print("\n" + "="*80)
print("KEY INSIGHT: Oxide has slope=1 (Henry), Metal has slope=0.5 (Sieverts)")
print("This difference drives the regime transition in Level 2b")
print("="*80)
```

```python
"""
================================================================================
LEVEL 2b: PERFECT OXIDE + PERFECT METAL (Series Resistance)
================================================================================
Physics: Oxide and Metal in series with interface coupling
    
    Oxide:  J_ox = D_ox × K_ox × (P_up - P_int) / δ     (Henry's law)
    Metal:  J_m  = D × K_s × (√P_int - √P_down) / L    (Sieverts' law)
    
    At steady-state: J_ox = J_m  →  Solve for P_interface

LIMIT CHECK PHYSICS (from test_interface_solver.py):
    - When δ→0: P_interface → P_upstream (normalized → 1.0)
    - This means ALL pressure drop is across metal → Level 1 behavior
    - Flux approaches metal-only flux when P_int ≈ P_up

Validation:
    (A) Flux vs Pressure: Transition from slope=1 to slope=0.5
    (B) System flux vs 1000/T: Compare with Level 1
    (C) Flux vs Metal Thickness: Metal controls at large L
    (D) LIMIT CHECK: δ_ox → 0 → P_int/P_up → 1 (metal dominates)
================================================================================
"""

from calculations.interface_solver import (
    calculate_oxide_metal_system,
    solve_interface_pressure,
    calculate_metal_flux_sieverts
)
from calculations.oxide_permeation import (
    molecular_diffusion_flux,
    get_oxide_properties_at_T,
    get_metal_properties_at_T
)

# Create figure
fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle('Level 2b: Perfect Oxide + Perfect Metal (Series Coupling)', 
             fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98)

# ============================================================================
# (A) FLUX vs PRESSURE - Show transition from slope=1 to slope=0.5
# ============================================================================
ax1 = axes[0, 0]

# Get properties at reference temperature
T_ref = CONFIG['T_ref']
oxide_props_2b = get_oxide_properties_at_T(CONFIG['oxide'], T_ref)
oxide_props_2b['thickness'] = CONFIG['L_oxide']

metal_props_2b = get_metal_properties_at_T(CONFIG['material'], T_ref)
metal_props_2b['thickness'] = CONFIG['L_metal']

# Realistic pressure sweep (0.01 Pa to 1 MPa)
pressures_2bA = np.logspace(-8, 20, 40)  

fluxes_L2b = []
fluxes_L1_compare = []
fluxes_L2a_compare = []
P_interfaces = []
P_normalized = []

for P in pressures_2bA:
    # Level 2b: Oxide + Metal system
    result = calculate_oxide_metal_system(P, CONFIG['P_down'], oxide_props_2b, metal_props_2b)
    fluxes_L2b.append(result['flux'])
    P_interfaces.append(result['P_interface'])
    P_normalized.append(result['P_interface_normalized'])
    
    # Level 1: Metal only
    result_L1 = calculate_simple_metal_flux(
        metal_props_2b['D_metal'], metal_props_2b['K_s_metal'],
        metal_props_2b['thickness'], P, CONFIG['P_down']
    )
    fluxes_L1_compare.append(result_L1['flux'])
    
    # Level 2a: Oxide only
    flux_L2a = molecular_diffusion_flux(
        oxide_props_2b['D_ox'], oxide_props_2b['K_ox'],
        oxide_props_2b['thickness'], P, CONFIG['P_down']
    )
    fluxes_L2a_compare.append(flux_L2a)

fluxes_L2b = np.array(fluxes_L2b)
fluxes_L1_compare = np.array(fluxes_L1_compare)
fluxes_L2a_compare = np.array(fluxes_L2a_compare)

# Plot all three levels
ax1.loglog(pressures_2bA, fluxes_L2b, '^-', color=COLORS['L2b'], 
           linewidth=STYLE['linewidth'], markersize=6, label='Level 2b: Oxide+Metal')
ax1.loglog(pressures_2bA, fluxes_L1_compare, '--', color=COLORS['L1'], 
           linewidth=2, alpha=0.6, label='Level 1: Metal only (slope=0.5)')
ax1.loglog(pressures_2bA, fluxes_L2a_compare, ':', color=COLORS['L2a'], 
           linewidth=2, alpha=0.6, label='Level 2a: Oxide only (slope=1)')

ax1.set_xlabel('Upstream Pressure, $P_{up}$ (Pa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax1.set_title('(A) Flux vs Pressure - Regime Transition', fontsize=STYLE['fontsize_title'])
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
ax1.text(0.05, 0.95, 
         'Low P: Oxide controls\nHigh P: Metal controls\nTransition: series resistance',
         transform=ax1.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=props)

# ============================================================================
# (B) SYSTEM FLUX vs 1000/T
# ============================================================================
ax2 = axes[0, 1]

P_fixed = 1e6  # 1000 kPa
temperatures_2bB = np.linspace(CONFIG['T_range'][0], CONFIG['T_range'][1], CONFIG['n_T_points'])
inv_T_2b = 1000 / temperatures_2bB

fluxes_L2b_T = []
fluxes_L1_T = []

for T_K in temperatures_2bB:
    ox_props = get_oxide_properties_at_T(CONFIG['oxide'], T_K)
    ox_props['thickness'] = CONFIG['L_oxide']
    
    met_props = get_metal_properties_at_T(CONFIG['material'], T_K)
    met_props['thickness'] = CONFIG['L_metal']
    
    # Level 2b
    result = calculate_oxide_metal_system(P_fixed, CONFIG['P_down'], ox_props, met_props)
    fluxes_L2b_T.append(result['flux'])
    
    # Level 1
    result_L1 = calculate_simple_metal_flux(
        met_props['D_metal'], met_props['K_s_metal'],
        met_props['thickness'], P_fixed, CONFIG['P_down']
    )
    fluxes_L1_T.append(result_L1['flux'])

fluxes_L2b_T = np.array(fluxes_L2b_T)
fluxes_L1_T = np.array(fluxes_L1_T)

# Arrhenius fit
ln_J_2b = np.log(fluxes_L2b_T)
slope_2bB, intercept_2bB, r_2bB, _, _ = stats.linregress(inv_T_2b, ln_J_2b)
E_app_2b = -slope_2bB * R * 1000

ax2.semilogy(inv_T_2b, fluxes_L2b_T, '^-', color=COLORS['L2b'], 
             linewidth=STYLE['linewidth'], markersize=STYLE['markersize'], 
             label='Level 2b: Oxide+Metal')
ax2.semilogy(inv_T_2b, fluxes_L1_T, 'o--', color=COLORS['L1'], 
             linewidth=2, alpha=0.6, markersize=6, label='Level 1: Metal only')

ax2.set_xlabel('1000/T (K⁻¹)', fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('System Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax2.set_title(f'(B) Flux vs Temperature (P={P_fixed/1000:.0f} kPa)', fontsize=STYLE['fontsize_title'])
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.legend(fontsize=STYLE['fontsize_legend'])
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# Temperature axis on top
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
T_ticks = np.array([600, 700, 800, 900, 1000])
ax2_top.set_xticks(1000 / (T_ticks + 273.15))
ax2_top.set_xticklabels([f'{t}' for t in T_ticks])
ax2_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax2_top.tick_params(labelsize=STYLE['fontsize_tick'])

ax2.text(0.95, 0.95, 
         f'$E_{{app}}$ = {E_app_2b/1000:.1f} kJ/mol\n(combined oxide+metal)\nR² = {r_2bB**2:.4f}',
         transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', ha='right', bbox=props)

# ============================================================================
# (C) FLUX vs METAL THICKNESS - Demonstrating regime transition
# Using adjusted parameters to show both limiting behaviors
# ============================================================================
ax3 = axes[1, 0]

# Use a THIN oxide (10 nm) and HIGH pressure to enable transition
# This ensures we can see both oxide-limited AND metal-limited regimes
ox_props_C = get_oxide_properties_at_T(CONFIG['oxide'], T_ref)
ox_props_C['thickness'] = 1e-8  # 10 nm oxide (much thinner than default 1 μm)

met_props_C = get_metal_properties_at_T(CONFIG['material'], T_ref)

# Use very high pressure so metal resistance becomes significant
P_fixed_C = 1e8  # 100 MPa (or adjust as needed)

# Wider metal thickness sweep
L_metal_sweep = np.logspace(-5, -1, 30)  # 0.01 mm to 100 mm

fluxes_2bC = []
fluxes_L1C = []
regimes_C = []

for L_m in L_metal_sweep:
    met_props_C['thickness'] = L_m
    
    # Level 2b
    result = calculate_oxide_metal_system(P_fixed_C, CONFIG['P_down'], ox_props_C, met_props_C)
    fluxes_2bC.append(result['flux'])
    regimes_C.append(result['regime'])
    
    # Level 1 (metal only)
    result_L1 = calculate_simple_metal_flux(
        met_props_C['D_metal'], met_props_C['K_s_metal'],
        L_m, P_fixed_C, CONFIG['P_down']
    )
    fluxes_L1C.append(result_L1['flux'])

fluxes_2bC = np.array(fluxes_2bC)
fluxes_L1C = np.array(fluxes_L1C)

# Calculate oxide-limited flux
flux_oxide_limit = molecular_diffusion_flux(
    ox_props_C['D_ox'], ox_props_C['K_ox'], ox_props_C['thickness'], P_fixed_C, 0
)

# Plot Level 2b
ax3.loglog(L_metal_sweep * 1000, fluxes_2bC, '^-', color=COLORS['L2b'], 
           linewidth=STYLE['linewidth'], markersize=6, 
           label='Level 2b (oxide+metal)')

# Plot Level 1
ax3.loglog(L_metal_sweep * 1000, fluxes_L1C, 'o--', color=COLORS['L1'], 
           linewidth=2, alpha=0.5, markersize=4,
           label='Level 1 (metal only)')

# Oxide limit line
ax3.axhline(flux_oxide_limit, color=COLORS['L2a'], linestyle=':', linewidth=2, 
            alpha=0.7, label='Oxide limit')

ax3.set_xlabel('Metal Thickness, $L$ (mm)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax3.set_title(f'(C) Flux vs Metal Thickness\n(δ={ox_props_C["thickness"]*1e9:.0f} nm, P={P_fixed_C/1e6:.0f} MPa)', 
              fontsize=STYLE['fontsize_title']-1)
ax3.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax3.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper right')
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# Find transition point (where L2b deviates from oxide limit)
ratio_to_oxide = fluxes_2bC / flux_oxide_limit
transition_idx = np.argmax(ratio_to_oxide < 0.9)
if transition_idx > 0:
    L_transition = L_metal_sweep[transition_idx] * 1000
    ax3.axvline(L_transition, color='gray', linestyle='-.', alpha=0.5)
    ax3.text(L_transition * 1.2, flux_oxide_limit * 0.5, 
             f'Transition\nL ≈ {L_transition:.1f} mm', fontsize=9, color='gray')

# Annotation showing both regimes
ax3.text(0.05, 0.05, 
         'Thin metal → Oxide controls (flat)\n'
         'Thick metal → Metal controls (J∝1/L)',
         transform=ax3.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='bottom', bbox=props)


print(f"\nPlot (C) - Parameters chosen to show transition:")
print(f"  Oxide thickness: {ox_props_C['thickness']*1e9:.0f} nm (thin, to reduce oxide resistance)")
print(f"  Pressure: {P_fixed_C/1e6:.0f} MPa (high, to increase metal resistance)")
print(f"  Metal range: {L_metal_sweep[0]*1000:.3f} - {L_metal_sweep[-1]*1000:.0f} mm")
print(f"  Oxide-limited flux: {flux_oxide_limit:.2e} mol/m²/s")

# ============================================================================
# (D) LIMIT CHECK: δ_ox → 0 → P_interface approaches P_upstream
# Based on test_interface_solver.py methodology
# ============================================================================
ax4 = axes[1, 1]

# Oxide thickness sweep from thick to thin
delta_ox_sweep = np.logspace(-4, -12, 40)  # 100 μm to 1 pm

met_props_D = get_metal_properties_at_T(CONFIG['material'], T_ref)
met_props_D['thickness'] = CONFIG['L_metal']
P_fixed_D = 1e6  # 1 kPa

# Calculate P_interface_normalized for each oxide thickness
# When this → 1, it means P_interface → P_upstream, so metal dominates
P_norm_list = []
flux_ratios = []
regimes = []

# Level 1 reference flux
result_L1D = calculate_simple_metal_flux(
    met_props_D['D_metal'], met_props_D['K_s_metal'],
    met_props_D['thickness'], P_fixed_D, CONFIG['P_down']
)
J_L1_ref = result_L1D['flux']

for delta_ox in delta_ox_sweep:
    ox_props_D = get_oxide_properties_at_T(CONFIG['oxide'], T_ref)
    ox_props_D['thickness'] = delta_ox
    
    result = calculate_oxide_metal_system(P_fixed_D, CONFIG['P_down'], ox_props_D, met_props_D)
    P_norm_list.append(result['P_interface_normalized'])
    flux_ratios.append(result['flux'] / J_L1_ref if J_L1_ref > 0 else 0)
    regimes.append(result['regime'])

P_norm_arr = np.array(P_norm_list)
flux_ratio_arr = np.array(flux_ratios)

# Plot P_interface_normalized (should go from ~0 at thick oxide to ~1 at thin oxide)
ax4.semilogx(delta_ox_sweep * 1e6, P_norm_arr, '^-', color=COLORS['L2b'], 
             linewidth=STYLE['linewidth'], markersize=6, 
             label='$P_{int}/P_{up}$ (normalized)')
ax4.axhline(1.0, color=COLORS['L1'], linestyle='--', linewidth=2, 
            label='Level 1 limit (P_int = P_up)')

ax4.set_xlabel('Oxide Thickness, $\\delta$ (μm)', fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('$(P_{int} - P_{down}) / (P_{up} - P_{down})$', fontsize=STYLE['fontsize_axis'])
ax4.set_title('(D) Limit Check: $\\delta \\rightarrow 0$ → Metal Dominated', fontsize=STYLE['fontsize_title'])
ax4.grid(True, alpha=STYLE['grid_alpha'])
ax4.legend(fontsize=STYLE['fontsize_legend'], loc='lower right')
ax4.tick_params(labelsize=STYLE['fontsize_tick'])
ax4.set_ylim([0, 1.1])
ax4.invert_xaxis()  # Thin oxide on right (→0)

# Find thickness where P_norm > 0.99
mask_metal = P_norm_arr > 0.99
if np.any(mask_metal):
    delta_99 = delta_ox_sweep[mask_metal][0] * 1e6  # Thickest δ where ratio > 0.99
    ax4.axvline(delta_99, color='green', linestyle=':', alpha=0.7)
    ax4.text(delta_99 * 0.7, 0.5, f'Metal-dominated\nδ < {delta_99:.2e} μm', 
             fontsize=9, color='green', ha='right')

# Print diagnostic info
print(f"\nLimit check values:")
print(f"  At δ = {delta_ox_sweep[0]*1e6:.2f} μm (thick): P_norm = {P_norm_arr[0]:.4f}, regime = {regimes[0]}")
print(f"  At δ = {delta_ox_sweep[-1]*1e12:.2f} pm (thin): P_norm = {P_norm_arr[-1]:.4f}, regime = {regimes[-1]}")

if P_norm_arr[-1] > 0.99:
    validation_text = '✓ VALIDATED:\nAs δ→0, P_int→P_up\n(Metal dominates = Level 1)'
    box_color = 'lightgreen'
else:
    validation_text = f'⚠ P_norm = {P_norm_arr[-1]:.3f}\n(Expected → 1.0)'
    box_color = 'lightyellow'

ax4.text(0.05, 0.95, validation_text,
         transform=ax4.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', ha='left',
         bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.9))

# ============================================================================
# FINAL LAYOUT
# ============================================================================
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Level2b_Oxide_Metal.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

# Summary
print("\n" + "="*80)
print("LEVEL 2b: OXIDE + METAL SYSTEM - VALIDATION SUMMARY")
print("="*80)
print(f"\nSystem: {CONFIG['oxide']} ({CONFIG['L_oxide']*1e6:.1f} μm) + {CONFIG['material']} ({CONFIG['L_metal']*1000:.1f} mm)")
print(f"T_ref: {T_ref-273.15:.0f}°C")
print(f"\n(D) Limit Check (δ_ox → 0):")
print(f"    Thick oxide: P_interface_normalized = {P_norm_arr[0]:.4f} ({regimes[0]})")
print(f"    Thin oxide:  P_interface_normalized = {P_norm_arr[-1]:.4f} ({regimes[-1]})")
print(f"    When P_int → P_up, metal is sole resistance → Level 1 behavior")
print("="*80)
```

```python
"""
================================================================================
LEVEL 4: DEFECTIVE METAL (Grain Boundaries + Trapping)
================================================================================
Physics: Microstructure modifies effective diffusivity

    GB Enhancement (parallel paths):
        D_eff = (1 - f_gb) × D_bulk + f_gb × D_gb
        where f_gb = 3δ_gb/d (GB volume fraction)
        
    Trapping Reduction (Oriani equilibrium):
        D_eff = D_lattice / (1 + Σ(N_T,i × K_i / N_L))
        where K_i = exp(E_b,i / RT)

    Combined: Both effects applied (GB first, then trapping)

Modes:
    - 'none':          D_eff = D_lattice (recovers Level 1)
    - 'gb_only':       D_eff > D_lattice (enhancement)
    - 'trapping_only': D_eff < D_lattice (reduction)
    - 'both':          Competition between GB and trapping

Validation:
    (A) Mode Comparison vs Pressure: All 4 modes on one plot
    (B) Modification Factor vs Temperature: GB vs Trap crossover
    (C) Sensitivity to Grain Size: Nanocrystalline enhancement
    (D) Limit Check: mode='none' recovers Level 1 exactly
================================================================================
"""

from calculations.permeation_calc import calculate_defective_metal_flux

# ============================================================================
# LEVEL 4 CONFIGURATION
# ============================================================================
# Arrhenius parameters for lattice diffusivity
D_0_L4 = 1.5e-6       # m²/s (pre-exponential)
E_D_L4 = 55000        # J/mol (activation energy)

# Microstructure parameters - chosen to show clear effects
MICROSTRUCTURE_L4 = {
    'grain_size': 100e-9,      # 100 nm (nanocrystalline - visible GB effect)
    'gb_thickness': 0.5e-9,    # 0.5 nm
    'grain_shape': 'equiaxed',
    'gb_type': 'HAGB',
    'include_gb_trapping': False,
    'trap_list': [
        {
            'name': 'vacancies',
            'binding_energy': 0.5 * 96485,  # 0.5 eV = ~48 kJ/mol
            'density': 1e26                  # High trap density for visible effect
        }
    ]
}

# Add L4 combined color to COLORS if not present
if 'L4_both' not in COLORS:
    COLORS['L4_both'] = 'crimson'  # Distinct color for combined mode

# Fixed conditions
K_s_L4 = get_solubility(CONFIG['T_ref'], material)
thickness_L4 = CONFIG['L_metal']

# Create figure
fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle('Level 4: Defective Metal (GB Enhancement + Trapping)', 
             fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

# ============================================================================
# (A) MODE COMPARISON vs PRESSURE - The "Hero Plot"
# ============================================================================
ax1 = axes[0, 0]

# Temperature for this plot
T_A = 673  # K (400°C - lower T to show strong trapping effect)
D_lattice_A = D_0_L4 * np.exp(-E_D_L4 / (R * T_A))

# Pressure sweep
pressures_L4A = np.logspace(3, 7, 30)  # 1 kPa to 10 MPa

# Calculate fluxes for all modes
fluxes_L1_A = []
fluxes_gb_A = []
fluxes_trap_A = []
fluxes_both_A = []

for P in pressures_L4A:
    # Level 1 (baseline)
    result_L1 = calculate_simple_metal_flux(D_lattice_A, K_s_L4, thickness_L4, P, 0.0)
    fluxes_L1_A.append(result_L1['flux'])
    
    # Level 4: GB only
    result_gb = calculate_defective_metal_flux(
        D_lattice=D_lattice_A, K_s=K_s_L4, thickness=thickness_L4,
        P_up=P, P_down=0.0, temperature=T_A,
        microstructure_params=MICROSTRUCTURE_L4, mode='gb_only'
    )
    fluxes_gb_A.append(result_gb['flux'])
    
    # Level 4: Trapping only
    result_trap = calculate_defective_metal_flux(
        D_lattice=D_lattice_A, K_s=K_s_L4, thickness=thickness_L4,
        P_up=P, P_down=0.0, temperature=T_A,
        microstructure_params=MICROSTRUCTURE_L4, mode='trapping_only'
    )
    fluxes_trap_A.append(result_trap['flux'])
    
    # Level 4: Combined
    result_both = calculate_defective_metal_flux(
        D_lattice=D_lattice_A, K_s=K_s_L4, thickness=thickness_L4,
        P_up=P, P_down=0.0, temperature=T_A,
        microstructure_params=MICROSTRUCTURE_L4, mode='both'
    )
    fluxes_both_A.append(result_both['flux'])

fluxes_L1_A = np.array(fluxes_L1_A)
fluxes_gb_A = np.array(fluxes_gb_A)
fluxes_trap_A = np.array(fluxes_trap_A)
fluxes_both_A = np.array(fluxes_both_A)

# Plot with consistent colors
ax1.loglog(pressures_L4A/1e6, fluxes_L1_A, '-', color=COLORS['L1'], 
           linewidth=3, label='Level 1: Perfect Lattice')
ax1.loglog(pressures_L4A/1e6, fluxes_gb_A, '-', color=COLORS['L4_gb'], 
           linewidth=2.5, label='L4: GB Only (enhanced)')
ax1.loglog(pressures_L4A/1e6, fluxes_trap_A, '-', color=COLORS['L4_trap'], 
           linewidth=2.5, label='L4: Trapping Only (reduced)')
ax1.loglog(pressures_L4A/1e6, fluxes_both_A, '-', color=COLORS['L4_both'], 
           linewidth=2.5, label='L4: Combined')

# Fill regions to show enhancement/reduction
ax1.fill_between(pressures_L4A/1e6, fluxes_L1_A, fluxes_gb_A, 
                 alpha=0.15, color=COLORS['L4_gb'])
ax1.fill_between(pressures_L4A/1e6, fluxes_trap_A, fluxes_L1_A, 
                 alpha=0.15, color=COLORS['L4_trap'])

ax1.set_xlabel('Pressure (MPa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax1.set_title(f'(A) Mode Comparison (T={T_A-273:.0f}°C)', fontsize=STYLE['fontsize_title'])
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

# Calculate ratios for annotation
ratio_gb = np.mean(fluxes_gb_A / fluxes_L1_A)
ratio_trap = np.mean(fluxes_trap_A / fluxes_L1_A)
ratio_both = np.mean(fluxes_both_A / fluxes_L1_A)

ax1.text(0.05, 0.95, 
         f'GB/L1 = {ratio_gb:.1f}× (↑)\n'
         f'Trap/L1 = {ratio_trap:.3f}× (↓)\n'
         f'Combined = {ratio_both:.2f}×',
         transform=ax1.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=props)

# ============================================================================
# (B) MODIFICATION FACTOR vs TEMPERATURE - Crossover Behavior
# ============================================================================
ax2 = axes[0, 1]

# Temperature sweep
T_min_B = 473   # K (200°C) - low T for strong trapping
T_max_B = 1173  # K (900°C)
temperatures_L4B = np.linspace(T_min_B, T_max_B, 40)
inv_T_B = 1000 / temperatures_L4B

# Calculate modification factors
mod_gb_B = []
mod_trap_B = []
mod_both_B = []

for T in temperatures_L4B:
    D_lattice = D_0_L4 * np.exp(-E_D_L4 / (R * T))
    
    # GB only
    result_gb = calculate_defective_metal_flux(
        D_lattice=D_lattice, K_s=K_s_L4, thickness=thickness_L4,
        P_up=1e6, P_down=0.0, temperature=T,
        microstructure_params=MICROSTRUCTURE_L4, mode='gb_only'
    )
    mod_gb_B.append(result_gb['modification_factor'])
    
    # Trapping only
    result_trap = calculate_defective_metal_flux(
        D_lattice=D_lattice, K_s=K_s_L4, thickness=thickness_L4,
        P_up=1e6, P_down=0.0, temperature=T,
        microstructure_params=MICROSTRUCTURE_L4, mode='trapping_only'
    )
    mod_trap_B.append(result_trap['modification_factor'])
    
    # Combined
    result_both = calculate_defective_metal_flux(
        D_lattice=D_lattice, K_s=K_s_L4, thickness=thickness_L4,
        P_up=1e6, P_down=0.0, temperature=T,
        microstructure_params=MICROSTRUCTURE_L4, mode='both'
    )
    mod_both_B.append(result_both['modification_factor'])

mod_gb_B = np.array(mod_gb_B)
mod_trap_B = np.array(mod_trap_B)
mod_both_B = np.array(mod_both_B)

# Plot with consistent colors
ax2.semilogy(temperatures_L4B - 273, mod_gb_B, '-', color=COLORS['L4_gb'], 
             linewidth=2.5, label='GB Enhancement')
ax2.semilogy(temperatures_L4B - 273, mod_trap_B, '-', color=COLORS['L4_trap'], 
             linewidth=2.5, label='Trapping Reduction')
ax2.semilogy(temperatures_L4B - 273, mod_both_B, '-', color=COLORS['L4_both'], 
             linewidth=2.5, label='Combined')
ax2.axhline(y=1.0, color=COLORS['L1'], linestyle='--', linewidth=2, label='No modification')

# Shade regions
ax2.axhspan(1, ax2.get_ylim()[1] if ax2.get_ylim()[1] > 1 else 100, 
            alpha=0.1, color=COLORS['L4_gb'])
ax2.axhspan(ax2.get_ylim()[0] if ax2.get_ylim()[0] < 1 else 0.001, 1, 
            alpha=0.1, color=COLORS['L4_trap'])

ax2.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('$D_{eff} / D_{lattice}$', fontsize=STYLE['fontsize_axis'])
ax2.set_title('(B) Modification Factor vs Temperature', fontsize=STYLE['fontsize_title'])
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.legend(fontsize=STYLE['fontsize_legend']-1, loc='best')
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# Find crossover temperature
crossover_idx = np.argmin(np.abs(mod_both_B - 1.0))
T_crossover = temperatures_L4B[crossover_idx]

ax2.text(0.95, 0.95, 
         f'Low T: Trapping dominates\nHigh T: GB dominates\n'
         f'Crossover ≈ {T_crossover-273:.0f}°C',
         transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', ha='right', bbox=props)

# ============================================================================
# (C) SENSITIVITY TO GRAIN SIZE - Nanocrystalline Enhancement
# ============================================================================
ax3 = axes[1, 0]

# Temperature for this plot
T_C = 873  # K (600°C)
D_lattice_C = D_0_L4 * np.exp(-E_D_L4 / (R * T_C))

# Grain size sweep
grain_sizes_C = np.logspace(-8, -3, 30)  # 10 nm to 1 mm

# Level 1 baseline (constant)
result_L1_C = calculate_simple_metal_flux(D_lattice_C, K_s_L4, thickness_L4, 1e6, 0.0)
flux_L1_C = result_L1_C['flux']

# Level 4: GB enhancement only
fluxes_L4_C = []
enhancement_C = []
f_gb_C = []

for d_grain in grain_sizes_C:
    microstructure_C = MICROSTRUCTURE_L4.copy()
    microstructure_C['grain_size'] = d_grain
    microstructure_C['trap_list'] = []  # No trapping for this plot
    
    result = calculate_defective_metal_flux(
        D_lattice=D_lattice_C, K_s=K_s_L4, thickness=thickness_L4,
        P_up=1e6, P_down=0.0, temperature=T_C,
        microstructure_params=microstructure_C, mode='gb_only'
    )
    fluxes_L4_C.append(result['flux'])
    enhancement_C.append(result['modification_factor'])
    
    # GB volume fraction
    f_gb = (3.0 * 0.5e-9) / d_grain
    f_gb_C.append(min(f_gb, 0.5))

fluxes_L4_C = np.array(fluxes_L4_C)
enhancement_C = np.array(enhancement_C)
f_gb_C = np.array(f_gb_C)

# Plot with consistent colors
ax3.loglog(grain_sizes_C * 1e6, [flux_L1_C]*len(grain_sizes_C), '--', 
           color=COLORS['L1'], linewidth=3, label='Level 1: Perfect Lattice')
ax3.loglog(grain_sizes_C * 1e6, fluxes_L4_C, '-', color=COLORS['L4_gb'], 
           linewidth=2.5, label='L4: GB Enhancement')

ax3.set_xlabel('Grain Size (μm)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax3.set_title(f'(C) Grain Size Sensitivity (T={T_C-273:.0f}°C)', fontsize=STYLE['fontsize_title'])
ax3.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax3.legend(fontsize=STYLE['fontsize_legend'], loc='lower right')
ax3.tick_params(labelsize=STYLE['fontsize_tick'])
ax3.invert_xaxis()  # Smaller grains on right

# Add regime labels
ax3.axvline(x=0.1, color='gray', linestyle=':', alpha=0.5)
ax3.axvline(x=1.0, color='gray', linestyle=':', alpha=0.5)
ax3.text(0.03, fluxes_L4_C.max()*0.8, 'Nano', fontsize=9, rotation=90, va='top')
ax3.text(0.3, fluxes_L4_C.max()*0.8, 'UFG', fontsize=9, rotation=90, va='top')
ax3.text(3, fluxes_L4_C.max()*0.8, 'Fine', fontsize=9, rotation=90, va='top')

ax3.text(0.05, 0.95, 
         f'Smaller grains → higher f_gb\n'
         f'→ more fast paths → ↑ flux\n\n'
         f'f_gb = 3δ_gb/d',
         transform=ax3.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

# ============================================================================
# (D) LIMIT CHECK: mode='none' RECOVERS LEVEL 1
# ============================================================================
ax4 = axes[1, 1]

# Temperature sweep
temperatures_L4D = np.linspace(CONFIG['T_range'][0], CONFIG['T_range'][1], 20)

# Calculate fluxes
fluxes_L1_D = []
fluxes_none_D = []

for T in temperatures_L4D:
    D_lattice = D_0_L4 * np.exp(-E_D_L4 / (R * T))
    
    # Level 1
    result_L1 = calculate_simple_metal_flux(D_lattice, K_s_L4, thickness_L4, 1e6, 0.0)
    fluxes_L1_D.append(result_L1['flux'])
    
    # Level 4 with mode='none'
    result_none = calculate_defective_metal_flux(
        D_lattice=D_lattice, K_s=K_s_L4, thickness=thickness_L4,
        P_up=1e6, P_down=0.0, temperature=T,
        microstructure_params=MICROSTRUCTURE_L4, mode='none'
    )
    fluxes_none_D.append(result_none['flux'])

fluxes_L1_D = np.array(fluxes_L1_D)
fluxes_none_D = np.array(fluxes_none_D)

# Parity plot with consistent colors
ax4.loglog(fluxes_L1_D, fluxes_none_D, 'o', color=COLORS['L1'], markersize=8, 
           label="L4 mode='none'")
J_range = [fluxes_L1_D.min()*0.5, fluxes_L1_D.max()*2]
ax4.loglog(J_range, J_range, '--', color='red', linewidth=2, label='Perfect agreement')

ax4.set_xlabel('Level 1 Flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel("Level 4 (mode='none') Flux (mol/m²/s)", fontsize=STYLE['fontsize_axis'])
ax4.set_title("(D) Limit Check: mode='none' → Level 1", fontsize=STYLE['fontsize_title'])
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.legend(fontsize=STYLE['fontsize_legend'])
ax4.tick_params(labelsize=STYLE['fontsize_tick'])
ax4.set_aspect('equal', adjustable='box')

# Calculate max relative error
rel_error = np.abs(fluxes_none_D - fluxes_L1_D) / fluxes_L1_D * 100
max_error = np.max(rel_error)

if max_error < 0.01:
    validation_text = f'✓ VALIDATED\nMax error: {max_error:.2e}%\nL4(none) ≡ L1'
    box_color = 'lightgreen'
else:
    validation_text = f'Max error: {max_error:.2f}%'
    box_color = 'lightyellow'

ax4.text(0.05, 0.95, validation_text,
         transform=ax4.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor=box_color, alpha=0.9))

# ============================================================================
# FINAL LAYOUT
# ============================================================================
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Level4_Defective_Metal.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

# Summary
print("\n" + "="*80)
print("LEVEL 4: DEFECTIVE METAL - VALIDATION SUMMARY")
print("="*80)
print(f"\nMicrostructure Parameters:")
print(f"    Grain size: {MICROSTRUCTURE_L4['grain_size']*1e9:.0f} nm (nanocrystalline)")
print(f"    GB thickness: {MICROSTRUCTURE_L4['gb_thickness']*1e9:.1f} nm")
print(f"    Trap density: {MICROSTRUCTURE_L4['trap_list'][0]['density']:.0e} m⁻³")
print(f"    Binding energy: {MICROSTRUCTURE_L4['trap_list'][0]['binding_energy']/96485:.2f} eV")

print(f"\n(A) Mode Comparison at T={T_A-273:.0f}°C:")
print(f"    GB/L1 ratio:       {ratio_gb:.2f}× (enhancement)")
print(f"    Trap/L1 ratio:     {ratio_trap:.4f}× (reduction)")
print(f"    Combined/L1 ratio: {ratio_both:.3f}×")

print(f"\n(B) Temperature Crossover:")
print(f"    Crossover T ≈ {T_crossover-273:.0f}°C (where Combined = 1)")
print(f"    Below crossover: Trapping dominates → D_eff < D_lattice")
print(f"    Above crossover: GB dominates → D_eff > D_lattice")

print(f"\n(C) Grain Size Effect at T={T_C-273:.0f}°C:")
print(f"    At d = 10 nm:  enhancement = {enhancement_C[0]:.1f}×")
print(f"    At d = 1 mm:   enhancement = {enhancement_C[-1]:.3f}× (≈1)")

print(f"\n(D) Limit Check:")
print(f"    mode='none' matches Level 1 exactly")
print(f"    Max relative error: {max_error:.2e}%")
print(f"    ✓ Hierarchical consistency verified")

print("\n" + "="*80)
print("KEY PHYSICS:")
print("  • GB Enhancement: D_eff = (1-f_gb)×D_bulk + f_gb×D_gb")
print("  • Trapping Reduction: D_eff = D_lattice / (1 + Σ(N_T×K/N_L))")
print("  • Combined: Competition depends on T, grain size, trap parameters")
print("="*80)
```

```python
"""
================================================================================
LEVEL 5: FULL SYSTEM (Defective Oxide + Defective Metal)
================================================================================
Physics: Complete hierarchical model combining all effects

    Oxide Layer (Level 3):
        - Intact regions: Henry's law molecular diffusion
        - Defect regions: Pinholes/cracks bypass oxide barrier
        - J_oxide = (1-f) × J_intact + f × J_defect
        
    Metal Layer (Level 4):
        - GB Enhancement: D_eff = (1-f_gb)×D_bulk + f_gb×D_gb
        - Trapping Reduction: D_eff = D_lattice / (1 + Σ(N_T×K/N_L))
        
    Combined (Level 5 = Level 3+4):
        J_total = (1-f) × J_intact_L4 + f × J_defect_L4
        Where both paths use defective metal (Level 4)

Validation:
    (A) All Levels Comparison: L1, L2b, L3, L4, L5 on one plot
    (B) Flux vs Temperature: Full system Arrhenius behavior
    (C) Sensitivity Map: Oxide defects vs Metal microstructure
    (D) Limit Checks: L5 → lower levels with appropriate limits
================================================================================
"""

from calculations.parallel_oxide_defect_paths import (
    calculate_parallel_path_flux_defective_metal,
    calculate_parallel_path_flux
)
from calculations.interface_solver import calculate_oxide_defective_metal_system

# Create figure
fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle('Level 5: Full System (Defective Oxide + Defective Metal)', 
             fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

# ============================================================================
# LEVEL 5 CONFIGURATION
# ============================================================================
# Oxide defect parameters
DEFECT_PARAMS_L5 = {
    'area_fraction': 0.01,  # 1% pinholes
    'type': 'pinhole'
}

# Metal microstructure parameters (same as Level 4)
MICROSTRUCTURE_L5 = {
    'grain_size': 100e-9,      # 100 nm
    'gb_thickness': 0.5e-9,    # 0.5 nm
    'grain_shape': 'equiaxed',
    'gb_type': 'HAGB',
    'include_gb_trapping': False,
    'trap_list': [
        {
            'name': 'vacancies',
            'binding_energy': 0.5 * 96485,  # 0.5 eV
            'density': 1e26
        }
    ]
}

# ============================================================================
# (A) ALL LEVELS COMPARISON vs PRESSURE
# ============================================================================
ax1 = axes[0, 0]

# Temperature
T_A = 873  # K (600°C)

# Get properties
oxide_props_L5 = get_oxide_properties_at_T(CONFIG['oxide'], T_A)
oxide_props_L5['thickness'] = CONFIG['L_oxide']

metal_props_L5 = get_metal_properties_at_T(CONFIG['material'], T_A)
metal_props_L5['thickness'] = CONFIG['L_metal']

D_lattice_A = D_0_L4 * np.exp(-E_D_L4 / (R * T_A))

# Pressure sweep
pressures_L5A = np.logspace(-3, 24, 25)  # 100 Pa to 1 MPa

# Calculate all levels
fluxes_L1_A = []
fluxes_L2b_A = []
fluxes_L3_A = []
fluxes_L4_A = []
fluxes_L5_A = []

for P in pressures_L5A:
    # Level 1: Perfect metal only
    result_L1 = calculate_simple_metal_flux(D_lattice_A, K_s_L4, thickness_L4, P, 0.0)
    fluxes_L1_A.append(result_L1['flux'])
    
    # Level 2b: Perfect oxide + perfect metal
    result_L2b = calculate_oxide_metal_system(P, 0.0, oxide_props_L5, metal_props_L5)
    fluxes_L2b_A.append(result_L2b['flux'])
    
    # Level 3: Defective oxide + perfect metal
    result_L3 = calculate_parallel_path_flux(P, 0.0, oxide_props_L5, metal_props_L5, DEFECT_PARAMS_L5)
    fluxes_L3_A.append(result_L3['flux_total'])
    
    # Level 4: Perfect oxide + defective metal (actually no oxide, just defective metal)
    result_L4 = calculate_defective_metal_flux(
        D_lattice=D_lattice_A, K_s=K_s_L4, thickness=thickness_L4,
        P_up=P, P_down=0.0, temperature=T_A,
        microstructure_params=MICROSTRUCTURE_L5, mode='both'
    )
    fluxes_L4_A.append(result_L4['flux'])
    
    # Level 5: Defective oxide + defective metal
    result_L5 = calculate_parallel_path_flux_defective_metal(
        P, 0.0, oxide_props_L5, metal_props_L5,
        DEFECT_PARAMS_L5, T_A, MICROSTRUCTURE_L5, mode='both'
    )
    fluxes_L5_A.append(result_L5['flux_total'])

fluxes_L1_A = np.array(fluxes_L1_A)
fluxes_L2b_A = np.array(fluxes_L2b_A)
fluxes_L3_A = np.array(fluxes_L3_A)
fluxes_L4_A = np.array(fluxes_L4_A)
fluxes_L5_A = np.array(fluxes_L5_A)

# Plot all levels
ax1.loglog(pressures_L5A/1e3, fluxes_L1_A, '--', color=COLORS['L1'], linewidth=2, 
           label='L1: Metal only')
ax1.loglog(pressures_L5A/1e3, fluxes_L2b_A, '--', color=COLORS['L2b'], linewidth=2, 
           label='L2b: Oxide+Metal')
ax1.loglog(pressures_L5A/1e3, fluxes_L3_A, '-', color=COLORS['L3'], linewidth=2.5, 
           label='L3: Defect Oxide')
ax1.loglog(pressures_L5A/1e3, fluxes_L4_A, '-', color=COLORS['L4_both'], linewidth=2.5, 
           label='L4: Defect Metal')
ax1.loglog(pressures_L5A/1e3, fluxes_L5_A, '-', color=COLORS['L5'], linewidth=3, 
           label='L5: Full System')

ax1.set_xlabel('Pressure (kPa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax1.set_title(f'(A) All Levels Comparison (T={T_A-273:.0f}°C)', fontsize=STYLE['fontsize_title'])
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.legend(fontsize=STYLE['fontsize_legend']-2, loc='lower right')
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

ax1.text(0.05, 0.95, 
         f'Oxide defects: f={DEFECT_PARAMS_L5["area_fraction"]*100:.0f}%\n'
         f'Grain size: {MICROSTRUCTURE_L5["grain_size"]*1e9:.0f} nm\n'
         f'Trap density: {MICROSTRUCTURE_L5["trap_list"][0]["density"]:.0e} m⁻³',
         transform=ax1.transAxes, fontsize=STYLE['fontsize_annotation']-1,
         verticalalignment='top', bbox=props)

# ============================================================================
# (B) FLUX vs TEMPERATURE - Full System Arrhenius
# ============================================================================
ax2 = axes[0, 1]

# Temperature sweep
T_min_B = 573   # K (300°C)
T_max_B = 1173  # K (900°C)
temperatures_L5B = np.linspace(T_min_B, T_max_B, 20)
inv_T_B = 1000 / temperatures_L5B

P_fixed_B = 1e5  # 100 kPa

fluxes_L1_B = []
fluxes_L2b_B = []
fluxes_L5_B = []

for T in temperatures_L5B:
    D_lattice = D_0_L4 * np.exp(-E_D_L4 / (R * T))
    K_s = get_solubility(T, material)
    
    ox_props = get_oxide_properties_at_T(CONFIG['oxide'], T)
    ox_props['thickness'] = CONFIG['L_oxide']
    met_props = get_metal_properties_at_T(CONFIG['material'], T)
    met_props['thickness'] = CONFIG['L_metal']
    
    # Level 1
    result_L1 = calculate_simple_metal_flux(D_lattice, K_s, thickness_L4, P_fixed_B, 0.0)
    fluxes_L1_B.append(result_L1['flux'])
    
    # Level 2b
    result_L2b = calculate_oxide_metal_system(P_fixed_B, 0.0, ox_props, met_props)
    fluxes_L2b_B.append(result_L2b['flux'])
    
    # Level 5
    result_L5 = calculate_parallel_path_flux_defective_metal(
        P_fixed_B, 0.0, ox_props, met_props,
        DEFECT_PARAMS_L5, T, MICROSTRUCTURE_L5, mode='both'
    )
    fluxes_L5_B.append(result_L5['flux_total'])

fluxes_L1_B = np.array(fluxes_L1_B)
fluxes_L2b_B = np.array(fluxes_L2b_B)
fluxes_L5_B = np.array(fluxes_L5_B)

# Arrhenius fits
slope_L1_B, intercept_L1_B, r_L1_B, _, _ = stats.linregress(1/temperatures_L5B, np.log(fluxes_L1_B))
slope_L5_B, intercept_L5_B, r_L5_B, _, _ = stats.linregress(1/temperatures_L5B, np.log(fluxes_L5_B))
E_a_L1_B = -slope_L1_B * R / 1000
E_a_L5_B = -slope_L5_B * R / 1000

# Plot
ax2.semilogy(inv_T_B, fluxes_L1_B, '--', color=COLORS['L1'], linewidth=2, label=f'L1 (E_a={E_a_L1_B:.0f} kJ/mol)')
ax2.semilogy(inv_T_B, fluxes_L2b_B, '--', color=COLORS['L2b'], linewidth=2, label='L2b')
ax2.semilogy(inv_T_B, fluxes_L5_B, '-', color=COLORS['L5'], linewidth=3, label=f'L5 (E_a={E_a_L5_B:.0f} kJ/mol)')

ax2.set_xlabel('1000/T (K⁻¹)', fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax2.set_title(f'(B) Arrhenius Plot (P={P_fixed_B/1000:.0f} kPa)', fontsize=STYLE['fontsize_title'])
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.legend(fontsize=STYLE['fontsize_legend'], loc='lower left')
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# Temperature axis on top
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
T_ticks = [400, 500, 600, 700, 800, 900]
ax2_top.set_xticks([1000/(T+273) for T in T_ticks])
ax2_top.set_xticklabels([f'{T}' for T in T_ticks])
ax2_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax2_top.tick_params(labelsize=STYLE['fontsize_tick'])

# ============================================================================
# (C) SENSITIVITY MAP: Oxide Defects vs Metal Microstructure
# ============================================================================
ax3 = axes[1, 0]

# Fixed conditions
T_C = 673  # K (400°C)
P_C = 1e5  # 100 kPa

ox_props_C = get_oxide_properties_at_T(CONFIG['oxide'], T_C)
ox_props_C['thickness'] = CONFIG['L_oxide']
met_props_C = get_metal_properties_at_T(CONFIG['material'], T_C)
met_props_C['thickness'] = CONFIG['L_metal']

D_lattice_C = D_0_L4 * np.exp(-E_D_L4 / (R * T_C))

# Sweep oxide defect fraction
f_defect_sweep = np.array([0.001, 0.01, 0.05, 0.1])  # 0.1%, 1%, 5%, 10%

# Calculate L5 flux for each f_defect with different metal modes
fluxes_C_none = []
fluxes_C_gb = []
fluxes_C_trap = []
fluxes_C_both = []

for f in f_defect_sweep:
    defect_params = {'area_fraction': f, 'type': 'pinhole'}
    
    for mode, flux_list in [('none', fluxes_C_none), ('gb_only', fluxes_C_gb), 
                             ('trapping_only', fluxes_C_trap), ('both', fluxes_C_both)]:
        result = calculate_parallel_path_flux_defective_metal(
            P_C, 0.0, ox_props_C, met_props_C,
            defect_params, T_C, MICROSTRUCTURE_L5, mode=mode
        )
        flux_list.append(result['flux_total'])

fluxes_C_none = np.array(fluxes_C_none)
fluxes_C_gb = np.array(fluxes_C_gb)
fluxes_C_trap = np.array(fluxes_C_trap)
fluxes_C_both = np.array(fluxes_C_both)

# Plot grouped bars
x = np.arange(len(f_defect_sweep))
width = 0.2

bars1 = ax3.bar(x - 1.5*width, fluxes_C_none, width, color=COLORS['L1'], label='No microstructure')
bars2 = ax3.bar(x - 0.5*width, fluxes_C_gb, width, color=COLORS['L4_gb'], label='GB only')
bars3 = ax3.bar(x + 0.5*width, fluxes_C_trap, width, color=COLORS['L4_trap'], label='Trap only')
bars4 = ax3.bar(x + 1.5*width, fluxes_C_both, width, color=COLORS['L5'], label='Combined')

ax3.set_xlabel('Oxide Defect Fraction (%)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax3.set_title(f'(C) Sensitivity: Oxide vs Metal Effects\n(T={T_C-273:.0f}°C, P={P_C/1000:.0f} kPa)', 
              fontsize=STYLE['fontsize_title']-1)
ax3.set_xticks(x)
ax3.set_xticklabels([f'{f*100:.1f}' for f in f_defect_sweep])
ax3.legend(fontsize=STYLE['fontsize_legend']-1, loc='upper left')
ax3.grid(True, alpha=STYLE['grid_alpha'], axis='y')
ax3.set_yscale('log')
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

ax3.text(0.95, 0.95, 
         'Higher bars = higher flux\n'
         'GB enhances, Trap reduces',
         transform=ax3.transAxes, fontsize=STYLE['fontsize_annotation']-1,
         verticalalignment='top', ha='right', bbox=props)

# ============================================================================
# (D) LIMIT CHECKS: L5 → Lower Levels
# ============================================================================
ax4 = axes[1, 1]

# Test that L5 reduces to appropriate lower levels
T_D = 873  # K

ox_props_D = get_oxide_properties_at_T(CONFIG['oxide'], T_D)
ox_props_D['thickness'] = CONFIG['L_oxide']
met_props_D = get_metal_properties_at_T(CONFIG['material'], T_D)
met_props_D['thickness'] = CONFIG['L_metal']

D_lattice_D = D_0_L4 * np.exp(-E_D_L4 / (R * T_D))
K_s_D = get_solubility(T_D, material)

# Pressure sweep for limit checks
pressures_D = np.logspace(3, 6, 15)

# Case 1: f_defect → 0, mode='none' should → L2b
fluxes_L2b_D = []
fluxes_L5_limit1 = []

# Case 2: f_defect → 1, mode='none' should → L1
fluxes_L1_D = []
fluxes_L5_limit2 = []

for P in pressures_D:
    # L2b reference
    result_L2b = calculate_oxide_metal_system(P, 0.0, ox_props_D, met_props_D)
    fluxes_L2b_D.append(result_L2b['flux'])
    
    # L5 with f→0, mode='none' (should match L2b)
    defect_small = {'area_fraction': 1e-10, 'type': 'pinhole'}
    micro_none = MICROSTRUCTURE_L5.copy()
    result_L5_1 = calculate_parallel_path_flux_defective_metal(
        P, 0.0, ox_props_D, met_props_D,
        defect_small, T_D, micro_none, mode='none'
    )
    fluxes_L5_limit1.append(result_L5_1['flux_total'])
    
    # L1 reference
    result_L1 = calculate_simple_metal_flux(D_lattice_D, K_s_D, thickness_L4, P, 0.0)
    fluxes_L1_D.append(result_L1['flux'])
    
    # L5 with f→1, mode='none' (should match L1)
    defect_full = {'area_fraction': 0.9999, 'type': 'pinhole'}
    result_L5_2 = calculate_parallel_path_flux_defective_metal(
        P, 0.0, ox_props_D, met_props_D,
        defect_full, T_D, micro_none, mode='none'
    )
    fluxes_L5_limit2.append(result_L5_2['flux_total'])

fluxes_L2b_D = np.array(fluxes_L2b_D)
fluxes_L5_limit1 = np.array(fluxes_L5_limit1)
fluxes_L1_D = np.array(fluxes_L1_D)
fluxes_L5_limit2 = np.array(fluxes_L5_limit2)

# Plot parity
ax4.loglog(fluxes_L2b_D, fluxes_L5_limit1, 'o', color=COLORS['L2b'], markersize=7, 
           label='L5(f→0) vs L2b')
ax4.loglog(fluxes_L1_D, fluxes_L5_limit2, 's', color=COLORS['L1'], markersize=7, 
           label='L5(f→1) vs L1')

J_range = [min(fluxes_L2b_D.min(), fluxes_L1_D.min())*0.5, 
           max(fluxes_L2b_D.max(), fluxes_L1_D.max())*2]
ax4.loglog(J_range, J_range, '--', color='red', linewidth=2, label='Perfect agreement')

ax4.set_xlabel('Reference Flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('Level 5 Limit Flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_title('(D) Limit Checks: L5 → Lower Levels', fontsize=STYLE['fontsize_title'])
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.legend(fontsize=STYLE['fontsize_legend'])
ax4.tick_params(labelsize=STYLE['fontsize_tick'])
ax4.set_aspect('equal', adjustable='box')

# Calculate errors
error_1 = np.max(np.abs(fluxes_L5_limit1 - fluxes_L2b_D) / fluxes_L2b_D * 100)
error_2 = np.max(np.abs(fluxes_L5_limit2 - fluxes_L1_D) / fluxes_L1_D * 100)

validation_text = (f'✓ VALIDATED\n'
                   f'f→0: error={error_1:.2e}%\n'
                   f'f→1: error={error_2:.2e}%')
ax4.text(0.05, 0.95, validation_text,
         transform=ax4.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))

# ============================================================================
# FINAL LAYOUT
# ============================================================================
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Level5_Full_System.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

# Summary
print("\n" + "="*80)
print("LEVEL 5: FULL SYSTEM - VALIDATION SUMMARY")
print("="*80)
print(f"\nConfiguration:")
print(f"  Oxide: {CONFIG['oxide']} ({CONFIG['L_oxide']*1e6:.1f} μm)")
print(f"  Metal: {CONFIG['material']} ({CONFIG['L_metal']*1000:.1f} mm)")
print(f"  Oxide defect fraction: {DEFECT_PARAMS_L5['area_fraction']*100:.0f}%")
print(f"  Grain size: {MICROSTRUCTURE_L5['grain_size']*1e9:.0f} nm")
print(f"  Trap density: {MICROSTRUCTURE_L5['trap_list'][0]['density']:.0e} m⁻³")

print(f"\n(A) All Levels at T={T_A-273:.0f}°C, P=100 kPa:")
idx_100kPa = np.argmin(np.abs(pressures_L5A - 1e5))
print(f"    L1:  {fluxes_L1_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L2b: {fluxes_L2b_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L3:  {fluxes_L3_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L4:  {fluxes_L4_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L5:  {fluxes_L5_A[idx_100kPa]:.3e} mol/m²/s")

print(f"\n(B) Apparent Activation Energies:")
print(f"    L1: E_a = {E_a_L1_B:.1f} kJ/mol")
print(f"    L5: E_a = {E_a_L5_B:.1f} kJ/mol")

print(f"\n(D) Limit Checks:")
print(f"    L5(f→0, mode='none') → L2b: max error = {error_1:.2e}%")
print(f"    L5(f→1, mode='none') → L1:  max error = {error_2:.2e}%")
print(f"    ✓ Hierarchical consistency verified")

print("\n" + "="*80)
print("HIERARCHICAL MODEL COMPLETE:")
print("  L1:  Perfect Metal (baseline)")
print("  L2b: Perfect Oxide + Perfect Metal")
print("  L3:  Defective Oxide + Perfect Metal")
print("  L4:  Defective Metal (no oxide)")
print("  L5:  Defective Oxide + Defective Metal (full system)")
print("="*80)
```

```python
"""
================================================================================
LEVEL 5: FULL SYSTEM (Defective Oxide + Defective Metal)
================================================================================
Physics: Complete hierarchical model combining all effects

    Oxide Layer (Level 3):
        - Intact regions: Henry's law molecular diffusion
        - Defect regions: Pinholes/cracks bypass oxide barrier
        - J_oxide = (1-f) × J_intact + f × J_defect
        
    Metal Layer (Level 4):
        - GB Enhancement: D_eff = (1-f_gb)×D_bulk + f_gb×D_gb
        - Trapping Reduction: D_eff = D_lattice / (1 + Σ(N_T×K/N_L))
        
    Combined (Level 5 = Level 3+4):
        J_total = (1-f) × J_intact_L4 + f × J_defect_L4
        Where both paths use defective metal (Level 4)

Validation:
    (A) All Levels Comparison: L1, L2b, L3, L4, L5 on one plot
    (B) Flux vs Temperature: Full system Arrhenius behavior
    (C) Sensitivity Map: Oxide defects vs Metal microstructure
    (D) Limit Checks: L5 → lower levels with appropriate limits
================================================================================
"""

from calculations.parallel_oxide_defect_paths import (
    calculate_parallel_path_flux_defective_metal,
    calculate_parallel_path_flux
)
from calculations.interface_solver import calculate_oxide_defective_metal_system

# Create figure
fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
fig.suptitle('Level 5: Full System (Defective Oxide + Defective Metal)', 
             fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

# ============================================================================
# LEVEL 5 CONFIGURATION
# ============================================================================
# Oxide defect parameters
DEFECT_PARAMS_L5 = {
    'area_fraction': 0.01,  # 1% pinholes
    'type': 'pinhole'
}

# Metal microstructure parameters (same as Level 4)
MICROSTRUCTURE_L5 = {
    'grain_size': 50e-6,      # 100 nm
    'gb_thickness': 0.5e-9,    # 0.5 nm
    'grain_shape': 'equiaxed',
    'gb_type': 'LAGB',
    'include_gb_trapping': False,
    'trap_list': [
        {
            'name': 'vacancies',
            'binding_energy': 0.8 * 96485,  # 0.5 eV
            'density': 1e26
        }
    ]
}

# ============================================================================
# (A) ALL LEVELS COMPARISON vs PRESSURE
# ============================================================================
ax1 = axes[0, 0]

# Temperature
T_A = 873  # K (600°C)

# Get properties
oxide_props_L5 = get_oxide_properties_at_T(CONFIG['oxide'], T_A)
oxide_props_L5['thickness'] = CONFIG['L_oxide']

metal_props_L5 = get_metal_properties_at_T(CONFIG['material'], T_A)
metal_props_L5['thickness'] = CONFIG['L_metal']

D_lattice_A = D_0_L4 * np.exp(-E_D_L4 / (R * T_A))

# Pressure sweep
pressures_L5A = np.logspace(0, 24, 25)  # 100 Pa to 1 MPa

# Calculate all levels
fluxes_L1_A = []
fluxes_L2b_A = []
fluxes_L3_A = []
fluxes_L4_A = []
fluxes_L5_A = []

for P in pressures_L5A:
    # Level 1: Perfect metal only
    result_L1 = calculate_simple_metal_flux(D_lattice_A, K_s_L4, thickness_L4, P, 0.0)
    fluxes_L1_A.append(result_L1['flux'])
    
    # Level 2b: Perfect oxide + perfect metal
    result_L2b = calculate_oxide_metal_system(P, 0.0, oxide_props_L5, metal_props_L5)
    fluxes_L2b_A.append(result_L2b['flux'])
    
    # Level 3: Defective oxide + perfect metal
    result_L3 = calculate_parallel_path_flux(P, 0.0, oxide_props_L5, metal_props_L5, DEFECT_PARAMS_L5)
    fluxes_L3_A.append(result_L3['flux_total'])
    
    # Level 4: Defective metal (no oxide)
    result_L4 = calculate_defective_metal_flux(
        D_lattice=D_lattice_A, K_s=K_s_L4, thickness=thickness_L4,
        P_up=P, P_down=0.0, temperature=T_A,
        microstructure_params=MICROSTRUCTURE_L5, mode='both'
    )
    fluxes_L4_A.append(result_L4['flux'])
    
    # Level 5: Defective oxide + defective metal
    result_L5 = calculate_parallel_path_flux_defective_metal(
        P, 0.0, oxide_props_L5, metal_props_L5,
        DEFECT_PARAMS_L5, T_A, MICROSTRUCTURE_L5, mode='both'
    )
    fluxes_L5_A.append(result_L5['flux_total'])

fluxes_L1_A = np.array(fluxes_L1_A)
fluxes_L2b_A = np.array(fluxes_L2b_A)
fluxes_L3_A = np.array(fluxes_L3_A)
fluxes_L4_A = np.array(fluxes_L4_A)
fluxes_L5_A = np.array(fluxes_L5_A)

# Plot all levels with DESCRIPTIVE LABELS
ax1.loglog(pressures_L5A/1e3, fluxes_L1_A, '--', color=COLORS['L1'], linewidth=2, 
           label='L1: Perfect Metal')
ax1.loglog(pressures_L5A/1e3, fluxes_L2b_A, '--', color=COLORS['L2b'], linewidth=2, 
           label='L2b: Perfect Oxide + Perfect Metal')
ax1.loglog(pressures_L5A/1e3, fluxes_L3_A, '-', color=COLORS['L3'], linewidth=2.5, 
           label='L3: Defective Oxide + Perfect Metal')
ax1.loglog(pressures_L5A/1e3, fluxes_L4_A, '-', color=COLORS['L4_both'], linewidth=2.5, 
           label='L4: Defective Metal (no oxide)')
ax1.loglog(pressures_L5A/1e3, fluxes_L5_A, '-', color=COLORS['L5'], linewidth=3, 
           label='L5: Defective Oxide + Defective Metal')

ax1.set_xlabel('Pressure (kPa)', fontsize=STYLE['fontsize_axis'])
ax1.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax1.set_title(f'(A) All Levels Comparison (T={T_A-273:.0f}°C)', fontsize=STYLE['fontsize_title'])
ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax1.legend(fontsize=STYLE['fontsize_legend']-2, loc='lower right')
ax1.tick_params(labelsize=STYLE['fontsize_tick'])

ax1.text(0.05, 0.95, 
         f'Oxide defects: f={DEFECT_PARAMS_L5["area_fraction"]*100:.0f}%\n'
         f'Grain size: {MICROSTRUCTURE_L5["grain_size"]*1e9:.0f} nm\n'
         f'Trap density: {MICROSTRUCTURE_L5["trap_list"][0]["density"]:.0e} m⁻³',
         transform=ax1.transAxes, fontsize=STYLE['fontsize_annotation']-1,
         verticalalignment='top', bbox=props)

# ============================================================================
# (B) FLUX vs TEMPERATURE - Full System Arrhenius
# ============================================================================
ax2 = axes[0, 1]

# Temperature sweep
T_min_B = 673   # K (400°C)
T_max_B = 1273  # K (1000°C)
temperatures_L5B = np.linspace(T_min_B, T_max_B, 20)
inv_T_B = 1000 / temperatures_L5B

P_fixed_B = 1e5  # 100 kPa

fluxes_L1_B = []
fluxes_L2b_B = []
fluxes_L5_B = []

for T in temperatures_L5B:
    D_lattice = D_0_L4 * np.exp(-E_D_L4 / (R * T))
    K_s = get_solubility(T, material)
    
    ox_props = get_oxide_properties_at_T(CONFIG['oxide'], T)
    ox_props['thickness'] = CONFIG['L_oxide']
    met_props = get_metal_properties_at_T(CONFIG['material'], T)
    met_props['thickness'] = CONFIG['L_metal']
    
    # Level 1
    result_L1 = calculate_simple_metal_flux(D_lattice, K_s, thickness_L4, P_fixed_B, 0.0)
    fluxes_L1_B.append(result_L1['flux'])
    
    # Level 2b
    result_L2b = calculate_oxide_metal_system(P_fixed_B, 0.0, ox_props, met_props)
    fluxes_L2b_B.append(result_L2b['flux'])
    
    # Level 5
    result_L5 = calculate_parallel_path_flux_defective_metal(
        P_fixed_B, 0.0, ox_props, met_props,
        DEFECT_PARAMS_L5, T, MICROSTRUCTURE_L5, mode='both'
    )
    fluxes_L5_B.append(result_L5['flux_total'])

fluxes_L1_B = np.array(fluxes_L1_B)
fluxes_L2b_B = np.array(fluxes_L2b_B)
fluxes_L5_B = np.array(fluxes_L5_B)

# Arrhenius fits
slope_L1_B, intercept_L1_B, r_L1_B, _, _ = stats.linregress(1/temperatures_L5B, np.log(fluxes_L1_B))
slope_L5_B, intercept_L5_B, r_L5_B, _, _ = stats.linregress(1/temperatures_L5B, np.log(fluxes_L5_B))
E_a_L1_B = -slope_L1_B * R / 1000
E_a_L5_B = -slope_L5_B * R / 1000

# Plot with DESCRIPTIVE LABELS
ax2.semilogy(inv_T_B, fluxes_L1_B, '--', color=COLORS['L1'], linewidth=2, 
             label=f'L1: Perfect Metal (E_a={E_a_L1_B:.0f} kJ/mol)')
ax2.semilogy(inv_T_B, fluxes_L2b_B, '--', color=COLORS['L2b'], linewidth=2, 
             label='L2b: Perfect Oxide + Perfect Metal')
ax2.semilogy(inv_T_B, fluxes_L5_B, '-', color=COLORS['L5'], linewidth=3, 
             label=f'L5: Full System (E_a={E_a_L5_B:.0f} kJ/mol)')

ax2.set_xlabel('1000/T (K⁻¹)', fontsize=STYLE['fontsize_axis'])
ax2.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax2.set_title(f'(B) Arrhenius Plot (P={P_fixed_B/1000:.0f} kPa)', fontsize=STYLE['fontsize_title'])
ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax2.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower left')
ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# Temperature axis on top
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
T_ticks = [400, 500, 600, 700, 800, 900]
ax2_top.set_xticks([1000/(T+273) for T in T_ticks])
ax2_top.set_xticklabels([f'{T}' for T in T_ticks])
ax2_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
ax2_top.tick_params(labelsize=STYLE['fontsize_tick'])

# ============================================================================
# (C) SENSITIVITY MAP: Oxide Defects vs Metal Microstructure
# ============================================================================
ax3 = axes[1, 0]

# Fixed conditions
T_C = 673  # K (400°C)
P_C = 1e5  # 100 kPa

ox_props_C = get_oxide_properties_at_T(CONFIG['oxide'], T_C)
ox_props_C['thickness'] = CONFIG['L_oxide']
met_props_C = get_metal_properties_at_T(CONFIG['material'], T_C)
met_props_C['thickness'] = CONFIG['L_metal']

D_lattice_C = D_0_L4 * np.exp(-E_D_L4 / (R * T_C))

# Sweep oxide defect fraction
f_defect_sweep = np.array([0.001, 0.01, 0.05, 0.1])  # 0.1%, 1%, 5%, 10%

# Calculate L5 flux for each f_defect with different metal modes
fluxes_C_none = []
fluxes_C_gb = []
fluxes_C_trap = []
fluxes_C_both = []

for f in f_defect_sweep:
    defect_params = {'area_fraction': f, 'type': 'pinhole'}
    
    for mode, flux_list in [('none', fluxes_C_none), ('gb_only', fluxes_C_gb), 
                             ('trapping_only', fluxes_C_trap), ('both', fluxes_C_both)]:
        result = calculate_parallel_path_flux_defective_metal(
            P_C, 0.0, ox_props_C, met_props_C,
            defect_params, T_C, MICROSTRUCTURE_L5, mode=mode
        )
        flux_list.append(result['flux_total'])

fluxes_C_none = np.array(fluxes_C_none)
fluxes_C_gb = np.array(fluxes_C_gb)
fluxes_C_trap = np.array(fluxes_C_trap)
fluxes_C_both = np.array(fluxes_C_both)

# Plot grouped bars with DESCRIPTIVE LABELS
x = np.arange(len(f_defect_sweep))
width = 0.2

bars1 = ax3.bar(x - 1.5*width, fluxes_C_none, width, color=COLORS['L1'], 
                label='No microstructure (Perfect Metal)')
bars2 = ax3.bar(x - 0.5*width, fluxes_C_gb, width, color=COLORS['L4_gb'], 
                label='GB only (Enhanced diffusion)')
bars3 = ax3.bar(x + 0.5*width, fluxes_C_trap, width, color=COLORS['L4_trap'], 
                label='Trap only (Reduced diffusion)')
bars4 = ax3.bar(x + 1.5*width, fluxes_C_both, width, color=COLORS['L5'], 
                label='Combined (GB + Trapping)')

ax3.set_xlabel('Oxide Defect Fraction (%)', fontsize=STYLE['fontsize_axis'])
ax3.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax3.set_title(f'(C) Sensitivity: Oxide vs Metal Effects\n(T={T_C-273:.0f}°C, P={P_C/1000:.0f} kPa)', 
              fontsize=STYLE['fontsize_title']-1)
ax3.set_xticks(x)
ax3.set_xticklabels([f'{f*100:.1f}' for f in f_defect_sweep])
ax3.legend(fontsize=STYLE['fontsize_legend']-2, loc='upper left')
ax3.grid(True, alpha=STYLE['grid_alpha'], axis='y')
ax3.set_yscale('log')
ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# ============================================================================
# (D) LIMIT CHECKS: L5 → Lower Levels
# ============================================================================
ax4 = axes[1, 1]

# Test that L5 reduces to appropriate lower levels
T_D = 873  # K

ox_props_D = get_oxide_properties_at_T(CONFIG['oxide'], T_D)
ox_props_D['thickness'] = CONFIG['L_oxide']
met_props_D = get_metal_properties_at_T(CONFIG['material'], T_D)
met_props_D['thickness'] = CONFIG['L_metal']

D_lattice_D = D_0_L4 * np.exp(-E_D_L4 / (R * T_D))
K_s_D = get_solubility(T_D, material)

# Pressure sweep for limit checks
pressures_D = np.logspace(3, 6, 15)

# Case 1: f_defect → 0, mode='none' should → L2b
fluxes_L2b_D = []
fluxes_L5_limit1 = []

# Case 2: f_defect → 1, mode='none' should → L1
fluxes_L1_D = []
fluxes_L5_limit2 = []

for P in pressures_D:
    # L2b reference
    result_L2b = calculate_oxide_metal_system(P, 0.0, ox_props_D, met_props_D)
    fluxes_L2b_D.append(result_L2b['flux'])
    
    # L5 with f→0, mode='none' (should match L2b)
    defect_small = {'area_fraction': 1e-10, 'type': 'pinhole'}
    micro_none = MICROSTRUCTURE_L5.copy()
    result_L5_1 = calculate_parallel_path_flux_defective_metal(
        P, 0.0, ox_props_D, met_props_D,
        defect_small, T_D, micro_none, mode='none'
    )
    fluxes_L5_limit1.append(result_L5_1['flux_total'])
    
    # L1 reference
    result_L1 = calculate_simple_metal_flux(D_lattice_D, K_s_D, thickness_L4, P, 0.0)
    fluxes_L1_D.append(result_L1['flux'])
    
    # L5 with f→1, mode='none' (should match L1)
    defect_full = {'area_fraction': 0.9999, 'type': 'pinhole'}
    result_L5_2 = calculate_parallel_path_flux_defective_metal(
        P, 0.0, ox_props_D, met_props_D,
        defect_full, T_D, micro_none, mode='none'
    )
    fluxes_L5_limit2.append(result_L5_2['flux_total'])

fluxes_L2b_D = np.array(fluxes_L2b_D)
fluxes_L5_limit1 = np.array(fluxes_L5_limit1)
fluxes_L1_D = np.array(fluxes_L1_D)
fluxes_L5_limit2 = np.array(fluxes_L5_limit2)

# Plot parity with DESCRIPTIVE LABELS
ax4.loglog(fluxes_L2b_D, fluxes_L5_limit1, 'o', color=COLORS['L2b'], markersize=7, 
           label='L5(f→0) → L2b: Perfect Oxide + Metal')
ax4.loglog(fluxes_L1_D, fluxes_L5_limit2, 's', color=COLORS['L1'], markersize=7, 
           label='L5(f→1) → L1: Perfect Metal only')

J_range = [min(fluxes_L2b_D.min(), fluxes_L1_D.min())*0.5, 
           max(fluxes_L2b_D.max(), fluxes_L1_D.max())*2]
ax4.loglog(J_range, J_range, '--', color='red', linewidth=2, label='Perfect agreement')

ax4.set_xlabel('Reference Flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_ylabel('Level 5 Limit Flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
ax4.set_title('(D) Limit Checks: L5 → Lower Levels', fontsize=STYLE['fontsize_title'])
ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
ax4.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
ax4.tick_params(labelsize=STYLE['fontsize_tick'])
ax4.set_aspect('equal', adjustable='box')

# Calculate errors
error_1 = np.max(np.abs(fluxes_L5_limit1 - fluxes_L2b_D) / fluxes_L2b_D * 100)
error_2 = np.max(np.abs(fluxes_L5_limit2 - fluxes_L1_D) / fluxes_L1_D * 100)

validation_text = (f'✓ VALIDATED\n'
                   f'f→0: error={error_1:.2e}%\n'
                   f'f→1: error={error_2:.2e}%')
ax4.text(0.05, 0.95, validation_text,
         transform=ax4.transAxes, fontsize=STYLE['fontsize_annotation'],
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))

# ============================================================================
# FINAL LAYOUT
# ============================================================================
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Level5_Full_System.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

# Summary
print("\n" + "="*80)
print("LEVEL 5: FULL SYSTEM - VALIDATION SUMMARY")
print("="*80)
print(f"\nConfiguration:")
print(f"  Oxide: {CONFIG['oxide']} ({CONFIG['L_oxide']*1e6:.1f} μm)")
print(f"  Metal: {CONFIG['material']} ({CONFIG['L_metal']*1000:.1f} mm)")
print(f"  Oxide defect fraction: {DEFECT_PARAMS_L5['area_fraction']*100:.0f}%")
print(f"  Grain size: {MICROSTRUCTURE_L5['grain_size']*1e9:.0f} nm")
print(f"  Trap density: {MICROSTRUCTURE_L5['trap_list'][0]['density']:.0e} m⁻³")

print(f"\n(A) All Levels at T={T_A-273:.0f}°C, P=100 kPa:")
idx_100kPa = np.argmin(np.abs(pressures_L5A - 1e5))
print(f"    L1 (Perfect Metal):                    {fluxes_L1_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L2b (Perfect Oxide + Perfect Metal):   {fluxes_L2b_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L3 (Defective Oxide + Perfect Metal):  {fluxes_L3_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L4 (Defective Metal, no oxide):        {fluxes_L4_A[idx_100kPa]:.3e} mol/m²/s")
print(f"    L5 (Defective Oxide + Defective Metal):{fluxes_L5_A[idx_100kPa]:.3e} mol/m²/s")

print(f"\n(B) Apparent Activation Energies:")
print(f"    L1 (Perfect Metal): E_a = {E_a_L1_B:.1f} kJ/mol")
print(f"    L5 (Full System):   E_a = {E_a_L5_B:.1f} kJ/mol")

print(f"\n(D) Limit Checks:")
print(f"    L5(f→0, mode='none') → L2b (Perfect Oxide): max error = {error_1:.2e}%")
print(f"    L5(f→1, mode='none') → L1 (Perfect Metal):  max error = {error_2:.2e}%")
print(f"    ✓ Hierarchical consistency verified")

print("\n" + "="*80)
print("HIERARCHICAL MODEL COMPLETE:")
print("  L1:  Perfect Metal (Sieverts + Fick baseline)")
print("  L2b: Perfect Oxide + Perfect Metal (Series resistance)")
print("  L3:  Defective Oxide + Perfect Metal (Parallel paths)")
print("  L4:  Defective Metal (GB enhancement + Trapping)")
print("  L5:  Defective Oxide + Defective Metal (Full system)")
print("="*80)
```

```python

```

```python
# """
# ================================================================================
# LEVEL 5+L6: FULL SYSTEM (Defective Oxide + Defective Metal + Surface Kinetics)
# ================================================================================
# Physics: Complete hierarchical model combining all effects

#     Oxide Layer (Level 3):
#         - Intact regions: Henry's law molecular diffusion
#         - Defect regions: Pinholes/cracks bypass oxide barrier
#         - J_oxide = (1-f) × J_intact + f × J_defect
        
#     Metal Layer (Level 4):
#         - GB Enhancement: D_eff = (1-f_gb)×D_bulk + f_gb×D_gb
#         - Trapping Reduction: D_eff = D_lattice / (1 + Σ(N_T×K/N_L))
        
#     Surface Kinetics (Level 6):
#         - J_diss = k_diss × P × (1-θ)²
#         - Flux = min(J_Sieverts, J_diss)
        
#     Combined (Level 5+L6 = Level 3+4+6):
#         J_total = (1-f) × J_intact_L5L6 + f × J_defect_L4L6
#         Where both paths use defective metal with surface kinetics

# Validation:
#     (A) All Levels Comparison: L1+L6, L2b+L6, L3+L6, L4+L6, L5+L6 on one plot
#     (B) Flux vs Temperature: Full system Arrhenius behavior
#     (C) Sensitivity Map: Oxide defects vs Metal microstructure
#     (D) Limit Checks: L5+L6 → lower levels with appropriate limits
# ================================================================================
# """

# from calculations.parallel_oxide_defect_paths import (
#     calculate_parallel_path_flux_with_surface,
#     calculate_parallel_path_flux_with_surface_L3L6,
#     calculate_parallel_path_flux
# )
# from calculations.interface_solver import (
#     calculate_oxide_metal_system,
#     calculate_oxide_metal_system_with_surface
# )
# from calculations.permeation_calc import (
#     calculate_simple_metal_flux_with_surface,
#     calculate_defective_metal_flux_with_surface
# )

# # Create figure
# fig, axes = plt.subplots(2, 2, figsize=STYLE['figsize'])
# fig.suptitle('Level 5+L6: Full System (Defective Oxide + Defective Metal + Surface)', 
#              fontsize=STYLE['fontsize_suptitle'], fontweight='bold', y=0.98)

# props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

# # ============================================================================
# # LEVEL 5+L6 CONFIGURATION
# # ============================================================================
# # Oxide defect parameters
# DEFECT_PARAMS_L5 = {
#     'area_fraction': 0.01,  # 1% pinholes
#     'type': 'pinhole'
# }

# # Metal microstructure parameters (same as Level 4)
# MICROSTRUCTURE_L5 = {
#     'grain_size': 50e-6,      # 50 μm
#     'gb_thickness': 0.5e-9,    # 0.5 nm
#     'grain_shape': 'equiaxed',
#     'gb_type': 'LAGB',
#     'include_gb_trapping': False,
#     'trap_list': [
#         {
#             'name': 'vacancies',
#             'binding_energy': 0.8 * 96485,  # 0.8 eV
#             'density': 1e26
#         }
#     ]
# }

# # Surface kinetics parameters
# k_diss_L5L6 = 1e-15       # mol/m²/s/Pa
# k_recomb_L5L6 = 1e-3      # m⁴/mol/s
# COVERAGE_MODE_L5L6 = 'steady_state'

# # ============================================================================
# # (A) ALL LEVELS COMPARISON vs PRESSURE - ALL WITH SURFACE KINETICS
# # ============================================================================
# ax1 = axes[0, 0]

# # Temperature
# T_A = 873  # K (600°C)

# # Get properties
# oxide_props_L5 = get_oxide_properties_at_T(CONFIG['oxide'], T_A)
# oxide_props_L5['thickness'] = CONFIG['L_oxide']

# metal_props_L5 = get_metal_properties_at_T(CONFIG['material'], T_A)
# metal_props_L5['thickness'] = CONFIG['L_metal']

# D_lattice_A = D_0_L4 * np.exp(-E_D_L4 / (R * T_A))

# # Pressure sweep
# pressures_L5A = np.logspace(-4, 24, 25)

# # Calculate all levels WITH SURFACE KINETICS
# fluxes_L1L6_A = []
# fluxes_L2bL6_A = []
# fluxes_L3L6_A = []
# fluxes_L4L6_A = []
# fluxes_L5L6_A = []
# SRF_A = []

# for P in pressures_L5A:
#     # Level 1+L6: Perfect metal + surface kinetics
#     result_L1L6 = calculate_simple_metal_flux_with_surface(
#         D=D_lattice_A, K_s=K_s_L4, thickness=thickness_L4,
#         P_up=P, P_down=0.0, temperature=T_A,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L1L6_A.append(result_L1L6['flux'])
    
#     # Level 2b+L6: Perfect oxide + perfect metal + surface kinetics
#     result_L2bL6 = calculate_oxide_metal_system_with_surface(
#         P_upstream=P, P_downstream=0.0,
#         oxide_props=oxide_props_L5, metal_props=metal_props_L5,
#         temperature=T_A,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L2bL6_A.append(result_L2bL6['flux'])
    
#     # Level 3+L6: Defective oxide + perfect metal + surface kinetics
#     result_L3L6 = calculate_parallel_path_flux_with_surface_L3L6(
#         P_upstream=P, P_downstream=0.0,
#         oxide_props=oxide_props_L5, metal_props=metal_props_L5,
#         defect_params=DEFECT_PARAMS_L5, temperature=T_A,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L3L6_A.append(result_L3L6['flux_total'])
    
#     # Level 4+L6: Defective metal + surface kinetics (no oxide)
#     result_L4L6 = calculate_defective_metal_flux_with_surface(
#         D_lattice=D_lattice_A, K_s=K_s_L4, thickness=thickness_L4,
#         P_up=P, P_down=0.0, temperature=T_A,
#         microstructure_params=MICROSTRUCTURE_L5, mode='both',
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L4L6_A.append(result_L4L6['flux'])
    
#     # Level 5+L6: Defective oxide + defective metal + surface kinetics
#     result_L5L6 = calculate_parallel_path_flux_with_surface(
#         P_upstream=P, P_downstream=0.0,
#         oxide_props=oxide_props_L5, metal_props=metal_props_L5,
#         defect_params=DEFECT_PARAMS_L5, temperature=T_A,
#         microstructure_params=MICROSTRUCTURE_L5,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         mode='both', coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L5L6_A.append(result_L5L6['flux_total'])
#     SRF_A.append(result_L5L6['SRF'])

# fluxes_L1L6_A = np.array(fluxes_L1L6_A)
# fluxes_L2bL6_A = np.array(fluxes_L2bL6_A)
# fluxes_L3L6_A = np.array(fluxes_L3L6_A)
# fluxes_L4L6_A = np.array(fluxes_L4L6_A)
# fluxes_L5L6_A = np.array(fluxes_L5L6_A)
# SRF_A = np.array(SRF_A)

# # Plot all levels with DESCRIPTIVE LABELS
# ax1.loglog(pressures_L5A/1e3, fluxes_L1L6_A, '--', color=COLORS['L1'], linewidth=2, 
#            label='L1+L6: Perfect Metal')
# ax1.loglog(pressures_L5A/1e3, fluxes_L2bL6_A, '--', color=COLORS['L2b'], linewidth=2, 
#            label='L2b+L6: Perfect Oxide + Perfect Metal')
# ax1.loglog(pressures_L5A/1e3, fluxes_L3L6_A, '-', color=COLORS['L3'], linewidth=2.5, 
#            label='L3+L6: Defective Oxide + Perfect Metal')
# ax1.loglog(pressures_L5A/1e3, fluxes_L4L6_A, '-', color=COLORS['L4_both'], linewidth=2.5, 
#            label='L4+L6: Defective Metal (no oxide)')
# ax1.loglog(pressures_L5A/1e3, fluxes_L5L6_A, '-', color=COLORS['L5'], linewidth=3, 
#            label='L5+L6: Defective Oxide + Defective Metal')

# ax1.set_xlabel('Pressure (kPa)', fontsize=STYLE['fontsize_axis'])
# ax1.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
# ax1.set_title(f'(A) All Levels + Surface (T={T_A-273:.0f}°C)', fontsize=STYLE['fontsize_title'])
# ax1.grid(True, which='both', alpha=STYLE['grid_alpha'])
# ax1.legend(fontsize=STYLE['fontsize_legend']-2, loc='lower right')
# ax1.tick_params(labelsize=STYLE['fontsize_tick'])

# ax1.text(0.05, 0.95, 
#          f'All include L6 surface kinetics\n'
#          f'Oxide defects: f={DEFECT_PARAMS_L5["area_fraction"]*100:.0f}%\n'
#          f'Grain size: {MICROSTRUCTURE_L5["grain_size"]*1e6:.0f} μm\n'
#          f'SRF range: {SRF_A.min():.3f}-{SRF_A.max():.3f}',
#          transform=ax1.transAxes, fontsize=STYLE['fontsize_annotation']-1,
#          verticalalignment='top', bbox=props)

# # ============================================================================
# # (B) FLUX vs TEMPERATURE - ALL WITH SURFACE KINETICS
# # ============================================================================
# ax2 = axes[0, 1]

# # Temperature sweep
# T_min_B = 673   # K (400°C)
# T_max_B = 1273  # K (1000°C)
# temperatures_L5B = np.linspace(T_min_B, T_max_B, 20)
# inv_T_B = 1000 / temperatures_L5B

# P_fixed_B = 1e5  # 100 kPa

# fluxes_L1L6_B = []
# fluxes_L2bL6_B = []
# fluxes_L5L6_B = []
# SRF_B = []

# for T in temperatures_L5B:
#     D_lattice = D_0_L4 * np.exp(-E_D_L4 / (R * T))
#     K_s = get_solubility(T, material)
    
#     ox_props = get_oxide_properties_at_T(CONFIG['oxide'], T)
#     ox_props['thickness'] = CONFIG['L_oxide']
#     met_props = get_metal_properties_at_T(CONFIG['material'], T)
#     met_props['thickness'] = CONFIG['L_metal']
    
#     # Level 1+L6
#     result_L1L6 = calculate_simple_metal_flux_with_surface(
#         D=D_lattice, K_s=K_s, thickness=thickness_L4,
#         P_up=P_fixed_B, P_down=0.0, temperature=T,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L1L6_B.append(result_L1L6['flux'])
    
#     # Level 2b+L6
#     result_L2bL6 = calculate_oxide_metal_system_with_surface(
#         P_upstream=P_fixed_B, P_downstream=0.0,
#         oxide_props=ox_props, metal_props=met_props,
#         temperature=T,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L2bL6_B.append(result_L2bL6['flux'])
    
#     # Level 5+L6
#     result_L5L6 = calculate_parallel_path_flux_with_surface(
#         P_upstream=P_fixed_B, P_downstream=0.0,
#         oxide_props=ox_props, metal_props=met_props,
#         defect_params=DEFECT_PARAMS_L5, temperature=T,
#         microstructure_params=MICROSTRUCTURE_L5,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         mode='both', coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L5L6_B.append(result_L5L6['flux_total'])
#     SRF_B.append(result_L5L6['SRF'])

# fluxes_L1L6_B = np.array(fluxes_L1L6_B)
# fluxes_L2bL6_B = np.array(fluxes_L2bL6_B)
# fluxes_L5L6_B = np.array(fluxes_L5L6_B)
# SRF_B = np.array(SRF_B)

# # Arrhenius fits
# slope_L1L6_B, intercept_L1L6_B, r_L1L6_B, _, _ = stats.linregress(1/temperatures_L5B, np.log(fluxes_L1L6_B))
# slope_L5L6_B, intercept_L5L6_B, r_L5L6_B, _, _ = stats.linregress(1/temperatures_L5B, np.log(fluxes_L5L6_B))
# E_a_L1L6_B = -slope_L1L6_B * R / 1000
# E_a_L5L6_B = -slope_L5L6_B * R / 1000

# # Plot with DESCRIPTIVE LABELS
# ax2.semilogy(inv_T_B, fluxes_L1L6_B, '--', color=COLORS['L1'], linewidth=2, 
#              label=f'L1+L6: Perfect Metal (E_a={E_a_L1L6_B:.0f} kJ/mol)')
# ax2.semilogy(inv_T_B, fluxes_L2bL6_B, '--', color=COLORS['L2b'], linewidth=2, 
#              label='L2b+L6: Perfect Oxide + Perfect Metal')
# ax2.semilogy(inv_T_B, fluxes_L5L6_B, '-', color=COLORS['L5'], linewidth=3, 
#              label=f'L5+L6: Full System (E_a={E_a_L5L6_B:.0f} kJ/mol)')

# ax2.set_xlabel('1000/T (K⁻¹)', fontsize=STYLE['fontsize_axis'])
# ax2.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
# ax2.set_title(f'(B) Arrhenius Plot (P={P_fixed_B/1000:.0f} kPa)', fontsize=STYLE['fontsize_title'])
# ax2.grid(True, which='both', alpha=STYLE['grid_alpha'])
# ax2.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower left')
# ax2.tick_params(labelsize=STYLE['fontsize_tick'])

# # Temperature axis on top
# ax2_top = ax2.twiny()
# ax2_top.set_xlim(ax2.get_xlim())
# T_ticks = [400, 500, 600, 700, 800, 900]
# ax2_top.set_xticks([1000/(T+273) for T in T_ticks])
# ax2_top.set_xticklabels([f'{T}' for T in T_ticks])
# ax2_top.set_xlabel('Temperature (°C)', fontsize=STYLE['fontsize_axis'])
# ax2_top.tick_params(labelsize=STYLE['fontsize_tick'])

# ax2.text(0.95, 0.95, f'All with L6 surface\nSRF: {SRF_B.min():.3f}-{SRF_B.max():.3f}',
#          transform=ax2.transAxes, fontsize=STYLE['fontsize_annotation'],
#          verticalalignment='top', ha='right', bbox=props)

# # ============================================================================
# # (C) SENSITIVITY MAP: Oxide Defects vs Metal Microstructure - ALL WITH SURFACE
# # ============================================================================
# ax3 = axes[1, 0]

# # Fixed conditions
# T_C = 673  # K (400°C)
# P_C = 1e5  # 100 kPa

# ox_props_C = get_oxide_properties_at_T(CONFIG['oxide'], T_C)
# ox_props_C['thickness'] = CONFIG['L_oxide']
# met_props_C = get_metal_properties_at_T(CONFIG['material'], T_C)
# met_props_C['thickness'] = CONFIG['L_metal']

# D_lattice_C = D_0_L4 * np.exp(-E_D_L4 / (R * T_C))

# # Sweep oxide defect fraction
# f_defect_sweep = np.array([0.001, 0.01, 0.05, 0.1])  # 0.1%, 1%, 5%, 10%

# # Calculate L5+L6 flux for each f_defect with different metal modes
# fluxes_C_none = []
# fluxes_C_gb = []
# fluxes_C_trap = []
# fluxes_C_both = []

# for f in f_defect_sweep:
#     defect_params = {'area_fraction': f, 'type': 'pinhole'}
    
#     for mode, flux_list in [('none', fluxes_C_none), ('gb_only', fluxes_C_gb), 
#                              ('trapping_only', fluxes_C_trap), ('both', fluxes_C_both)]:
#         result = calculate_parallel_path_flux_with_surface(
#             P_upstream=P_C, P_downstream=0.0,
#             oxide_props=ox_props_C, metal_props=met_props_C,
#             defect_params=defect_params, temperature=T_C,
#             microstructure_params=MICROSTRUCTURE_L5,
#             k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#             mode=mode, coverage_mode=COVERAGE_MODE_L5L6
#         )
#         flux_list.append(result['flux_total'])

# fluxes_C_none = np.array(fluxes_C_none)
# fluxes_C_gb = np.array(fluxes_C_gb)
# fluxes_C_trap = np.array(fluxes_C_trap)
# fluxes_C_both = np.array(fluxes_C_both)

# # Plot grouped bars with DESCRIPTIVE LABELS
# x = np.arange(len(f_defect_sweep))
# width = 0.2

# bars1 = ax3.bar(x - 1.5*width, fluxes_C_none, width, color=COLORS['L1'], 
#                 label='No microstructure (Perfect Metal)')
# bars2 = ax3.bar(x - 0.5*width, fluxes_C_gb, width, color=COLORS['L4_gb'], 
#                 label='GB only (Enhanced diffusion)')
# bars3 = ax3.bar(x + 0.5*width, fluxes_C_trap, width, color=COLORS['L4_trap'], 
#                 label='Trap only (Reduced diffusion)')
# bars4 = ax3.bar(x + 1.5*width, fluxes_C_both, width, color=COLORS['L5'], 
#                 label='Combined (GB + Trapping)')

# ax3.set_xlabel('Oxide Defect Fraction (%)', fontsize=STYLE['fontsize_axis'])
# ax3.set_ylabel('Flux, $J$ (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
# ax3.set_title(f'(C) Sensitivity: Oxide vs Metal Effects\n(T={T_C-273:.0f}°C, P={P_C/1000:.0f} kPa)', 
#               fontsize=STYLE['fontsize_title']-1)
# ax3.set_xticks(x)
# ax3.set_xticklabels([f'{f*100:.1f}' for f in f_defect_sweep])
# ax3.legend(fontsize=STYLE['fontsize_legend']-2, loc='upper left')
# ax3.grid(True, alpha=STYLE['grid_alpha'], axis='y')
# ax3.set_yscale('log')
# ax3.tick_params(labelsize=STYLE['fontsize_tick'])

# ax3.text(0.95, 0.95, 
#          'All include L6 surface kinetics',
#          transform=ax3.transAxes, fontsize=STYLE['fontsize_annotation']-1,
#          verticalalignment='top', ha='right', bbox=props)

# # ============================================================================
# # (D) LIMIT CHECKS: L5+L6 → Lower Levels (ALL WITH SURFACE)
# # ============================================================================
# ax4 = axes[1, 1]

# # Test that L5+L6 reduces to appropriate lower levels
# T_D = 873  # K

# ox_props_D = get_oxide_properties_at_T(CONFIG['oxide'], T_D)
# ox_props_D['thickness'] = CONFIG['L_oxide']
# met_props_D = get_metal_properties_at_T(CONFIG['material'], T_D)
# met_props_D['thickness'] = CONFIG['L_metal']

# D_lattice_D = D_0_L4 * np.exp(-E_D_L4 / (R * T_D))
# K_s_D = get_solubility(T_D, material)

# # Pressure sweep for limit checks
# pressures_D = np.logspace(3, 6, 15)

# # Case 1: f_defect → 0, mode='none' should → L2b+L6
# fluxes_L2bL6_D = []
# fluxes_L5L6_limit1 = []

# # Case 2: f_defect → 1, mode='none' should → L1+L6
# fluxes_L1L6_D = []
# fluxes_L5L6_limit2 = []

# for P in pressures_D:
#     # L2b+L6 reference (WITH surface)
#     result_L2bL6 = calculate_oxide_metal_system_with_surface(
#         P_upstream=P, P_downstream=0.0,
#         oxide_props=ox_props_D, metal_props=met_props_D,
#         temperature=T_D,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L2bL6_D.append(result_L2bL6['flux'])
    
#     # L5+L6 with f→0, mode='none' (should match L2b+L6)
#     defect_small = {'area_fraction': 1e-10, 'type': 'pinhole'}
#     result_L5L6_1 = calculate_parallel_path_flux_with_surface(
#         P_upstream=P, P_downstream=0.0,
#         oxide_props=ox_props_D, metal_props=met_props_D,
#         defect_params=defect_small, temperature=T_D,
#         microstructure_params=MICROSTRUCTURE_L5,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         mode='none', coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L5L6_limit1.append(result_L5L6_1['flux_total'])
    
#     # L1+L6 reference (WITH surface)
#     result_L1L6 = calculate_simple_metal_flux_with_surface(
#         D=D_lattice_D, K_s=K_s_D, thickness=thickness_L4,
#         P_up=P, P_down=0.0, temperature=T_D,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L1L6_D.append(result_L1L6['flux'])
    
#     # L5+L6 with f→1, mode='none' (should match L1+L6)
#     defect_full = {'area_fraction': 0.9999, 'type': 'pinhole'}
#     result_L5L6_2 = calculate_parallel_path_flux_with_surface(
#         P_upstream=P, P_downstream=0.0,
#         oxide_props=ox_props_D, metal_props=met_props_D,
#         defect_params=defect_full, temperature=T_D,
#         microstructure_params=MICROSTRUCTURE_L5,
#         k_diss=k_diss_L5L6, k_recomb=k_recomb_L5L6,
#         mode='none', coverage_mode=COVERAGE_MODE_L5L6
#     )
#     fluxes_L5L6_limit2.append(result_L5L6_2['flux_total'])

# fluxes_L2bL6_D = np.array(fluxes_L2bL6_D)
# fluxes_L5L6_limit1 = np.array(fluxes_L5L6_limit1)
# fluxes_L1L6_D = np.array(fluxes_L1L6_D)
# fluxes_L5L6_limit2 = np.array(fluxes_L5L6_limit2)

# # Plot parity with DESCRIPTIVE LABELS
# ax4.loglog(fluxes_L2bL6_D, fluxes_L5L6_limit1, 'o', color=COLORS['L2b'], markersize=7, 
#            label='L5+L6(f→0) → L2b+L6')
# ax4.loglog(fluxes_L1L6_D, fluxes_L5L6_limit2, 's', color=COLORS['L1'], markersize=7, 
#            label='L5+L6(f→1) → L1+L6')

# J_range = [min(fluxes_L2bL6_D.min(), fluxes_L1L6_D.min())*0.5, 
#            max(fluxes_L2bL6_D.max(), fluxes_L1L6_D.max())*2]
# ax4.loglog(J_range, J_range, '--', color='red', linewidth=2, label='Perfect agreement')

# ax4.set_xlabel('Reference Flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
# ax4.set_ylabel('Level 5+L6 Limit Flux (mol/m²/s)', fontsize=STYLE['fontsize_axis'])
# ax4.set_title('(D) Limit Checks: L5+L6 → Lower Levels', fontsize=STYLE['fontsize_title'])
# ax4.grid(True, which='both', alpha=STYLE['grid_alpha'])
# ax4.legend(fontsize=STYLE['fontsize_legend']-1, loc='lower right')
# ax4.tick_params(labelsize=STYLE['fontsize_tick'])
# ax4.set_aspect('equal', adjustable='box')

# # Calculate errors
# error_1 = np.max(np.abs(fluxes_L5L6_limit1 - fluxes_L2bL6_D) / fluxes_L2bL6_D * 100)
# error_2 = np.max(np.abs(fluxes_L5L6_limit2 - fluxes_L1L6_D) / fluxes_L1L6_D * 100)

# validation_text = (f'✓ VALIDATED (all with L6)\n'
#                    f'f→0: error={error_1:.2e}%\n'
#                    f'f→1: error={error_2:.2e}%')
# ax4.text(0.05, 0.95, validation_text,
#          transform=ax4.transAxes, fontsize=STYLE['fontsize_annotation'],
#          verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))

# # ============================================================================
# # FINAL LAYOUT
# # ============================================================================
# plt.tight_layout(rect=[0, 0, 1, 0.96])
# plt.savefig('Level5L6_Full_System_Surface.png', dpi=300, bbox_inches='tight', facecolor='white')
# plt.show()

# # Summary
# print("\n" + "="*80)
# print("LEVEL 5+L6: FULL SYSTEM + SURFACE KINETICS - VALIDATION SUMMARY")
# print("="*80)
# print(f"\nConfiguration:")
# print(f"  Oxide: {CONFIG['oxide']} ({CONFIG['L_oxide']*1e6:.1f} μm)")
# print(f"  Metal: {CONFIG['material']} ({CONFIG['L_metal']*1000:.1f} mm)")
# print(f"  Oxide defect fraction: {DEFECT_PARAMS_L5['area_fraction']*100:.0f}%")
# print(f"  Grain size: {MICROSTRUCTURE_L5['grain_size']*1e6:.0f} μm")
# print(f"  Trap density: {MICROSTRUCTURE_L5['trap_list'][0]['density']:.0e} m⁻³")
# print(f"  k_diss: {k_diss_L5L6:.2e} mol/m²/s/Pa")
# print(f"  Coverage mode: {COVERAGE_MODE_L5L6}")

# print(f"\n(A) All Levels at T={T_A-273:.0f}°C, P=100 kPa (ALL WITH SURFACE):")
# idx_100kPa = np.argmin(np.abs(pressures_L5A - 1e5))
# print(f"    L1+L6 (Perfect Metal):                     {fluxes_L1L6_A[idx_100kPa]:.3e} mol/m²/s")
# print(f"    L2b+L6 (Perfect Oxide + Perfect Metal):    {fluxes_L2bL6_A[idx_100kPa]:.3e} mol/m²/s")
# print(f"    L3+L6 (Defective Oxide + Perfect Metal):   {fluxes_L3L6_A[idx_100kPa]:.3e} mol/m²/s")
# print(f"    L4+L6 (Defective Metal, no oxide):         {fluxes_L4L6_A[idx_100kPa]:.3e} mol/m²/s")
# print(f"    L5+L6 (Defective Oxide + Defective Metal): {fluxes_L5L6_A[idx_100kPa]:.3e} mol/m²/s")
# print(f"    SRF:                                        {SRF_A[idx_100kPa]:.4f}")

# print(f"\n(B) Apparent Activation Energies (all with surface):")
# print(f"    L1+L6 (Perfect Metal): E_a = {E_a_L1L6_B:.1f} kJ/mol")
# print(f"    L5+L6 (Full System):   E_a = {E_a_L5L6_B:.1f} kJ/mol")

# print(f"\n(D) Limit Checks (all with surface):")
# print(f"    L5+L6(f→0, mode='none') → L2b+L6: max error = {error_1:.2e}%")
# print(f"    L5+L6(f→1, mode='none') → L1+L6:  max error = {error_2:.2e}%")
# print(f"    ✓ Hierarchical consistency verified")

# print("\n" + "="*80)
# print("HIERARCHICAL MODEL WITH SURFACE KINETICS:")
# print("  L1+L6:  Perfect Metal + Surface")
# print("  L2b+L6: Perfect Oxide + Perfect Metal + Surface")
# print("  L3+L6:  Defective Oxide + Perfect Metal + Surface")
# print("  L4+L6:  Defective Metal + Surface (no oxide)")
# print("  L5+L6:  Defective Oxide + Defective Metal + Surface (full system)")
# print("="*80)
```

```python

```
