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
import numpy as np
from scipy.optimize import brentq

# Import data dictionaries
from data.surface_kinetics_data import SURFACE_KINETICS, get_surface_kinetics
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS
R = 8.314  # J/(mol·K)
```

# L2+L6: Perfect Oxide + Prfect Metal with Surface Kinetics at Gas-Oxide Interface

<!-- #region -->
### The Three Fluxes

**1. Surface Dissociation Flux ($J_{surface}$):**

$$J_{surface} = k_{diss} \cdot P_{up} \cdot (1 - \theta_{ox})^2 - k_{recomb} \cdot \theta_{ox}^2$$

This is the net rate of H atoms entering the oxide from the gas phase.

At equilibrium:

$$k_{diss} \cdot P \cdot (1 - \theta)^2 = k_{recomb} \cdot \theta^2$$

Rearranging gives:

$$\frac{\theta}{1 - \theta} = \sqrt{\frac{k_{diss}}{k_{recomb}} \cdot P}$$

So if you define:

$$K_{eq} \equiv \frac{k_{diss}}{k_{recomb}}$$

then:

$$\theta = \frac{\sqrt{K_{eq} \cdot P}}{1 + \sqrt{K_{eq} \cdot P}}$$

Now, if you input $k_{diss}$ and $K_{eq}$, and compute:

$$k_{recomb} = \frac{k_{diss}}{K_{eq}}$$

everything remains mathematically and thermodynamically consistent.

To calculate $C_up$ 

From 
$$\frac{\theta}{1 - \theta} = \sqrt{\frac{k_{diss}}{k_{recomb}} \cdot P}$$

where can make P subject: 

$$ \sqrt(P) = \frac{1}{\sqrt(K_eq) } \cdot \frac{\theta}{1-\theta}$$

**2. Oxide Diffusion Flux ($J_{oxide}$):**

$$J_{oxide} = \frac{D_{ox}}{L_{ox}} \left( C_{ox,up} - C_{ox,int} \right)$$

Where:
- $C_{ox,up} = \frac{K_{ox}}{\sqrt{K_{eq}}} \cdot \frac{\theta_{ox}}{1 - \theta_{ox}}$ (from surface coverage)
- $C_{ox,int} = K_{ox} \cdot \sqrt{P_{int}}$ (Sieverts-like at oxide-metal interface)


**3. Metal Diffusion Flux ($J_{metal}$):**

$$J_{metal} = \frac{D_m \cdot K_{s,m}}{L_m} \left( \sqrt{P_{int}} - \sqrt{P_{down}} \right)$$


### Steady-State Condition

At steady state, all three fluxes must be equal:

$$\boxed{J_{surface} = J_{oxide} = J_{metal}}$$

This gives us **two equations** for **two unknowns**:
- $\theta_{ox}$ (surface coverage on oxide)
- $\sqrt{P_{int}}$ (interface pressure between oxide and metal)


### The Two Equations to Solve

**Equation 1: Surface = Oxide flux**

$$k_{diss} \cdot P_{up} \cdot (1 - \theta_{ox})^2 - k_{recomb} \cdot \theta_{ox}^2 = \frac{D_{ox}}{L_{ox}} \left( \frac{K_{ox}}{\sqrt{K_{eq}}} \cdot \frac{\theta_{ox}}{1 - \theta_{ox}} - K_{ox} \cdot \sqrt{P_{int}} \right)$$

**Equation 2: Oxide = Metal flux**

$$\frac{D_{ox}}{L_{ox}} \left( \frac{K_{ox}}{\sqrt{K_{eq}}} \cdot \frac{\theta_{ox}}{1 - \theta_{ox}} - K_{ox} \cdot \sqrt{P_{int}} \right) = \frac{D_m \cdot K_{s,m}}{L_m} \left( \sqrt{P_{int}} - \sqrt{P_{down}} \right)$$


### Solution Strategy

#### Step 1: Define Permeances and Functions

Let:
- $\alpha = \frac{D_{ox} \cdot K_{ox}}{L_{ox}}$ (oxide permeance)
- $\beta = \frac{D_m \cdot K_{s,m}}{L_m}$ (metal permeance)
- $g(\theta) = \frac{\theta}{(1 - \theta) \sqrt{K_{eq}}}$ (concentration function)

#### Step 2: Solve for $\sqrt{P_{int}}$ from Equation 2

From flux balance at oxide-metal interface:

$$\alpha \left( g(\theta_{ox}) - \sqrt{P_{int}} \right) = \beta \left( \sqrt{P_{int}} - \sqrt{P_{down}} \right)$$

Solving for $\sqrt{P_{int}}$:

$$\boxed{\sqrt{P_{int}} = \frac{\alpha \cdot g(\theta_{ox}) + \beta \cdot \sqrt{P_{down}}}{\alpha + \beta}}$$

#### Step 3: Substitute into Equation 1

Replace $\sqrt{P_{int}}(\theta_{ox})$ in Equation 1 to get a single equation in $\theta_{ox}$:

$$k_{diss} \cdot P_{up} \cdot (1 - \theta_{ox})^2 - k_{recomb} \cdot \theta_{ox}^2 = \alpha \left( g(\theta_{ox}) - \sqrt{P_{int}}(\theta_{ox}) \right)$$

#### Step 4: Numerical Root-Finding

Solve for $\theta_{ox}$ using `brentq` on the residual:

$$f(\theta_{ox}) = J_{surface}(\theta_{ox}) - J_{oxide}(\theta_{ox}) = 0$$

#### Step 5: Calculate Results

Once $\theta_{ox}$ is found:
1. Calculate $\sqrt{P_{int}}$ from the analytical expression
2. Calculate steady-state flux: $J_{ss} = \beta \left( \sqrt{P_{int}} - \sqrt{P_{down}} \right)$


### Rate-Limiting Analysis

The system can be characterized by resistances:

| Step | Resistance | Physical Meaning |
|------|------------|------------------|
| Surface | $R_{surface} \approx \frac{1}{k_{diss} \cdot P_{up}}$ | Dissociation kinetics |
| Oxide | $R_{oxide} = \frac{1}{\alpha} = \frac{L_{ox}}{D_{ox} \cdot K_{ox}}$ | Oxide diffusion |
| Metal | $R_{metal} = \frac{1}{\beta} = \frac{L_m}{D_m \cdot K_{s,m}}$ | Metal diffusion |

**Fractional contributions:**

$$f_i = \frac{R_i}{R_{total}}, \quad R_{total} = R_{surface} + R_{oxide} + R_{metal}$$

**Rate-limiting identification:**
- $f_{surface} > 0.5$ → Surface-limited (slope ≈ 1 in log J vs log P)
- $f_{oxide} > 0.5$ → Oxide-limited (slope ≈ 0.5)
- $f_{metal} > 0.5$ → Metal-limited (slope ≈ 0.5)
<!-- #endregion -->

```python
#Get all the known parameters
def get_all_properties(oxide_name, metal_name, temperature_K):
    """Get all properties needed for system calculation:
    - Surface kinetics
    - Oxide transport properties
    - Metal transport properties
    """

    # 1. Surface kinetics from oxide
    surf_kin = get_surface_kinetics(oxide_name, temperature_K)
    
    # 2. Oxide transport properties
    oxide_data = OXIDE_PROPERTIES[oxide_name]
    T_ref_ox = oxide_data['T_ref']
    D_ox = oxide_data['D_ox_ref'] * np.exp((-oxide_data['E_D_ox'] / R) * (1/temperature_K - 1/T_ref_ox))
    K_ox = oxide_data['K_ox_ref'] * np.exp((-oxide_data['H_sol_ox'] / R) * (1/temperature_K - 1/T_ref_ox))
    L_ox = oxide_data['thickness']
    
    # 3. Metal transport properties
    metal_data = MATERIALS[metal_name]
    T_ref_m = metal_data['T_ref']
    D_m = metal_data['D_ref'] * np.exp((-metal_data['E_D'] / R) * (1/temperature_K - 1/T_ref_m))
    K_s_m = metal_data['K_s_ref'] * np.exp((-metal_data['H_s'] / R) * (1/temperature_K - 1/T_ref_m))
    
    
    return {
        # Surface kinetics
        'k_diss': surf_kin['k_diss'],
        'k_recomb': surf_kin['k_recomb'],
        'K_eq': surf_kin['K_eq'],
        # Oxide
        'D_ox': D_ox, 'K_ox': K_ox, 'L_ox': L_ox,
        # Metal
        'D_m': D_m, 'K_s_m': K_s_m
    }
```

```python
#Solution step 1: Solve for sqrt(P_int)
def compute_permeances(props,L_m):
    """Compute oxide and metal permeances."""
    alpha = props['D_ox'] * props['K_ox'] / props['L_ox']  # Oxide permeance
    beta = props['D_m'] * props['K_s_m'] / L_m            # Metal permeance
    return alpha, beta

def g_theta(theta, K_eq):
    """
    Concentration function: g(θ) = θ / ((1-θ) × √K_eq)
    
    This converts surface coverage to effective √P at oxide surface.
    """
    if theta >= 1.0:
        return np.inf
    return theta / ((1.0 - theta) * np.sqrt(K_eq))

def sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down):
    """
    Solve for √P_int analytically from flux balance (Eq 2 = Eq 3).
    
    √P_int = (α × g(θ) + β × √P_down) / (α + β)
    """
    g = g_theta(theta, K_eq)
    sqrt_P_down = np.sqrt(P_down)
    return (alpha * g + beta * sqrt_P_down) / (alpha + beta)
```

```python
#Solution step 2: Substitute sqrt(P_int) into eqaution 1 and solve numerically to get theta
def surface_flux(theta, P_up, k_diss, K_eq):
    k_recomb = k_diss/K_eq
    J_surface = k_diss * P_up * (1 - theta)**2 - k_recomb * theta**2
    return J_surface

def oxide_flux(theta, alpha, beta, K_eq, P_down):
    g = g_theta(theta, K_eq)
    sqrt_P_int = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
    J_oxide = alpha * (g - sqrt_P_int)
    return J_oxide 

def metal_flux(theta, alpha, beta, K_eq, P_down):
    sqrt_P_int = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
    sqrt_P_down = np.sqrt(P_down)
    J_metal = beta * (sqrt_P_int - sqrt_P_down)
    return J_metal

def surface_flux_residual(theta, P_up, alpha, beta, P_down, k_diss, K_eq):
    """
    Residual for root-finding: J_surface - J_oxide = 0
    
    We substitute √P_int(θ) into J_oxide and solve for θ.
    """
    # Surface flux
    J_surface = surface_flux(theta, P_up, k_diss, K_eq)
    
    # Oxide flux (using g(θ) - √P_int)
    J_oxide = oxide_flux(theta, alpha, beta, K_eq, P_down)
    
    return J_surface - J_oxide

def solve_steady_state_flux(P_up, P_down, L_m, oxide_name, metal_name, temperature_K):
    """
    Solve coupled system with direct parameter input (no data file lookup).
    Includes rate-limiting analysis.
    """
    # Compute permeances
    props = get_all_properties(oxide_name, metal_name, temperature_K)
    alpha = props['D_ox'] * props['K_ox'] / props['L_ox']  # Oxide permeance
    beta = props['D_m'] * props['K_s_m'] / L_m   # Metal permeance

    k_diss  = props['k_diss']
    K_eq    = props['K_eq']

    # Solve for θ using brentq
    theta_ss = brentq(
        surface_flux_residual,
        1e-10, 1.0 - 1e-10,
        args=(P_up, alpha, beta, P_down, k_diss, K_eq)
    )
    
    # Calculate √P_int from θ
    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
    P_int = sqrt_P_int**2
    
    # Calculate steady-state flux
    J_ss = metal_flux(theta_ss, alpha, beta, K_eq, P_down)
    
    # Verify fluxes
    J_surf = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_ox = oxide_flux(theta_ss, alpha, beta, K_eq, P_down)
    

    # Surface resistance (linearized approximation)
    R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    
    # Oxide resistance
    R_oxide = 1.0 / alpha   # = L_ox / (D_ox × K_ox)
    
    # Metal resistance
    
    R_metal = 1.0 / beta
    
    # Total resistance
    R_total = R_surface + R_oxide + R_metal
    
    # Fractional contributions
    f_surface = R_surface / R_total if R_total > 0 else 0
    f_oxide = R_oxide / R_total if R_total > 0 else 0
    f_metal = R_metal / R_total if R_total > 0 else 0
    

    if f_surface > 0.5:
        rate_limiting = 'surface (dissociation)'
    elif f_oxide > 0.5:
        rate_limiting = 'oxide (diffusion)'
    else:
        rate_limiting = 'metal (diffusion)'

   
    
    return {
        'theta': theta_ss,
        'P_int': P_int,
        'J_ss': J_ss,
        'J_surface': J_surf,
        'J_oxide': J_ox,
        'J_metal': J_ss,
        'alpha': alpha,
        'beta': beta,
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface': R_surface,
            'R_oxide': R_oxide,
            'R_metal': R_metal,
            'R_total': R_total,
            'fraction_surface': f_surface,
            'fraction_oxide': f_oxide,
            'fraction_metal': f_metal,
        },
    }
```

```python
def solve_steady_state_flux_direct(P_up, P_down, L_m, k_diss, K_eq, D_ox, K_ox, L_ox, D_m, K_s_m):
    """
    Solve coupled system with direct parameter input (no data file lookup).
    Includes rate-limiting analysis.
    """
    # Compute permeances
    alpha = D_ox * K_ox / L_ox  # Oxide permeance
    beta = D_m * K_s_m / L_m    # Metal permeance
    
    # Solve for θ using brentq
    theta_ss = brentq(
        surface_flux_residual,
        1e-10, 1.0 - 1e-10,
        args=(P_up, alpha, beta, P_down, k_diss, K_eq)
    )
    
    # Calculate √P_int from θ
    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
    P_int = sqrt_P_int**2
    
    # Calculate steady-state flux
    J_ss = metal_flux(theta_ss, alpha, beta, K_eq, P_down)
    
    # Verify fluxes
    J_surf = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_ox = oxide_flux(theta_ss, alpha, beta, K_eq, P_down)
    
    # Surface resistance (linearized approximation)
    R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    
    # Oxide resistance
    R_oxide = 1.0 / alpha   # = L_ox / (D_ox × K_ox)
    
    # Metal resistance
    R_metal = 1.0 / beta
    
    # Total resistance
    R_total = R_surface + R_oxide + R_metal
    
    # Fractional contributions
    f_surface = R_surface / R_total if R_total > 0 else 0
    f_oxide = R_oxide / R_total if R_total > 0 else 0
    f_metal = R_metal / R_total if R_total > 0 else 0
    
    # FIX: Add 'mixed' case when no single resistance dominates
    if f_surface > 0.5:
        rate_limiting = 'surface'
    elif f_oxide > 0.5:
        rate_limiting = 'oxide'
    elif f_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'
 
    return {
        'theta': theta_ss,
        'P_int': P_int,
        'J_ss': J_ss,
        'J_surface': J_surf,
        'J_oxide': J_ox,
        'J_metal': J_ss,
        'alpha': alpha,
        'beta': beta,
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface': R_surface,
            'R_oxide': R_oxide,
            'R_metal': R_metal,
            'R_total': R_total,
            'fraction_surface': f_surface,
            'fraction_oxide': f_oxide,
            'fraction_metal': f_metal,
        },
    }
```

```python
import ipywidgets as widgets
from ipywidgets import interact
import matplotlib.pyplot as plt
import pandas as pd
from itertools import groupby
from operator import itemgetter

@interact(
    # Operating conditions
    P_down=widgets.FloatLogSlider(value=1e-10, base=10, min=-2, max=4, step=0.5, description='P_down (Pa)'),
    L_m=widgets.FloatLogSlider(value=1e-3, base=10, min=-4, max=-1, step=0.5, description='L_m (m)'),
    # Surface kinetics
    k_diss=widgets.FloatLogSlider(value=1e-15, base=10, min=-18, max=-3, step=0.5, description='k_diss'),
    K_eq=widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-1, step=0.5, description='K_eq'),
    # Oxide properties
    D_ox=widgets.FloatLogSlider(value=1e-11, base=10, min=-18, max=-5, step=0.5, description='D_ox (m²/s)'),
    K_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-14, max=-4, step=0.5, description='K_ox'),
    L_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-8, max=-4, step=0.5, description='L_ox (m)'),
    # Metal properties
    D_m=widgets.FloatLogSlider(value=1.0e-12, base=10, min=-13, max=-6, step=0.5, description='D_m (m²/s)'),
    K_s_m=widgets.FloatLogSlider(value=3.16e-4, base=10, min=-6, max=0, step=0.5, description='K_s_m'),
)
def interactive_L26_solver(P_down, L_m, k_diss, K_eq,
                           D_ox, K_ox, L_ox, D_m, K_s_m):
    """
    L2+L6 solver with single-pass loop: plot arrays and DataFrame rows
    built together from the same computed values.
    """

    P_up_range = np.logspace(0, 10, 40)  # 1 Pa to 10 GPa

    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    J_system = []
    theta_values = []
    f_surface_list = []
    f_oxide_list = []
    f_metal_list = []
    rows = []

    for P_up in P_up_range:
        try:
            r = solve_steady_state_flux_direct(
                P_up, P_down, L_m, k_diss, K_eq,
                D_ox, K_ox, L_ox, D_m, K_s_m
            )
            
            # Extract values once
            J_ss = r['J_ss']
            theta = r['theta']
            f_s = r['resistances']['fraction_surface']
            f_o = r['resistances']['fraction_oxide']
            f_m = r['resistances']['fraction_metal']
            
            # Derive rate-limiting label ONCE (same logic as function)
            if   f_s > 0.5: rate_lim = 'surface'
            elif f_o > 0.5: rate_lim = 'oxide'
            elif f_m > 0.5: rate_lim = 'metal'
            else:           rate_lim = 'mixed'
            
            # Append to plot arrays
            J_system.append(J_ss)
            theta_values.append(theta)
            f_surface_list.append(f_s)
            f_oxide_list.append(f_o)
            f_metal_list.append(f_m)
            
            # Append to DataFrame rows (same values, no recomputation)
            rows.append({
                "P_up (Pa)":           P_up,
                "P_int (Pa)":          r["P_int"],
                "J_ss (mol/m²/s)":     J_ss,
                "θ_surface":           theta,
                "f_surface (%)":       f_s * 100,
                "f_oxide (%)":         f_o * 100,
                "f_metal (%)":         f_m * 100,
                "Rate-Limiting":       rate_lim.upper(),
                "α_oxide":             r["alpha"],
                "β_metal":             r["beta"],
            })
            
        except Exception as e:
            J_system.append(np.nan)
            theta_values.append(np.nan)
            f_surface_list.append(np.nan)
            f_oxide_list.append(np.nan)
            f_metal_list.append(np.nan)
            rows.append({"P_up (Pa)": P_up, "Rate-Limiting": "ERROR", "Error": str(e)})

    # Convert to arrays
    J_system = np.array(J_system)
    theta_values = np.array(theta_values)
    f_surface = np.array(f_surface_list)
    f_oxide = np.array(f_oxide_list)
    f_metal = np.array(f_metal_list)

    # =========================================================================
    # Rate-limiting array for plot shading — 4-branch classification
    # =========================================================================
    rate_limiting_arr = np.where(
        f_surface > 0.5, 'surface',
        np.where(f_oxide > 0.5, 'oxide',
        np.where(f_metal > 0.5, 'metal',
                 'mixed'))
    )

    # =========================================================================
    # Create figure with 2 subplots
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # ===== Plot 1: Flux vs Pressure =====
    valid_idx = ~np.isnan(J_system)
    ax1.loglog(P_up_range, J_system, 'k-', linewidth=2.5, label='L2+L6 Model')

    if np.any(valid_idx):
        # Slope = 1 reference at the beginning
        P_ref1 = P_up_range[0]
        J_ref1 = J_system[valid_idx][0]
        J_slope1 = J_ref1 * (P_up_range / P_ref1) ** 1.0
        ax1.loglog(P_up_range, J_slope1, 'r--', linewidth=1.5, alpha=0.5, label='Slope = 1 (surface)')

        # Slope = 0.5 reference at the end
        P_ref05 = P_up_range[-1]
        J_ref05 = J_system[valid_idx][-1]
        J_slope05 = J_ref05 * (P_up_range / P_ref05) ** 0.5
        ax1.loglog(P_up_range, J_slope05, 'g--', linewidth=1.5, alpha=0.5, label='Slope = 0.5 (diffusion)')

    # FIX: 4 regions including 'mixed'
    regions = [
        {'mask': rate_limiting_arr == 'surface', 'color': 'red',    'label': 'Surface-limited'},
        {'mask': rate_limiting_arr == 'oxide',   'color': 'orange', 'label': 'Oxide-limited'},
        {'mask': rate_limiting_arr == 'metal',   'color': 'blue',   'label': 'Metal-limited'},
        {'mask': rate_limiting_arr == 'mixed',   'color': 'green',  'label': 'Mixed'},
    ]

    for region in regions:
        mask = region['mask'] & valid_idx
        if np.any(mask):
            idxs = np.where(mask)[0]
            for k, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
                group = list(map(itemgetter(1), g))
                if len(group) > 2:
                    P_seg = P_up_range[group]
                    J_seg = J_system[group]
                    ax1.loglog(P_seg, J_seg, color=region['color'], linewidth=4, alpha=0.7)
                    # Slope calculation
                    logP = np.log10(P_seg)
                    logJ = np.log10(np.abs(J_seg))
                    slope, _ = np.polyfit(logP, logJ, 1)
                    # Annotate
                    mid = len(group) // 2
                    ax1.text(P_seg[mid], J_seg[mid], f"{region['label']}\nSlope={slope:.2f}",
                             color=region['color'], fontsize=10, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                    
    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L2+L6: Flux vs Pressure', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    # ===== Plot 2: Rate-Limiting Fractions =====
    ax2.semilogx(P_up_range, f_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, f_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, f_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
    ax2.axhline(50, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
    
    ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax2.set_ylabel('Resistance Fraction (%)', fontsize=12)
    ax2.set_title('Rate-Limiting Step Analysis', fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')

    plt.tight_layout()
    plt.show()

    # =========================================================================
    # DataFrame — already built inside the loop, nothing to recompute
    # =========================================================================
    df = pd.DataFrame(rows)
    display(df)
```

```python
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS

@interact(
    # === GROUP 1: Material Selection ===
    oxide_name=widgets.Dropdown(options=list(OXIDE_PROPERTIES.keys()), description='Oxide'),
    metal_name=widgets.Dropdown(options=list(MATERIALS.keys()), description='Metal'),
    # === Everything else stays exactly as Reference 3 ===
    P_down=widgets.FloatLogSlider(value=1e-10, base=10, min=-2, max=4, step=0.5, description='P_down (Pa)'),
    L_m=widgets.FloatLogSlider(value=1e-3, base=10, min=-4, max=-1, step=0.5, description='L_m (m)'),
    k_diss=widgets.FloatLogSlider(value=1e-15, base=10, min=-18, max=-3, step=0.5, description='k_diss'),
    K_eq=widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-1, step=0.5, description='K_eq'),
    D_ox=widgets.FloatLogSlider(value=1e-11, base=10, min=-18, max=-5, step=0.5, description='D_ox (m²/s)'),
    K_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-14, max=-4, step=0.5, description='K_ox'),
    L_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-8, max=-4, step=0.5, description='L_ox (m)'),
    D_m=widgets.FloatLogSlider(value=1.0e-12, base=10, min=-13, max=-6, step=0.5, description='D_m (m²/s)'),
    K_s_m=widgets.FloatLogSlider(value=3.16e-4, base=10, min=-6, max=0, step=0.5, description='K_s_m'),
)
def interactive_L26_solver(oxide_name, metal_name,       # <-- NEW
                           P_down, L_m, k_diss, K_eq,
                           D_ox, K_ox, L_ox, D_m, K_s_m):

    P_up_range = np.logspace(0, 10, 40)

    J_system = []
    theta_values = []
    f_surface_list = []
    f_oxide_list = []
    f_metal_list = []

    for P_up in P_up_range:
        try:
            result = solve_steady_state_flux_direct(
                P_up, P_down, L_m, k_diss, K_eq,
                D_ox, K_ox, L_ox, D_m, K_s_m
            )
            J_system.append(result['J_ss'])
            theta_values.append(result['theta'])
            f_surface_list.append(result['resistances']['fraction_surface'])
            f_oxide_list.append(result['resistances']['fraction_oxide'])
            f_metal_list.append(result['resistances']['fraction_metal'])
        except:
            J_system.append(np.nan)
            theta_values.append(np.nan)
            f_surface_list.append(np.nan)
            f_oxide_list.append(np.nan)
            f_metal_list.append(np.nan)

    J_system = np.array(J_system)
    theta_values = np.array(theta_values)
    f_surface = np.array(f_surface_list)
    f_oxide = np.array(f_oxide_list)
    f_metal = np.array(f_metal_list)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    valid_idx = ~np.isnan(J_system)
    ax1.loglog(P_up_range, J_system, 'k-', linewidth=2.5, label='L2+L6 Model')

    if np.any(valid_idx):
        P_ref1 = P_up_range[0]
        J_ref1 = J_system[0] if not np.isnan(J_system[0]) else J_system[valid_idx][0]
        J_slope1 = J_ref1 * (P_up_range / P_ref1) ** 1.0
        ax1.loglog(P_up_range, J_slope1, 'r--', linewidth=1.5, alpha=0.5, label='Slope = 1 (surface)')

        P_ref05 = P_up_range[-1]
        J_ref05 = J_system[-1] if not np.isnan(J_system[-1]) else J_system[valid_idx][-1]
        J_slope05 = J_ref05 * (P_up_range / P_ref05) ** 0.5
        ax1.loglog(P_up_range, J_slope05, 'g--', linewidth=1.5, alpha=0.5, label='Slope = 0.5 (diffusion)')

    rate_limiting_arr = np.where(f_surface > 0.5, 'surface', np.where(f_oxide > 0.5, 'oxide', 'metal'))
    regions = [
        {'mask': rate_limiting_arr == 'surface', 'color': 'red',    'label': 'Surface-limited'},
        {'mask': rate_limiting_arr == 'metal',   'color': 'blue',   'label': 'Metal-limited'},
        {'mask': rate_limiting_arr == 'oxide',   'color': 'orange', 'label': 'Oxide-limited'},
    ]
    for region in regions:
        mask = region['mask'] & valid_idx
        if np.any(mask):
            from itertools import groupby
            from operator import itemgetter
            idxs = np.where(mask)[0]
            for k, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
                group = list(map(itemgetter(1), g))
                if len(group) > 2:
                    P_seg = P_up_range[group]
                    J_seg = J_system[group]
                    ax1.loglog(P_seg, J_seg, color=region['color'], linewidth=4, alpha=0.7)
                    slope, _ = np.polyfit(np.log10(P_seg), np.log10(np.abs(J_seg)), 1)
                    mid = len(group) // 2
                    ax1.text(P_seg[mid], J_seg[mid], f"{region['label']}\nSlope={slope:.2f}",
                             color=region['color'], fontsize=10, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L2+L6: Flux vs Pressure', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    ax2.semilogx(P_up_range, f_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, f_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, f_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
    ax2.axhline(50, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
    ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax2.set_ylabel('Resistance Fraction (%)', fontsize=12)
    ax2.set_title('Rate-Limiting Step Analysis', fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')

    plt.tight_layout()
    plt.show()

    import pandas as pd
    rows = []
    for P_up in P_up_range:
        try:
            r = solve_steady_state_flux_direct(P_up, P_down, L_m, k_diss, K_eq, D_ox, K_ox, L_ox, D_m, K_s_m)
            rows.append({
                "P_up (Pa)":       P_up,
                "P_int (Pa)":      r["P_int"],
                "J_ss (mol/m²/s)": r["J_ss"],
                "θ_surface":       r["theta"],
                "f_surface (%)":   r["resistances"]["fraction_surface"] * 100,
                "f_oxide (%)":     r["resistances"]["fraction_oxide"]   * 100,
                "f_metal (%)":     r["resistances"]["fraction_metal"]   * 100,
                "Rate-Limiting":   r["rate_limiting"].upper(),
                "α_oxide":         r["alpha"],
                "β_metal":         r["beta"],
            })
        except:
            rows.append({"P_up (Pa)": P_up})
    display(pd.DataFrame(rows))
```

```python
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS
from data.surface_kinetics_data import get_surface_kinetics

R_gas = 8.314  # J/(mol·K)

# === GROUP 1: Dropdowns ===
oxide_dd = widgets.Dropdown(options=list(OXIDE_PROPERTIES.keys()), description='Oxide')
metal_dd = widgets.Dropdown(options=list(MATERIALS.keys()),        description='Metal')

# === GROUP 2: Arrhenius sliders — initial values from first entry in database ===
_ox = OXIDE_PROPERTIES[oxide_dd.value]
_mt = MATERIALS[metal_dd.value]

T_K_slider     = widgets.IntSlider(value=773, min=300, max=1500, step=25, description='T (K)')
P_down_slider  = widgets.FloatLogSlider(value=1e-10, base=10, min=-2, max=4,  step=0.5, description='P_down (Pa)')
L_m_slider     = widgets.FloatLogSlider(value=1e-3,  base=10, min=-4, max=-1, step=0.5, description='L_m (m)')

D_ox_ref_slider  = widgets.FloatLogSlider(value=_ox['D_ox_ref'], base=10, min=-18, max=-5,  step=0.5, description='D_ox_ref')
E_D_ox_slider    = widgets.FloatSlider(value=_ox['E_D_ox'],  min=0, max=200000, step=1000, description='E_D_ox (J/mol)')
T_ref_ox_slider  = widgets.IntSlider(value=_ox['T_ref'],     min=300, max=1500, step=25,   description='T_ref_ox (K)')
K_ox_ref_slider  = widgets.FloatLogSlider(value=_ox['K_ox_ref'], base=10, min=-14, max=-2, step=0.5, description='K_ox_ref')
H_sol_ox_slider  = widgets.FloatSlider(value=_ox['H_sol_ox'], min=-100000, max=100000, step=1000, description='H_sol_ox (J/mol)')

D_ref_slider     = widgets.FloatLogSlider(value=_mt['D_ref'],   base=10, min=-18, max=-5,  step=0.5, description='D_ref')
E_D_slider       = widgets.FloatSlider(value=_mt['E_D'],    min=0, max=200000, step=1000,  description='E_D (J/mol)')
T_ref_m_slider   = widgets.IntSlider(value=_mt['T_ref'],    min=300, max=1500, step=25,    description='T_ref_m (K)')
K_s_ref_slider   = widgets.FloatLogSlider(value=_mt['K_s_ref'], base=10, min=-8, max=2,    step=0.5, description='K_s_ref')
H_s_slider       = widgets.FloatSlider(value=_mt['H_s'],    min=-100000, max=100000, step=1000, description='H_s (J/mol)')

# === GROUP 3: Read-only displays ===
out_D_ox   = widgets.FloatText(value=0, description='D_ox (m²/s)', disabled=True, style={'description_width': 'initial'})
out_K_ox   = widgets.FloatText(value=0, description='K_ox',        disabled=True, style={'description_width': 'initial'})
out_D_m    = widgets.FloatText(value=0, description='D_m (m²/s)',  disabled=True, style={'description_width': 'initial'})
out_K_s_m  = widgets.FloatText(value=0, description='K_s_m',       disabled=True, style={'description_width': 'initial'})
out_k_diss = widgets.FloatText(value=0, description='k_diss',      disabled=True, style={'description_width': 'initial'})
out_K_eq   = widgets.FloatText(value=0, description='K_eq',        disabled=True, style={'description_width': 'initial'})

# === Wire dropdowns → update Group 2 sliders from database ===
def update_arrhenius_from_dropdown(change):
    _ox = OXIDE_PROPERTIES[oxide_dd.value]
    _mt = MATERIALS[metal_dd.value]
    # Oxide
    D_ox_ref_slider.value = _ox['D_ox_ref']
    E_D_ox_slider.value   = _ox['E_D_ox']
    T_ref_ox_slider.value = _ox['T_ref']
    K_ox_ref_slider.value = _ox['K_ox_ref']
    H_sol_ox_slider.value = _ox['H_sol_ox']
    # Metal
    D_ref_slider.value    = _mt['D_ref']
    E_D_slider.value      = _mt['E_D']
    T_ref_m_slider.value  = _mt['T_ref']
    K_s_ref_slider.value  = _mt['K_s_ref']
    H_s_slider.value      = _mt['H_s']

oxide_dd.observe(update_arrhenius_from_dropdown, names='value')
metal_dd.observe(update_arrhenius_from_dropdown, names='value')

# Display Group 3 panel
display(widgets.HTML("<b>── Computed Parameters at T (read-only) ──</b>"))
display(widgets.HBox([out_D_ox, out_K_ox, out_D_m, out_K_s_m, out_k_diss, out_K_eq]))

@interact(
    oxide_name = oxide_dd,
    metal_name = metal_dd,
    T_K        = T_K_slider,
    P_down     = P_down_slider,
    L_m        = L_m_slider,
    D_ox_ref   = D_ox_ref_slider,
    E_D_ox     = E_D_ox_slider,
    T_ref_ox   = T_ref_ox_slider,
    K_ox_ref   = K_ox_ref_slider,
    H_sol_ox   = H_sol_ox_slider,
    D_ref      = D_ref_slider,
    E_D        = E_D_slider,
    T_ref_m    = T_ref_m_slider,
    K_s_ref    = K_s_ref_slider,
    H_s        = H_s_slider,
)
def interactive_L26_solver(oxide_name, metal_name,
                           T_K, P_down, L_m,
                           D_ox_ref, E_D_ox, T_ref_ox, K_ox_ref, H_sol_ox,
                           D_ref, E_D, T_ref_m, K_s_ref, H_s):

    # === Compute parameters from Arrhenius ===
    D_ox  = D_ox_ref * np.exp((-E_D_ox   / R_gas) * (1/T_K - 1/T_ref_ox))
    K_ox  = K_ox_ref * np.exp((-H_sol_ox / R_gas) * (1/T_K - 1/T_ref_ox))
    D_m   = D_ref    * np.exp((-E_D      / R_gas) * (1/T_K - 1/T_ref_m))
    K_s_m = K_s_ref  * np.exp((-H_s      / R_gas) * (1/T_K - 1/T_ref_m))
    L_ox  = OXIDE_PROPERTIES[oxide_name]['thickness']

    surf   = get_surface_kinetics(oxide_name, T_K)
    k_diss = surf['k_diss']
    K_eq   = surf['K_eq']

    # === Update Group 3 read-only displays ===
    out_D_ox.value   = D_ox
    out_K_ox.value   = K_ox
    out_D_m.value    = D_m
    out_K_s_m.value  = K_s_m
    out_k_diss.value = k_diss
    out_K_eq.value   = K_eq

    # === Main loop (unchanged) ===
    P_up_range = np.logspace(0, 10, 40)
    J_system, theta_values = [], []
    f_surface_list, f_oxide_list, f_metal_list = [], [], []

    for P_up in P_up_range:
        try:
            result = solve_steady_state_flux_direct(
                P_up, P_down, L_m, k_diss, K_eq,
                D_ox, K_ox, L_ox, D_m, K_s_m
            )
            J_system.append(result['J_ss'])
            theta_values.append(result['theta'])
            f_surface_list.append(result['resistances']['fraction_surface'])
            f_oxide_list.append(result['resistances']['fraction_oxide'])
            f_metal_list.append(result['resistances']['fraction_metal'])
        except:
            J_system.append(np.nan)
            theta_values.append(np.nan)
            f_surface_list.append(np.nan)
            f_oxide_list.append(np.nan)
            f_metal_list.append(np.nan)

    J_system  = np.array(J_system)
    f_surface = np.array(f_surface_list)
    f_oxide   = np.array(f_oxide_list)
    f_metal   = np.array(f_metal_list)

    # === Plots (unchanged) ===
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    valid_idx = ~np.isnan(J_system)
    ax1.loglog(P_up_range, J_system, 'k-', linewidth=2.5, label='L2+L6 Model')

    if np.any(valid_idx):
        P_ref1 = P_up_range[0]
        J_ref1 = J_system[0] if not np.isnan(J_system[0]) else J_system[valid_idx][0]
        ax1.loglog(P_up_range, J_ref1 * (P_up_range / P_ref1)**1.0, 'r--', linewidth=1.5, alpha=0.5, label='Slope = 1 (surface)')
        P_ref05 = P_up_range[-1]
        J_ref05 = J_system[-1] if not np.isnan(J_system[-1]) else J_system[valid_idx][-1]
        ax1.loglog(P_up_range, J_ref05 * (P_up_range / P_ref05)**0.5, 'g--', linewidth=1.5, alpha=0.5, label='Slope = 0.5 (diffusion)')

    rate_limiting_arr = np.where(f_surface > 0.8, 'surface', np.where(f_oxide > 0.8, 'oxide', 'metal'))
    regions = [
        {'mask': rate_limiting_arr == 'surface', 'color': 'red',    'label': 'Surface-limited'},
        {'mask': rate_limiting_arr == 'metal',   'color': 'blue',   'label': 'Metal-limited'},
        {'mask': rate_limiting_arr == 'oxide',   'color': 'orange', 'label': 'Oxide-limited'},
    ]
    for region in regions:
        mask = region['mask'] & valid_idx
        if np.any(mask):
            from itertools import groupby
            from operator import itemgetter
            idxs = np.where(mask)[0]
            for k, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
                group = list(map(itemgetter(1), g))
                if len(group) > 2:
                    P_seg, J_seg = P_up_range[group], J_system[group]
                    ax1.loglog(P_seg, J_seg, color=region['color'], linewidth=4, alpha=0.7)
                    slope, _ = np.polyfit(np.log10(P_seg), np.log10(np.abs(J_seg)), 1)
                    mid = len(group) // 2
                    ax1.text(P_seg[mid], J_seg[mid], f"{region['label']}\nSlope={slope:.2f}",
                             color=region['color'], fontsize=10, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title(f'L2+L6: Flux vs Pressure  |  {oxide_name} / {metal_name}  |  T={T_K} K', fontsize=13)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    ax2.semilogx(P_up_range, f_surface * 100, 'r-',          linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, f_oxide   * 100, color='orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, f_metal   * 100, 'b-',          linewidth=2, label='Metal (diffusion)')
    ax2.axhline(50, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
    ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax2.set_ylabel('Resistance Fraction (%)', fontsize=12)
    ax2.set_title('Rate-Limiting Step Analysis', fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')

    plt.tight_layout()
    plt.show()

    # === Table (unchanged) ===
    import pandas as pd
    rows = []
    for P_up in P_up_range:
        try:
            r = solve_steady_state_flux_direct(P_up, P_down, L_m, k_diss, K_eq, D_ox, K_ox, L_ox, D_m, K_s_m)
            rows.append({
                "P_up (Pa)":       P_up,
                "P_int (Pa)":      r["P_int"],
                "J_ss (mol/m²/s)": r["J_ss"],
                "θ_surface":       r["theta"],
                "f_surface (%)":   r["resistances"]["fraction_surface"] * 100,
                "f_oxide (%)":     r["resistances"]["fraction_oxide"]   * 100,
                "f_metal (%)":     r["resistances"]["fraction_metal"]   * 100,
                "Rate-Limiting":   r["rate_limiting"].upper(),
                "α_oxide":         r["alpha"],
                "β_metal":         r["beta"],
            })
        except:
            rows.append({"P_up (Pa)": P_up})
    display(pd.DataFrame(rows))
```

### Perfect Oxide + Defective Metal Flux with Surface Chemistry

```python
# =============================================================================
# LEVEL 4 + L6: Perfect Oxide + Defective Metal Flux with Surface Chemistry
# =============================================================================

def calculate_defective_metal_flux_L6(
    P_up, P_down, thickness, temperature,
    # Surface kinetics parameters
    k_diss, K_eq,
    # Oxide parameters  
    D_ox, K_ox, L_ox,
    # Metal parameters
    D_lattice, K_s_m,
    # Microstructure parameters
    microstructure_params,
    lattice_density=1.06e29,
    method='average', n_points=50, mode='both',
    max_iterations=15, tolerance=1e-6
):
    """
    Calculate hydrogen permeation flux through defective metal with surface chemistry.
    
    Combines:
    - L6: Surface kinetics at gas-oxide interface (θ-based)
    - L2: Oxide diffusion (atomic H, √P dependence)
    - L4: Defective metal diffusion (GB enhancement + trapping)
    
    The key difference from calculate_defective_metal_flux():
    - Upstream concentration is NOT from Sieverts' law
    - Instead, solve coupled system: J_surface = J_oxide = J_metal
    - Surface coverage θ determines subsurface concentration
    - D_eff and θ are solved self-consistently via iteration
    
    Parameters
    ----------
    P_up : float
        Upstream H2 pressure [Pa]
    P_down : float
        Downstream H2 pressure [Pa]
    thickness : float
        Metal thickness [m]
    temperature : float
        Temperature [K]
    k_diss : float
        Dissociation rate constant [mol/m²/s/Pa]
    K_eq : float
        Surface equilibrium constant [Pa⁻¹]
    D_ox : float
        Oxide diffusivity [m²/s]
    K_ox : float
        Oxide solubility [mol/m³/Pa^0.5]
    L_ox : float
        Oxide thickness [m]
    D_lattice : float
        Metal lattice diffusivity [m²/s]
    K_s_m : float
        Metal solubility [mol/m³/Pa^0.5]
    microstructure_params : dict
        Microstructure specification (same as calculate_defective_metal_flux)
    lattice_density : float, optional
        Lattice site density [sites/m³], default 1.06e29
    method : str, optional
        Method for averaging D_eff: 'average', 'harmonic', 'inlet', 'outlet'
    n_points : int, optional
        Number of points for concentration profile discretization
    mode : str, optional
        Microstructure mode: 'both', 'gb_only', 'trapping_only', 'none'
    max_iterations : int, optional
        Maximum iterations for self-consistent D_eff convergence
    tolerance : float, optional
        Relative tolerance for D_eff convergence
    
    Returns
    -------
    dict
        Dictionary containing flux, concentrations, θ, P_int, D_eff, etc.
    """
    from calculations.defective_metal import combined_microstructure_model
    
    # =========================================================================
    # Input Validation
    # =========================================================================
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    if thickness <= 0:
        raise ValueError(f"Thickness must be positive: {thickness} m")
    if D_lattice <= 0:
        raise ValueError(f"Diffusion coefficient must be positive")
    if k_diss <= 0 or K_eq <= 0:
        raise ValueError("Surface kinetics parameters must be positive")
    
    # =========================================================================
    # Step 1: Compute oxide permeance (fixed)
    # =========================================================================
    alpha = D_ox * K_ox / L_ox  # Oxide permeance
    
    # =========================================================================
    # Step 2: Iterate to find self-consistent D_eff and theta_ss
    # =========================================================================
    # The key insight: theta_ss depends on beta, and D_eff depends on the
    # concentration profile which depends on theta_ss. We must iterate.
    
    D_eff = D_lattice  # Initial guess
    sqrt_P_down = np.sqrt(max(P_down, 0))
    
    convergence_history = []
    
    for iteration in range(max_iterations):
        # Current metal permeance with current D_eff
        beta = D_eff * K_s_m / thickness
        
        # Solve for θ with current beta
        try:
            theta_ss = brentq(
                surface_flux_residual,
                1e-10, 1.0 - 1e-10,
                args=(P_up, alpha, beta, P_down, k_diss, K_eq)
            )
        except ValueError as e:
            return {
                'flux': np.nan,
                'error': f'Failed to solve for theta at iteration {iteration}: {str(e)}',
                'convergence_history': convergence_history
            }
        
        # Calculate sqrt_P_int consistent with this theta and beta
        sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
        
        # Build concentration profile through metal
        x_array = np.linspace(0, thickness, n_points)
        sqrt_P_array = sqrt_P_int - (sqrt_P_int - sqrt_P_down) * x_array / thickness
        C_array = K_s_m * sqrt_P_array
        
        # Calculate D_eff at each position using microstructure model
        D_array = np.zeros(n_points)
        theta_trap_array = np.zeros(n_points)
        gb_factor_array = np.zeros(n_points)
        
        for i, C_local in enumerate(C_array):
            C_local = max(C_local, 1e-20)
            
            result_i = combined_microstructure_model(
                D_lattice=D_lattice,
                temperature=temperature,
                microstructure_params=microstructure_params,
                lattice_concentration=C_local,
                lattice_density=lattice_density,
                mode=mode
            )
            
            D_array[i] = result_i['D_eff']
            
            # Extract trapping info
            if 'trapping' in result_i and result_i['trapping'] is not None:
                theta_trap_array[i] = result_i['trapping'].get('theta_total', 0.0)
            elif 'theta_total' in result_i:
                theta_trap_array[i] = result_i['theta_total']
            else:
                theta_trap_array[i] = 0.0
                
            # Extract GB enhancement info
            if 'gb_enhancement' in result_i and result_i['gb_enhancement'] is not None:
                gb_factor_array[i] = result_i['gb_enhancement'].get('factor', 1.0)
            elif 'gb_enhancement_factor' in result_i:
                gb_factor_array[i] = result_i['gb_enhancement_factor']
            else:
                gb_factor_array[i] = 1.0
        
        # Calculate new D_eff based on method
        if method == 'average':
            D_eff_new = np.mean(D_array)
        elif method == 'harmonic':
            D_eff_new = len(D_array) / np.sum(1.0 / D_array)
        elif method == 'inlet':
            D_eff_new = D_array[0]
        elif method == 'outlet':
            D_eff_new = D_array[-1]
        else:
            D_eff_new = np.mean(D_array)
        
        # Track convergence
        convergence_history.append({
            'iteration': iteration,
            'D_eff': D_eff_new,
            'theta': theta_ss,
            'relative_change': abs(D_eff_new - D_eff) / D_eff if D_eff > 0 else np.inf
        })
        
        # Check convergence
        if abs(D_eff_new - D_eff) / D_eff < tolerance:
            D_eff = D_eff_new
            break
        
        D_eff = D_eff_new
    
    # =========================================================================
    # Step 3: Final calculations with converged D_eff and theta_ss
    # =========================================================================
    
    # Final beta_eff
    beta_eff = D_eff * K_s_m / thickness
    
    # Final P_int (recalculate with converged values)
    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta_eff, K_eq, P_down)
    P_int = sqrt_P_int**2
    
    # Concentrations at metal boundaries
    C_up = K_s_m * sqrt_P_int  # mol/m³ at metal inlet (from P_int)
    C_down = K_s_m * sqrt_P_down if P_down > 0 else 0.0
    
    # =========================================================================
    # Step 4: Calculate flux through metal (self-consistent)
    # =========================================================================
    
    # Metal flux using effective diffusivity - now theta_ss and beta_eff are consistent
    J_metal = metal_flux(theta_ss, alpha, beta_eff, K_eq, P_down)
    
    # Verify with surface and oxide fluxes
    J_surface = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_oxide = oxide_flux(theta_ss, alpha, beta_eff, K_eq, P_down)
    
    # =========================================================================
    # Step 5: Diagnostic information
    # =========================================================================
    
    modification_factor = D_eff / D_lattice
    avg_gb_factor = np.mean(gb_factor_array)
    avg_theta_trap = np.mean(theta_trap_array)
    trap_reduction = 1.0 / (1.0 + avg_theta_trap) if avg_theta_trap > 0 else 1.0
    
    # Determine dominant effect
    if avg_gb_factor > 1.5 and trap_reduction > 0.5:
        dominant_effect = 'gb_enhancement'
    elif avg_gb_factor < 1.5 and trap_reduction < 0.5:
        dominant_effect = 'trapping'
    else:
        dominant_effect = 'balanced'

    # =========================================================================
    # Step 6: Rate-Limiting Analysis
    # =========================================================================

    # Surface resistance (linearized approximation)
    R_surface_approx = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0

    # Oxide resistance
    R_oxide = 1.0 / alpha  # = L_ox / (D_ox × K_ox)

    # Metal resistance (with microstructure)
    R_metal = 1.0 / beta_eff  # = thickness / (D_eff × K_s_m)

    # Total resistance
    R_total = R_surface_approx + R_oxide + R_metal

    # Fractional contributions
    f_surface = R_surface_approx / R_total if R_total > 0 else 0
    f_oxide = R_oxide / R_total if R_total > 0 else 0
    f_metal = R_metal / R_total if R_total > 0 else 0

    # Determine rate-limiting step (>50% of resistance)
    if f_surface > 0.5:
        rate_limiting = 'surface'
    elif f_oxide > 0.5:
        rate_limiting = 'oxide'
    elif f_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'

    # =========================================================================
    # Step 7: Return Results
    # =========================================================================
    
    return {
        # Primary outputs
        'flux': J_metal,
        'J_surface': J_surface,
        'J_oxide': J_oxide,
        'J_metal': J_metal,
        
        # Surface chemistry outputs (L6)
        'theta_surface': theta_ss,  # Surface coverage at gas-oxide interface
        'P_int': P_int,             # Pressure at oxide-metal interface
        
        # Concentration outputs
        'C_up': C_up,               # Concentration at metal inlet (from P_int)
        'C_down': C_down,           # Concentration at metal outlet
        
        # Permeances
        'alpha': alpha,             # Oxide permeance
        'beta_lattice': D_lattice * K_s_m / thickness,  # Lattice-only metal permeance
        'beta_eff': beta_eff,       # Effective metal permeance (with microstructure)
        
        # Level 4 outputs (microstructure)
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'modification_factor': modification_factor,
        
        # Microstructure details
        'microstructure_details': {
            'average_theta_trap': avg_theta_trap,
            'average_gb_enhancement': avg_gb_factor,
            'trap_reduction_factor': trap_reduction,
            'dominant_effect': dominant_effect,
            'method_used': method,
        },
        
        # Flux verification
        'flux_balance': {
            'J_surface': J_surface,
            'J_oxide': J_oxide,
            'J_metal': J_metal,
            'balanced': np.allclose(J_surface, J_oxide, rtol=1e-4) and np.allclose(J_oxide, J_metal, rtol=1e-4)
        },
        
        # Profile data
        'profiles': {
            'x': x_array,
            'D': D_array,
            'C': C_array,
            'theta_trap': theta_trap_array,
            'gb_factor': gb_factor_array
        },

        # Rate-limiting analysis
        'rate_limiting': rate_limiting,
        'resistances': {
            'R_surface': R_surface_approx,
            'R_oxide': R_oxide,
            'R_metal': R_metal,
            'R_total': R_total,
            'fraction_surface': f_surface,
            'fraction_oxide': f_oxide,
            'fraction_metal': f_metal,
        },
        
        # Convergence info
        'convergence': {
            'iterations': len(convergence_history),
            'converged': len(convergence_history) < max_iterations,
            'history': convergence_history
        },
        
        # Units
        'units': {
            'flux': 'mol/m²/s',
            'concentration': 'mol/m³',
            'pressure': 'Pa',
            'diffusivity': 'm²/s',
            'theta': 'dimensionless'
        }
    }
```

```python
# =============================================================================
# Interactive Widget for L6 + L4 Combined Model
# Perfect Oxide + Defective Metal Flux with Surface Chemistry
# =============================================================================

import pandas as pd
from itertools import groupby
from operator import itemgetter

# Define default microstructure (can be modified in the widget)
default_microstructure = {
    'grain_size': 500e-6,
    'grain_shape': 'equiaxed',
    'gb_type': 'LAGB',
    'trap_list': [
        {'name': 'dislocations', 'density': 1e15, 'binding_energy': 27e3}
    ]
}

@interact(
    # Operating conditions
    P_down=widgets.FloatLogSlider(value=1e-10, base=10, min=-2, max=4, step=0.5, description='P_down (Pa)'),
    thickness=widgets.FloatLogSlider(value=1e-3, base=10, min=-4, max=-1, step=0.5, description='L_m (m)'),
    temperature=widgets.IntSlider(value=973, min=573, max=1273, step=50, description='T (K)'),
    # Surface kinetics
    k_diss=widgets.FloatLogSlider(value=1e-11, base=10, min=-18, max=-3, step=0.5, description='k_diss'),
    K_eq=widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-1, step=0.5, description='K_eq'),
    # Oxide properties
    D_ox=widgets.FloatLogSlider(value=1e-9, base=10, min=-18, max=-5, step=0.5, description='D_ox (m²/s)'),
    K_ox=widgets.FloatLogSlider(value=3.16e-6, base=10, min=-14, max=-4, step=0.5, description='K_ox'),
    L_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-8, max=-4, step=0.5, description='L_ox (m)'),
    # Metal properties
    D_lattice=widgets.FloatLogSlider(value=1.0e-10, base=10, min=-12, max=-6, step=0.5, description='D_lattice'),
    K_s_m=widgets.FloatLogSlider(value=3.16e-2, base=10, min=-4, max=0, step=0.5, description='K_s_m'),
    # Microstructure
    grain_size=widgets.FloatLogSlider(value=31e-6, base=10, min=-6, max=-3, step=0.5, description='Grain (m)'),
    trap_density=widgets.FloatLogSlider(value=3.16e15, base=10, min=12, max=18, step=0.5, description='ρ_trap (m⁻²)'),
)
def interactive_L246_solver(P_down, thickness, temperature, k_diss, K_eq,
                            D_ox, K_ox, L_ox, D_lattice, K_s_m,
                            grain_size, trap_density):
    """
    L2+L4+L6 solver with single-pass loop: plot arrays and DataFrame rows
    built together from the same computed values.
    """
    
    # Build microstructure dict from widget inputs
    microstructure = {
        'grain_size': grain_size,
        'grain_shape': 'equiaxed',
        'gb_type': 'LAGB',
        'trap_list': [
            {'name': 'dislocations', 'density': trap_density, 'binding_energy': 27e3}
        ]
    }
    
    P_up_range = np.logspace(0, 12, 40)  # 1 Pa to 1 TPa

    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    J_system = []
    theta_values = []
    f_surface_list = []
    f_oxide_list = []
    f_metal_list = []
    rows = []

    for P_up in P_up_range:
        try:
            r = calculate_defective_metal_flux_L6(
                P_up=P_up, P_down=P_down, thickness=thickness, temperature=temperature,
                k_diss=k_diss, K_eq=K_eq,
                D_ox=D_ox, K_ox=K_ox, L_ox=L_ox,
                D_lattice=D_lattice, K_s_m=K_s_m,
                microstructure_params=microstructure
            )
            
            # Extract values once
            J_ss = r['flux']
            theta = r['theta_surface']
            f_s = r['resistances']['fraction_surface']
            f_o = r['resistances']['fraction_oxide']
            f_m = r['resistances']['fraction_metal']
            
            # Derive rate-limiting label ONCE (same logic as function)
            if   f_s > 0.5: rate_lim = 'surface'
            elif f_o > 0.5: rate_lim = 'oxide'
            elif f_m > 0.5: rate_lim = 'metal'
            else:           rate_lim = 'mixed'
            
            # Append to plot arrays
            J_system.append(J_ss)
            theta_values.append(theta)
            f_surface_list.append(f_s)
            f_oxide_list.append(f_o)
            f_metal_list.append(f_m)
            
            # Append to DataFrame rows (same values, no recomputation)
            rows.append({
                "P_up (Pa)":           P_up,
                "P_int (Pa)":          r["P_int"],
                "J_ss (mol/m²/s)":     J_ss,
                "θ_surface":           theta,
                "D_eff (m²/s)":        r["D_eff"],
                "D_eff/D_lattice":     r["modification_factor"],
                "f_surface (%)":       f_s * 100,
                "f_oxide (%)":         f_o * 100,
                "f_metal (%)":         f_m * 100,
                "Rate-Limiting":       rate_lim.upper(),
                "α_oxide":             r["alpha"],
                "β_eff":               r["beta_eff"],
            })
            
        except Exception as e:
            J_system.append(np.nan)
            theta_values.append(np.nan)
            f_surface_list.append(np.nan)
            f_oxide_list.append(np.nan)
            f_metal_list.append(np.nan)
            rows.append({"P_up (Pa)": P_up, "Rate-Limiting": "ERROR", "Error": str(e)})

    # Convert to arrays
    J_system = np.array(J_system)
    theta_values = np.array(theta_values)
    f_surface = np.array(f_surface_list)
    f_oxide = np.array(f_oxide_list)
    f_metal = np.array(f_metal_list)

    # =========================================================================
    # Rate-limiting array for plot shading — 4-branch classification
    # =========================================================================
    rate_limiting_arr = np.where(
        f_surface > 0.5, 'surface',
        np.where(f_oxide > 0.5, 'oxide',
        np.where(f_metal > 0.5, 'metal',
                 'mixed'))
    )




# Diagnostic: verify Plot 1 shading matches Plot 2 crossings
    print("\n=== Rate-Limiting Consistency Check ===")
    for i in [0, len(P_up_range)//2, -1]:  # Check start, middle, end
        P = P_up_range[i]
        fs, fo, fm = f_surface[i]*100, f_oxide[i]*100, f_metal[i]*100
        rl = rate_limiting_arr[i]
        print(f"P={P:.1e}: f_surf={fs:.1f}%, f_ox={fo:.1f}%, f_metal={fm:.1f}% → {rl}")
        
        # Verify classification matches
        if fs > 50 and rl != 'surface':
            print(f"  ⚠️ MISMATCH: f_surface > 50% but rate_limiting={rl}")
        if fo > 50 and rl != 'oxide':
            print(f"  ⚠️ MISMATCH: f_oxide > 50% but rate_limiting={rl}")
        if fm > 50 and rl != 'metal':
            print(f"  ⚠️ MISMATCH: f_metal > 50% but rate_limiting={rl}")

    # =========================================================================
    # Create figure with 2 subplots
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # ===== Plot 1: Flux vs Pressure =====
    valid_idx = ~np.isnan(J_system)
    ax1.loglog(P_up_range, J_system, 'k-', linewidth=2.5, label='L2+L4+L6 Model')

    if np.any(valid_idx):
        # Slope = 1 reference at the beginning
        P_ref1 = P_up_range[0]
        J_ref1 = J_system[valid_idx][0]
        J_slope1 = J_ref1 * (P_up_range / P_ref1) ** 1.0
        ax1.loglog(P_up_range, J_slope1, 'r--', linewidth=1.5, alpha=0.5, label='Slope = 1 (surface)')

        # Slope = 0.5 reference at the end
        P_ref05 = P_up_range[-1]
        J_ref05 = J_system[valid_idx][-1]
        J_slope05 = J_ref05 * (P_up_range / P_ref05) ** 0.5
        ax1.loglog(P_up_range, J_slope05, 'g--', linewidth=1.5, alpha=0.5, label='Slope = 0.5 (diffusion)')

    # FIX: 4 regions including 'mixed'
    regions = [
        {'mask': rate_limiting_arr == 'surface', 'color': 'red',    'label': 'Surface-limited'},
        {'mask': rate_limiting_arr == 'oxide',   'color': 'orange', 'label': 'Oxide-limited'},
        {'mask': rate_limiting_arr == 'metal',   'color': 'blue',   'label': 'Metal-limited'},
        {'mask': rate_limiting_arr == 'mixed',   'color': 'green',  'label': 'Mixed'},
    ]

    for region in regions:
        mask = region['mask'] & valid_idx
        if np.any(mask):
            idxs = np.where(mask)[0]
            for k, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
                group = list(map(itemgetter(1), g))
                if len(group) > 2:
                    P_seg = P_up_range[group]
                    J_seg = J_system[group]
                    ax1.loglog(P_seg, J_seg, color=region['color'], linewidth=4, alpha=0.7)
                    # Slope calculation
                    logP = np.log10(P_seg)
                    logJ = np.log10(np.abs(J_seg))
                    slope, _ = np.polyfit(logP, logJ, 1)
                    # Annotate
                    mid = len(group) // 2
                    ax1.text(P_seg[mid], J_seg[mid], f"{region['label']}\nSlope={slope:.2f}",
                             color=region['color'], fontsize=10, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L2+L4+L6: Flux vs Pressure', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    # ===== Plot 2: Rate-Limiting Fractions =====
    ax2.semilogx(P_up_range, f_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, f_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, f_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
    ax2.axhline(50, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
    
    ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax2.set_ylabel('Resistance Fraction (%)', fontsize=12)
    ax2.set_title('Rate-Limiting Step Analysis', fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')

    plt.tight_layout()
    plt.show()

    # =========================================================================
    # DataFrame — already built inside the loop, nothing to recompute
    # =========================================================================
    df = pd.DataFrame(rows)
    display(df)
```

## Defective Oxide + Perfect Metal


### Plan: Parallel Path Model with Surface Chemistry (L3 + L6)

### Overview

Combine the **Strehlow & Savage (1974)** parallel path model for defective oxides with the **L6 surface kinetics** model. This creates a comprehensive framework where:

- **Intact oxide regions**: Full L2+L6 coupled model (surface kinetics + oxide diffusion + metal diffusion)
- **Defect regions**: Modified transport depending on defect type


<!-- #region -->

### Step 1: Define the Parallel Path Framework

#### Total Flux (Area-Weighted Sum)

$$\boxed{J_{total} = f_{intact} \cdot J_{intact} + f_{defect} \cdot J_{defect}}$$

Where:
- $f_{intact} = 1 - f_{defect}$ = area fraction of intact oxide
- $f_{defect}$ = area fraction with defects (pinholes, cracks, grain boundaries)
- $J_{intact}$ = flux through intact oxide+metal (from L2+L6 model)
- $J_{defect}$ = flux through defect regions

### Step 2: Intact Oxide Path (L2 + L6)

#### Surface Kinetics (L6) at Gas-Oxide Interface

$$J_{surface} = k_{diss} \cdot P_{up} \cdot (1 - \theta)^2 - k_{recomb} \cdot \theta^2$$

Where $k_{recomb} = k_{diss} / K_{eq}$

#### Oxide Diffusion (L2)

$$J_{oxide} = \alpha_{intact} \cdot \left( g(\theta) - \sqrt{P_{int}} \right)$$

Where:
- $\alpha_{intact} = \frac{D_{ox} \cdot K_{ox}}{L_{ox}}$ (intact oxide permeance)
- $g(\theta) = \frac{\theta}{(1 - \theta) \cdot \sqrt{K_{eq}}}$

#### Metal Diffusion

$$J_{metal} = \beta \cdot \left( \sqrt{P_{int}} - \sqrt{P_{down}} \right)$$

Where $\beta = \frac{D_m \cdot K_{s,m}}{L_m}$

#### Solve for $\theta_{intact}$ and $P_{int,intact}$

From flux balance $J_{oxide} = J_{metal}$:

$$\boxed{\sqrt{P_{int}} = \frac{\alpha_{intact} \cdot g(\theta) + \beta \cdot \sqrt{P_{down}}}{\alpha_{intact} + \beta}}$$

Then solve $J_{surface} = J_{oxide}$ for $\theta$ using root-finding.

**Result:** $J_{intact} = J_{metal}(\theta_{intact})$

### Step 3: Defect Paths

#### 3a. Pinhole (No Oxide Barrier)

At pinholes, gas directly contacts metal. Two options:

**Option A: Simple Sieverts' Law** (fast surface kinetics)

$$J_{pinhole} = \beta \cdot \left( \sqrt{P_{up}} - \sqrt{P_{down}} \right)$$

**Option B: With Surface Kinetics at Metal Surface**

Solve coupled system with $\alpha = \infty$ (no oxide resistance):

$$J_{surface,ph} = k_{diss,metal} \cdot P_{up} \cdot (1 - \theta_{ph})^2 - k_{recomb,metal} \cdot \theta_{ph}^2$$

$$J_{metal,ph} = \beta \cdot \left( \sqrt{P_{int,ph}} - \sqrt{P_{down}} \right)$$

Where $P_{int,ph}$ is determined by surface equilibrium at the metal.


#### 3b. Crack (Thin Oxide Layer)

Cracks have reduced oxide thickness: $L_{crack} = \gamma \cdot L_{ox}$ where $\gamma < 1$

**Modified oxide permeance:**

$$\alpha_{crack} = \frac{D_{ox} \cdot K_{ox}}{\gamma \cdot L_{ox}} = \frac{\alpha_{intact}}{\gamma}$$

Solve the same L2+L6 system with $\alpha_{crack}$:

$$\sqrt{P_{int,crack}} = \frac{\alpha_{crack} \cdot g(\theta_{crack}) + \beta \cdot \sqrt{P_{down}}}{\alpha_{crack} + \beta}$$

Solve $J_{surface} = J_{oxide,crack}$ for $\theta_{crack}$

**Result:** $J_{crack}$


#### 3c. Oxide Grain Boundary (Enhanced Diffusion)

GB regions have enhanced oxide diffusivity: $D_{ox,gb} = \delta \cdot D_{ox}$ where $\delta > 1$

**Modified oxide permeance:**

$$\alpha_{gb} = \frac{\delta \cdot D_{ox} \cdot K_{ox}}{L_{ox}} = \delta \cdot \alpha_{intact}$$

Solve L2+L6 with $\alpha_{gb}$:

$$\sqrt{P_{int,gb}} = \frac{\alpha_{gb} \cdot g(\theta_{gb}) + \beta \cdot \sqrt{P_{down}}}{\alpha_{gb} + \beta}$$

**Result:** $J_{gb}$


### Step 4: Mixed Defects (Composite)

If defect region contains multiple types:

$$J_{defect} = \frac{f_{ph} \cdot J_{pinhole} + f_{crack} \cdot J_{crack} + f_{gb} \cdot J_{gb}}{f_{ph} + f_{crack} + f_{gb}}$$

Where $f_{ph} + f_{crack} + f_{gb} = f_{defect}$

Or more directly:

$$J_{total} = f_{intact} \cdot J_{intact} + f_{ph} \cdot J_{pinhole} + f_{crack} \cdot J_{crack} + f_{gb} \cdot J_{gb}$$


### Step 5: Total System Flux

$$\boxed{J_{total} = (1 - f_{defect}) \cdot J_{intact} + f_{defect} \cdot J_{defect}}$$


### Step 6: Rate-Limiting Analysis for Each Path

For each path, calculate resistances:

| Resistance | Formula | Physical Meaning |
|------------|---------|------------------|
| **Surface** | $R_{surface} = \frac{1}{k_{diss} \cdot P_{up}}$ (if $\theta < 0.9$) | Dissociation kinetics |
| **Oxide** | $R_{oxide} = \frac{1}{\alpha} = \frac{L_{ox}}{D_{ox} \cdot K_{ox}}$ | Oxide diffusion |
| **Metal** | $R_{metal} = \frac{1}{\beta} = \frac{L_m}{D_m \cdot K_{s,m}}$ | Metal diffusion |

**Total Resistance:**

$$R_{total} = R_{surface} + R_{oxide} + R_{metal}$$

**Fractional Contributions:**

$$f_i = \frac{R_i}{R_{total}} \quad \text{for } i \in \{\text{surface, oxide, metal}\}$$

**Rate-Limiting Identification:**

| Condition | Rate-Limiting Step | Pressure Dependence |
|-----------|-------------------|---------------------|
| $f_{surface} > 0.5$ | Surface dissociation | $J \propto P$ (slope ≈ 1) |
| $f_{oxide} > 0.5$ | Oxide diffusion | $J \propto \sqrt{P}$ (slope ≈ 0.5) |
| $f_{metal} > 0.5$ | Metal diffusion | $J \propto \sqrt{P}$ (slope ≈ 0.5) |


### Step 7: Enhancement Factor

Compare to perfect oxide:

$$\boxed{\text{Enhancement Factor} = \frac{J_{total}}{J_{intact}}}$$

**Physical Interpretation:**
- Enhancement = 1: No effect of defects
- Enhancement > 1: Defects increase permeation
- Enhancement >> 1: Defect-dominated transport


### Summary Table: Path Characteristics

| Path | Oxide Permeance | Key Parameter | Physical Picture |
|------|-----------------|---------------|------------------|
| **Intact** | $\alpha_{intact} = \frac{D_{ox} K_{ox}}{L_{ox}}$ | — | Full oxide barrier |
| **Crack** | $\alpha_{crack} = \frac{\alpha_{intact}}{\gamma}$ | $\gamma < 1$ | Thinner oxide |
| **GB** | $\alpha_{gb} = \delta \cdot \alpha_{intact}$ | $\delta > 1$ | Fast diffusion path |
| **Pinhole** | $\alpha \to \infty$ | — | No oxide barrier |
<!-- #endregion -->

# Defective Oxide + Perfect Metal (L3 +L6)


### STEP 1: single path flux calculation (L3 + L6)

```python
# =============================================================================
# STEP 1: Parallel Path Framework (L3 + L6) # For single path flux calculation, we can reuse the L2+L6 model but allow α to vary for each path type.
# =============================================================================
"""
Parallel Path Model Foundation
------------------------------
Based on Strehlow & Savage (1974), the total flux through a defective oxide is:

    J_total = f_intact × J_intact + f_defect × J_defect

Where:
    f_intact = (1 - f_defect) = area fraction of intact oxide
    f_defect = area fraction of defects (pinholes, cracks, grain boundaries)
    
Each path uses the L2+L6 coupled model with modified oxide permeance:
    - Intact: α_intact = D_ox × K_ox / L_ox  (full oxide barrier)
    - Pinhole: α → ∞ (no oxide, direct metal exposure)
    - Crack: α_crack = D_ox × K_ox / (γ × L_ox) where γ < 1 (thinner oxide)
    - GB: α_gb = δ × D_ox × K_ox / L_ox where δ > 1 (enhanced diffusivity)

The key insight: Each parallel path still satisfies:
    J_surface = J_oxide = J_metal (steady-state for that path)
    
But they have different α values, leading to different θ and P_int for each path.
"""
# For single path flux calculation, we can reuse the L2+L6 model but allow α to vary for each path type.
def calculate_path_flux_L6(
    P_up, P_down, L_m,
    k_diss, K_eq,
    alpha,  # Effective oxide permeance for this path (use np.inf for pinhole)
    D_m, K_s_m,
    path_type='intact'
):
    """
    Calculate steady-state flux through a single path (intact or defect).
    
    This is the core building block for the parallel path model.
    Each path has its own oxide permeance α, leading to different θ and P_int.
    
    For pinhole (α → ∞): Uses analytical limit where √P_int = g(θ) directly.
    
    Parameters
    ----------
    P_up : float
        Upstream H2 pressure [Pa]
    P_down : float
        Downstream H2 pressure [Pa]
    L_m : float
        Metal thickness [m]
    k_diss : float
        Surface dissociation rate constant [mol/m²/s/Pa]
    K_eq : float
        Surface equilibrium constant [Pa⁻¹]
    alpha : float or np.inf
        Oxide permeance for this path [mol/m²/s/Pa^0.5]
        - Intact: D_ox × K_ox / L_ox
        - Crack: D_ox × K_ox / (γ × L_ox) where γ < 1
        - GB: δ × D_ox × K_ox / L_ox where δ > 1
        - Pinhole: np.inf (no oxide barrier)
    D_m : float
        Metal diffusivity [m²/s]
    K_s_m : float
        Metal Sieverts constant/Solubility [mol/m³/Pa^0.5]
    path_type : str
        'intact', 'pinhole', 'crack', or 'grain_boundary' (for diagnostics)
    
    Returns
    -------
    dict
        Contains flux, theta, P_int, resistances, and rate-limiting info
    """
    
    # Metal permeance
    beta = D_m * K_s_m / L_m
    sqrt_P_down = np.sqrt(max(P_down, 0))
    
    # Check if pinhole (no oxide) - handle α → ∞ limit analytically
    is_pinhole = (path_type == 'pinhole' or alpha == np.inf or alpha > 1e10)
    
    # =========================================================================
    # Define residual function that handles both cases (pinhole = no oxide resistance and other cases with finite α and oxide resistance)
    # =========================================================================
    def residual(theta):
        if theta <= 0 or theta >= 1:
            return np.inf
        
        # Surface flux (same for both cases)
        J_surf = surface_flux(theta, P_up, k_diss, K_eq)
        
        # √P_int depends on whether we have oxide
        if is_pinhole:
            # Limit as α → ∞: √P_int = g(θ)
            sqrt_P_int = g_theta(theta, K_eq)
        else:
            # Standard case with finite α
            sqrt_P_int = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
        
        # Metal flux
        J_metal = beta * (sqrt_P_int - sqrt_P_down)
        #J_metal = metal_flux(theta_ss, alpha, beta, K_eq, P_down)
        
        return J_surf - J_metal
    
    # =========================================================================
    # Solve for θ
    # =========================================================================
    try:
        theta_ss = brentq(residual, 1e-10, 1.0 - 1e-10)
    except ValueError as e:
        return {
            'flux': np.nan,
            'theta': np.nan,
            'P_int': np.nan,
            'path_type': path_type,
            'alpha': np.inf if is_pinhole else alpha,
            'beta': beta,
            'error': str(e)
        }
    
    # =========================================================================
    # Calculate results based on path type
    # =========================================================================
    if is_pinhole:
        # Pinhole: √P_int = g(θ), no oxide resistance
        sqrt_P_int = g_theta(theta_ss, K_eq)
        R_oxide = 0.0
        J_ox = np.nan  # No oxide flux for pinhole
        alpha_out = np.inf
    else:
        # Standard case with finite α
        sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
        R_oxide = 1.0 / alpha
        J_ox = oxide_flux(theta_ss, alpha, beta, K_eq, P_down)
        alpha_out = alpha
    
    P_int = sqrt_P_int**2
    J_ss = beta * (sqrt_P_int - sqrt_P_down)
    J_surf = surface_flux(theta_ss, P_up, k_diss, K_eq)
    
    # =========================================================================
    # Rate-limiting analysis
    # =========================================================================
    R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    R_metal = 1.0 / beta
    R_total = R_surface + R_oxide + R_metal
    
    f_surface = R_surface / R_total if R_total > 0 else 0
    f_oxide = R_oxide / R_total if R_total > 0 else 0
    f_metal = R_metal / R_total if R_total > 0 else 0
    
    # Determine rate-limiting step with 'mixed' case
    if f_surface > 0.5:
        rate_limiting = 'surface'
    elif f_oxide > 0.5:
        rate_limiting = 'oxide'
    elif f_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'
    
    # Flux balance check
    if is_pinhole:
        flux_balanced = np.allclose(J_surf, J_ss, rtol=1e-6)
    else:
        flux_balanced = np.allclose(J_surf, J_ox, rtol=1e-6) and np.allclose(J_ox, J_ss, rtol=1e-6)
    
    return {
        'flux': J_ss,
        'theta': theta_ss,
        'P_int': P_int,
        'path_type': path_type,
        'alpha': alpha_out,
        'beta': beta,
        'flux_balance': {
            'J_surface': J_surf,
            'J_oxide': J_ox,
            'J_metal': J_ss,
            'balanced': flux_balanced 
        },
        'resistances': {
            'R_surface': R_surface,
            'R_oxide': R_oxide,
            'R_metal': R_metal,
            'R_total': R_total,
            'f_surface': f_surface,
            'f_oxide': f_oxide,
            'f_metal': f_metal,
        },
        'rate_limiting': rate_limiting
    }


# For Parallel Path Model: Combine intact and defect paths using area fractions
def calculate_parallel_path_flux_L6(
    P_up, P_down, L_m,
    k_diss, K_eq,
    D_ox, K_ox, L_ox,
    D_m, K_s_m,
    defect_area_fraction,
    defect_type='crack',
    thickness_factor=0.1,      # γ for cracks (L_defect = γ × L_ox)
    diffusivity_factor=100.0,  # δ for GB (D_defect = δ × D_ox)
):
    """
    Calculate total flux through parallel intact + defect paths.
    
    Implements Strehlow & Savage (1974) parallel path model with L6 surface kinetics.
    
    J_total = (1 - f_defect) × J_intact + f_defect × J_defect
    
    Parameters
    ----------
    P_up, P_down : float
        Upstream and downstream H2 pressures [Pa]
    L_m : float
        Metal thickness [m]
    k_diss, K_eq : float
        Surface kinetics parameters
    D_ox, K_ox, L_ox : float
        Intact oxide properties
    D_m, K_s_m : float
        Metal properties
    defect_area_fraction : float
        f_defect: fraction of surface area with defects (0 to 1)
    defect_type : str
        'pinhole' : No oxide (α → ∞)
        'crack'   : Thin oxide (α_crack = α_intact / γ)
        'grain_boundary' : Enhanced diffusion (α_gb = δ × α_intact)
    thickness_factor : float
        γ: For cracks, L_defect = γ × L_ox (default 0.1 = 10× thinner)
    diffusivity_factor : float
        δ: For GB, D_defect = δ × D_ox (default 100 = 100× faster)
    
    Returns
    -------
    dict
        Total flux, individual path results, and enhancement factor
    """
    
    # Validate inputs
    if not 0 <= defect_area_fraction <= 1:
        raise ValueError(f"defect_area_fraction must be between 0 and 1, got {defect_area_fraction}")
    
    # Calculate intact oxide permeance
    alpha_intact = D_ox * K_ox / L_ox
    
    # Calculate defect permeance based on type
    if defect_type == 'pinhole':
        # No oxide barrier - use np.inf to trigger analytical limit
        alpha_defect = np.inf
    elif defect_type == 'crack':
        # Thinner oxide at crack: L_defect = γ × L_ox
        # α_crack = D_ox × K_ox / (γ × L_ox) = α_intact / γ
        alpha_defect = alpha_intact / thickness_factor
    elif defect_type == 'grain_boundary':
        # Enhanced diffusion at GB: D_defect = δ × D_ox
        # α_gb = δ × D_ox × K_ox / L_ox = δ × α_intact
        alpha_defect = diffusivity_factor * alpha_intact
    else:
        raise ValueError(f"Unknown defect_type: {defect_type}. Use 'pinhole', 'crack', or 'grain_boundary'")
    
    # =========================================================================
    # Calculate flux through intact path
    # =========================================================================
    intact_result = calculate_path_flux_L6(
        P_up, P_down, L_m,
        k_diss, K_eq,
        alpha_intact,
        D_m, K_s_m,
        path_type='intact'
    )
    
    # =========================================================================
    # Calculate flux through defect path
    # =========================================================================
    defect_result = calculate_path_flux_L6(
        P_up, P_down, L_m,
        k_diss, K_eq,
        alpha_defect,
        D_m, K_s_m,
        path_type=defect_type
    )
    
    # =========================================================================
    # Area-weighted total flux
    # =========================================================================
    f_intact = 1.0 - defect_area_fraction
    f_defect = defect_area_fraction
    
    J_intact = intact_result['flux']
    J_defect = defect_result['flux']
    
    J_total = f_intact * J_intact + f_defect * J_defect
    
    # Enhancement factor compared to intact-only case
    enhancement_factor = J_total / J_intact if J_intact > 0 else np.inf
    
    # Determine which path dominates the total flux
    flux_from_intact = f_intact * J_intact
    flux_from_defect = f_defect * J_defect
    
    if flux_from_defect > 0.5 * J_total:
        dominant_path = 'defect'
    elif flux_from_intact > 0.5 * J_total:
        dominant_path = 'intact'
    else:
        dominant_path = 'mixed'
    
    return {
        # Total system results
        'J_total': J_total,
        'enhancement_factor': enhancement_factor,
        'dominant_path': dominant_path,
        
        # Area fractions
        'f_intact': f_intact,
        'f_defect': f_defect,
        
        # Flux contributions
        'flux_from_intact': flux_from_intact,
        'flux_from_defect': flux_from_defect,
        'fraction_from_defect': flux_from_defect / J_total if J_total > 0 else 0,
        
        # Individual path results
        'intact_path': intact_result,
        'defect_path': defect_result,
        
        # Permeances for reference
        'alpha_intact': alpha_intact,
        'alpha_defect': alpha_defect,
        'alpha_ratio': alpha_defect / alpha_intact,
        
        # Input parameters (for verification)
        'defect_type': defect_type,
        'thickness_factor': thickness_factor if defect_type == 'crack' else None,
        'diffusivity_factor': diffusivity_factor if defect_type == 'grain_boundary' else None,
    }
```

### STEP 2: Multiple Defect Types (Mixed Defects) + System Rate-Limiting

```python
# =============================================================================
# STEP 2: Multiple Defect Types (Mixed Defects) + System Rate-Limiting
# =============================================================================
"""
Step 2: Mixed Defect Model
--------------------------
Real oxide films often have multiple defect types simultaneously:
- Pinholes (complete oxide failure)
- Cracks (thin oxide regions)  
- Grain boundaries (enhanced diffusion paths)

The total flux becomes:

    J_total = f_intact × J_intact + Σᵢ (f_defect,i × J_defect,i)

Where the sum is over all defect types i.

Constraint: f_intact + Σᵢ f_defect,i = 1

This step also adds:
1. System-level rate-limiting analysis (which path dominates?)
2. Effective system resistance
"""

def calculate_mixed_defect_flux_L6(
    P_up, P_down, L_m,
    k_diss, K_eq,
    D_ox, K_ox, L_ox,
    D_m, K_s_m,
    defect_config,
):
    """
    Calculate flux through oxide with multiple defect types.
    
    Parameters
    ----------
    P_up, P_down : float
        Upstream and downstream H2 pressures [Pa]
    L_m : float
        Metal thickness [m]
    k_diss, K_eq : float
        Surface kinetics parameters
    D_ox, K_ox, L_ox : float
        Intact oxide properties
    D_m, K_s_m : float
        Metal properties
    defect_config : dict
        Configuration for each defect type. Example:
        {
            'pinhole': {
                'area_fraction': 0.001,  # 0.1% pinholes
            },
            'crack': {
                'area_fraction': 0.005,  # 0.5% cracks
                'thickness_factor': 0.1, # γ = 0.1 (10× thinner)
            },
            'grain_boundary': {
                'area_fraction': 0.02,   # 2% GB area
                'diffusivity_factor': 100, # δ = 100 (100× faster)
            }
        }
    
    Returns
    -------
    dict
        Total flux, individual path contributions, system analysis
    """
    
    # =========================================================================
    # Validate defect configuration
    # =========================================================================
    valid_defect_types = ['pinhole', 'crack', 'grain_boundary']
    
    total_defect_fraction = 0.0
    for defect_type, config in defect_config.items():
        if defect_type not in valid_defect_types:
            raise ValueError(f"Unknown defect type: {defect_type}. Valid types: {valid_defect_types}")
        total_defect_fraction += config.get('area_fraction', 0.0)
    
    if total_defect_fraction > 1.0:
        raise ValueError(f"Total defect area fraction ({total_defect_fraction:.3f}) exceeds 1.0")
    
    f_intact = 1.0 - total_defect_fraction
    
    # =========================================================================
    # Calculate intact oxide permeance
    # =========================================================================
    alpha_intact = D_ox * K_ox / L_ox
    
    # =========================================================================
    # Calculate flux through intact path
    # =========================================================================
    intact_result = calculate_path_flux_L6(
        P_up, P_down, L_m,
        k_diss, K_eq,
        alpha_intact,
        D_m, K_s_m,
        path_type='intact'
    )
    
    # =========================================================================
    # Calculate flux through each defect path
    # =========================================================================
    defect_results = {}
    
    for defect_type, config in defect_config.items():
        f_defect = config.get('area_fraction', 0.0)
        
        if f_defect <= 0:
            continue
        
        # Calculate defect permeance
        if defect_type == 'pinhole':
            alpha_defect = np.inf  # No oxide barrier
        elif defect_type == 'crack':
            gamma = config.get('thickness_factor', 0.1)
            alpha_defect = alpha_intact / gamma
        elif defect_type == 'grain_boundary':
            delta = config.get('diffusivity_factor', 100.0)
            alpha_defect = delta * alpha_intact
        
        # Calculate flux through this defect type
        path_result = calculate_path_flux_L6(
            P_up, P_down, L_m,
            k_diss, K_eq,
            alpha_defect,
            D_m, K_s_m,
            path_type=defect_type
        )
        
        defect_results[defect_type] = {
            'area_fraction': f_defect,
            'alpha': alpha_defect,
            'alpha_ratio': alpha_defect / alpha_intact,
            'path_result': path_result,
            'flux_contribution': f_defect * path_result['flux'],
        }
    
    # =========================================================================
    # Calculate total flux (area-weighted sum)
    # =========================================================================
    J_intact_contribution = f_intact * intact_result['flux']
    J_total = J_intact_contribution
    
    for defect_type, data in defect_results.items():
        J_total += data['flux_contribution']
    
    # =========================================================================
    # Flux breakdown analysis
    # =========================================================================
    flux_breakdown = {
        'intact': {
            'area_fraction': f_intact,
            'flux': intact_result['flux'],
            'contribution': J_intact_contribution,
            'fraction_of_total': J_intact_contribution / J_total if J_total > 0 else 0,
        }
    }
    
    for defect_type, data in defect_results.items():
        flux_breakdown[defect_type] = {
            'area_fraction': data['area_fraction'],
            'flux': data['path_result']['flux'],
            'contribution': data['flux_contribution'],
            'fraction_of_total': data['flux_contribution'] / J_total if J_total > 0 else 0,
        }
    
    # =========================================================================
    # Determine dominant path
    # =========================================================================
    max_contribution = J_intact_contribution
    dominant_path = 'intact'
    
    for defect_type, data in defect_results.items():
        if data['flux_contribution'] > max_contribution:
            max_contribution = data['flux_contribution']
            dominant_path = defect_type
    
    # Check if mixed (no single path > 70%)
    dominant_fraction = max_contribution / J_total if J_total > 0 else 0
    if dominant_fraction < 0.7:
        dominant_path = 'mixed'
    
    # =========================================================================
    # System-level rate-limiting analysis
    # =========================================================================
    # For the dominant path, what's the rate-limiting step?
    if dominant_path == 'intact':
        dominant_rate_limiting = intact_result['rate_limiting']
        dominant_resistances = intact_result['resistances']
    elif dominant_path in defect_results:
        dominant_rate_limiting = defect_results[dominant_path]['path_result']['rate_limiting']
        dominant_resistances = defect_results[dominant_path]['path_result']['resistances']
    else:
        # Mixed - report the path with highest flux contribution
        dominant_rate_limiting = 'mixed (multiple paths)'
        dominant_resistances = None
    
    # Enhancement factor vs intact-only
    enhancement_factor = J_total / intact_result['flux'] if intact_result['flux'] > 0 else np.inf
    
    # =========================================================================
    # Return comprehensive results
    # =========================================================================
    return {
        # Primary outputs
        'J_total': J_total,
        'enhancement_factor': enhancement_factor,
        'dominant_path': dominant_path,
        'dominant_fraction': dominant_fraction,
        
        # Flux breakdown
        'flux_breakdown': flux_breakdown,
        
        # Area fractions
        'f_intact': f_intact,
        'total_defect_fraction': total_defect_fraction,
        
        # Individual path results
        'intact_path': intact_result,
        'defect_paths': defect_results,
        
        # System rate-limiting
        'system_rate_limiting': dominant_rate_limiting,
        'dominant_resistances': dominant_resistances,
        
        # Reference values
        'alpha_intact': alpha_intact,
        
        # Input config (for verification)
        'defect_config': defect_config,
    }
```

```python
# =============================================================================
# Interactive Widget for L3 + L6 Model (Defective Oxide + Perfect Metal)
# Single-pass loop: plot arrays and DataFrame rows built together
# =============================================================================

import pandas as pd
from itertools import groupby
from operator import itemgetter

@interact(
    # Operating conditions
    P_down=widgets.FloatLogSlider(value=1e0, base=10, min=-2, max=4, step=0.5, description='P_down (Pa)'),
    L_m=widgets.FloatLogSlider(value=1e-3, base=10, min=-4, max=-1, step=0.5, description='L_m (m)'),
    # Surface kinetics
    k_diss=widgets.FloatLogSlider(value=1e-18, base=10, min=-18, max=-10, step=0.5, description='k_diss'),
    K_eq=widgets.FloatLogSlider(value=1e-5, base=10, min=-15, max=-5, step=0.5, description='K_eq'),
    # Oxide properties
    D_ox=widgets.FloatLogSlider(value=1e-11, base=10, min=-22, max=-9, step=0.5, description='D_ox (m²/s)'),
    K_ox=widgets.FloatLogSlider(value=1e-10, base=10, min=-12, max=0, step=0.5, description='K_ox'),
    L_ox=widgets.FloatLogSlider(value=1e-7, base=10, min=-8, max=-4, step=0.5, description='L_ox (m)'),
    # Metal properties
    D_m=widgets.FloatLogSlider(value=1e-12, base=10, min=-12, max=-6, step=0.5, description='D_m (m²/s)'),
    K_s_m=widgets.FloatLogSlider(value=0.000316, base=10, min=-6, max=0, step=0.5, description='K_s_m'),
    # Defect parameters
    f_pinhole=widgets.FloatSlider(value=0.1, min=0, max=5, step=0.1, description='Pinhole %'),
    f_crack=widgets.FloatSlider(value=0.5, min=0, max=10, step=0.5, description='Crack %'),
    f_gb=widgets.FloatSlider(value=0.2, min=0, max=20, step=1.0, description='GB %'),
    gamma=widgets.FloatLogSlider(value=0.1, base=10, min=-2, max=0, step=0.5, description='γ (crack)'),
    delta=widgets.FloatLogSlider(value=100, base=10, min=0, max=4, step=0.5, description='δ (GB)'),
)
def interactive_L36_solver(P_down, L_m, k_diss, K_eq,
                           D_ox, K_ox, L_ox, D_m, K_s_m,
                           f_pinhole, f_crack, f_gb, gamma, delta):
    """
    Interactive solver for L3+L6 parallel path model.
    Single-pass loop: plot arrays and DataFrame rows built together from the
    same computed values, so Rate-Lim is identical in both by construction.
    """
    
    # =========================================================================
    # Setup
    # =========================================================================
    # Check total defect fraction
    total_defect = f_pinhole + f_crack + f_gb
    if total_defect > 1.0:
        print(f"⚠️ Total defect fraction ({total_defect*100:.1f}%) exceeds 100%!")
        return
    
    # Build defect config
    defect_config = {}
    if f_pinhole > 0:
        defect_config['pinhole'] = {'area_fraction': f_pinhole}
    if f_crack > 0:
        defect_config['crack'] = {'area_fraction': f_crack, 'thickness_factor': gamma}
    if f_gb > 0:
        defect_config['grain_boundary'] = {'area_fraction': f_gb, 'diffusivity_factor': delta}
    
    # Pressure range for sweep
    P_up_range = np.logspace(2, 12, 40)
    alpha_intact = D_ox * K_ox / L_ox
    
    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    plot_data = {k: [] for k in [
        'J_total', 'J_intact', 'enhancement', 'theta_intact',
        'frac_intact', 'frac_pinhole', 'frac_crack', 'frac_gb',
        'f_surface', 'f_oxide', 'f_metal',
    ]}
    rows = []
    
    for P_up in P_up_range:
        try:
            if defect_config:
                r = calculate_mixed_defect_flux_L6(
                    P_up, P_down, L_m, k_diss, K_eq,
                    D_ox, K_ox, L_ox, D_m, K_s_m,
                    defect_config=defect_config
                )
                
                # Extract values ONCE
                J_tot = r['J_total']
                J_int = r['intact_path']['flux']
                theta = r['intact_path']['theta']
                P_int = r['intact_path']['P_int']
                enhancement = r['enhancement_factor']
                
                # Flux fractions
                frac_intact = r['flux_breakdown']['intact']['fraction_of_total']
                frac_pinhole = r['flux_breakdown'].get('pinhole', {}).get('fraction_of_total', 0)
                frac_crack = r['flux_breakdown'].get('crack', {}).get('fraction_of_total', 0)
                frac_gb = r['flux_breakdown'].get('grain_boundary', {}).get('fraction_of_total', 0)
                
                # Flux-weighted resistances across ALL paths (like L346 model)
                ws = wo = wm = 0.0
                J_intact_contrib = frac_intact * J_tot
                
                # Intact path contribution
                res_intact = r['intact_path']['resistances']
                w = J_intact_contrib / J_tot if J_tot > 0 else 0.0
                ws += w * res_intact['f_surface']
                wo += w * res_intact['f_oxide']
                wm += w * res_intact['f_metal']
                
                # Defect path contributions
                for dt, data in r['defect_paths'].items():
                    pr = data['path_result']
                    J_c = data['flux_contribution']
                    w = J_c / J_tot if J_tot > 0 else 0.0
                    ws += w * pr['resistances']['f_surface']
                    wo += w * pr['resistances']['f_oxide']
                    wm += w * pr['resistances']['f_metal']
                
                # Rate-limiting — single evaluation shared by plot and DataFrame
                if   ws > 0.5: rate_lim = 'surface'
                elif wo > 0.5: rate_lim = 'oxide'
                elif wm > 0.5: rate_lim = 'metal'
                else:          rate_lim = 'mixed'
                
                # Append to plot arrays
                plot_data['J_total'].append(J_tot)
                plot_data['J_intact'].append(J_int)
                plot_data['enhancement'].append(enhancement)
                plot_data['theta_intact'].append(theta)
                plot_data['frac_intact'].append(frac_intact)
                plot_data['frac_pinhole'].append(frac_pinhole)
                plot_data['frac_crack'].append(frac_crack)
                plot_data['frac_gb'].append(frac_gb)
                plot_data['f_surface'].append(ws)
                plot_data['f_oxide'].append(wo)
                plot_data['f_metal'].append(wm)
                
                # Append to DataFrame rows (same values, no recomputation)
                rows.append({
                    "P_up (Pa)":           P_up,
                    "J_total (mol/m²/s)":  J_tot,
                    "J_intact (mol/m²/s)": J_int,
                    "Enhancement":         enhancement,
                    "θ_intact":            theta,
                    "P_int_intact (Pa)":   P_int,
                    "f_intact (%)":        frac_intact * 100,
                    "f_pinhole (%)":       frac_pinhole * 100,
                    "f_crack (%)":         frac_crack * 100,
                    "f_gb (%)":            frac_gb * 100,
                    "Dominant Path":       r["dominant_path"].upper(),
                    "Rate-Limiting":       rate_lim.upper(),
                    "α_intact":            r["alpha_intact"],
                })
                
            else:
                # No defects - intact only
                r = calculate_path_flux_L6(
                    P_up, P_down, L_m, k_diss, K_eq,
                    alpha_intact, D_m, K_s_m, path_type='intact'
                )
                
                # Extract values ONCE
                J_ss = r['flux']
                theta = r['theta']
                P_int = r['P_int']
                ws = r['resistances']['f_surface']
                wo = r['resistances']['f_oxide']
                wm = r['resistances']['f_metal']
                
                # Rate-limiting — single evaluation
                if   ws > 0.5: rate_lim = 'surface'
                elif wo > 0.5: rate_lim = 'oxide'
                elif wm > 0.5: rate_lim = 'metal'
                else:          rate_lim = 'mixed'
                
                # Append to plot arrays
                plot_data['J_total'].append(J_ss)
                plot_data['J_intact'].append(J_ss)
                plot_data['enhancement'].append(1.0)
                plot_data['theta_intact'].append(theta)
                plot_data['frac_intact'].append(1.0)
                plot_data['frac_pinhole'].append(0)
                plot_data['frac_crack'].append(0)
                plot_data['frac_gb'].append(0)
                plot_data['f_surface'].append(ws)
                plot_data['f_oxide'].append(wo)
                plot_data['f_metal'].append(wm)
                
                # Append to DataFrame rows
                rows.append({
                    "P_up (Pa)":           P_up,
                    "J_total (mol/m²/s)":  J_ss,
                    "J_intact (mol/m²/s)": J_ss,
                    "Enhancement":         1.0,
                    "θ_intact":            theta,
                    "P_int_intact (Pa)":   P_int,
                    "f_intact (%)":        100.0,
                    "f_pinhole (%)":       0.0,
                    "f_crack (%)":         0.0,
                    "f_gb (%)":            0.0,
                    "Dominant Path":       "INTACT",
                    "Rate-Limiting":       rate_lim.upper(),
                    "α_intact":            r["alpha"],
                })
                
        except Exception as e:
            for k in plot_data:
                plot_data[k].append(np.nan)
            rows.append({"P_up (Pa)": P_up, "Rate-Limiting": "ERROR", "Error": str(e)})
    
    # =========================================================================
    # Convert to arrays
    # =========================================================================
    J_total = np.array(plot_data['J_total'])
    J_intact = np.array(plot_data['J_intact'])
    f_surface = np.array(plot_data['f_surface'])
    f_oxide = np.array(plot_data['f_oxide'])
    f_metal = np.array(plot_data['f_metal'])
    frac_intact = np.array(plot_data['frac_intact'])
    frac_pinhole = np.array(plot_data['frac_pinhole'])
    frac_crack = np.array(plot_data['frac_crack'])
    frac_gb_arr = np.array(plot_data['frac_gb'])
    
    # Vectorised rate_limiting_arr for plot shading — derived from the same
    # f_surface/f_oxide/f_metal arrays written inside the loop
    rate_limiting_arr = np.where(
        f_surface > 0.5, 'surface',
        np.where(f_oxide > 0.5, 'oxide',
        np.where(f_metal > 0.5, 'metal',
                 'mixed'))
    )
    
    # =========================================================================
    # Create figure with 2 subplots
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # ===== Plot 1: Flux vs Pressure =====
    valid_idx = ~np.isnan(J_total)
    ax1.loglog(P_up_range, J_total, 'k-', linewidth=2.5, label='L3+L6 Model (Total)')
    ax1.loglog(P_up_range, J_intact, 'b--', linewidth=1.5, alpha=0.7, label='Intact only')
    
    if np.any(valid_idx):
        # Slope = 1 reference at the beginning
        P_ref1 = P_up_range[0]
        J_ref1 = J_total[valid_idx][0]
        J_slope1 = J_ref1 * (P_up_range / P_ref1) ** 1.0
        ax1.loglog(P_up_range, J_slope1, 'r--', linewidth=1.5, alpha=0.5, label='Slope = 1 (surface)')

        # Slope = 0.5 reference at the end
        P_ref05 = P_up_range[-1]
        J_ref05 = J_total[valid_idx][-1]
        J_slope05 = J_ref05 * (P_up_range / P_ref05) ** 0.5
        ax1.loglog(P_up_range, J_slope05, 'g--', linewidth=1.5, alpha=0.5, label='Slope = 0.5 (diffusion)')
    
    # Rate-limiting regions with 4-branch classification
    regions = [
        {'mask': rate_limiting_arr == 'surface', 'color': 'red',    'label': 'Surface-limited'},
        {'mask': rate_limiting_arr == 'oxide',   'color': 'orange', 'label': 'Oxide-limited'},
        {'mask': rate_limiting_arr == 'metal',   'color': 'blue',   'label': 'Metal-limited'},
        {'mask': rate_limiting_arr == 'mixed',   'color': 'green',  'label': 'Mixed'},
    ]

    for region in regions:
        mask = region['mask'] & valid_idx
        if np.any(mask):
            idxs = np.where(mask)[0]
            for k, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
                group = list(map(itemgetter(1), g))
                if len(group) > 2:
                    P_seg = P_up_range[group]
                    J_seg = J_total[group]
                    ax1.loglog(P_seg, J_seg, color=region['color'], linewidth=4, alpha=0.7)
                    slope, _ = np.polyfit(np.log10(P_seg), np.log10(np.abs(J_seg)), 1)
                    mid = len(group) // 2
                    ax1.text(P_seg[mid], J_seg[mid], f"{region['label']}\nSlope={slope:.2f}",
                             color=region['color'], fontsize=10, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L3+L6: Flux vs Pressure (Defective Oxide)', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    # ===== Plot 2: Rate-Limiting Fractions =====
    ax2.semilogx(P_up_range, f_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, f_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, f_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
    ax2.axhline(50, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
    
    ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax2.set_ylabel('Resistance Fraction (%)', fontsize=12)
    ax2.set_title('Rate-Limiting Step Analysis', fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')

    plt.tight_layout()
    plt.show()
    
    # =========================================================================
    # DataFrame — already built inside the loop, nothing to recompute
    # =========================================================================
    df = pd.DataFrame(rows)
    display(df)
```




# Defective Oxide + Defective Metal + Surface Chemistry

```python
# =============================================================================
# FULL MODEL: L3 + L4 + L6
# Defective Oxide + Defective Metal + Surface Chemistry
# =============================================================================
"""
Complete Hydrogen Permeation Model
===================================

This model combines three physics layers:

L6 - Surface Chemistry:
    - Dissociation/recombination kinetics at gas-oxide interface
    - Surface coverage θ determines boundary condition
    - J_surface = k_diss × P × (1-θ)² - k_recomb × θ²

L3 - Defective Oxide (Parallel Paths):
    - Intact oxide: α_intact = D_ox × K_ox / L_ox
    - Pinholes: α → ∞ (no oxide barrier)
    - Cracks: α_crack = α_intact / γ (thinner oxide)
    - Grain boundaries: α_gb = δ × α_intact (enhanced diffusion)
    
L4 - Defective Metal (Microstructure):
    - Grain boundary enhancement: D_eff = f_gb × D_lattice
    - Trapping reduction: D_eff = D_lattice / (1 + θ_trap)
    - Position-dependent effects through concentration profile

Total flux = Σᵢ (fᵢ × Jᵢ) where i = {intact, pinhole, crack, grain_boundary}
"""

def calculate_path_flux_L346_v2(
    P_up, P_down, L_m, temperature,
    k_diss, K_eq,
    alpha,  # Oxide permeance for this path (np.inf for pinhole)
    D_lattice, K_s_m,
    microstructure_params,
    path_type='intact',
    lattice_density=1.06e29,
    method='average',
    n_points=20,
    mode='both',
    max_iterations=10,
    tolerance=1e-5
):
    """
    Calculate flux through a single oxide path with defective metal.
    
    This is the core building block that combines:
    - L6: Surface kinetics (θ-based)
    - Variable α for oxide path type
    - L4: Microstructure effects in metal
    
    Parameters
    ----------
    P_up, P_down : float
        Upstream and downstream H2 pressures [Pa]
    L_m : float
        Metal thickness [m]
    temperature : float
        Temperature [K]
    k_diss : float
        Surface dissociation rate constant [mol/m²/s/Pa]
    K_eq : float
        Surface equilibrium constant [Pa⁻¹]
    alpha : float or np.inf
        Oxide permeance for this path [mol/m²/s/Pa^0.5]
    D_lattice : float
        Metal lattice diffusivity [m²/s]
    K_s_m : float
        Metal Sieverts constant [mol/m³/Pa^0.5]
    microstructure_params : dict
        Metal microstructure specification
    path_type : str
        'intact', 'pinhole', 'crack', or 'grain_boundary'
    
    Returns
    -------
    dict
        Flux, theta, P_int, D_eff, resistances, rate-limiting info
    """
    from calculations.defective_metal import combined_microstructure_model
    
    # =========================================================================
    # Handle pinhole case (α → ∞)
    # =========================================================================
    is_pinhole = (path_type == 'pinhole' or alpha == np.inf or alpha > 1e10)
    sqrt_P_down = np.sqrt(max(P_down, 0))
    
    # =========================================================================
    # Iterative solution: D_eff and θ are coupled
    # =========================================================================
    D_eff = D_lattice  # Initial guess
    convergence_history = []
    
    for iteration in range(max_iterations):
        beta = D_eff * K_s_m / L_m
        
        # ---------------------------------------------------------------------
        # Define residual for θ (surface flux = metal flux)
        # ---------------------------------------------------------------------
        def residual(theta):
            if theta <= 0 or theta >= 1:
                return np.inf
            
            J_surf = surface_flux(theta, P_up, k_diss, K_eq)
            
            if is_pinhole:
                # No oxide: √P_int = g(θ) directly from surface equilibrium
                sqrt_P_int = g_theta(theta, K_eq)
            else:
                # With oxide: solve flux balance
                sqrt_P_int = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
            
            J_metal = beta * (sqrt_P_int - sqrt_P_down)
            return J_surf - J_metal
        
        # ---------------------------------------------------------------------
        # Solve for θ
        # ---------------------------------------------------------------------
        try:
            theta_ss = brentq(residual, 1e-12, 1.0 - 1e-12)
        except ValueError as e:
            return {
                'flux': np.nan,
                'theta': np.nan,
                'P_int': np.nan,
                'D_eff': np.nan,
                'path_type': path_type,
                'error': f'Failed to solve for theta: {str(e)}',
                'iteration': iteration
            }
        
        # ---------------------------------------------------------------------
        # Calculate P_int from θ
        # ---------------------------------------------------------------------
        if is_pinhole:
            sqrt_P_int = g_theta(theta_ss, K_eq)
        else:
            sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
        
        # ---------------------------------------------------------------------
        # Calculate D_eff through metal using microstructure model
        # ---------------------------------------------------------------------
        x_array = np.linspace(0, L_m, n_points)
        sqrt_P_array = sqrt_P_int - (sqrt_P_int - sqrt_P_down) * x_array / L_m
        C_array = K_s_m * sqrt_P_array
        
        D_array = np.zeros(n_points)
        theta_trap_array = np.zeros(n_points)
        gb_factor_array = np.zeros(n_points)
        
        for i, C_local in enumerate(C_array):
            C_local = max(C_local, 1e-20)
            
            result_i = combined_microstructure_model(
                D_lattice=D_lattice,
                temperature=temperature,
                microstructure_params=microstructure_params,
                lattice_concentration=C_local,
                lattice_density=lattice_density,
                mode=mode
            )
            
            D_array[i] = result_i['D_eff']
            
            # Extract trapping info
            if 'trapping' in result_i and result_i['trapping'] is not None:
                theta_trap_array[i] = result_i['trapping'].get('theta_total', 0.0)
            elif 'theta_total' in result_i:
                theta_trap_array[i] = result_i['theta_total']
            
            # Extract GB enhancement
            if 'gb_enhancement' in result_i and result_i['gb_enhancement'] is not None:
                gb_factor_array[i] = result_i['gb_enhancement'].get('factor', 1.0)
            elif 'gb_enhancement_factor' in result_i:
                gb_factor_array[i] = result_i['gb_enhancement_factor']
            else:
                gb_factor_array[i] = 1.0
        
        # Average D_eff
        if method == 'average':
            D_eff_new = np.mean(D_array)
        elif method == 'harmonic':
            D_eff_new = len(D_array) / np.sum(1.0 / D_array)
        elif method == 'inlet':
            D_eff_new = D_array[0]
        elif method == 'outlet':
            D_eff_new = D_array[-1]
        else:
            D_eff_new = np.mean(D_array)
        
        # Track convergence
        rel_change = abs(D_eff_new - D_eff) / D_eff if D_eff > 0 else np.inf
        convergence_history.append({
            'iteration': iteration,
            'D_eff': D_eff_new,
            'theta': theta_ss,
            'relative_change': rel_change
        })
        
        # Check convergence
        if rel_change < tolerance:
            D_eff = D_eff_new
            break
        
        D_eff = D_eff_new
    
    # =========================================================================
    # Final calculations with converged values
    # =========================================================================
    beta_eff = D_eff * K_s_m / L_m
    P_int = sqrt_P_int**2
    
    # Fluxes
    J_metal = beta_eff * (sqrt_P_int - sqrt_P_down)
    J_surface = surface_flux(theta_ss, P_up, k_diss, K_eq)
    
    if is_pinhole:
        J_oxide = np.nan  # No oxide
        R_oxide = 0.0
        alpha_out = np.inf
    else:
        J_oxide = oxide_flux(theta_ss, alpha, beta_eff, K_eq, P_down)
        R_oxide = 1.0 / alpha
        alpha_out = alpha
    
    # =========================================================================
    # Rate-limiting analysis
    # =========================================================================
    R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    R_metal = 1.0 / beta_eff
    R_total = R_surface + R_oxide + R_metal
    
    f_surface = R_surface / R_total if R_total > 0 else 0
    f_oxide = R_oxide / R_total if R_total > 0 else 0
    f_metal = R_metal / R_total if R_total > 0 else 0
    
    if f_surface > 0.5:
        rate_limiting = 'surface'
    elif f_oxide > 0.5:
        rate_limiting = 'oxide'
    elif f_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'
    
    # Microstructure diagnostics
    modification_factor = D_eff / D_lattice
    avg_gb_factor = np.mean(gb_factor_array)
    avg_theta_trap = np.mean(theta_trap_array)
    
    return {
        'flux': J_metal,
        'theta': theta_ss,
        'P_int': P_int,
        'path_type': path_type,
        
        # Permeances
        'alpha': alpha_out,
        'beta_lattice': D_lattice * K_s_m / L_m,
        'beta_eff': beta_eff,
        
        # Microstructure
        'D_eff': D_eff,
        'D_lattice': D_lattice,
        'modification_factor': modification_factor,
        'microstructure': {
            'avg_gb_factor': avg_gb_factor,
            'avg_theta_trap': avg_theta_trap,
        },
        
        # Flux verification
        'flux_balance': {
            'J_surface': J_surface,
            'J_oxide': J_oxide,
            'J_metal': J_metal,
        },
        
        # Rate-limiting
        'resistances': {
            'R_surface': R_surface,
            'R_oxide': R_oxide,
            'R_metal': R_metal,
            'R_total': R_total,
            'f_surface': f_surface,
            'f_oxide': f_oxide,
            'f_metal': f_metal,
        },
        'rate_limiting': rate_limiting,
        
        # Convergence
        'convergence': {
            'iterations': len(convergence_history),
            'converged': len(convergence_history) < max_iterations,
            'history': convergence_history
        },
        
        # Profiles
        'profiles': {
            'x': x_array,
            'D': D_array,
            'C': C_array,
            'theta_trap': theta_trap_array,
            'gb_factor': gb_factor_array
        }
    }


def calculate_full_model_flux_L346_v2(
    P_up, P_down, L_m, temperature,
    k_diss, K_eq,
    D_ox, K_ox, L_ox,
    D_lattice, K_s_m,
    microstructure_params,
    defect_config,
    lattice_density=1.06e29,
    method='average',
    n_points=20,
    mode='both'
):
    """
    Calculate total flux with defective oxide AND defective metal.
    
    Full L3 + L4 + L6 model combining:
    - L6: Surface kinetics at gas-oxide interface
    - L3: Parallel oxide paths (intact + defects)
    - L4: Defective metal microstructure
    
    Parameters
    ----------
    P_up, P_down : float
        Upstream and downstream pressures [Pa]
    L_m : float
        Metal thickness [m]
    temperature : float
        Temperature [K]
    k_diss, K_eq : float
        Surface kinetics parameters
    D_ox, K_ox, L_ox : float
        Intact oxide properties
    D_lattice, K_s_m : float
        Metal lattice properties
    microstructure_params : dict
        Metal microstructure (grain size, traps, etc.)
    defect_config : dict
        Oxide defect configuration:
        {
            'pinhole': {'area_fraction': 0.001},
            'crack': {'area_fraction': 0.005, 'thickness_factor': 0.1},
            'grain_boundary': {'area_fraction': 0.02, 'diffusivity_factor': 100}
        }
    
    Returns
    -------
    dict
        Total flux, path contributions, microstructure effects, rate-limiting
    """
    
    # =========================================================================
    # Validate defect configuration
    # =========================================================================
    valid_defect_types = ['pinhole', 'crack', 'grain_boundary']
    
    total_defect_fraction = 0.0
    for defect_type, config in defect_config.items():
        if defect_type not in valid_defect_types:
            raise ValueError(f"Unknown defect type: {defect_type}")
        total_defect_fraction += config.get('area_fraction', 0.0)
    
    if total_defect_fraction > 1.0:
        raise ValueError(f"Total defect fraction ({total_defect_fraction:.3f}) exceeds 1.0")
    
    f_intact = 1.0 - total_defect_fraction
    
    # =========================================================================
    # Calculate intact oxide permeance
    # =========================================================================
    alpha_intact = D_ox * K_ox / L_ox
    
    # =========================================================================
    # Calculate flux through intact path
    # =========================================================================
    intact_result = calculate_path_flux_L346_v2(
        P_up, P_down, L_m, temperature,
        k_diss, K_eq,
        alpha_intact,
        D_lattice, K_s_m,
        microstructure_params,
        path_type='intact',
        lattice_density=lattice_density,
        method=method,
        n_points=n_points,
        mode=mode
    )
    
    # =========================================================================
    # Calculate flux through each defect path
    # =========================================================================
    defect_results = {}
    
    for defect_type, config in defect_config.items():
        f_defect = config.get('area_fraction', 0.0)
        
        if f_defect <= 0:
            continue
        
        # Calculate defect permeance based on type
        if defect_type == 'pinhole':
            alpha_defect = np.inf  # No oxide barrier
        elif defect_type == 'crack':
            gamma = config.get('thickness_factor', 0.1)
            alpha_defect = alpha_intact / gamma
        elif defect_type == 'grain_boundary':
            delta = config.get('diffusivity_factor', 100.0)
            alpha_defect = delta * alpha_intact
        
        # Calculate flux through this defect path
        path_result = calculate_path_flux_L346_v2(
            P_up, P_down, L_m, temperature,
            k_diss, K_eq,
            alpha_defect,
            D_lattice, K_s_m,
            microstructure_params,
            path_type=defect_type,
            lattice_density=lattice_density,
            method=method,
            n_points=n_points,
            mode=mode
        )
        
        defect_results[defect_type] = {
            'area_fraction': f_defect,
            'alpha': alpha_defect,
            'alpha_ratio': alpha_defect / alpha_intact if alpha_defect != np.inf else np.inf,
            'path_result': path_result,
            'flux_contribution': f_defect * path_result['flux'],
        }
    
    # =========================================================================
    # Calculate total flux (area-weighted sum)
    # =========================================================================
    J_intact_contribution = f_intact * intact_result['flux']
    J_total = J_intact_contribution
    
    for defect_type, data in defect_results.items():
        J_total += data['flux_contribution']
    
    # =========================================================================
    # Flux breakdown analysis
    # =========================================================================
    flux_breakdown = {
        'intact': {
            'area_fraction': f_intact,
            'flux': intact_result['flux'],
            'contribution': J_intact_contribution,
            'fraction_of_total': J_intact_contribution / J_total if J_total > 0 else 0,
            'D_eff': intact_result['D_eff'],
            'modification_factor': intact_result['modification_factor'],
            'theta': intact_result['theta'],
            'P_int': intact_result['P_int'],
        }
    }
    
    for defect_type, data in defect_results.items():
        pr = data['path_result']
        flux_breakdown[defect_type] = {
            'area_fraction': data['area_fraction'],
            'flux': pr['flux'],
            'contribution': data['flux_contribution'],
            'fraction_of_total': data['flux_contribution'] / J_total if J_total > 0 else 0,
            'D_eff': pr['D_eff'],
            'modification_factor': pr['modification_factor'],
            'theta': pr['theta'],
            'P_int': pr['P_int'],
        }
    
    # =========================================================================
    # Determine dominant path
    # =========================================================================
    max_contribution = J_intact_contribution
    dominant_path = 'intact'
    
    for defect_type, data in defect_results.items():
        if data['flux_contribution'] > max_contribution:
            max_contribution = data['flux_contribution']
            dominant_path = defect_type
    
    dominant_fraction = max_contribution / J_total if J_total > 0 else 0
    if dominant_fraction < 0.7:
        dominant_path = 'mixed'
    
    # =========================================================================
    # Enhancement factors
    # =========================================================================
    enhancement_vs_intact = J_total / intact_result['flux'] if intact_result['flux'] > 0 else np.inf
    
    # =========================================================================
    # System-level rate-limiting (from dominant path)
    # =========================================================================
    if dominant_path == 'intact':
        system_rate_limiting = intact_result['rate_limiting']
        system_resistances = intact_result['resistances']
    elif dominant_path in defect_results:
        system_rate_limiting = defect_results[dominant_path]['path_result']['rate_limiting']
        system_resistances = defect_results[dominant_path]['path_result']['resistances']
    else:
        system_rate_limiting = 'mixed'
        system_resistances = None
    
    # =========================================================================
    # Average D_eff across all paths (flux-weighted)
    # =========================================================================
    D_eff_avg = 0.0
    for path, data in flux_breakdown.items():
        D_eff_avg += data['fraction_of_total'] * data['D_eff']
    
    # =========================================================================
    # Return comprehensive results
    # =========================================================================
    return {
        # Primary outputs
        'J_total': J_total,
        'enhancement_vs_intact': enhancement_vs_intact,
        'dominant_path': dominant_path,
        'dominant_fraction': dominant_fraction,
        
        # Flux breakdown
        'flux_breakdown': flux_breakdown,
        
        # Area fractions
        'f_intact': f_intact,
        'total_defect_fraction': total_defect_fraction,
        
        # Individual path results
        'intact_path': intact_result,
        'defect_paths': defect_results,
        
        # Microstructure summary (flux-weighted average)
        'D_eff_avg': D_eff_avg,
        'D_lattice': D_lattice,
        'overall_modification_factor': D_eff_avg / D_lattice,
        
        # System analysis
        'system_rate_limiting': system_rate_limiting,
        'system_resistances': system_resistances,
        
        # Reference values
        'alpha_intact': alpha_intact,
        
        # Input configs (for verification)
        'defect_config': defect_config,
        'microstructure_params': microstructure_params,
    }
```

```python
# =============================================================================
# Interactive Widget for Full L3 + L4 + L6 Model (v2)
# Defective Oxide + Defective Metal + Surface Chemistry
# =============================================================================

import pandas as pd
from itertools import groupby
from operator import itemgetter

@interact(
    P_down      = widgets.FloatLogSlider(value=1e0,   base=10, min=-2,  max=4,  step=0.5, description='P_down (Pa)'),
    L_m         = widgets.FloatLogSlider(value=1e-3,  base=10, min=-4,  max=-1, step=0.5, description='L_m (m)'),
    temperature = widgets.IntSlider(     value=773,            min=473, max=1273, step=50, description='T (K)'),
    k_diss      = widgets.FloatLogSlider(value=1e-15, base=10, min=-18, max=-10, step=0.5, description='k_diss'),
    K_eq        = widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-5,  step=0.5, description='K_eq'),
    D_ox        = widgets.FloatLogSlider(value=1e-18, base=10, min=-22, max=-12, step=0.5, description='D_ox (m²/s)'),
    K_ox        = widgets.FloatLogSlider(value=1e-3,  base=10, min=-6,  max=0,   step=0.5, description='K_ox'),
    L_ox        = widgets.FloatLogSlider(value=1e-6,  base=10, min=-8,  max=-4,  step=0.5, description='L_ox (m)'),
    D_lattice   = widgets.FloatLogSlider(value=1e-10, base=10, min=-14, max=-6,  step=0.5, description='D_lattice'),
    K_s_m       = widgets.FloatLogSlider(value=0.03,  base=10, min=-4,  max=0,   step=0.5, description='K_s_m'),
    f_pinhole   = widgets.FloatSlider(   value=0.1,            min=0,   max=2,   step=0.05, description='Pinhole %'),
    f_crack     = widgets.FloatSlider(   value=0.5,            min=0,   max=5,   step=0.1,  description='Crack %'),
    f_gb_ox     = widgets.FloatSlider(   value=0.2,            min=0,   max=10,  step=0.5,  description='Oxide GB %'),
    gamma       = widgets.FloatLogSlider(value=0.1,   base=10, min=-2,  max=0,   step=0.5, description='γ (crack)'),
    delta       = widgets.FloatLogSlider(value=100,   base=10, min=0,   max=4,   step=0.5, description='δ (ox GB)'),
    grain_size  = widgets.FloatLogSlider(value=50e-6, base=10, min=-6,  max=-3,  step=0.5, description='Grain (m)'),
    trap_density= widgets.FloatLogSlider(value=1e14,  base=10, min=12,  max=17,  step=0.5, description='ρ_trap (m⁻²)'),
    E_trap      = widgets.FloatSlider(   value=27,             min=15,  max=60,  step=5,   description='E_trap (kJ/mol)'),
)
def interactive_L346_full_model(P_down, L_m, temperature, k_diss, K_eq,
                                 D_ox, K_ox, L_ox, D_lattice, K_s_m,
                                 f_pinhole, f_crack, f_gb_ox, gamma, delta,
                                 grain_size, trap_density, E_trap):
    """
    Interactive widget for the full L3+L4+L6 model.
    Defective Oxide + Defective Metal + Surface Chemistry.
    Single-pass loop: plot arrays and DataFrame rows built together from the
    same computed values, so Rate-Lim is identical in both by construction.
    """

    # =========================================================================
    # Setup
    # =========================================================================
    f_pinhole = f_pinhole / 100
    f_crack   = f_crack   / 100
    f_gb_ox   = f_gb_ox   / 100

    total_defect = f_pinhole + f_crack + f_gb_ox
    if total_defect > 1.0:
        print(f"⚠️ Total defect fraction ({total_defect*100:.1f}%) exceeds 100%!")
        return

    defect_config = {}
    if f_pinhole > 0:
        defect_config['pinhole'] = {'area_fraction': f_pinhole}
    if f_crack > 0:
        defect_config['crack'] = {'area_fraction': f_crack, 'thickness_factor': gamma}
    if f_gb_ox > 0:
        defect_config['grain_boundary'] = {'area_fraction': f_gb_ox, 'diffusivity_factor': delta}

    microstructure_params = {
        'grain_size':  grain_size,
        'grain_shape': 'equiaxed',
        'gb_type':     'HAGB',
        'trap_list': [
            {'name': 'dislocations', 'density': trap_density, 'binding_energy': E_trap * 1000}
        ] if trap_density > 1e12 else []
    }

    microstructure_lattice = {
        'grain_size':  1e-3,
        'grain_shape': 'equiaxed',
        'gb_type':     'LAGB',
        'trap_list':   []
    }

    P_up_range   = np.logspace(1, 8, 40)
    alpha_intact = D_ox * K_ox / L_ox

    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    plot_data = {k: [] for k in [
        'J_full', 'J_no_micro', 'J_no_defects', 'J_baseline',
        'frac_intact', 'frac_pinhole', 'frac_crack', 'frac_gb',
        'D_mod_factor', 'theta_intact',
        'f_surface', 'f_oxide', 'f_metal',
    ]}
    rows = []

    for P_up in P_up_range:
        try:
            # -----------------------------------------------------------------
            # Primary computation
            # -----------------------------------------------------------------
            if defect_config:
                r = calculate_full_model_flux_L346_v2(
                    P_up, P_down, L_m, temperature,
                    k_diss, K_eq, D_ox, K_ox, L_ox,
                    D_lattice, K_s_m,
                    microstructure_params, defect_config,
                    n_points=10
                )

                # Flux-weighted resistances across ALL paths
                J_tot = r['J_total']
                ws = wo = wm = 0.0
                J_intact_contrib = r['flux_breakdown']['intact']['fraction_of_total'] * J_tot
                for J_c, res in (
                    [(J_intact_contrib, r['intact_path']['resistances'])] +
                    [(d['flux_contribution'], d['path_result']['resistances'])
                     for d in r['defect_paths'].values()]
                ):
                    w   = J_c / J_tot if J_tot > 0 else 0.0
                    ws += w * res['f_surface']
                    wo += w * res['f_oxide']
                    wm += w * res['f_metal']

                # Rate-limiting — single evaluation shared by plot and DataFrame
                if   ws > 0.5: rate_lim = 'surface'
                elif wo > 0.5: rate_lim = 'oxide'
                elif wm > 0.5: rate_lim = 'metal'
                else:          rate_lim = 'mixed'

                # Plot scalars
                plot_data['J_full'].append(r['J_total'])
                plot_data['D_mod_factor'].append(r['overall_modification_factor'])
                plot_data['theta_intact'].append(r['intact_path']['theta'])
                plot_data['frac_intact'].append(r['flux_breakdown']['intact']['fraction_of_total'])
                plot_data['frac_pinhole'].append(r['flux_breakdown'].get('pinhole', {}).get('fraction_of_total', 0))
                plot_data['frac_crack'].append(r['flux_breakdown'].get('crack', {}).get('fraction_of_total', 0))
                plot_data['frac_gb'].append(r['flux_breakdown'].get('grain_boundary', {}).get('fraction_of_total', 0))
                plot_data['f_surface'].append(ws)
                plot_data['f_oxide'].append(wo)
                plot_data['f_metal'].append(wm)

                # DataFrame row
                rows.append({
                    "P_up (Pa)":          P_up,
                    "J_total (mol/m²/s)": r['J_total'],
                    "J_intact":           r['intact_path']['flux'],
                    "Enhancement":        r['enhancement_vs_intact'],
                    "θ_intact":           r['intact_path']['theta'],
                    "P_int (Pa)":         r['intact_path']['P_int'],
                    "D_eff/D_L":          r['overall_modification_factor'],
                    "f_intact (%)":       r['flux_breakdown']['intact']['fraction_of_total'] * 100,
                    "f_pinhole (%)":      r['flux_breakdown'].get('pinhole', {}).get('fraction_of_total', 0) * 100,
                    "f_crack (%)":        r['flux_breakdown'].get('crack', {}).get('fraction_of_total', 0) * 100,
                    "f_gb (%)":           r['flux_breakdown'].get('grain_boundary', {}).get('fraction_of_total', 0) * 100,
                    "Dominant":           r['dominant_path'].upper(),
                    "Rate-Lim":           rate_lim.upper(),
                })

                r_no_micro = calculate_full_model_flux_L346_v2(
                    P_up, P_down, L_m, temperature,
                    k_diss, K_eq, D_ox, K_ox, L_ox,
                    D_lattice, K_s_m,
                    microstructure_lattice, defect_config,
                    n_points=10
                )
                plot_data['J_no_micro'].append(r_no_micro['J_total'])

            else:
                # No defects — single path
                r = calculate_path_flux_L346_v2(
                    P_up, P_down, L_m, temperature,
                    k_diss, K_eq, alpha_intact,
                    D_lattice, K_s_m, microstructure_params,
                    n_points=10
                )

                ws = r['resistances']['f_surface']
                wo = r['resistances']['f_oxide']
                wm = r['resistances']['f_metal']

                if   ws > 0.5: rate_lim = 'surface'
                elif wo > 0.5: rate_lim = 'oxide'
                elif wm > 0.5: rate_lim = 'metal'
                else:          rate_lim = 'mixed'

                plot_data['J_full'].append(r['flux'])
                plot_data['D_mod_factor'].append(r['modification_factor'])
                plot_data['theta_intact'].append(r['theta'])
                plot_data['frac_intact'].append(1.0)
                plot_data['frac_pinhole'].append(0)
                plot_data['frac_crack'].append(0)
                plot_data['frac_gb'].append(0)
                plot_data['f_surface'].append(ws)
                plot_data['f_oxide'].append(wo)
                plot_data['f_metal'].append(wm)

                rows.append({
                    "P_up (Pa)":          P_up,
                    "J_total (mol/m²/s)": r['flux'],
                    "J_intact":           r['flux'],
                    "Enhancement":        1.0,
                    "θ_intact":           r['theta'],
                    "P_int (Pa)":         r['P_int'],
                    "D_eff/D_L":          r['modification_factor'],
                    "f_intact (%)":       100.0,
                    "f_pinhole (%)":      0.0,
                    "f_crack (%)":        0.0,
                    "f_gb (%)":           0.0,
                    "Dominant":           "INTACT",
                    "Rate-Lim":           rate_lim.upper(),
                })

                r_no_micro = calculate_path_flux_L346_v2(
                    P_up, P_down, L_m, temperature,
                    k_diss, K_eq, alpha_intact,
                    D_lattice, K_s_m, microstructure_lattice,
                    n_points=10
                )
                plot_data['J_no_micro'].append(r_no_micro['flux'])

            # Comparison baselines (shared between both branches)
            r_nd = calculate_path_flux_L346_v2(
                P_up, P_down, L_m, temperature,
                k_diss, K_eq, alpha_intact,
                D_lattice, K_s_m, microstructure_params,
                n_points=10
            )
            plot_data['J_no_defects'].append(r_nd['flux'])

            r_bl = calculate_path_flux_L346_v2(
                P_up, P_down, L_m, temperature,
                k_diss, K_eq, alpha_intact,
                D_lattice, K_s_m, microstructure_lattice,
                n_points=10
            )
            plot_data['J_baseline'].append(r_bl['flux'])

        except Exception as e:
            for k in plot_data:
                plot_data[k].append(np.nan)
            rows.append({"P_up (Pa)": P_up, "Rate-Lim": "ERROR", "Error": str(e)})

    # =========================================================================
    # Convert to arrays
    # =========================================================================
    J_full       = np.array(plot_data['J_full'])
    J_no_micro   = np.array(plot_data['J_no_micro'])
    J_no_defects = np.array(plot_data['J_no_defects'])
    J_baseline   = np.array(plot_data['J_baseline'])
    f_surface    = np.array(plot_data['f_surface'])
    f_oxide      = np.array(plot_data['f_oxide'])
    f_metal      = np.array(plot_data['f_metal'])
    frac_intact  = np.array(plot_data['frac_intact'])
    frac_pinhole = np.array(plot_data['frac_pinhole'])
    frac_crack   = np.array(plot_data['frac_crack'])
    frac_gb      = np.array(plot_data['frac_gb'])
    D_mod_factor = np.array(plot_data['D_mod_factor'])

    # Vectorised rate_limiting_arr for plot shading — derived from the same
    # f_surface/f_oxide/f_metal arrays written inside the loop, so it is
    # identical to each row's Rate-Lim by construction.
    rate_limiting_arr = np.where(
        f_surface > 0.5, 'surface',
        np.where(f_oxide > 0.5, 'oxide',
        np.where(f_metal > 0.5, 'metal',
                 'mixed'))
    )

    # =========================================================================
    # Plot 1: Flux vs Pressure with Rate-Limiting Regions (Top Left)
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax1 = axes[0, 0]
    valid_idx = ~np.isnan(J_full)

    ax1.loglog(P_up_range, J_full,     'k-',   linewidth=2.5, label='Full L3+L4+L6')
    ax1.loglog(P_up_range, J_baseline, 'gray', linewidth=1.5, linestyle=':', alpha=0.7, label='Baseline (perfect)')

    if np.any(valid_idx):
        P_ref1  = P_up_range[0]
        J_ref1  = J_full[valid_idx][0]
        ax1.loglog(P_up_range, J_ref1 * (P_up_range / P_ref1) ** 1.0,
                   'r--', linewidth=1.5, alpha=0.5, label='Slope = 1')

        P_ref05 = P_up_range[-1]
        J_ref05 = J_full[valid_idx][-1]
        ax1.loglog(P_up_range, J_ref05 * (P_up_range / P_ref05) ** 0.5,
                   'g--', linewidth=1.5, alpha=0.5, label='Slope = 0.5')

    regions = [
        {'mask': rate_limiting_arr == 'surface', 'color': 'red',    'label': 'Surface-limited'},
        {'mask': rate_limiting_arr == 'oxide',   'color': 'orange', 'label': 'Oxide-limited'},
        {'mask': rate_limiting_arr == 'metal',   'color': 'blue',   'label': 'Metal-limited'},
        {'mask': rate_limiting_arr == 'mixed',   'color': 'green',  'label': 'Mixed'},
    ]

    for region in regions:
        mask = region['mask'] & valid_idx
        if np.any(mask):
            idxs = np.where(mask)[0]
            for k, g in groupby(enumerate(idxs), lambda x: x[0] - x[1]):
                group = list(map(itemgetter(1), g))
                if len(group) > 2:
                    P_seg = P_up_range[group]
                    J_seg = J_full[group]
                    ax1.loglog(P_seg, J_seg, color=region['color'], linewidth=4, alpha=0.7)
                    slope, _ = np.polyfit(np.log10(P_seg), np.log10(np.abs(J_seg)), 1)
                    mid = len(group) // 2
                    ax1.text(P_seg[mid], J_seg[mid],
                             f"{region['label']}\nSlope={slope:.2f}",
                             color=region['color'], fontsize=9, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=11)
    ax1.set_ylabel('Flux $J$ (mol/m²/s)', fontsize=11)
    ax1.set_title('Flux vs Pressure (Full L3+L4+L6)', fontsize=12, fontweight='bold')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left', fontsize=8)

    # =========================================================================
    # Plot 2: Rate-Limiting Fractions (Top Right)
    # =========================================================================
    ax2 = axes[0, 1]
    ax2.semilogx(P_up_range, f_surface * 100, 'r-',      linewidth=2, label='Surface')
    ax2.semilogx(P_up_range, f_oxide   * 100, 'orange',  linewidth=2, label='Oxide')
    ax2.semilogx(P_up_range, f_metal   * 100, 'b-',      linewidth=2, label='Metal')
    ax2.axhline(50, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
    ax2.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=11)
    ax2.set_ylabel('Resistance Fraction (%)', fontsize=11)
    ax2.set_title('Rate-Limiting Step Analysis', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=9)

    # =========================================================================
    # Plot 3: Flux Comparison (Bottom Left)
    # =========================================================================
    ax3 = axes[1, 0]
    ax3.loglog(P_up_range[valid_idx], J_full[valid_idx],       'b-',  linewidth=2.5, label='Full L3+L4+L6')
    ax3.loglog(P_up_range[valid_idx], J_no_micro[valid_idx],   'g--', linewidth=2,   label='L3+L6 (no microstructure)')
    ax3.loglog(P_up_range[valid_idx], J_no_defects[valid_idx], 'r--', linewidth=2,   label='L4+L6 (no oxide defects)')
    ax3.loglog(P_up_range[valid_idx], J_baseline[valid_idx],   'k:',  linewidth=1.5, label='Baseline (perfect)')
    ax3.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=11)
    ax3.set_ylabel('Flux $J$ (mol/m²/s)', fontsize=11)
    ax3.set_title('Effect of Defects vs Microstructure', fontsize=12, fontweight='bold')
    ax3.grid(True, which='both', alpha=0.3)
    ax3.legend(loc='upper left', fontsize=9)

    info_text = (f'Oxide Defects:\n'
                 f'  Pinhole: {f_pinhole:.2f}%\n'
                 f'  Crack: {f_crack:.2f}% (γ={gamma:.2f})\n'
                 f'  GB: {f_gb_ox:.2f}% (δ={delta:.0f})\n'
                 f'Metal Microstructure:\n'
                 f'  Grain: {grain_size*1e6:.1f} μm\n'
                 f'  Traps: {trap_density:.1e} m⁻²')
    ax3.text(0.98, 0.02, info_text, transform=ax3.transAxes, fontsize=8,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # =========================================================================
    # Plot 4: Flux Contribution by Path (Bottom Right)
    # =========================================================================
    ax4 = axes[1, 1]
    ax4.semilogx(P_up_range, frac_intact * 100, 'k-', linewidth=2, label='Intact')
    if f_pinhole > 0:
        ax4.semilogx(P_up_range, frac_pinhole * 100, 'r-',      linewidth=2, label='Pinhole')
    if f_crack > 0:
        ax4.semilogx(P_up_range, frac_crack   * 100, 'orange',  linewidth=2, label='Crack')
    if f_gb_ox > 0:
        ax4.semilogx(P_up_range, frac_gb      * 100, 'g-',      linewidth=2, label='Oxide GB')
    ax4.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=11)
    ax4.set_ylabel('Flux Contribution (%)', fontsize=11)
    ax4.set_title('Which Oxide Path Dominates?', fontsize=12, fontweight='bold')
    ax4.set_ylim(0, 100)
    ax4.grid(True, alpha=0.3)
    ax4.legend(loc='best', fontsize=9)

    plt.tight_layout()
    plt.show()

    # =========================================================================
    # DataFrame — already built inside the loop, nothing to recompute
    # =========================================================================
    df = pd.DataFrame(rows)
    display(df)
```

```python

```

```python

```

```python

```

```python

```

```python

```
