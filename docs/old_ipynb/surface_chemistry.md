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
    display_name: rmg_rms_env
    language: python
    name: python3
---

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
import numpy as np
from scipy.optimize import brentq

# Import data dictionaries
from data.surface_kinetics_data import SURFACE_KINETICS, get_surface_kinetics
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS
R = 8.314  # J/(mol·K)
```

# L2+L6: Perfect Oxide + Perfect Metal with Surface Kinetics at Gas-Oxide Interface

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

$$ \sqrt{P} = \frac{1}{\sqrt{K_{eq}}} \cdot \frac{\theta}{1-\theta}$$

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

$$fraction_i = \frac{R_i}{R_{total}}, \quad R_{total} = R_{surface} + R_{oxide} + R_{metal}$$

**Rate-limiting identification:**
- $fraction_{surface} > 0.5$ → Surface-limited (slope ≈ 1 in log J vs log P)
- $fraction_{oxide} > 0.5$ → Oxide-limited (slope ≈ 0.5)
- $fraction_{metal} > 0.5$ → Metal-limited (slope ≈ 0.5)
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
    surfraction_kin = get_surface_kinetics(oxide_name, temperature_K)
    
    # 2. Oxide transport properties
    oxide_data = OXIDE_PROPERTIES[oxide_name]
    T_refraction_ox = oxide_data['T_ref']
    D_ox = oxide_data['D_ox_ref'] * np.exp((-oxide_data['E_D_ox'] / R) * (1/temperature_K - 1/T_refraction_ox))
    K_ox = oxide_data['K_ox_ref'] * np.exp((-oxide_data['H_sol_ox'] / R) * (1/temperature_K - 1/T_refraction_ox))
    L_ox = oxide_data['thickness']
    
    # 3. Metal transport properties
    metal_data = MATERIALS[metal_name]
    T_refraction_m = metal_data['T_ref']
    D_m = metal_data['D_ref'] * np.exp((-metal_data['E_D'] / R) * (1/temperature_K - 1/T_refraction_m))
    K_s_m = metal_data['K_s_ref'] * np.exp((-metal_data['H_s'] / R) * (1/temperature_K - 1/T_refraction_m))
    
    
    return {
        # Surface kinetics
        'k_diss': surfraction_kin['k_diss'],
        'k_recomb': surfraction_kin['k_recomb'],
        'K_eq': surfraction_kin['K_eq'],
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

#### PATCH

```python
# ...existing code...

# Add once near other helper functions (e.g., after g_theta / sqrt_P_int_from_theta)
def smooth_surface_resistance(k_diss_eff, P_up, theta, eps_theta=1e-12, eps_rate=1e-60):
    """
    Smooth surface resistance model (no hard cutoff):
        R_surface ~ 1 / (k_diss * P_up * (1 - theta)^2)
    """
    theta_clamped = min(max(theta, 0.0), 1.0 - eps_theta)
    vacant = max(1.0 - theta_clamped, eps_theta)
    rate = k_diss_eff * max(P_up, eps_theta) * (vacant ** 2)
    return 1.0 / max(rate, eps_rate)

# # ...existing code...

# # In solve_steady_state_flux(...), replace:
# # R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0

# # ...existing code...
# R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
# # ...existing code...


# # In calculate_defective_metal_flux_L6(...), replace:
# # R_surface_approx = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0

# # ...existing code...
# R_surface_approx = smooth_surface_resistance(k_diss, P_up, theta_ss)
# # ...existing code...

# # In calculate_defective_metal_flux_L6(...), replace:
# # R_surface_approx = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0

# # ...existing code...
# R_surface_approx = smooth_surface_resistance(k_diss, P_up, theta_ss)
# # ...existing code...

# # In calculate_path_flux_L6(...), finite-kinetics pinhole branch replace:
# # R_surface = 1.0 / (k_diss_eff * P_up) if theta_ss < 0.9 else 0.0

# # ...existing code...
# R_surface = smooth_surface_resistance(k_diss_eff, P_up, theta_ss)
# R_oxide = 0.0
# R_metal = 1.0 / beta
# # ...existing code...

# # In calculate_path_flux_L6(...), non-pinhole branch replace:
# # R_surface = 1.0 / (k_diss * P_up) if theta_ss < 1 else 0.0

# # ...existing code...
# R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
# R_oxide = 1.0 / alpha
# R_metal = 1.0 / beta
# # ...existing code...


# # In calculate_path_flux_L346_v2(...), finite-kinetics pinhole branch replace:
# # R_surface = 1.0 / (k_diss_eff * P_up) if theta_ss < 1 else 0.0

# # ...existing code...
# R_surface = smooth_surface_resistance(k_diss_eff, P_up, theta_ss)
# # ...existing code...

# # In calculate_path_flux_L346_v2(...), non-pinhole branch replace:
# # R_surface = 1.0 / (k_diss * P_up) if theta_ss < 1 else 0.0

# # ...existing code...
# R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss)
# # ...existing code...


# # Optional diagnostics: in interactive_L36_solver(...), add to DataFrame rows in BOTH branches:

# # ...existing code...
# rows.append({
#     "P_up (Pa)":           P_up,
#     # ...existing keys...
#     "fraction_surface (%)": ws * 100,
#     "fraction_oxide (%)":   wo * 100,
#     "fraction_metal (%)":   wm * 100,
#     "pinhole_flux_fraction (%)": frac_pinhole * 100,
#     # ...existing keys...
# })
# # ...existing code...
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
    # R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss) #PATCH
    
    # Oxide resistance
    R_oxide = 1.0 / alpha   # = L_ox / (D_ox × K_ox)
    
    # Metal resistance
    
    R_metal = 1.0 / beta
    
    # Total resistance
    R_total = R_surface + R_oxide + R_metal
    
    # Fractional contributions
    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide = R_oxide / R_total if R_total > 0 else 0
    fraction_metal = R_metal / R_total if R_total > 0 else 0
    

    if fraction_surface > 0.5:
        rate_limiting = 'surface (dissociation)'
    elif fraction_oxide > 0.5:
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
            'fraction_surface': fraction_surface,
            'fraction_oxide': fraction_oxide,
            'fraction_metal': fraction_metal,
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
    # R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss) #PATCH
    
    # Oxide resistance
    R_oxide = 1.0 / alpha   # = L_ox / (D_ox × K_ox)
    
    # Metal resistance
    R_metal = 1.0 / beta
    
    # Total resistance
    R_total = R_surface + R_oxide + R_metal
    
    # Fractional contributions
    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide = R_oxide / R_total if R_total > 0 else 0
    fraction_metal = R_metal / R_total if R_total > 0 else 0
    
    # FIX: Add 'mixed' case when no single resistance dominates
    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
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
            'fraction_surface': fraction_surface,
            'fraction_oxide': fraction_oxide,
            'fraction_metal': fraction_metal,
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
    # P_down=widgets.FloatLogSlider(value=0, base=10, min=0, max=4, step=0.5, description='P_down (Pa)'),
    P_down=widgets.FloatSlider(value=0, min=0, max=10000, step=100, description='P_down (Pa)'),
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

    P_up_range = np.logspace(-3, 12, 40)  # 1 Pa to 10 GPa

    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    J_system = []
    theta_values = []
    fraction_surface_list = []
    fraction_oxide_list = []
    fraction_metal_list = []
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
            fraction_s = r['resistances']['fraction_surface']
            fraction_o = r['resistances']['fraction_oxide']
            fraction_m = r['resistances']['fraction_metal']
            
            # Derive rate-limiting label ONCE (same logic as function)
            if   fraction_s > 0.5: rate_lim = 'surface'
            elif fraction_o > 0.5: rate_lim = 'oxide'
            elif fraction_m > 0.5: rate_lim = 'metal'
            else:           rate_lim = 'mixed'
            
            # Append to plot arrays
            J_system.append(J_ss)
            theta_values.append(theta)
            fraction_surface_list.append(fraction_s)
            fraction_oxide_list.append(fraction_o)
            fraction_metal_list.append(fraction_m)
            
            # Append to DataFrame rows (same values, no recomputation)
            rows.append({
                "P_up (Pa)":           P_up,
                "P_int (Pa)":          r["P_int"],
                "J_ss (mol/m²/s)":     J_ss,
                "θ_surface":           theta,
                "fraction_surface (%)":       fraction_s * 100,
                "fraction_oxide (%)":         fraction_o * 100,
                "fraction_metal (%)":         fraction_m * 100,
                "Rate-Limiting":       rate_lim.upper(),
                "α_oxide":             r["alpha"],
                "β_metal":             r["beta"],
            })
            
        except Exception as e:
            J_system.append(np.nan)
            theta_values.append(np.nan)
            fraction_surface_list.append(np.nan)
            fraction_oxide_list.append(np.nan)
            fraction_metal_list.append(np.nan)
            rows.append({"P_up (Pa)": P_up, "Rate-Limiting": "ERROR", "Error": str(e)})

    # Convert to arrays
    J_system = np.array(J_system)
    theta_values = np.array(theta_values)
    fraction_surface = np.array(fraction_surface_list)
    fraction_oxide = np.array(fraction_oxide_list)
    fraction_metal = np.array(fraction_metal_list)

    # =========================================================================
    # Rate-limiting array for plot shading — 4-branch classification
    # =========================================================================
    rate_limiting_arr = np.where(
        fraction_surface > 0.5, 'surface',
        np.where(fraction_oxide > 0.5, 'oxide',
        np.where(fraction_metal > 0.5, 'metal',
                 'mixed'))
    )

    # =========================================================================
    # Create figure with 2 subplots
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # ===== Plot 1: Flux vs Pressure =====
    valid_idx = ~np.isnan(J_system)
    ax1.loglog(P_up_range, J_system, 'k-', linewidth=2.5, label='L2+L6 Model')
    # Slope calculation
    logP_net = np.log10(P_up_range)
    logJ_net = np.log10(np.abs(J_system))
    slope_net, _ = np.polyfit(logP_net, logJ_net, 1)

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
                    
    ax1.text(
    0.98, 0.02,
    f'Net_slope = {slope_net:.2f}',
    transform=ax1.transAxes,
    ha='right', va='bottom',
    fontsize=10,
    fontweight='bold',
    color='black',
    bbox=dict(boxstyle='square', fc='wheat', ec='gray', alpha=1)
)
    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L2+L6: Flux vs Pressure', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    # ===== Plot 2: Rate-Limiting Fractions =====
    ax2.semilogx(P_up_range, fraction_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, fraction_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, fraction_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
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
    P_down=widgets.FloatSlider(value=0, min=0, max=10000, step=100, description='P_down (Pa)'),
    L_m=widgets.FloatLogSlider(value=1e-3, base=10, min=-4, max=-1, step=0.5, description='L_m (m)'),
    k_diss=widgets.FloatLogSlider(value=1e-15, base=10, min=-18, max=-3, step=0.5, description='k_diss'),
    K_eq=widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-1, step=0.5, description='K_eq'),
    D_ox=widgets.FloatLogSlider(value=1e-11, base=10, min=-18, max=-5, step=0.5, description='D_ox (m²/s)'),
    K_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-14, max=-4, step=0.5, description='K_ox'),
    L_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-8, max=-4, step=0.5, description='L_ox (m)'),
    D_m=widgets.FloatLogSlider(value=1.0e-12, base=10, min=-13, max=-6, step=0.5, description='D_m (m²/s)'),
    K_s_m=widgets.FloatLogSlider(value=3.16e-4, base=10, min=-6, max=0, step=0.5, description='K_s_m'),
)
def interactive_L26_solver(oxide_name, metal_name,       
                           P_down, L_m, k_diss, K_eq,
                           D_ox, K_ox, L_ox, D_m, K_s_m):

    P_up_range = np.logspace(-3, 12, 40)

    J_system = []
    theta_values = []
    fraction_surface_list = []
    fraction_oxide_list = []
    fraction_metal_list = []

    for P_up in P_up_range:
        try:
            result = solve_steady_state_flux_direct(
                P_up, P_down, L_m, k_diss, K_eq,
                D_ox, K_ox, L_ox, D_m, K_s_m
            )
            J_system.append(result['J_ss'])
            theta_values.append(result['theta'])
            fraction_surface_list.append(result['resistances']['fraction_surface'])
            fraction_oxide_list.append(result['resistances']['fraction_oxide'])
            fraction_metal_list.append(result['resistances']['fraction_metal'])
        except:
            J_system.append(np.nan)
            theta_values.append(np.nan)
            fraction_surface_list.append(np.nan)
            fraction_oxide_list.append(np.nan)
            fraction_metal_list.append(np.nan)

    J_system = np.array(J_system)
    theta_values = np.array(theta_values)
    fraction_surface = np.array(fraction_surface_list)
    fraction_oxide = np.array(fraction_oxide_list)
    fraction_metal = np.array(fraction_metal_list)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    valid_idx = ~np.isnan(J_system)
    ax1.loglog(P_up_range, J_system, 'k-', linewidth=2.5, label='L2+L6 Model')
    # Slope calculation
    logP_net = np.log10(P_up_range)
    logJ_net = np.log10(np.abs(J_system))
    slope_net, _ = np.polyfit(logP_net, logJ_net, 1)

    if np.any(valid_idx):
        P_ref1 = P_up_range[0]
        J_ref1 = J_system[0] if not np.isnan(J_system[0]) else J_system[valid_idx][0]
        J_slope1 = J_ref1 * (P_up_range / P_ref1) ** 1.0
        ax1.loglog(P_up_range, J_slope1, 'r--', linewidth=1.5, alpha=0.5, label='Slope = 1 (surface)')

        P_ref05 = P_up_range[-1]
        J_ref05 = J_system[-1] if not np.isnan(J_system[-1]) else J_system[valid_idx][-1]
        J_slope05 = J_ref05 * (P_up_range / P_ref05) ** 0.5
        ax1.loglog(P_up_range, J_slope05, 'g--', linewidth=1.5, alpha=0.5, label='Slope = 0.5 (diffusion)')

    rate_limiting_arr = np.where(fraction_surface > 0.5, 'surface', np.where(fraction_oxide > 0.5, 'oxide', 'metal'))
    regions = [
        {'mask': rate_limiting_arr == 'surface', 'color': 'red',    'label': 'Surface-limited'},
        {'mask': rate_limiting_arr == 'metal',   'color': 'blue',   'label': 'Metal-limited'},
        {'mask': rate_limiting_arr == 'oxide',   'color': 'orange', 'label': 'Oxide-limited'},
        {'mask': rate_limiting_arr == 'mixed',   'color': 'green',  'label': 'Mixed'},
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

    ax1.text(
    0.98, 0.02,
    f'Net_slope = {slope_net:.2f}',
    transform=ax1.transAxes,
    ha='right', va='bottom',
    fontsize=10,
    fontweight='bold',
    color='black',
    bbox=dict(boxstyle='square', fc='wheat', ec='gray', alpha=1)
)
    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L2+L6: Flux vs Pressure', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    ax2.semilogx(P_up_range, fraction_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, fraction_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, fraction_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
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
                "fraction_surface (%)":   r["resistances"]["fraction_surface"] * 100,
                "fraction_oxide (%)":     r["resistances"]["fraction_oxide"]   * 100,
                "fraction_metal (%)":     r["resistances"]["fraction_metal"]   * 100,
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
    #R_surface_approx = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    R_surface_approx = smooth_surface_resistance(k_diss, P_up, theta_ss)

    # Oxide resistance
    R_oxide = 1.0 / alpha  # = L_ox / (D_ox × K_ox)

    # Metal resistance (with microstructure)
    R_metal = 1.0 / beta_eff  # = thickness / (D_eff × K_s_m)

    # Total resistance
    R_total = R_surface_approx + R_oxide + R_metal

    # Fractional contributions
    fraction_surface = R_surface_approx / R_total if R_total > 0 else 0
    fraction_oxide = R_oxide / R_total if R_total > 0 else 0
    fraction_metal = R_metal / R_total if R_total > 0 else 0

    # Determine rate-limiting step (>50% of resistance)
    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
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
            'fraction_surface': fraction_surface,
            'fraction_oxide': fraction_oxide,
            'fraction_metal': fraction_metal,
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
    P_down=widgets.FloatSlider(value=0, min=0, max=10000, step=100, description='P_down (Pa)'),
    thickness=widgets.FloatLogSlider(value=1e-3, base=10, min=-4, max=-1, step=0.5, description='L_m (m)'),
    temperature=widgets.IntSlider(value=973, min=573, max=1273, step=50, description='T (K)'),
    # Surface kinetics
    k_diss=widgets.FloatLogSlider(value=1e-15, base=10, min=-18, max=-3, step=0.5, description='k_diss'),
    K_eq=widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-1, step=0.5, description='K_eq'),
    # Oxide properties
    D_ox=widgets.FloatLogSlider(value=1e-11, base=10, min=-18, max=-5, step=0.5, description='D_ox (m²/s)'),
    K_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-14, max=-4, step=0.5, description='K_ox'),
    L_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-8, max=-4, step=0.5, description='L_ox (m)'),
    # Metal properties
    D_lattice=widgets.FloatLogSlider(value=1.0e-12, base=10, min=-13, max=-6, step=0.5, description='D_lattice'),
    K_s_m=widgets.FloatLogSlider(value=3.16e-4, base=10, min=-6, max=0, step=0.5, description='K_s_m'),
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
    
    P_up_range = np.logspace(-3, 12, 100)  # 1 Pa to 1 TPa

    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    J_system = []
    theta_values = []
    fraction_surface_list = []
    fraction_oxide_list = []
    fraction_metal_list = []
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
            fraction_s = r['resistances']['fraction_surface']
            fraction_o = r['resistances']['fraction_oxide']
            fraction_m = r['resistances']['fraction_metal']
            
            # Derive rate-limiting label ONCE (same logic as function)
            if   fraction_s > 0.5: rate_lim = 'surface'
            elif fraction_o > 0.5: rate_lim = 'oxide'
            elif fraction_m > 0.5: rate_lim = 'metal'
            else:           rate_lim = 'mixed'
            
            # Append to plot arrays
            J_system.append(J_ss)
            theta_values.append(theta)
            fraction_surface_list.append(fraction_s)
            fraction_oxide_list.append(fraction_o)
            fraction_metal_list.append(fraction_m)
            
            # Append to DataFrame rows (same values, no recomputation)
            rows.append({
                "P_up (Pa)":           P_up,
                "P_int (Pa)":          r["P_int"],
                "J_ss (mol/m²/s)":     J_ss,
                "θ_surface":           theta,
                "D_eff (m²/s)":        r["D_eff"],
                "D_eff/D_lattice":     r["modification_factor"],
                "fraction_surface (%)":       fraction_s * 100,
                "fraction_oxide (%)":         fraction_o * 100,
                "fraction_metal (%)":         fraction_m * 100,
                "Rate-Limiting":       rate_lim.upper(),
                "α_oxide":             r["alpha"],
                "β_eff":               r["beta_eff"],
            })
            
        except Exception as e:
            J_system.append(np.nan)
            theta_values.append(np.nan)
            fraction_surface_list.append(np.nan)
            fraction_oxide_list.append(np.nan)
            fraction_metal_list.append(np.nan)
            rows.append({"P_up (Pa)": P_up, "Rate-Limiting": "ERROR", "Error": str(e)})

    # Convert to arrays
    J_system = np.array(J_system)
    theta_values = np.array(theta_values)
    fraction_surface = np.array(fraction_surface_list)
    fraction_oxide = np.array(fraction_oxide_list)
    fraction_metal = np.array(fraction_metal_list)

    # =========================================================================
    # Rate-limiting array for plot shading — 4-branch classification
    # =========================================================================
    rate_limiting_arr = np.where(
        fraction_surface > 0.5, 'surface',
        np.where(fraction_oxide > 0.5, 'oxide',
        np.where(fraction_metal > 0.5, 'metal',
                 'mixed'))
    )




# Diagnostic: verify Plot 1 shading matches Plot 2 crossings
    print("\n=== Rate-Limiting Consistency Check ===")
    for i in [0, len(P_up_range)//2, -1]:  # Check start, middle, end
        P = P_up_range[i]
        fs, fo, fm = fraction_surface[i]*100, fraction_oxide[i]*100, fraction_metal[i]*100
        rl = rate_limiting_arr[i]
        print(f"P={P:.1e}: fraction_surf={fs:.1f}%, fraction_ox={fo:.1f}%, fraction_metal={fm:.1f}% → {rl}")
        
        # Verify classification matches
        if fs > 50 and rl != 'surface':
            print(f"  ⚠️ MISMATCH: fraction_surface > 50% but rate_limiting={rl}")
        if fo > 50 and rl != 'oxide':
            print(f"  ⚠️ MISMATCH: fraction_oxide > 50% but rate_limiting={rl}")
        if fm > 50 and rl != 'metal':
            print(f"  ⚠️ MISMATCH: fraction_metal > 50% but rate_limiting={rl}")

    # =========================================================================
    # Create figure with 2 subplots
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # ===== Plot 1: Flux vs Pressure =====
    valid_idx = ~np.isnan(J_system)
    ax1.loglog(P_up_range, J_system, 'k-', linewidth=2.5, label='L2+L4+L6 Model')
    # Slope calculation
    logP_net = np.log10(P_up_range)
    logJ_net = np.log10(np.abs(J_system))
    slope_net, _ = np.polyfit(logP_net, logJ_net, 1)

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

    ax1.text(0.98, 0.02, f'Net_slope = {slope_net:.2f}', transform=ax1.transAxes,ha='right', va='bottom',
    fontsize=10,fontweight='bold',color='black',bbox=dict(boxstyle='square', fc='wheat', ec='gray', alpha=1))
    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L2+L4+L6: Flux vs Pressure', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    # ===== Plot 2: Rate-Limiting Fractions =====
    ax2.semilogx(P_up_range, fraction_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, fraction_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, fraction_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
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

$$\boxed{J_{total} = fraction_{intact} \cdot J_{intact} + fraction_{defect} \cdot J_{defect}}$$

Where:
- $fraction_{intact} = 1 - fraction_{defect}$ = area fraction of intact oxide
- $fraction_{defect}$ = area fraction with defects (pinholes, cracks, grain boundaries)
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

$$J_{defect} = \frac{fraction_{ph} \cdot J_{pinhole} + fraction_{crack} \cdot J_{crack} + fraction_{gb} \cdot J_{gb}}{fraction_{ph} + fraction_{crack} + fraction_{gb}}$$

Where $fraction_{ph} + fraction_{crack} + fraction_{gb} = fraction_{defect}$

Or more directly:

$$J_{total} = fraction_{intact} \cdot J_{intact} + fraction_{ph} \cdot J_{pinhole} + fraction_{crack} \cdot J_{crack} + fraction_{gb} \cdot J_{gb}$$


### Step 5: Total System Flux

$$\boxed{J_{total} = (1 - fraction_{defect}) \cdot J_{intact} + fraction_{defect} \cdot J_{defect}}$$


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

$$fraction_i = \frac{R_i}{R_{total}} \quad \text{for } i \in \{\text{surface, oxide, metal}\}$$

**Rate-Limiting Identification:**

| Condition | Rate-Limiting Step | Pressure Dependence |
|-----------|-------------------|---------------------|
| $fraction_{surface} > 0.5$ | Surface dissociation | $J \propto P$ (slope ≈ 1) |
| $fraction_{oxide} > 0.5$ | Oxide diffusion | $J \propto \sqrt{P}$ (slope ≈ 0.5) |
| $fraction_{metal} > 0.5$ | Metal diffusion | $J \propto \sqrt{P}$ (slope ≈ 0.5) |


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

    J_total = fraction_intact × J_intact + fraction_defect × J_defect

Where:
    fraction_intact = (1 - fraction_defect) = area fraction of intact oxide
    fraction_defect = area fraction of defects (pinholes, cracks, grain boundaries)
    
Each path uses the L2+L6 coupled model with modified oxide permeance:
    - Intact: α_intact = D_ox × K_ox / L_ox  (full oxide barrier)
    - Pinhole: α → ∞ (no oxide, direct metal exposure)
    - Crack: α_crack = D_ox × K_ox / (γ × L_ox) where γ < 1 (thinner oxide)
    - GB: α_gb = δ × D_ox × K_ox / L_ox where δ > 1 (enhanced diffusivity)

The key insight: Each parallel path still satisfies:
    J_surface = J_oxide = J_metal (steady-state for that path)
"""
def calculate_path_flux_L6(
    P_up, P_down, L_m,
    k_diss, K_eq,
    alpha,  # Effective oxide permeance for this path (use np.inf for pinhole)
    D_m, K_s_m,
    path_type='intact',
    # NEW: Metal surface kinetics for pinhole path
    k_diss_metal=None,  # If None, use Sieverts' law limit for pinhole
    K_eq_metal=None,
):
    """
    Calculate steady-state flux through a single path (intact or defect).
    
    For pinhole paths:
    - If k_diss_metal and K_eq_metal are provided: use metal surface kinetics
    - If not provided: assume fast kinetics (Sieverts' law limit)
    
    Parameters
    ----------
    ... (existing parameters) ...
    k_diss_metal : float, optional
        Metal surface dissociation rate [mol/m²/s/Pa]. For pinhole only.
    K_eq_metal : float, optional
        Metal surface equilibrium constant [Pa⁻¹]. For pinhole only.
    """
    
    # Metal permeance
    beta = D_m * K_s_m / L_m
    sqrt_P_down = np.sqrt(max(P_down, 0))
    sqrt_P_up = np.sqrt(max(P_up, 0))
    
    # Check if pinhole (no oxide)
    is_pinhole = (path_type == 'pinhole' or alpha == np.inf or alpha > 1e10)
    
    # =========================================================================
    # PINHOLE PATH: Use metal kinetics or Sieverts' limit
    # =========================================================================
    if is_pinhole:
        # Determine which kinetics to use for pinhole
        if k_diss_metal is not None and K_eq_metal is not None:
            # Use metal surface kinetics
            k_diss_eff = k_diss_metal
            K_eq_eff = K_eq_metal
            use_sieverts_limit = False
        else:
            # Assume fast metal kinetics → Sieverts' law limit
            # In this limit, √P_int = √P_up (no surface resistance)
            use_sieverts_limit = True
        
        if use_sieverts_limit:
            # =========================================================
            # Sieverts' Law Limit (fast metal surface kinetics)
            # =========================================================
            # No surface resistance: √P_int = √P_up
            sqrt_P_int = sqrt_P_up
            P_int = P_up
            
            # Flux is purely metal-diffusion limited
            J_ss = beta * (sqrt_P_int - sqrt_P_down)
            
            # Surface coverage is at equilibrium (θ → 0 for fast kinetics)
            theta_ss = 0.0  # Or could estimate from K_eq_metal if provided
            
            # No surface flux calculation needed
            J_surf = J_ss  # By definition at steady state
            J_ox = np.nan  # No oxide
            
            # Resistances
            R_surface = 0.0  # Fast kinetics = no resistance
            R_oxide = 0.0    # No oxide
            R_metal = 1.0 / beta
            R_total = R_metal
            
            fraction_surface = 0.0
            fraction_oxide = 0.0
            fraction_metal = 1.0
            rate_limiting = 'metal'
            
        else:
            # =========================================================
            # Metal Surface Kinetics (finite k_diss_metal)
            # =========================================================
            def residual(theta):
                if theta <= 0 or theta >= 1:
                    return np.inf
                
                # Surface flux with METAL kinetics
                J_surf = surface_flux(theta, P_up, k_diss_eff, K_eq_eff)
                
                # At pinhole: √P_int = g(θ) with METAL K_eq
                sqrt_P_int = g_theta(theta, K_eq_eff)
                
                # Metal flux
                J_metal = beta * (sqrt_P_int - sqrt_P_down)
                
                return J_surf - J_metal
            
            try:
                theta_ss = brentq(residual, 1e-10, 1.0 - 1e-10)
            except ValueError as e:
                return {
                    'flux': np.nan, 'theta': np.nan, 'P_int': np.nan,
                    'path_type': path_type, 'alpha': np.inf, 'beta': beta,
                    'error': str(e)
                }
            
            sqrt_P_int = g_theta(theta_ss, K_eq_eff)
            P_int = sqrt_P_int**2
            J_ss = beta * (sqrt_P_int - sqrt_P_down)
            J_surf = surface_flux(theta_ss, P_up, k_diss_eff, K_eq_eff)
            J_ox = np.nan  # No oxide
            
            # Resistances with METAL kinetics
            # R_surface = 1.0 / (k_diss_eff * P_up) if theta_ss < 0.9 else 0.0
            R_surface = smooth_surface_resistance(k_diss_eff, P_up, theta_ss) #PATCH
            R_oxide = 0.0  # No oxide
            R_metal = 1.0 / beta
            R_total = R_surface + R_metal
            
            fraction_surface = R_surface / R_total if R_total > 0 else 0
            fraction_oxide = 0.0
            fraction_metal = R_metal / R_total if R_total > 0 else 0
            
            if fraction_surface > 0.5:
                rate_limiting = 'surface'
            elif fraction_metal > 0.5:
                rate_limiting = 'metal'
            else:
                rate_limiting = 'mixed'
        
        # Return pinhole result
        return {
            'flux': J_ss,
            'theta': theta_ss,
            'P_int': P_int,
            'path_type': path_type,
            'alpha': np.inf,
            'beta': beta,
            'kinetics_used': 'sieverts' if use_sieverts_limit else 'metal',
            'flux_balance': {
                'J_surface': J_surf,
                'J_oxide': J_ox,
                'J_metal': J_ss,
                'balanced': True
            },
            'resistances': {
                'R_surface': R_surface,
                'R_oxide': R_oxide,
                'R_metal': R_metal,
                'R_total': R_total,
                'fraction_surface': fraction_surface,
                'fraction_oxide': fraction_oxide,
                'fraction_metal': fraction_metal,
            },
            'rate_limiting': rate_limiting
        }
    
    # =========================================================================
    # NON-PINHOLE PATHS (intact, crack, GB): Use oxide kinetics (unchanged)
    # =========================================================================
    def residual(theta):
        if theta <= 0 or theta >= 1:
            return np.inf
        
        # Surface flux with OXIDE kinetics
        J_surf = surface_flux(theta, P_up, k_diss, K_eq)
        
        # Standard case with finite α
        sqrt_P_int = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
        
        # Metal flux
        J_metal = beta * (sqrt_P_int - sqrt_P_down)
        
        return J_surf - J_metal
    
    try:
        theta_ss = brentq(residual, 1e-10, 1.0 - 1e-10)
    except ValueError as e:
        return {
            'flux': np.nan, 'theta': np.nan, 'P_int': np.nan,
            'path_type': path_type, 'alpha': alpha, 'beta': beta,
            'error': str(e)
        }
    
    sqrt_P_int = sqrt_P_int_from_theta(theta_ss, alpha, beta, K_eq, P_down)
    P_int = sqrt_P_int**2
    J_ss = beta * (sqrt_P_int - sqrt_P_down)
    J_surf = surface_flux(theta_ss, P_up, k_diss, K_eq)
    J_ox = oxide_flux(theta_ss, alpha, beta, K_eq, P_down)
    
    #R_surface = 1.0 / (k_diss * P_up) if theta_ss < 0.9 else 0.0
    R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss) #PATCH
    R_oxide = 1.0 / alpha
    R_metal = 1.0 / beta
    R_total = R_surface + R_oxide + R_metal
    
    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide = R_oxide / R_total if R_total > 0 else 0
    fraction_metal = R_metal / R_total if R_total > 0 else 0
    
    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
        rate_limiting = 'metal'
    else:
        rate_limiting = 'mixed'
    
    return {
        'flux': J_ss,
        'theta': theta_ss,
        'P_int': P_int,
        'path_type': path_type,
        'alpha': alpha,
        'beta': beta,
        'kinetics_used': 'oxide',
        'flux_balance': {
            'J_surface': J_surf,
            'J_oxide': J_ox,
            'J_metal': J_ss,
            'balanced': np.allclose(J_surf, J_ox, rtol=1e-6) and np.allclose(J_ox, J_ss, rtol=1e-6)
        },
        'resistances': {
            'R_surface': R_surface,
            'R_oxide': R_oxide,
            'R_metal': R_metal,
            'R_total': R_total,
            'fraction_surface': fraction_surface,
            'fraction_oxide': fraction_oxide,
            'fraction_metal': fraction_metal,
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
    k_diss_metal=None,
    K_eq_metal=None,
):
    """
    Calculate total flux through parallel intact + defect paths.
    
    Implements Strehlow & Savage (1974) parallel path model with L6 surface kinetics.
    
    J_total = (1 - fraction_defect) × J_intact + fraction_defect × J_defect
    
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
        fraction_defect: fraction of surface area with defects (0 to 1)
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
        path_type=defect_type,
        k_diss_metal=k_diss_metal,  # Pass metal kinetics
        K_eq_metal=K_eq_metal
    )
    
    # =========================================================================
    # Area-weighted total flux
    # =========================================================================
    fraction_intact = 1.0 - defect_area_fraction
    fraction_defect = defect_area_fraction
    
    J_intact = intact_result['flux']
    J_defect = defect_result['flux']
    
    J_total = fraction_intact * J_intact + fraction_defect * J_defect
    
    # Enhancement factor compared to intact-only case
    enhancement_factor = J_total / J_intact if J_intact > 0 else np.inf
    
    # Determine which path dominates the total flux
    flux_from_intact = fraction_intact * J_intact
    flux_from_defect = fraction_defect * J_defect
    
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
        'fraction_intact': fraction_intact,
        'fraction_defect': fraction_defect,
        
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
        'alpha_ratio': alpha_defect / alpha_intact if alpha_defect != np.inf else np.inf,
        
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

    J_total = fraction_intact × J_intact + Σᵢ (fraction_defect,i × J_defect,i)

Where the sum is over all defect types i.

Constraint: fraction_intact + Σᵢ fraction_defect,i = 1

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
    k_diss_metal=None,  # Metal kinetics for pinhole
    K_eq_metal=None,    # Metal kinetics for pinhole
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
        Surface kinetics parameters (oxide surface)
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
    k_diss_metal : float, optional
        Metal surface dissociation rate [mol/m²/s/Pa]. For pinhole only.
    K_eq_metal : float, optional
        Metal surface equilibrium constant [Pa⁻¹]. For pinhole only.
    
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
    
    fraction_intact = 1.0 - total_defect_fraction
    
    # =========================================================================
    # Calculate intact oxide permeance
    # =========================================================================
    alpha_intact = D_ox * K_ox / L_ox
    
    # =========================================================================
    # Calculate flux through intact path
    # =========================================================================
    try:
        intact_result = calculate_path_flux_L6(
            P_up, P_down, L_m,
            k_diss, K_eq,
            alpha_intact,
            D_m, K_s_m,
            path_type='intact'
        )
        
        if 'error' in intact_result:
            return {
                'J_total': np.nan,
                'error': f"Intact path failed: {intact_result['error']}",
            }
    except Exception as e:
        return {
            'J_total': np.nan,
            'error': f"Intact path calculation failed: {e}",
        }
    
    # =========================================================================
    # Calculate flux through each defect path
    # =========================================================================
    defect_results = {}
    
    for defect_type, config in defect_config.items():
        fraction_defect = config.get('area_fraction', 0.0)
        
        if fraction_defect <= 0:
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
        else:
            continue  # Skip unknown types (already validated above)
        
        # Calculate flux through this defect type with error handling
        try:
            path_result = calculate_path_flux_L6(
                P_up, P_down, L_m,
                k_diss, K_eq,
                alpha_defect,
                D_m, K_s_m,
                path_type=defect_type,
                k_diss_metal=k_diss_metal if defect_type == 'pinhole' else None,
                K_eq_metal=K_eq_metal if defect_type == 'pinhole' else None
            )
            
            # Check for error in result
            if 'error' in path_result:
                print(f"Warning: {defect_type} path failed: {path_result['error']}")
                continue
                
        except Exception as e:
            print(f"Warning: {defect_type} path calculation failed: {e}")
            continue
        
        defect_results[defect_type] = {
            'area_fraction': fraction_defect,
            'alpha': alpha_defect,
            'alpha_ratio': alpha_defect / alpha_intact if alpha_defect != np.inf else np.inf,
            'path_result': path_result,
            'flux_contribution': fraction_defect * path_result['flux'],
        }
    
    # =========================================================================
    # Calculate total flux (area-weighted sum)
    # =========================================================================
    J_intact_contribution = fraction_intact * intact_result['flux']
    J_total = J_intact_contribution
    
    for defect_type, data in defect_results.items():
        J_total += data['flux_contribution']
    
    # =========================================================================
    # Flux breakdown analysis
    # =========================================================================
    flux_breakdown = {
        'intact': {
            'area_fraction': fraction_intact,
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
        'fraction_intact': fraction_intact,
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
        
        # Units dictionary
        'units': {
            'flux': 'mol/m²/s',
            'pressure': 'Pa',
            'permeance': 'mol/m²/s/Pa^0.5',
            'area_fraction': 'dimensionless',
        }
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
    P_down=widgets.FloatSlider(value=0,  min=0, max=10000, step=100, description='P_down (Pa)'),
    L_m=widgets.FloatLogSlider(value=1e-3, base=10, min=-4, max=-1, step=0.5, description='L_m (m)'),
    # Surface kinetics
    k_diss=widgets.FloatLogSlider(value=1e-15, base=10, min=-18, max=-10, step=0.5, description='k_diss'),
    K_eq=widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-5, step=0.5, description='K_eq'),
    # ADD: Metal surface kinetics for pinhole
    k_diss_metal=widgets.FloatLogSlider(value=1e-12, base=10, min=-15, max=-8, step=0.5, description='k_diss_metal'),
    K_eq_metal=widgets.FloatLogSlider(value=1e-8, base=10, min=-12, max=-4, step=0.5, description='K_eq_metal'),
    use_sieverts_pinhole=widgets.Checkbox(value=False, description='Sieverts limit for pinhole'),
    # Oxide properties
    D_ox=widgets.FloatLogSlider(value=1e-11, base=10, min=-22, max=-9, step=0.5, description='D_ox (m²/s)'),
    K_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-12, max=0, step=0.5, description='K_ox'),
    L_ox=widgets.FloatLogSlider(value=1e-6, base=10, min=-8, max=-4, step=0.5, description='L_ox (m)'),
    # Metal properties
    D_m=widgets.FloatLogSlider(value=1e-12, base=10, min=-12, max=-6, step=0.5, description='D_m (m²/s)'),
    K_s_m=widgets.FloatLogSlider(value=3.16e-4, base=10, min=-6, max=0, step=0.5, description='K_s_m'),
    # Defect parameters
    fraction_pinhole=widgets.FloatSlider(value=0.01, min=0, max=5, step=0.1, description='Pinhole %'),
    fraction_crack=widgets.FloatSlider(value=0.05, min=0, max=10, step=0.5, description='Crack %'),
    fraction_gb=widgets.FloatSlider(value=0.02, min=0, max=20, step=1.0, description='GB %'),
    gamma=widgets.FloatLogSlider(value=0.1, base=10, min=-2, max=0, step=0.5, description='γ (crack)'),
    delta=widgets.FloatLogSlider(value=100, base=10, min=0, max=4, step=0.5, description='δ (GB)'),
)
def interactive_L36_solver(P_down, L_m, k_diss, K_eq,k_diss_metal, K_eq_metal, use_sieverts_pinhole,
                           D_ox, K_ox, L_ox, D_m, K_s_m,
                           fraction_pinhole, fraction_crack, fraction_gb, gamma, delta):
    """
    Interactive solver for L3+L6 parallel path model.
    Single-pass loop: plot arrays and DataFrame rows built together from the
    same computed values, so Rate-Lim is identical in both by construction.
    """
    
    # =========================================================================
    # Setup
    # =========================================================================
    # Convert percentages to fractions
    fraction_pinhole_frac = fraction_pinhole / 100
    fraction_crack_frac = fraction_crack / 100
    fraction_gb_frac = fraction_gb / 100
    
    # Check total defect fraction
    total_defect = fraction_pinhole_frac + fraction_crack_frac + fraction_gb_frac
    if total_defect > 1.0:
        print(f"⚠️ Total defect fraction ({total_defect*100:.1f}%) exceeds 100%!")
        return
    
    # Build defect config with converted fractions
    defect_config = {}
    if fraction_pinhole_frac > 0:
        defect_config['pinhole'] = {'area_fraction': fraction_pinhole_frac}
    if fraction_crack_frac > 0:
        defect_config['crack'] = {'area_fraction': fraction_crack_frac, 'thickness_factor': gamma}
    if fraction_gb_frac > 0:
        defect_config['grain_boundary'] = {'area_fraction': fraction_gb_frac, 'diffusivity_factor': delta}
    

    # Pressure range for sweep
    P_up_range = np.logspace(0, 12, 50)
    alpha_intact = D_ox * K_ox / L_ox


        
    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    plot_data = {k: [] for k in [
        'J_total', 'J_intact', 'enhancement', 'theta_intact',
        'frac_intact', 'frac_pinhole', 'frac_crack', 'frac_gb',
        'fraction_surface', 'fraction_oxide', 'fraction_metal',
    ]}
    rows = []
    
    for P_up in P_up_range:
        try:
            # Pass metal kinetics to calculation
            if defect_config:
                r = calculate_mixed_defect_flux_L6(
                    P_up, P_down, L_m, k_diss, K_eq,
                    D_ox, K_ox, L_ox, D_m, K_s_m,
                    defect_config=defect_config,
                    k_diss_metal=None if use_sieverts_pinhole else k_diss_metal,
                    K_eq_metal=None if use_sieverts_pinhole else K_eq_metal,
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
                ws += w * res_intact['fraction_surface']
                wo += w * res_intact['fraction_oxide']
                wm += w * res_intact['fraction_metal']
                
                # Defect path contributions
                for dt, data in r['defect_paths'].items():
                    pr = data['path_result']
                    J_c = data['flux_contribution']
                    w = J_c / J_tot if J_tot > 0 else 0.0
                    ws += w * pr['resistances']['fraction_surface']
                    wo += w * pr['resistances']['fraction_oxide']
                    wm += w * pr['resistances']['fraction_metal']
                
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
                plot_data['fraction_surface'].append(ws)
                plot_data['fraction_oxide'].append(wo)
                plot_data['fraction_metal'].append(wm)
                
                # Append to DataFrame rows (same values, no recomputation)
                rows.append({
                    "P_up (Pa)":           P_up,
                    "J_total (mol/m²/s)":  J_tot,
                    "J_intact (mol/m²/s)": J_int,
                    "Enhancement":         enhancement,
                    "θ_intact":            theta,
                    "P_int_intact (Pa)":   P_int,
                    "fraction_surface (%)": ws * 100,
                    "fraction_oxide (%)":   wo * 100,
                    "fraction_metal (%)":   wm * 100,
                    "pinhole_flux_fraction (%)": frac_pinhole * 100,
                    "fraction_intact (%)":        frac_intact * 100,
                    "fraction_pinhole (%)":       frac_pinhole * 100,
                    "fraction_crack (%)":         frac_crack * 100,
                    "fraction_gb (%)":            frac_gb * 100,
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
                ws = r['resistances']['fraction_surface']
                wo = r['resistances']['fraction_oxide']
                wm = r['resistances']['fraction_metal']
                
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
                plot_data['fraction_surface'].append(ws)
                plot_data['fraction_oxide'].append(wo)
                plot_data['fraction_metal'].append(wm)
                
                # Append to DataFrame rows
                rows.append({
                    "P_up (Pa)":           P_up,
                    "J_total (mol/m²/s)":  J_ss,
                    "J_intact (mol/m²/s)": J_ss,
                    "Enhancement":         1.0,
                    "θ_intact":            theta,
                    "P_int_intact (Pa)":   P_int,
                    "fraction_intact (%)":        100.0,
                    "fraction_pinhole (%)":       0.0,
                    "fraction_crack (%)":         0.0,
                    "fraction_gb (%)":            0.0,
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
    fraction_surface = np.array(plot_data['fraction_surface'])
    fraction_oxide = np.array(plot_data['fraction_oxide'])
    fraction_metal = np.array(plot_data['fraction_metal'])
    frac_intact = np.array(plot_data['frac_intact'])
    frac_pinhole = np.array(plot_data['frac_pinhole'])
    frac_crack = np.array(plot_data['frac_crack'])
    frac_gb_arr = np.array(plot_data['frac_gb'])
    
    # Vectorised rate_limiting_arr for plot shading — derived from the same
    # fraction_surface/fraction_oxide/fraction_metal arrays written inside the loop
    rate_limiting_arr = np.where(
        fraction_surface > 0.5, 'surface',
        np.where(fraction_oxide > 0.5, 'oxide',
        np.where(fraction_metal > 0.5, 'metal',
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
    #Net slope
    logP_net = np.log10(P_up_range)
    logJ_net = np.log10(np.abs(J_total))
    slope_net, _ = np.polyfit(logP_net, logJ_net, 1)
    # Intact path slope
    logP_intact = np.log10(P_up_range)
    logJ_intact = np.log10(np.abs(J_intact))
    slope_intact, _ = np.polyfit(logP_intact, logJ_intact, 1)
    
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



    ax1.text(0.98, 0.02, f'Net_slope = {slope_net:.2f}', transform=ax1.transAxes,ha='right', va='bottom',
    fontsize=10,fontweight='bold',color='black',bbox=dict(boxstyle='square', fc='wheat', ec='gray', alpha=1))
    ax1.text(0.98, 0.10, f'Intact_slope = {slope_intact:.2f}', transform=ax1.transAxes,ha='right', va='bottom',
    fontsize=10,fontweight='bold',color='blue',bbox=dict(boxstyle='square', fc='lightblue', ec='blue', alpha=1))
    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=12)
    ax1.set_ylabel('Steady-State Flux $J_{ss}$ (mol/m²/s)', fontsize=12)
    ax1.set_title('L3+L6: Flux vs Pressure (Defective Oxide)', fontsize=14)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')

    # ===== Plot 2: Rate-Limiting Fractions =====
    ax2.semilogx(P_up_range, fraction_surface * 100, 'r-', linewidth=2, label='Surface (dissociation)')
    ax2.semilogx(P_up_range, fraction_oxide * 100, 'orange', linewidth=2, label='Oxide (diffusion)')
    ax2.semilogx(P_up_range, fraction_metal * 100, 'b-', linewidth=2, label='Metal (diffusion)')
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
    - Grain boundary enhancement: D_eff = fraction_gb × D_lattice
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
    tolerance=1e-5,
    # NEW: Metal surface kinetics for pinhole
    k_diss_metal=None,
    K_eq_metal=None,
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
    sqrt_P_up = np.sqrt(max(P_up, 0))
    
    # Determine kinetics for pinhole
    if is_pinhole:
        if k_diss_metal is not None and K_eq_metal is not None:
            k_diss_eff = k_diss_metal
            K_eq_eff = K_eq_metal
            use_sieverts_limit = False
        else:
            use_sieverts_limit = True
    else:
        k_diss_eff = k_diss
        K_eq_eff = K_eq
        use_sieverts_limit = False
    # =========================================================================
    # Iterative solution: D_eff and θ are coupled
    # =========================================================================
    D_eff = D_lattice  # Initial guess
    convergence_history = []
    
    for iteration in range(max_iterations):
        beta = D_eff * K_s_m / L_m
        
        # ---------------------------------------------------------------------
        # Handle Sieverts' limit for pinhole (fast metal kinetics)
        # ---------------------------------------------------------------------
        if is_pinhole and use_sieverts_limit:
            # No surface resistance: √P_int = √P_up
            sqrt_P_int = sqrt_P_up
            theta_ss = 0.0  # Fast kinetics = no coverage
            
        else:
            # ---------------------------------------------------------------------
            # Define residual for θ (surface flux = metal flux)
            # ---------------------------------------------------------------------
            def residual(theta):
                if theta <= 0 or theta >= 1:
                    return np.inf
                
                J_surf = surface_flux(theta, P_up, k_diss_eff, K_eq_eff)
                
                if is_pinhole:
                    # No oxide: √P_int = g(θ) with appropriate K_eq
                    sqrt_P_int_local = g_theta(theta, K_eq_eff)
                else:
                    # With oxide: solve flux balance
                    sqrt_P_int_local = sqrt_P_int_from_theta(theta, alpha, beta, K_eq, P_down)
                
                J_metal = beta * (sqrt_P_int_local - sqrt_P_down)
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
            
            # Calculate P_int from θ
            if is_pinhole:
                sqrt_P_int = g_theta(theta_ss, K_eq_eff)
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
    
    if is_pinhole and use_sieverts_limit:
        J_surface = J_metal  # By definition at steady state
        J_oxide = np.nan
        R_oxide = 0.0
        R_surface = 0.0  # Fast kinetics = no resistance
        alpha_out = np.inf
        kinetics_used = 'sieverts'
    elif is_pinhole:
        J_surface = surface_flux(theta_ss, P_up, k_diss_eff, K_eq_eff)
        J_oxide = np.nan
        R_oxide = 0.0
        #R_surface = 1.0 / (k_diss_eff * P_up) if theta_ss < 1 else 0.0
        R_surface = smooth_surface_resistance(k_diss_eff, P_up, theta_ss) #PATCH
        alpha_out = np.inf
        kinetics_used = 'metal'
    else:
        J_surface = surface_flux(theta_ss, P_up, k_diss, K_eq)
        J_oxide = oxide_flux(theta_ss, alpha, beta_eff, K_eq, P_down)
        R_oxide = 1.0 / alpha
        #R_surface = 1.0 / (k_diss * P_up) if theta_ss < 1 else 0.0
        R_surface = smooth_surface_resistance(k_diss, P_up, theta_ss) #PATCH
        alpha_out = alpha
        kinetics_used = 'oxide'
    
    # =========================================================================
    # Rate-limiting analysis
    # =========================================================================
    R_metal = 1.0 / beta_eff
    R_total = R_surface + R_oxide + R_metal
    
    fraction_surface = R_surface / R_total if R_total > 0 else 0
    fraction_oxide = R_oxide / R_total if R_total > 0 else 0
    fraction_metal = R_metal / R_total if R_total > 0 else 0
    
    if fraction_surface > 0.5:
        rate_limiting = 'surface'
    elif fraction_oxide > 0.5:
        rate_limiting = 'oxide'
    elif fraction_metal > 0.5:
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
        'kinetics_used': kinetics_used,
        
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
            'fraction_surface': fraction_surface,
            'fraction_oxide': fraction_oxide,
            'fraction_metal': fraction_metal,
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
    mode='both',
    # NEW: Metal surface kinetics for pinhole
    k_diss_metal=None,
    K_eq_metal=None,
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
        Surface kinetics parameters (oxide surface)
    D_ox, K_ox, L_ox : float
        Intact oxide properties
    D_lattice, K_s_m : float
        Metal lattice properties
    microstructure_params : dict
        Metal microstructure (grain size, traps, etc.)
    defect_config : dict
        Oxide defect configuration
    k_diss_metal : float, optional
        Metal surface dissociation rate [mol/m²/s/Pa]. For pinhole only.
    K_eq_metal : float, optional
        Metal surface equilibrium constant [Pa⁻¹]. For pinhole only.
    
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
    
    fraction_intact = 1.0 - total_defect_fraction
    
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
        # No metal kinetics for intact path (uses oxide kinetics)
    )
    
    # =========================================================================
    # Calculate flux through each defect path
    # =========================================================================
    defect_results = {}
    
    for defect_type, config in defect_config.items():
        fraction_defect = config.get('area_fraction', 0.0)
        
        if fraction_defect <= 0:
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
        # Pass metal kinetics ONLY for pinhole paths
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
            mode=mode,
            # NEW: Pass metal kinetics for pinhole only
            k_diss_metal=k_diss_metal if defect_type == 'pinhole' else None,
            K_eq_metal=K_eq_metal if defect_type == 'pinhole' else None,
        )
        
        defect_results[defect_type] = {
            'area_fraction': fraction_defect,
            'alpha': alpha_defect,
            'alpha_ratio': alpha_defect / alpha_intact if alpha_defect != np.inf else np.inf,
            'path_result': path_result,
            'flux_contribution': fraction_defect * path_result['flux'],
        }
    
    # =========================================================================
    # Calculate total flux (area-weighted sum)
    # =========================================================================
    J_intact_contribution = fraction_intact * intact_result['flux']
    J_total = J_intact_contribution
    
    for defect_type, data in defect_results.items():
        J_total += data['flux_contribution']
    
    # =========================================================================
    # Flux breakdown analysis
    # =========================================================================
    flux_breakdown = {
        'intact': {
            'area_fraction': fraction_intact,
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
            'kinetics_used': pr.get('kinetics_used', 'oxide'),  # NEW: Track kinetics type
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
        'fraction_intact': fraction_intact,
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
    # Operating conditions
    P_down      = widgets.FloatSlider(value=0,   min=0,  max=10000,  step=100, description='P_down (Pa)'),
    L_m         = widgets.FloatLogSlider(value=1e-3,  base=10, min=-4,  max=-1, step=0.5, description='L_m (m)'),
    temperature = widgets.IntSlider(     value=773,            min=473, max=1273, step=50, description='T (K)'),
    # Surface kinetics (oxide)
    k_diss      = widgets.FloatLogSlider(value=1e-15, base=10, min=-18, max=-10, step=0.5, description='k_diss'),
    K_eq        = widgets.FloatLogSlider(value=1e-10, base=10, min=-15, max=-5,  step=0.5, description='K_eq'),
    # Metal surface kinetics for pinhole (NEW)
    k_diss_metal    = widgets.FloatLogSlider(value=1e-12, base=10, min=-15, max=-8, step=0.5, description='k_diss_metal'),
    K_eq_metal      = widgets.FloatLogSlider(value=1e-8,  base=10, min=-12, max=-4, step=0.5, description='K_eq_metal'),
    use_sieverts_pinhole = widgets.Checkbox(value=False, description='Sieverts limit for pinhole'),
    # Oxide properties
    D_ox        = widgets.FloatLogSlider(value=1e-11, base=10, min=-22, max=-12, step=0.5, description='D_ox (m²/s)'),
    K_ox        = widgets.FloatLogSlider(value=1e-6,  base=10, min=-6,  max=0,   step=0.5, description='K_ox'),
    L_ox        = widgets.FloatLogSlider(value=1e-6,  base=10, min=-8,  max=-4,  step=0.5, description='L_ox (m)'),
    # Metal properties
    D_lattice   = widgets.FloatLogSlider(value=1e-12, base=10, min=-14, max=-6,  step=0.5, description='D_lattice'),
    K_s_m       = widgets.FloatLogSlider(value=3.16e-4,  base=10, min=-4,  max=0,   step=0.5, description='K_s_m'),
    # Oxide defect parameters
    fraction_pinhole   = widgets.FloatSlider(value=0.01, min=0, max=2,  step=0.05, description='Pinhole %'),
    fraction_crack     = widgets.FloatSlider(value=0.05, min=0, max=5,  step=0.1,  description='Crack %'),
    fraction_gb_ox     = widgets.FloatSlider(value=0.02, min=0, max=10, step=0.5,  description='Oxide GB %'),
    gamma       = widgets.FloatLogSlider(value=0.1,   base=10, min=-2,  max=0,   step=0.5, description='γ (crack)'),
    delta       = widgets.FloatLogSlider(value=100,   base=10, min=0,   max=4,   step=0.5, description='δ (ox GB)'),
    # Metal microstructure parameters
    grain_size  = widgets.FloatLogSlider(value=31e-6, base=10, min=-6,  max=-3,  step=0.5, description='Grain (m)'),
    trap_density= widgets.FloatLogSlider(value=3.16e15,  base=10, min=12,  max=17,  step=0.5, description='ρ_trap (m⁻²)'),
    E_trap      = widgets.FloatSlider(   value=27,             min=15,  max=60,  step=5,   description='E_trap (kJ/mol)'),
)
def interactive_L346_full_model(P_down, L_m, temperature, k_diss, K_eq,
                                 k_diss_metal, K_eq_metal, use_sieverts_pinhole,
                                 D_ox, K_ox, L_ox, D_lattice, K_s_m,
                                 fraction_pinhole, fraction_crack, fraction_gb_ox, gamma, delta,
                                 grain_size, trap_density, E_trap):
    """
    Interactive widget for the full L3+L4+L6 model.
    Defective Oxide + Defective Metal + Surface Chemistry.
    
    Now includes:
    - Metal surface kinetics option for pinhole paths
    - Consistent handling with L36 solver
    """

    # =========================================================================
    # Setup
    # =========================================================================
    # Convert percentages to fractions
    fraction_pinhole_frac = fraction_pinhole / 100
    fraction_crack_frac   = fraction_crack   / 100
    fraction_gb_ox_frac   = fraction_gb_ox   / 100

    total_defect = fraction_pinhole_frac + fraction_crack_frac + fraction_gb_ox_frac
    if total_defect > 1.0:
        print(f"⚠️ Total defect fraction ({total_defect*100:.1f}%) exceeds 100%!")
        return

    # Build defect config
    defect_config = {}
    if fraction_pinhole_frac > 0:
        defect_config['pinhole'] = {'area_fraction': fraction_pinhole_frac}
    if fraction_crack_frac > 0:
        defect_config['crack'] = {'area_fraction': fraction_crack_frac, 'thickness_factor': gamma}
    if fraction_gb_ox_frac > 0:
        defect_config['grain_boundary'] = {'area_fraction': fraction_gb_ox_frac, 'diffusivity_factor': delta}

    # Metal microstructure config
    microstructure_params = {
        'grain_size':  grain_size,
        'grain_shape': 'equiaxed',
        'gb_type':     'HAGB',
        'trap_list': [
            {'name': 'dislocations', 'density': trap_density, 'binding_energy': E_trap * 1000}
        ] if trap_density > 1e12 else []
    }

    # Baseline microstructure (no defects)
    microstructure_lattice = {
        'grain_size':  1e-3,
        'grain_shape': 'equiaxed',
        'gb_type':     'LAGB',
        'trap_list':   []
    }

    # Metal kinetics for pinhole (None = Sieverts limit)
    k_diss_metal_eff = None if use_sieverts_pinhole else k_diss_metal
    K_eq_metal_eff = None if use_sieverts_pinhole else K_eq_metal

    P_up_range   = np.logspace(1, 8, 40)
    alpha_intact = D_ox * K_ox / L_ox

    # =========================================================================
    # Single-pass loop — builds plot arrays and DataFrame rows together
    # =========================================================================
    plot_data = {k: [] for k in [
        'J_full', 'J_no_micro', 'J_no_defects', 'J_baseline',
        'frac_intact', 'frac_pinhole', 'frac_crack', 'frac_gb',
        'D_mod_factor', 'theta_intact',
        'fraction_surface', 'fraction_oxide', 'fraction_metal',
    ]}
    rows = []

    for P_up in P_up_range:
        try:
            # -----------------------------------------------------------------
            # Primary computation with defects
            # -----------------------------------------------------------------
            if defect_config:
                r = calculate_full_model_flux_L346_v2(
                    P_up, P_down, L_m, temperature,
                    k_diss, K_eq, D_ox, K_ox, L_ox,
                    D_lattice, K_s_m,
                    microstructure_params, defect_config,
                    n_points=10,
                    k_diss_metal=k_diss_metal_eff,
                    K_eq_metal=K_eq_metal_eff,
                )

                # Flux-weighted resistances across ALL paths
                J_tot = r['J_total']
                ws = wo = wm = 0.0
                
                # Intact path contribution
                frac_intact_val = r['flux_breakdown']['intact']['fraction_of_total']
                J_intact_contrib = frac_intact_val * J_tot
                res_intact = r['intact_path']['resistances']
                w = J_intact_contrib / J_tot if J_tot > 0 else 0.0
                ws += w * res_intact['fraction_surface']
                wo += w * res_intact['fraction_oxide']
                wm += w * res_intact['fraction_metal']
                
                # Defect path contributions
                for dt, data in r['defect_paths'].items():
                    pr = data['path_result']
                    J_c = data['flux_contribution']
                    w = J_c / J_tot if J_tot > 0 else 0.0
                    ws += w * pr['resistances']['fraction_surface']
                    wo += w * pr['resistances']['fraction_oxide']
                    wm += w * pr['resistances']['fraction_metal']

                # Rate-limiting — single evaluation shared by plot and DataFrame
                if   ws > 0.5: rate_lim = 'surface'
                elif wo > 0.5: rate_lim = 'oxide'
                elif wm > 0.5: rate_lim = 'metal'
                else:          rate_lim = 'mixed'

                # Plot scalars
                plot_data['J_full'].append(r['J_total'])
                plot_data['D_mod_factor'].append(r['overall_modification_factor'])
                plot_data['theta_intact'].append(r['intact_path']['theta'])
                plot_data['frac_intact'].append(frac_intact_val)
                plot_data['frac_pinhole'].append(r['flux_breakdown'].get('pinhole', {}).get('fraction_of_total', 0))
                plot_data['frac_crack'].append(r['flux_breakdown'].get('crack', {}).get('fraction_of_total', 0))
                plot_data['frac_gb'].append(r['flux_breakdown'].get('grain_boundary', {}).get('fraction_of_total', 0))
                plot_data['fraction_surface'].append(ws)
                plot_data['fraction_oxide'].append(wo)
                plot_data['fraction_metal'].append(wm)

                # DataFrame row
                rows.append({
                    "P_up (Pa)":          P_up,
                    "J_total (mol/m²/s)": r['J_total'],
                    "J_intact":           r['intact_path']['flux'],
                    "Enhancement":        r['enhancement_vs_intact'],
                    "θ_intact":           r['intact_path']['theta'],
                    "P_int (Pa)":         r['intact_path']['P_int'],
                    "D_eff/D_L":          r['overall_modification_factor'],
                    "fraction_intact (%)":       frac_intact_val * 100,
                    "fraction_pinhole (%)":      r['flux_breakdown'].get('pinhole', {}).get('fraction_of_total', 0) * 100,
                    "fraction_crack (%)":        r['flux_breakdown'].get('crack', {}).get('fraction_of_total', 0) * 100,
                    "fraction_gb (%)":           r['flux_breakdown'].get('grain_boundary', {}).get('fraction_of_total', 0) * 100,
                    "Dominant":           r['dominant_path'].upper(),
                    "Rate-Lim":           rate_lim.upper(),
                })

                # Comparison: no microstructure
                r_no_micro = calculate_full_model_flux_L346_v2(
                    P_up, P_down, L_m, temperature,
                    k_diss, K_eq, D_ox, K_ox, L_ox,
                    D_lattice, K_s_m,
                    microstructure_lattice, defect_config,
                    n_points=10,
                    k_diss_metal=k_diss_metal_eff,
                    K_eq_metal=K_eq_metal_eff,
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

                ws = r['resistances']['fraction_surface']
                wo = r['resistances']['fraction_oxide']
                wm = r['resistances']['fraction_metal']

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
                plot_data['fraction_surface'].append(ws)
                plot_data['fraction_oxide'].append(wo)
                plot_data['fraction_metal'].append(wm)

                rows.append({
                    "P_up (Pa)":          P_up,
                    "J_total (mol/m²/s)": r['flux'],
                    "J_intact":           r['flux'],
                    "Enhancement":        1.0,
                    "θ_intact":           r['theta'],
                    "P_int (Pa)":         r['P_int'],
                    "D_eff/D_L":          r['modification_factor'],
                    "fraction_intact (%)":       100.0,
                    "fraction_pinhole (%)":      0.0,
                    "fraction_crack (%)":        0.0,
                    "fraction_gb (%)":           0.0,
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
    fraction_surface    = np.array(plot_data['fraction_surface'])
    fraction_oxide      = np.array(plot_data['fraction_oxide'])
    fraction_metal      = np.array(plot_data['fraction_metal'])
    frac_intact  = np.array(plot_data['frac_intact'])
    frac_pinhole = np.array(plot_data['frac_pinhole'])
    frac_crack   = np.array(plot_data['frac_crack'])
    frac_gb      = np.array(plot_data['frac_gb'])
    D_mod_factor = np.array(plot_data['D_mod_factor'])

    # Vectorised rate_limiting_arr for plot shading
    rate_limiting_arr = np.where(
        fraction_surface > 0.5, 'surface',
        np.where(fraction_oxide > 0.5, 'oxide',
        np.where(fraction_metal > 0.5, 'metal',
                 'mixed'))
    )

    # =========================================================================
    # Plot 1: Flux vs Pressure with Rate-Limiting Regions (Top Left)
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax1 = axes[0, 0]
    valid_idx = ~np.isnan(J_full)

    ax1.loglog(P_up_range, J_full,     'k-',   linewidth=2.5, label='Full L3+L4+L6')
    ax1.loglog(P_up_range, J_baseline, 'b', linewidth=1.5, linestyle=':', alpha=0.7, label='Baseline (perfect)')
    #slope_net
    logP_net = np.log10(P_up_range)
    logJ_net = np.log10(np.abs(J_full))
    slope_net, _ = np.polyfit(logP_net, logJ_net, 1)

    #slope_baseline
    logP_base = np.log10(P_up_range)
    logJ_base = np.log10(np.abs(J_baseline))
    slope_base, _ = np.polyfit(logP_base, logJ_base, 1)

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

    ax1.text(0.75, 0.05, f"Slope (net) = {slope_net:.2f}", transform=ax1.transAxes,
             fontsize=10, fontweight ='bold', verticalalignment='bottom', horizontalalignment='left',
             bbox=dict(boxstyle='square', facecolor='wheat', alpha=0.8))

    ax1.text(0.75, 0.15, f"Slope (net) = {slope_base:.2f}", transform=ax1.transAxes,
             fontsize=10, fontweight ='bold', verticalalignment='bottom', horizontalalignment='left',
             bbox=dict(boxstyle='square', facecolor='wheat', alpha=0.8))
    ax1.set_xlabel('Upstream Pressure $P_{up}$ (Pa)', fontsize=11)
    ax1.set_ylabel('Flux $J$ (mol/m²/s)', fontsize=11)
    ax1.set_title('Flux vs Pressure (Full L3+L4+L6)', fontsize=12, fontweight='bold')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left', fontsize=8)

    # =========================================================================
    # Plot 2: Rate-Limiting Fractions (Top Right)
    # =========================================================================
    ax2 = axes[0, 1]
    ax2.semilogx(P_up_range, fraction_surface * 100, 'r-',      linewidth=2, label='Surface')
    ax2.semilogx(P_up_range, fraction_oxide   * 100, 'orange',  linewidth=2, label='Oxide')
    ax2.semilogx(P_up_range, fraction_metal   * 100, 'b-',      linewidth=2, label='Metal')
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
                 f'  Pinhole: {fraction_pinhole:.2f}%\n'
                 f'  Crack: {fraction_crack:.2f}% (γ={gamma:.2f})\n'
                 f'  GB: {fraction_gb_ox:.2f}% (δ={delta:.0f})\n'
                 f'Metal Microstructure:\n'
                 f'  Grain: {grain_size*1e6:.1f} μm\n'
                 f'  Traps: {trap_density:.1e} m⁻²\n'
                 f'Pinhole kinetics: {"Sieverts" if use_sieverts_pinhole else "Metal surface"}')
    ax3.text(0.98, 0.02, info_text, transform=ax3.transAxes, fontsize=8,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # =========================================================================
    # Plot 4: Flux Contribution by Path (Bottom Right)
    # =========================================================================
    ax4 = axes[1, 1]
    ax4.semilogx(P_up_range, frac_intact * 100, 'k-', linewidth=2, label='Intact')
    if fraction_pinhole > 0:
        ax4.semilogx(P_up_range, frac_pinhole * 100, 'r-',      linewidth=2, label='Pinhole')
    if fraction_crack > 0:
        ax4.semilogx(P_up_range, frac_crack   * 100, 'orange',  linewidth=2, label='Crack')
    if fraction_gb_ox > 0:
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
