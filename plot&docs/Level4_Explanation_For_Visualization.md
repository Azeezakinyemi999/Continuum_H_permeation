# ​Level 4: Defective Metal Model - Explanation for Visualization

## Overview

Level 4 models how **microstructural features in the metal** affect hydrogen transport. While Level 1 treats metal as a perfect crystal with uniform diffusion, real metals have:

1.  **Grain boundaries** - fast diffusion highways
2.  **Traps** - sites that capture and hold hydrogen (vacancies, dislocations, precipitates)

These two effects work in **opposite directions**:

-   Grain boundaries **enhance** diffusion (D↑)
-   Traps **reduce** effective diffusion (D↓)

---

## Physical Picture

```
Level 1 (Perfect Metal):┌─────────────────────────────────────────┐│  ○ → → → → → → → → → → → → → → → → ○  │  Uniform diffusion│  ○ → → → → → → → → → → → → → → → → ○  │  D = D_lattice│  ○ → → → → → → → → → → → → → → → → ○  │└─────────────────────────────────────────┘Level 4 (Defective Metal):┌─────────────────────────────────────────┐│  ○ → → → ║→→→→→║ → → → → ○ → → → → ○  │  ║ = grain boundary (fast)│  ○ → ●●● ║→→→→→║ → → ●●● → → → → → ○  │  ●●● = trap cluster (slow)│  ○ → → → ║→→→→→║ → → → → → → ●●● → ○  │└─────────────────────────────────────────┘             ↑ Fast                ↑ Slow         GB highway           Trapped H
```

---

## Part A: Grain Boundary Enhancement

### The Physics

Grain boundaries are disordered regions between crystalline grains. The atomic disorder creates:

-   More free volume
-   Lower activation energy for diffusion
-   "Fast pipes" for hydrogen transport

### Key Parameters

Parameter

Symbol

Typical Value

Description

Grain size

d

10-100 μm

Average grain diameter

GB thickness

δ

0.5-1 nm

Width of disordered region

Enhancement factor

α

10-1000×

D_gb/D_bulk ratio

### The Math (Hart Equation)

The effective diffusivity is a volume-weighted average:

```
D_eff = (1 - f_gb) × D_bulk + f_gb × D_gb
```

Where:

-   **f_gb** = grain boundary volume fraction = 3δ/d (for cubic grains)
-   **D_gb** = grain boundary diffusivity = α × D_bulk

For small f_gb (typically ~10⁻⁴):

```
D_eff ≈ D_bulk × (1 + f_gb × (α - 1))D_eff ≈ D_bulk × (1 + 3δα/d)
```

### GB Enhancement Factor

The code calculates:

```python
enhancement = 1 + f_gb * (gb_diffusivity_ratio - 1)# Where:# f_gb = 3 * gb_thickness / grain_size# gb_diffusivity_ratio = exp((E_bulk - E_gb) / RT)
```

### Temperature Dependence

The enhancement ratio α depends on temperature:

```
α = D_gb/D_bulk = (D0_gb/D0_bulk) × exp((E_bulk - E_gb)/RT)
```

Since E_gb < E_bulk:

-   At **high T**: α is smaller (bulk catches up)
-   At **low T**: α is larger (GB advantage is bigger)

**Visualization implication**: GB enhancement matters more at low temperatures!

---

## Part B: Trapping (Oriani Model)

### The Physics

Traps are sites where hydrogen has **lower energy** than in normal lattice sites:

-   Hydrogen "falls into" the trap
-   Must overcome binding energy to escape
-   Creates a **delayed diffusion** effect

### Types of Traps

Trap Type

Binding Energy

Density

Character

Dislocations

20-30 kJ/mol

10²³-10²⁴ /m³

Weak, reversible

Vacancies

40-50 kJ/mol

10²²-10²³ /m³

Moderate

Interfaces

50-70 kJ/mol

Depends on microstructure

Strong

Precipitates

60-100 kJ/mol

10²¹-10²³ /m³

Very strong

### The Oriani Equilibrium Model (1970)

**Key assumption**: Hydrogen in traps is in local equilibrium with hydrogen in lattice sites.

The equilibrium constant for trap ↔ lattice exchange:

```
K = exp(E_b / RT)
```

Where E_b = binding energy (positive value = trap is lower energy)

### Trap Occupancy

The fraction of trap sites occupied:

```
θ_T = θ_L × K / (1 + θ_L × (K - 1))
```

Where:

-   θ_L = lattice site occupancy (typically << 1)
-   θ_T = trap site occupancy
-   K = equilibrium constant

For dilute hydrogen (θ_L << 1):

```
θ_T ≈ θ_L × K / (1 + θ_L × K)
```

### Effective Diffusivity with Trapping

The Oriani formula for effective diffusivity:

```
D_eff = D_lattice / (1 + Σ(N_T,i × K_i / N_L))
```

Where:

-   N_T,i = density of trap type i (traps/m³)
-   K_i = equilibrium constant for trap type i
-   N_L = density of lattice sites (~10²⁹ /m³ for metals)

### Simplification for Single Trap Type

```
D_eff = D_lattice / (1 + N_T × K / N_L)Let β = N_T × K / N_L (dimensionless trapping parameter)D_eff = D_lattice / (1 + β)
```

### Temperature Dependence of Trapping

Since K = exp(E_b/RT):

-   **High T**: K is smaller → less trapping → D_eff → D_lattice
-   **Low T**: K is larger → more trapping → D_eff << D_lattice

**Visualization implication**: Trapping matters more at low temperatures!

---

## Combined Model (Level 4 Full)

### How the Effects Combine

The code applies both effects sequentially:

```python
# Step 1: Calculate GB-enhanced diffusivityD_gb_enhanced = D_lattice * gb_enhancement_factor# Step 2: Apply trapping reductionD_eff = D_gb_enhanced / (1 + trapping_factor)
```

This gives:

```
D_eff = D_lattice × (1 + f_gb × (α - 1)) / (1 + N_T × K / N_L)
```

### The Modification Factor

The code reports a "modification factor" M:

```
M = D_eff / D_lattice = GB_enhancement / (1 + trapping_factor)
```

-   M > 1: GB enhancement dominates → faster than perfect crystal
-   M < 1: Trapping dominates → slower than perfect crystal
-   M = 1: Effects cancel (or no defects)

---

## Visualization Suggestions

### Panel 1: GB Enhancement Only

```
Perfect Crystal          vs          With Grain Boundaries┌────────────────┐                   ┌────────────────┐│ ○ → → → → → ○ │                   │ ○ →║→→→║→ → ○ ││ D = D_lattice  │                   │ D = D_eff > D_L│└────────────────┘                   └────────────────┘
```

-   Show flux increase
-   Label: "GB Enhancement: D_eff/D_L = 1 + 3δα/d"

### Panel 2: Trapping Only

```
Perfect Crystal          vs          With Traps┌────────────────┐                   ┌────────────────┐│ ○ → → → → → ○ │                   │ ○ → ●● → → → ○││ D = D_lattice  │                   │ D = D_eff < D_L│└────────────────┘                   └────────────────┘
```

-   Show flux decrease
-   Label: "Trapping: D_eff/D_L = 1/(1 + N_T×K/N_L)"

### Panel 3: Combined Effects

```
┌──────────────────────────────────┐│   ○ → →║→→→║→ ●● →║→→→║→ → ○   ││                                  ││   GB makes it faster             ││   Traps make it slower           ││   Net effect depends on T        │└──────────────────────────────────┘
```

### Panel 4: Temperature Dependence Plot

```
log(D_eff/D_L)     │     │    ╱ GB-only (D > D_L)     │   ╱   0 │──────────────────── Level 1 (D = D_L)     │          ╲     │           ╲ Trap-only (D < D_L)     │            ╲     │──────────────────────────→ 1000/T (K⁻¹)         High T        Low T
```

---

## Key Equations Summary

Effect

Formula

Temperature Dependence

GB enhancement

D_gb/D_bulk = α = exp((E_bulk - E_gb)/RT)

↑ at low T

GB volume fraction

f_gb = 3δ/d

Independent of T

Trap equilibrium

K = exp(E_b/RT)

↑ at low T

Trapping factor

β = N_T × K / N_L

↑ at low T

Combined D_eff

D_L × (1 + f_gb(α-1)) / (1 + β)

Competing effects

---

## Code References

### Main Functions in `defective_metal.py`:

1.  **`gb_enhancement_factor()`** - Calculates the GB enhancement multiplier
2.  **`trap_occupancy()`** - Calculates Oriani equilibrium trap filling
3.  **`calculate_effective_diffusivity_trapping()`** - Applies trapping reduction
4.  **`calculate_gb_enhanced_diffusivity()`** - Applies GB enhancement
5.  **`combined_microstructure_model()`** - Combines both effects

### Main Function in `permeation_calc.py`:

**`calculate_defective_metal_flux()`** - Uses the combined model to get flux:

```python
J = D_eff × (C_up - C_down) / L
```

---

## Physical Intuition Checklist

✅ Grain boundaries are "fast pipes" - more free volume, lower activation energy✅ Smaller grains = more GB area = more enhancement✅ Traps have lower energy than lattice - hydrogen prefers them✅ Stronger traps (higher E_b) = more trapping effect✅ Both effects are stronger at low temperature✅ Net effect depends on microstructure and temperature✅ M > 1 means faster than ideal (GB wins)✅ M < 1 means slower than ideal (traps win)