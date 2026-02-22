# Level 6: Surface Kinetics Implementation Documentation

## Table of Contents
1. [Executive Summary](#1-executive-summary)
2. [Physical Theory and Motivation](#2-physical-theory-and-motivation)
3. [Mathematical Formulation](#3-mathematical-formulation)
4. [Implementation Details](#4-implementation-details)
5. [Code Architecture](#5-code-architecture)
6. [Coverage Mode Analysis](#6-coverage-mode-analysis)
7. [Validation and Testing](#7-validation-and-testing)
8. [Integration with Model Hierarchy](#8-integration-with-model-hierarchy)
9. [References](#9-references)

---

## 1. Executive Summary

### What Was Done
Level 6 implements **surface kinetics for hydrogen dissociative adsorption** as the rate-limiting step in hydrogen permeation. This adds a physics layer capturing the molecular-to-atomic transition at metal surfaces.

### Why It Matters
In many real systems, especially at:
- Low temperatures (reduced thermal activation)
- Low pressures (fewer molecular collisions)
- High diffusion rates (bulk transport faster than surface supply)
- Oxidized or contaminated surfaces (reduced active sites)

...the **surface dissociation step** can become rate-limiting rather than bulk diffusion.

### Key Outcome
Three coverage modes implemented:
| Mode | Use Case | Physical Meaning |
|------|----------|------------------|
| `equilibrium` | Da >> 1 | Fast kinetics → Sieverts' Law recovered |
| `steady_state` | Da ~ 1 or Da << 1 | Flux = min(J_diss, J_Sieverts) |
| `forced` | Parametric studies | User-specified coverage for "what if" analysis |

### Critical Physics Insight
**At Langmuir equilibrium, the model exactly recovers Sieverts' Law** — this is mathematically exact, not an approximation. Surface limitation only manifests when coverage deviates from equilibrium.

---

## 2. Physical Theory and Motivation

### 2.1 The Dissociation Problem

Hydrogen molecules (H₂) cannot directly dissolve in metals. They must first **dissociate** into atomic hydrogen (H) at the surface:

$$\text{H}_2(\text{gas}) + 2S \xrightarrow{k_{diss}} 2\text{H(ads)}$$

where $S$ represents an empty surface adsorption site.

The reverse process is **recombinative desorption**:

$$2\text{H(ads)} \xrightarrow{k_{recomb}} \text{H}_2(\text{gas}) + 2S$$

### 2.2 Surface Coverage θ

The **surface coverage** θ represents the fraction of surface sites occupied by adsorbed hydrogen:

$$\theta = \frac{N_{occupied}}{N_{total}} \quad (0 \leq \theta \leq 1)$$

Coverage affects the available sites for new dissociation events. If θ → 1 (surface saturated), no more molecules can adsorb.

### 2.3 Rate Laws

**Dissociation flux** (requires two adjacent empty sites):
$$J_{diss} = k_{diss} \times P \times (1-\theta)^2$$

**Recombination flux** (requires two adjacent occupied sites):
$$J_{recomb} = k_{recomb} \times \theta^2$$

### 2.4 Why Level 1-5 Are Insufficient

Levels 1-5 assume:
- Instantaneous equilibration at the surface
- Sieverts' Law always applies: $C_{surface} = K_s \sqrt{P}$

This fails when:
- Surface kinetics are slow relative to bulk diffusion
- Surface contamination reduces active site density
- Low temperatures suppress dissociation activation

---

## 3. Mathematical Formulation

### 3.1 Langmuir Isotherm (Equilibrium Coverage)

At equilibrium, dissociation and recombination rates are equal:
$$k_{diss} \times P \times (1-\theta)^2 = k_{recomb} \times \theta^2$$

Defining the equilibrium constant:
$$K_{eq} = \frac{k_{diss}}{k_{recomb}}$$

Solving for θ:
$$\frac{\theta}{1-\theta} = \sqrt{K_{eq} \times P}$$

Therefore:
$$\boxed{\theta = \frac{\sqrt{K_{eq} \times P}}{1 + \sqrt{K_{eq} \times P}}}$$

**Coverage regimes:**
| Condition | Regime | θ approximation |
|-----------|--------|-----------------|
| $K_{eq} \times P \ll 1$ | Low coverage | $\theta \approx \sqrt{K_{eq} \times P}$ |
| $K_{eq} \times P \gg 1$ | High coverage | $\theta \approx 1$ |
| $K_{eq} \times P \sim 1$ | Intermediate | Use full formula |

### 3.2 Subsurface Concentration

The subsurface concentration (just below the surface) relates to coverage through:

$$C_{surface} = K_s \times \frac{\theta}{1-\theta} \times \frac{1}{\sqrt{K_{eq}}}$$

### 3.3 Critical Identity: Langmuir → Sieverts

**At Langmuir equilibrium**, substituting $\theta/(1-\theta) = \sqrt{K_{eq} \times P}$:

$$C_{surface} = K_s \times \sqrt{K_{eq} \times P} \times \frac{1}{\sqrt{K_{eq}}} = K_s \times \sqrt{P}$$

This **is Sieverts' Law**! This proves that at equilibrium coverage, surface kinetics parameters drop out entirely.

### 3.4 Damköhler Number

The **Damköhler number** compares surface reaction rate to bulk transport rate:

$$\boxed{Da = \frac{k_{diss} \times K_s \times L}{D}}$$

| Da Range | Regime | Physical Meaning |
|----------|--------|------------------|
| Da >> 1 (>10) | Diffusion-limited | Fast dissociation; Sieverts applies |
| Da << 1 (<0.1) | Surface-limited | Slow dissociation controls flux |
| Da ~ 1 | Mixed | Comparable resistances |

### 3.5 Temperature Dependence (Arrhenius)

Rate constants follow reference-temperature Arrhenius form:

$$k_{diss}(T) = k_{diss,ref} \times \exp\left(\frac{-E_{diss}}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right)$$

$$k_{recomb}(T) = k_{recomb,ref} \times \exp\left(\frac{-E_{recomb}}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right)$$

**Typical values (Incoloy 800, T_ref = 1073.15 K):**
| Parameter | Value | Units |
|-----------|-------|-------|
| $k_{diss,ref}$ | 1.70×10⁻³ | m⁴/(mol·s) |
| $k_{recomb,ref}$ | 1.55×10⁻¹¹ | m²/s |
| $E_{diss}$ | 15,000 | J/mol |
| $E_{recomb}$ | 80,000 | J/mol |

---

## 4. Implementation Details

### 4.1 Core Functions

#### `surface_equilibrium_coverage(pressure, k_diss, k_recomb)`
**Location:** [calculations/surface_kinetics.py](calculations/surface_kinetics.py#L286-L320)

```
Input: pressure (Pa), k_diss (mol/(m²·s·Pa)), k_recomb (m⁴/(mol·s))
Output: {theta, K_eq, sqrt_K_eq_P, regime}
```

Implements:
$$\theta = \frac{\sqrt{K_{eq} \times P}}{1 + \sqrt{K_{eq} \times P}}$$

#### `dissociation_flux(pressure, theta, k_diss)`
**Location:** [calculations/surface_kinetics.py](calculations/surface_kinetics.py#L323-L342)

Implements:
$$J_{diss} = k_{diss} \times P \times (1-\theta)^2$$

#### `calculate_damkohler_number(k_diss, D, K_s, thickness)`
**Location:** [calculations/surface_kinetics.py](calculations/surface_kinetics.py#L345-L365)

Implements:
$$Da = \frac{k_{diss} \times K_s \times L}{D}$$

Returns regime classification.

#### `calculate_surface_limited_flux(..., coverage_mode)`
**Location:** [calculations/surface_kinetics.py](calculations/surface_kinetics.py#L540-L822)

**Main function** integrating all surface kinetics. Supports three coverage modes.

### 4.2 Data Module

**Location:** [data/surface_kinetics_data.py](data/surface_kinetics_data.py)

Contains:
- `SURFACE_KINETICS` dictionary with material parameters
- `SURFACE_SITE_DENSITY` for crystallographic surfaces
- `get_surface_kinetics(material_name, temperature)` function

**Supported materials:**
- Incoloy800 (primary)
- Fe_alpha (BCC iron)
- Ni (FCC nickel)
- SS316L (austenitic stainless)

---

## 5. Code Architecture

### 5.1 File Structure

```
calculations/
├── surface_kinetics.py          # Core L6 functions
├── permeation_calc.py          # Level 1 base (imports from L6)
├── interface_solver.py         # L2+L6 integration
└── parallel_oxide_defect_paths.py  # L3+L6 integration

data/
├── surface_kinetics_data.py    # Material kinetics parameters
└── material_data.py            # Bulk properties

validation/
└── test_surface_kinetics.py    # Comprehensive unit tests
```

### 5.2 Function Call Hierarchy

```
calculate_oxide_metal_system_with_surface() [interface_solver.py]
    ├── solve_interface_pressure_with_surface()
    │   ├── flux_balance_equation_with_surface()
    │   │   ├── molecular_diffusion_flux()          [L2: oxide]
    │   │   └── calculate_metal_flux_with_surface() [L6: metal]
    │   │       └── calculate_surface_limited_flux()
    │   │           ├── surface_equilibrium_coverage()
    │   │           ├── dissociation_flux()
    │   │           └── calculate_damkohler_number()
    │   └── get_surface_kinetics()                  [data module]
    └── calculate_damkohler_number()
```

### 5.3 Integration Points

| Level | Integration File | Key Function |
|-------|-----------------|--------------|
| L2+L6 | interface_solver.py | `solve_interface_pressure_with_surface()` |
| L3+L6 | parallel_oxide_defect_paths.py | `calculate_parallel_path_flux_with_surface()` |
| L5+L6 | sensitivity_level1.py | `level5L6_model_wrapper()` |

---

## 6. Coverage Mode Analysis

### 6.1 Mode: `equilibrium` (Default)

**Usage:** When surface kinetics are fast (Da >> 1)

**Implementation:**
1. Calculate θ from Langmuir isotherm
2. Calculate $C_{surface} = K_s \sqrt{P}$ (Sieverts)
3. Apply Fick's Law for bulk flux

**Result:** Exactly recovers Level 1 (Sieverts' Law)

**When to use:**
- High temperatures (> 800°C for most metals)
- Clean, unoxidized surfaces
- Thick membranes (high L in Da denominator)

### 6.2 Mode: `steady_state`

**Usage:** When surface and bulk resistances are comparable (Da ~ 1)

**Implementation:**
1. Calculate θ from Langmuir (for reference)
2. Calculate both J_diss and J_Sieverts
3. Flux = min(J_diss, J_Sieverts)

**Physics:** The slower process limits overall flux.

**When to use:**
- Mixed regime (0.1 < Da < 10)
- Contaminated surfaces
- Low-temperature operation

### 6.3 Mode: `forced`

**Usage:** Parametric studies, sensitivity analysis

**Implementation:**
1. User specifies θ directly (0 ≤ θ < 1)
2. Calculate $C_{surface} = K_s \times \frac{\theta}{1-\theta} \times \frac{1}{\sqrt{K_{eq}}}$
3. Apply Fick's Law

**Applications:**
- "What if" scenarios
- Matching experimental θ measurements
- Sensitivity analysis of coverage effects

---

## 7. Validation and Testing

### 7.1 Test File
**Location:** [validation/test_surface_kinetics.py](validation/test_surface_kinetics.py)

### 7.2 Key Test Cases

| Test | Validates | Expected Result |
|------|-----------|-----------------|
| `test_equilibrium_coverage_zero_pressure` | θ = 0 at P = 0 | ✓ |
| `test_equilibrium_coverage_langmuir_limit` | θ → 1 at high K_eq×P | ✓ |
| `test_dissociation_flux_full_coverage` | J_diss = 0 at θ = 1 | ✓ |
| `test_fast_kinetics_approaches_sieverts` | SRF ≈ 1.0 at equilibrium | ✓ |
| `test_equilibrium_coverage_recovers_sieverts` | Langmuir → Sieverts identity | ✓ |
| `test_slow_kinetics_reduces_flux` | SRF < 1 with forced low θ | ✓ |

### 7.3 Critical Physics Validation

The test suite proves the fundamental identity:

$$\theta_{Langmuir} = \frac{\sqrt{K_{eq} P}}{1 + \sqrt{K_{eq} P}} \Rightarrow C_{surface} = K_s \sqrt{P}$$

Tests sweep k_diss over 20 orders of magnitude and verify SRF = 1.0 at Langmuir equilibrium for all values.

---

## 8. Integration with Model Hierarchy

### 8.1 Complete Model Stack (L5+L6)

```
Level 6: Surface Kinetics (dissociation rate-limiting)
    ↓
Level 5: Complete hierarchical model
    ├── Level 4: Metal microstructure (GB + trapping)
    ├── Level 3: Oxide defects (parallel paths)
    ├── Level 2: Oxide barrier layer
    └── Level 1: Base metal (Sieverts + Fick)
```

### 8.2 Sensitivity Analysis Integration

**Location:** [validation/sensitivity_level1.py](validation/sensitivity_level1.py#L1627)

The `level5L6_model_wrapper()` function:
- Accepts all L1-L4 parameters plus L6 surface kinetics
- Returns standard outputs plus: SRF, theta, Da
- Supports Morris and Sobol sensitivity analysis

**Key L6 parameters for sensitivity:**
- `E_diss`: Dissociation activation energy
- `E_recomb`: Recombination activation energy
- `coverage_mode`: 'equilibrium', 'steady_state', or 'forced'

### 8.3 Surface Reduction Factor (SRF)

The **Surface Reduction Factor** quantifies L6 impact:

$$SRF = \frac{J_{L6}}{J_{Sieverts}}$$

| SRF | Interpretation |
|-----|----------------|
| 1.0 | No surface effect (Sieverts recovered) |
| < 1 | Surface-limited flux reduction |
| > 1 | Non-physical (θ > θ_equilibrium) |

---

## 9. References

### 9.1 Primary Literature

1. **Pick, M.A. & Sonnenberg, K. (1985)**  
   "A model for atomic hydrogen-metal interactions - Application to recycling, recombination and permeation."  
   *J. Nucl. Mater.* 131, 208-220.  
   DOI: [10.1016/0022-3115(85)90459-3](https://doi.org/10.1016/0022-3115(85)90459-3)  
   *→ Foundational paper for Langmuir isotherm + permeation coupling*

2. **Baskes, M.I. (1980)**  
   "A calculation of the surface recombination rate constant for hydrogen isotopes on metals."  
   *J. Nucl. Mater.* 92, 318-324.  
   DOI: [10.1016/0022-3115(80)90117-8](https://doi.org/10.1016/0022-3115(80)90117-8)  
   *→ Recombination rate constants for various metals*

3. **Andrew, P.L. & Haasz, A.A. (1992)**  
   "Models for hydrogen permeation in metals."  
   *J. Appl. Phys.* 72, 2749-2757.  
   DOI: [10.1063/1.351526](https://doi.org/10.1063/1.351526)  
   *→ Review of permeation models including surface effects*

4. **Causey, R.A. (2002)**  
   "Hydrogen isotope retention and recycling in fusion reactor plasma-facing components."  
   *J. Nucl. Mater.* 300, 91-117.  
   DOI: [10.1016/S0022-3115(01)00732-2](https://doi.org/10.1016/S0022-3115(01)00732-2)  
   *→ Comprehensive review with kinetics data*

5. **Wampler, W.R. (1986)**  
   "Surface recombination of hydrogen on clean nickel."  
   *Appl. Phys. Lett.* 48, 405-407.  
   DOI: [10.1063/1.96521](https://doi.org/10.1063/1.96521)  
   *→ Experimental Ni kinetics data*

### 9.2 Code Files

| File | Purpose |
|------|---------|
| [calculations/surface_kinetics.py](calculations/surface_kinetics.py) | Core L6 functions |
| [data/surface_kinetics_data.py](data/surface_kinetics_data.py) | Material parameters |
| [calculations/interface_solver.py](calculations/interface_solver.py) | L2+L6 integration |
| [validation/test_surface_kinetics.py](validation/test_surface_kinetics.py) | Unit tests |
| [validation/sensitivity_level1.py](validation/sensitivity_level1.py) | Sensitivity analysis |

---

## Appendix A: Equation Derivation Details

### A.1 Langmuir Isotherm Derivation

Starting from rate equality at equilibrium:
$$k_{diss} \cdot P \cdot (1-\theta)^2 = k_{recomb} \cdot \theta^2$$

Rearranging:
$$\frac{(1-\theta)^2}{\theta^2} = \frac{k_{recomb}}{k_{diss} \cdot P} = \frac{1}{K_{eq} \cdot P}$$

Taking square root:
$$\frac{1-\theta}{\theta} = \frac{1}{\sqrt{K_{eq} \cdot P}}$$

Inverting:
$$\frac{\theta}{1-\theta} = \sqrt{K_{eq} \cdot P}$$

Solving for θ:
$$\theta = (1-\theta) \cdot \sqrt{K_{eq} \cdot P}$$
$$\theta = \sqrt{K_{eq} \cdot P} - \theta \cdot \sqrt{K_{eq} \cdot P}$$
$$\theta \cdot (1 + \sqrt{K_{eq} \cdot P}) = \sqrt{K_{eq} \cdot P}$$
$$\boxed{\theta = \frac{\sqrt{K_{eq} \cdot P}}{1 + \sqrt{K_{eq} \cdot P}}}$$

### A.2 Sieverts Recovery Proof

Given:
$$C_{surface} = K_s \cdot \frac{\theta}{1-\theta} \cdot \frac{1}{\sqrt{K_{eq}}}$$

At Langmuir equilibrium:
$$\frac{\theta}{1-\theta} = \sqrt{K_{eq} \cdot P}$$

Substituting:
$$C_{surface} = K_s \cdot \sqrt{K_{eq} \cdot P} \cdot \frac{1}{\sqrt{K_{eq}}}$$
$$C_{surface} = K_s \cdot \sqrt{\frac{K_{eq} \cdot P}{K_{eq}}}$$
$$\boxed{C_{surface} = K_s \cdot \sqrt{P}}$$

**Q.E.D.** — This is exactly Sieverts' Law.

### A.3 Damköhler Number Physical Interpretation

The Damköhler number represents:
$$Da = \frac{\text{Surface reaction rate}}{\text{Bulk transport rate}}$$

Derivation:
- Surface flux scale: $J_{surf} \sim k_{diss} \cdot K_s \cdot L \cdot P^{0.5}$
- Bulk flux scale: $J_{bulk} \sim D \cdot K_s \cdot P^{0.5} / L$

Ratio:
$$Da = \frac{J_{surf}}{J_{bulk}} = \frac{k_{diss} \cdot K_s \cdot L}{D}$$

---

*Document generated for Hierarchical Hydrogen Permeation Model*  
*Level 6: Surface Kinetics - Phases 1-3 Complete*  
*Last updated: 2025*
