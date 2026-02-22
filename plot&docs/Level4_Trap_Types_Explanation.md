# Level 4: Trap Types Explanation for Visualization

## Overview

The model implements **4 distinct trap types**, each with its own:
- **Binding energy** (E_b) — how strongly it holds hydrogen
- **Density** — how many trap sites exist per unit volume
- **Physical origin** — what creates this trap type

All traps use the **same Oriani equilibrium math**, but differ in their parameters.

---

## The 4 Trap Types

### 1. Dislocations (Weak, Reversible Traps)

| Property | Value | Notes |
|----------|-------|-------|
| Binding energy | **15 kJ/mol** | Lowest — weakest trap |
| Capture radius | 1 nm | Size of strain field |
| Typical density | 10¹⁴ – 10¹⁶ m⁻³ | Depends on cold work |

**Physical picture:**
```
   _______________
  |   ⊥   ⊥   ⊥   |  ⊥ = edge dislocation
  |_______________|
  
  H atoms sit in the tensile strain field below the dislocation line
```

**How density varies:**
- Annealed: 10¹⁴ m⁻³
- 10% cold work: 5×10¹⁵ m⁻³
- 20% cold work: 10¹⁶ m⁻³

**Key physics:** Dislocations are **weak, reversible traps**. Hydrogen exchanges rapidly with lattice. Important for materials with plastic deformation history.

---

### 2. Grain Boundaries (Moderate Traps + Fast Diffusion Paths)

| Property | Value | Notes |
|----------|-------|-------|
| Binding energy | **20 kJ/mol** | Moderate strength |
| Sites per area | 10¹⁹ m⁻² | Per unit GB area |
| GB thickness | 0.5 nm | Width of disordered region |

**Physical picture:**
```
  Grain 1 ║ Grain 2
          ║
  ○ → → → ║⇒⇒⇒⇒⇒ → → ○
          ║
          ║ ← disordered GB region
```

**Dual role:**
- **Trap:** H binds to excess free volume at GB
- **Fast path:** D_gb >> D_bulk (this is the Level 4 enhancement effect)

**How density is calculated:**
The code converts grain size to GB trap density:
```
N_T,GB = (sites_per_area × GB_area) / volume
       = (10¹⁹ m⁻²) × (3/d) for cubic grains
```
Smaller grains → more GB area → more traps

---

### 3. Vacancies (Moderate Traps, Temperature-Dependent)

| Property | Value | Notes |
|----------|-------|-------|
| Binding energy | **20 kJ/mol** | Similar to GB |
| Formation energy | 140 kJ/mol | Cost to create vacancy |
| Max occupancy | 6 H/vacancy | Multiple H can occupy |

**Physical picture:**
```
  ○ ○ ○ ○ ○
  ○ ○ □ ○ ○   □ = vacancy (missing atom)
  ○ ○ ○ ○ ○   H sits in the empty space
```

**How density varies with temperature:**
The code calculates thermal equilibrium vacancy concentration:
```python
C_v = N_sites × exp(-E_formation / RT)
```

| Temperature | Equilibrium Concentration |
|-------------|--------------------------|
| 600°C | 10²⁰ m⁻³ |
| 800°C | 10²¹ m⁻³ |
| 1000°C | 5×10²² m⁻³ |
| Quenched from 1000°C | 10²³ m⁻³ |

**Key physics:** Vacancy concentration is **temperature-dependent**! Higher T → more vacancies → more trapping (but also higher K → weaker binding effect).

---

### 4. Carbide Precipitates (Strong, Irreversible Traps)

| Property | Value | Notes |
|----------|-------|-------|
| Binding energy | **12 kJ/mol*** | At carbide/matrix interface |
| Interface area | ~100 m²/m³ | Specific interface area |
| Typical density | 10²⁰ – 10²² m⁻³ | Depends on aging |

*Note: The data file shows 12 kJ/mol for testing; literature values are 60-90 kJ/mol for real precipitates.

**Physical picture:**
```
  ┌─────────────────────┐
  │ Matrix   ████ Carbide│
  │         ████        │
  │    H →→ ████ ← H    │  H binds at interface
  │         ████        │
  └─────────────────────┘
```

**How density varies with heat treatment:**
- Solution treated: 10²⁰ m⁻³
- Aged 700°C/100h: 5×10²¹ m⁻³
- Aged 700°C/1000h: 10²² m⁻³

**Key physics:** Carbides are **strong traps** — hydrogen doesn't easily escape. Important for aged/heat-treated materials.

---

## How the Oriani Math Works for Each Trap Type

### The Universal Formula

All trap types use the **SAME core Oriani formula**:

```
D_eff = D_lattice / (1 + Σᵢ βᵢ)

where βᵢ = (N_T,i / N_L) × Kᵢ = (N_T,i / N_L) × exp(E_b,i / RT)
```

**What differs between trap types:**
1. **E_b** (binding energy) → different for each trap type
2. **N_T** (density) → calculated differently for each trap type
3. **Physical interpretation** → different trapping mechanism

---

### Trap-Specific Implementation Details

#### 1. Dislocations — Direct Input

```python
# In trap_list:
{'name': 'dislocations', 'binding_energy': 15e3, 'density': N_T_disl}

# Oriani contribution:
K_disl = exp(15000 / (R * T))
β_disl = (N_T_disl / N_L) * K_disl
```

**N_T source:** Directly specified from material condition
- Annealed: `N_T = 1e14 m⁻³`
- Cold-worked 20%: `N_T = 1e16 m⁻³`

**Key feature:** Density is a **fixed input** based on processing history.

---

#### 2. Grain Boundaries — Calculated from Geometry

```python
# First, calculate GB trap density from grain size:
gb_density_result = grain_boundary_density(
    grain_size=d,           # e.g., 50e-6 m
    gb_thickness=0.5e-9,    # 0.5 nm
    sites_per_area=1e19,    # trap sites/m²
    grain_shape='equiaxed'
)
N_T_gb = gb_density_result['trap_density']  # = (3/d) × 1e19

# Then apply Oriani:
K_gb = exp(20000 / (R * T))
β_gb = (N_T_gb / N_L) * K_gb
```

**N_T source:** Calculated via stereology formula:
```
N_T,GB = S_v × sites_per_area = (3/d) × 1e19
```
where `S_v = 3/d` is the GB surface area per volume.

**Key feature:** Density **depends on microstructure** (grain size).

---

#### 3. Vacancies — Temperature-Dependent Density

```python
# First, calculate thermal vacancy concentration:
vac_result = vacancy_concentration(
    temperature=T,
    material='Incoloy800',
    condition='equilibrium'  # or 'quenched'
)
N_T_vac = vac_result['concentration']  # = N_sites × exp(-E_f/RT)

# Then apply Oriani:
K_vac = exp(20000 / (R * T))
β_vac = (N_T_vac / N_L) * K_vac
```

**N_T source:** Thermodynamic equilibrium:
```
C_v = N_sites × exp(-E_formation / RT)
```
where `E_formation = 140 kJ/mol` for Incoloy 800.

**Key feature:** Density **depends on temperature** (self-consistently)!

This creates a **double temperature dependence**:
- K increases at low T (stronger binding)
- N_T decreases at low T (fewer vacancies)
- Net effect depends on which dominates

---

#### 4. Carbide Precipitates — Heat Treatment Dependent

```python
# In trap_list:
{'name': 'carbide_precipitates', 'binding_energy': 60e3, 'density': N_T_carb}

# Oriani contribution:
K_carb = exp(60000 / (R * T))
β_carb = (N_T_carb / N_L) * K_carb
```

**N_T source:** Directly specified from heat treatment
- Solution treated: `N_T = 1e20 m⁻³`
- Aged 700°C/100h: `N_T = 5e21 m⁻³`
- Aged 700°C/1000h: `N_T = 1e22 m⁻³`

**Key feature:** Very high E_b makes K **extremely large**, so even moderate N_T gives strong trapping.

---

### Code Implementation in `calculate_effective_diffusivity_trapping()`

The actual loop that handles all trap types identically:

```python
R_gas = 8.314  # J/mol/K
trapping_term = 0.0

for trap in trap_list:
    # Same formula for ALL trap types:
    K_i = np.exp(trap['binding_energy'] / (R_gas * temperature))
    trapping_contribution_i = (trap['density'] / lattice_density) * K_i
    trapping_term += trapping_contribution_i

# Final effective diffusivity
D_eff = D_lattice / (1.0 + trapping_term)
```

**The key insight:** The Oriani math is **trap-agnostic**. What differs is how you **determine the inputs** (E_b and N_T) for each trap type.

---

### Example Calculation at T = 800 K:

| Trap Type | E_b (kJ/mol) | K = exp(E_b/RT) | N_T (m⁻³) | β = N_T×K/N_L |
|-----------|--------------|-----------------|-----------|---------------|
| Dislocations | 15 | 6.0 | 10¹⁵ | 5.7×10⁻¹⁴ |
| Grain boundaries | 20 | 20.2 | 10²³ | 1.9×10⁻⁵ |
| Vacancies | 20 | 20.2 | 10²¹ | 1.9×10⁻⁷ |
| Carbides | 60 | 8,100 | 10²¹ | 7.6×10⁻⁵ |

With N_L = 1.06×10²⁹ m⁻³ (lattice site density for FCC)

**Total trapping term:** Σβ ≈ 9.5×10⁻⁵

**Reduction factor:** D_eff/D_L = 1/(1 + Σβ) ≈ 0.99990

At **lower temperatures** or **higher trap densities**, the effect becomes much stronger!

---

## Code Structure Summary

```python
# In data/microstruture_parameters.py:
TRAP_PROPERTIES = {
    'dislocations': {'binding_energy': 15e3, 'density_range': {...}},
    'grain_boundaries': {'binding_energy': 20e3, 'sites_per_area': 1e19},
    'vacancies': {'binding_energy': 20e3, 'formation_energy': 140e3},
    'carbide_precipitates': {'binding_energy': 12e3, 'density_range': {...}}
}

# In calculations/defective_metal.py:
def calculate_effective_diffusivity_trapping(D_lattice, T, trap_list, ...):
    """
    trap_list = [
        {'name': 'dislocations', 'binding_energy': 15e3, 'density': 1e15},
        {'name': 'vacancies', 'binding_energy': 20e3, 'density': 1e21},
        ...
    ]
    """
    trapping_term = 0
    for trap in trap_list:
        K = exp(trap['binding_energy'] / (R * T))
        β = (trap['density'] / N_L) * K
        trapping_term += β
    
    D_eff = D_lattice / (1 + trapping_term)
    return D_eff
```

---

## Key Differences Between Trap Types

| Trap Type | Binding Strength | Density Control | Reversibility | Dominant When? |
|-----------|-----------------|-----------------|---------------|----------------|
| Dislocations | Weak (15 kJ/mol) | Cold work | Reversible | High strain |
| Grain Boundaries | Moderate (20 kJ/mol) | Grain size | Semi-reversible | Fine grains |
| Vacancies | Moderate (20 kJ/mol) | Temperature | Reversible | High T or quenched |
| Carbides | Strong (60-90 kJ/mol) | Heat treatment | Irreversible | Aged materials |

---

## Temperature Dependence

The equilibrium constant K = exp(E_b/RT) makes trapping **stronger at low temperatures**:

| Trap Type | K at 500 K | K at 800 K | K at 1000 K |
|-----------|------------|------------|-------------|
| Dislocations (15 kJ/mol) | 37 | 6.0 | 3.0 |
| GB/Vacancies (20 kJ/mol) | 122 | 20 | 11 |
| Carbides (60 kJ/mol)* | 1.7×10⁶ | 8,100 | 1,400 |

*Using literature value of 60 kJ/mol

**Implication:** At low T, strong traps (carbides) dominate. At high T, all traps are less effective.

---

## Physical Intuition Checklist

✅ Dislocations = strain field traps → weak, density depends on cold work
✅ Grain boundaries = disorder traps → dual role (trap + fast path)
✅ Vacancies = empty site traps → temperature-dependent concentration
✅ Carbides = interface traps → strong, density depends on aging
✅ All traps use same Oriani math: D_eff = D_L / (1 + Σ N_T×K/N_L)
✅ Lower T → larger K → stronger trapping
✅ Competing effects: more traps slow diffusion, but effect weakens at high T
