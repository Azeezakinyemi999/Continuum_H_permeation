# Level 6 Roadmap: Phase 4 and Phase 5 Implementation

## Table of Contents
1. [Overview](#1-overview)
2. [Current State (Phases 1-3)](#2-current-state-phases-1-3)
3. [Phase 4: Co-Adsorbate Competition](#3-phase-4-co-adsorbate-competition)
4. [Phase 5: Dynamic Surface Reconstruction](#4-phase-5-dynamic-surface-reconstruction)
5. [Implementation Timeline](#5-implementation-timeline)
6. [Appendix: Equation Derivations](#appendix-equation-derivations)

---

## 1. Overview

### Complexity Progression

| Phase | Name | Status | Key Physics |
|-------|------|--------|-------------|
| 1 | Single-step kinetics | ✅ Complete | J = k × P |
| 2 | Dual-step (diss/recomb) | ✅ Complete | J_diss = k_diss × P × (1-θ)² |
| 3 | Langmuir coverage | ✅ Complete | θ = √(K_eq×P)/(1+√(K_eq×P)) |
| **4** | **Co-adsorbate competition** | 🔲 Planned | Multi-species Langmuir |
| **5** | **Dynamic reconstruction** | 🔲 Planned | Surface structure changes |

### Why These Phases Matter

**Phase 4 (Co-adsorbates):** Real industrial environments contain multiple gas species (H₂O, CO, CO₂, O₂, N₂) that compete for surface sites, reducing H₂ dissociation and permeation.

**Phase 5 (Reconstruction):** At high coverages or elevated temperatures, metal surfaces can reconstruct (change atomic arrangement), modifying adsorption energetics and kinetics.

---

## 2. Current State (Phases 1-3)

### Implemented Equations

**Langmuir isotherm (single species):**
$$\theta_H = \frac{\sqrt{K_{eq} P_{H_2}}}{1 + \sqrt{K_{eq} P_{H_2}}}$$

**Dissociation flux:**
$$J_{diss} = k_{diss} \times P_{H_2} \times (1-\theta_H)^2$$

**Damköhler number:**
$$Da = \frac{k_{diss} \times K_s \times L}{D}$$

### Key Files

| File | Content |
|------|---------|
| `calculations/surface_kinetics.py` | Core L6 functions |
| `data/surface_kinetics_data.py` | Material kinetics data |
| `calculations/interface_solver.py` | L2+L6 integration |

### Limitations of Current Implementation

1. **Single adsorbate only** — ignores competitive adsorption from H₂O, O₂, CO
2. **Static surface** — assumes fixed surface structure regardless of conditions
3. **No impurity effects** — ignores surface segregation of bulk impurities (S, C, N)
4. **No dynamic coverage** — no transient coverage evolution

---

## 3. Phase 4: Co-Adsorbate Competition

### 3.1 Physical Motivation

In real industrial environments, multiple gas species compete for the same surface sites:

| Species | Effect on H₂ Permeation | Typical Industrial Source |
|---------|------------------------|---------------------------|
| H₂O | Strong site blocking | Steam, humidity |
| O₂ | Oxide formation, strong blocking | Air ingress, steam |
| CO | Moderate blocking | Combustion products |
| CO₂ | Weak blocking (low sticking) | Combustion products |
| N₂ | Negligible (physisorption only) | Air |
| S-compounds | Severe poisoning (ppm levels) | H₂S, SO₂ |

### 3.2 Mathematical Formulation

#### 3.2.1 Multi-Species Langmuir Isotherm

For N competing species, the coverage of species i is:

$$\boxed{\theta_i = \frac{K_i P_i}{1 + \sum_{j=1}^{N} K_j P_j}}$$

For hydrogen (dissociative adsorption, occupies 2 sites):

$$\theta_H = \frac{\sqrt{K_H P_{H_2}}}{1 + \sqrt{K_H P_{H_2}} + \sum_{j \neq H} K_j P_j}$$

#### 3.2.2 Effective Site Availability

The fraction of sites available for H₂ dissociation:

$$f_{available} = (1 - \theta_H - \theta_{block})^2$$

where:
$$\theta_{block} = \sum_{j \neq H} \theta_j$$

#### 3.2.3 Modified Dissociation Flux

$$\boxed{J_{diss} = k_{diss} \times P_{H_2} \times (1 - \theta_H - \theta_{block})^2}$$

#### 3.2.4 Blocking Reduction Factor (BRF)

Define a **Blocking Reduction Factor** to quantify co-adsorbate impact:

$$BRF = \frac{J_{diss,blocked}}{J_{diss,clean}} = \left(\frac{1 - \theta_H - \theta_{block}}{1 - \theta_H}\right)^2$$

| BRF | Interpretation |
|-----|----------------|
| 1.0 | No blocking (clean surface) |
| 0.5 | Moderate blocking |
| < 0.1 | Severe blocking |
| 0 | Complete poisoning |

### 3.3 Implementation Plan

#### 3.3.1 Data Structure Updates

**New file: `data/coadsorbate_data.py`**

```python
COADSORBATES = {
    'H2O': {
        'K_ads_ref': 1e-5,      # Adsorption equilibrium constant at T_ref
        'T_ref': 1073.15,       # Reference temperature [K]
        'E_ads': 40000,         # Adsorption enthalpy [J/mol]
        'sticking_coeff': 0.1,  # Sticking probability
        'sites_occupied': 1,    # Number of sites per molecule
        'desorption_type': 'molecular',  # 'molecular' or 'dissociative'
        'notes': 'Moderate blocking, desorbs at ~400°C'
    },
    'O2': {
        'K_ads_ref': 1e-3,
        'T_ref': 1073.15,
        'E_ads': 80000,         # Strong binding (forms oxide!)
        'sticking_coeff': 0.5,
        'sites_occupied': 2,    # Dissociative adsorption
        'desorption_type': 'dissociative',
        'notes': 'Can oxidize surface, very strong blocking'
    },
    'CO': {
        'K_ads_ref': 1e-6,
        'T_ref': 1073.15,
        'E_ads': 30000,
        'sticking_coeff': 0.05,
        'sites_occupied': 1,
        'desorption_type': 'molecular',
        'notes': 'Moderate blocking, more significant at lower T'
    },
    'H2S': {
        'K_ads_ref': 1e-2,
        'T_ref': 1073.15,
        'E_ads': 100000,        # Extremely strong (poisoning!)
        'sticking_coeff': 0.8,
        'sites_occupied': 1,
        'desorption_type': 'dissociative',
        'notes': 'SEVERE poisoning even at ppm levels!'
    }
}
```

#### 3.3.2 New Functions

**File: `calculations/surface_kinetics.py` (extensions)**

```python
def competitive_langmuir_coverage(P_H2, P_coadsorbates, k_diss, k_recomb, 
                                    coadsorbate_K, temperature):
    """
    Calculate H₂ coverage with competing species.
    
    Parameters
    ----------
    P_H2 : float
        Hydrogen partial pressure [Pa]
    P_coadsorbates : dict
        {species_name: partial_pressure} for blocking species
    k_diss, k_recomb : float
        H₂ dissociation/recombination kinetics
    coadsorbate_K : dict
        {species_name: adsorption_equilibrium_constant}
    temperature : float
        Temperature [K]
    
    Returns
    -------
    dict
        theta_H, theta_block, f_available, BRF
    """
    pass


def calculate_surface_flux_with_blocking(D, K_s, thickness, P_H2, P_coadsorbates,
                                          temperature, k_diss, k_recomb,
                                          coadsorbate_K, coverage_mode='equilibrium'):
    """
    Main L6 function extended for co-adsorbate competition.
    """
    pass
```

#### 3.3.3 Testing Strategy

| Test Case | Expected Result |
|-----------|-----------------|
| Pure H₂ (no blockers) | BRF = 1.0, recovers Phase 3 |
| H₂ + 10% H₂O | BRF < 1, reduced flux |
| H₂ + ppm H₂S | BRF << 1, severe poisoning |
| High T limit | Blockers desorb, BRF → 1 |
| Low T limit | Blockers dominate, BRF → 0 |

### 3.4 References for Phase 4

1. **Campbell, C.T. (1985)**
   "Atomic and molecular oxygen adsorption on Ni(111)"
   *Surf. Sci.* 157, 43-60.
   *→ O₂ adsorption kinetics on Ni*

2. **Ferrin, P. & Mavrikakis, M. (2009)**
   "Structure Sensitivity of Methanol Electrooxidation on Transition Metals"
   *J. Am. Chem. Soc.* 131, 14381-14389.
   *→ CO adsorption data*

3. **Oudar, J. (1980)**
   "Sulfur Adsorption and Poisoning of Metallic Catalysts"
   *Catal. Rev. Sci. Eng.* 22, 171-195.
   *→ Foundational paper on S poisoning*

4. **Mims, C.A. et al. (1993)**
   "Hydrogen transport through metal membranes"
   *J. Phys. Chem.* 97, 2483-2489.
   *→ H₂O effect on Pd membranes*

---

## 4. Phase 5: Dynamic Surface Reconstruction

### 4.1 Physical Motivation

Metal surfaces are not rigid lattices. Under certain conditions, atoms can rearrange:

| Phenomenon | Trigger | Effect on H₂ |
|------------|---------|--------------|
| **Coverage-induced reconstruction** | θ > θ_crit | New adsorption sites |
| **Thermal reconstruction** | T > T_crit | Different surface facets |
| **Adsorbate-induced** | Specific species | Lifting of reconstruction |
| **Subsurface hydrogen** | High θ_H | H penetrates first layer |

### 4.2 Mathematical Formulation

#### 4.2.1 Two-State Surface Model

Assume the surface exists in two states:
- **State A:** Unreconstructed (default at low θ, high T)
- **State B:** Reconstructed (favored at high θ or low T)

The fraction in state B:

$$f_B = \frac{1}{1 + \exp\left(\frac{\Delta G_{A \to B} + n_{H}\Delta\mu_H}{k_B T}\right)}$$

where:
- $\Delta G_{A \to B}$ = Free energy change for reconstruction
- $n_H$ = Number of H atoms involved in reconstruction
- $\Delta\mu_H$ = Chemical potential change of H upon reconstruction

#### 4.2.2 State-Dependent Kinetics

Each state has different kinetics:

| State | k_diss | k_recomb | E_diss | E_recomb |
|-------|--------|----------|--------|----------|
| A (unreconstructed) | $k_{diss}^A$ | $k_{recomb}^A$ | $E_{diss}^A$ | $E_{recomb}^A$ |
| B (reconstructed) | $k_{diss}^B$ | $k_{recomb}^B$ | $E_{diss}^B$ | $E_{recomb}^B$ |

Effective rate constants:

$$\boxed{k_{diss}^{eff} = (1 - f_B) \cdot k_{diss}^A + f_B \cdot k_{diss}^B}$$

$$\boxed{k_{recomb}^{eff} = (1 - f_B) \cdot k_{recomb}^A + f_B \cdot k_{recomb}^B}$$

#### 4.2.3 Coupled Coverage-Reconstruction Equations

The system becomes coupled because:
1. θ depends on kinetics → θ = f(k_diss, k_recomb)
2. Kinetics depend on f_B → k = f(f_B)
3. f_B depends on θ → f_B = f(θ)

This requires **self-consistent iteration**:

```
Initialize: θ⁰ = θ_Langmuir(P, k_diss⁰, k_recomb⁰)
Iterate until convergence:
    1. f_B^n = calculate_reconstruction_fraction(θ^n, T)
    2. k_diss^n, k_recomb^n = get_effective_kinetics(f_B^n)
    3. θ^{n+1} = θ_Langmuir(P, k_diss^n, k_recomb^n)
    4. Check: |θ^{n+1} - θ^n| < tolerance?
```

#### 4.2.4 Surface Reconstruction Factor (SuRF)

Quantify reconstruction impact:

$$SuRF = \frac{J_{with\ reconstruction}}{J_{without\ reconstruction}}$$

### 4.3 Implementation Plan

#### 4.3.1 Data Structure Updates

**New file: `data/surface_reconstruction_data.py`**

```python
RECONSTRUCTION_PARAMS = {
    'Ni_111': {
        'delta_G_AB': -5000,    # Free energy of reconstruction [J/mol]
        'theta_critical': 0.5,  # Coverage threshold for reconstruction
        'T_reconstruction': 600,  # Temperature onset [K]
        
        # State A (unreconstructed) parameters
        'k_diss_A_ref': 2.14e-4,
        'k_recomb_A_ref': 1.23e-12,
        'E_diss_A': 20000,
        'E_recomb_A': 85000,
        
        # State B (reconstructed) parameters
        'k_diss_B_ref': 5.0e-5,   # Often slower (fewer step sites)
        'k_recomb_B_ref': 3.0e-12,
        'E_diss_B': 25000,
        'E_recomb_B': 75000,
        
        'reference': 'Christmann (1988)',
        'notes': 'Ni(111) reconstructs above θ_H ≈ 0.5'
    },
    'Fe_110': {
        'delta_G_AB': -3000,
        'theta_critical': 0.7,
        'T_reconstruction': 500,
        
        'k_diss_A_ref': 1.65e-2,
        'k_recomb_A_ref': 2.17e-12,
        'E_diss_A': 10000,
        'E_recomb_A': 70000,
        
        'k_diss_B_ref': 8.0e-3,
        'k_recomb_B_ref': 4.0e-12,
        'E_diss_B': 12000,
        'E_recomb_B': 65000,
        
        'reference': 'Estrup (1970)',
        'notes': 'BCC Fe less prone to reconstruction'
    }
}
```

#### 4.3.2 New Functions

**File: `calculations/surface_kinetics.py` (extensions)**

```python
def calculate_reconstruction_fraction(theta, temperature, material_name):
    """
    Calculate fraction of surface in reconstructed state.
    
    Uses two-state model with Boltzmann statistics.
    
    Returns
    -------
    dict
        f_B, delta_G, reconstruction_active (bool)
    """
    pass


def get_effective_kinetics_with_reconstruction(f_B, params_A, params_B, temperature):
    """
    Calculate effective k_diss, k_recomb from weighted average.
    """
    pass


def solve_coverage_with_reconstruction(P_H2, temperature, material_name,
                                        max_iter=100, tol=1e-8):
    """
    Self-consistent solver for coupled coverage-reconstruction problem.
    
    Iterates until coverage and reconstruction fraction converge.
    
    Returns
    -------
    dict
        theta, f_B, k_eff_diss, k_eff_recomb, iterations, converged
    """
    pass


def calculate_surface_flux_with_reconstruction(D, K_s, thickness, P_up, P_down,
                                                temperature, material_name):
    """
    Main L6 function extended for dynamic surface reconstruction.
    """
    pass
```

#### 4.3.3 Testing Strategy

| Test Case | Expected Result |
|-----------|-----------------|
| Low θ (θ < θ_crit) | f_B ≈ 0, recovers Phase 3 |
| High θ (θ > θ_crit) | f_B > 0, modified kinetics |
| High T | Reconstruction suppressed |
| Low T + high θ | Maximum reconstruction |
| State A = State B | SuRF = 1 (no effect) |

### 4.4 References for Phase 5

1. **Christmann, K. (1988)**
   "Interaction of hydrogen with solid surfaces"
   *Surf. Sci. Rep.* 9, 1-163.
   DOI: [10.1016/0167-5729(88)90009-X](https://doi.org/10.1016/0167-5729(88)90009-X)
   *→ Comprehensive review of H/metal surface interactions*

2. **Estrup, P.J. & McRae, E.G. (1971)**
   "Surface Studies by Low-Energy Electron Diffraction"
   *Surf. Sci.* 25, 1-52.
   *→ LEED studies of surface reconstruction*

3. **Somorjai, G.A. & Li, Y. (2010)**
   "Introduction to Surface Chemistry and Catalysis"
   Wiley, 2nd Edition.
   *→ Textbook coverage of reconstruction phenomena*

4. **Besenbacher, F. et al. (1990)**
   "STM studies of metal surfaces"
   *Rep. Prog. Phys.* 53, 1253-1295.
   *→ Direct imaging of reconstruction*

---

## 5. Implementation Timeline

### Recommended Sequence

```
Week 1-2: Phase 4 Foundation
├── Create coadsorbate_data.py
├── Implement competitive_langmuir_coverage()
├── Write unit tests for multi-species
└── Validate against single-species limit

Week 3-4: Phase 4 Integration
├── Extend calculate_surface_flux_with_blocking()
├── Add BRF to output dictionary
├── Integrate with interface_solver.py
└── Sensitivity analysis on P_coadsorbates

Week 5-6: Phase 5 Foundation
├── Create surface_reconstruction_data.py
├── Implement calculate_reconstruction_fraction()
├── Implement self-consistent solver
└── Unit tests for two-state model

Week 7-8: Phase 5 Integration
├── Extend main surface flux function
├── Add SuRF to output dictionary
├── Integrate with model hierarchy
└── Combined Phase 4+5 testing

Week 9-10: Validation & Documentation
├── Experimental comparison (if data available)
├── Parameter sensitivity analysis
├── Update documentation
└── Code review and cleanup
```

### Priority Assessment

| Phase | Complexity | Industrial Relevance | Priority |
|-------|------------|---------------------|----------|
| Phase 4 | Medium | **Very High** (real environments have impurities) | **High** |
| Phase 5 | High | Medium (significant at extreme conditions) | Medium |

**Recommendation:** Implement Phase 4 first. It has:
- Lower implementation complexity
- Higher industrial relevance
- More available experimental data for validation

---

## Appendix: Equation Derivations

### A.1 Multi-Species Langmuir Derivation

For species i adsorbing molecularly:
$$\frac{d\theta_i}{dt} = k_{ads,i} P_i (1 - \theta_{total}) - k_{des,i} \theta_i$$

At equilibrium ($d\theta_i/dt = 0$):
$$\theta_i = K_i P_i (1 - \theta_{total})$$

where $K_i = k_{ads,i}/k_{des,i}$

Summing over all species:
$$\theta_{total} = \sum_i K_i P_i (1 - \theta_{total})$$
$$\theta_{total} = (1 - \theta_{total}) \sum_i K_i P_i$$
$$\theta_{total} \left(1 + \sum_i K_i P_i\right) = \sum_i K_i P_i$$
$$1 - \theta_{total} = \frac{1}{1 + \sum_i K_i P_i}$$

Therefore:
$$\boxed{\theta_i = \frac{K_i P_i}{1 + \sum_j K_j P_j}}$$

### A.2 Two-State Reconstruction Model

The free energy difference between states:
$$\Delta G = \Delta G^0 + n_H \Delta\mu_H$$

where:
- $\Delta G^0$ = intrinsic reconstruction free energy
- $\Delta\mu_H = k_B T \ln(\theta/(1-\theta))$ = H chemical potential on surface

The Boltzmann probability of state B:
$$\frac{P_B}{P_A} = \exp\left(-\frac{\Delta G}{k_B T}\right)$$

Since $P_A + P_B = 1$:
$$f_B = \frac{1}{1 + \exp(\Delta G / k_B T)}$$

This is a sigmoidal function of θ (through $\Delta\mu_H$).

### A.3 Self-Consistent Iteration Convergence

The iteration:
$$\theta^{n+1} = F(\theta^n)$$

converges if:
$$\left|\frac{dF}{d\theta}\right|_{\theta^*} < 1$$

where $\theta^*$ is the fixed point.

For most physical systems with reasonable reconstruction parameters, this condition is satisfied because:
1. f_B is a smooth sigmoid
2. Kinetics change gradually
3. θ depends weakly on kinetics perturbations

Under-relaxation can be used if convergence is slow:
$$\theta^{n+1} = \alpha \cdot F(\theta^n) + (1-\alpha) \cdot \theta^n$$

with $0 < \alpha < 1$.

---

*Roadmap document for Hierarchical Hydrogen Permeation Model*  
*Level 6: Surface Kinetics - Phases 4 & 5 Planning*  
*Last updated: 2025*
