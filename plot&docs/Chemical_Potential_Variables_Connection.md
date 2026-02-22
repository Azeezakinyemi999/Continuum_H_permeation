# How Your Variables Connect to Chemical Potential

## 1. Gas Phase (H₂) — The Starting Point

The chemical potential of molecular hydrogen in the gas phase:

**μ_H₂ = μ°_H₂ + RT·ln(P/P₀)**

**Where:**
- μ°_H₂ = standard chemical potential at reference pressure P₀ (usually 1 atm)
- P = actual H₂ pressure
- R = gas constant (8.314 J/mol·K)
- T = absolute temperature (K)

**Physical meaning:** This describes the "thermodynamic availability" of H₂ molecules to participate in reactions. Higher pressure → higher chemical potential → greater driving force for processes that consume H₂.

---

## 2. Surface Coverage (θ) — The Langmuir Connection

### 2.1 Chemical Potential of Adsorbed Hydrogen

For atomic hydrogen adsorbed on the surface, each H atom has chemical potential:

**μ_H,ads = ½·μ_H₂,virtual = ½·[μ°_H₂ + RT·ln(P_virtual)]**

**Key insight:** We use a "virtual pressure" concept because adsorbed atoms aren't actually gas molecules, but we can define an equivalent pressure that would give the same chemical potential.

### 2.2 Statistical Mechanics of Coverage

From statistical thermodynamics, adsorbed H atoms on a surface have chemical potential:

**μ_H,ads = μ°_H,ads + RT·ln[θ/(1−θ)]**

**Where:**
- θ = fractional coverage (occupied sites / total sites)
- θ/(1−θ) = "activity" of adsorbed atoms (accounts for site saturation)
- μ°_H,ads = standard chemical potential of adsorbed state

**Physical meaning:** The θ/(1−θ) term captures:
- As θ → 0: many empty sites, low chemical potential
- As θ → 1: few empty sites, chemical potential rises sharply (site blocking)

### 2.3 Derivation of Langmuir Isotherm from Equilibrium

At **equilibrium**, chemical potentials must balance:

**μ_H,ads = ½·μ_H₂,gas**

Substituting expressions:

**μ°_H,ads + RT·ln[θ/(1−θ)] = ½·[μ°_H₂ + RT·ln(P)]**

Rearranging:

**ln[θ/(1−θ)] = ½·ln(P) + (μ°_H₂/2 − μ°_H,ads)/(RT)**

Define the equilibrium constant:

**K_eq = exp[(μ°_H₂/2 − μ°_H,ads)/(RT)] = K_eq,ref·exp[(−ΔH_ads/R)·(1/T − 1/T_ref)]**

This gives your Langmuir isotherm:

┌─────────────────────────────┐
│  **θ/(1−θ) = √(K_eq·P)**   │
└─────────────────────────────┘

Or equivalently:

**θ = √(K_eq·P) / [1 + √(K_eq·P)]**

### 2.4 The Virtual Pressure Concept

We can invert the Langmuir isotherm to find the **virtual pressure** corresponding to any coverage:

**P_virtual(θ) = (1/K_eq)·[θ/(1−θ)]²**

**Physical interpretation:**
- P_virtual = the gas pressure that would be in equilibrium with coverage θ
- If actual P > P_virtual: net adsorption occurs
- If actual P < P_virtual: net desorption occurs
- The difference drives surface kinetics

---

## 3. Subsurface Concentration — The Sieverts Connection

### 3.1 Your Equation Derived

Your concentration equation:

**C_ox,up = (K_ox/√K_eq)·[θ/(1−θ)]**

**Derivation from chemical potential:**

At the surface-subsurface interface, dissolved hydrogen is in equilibrium with adsorbed hydrogen:

**μ_H,dissolved = μ_H,ads**

For dissolved H in the oxide (Sieverts' law form):

**μ_H,dissolved = μ°_H,dissolved + RT·ln(C/C₀)**

Setting equal to adsorbed H:

**μ°_H,dissolved + RT·ln(C/C₀) = μ°_H,ads + RT·ln[θ/(1−θ)]**

Solving for C:

**C = C₀·exp[(μ°_H,ads − μ°_H,dissolved)/(RT)]·[θ/(1−θ)]**

The prefactor contains the equilibrium constants:

**C_ox,up = (K_ox/√K_eq)·[θ/(1−θ)]**

### 3.2 Equivalence to Effective Pressure Form

Your equation is mathematically equivalent to:

**C_ox,up = K_ox·√P_effective**

**Where the effective pressure is:**

**P_effective = P_virtual(θ) = (1/K_eq)·[θ/(1−θ)]²**

**Proof:**

K_ox·√P_effective = K_ox·√{(1/K_eq)·[θ/(1−θ)]²} = (K_ox/√K_eq)·[θ/(1−θ)]  ✓

### 3.3 Physical Meaning of P_effective

| Condition | P_effective vs P_gas | Physical Interpretation |
|-----------|----------------------|-------------------------|
| Equilibrium | P_eff = P_gas | Surface kinetics not rate-limiting |
| Surface-limited (upstream) | P_eff < P_gas | Adsorption can't keep up; surface "starved" |
| Surface-limited (downstream) | P_eff > P_gas | Desorption can't keep up; H backs up |

---

## 4. Rate of Dissociation — Departure from Equilibrium

### 4.1 Net Dissociation Rate

Your rate expressions:

**r_diss = k_diss·P·(1−θ)²**

**r_recomb = k_recomb·θ²**

The **net rate** (flux through surface):

**J_surface = r_diss − r_recomb = k_diss·P·(1−θ)² − k_recomb·θ²**

### 4.2 Connection to Chemical Potential Driving Force

The net rate can be written as:

**J_surface = k_recomb·(1−θ)²·[K_eq·P − (θ/(1−θ))²]**

Using virtual pressure:

**J_surface = k_recomb·(1−θ)²·K_eq·[P − P_virtual(θ)]**

**Or in chemical potential terms:**

**J_surface = L_surf·[μ_H₂,gas − μ_H₂,virtual(θ)]**

Where L_surf is the surface kinetic coefficient (Onsager coefficient).

**Physical meaning:** The surface flux is driven by the chemical potential difference between:
- The actual gas phase (μ_H₂,gas)
- The "virtual gas" in equilibrium with current coverage (μ_H₂,virtual)

---

## 5. Equilibrium Constant — The Thermodynamic Anchor

### 5.1 Temperature Dependence

**K_eq(T) = K_eq,ref·exp[(−ΔH_ads/R)·(1/T − 1/T_ref)]**

**Where:**
- K_eq,ref = equilibrium constant at T_ref = 1073.15 K
- ΔH_ads = enthalpy of adsorption (negative for exothermic)

### 5.2 Relation to Free Energy

**K_eq = exp(−ΔG°_ads/(RT)) = exp[(−ΔH°_ads + T·ΔS°_ads)/(RT)]**

### 5.3 Microscopic Reversibility (Detailed Balance)

At equilibrium:

**k_diss·P_eq·(1−θ_eq)² = k_recomb·θ_eq²**

This gives:

**K_eq = k_diss/k_recomb = (k_diss,ref/k_recomb,ref)·exp[(−(E_diss − E_recomb)/R)·(1/T − 1/T_ref)]**

**Note:** E_diss − E_recomb = ΔH_ads (thermodynamic consistency!)

---

## 6. Complete Picture — Chemical Potential Profile

```
μ_H₂,gas (upstream)     ← Driving force for entire system →     μ_H₂,gas (downstream)
     |                                                                    |
     ↓ (if surface-limited)                                              ↑
μ_H₂,virtual (θ_up)     ← Gradient through oxide →              μ_H₂,virtual (θ_down)
     |                                                                    |
     ↓ (instantaneous equilibrium)                                       ↑
μ_H,dissolved (C_ox,up) ← Fick's diffusion →                    μ_H,dissolved (C_ox,down)
```

**Your model captures this by:**

| Variable | Chemical Potential Equivalent |
|----------|------------------------------|
| P | μ_H₂,gas = μ° + RT·ln(P/P₀) |
| θ | μ_H,ads = μ°_ads + RT·ln[θ/(1−θ)] |
| C_ox | μ_H,dissolved = μ°_diss + RT·ln(C/C₀) |
| k_diss, k_recomb | Kinetic coefficients in J = L·∇μ |
| K_eq | exp(−ΔG°/(RT)) — sets equilibrium coverage |

---

## 7. Summary: Your Driving Force

For your coupled L2+L6 model:

┌──────────────────────────────────────────────────────────────────────┐
│  **Overall: Δμ_total = (RT/2)·ln(P_up/P_down) = RT·ln(√P_up/√P_down)**  │
└──────────────────────────────────────────────────────────────────────┘

**Decomposed through your variables:**

**Δμ_total = Δμ_surface,up + Δμ_oxide + Δμ_surface,down**
             (if surface-    (∝ ∇C_ox)   (if surface-
              limited)                     limited)

**Practical driving force in your code:** The √P gradient works because it's proportional to the chemical potential gradient for Sieverts' law systems.
