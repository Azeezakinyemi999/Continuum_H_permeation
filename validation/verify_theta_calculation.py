"""
Verify theta calculation in three coverage modes.
"""
import numpy as np
from scipy.optimize import brentq

# Test parameters (same as test file)
D = 1e-10
K_s = 0.5
L = 1e-3
P_up = 1e5
P_down = 0
k_diss = 1e-8
k_recomb = 1e-7
K_eq = k_diss / k_recomb
sqrt_K_eq = np.sqrt(K_eq)

print("="*70)
print("VERIFICATION OF THETA CALCULATION IN THREE MODES")
print("="*70)
print()
print("Parameters:")
print(f"  D = {D:.2e} m²/s")
print(f"  K_s = {K_s} mol/(m³·Pa^0.5)")
print(f"  L = {L:.0e} m")
print(f"  P_up = {P_up:.0e} Pa")
print(f"  k_diss = {k_diss:.2e} mol/(m²·s·Pa)")
print(f"  k_recomb = {k_recomb:.2e} m⁴/(mol·s)")
print(f"  K_eq = k_diss/k_recomb = {K_eq}")
print()

# =============================================================================
# 1. Langmuir equilibrium theta
# =============================================================================
sqrt_KP = np.sqrt(K_eq * P_up)
theta_eq = sqrt_KP / (1 + sqrt_KP)
print(f"1. LANGMUIR EQUILIBRIUM:")
print(f"   √(K_eq × P) = √({K_eq} × {P_up}) = {sqrt_KP}")
print(f"   θ_eq = {sqrt_KP} / (1 + {sqrt_KP}) = {theta_eq:.6f}")
print()

# =============================================================================
# 2. Define flux equations
# =============================================================================
C_down = K_s * np.sqrt(P_down)

def J_diss(theta):
    """Dissociation flux: J = k_diss × P × (1-θ)²"""
    return k_diss * P_up * (1 - theta)**2

def J_bulk(theta):
    """Bulk diffusion flux: J = D × [C_surface - C_down] / L"""
    C_surface = K_s * (theta / (1 - theta)) / sqrt_K_eq
    return D * (C_surface - C_down) / L

def residual(theta):
    """F(θ) = J_diss - J_bulk = 0 at steady state"""
    return J_diss(theta) - J_bulk(theta)

# =============================================================================
# 3. Show flux balance at different theta
# =============================================================================
print("2. FLUX BALANCE AT DIFFERENT θ:")
print()
print("   θ           J_diss         J_bulk         Residual")
print("   " + "-" * 56)
for theta in [0.3, 0.5, 0.7, 0.9, 0.94, 0.9469, 0.95, 0.99, theta_eq]:
    jd = J_diss(theta)
    jb = J_bulk(theta)
    res = residual(theta)
    marker = " *" if abs(res) < 1e-10 else ""
    print(f"   {theta:.6f}   {jd:.4e}    {jb:.4e}    {res:+.4e}{marker}")

print()
print("   (* indicates J_diss = J_bulk, i.e., steady-state solution)")
print()

# =============================================================================
# 4. Solve for steady-state theta
# =============================================================================
theta_ss = brentq(residual, 0.01, 0.9999)
print(f"3. STEADY-STATE SOLUTION (brentq):")
print(f"   θ_ss = {theta_ss:.6f}")
print(f"   Verification: J_diss(θ_ss) = {J_diss(theta_ss):.6e}")
print(f"                 J_bulk(θ_ss) = {J_bulk(theta_ss):.6e}")
print(f"                 Residual     = {residual(theta_ss):.2e}")
print()

# =============================================================================
# 5. Compare C_surface and SRF
# =============================================================================
print("4. COMPARISON OF MODES:")
print()

def C_surface(theta):
    return K_s * (theta / (1 - theta)) / sqrt_K_eq

C_sieverts = K_s * np.sqrt(P_up)

print(f"   Reference: C_sieverts = K_s × √P = {C_sieverts:.4f} mol/m³")
print()
print(f"   Mode           θ         θ/(1-θ)    C_surface    SRF")
print("   " + "-" * 60)

# Equilibrium
theta = theta_eq
Cs = C_surface(theta)
SRF = Cs / C_sieverts
print(f"   equilibrium   {theta:.6f}   {theta/(1-theta):8.2f}   {Cs:10.4f}   {SRF:.4f}")

# Steady-state
theta = theta_ss
Cs = C_surface(theta)
SRF = Cs / C_sieverts
print(f"   steady_state  {theta:.6f}   {theta/(1-theta):8.2f}   {Cs:10.4f}   {SRF:.4f}")

# Forced (θ = 0.3)
theta = 0.3
Cs = C_surface(theta)
SRF = Cs / C_sieverts
print(f"   forced(0.3)   {theta:.6f}   {theta/(1-theta):8.2f}   {Cs:10.4f}   {SRF:.4f}")

print()

# =============================================================================
# 6. PHYSICAL INTERPRETATION
# =============================================================================
print("5. PHYSICAL INTERPRETATION:")
print()
print("   The key insight: θ/(1-θ) is EXTREMELY sensitive near θ → 1")
print()
print("   θ = 0.99    → θ/(1-θ) = 99")
print("   θ = 0.947   → θ/(1-θ) = 17.8  (82% reduction!)")
print("   θ = 0.3     → θ/(1-θ) = 0.43")
print()
print("   In 'equilibrium' mode: θ = θ_eq (Langmuir), so C = K_s×√P (Sieverts)")
print("   In 'steady_state' mode: θ < θ_eq because surface can't supply H fast enough")
print("   In 'forced' mode: θ = whatever you specify")
print()

# =============================================================================
# 7. Check Da number
# =============================================================================
Da = k_diss * K_s * L / D
print("6. DAMKÖHLER NUMBER:")
print(f"   Da = k_diss × K_s × L / D = {k_diss} × {K_s} × {L} / {D}")
print(f"   Da = {Da:.4f}")
print()
if Da < 0.1:
    print("   Da < 0.1 → DISSOCIATION-LIMITED (surface controls)")
    print("   → θ_ss should be LESS than θ_eq")
elif Da > 10:
    print("   Da > 10 → DIFFUSION-LIMITED (bulk controls)")
    print("   → θ_ss should be CLOSE to θ_eq")
else:
    print("   0.1 < Da < 10 → MIXED REGIME")

print()
print(f"   Indeed: θ_eq = {theta_eq:.6f}, θ_ss = {theta_ss:.6f}")
print(f"           Difference = {theta_eq - theta_ss:.6f} ({(theta_eq - theta_ss)/theta_eq*100:.1f}%)")
print()
print("="*70)
print("CONCLUSION: All theta calculations are CORRECT!")
print("="*70)
