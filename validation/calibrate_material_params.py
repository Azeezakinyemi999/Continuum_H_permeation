"""
Calibrate material parameters to match JAERI experimental data exactly.
Uses Arrhenius fit to extract P_0 and E_p directly from data.
"""

import numpy as np
from scipy.optimize import curve_fit
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.experimental_data import get_experimental_data, convert_to_SI

R = 8.314  # J/mol/K


def arrhenius_model(T, P_0, E_p):
    """Permeability = P_0 * exp(-E_p / RT)"""
    return P_0 * np.exp(-E_p / (R * T))


def calibrate_to_JAERI():
    """Fit Arrhenius parameters to JAERI experimental data."""
    
    print("=" * 70)
    print("MATERIAL PARAMETER CALIBRATION")
    print("=" * 70)
    
    # Get experimental data
    exp_data_raw = get_experimental_data('Incoloy 800', 'JAERI')
    exp_data = convert_to_SI(exp_data_raw)
    
    T_exp = exp_data['temperatures_K']
    P_exp = exp_data['permeabilities']
    
    print(f"\nExperimental data: {len(T_exp)} points")
    print(f"Temperature range: {T_exp.min():.1f} K to {T_exp.max():.1f} K")
    print(f"                   ({T_exp.min()-273:.0f}°C to {T_exp.max()-273:.0f}°C)")
    print(f"Permeability range: {P_exp.min():.3e} to {P_exp.max():.3e} mol/m/s/Pa^0.5")
    
    # Method 1: Linear fit to ln(P) vs 1/T
    print("\n" + "-" * 70)
    print("METHOD 1: Linear Arrhenius fit")
    print("-" * 70)
    
    ln_P = np.log(P_exp)
    inv_T = 1.0 / T_exp
    
    # ln(P) = ln(P_0) - E_p/R * (1/T)
    slope, intercept = np.polyfit(inv_T, ln_P, 1)
    
    E_p_fit1 = -slope * R
    P_0_fit1 = np.exp(intercept)
    
    print(f"  ln(P) = {intercept:.4f} + ({slope:.1f}) × (1/T)")
    print(f"  P_0 = {P_0_fit1:.4e} mol/m/s/Pa^0.5")
    print(f"  E_p = {E_p_fit1:.0f} J/mol ({E_p_fit1/1000:.1f} kJ/mol)")
    
    # Calculate R² for this fit
    P_model1 = P_0_fit1 * np.exp(-E_p_fit1 / (R * T_exp))
    ss_res = np.sum((np.log(P_exp) - np.log(P_model1))**2)
    ss_tot = np.sum((np.log(P_exp) - np.mean(np.log(P_exp)))**2)
    r2_fit1 = 1 - ss_res / ss_tot
    
    print(f"  R² (log scale) = {r2_fit1:.6f}")
    
    # Method 2: Non-linear curve fit
    print("\n" + "-" * 70)
    print("METHOD 2: Non-linear curve fit")
    print("-" * 70)
    
    try:
        popt, pcov = curve_fit(
            arrhenius_model, T_exp, P_exp,
            p0=[P_0_fit1, E_p_fit1],
            bounds=([1e-10, 10000], [1e0, 200000])
        )
        P_0_fit2, E_p_fit2 = popt
        
        print(f"  P_0 = {P_0_fit2:.4e} mol/m/s/Pa^0.5")
        print(f"  E_p = {E_p_fit2:.0f} J/mol ({E_p_fit2/1000:.1f} kJ/mol)")
        
        P_model2 = arrhenius_model(T_exp, P_0_fit2, E_p_fit2)
        mape2 = np.mean(np.abs(P_model2 - P_exp) / P_exp) * 100
        print(f"  MAPE = {mape2:.1f}%")
    except Exception as e:
        print(f"  Non-linear fit failed: {e}")
        P_0_fit2, E_p_fit2 = P_0_fit1, E_p_fit1
    
    # Final recommended parameters
    print("\n" + "=" * 70)
    print("RECOMMENDED PARAMETERS FOR material_data.py")
    print("=" * 70)
    
    # Use linear fit (more robust)
    P_0_final = P_0_fit1
    E_p_final = E_p_fit1
    
    # Decompose into D and K_s
    # P = D * K_s, E_p = E_D + H_s (if H_s > 0) or E_p = E_D - |H_s| (if H_s < 0)
    # 
    # Keep E_D = 54000 J/mol (standard literature value)
    # Then H_s = E_p - E_D
    
    E_D = 54000  # J/mol (keep from literature)
    H_s = E_p_final - E_D
    
    # P_0 = D_0 * K_s0
    # Keep D_0 reasonable: ~1e-7 to 1e-6 m²/s for metals
    D_0 = 6.4e-7  # m²/s
    K_s0 = P_0_final / D_0
    
    print(f"""
'Incoloy800': {{
    # Calibrated to JAERI-Tech 2002-090 experimental data
    # Arrhenius fit: P = {P_0_final:.4e} × exp(-{E_p_final:.0f} / RT)
    # R² = {r2_fit1:.4f}
    
    'D_0': {D_0:.2e},       # m²/s
    'E_D': {E_D},           # J/mol
    'K_s0': {K_s0:.2e},     # mol/m³/Pa^0.5
    'H_s': {H_s:.0f},       # J/mol
    
    'reference': 'Calibrated to JAERI-Tech 2002-090',
    'temp_range': [600, 1000],
}}
""")
    
    # Verify the fit
    print("-" * 70)
    print("VERIFICATION")
    print("-" * 70)
    
    print(f"\n{'T (°C)':<10} | {'T (K)':<10} | {'P_exp':<14} | {'P_model':<14} | {'Error %':<10}")
    print("-" * 65)
    
    for i in range(0, len(T_exp), 3):  # Every 3rd point
        T_K = T_exp[i]
        P_e = P_exp[i]
        
        D = D_0 * np.exp(-E_D / (R * T_K))
        K_s = K_s0 * np.exp(-H_s / (R * T_K))
        P_m = D * K_s
        
        error = (P_m - P_e) / P_e * 100
        print(f"{T_K-273:<10.0f} | {T_K:<10.1f} | {P_e:<14.4e} | {P_m:<14.4e} | {error:<+10.1f}")
    
    # Calculate final metrics
    P_model_final = []
    for T_K in T_exp:
        D = D_0 * np.exp(-E_D / (R * T_K))
        K_s = K_s0 * np.exp(-H_s / (R * T_K))
        P_model_final.append(D * K_s)
    P_model_final = np.array(P_model_final)
    
    mape_final = np.mean(np.abs(P_model_final - P_exp) / P_exp) * 100
    
    ss_res = np.sum((np.log(P_exp) - np.log(P_model_final))**2)
    ss_tot = np.sum((np.log(P_exp) - np.mean(np.log(P_exp)))**2)
    r2_final = 1 - ss_res / ss_tot
    
    print(f"\nFinal metrics:")
    print(f"  R² (log scale) = {r2_final:.4f}")
    print(f"  MAPE = {mape_final:.1f}%")
    
    return {
        'P_0': P_0_final,
        'E_p': E_p_final,
        'D_0': D_0,
        'E_D': E_D,
        'K_s0': K_s0,
        'H_s': H_s,
        'R2': r2_final,
        'MAPE': mape_final
    }


if __name__ == "__main__":
    results = calibrate_to_JAERI()