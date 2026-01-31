# """
# Surface Kinetics Model for Hydrogen Permeation (Level 6)

# This module implements surface chemistry effects on hydrogen permeation,
# incorporating dissociative adsorption and recombinative desorption with
# Langmuir coverage dependence.

# Theory:
# -------
# At metal surfaces, hydrogen transport involves three sequential steps:
# 1. Dissociative adsorption at upstream surface
# 2. Bulk diffusion through the metal
# 3. Recombinative desorption at downstream surface

# When surface kinetics are fast (Damköhler number >> 1), equilibrium is
# established and Sieverts' Law applies (Level 1). When surface kinetics 
# are slow (Da << 1), surface processes become rate-limiting.

# The Pick & Sonnenberg (1985) model treats these as resistances in series:
#     R_total = R_diss + R_bulk + R_recomb

# With Langmuir coverage dependence:
#     J_diss = k_diss × P × (1 - θ_up)²      [dissociation rate]
#     J_recomb = k_recomb × θ_down²           [recombination rate]

# At steady state, all fluxes are equal:
#     J = J_diss = J_bulk = J_recomb

# This requires solving a nonlinear system because:
# - θ depends on the flux (mass balance)
# - Flux depends on θ (rate equations)

# Mathematical Framework:
# -----------------------
# Steady-state mass balance at upstream surface:
#     J = k_diss × P_up × (1 - θ_up)² - k_des × θ_up²
    
# where k_des is the desorption rate (related to k_recomb).

# Langmuir isotherm (equilibrium limit):
#     θ/(1-θ) = √(K_eq × P)
#     K_eq = k_diss/k_des

# Damköhler numbers for regime classification:
#     Da_diss = k_diss × K_s × thickness / D      [upstream]
#     Da_recomb = k_recomb × K_s × thickness / D  [downstream]

# Module Structure:
# -----------------
# Core Functions:
# 1. surface_equilibrium_coverage() - Langmuir isotherm θ(P)
# 2. dissociation_flux() - J = k_diss × P × (1-θ)²
# 3. recombination_flux() - J = k_recomb × θ²
# 4. calculate_damkohler_numbers() - Regime classification
# 5. solve_steady_state_coverage() - Self-consistent solution for J and θ
# 6. calculate_surface_limited_flux() - Main interface function (Level 6)

# References:
# -----------
# 1. Pick, M.A. & Sonnenberg, K. (1985). "A model for atomic hydrogen-metal 
#    interactions." J. Nucl. Mater. 131, 208-220.

# 2. Baskes, M.I. (1980). "A calculation of the surface recombination rate 
#    constant for hydrogen isotopes on metals." J. Nucl. Mater. 92, 318-324.

# 3. Andrew, P.L. & Haasz, A.A. (1992). "Models for hydrogen permeation in 
#    metals." J. Appl. Phys. 72, 2749-2757.

# 4. Causey, R.A. (2002). "Hydrogen isotope retention and recycling in fusion 
#    reactor plasma-facing components." J. Nucl. Mater. 300, 91-117.
# """

# import numpy as np
# from scipy.optimize import least_squares

# # Import existing functions from permeation_calc
# from calculations.permeation_calc import sieverts_concentration, fick_flux

# # Gas constant
# R = 8.314  # J/(mol·K)


# def surface_equilibrium_coverage(pressure, k_diss, k_recomb):
#     """
#     Calculate equilibrium surface coverage using Langmuir isotherm.
    
#     Theory:
#     -------
#     At equilibrium, adsorption and desorption rates balance:
#         k_diss × P × (1-θ)² = k_recomb × θ²
    
#     Rearranging:
#         θ/(1-θ) = √(K_eq × P)
        
#     where K_eq = k_diss/k_recomb is the equilibrium constant.
    
#     Solving for θ:
#         θ = √(K_eq × P) / (1 + √(K_eq × P))
    
#     Parameters
#     ----------
#     pressure : float
#         Hydrogen partial pressure [Pa]
#     k_diss : float
#         Dissociation rate constant [m⁴/(mol·s)]
#     k_recomb : float
#         Recombination rate constant [m²/s]
    
#     Returns
#     -------
#     dict
#         Dictionary containing:
#         - 'theta': Surface coverage fraction (0 ≤ θ ≤ 1) [-]
#         - 'K_eq': Equilibrium constant [m²/(mol·Pa)]
#         - 'sqrt_K_eq_P': Dimensionless parameter √(K_eq×P) [-]
#         - 'regime': 'low_coverage', 'intermediate', or 'high_coverage'
    
#     Raises
#     ------
#     ValueError
#         If pressure < 0
#         If k_diss ≤ 0 or k_recomb ≤ 0
    
#     Example
#     -------
#     >>> result = surface_equilibrium_coverage(
#     ...     pressure=1e5,  # Pa
#     ...     k_diss=1e-2,
#     ...     k_recomb=1e-7
#     ... )
#     >>> print(f"Coverage: {result['theta']:.3f}")
#     """
#     # Input validation
#     if pressure < 0:
#         raise ValueError(f"Pressure must be non-negative, got {pressure} Pa")
#     if k_diss <= 0:
#         raise ValueError(f"k_diss must be positive, got {k_diss}")
#     if k_recomb <= 0:
#         raise ValueError(f"k_recomb must be positive, got {k_recomb}")
    
#     # Handle zero pressure case
#     if pressure == 0:
#         return {
#             'theta': 0.0,
#             'K_eq': k_diss / k_recomb,
#             'sqrt_K_eq_P': 0.0,
#             'regime': 'low_coverage'
#         }
    
#     # Equilibrium constant
#     K_eq = k_diss / k_recomb
    
#     # Dimensionless parameter
#     sqrt_K_eq_P = np.sqrt(K_eq * pressure)
    
#     # Langmuir isotherm solution
#     theta = sqrt_K_eq_P / (1.0 + sqrt_K_eq_P)
    
#     # Classify coverage regime
#     if theta < 0.1:
#         regime = 'low_coverage'
#     elif theta > 0.9:
#         regime = 'high_coverage'
#     else:
#         regime = 'intermediate'
    
#     return {
#         'theta': theta,
#         'K_eq': K_eq,
#         'sqrt_K_eq_P': sqrt_K_eq_P,
#         'regime': regime
#     }


# def dissociation_flux(pressure, theta, k_diss):
#     """
#     Calculate dissociative adsorption flux with Langmuir blocking.
    
#     Theory:
#     -------
#     Dissociation requires two adjacent empty sites:
#         H₂(g) + 2S → 2H(ads)
    
#     The probability of finding two empty sites is (1-θ)².
    
#     Rate equation:
#         J_diss = k_diss × P × (1 - θ)²
    
#     Limiting cases:
#     - θ → 0: J_diss → k_diss × P (maximum rate)
#     - θ → 1: J_diss → 0 (surface blocked)
    
#     Parameters
#     ----------
#     pressure : float
#         Hydrogen partial pressure [Pa]
#     theta : float
#         Surface coverage fraction [-]
#     k_diss : float
#         Dissociation rate constant [m⁴/(mol·s)]
    
#     Returns
#     -------
#     float
#         Dissociation flux [mol/(m²·s)]
    
#     Raises
#     ------
#     ValueError
#         If theta not in [0, 1]
#         If pressure < 0
#     """
#     # Input validation
#     if pressure < 0:
#         raise ValueError(f"Pressure must be non-negative, got {pressure}")
#     if theta < 0 or theta > 1:
#         raise ValueError(f"Coverage must be in [0, 1], got {theta}")
    
#     # Langmuir blocking factor
#     blocking = (1.0 - theta) ** 2
    
#     # Dissociation flux
#     flux = k_diss * pressure * blocking
    
#     return flux


# def recombination_flux(theta, k_recomb, N_surf):
#     """
#     Calculate recombinative desorption flux.
    
#     Theory:
#     -------
#     Recombination requires two adjacent occupied sites:
#         2H(ads) → H₂(g) + 2S
    
#     The probability of finding two adjacent H atoms is proportional to θ².
    
#     Rate equation:
#         J_recomb = k_recomb × (N_surf × θ / N_A)²
    
#     where N_surf is the surface site density [m⁻²] and N_A is Avogadro's number.
    
#     Parameters
#     ----------
#     theta : float
#         Surface coverage fraction [-]
#     k_recomb : float
#         Recombination rate constant [m²/s]
#     N_surf : float
#         Surface site density [m⁻²]
    
#     Returns
#     -------
#     float
#         Recombination flux [mol/(m²·s)]
    
#     Raises
#     ------
#     ValueError
#         If theta not in [0, 1]
#     """
#     # Input validation
#     if theta < 0 or theta > 1:
#         raise ValueError(f"Coverage must be in [0, 1], got {theta}")
    
#     # Convert to surface concentration [mol/m²]
#     N_A = 6.022e23  # Avogadro's number
#     C_surf = N_surf * theta / N_A  # mol/m²
    
#     # Recombination flux (second-order in surface concentration)
#     flux = k_recomb * C_surf ** 2
    
#     return flux


# def calculate_damkohler_numbers(k_diss, k_recomb, D, K_s, thickness, N_surf):
#     """
#     Calculate Damköhler numbers for regime identification.
    
#     Theory:
#     -------
#     The Damköhler number compares surface reaction rate to bulk diffusion rate:
#         Da = (characteristic reaction rate) / (characteristic diffusion rate)
    
#     For hydrogen permeation:
#         Da >> 1: Fast surface kinetics → Sieverts' Law applies (Level 1)
#         Da << 1: Slow surface kinetics → Surface-limited regime (Level 6)
#         Da ~ 1: Mixed control → Full model needed
    
#     Definitions used:
#         Da_diss = k_diss × K_s × thickness / D
#         Da_recomb = k_recomb × C_surf_ref × thickness / D
        
#     where C_surf_ref = N_surf / N_A is a reference surface concentration.
    
#     Parameters
#     ----------
#     k_diss : float
#         Dissociation rate constant [m⁴/(mol·s)]
#     k_recomb : float
#         Recombination rate constant [m²/s]
#     D : float
#         Bulk diffusivity [m²/s]
#     K_s : float
#         Sieverts' constant [mol/(m³·Pa^0.5)]
#     thickness : float
#         Membrane thickness [m]
#     N_surf : float
#         Surface site density [m⁻²]
    
#     Returns
#     -------
#     dict
#         Dictionary containing:
#         - 'Da_diss': Upstream Damköhler number [-]
#         - 'Da_recomb': Downstream Damköhler number [-]
#         - 'regime': 'diffusion-limited', 'dissociation-limited', 
#                    'recombination-limited', 'surface-limited', or 'mixed'
#         - 'description': Human-readable regime description
    
#     Example
#     -------
#     >>> result = calculate_damkohler_numbers(
#     ...     k_diss=1e-2, k_recomb=1e-7,
#     ...     D=1e-10, K_s=0.5, thickness=1e-3, N_surf=1e19
#     ... )
#     >>> print(f"Regime: {result['regime']}")
#     """
#     # Input validation
#     if D <= 0:
#         raise ValueError(f"Diffusivity must be positive, got {D}")
#     if thickness <= 0:
#         raise ValueError(f"Thickness must be positive, got {thickness}")
    
#     # Upstream Damköhler (dissociation vs diffusion)
#     Da_diss = k_diss * K_s * thickness / D
    
#     # Downstream Damköhler (recombination vs diffusion)
#     N_A = 6.022e23
#     C_surf_ref = N_surf / N_A  # Reference surface concentration [mol/m²]
#     Da_recomb = k_recomb * C_surf_ref * thickness / D
    
#     # Regime classification
#     threshold_fast = 10.0   # Da > 10 → equilibrium
#     threshold_slow = 0.1    # Da < 0.1 → surface-limited
    
#     if Da_diss > threshold_fast and Da_recomb > threshold_fast:
#         regime = 'diffusion-limited'
#         description = "Fast surface kinetics; Sieverts' Law applies (Level 1)"
#     elif Da_diss < threshold_slow and Da_recomb > threshold_fast:
#         regime = 'dissociation-limited'
#         description = "Slow upstream dissociation controls flux"
#     elif Da_recomb < threshold_slow and Da_diss > threshold_fast:
#         regime = 'recombination-limited'
#         description = "Slow downstream recombination controls flux"
#     elif Da_diss < threshold_slow and Da_recomb < threshold_slow:
#         regime = 'surface-limited'
#         description = "Both surface processes are slow"
#     else:
#         regime = 'mixed'
#         description = "Comparable surface and bulk resistances"
    
#     return {
#         'Da_diss': Da_diss,
#         'Da_recomb': Da_recomb,
#         'regime': regime,
#         'description': description
#     }


# # def solve_steady_state_coverage(P_up, P_down, D, K_s, thickness, 
# #                                  k_diss, k_recomb, N_surf):
# #     """
# #     Solve for steady-state surface coverages and flux self-consistently.
    
# #     Theory:
# #     -------
# #     At steady state, the flux through each step must be equal:
# #         J = J_diss = J_bulk = J_recomb
    
# #     This gives coupled equations:
    
# #     1. Upstream surface balance:
# #        J = k_diss × P_up × (1 - θ_up)² - k_recomb × (N_surf × θ_up / N_A)²
       
# #     2. Bulk diffusion (using effective pressures from coverage):
# #        J = D × (C_up - C_down) / thickness
# #        where C = K_s × √P_eff and P_eff = [θ/(1-θ)]² / K_eq
       
# #     3. Downstream surface balance:
# #        J = k_recomb × (N_surf × θ_down / N_A)²
# #        (assuming P_down contribution is small)
    
# #     Parameters
# #     ----------
# #     P_up : float
# #         Upstream hydrogen pressure [Pa]
# #     P_down : float
# #         Downstream hydrogen pressure [Pa]
# #     D : float
# #         Bulk diffusivity [m²/s]
# #     K_s : float
# #         Sieverts' constant [mol/(m³·Pa^0.5)]
# #     thickness : float
# #         Membrane thickness [m]
# #     k_diss : float
# #         Dissociation rate constant [m⁴/(mol·s)]
# #     k_recomb : float
# #         Recombination rate constant [m²/s]
# #     N_surf : float
# #         Surface site density [m⁻²]
    
# #     Returns
# #     -------
# #     dict
# #         Dictionary containing:
# #         - 'flux': Steady-state flux [mol/(m²·s)]
# #         - 'theta_up': Upstream surface coverage [-]
# #         - 'theta_down': Downstream surface coverage [-]
# #         - 'C_up': Upstream subsurface concentration [mol/m³]
# #         - 'C_down': Downstream subsurface concentration [mol/m³]
# #         - 'converged': Whether solution converged
# #         - 'iterations': Number of iterations
# #         - 'flux_sieverts': Sieverts' Law flux for comparison [mol/(m²·s)]
# #     """
# #     # Input validation
# #     if P_up < 0 or P_down < 0:
# #         raise ValueError("Pressures must be non-negative")
    
# #     N_A = 6.022e23
    
# #     # Sieverts' Law flux for comparison (Level 1)
# #     C_up_sieverts = sieverts_concentration(K_s, P_up)
# #     C_down_sieverts = sieverts_concentration(K_s, P_down)
# #     flux_sieverts = fick_flux(D, C_up_sieverts, C_down_sieverts, thickness)
    
# #     # Handle zero flux case
# #     if P_up == 0 or (P_up == P_down):
# #         return {
# #             'flux': 0.0,
# #             'theta_up': 0.0,
# #             'theta_down': 0.0,
# #             'C_up': C_up_sieverts,
# #             'C_down': C_down_sieverts,
# #             'converged': True,
# #             'iterations': 0,
# #             'flux_sieverts': flux_sieverts
# #         }
    
# #     # Equilibrium constant
# #     K_eq = k_diss / k_recomb
    
# #     # Initial guess: equilibrium coverages
# #     theta_up_eq = surface_equilibrium_coverage(P_up, k_diss, k_recomb)['theta']
# #     theta_down_eq = surface_equilibrium_coverage(P_down, k_diss, k_recomb)['theta']
    
# #     # Define the system of equations
# #     def equations(x):
# #         theta_up, theta_down = x
        
# #         # Ensure physical bounds
# #         theta_up = np.clip(theta_up, 1e-10, 1.0 - 1e-10)
# #         theta_down = np.clip(theta_down, 1e-10, 1.0 - 1e-10)
        
# #         # Upstream: net flux = dissociation in - recombination out
# #         J_diss_in = dissociation_flux(P_up, theta_up, k_diss)
# #         J_recomb_up = recombination_flux(theta_up, k_recomb, N_surf)
# #         J_net_up = J_diss_in - J_recomb_up
        
# #         # Downstream: net flux = recombination out - dissociation in
# #         J_diss_down = dissociation_flux(P_down, theta_down, k_diss)
# #         J_recomb_down = recombination_flux(theta_down, k_recomb, N_surf)
# #         J_net_down = J_recomb_down - J_diss_down
        
# #         # Bulk diffusion flux using effective pressures from coverage
# #         # From Langmuir: θ/(1-θ) = √(K_eq × P_eff) → P_eff = [θ/(1-θ)]² / K_eq
# #         P_eff_up = (theta_up / (1 - theta_up)) ** 2 / K_eq if theta_up < 1 else P_up
# #         P_eff_down = (theta_down / (1 - theta_down)) ** 2 / K_eq if theta_down < 1 else P_down
        
# #         # Use sieverts_concentration for consistency
# #         C_up = sieverts_concentration(K_s, max(P_eff_up, 0))
# #         C_down = sieverts_concentration(K_s, max(P_eff_down, 0))
# #         J_bulk = fick_flux(D, C_up, C_down, thickness)
        
# #         # Residuals: all fluxes should be equal
# #         eq1 = J_net_up - J_bulk      # Upstream balance
# #         eq2 = J_bulk - J_net_down    # Downstream balance
        
# #         return [eq1, eq2]
    
# #     # Solve the system
# #     x0 = [theta_up_eq, max(theta_down_eq, 0.01)]
    
# #     try:
# #         solution, info, ier, mesg = fsolve(equations, x0, full_output=True)
# #         converged = (ier == 1)
# #         theta_up_sol = np.clip(solution[0], 0, 1)
# #         theta_down_sol = np.clip(solution[1], 0, 1)
# #         iterations = info.get('nfev', 0)
# #     except Exception:
# #         # Fall back to equilibrium values
# #         theta_up_sol = theta_up_eq
# #         theta_down_sol = theta_down_eq
# #         converged = False
# #         iterations = 0
    
# #     # Calculate final flux from downstream recombination
# #     flux = recombination_flux(theta_down_sol, k_recomb, N_surf)
    
# #     # Calculate subsurface concentrations
# #     P_eff_up = (theta_up_sol / (1 - theta_up_sol)) ** 2 / K_eq if theta_up_sol < 1 else P_up
# #     P_eff_down = (theta_down_sol / (1 - theta_down_sol)) ** 2 / K_eq if theta_down_sol < 1 else P_down
# #     C_up_final = sieverts_concentration(K_s, max(P_eff_up, 0))
# #     C_down_final = sieverts_concentration(K_s, max(P_eff_down, 0))
    
# #     return {
# #         'flux': flux,
# #         'theta_up': theta_up_sol,
# #         'theta_down': theta_down_sol,
# #         'C_up': C_up_final,
# #         'C_down': C_down_final,
# #         'converged': converged,
# #         'iterations': iterations,
# #         'flux_sieverts': flux_sieverts
# #     }

# # def solve_steady_state_coverage(P_up, P_down, D, K_s, thickness, 
# #                                  k_diss, k_recomb, N_surf):
# #     """
# #     Solve for steady-state surface coverages and flux self-consistently.
    
# #     Theory:
# #     -------
# #     At steady state, the flux through each step must be equal:
# #         J = J_diss = J_bulk = J_recomb
    
# #     This gives coupled equations:
    
# #     1. Upstream surface balance:
# #        J = k_diss × P_up × (1 - θ_up)² - k_recomb × (N_surf × θ_up / N_A)²
       
# #     2. Bulk diffusion (using effective pressures from coverage):
# #        J = D × (C_up - C_down) / thickness
# #        where C = K_s × √P_eff and P_eff = [θ/(1-θ)]² / K_eq
       
# #     3. Downstream surface balance:
# #        J = k_recomb × (N_surf × θ_down / N_A)²
# #        (assuming P_down contribution is small)
    
# #     Parameters
# #     ----------
# #     P_up : float
# #         Upstream hydrogen pressure [Pa]
# #     P_down : float
# #         Downstream hydrogen pressure [Pa]
# #     D : float
# #         Bulk diffusivity [m²/s]
# #     K_s : float
# #         Sieverts' constant [mol/(m³·Pa^0.5)]
# #     thickness : float
# #         Membrane thickness [m]
# #     k_diss : float
# #         Dissociation rate constant [m⁴/(mol·s)]
# #     k_recomb : float
# #         Recombination rate constant [m²/s]
# #     N_surf : float
# #         Surface site density [m⁻²]
    
# #     Returns
# #     -------
# #     dict
# #         Dictionary containing:
# #         - 'flux': Steady-state flux [mol/(m²·s)]
# #         - 'theta_up': Upstream surface coverage [-]
# #         - 'theta_down': Downstream surface coverage [-]
# #         - 'C_up': Upstream subsurface concentration [mol/m³]
# #         - 'C_down': Downstream subsurface concentration [mol/m³]
# #         - 'converged': Whether solution converged
# #         - 'iterations': Number of iterations
# #         - 'flux_sieverts': Sieverts' Law flux for comparison [mol/(m²·s)]
# #     """
# #     # Input validation
# #     if P_up < 0 or P_down < 0:
# #         raise ValueError("Pressures must be non-negative")
    
# #     N_A = 6.022e23
    
# #     # Sieverts' Law flux for comparison (Level 1)
# #     C_up_sieverts = sieverts_concentration(K_s, P_up)
# #     C_down_sieverts = sieverts_concentration(K_s, P_down)
# #     flux_sieverts = fick_flux(D, C_up_sieverts, C_down_sieverts, thickness)
    
# #     # Handle zero flux case
# #     if P_up == 0 or (P_up == P_down):
# #         return {
# #             'flux': 0.0,
# #             'theta_up': 0.0,
# #             'theta_down': 0.0,
# #             'C_up': C_up_sieverts,
# #             'C_down': C_down_sieverts,
# #             'converged': True,
# #             'iterations': 0,
# #             'flux_sieverts': flux_sieverts
# #         }
    
# #     # Equilibrium constant
# #     K_eq = k_diss / k_recomb
    
# #     # Initial guess: equilibrium coverages
# #     theta_up_eq = surface_equilibrium_coverage(P_up, k_diss, k_recomb)['theta']
    
# #     # For P_down = 0, use a small initial guess for theta_down
# #     if P_down == 0:
# #         theta_down_eq = 0.01  # Small but non-zero
# #     else:
# #         theta_down_eq = surface_equilibrium_coverage(P_down, k_diss, k_recomb)['theta']
    
# #     # Ensure reasonable initial guesses
# #     theta_up_init = np.clip(theta_up_eq, 0.01, 0.99)
# #     theta_down_init = np.clip(theta_down_eq, 0.001, 0.5)  # Downstream typically lower
    
# #     # Define the residual function for least_squares
# #     # No need for manual clipping - least_squares handles bounds natively
# #     def residuals(x):
# #         theta_up, theta_down = x
        
# #         # Upstream: net flux = dissociation in - recombination out
# #         J_diss_in = dissociation_flux(P_up, theta_up, k_diss)
# #         J_recomb_up = recombination_flux(theta_up, k_recomb, N_surf)
# #         J_net_up = J_diss_in - J_recomb_up
        
# #         # Downstream: net flux out = recombination out - dissociation in
# #         J_diss_down = dissociation_flux(P_down, theta_down, k_diss)
# #         J_recomb_down = recombination_flux(theta_down, k_recomb, N_surf)
# #         J_net_down = J_recomb_down - J_diss_down
        
# #         # Bulk diffusion flux using effective pressures from coverage
# #         # From Langmuir: θ/(1-θ) = √(K_eq × P_eff) → P_eff = [θ/(1-θ)]² / K_eq
# #         P_eff_up = (theta_up / (1 - theta_up)) ** 2 / K_eq
# #         P_eff_down = (theta_down / (1 - theta_down)) ** 2 / K_eq
        
# #         # Use sieverts_concentration for consistency
# #         C_up = sieverts_concentration(K_s, max(P_eff_up, 0))
# #         C_down = sieverts_concentration(K_s, max(P_eff_down, 0))
# #         J_bulk = fick_flux(D, C_up, C_down, thickness)
        
# #         # Residuals: all fluxes should be equal
# #         # Normalize residuals by flux_sieverts for better conditioning
# #         scale = max(flux_sieverts, 1e-20)
# #         eq1 = (J_net_up - J_bulk) / scale
# #         eq2 = (J_bulk - J_net_down) / scale
        
# #         return [eq1, eq2]
    
# #     # Use least_squares with Trust Region Reflective (trf) method
# #     # Key advantages over fsolve:
# #     # 1. Native bound constraints (no need for clipping hacks)
# #     # 2. More robust for ill-conditioned problems
# #     # 3. Better convergence near boundaries (θ → 0 or θ → 1)
    
# #     # Physical bounds: 0 < θ < 1
# #     bounds = ([1e-10, 1e-10], [1.0 - 1e-10, 1.0 - 1e-10])
    
# #     # Initial guess
# #     x0 = [theta_up_init, theta_down_init]
    
# #     try:
# #         result = least_squares(
# #             residuals,
# #             x0,
# #             bounds=bounds,
# #             method='trf',  # Trust Region Reflective - robust for bounded problems
# #             ftol=1e-10,
# #             xtol=1e-10,
# #             gtol=1e-10,
# #             max_nfev=1000
# #         )
        
# #         converged = result.success
# #         theta_up_sol = result.x[0]
# #         theta_down_sol = result.x[1]
# #         iterations = result.nfev
        
# #     except Exception:
# #         # Fall back to equilibrium values if solver fails completely
# #         converged = False
# #         theta_up_sol = theta_up_eq
# #         theta_down_sol = max(theta_down_eq, 0.001)
# #         iterations = 0
    
# #     # Calculate final flux from downstream recombination
# #     flux = recombination_flux(theta_down_sol, k_recomb, N_surf)
    
# #     # Calculate subsurface concentrations
# #     if theta_up_sol < 1 - 1e-10:
# #         P_eff_up = (theta_up_sol / (1 - theta_up_sol)) ** 2 / K_eq
# #     else:
# #         P_eff_up = P_up
        
# #     if theta_down_sol < 1 - 1e-10:
# #         P_eff_down = (theta_down_sol / (1 - theta_down_sol)) ** 2 / K_eq
# #     else:
# #         P_eff_down = P_down
        
# #     C_up_final = sieverts_concentration(K_s, max(P_eff_up, 0))
# #     C_down_final = sieverts_concentration(K_s, max(P_eff_down, 0))
    
# #     return {
# #         'flux': flux,
# #         'theta_up': theta_up_sol,
# #         'theta_down': theta_down_sol,
# #         'C_up': C_up_final,
# #         'C_down': C_down_final,
# #         'converged': converged,
# #         'iterations': iterations,
# #         'flux_sieverts': flux_sieverts
# #     }

# def solve_steady_state_coverage(P_up, P_down, D, K_s, thickness, 
#                                  k_diss, k_recomb, N_surf):
#     """
#     Solve for steady-state surface coverages and flux self-consistently.
    
#     Theory:
#     -------
#     At steady state, the flux through each step must be equal:
#         J = J_diss = J_bulk = J_recomb
    
#     Surface mass balances:
#     - Upstream: J = k_diss × P_up × (1-θ_up)² - k_recomb × (C_surf_up)²
#     - Downstream: J = k_recomb × (C_surf_down)² - k_diss × P_down × (1-θ_down)²
    
#     Bulk diffusion connects subsurface concentrations via Fick's law.
#     The subsurface concentration relates to surface coverage through 
#     local equilibrium between surface and immediate subsurface.
    
#     Parameters
#     ----------
#     (same as before)
    
#     Returns
#     -------
#     (same as before)
#     """
#     # Input validation
#     if P_up < 0 or P_down < 0:
#         raise ValueError("Pressures must be non-negative")
    
#     N_A = 6.022e23
    
#     # Sieverts' Law flux for comparison (Level 1)
#     C_up_sieverts = sieverts_concentration(K_s, P_up)
#     C_down_sieverts = sieverts_concentration(K_s, P_down)
#     flux_sieverts = fick_flux(D, C_up_sieverts, C_down_sieverts, thickness)
    
#     # Handle zero flux case
#     if P_up == 0 or (P_up == P_down):
#         return {
#             'flux': 0.0,
#             'theta_up': 0.0,
#             'theta_down': 0.0,
#             'C_up': C_up_sieverts,
#             'C_down': C_down_sieverts,
#             'converged': True,
#             'iterations': 0,
#             'flux_sieverts': flux_sieverts
#         }
    
#     # Equilibrium constant
#     K_eq = k_diss / k_recomb
    
#     # Check Damköhler numbers for regime
#     C_surf_ref = N_surf / N_A
#     Da_diss = k_diss * K_s * thickness / D
#     Da_recomb = k_recomb * C_surf_ref * thickness / D
    
#     # For very fast kinetics (Da >> 1), skip solver and use Sieverts directly
#     if Da_diss > 100 and Da_recomb > 100:
#         # Equilibrium coverages
#         theta_up_eq = surface_equilibrium_coverage(P_up, k_diss, k_recomb)['theta']
#         theta_down_eq = surface_equilibrium_coverage(P_down, k_diss, k_recomb)['theta']
        
#         return {
#             'flux': flux_sieverts,
#             'theta_up': theta_up_eq,
#             'theta_down': theta_down_eq if P_down > 0 else 0.0,
#             'C_up': C_up_sieverts,
#             'C_down': C_down_sieverts,
#             'converged': True,
#             'iterations': 0,
#             'flux_sieverts': flux_sieverts
#         }
    
#     # Initial guess: equilibrium coverages
#     theta_up_eq = surface_equilibrium_coverage(P_up, k_diss, k_recomb)['theta']
    
#     # For P_down = 0, use a small initial guess for theta_down
#     if P_down == 0:
#         theta_down_eq = 0.01
#     else:
#         theta_down_eq = surface_equilibrium_coverage(P_down, k_diss, k_recomb)['theta']
    
#     # Ensure reasonable initial guesses
#     theta_up_init = np.clip(theta_up_eq, 0.01, 0.99)
#     theta_down_init = np.clip(theta_down_eq, 0.001, 0.5)
    
#     def residuals(x):
#         theta_up, theta_down = x
        
#         # Surface concentration [mol/m²]
#         C_surf_up = N_surf * theta_up / N_A
#         C_surf_down = N_surf * theta_down / N_A
        
#         # Upstream surface: net adsorption flux
#         J_diss_in = k_diss * P_up * (1 - theta_up) ** 2
#         J_recomb_up = k_recomb * C_surf_up ** 2
#         J_net_up = J_diss_in - J_recomb_up
        
#         # Downstream surface: net desorption flux
#         J_diss_down = k_diss * P_down * (1 - theta_down) ** 2
#         J_recomb_down = k_recomb * C_surf_down ** 2
#         J_net_down = J_recomb_down - J_diss_down
        
#         # Subsurface concentrations from local equilibrium with surface
#         # At equilibrium: C_subsurface = K_s × √P where P relates to θ via Langmuir
#         # For non-equilibrium, use surface coverage to estimate subsurface C
#         # C_subsurface ≈ (C_surf / C_surf_sat) × C_sieverts_max
#         # Simpler approach: interpolate between 0 and Sieverts value based on θ
        
#         # Use effective pressure from coverage (bounded to avoid overflow)
#         if theta_up < 0.9999:
#             sqrt_K_eq_P_up = theta_up / (1 - theta_up)
#             P_eff_up = min(sqrt_K_eq_P_up ** 2 / K_eq, P_up * 10)  # Cap at 10× P_up
#         else:
#             P_eff_up = P_up
            
#         if theta_down < 0.9999:
#             sqrt_K_eq_P_down = theta_down / (1 - theta_down)
#             P_eff_down = min(sqrt_K_eq_P_down ** 2 / K_eq, P_up * 10)
#         else:
#             P_eff_down = P_down if P_down > 0 else 0
        
#         C_up = sieverts_concentration(K_s, max(P_eff_up, 0))
#         C_down = sieverts_concentration(K_s, max(P_eff_down, 0))
#         J_bulk = fick_flux(D, C_up, C_down, thickness)
        
#         # Residuals normalized by characteristic flux
#         scale = max(flux_sieverts, 1e-20)
#         eq1 = (J_net_up - J_bulk) / scale
#         eq2 = (J_bulk - J_net_down) / scale
        
#         return [eq1, eq2]
    
#     # Bounds
#     bounds = ([1e-10, 1e-10], [1.0 - 1e-10, 1.0 - 1e-10])
#     x0 = [theta_up_init, theta_down_init]
    
#     try:
#         result = least_squares(
#             residuals,
#             x0,
#             bounds=bounds,
#             method='trf',
#             ftol=1e-10,
#             xtol=1e-10,
#             gtol=1e-10,
#             max_nfev=1000
#         )
        
#         converged = result.success
#         theta_up_sol = result.x[0]
#         theta_down_sol = result.x[1]
#         iterations = result.nfev
        
#     except Exception:
#         converged = False
#         theta_up_sol = theta_up_eq
#         theta_down_sol = max(theta_down_eq, 0.001)
#         iterations = 0
    
#     # Calculate final flux from downstream recombination
#     C_surf_down = N_surf * theta_down_sol / N_A
#     flux = k_recomb * C_surf_down ** 2
    
#     # Ensure flux doesn't exceed Sieverts flux (physical constraint)
#     flux = min(flux, flux_sieverts)
    
#     # Calculate subsurface concentrations
#     if theta_up_sol < 1 - 1e-10:
#         sqrt_K_eq_P = theta_up_sol / (1 - theta_up_sol)
#         P_eff_up = min(sqrt_K_eq_P ** 2 / K_eq, P_up * 10)
#     else:
#         P_eff_up = P_up
        
#     if theta_down_sol < 1 - 1e-10:
#         sqrt_K_eq_P = theta_down_sol / (1 - theta_down_sol)
#         P_eff_down = min(sqrt_K_eq_P ** 2 / K_eq, P_up * 10)
#     else:
#         P_eff_down = P_down
        
#     C_up_final = sieverts_concentration(K_s, max(P_eff_up, 0))
#     C_down_final = sieverts_concentration(K_s, max(P_eff_down, 0))
    
#     return {
#         'flux': flux,
#         'theta_up': theta_up_sol,
#         'theta_down': theta_down_sol,
#         'C_up': C_up_final,
#         'C_down': C_down_final,
#         'converged': converged,
#         'iterations': iterations,
#         'flux_sieverts': flux_sieverts
#     }

# def calculate_surface_limited_flux(D, K_s, thickness, P_up, P_down,
#                                     temperature, material_name=None,
#                                     k_diss=None, k_recomb=None, N_surf=None,
#                                     theta_up=None, theta_down=None):
#     """
#     Calculate hydrogen permeation flux including surface kinetics with coverage.
    
#     This is the Level 6 equivalent of calculate_simple_metal_flux().
#     Combines bulk diffusion with Langmuir-type surface kinetics.
    
#     Theory:
#     -------
#     The flux is determined by the slowest step in the series:
#     1. Dissociative adsorption: J_diss = k_diss × P × (1-θ)²
#     2. Bulk diffusion: J_bulk = D × (C_up - C_down) / thickness
#     3. Recombinative desorption: J_recomb = k_recomb × θ²
    
#     At steady state: J_diss = J_bulk = J_recomb = J
    
#     The surface coverage θ can be:
#     - Calculated self-consistently (default)
#     - Provided directly via theta_up/theta_down overrides
    
#     Parameters
#     ----------
#     D : float
#         Diffusion coefficient [m²/s]
#     K_s : float
#         Solubility constant [mol/m³/Pa^0.5]
#     thickness : float
#         Metal thickness [m]
#     P_up : float
#         Upstream hydrogen partial pressure [Pa]
#     P_down : float
#         Downstream hydrogen partial pressure [Pa]
#     temperature : float
#         Temperature [K]
#     material_name : str, optional
#         Material name for property lookup (e.g., 'Incoloy800').
#         Required if k_diss, k_recomb, or N_surf not provided.
#     k_diss : float, optional
#         Override dissociation rate constant [m⁴/(mol·s)]
#     k_recomb : float, optional
#         Override recombination rate constant [m²/s]
#     N_surf : float, optional
#         Override surface site density [m⁻²]
#     theta_up : float, optional
#         Override upstream surface coverage [-]. If provided with theta_down,
#         skips self-consistent solver and uses fixed coverages.
#     theta_down : float, optional
#         Override downstream surface coverage [-]. If provided with theta_up,
#         skips self-consistent solver and uses fixed coverages.
    
#     Returns
#     -------
#     dict
#         Dictionary containing:
#         - 'flux': Permeation flux [mol/(m²·s)]
#         - 'C_up': Upstream subsurface concentration [mol/m³]
#         - 'C_down': Downstream subsurface concentration [mol/m³]
#         - 'permeability': Effective permeability [mol/m/s/Pa^0.5]
#         - 'flux_sieverts': Sieverts' Law flux for comparison [mol/(m²·s)]
#         - 'surface_reduction_factor': flux/flux_sieverts [-]
#         - 'theta_up': Upstream surface coverage [-]
#         - 'theta_down': Downstream surface coverage [-]
#         - 'damkohler': Dict with Da_diss, Da_recomb, regime
#         - 'surface_kinetics': Dict with k_diss, k_recomb, N_surf used
#         - 'converged': Whether solution converged (True if coverages provided)
#         - 'coverage_mode': 'calculated' or 'fixed'
#         - 'units': Dictionary of units
    
#     Raises
#     ------
#     ValueError
#         If material_name not provided and k_diss/k_recomb/N_surf missing
#         If only one of theta_up/theta_down provided (need both or neither)
#         If theta values not in [0, 1]
    
#     Example
#     -------
#     >>> # Self-consistent solution (default)
#     >>> result = calculate_surface_limited_flux(
#     ...     D=1e-10, K_s=0.5, thickness=1e-3,
#     ...     P_up=1e5, P_down=0,
#     ...     temperature=1073,
#     ...     material_name='Incoloy800'
#     ... )
    
#     >>> # Fixed coverage (skip solver)
#     >>> result = calculate_surface_limited_flux(
#     ...     D=1e-10, K_s=0.5, thickness=1e-3,
#     ...     P_up=1e5, P_down=0,
#     ...     temperature=1073,
#     ...     k_diss=1e-2, k_recomb=1e-7, N_surf=1e19,
#     ...     theta_up=0.8, theta_down=0.1
#     ... )
#     """
#     # Import surface kinetics data
#     from data.surface_kinetics_data import SURFACE_KINETICS
    
#     # Input validation
#     if P_up < 0 or P_down < 0:
#         raise ValueError("Pressures must be non-negative")
#     if thickness <= 0:
#         raise ValueError(f"Thickness must be positive: {thickness} m")
#     if D <= 0:
#         raise ValueError(f"Diffusion coefficient must be positive: {D} m²/s")
#     if K_s <= 0:
#         raise ValueError(f"Solubility constant must be positive: {K_s}")
#     if temperature <= 0:
#         raise ValueError(f"Temperature must be positive: {temperature} K")
    
#     # Validate coverage inputs (must provide both or neither)
#     if (theta_up is None) != (theta_down is None):
#         raise ValueError("Must provide both theta_up and theta_down, or neither")
    
#     if theta_up is not None:
#         if not (0 <= theta_up <= 1):
#             raise ValueError(f"theta_up must be in [0, 1], got {theta_up}")
#         if not (0 <= theta_down <= 1):
#             raise ValueError(f"theta_down must be in [0, 1], got {theta_down}")
    
#     # Warning for reverse flow
#     if P_down > P_up:
#         print(f"Warning: Downstream pressure ({P_down} Pa) > Upstream pressure ({P_up} Pa)")
#         print("Flux will be negative (backward flow)")
    
#     # Get surface kinetics parameters
#     if k_diss is None or k_recomb is None or N_surf is None:
#         if material_name is None:
#             raise ValueError("Must provide material_name or all of k_diss, k_recomb, N_surf")
#         if material_name not in SURFACE_KINETICS:
#             raise ValueError(f"Material '{material_name}' not found in SURFACE_KINETICS database")
#         surf = SURFACE_KINETICS[material_name]
        
#         if k_diss is None:
#             k_diss = surf['k_diss_0'] * np.exp(-surf['E_diss'] / (R * temperature))
#         if k_recomb is None:
#             k_recomb = surf['k_recomb_0'] * np.exp(-surf['E_recomb'] / (R * temperature))
#         if N_surf is None:
#             N_surf = surf['N_surf']
    
#     # Calculate Sieverts' Law flux for comparison (Level 1)
#     C_up_sieverts = sieverts_concentration(K_s, P_up)
#     C_down_sieverts = sieverts_concentration(K_s, P_down)
#     flux_sieverts = fick_flux(D, C_up_sieverts, C_down_sieverts, thickness)
    
#     # Determine coverage mode
#     if theta_up is not None and theta_down is not None:
#         # =====================================================================
#         # FIXED COVERAGE MODE: Use provided θ values directly
#         # =====================================================================
#         coverage_mode = 'fixed'
#         converged = True
#         iterations = 0
        
#         # Calculate flux from downstream recombination
#         flux = recombination_flux(theta_down, k_recomb, N_surf)
        
#         # Calculate subsurface concentrations from coverage
#         # From Langmuir: θ/(1-θ) = √(K_eq × P_eff) → P_eff = [θ/(1-θ)]² / K_eq
#         K_eq = k_diss / k_recomb
        
#         if theta_up < 1:
#             P_eff_up = (theta_up / (1 - theta_up)) ** 2 / K_eq
#         else:
#             P_eff_up = P_up
            
#         if theta_down < 1:
#             P_eff_down = (theta_down / (1 - theta_down)) ** 2 / K_eq
#         else:
#             P_eff_down = P_down
        
#         C_up = sieverts_concentration(K_s, max(P_eff_up, 0))
#         C_down = sieverts_concentration(K_s, max(P_eff_down, 0))
        
#     else:
#         # =====================================================================
#         # CALCULATED COVERAGE MODE: Solve self-consistently
#         # =====================================================================
#         coverage_mode = 'calculated'
        
#         solution = solve_steady_state_coverage(
#             P_up=P_up,
#             P_down=P_down,
#             D=D,
#             K_s=K_s,
#             thickness=thickness,
#             k_diss=k_diss,
#             k_recomb=k_recomb,
#             N_surf=N_surf
#         )
        
#         flux = solution['flux']
#         theta_up = solution['theta_up']
#         theta_down = solution['theta_down']
#         C_up = solution['C_up']
#         C_down = solution['C_down']
#         converged = solution['converged']
#         iterations = solution['iterations']
    
#     # Calculate Damköhler numbers
#     damkohler = calculate_damkohler_numbers(
#         k_diss=k_diss,
#         k_recomb=k_recomb,
#         D=D,
#         K_s=K_s,
#         thickness=thickness,
#         N_surf=N_surf
#     )
    
#     # Surface reduction factor
#     surface_reduction_factor = flux / flux_sieverts if flux_sieverts > 0 else 0.0
    
#     # Effective permeability (based on actual flux)
#     delta_sqrt_P = np.sqrt(P_up) - np.sqrt(P_down)
#     if delta_sqrt_P > 0:
#         permeability = flux * thickness / delta_sqrt_P
#     else:
#         permeability = D * K_s  # Fall back to bulk value
    
#     return {
#         # Primary outputs (compatible with calculate_simple_metal_flux)
#         'flux': flux,
#         'C_up': C_up,
#         'C_down': C_down,
#         'permeability': permeability,
        
#         # Level 6 specific outputs
#         'flux_sieverts': flux_sieverts,
#         'surface_reduction_factor': surface_reduction_factor,
#         'theta_up': theta_up,
#         'theta_down': theta_down,
#         'damkohler': damkohler,
        
#         # Surface kinetics parameters used
#         'surface_kinetics': {
#             'k_diss': k_diss,
#             'k_recomb': k_recomb,
#             'N_surf': N_surf,
#             'temperature': temperature
#         },
        
#         # Convergence info
#         'converged': converged,
#         'iterations': iterations,
#         'coverage_mode': coverage_mode,
        
#         # Units
#         'units': {
#             'flux': 'mol/m²/s',
#             'concentration': 'mol/m³',
#             'permeability': 'mol/m/s/Pa^0.5',
#             'coverage': 'dimensionless (0-1)',
#             'k_diss': 'm⁴/(mol·s)',
#             'k_recomb': 'm²/s',
#             'N_surf': 'm⁻²'
#         }
#     }

# # # Self-consistent (default)
# # result = calculate_surface_limited_flux(D=1e-10, K_s=0.5, thickness=1e-3,
# #                                          P_up=1e5, P_down=0, temperature=1073,
# #                                          material_name='Incoloy800')

# # # Fixed coverage with material lookup
# # result = calculate_surface_limited_flux(D=1e-10, K_s=0.5, thickness=1e-3,
# #                                          P_up=1e5, P_down=0, temperature=1073,
# #                                          material_name='Incoloy800',
# #                                          theta_up=0.8, theta_down=0.1)

# # # Fixed coverage with explicit kinetics (no material lookup needed)
# # result = calculate_surface_limited_flux(D=1e-10, K_s=0.5, thickness=1e-3,
# #                                          P_up=1e5, P_down=0, temperature=1073,
# #                                          k_diss=1e-2, k_recomb=1e-7, N_surf=1e19,
# #                                          theta_up=0.8, theta_down=0.1)




"""
Surface Kinetics Model for Hydrogen Permeation (Level 6) - Dissociation Only

Simplified model where only upstream dissociation can be rate-limiting.
Downstream recombination is assumed fast (instantaneous).

Theory:
-------
J = min(J_diss, J_Sieverts)

where:
    J_diss = k_diss × P × (1-θ)²
    θ = √(K_eq × P) / (1 + √(K_eq × P))  [Langmuir isotherm]
    K_eq = k_diss / k_recomb

Damköhler number:
    Da = k_diss × K_s × L / D
    Da >> 1: Sieverts' Law (diffusion-limited)
    Da << 1: Surface-limited (dissociation controls)

References:
-----------
1. Pick, M.A. & Sonnenberg, K. (1985). J. Nucl. Mater. 131, 208-220.
2. Baskes, M.I. (1980). J. Nucl. Mater. 92, 318-324.
"""

import numpy as np

# Import existing functions
from calculations.permeation_calc import sieverts_concentration, fick_flux

# Gas constant
R = 8.314  # J/(mol·K)


def surface_equilibrium_coverage(pressure, k_diss, k_recomb):
    """
    Calculate equilibrium surface coverage using Langmuir isotherm.
    
    θ = √(K_eq × P) / (1 + √(K_eq × P))
    
    Parameters
    ----------
    pressure : float
        Hydrogen partial pressure [Pa]
    k_diss : float
        Dissociation rate constant [m/s] or appropriate units
    k_recomb : float
        Recombination rate constant [m⁴/(mol·s)] or appropriate units
    
    Returns
    -------
    dict with 'theta', 'K_eq', 'sqrt_K_eq_P', 'regime'
    """
    if pressure < 0:
        raise ValueError(f"Pressure must be non-negative, got {pressure}")
    if k_diss <= 0 or k_recomb <= 0:
        raise ValueError("Rate constants must be positive")
    
    if pressure == 0:
        return {'theta': 0.0, 'K_eq': k_diss/k_recomb, 'sqrt_K_eq_P': 0.0, 'regime': 'low_coverage'}
    
    K_eq = k_diss / k_recomb
    sqrt_K_eq_P = np.sqrt(K_eq * pressure)
    theta = sqrt_K_eq_P / (1.0 + sqrt_K_eq_P)
    
    if theta < 0.1:
        regime = 'low_coverage'
    elif theta > 0.9:
        regime = 'high_coverage'
    else:
        regime = 'intermediate'
    
    return {'theta': theta, 'K_eq': K_eq, 'sqrt_K_eq_P': sqrt_K_eq_P, 'regime': regime}


def dissociation_flux(pressure, theta, k_diss):
    """
    Calculate dissociative adsorption flux.
    
    J_diss = k_diss × P × (1-θ)²
    
    Parameters
    ----------
    pressure : float
        Hydrogen partial pressure [Pa]
    theta : float
        Surface coverage fraction (0 ≤ θ ≤ 1)
    k_diss : float
        Dissociation rate constant
    
    Returns
    -------
    float
        Dissociation flux [mol/(m²·s)]
    """
    blocking_factor = (1.0 - theta) ** 2
    return k_diss * pressure * blocking_factor


def calculate_damkohler_number(k_diss, D, K_s, thickness):
    """
    Calculate Damköhler number for dissociation.
    
    Da = k_diss × K_s × L / D
    
    Da >> 1: Diffusion-limited (Sieverts applies)
    Da << 1: Dissociation-limited (surface controls)
    
    Returns
    -------
    dict with 'Da', 'regime', 'description'
    """
    Da = k_diss * K_s * thickness / D
    
    if Da > 10:
        regime = 'diffusion-limited'
        description = "Fast dissociation; Sieverts' Law applies"
    elif Da < 0.1:
        regime = 'dissociation-limited'
        description = "Slow dissociation controls flux"
    else:
        regime = 'mixed'
        description = "Comparable surface and bulk resistances"
    
    return {'Da': Da, 'regime': regime, 'description': description}


def calculate_surface_limited_flux(D, K_s, thickness, P_up, P_down, temperature,
                                    k_diss=None, k_recomb=None,
                                    material_name=None, theta_fixed=None):
    """
    Calculate hydrogen flux with dissociation-limited surface kinetics.
    
    Simple model: J = min(J_diss, J_Sieverts)
    
    Parameters
    ----------
    D : float
        Bulk diffusivity [m²/s]
    K_s : float
        Sieverts' constant [mol/(m³·Pa^0.5)]
    thickness : float
        Membrane thickness [m]
    P_up : float
        Upstream pressure [Pa]
    P_down : float
        Downstream pressure [Pa]
    temperature : float
        Temperature [K]
    k_diss : float, optional
        Dissociation rate constant
    k_recomb : float, optional
        Recombination rate constant (for K_eq calculation)
    material_name : str, optional
        Material name to look up kinetics
    theta_fixed : float, optional
        Fixed coverage (skip Langmuir calculation)
    
    Returns
    -------
    dict with flux, theta, damkohler info, surface_reduction_factor
    """
    # Input validation
    if P_up < 0 or P_down < 0:
        raise ValueError("Pressures must be non-negative")
    
    # Get kinetics parameters
    if k_diss is None or k_recomb is None:
        if material_name is None:
            raise ValueError("Must provide k_diss/k_recomb or material_name")
        from data.surface_kinetics_data import get_surface_kinetics
        kinetics = get_surface_kinetics(material_name, temperature)
        k_diss = kinetics['k_diss']
        k_recomb = kinetics['k_recomb']
    
    # Sieverts flux (Level 1 reference)
    C_up = sieverts_concentration(K_s, P_up)
    C_down = sieverts_concentration(K_s, P_down)
    flux_sieverts = fick_flux(D, C_up, C_down, thickness)
    
    # Damköhler number
    Da_info = calculate_damkohler_number(k_diss, D, K_s, thickness)
    
    # Surface coverage
    if theta_fixed is not None:
        theta = theta_fixed
        coverage_mode = 'fixed'
    else:
        coverage_result = surface_equilibrium_coverage(P_up, k_diss, k_recomb)
        theta = coverage_result['theta']
        coverage_mode = 'calculated'
    
    # Dissociation flux
    J_diss = dissociation_flux(P_up, theta, k_diss)
    
    # Final flux: minimum of dissociation and Sieverts
    flux = min(J_diss, flux_sieverts)
    
    # Surface reduction factor
    SRF = flux / flux_sieverts if flux_sieverts > 0 else 1.0
    
    return {
        'flux': flux,
        'flux_sieverts': flux_sieverts,
        'flux_dissociation': J_diss,
        'theta': theta,
        'coverage_mode': coverage_mode,
        'damkohler': Da_info,
        'surface_reduction_factor': SRF,
        'converged': True,
        'K_eq': k_diss / k_recomb
    }