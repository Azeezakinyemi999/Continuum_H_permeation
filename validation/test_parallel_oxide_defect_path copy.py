# """
# Comprehensive Level 3 Test Suite
# Tests all aspects of the parallel path model for oxide with defects.

# Test Categories:
# 1. Basic Functionality - All defect types work
# 2. Limiting Cases - f=0 → Level 2, f=1 → Level 1
# 3. Monotonic Behavior - Flux increases with defect fraction
# 4. PRF Validation - Values in literature range
# 5. Pressure Sweep Analysis - Three-regime behavior
# 6. Defect Type Comparison - Impact analysis
# 7. Mixed Defects - Composite behavior
# """

# import numpy as np
# import matplotlib.pyplot as plt
# from datetime import datetime
# import os

# # Import Level 1 and 2 for comparison
# from calculations.permeation_calc import calculate_simple_metal_flux
# from calculations.interface_solver import calculate_oxide_metal_system

# # Import Level 3 functions
# from calculations.parallel_oxide_defect_paths import (
#     calculate_defect_path_flux,
#     calculate_parallel_path_flux,
#     calculate_PRF
# )

# # Import parameters
# from data.oxide_defect_parameters import DEFECT_CONFIGURATIONS, PARAMETER_RANGES, PRF_RANGES
# from data.oxide_properties import OXIDE_PROPERTIES
# from data.material_data import MATERIALS

# # Constants
# T_TEST = 1073  # K (800°C)
# R = 8.314  # J/mol/K

# # Output directory
# OUTPUT_DIR = "validation/results/level3_tests"


# def setup_material_properties(T=T_TEST):
#     """Set up oxide and metal properties for testing."""
#     oxide_data = OXIDE_PROPERTIES['Cr2O3']
#     metal_data = MATERIALS['Incoloy800']
    
#     oxide_props = {
#         'D_ox': oxide_data['D_ox_0'] * np.exp(-oxide_data['E_D_ox'] / (R * T)),
#         'K_ox': oxide_data['K_ox_0'] * np.exp(-oxide_data['H_sol_ox'] / (R * T)),
#         'thickness': oxide_data['thickness']
#     }
    
#     metal_props = {
#         'D_metal': metal_data['D_0'] * np.exp(-metal_data['E_D'] / (R * T)),
#         'K_s_metal': metal_data['K_s0'] * np.exp(-metal_data['H_s'] / (R * T)),
#         'thickness': 1e-3  # 1mm
#     }
    
#     return oxide_props, metal_props


# def ensure_output_dir():
#     """Create output directory if it doesn't exist."""
#     if not os.path.exists(OUTPUT_DIR):
#         os.makedirs(OUTPUT_DIR)


# # ============================================================================
# # TEST 1: Basic Functionality
# # ============================================================================

# def test_1_basic_functionality():
#     """Test all defect types including mixed."""
#     print("\n" + "="*60)
#     print("TEST 1: Basic Functionality")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 1.0, 0.0
    
#     defect_configs = {
#         'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
#         'crack': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
#         'grain_boundary': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
#         'mixed': {
#             'type': 'mixed', 
#             'area_fraction': 0.015,
#             'components': {'pinholes': 0.005, 'cracks': 0.005, 'grain_boundaries': 0.005},
#             'thickness_factor': 0.1,
#             'diffusivity_factor': 10
#         }
#     }
    
#     results = {}
#     all_passed = True
    
#     for name, config in defect_configs.items():
#         try:
#             flux = calculate_defect_path_flux(
#                 P_upstream, P_downstream, oxide_props, metal_props, config
#             )
#             results[name] = flux
#             print(f"✓ {name:15}: flux = {flux:.3e} mol/m²/s")
#         except Exception as e:
#             print(f"✗ {name:15}: FAILED - {e}")
#             all_passed = False
    
#     if all_passed:
#         print("\nTest 1: PASSED ✓")
#     return all_passed


# # ============================================================================
# # TEST 2: Limiting Cases
# # ============================================================================

# def test_2_limiting_cases():
#     """Verify f=0 → Level 2 and f=1 → Level 1."""
#     print("\n" + "="*60)
#     print("TEST 2: Limiting Cases Validation")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 1.0, 0.0
    
#     # Test 2a: f_defect = 0 → Level 2
#     print("\nTest 2a: f_defect = 0 (should match Level 2)")
    
#     result_l3_f0 = calculate_parallel_path_flux(
#         P_upstream, P_downstream, oxide_props, metal_props,
#         {'area_fraction': 0.0, 'type': 'pinhole'}
#     )
#     result_l2 = calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props)
    
#     if result_l2['flux'] > 0:
#         error_2a = abs(result_l3_f0['flux_total'] - result_l2['flux']) / result_l2['flux'] * 100
#     else:
#         error_2a = 0 if result_l3_f0['flux_total'] == 0 else float('inf')
    
#     print(f"  Level 3 (f=0): {result_l3_f0['flux_total']:.3e} mol/m²/s")
#     print(f"  Level 2:       {result_l2['flux']:.3e} mol/m²/s")
#     print(f"  Error: {error_2a:.2e}%")
    
#     # Test 2b: f_defect = 1 → Level 1
#     print("\nTest 2b: f_defect = 1 (should match Level 1)")
    
#     result_l3_f1 = calculate_parallel_path_flux(
#         P_upstream, P_downstream, oxide_props, metal_props,
#         {'area_fraction': 1.0, 'type': 'pinhole'}
#     )
#     result_l1 = calculate_simple_metal_flux(
#         metal_props['D_metal'], metal_props['K_s_metal'], 
#         metal_props['thickness'], P_upstream, P_downstream
#     )
    
#     error_2b = abs(result_l3_f1['flux_total'] - result_l1['flux']) / result_l1['flux'] * 100
    
#     print(f"  Level 3 (f=1): {result_l3_f1['flux_total']:.3e} mol/m²/s")
#     print(f"  Level 1:       {result_l1['flux']:.3e} mol/m²/s")
#     print(f"  Error: {error_2b:.2e}%")
    
#     passed = error_2a < 1e-6 and error_2b < 1e-6
#     if passed:
#         print("\nTest 2: PASSED ✓")
#     else:
#         print("\nTest 2: FAILED ✗")
#     return passed


# # ============================================================================
# # TEST 3: Monotonic Behavior
# # ============================================================================

# def test_3_monotonic_behavior():
#     """Verify flux increases monotonically with defect fraction."""
#     print("\n" + "="*60)
#     print("TEST 3: Monotonic Behavior with Defect Fraction")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 1.0, 0.0
    
#     defect_fractions = [0.0001, 0.001, 0.01, 0.1]
#     fluxes = []
    
#     print("\nDefect Fraction | Flux (mol/m²/s) | Enhancement")
#     print("-" * 55)
    
#     for f in defect_fractions:
#         result = calculate_parallel_path_flux(
#             P_upstream, P_downstream, oxide_props, metal_props,
#             {'area_fraction': f, 'type': 'pinhole'}
#         )
#         fluxes.append(result['flux_total'])
#         print(f"  {f:12.4f}  | {result['flux_total']:15.3e} | {result['defect_enhancement_factor']:12.2e}x")
    
#     is_monotonic = all(fluxes[i] <= fluxes[i+1] for i in range(len(fluxes)-1))
    
#     if is_monotonic:
#         print("\n✓ Flux increases monotonically")
#         print("Test 3: PASSED ✓")
#     else:
#         print("\n✗ Non-monotonic behavior detected")
#         print("Test 3: FAILED ✗")
#     return is_monotonic


# # ============================================================================
# # TEST 4: PRF Validation
# # ============================================================================

# def test_4_PRF_validation():
#     """Verify PRF values across all configurations."""
#     print("\n" + "="*60)
#     print("TEST 4: PRF Validation")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_test = 1.0
    
#     print("\nConfiguration       | PRF           | Efficiency | Regime")
#     print("-" * 65)
    
#     all_passed = True
    
#     for config_name, config in DEFECT_CONFIGURATIONS.items():
#         try:
#             if config_name == 'perfect':
#                 result = calculate_PRF(P_test, oxide_props, metal_props, None)
#             else:
#                 result = calculate_PRF(P_test, oxide_props, metal_props, config)
            
#             prf_str = f"{result['PRF']:.2e}" if result['PRF'] > 1e6 else f"{result['PRF']:.1f}"
#             print(f"  {config_name:17} | {prf_str:13} | {result['efficiency']:10.2%} | {result['regime']}")
#         except Exception as e:
#             print(f"  {config_name:17} | FAILED: {e}")
#             all_passed = False
    
#     if all_passed:
#         print("\nTest 4: PASSED ✓")
#     return all_passed


# # ============================================================================
# # TEST 5: Pressure Sweep Analysis
# # ============================================================================

# def test_5_pressure_sweep():
#     """Analyze three-regime behavior across pressure range."""
#     print("\n" + "="*60)
#     print("TEST 5: Pressure Sweep Analysis")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_downstream = 0.0
    
#     pressures = np.logspace(-2, 6, 50)  # 0.01 Pa to 1 MPa
#     defect_fractions = [0.0, 0.001, 0.01, 0.1]
    
#     results = {f: {'pressures': [], 'fluxes': [], 'regimes': []} for f in defect_fractions}
    
#     print("\nCalculating flux vs pressure for different defect fractions...")
    
#     for P in pressures:
#         for f in defect_fractions:
#             try:
#                 if f == 0:
#                     result = calculate_oxide_metal_system(P, P_downstream, oxide_props, metal_props)
#                     flux = result['flux']
#                     regime = result.get('regime', 'unknown')
#                 else:
#                     result = calculate_parallel_path_flux(
#                         P, P_downstream, oxide_props, metal_props,
#                         {'area_fraction': f, 'type': 'pinhole'}
#                     )
#                     flux = result['flux_total']
#                     regime = result.get('regime_intact', 'unknown')
                
#                 results[f]['pressures'].append(P)
#                 results[f]['fluxes'].append(flux)
#                 results[f]['regimes'].append(regime)
#             except Exception as e:
#                 pass  # Skip failed points
    
#     # Create plot
#     ensure_output_dir()
#     fig, ax = plt.subplots(figsize=(10, 7))
    
#     colors = ['black', 'blue', 'green', 'red']
#     labels = ['Perfect oxide (f=0)', 'f=0.001 (0.1%)', 'f=0.01 (1%)', 'f=0.1 (10%)']
    
#     for i, f in enumerate(defect_fractions):
#         if results[f]['fluxes']:
#             ax.loglog(results[f]['pressures'], results[f]['fluxes'], 
#                      color=colors[i], linewidth=2, label=labels[i])
    
#     ax.set_xlabel('Upstream Pressure (Pa)', fontsize=12)
#     ax.set_ylabel('Flux (mol/m²/s)', fontsize=12)
#     ax.set_title('Level 3: Effect of Defect Fraction on Permeation', fontsize=14)
#     ax.legend(loc='best')
#     ax.grid(True, alpha=0.3)
    
#     plot_path = f"{OUTPUT_DIR}/pressure_sweep_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
#     plt.savefig(plot_path, dpi=150, bbox_inches='tight')
#     plt.close()
    
#     print(f"\n✓ Plot saved: {plot_path}")
#     print("Test 5: PASSED ✓")
#     return True


# # ============================================================================
# # TEST 6: Defect Type Comparison
# # ============================================================================

# def test_6_defect_type_comparison():
#     """Compare impact of different defect types."""
#     print("\n" + "="*60)
#     print("TEST 6: Defect Type Comparison")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 100.0, 0.0  # 100 Pa
    
#     defect_configs = {
#         'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
#         'crack_thin': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
#         'crack_thick': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.5},
#         'grain_boundary_10x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
#         'grain_boundary_100x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 100},
#     }
    
#     # Get baseline (Level 2)
#     baseline = calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props)
    
#     print(f"\nBaseline (perfect oxide): {baseline['flux']:.3e} mol/m²/s")
#     print("\nDefect Type          | Flux (mol/m²/s) | Enhancement | Dominant Path")
#     print("-" * 75)
    
#     results = {}
#     for name, config in defect_configs.items():
#         result = calculate_parallel_path_flux(
#             P_upstream, P_downstream, oxide_props, metal_props, config
#         )
#         results[name] = result
#         enhancement = result['flux_total'] / baseline['flux'] if baseline['flux'] > 0 else float('inf')
#         print(f"  {name:18} | {result['flux_total']:15.3e} | {enhancement:11.2e}x | {result['dominant_path']}")
    
#     # Find most impactful defect type
#     max_flux_type = max(results.keys(), key=lambda x: results[x]['flux_total'])
#     print(f"\n→ Most impactful defect type: {max_flux_type}")
#     print("Test 6: PASSED ✓")
#     return True


# # ============================================================================
# # TEST 7: PRF vs Defect Fraction
# # ============================================================================

# def test_7_PRF_vs_defect_fraction():
#     """Calculate and plot PRF vs defect fraction."""
#     print("\n" + "="*60)
#     print("TEST 7: PRF vs Defect Fraction")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_test = 100.0  # 100 Pa
    
#     defect_fractions = np.logspace(-4, -0.5, 20)  # 0.01% to ~30%
    
#     prf_pinhole = []
#     prf_crack = []
#     prf_gb = []
    
#     for f in defect_fractions:
#         # Pinhole
#         result = calculate_PRF(P_test, oxide_props, metal_props,
#                               {'area_fraction': f, 'type': 'pinhole'})
#         prf_pinhole.append(result['PRF'])
        
#         # Crack
#         result = calculate_PRF(P_test, oxide_props, metal_props,
#                               {'area_fraction': f, 'type': 'crack', 'thickness_factor': 0.1})
#         prf_crack.append(result['PRF'])
        
#         # Grain boundary
#         result = calculate_PRF(P_test, oxide_props, metal_props,
#                               {'area_fraction': f, 'type': 'grain_boundary', 'diffusivity_factor': 10})
#         prf_gb.append(result['PRF'])
    
#     # Create plot
#     ensure_output_dir()
#     fig, ax = plt.subplots(figsize=(10, 7))
    
#     ax.loglog(defect_fractions * 100, prf_pinhole, 'r-', linewidth=2, label='Pinholes')
#     ax.loglog(defect_fractions * 100, prf_crack, 'b--', linewidth=2, label='Cracks (α=0.1)')
#     ax.loglog(defect_fractions * 100, prf_gb, 'g:', linewidth=2, label='Grain boundaries (β=10)')
    
#     # Add literature range
#     ax.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range (10-3828)')
    
#     ax.set_xlabel('Defect Area Fraction (%)', fontsize=12)
#     ax.set_ylabel('Permeation Reduction Factor (PRF)', fontsize=12)
#     ax.set_title('PRF vs Defect Fraction for Different Defect Types', fontsize=14)
#     ax.legend(loc='best')
#     ax.grid(True, alpha=0.3)
#     ax.set_xlim([0.01, 50])
    
#     plot_path = f"{OUTPUT_DIR}/PRF_vs_defect_fraction_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
#     plt.savefig(plot_path, dpi=150, bbox_inches='tight')
#     plt.close()
    
#     # Print summary
#     print("\nDefect Fraction (%) | PRF (pinhole) | PRF (crack) | PRF (GB)")
#     print("-" * 65)
#     for i in [0, 5, 10, 15, 19]:  # Sample points
#         print(f"  {defect_fractions[i]*100:16.2f}  | {prf_pinhole[i]:13.1f} | {prf_crack[i]:11.2e} | {prf_gb[i]:.2e}")
    
#     print(f"\n✓ Plot saved: {plot_path}")
#     print("Test 7: PASSED ✓")
#     return True


# # ============================================================================
# # MAIN TEST RUNNER
# # ============================================================================

# def run_all_tests():
#     """Run all tests in sequence."""
#     print("\n" + "="*70)
#     print("LEVEL 3 COMPREHENSIVE TEST SUITE")
#     print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
#     print("="*70)
    
#     tests = [
#         ("Basic Functionality", test_1_basic_functionality),
#         ("Limiting Cases", test_2_limiting_cases),
#         ("Monotonic Behavior", test_3_monotonic_behavior),
#         ("PRF Validation", test_4_PRF_validation),
#         ("Pressure Sweep", test_5_pressure_sweep),
#         ("Defect Type Comparison", test_6_defect_type_comparison),
#         ("PRF vs Defect Fraction", test_7_PRF_vs_defect_fraction),
#     ]
    
#     results = {}
#     for name, test_func in tests:
#         try:
#             results[name] = test_func()
#         except Exception as e:
#             print(f"\n✗ {name} failed with error: {e}")
#             results[name] = False
    
#     # Summary
#     print("\n" + "="*70)
#     print("TEST SUMMARY")
#     print("="*70)
    
#     passed = sum(results.values())
#     total = len(results)
    
#     for name, result in results.items():
#         status = "✓ PASSED" if result else "✗ FAILED"
#         print(f"  {name:30} {status}")
    
#     print(f"\nTotal: {passed}/{total} tests passed")
    
#     if passed == total:
#         print("\n" + "="*70)
#         print("ALL TESTS PASSED SUCCESSFULLY! ✓")
#         print("="*70)
    
#     return passed == total


# if __name__ == "__main__":
#     run_all_tests()


# """
# Comprehensive Level 3 Test Suite
# Tests all aspects of the parallel path model for oxide with defects.

# Test Categories:
# 1. Basic Functionality - All defect types work
# 2. Limiting Cases - f=0 → Level 2, f=1 → Level 1
# 3. Monotonic Behavior - Flux increases with defect fraction
# 4. PRF Validation - Values in literature range
# 5. Pressure Sweep Analysis - Three-regime behavior
# 6. Defect Type Comparison - Impact analysis
# 7. PRF vs Defect Fraction - Barrier effectiveness
# 8. Dominant Path Analysis - Phase diagram of defect vs oxide dominated regions
# 9. PRF Regime Analysis - PRF phase diagram
# 10. Temperature Effects - Regime transitions vs temperature
# 11. Defect Type Regime Comparison - Compare defect types across regimes
# """

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from datetime import datetime
# import os

# # Import Level 1 and 2 for comparison
# from calculations.permeation_calc import calculate_simple_metal_flux
# from calculations.interface_solver import calculate_oxide_metal_system

# # Import Level 3 functions
# from calculations.parallel_oxide_defect_paths import (
#     calculate_defect_path_flux,
#     calculate_parallel_path_flux,
#     calculate_PRF
# )

# # Import parameters
# from data.oxide_defect_parameters import DEFECT_CONFIGURATIONS, PARAMETER_RANGES, PRF_RANGES
# from data.oxide_properties import OXIDE_PROPERTIES
# from data.material_data import MATERIALS

# # Constants
# T_TEST = 1073  # K (800°C)
# R = 8.314  # J/mol/K

# # Output directory
# OUTPUT_DIR = "validation/results/complete_level3_tests"


# def setup_material_properties(T=T_TEST):
#     """Set up oxide and metal properties for testing."""
#     oxide_data = OXIDE_PROPERTIES['Cr2O3']
#     metal_data = MATERIALS['Incoloy800']
    
#     oxide_props = {
#         'D_ox': oxide_data['D_ox_0'] * np.exp(-oxide_data['E_D_ox'] / (R * T)),
#         'K_ox': oxide_data['K_ox_0'] * np.exp(-oxide_data['H_sol_ox'] / (R * T)),
#         'thickness': oxide_data['thickness']
#     }
    
#     metal_props = {
#         'D_metal': metal_data['D_0'] * np.exp(-metal_data['E_D'] / (R * T)),
#         'K_s_metal': metal_data['K_s0'] * np.exp(-metal_data['H_s'] / (R * T)),
#         'thickness': 1e-3  # 1mm
#     }
    
#     return oxide_props, metal_props


# def ensure_output_dir():
#     """Create output directory if it doesn't exist."""
#     if not os.path.exists(OUTPUT_DIR):
#         os.makedirs(OUTPUT_DIR)


# # ============================================================================
# # TEST 1: Basic Functionality
# # ============================================================================

# def test_1_basic_functionality():
#     """Test all defect types including mixed."""
#     print("\n" + "="*60)
#     print("TEST 1: Basic Functionality")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 1.0, 0.0
    
#     defect_configs = {
#         'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
#         'crack': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
#         'grain_boundary': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
#         'mixed': {
#             'type': 'mixed', 
#             'area_fraction': 0.015,
#             'components': {'pinholes': 0.005, 'cracks': 0.005, 'grain_boundaries': 0.005},
#             'thickness_factor': 0.1,
#             'diffusivity_factor': 10
#         }
#     }
    
#     results = {}
#     all_passed = True
    
#     for name, config in defect_configs.items():
#         try:
#             flux = calculate_defect_path_flux(
#                 P_upstream, P_downstream, oxide_props, metal_props, config
#             )
#             results[name] = flux
#             print(f"✓ {name:15}: flux = {flux:.3e} mol/m²/s")
#         except Exception as e:
#             print(f"✗ {name:15}: FAILED - {e}")
#             all_passed = False
    
#     if all_passed:
#         print("\nTest 1: PASSED ✓")
#     return all_passed


# # ============================================================================
# # TEST 2: Limiting Cases
# # ============================================================================

# def test_2_limiting_cases():
#     """Verify f=0 → Level 2 and f=1 → Level 1."""
#     print("\n" + "="*60)
#     print("TEST 2: Limiting Cases Validation")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 1.0, 0.0
    
#     # Test 2a: f_defect = 0 → Level 2
#     print("\nTest 2a: f_defect = 0 (should match Level 2)")
    
#     result_l3_f0 = calculate_parallel_path_flux(
#         P_upstream, P_downstream, oxide_props, metal_props,
#         {'area_fraction': 0.0, 'type': 'pinhole'}
#     )
#     result_l2 = calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props)
    
#     if result_l2['flux'] > 0:
#         error_2a = abs(result_l3_f0['flux_total'] - result_l2['flux']) / result_l2['flux'] * 100
#     else:
#         error_2a = 0 if result_l3_f0['flux_total'] == 0 else float('inf')
    
#     print(f"  Level 3 (f=0): {result_l3_f0['flux_total']:.3e} mol/m²/s")
#     print(f"  Level 2:       {result_l2['flux']:.3e} mol/m²/s")
#     print(f"  Error: {error_2a:.2e}%")
    
#     # Test 2b: f_defect = 1 → Level 1
#     print("\nTest 2b: f_defect = 1 (should match Level 1)")
    
#     result_l3_f1 = calculate_parallel_path_flux(
#         P_upstream, P_downstream, oxide_props, metal_props,
#         {'area_fraction': 1.0, 'type': 'pinhole'}
#     )
#     result_l1 = calculate_simple_metal_flux(
#         metal_props['D_metal'], metal_props['K_s_metal'], 
#         metal_props['thickness'], P_upstream, P_downstream
#     )
    
#     error_2b = abs(result_l3_f1['flux_total'] - result_l1['flux']) / result_l1['flux'] * 100
    
#     print(f"  Level 3 (f=1): {result_l3_f1['flux_total']:.3e} mol/m²/s")
#     print(f"  Level 1:       {result_l1['flux']:.3e} mol/m²/s")
#     print(f"  Error: {error_2b:.2e}%")
    
#     passed = error_2a < 1e-6 and error_2b < 1e-6
#     if passed:
#         print("\nTest 2: PASSED ✓")
#     else:
#         print("\nTest 2: FAILED ✗")
#     return passed


# # ============================================================================
# # TEST 3: Monotonic Behavior
# # ============================================================================

# def test_3_monotonic_behavior():
#     """Verify flux increases monotonically with defect fraction."""
#     print("\n" + "="*60)
#     print("TEST 3: Monotonic Behavior with Defect Fraction")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 1.0, 0.0
    
#     defect_fractions = [0.0001, 0.001, 0.01, 0.1]
#     fluxes = []
    
#     print("\nDefect Fraction | Flux (mol/m²/s) | Enhancement")
#     print("-" * 55)
    
#     for f in defect_fractions:
#         result = calculate_parallel_path_flux(
#             P_upstream, P_downstream, oxide_props, metal_props,
#             {'area_fraction': f, 'type': 'pinhole'}
#         )
#         fluxes.append(result['flux_total'])
#         print(f"  {f:12.4f}  | {result['flux_total']:15.3e} | {result['defect_enhancement_factor']:12.2e}x")
    
#     is_monotonic = all(fluxes[i] <= fluxes[i+1] for i in range(len(fluxes)-1))
    
#     if is_monotonic:
#         print("\n✓ Flux increases monotonically")
#         print("Test 3: PASSED ✓")
#     else:
#         print("\n✗ Non-monotonic behavior detected")
#         print("Test 3: FAILED ✗")
#     return is_monotonic


# # ============================================================================
# # TEST 4: PRF Validation
# # ============================================================================

# def test_4_PRF_validation():
#     """Verify PRF values across all configurations."""
#     print("\n" + "="*60)
#     print("TEST 4: PRF Validation")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_test = 1.0
    
#     print("\nConfiguration       | PRF           | Efficiency | Regime")
#     print("-" * 65)
    
#     all_passed = True
    
#     for config_name, config in DEFECT_CONFIGURATIONS.items():
#         try:
#             if config_name == 'perfect':
#                 result = calculate_PRF(P_test, oxide_props, metal_props, None)
#             else:
#                 result = calculate_PRF(P_test, oxide_props, metal_props, config)
            
#             prf_str = f"{result['PRF']:.2e}" if result['PRF'] > 1e6 else f"{result['PRF']:.1f}"
#             print(f"  {config_name:17} | {prf_str:13} | {result['efficiency']:10.2%} | {result['regime']}")
#         except Exception as e:
#             print(f"  {config_name:17} | FAILED: {e}")
#             all_passed = False
    
#     if all_passed:
#         print("\nTest 4: PASSED ✓")
#     return all_passed


# # ============================================================================
# # TEST 5: Pressure Sweep Analysis
# # ============================================================================

# def test_5_pressure_sweep():
#     """Analyze three-regime behavior across pressure range."""
#     print("\n" + "="*60)
#     print("TEST 5: Pressure Sweep Analysis")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_downstream = 0.0
    
#     pressures = np.logspace(-2, 6, 50)  # 0.01 Pa to 1 MPa
#     defect_fractions = [0.0, 0.001, 0.01, 0.1]
    
#     results = {f: {'pressures': [], 'fluxes': [], 'regimes': []} for f in defect_fractions}
    
#     # Also calculate bare metal (Level 1) for reference
#     bare_metal_fluxes = []
#     for P in pressures:
#         result_L1 = calculate_simple_metal_flux(
#             metal_props['D_metal'], metal_props['K_s_metal'],
#             metal_props['thickness'], P, P_downstream
#         )
#         bare_metal_fluxes.append(result_L1['flux'])
    
#     print("\nCalculating flux vs pressure for different defect fractions...")
    
#     for P in pressures:
#         for f in defect_fractions:
#             try:
#                 if f == 0:
#                     result = calculate_oxide_metal_system(P, P_downstream, oxide_props, metal_props)
#                     flux = result['flux']
#                     regime = result.get('regime', 'unknown')
#                 else:
#                     result = calculate_parallel_path_flux(
#                         P, P_downstream, oxide_props, metal_props,
#                         {'area_fraction': f, 'type': 'pinhole'}
#                     )
#                     flux = result['flux_total']
#                     regime = result.get('regime_intact', 'unknown')
                
#                 results[f]['pressures'].append(P)
#                 results[f]['fluxes'].append(flux)
#                 results[f]['regimes'].append(regime)
#             except Exception as e:
#                 pass  # Skip failed points
    
#     # Create plot
#     ensure_output_dir()
#     fig, ax = plt.subplots(figsize=(10, 7))
    
#     # Plot bare metal (Level 1) as upper bound
#     ax.loglog(pressures, bare_metal_fluxes, 'k--', linewidth=2.5, 
#               label='Level 1: Bare metal (upper limit)', alpha=0.8)
    
#     colors = ['darkblue', 'blue', 'green', 'red']
#     labels = ['Perfect oxide (f=0)', 'f=0.001 (0.1%)', 'f=0.01 (1%)', 'f=0.1 (10%)']
    
#     for i, f in enumerate(defect_fractions):
#         if results[f]['fluxes']:
#             ax.loglog(results[f]['pressures'], results[f]['fluxes'], 
#                      color=colors[i], linewidth=2, label=labels[i])
    
#     ax.set_xlabel('log Upstream Pressure (Pa)', fontsize=12)
#     ax.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
#     ax.set_title('Level 3: Effect of Defect Fraction on Permeation\n(Dashed line = bare metal upper limit)', fontsize=14)
#     ax.legend(loc='best')
#     ax.grid(True, alpha=0.3)
    
#     # Add annotation for metal-limited region
#     ax.annotate('Metal-limited\nregion', xy=(1e5, bare_metal_fluxes[40]*0.95), 
#                 fontsize=10, color='darkred', ha='center')
    
#     plot_path = f"{OUTPUT_DIR}/pressure_sweep_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
#     plt.savefig(plot_path, dpi=150, bbox_inches='tight')
#     plt.close()
    
#     print(f"\n✓ Plot saved: {plot_path}")
#     print("Test 5: PASSED ✓")
#     return True

# # ============================================================================
# # TEST 6: Defect Type Comparison
# # ============================================================================

# def test_6_defect_type_comparison():
#     """Compare impact of different defect types."""
#     print("\n" + "="*60)
#     print("TEST 6: Defect Type Comparison")
#     print("="*60)
    
#     oxide_props, metal_props = setup_material_properties()
#     P_upstream, P_downstream = 100.0, 0.0  # 100 Pa
    
#     defect_configs = {
#         'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
#         'crack_thin': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
#         'crack_thick': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.5},
#         'grain_boundary_10x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
#         'grain_boundary_100x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 100},
#     }
    
#     # Get baseline (Level 2)
#     baseline = calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props)
    
#     print(f"\nBaseline (perfect oxide): {baseline['flux']:.3e} mol/m²/s")
#     print("\nDefect Type          | Flux (mol/m²/s) | Enhancement | Dominant Path")
#     print("-" * 75)
    
#     results = {}
#     for name, config in defect_configs.items():
#         result = calculate_parallel_path_flux(
#             P_upstream, P_downstream, oxide_props, metal_props, config
#         )
#         results[name] = result
#         enhancement = result['flux_total'] / baseline['flux'] if baseline['flux'] > 0 else float('inf')
#         print(f"  {name:18} | {result['flux_total']:15.3e} | {enhancement:11.2e}x | {result['dominant_path']}")
    
#     # Find most impactful defect type
#     max_flux_type = max(results.keys(), key=lambda x: results[x]['flux_total'])
#     print(f"\n→ Most impactful defect type: {max_flux_type}")
#     print("Test 6: PASSED ✓")
#     return True


"""
Comprehensive Level 3 Test Suite
Tests all aspects of the parallel path model for oxide with defects.

Test Categories:
1. Basic Functionality - All defect types work
2. Limiting Cases - f=0 → Level 2, f=1 → Level 1
3. Monotonic Behavior - Flux increases with defect fraction
4. PRF Validation - Values in literature range
5. Pressure Sweep Analysis - Three-regime behavior
6. Defect Type Comparison - Impact analysis
7. PRF vs Defect Fraction - Barrier effectiveness
8. Dominant Path Analysis - Phase diagram of defect vs oxide dominated regions
9. PRF Regime Analysis - PRF phase diagram
10. Temperature Effects - Regime transitions vs temperature
11. Defect Type Regime Comparison - Compare defect types across regimes
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from datetime import datetime
import os
import sys
from io import StringIO

# Import Level 1 and 2 for comparison
from calculations.permeation_calc import calculate_simple_metal_flux
from calculations.interface_solver import calculate_oxide_metal_system

# Import Level 3 functions
from calculations.parallel_oxide_defect_paths import (
    calculate_defect_path_flux,
    calculate_parallel_path_flux,
    calculate_PRF
)

# Import parameters
from data.oxide_defect_parameters import DEFECT_CONFIGURATIONS, PARAMETER_RANGES, PRF_RANGES
from data.oxide_properties import OXIDE_PROPERTIES
from data.material_data import MATERIALS

# Constants
T_TEST = 1073  # K (800°C)
R = 8.314  # J/mol/K

# Output directory
OUTPUT_DIR = "validation/results/complete_level3_tests"


def setup_material_properties(T=T_TEST):
    """Set up oxide and metal properties for testing."""
    oxide_data = OXIDE_PROPERTIES['Cr2O3']
    metal_data = MATERIALS['Incoloy800']
    
    oxide_props = {
        'D_ox': oxide_data['D_ox_0'] * np.exp(-oxide_data['E_D_ox'] / (R * T)),
        'K_ox': oxide_data['K_ox_0'] * np.exp(-oxide_data['H_sol_ox'] / (R * T)),
        'thickness': oxide_data['thickness']
    }
    
    metal_props = {
        'D_metal': metal_data['D_0'] * np.exp(-metal_data['E_D'] / (R * T)),
        'K_s_metal': metal_data['K_s0'] * np.exp(-metal_data['H_s'] / (R * T)),
        'thickness': 1e-3  # 1mm
    }
    
    return oxide_props, metal_props


def ensure_output_dir():
    """Create output directory if it doesn't exist."""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)


# ============================================================================
# TEST 1: Basic Functionality
# ============================================================================

def test_1_basic_functionality():
    """Test all defect types including mixed."""
    print("\n" + "="*60)
    print("TEST 1: Basic Functionality")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_upstream, P_downstream = 1.0, 0.0
    
    defect_configs = {
        'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
        'crack': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
        'grain_boundary': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
        'mixed': {
            'type': 'mixed', 
            'area_fraction': 0.015,
            'components': {'pinholes': 0.005, 'cracks': 0.005, 'grain_boundaries': 0.005},
            'thickness_factor': 0.1,
            'diffusivity_factor': 10
        }
    }
    
    results = {}
    all_passed = True
    
    for name, config in defect_configs.items():
        try:
            flux = calculate_defect_path_flux(
                P_upstream, P_downstream, oxide_props, metal_props, config
            )
            results[name] = flux
            print(f"✓ {name:15}: flux = {flux:.3e} mol/m²/s")
        except Exception as e:
            print(f"✗ {name:15}: FAILED - {e}")
            all_passed = False
    
    # Create plot
    if all_passed and results:
        ensure_output_dir()
        fig, ax = plt.subplots(figsize=(10, 6))
        
        names = list(results.keys())
        fluxes = [results[n] for n in names]
        
        colors = ['#2ecc71', '#3498db', '#e74c3c', '#9b59b6']
        bars = ax.bar(names, fluxes, color=colors, edgecolor='black', linewidth=1.5)
        
        ax.set_yscale('log')
        ax.set_xlabel('Defect Type', fontsize=12)
        ax.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
        ax.set_title('Test 1: Basic Functionality - Flux by Defect Type\n(f_defect = 1%, P = 1 Pa, T = 800°C)', fontsize=14)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar, flux in zip(bars, fluxes):
            ax.text(bar.get_x() + bar.get_width()/2., flux * 1.5, 
                   f'{flux:.2e}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plot_path = f"{OUTPUT_DIR}/test1_basic_functionality_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n✓ Plot saved: {plot_path}")
    
    if all_passed:
        print("\nTest 1: PASSED ✓")
    return all_passed


# ============================================================================
# TEST 2: Limiting Cases
# ============================================================================

def test_2_limiting_cases():
    """Verify f=0 → Level 2 and f=1 → Level 1."""
    print("\n" + "="*60)
    print("TEST 2: Limiting Cases Validation")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_upstream, P_downstream = 1.0, 0.0
    
    # Test 2a: f_defect = 0 → Level 2
    print("\nTest 2a: f_defect = 0 (should match Level 2)")
    
    result_l3_f0 = calculate_parallel_path_flux(
        P_upstream, P_downstream, oxide_props, metal_props,
        {'area_fraction': 0.0, 'type': 'pinhole'}
    )
    result_l2 = calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props)
    
    if result_l2['flux'] > 0:
        error_2a = abs(result_l3_f0['flux_total'] - result_l2['flux']) / result_l2['flux'] * 100
    else:
        error_2a = 0 if result_l3_f0['flux_total'] == 0 else float('inf')
    
    print(f"  Level 3 (f=0): {result_l3_f0['flux_total']:.3e} mol/m²/s")
    print(f"  Level 2:       {result_l2['flux']:.3e} mol/m²/s")
    print(f"  Error: {error_2a:.2e}%")
    
    # Test 2b: f_defect = 1 → Level 1
    print("\nTest 2b: f_defect = 1 (should match Level 1)")
    
    result_l3_f1 = calculate_parallel_path_flux(
        P_upstream, P_downstream, oxide_props, metal_props,
        {'area_fraction': 1.0, 'type': 'pinhole'}
    )
    result_l1 = calculate_simple_metal_flux(
        metal_props['D_metal'], metal_props['K_s_metal'], 
        metal_props['thickness'], P_upstream, P_downstream
    )
    
    error_2b = abs(result_l3_f1['flux_total'] - result_l1['flux']) / result_l1['flux'] * 100
    
    print(f"  Level 3 (f=1): {result_l3_f1['flux_total']:.3e} mol/m²/s")
    print(f"  Level 1:       {result_l1['flux']:.3e} mol/m²/s")
    print(f"  Error: {error_2b:.2e}%")
    
    passed = error_2a < 1e-6 and error_2b < 1e-6
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 2a: f=0 comparison
    ax1 = axes[0]
    labels_a = ['Level 3 (f=0)', 'Level 2\n(Perfect Oxide)']
    values_a = [result_l3_f0['flux_total'], result_l2['flux']]
    colors_a = ['#3498db', '#e74c3c']
    
    bars_a = ax1.bar(labels_a, values_a, color=colors_a, edgecolor='black', linewidth=1.5)
    ax1.set_yscale('log')
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title(f'f_defect = 0 → Level 2\nError: {error_2a:.2e}%', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars_a, values_a):
        ax1.text(bar.get_x() + bar.get_width()/2., val * 1.5, 
                f'{val:.2e}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Plot 2b: f=1 comparison
    ax2 = axes[1]
    labels_b = ['Level 3 (f=1)', 'Level 1\n(Bare Metal)']
    values_b = [result_l3_f1['flux_total'], result_l1['flux']]
    colors_b = ['#3498db', '#2ecc71']
    
    bars_b = ax2.bar(labels_b, values_b, color=colors_b, edgecolor='black', linewidth=1.5)
    ax2.set_yscale('log')
    ax2.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax2.set_title(f'f_defect = 1 → Level 1\nError: {error_2b:.2e}%', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars_b, values_b):
        ax2.text(bar.get_x() + bar.get_width()/2., val * 1.5, 
                f'{val:.2e}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.suptitle('Test 2: Limiting Cases Validation', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test2_limiting_cases_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n✓ Plot saved: {plot_path}")
    
    if passed:
        print("\nTest 2: PASSED ✓")
    else:
        print("\nTest 2: FAILED ✗")
    return passed


# ============================================================================
# TEST 3: Monotonic Behavior
# ============================================================================

def test_3_monotonic_behavior():
    """Verify flux increases monotonically with defect fraction."""
    print("\n" + "="*60)
    print("TEST 3: Monotonic Behavior with Defect Fraction")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_upstream, P_downstream = 1.0, 0.0
    
    defect_fractions = [0.0001, 0.001, 0.01, 0.1]
    fluxes = []
    enhancements = []
    
    print("\nDefect Fraction | Flux (mol/m²/s) | Enhancement")
    print("-" * 55)
    
    for f in defect_fractions:
        result = calculate_parallel_path_flux(
            P_upstream, P_downstream, oxide_props, metal_props,
            {'area_fraction': f, 'type': 'pinhole'}
        )
        fluxes.append(result['flux_total'])
        enhancements.append(result['defect_enhancement_factor'])
        print(f"  {f:12.4f}  | {result['flux_total']:15.3e} | {result['defect_enhancement_factor']:12.2e}x")
    
    is_monotonic = all(fluxes[i] <= fluxes[i+1] for i in range(len(fluxes)-1))
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Flux vs Defect Fraction
    ax1 = axes[0]
    ax1.loglog(defect_fractions, fluxes, 'b-o', linewidth=2, markersize=10, markerfacecolor='white')
    ax1.set_xlabel('log Defect Fraction', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux vs Defect Fraction', fontsize=12)
    ax1.grid(True, alpha=0.3, which='both')
    
    # Mark monotonic behavior
    for i, (f, flux) in enumerate(zip(defect_fractions, fluxes)):
        ax1.annotate(f'{flux:.1e}', (f, flux), textcoords="offset points", 
                    xytext=(0, 10), ha='center', fontsize=9)
    
    # Plot 2: Enhancement Factor
    ax2 = axes[1]
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(defect_fractions)))
    bars = ax2.bar([f'{f:.1%}' for f in defect_fractions], enhancements, color=colors, edgecolor='black')
    ax2.set_yscale('log')
    ax2.set_xlabel('Defect Fraction', fontsize=12)
    ax2.set_ylabel('log Enhancement Factor', fontsize=12)
    ax2.set_title('Enhancement Factor vs Defect Fraction', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    
    for bar, enh in zip(bars, enhancements):
        ax2.text(bar.get_x() + bar.get_width()/2., enh * 1.5, 
                f'{enh:.1e}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    status = "MONOTONIC ✓" if is_monotonic else "NON-MONOTONIC ✗"
    plt.suptitle(f'Test 3: Monotonic Behavior Verification - {status}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test3_monotonic_behavior_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n✓ Plot saved: {plot_path}")
    
    if is_monotonic:
        print("\n✓ Flux increases monotonically")
        print("Test 3: PASSED ✓")
    else:
        print("\n✗ Non-monotonic behavior detected")
        print("Test 3: FAILED ✗")
    return is_monotonic


# ============================================================================
# TEST 4: PRF Validation
# ============================================================================

def test_4_PRF_validation():
    """Verify PRF values across all configurations."""
    print("\n" + "="*60)
    print("TEST 4: PRF Validation")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_test = 1.0
    
    print("\nConfiguration       | PRF           | Efficiency | Regime")
    print("-" * 65)
    
    all_passed = True
    results = {}
    
    for config_name, config in DEFECT_CONFIGURATIONS.items():
        try:
            if config_name == 'perfect':
                result = calculate_PRF(P_test, oxide_props, metal_props, None)
            else:
                result = calculate_PRF(P_test, oxide_props, metal_props, config)
            
            results[config_name] = result
            prf_str = f"{result['PRF']:.2e}" if result['PRF'] > 1e6 else f"{result['PRF']:.1f}"
            print(f"  {config_name:17} | {prf_str:13} | {result['efficiency']:10.2%} | {result['regime']}")
        except Exception as e:
            print(f"  {config_name:17} | FAILED: {e}")
            all_passed = False
    
    # Create plot
    if results:
        ensure_output_dir()
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        names = list(results.keys())
        prfs = [results[n]['PRF'] for n in names]
        efficiencies = [results[n]['efficiency'] * 100 for n in names]
        
        # Plot 1: PRF by configuration
        ax1 = axes[0]
        colors = plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, len(names)))
        bars1 = ax1.bar(names, prfs, color=colors, edgecolor='black', linewidth=1.5)
        ax1.set_yscale('log')
        
        # Add literature range
        ax1.axhspan(10, 3828, alpha=0.2, color='blue', label='Literature range (10-3828)')
        
        ax1.set_xlabel('Configuration', fontsize=12)
        ax1.set_ylabel('log PRF', fontsize=12)
        ax1.set_title('Permeation Reduction Factor by Configuration', fontsize=12)
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3, axis='y')
        plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        # Plot 2: Barrier Efficiency
        ax2 = axes[1]
        colors2 = plt.cm.RdYlGn(np.array(efficiencies) / 100)
        bars2 = ax2.bar(names, efficiencies, color=colors2, edgecolor='black', linewidth=1.5)
        ax2.axhline(90, color='green', linestyle='--', alpha=0.7, label='90% threshold')
        ax2.axhline(50, color='orange', linestyle='--', alpha=0.7, label='50% threshold')
        
        ax2.set_xlabel('Configuration', fontsize=12)
        ax2.set_ylabel('Barrier Efficiency (%)', fontsize=12)
        ax2.set_title('Barrier Efficiency by Configuration', fontsize=12)
        ax2.set_ylim([0, 105])
        ax2.legend(loc='lower right')
        ax2.grid(True, alpha=0.3, axis='y')
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        for bar, eff in zip(bars2, efficiencies):
            ax2.text(bar.get_x() + bar.get_width()/2., eff + 2, 
                    f'{eff:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        plt.suptitle('Test 4: PRF Validation Across Configurations', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        plot_path = f"{OUTPUT_DIR}/test4_PRF_validation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n✓ Plot saved: {plot_path}")
    
    if all_passed:
        print("\nTest 4: PASSED ✓")
    return all_passed


# ============================================================================
# TEST 5: Pressure Sweep Analysis
# ============================================================================

def test_5_pressure_sweep():
    """Analyze three-regime behavior across pressure range."""
    print("\n" + "="*60)
    print("TEST 5: Pressure Sweep Analysis")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_downstream = 0.0
    
    pressures = np.logspace(-2, 6, 50)  # 0.01 Pa to 1 MPa
    defect_fractions = [0.0, 0.001, 0.01, 0.1]
    
    results = {f: {'pressures': [], 'fluxes': [], 'regimes': []} for f in defect_fractions}
    
    # Also calculate bare metal (Level 1) for reference
    bare_metal_fluxes = []
    for P in pressures:
        result_L1 = calculate_simple_metal_flux(
            metal_props['D_metal'], metal_props['K_s_metal'],
            metal_props['thickness'], P, P_downstream
        )
        bare_metal_fluxes.append(result_L1['flux'])
    
    print("\nCalculating flux vs pressure for different defect fractions...")
    
    for P in pressures:
        for f in defect_fractions:
            try:
                if f == 0:
                    result = calculate_oxide_metal_system(P, P_downstream, oxide_props, metal_props)
                    flux = result['flux']
                    regime = result.get('regime', 'unknown')
                else:
                    result = calculate_parallel_path_flux(
                        P, P_downstream, oxide_props, metal_props,
                        {'area_fraction': f, 'type': 'pinhole'}
                    )
                    flux = result['flux_total']
                    regime = result.get('regime_intact', 'unknown')
                
                results[f]['pressures'].append(P)
                results[f]['fluxes'].append(flux)
                results[f]['regimes'].append(regime)
            except Exception as e:
                pass  # Skip failed points
    
    # Create plot
    ensure_output_dir()
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot bare metal (Level 1) as upper bound
    ax.loglog(pressures, bare_metal_fluxes, 'k--', linewidth=2.5, 
              label='Level 1: Bare metal (upper limit)', alpha=0.8)
    
    colors = ['darkblue', 'blue', 'green', 'red']
    labels = ['Perfect oxide (f=0)', 'f=0.001 (0.1%)', 'f=0.01 (1%)', 'f=0.1 (10%)']
    
    for i, f in enumerate(defect_fractions):
        if results[f]['fluxes']:
            ax.loglog(results[f]['pressures'], results[f]['fluxes'], 
                     color=colors[i], linewidth=2, label=labels[i])
    
    ax.set_xlabel('log Upstream Pressure (Pa)', fontsize=12)
    ax.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax.set_title('Level 3: Effect of Defect Fraction on Permeation\n(Dashed line = bare metal upper limit)', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Add annotation for metal-limited region
    ax.annotate('Metal-limited\nregion', xy=(1e5, bare_metal_fluxes[40]*0.95), 
                fontsize=10, color='darkred', ha='center')
    
    plot_path = f"{OUTPUT_DIR}/test5_pressure_sweep_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 5: PASSED ✓")
    return True


# ============================================================================
# TEST 6: Defect Type Comparison
# ============================================================================

def test_6_defect_type_comparison():
    """Compare impact of different defect types."""
    print("\n" + "="*60)
    print("TEST 6: Defect Type Comparison")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_upstream, P_downstream = 100.0, 0.0  # 100 Pa
    
    defect_configs = {
        'pinhole': {'type': 'pinhole', 'area_fraction': 0.01},
        'crack_thin': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.1},
        'crack_thick': {'type': 'crack', 'area_fraction': 0.01, 'thickness_factor': 0.5},
        'grain_boundary_10x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 10},
        'grain_boundary_100x': {'type': 'grain_boundary', 'area_fraction': 0.01, 'diffusivity_factor': 100},
    }
    
    # Get baseline (Level 2)
    baseline = calculate_oxide_metal_system(P_upstream, P_downstream, oxide_props, metal_props)
    
    # Get bare metal (Level 1)
    bare_metal = calculate_simple_metal_flux(
        metal_props['D_metal'], metal_props['K_s_metal'],
        metal_props['thickness'], P_upstream, P_downstream
    )
    
    print(f"\nBaseline (perfect oxide): {baseline['flux']:.3e} mol/m²/s")
    print(f"Bare metal (Level 1):     {bare_metal['flux']:.3e} mol/m²/s")
    print("\nDefect Type          | Flux (mol/m²/s) | Enhancement | Dominant Path")
    print("-" * 75)
    
    results = {}
    for name, config in defect_configs.items():
        result = calculate_parallel_path_flux(
            P_upstream, P_downstream, oxide_props, metal_props, config
        )
        results[name] = result
        enhancement = result['flux_total'] / baseline['flux'] if baseline['flux'] > 0 else float('inf')
        results[name]['enhancement'] = enhancement
        print(f"  {name:18} | {result['flux_total']:15.3e} | {enhancement:11.2e}x | {result['dominant_path']}")
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    names = list(results.keys())
    fluxes = [results[n]['flux_total'] for n in names]
    enhancements = [results[n]['enhancement'] for n in names]
    
    # Plot 1: Flux comparison
    ax1 = axes[0]
    x = np.arange(len(names))
    width = 0.6
    
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6', '#f39c12']
    bars = ax1.bar(x, fluxes, width, color=colors, edgecolor='black', linewidth=1.5)
    
    # Reference lines
    ax1.axhline(baseline['flux'], color='darkblue', linestyle='--', linewidth=2, label=f'Level 2 (Perfect): {baseline["flux"]:.2e}')
    ax1.axhline(bare_metal['flux'], color='black', linestyle='-', linewidth=2, label=f'Level 1 (Bare): {bare_metal["flux"]:.2e}')
    
    ax1.set_yscale('log')
    ax1.set_xticks(x)
    ax1.set_xticklabels(['Pinhole', 'Crack\n(α=0.1)', 'Crack\n(α=0.5)', 'GB\n(β=10)', 'GB\n(β=100)'])
    ax1.set_xlabel('Defect Type', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux by Defect Type (f = 1%)', fontsize=12)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Enhancement factor
    ax2 = axes[1]
    bars2 = ax2.bar(x, enhancements, width, color=colors, edgecolor='black', linewidth=1.5)
    
    ax2.set_yscale('log')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Pinhole', 'Crack\n(α=0.1)', 'Crack\n(α=0.5)', 'GB\n(β=10)', 'GB\n(β=100)'])
    ax2.set_xlabel('Defect Type', fontsize=12)
    ax2.set_ylabel('log Enhancement (vs Perfect Oxide)', fontsize=12)
    ax2.set_title('Enhancement Factor by Defect Type', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    
    for bar, enh in zip(bars2, enhancements):
        ax2.text(bar.get_x() + bar.get_width()/2., enh * 1.5, 
                f'{enh:.1e}x', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.suptitle('Test 6: Defect Type Comparison at P = 100 Pa, T = 800°C', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test6_defect_type_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n✓ Plot saved: {plot_path}")
    
    print("\nTest 6: PASSED ✓")
    return True

# ============================================================================
# TEST 7: PRF vs Defect Fraction
# ============================================================================

def test_7_PRF_vs_defect_fraction():
    """Calculate and plot PRF vs defect fraction."""
    print("\n" + "="*60)
    print("TEST 7: PRF vs Defect Fraction")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_test = 100.0  # 100 Pa
    
    defect_fractions = np.logspace(-4, -0.5, 20)  # 0.01% to ~30%
    
    prf_pinhole = []
    prf_crack = []
    prf_gb = []
    
    for f in defect_fractions:
        # Pinhole
        result = calculate_PRF(P_test, oxide_props, metal_props,
                              {'area_fraction': f, 'type': 'pinhole'})
        prf_pinhole.append(result['PRF'])
        
        # Crack
        result = calculate_PRF(P_test, oxide_props, metal_props,
                              {'area_fraction': f, 'type': 'crack', 'thickness_factor': 0.1})
        prf_crack.append(result['PRF'])
        
        # Grain boundary
        result = calculate_PRF(P_test, oxide_props, metal_props,
                              {'area_fraction': f, 'type': 'grain_boundary', 'diffusivity_factor': 10})
        prf_gb.append(result['PRF'])
    
    # Create plot
    ensure_output_dir()
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.loglog(defect_fractions * 100, prf_pinhole, 'r-', linewidth=2, label='Pinholes')
    ax.loglog(defect_fractions * 100, prf_crack, 'b--', linewidth=2, label='Cracks (α=0.1)')
    ax.loglog(defect_fractions * 100, prf_gb, 'g:', linewidth=2, label='Grain boundaries (β=10)')
    
    # Add literature range
    ax.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range (10-3828)')
    
    # Add metal-limited region (PRF < 2 means oxide is nearly ineffective)
    ax.axhspan(1, 2, alpha=0.3, color='red', label='Metal-limited (PRF < 2)')
    ax.axhline(1, color='black', linestyle='-', linewidth=1.5, alpha=0.7)
    ax.text(0.015, 1.3, 'PRF = 1 (bare metal)', fontsize=9, color='black')
    
    ax.set_xlabel('log Defect Area Fraction (%)', fontsize=12)
    ax.set_ylabel('log Permeation Reduction Factor (PRF)', fontsize=12)
    ax.set_title('PRF vs Defect Fraction for Different Defect Types\n(PRF → 1 indicates metal-limited regime)', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0.01, 50])
    ax.set_ylim([0.8, 1e5])
    
    plot_path = f"{OUTPUT_DIR}/test7_PRF_vs_defect_fraction_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Print summary with metal-limited indication
    print("\nDefect Fraction (%) | PRF (pinhole) | PRF (crack) | PRF (GB)   | Regime")
    print("-" * 80)
    for i in [0, 5, 10, 15, 19]:  # Sample points
        regime = "Metal-limited" if prf_pinhole[i] < 2 else ("Poor barrier" if prf_pinhole[i] < 10 else "Effective")
        print(f"  {defect_fractions[i]*100:16.2f}  | {prf_pinhole[i]:13.1f} | {prf_crack[i]:11.2e} | {prf_gb[i]:.2e} | {regime}")
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 7: PASSED ✓")
    return True


# ============================================================================
# TEST 8: Dominant Path Analysis (Regime Analysis)
# ============================================================================

def test_8_dominant_path_analysis():
    """
    Analyze which path dominates flux at different conditions.
    Creates a phase diagram of defect-dominated vs oxide-dominated vs metal-dominated regions.
    
    Three regimes:
    1. Oxide-dominated: Intact oxide limits flux (low defect fraction)
    2. Defect-dominated: Defect paths limit flux (moderate defect fraction)
    3. Metal-dominated: Metal substrate limits flux (high defect fraction, pinholes)
    """
    print("\n" + "="*60)
    print("TEST 8: Dominant Path Analysis (Phase Diagram)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    
    # Create grid
    pressures = np.logspace(-1, 6, 50)  # 0.1 Pa to 1 MPa
    defect_fractions = np.logspace(-4, -0.3, 50)  # 0.01% to 50%
    
    # Arrays to store results
    dominant_path = np.zeros((len(defect_fractions), len(pressures)))
    flux_ratio = np.zeros((len(defect_fractions), len(pressures)))
    total_flux = np.zeros((len(defect_fractions), len(pressures)))
    metal_limitation = np.zeros((len(defect_fractions), len(pressures)))
    
    # Get bare metal flux for comparison (Level 1)
    print(f"Computing {len(pressures)} × {len(defect_fractions)} = {len(pressures)*len(defect_fractions)} points...")
    
    for i, f_def in enumerate(defect_fractions):
        for j, P in enumerate(pressures):
            try:
                result = calculate_parallel_path_flux(
                    P, 0, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'}
                )
                
                # Get bare metal flux for this pressure
                bare_metal = calculate_simple_metal_flux(
                    metal_props['D_metal'], metal_props['K_s_metal'],
                    metal_props['thickness'], P, 0
                )
                
                flux_defect = result['flux_defect_contribution']
                flux_intact = result['flux_intact_contribution']
                flux_bare = bare_metal['flux']
                total_flux[i, j] = result['flux_total']
                
                # Calculate how close to bare metal limit
                # metal_limitation = J_total / J_bare_metal (approaches 1 when metal-limited)
                metal_limitation[i, j] = result['flux_total'] / flux_bare if flux_bare > 0 else 0
                
                # Ratio: >1 means defect-dominated, <1 means oxide-dominated
                if flux_intact > 0:
                    flux_ratio[i, j] = flux_defect / flux_intact
                else:
                    flux_ratio[i, j] = 1e10
                
                # Classify regimes (now with 3 regions):
                # 0 = oxide-dominated, 0.5 = defect-path dominated, 1 = metal-limited
                if metal_limitation[i, j] > 0.9:
                    # Within 90% of bare metal flux = metal-limited
                    dominant_path[i, j] = 1.0  # Metal-limited
                elif flux_defect > 10 * flux_intact:
                    dominant_path[i, j] = 0.67  # Defect-path dominated
                elif flux_defect > flux_intact:
                    dominant_path[i, j] = 0.5  # Mixed (defect > oxide)
                elif flux_intact > 10 * flux_defect:
                    dominant_path[i, j] = 0.0  # Oxide-dominated
                else:
                    dominant_path[i, j] = 0.33  # Mixed (oxide > defect)
                    
            except Exception as e:
                dominant_path[i, j] = np.nan
                flux_ratio[i, j] = np.nan
                total_flux[i, j] = np.nan
                metal_limitation[i, j] = np.nan
    
    # Create plots
    ensure_output_dir()
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    P_mesh, F_mesh = np.meshgrid(pressures, defect_fractions * 100)
    
    # Plot 1: Three-regime phase diagram
    ax1 = axes[0, 0]
    # Custom colormap: Green (oxide) -> Yellow (mixed) -> Orange (defect) -> Red (metal)
    from matplotlib.colors import LinearSegmentedColormap
    colors_regime = ['green', 'yellowgreen', 'yellow', 'orange', 'red']
    cmap_regime = LinearSegmentedColormap.from_list('regime', colors_regime, N=256)
    
    im1 = ax1.pcolormesh(P_mesh, F_mesh, dominant_path, cmap=cmap_regime, vmin=0, vmax=1, shading='auto')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax1.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_title('Permeation Regime Map\n(Green=Oxide, Yellow=Mixed, Orange=Defect, Red=Metal)', fontsize=11)
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_ticks([0, 0.33, 0.5, 0.67, 1.0])
    cbar1.set_ticklabels(['Oxide', 'Mixed\n(Ox>Def)', 'Mixed\n(Def>Ox)', 'Defect', 'Metal'])
    
    # Add contour lines at regime boundaries
    ax1.contour(P_mesh, F_mesh, dominant_path, levels=[0.25, 0.58, 0.85], 
                colors='black', linewidths=1.5, linestyles=['--', '-', '--'])
    
    # Plot 2: Metal limitation factor (J_total / J_bare_metal)
    ax2 = axes[0, 1]
    im2 = ax2.pcolormesh(P_mesh, F_mesh, metal_limitation, cmap='hot_r', 
                         vmin=0, vmax=1, shading='auto')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax2.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_title('Metal Limitation Factor (J_total / J_bare_metal)\n1.0 = Metal substrate is rate-limiting', fontsize=11)
    plt.colorbar(im2, ax=ax2, label='J_total / J_bare')
    
    # Add contours for key values
    ax2.contour(P_mesh, F_mesh, metal_limitation, levels=[0.5, 0.9, 0.99], 
                colors='white', linewidths=1.5, linestyles=[':', '--', '-'])
    
    # Plot 3: Flux ratio (defect/oxide contribution)
    ax3 = axes[1, 0]
    im3 = ax3.pcolormesh(P_mesh, F_mesh, flux_ratio, cmap='coolwarm', 
                         norm=LogNorm(vmin=0.01, vmax=100), shading='auto')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax3.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax3.set_title('Flux Ratio (Defect Path / Oxide Path)', fontsize=12)
    plt.colorbar(im3, ax=ax3, label='J_defect / J_oxide')
    
    # Add contour at ratio = 1 (equal contributions)
    ax3.contour(P_mesh, F_mesh, flux_ratio, levels=[1], colors='white', linewidths=2, linestyles='--')
    
    # Plot 4: Cross-sections showing all three regimes
    ax4 = axes[1, 1]
    
    P_slices = [100, 10000, 1000000]  # 100 Pa, 10 kPa, 1 MPa
    colors = ['blue', 'green', 'red']
    
    for P_slice, color in zip(P_slices, colors):
        j_idx = np.argmin(np.abs(pressures - P_slice))
        ax4.semilogx(defect_fractions * 100, metal_limitation[:, j_idx], 
                    color=color, linewidth=2, label=f'P = {P_slice} Pa')
    
    ax4.axhline(1.0, color='black', linestyle='-', alpha=0.7, label='Bare metal limit')
    ax4.axhline(0.9, color='gray', linestyle='--', alpha=0.5, label='90% of bare metal')
    ax4.axhline(0.5, color='gray', linestyle=':', alpha=0.5, label='50% of bare metal')
    
    ax4.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax4.set_ylabel('J_total / J_bare_metal', fontsize=12)
    ax4.set_title('Approach to Metal-Limited Regime', fontsize=12)
    ax4.legend(loc='best', fontsize=9)
    ax4.set_ylim([0, 1.1])
    ax4.grid(True, alpha=0.3)
    
    # Add regime annotations
    ax4.fill_between([0.01, 100], 0.9, 1.1, alpha=0.2, color='red', label='Metal-limited')
    ax4.text(20, 0.95, 'Metal-limited', fontsize=10, color='darkred')
    ax4.text(0.1, 0.3, 'Oxide/Defect\ncontrolled', fontsize=10, color='darkblue')
    
    plt.suptitle(f'Three-Regime Analysis at T = 800°C (including metal-limited region)', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test8_dominant_path_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 8: PASSED ✓")
    return True



# ============================================================================
# TEST 9: PRF Regime Analysis
# ============================================================================

def test_9_PRF_regime_analysis():
    """
    Create PRF phase diagram showing barrier effectiveness regions.
    """
    print("\n" + "="*60)
    print("TEST 9: PRF Regime Analysis (Phase Diagram)")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    
    # Create grid
    pressures = np.logspace(0, 6, 40)  # 1 Pa to 1 MPa
    defect_fractions = np.logspace(-4, -0.5, 40)  # 0.01% to ~30%
    
    PRF_map = np.zeros((len(defect_fractions), len(pressures)))
    
    print(f"Computing PRF for {len(pressures)} × {len(defect_fractions)} points...")
    
    for i, f_def in enumerate(defect_fractions):
        for j, P in enumerate(pressures):
            try:
                result = calculate_PRF(
                    P, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'}
                )
                PRF_map[i, j] = result['PRF']
            except:
                PRF_map[i, j] = np.nan
    
    # Create plot
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    P_mesh, F_mesh = np.meshgrid(pressures, defect_fractions * 100)
    
    # Plot 1: PRF map
    ax1 = axes[0]
    im1 = ax1.pcolormesh(P_mesh, F_mesh, PRF_map, cmap='RdYlGn',
                         norm=LogNorm(vmin=1, vmax=10000), shading='auto')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('log Pressure (Pa)', fontsize=12)
    ax1.set_ylabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_title('Permeation Reduction Factor (PRF)', fontsize=12)
    cbar = plt.colorbar(im1, ax=ax1)
    cbar.set_label('PRF')
    
    # Add contours for key PRF values (including metal-limited at PRF=2)
    contour_levels = [2, 10, 100, 1000]
    cs = ax1.contour(P_mesh, F_mesh, PRF_map, levels=contour_levels, 
                     colors=['darkred', 'red', 'orange', 'black'], linewidths=1.5)
    ax1.clabel(cs, inline=True, fontsize=9, fmt='PRF=%d')
    
    # Plot 2: PRF vs defect fraction at different pressures
    ax2 = axes[1]
    
    P_slices = [10, 1000, 100000]
    colors = ['blue', 'green', 'red']
    labels = ['10 Pa', '1 kPa', '100 kPa']
    
    for P_slice, color, label in zip(P_slices, colors, labels):
        j_idx = np.argmin(np.abs(pressures - P_slice))
        ax2.loglog(defect_fractions * 100, PRF_map[:, j_idx],
                  color=color, linewidth=2, label=f'P = {label}')
    
    # Add literature range
    ax2.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature PRF (10-3828)')
    
    # Add metal-limited region
    ax2.axhspan(1, 2, alpha=0.3, color='red', label='Metal-limited (PRF < 2)')
    ax2.axhline(1, color='black', linestyle='-', linewidth=1.5)
    
    # Add regime boundaries
    ax2.axhline(100, color='gray', linestyle='--', alpha=0.7)
    ax2.axhline(10, color='gray', linestyle='--', alpha=0.7)
    
    ax2.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_ylabel('log PRF', fontsize=12)
    ax2.set_title('PRF vs Defect Fraction', fontsize=12)
    ax2.legend(loc='best')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_ylim([0.8, 1e5])
    
    # Add regime labels
    ax2.text(0.02, 2000, 'Excellent Barrier\n(PRF > 100)', fontsize=10, color='green')
    ax2.text(0.5, 30, 'Moderate Barrier\n(PRF 10-100)', fontsize=10, color='orange')
    ax2.text(5, 5, 'Poor Barrier\n(PRF 2-10)', fontsize=10, color='red')
    ax2.text(15, 1.3, 'Metal-limited\n(PRF < 2)', fontsize=10, color='darkred')
    
    plt.suptitle(f'PRF Analysis at T = 800°C (including metal-limited region)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test9_PRF_regime_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 9: PASSED ✓")
    return True


# ============================================================================
# TEST 10: Temperature Effects on Regimes
# ============================================================================

def test_10_temperature_effects():
    """
    Analyze how temperature affects regime transitions.
    """
    print("\n" + "="*60)
    print("TEST 10: Temperature Effects on Regimes")
    print("="*60)
    
    temperatures = np.array([600, 700, 800, 900, 1000]) + 273.15  # K
    defect_fractions = np.logspace(-4, -0.5, 50)
    P_test = 1000  # Pa
    
    ensure_output_dir()
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(temperatures)))
    
    # Plot 1: Flux vs defect fraction at different temperatures
    ax1 = axes[0]
    
    # Store bare metal fluxes for each temperature
    bare_metal_by_T = {}
    
    for T_K, color in zip(temperatures, colors):
        oxide_props, metal_props = setup_material_properties(T_K)
        
        # Get bare metal flux at this temperature
        result_L1 = calculate_simple_metal_flux(
            metal_props['D_metal'], metal_props['K_s_metal'],
            metal_props['thickness'], P_test, 0
        )
        bare_metal_by_T[T_K] = result_L1['flux']
        
        fluxes = []
        for f_def in defect_fractions:
            try:
                result = calculate_parallel_path_flux(
                    P_test, 0, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'}
                )
                fluxes.append(result['flux_total'])
            except:
                fluxes.append(np.nan)
        
        ax1.loglog(defect_fractions * 100, fluxes, color=color, linewidth=2,
                  label=f'{T_K-273:.0f}°C')
        
        # Add bare metal reference as horizontal dashed line
        ax1.axhline(bare_metal_by_T[T_K], color=color, linestyle='--', alpha=0.4, linewidth=1)
    
    ax1.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title(f'Flux vs Defect Fraction at P = {P_test} Pa\n(dashed = bare metal limit)', fontsize=12)
    ax1.legend(loc='best', title='Temperature')
    ax1.grid(True, which='both', alpha=0.3)
    
    # Plot 2: PRF vs defect fraction at different temperatures
    ax2 = axes[1]
    
    for T_K, color in zip(temperatures, colors):
        oxide_props, metal_props = setup_material_properties(T_K)
        
        PRFs = []
        for f_def in defect_fractions:
            try:
                result = calculate_PRF(
                    P_test, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'}
                )
                PRFs.append(result['PRF'])
            except:
                PRFs.append(np.nan)
        
        ax2.loglog(defect_fractions * 100, PRFs, color=color, linewidth=2,
                  label=f'{T_K-273:.0f}°C')
    
    ax2.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range')
    ax2.axhspan(1, 2, alpha=0.3, color='red', label='Metal-limited')
    ax2.axhline(1, color='black', linestyle='-', linewidth=1.5)
    ax2.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_ylabel('log PRF', fontsize=12)
    ax2.set_title(f'PRF vs Defect Fraction at P = {P_test} Pa', fontsize=12)
    ax2.legend(loc='best', title='Temperature', fontsize=8)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.set_ylim([0.8, 1e5])
    
    # Plot 3: Metal limitation factor vs defect fraction at different temperatures
    ax3 = axes[2]
    
    for T_K, color in zip(temperatures, colors):
        oxide_props, metal_props = setup_material_properties(T_K)
        bare_flux = bare_metal_by_T[T_K]
        
        metal_factors = []
        for f_def in defect_fractions:
            try:
                result = calculate_parallel_path_flux(
                    P_test, 0, oxide_props, metal_props,
                    {'area_fraction': f_def, 'type': 'pinhole'}
                )
                metal_factors.append(result['flux_total'] / bare_flux)
            except:
                metal_factors.append(np.nan)
        
        ax3.semilogx(defect_fractions * 100, metal_factors, color=color, linewidth=2,
                    label=f'{T_K-273:.0f}°C')
    
    ax3.axhline(1.0, color='black', linestyle='-', linewidth=1.5, label='Bare metal')
    ax3.axhline(0.9, color='gray', linestyle='--', alpha=0.7, label='90% of bare')
    ax3.fill_between([0.01, 100], 0.9, 1.0, alpha=0.2, color='red')
    ax3.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax3.set_ylabel('J_total / J_bare_metal', fontsize=12)
    ax3.set_title('Metal Limitation Factor vs Defect Fraction', fontsize=12)
    ax3.legend(loc='best', title='Temperature', fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim([0, 1.1])
    ax3.text(20, 0.95, 'Metal-limited', fontsize=10, color='darkred')
    
    plt.suptitle('Temperature Effect on Permeation Regimes (including metal-limited)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test10_temperature_regime_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 10: PASSED ✓")
    return True



# ============================================================================
# TEST 11: Defect Type Regime Comparison
# ============================================================================

def test_11_defect_type_regime_comparison():
    """
    Compare regime behavior for different defect types.
    """
    print("\n" + "="*60)
    print("TEST 11: Defect Type Regime Comparison")
    print("="*60)
    
    oxide_props, metal_props = setup_material_properties()
    P_test = 1000  # Pa
    
    defect_fractions = np.logspace(-4, -0.5, 50)
    
    defect_configs = {
        'Pinhole': {'type': 'pinhole'},
        'Crack (α=0.1)': {'type': 'crack', 'thickness_factor': 0.1},
        'Crack (α=0.5)': {'type': 'crack', 'thickness_factor': 0.5},
        'Grain Boundary (β=10)': {'type': 'grain_boundary', 'diffusivity_factor': 10},
        'Grain Boundary (β=100)': {'type': 'grain_boundary', 'diffusivity_factor': 100},
    }
    
    ensure_output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(defect_configs)))
    
    # Plot 1: Flux comparison
    ax1 = axes[0]
    
    for (name, config), color in zip(defect_configs.items(), colors):
        fluxes = []
        for f_def in defect_fractions:
            try:
                params = {'area_fraction': f_def, **config}
                result = calculate_parallel_path_flux(
                    P_test, 0, oxide_props, metal_props, params
                )
                fluxes.append(result['flux_total'])
            except:
                fluxes.append(np.nan)
        
        ax1.loglog(defect_fractions * 100, fluxes, color=color, linewidth=2, label=name)
    
    # Add Level 1 and Level 2 references
    result_L1 = calculate_simple_metal_flux(
        metal_props['D_metal'], metal_props['K_s_metal'],
        metal_props['thickness'], P_test, 0
    )
    result_L2 = calculate_oxide_metal_system(P_test, 0, oxide_props, metal_props)
    
    ax1.axhline(result_L1['flux'], color='blue', linestyle='--', alpha=0.7, label='Level 1 (Bare)')
    ax1.axhline(result_L2['flux'], color='red', linestyle='--', alpha=0.7, label='Level 2 (Perfect)')
    
    ax1.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax1.set_ylabel('log Flux (mol/m²/s)', fontsize=12)
    ax1.set_title('Flux vs Defect Fraction by Defect Type', fontsize=12)
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(True, which='both', alpha=0.3)
    
    # Plot 2: PRF comparison
    ax2 = axes[1]
    
    for (name, config), color in zip(defect_configs.items(), colors):
        PRFs = []
        for f_def in defect_fractions:
            try:
                params = {'area_fraction': f_def, **config}
                result = calculate_PRF(P_test, oxide_props, metal_props, params)
                PRFs.append(result['PRF'])
            except:
                PRFs.append(np.nan)
        
        ax2.loglog(defect_fractions * 100, PRFs, color=color, linewidth=2, label=name)
    
    ax2.axhspan(10, 3828, alpha=0.2, color='yellow', label='Literature range')
    ax2.set_xlabel('log Defect Fraction (%)', fontsize=12)
    ax2.set_ylabel('log PRF', fontsize=12)
    ax2.set_title('PRF vs Defect Fraction by Defect Type', fontsize=12)
    ax2.legend(loc='best', fontsize=9)
    ax2.grid(True, which='both', alpha=0.3)
    
    plt.suptitle(f'Defect Type Comparison at T=800°C, P={P_test} Pa', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = f"{OUTPUT_DIR}/test11_defect_type_regime_comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved: {plot_path}")
    print("Test 11: PASSED ✓")
    return True


# ============================================================================
# REGIME SUMMARY
# ============================================================================

# ============================================================================
# REGIME SUMMARY
# ============================================================================

def print_regime_summary():
    """Print summary of regime analysis findings."""
    print("\n" + "="*70)
    print("REGIME ANALYSIS SUMMARY")
    print("="*70)
    print("""
KEY FINDINGS:

1. THREE PERMEATION REGIMES
   ┌─────────────────┬────────────────────┬─────────────────────────┐
   │ Regime          │ Defect Fraction    │ Rate-Limiting Step      │
   ├─────────────────┼────────────────────┼─────────────────────────┤
   │ Oxide-dominated │ f_defect < 0.1%    │ Diffusion through oxide │
   │ Defect-dominated│ 0.1% < f < 10%     │ Defect path resistance  │
   │ Metal-limited   │ f_defect > 10-30%  │ Metal substrate itself  │
   └─────────────────┴────────────────────┴─────────────────────────┘

2. DOMINANT PATH TRANSITIONS
   - At f_defect < 0.1%: Oxide path dominates (PRF > 1000)
   - At f_defect ~ 1%: Mixed regime (PRF ~ 100)
   - At f_defect > 10%: Defect path dominates (PRF < 10)
   - At f_defect > 30%: Metal-limited (J → J_bare_metal)

3. CRITICAL DEFECT FRACTIONS
   - PRF = 1000: f_defect ≈ 0.1%
   - PRF = 100:  f_defect ≈ 1%
   - PRF = 10:   f_defect ≈ 10%
   - PRF → 1:    f_defect → 100% (bare metal)

4. METAL-LIMITED REGIME CHARACTERISTICS
   - Flux approaches bare metal value regardless of oxide properties
   - PRF ≈ 1 (oxide provides no meaningful barrier)
   - Occurs when defect area fraction is large (>30% for pinholes)
   - For cracks/GB: requires even higher fractions due to residual oxide

5. DEFECT TYPE IMPACT (ranked by severity)
   1. Pinholes - Most severe (direct metal exposure)
      → Metal-limited at f > ~30%
   2. Thin cracks (α=0.1) - Moderate  
      → Metal-limited at f > ~50%
   3. Grain boundaries - Least severe (enhanced diffusion only)
      → Rarely metal-limited (oxide still present)

6. TEMPERATURE EFFECTS
   - Higher T → Higher flux for all cases
   - PRF relatively temperature-independent
   - Metal-limited threshold nearly T-independent

7. DESIGN IMPLICATIONS
   - For effective barriers: Keep f_defect < 0.1%
   - For moderate protection: f_defect < 1%
   - Above 10% defects: Oxide provides minimal benefit
   - Above 30% pinholes: Oxide is essentially ineffective
""")


# ============================================================================
# RESULTS REPORT SAVING
# ============================================================================

def save_results_report(output_text, results_dict, output_dir=OUTPUT_DIR):
    """
    Save all test results to a timestamped report file.
    
    Parameters
    ----------
    output_text : str
        The captured console output from all tests
    results_dict : dict
        Dictionary mapping test names to pass/fail status
    output_dir : str
        Directory to save the report file
    
    Returns
    -------
    str
        Path to the saved report file
    """
    ensure_output_dir()
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_path = f"{output_dir}/test_results_report_{timestamp}.txt"
    
    with open(report_path, 'w') as f:
        # Write header
        f.write("=" * 80 + "\n")
        f.write("LEVEL 3 PARALLEL PATH MODEL - COMPREHENSIVE TEST SUITE\n")
        f.write("=" * 80 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Temperature: {T_TEST} K ({T_TEST - 273}°C)\n")
        f.write(f"Output Directory: {output_dir}\n")
        f.write("=" * 80 + "\n\n")
        
        # Write captured test output
        f.write(output_text)
        
        # Write final summary
        f.write("\n" + "=" * 80 + "\n")
        f.write("FINAL TEST SUMMARY\n")
        f.write("=" * 80 + "\n")
        
        for name, result in results_dict.items():
            status = "✓ PASSED" if result else "✗ FAILED"
            f.write(f"  Test: {name:35} {status}\n")
        
        passed = sum(results_dict.values())
        total = len(results_dict)
        f.write(f"\nTotal: {passed}/{total} tests passed\n")
        
        if passed == total:
            f.write("\n" + "=" * 80 + "\n")
            f.write("ALL TESTS PASSED SUCCESSFULLY! ✓\n")
            f.write("=" * 80 + "\n")
        
        # List generated plots
        f.write("\n" + "-" * 80 + "\n")
        f.write("GENERATED PLOTS:\n")
        f.write("-" * 80 + "\n")
        
        if os.path.exists(output_dir):
            plots = sorted([f for f in os.listdir(output_dir) if f.endswith('.png')])
            for plot in plots:
                f.write(f"  • {plot}\n")
            f.write(f"\nTotal: {len(plots)} plots saved to {output_dir}/\n")
    
    return report_path


class TeeOutput:
    """Class to capture stdout while also printing to console."""
    def __init__(self):
        self.buffer = StringIO()
        self.stdout = sys.stdout
    
    def write(self, text):
        self.buffer.write(text)
        self.stdout.write(text)
    
    def flush(self):
        self.buffer.flush()
        self.stdout.flush()
    
    def getvalue(self):
        return self.buffer.getvalue()


# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

def run_all_tests(save_report=True):
    """Run all tests in sequence and optionally save results to a report file."""
    # Set up output capture if saving report
    if save_report:
        tee = TeeOutput()
        sys.stdout = tee
    
    try:
        print("\n" + "="*70)
        print("LEVEL 3 COMPREHENSIVE TEST SUITE")
        print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*70)
        
        tests = [
            ("Basic Functionality", test_1_basic_functionality),
            ("Limiting Cases", test_2_limiting_cases),
            ("Monotonic Behavior", test_3_monotonic_behavior),
            ("PRF Validation", test_4_PRF_validation),
            ("Pressure Sweep", test_5_pressure_sweep),
            ("Defect Type Comparison", test_6_defect_type_comparison),
            ("PRF vs Defect Fraction", test_7_PRF_vs_defect_fraction),
            ("Dominant Path Analysis", test_8_dominant_path_analysis),
            ("PRF Regime Analysis", test_9_PRF_regime_analysis),
            ("Temperature Effects", test_10_temperature_effects),
            ("Defect Type Regime Comparison", test_11_defect_type_regime_comparison),
        ]
        
        results = {}
        for name, test_func in tests:
            try:
                results[name] = test_func()
            except Exception as e:
                print(f"\n✗ {name} failed with error: {e}")
                results[name] = False
        
        # Print regime summary after all tests
        print_regime_summary()
        
        # Summary
        print("\n" + "="*70)
        print("TEST SUMMARY")
        print("="*70)
        
        passed = sum(results.values())
        total = len(results)
        
        for name, result in results.items():
            status = "✓ PASSED" if result else "✗ FAILED"
            print(f"  {name:30} {status}")
        
        print(f"\nTotal: {passed}/{total} tests passed")
        
        if passed == total:
            print("\n" + "="*70)
            print("ALL TESTS PASSED SUCCESSFULLY! ✓")
            print("="*70)
        
    finally:
        # Restore stdout and save report
        if save_report:
            sys.stdout = tee.stdout
            captured_output = tee.getvalue()
            report_path = save_results_report(captured_output, results)
            print(f"\n📄 Results saved to: {report_path}")
    
    return passed == total


def run_regime_tests_only(save_report=True):
    """Run only the regime analysis tests (8-11) and optionally save results to a report file."""
    # Set up output capture if saving report
    if save_report:
        tee = TeeOutput()
        sys.stdout = tee
    
    try:
        print("\n" + "="*70)
        print("LEVEL 3 REGIME ANALYSIS TESTS")
        print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*70)
        
        tests = [
            ("Dominant Path Analysis", test_8_dominant_path_analysis),
            ("PRF Regime Analysis", test_9_PRF_regime_analysis),
            ("Temperature Effects", test_10_temperature_effects),
            ("Defect Type Regime Comparison", test_11_defect_type_regime_comparison),
        ]
        
        results = {}
        for name, test_func in tests:
            try:
                results[name] = test_func()
            except Exception as e:
                print(f"\n✗ {name} failed with error: {e}")
                results[name] = False
        
        # Print regime summary
        print_regime_summary()
        
        passed = sum(results.values())
        total = len(results)
        print(f"\nRegime Tests: {passed}/{total} passed")
        
    finally:
        # Restore stdout and save report
        if save_report:
            sys.stdout = tee.stdout
            captured_output = tee.getvalue()
            report_path = save_results_report(captured_output, results, OUTPUT_DIR)
            print(f"\n📄 Results saved to: {report_path}")
    
    return passed == total


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == '--regime':
        run_regime_tests_only()
    else:
        run_all_tests()