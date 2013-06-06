
source "../../configuration.tcl"

$config set "clc_polarization_start" 1
$config set "clc_polarization_end"   1
$config set "clc_include_exchange_shifts" 1
$config set "clc_use_screening" 1
$config set "clc_bands_extrapolate_kmax_n_factor" 3
$config set "clc_gmres_tolerance" 1.0e-9
$config set "clc_solve_self_consistent_lambda_k" 2
$config set "clc_radial_angular_integration_well_num_points" 60

set num_k_values   64
set kmax           2.4
set nitrides_use_vurgaftman_2007 1
set slc_delta_k 0.015

set wurzite_strain_axis_1  1
set wurzite_strain_axis_2  2
set wurzite_strain_axis_3  0

set well_transverse_direction [Vector3D -args 1.0 1.0 0.0];
set well_quantized_direction  [Vector3D -args 0.0 0.0 1.0];

set bulk_transverse_direction [Vector3D -args 1.0 1.0 1.0];

set clc_static_dielectric_constant 10.3
set slc_refractive_index 2.27
set slc_homogenous_broadening [expr {6.6 / $constants_hbar / 1000.0}]
#set slc_homogenous_broadening [expr {9.0 / $constants_hbar / 1000.0}]
#set slc_homogenous_broadening [expr {4.0 / $constants_hbar / 1000.0}]

set slc_omega_range_mode        1
set slc_omega_min               [expr { 2.5 / $constants_hbar } ]
set slc_omega_max               [expr { 3.0 / $constants_hbar } ]
#set slc_omega_min               [expr { 2.7 / $constants_hbar } ]
#set slc_omega_max               [expr { 3.4 / $constants_hbar } ]

set coulomb_q_min               1.0e-3
set coulomb_q_max               [expr {$kmax * 2}]
set coulomb_q_num               166
set coulomb_progression         3.0

set slc_omega_num               300

$config set "kpmatrix_wurtzite_include_spin_orbit_interaction" 0
$config set "kpmatrix_ignore_second_order_strain_dependence" 1
$config set "kpmatrix_ignore_first_order_strain_dependence" 1
$config set "kpmatrix_include_spin_orbit_strain_dependence" 0
$config set "assembly_parallel_energy_guess_vb_poly_coeff_a" 0.0
$config set "assembly_parallel_energy_guess_vb_poly_coeff_b" 0.0
$config set "slc_selected_lineshape_function" 2
