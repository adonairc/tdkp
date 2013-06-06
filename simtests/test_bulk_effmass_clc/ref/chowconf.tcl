

$config set "clc_polarization_start" 1
$config set "clc_include_exchange_shifts" 1
$config set "clc_use_screening" 1
$config set "clc_bands_extrapolate_kmax_n_factor" 3
$config set "clc_gmres_tolerance" 1.0e-9
$config set "clc_solve_self_consistent_lambda_k" 2
$config set "output_clc_screening_epsilon" 1
$config set "output_clc_exchange_shifts"  1

set num_k_values   100
set kmax           0.9
set nitrides_use_vurgaftman_2007 1
set slc_delta_k 0.005

set wurzite_strain_axis_1  1
set wurzite_strain_axis_2  2
set wurzite_strain_axis_3  0

set well_transverse_direction [Vector3D -args 1.0 1.0 0.0];
set well_quantized_direction  [Vector3D -args 0.0 0.0 1.0];

set bulk_transverse_direction [Vector3D -args 1.0 1.0 1.0];

set clc_static_dielectric_constant 12.0
set slc_refractive_index 3.6
set slc_homogenous_broadening [expr {6.6 / $constants_hbar / 1000.0}]

# (0.473)^2 * 2 / m0 * (m0*m0*(1.5/hbar^2))
set effmass_optical_matrix_param 8.8082

set slc_omega_range_mode        1
set slc_omega_min               [expr { 1.3 / $constants_hbar } ]
set slc_omega_max               [expr { 1.8 / $constants_hbar } ]

set coulomb_q_min               1.0e-3
set coulomb_q_max               [expr {$kmax * 2}]
set coulomb_q_num               166
set coulomb_progression         3.0

set slc_omega_num               300

