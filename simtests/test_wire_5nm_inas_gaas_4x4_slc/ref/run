#####################################################
# kpXxX calculation for the quantum well
#####################################################

set ONLY_PROC_DEFINITIONS 1
source "../../../input/default_config.tcl"
source "../../../input/tdkp_helpers.tcl"
source "../../../input/matrix_elements.tcl"

set gridfile     n1_wurst_tensor.grd
set numsol       16
set num_k_values 20

# k space definition
set kmax [expr {(2.0 * $constants_pi / 0.5) * 0.20}]
set kmin 0

# setup
set config [get_configuration_object];
set logmsg [get_logger_object];

# be extremly chatty
$logmsg set_level $LOG_INFO_DEVEL2;

# set configuration options
source "../../configuration.tcl"

# set up material database
set material_database [MaterialDatabase];

# create i/o parser
set parser [InputParser]

# read geometry/grid file
set geometry [$parser read_geometry [DFISEGridReader -args $gridfile]]
$geometry set_materials $material_database



if {1 == 1} {

# solve for electrons
set problem [EffectiveMass -args $geometry $material_database]


$problem set_solution_type $electrons
$problem solve $numsol
set cb_bands [$problem get_bandstructure];
BandstructureDomaincomplex -this $cb_bands

# solve for holes
set problem [KP4x41D2D -args $geometry $material_database]

$problem set_energy_guess 0 0.0
$problem set_axes $ez $ex $ey
$problem solve $numsol $kmin $kmax $num_k_values
# get bandstructure object
set vb_bands [$problem get_bandstructure 0];
BandstructureDomaincomplex -this $vb_bands

#$cb_bands write_binary "cb_bands.dat"
#$vb_bands write_binary "vb_bands.dat"

} else {
    set vb_bands [BandstructureDomaincomplex -args "vb_bands.dat"]
    set cb_bands [BandstructureDomaincomplex -args "cb_bands.dat"]
    $material_database load_material "GaAs"
    $material_database load_material "InAs"
}

# get bulk material
set bulk  [$material_database get_material "GaAs"]
set quant [$material_database get_material "InAs"]

# create cb dispersion
set eff_tr [$quant get "electron_effective_mass_transverse"]
set domain [$vb_bands get_domain]
set cb_bands_disp [create_effmass_dispersion  $eff_tr $domain $cb_bands]

# get num valid and bound
set cb_num_valid [find_num_valid_states $cb_bands $electrons]
set vb_num_valid [find_num_valid_states $vb_bands $holes]
set cb_num_bound [find_num_bound_states $cb_bands [$bulk get "conduction_band_edge"] $electrons]
set vb_num_bound [find_num_bound_states $vb_bands [$bulk get "valence_band_edge"] $holes]
set cb_num_bands [min $cb_num_valid $cb_num_bound]
set vb_num_bands [min $vb_num_valid $vb_num_bound]
puts "cb: valid = ${cb_num_valid}, bound = ${cb_num_bound}, taking = ${cb_num_bands}"
puts "vb: valid = ${vb_num_valid}, bound = ${vb_num_bound}, taking = ${vb_num_bands}"

# store bands into ascii file (for testing)
$parser write_ascii_cplx2real $cb_bands_disp "cb_bands_ascii.dat"
$parser write_ascii_cplx2real $vb_bands      "vb_bands_ascii.dat"

# create momentum operator and matrix elements class
set momentum_operator [MomentumOperator4x4 -args $geometry $material_database]
set matrix_elements   [MatrixElements -args $momentum_operator]

# calculate matrix elements
$matrix_elements calculate $cb_num_bands $vb_num_bands $cb_bands $vb_bands

# find degenerate transitions and dump file transitions
set transitions [matrix_elements_find_degenerate_transitions $matrix_elements $cb_bands $vb_bands]
matrix_elements_write_degenerate_to_file $matrix_elements $transitions "transitions.dat" "transition_matrix.dat"

# build slc
set gaynous_broadening 1.0e13
set refractive_index 3.5
set well_width 5.0
set cb_fermi 0.7
set vb_fermi 0.0
set temp 300.0
set slc [SLC -args $cb_bands_disp $vb_bands $matrix_elements $gaynous_broadening $refractive_index $well_width]

$slc calculate $cb_fermi $vb_fermi $temp

$parser write_ascii_double $slc "luminescence.dat"

quit


