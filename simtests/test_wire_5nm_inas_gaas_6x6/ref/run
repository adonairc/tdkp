#####################################################
# kpXxX calculation for the quantum well
#####################################################


set gridfile     n1_0_1_both_msh.grd
set kpmethod     kp6x6
set numsol       16
set num_k_values 20
set ordering     "foreman"


# k space definition
set kmax [expr {(2.0 * $constants_pi / 0.5) * 0.12}]
set kmin 0

# setup 
set config [get_configuration_object];
set logmsg [get_logger_object];

# be extremly chatty
$logmsg set_level $LOG_INFO_DEVEL2;

# set configuration options
$config set "output_eigensolver_statistics" 		         0
$config set "assembly_build_nonsymmetric_matrices" 	         1
$config set "assembly_check_matrix_for_symmetry"   	         0
$config set "assembly_check_matrix_for_symmetry_tolerance" 1.0e-10
$config set "assembly_track_nonhermitian_nodes"                  0
$config set "assembly_save_matrices_to_file"		         0
$config set "output_kpmatrix_dump_to_file"                       0

if { $ordering == "symmetric" } {
  $config set "kpmatrix_disable_foreman_enable_symmetrization"     1
} else {
  $config set "kpmatrix_disable_foreman_enable_symmetrization"     0
}
# set up material database
set material_database [MaterialDatabase];

# create i/o parser
set parser [InputParser]

# read geometry/grid file
set geometry [$parser read_geometry_from_file $gridfile]

# create problem
switch $kpmethod {
    effmass {
	set problem [EffectiveMass1D2D]
        $problem set_geometry $geometry
    }
    kp4x4 {
 	set problem [KP4x41D2D]
        $problem set_geometry $geometry
        $problem set_axes $ez $ex $ey
    }
    kp6x6 {
	set problem [KP6x61D2D]
	$problem set_geometry $geometry
        $problem set_axes $ez $ex $ey
    }
    kp8x8h {
	set problem [KP8x81D2D]	
	$problem set_geometry $geometry
	$problem set_solution_type $holes
        $problem set_axes $ez $ex $ey
        $problem set_energy_barrier 0.2
    }
    kp8x8e {
	set problem [KP8x81D2D]
	$problem set_geometry $geometry
	$problem set_solution_type $electrons     
        $problem set_axes $ez $ex $ey
        $problem set_energy_barrier 0.2
    }
    default {
	puts " **error** unknown solution method ${kpmethod}"
	exit
    }
}

# add geometry and material
$problem set_material_db $material_database

# solve problem and display some information
if { $kpmethod == "effmass" } {
  $problem solve $numsol
} else {
  $problem solve $numsol $kmin $kmax $num_k_values
}

# get bandstructure object
set bandstructure [$problem get_bandstructure 0];
Bandstructurecomplex -this $bandstructure

# write xy plot file
set xydatafile "bandstructure_${kpmethod}.dat"
$parser write_ascii_cplx $bandstructure $xydatafile
$parser write_ascii_cplx2real  $bandstructure "real_${xydatafile}"

