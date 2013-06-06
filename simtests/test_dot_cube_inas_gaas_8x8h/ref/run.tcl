#####################################################
# kpXxX calculation for the infinite quantum dot
#####################################################


set kpmethod kp8x8h
set numsol 16
set gridfile dot_1483_verts.grd

# check if gridfile exists

# setup 
set config [get_configuration_object];
set logmsg [get_logger_object];

# be extremly chatty
$logmsg set_level $LOG_INFO_DEVEL2;

# set configuration options
$config set "output_eigensolver_statistics" 		0
$config set "assembly_build_nonsymmetric_matrices" 	1
$config set "assembly_check_matrix_for_symmetry"   	0
$config set "assembly_save_matrices_to_file"		0

# set up material database
set material_database [MaterialDatabase];

# create i/o parser
set parser [InputParser]

# read geometry/grid file
set geometry [$parser read_geometry_from_file $gridfile]

# create problem
switch $kpmethod {
    effmass {
	set problem [EffectiveMass3D]
        $problem set_geometry $geometry
    }
    kp4x4 {
 	set problem [KP4x43D]
        $problem set_geometry $geometry
        $problem set_axes $ex $ey $ez
    }
    kp6x6 {
	set problem [KP6x63D]
	$problem set_geometry $geometry
        $problem set_axes $ex $ey $ez
    }
    kp8x8h {
	set problem [KP8x83D]	
	$problem set_geometry $geometry
	$problem set_solution_type $holes
        $problem set_axes $ex $ey $ez
    }
    kp8x8e {
	set problem [KP8x83D]
	$problem set_geometry $geometry
	$problem set_solution_type $electrons
        $problem set_axes $ex $ey $ez
    }
    default {
	puts " **error** unknown solution method ${kpmethod}"
	exit
    }
}

# add geometry and material
$problem set_material_db $material_database

# solve problem and display some information
$problem solve $numsol

# store (not all) energies in file
set fo [open "energies.dat" w]
for {set ii 0} {$ii < [expr {$numsol * 0.75}]} {incr ii} {
  set re [complex2real [$problem get_eigenvalue $ii]]
  set im [complex2imag [$problem get_eigenvalue $ii]]
  set str "${re} "
  puts $fo $str
}

close $fo

