#####################################################
# kpXxX calculation for the infinite quantum dot
#####################################################


set kpmethod kp6x6
set numsol 16
set gridfile dot_1483_verts.asc.gz

# check if gridfile exists

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
set geometry [$parser read_geometry [AsciiGridReader -args $gridfile]]
$geometry set_materials $material_database

# create problem
switch $kpmethod {
    effmass {
	set problem [EffectiveMass -args $geometry $material_database]

    }
    kp4x4 {
 	set problem [KP4x43D -args $geometry $material_database]

        $problem set_axes $ex $ey $ez
    }
    kp6x6 {
	set problem [KP6x63D -args $geometry $material_database]

        $problem set_axes $ex $ey $ez
    }
    kp8x8h {
	set problem [KP8x83D -args $geometry $material_database]

	$problem set_solution_type $holes
        $problem set_axes $ex $ey $ez
    }
    kp8x8e {
	set problem [KP8x83D -args $geometry $material_database]

	$problem set_solution_type $electrons
        $problem set_axes $ex $ey $ez
    }
    default {
	puts " **error** unknown solution method ${kpmethod}"
	exit
    }
}

# add geometry and material


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
quit

