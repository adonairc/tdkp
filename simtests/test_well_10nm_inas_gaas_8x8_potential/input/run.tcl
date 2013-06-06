#####################################################
# kpXxX calculation for the quantum well
#####################################################

set gridfile 		well_center_kp.asc.gz

set numsol   		14
set num_k_values 	16

# k space definition
set kmax 1.5
set kmin 0

# transversal direction [010]
set kt_dir $ez
# quantization dir      [100]
set kq_dir $ex

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
$geometry rescale_node_coordinates 1000.0
$geometry prepare

set problem [KP8x81D2D -args $geometry $material_database]
# add geometry and material



foreach {idx} {27 37 41 44} {
   source ./calc.tcl
}



quit