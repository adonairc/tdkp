#####################################################
# kpXxX calculation for the quantum well
#####################################################

set gridfile 		well_center_kp.grd

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
$config set "output_eigensolver_statistics" 		                 0
$config set "assembly_build_nonsymmetric_matrices" 	           1
$config set "assembly_check_matrix_for_symmetry"   	           0
$config set "assembly_check_matrix_for_symmetry_tolerance" 1.0e-10
$config set "assembly_track_nonhermitian_nodes"                  0
$config set "assembly_save_matrices_to_file"		                 0
$config set "output_kpmatrix_dump_to_file"                       0

# set up material database
set material_database [MaterialDatabase];

# create i/o parser
set parser [InputParser]

# read geometry/grid file
set geometry [$parser read_geometry_from_file $gridfile]
$geometry rescale_vertex_coordinates 1000.0
$geometry prepare

set problem [KP8x81D2D]  
# add geometry and material
$problem set_geometry $geometry
$problem set_material_db $material_database

foreach {idx} {27 37 41 44} {
   source ./calc.tcl
}



quit