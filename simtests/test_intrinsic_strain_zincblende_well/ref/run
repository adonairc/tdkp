
set gridfile "well.grd"

# setup
set config [get_configuration_object];
set logmsg [get_logger_object];


# be extremly chatty
$logmsg set_level $LOG_INFO_DEVEL2;

# set configuration options
$config set "output_eigensolver_statistics"                      0
$config set "assembly_build_nonsymmetric_matrices"               1
$config set "assembly_check_matrix_for_symmetry"                 1
$config set "assembly_check_matrix_for_symmetry_tolerance" 1.0e-10
$config set "assembly_track_nonhermitian_nodes"                  0
$config set "assembly_save_matrices_to_file"                     0
$config set "output_kpmatrix_dump_to_file"                       0


# set up material database
set material_database [MaterialDatabase];

# create i/o parser
set parser [InputParser]

# read geometry/grid file
set geometry [$parser read_geometry_from_file $gridfile]

# let right side float 
[$geometry get_vertex 384] set_index_internal 383

set problem [IntrinsicStrain]
$problem set_geometry $geometry
# add geometry and material
$problem set_material_db $material_database

$problem set_reference_lattice_constant 5.6503

# define x-axis to be unstrained
$problem set_strained_axes 0 0

$problem solve 1

set field [$problem get_strain_field]

# produce strain vertex file
$parser write_ascii_double $geometry $field "element_strain.dat"



       
