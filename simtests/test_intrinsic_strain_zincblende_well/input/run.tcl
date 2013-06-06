
set gridfile "well.asc.gz"

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

# let right side float
[$geometry get_node 384] set_index_internal 383

set problem [IntrinsicStrain -args $geometry $material_database]

# add geometry and material


$problem set_reference_lattice_constant 5.6503

# define x-axis to be unstrained
$problem set_strained_axes 0 0

$problem solve 1

set field [$problem get_strain_field]

# produce strain vertex file
$parser write_ascii_double $geometry $field "element_strain.dat"


quit


