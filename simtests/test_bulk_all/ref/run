###########################################################
# calcualte kp8x8 bulk bandstructure
###########################################################

# setup 
set config [get_configuration_object];
set logmsg [get_logger_object];

set material_name "GaAs"

# be extremly chatty
$logmsg set_level $LOG_INFO_DEVEL2;

# set up material database
set material_database [MaterialDatabase];
# load material
$material_database load_material $material_name;

# create i/o parser
set parser [InputParser];

# create rotation matrix
set rotation [list 1 0 0 0 1 0 0 0 1]
set rotation_matrix [RotationMatrix -args 3 3]

for {set ii 0} {$ii < 3} {incr ii} {
  for {set jj 0} {$jj < 3} {incr jj} {
    $rotation_matrix set $ii $jj [lindex $rotation [expr {$ii * 3 + $jj}]]
  }
}  


# create instance of kp matrix
set kpmatrix8x8 [new_KPMatrix8x8EndersForeman];
set kpmatrix4x4 [new_KPMatrix4x4EndersForeman];
set kpmatrix6x6 [new_KPMatrix6x6EndersForeman];

# set material to kp matrix
$kpmatrix4x4 set_material [$material_database get_material $material_name]
$kpmatrix4x4 set_rotation $rotation_matrix
# set material to kp matrix
$kpmatrix6x6 set_material [$material_database get_material $material_name]
$kpmatrix6x6 set_rotation $rotation_matrix
# set material to kp matrix
$kpmatrix8x8 set_material [$material_database get_material $material_name]
$kpmatrix8x8 set_rotation $rotation_matrix

# set direction of bandstructure
set direction [Vector3D -args 0.0 0.0 1.0];

# set bandstructure range
set kmin 0;
set kmax [expr {0.15 * (2.0 * $constants_pi / 0.5)}]; 
set nk   102;

# create bulk solver
set solver [BulkBandstructureSolver -args $kpmatrix4x4]
# get bands
$solver solve $direction $kmin $kmax $nk
# write file bands.dat
$parser write_ascii_cplx2real [$solver get_bandstructure 0] "bands4x4.dat";


# create bulk solver
set solver [BulkBandstructureSolver -args $kpmatrix6x6]
# get bands
$solver solve $direction $kmin $kmax $nk
# write file bands.dat
$parser write_ascii_cplx2real [$solver get_bandstructure 0] "bands6x6.dat";


# create bulk solver
set solver [BulkBandstructureSolver -args $kpmatrix8x8]
# get bands
$solver solve $direction $kmin $kmax $nk
# write file bands.dat
$parser write_ascii_cplx2real [$solver get_bandstructure 0] "bands8x8.dat";


# calculate density of states
set dos_calculator [DensityOfStates -args 3]
$dos_calculator set_bandstructure 8 [$solver get_bandstructure 0]
$dos_calculator calculate

set totaldos [$dos_calculator get_full_dos]

set estart -0.2
set estop  1.6
set epoint 100
set edelta [expr {($estop - $estart) / ($epoint - 1)}]

set fo [open "dos001.dat" w]
for {set ii 0} {$ii < $epoint} {incr ii} {
  set energy [expr {$estart + $ii * $edelta}]
  set tdos   [$totaldos EvalAt $energy]
  puts $fo "$energy \t $tdos"
}

close $fo

