
# read in potential energy field
set potfield [PotentialEnergyField -args [$geometry get_num_elements]]
#  Slurp up the data file
set fp [open "potential${idx}.pot" r]
set data [read $fp]
close $fp
#  Process data file
set data [split $data "\n"]
set ii -1
set cnt 0
set avg 0.0
set nelem [$geometry get_num_elements]
foreach line $data {
   if {$ii >= 0 && $line != ""} {
      $potfield set_element_value $ii 0 [expr {($line + $last_line) / 2.0}]
      set avg [expr {$avg + $line}]
      incr cnt
   }
   set last_line $line
   incr ii
}

set avg [expr {$avg / $cnt}]



$problem set_field $potfield
$problem set_energy_barrier [expr {-5.1 + $avg}]

$problem set_solution_type $electrons     
set offset [expr {$constants_hbar * $constants_hbar * $kmax * $kmax / (2 * $constants_m0 *0.067 * ($num_k_values - 1))}]
#$problem set_target_energy_offset $offset
$problem set_energy_guess 0 [expr {-4.8 + $avg}]
$problem set_axes $kt_dir $kq_dir

$problem solve $numsol $kmin $kmax $num_k_values

# get bandstructure object
set bandstructure [$problem get_bandstructure 0];
BandstructureSingleDispersionComplex -this $bandstructure
#Bandstructurecomplex -this $bandstructure

# write xy plot file
set xydatafile "bandstructure_electrons_${idx}.dat"
$parser write_ascii_cplx       $bandstructure $xydatafile
$parser write_ascii_cplx2real  $bandstructure "real_${xydatafile}"
$problem delete_solutions

$problem set_solution_type $holes
$problem set_energy_guess 0 [expr {-5.20 + $avg}]

$problem solve $numsol $kmin $kmax $num_k_values

# get bandstructure object
set bandstructure [$problem get_bandstructure 0];
BandstructureSingleDispersionComplex -this $bandstructure
#Bandstructurecomplex -this $bandstructure

# write xy plot file
set xydatafile "bandstructure_holes_${idx}.dat"
$parser write_ascii_cplx       $bandstructure $xydatafile
$parser write_ascii_cplx2real  $bandstructure "real_${xydatafile}"
$problem delete_solutions
