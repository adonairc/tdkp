#####################################################
# kpXxX calculation for the quantum well
#####################################################

# --------------------------------------------------------------------
# check matrix elements bands and build list of degenerate transitions
# so basically we return a list of transitions, where each transition
# has two lists of cb and vb indices
# --------------------------------------------------------------------
proc matrix_elements_find_degenerate_transitions {matrix_elements cb_bands vb_bands} {

    set cb_num [$matrix_elements get_num_cb_bands]
    set vb_num [$matrix_elements get_num_vb_bands]

    global band_degeneracy_threshold

    set transitions  [list]
    set new_vb_trans [list]
    set new_cb_trans [list]

    for {set cc 0} {$cc < $cb_num} {incr cc} {
        # collect nondegenerate cb bands
        set new_cb_trans [list $cc]
     #   incr cc
     #   lappend new_cb_trans $cc
        while {$cc < [expr {$cb_num-1}] && [check_degenerate_bands $cb_bands $cc [expr {$cc + 1}] $band_degeneracy_threshold] == 1} {
            incr cc
            lappend new_cb_trans $cc
        }
        for {set vv 0} {$vv < $vb_num} {incr vv} {
            lappend new_vb_trans $vv
      #      incr vv
      #      lappend new_vb_trans $vv
            # if last or next energy differs, store transition and reset
            if {[expr {$vv + 1 == $vb_num}] || [check_degenerate_bands $vb_bands $vv [expr {$vv + 1}] $band_degeneracy_threshold] == 0} {
                set tmp [list $new_cb_trans $new_vb_trans]
                lappend transitions $tmp
                set new_vb_trans [list]
            }
        }
    }
    return $transitions
}

# --------------------------------------------------------------------
# compare two dispersion relations if they are equal (degenerate bands)
# (equal -> return == 1, else return 0)
# --------------------------------------------------------------------
proc check_degenerate_bands {bands idx_1 idx_2 tol} {
    set numk [$bands get_number_of_k_values]
    for {set kk 0} {$kk < $numk} {incr kk} {
        set e1 [complex2real [$bands get_energy $kk $idx_1]]
        set e2 [complex2real [$bands get_energy $kk $idx_2]]
        if {[expr {abs($e1 - $e2)}] > $tol} {
            return 0
        }
    }
    return 1
}

# --------------------------------------------------------------------
# plot transition matrix
# transition matrix is a file with 3 matrices where
# for each ii is the degenerate vb idx and jj is the degenerate cb idx
# the values of the matrix are the strength of transitions at k = 0
# normalized such that 100 is the maximum
# --------------------------------------------------------------------
proc write_transition_matrix_file {transition_values transitions matrixfile} {

    set max_cb 0
    set max_vb 0
    # find maximum transition value
    set max_values {0.0 0.0 0.0}
    for {set ii 0} {$ii < [llength $transition_values]} {incr ii} {
        set pp [lindex $transition_values $ii]
        for {set jj 0} {$jj < 3} {incr jj} {
            if { [lindex $max_values $jj] < [lindex $pp $jj] } {
                set max_values [lreplace $max_values $jj $jj [lindex $pp $jj]]
            }
        }
    }
    set fo [open $matrixfile w]
    # for every direction
    for {set pp 0} {$pp < 3} {incr pp} {
        set line ""
        set last_cb 0
        set max_value [lindex $max_values $pp]
        for {set tt 0} {$tt < [llength $transitions]} {incr tt} {
            set cb_trans [lindex [lindex $transitions $tt] 0]
            set vb_trans [lindex [lindex $transitions $tt] 1]
            # new line?
            if { [lindex $cb_trans 0] != $last_cb } {
                set line "${line} \n"
                set last_cb [lindex $cb_trans 0]
            }
            set value [lindex [lindex $transition_values $tt] $pp]
            set value [expr {$value / $max_value * 100}]
            set value [format "%3.0f" $value]
            set line "${line} ${value} "
        }
        puts $fo "# ----------- ${pp}, max 100 = ${max_value} -------------------------"
        puts $fo $line
        puts $fo "\n"
    }

}


# --------------------------------------------------------------------
# store transition matrix element (sum over degeneracies)
# --------------------------------------------------------------------
proc matrix_elements_write_degenerate_to_file {matrix_elements transitions filename matrixfile} {
    set fo      [open $filename w]
    set ntrans  [llength $transitions]
    set numk    [$matrix_elements get_num_k_values]
    set transk0 [list]
    for {set kk 0} {$kk < $numk} {incr kk} {
        set kval [[[$matrix_elements get_domain] get_point $kk] get_coord_abs]
        set outstr "${kval} "
        for {set tt 0} {$tt < $ntrans} {incr tt} {
            set cb_trans [lindex [lindex $transitions $tt] 0]
            set vb_trans [lindex [lindex $transitions $tt] 1]
            set trans_store [list]
            for {set dr 0} {$dr < 3} {incr dr} {
                set val 0
                foreach cc $cb_trans {
                    foreach vv $vb_trans {
                        set vadd [eval_dref [$matrix_elements get_abs_square $cc $vv $dr $kk]]
                        set val [expr {$val + $vadd}]
                    }
                }
                lappend trans_store $val
                set outstr "${outstr} \t$val"
            }
            if {$kk == 0} {
                lappend transk0 $trans_store
            }
        }
        puts $fo $outstr
    }
    close $fo
    puts "-> wrote $ntrans transitions to file $filename with the following transitions"
    for {set tt 0} {$tt < $ntrans} {incr tt} {
        set cstr ""
        set vstr ""
        set cb_trans [lindex [lindex $transitions $tt] 0]
        set vb_trans [lindex [lindex $transitions $tt] 1]
        foreach cc $cb_trans {
            if { $cstr == "" } {
                set cstr "cb $cc"
            } else {
                set cstr "$cstr, $cc"
            }
        }
        set vstr ""
        foreach vv $vb_trans {
            if { $vstr == "" } {
                set vstr "vb $vv"
            } else {
                set vstr "$vstr, $vv"
            }
        }
        puts "     $cstr -> $vstr"
    }
    write_transition_matrix_file $transk0 $transitions $matrixfile

}

source "../../../input/tdkp_helpers.tcl"

set gridfile     well384.asc.gz
set numsol       16
set num_k_values 20

set band_degeneracy_threshold 1.0e-3

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
set geometry [$parser read_geometry [AsciiGridReader -args $gridfile]]
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


set v110 [Vector3D -argc 1.0 1.0 0.0]
$problem set_axes $v110 $ez
$problem set_energy_guess 0 0.0
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

# build nonsymmetric matrices (momentum operator requires that!)
$config set "assembly_build_nonsymmetric_matrices" 1.0
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


