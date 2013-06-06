# ------------------------------------------------------------
#
# This file is part of tdkp, a simulation tool for nanostrutctures
# of optoelectronics developed at ETH Zurich
#
# (C) 2005-2009 Ratko G. Veprek, ETH Zurich, veprek@iis.ee.ethz.ch
#
# 1) As of 18.6.2009 this code is property of ETH Zurich and must not be
# transferred, modified or used by third parties without appropriate
# licenses issued by authorized agents of ETH Zurich.
#
# 2) Violation of this will result in judicial action according to civil
# and penal law.
#
# 3) Any claim of authorship other than by the author himself is
# strictly forbidden.
#
# 4) The source code must retain the copyright notice, this list
# of conditions and the following disclaimer.
#
# THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ------------------------------------------------------------

#########################################################
#                                                       #
# package matrix_elements.tcl                           #
#                                                       #
# calculate momentum matrix elements                    #
#                                                       #
#########################################################


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
    global LOG_INFO
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
    logger_emit "TCLmatrix_elements: -> wrote $ntrans transitions to file $filename with the following transitions" $LOG_INFO
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

######################################################
# calculation                                        #
######################################################
if { ![info exists ONLY_PROC_DEFINITIONS] } {


    ######################################################
    # create momentum operator and matrix elements class #
    ######################################################
    if { $dimension == 0 } {

        if { $options(-model) == "kp4x4" } {
            set momentum_operator [MomentumOperatorBulk4x4 -args $material]
        } elseif { $options(-model) == "kp6x6" } {
            set momentum_operator [MomentumOperatorBulk6x6 -args $material]
        } elseif { $options(-model) == "kp8x8" } {
            set kpmatrix   [KPMatrix8x8EndersForeman -args $material]
            $kpmatrix calculate
            set momentum_operator [MomentumOperatorBulk8x8 -args $kpmatrix]
        } elseif { $options(-model) == "kp6x6WZ" } {
            set momentum_operator [MomentumOperatorBulk6x6WZ -args $material]
        } elseif { $options(-model) == "kp8x8WZ" } {
            set kpmatrix   [KPMatrix8x8Wurtzite -args $material]
            $kpmatrix calculate
            set momentum_operator [MomentumOperatorBulk8x8 -args $kpmatrix]
        } elseif { $options(-model) == "effmass" } {
            set momentum_operator [MomentumOperatorBulkEffectiveMass -args $effmass_optical_matrix_param]
        } else {
            show_bands_tcl_options "invalid kp model requested: $options(-model)"
            exit
        }

    } else {

        if { $options(-model) == "kp4x4" } {
            set momentum_operator [MomentumOperator4x4 -args $geometry $material_database]
        } elseif { $options(-model) == "kp6x6" } {
            set momentum_operator [MomentumOperator6x6 -args $geometry $material_database]
        } elseif { $options(-model) == "kp8x8" } {
            set momentum_operator [MomentumOperator8x8 -args $geometry $material_database]
        } elseif { $options(-model) == "kp6x6WZ" } {
            set momentum_operator [MomentumOperator6x6WZ -args $geometry $material_database]
        } elseif { $options(-model) == "kp8x8WZ" } {
            set momentum_operator [MomentumOperator8x8WZ -args $geometry $material_database]
        } elseif { $options(-model) == "effmass" } {
            set momentum_operator [MomentumOperatorEffectiveMass -args $geometry $material_database]
        } else {
            show_bands_tcl_options "invalid kp model requested: $options(-model)"
            exit
        }

    }
    set matrix_elements [MatrixElements -args $momentum_operator]

    ##################################
    # load matrix elements from file #
    ##################################
    if { [info exists options(-load_matrix_elements)] } {

        logger_emit "TCLmatrix_elements: loading matrix elements from binary file"
        $matrix_elements read_binary "$result_directory/$binary_matrix_elements_file"

    #####################################################################
    # calculate matrix elements if requested of if optics are requested #
    #####################################################################
    } elseif { [info exists options(-matrix_elements)] } {

        if { $matrix_elements_bound_only == 1 } {
            $matrix_elements calculate $num_cb_bound $num_vb_bound $cb_bands $vb_bands
        } else {
            $matrix_elements calculate $cb_bands $vb_bands
        }

        ######################################
        # store matrix elements if requested #
        ######################################
        if { $binary_write_matrix_elements == 1 } {
			$matrix_elements write_binary "$result_directory/$binary_matrix_elements_file"            
        }


    }


	###################################################
	# plot matrix elements                            #
	###################################################
	if { [plot_requested "matrix_elements"] } {
		plot_to_file "matrix_elements" $matrix_elements "$result_directory/$plot_matrix_elements_file"
	}

    ###################################################
    # plot transitions of degenerate bands on request #
    ###################################################
    if { [plot_requested "transitions"] && [info exists matrix_elements] } {
        set transitions [matrix_elements_find_degenerate_transitions $matrix_elements $cb_bands $vb_bands]
        matrix_elements_write_degenerate_to_file $matrix_elements $transitions "$result_directory/$plot_transitions_file" "$result_directory/$plot_transition_strength_file"
    }
}
