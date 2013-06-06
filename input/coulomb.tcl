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
# package scl.tcl                                       #
#                                                       #
# calculate simple luminescence                         #
#                                                       #
#########################################################


# ------------------------------------------------------------
# prepare coulomb object with states at k = 0!
# ------------------------------------------------------------
proc prepare_coulomb_object { variable geometry material_database cb_bands vb_bands } {

    upvar $variable coulomb
    global coulomb_q_min
    global coulomb_q_max
    global coulomb_q_num
    global coulomb_progression

    # create object
    set coulomb [CoulombMatrixElementQuantized -args $geometry $material_database [$vb_bands get_basis_size]]
    # register cb wavefunctions
    for {set ii 0} {$ii < [$cb_bands get_number_of_bands]} {incr ii} {
        $coulomb add_wavefunction [$cb_bands get_eigensolution_ptr 0 $ii] [format "cb%d" $ii]
    }
    # register vb wavefunctions
    for {set ii 0} {$ii < [$vb_bands get_number_of_bands]} {incr ii} {
        $coulomb add_wavefunction [$vb_bands get_eigensolution_ptr 0 $ii] [format "vb%d" $ii]
    }
    $coulomb set_q_range $coulomb_q_min $coulomb_q_max $coulomb_q_num $coulomb_progression
}

# ------------------------------------------------------------
# prepare bulk coulomb object with states at k = 0!
# ------------------------------------------------------------
proc prepare_coulomb_object_bulk { variable cb_bands vb_bands } {

    upvar $variable coulomb
    global coulomb_q_min
    global coulomb_q_max
    global coulomb_q_num
    global coulomb_progression

    # create object
    set coulomb [CoulombMatrixElementBulk]
    # register cb wavefunctions
    for {set ii 0} {$ii < [$cb_bands get_number_of_bands]} {incr ii} {
        $coulomb add_wavefunction [$cb_bands get_eigensolution_ptr 0 $ii] [format "cb%d" $ii]
    }
    # register vb wavefunctions
    for {set ii 0} {$ii < [$vb_bands get_number_of_bands]} {incr ii} {
        $coulomb add_wavefunction [$vb_bands get_eigensolution_ptr 0 $ii] [format "vb%d" $ii]
    }

    $coulomb set_q_range $coulomb_q_min $coulomb_q_max $coulomb_q_num $coulomb_progression
}

# ------------------------------------------------------------
# calculate diagonal coulomb matrix elements
# ------------------------------------------------------------
proc calculate_coulomb_diag { variable cb_bands vb_bands } {

    upvar $variable coulomb

    # register diagonal objects
    set num_cb  [$cb_bands get_number_of_bands]
    set num_vb  [$vb_bands get_number_of_bands]
    set num_tot [expr {$num_cb + $num_vb}]

    for {set ii 0} {$ii < $num_tot} {incr ii} {
        for {set jj 0} {$jj < $num_tot} {incr jj} {
            $coulomb request_matrix_element $ii $jj $jj $ii
        }
    }

    $coulomb calculate

}

# ------------------------------------------------------------
# calculate ALL coulomb matrix elements
# ------------------------------------------------------------
proc calculate_coulomb_all { variable cb_bands vb_bands } {

    upvar $variable coulomb

    # register diagonal objects
    set num_cb  [$cb_bands get_number_of_bands]
    set num_vb  [$vb_bands get_number_of_bands]
    set num_tot [expr {$num_cb + $num_vb}]

    for {set ii 0} {$ii < $num_tot} {incr ii} {
        for {set jj 0} {$jj < $num_tot} {incr jj} {
            for {set kk 0} {$kk < $num_tot} {incr kk} {
                for {set ll 0} {$ll < $num_tot} {incr ll} {
                    $coulomb request_matrix_element $ii $jj $kk $ll
                }
            }
        }
    }

    $coulomb calculate

}

# ------------------------------------------------------------
# calculate SELECTED coulomb matrix elements
# ------------------------------------------------------------
proc calculate_coulomb_selected { variable cb_bands vb_bands my_selection} {

    upvar $variable coulomb
    set sp [split $my_selection ":"]
    # register requested object
    foreach request $sp {
        set cc [split $request "-"]
        if { [llength $cc] == 4 } {
            $coulomb request_matrix_element [lindex $cc 0] [lindex $cc 1] [lindex $cc 2] [lindex $cc 3]
        } else {
            puts " ** error **: can not parse coulomb integral request $request"
        }
    }
    $coulomb calculate

}

# -----------------------------------------------------------
# executed code
# -----------------------------------------------------------

# -----------------------------------------------------------
# load coulomb matrix elements if they exist
# -----------------------------------------------------------
if { [info exists options(-load_coulomb)] && [file exists "$result_directory/$coulomb_binary_file"] } {
    set coulomb_data [CoulombMatrixElementData -args "$result_directory/$coulomb_binary_file"]
} else {
    # -----------------------------------------------------------
    # bulk is different
    # -----------------------------------------------------------
    if { [info exists options(-bulk)] } {
        prepare_coulomb_object_bulk coulomb $cb_bands $vb_bands
    } else {
        $geometry set_boundary_conditions [BCIncludeAll -args $geometry]
        # create object and add wavefunctions
        prepare_coulomb_object coulomb $geometry $material_database $cb_bands $vb_bands
    }
    # -----------------------------------------------------------
    # fresh calculate coulomb matrix elements
    # -----------------------------------------------------------
    set sp [split $options(-coulomb) ":"]
    if { $options(-coulomb) == "all" } {
        # create object and add wavefunctions
        calculate_coulomb_all  coulomb $cb_bands $vb_bands
    } elseif { [llength $sp] > 1 } {
        calculate_coulomb_selected coulomb $cb_bands $vb_bands $options(-coulomb)
    } else {
        calculate_coulomb_diag coulomb $cb_bands $vb_bands
    }
    set coulomb_data [$coulomb get_data_object]

    # -----------------------------------------------------------
    # write data to cache
    # -----------------------------------------------------------
	$coulomb_data write_binary "$result_directory/$coulomb_binary_file"
}
if { [plot_requested "coulomb"] } {
	plot_to_file "coulomb" $coulomb_data "$result_directory/$plot_coulomb_file"
}


