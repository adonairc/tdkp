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
# package bandedges.tcl                                 #
#                                                       #
# calculate bandedges                                   #
#                                                       #
#########################################################


########################
# create problem class #
########################
if { ![info exists vb_problem] } {
    source "$sourcepath/model_init.tcl"
}

#########################################
# pure effective mass                   #
#########################################
source "$sourcepath/bandedges_calculate.tcl"

#############################################
# calculate max bandedge from problem class #
#############################################
set node_bandedges [StdNodeDatadouble -args $geometry $bandedges]
set elem_bandedges $bandedges

# delete old data
if { [info exists unbound_cb_edge] } {
    unset unbound_cb_edge
    unset unbound_vb_edge
    unset min_cb_edge
    unset min_vb_edge
}

for {set nn 0} {$nn < [$geometry get_num_nodes]} {incr nn} {
    set cb_val [eval_dref [$node_bandedges get_node_value $nn $cb_idx]]
    set vb_val [eval_dref [$node_bandedges get_node_value $nn $vb_idx]]
    set node [$geometry get_node $nn]
    # set unbound edge
    if { [$node get_location] == $location_exterior } {
        if { ![info exists unbound_cb_edge] } {
            set unbound_cb_edge $cb_val
            set unbound_vb_edge $vb_val
        }
        set unbound_cb_edge [min $unbound_cb_edge $cb_val]
        set unbound_vb_edge [max $unbound_vb_edge $vb_val]
    }
    # set min edge
    if { ![info exists min_cb_edge] } {
        set min_cb_edge $cb_val
        set min_vb_edge $vb_val
    }
    set min_cb_edge [min $cb_val $min_cb_edge]
    set min_vb_edge [max $vb_val $min_vb_edge]
}
if { $use_user_defined_unbound_edges == 1 } {
    set unbound_cb_edge $user_defined_unbound_cb_edge
    set unbound_vb_edge $user_defined_unbound_vb_edge
    logger_emit [format "TCLbandedges: using user supplied unbound edges"]
}
logger_emit [format "TCLbandedges: lowest cb edge = %g, highest vb edge = %g" $min_cb_edge $min_vb_edge]
logger_emit [format "TCLbandedges: unbound cb edge = %g, unbound vb edge = %g" $unbound_cb_edge $unbound_vb_edge]

###########################################
# write bandedges to output files         #
# (binary, ascii, etc.)                   #
###########################################
if { [plot_requested "bandedges"] } {
	plot_to_file "bandedges" $bandedges "$result_directory/$plot_bandedges_file"
}

###########################################
# calculate delta_strain                  #
###########################################
if { [info exists options(-delta_strain)] } {
    ###########################################
    # remove strain                           #
    ###########################################
    if { [info exists em_problem] } {
        $em_problem remove_strain_field
        $em_problem remove_potential_energy_field
    }
    if { [info exists vb_problem] } {
        $vb_problem remove_strain_field
        $vb_problem remove_potential_energy_field
    }
    source "$sourcepath/bandedges_calculate.tcl"
    set unstr_edges $bandedges

    ##########################################
    # calculate delta bandedges              #
    ##########################################
    set delta_strain [StdElementDatadouble -args 2 [$geometry get_num_elements]]
    $delta_strain set_identifier 0 "delta_cb_edge"
    $delta_strain set_identifier 1 "delta_vb_edge"
    for {set ii 0} {$ii < [$geometry get_num_elements]} {incr ii} {
        if { [[$geometry get_element $ii] enabled] } {
            set cb_val_s [eval_dref [$elem_bandedges get_element_value $ii $cb_idx]]
            set vb_val_s [eval_dref [$elem_bandedges get_element_value $ii $vb_idx]]
            set cb_val   [eval_dref [$unstr_edges get_element_value $ii $cb_idx]]
            set vb_val   [eval_dref [$unstr_edges get_element_value $ii $vb_idx]]
            $delta_strain set_element_value $ii 0 [expr {$cb_val_s - $cb_val}]
            $delta_strain set_element_value $ii 1 [expr {$vb_val_s - $vb_val}]
        }
    }
	#######################################################
	# write them into a binary file (which is then used   #
	# in aqua                                             #
	# this guy here is stored element wise, thats why     #
	# we don't use the plot_to_file procedure             #
	#######################################################
	set writer [TCLDataIOBinary -args [BinaryDataIO]]
	$writer write $delta_strain "$result_directory/$binary_delta_bandedges_file"

    #######################################################
    # check material database for Gas materials           #
    #######################################################
    set disabled_materials [list]
    for {set ii 0} {$ii < [$material_database get_num_materials]} {incr ii} {
        set matname [string2char [$material_database get_material_name $ii]]
        if { $matname == "Contact" || $matname == "Gas" } {
            logger_emit "TCLbandedges: ignoring material $matname for calculating delta bandedge" $LOG_WARN
            lappend disabled_materials $ii
        }
    }

    #######################################################
    # average to nodal data, but exclude contacts and air #
    #######################################################
    if { [llength $disabled_materials] > 0 } {
        set node_delta [StdNodeDatadouble -args 2 [$geometry get_num_nodes]]
        $node_delta set_identifier 0 "delta_cb_edge"
        $node_delta set_identifier 1 "delta_vb_edge"
        set contrib_counter [list]
        for {set ii 0} {$ii < [$geometry get_num_nodes]} {incr ii} {
            lappend contrib_counter 0
        }
        set skipped 0
        for {set ii 0} {$ii < [$geometry get_num_elements]} {incr ii} {
            set elem [$geometry get_element $ii]
            if { [$elem enabled] && ![exists_in_list [[$elem get_material] get_id] $disabled_materials] } {
                set cb_val [eval_dref [$delta_strain get_element_value $ii 0]]
                set vb_val [eval_dref [$delta_strain get_element_value $ii 1]]
                for {set nn 0} {$nn < [$elem get_num_nodes]} {incr nn} {
                    set node [$elem get_node $nn]
                    set nidx [$node get_index_global]
                    $node_delta set_node_value $nidx 0 [expr {[eval_dref [$node_delta get_node_value $nidx 0]] + $cb_val}]
                    $node_delta set_node_value $nidx 1 [expr {[eval_dref [$node_delta get_node_value $nidx 1]] + $vb_val}]
                    lset contrib_counter $nidx [expr {[lindex $contrib_counter $nidx] + 1}]
                }
            } else {
                if { [$elem enabled] } {
                    set skipped [expr {$skipped + 1}]
                }
            }
        }
        logger_emit "TCLbandedges: ignored $skipped noncontact-elements for delta bandedge calculation" $LOG_INFO_DEVEL2
        for {set ii 0} {$ii < [$geometry get_num_nodes]} {incr ii} {
            set div [lindex $contrib_counter $ii]
            if { $div > 0 } {
                $node_delta set_node_value $ii 0 [expr {[eval_dref [$node_delta get_node_value $ii 0]] / $div}]
                $node_delta set_node_value $ii 1 [expr {[eval_dref [$node_delta get_node_value $ii 1]] / $div}]
            } else {
                $node_delta set_node_value $ii 0 0
                $node_delta set_node_value $ii 1 0
            }
        }
    } else {
        # old code, did not work with Gas surrounding (as surrounding gas has been included in the averaging, but band edge within gas is bullshit anyway ...)
        set node_delta [StdNodeDatadouble -args $geometry $delta_strain]
    }

	if { [plot_requested "delta_bandedges" ] } {
		plot_to_file "delta_bandedges" $node_delta "$result_directory/$plot_delta_bandedges_file"
	}

}


