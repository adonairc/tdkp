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
# package piezo.tcl                                     #
#                                                       #
# calculate surface and volume charges stemming from    #
# spontaneous and piezo-electric effect                 #
#                                                       #
#########################################################

proc init_polarization_calculator { pc geometry material_database } {
    upvar $pc polarization_calculator
    global options
    $geometry prepare_boundaries
    if { $options(-model) == "kp6x6WZ" || $options(-model) == "kp8x8WZ" } {
        set polarization_calculator [PolarizationChargeWurtzite -args $geometry $material_database]
    } else {
        set polarization_calculator [PolarizationChargeZincBlende -args $geometry $material_database]
    }
}

if { ![info exists options(-load_piezo)] && ([info exists options(-nocache)] || ![file exists "$result_directory/$binary_pol_surf_charge_file"] || ![file exists "$result_directory/$binary_pol_vol_charge_file"]) } {

    ###########################################
    # prepare geometry for boundary integrals #
    # and init charge calculator              #
    ###########################################
    init_polarization_calculator polarization_calculator $geometry $material_database

    ###############################################
    # check for correct rotation matrix           #
    # if it does not exist, source model_init.tcl #
    # which create a tdkp problem class to pass   #
    # the right rotation matrix                   #
    ###############################################
    if { ![info exists rotation_matrix] } {
        source "$sourcepath/model_init.tcl"
        if { ![info exists rotation_matrix] } {
            logger_emit "TCLpiezo: i source the model_init.tcl script and expected now an existing rotation matrix. but i didn't get it. can you fix me?" $LOG_ERROR
            quit
        }
    }
    $polarization_calculator set_rotation $rotation_matrix

    ###############################################
    # are strains loaded?                         #
    ###############################################
    if { [info exists strain] } {
        $polarization_calculator set_strains $strain
    }

    ##################################
    # calculate, get values and plot #
    ##################################
    $polarization_calculator calculate
    set polarization_surface_charge [$polarization_calculator get_surface_charge]
    set polarization_volume_charge  [$polarization_calculator get_volume_charge]
    $parser write_binary $polarization_surface_charge "$result_directory/$binary_pol_surf_charge_file"
    $parser write_binary $polarization_volume_charge "$result_directory/$binary_pol_vol_charge_file"

} else {

    #################################
    # quit if no charge files exist #
    #################################
    if { [info exists options(-load_piezo)] && (![file exists "$result_directory/$binary_pol_surf_charge_file"] || ![file exists "$result_directory/$binary_pol_vol_charge_file"]) } {
        logger_emit "TCLpiezo: -load_piezo requested, but piezo charge data files are not available!" $LOG_ERROR
        quit
    }

    ############################
    # simpliy load the charges #
    ############################
    set polarization_surface_charge [StdNodeDatadouble]
    $parser read_binary $polarization_surface_charge "$result_directory/$binary_pol_surf_charge_file"
    set polarization_volume_charge  [StdElementDatadouble]
    $parser read_binary $polarization_volume_charge "$result_directory/$binary_pol_vol_charge_file"
}

################
# plot charges #
################
if { [plot_requested "piezo"] } {

    #############################################
    # make sure we have such a guy for plotting #
    #############################################
    if { ![info exists polarization_calculator] } {
        init_polarization_calculator polarization_calculator $geometry $material_database
    }

    #########################################
    # map surface charges to volume charges #
    #########################################
    set surf_vol_charge [StdElementDatadouble]
    $polarization_calculator map_surface_charge_to_volume_charge $polarization_surface_charge $surf_vol_charge

    #####################################
    # prepare new element output object #
    #####################################
    set tmp_charge [StdElementDatadouble -args 3 [$geometry get_num_elements]]
    $tmp_charge set_identifier 0 "total_charge_density"
    $tmp_charge set_identifier 1 "volume_charge_density"
    $tmp_charge set_identifier 2 "surface_charge_density"
    for {set ii 0} {$ii < [$geometry get_num_elements]} {incr ii} {
        set vc [eval_dref [$polarization_volume_charge get_element_value $ii 0]]
        set sc [eval_dref [$surf_vol_charge get_element_value $ii 0]]
        $tmp_charge set_element_value $ii 0 [expr {$vc + $sc}]
        $tmp_charge set_element_value $ii 1 $vc
        $tmp_charge set_element_value $ii 2 $sc
    }
	plot_to_file "piezo" $tmp_charge "$result_directory/$plot_all_charge_file"
}

#########################################################
# solve poisson's equation                              #
#########################################################
if { [info exists options(-nocache)] || ![file exists "$result_directory/$binary_element_potential"] } {

    # set poisson boundary conditions
    if { [info exists bc_poisson] } {
        # user defined
        logger_emit "TCLpiezo: PoissonEquation uses user defined boundary conditions" $LOG_INFO_DEVEL2
        $geometry set_boundary_conditions $bc_poisson
    } else {
        # standard, we take dirichlet
        logger_emit "TCLpiezo: PoissonEquation uses Dirichlet boundary conditions everywhere" $LOG_INFO_DEVEL2
        $geometry set_boundary_conditions [BCDirichlet -args $geometry]
    }
    set poisson [PoissonEquation -args $geometry $material_database]
    if { [info exists polarization_surface_charge] } {
        $poisson set_surface_charge_density $polarization_surface_charge
    }
    if { [info exists polarization_volume_charge] } {
        $poisson set_element_charge_density $polarization_volume_charge
    }
    $poisson solve 1
    set potential [$poisson get_solution]

    if { [plot_requested "potential"] } {
		plot_to_file "potential" $potential "$result_directory/${plot_potential_file}"    }

    #############################################
    # build potential field for bandstructure   #
    #############################################
    set potential_field [PotentialEnergyField -args [$geometry get_num_nodes]]
    for {set ii 0} {$ii < [$potential_field get_length]} {incr ii} {
        set pt [eval_dref [$potential get_node_value $ii 0]]
        set pt [expr {$pt * $constants_ec * (-1.0)}]
        $potential_field set_node_value $ii 0 $pt
    }
    $parser write_binary $potential_field "$result_directory/$binary_element_potential"
} else {
    set potential_field [PotentialEnergyField -args [$geometry get_num_nodes]]
    $parser read_binary $potential_field "$result_directory/$binary_element_potential"
    logger_emit "TCLpiezo: loaded potential from file $result_directory/$binary_element_potential"
    if { [plot_requested "potential"] } {
		logger_emit "TCLpiezo: you requested a plot of the electrostatic potential which is not possible as it is loaded from the binary element potential" $LOG_WARN
	}
}

if { [info exists vb_problem] } {
    unset vb_problem
}
if { [info exists em_problem] } {
    unset em_problem
}
if { [info exists cb_problem] } {
    unset cb_problem
}
