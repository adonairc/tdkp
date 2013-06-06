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
# package strain.tcl                                    #
#                                                       #
# calculate intrinsic strain in device                  #
#                                                       #
#########################################################


##############################
# fresh evaluation of strain #
##############################
if { [info exists options(-nocache)] || ![file exists "$result_directory/$binary_strain_file"] } {

    #######################################
    # set appropriate boundary conditions #
    #######################################
    if { [info exists bc_strain] } {
        $geometry set_boundary_conditions $bc_strain
        $geometry prepare
    } else {
        if { $dimension ==  1 } {
            $geometry set_boundary_conditions [BCWellStrainFloatRight -args $geometry]
        } else {
            # free floating?
            $geometry set_boundary_conditions [BCIncludeAll -args $geometry]
        }
    }

    ###################################
    # different treatment for wurzite #
    ###################################
    if { $options(-model) == "kp6x6WZ" || $options(-model) == "kp8x8WZ" } {
        set problem [IntrinsicStrainWurzite -args $geometry $material_database]
        # set axis permutation
        $problem set_axes $wurzite_strain_axis_1 $wurzite_strain_axis_2 $wurzite_strain_axis_3
        # set reference lattice constants
        if { $options(-strain) == 1 } {
            set ref_material [$material_database get_material [find_reference_lattice_constant_material $geometry $material_database "lattice_constant_a"]]
        } else {
            set ref_material [$material_database get_material $options(-strain)]
        }
        $problem set_reference_lattice_constant [$ref_material get "lattice_constant_a"] [$ref_material get "lattice_constant_c"]
    } else {
        set problem [IntrinsicStrain -args $geometry $material_database]
        # user supplied reference const (lattice const of 1 is very unusual ... so i hope nobody will use that ever)
        if { $options(-strain) == 1 } {
            set ref_material [$material_database get_material [find_reference_lattice_constant_material $geometry $material_database "lattice_constant"]]
        } else {
            set ref_material [$material_database get_material $options(-strain)]
        }
        set ref_lattice_const [$ref_material get "lattice_constant"]
        logger_emit "TCLStrain: using $ref_lattice_const as reference lattice constant"
        $problem set_reference_lattice_constant $ref_lattice_const
    }
    #############################################
    # general case: assume biaxial strain in QW #
    #############################################
    if { $dimension == 1 } {
        # define x-axis (quantized direction) to be unstrained (bi-axial strain)
        $problem set_strained_axes 0 0
    }
    #############################################
    # set user defined unstrained axes          #
    #############################################
    if { $axis_x_strained == 0 } {
        $problem set_strained_axes 0 0
    }
    if { $axis_y_strained == 0 } {
        $problem set_strained_axes 1 0
    }
    if { $axis_z_strained == 0 } {
        $problem set_strained_axes 2 0
    }

    $problem solve 1
    set strain [$problem get_strain_field]
	# todo: change this into the general output wrapper!
    $strain write_binary "$result_directory/$binary_strain_file"

    ############################################
    # plot displacement                        #
    ############################################	
	if { [plot_requested "displacement"] } {
        set displacement [$problem get_solution]
		plot_to_file "displacement" $displacement "$result_directory/$plot_displacement_file"
    }
    $geometry set_boundary_conditions [BCDirichlet -args $geometry]

} else {
    ############################################
    # load from binary file                    #
    ############################################
    logger_emit "TCLStrain: loading strain data from $result_directory/$binary_strain_file" $LOG_INFO_DEVEL2
    set strain [StrainField -args "$result_directory/$binary_strain_file"]
    ##################################################
	# emit warning if displacement plot is requested #
    ##################################################
	if { [plot_requested "displacement"] } {
     	logger_emit "TCLStrain: can't plot displacement for cached strain data file" $LOG_WARN
	}
}

############################################
# plot strain on request
############################################
if { [plot_requested "strain"] } {
	plot_to_file "strain" $strain "$result_directory/$plot_strain_file"
}

#########################################
# compute hydrostatic strain on request #
#########################################
if { [plot_requested "hydrostatic_strain"] } {
    set hydro [StdElementDatadouble -args 2 [$strain get_length]]
    $hydro set_identifier 0 "hydrostatic_strain"
    $hydro set_identifier 1 "deviatoric_strain"
    for {set ii 0} {$ii < [$strain get_length]} {incr ii} {
        set tensor [$strain get $ii]
        set exx [eval_dref [$tensor get $iexx]]
        set eyy [eval_dref [$tensor get $ieyy]]
        set ezz [eval_dref [$tensor get $iezz]]
        set trE [expr {$exx + $eyy + $ezz}]
        $hydro set_element_value $ii 0 $trE
        #####################
        # deviatoric strain #
        #####################
        set exy [eval_dref [$tensor get $iexy]]
        set exz [eval_dref [$tensor get $iexz]]
        set eyz [eval_dref [$tensor get $ieyz]]
        set exxr [expr {$exx - 1.0 / 3.0 * $trE}]
        set eyyr [expr {$eyy - 1.0 / 3.0 * $trE}]
        set ezzr [expr {$ezz - 1.0 / 3.0 * $trE}]
        set dev [expr {$exxr*$exxr + $eyyr * $eyyr + $ezzr * $ezzr + 2.0 * $exy * $exy + 2.0 * $exz * $exz + 2.0 * $eyz * $eyz}]
        $hydro set_element_value $ii 1 $dev
    }
	plot_to_file "hydrostatic_strain" $hydro "$result_directory/$plot_hydrostatic_strain_file"
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
