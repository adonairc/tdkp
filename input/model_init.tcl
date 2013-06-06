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
# package model_init.tcl                                #
#                                                       #
# create and initialize problem class objects           #
#                                                       #
#########################################################

######################
# delete old objects #
######################
set cb_effmass 0
set vb_effmass 0

if { [info exists vb_problem] } {
    $vb_problem -delete
}
if { [info exists cb_problem] } {
    $cb_problem -delete
}
if { [info exists em_problem] } {
    $em_problem -delete
}

############################################
# reset rotation matrix to identity matrix #
############################################
set rotation_matrix $identity_matrix

if { $dimension < 3 && $dimension > 0 } {

    if { $options(-model) == "kp4x4" } {
        set vb_problem  [KP4x41D2D -args $geometry $material_database]
        set cb_effmass 1
    } elseif { $options(-model) == "kp6x6" } {
        set vb_problem  [KP6x61D2D -args $geometry $material_database]
        set cb_effmass 1
    } elseif { $options(-model) == "kp8x8" } {
        set vb_problem  [KP8x81D2D -args $geometry $material_database]
    } elseif { $options(-model) == "kp6x6WZ" } {
        set vb_problem  [KP6x61D2DWZ -args $geometry $material_database]
        set cb_effmass 1
    } elseif { $options(-model) == "kp8x8WZ" } {
        set vb_problem  [KP8x81D2DWZ -args $geometry $material_database]
    } elseif { $options(-model) == "effmass" } {
        set cb_effmass 1
        set vb_effmass 1
    } else {
        show_bands_tcl_options "invalid kp model requested: $options(-model)"
        exit
    }

    #############################################
    # create effective mass object if necessary #
    #############################################
    if { $cb_effmass == 1 || $vb_effmass == 1 } {
        set em_problem [EffectiveMass -args $geometry $material_database]
        if { [info exists strain] } {
            if {$options(-model) == "kp6x6WZ"} {
                ###########################################
                # set hydrostatic potential field names   #
                ###########################################
                # initially set all to aa type
                for {set ii 0} {$ii < 3} {incr ii} {
                    $em_problem set_hydro_strain_potential_field_name_cb $ii "strain_potential_ac_aa"
                    $em_problem set_hydro_strain_potential_field_name_vb $ii "strain_potential_av_aa"
                }
                $em_problem set_hydro_strain_potential_field_name_cb $wurzite_strain_axis_3 "strain_potential_ac_cc"
                $em_problem set_hydro_strain_potential_field_name_vb $wurzite_strain_axis_3 "strain_potential_av_cc"
                logger_emit "you are using wurtzite with effective mass. please check that the wurtzite strain axes give the correct axis permutation for the considered system." $LOG_WARN
            }
            $em_problem set_field $strain
        }
        if { [info exists potential_field] } {
            $em_problem set_field $potential_field
        }
    }

    #################################
    # set axes and strain to system #
    #################################
    if { [info exists vb_problem] } {
        if {$dimension == 1} {
            $vb_problem set_axes $well_transverse_direction $well_quantized_direction
        } else {
            $vb_problem set_axes $wire_transverse_direction_z $wire_quantized_direction_x $wire_quantized_direction_y
        }

        ###############################################
        # catch rotation matrix from vb problem class #
        ###############################################
        set rotation_matrix [$vb_problem get_rotation_matrix]

        ###############################################
        # set strains and potential field             #
        ###############################################
        if { [info exists strain] } {
            $vb_problem set_field $strain
        }
        if { [info exists potential_field] } {
            $vb_problem set_field $potential_field
        }
    }

} elseif { $dimension == 3 } {

    if { $options(-model) == "kp4x4" } {
        set vb_problem  [KP4x43D -args $geometry $material_database]
        set cb_effmass 1
    } elseif { $options(-model) == "kp6x6" } {
        set vb_problem  [KP6x63D -args $geometry $material_database]
        set cb_effmass 1
    } elseif { $options(-model) == "kp8x8" } {
        set vb_problem  [KP8x83D -args $geometry $material_database]
    } elseif { $options(-model) == "kp6x6WZ" } {
        set vb_problem  [KP6x63DWZ -args $geometry $material_database]
        set cb_effmass 1
    } elseif { $options(-model) == "kp8x8WZ" } {
        set vb_problem  [KP8x83DWZ -args $geometry $material_database]
    } elseif { $options(-model) == "effmass" } {
        set cb_effmass 1
        set vb_effmass 1
    } else {
        show_bands_tcl_options "invalid kp model requested: $options(-model)"
        exit
    }

    ######################################
    # create effmass object if necessary #
    ######################################
    if { $cb_effmass == 1 || $vb_effmass == 1 } {
        set em_problem [EffectiveMass -args $geometry $material_database]
        if { [info exists strain] } {
            if {$options(-model) == "kp6x6WZ"} {
                ###########################################
                # set hydrostatic potential field names   #
                ###########################################
                # initially set all to aa type
                for {set ii 0} {$ii < 3} {incr ii} {
                    $em_problem set_hydro_strain_potential_field_name_cb $ii "strain_potential_ac_aa"
                    $em_problem set_hydro_strain_potential_field_name_vb $ii "strain_potential_av_aa"
                }
                $em_problem set_hydro_strain_potential_field_name_cb $wurzite_strain_axis_3 "strain_potential_ac_cc"
                $em_problem set_hydro_strain_potential_field_name_vb $wurzite_strain_axis_3 "strain_potential_av_cc"
                logger_emit "you are using wurtzite with effective mass. please check that the wurtzite strain axes give the correct axis permutation for the considered system." $LOG_WARN
            }
            $em_problem set_field $strain
        }
        if { [info exists potential_field] } {
            $em_problem set_field $potential_field
        }
    }


    ######################
    # set axes to system #
    ######################
    if { [info exists vb_problem] } {
        $vb_problem set_axes $dot_direction_x $dot_direction_y $dot_direction_z
		set rotation_matrix [$vb_problem get_rotation_matrix]
    }
    #########################################
    # rotate strain for nitride problems    #
    #########################################
    if {[info exists strain] && ($options(-model) == "kp6x6WZ" || $options(-model) == "kp8x8WZ")} {
        if {[info exists strain_is_rotated]} {
            logger_emit "warning, it looks like the strain has been rotated already and is now rotated again!" $LOG_WARN
            exec sleep 5
        }
        set strain_is_rotated 1
        ##############################################################
        # we calculate the strain in a system with permuted axes e'' #
        # e'' = P e P^{t}                                            #
        # and the system of the envelopes is given in terms of       #
        # e' = R e R^{t}                                             #
        # therefore the right transformation is given by             #
        # e' = R P^{t} e P R^{t}                                     #
        ##############################################################
        set perm [RotationMatrix -args 3 3]
        $perm set $wurzite_strain_axis_1 0 1
        $perm set $wurzite_strain_axis_2 1 1
        $perm set $wurzite_strain_axis_3 2 1
        set perm [$perm get_transpose]
        set rot [$vb_problem get_rotation_matrix]
        set perm [$rot mm $perm]
        $strain rotate $perm
    }
    if { [info exists strain] } {
        $vb_problem set_field $strain
    }
    if { [info exists potential_field] } {
        $vb_problem set_field $potential_field
    }

}

