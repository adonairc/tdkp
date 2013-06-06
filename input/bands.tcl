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
# package bands.tcl                                     #
#                                                       #
# - calculate bandstructure for valence and conduction  #
#   bands                                               #
# - determine band edges                                #
# - determine bound states                              #
# - saves bandstructure to binary files                 #
# - create dispersive band structures                   #
#                                                       #
#########################################################


####################################
# set energy guess list to problem #
####################################
proc set_energy_guess_list { problem guess } {
    for {set ii 0} {$ii < [llength $guess]} {incr ii} {
        set en [lindex $guess $ii]
        $problem set_energy_guess $ii $en
    }
}

################################
# load bandstructure from file #
################################
if { [info exists options(-load_bands)] } {

    #########################################################
    # check bandstructure if it agrees with the given model #
    #########################################################
    set expect_cb 1
    set expect_vb 1
    if { $options(-model) == "kp4x4" } {
        set expect_vb 4
    } elseif { $options(-model) == "kp6x6" } {
        set expect_vb 6
    } elseif { $options(-model) == "kp8x8" } {
        set expect_cb 8
        set expect_vb 8
    } elseif { $options(-model) == "kp6x6WZ" } {
        set expect_vb 6
    } elseif { $options(-model) == "kp8x8WZ" } {
        set expect_cb 8
        set expect_vb 8
    } elseif { $options(-model) == "effmass" } {
        set expect_vb 1
        set expect_cb 1
    } else {
        show_bands_tcl_options "TCLbands: invalid kp model requested $options(-model)"
        exit
    }

    if { $num_vb_sol > 0 } {
        set vb_bands [BandstructureDomaincomplex -args "$result_directory/$vb_bands_binary_file"]
        set num_vb_bound [$vb_bands get_number_of_bands]
        set vb_basis_size [$vb_bands get_basis_size]
        if { $vb_basis_size != $expect_vb } {
            logger_emit "TCLbands: the basis size in the binary vb file is $vb_basis_size instead of $expect_vb" $LOG_ERROR
            exit
        }
        ########################################################
        # if fewer bands are desired than stored, extract them #
        ########################################################
        if { $num_vb_sol < $num_vb_bound } {
            set num_vb_bound $num_vb_sol
            set vb_bands_old $vb_bands
            set vb_bands     [$vb_bands_old extract_bands $num_vb_bound 0]
            $vb_bands_old -delete
            unset vb_bands_old
        }
        logger_emit "TCLbands: loaded $num_vb_bound valence bands"
    }
    if { $num_cb_sol > 0 } {
        set cb_bands [BandstructureDomaincomplex -args "$result_directory/$cb_bands_binary_file"]
        set num_cb_bound [$cb_bands get_number_of_bands]
        set cb_basis_size [$cb_bands get_basis_size]
        if { $cb_basis_size != $expect_cb } {
            logger_emit "TCLbands: the basis size in the binary cb file is $cb_basis_size instead of $expect_cb" $LOG_ERROR
            exit
        }
        ########################################################
        # if fewer bands are desired than stored, extract them #
        ########################################################
        if { $num_cb_sol < $num_cb_bound } {
            set num_cb_bound $num_cb_sol
            set cb_bands_old $cb_bands
            set cb_bands     [$cb_bands_old extract_bands $num_cb_bound 0]
            $cb_bands_old -delete
            unset cb_bands_old
        }
        logger_emit "TCLbands: loaded $num_cb_bound conduction bands"
    }


##############################################
# calculate bulk bandstructures if requested #
##############################################
} elseif { [info exists options(-bulk)] && [info exists options(-bands)] } {

    set cb_effmass 0
    set vb_effmass 0
    set cb_offset  -1

    ################################
    # create appropriate kp matrix #
    ################################
    if { $options(-model) == "kp4x4" } {
        set kpmatrix   [KPMatrix4x4EndersForeman]
        set cb_effmass 1
    } elseif { $options(-model) == "kp6x6" } {
        set kpmatrix   [KPMatrix6x6EndersForeman]
        set cb_effmass 1
    } elseif { $options(-model) == "kp8x8" } {
        set kpmatrix   [KPMatrix8x8EndersForeman]
        set cb_offset 6
    } elseif { $options(-model) == "kp6x6WZ" } {
        set kpmatrix   [KPMatrix6x6Wurtzite]
        set cb_effmass 1
    } elseif { $options(-model) == "kp8x8WZ" } {
        set kpmatrix   [KPMatrix8x8Wurtzite]
        set cb_offset 6
    } elseif { $options(-model) == "effmass" } {
        set cb_effmass 1
        set vb_effmass 1
    } else {
        show_bands_tcl_options "invalid kp model requested $options(-model)"
        exit
    }

    ##############################################
    # set strain potentials for hydrostatic case #
    ##############################################
    set effmass_cb_strain_shift 0.0
    set effmass_vb_strain_shift 0.0
    if { [info exists bulk_strain] } {
        set bulk_mat [$material_database get_material $options(-bulk)]
        if { $options(-model) == "kp6x6WZ" } {
            set hydro_ac($wurzite_strain_axis_1) [$material get "strain_potential_ac_aa"]
            set hydro_ac($wurzite_strain_axis_2) [$material get "strain_potential_ac_aa"]
            set hydro_ac($wurzite_strain_axis_3) [$material get "strain_potential_ac_cc"]
        } elseif { $cb_effmass == 1 } {
            set hydro_ac(0) [$material get "strain_potential_ac"]
            set hydro_ac(1) [$material get "strain_potential_ac"]
            set hydro_ac(2) [$material get "strain_potential_ac"]
        }
        if { $vb_effmass == 1 } {
            set hydro_av(0) [$material get "strain_potential_av"]
            set hydro_av(1) [$material get "strain_potential_av"]
            set hydro_av(2) [$material get "strain_potential_av"]
        }

        for {set ii 0} {$ii < 3} {incr ii} {
            set eii [$bulk_strain get $ii $ii]
            if { [info exists hydro_ac] } {
                set effmass_cb_strain_shift [expr {$effmass_cb_strain_shift + $eii * $hydro_ac($ii)}]
            }
            if { [info exists hydro_av] } {
                set effmass_vb_strain_shift [expr {$effmass_vb_strain_shift + $eii * $hydro_av($ii)}]
            }
        }
    }

    ######################
    # solve valence band #
    ######################
    set domain [DomainMaster]
    create_3D_domain_radial $domain $bulk_transverse_direction $kmin $kmax $num_k_values
    set sp [DomainNodeSingularPoint -args 1.0]
    set singular_domain [DomainMaster -args $sp]
    if { $vb_effmass == 0 } {
        $kpmatrix set_material [$material_database get_material $options(-bulk)]
        ######################
        # set strain         #
        ######################
        if {[info exists bulk_strain] && [info exists kpmatrix]} {
            logger_emit "setting strains to kp matrix" $LOG_INFO
            $kpmatrix calculate
            $kpmatrix set_strains $bulk_strain
        }
        set solver [BulkBandstructureSolver -args $kpmatrix]
        ############################
        # use user defined domain? #
        ############################
        if { [info exists user_defined_bands_domain] } {
            $solver solve $user_defined_bands_domain
        } else {
            $solver solve $bulk_transverse_direction $kmin $kmax $num_k_values
        }
        set all_bands [$solver get_bandstructure]
        if { $cb_offset > -1 } {
            set num_bands [$all_bands get_number_of_bands]
            set vb_bands [$all_bands extract_bands $cb_offset 0]
        } else {
            set vb_bands $all_bands
        }



    ########################
    # effmass valence band #
    ########################
    } else {
        set vb_bands [BandstructureDomaincomplex -args 1 1 1 $singular_domain]
        set es [EigenSolutioncomplex -args 1 1]
        $es set_node_value 0 0 [double2complex 1.0 0.0]
        set energy     [expr {[$material get "valence_band_edge"] + $effmass_vb_strain_shift}]
        $es set_energy [double2complex $energy 0.0]
        # add strain?
        $vb_bands add_eigensolution 0 0 $es
    }

    ########################
    # extract cb band data #
    ########################
    if { $cb_effmass == 0 } {
        set cb_bands [$all_bands extract_bands [expr {$num_bands - $cb_offset}] $cb_offset]
    } else {
        set cb_bands [BandstructureDomaincomplex -args 1 1 1 $singular_domain]
        set es [EigenSolutioncomplex -args 1 1]
        $es set_node_value 0 0 [double2complex 1.0 0.0]
        set energy     [expr {[$material get "conduction_band_edge"] + $effmass_cb_strain_shift}]
        $es set_energy [double2complex $energy 0.0]
        $cb_bands add_eigensolution 0 0 $es
    }

    #########################################
    # create decross object for regrouping  #
    # obtained wavefunctions such that spin #
    # degenerate bands are composed and     #
    # have same phase factor. needed to     #
    # create continuous matrix elements     #
    # which can be interpolated properly    #
    #########################################
    if { [info exists decross_bands] && $decross_bands == 1 } {
        set decross [DeCross]
    }
    #############################################################
    # decross if requested (should only be necessary for kp 8x8 #
    #############################################################
    if { [info exists decross] && [info exists cb_bands] && [$cb_bands get_basis_size] > 1 } {
        $decross resort $cb_bands
    }
    if { [info exists decross] && [info exists vb_bands] && [$vb_bands get_basis_size] > 1 } {
        $decross resort $vb_bands
    }
    if { [info exists decross] } {
        $decross -delete
        unset decross
    }

    ############################
    # remove some bands?       #
    ############################
    if { [$vb_bands get_number_of_bands] != $num_vb_sol } {
        # bands are ordered from lowest to highest
        set offset [expr {[$vb_bands get_number_of_bands] - $num_vb_sol}]
        if { $offset >= 0 } {
            set vb_bands [$vb_bands extract_bands $num_vb_sol $offset]
        }
    }
    if { [$cb_bands get_number_of_bands] != $num_cb_sol } {
        # bands are ordered from lowest to highest
        if { [$cb_bands get_number_of_bands] > $num_cb_sol } {
            set cb_bands [$cb_bands extract_bands $num_cb_sol 0]
        }
    }


    ##############################
    # store band data on request #
    ##############################
    if { $binary_write_bands == 1 } {
        if { [info exists cb_bands] } {
            $cb_bands write_binary "$result_directory/$cb_bands_binary_file"
        }
        if { [info exists vb_bands] } {
            $vb_bands write_binary "$result_directory/$vb_bands_binary_file"
        }
    }

    set num_cb_bound [$cb_bands get_number_of_bands]
    set num_vb_bound [$vb_bands get_number_of_bands]

    ##################
    # delete objects #
    ##################
    if { [info exists kpmatrix] } {
        $kpmatrix -delete
        unset kpmatrix
    }
    $domain -delete
    $singular_domain -delete
    unset domain
    unset singular_domain

##############################
# calculate quantized system #
##############################
} elseif { [info exists options(-bands)] } {

    #######################################
    # set boundary conditions             #
    #######################################
    if { [info exists bc_bands] } {
        $geometry set_boundary_conditions $bc_bands
    } else {
		# ignore gas?
		if { $disable_gas_for_bandstructure_calc == 1 } {
			$geometry set_boundary_conditions [BCIgnoreGas -args [BCDirichlet -args $geometry] $material_database]
		} else {
        	$geometry set_boundary_conditions [BCDirichlet -args $geometry]
		}
    }

    ######################################
    # apply graph reordering             #
    ######################################
    if { $graph_reordering_using_metis == 1 } {
		GraphReordering_apply_reordering $geometry
    }

    ######################################
    # set energy guess                   #
    ######################################
    if { [info exists options(-energy_guess)] } {
        set sp [split $options(-energy_guess) ":"]
        if { [llength $sp] != 2 } {
            logger_emit [format "TCLbands: invalid energy guess supplied: %s" $options(-energy_guess)] $LOG_WARN
            exit
        }
        set cb_energy_guess [lindex $sp 0]
        set vb_energy_guess [lindex $sp 1]
    } else {
        set cb_energy_guess [expr {($unbound_cb_edge - $min_cb_edge) / 4.0 + $min_cb_edge}]
        set vb_energy_guess $min_vb_edge
    }

    #######################################
    # solve kp problem for 1D/2D problems #
    #######################################
    if { $dimension < 3 } {

        ###############################
        # setup problem classes       #
        ###############################
        if { ![info exists vb_problem] } {
            source "$sourcepath/model_init.tcl"
        }

        #########################################
        # create energy guess for parallel eval #
        #########################################
        if {[info exists options(-par)]} {
            $config set "desired_eigenvalue_solver" 4
            $config set "assembly_build_nonsymmetric_matrices" 1.0
        }

        ###############################
        # solve conduction band first #
        ###############################
        if {$num_cb_sol > 0} {
            ##################
            # effective mass #
            ##################
            if {$cb_effmass == 1} {
                $em_problem set_solution_type $electrons
                $em_problem set_energy_guess  $cb_energy_guess
                ##################################
                # PML or standard?               #
                ##################################
                if { $enable_pml == 1 } {
                    schroedinger_pml_factory em_pml $options(-model) $em_problem
                    $em_pml solve $num_cb_sol
                    if { $pml_print_bands_orderered_by_tau == 1 } {
                        $em_pml display_solution_info
                    }
                    if { $pml_resort_according_to_lifetime == 1 } {
                        set cb_bands [BandstructureDomaincomplex]
                        extract_stable_bands $pml_threshold_stable_states_max_imag_value [$em_problem get_bandstructure] $cb_bands
                    } else {
                        set cb_bands [BandstructureDomaincomplex -args [$em_problem get_bandstructure]];
                    }
                } else {
                    $em_problem solve $num_cb_sol
                    set cb_bands [BandstructureDomaincomplex -args [$em_problem get_bandstructure]];
                }
                $em_problem delete_solutions
            ###################################
            # some more complex bandstructure #
            ###################################
            } else {
                $vb_problem set_solution_type  $electrons
                $vb_problem set_energy_barrier [expr {($min_cb_edge + $min_vb_edge) / 2.0}]
                $vb_problem set_energy_guess 0 $cb_energy_guess
                ##################################
                # PML or standard?               #
                ##################################
                if { $enable_pml == 1 } {
                    schroedinger_pml_factory vb_pml $options(-model) $vb_problem
                    $vb_pml solve $num_cb_sol $kmin $kmax $num_k_values
                    if { $pml_print_bands_orderered_by_tau == 1 } {
                        $vb_pml display_solution_info
                    }
                    if { $pml_resort_according_to_lifetime == 1 } {
                        set cb_bands [BandstructureDomaincomplex]
                        extract_stable_bands $pml_threshold_stable_states_max_imag_value [$vb_problem get_bandstructure 0] $cb_bands
                    } else {
                        set cb_bands [BandstructureDomaincomplex -args [$vb_problem get_bandstructure 0]];
                    }
                } else {
                    $vb_problem solve $num_cb_sol $kmin $kmax $num_k_values
                    set cb_bands [BandstructureDomaincomplex -args [$vb_problem get_bandstructure]];
                }
                $vb_problem remove_all_energy_guesses
                set cb_bands [$vb_problem get_bandstructure 0];
            }


            ######################################################
            # check bound / valid bands and extract if requested #
            ######################################################
            set num_cb_bound [find_num_bound_states $cb_bands $unbound_cb_edge $electrons]
            set num_cb_valid [find_num_valid_states $cb_bands $electrons]
            set num_to_keep  [$cb_bands get_number_of_bands]
            if { $binary_keep_only_valid == 1} {
                set num_to_keep $num_cb_valid
            }
            if { $binary_keep_only_bound == 1 } {
                set num_to_keep [min $num_cb_bound $num_to_keep]
            }
            if { $num_to_keep == 0 } {
                logger_emit "TCLbands: all cb bands are invalid or unbound\nbound: $num_cb_bound, valid: $num_cb_valid\n-> won't extract any bands. inspect your bandstructure!" $LOG_ERROR
            } elseif { [$cb_bands get_number_of_bands] > $num_to_keep } {
                set cb_bands [$cb_bands extract_bands $num_to_keep 0]
            }

        }

        ##########################
        # solve for valence band #
        ##########################
        if {$num_vb_sol > 0} {
            if {$vb_effmass == 0} {
                $vb_problem set_solution_type $holes
                $vb_problem set_energy_guess  0 $vb_energy_guess
                ##################################
                # PML or standard?               #
                ##################################
                if { $enable_pml == 1 } {
                    if { ![info exists vb_pml] } {
                        schroedinger_pml_factory vb_pml $options(-model) $vb_problem
                    }
                    $vb_pml solve $num_vb_sol $kmin $kmax $num_k_values
                    if { $pml_print_bands_orderered_by_tau == 1 } {
                        $vb_pml display_solution_info
                    }
                    if {$cb_effmass == 1 || $num_cb_sol == 0} {
                        set vb_bands_tmp [$vb_problem get_bandstructure 0];
                    } else {
                        set vb_bands_tmp [$vb_problem get_bandstructure 1];
                    }
                    if { $pml_resort_according_to_lifetime == 1 } {
                        set vb_bands [BandstructureDomaincomplex]
                        extract_stable_bands $pml_threshold_stable_states_max_imag_value $vb_bands_tmp $vb_bands
                        $vb_bands_tmp -delete
                    } else {
                        set vb_bands $vb_bands_tmp
                    }
                    unset vb_bands_tmp
                } else {
                    $vb_problem solve $num_vb_sol $kmin $kmax $num_k_values
                    if {$cb_effmass == 1 || $num_cb_sol == 0} {
                        set vb_bands [$vb_problem get_bandstructure 0];
                    } else {
                        set vb_bands [$vb_problem get_bandstructure 1];
                    }
                }
            ##########################
            # solve eff mass problem #
            ##########################
            } else {
                $em_problem set_solution_type $holes
                $em_problem set_energy_guess  $vb_energy_guess
                ##################################
                # PML or standard?               #
                ##################################
                if { $enable_pml == 1 } {
                    if { ![info exists em_pml] } {
                        schroedinger_pml_factory em_pml $options(-model) $em_problem
                    }
                    $em_pml solve $num_vb_sol
                    if { $pml_print_bands_orderered_by_tau == 1 } {
                        $em_pml display_solution_info
                    }
                    if { $pml_resort_according_to_lifetime == 1 } {
                        set vb_bands [BandstructureDomaincomplex]
                        extract_stable_bands $pml_threshold_stable_states_max_imag_value [$em_problem get_bandstructure] $vb_bands
                    } else {
                        set vb_bands [BandstructureDomaincomplex -args [$em_problem get_bandstructure]];
                    }
                } else {
                    $em_problem solve $num_vb_sol
                    set vb_bands [$em_problem get_bandstructure];
                }
            }

            #############################
            # check bound / valid bands #
            #############################
            set num_vb_bound [find_num_bound_states $vb_bands $unbound_vb_edge $holes]
            set num_vb_valid [find_num_valid_states $vb_bands $holes]
            set num_to_keep  [$vb_bands get_number_of_bands]
            if { $binary_keep_only_valid == 1 } {
                set num_to_keep $num_vb_valid
            }
            if { $binary_keep_only_bound == 1 } {
                set num_to_keep [min $num_vb_bound $num_to_keep]
            }
            if { $num_to_keep == 0 } {
                logger_emit "TCLbands: all vb bands are invalid or unbound\nbound: $num_vb_bound, valid: $num_vb_valid\n-> won't extract any bands. inspect your bandstructure!" $LOG_ERROR
            } elseif { [$vb_bands get_number_of_bands] > $num_to_keep } {
                set vb_bands [$vb_bands extract_bands $num_to_keep 0]
            }
        }

        #########################################
        # delete PML objects                    #
        #########################################
        if { [info exists em_pml] } {
            $em_pml -delete
            unset em_pml
        }
        if { [info exists vb_pml] } {
            $vb_pml -delete
            unset vb_pml
        }

        #########################################
        # create decross object for regrouping  #
        # obtained wavefunctions such that spin #
        # degenerate bands are composed and     #
        # have same phase factor. needed to     #
        # create continuous matrix elements     #
        # which can be interpolated properly    #
        # attention: if you use PML's, your     #
        # states may not really be orthogonal!  #
        # then decrossing would fail terribly.  #
        #########################################
        if { [info exists decross_bands] && $decross_bands == 1  && $num_k_values > 1 && $enable_pml == 0 } {
            set decross [DeCross -args $geometry $material_database]
        }
        ############################
        # decross if requested     #
        ############################
        if { [info exists decross] && [info exists cb_bands] && [$cb_bands get_basis_size] > 1 && $num_k_values > 1} {
            if { [expr {[$cb_bands get_number_of_bands] % 2}] == 0 } {
                $decross resort $cb_bands
            } else {
                set tcb [$cb_bands get_number_of_bands]
                logger_emit "TCLbands: there are $tcb cb bands. i guess the are not spin degenerate. skipping decross" $LOG_WARN
            }
        }
        if { [info exists decross] && [info exists vb_bands] && [$vb_bands get_basis_size] > 1 && $num_k_values > 1} {
            if { [expr {[$vb_bands get_number_of_bands] % 2}] == 0 } {
                $decross resort $vb_bands
            } else {
                set tvb [$vb_bands get_number_of_bands]
                logger_emit "TCLbands: there are $tcb vb bands. i guess the are not spin degenerate. skipping decross" $LOG_WARN
            }
        }
        if { [info exists decross] } {
            $decross -delete
            unset decross
        }

        ############################
        # store bands if requested #
        ############################
        if { $binary_write_bands && [info exists cb_bands] } {
            $cb_bands write_binary "$result_directory/$cb_bands_binary_file"
        }
        if { $binary_write_bands && [info exists vb_bands] } {
            $vb_bands write_binary "$result_directory/$vb_bands_binary_file"
        }

        if { $num_cb_sol > 0 && $num_vb_sol > 0 } {
            logger_emit "TCLbands: cb: $num_cb_bound bound / $num_cb_valid valid subband, vb: $num_vb_bound bound / $num_vb_valid valid subband"
        } elseif { $num_cb_sol > 0 } {
            logger_emit "TCLbands: cb: $num_cb_bound bound / $num_cb_valid valid subband, vb: not calculated"
        } else {
            logger_emit "TCLbands: cb: not calculated, vb: $num_vb_bound bound / $num_vb_valid valid subband"
        }

    ####################################
    # solve kp problem for 3D problems #
    ####################################
    } else {
	
        ###############################
        # setup problem class         #
        ###############################
        if { ![info exists vb_problem] } {
            source "$sourcepath/model_init.tcl"
        }

        ###############################
        # solve conduction band first #
        ###############################
        if {$num_cb_sol > 0} {
            ##################
            # effective mass #
            ##################
            if {$cb_effmass == 1} {
                $em_problem set_energy_guess $cb_energy_guess
                set problem $em_problem
            ###################################
            # some more complex bandstructure #
            ###################################
            } else {
                $vb_problem set_energy_guess $cb_energy_guess
                set problem $vb_problem
            }
            $problem set_solution_type     $electrons
            $problem solve                 $num_cb_sol
            $problem drop_energy_guess

            # must make a copy of that
            set cb_bands [BandstructureDomaincomplex -args [$problem get_bandstructure]];
            $problem delete_solutions
            $problem set_solution_type     $holes
            $problem set_energy_guess      $vb_energy_guess

            #########################################################
            # there is now barrier sorting for 3D. so we do it here #
            # -> find states below the min cb edge (should be vb)   #
            #########################################################
            set cb_states_start [find_states_above_energy $cb_bands $min_cb_edge $electrons]

            ######################################################
            # check bound / valid bands and extract if requested #
            ######################################################
            set num_cb_bound [expr {[find_num_bound_states $cb_bands $unbound_cb_edge $electrons] - $cb_states_start}]
            set num_to_keep  [expr {[$cb_bands get_number_of_bands] - $cb_states_start}]
            if { $binary_keep_only_bound == 1 } {
                set num_to_keep $num_cb_bound
            }
            if { $num_to_keep == 0 } {
                logger_emit "TCLbands: all states are unbound\nwon't extract any states. inspect your bandstructure!" $LOG_ERROR
            } elseif { [$cb_bands get_number_of_bands] > $num_to_keep } {
                set cb_bands [$cb_bands extract_bands $num_to_keep $cb_states_start]
            }
        }

        ############################
        # store bands if requested #
        ############################
        if { $binary_write_bands && [info exists cb_bands] } {
            $cb_bands write_binary "$result_directory/$cb_bands_binary_file"
        }

        ##########################
        # solve for valence band #
        ##########################
        if {$num_vb_sol > 0} {
            if {$vb_effmass == 0} {
                set problem $vb_problem
            ##########################
            # solve eff mass problem #
            ##########################
            } else {
                set problem $em_problem
            }
            $problem set_energy_guess      $vb_energy_guess
            $problem solve $num_vb_sol
            set vb_bands [$problem get_bandstructure];

            #########################################################
            # there is now barrier sorting for 3D. so we do it here #
            # -> find states above the min vb edge (should be cb)   #
            #########################################################
            set vb_states_start [find_states_above_energy $vb_bands $min_vb_edge $holes]

            ######################################################
            # check bound / valid bands and extract if requested #
            ######################################################
            set num_vb_bound [expr {[find_num_bound_states $vb_bands $unbound_vb_edge $holes] - $vb_states_start}]
            set num_to_keep  [expr {[$vb_bands get_number_of_bands] - $vb_states_start}]
            if { $binary_keep_only_bound == 1 } {
                set num_to_keep $num_vb_bound
            }
            if { $num_to_keep == 0 } {
                logger_emit "TCLbands: all states are unbound\n-> won't extract any states. inspect your bandstructure!" $LOG_ERROR
            } elseif { [$vb_bands get_number_of_bands] > $num_to_keep } {
                set vb_bands [$vb_bands extract_bands $num_to_keep $vb_states_start]
            }
        }
        ############################
        # store bands if requested #
        ############################
        if { $binary_write_bands && [info exists vb_bands] } {
            $vb_bands write_binary "$result_directory/$vb_bands_binary_file"
        }

        if { $num_cb_sol > 0 && $num_vb_sol > 0 } {
            logger_emit "TCLbands: cb: $num_cb_bound bound state, vb: $num_vb_bound bound state"
        } elseif { $num_cb_sol > 0 } {
            logger_emit "TCLbands: cb: $num_cb_bound bound state, vb: not calculated"
        } else {
            logger_emit "TCLbands: cb: not calculated, vb: $num_vb_bound bound state"
        }
    }


}

######################################################################
# if we have effective mass subbands, create dispersive band objects #
######################################################################
if { ([info exists vb_bands] || [info exists cb_bands]) && $dimension < 3 } {

    if { [info exists cb_bands_disp] } {
        unset cb_bands_disp
    }
    if { [info exists vb_bands_disp] } {
        unset vb_bands_disp
    }

    ################
    # build domain #
    ################
    set create_new_domain 1
    if { $num_vb_sol > 0 } {
        if { [$vb_bands get_basis_size] > 1 } {
            set domain [$vb_bands get_domain]
            set create_new_domain 0
        }
    }
    if { $num_cb_sol > 0 } {
        if { [$cb_bands get_basis_size] > 1 } {
            set domain [$cb_bands get_domain]
            set create_new_domain 0
        }
    }

    if { $create_new_domain == 1 } {
        set domain [DomainMaster]
        if { $dimension == 0 } {
            create_3D_domain_radial $domain $bulk_transverse_direction $kmin $kmax $num_k_values
        } elseif { $dimension == 1 } {
            create_2D_domain_radial $domain $kmin $kmax $num_k_values
        } else {
            create_1D_domain_wire_bands $domain $kmin $kmax $num_k_values
        }
    }
    #########################
    # which material to use #
    #########################
    if { $dimension == 0 } {
        set usemat $material
    } else {
        set usemat $quantized_material
    }
    if { $num_cb_sol > 0 } {
        if { [$cb_bands get_basis_size] == 1 } {
            logger_emit [format "TCLbands: cb transverse effective mass = %g" [eval_dref [$usemat get "electron_effective_mass_transverse"]]] $LOG_INFO_DEVEL2
            set cb_bands_disp [create_effmass_dispersion [$usemat get "electron_effective_mass_transverse"] $domain $cb_bands]
        }
    }
    if { $num_vb_sol > 0 } {
        if { [$vb_bands get_basis_size] == 1 } {
            logger_emit [format "TCLbands: vb transverse effective mass = %g" [eval_dref [$usemat get "hole_effective_mass_transverse"]]] $LOG_INFO_DEVEL2
            set vb_bands_disp [create_effmass_dispersion -[$usemat get "hole_effective_mass_transverse"] $domain $vb_bands]
        }
    }
}

################################################
# write dispersion relation / energies to file #
################################################
if { [plot_requested "bands"] && ([info exists vb_bands] || [info exists cb_bands]) } {
    if { $dimension < 3 } {
        if { $num_cb_sol > 0 } {
            if { [info exists cb_bands_disp] } {
                set cb_tmp $cb_bands_disp
            } else {
                set cb_tmp $cb_bands
            }
			plot_to_file "bands" $cb_tmp "$result_directory/$cb_bands_file" "real"
        }
        if { $num_vb_sol > 0 } {
            if { [info exists vb_bands_disp] } {
                set vb_tmp $vb_bands_disp
            } else {
                set vb_tmp $vb_bands
            }
			plot_to_file "bands" $vb_tmp "$result_directory/$vb_bands_file" "real"
        }
    } else {
        if { $num_cb_sol > 0 } {
            set f [open "$result_directory/${cb_bands_file}.dat" "w"]
            for {set ii 0} {$ii < [$cb_bands get_number_of_bands]} {incr ii} {
                set en [complex2real [$cb_bands get_energy 0 $ii]]
                puts $f "$en"
            }
            close $f
        }
        if { $num_vb_sol > 0 } {
            set f [open "$result_directory/${vb_bands_file}.dat" "w"]
            for {set ii 0} {$ii < [$vb_bands get_number_of_bands]} {incr ii} {
                set en [complex2real [$vb_bands get_energy 0 $ii]]
                puts $f "$en"
            }
            close $f
        }
    }
}

#######################################
# write probability densities to file #
#######################################
if { ([info exists vb_bands] || [info exists cb_bands]) && [plot_requested "probabilities"] && $dimension > 0 } {

    #######################################################################
    # determine for which k points we should plot the probability density #
    #######################################################################
    set sp [split $plot_probabilities_k_points ":"]
    if { [llength $sp] == 2 } {
        set plot_k_start [lindex $sp 0]
        set plot_k_end   [lindex $sp 1]
    } elseif { $plot_probabilities_k_points == "all" } {
        set plot_k_start 0
        set plot_k_end   10000000000
    } else {
        set plot_k_start $plot_probabilities_k_points
        set plot_k_end   [expr { $plot_probabilities_k_points + 1 }]
    }

    if { $num_cb_sol > 0 } {
        for {set bb 0} {$bb < [$cb_bands get_number_of_bands]} {incr bb} {
            for {set ii $plot_k_start} {$ii < [$cb_bands get_number_of_k_values] && $ii < $plot_k_end} {incr ii} {
                set prob [$cb_bands get_probability $ii $bb]
				set fname       [format $cb_probability_file $bb $ii]
				plot_to_file "probabilities" $prob "$result_directory/$fname"
                delete_EigenSolutiondouble $prob
                unset prob
            }
        }
    }
    if { $num_vb_sol > 0 } {
        for {set bb 0} {$bb < [$vb_bands get_number_of_bands]} {incr bb} {
            for {set ii $plot_k_start} {$ii < [$vb_bands get_number_of_k_values] && $ii < $plot_k_end} {incr ii} {
                set prob [$vb_bands get_probability $ii $bb]
				set fname       [format $vb_probability_file $bb $ii]
				plot_to_file "probabilities" $prob "$result_directory/$fname"
                delete_EigenSolutiondouble $prob
                unset prob
            }
        }
    }
}


