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

# ---------------------------------------------------------------
#
# tdkp_helpers.tcl - collection of tcl functions that may help
# using the tcl interface
#
# (c) ratko veprek, iis, eth zuerich, nov. 2007
#
# procedures implemented:
# - find_num_bound_states{bands, max_edge, solution_type}
# - find_num_valid_states{bands, solution_type}
# - min {args} (args is unlimited number of parameters)
# - max {args}
# - matrix_elements_find_degenerate_transitions {matrix_elements, cb_bands, vb_bands}
# - check_degenerate_bands {bands, idx_1, idx_2, tol}
# - determine_minmax_edges {geometry material_database max_cb_edge max_vb_edge min_cb_edge min_vb_edge}
# - determine_quantized_material {geometry material_database variable}
# - determine_command_line_options { variable }
# - load_materials_from_file {geometry material_database}
# - exists_in_list { needle haystack }
# - find_reference_lattice_constant { geometry material_database }
# - check_and_process_user_input { options }
# - schroedinger_pml_factory { target model dimension }
# - plot_requested { optionname } 
# - plot_to_file { optionname plot_data basefilename }
# - process_command_line_plot_requests
# ---------------------------------------------------------------

# --------------------------------------------------------------------
# return the number of states that are initially below
# the specified edge
# --------------------------------------------------------------------
proc find_num_bound_states {bands max_edge solution_type} {

    set num_bands [$bands get_number_of_bands]
    global holes
    global electrons
    if {$solution_type == $holes} {
        set order -1.0
    } elseif {$solution_type == $electrons} {
        set order 1.0
    } else {
        puts "unknown solution type ${solution_type}"
        return 0
    }
    for {set ii 0} {$ii < $num_bands} {incr ii} {
        set energy [complex2real [$bands get_energy 0 $ii]]
        if {[expr {$energy * $order}] > [expr {$max_edge * $order}]} {
            return $ii;
        }
    }
    return $num_bands
}

# --------------------------------------------------------------------
# return the state idx where the valid states above the minedge
# begin
# --------------------------------------------------------------------
proc find_states_above_energy {bands min_edge solution_type} {

    global electrons

    for {set ii 0} {$ii < [$bands get_number_of_bands]} {incr ii} {
        set energy [complex2real [$bands get_energy 0 $ii]]
        if { $solution_type == $electrons } {
            if { $energy >= $min_edge } {
                return $ii
            }
        } else {
            if { $energy <= $min_edge } {
                return $ii
            }
        }
    }
    logger_output "no state above the band edge found!" $LOG_ERROR
    exit
}

# --------------------------------------------------------------------
# return the number of states that have a valid (non-jumping,
# induced by barrier sorting or bad sigma) dispersion
# --------------------------------------------------------------------
proc find_num_valid_states {bands solution_type} {
    global holes
    global electrons
    if {$solution_type == $holes} {
        set order -1.0
    } elseif {$solution_type == $electrons} {
        set order 1.0
    } else {
        puts "unknown solution type ${solution_type}"
        return 0
    }
    for {set ii 1} {$ii < [$bands get_number_of_bands]} {incr ii} {
        for {set jj 0} {$jj < [$bands get_number_of_k_values]} {incr jj} {
            set preii [expr {$ii - 1}]
            set e_pre_ii [complex2real [$bands get_energy $jj $preii]]
            set e_ii     [complex2real [$bands get_energy $jj $ii]]
            if {[expr {$e_pre_ii * $order}] > [expr {$e_ii * $order}]} {
                return $ii
            }
        }
    }
    return [$bands get_number_of_bands]
}

# --------------------------------------------------------------------
# return minimum in argument list
# --------------------------------------------------------------------
proc min {args} {
    set num [llength $args]
    if {$num == 0} {
        puts "no arguments passed to min"
        return 0
    }
    set ret [lindex $args 0]
    for {set ii 1} {$ii < $num} {incr ii} {
        set ta [lindex $args $ii]
        if {$ta < $ret} {
            set ret $ta
        }
    }
    return $ret
}

# --------------------------------------------------------------------
# return maximum in argument list
# --------------------------------------------------------------------
proc max {args} {
    set num [llength $args]
    if {$num == 0} {
        puts "no arguments passed to max"
        return 0
    }
    set ret [lindex $args 0]
    for {set ii 1} {$ii < $num} {incr ii} {
        set ta [lindex $args $ii]
        if {$ta > $ret} {
            set ret $ta
        }
    }
    return $ret
}


# --------------------------------------------------------------------
# determine upper and lower bounds for our bound particles
# --------------------------------------------------------------------
proc determine_deprecated_minmax_edges {geometry material_database max_cb_edge max_vb_edge min_cb_edge min_vb_edge} {

    # reference to passed variables
    upvar $max_cb_edge my_max_cb_edge
    upvar $max_vb_edge my_max_vb_edge
    upvar $min_cb_edge my_min_cb_edge
    upvar $min_vb_edge my_min_vb_edge

    for {set ii 0} {$ii < [$geometry get_num_regions]} {incr ii} {
        if { [[$geometry get_region $ii] enabled] == 1 } {
            set matname  [string2char [[$geometry get_region $ii] get_material_name]]
            set material [$material_database get_material $matname]
            set cb_edge  [$material get "conduction_band_edge"]
            set vb_edge  [$material get "valence_band_edge"]
            if {$ii == 0} {
                set my_max_cb_edge $cb_edge
                set my_max_vb_edge $vb_edge
                set my_min_cb_edge $cb_edge
                set my_min_vb_edge $vb_edge
            } else {
                set my_max_cb_edge [max $cb_edge $my_max_cb_edge]
                set my_max_vb_edge [min $vb_edge $my_max_vb_edge]
                set my_min_cb_edge [min $cb_edge $my_min_cb_edge]
                set my_min_vb_edge [max $vb_edge $my_min_vb_edge]
            }
        }
    }
}


# --------------------------------------------------------------------
# determine material with lowest cb edge
# --------------------------------------------------------------------
proc determine_quantized_material {geometry material_database variable} {

    upvar $variable quantized
    set quantized ""
    for {set ii 0} { $ii < [$geometry get_num_regions] } {incr ii} {
        if { [[$geometry get_region $ii] enabled] == 1 } {
            set matname [string2char [[$geometry get_region $ii] get_material_name]]
            if { $matname != "Gas" } {
                if { $quantized == "" } {
                    set quantized $matname
                } else {
                    if { [$material_database get $quantized "conduction_band_edge"] > [$material_database get $matname "conduction_band_edge"]} {
                        set quantized $matname
                    }
                }
            }
        }
    }
}


# --------------------------------------------------------------------
# parse command line arguments
# --------------------------------------------------------------------
proc determine_command_line_options { variable } {

    upvar $variable options

    for {set ii 0} { $ii < $::argc } {incr ii} {
        set argument [lindex $::argv $ii]
        # check if it is an option
        if { [string index $argument 0] == "-" } {
            # check if next items is an option or a value
            if { [expr {$ii + 1}] < $::argc } {
                set opt [lindex $::argv [expr {$ii + 1}]]
                if { [string index $opt 0] != "-" } {
                    set options($argument) $opt
                    incr ii
                } else {
                    # next is an value, so indicate that option is set
                    set options($argument) 1;
                }
            } else {
               set options($argument) 1;
            }
        } else {
            show_bands_tcl_options "can not parse $argument"
            exit
        }
    }
}

# --------------------------------------------------------------------
# load materials of geometry object into database
# --------------------------------------------------------------------
proc load_materials_from_file {geometry material_database} {
    for {set ii 0} { $ii < [$geometry get_num_regions] } {incr ii} {
        set matname [string2char [[$geometry get_region $ii] get_material_name]]
        if { ![$material_database material_exists $matname] } {
            $material_database load_material $matname
        }
    }
}

# --------------------------------------------------------------------
# check if a given item exists in a list
# --------------------------------------------------------------------
proc exists_in_list { needle haystack } {

    foreach p $haystack {
        if { $p == $needle } {
            return 1
        }
    }

    return 0

}

# --------------------------------------------------------------------
# determine the material for the reference lattice constant
# --------------------------------------------------------------------
proc find_reference_lattice_constant_material { geometry material_database lattice_constant_field } {

    set lattice_constants [list]
    set reference_idx -1

    # get lattice constant of every boundary material
    for {set ii 0} {$ii < [$material_database get_num_materials]} {incr ii} {
        if { [$geometry boundary_material $ii] && [[$material_database get_material $ii] valid_key $lattice_constant_field]} {
            lappend lattice_constants [[$material_database get_material $ii] get $lattice_constant_field]
            if { $reference_idx == -1 } {
                set reference_idx $ii
            }
        }
    }
    # check lattice constants
    if { [llength $lattice_constants] == 0 } {
        puts "-> ** error ** - no boundary materials available. can not determine reference lattice constant"
        exit
    }
    set reference [lindex $lattice_constants 0]
    set good 1
    for {set ii 1} {$ii < [llength $lattice_constants]} {incr ii} {
        set other [lindex $lattice_constants $ii]
        if { [expr {abs($other - $reference) > 0.0001}] } {
            set good 0
        }
    }
    if { $good == 1 } {
        return $reference_idx
    }
    puts "-> ** error ** - lattice constants of boundary materials differ. please define it by yourself."
    for {set ii 0} {$ii < [$material_database get_num_materials]} {incr ii} {
        if { [$geometry boundary_material $ii] } {
            set lattice_const [[$material_database get_material $ii] get $lattice_constant_field]
            set matname       [string2char [$material_database get_material_name $ii]]
            puts "   $matname = $lattice_const"
        }
    }
    exit

}

# ------------------------------------------------------------
# check user input to tow.tcl
# ------------------------------------------------------------
proc check_and_process_user_input { variable } {

    upvar $variable options

    global result_directory
    global binary_matrix_elements_file
    global coulomb_binary_file
    global vb_bands_binary_file

    # check if any calculation is requested
    set valid_option [list -bands -matrix_elements -slc -slc_fermi -fermi -strain -coulomb -iso -plot -tclscript -clc -clc_fermi -piezo]
    set valid 0
    foreach f $valid_option {
        if { [info exists options($f)] } {
            set valid 1
        }
    }
    if { $valid == 0 } {
        show_bands_tcl_options "please choose at least one simulation option"
        exit
    }
    # either bulk or a grid
    if { [info exists options(-bulk)] && [info exists options(-grid)] } {
        show_bands_tcl_options "-bulk and -grid are exclusive requests"
        exit
    } elseif { ![info exists options(-bulk)] && ![info exists options(-grid)] } {
        show_bands_tcl_options "-bulk or -grid must be set"
        exit
    }
    # model must be set
    if { ![info exists options(-model)] } {
        show_bands_tcl_options "setting the kp model via -model is mandatory"
        exit
    }
    # strain only exists for quantized stuff ...
    if { [info exists options(-strain)] && [info exists options(-bulk)] } {
        show_bands_tcl_options "switch strain is only implemented for quantized bandstructures!\nfor bulk strain, please create a StrainTensor in the config file named bulk_strain."
        exit
    }
    # piezo only exists for quantized stuff ...
    if { [info exists options(-piezo)] && [info exists options(-bulk)] } {
        show_bands_tcl_options "switch piezo is nonsense for bulk crystals!"
        exit
    }
    # delta strain requires strain
    if { [info exists options(-delta_strain)] && !([info exists options(-strain)] || [info exists options(-load_strain)]) } {
        show_bands_tcl_options "delta_strain requires strain calculation"
        exit
    }
    # slc_fermi and fermi are exclusive
    if { [info exists options(-slc_fermi)] && [info exists options(-fermi)] } {
        show_bands_tcl_options "slc_fermi and fermi are exclusive!"
        exit
    }
    # ensure that we don-t have load and calculate at the same time
    if { [info exists options(-load_bands)] && [info exists options(-bands)] } {
        show_bands_tcl_options "option bands and load_bands are exclusive"
        exit
    }
    if { [info exists options(-load_matrix_elements)] && [info exists options(-matrix_elements)] } {
        show_bands_tcl_options "option matrix_elements and load_matrix_elements are exclusive"
        exit
    }
    # determine what must be calculated
    if { [info exists options(-slc_fermi)] } {
        set options(-slc) 1
    }
    if { [info exists options(-clc_fermi)] } {
        set options(-clc) 1
    }
    # check if bands are needed
    set bands_needed 0
    if { [info exists options(-slc)] || [info exists options(-matrix_elements)] || [info exists options(-fermi)] || [info exists options(-slc_fermi)] || [info exists options(-coulomb)] || [info exists options(-iso)] || [info exists options(-coulomb)] || [info exists options(-clc)] } {
        set bands_needed 1
    }
    # check if momentum matrix elements are needed
    set matrix_elements_needed 0
    if { [info exists options(-slc)] || [info exists options(-clc)] } {
        set matrix_elements_needed 1
    }
    # check if coulomb matrix elements are needed
    set coulomb_matrix_elements_needed 0
    if { [info exists options(-clc)] } {
        set coulomb_matrix_elements_needed 1
    }
    #######################################################################
    # coulomb matrix elements -> all or diag: diag is n1 = n4 and n2 = n3 #
    #######################################################################
    if { [info exists options(-coulomb)] } {
        if { [info exists options(-load_coulomb)] } {
            show_bands_tcl_options "options -coulomb and -load_coulomb are exclusive"
            exit
        }
        if { ($options(-coulomb) != "all" && $options(-coulomb) != "diag" && [llength [split $options(-coulomb) ":"]] == 0) || $options(-coulomb) == 1 } {
            show_bands_tcl_options "invalid options for coulomb matrix elment"
            exit
        }
    }

    ###############################################
    # check if cached band structure is available #
    # and if bandstructure is needed              #
    ###############################################
    if { ![info exists options(-bands)] } {
        # calculate bands if
        if { ![info exists options(-nocache)] } {
            set bands_available [file exists "$result_directory/$vb_bands_binary_file"]
        } else {
            set bands_available 0
        }
        if { $bands_needed } {
            if { $bands_available } {
                set options(-load_bands) 1
            } else {
                set options(-bands) 1
            }
        }
    }
    #############################################
    # matrix elements?                          #
    #############################################
    if { ![info exists options(-matrix_elements)] } {
        if { ![info exists options(-nocache)] } {
            set matrix_element_available [file exists "$result_directory/${binary_matrix_elements_file}"]
        } else {
            set matrix_element_available 0
        }
        if { $matrix_elements_needed } {
            if { $matrix_element_available } {
                set options(-load_matrix_elements) 1
            } else {
                set options(-matrix_elements) 1
            }
        }
    }
    #############################################
    # fermi levels?                             #
    #############################################
    if { ![info exists options(-fermi)] && [info exists options(-slc)] } {
        set options(-fermi) $options(-slc)
    }
    if { ![info exists options(-fermi)] && [info exists options(-clc)] } {
        set options(-fermi) $options(-clc)
    }
    #############################################
    # diagonal coulomb matrix elements needed?  #
    #############################################
    if { $coulomb_matrix_elements_needed == 1 && ![info exists options(-coulomb)] && ![info exists options(-load_coulomb)] } {
        if { ![info exists options(-nocache)] } {
            set coulomb_available [file exists "$result_directory/$coulomb_binary_file"]
        } else {
            set coulomb_available 0
        }
        if { $coulomb_available == 1 } {
            set options(-load_coulomb) 1
        } else {
            set options(-coulomb) "diag"
        }
    }
}





# ------------------------------------------------------------
# create isosurface plots
# ------------------------------------------------------------
proc create_isosurface_plots { request geometry cb_bands vb_bands } {

    global result_directory

    if { [$geometry get_dimension] == 3 } {
        set sp [split $request ":"]
        set isosurf [IsoSurfaceGenerator]
        $isosurf set_geometry $geometry
        set ii 0
        foreach pl $sp {
            if { [regexp {^([cbvb]{2})([0-9]+)$} $pl matches band_type idx] } {
                if { $band_type == "cb" } {
                    set prob [$cb_bands get_probability 0 $idx]
                } else {
                    set prob [$vb_bands get_probability 0 $idx]
                }
                set maxvalue 0.0
                for {set jj 0} {$jj < [$prob get_length]} {incr jj} {
                    set maxvalue [max $maxvalue [eval_dref [$prob get_vertex_value $jj]]]
                }
                puts "maxvalue is $maxvalue"
                $isosurf set_data $prob
                if { $ii == 0 } {
                    $isosurf create_lists
                    set ii 1
               }
                set num 8
                set dn  [expr {$maxvalue / ($num + 2)}]
                for {set jj 0} {$jj < $num} {incr jj} {
       #             set isosurf [IsoSurfaceGenerator]
      #              $isosurf set_geometry $geometry
      #              $isosurf set_data $prob
       #             $isosurf create_lists
                    set n     [expr {($jj + 1.0) * $dn}]
                    set fname [format "$result_directory/%s_%d.pov" $pl $jj]
                    puts "  doing plot for $pl $n"
                    $isosurf create_surface_plot $fname $n
                }

            } else {
                puts " ** error **: invalid iso request $pl ..."
                exit
            }
        }
#        set isosurf [IsoSurfaceGenerator]
#        $isosurf set $geometry
#        $isosurf create_lists

    } else {
        puts " ** error **: isosurfaces can only be plotted in 3D ..."
        exit
    }
}

###################################################
# replace + and - from string (for plt file names #
###################################################
proc make_filename_plt_compatible { w } {

    regsub -all {\-} $w  "m" w1
    regsub -all {\+} $w1 "p" w2
    regsub -all {_}  $w2 "u" w3
    return $w3

}

###################################################
# create appropriate PML wrapper                  #
###################################################
proc schroedinger_pml_factory { target model base_problem } {

    upvar $target pml
    global electrons
    global holes
    global pml_alpha
    global pml_beta
    global pml_power

    set solution_type [$base_problem get_solution_type]

    #################################
    # easy catch for effmass        #
    #################################
    set effmass 0
    if { $model == "effmass" } {
        set effmass 1
    } elseif { $solution_type == $electrons } {
        if { $model != "kp8x8" && $model != "kp8x8WZ" } {
            set effmass 1
        }
    }
    if { $effmass == 1 } {
        set pml [EffectiveMassPML -args $base_problem]
    } else {
        #################################
        # create appropriate pml object #
        #################################
        if { [[$base_problem get_geometry] get_dimension] == 3 } {
            if { $model == "kp4x4" } {
                set pml [KP4x43DPML -args $base_problem]
            } elseif { $model == "kp6x6" } {
                set pml [KP6x63DPML -args $base_problem]
            } elseif { $model == "kp8x8" } {
                set pml [KP8x83DPML -args $base_problem]
            } elseif { $model == "kp6x6WZ" } {
                set pml [KP6x63DWZPML -args $base_problem]
            } elseif { $model == "kp8x8WZ" } {
                set pml [KP8x83DWZPML -args $base_problem]
            } else {
                logger_emit "not implemented model!" $LOG_ERROR
                exit
            }
        } else {
            if { $model == "kp4x4" } {
                set pml [KP4x41D2DPML -args $base_problem]
            } elseif { $model == "kp6x6" } {
                set pml [KP6x61D2DPML -args $base_problem]
            } elseif { $model == "kp8x8" } {
                set pml [KP8x81D2DPML -args $base_problem]
            } elseif { $model == "kp6x6WZ" } {
                set pml [KP6x61D2DWZPML -args $base_problem]
            } elseif { $model == "kp8x8WZ" } {
                set pml [KP8x81D2DWZPML -args $base_problem]
            } else {
                logger_emit "not implemented model!" $LOG_ERROR
                exit
            }
        }
    }
    $pml set_stretch_parameters $pml_alpha $pml_beta $pml_power
}

# ------------------------------------------------------
# processes command line plot requests and populate 
# internal plot_requests list
# ------------------------------------------------------
proc process_command_line_plot_requests { command_line_args } {

	global plot_request
	global plot_valid_options
	global plot_file_formats
	global plot_default_format
	global LOG_ERROR

	# is the output option 'all' contained in the command line?
	if { [lsearch $command_line_args "all"] != -1 } {
		set command_line_args $plot_valid_options
	}

	# modifiy requests such that requests on 'default values' are modified to
	# contain the default file format (easier to handle later on)
	foreach request $plot_request {
		set idx [lsearch $plot_valid_options $request]
		if { $idx != -1 } {
			set plot_request [lreplace $plot_request $idx 1 "${request}_${plot_default_format}"]
		}
	}
	
	# just build a list of all possible requests
	set valid_requests [list]
	foreach opt $plot_valid_options {
		# add default plot style
		lappend valid_requests $opt
		# add each combo
		foreach type $plot_file_formats {
			lappend valid_requests "${opt}_${type}"
		}
	}
	# test for every request if its valid
	foreach request $command_line_args {
		if { [lsearch $valid_requests $request] == -1 } {
			# this is an error ... 
			set out	"TCLplot: output request '${request}' is invalid. Valid options are: "
			foreach opt $plot_valid_options {
				set out "${out} ${opt}"
			}
			set out "${out}. This options can also be extended with a valid output type specifier: "
			foreach type $plot_file_formats {
				set out "${out} _${type}"
			}
			logger_emit "$out\n" $LOG_ERROR
			exit
		}
		# good option, so append default format in case that there was no default format requested
		set modified_request $request
		if { [lsearch $plot_valid_options $request] != -1 } {		
			set modified_request "${request}_${plot_default_format}"
		}
		# only append if this request is not already included 
		if { [lsearch $plot_request $modified_request] == -1 } {	
			lappend plot_request $modified_request
		}
	}
}

# ------------------------------------------------------
# test whether plot of type 'optionname' has been 
# requested by user (be it in default file format or
# in a specific file format)
# ------------------------------------------------------
proc plot_requested { optionname } {
	global plot_request
	global plot_valid_options
	global plot_file_formats
	global LOG_ERROR
	# iterative over plot_valid_options and check if optionname is valid
	set found 0
	foreach valid $plot_valid_options {
		if { $valid == $optionname } {
			set found 1
			break
		}
		foreach type $plot_file_formats {
			if { "${valid}_{$type}" == $optionname } {
				set found 1
				break
			}
		}
	}
	if { $found == 0 } {
		set out "TCLplot: the output option $optionname does not exist. Valid options are:"
		foreach opt $plot_valid_options {
			set out "${out} ${opt}"
		}
		logger_emit $out $LOG_ERROR
		exit 1
	} 
	# check if there is a 'default' file format request
	if { [lsearch $plot_request $optionname] != -1 } {
		return 1
	}
	# check if there is a plot request for it
	foreach type $plot_file_formats {
		if { [lsearch $plot_request "${optionname}_${type}"] != -1 } {
			return 1
		}
	}
	return 0
}

# ------------------------------------------------------
# processes plot requests and dispatches it to the
# correct output writer
# ------------------------------------------------------
proc plot_to_file { optionname plot_data basefilename { value_change ""} } {

	global plot_request
	global plot_valid_options
	global plot_file_formats
	global plot_default_format
	global __writer_med
	global geometry
	
	# this list takes the types for which we have to create outputs
	set plot_types [list]

	# determine which output types are requested		
	foreach type $plot_file_formats {
		# is type defined in request but has not yet been included in the type list?
		if { [lsearch $plot_request "${optionname}_${type}"] != -1 && [lsearch $plot_types $type] == -1} {

			# easy fix: in 0D, we don't have a geometry object, but we only need XY data plots -> we just allow ascii
			# and binary for the moment. for others, we need to change the constructors
			if { ![info exists geometry] && $type != "binary" } {
				if { [lsearch $plot_types "ascii" ] == -1 } {
					lappend plot_types "ascii"
				}

			# quick check: if med is requested as output but our grid is only 1D, then we 
			# automatically switch to ascii
			} elseif { [$geometry get_dimension] == 1 && $type == "med" } {
				if { [lsearch $plot_types "ascii" ] == -1 } {
					lappend plot_types "ascii"
				} 
			} else {
				lappend plot_types $type
			}
		}
	}

	# plot for each type
	foreach type $plot_types {
		# create writer object
		if { $type == "med" } {
			global options
			global result_directory
			# is there an old one?
			if { ![info exists __writer_med] } {
				set basefile [string range $options(-grid) 0 [expr {[string length $options(-grid)] - 5}]]
				# ensure that we don't overwrite old results
				set output_med "${result_directory}/${basefile}_out.med"
				set ii 0
				while { [file exists $output_med] } {
					set output_med "${result_directory}/${basefile}_out.${ii}.med"
					incr ii
				}				
				set __writer_med [TCLDataIOMED -args [MEDDataIO -args $geometry $output_med]]
			}
			set writer $__writer_med
		} elseif { $type == "sebise" } {
			set writer [TCLDataIOSebise -args [SebiseDataIO -args $geometry]]
		} elseif { $type == "ascii" } {
			if { [info exists geometry] } {
				set writer [TCLDataIOAscii -args [AsciiDataIO -args $geometry]]
			} else {
				set writer [TCLDataIOAscii -args [AsciiDataIO]]
			}
		} elseif { $type == "binary" } {
			set writer [TCLDataIOBinary -args [BinaryDataIO]]
		} else {
			global LOG_ERROR
			logger_emit "TCLplot: invalid output type requested: ${type}" $LOG_ERROR
			exit
		}
		# write data
		if { $value_change == "real" } {
			$writer write_real $plot_data $basefilename	
		} else {
			$writer write $plot_data $basefilename	
		}		
		unset writer
	}
}

