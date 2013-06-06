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
# package grid.tcl                                      #
#                                                       #
# setup grid, load materials                            #
#                                                       #
# requires:                                             #
#   options(-grid) grid filename                        #
#   or options(-bulk) bulk material name                #
#   $parser InputParser object                          #
#   $material_database MaterialDatabase object          #
#                                                       #
# sets:                                                 #
#   $geometry geometry object                           #
#   $dimension dimensionality of real space             #
#   $quantized_volume volume of the quantized region    #
#   $quantized_material name of quantized material      #
#    (required to obtain transverse effective mass)     #
#########################################################

#########################################################
# load tcl material database file                       #
#########################################################
if { $material_load_from_file == 0 } {
    source "$confpath/../input/material_database.tcl"
}

#########################################################
# load geometry and check dimensionality of the problem #
#########################################################
if { [info exists options(-grid)] } {

    if { [info exists options(-deg)] } {
        set polynomial_degree $options(-deg)
    } else {
        set polynomial_degree 1
    }
    ########################################
    # check type of grid file              #
    ########################################
    set fp [open $options(-grid) "r"]
    fconfigure $fp -buffering line
    gets $fp firstline
    close $fp
    # check if its a tdkp 1d grid file, everything else is assumed to be dfise.
    # (checking for dfise would be more complicated as these files can be gzipped)
    if { $firstline == "TDKP1DGRID" } {
        logger_emit "TCLgrid: $options(-grid) seems to be a tdkp 1d grid file" $LOG_INFO_DEVEL2
        set reader [OneDGridCreator -args $options(-grid)]
    } else {
		###################################################
		# get filename ending                             #
		###################################################
		set splitted_filename [split $options(-grid) "."]
		set nargs [llength $splitted_filename]
		set ending [lindex $splitted_filename [expr $nargs - 1]]
		set pre_ending [lindex $splitted_filename [expr $nargs - 2]]
		if { $ending == "med" } {
			set reader [MEDGridReader -args $options(-grid)]
			####################################################
			# check if there is some material to region file   #
			####################################################
			set matfile [format "%s.tcl" [string range $options(-grid) 0 [expr {[string length $options(-grid)] - 5}]]]
			if { [file exists $matfile] } {
				source $matfile
				for {set ii 0} {$ii < [$reader get_num_regions]} {incr ii} {
					set reg [string2char [$reader get_region_name $ii]]	
					if { [info exists med_materials($reg)] } {
						logger_emit "Grid: assigning $med_materials($reg) to region $reg" $LOG_INFO_DEVEL1
						$reader set_region_material_name $ii $med_materials($reg)
					} elseif { [info exists med_materials(default)]} {
						logger_emit "Grid: region $reg gets default material $med_materials(default)" $LOG_INFO_DEVEL1
						$reader set_region_material_name $ii $med_materials(default)
					} else {
						logger_emit "Grid: there is no material assigned to region $reg. You have to define it in a .tcl file (instead of .med) via assignments such as 'set med_materials(default) \"GaAs\"' or set med_materials(the_wire) \"GaAs\"'" $LOG_ERROR
						exit 1
					}
				}
			}
	    } elseif { $ending == "asc" || ($ending == "gz" && $pre_ending == "asc") } {
			set reader [AsciiGridReader -args $options(-grid)]
		} else {
	        set reader [SebiseGridReader -args $options(-grid)]
		}
    }

	####################################################
	# build the grid                                   #
	####################################################
    set geometry [$parser read_geometry $reader $polynomial_degree]
    $reader -delete
    unset reader

	####################################################
	# load the material database from file             #
	####################################################
    if { $material_load_from_file == 1 } {
        load_materials_from_file $geometry $material_database
    } else {
        load_materials_from_tcl_db $geometry $material_database $temperature
        if { $material_file_dump == 1 } {
            write_materials_to_files $material_database $result_directory
        }
    }
    $geometry set_materials $material_database
    set dimension [$geometry get_dimension]

    #######################################
    # scale node coordinates on request   #
    #######################################
    if { $scale_coordinates != 1 } {
        $geometry rescale_node_coordinates $scale_coordinates
    }

    ########################################
    # check if PMLs are present            #
    ########################################
    if { $enable_pml == 1 } {
        set pml_exists 0
        for {set ii 0} {$ii < [$geometry get_num_regions]} {incr ii} {
            set rname [string2char [[$geometry get_region $ii] get_name]]
            if { [string range $rname 0 2] == "PML" } {
                set pml_exists 1
            }
        }
        if { $pml_exists == 0 } {
            logger_emit "no PMLs have be defined therefore disabling PML calculation"
            set enable_pml 0
        }
    }

    ###########################
    # find quantized material #
    ###########################
    determine_quantized_material $geometry $material_database quantized_material
    logger_emit "TCLgrid: using $quantized_material as quantized material " $LOG_INFO_DEVEL1
    set quantized_material [$material_database get_material $quantized_material]

    ###############################
    # calculate quantized volume  #
    # or use regions from file to #
    # calculate it                #
    ###############################
    if { $quantized_volume_predef > 0.0 } {
        set quantized_volume $quantized_volume_predef
        logger_emit "TCLgrid: using user defined quantized volume: $quantized_volume"
    } else {
        if { $read_quantized_regions_from_file != "" && [file exists $read_quantized_regions_from_file]} {
            set quantized_calculator [UserDefinedQuantizedVolumeCalculator -args $geometry]
            set f [open $read_quantized_regions_from_file]
            while { [gets $f line] >= 0} {
                $quantized_calculator add $line
            }
        } else {
            set quantized_calculator [QuantizedVolumeCalculator -args $geometry]
        }
        set quantized_volume   [$quantized_calculator calculate_volume ]
        if { $quantized_volume == 0.0 } {
            logger_emit "TCLgrid: quantized volume is ZERO! means: i could not determine the quantized volume correctly. please write the names of the quantized regions into the file." $LOG_WARN
        } else {
            logger_emit "TCLgrid: quantized volume $quantized_volume"
        }
    }

    #########################################
    # use user defined boundary conditions? #
    # tcl file source, expect bc_bands or   #
    # bc_strain to be defined               #
    #########################################
    if { [info exists options(-bc)] } {
        source $options(-bc)
        if { ![info exists bc_bands] && ![info exists bc_strain] && ![info exists bc_poisson]  } {
            logger_emit "TCLgrid: user requested user-defined boundary conditions but didn't specify them" $LOG_ERROR
            quit
        }
    }

    #########################################
    # load potential on request             #
    #########################################
    if { [info exists options(-potential)] } {
        if { ![file exists $options(-potential)] } {
            logger_emit "TCLgrid: requested potential file $options(-potential) does not exist!" $LOG_ERROR
            quit
        } else {
            logger_emit "TCLgrid: reading potential from file $options(-potential)." $LOG_INFO_DEVEL1
        }
        ###################################################
        # different handling if pot file is binary (.bin) #
        ###################################################
        set potential_field [PotentialEnergyField -args [$geometry get_num_nodes]]
        if { "bin" == [lindex [split $options(-potential) "."] [expr {[llength [split $options(-potential) "."]] - 1}]] } {
            $parser read_binary $potential_field $options(-potential)
            if {[$geometry get_num_nodes] != [$potential_field get_length]} {
                set plength [$potential_field get_length]
                set nnodes  [$geometry get_num_nodes]
                logger_emit "TCLgrid: potential file $options(-potential) has length $plength, while num nodes is $nnodes!" $LOG_ERROR
                quit
            }
        } else {
            set f [open $options(-potential) "r"]
            set pot [list]
            while { [gets $f line] >= 0 } {
                lappend pot $line
            }
            close $f

            if { [llength $pot] == [$geometry get_num_nodes] } {
                for {set ii 0} {$ii < [llength $pot]} {incr ii} {
                    $potential_field set_node_value $ii 0 [lindex $pot $ii]
                }
            } else {
                set flength [llength $pot]
                set nnodes  [$geometry get_num_nodes]
                logger_emit "TCLgrid: requested potential file $options(-potential) has length $flength, while num nodes is $nnodes!" $LOG_ERROR
                quit
            }
        }
    }
} elseif { [info exists options(-bulk)] } {

    if { $material_load_from_file == 1 } {
        $material_database load_material $options(-bulk)
        set material [$material_database get_material $options(-bulk)]
    } else {
        # get material properties
        get_material_object mat $options(-bulk)
        set_temperature_dependent_properties mat $temperature
        # collect all material properties
        foreach p [array names mat] {
            set db(0,$p) $mat($p)
        }
        # set bandedge to zero
        set db(0,vb_edge) 0.0
        set db(0,cb_edge) $db(0,bandgap)

        set material [get_tdkp_material_object 0 db]
        $material_database add_material $options(-bulk) $material
        if { $material_tcl_dump == 1 } {
            puts "-> TCLgrid: Material $options(-bulk)"
            $material dump
        }
        if { $material_file_dump == 1 } {
            write_materials_to_files $material_database $result_directory
        }
    }
    set dimension 0
    set quantized_volume 1.0
}
