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

###############################################################
#  material_database.tcl                                      #
#                                                             #
#  collection of material properties                          #
#                                                             #
#  (c) ratko veprek, iis/ethz, jan. 2008                      #
#      veprek@iis.ee.ethz.ch                                  #
#                                                             #
#  purpose:                                                   #
#  - provide a commmon syntax / place for material paramters  #
#  - InAs / GaAs / AlGaAs parameters are taken from           #
#    Vurgaftman, Meyer and Ram-Mohan, Journal of Applied      #
#    Physics, Volume 89 (2001) [1]                            #
#  - old GaN / InN parameters are taken from                  #
#    Journal of Applied Physics, Volume 94 (2003) [2]         #
#  - new GaN / InN / AlN parameters are taken from sebi's     #
#    aqua simulator and are based on the values given by      #
#    vurgaftman in the corresponding book chapters of         #
#    Nitride Semiconductor Devices, Principles and Simulation #
#    edited by Joachim Piprek, wiley 2007 [3]                 #
#  - some parameters are also taken from sadao adachi,        #
#    properties of group IV, III-V and II-VI semiconductors,  #
#    wiley 2005 [4]                                           #
#                                                             #
###############################################################


#####################################################
# return material object of given matierial name    #
#####################################################
proc get_material_object { target matname } {

    upvar $target material
    global temperature
    set    update_temp_dependence 1

    # bare materials
    if { $matname == "InAs" } {
        param_InAs material
    } elseif { $matname == "GaAs" } {
        param_GaAs material
    } elseif { $matname == "AlAs" } {
        param_AlAs material
    } elseif { $matname == "InN" } {
        param_InN material
    } elseif { $matname == "GaN" } {
        param_GaN material
    } elseif { $matname == "AlN" } {
        param_AlN material
    } elseif { $matname == "InP" } {
        param_InP material
    } elseif { $matname == "Gas" } {
        set update_temp_dependence 0
        param_Gas material
    # ternary alloys?
    } else {
        # arsenides?
        if { [regexp {([AlGaIn]+)As([0-9\.]+)} $matname match ternary_name mf] } {
            if { $ternary_name == "AlGa" } {
                param_GaAs mat1
                param_AlAs mat2
                bowing_AlGaAs bowing
            } elseif { $ternary_name == "InGa" } {
                param_GaAs mat1
                param_InAs mat2
                bowing_InGaAs bowing
            } else {
                puts "-> ** error ** : unknown ternary material $matname $ternary_name $mf"
                exit
            }
        # ternary nitrides
        } elseif { [regexp {([AlInGa]+)N([0-9\.]+)} $matname match ternary_name mf] } {
            if { $ternary_name == "InGa" } {
                param_InN mat2
                param_GaN mat1
                bowing_InGaN bowing
            } elseif { $ternary_name == "AlGa" } {
                param_AlN mat2
                param_GaN mat1
                bowing_InGaN bowing
            } else {
                puts "-> ** error ** : unknown ternary material $matname $ternary_name $mf"
                exit
            }
        } else {
            puts "-> ** error ** : unknown material $matname"
            exit
        }
        # first, update bandgap etc, and then interpolate
        set_temperature_dependent_properties mat1 $temperature
        set_temperature_dependent_properties mat2 $temperature
        create_ternary_alloy material $mf mat1 mat2 bowing
        set material(identifier) $matname
        set update_temp_dependence 0
    }
    if { $update_temp_dependence == 1 } {
        set_temperature_dependent_properties material $temperature
    }

}

#####################################################
# check kane for the ellipticity of the equation    #
# system                                            #
#####################################################
proc zincblende_kp_non_ellipticity_ratio { L M N Nm } {

    set Np [expr {$N - $Nm}]
    lappend ev [expr {$M - $Nm}]
    lappend ev [expr {$M - $Nm}]
    lappend ev [expr {$M - $Nm}]
    lappend ev [expr {$M + $Nm}]
    lappend ev [expr {$M + $Nm}]
    lappend ev [expr {$M + $Nm}]
    lappend ev [expr {$L - $Np}]
    lappend ev [expr {$L - $Np}]
    lappend ev [expr {$L + 2.0 * $Np}]

    set plus  0.0
    set minus 0.0
    foreach p $ev {
        if { $p > 0.0 } {
            set plus [expr {$plus + $p}]
        } else {
            set minus [expr {$minus - $p}]
        }
    }
    if { $minus == 0.0 } {
        return 1.0
    } else {
        return [expr {$plus / $minus} ]
    }

}

#####################################################
# write material properties to files                #
#####################################################
proc write_materials_to_files { material_database resdir } {

    for {set ii 0} {$ii < [$material_database get_num_materials]} {incr ii} {
        set fname [format "%s/%s.mat" $resdir [string2char [$material_database get_material_name $ii]]]
        [$material_database get_material $ii] write_to_file $fname
    }
}

#####################################################
# calculate kp 8x8 luttinger params                 #
#####################################################
proc zincblende_kp_8x8_luttinger { target variable } {

    upvar $variable mat
    upvar $target   lut
    global constants_hbar
    global constants_m0

    set lut(Ac) [expr {1.0 / $mat(me_cb) - 2.0 / 3.0 * $mat(optical_matrix) / $mat(bandgap) - 1.0 / 3.0 * $mat(optical_matrix) / ($mat(bandgap) + $mat(spin_orbit))}]

    set lut(1)  [expr {$mat(lut1) - 1.0 / 3.0 * $mat(optical_matrix) / $mat(bandgap)} ]
    set lut(2)  [expr {$mat(lut2) - 1.0 / 6.0 * $mat(optical_matrix) / $mat(bandgap)} ]
    set lut(3)  [expr {$mat(lut3) - 1.0 / 6.0 * $mat(optical_matrix) / $mat(bandgap)} ]


}

#####################################################
# set optimal Nm for kp 8x8                         #
# function returns 3 values                         #
# best_Nm       <= optimal Nm for current param set #
# elliptic_lut3 <= reduced lut3 param for total     #
# elliptic equation system                          #
# elliptic_Nm <= Nm for reduced N                   #
# elliptic_optmat <= reduced optmat for elliptic eq #
# if L > 0 or M > 0, we have to reduce optmat       #
#####################################################
proc zincblende_kp_8x8_optimal_param { material } {

    upvar $material        mat
    global material_tcl_elliptic_threshold

    # save as we will modify it
    set optmat_orig $mat(optical_matrix)
    set lut3_orig   $mat(lut3)

    # this flag is set if Ac < 0 || L > 0 || N > 0 || M > 0 (even if non-ellipticity ratio is o.k.)
    set mat(rescale_mandatory) 0

    # loop until you are elliptic
    set loop        1
    set first_loop  1
    while { $loop == 1 } {

        # calculate kp 8x8 luttinger params
        zincblende_kp_8x8_luttinger lut8x8 mat
        # calculate kane params from 8x8 lut
        luttinger_to_kane kane $lut8x8(1) $lut8x8(2) $lut8x8(3)
        # get best Nm
        set best_Nm  [zincblende_kp_8x8_optimal_Nm $kane(L) $kane(M) $kane(N)]

        # check nonellipticity
        set ratio    [zincblende_kp_non_ellipticity_ratio $kane(L) $kane(M) $kane(N) $best_Nm]

        # save this if we are in the first loop
        if { $first_loop == 1 } {
            set mat(best_Nm_8x8) $best_Nm
            set mat(nonellipticity_original) $ratio
            set first_loop 0
            # BUT, if
            if {$lut8x8(Ac) < 0 || $kane(L) > 0 || $kane(M) > 0 || $kane(N) > 0} {
                set mat(rescale_mandatory) 1
            }
        }
        # quit if ratio is o.k
        if {$lut8x8(Ac) > 0.1 && $kane(L) < -0.1 && $kane(M) < -0.1 && $kane(N) < -0.1 && $ratio <= $material_tcl_elliptic_threshold } {
            set loop 0
        # modify material parameters
        } else {
            # first, bring these guys into their negative range
            if {$lut8x8(Ac) < 0 || $kane(L) > -0.1 || $kane(M) > -0.1 || $kane(N) > -0.1} {
                set mat(optical_matrix) [expr {$mat(optical_matrix) - 0.5}]
            # then reduce lut3 (~= N) by 5%
            } else {
                # and, reduce Ep slightly more, to decrease the required reduction of lut3
                set mat(optical_matrix) [expr {$mat(optical_matrix) - 0.5}]
                # Ep seems to be random, lut3 is non-controversial, so we should try to
                # keep this parameters
                set mat(lut3) [expr {$mat(lut3) * 0.99}]
                # prevent infinite loop
                if { $mat(lut3) < 0.01 } {
                    set loop 0
                    puts " ** warning ** - could not create elliptic params ... "
                }
            }
        }
    }

    # store elliptic parameters and reset lut3 and optmat to original
    set mat(elliptic_optical_matrix) $mat(optical_matrix)
    set mat(elliptic_lut3)           $mat(lut3)
    set mat(elliptic_Nm)             $best_Nm
    set mat(optical_matrix)          $optmat_orig
    set mat(lut3)                    $lut3_orig

	# also store eigenvalues
	set Np [expr {$kane(N) - $best_Nm}]
	set mat(elliptic_ev1)			 [expr {$kane(M) - $best_Nm}]
	set mat(elliptic_ev2)			 [expr {$kane(M) + $best_Nm}]
	set mat(elliptic_ev3)			 [expr {$kane(L) - $Np}]
	set mat(elliptic_ev4)			 [expr {$kane(L) + 2.0 * $Np}]

}

#####################################################
# determine best Nm for given set of kane params    #
#####################################################
proc zincblende_kp_8x8_optimal_Nm { L M N } {

    # foreman
    set Nm_final    [expr {$M - 1}]
    set best_ratio  0.0
    set num         100
    set Nm          [expr {- abs($M) - 5.0}]
    set end         [expr {abs($M) + 5.0}]
    set dN          [expr {($end - $Nm) / ($num - 1.0)}]
    set ii          0
    # check for best ratio
    while { $Nm < $end } {
        set ratio [zincblende_kp_non_ellipticity_ratio $L $M $N $Nm]
        if { $ii == 0 || $best_ratio > $ratio } {
            set best_ratio $ratio
            set Nm_final $Nm
            incr ii
        }
        set Nm [expr {$Nm + $dN}]
    }
    return $Nm_final

}


#####################################################
# create materials of geometry object and store     #
# them in the material database                     #
#####################################################
proc load_materials_from_tcl_db {geometry material_database temperature} {

    global material_tcl_dump
    global material_tcl_use_elliptic_params
    global material_tcl_use_optimal_Nm
    global material_tcl_elliptic_threshold
    global constants_hbar
    global constants_m0
    global LOG_INFO

    set material_names [list]
    # get all distinct materials
    for {set ii 0} { $ii < [$geometry get_num_regions] } {incr ii} {
        set matname [string2char [[$geometry get_region $ii] get_material_name]]
        # do we know this guy? and it shouldnt be "Contact", cause we ignore Contacts
        if { [exists_in_list $matname $material_names] == 0 && $matname != "Contact"} {
            lappend material_names $matname
        }
    }
    # now get all materials
    set top_material 0
    while { [lindex $material_names $top_material] == "Gas" && $top_material < [llength $material_names]} {
        incr top_material
    }
    for {set ii 0} { $ii < [llength $material_names] } {incr ii} {
        get_material_object mat [lindex $material_names $ii]
        if { $mat(crystal) == "Zincblende" } {
            zincblende_kp_8x8_optimal_param mat
        }
        # collect all material properties
        foreach p [array names mat] {
            set db($ii,$p) $mat($p)
        }
        # check for material with top valence band offset
        if {[lindex $material_names $ii] != "Gas" && $db($top_material,vbo) < $db($ii,vbo) } {
            set top_material $ii
        }
    }

    # set refernence vbo
    set vbo_ref $db($top_material,vbo)

    # calculate band edges and insert objects
    for {set ii 0} { $ii < [llength $material_names] } {incr ii} {
        set matname [lindex $material_names $ii]
        if { $matname != "Gas" } {
            # find relative vb offset
            set vbo_diff  [expr {$db($top_material,vbo) - $db($ii,vbo)}]
            set egT0_diff [expr {$db($ii,bandgapT0) - $db($top_material,bandgapT0)}]
            if { $egT0_diff >= 0 && $vbo_diff > 0} {
                set ratio    [expr {$vbo_diff / $egT0_diff} ]
                if { $ratio <= 1.0 } {
                    set Tgap_diff [expr {$db($ii,bandgapT0) - $db($ii,bandgap)}]
                    set db($ii,vb_edge) [expr {- $vbo_diff + $Tgap_diff * $ratio}]
                    set db($ii,cb_edge) [expr {$db($ii,vb_edge) + $db($ii,bandgap)}]
                } else {
                    logger_emit "valence band offset difference is bigger than bandgap difference!" $LOG_WARN
                    set db($ii,cb_edge) [expr {$db($ii,bandgap) - $vbo_diff}]
                    set db($ii,vb_edge) -$vbo_diff
                }
            } else {
                if { $egT0_diff < 0.0 || $vbo_diff != 0.0 } {
                    logger_emit "bandgap of material $ii is smaller than of top valence band edge material" $LOG_WARN
                }
                set db($ii,cb_edge) [expr {$db($ii,bandgap) - $vbo_diff}]
                set db($ii,vb_edge) -$vbo_diff
            }
            # adjust kp 8x8 params if requested
            if { $material_tcl_use_elliptic_params == 1 && $db($ii,crystal) == "Zincblende" } {
                # set elliptic parameters (if existing, else best_Nm is good enough)
                if { $db($ii,nonellipticity_original) > $material_tcl_elliptic_threshold || $db($ii,rescale_mandatory) == 1} {
                    set lut3_prct [format "%2.1f" [expr {(1.0 - ($db($ii,elliptic_lut3) / $db($ii,lut3))) * 100.0}]]
                    set ep_prct   [format "%2.1f" [expr {(1.0 - ($db($ii,elliptic_optical_matrix) / $db($ii,optical_matrix))) * 100.0}]]
                    logger_emit "$matname kp 8x8 ellipticity parameter resetting\n  optical_matrix: $db($ii,optical_matrix) to $db($ii,elliptic_optical_matrix) ($ep_prct %)\n  luttinger_par_3: $db($ii,lut3) to $db($ii,elliptic_lut3) ($lut3_prct %)\n  eigenvalues: $db($ii,elliptic_ev1), $db($ii,elliptic_ev2), $db($ii,elliptic_ev3), $db($ii,elliptic_ev4)" $LOG_INFO
                    set db($ii,best_Nm_8x8)    $db($ii,elliptic_Nm)
                    set db($ii,lut3)           $db($ii,elliptic_lut3)
                    set db($ii,optical_matrix) $db($ii,elliptic_optical_matrix)
                }
            }
            if { $material_tcl_use_optimal_Nm == 1 && $db($ii,crystal) == "Zincblende" } {
                # don't update optmat and lut3, just use best Nm
                set db($ii,Nm_8x8) $db($ii,best_Nm_8x8)
            }

            # Nm must be multplied by hbar^2/2m0
            # displayed kane params are in units of hbar. therefore switched
            # input of Nm into units of hbar^2/2m0. therefore, such a scaling
            # is not required here anymore
#            if { [info exists db($ii,Nm_8x8)] } {
#                set db($ii,Nm_8x8) [expr {$constants_hbar * $constants_hbar / (2.0 * $constants_m0) * $db($ii,Nm_8x8)}]
#            }
        }
        set obj [get_tdkp_material_object $ii db]
        $material_database add_material [lindex $material_names $ii] $obj

        if { $material_tcl_dump == 1 } {
            puts "-> Material $matname"
            $obj dump
        }
    }

}


#####################################################
# prepare material for given temperature            #
#####################################################
proc set_temperature_dependent_properties { material temperature } {

    upvar $material mat

    # create correct object
    if { $mat(crystal) == "Zincblende" } {
        set mat(bandgap)          [get_bandgap_at_T $mat(bandgapT0) $mat(alpha) $mat(beta) $temperature]
        set mat(lattice_constant) [expr {$mat(lattice_T300) + $mat(lattice_dT300) * ($temperature - 300.0)}]
        # preliminary set vb edge to vbo and cb edge to vb_edge + bandgap
        set mat(vb_edge)          $mat(vbo)
        set mat(cb_edge)          [expr {$mat(vbo) + $mat(bandgap)}]
        # set hole effective mass to heavy hole stuff
        set mat(mhh_z)            [expr {1.0 / ($mat(lut1) - 2.0 * $mat(lut2))}]
        set mat(mhh_tr)           [expr {1.0 / ($mat(lut1) + 2.0 * $mat(lut2))}]
        # electron effective mass
        set mat(me_cb_tr)         $mat(me_cb)
    } elseif { $mat(crystal) == "Wurtzite" } {
        set mat(bandgap)          [get_bandgap_at_T $mat(bandgapT0) $mat(alpha) $mat(beta) $temperature]
        # preliminary set vb edge to vbo and cb edge to vb_edge + bandgap
        set mat(vb_edge)          $mat(vbo)
        set mat(cb_edge)          [expr {$mat(vbo) + $mat(bandgap)}]
        # transverse effective masses (assume quantization is c-plane, transverse a-plane)
        set mat(mcb)              $mat(mcb_c)
        set mat(mcb_tr)           $mat(mcb_a)

    } else {
        puts "-> ** error **: unknown crystal"
        exit
    }
}


#####################################################
# create alloy of two materials                     #
# mole_frac = 1 -> mat1                             #
# mole_frac = 0 -> mat2                             #
# instead of luttinger parameters, we interpolate   #
# kane parameters linearly and take the luttinger   #
# params therefrom                                  #
#####################################################
proc create_ternary_alloy { target mole_frac mat1r mat2r bowing12r } {

    upvar $target    tern
    upvar $mat1r     mat1
    upvar $mat2r     mat2
    upvar $bowing12r bowing12

    set tern(mole_fraction) $mole_frac

    # overwrite bandgapT0 (x-dependent bowing)
    if { [info exists bowing12(bandgapT0_x)] } {
        set old_bandgapT0 $bowing12(bandgapT0)
        set bowing12(bandgapT0) [expr { $bowing12(bandgapT0) + $bowing12(bandgapT0_x) * $mole_frac }]
    }
    # overwrite bandgapT0 (x-dependent bowing)
    if { [info exists bowing12(bandgap_x)] } {
        set old_bandgap $bowing12(bandgap)
        set bowing12(bandgap) [expr { $bowing12(bandgap) + $bowing12(bandgap_x) * $mole_frac }]
    }

    # calculate kane params for zincblende
    if { $mat1(crystal) == "Zincblende" } {
        luttinger_to_kane kane1 $mat1(lut1) $mat1(lut2) $mat1(lut3)
        luttinger_to_kane kane2 $mat2(lut1) $mat2(lut2) $mat2(lut3)
        foreach p [array names kane1] {
            set mat1($p) $kane1($p)
            set mat2($p) $kane2($p)
        }
    }
    # first, standard treatment
    foreach p [array names mat1] {
        if { [info exists mat2($p)] } {
            set C 0.0
            if { [info exists bowing12($p)] } {
                set C $bowing12($p)
            }
            # only interpolate numeric values
            if { $p == "crystal" } {
                set tern($p) $mat1($p)
            } elseif { $p == "identifier" } {
                set tern($p) [format "%s%s%1.3f" $mat1($p) $mat2($p) $mole_frac]
            } else {
                set tern($p) [expr {$mat1($p) * (1.0 - $mole_frac) + $mat2($p) * $mole_frac - $mole_frac * (1.0 - $mole_frac) * $C } ]
            }
        }
    }

    # set linearly interpolated kane params to luttinger params
    if { $mat1(crystal) == "Zincblende" } {
        kane_to_luttinger lut $tern(L) $tern(M) $tern(N)
        foreach p [array names lut] {
            set tern($p) $lut($p)
        }
    }

    # set bandgapT0 back to its initial value (before bowing)
    if { [info exists old_bandgapT0] } {
        set bowing12(bandgapT0) $old_bandgapT0
    }
    # set bandgapT0 back to its initial value (before bowing)
    if { [info exists old_bandgap] } {
        set bowing12(bandgap) $old_bandgap
    }

}


#####################################################
# calculate kane parameters from luttinger params   #
# actually, its 2m0 / hbar * kane params            #
#####################################################
proc luttinger_to_kane { target lut1 lut2 lut3 } {

    upvar $target k

    set k(L)   [expr {- ($lut1 + 4.0 * $lut2) } ]
    set k(M)   [expr {- ($lut1 - 2.0 * $lut2) } ]
    set k(N)   [expr {- (6.0 * $lut3) } ]

}

#####################################################
# calculate luttinger parameters from kane          #
#####################################################
proc kane_to_luttinger { target L M N } {

    upvar $target k

    set k(lut1) [expr {1.0 / 3.0 * (- $L - 2.0 * $M) }]
    set k(lut2) [expr {1.0 / 6.0 * (- $L + $M) }]
    set k(lut3) [expr {- 1.0 / 6.0 * $N}]

}

#####################################################
# build material object from tcl associative array  #
#####################################################
proc get_tdkp_material_object { ii mat_ref } {

    upvar $mat_ref mat

    # create correct object
    if { $mat($ii,crystal) == "Zincblende" } {
        return [get_tdkp_zincblende_object $ii mat]
    } elseif { $mat($ii,crystal) == "Wurtzite" } {
        return [get_tdkp_wurtzite_object $ii mat]
    } elseif { $mat($ii,crystal) == "Gas" } {
        return [get_tdkp_gas_object $ii mat]
    } else {
        puts "-> ** error **: unknown crystal"
        exit
    }

}

#####################################################
# build zincblende material object from tcl         #
# associative array                                 #
#####################################################
proc get_tdkp_zincblende_object { ii mat_ref } {

    upvar $mat_ref mat

    # property container
    set obj [MaterialPropertyContainer_create_zincblende_material]

    # available params (to set)
    set par(bandgap)            "bandgap"
    set par(cb_edge)            "conduction_band_edge"
    set par(vb_edge)            "valence_band_edge"
    set par(spin_orbit)         "spin_orbit_splitting"
    set par(me_cb)              "electron_effective_mass"
    set par(me_cb_tr)           "electron_effective_mass_transverse"
    set par(lut1)               "luttinger_parameter_1"
    set par(lut2)               "luttinger_parameter_2"
    set par(lut3)               "luttinger_parameter_3"
    set par(optical_matrix)     "optical_matrix_element"
    set par(strain_ac)          "strain_potential_ac"
    set par(strain_av)          "strain_potential_av"
    set par(strain_b)           "strain_potential_b"
    set par(strain_d)           "strain_potential_d"
    set par(lattice_constant)   "lattice_constant"
    set par(elastic_c11)        "elastic_constant_C11"
    set par(elastic_c12)        "elastic_constant_C12"
    set par(elastic_c44)        "elastic_constant_C44"
    set par(piezo_e14)          "piezo_coefficient_e14"
    set par(Nm_8x8)             "foremans_kane_parameter_N_minus"
    set par(mhh_z)              "hole_effective_mass"
    set par(mhh_tr)             "hole_effective_mass_transverse"
    set par(el_permittivity)    "electrostatic_permittivity"

    # check bandgap
    if { [expr {abs($mat($ii,cb_edge) - ($mat($ii,vb_edge) + $mat($ii,bandgap)))}]  > 1.0e-6 } {
        puts "-> ** error ** - bandgap is not equal to cb_edge - vb_edge "
        exit
    }

    # set all available params
    foreach p [array names par] {
        if { [info exists mat($ii,$p)] } {
            $obj set $par($p) $mat($ii,$p)
        }
    }

    # test object
    if { ![$obj valid] } {
        puts "-> ** error ** - material object seems to be invalid "
        exit
    }

    return $obj

}

#####################################################
# build wurtzite material object from tcl           #
# associative array                                 #
#####################################################
proc get_tdkp_wurtzite_object { ii mat_ref } {

    upvar $mat_ref mat

    # property container
    set obj [MaterialPropertyContainer_create_wurtzite_material]

    # available params (to set)
    set par(spin_orbit_delta_1)     "spin_orbit_split_delta_1"
    set par(spin_orbit_delta_2)     "spin_orbit_split_delta_2"
    set par(spin_orbit_delta_3)     "spin_orbit_split_delta_3"
    set par(mcb_c)                  "electron_effective_mass_c_plane"
    set par(mcb_a)                  "electron_effective_mass_a_plane"
    set par(mcb)                    "electron_effective_mass"
    set par(mcb_tr)                 "electron_effective_mass_transverse"
    set par(cb_edge)                "conduction_band_edge"
    set par(vb_edge)                "valence_band_edge"
    set par(A1)                     "hole_effective_mass_A1"
    set par(A2)                     "hole_effective_mass_A2"
    set par(A3)                     "hole_effective_mass_A3"
    set par(A4)                     "hole_effective_mass_A4"
    set par(A5)                     "hole_effective_mass_A5"
    set par(A6)                     "hole_effective_mass_A6"
    set par(A7)                     "hole_spin_orbit_A7"
    set par(strain_D1)              "strain_potential_D1"
    set par(strain_D2)              "strain_potential_D2"
    set par(strain_D3)              "strain_potential_D3"
    set par(strain_D4)              "strain_potential_D4"
    set par(strain_D5)              "strain_potential_D5"
    set par(strain_D6)              "strain_potential_D6"
    set par(strain_ac_cc)           "strain_potential_ac_cc"
    set par(strain_ac_aa)           "strain_potential_ac_aa"
    set par(strain_av_cc)           "strain_potential_av_cc"
    set par(strain_av_aa)           "strain_potential_av_aa"
    set par(elastic_C11)            "elastic_constant_C11"
    set par(elastic_C12)            "elastic_constant_C12"
    set par(elastic_C13)            "elastic_constant_C13"
    set par(elastic_C33)            "elastic_constant_C33"
    set par(elastic_C44)            "elastic_constant_C44"
    set par(lattice_a_T300)         "lattice_constant_a"
    set par(lattice_c_T300)         "lattice_constant_c"
    # these guys here are optional
    set par(matrix_element_P1)      "optical_momentum_matrix_element_P1"
    set par(matrix_element_P2)      "optical_momentum_matrix_element_P2"

    set par(piezo_e31)              "piezo_coefficient_e31"
    set par(piezo_e33)              "piezo_coefficient_e33"
    set par(piezo_e15)              "piezo_coefficient_e15"
    set par(polarization_Psp)       "spontaneous_polarization_Psp_zz"

    set par(el_permittivity)        "electrostatic_permittivity"
    set par(el_permittivity_xx)     "electrostatic_permittivity_xx"
    set par(el_permittivity_yy)     "electrostatic_permittivity_yy"
    set par(el_permittivity_zz)     "electrostatic_permittivity_zz"

    # check bandgap
    if { [expr {abs($mat($ii,cb_edge) - ($mat($ii,vb_edge) + $mat($ii,bandgap)))}]  > 1.0e-6 } {
        puts "-> ** error ** - bandgap is not equal to cb_edge - vb_edge "
        exit
    }

    # set all available params
    foreach p [array names par] {
        if { [info exists mat($ii,$p)] } {
            $obj set $par($p) $mat($ii,$p)
        }
    }

    # test object
    if { ![$obj valid] } {
        puts "-> ** error ** - material object seems to be invalid "
        exit
    }

    return $obj

}

#####################################################
# build gas material object from tcl                #
# associative array                                 #
#####################################################
proc get_tdkp_gas_object { ii mat_ref } {

    upvar $mat_ref mat

    # property container
    set obj [MaterialPropertyContainer_create_gas_material]

    # available params (to set)
    set par(permittivity_e)         "electrostatic_permittivity"
    set par(elastic_C11)            "elastic_constant_C11"
    set par(elastic_C12)            "elastic_constant_C12"
    set par(elastic_C13)            "elastic_constant_C13"
    set par(elastic_C33)            "elastic_constant_C33"
    set par(elastic_C44)            "elastic_constant_C44"


    # set all available params
    foreach p [array names par] {
        if { [info exists mat($ii,$p)] } {
            $obj set $par($p) $mat($ii,$p)
        }
    }

    # test object
    if { ![$obj valid] } {
        puts "-> ** error ** - material object seems to be invalid "
        exit
    }

    return $obj

}


#####################################################
# Varshni form of bandgap on temperature dependece: #
# Eg(T) = Eg(T=0) - alpha * T^2 / (T + beta)        #
# see eq. (2.13) in [1]                             #
#####################################################
proc get_bandgap_at_T { bandgap0 alpha beta temperature } {
    if { $temperature < 0.001 } {
        return $bandgap0
    } else {
        return [expr {$bandgap0 - $alpha * $temperature * $temperature / ($temperature + $beta) } ]
    }
}


#####################################################
# return GaAs parameters                            #
# note: matrix element is smaller in order to       #
#       reduce non-ellipticity                      #
#####################################################
proc param_GaAs { variable } {

    upvar $variable mat

    set mat(identifier)         "GaAs"
    set mat(crystal)            "Zincblende"

    # bandedges, splittings and valence band offset
    set mat(bandgapT0)          1.519
    set mat(alpha)              0.540e-3
    set mat(beta)               204.0
    set mat(spin_orbit)         0.341
    # see fig. 10 in [1]
    set mat(vbo)                -0.80

    # effective mass parameterers
    set mat(me_cb)              0.067
    set mat(lut1)               6.98
    set mat(lut2)               2.06
    set mat(lut3)               2.93

    # problematic parameters ...
    # 28.8 eV would produce spurious solutions
    set mat(optical_matrix)     25.5

    # strain parameters (note in [1] a_gap = ac + av,
    # which means that av switched sign)
    set mat(strain_ac)          -7.17
    set mat(strain_av)          1.16
    set mat(strain_b)           -2.0
    set mat(strain_d)           -4.8


    # elastic and lattice constants
    set mat(lattice_T300)       5.65325
    set mat(lattice_dT300)      3.88e-5
    set mat(elastic_c11)        1221.0
    set mat(elastic_c12)        566.0
    set mat(elastic_c44)        600.0

    # piezo coefficient [4] (in tdkp's units)
    set mat(piezo_e14)          [expr { -0.16 / (1.60219e-19) / (1.0e9 * 1.0e9)}]

    # electrostatic permittivity [4]
    set mat(el_permittivity)    12.90

}

#####################################################
# return AlAs parameters                            #
#####################################################
proc param_AlAs { variable } {

    upvar $variable mat

    set mat(identifier)         "AlAs"
    set mat(crystal)            "Zincblende"

    # bandedges, splittings and valence band offset
    set mat(bandgapT0)          3.099
    set mat(alpha)              0.885e-3
    set mat(beta)               530.0
    set mat(spin_orbit)         0.28
    # see fig. 10 in [1] (valence band offset)
    set mat(vbo)                -1.33

    # effective mass parameterers
    set mat(me_cb)              0.15
    set mat(lut1)               3.76
    set mat(lut2)               0.82
    set mat(lut3)               1.42

    # problematic parameters ...
    set mat(optical_matrix)     21.1

    # strain parameters (note in [1] a_gap = ac + av,
    # which means that av switched sign (because a_gap = ac - av here))
    set mat(strain_ac)          -5.64
    set mat(strain_av)          2.47
    set mat(strain_b)           -2.3
    set mat(strain_d)           -3.4

    # elastic and lattice constants
    # these guys are used relative to each other, so the units just must be
    # consistent here ...
    set mat(lattice_T300)       5.6611
    set mat(lattice_dT300)      2.9e-5
    set mat(elastic_c11)        1250.0
    set mat(elastic_c12)        534.0
    set mat(elastic_c44)        542.0

    # electrostatic permittivity [4]
    set mat(el_permittivity)     10.06


}

#####################################################
# return InAs parameters                            #
#####################################################
proc param_InAs { variable } {

    upvar $variable mat

    set mat(identifier)         "InAs"
    set mat(crystal)            "Zincblende"

    # bandedges, splittings and valence band offset
    set mat(bandgapT0)          0.417
    set mat(alpha)              0.27e-3
    set mat(beta)               93.0
    set mat(spin_orbit)         0.39
    # see fig. 10 in [1] (valence band offset)
    set mat(vbo)                -0.59

    # effective mass parameterers
    # modified cb mass (using 0.0223, see chuang)
    set mat(me_cb)              0.0223
    set mat(lut1)               20.0
    set mat(lut2)               8.5
    set mat(lut3)               9.2


    # problematic parameters ...
    # modified (set to 19 to have elliptic equation system)
    set mat(optical_matrix)     21.0


    # strain parameters (note in [1] a_gap = ac + av,
    # which means that av switched sign (because a_gap = ac - av here))
    set mat(strain_ac)          -5.08
    set mat(strain_av)           1.00
    set mat(strain_b)           -1.8
    set mat(strain_d)           -3.6

    # elastic and lattice constants
    # these guys are used relative to each other, so the units just must be
    # consistent here ...
    set mat(lattice_T300)       6.0583
    set mat(lattice_dT300)      2.74e-5
    set mat(elastic_c11)        832.9
    set mat(elastic_c12)        452.6
    set mat(elastic_c44)        395.9

    # electrostatic permittivity [4]
    set mat(el_permittivity)    14.3

    # piezo coefficient [4] (in tdkp's units)
    set mat(piezo_e14)          [expr { 0.045 / (1.60219e-19) / (1.0e9 * 1.0e9)}]

}

#####################################################
# return InP parameters                            #
#####################################################
proc param_InP { variable } {

    upvar $variable mat

    set mat(identifier)         "InP"
    set mat(crystal)            "Zincblende"

    # bandedges, splittings and valence band offset
    set mat(bandgapT0)          1.4236
    set mat(alpha)              0.363e-3
    set mat(beta)               162.0
    set mat(spin_orbit)         0.108
    # see fig. 10 in [1] (valence band offset)
    set mat(vbo)                -0.94

    # effective mass parameterers
    set mat(me_cb)              0.0795
    set mat(lut1)               5.08
    set mat(lut2)               1.60
    set mat(lut3)               2.10

    # problematic parameters ...
    set mat(optical_matrix)     20.7

    # strain parameters (note in [1] a_gap = ac + av,
    # which means that av switched sign (because a_gap = ac - av here))
    set mat(strain_ac)          -6.0
    set mat(strain_av)           0.6
    set mat(strain_b)           -2.0
    set mat(strain_d)           -5.0

    # elastic and lattice constants
    # these guys are used relative to each other, so the units just must be
    # consistent here ...
    set mat(lattice_T300)       5.8697
    set mat(lattice_dT300)      2.79e-5
    set mat(elastic_c11)        1011.0
    set mat(elastic_c12)        561.0
    set mat(elastic_c44)        456.0

    # electrostatic permittivity [4]
#    set mat(el_permittivity)     ????


}


#####################################################
# bowing parameters for ternary alloy such that     #
# Al(x)Ga(1-x)As = x*Al + (1-x) Ga - x(1-x)*bowing  #
#####################################################
proc bowing_AlGaAs { variable } {

    upvar $variable bowing

    set bowing(bandgapT0)       -0.127
    set bowing(bandgap)         -0.127
    # x dependent bowing param (see table XII [1])
    set bowing(bandgapT0_x)     1.310
    set bowing(bandgap_x)       1.310

}

#####################################################
# bowing parameters for ternary alloy such that     #
# In(x)Ga(1-x)As = x*In + (1-x) Ga - x(1-x)*bowing  #
#####################################################
proc bowing_InGaAs { variable } {

    upvar $variable bowing

    set bowing(bandgapT0)       0.477
    set bowing(bandgap)         0.477
    set bowing(spin_orbit)      0.15
    set bowing(me_cb)           0.0091
    set bowing(optical_matrix)  -1.48
    set bowing(strain_ac)       2.61
    set bowing(vbo)             -0.38

}



# ---------------------------------------------------------------
#
# N I T R I D E S
#
# ---------------------------------------------------------------

#####################################################
# return GaN parameter [3]                          #
# and old ones from [2] depending on switch         #
#####################################################
proc param_GaN { variable } {

    upvar $variable mat
    global nitrides_use_vurgaftman_2007

    if { $nitrides_use_vurgaftman_2007 == 1 } {

        set mat(identifier)         "GaN"
        set mat(crystal)            "Wurtzite"

        # bandgaps
        set mat(bandgapT0)          3.510
        set mat(alpha)              0.914e-3
        set mat(beta)               825.0
        # delta_cr
        set mat(spin_orbit_delta_1) 0.010
        # delta 2 = delta 3 = delta_so / 3
        set mat(spin_orbit_delta_2) 0.006
        set mat(spin_orbit_delta_3) 0.006

        #################################################
        # according to marco tomamichel,                #
        # I. Vurgaftman and J. Meyer, unpublished (2006)#
        # give 0.8 as cb offset for InN / GaN           #
        # so we use vbo = 0 for InN and 0.8 offset      #
        # at T = 300 K to get vbo of GaN                #
        #################################################
        set mat(vbo)                -0.564

        # effective masses
        set mat(mcb_c)               0.20
        set mat(mcb_a)               0.21
        set mat(A1)                 -7.21
        set mat(A2)                 -0.44
        set mat(A3)                  6.68
        set mat(A4)                 -3.46
        set mat(A5)                 -3.40
        set mat(A6)                 -4.90
        set mat(A7)                 0.00937

        # strain potentials
        set mat(strain_D1)          -3.6
        set mat(strain_D2)           1.7
        set mat(strain_D3)           5.2
        set mat(strain_D4)          -2.7
        set mat(strain_D5)          -2.8
        set mat(strain_D6)          -4.3
	# very careful, credits go to fabio sacconi for pointing this out:
	# vurgaftman can be misunderstood. he names a1/a2 the interband deformation potentials
	# but the ones he quotes don't contain the crystal field splitting.
	# so in fact
	# 1) the real interband deformation potential is a1 - D3, a2 - D4 (where D3/D4 is the crystal field split)
        # 2) the hydrostatic response of wz in the vb is (D1 + D3) to ezz and (D2 + D4) to exx/eyy
        # 3) bandgap response is ag = ac - av, so ac = ag + av
        # 4) as a result, we have ac = a1 - D3 + (D1 + D3) = a1 + D1 and ac2 = a2 + D2
        # 5) check the literature, these params are not really known. so this here is nothing
        # more than an educated guess
        set a1 -7.1
        set mat(strain_ac_cc)       [expr {$a1 + $mat(strain_D1)}]
        set a2 -9.9
        set mat(strain_ac_aa)       [expr {$a2 + $mat(strain_D2)}]

        # approximate av_aa and av_cc
        set mat(strain_av_aa)       [expr {$mat(strain_D2) + $mat(strain_D4)}]
        set mat(strain_av_cc)       [expr {$mat(strain_D1) + $mat(strain_D3)}]


        # GPA (is unitless in the resulting equation anyway ...)
        set mat(elastic_C11)        390.0
        set mat(elastic_C12)        145.0
        set mat(elastic_C13)        106.0
        set mat(elastic_C33)        398.0
        set mat(elastic_C44)        105.0

        set mat(lattice_a_T300)     3.189
        set mat(lattice_c_T300)     5.185
        # piezo coefficents in electron charge / nm^2
        # (C / m^2 * (1/(1.6e-19)) * 1 / (1.0e9^2 nm^2)) = 6.241457
        # taken from bernadrini in piprek, nitride semiconductor devices
        set mat(piezo_e31)          -2.1096
        set mat(piezo_e33)           4.163051
        set mat(piezo_e15)          -1.04232331
        set mat(polarization_Psp)   [expr { -0.034 / (1.60219e-19) / (1.0e9 * 1.0e9)}]

        # directional permittivity [4]
        set mat(el_permittivity_xx) 9.6
        set mat(el_permittivity_yy) 9.6
        # c axis
        set mat(el_permittivity_zz) 10.6
        # average el permittivity
        set mat(el_permittivity)    [expr {($mat(el_permittivity_xx) + $mat(el_permittivity_yy) + $mat(el_permittivity_zz)) / 3.0}]

    } else {
        ##############################################
        # vurgaftman 2003 [2]                        #
        ##############################################
        set mat(identifier)         "GaN"
        set mat(crystal)            "Wurtzite"

        # bandgaps
        set mat(bandgapT0)          3.510
        set mat(alpha)              0.909e-3
        set mat(beta)               830.0
        # delta_cr
        set mat(spin_orbit_delta_1) 0.010
        # delta 2 = delta 3 = delta_so / 3
        set mat(spin_orbit_delta_2) 0.0056
        set mat(spin_orbit_delta_3) 0.0056

        #################################################
        # according to marco tomamichel,                #
        # I. Vurgaftman and J. Meyer, unpublished (2006)#
        # give 0.8 as cb offset for InN / GaN           #
        # so we use vbo = 0 for InN and 0.8 offset      #
        # at T = 300 K to get vbo of GaN                #
        #################################################
        set mat(vbo)                -0.546

        # effective masses
        set mat(mcb_c)              0.20
        set mat(mcb_a)              0.20
        set mat(A1)                 -7.21
        set mat(A2)                 -0.44
        set mat(A3)                 6.68
        set mat(A4)                 -3.46
        set mat(A5)                 -3.40
        set mat(A6)                 -4.90
        # attention, A7 is eV / Angstrom, therefore / 10 for eV / nanometer!
        set mat(A7)                 0.00937

        # strain potentials
        set mat(strain_D1)          -3.7
        set mat(strain_D2)           4.5
        set mat(strain_D3)           8.2
        set mat(strain_D4)          -4.1
        set mat(strain_D5)          -4.0
        set mat(strain_D6)          -5.5
        # a1 + D1
        set a1 -4.9
        set mat(strain_ac_cc)       [expr {$a1 + $mat(strain_D1)}]
        # a2 + D2 
        set a2 -11.3
        set mat(strain_ac_aa)       [expr {$a2 + $mat(strain_D2)}]

        # approximate av_aa and av_cc
        set mat(strain_av_aa)       [expr {$mat(strain_D2) + $mat(strain_D4)}]
        set mat(strain_av_cc)       [expr {$mat(strain_D1) + $mat(strain_D3)}]

        # GPA (is unitless anyway ...)
        set mat(elastic_C11)        390.0
        set mat(elastic_C12)        145.0
        set mat(elastic_C13)        106.0
        set mat(elastic_C33)        398.0
        set mat(elastic_C44)        105.0

        set mat(lattice_a_T300)     3.189
        set mat(lattice_c_T300)     5.185

    }

}

#####################################################
# return InN parameter [3]                          #
#####################################################
proc param_InN { variable } {

    upvar $variable mat

    global nitrides_use_vurgaftman_2007

    if { $nitrides_use_vurgaftman_2007 == 1 } {

        set mat(identifier)         "InN"
        set mat(crystal)            "Wurtzite"

        # bandgaps
        set mat(bandgapT0)          0.69
        set mat(alpha)              4.14e-4
        set mat(beta)               154.0
        # delta_cr
        set mat(spin_orbit_delta_1) 0.024
        # delta 2 = delta 3 = delta_so / 3
        set mat(spin_orbit_delta_2) 0.0017
        set mat(spin_orbit_delta_3) 0.0017

        #################################################
        # according to marco tomamichel,                #
        # I. Vurgaftman and J. Meyer, unpublished (2006)#
        # give 0.8 as cb offset for InN / GaN           #
        # so we use vbo = 0 for InN and 0.8 offset      #
        # at T = 300 K to get vbo of GaN                #
        #################################################
        set mat(vbo)                0.0

        # effective masses
        set mat(mcb_c)              0.07
        set mat(mcb_a)              0.07
        set mat(A1)                 -8.21
        set mat(A2)                 -0.68
        set mat(A3)                  7.57
        set mat(A4)                 -5.23
        set mat(A5)                 -5.11
        set mat(A6)                 -5.96
        set mat(A7)                  0.0

        # strain potentials
        set mat(strain_D1)          -3.6
        set mat(strain_D2)           1.7
        set mat(strain_D3)           5.2
        set mat(strain_D4)          -2.7
        set mat(strain_D5)          -2.8
        set mat(strain_D6)          -4.3
        # a1 + D1 
        set a1 -4.2
        set mat(strain_ac_cc)       [expr {$a1 + $mat(strain_D1)}]
        # a2 + (D2 + D4)
        set a2 -4.2
        set mat(strain_ac_aa)       [expr {$a2 + $mat(strain_D2)}]

        # approximate av_aa and av_cc
        set mat(strain_av_aa)       [expr {$mat(strain_D2) + $mat(strain_D4)}]
        set mat(strain_av_cc)       [expr {$mat(strain_D1) + $mat(strain_D3)}]

        # GPA (is unitless anyway ...)
        set mat(elastic_C11)        223.0
        set mat(elastic_C12)        115.0
        set mat(elastic_C13)        92.0
        set mat(elastic_C33)        224.0
        set mat(elastic_C44)        48.0

        set mat(lattice_a_T300)     3.545
        set mat(lattice_c_T300)     5.703

        # piezo coefficents in electron charge / nm^2
        # (C / m^2 * (1/(1.6e-19)) * 1 / (1.0e9^2 nm^2)) = 6.241457
        # taken from bernadrini in piprek, nitride semiconductor devices
        set mat(piezo_e31)          -2.57148
        set mat(piezo_e33)           5.0867
        set mat(piezo_e15)          -0.699
		set mat(polarization_Psp)   [expr { -0.042 / (1.60219e-19) / (1.0e9 * 1.0e9)}]

        # directional permittivity [4]
        set mat(el_permittivity_xx) 13.1
        set mat(el_permittivity_yy) 13.1
        # c axis
        set mat(el_permittivity_zz) 14.4
        # average el permittivity
        set mat(el_permittivity)    [expr {($mat(el_permittivity_xx) + $mat(el_permittivity_yy) + $mat(el_permittivity_zz)) / 3.0}]

    } else {

        ##############################################
        # vurgaftman 2003 [2]                        #
        ##############################################
        set mat(identifier)         "InN"
        set mat(crystal)            "Wurtzite"

        # bandgaps
        set mat(bandgapT0)          0.78
        set mat(alpha)              0.245e-3
        set mat(beta)               624.0
        # delta_cr
        set mat(spin_orbit_delta_1) 0.040
        # delta 2 = delta 3 = delta_so / 3
        set mat(spin_orbit_delta_2) 0.00167
        set mat(spin_orbit_delta_3) 0.00167

        #################################################
        # according to marco tomamichel,                #
        # I. Vurgaftman and J. Meyer, unpublished (2006)#
        # give 0.8 as cb offset for InN / GaN           #
        # so we use vbo = 0 for InN and 0.8 offset      #
        # at T = 300 K to get vbo of GaN                #
        #################################################
        set mat(vbo)                0.0

        # effective masses
        set mat(mcb_c)              0.07
        set mat(mcb_a)              0.07
        set mat(A1)                 -8.21
        set mat(A2)                 -0.68
        set mat(A3)                 7.57
        set mat(A4)                 -5.23
        set mat(A5)                 -5.11
        set mat(A6)                 -5.96
        set mat(A7)                 0.0

        # strain potentials
        set mat(strain_D1)          -3.7
        set mat(strain_D2)          4.5
        set mat(strain_D3)          8.2
        set mat(strain_D4)          -4.1
        set mat(strain_D5)          -4.0
        set mat(strain_D6)          -5.5
        # a1 + D1 
        set a1  -3.5
        set mat(strain_ac_cc)       [expr {$a1 + $mat(strain_D1)}]
        # a2 + D2 
        set a2  -3.5
        set mat(strain_ac_aa)       [expr {$a2 + $mat(strain_D2)}]

        # approximate av_aa and av_cc
        set mat(strain_av_aa)       [expr {$mat(strain_D2) + $mat(strain_D4)}]
        set mat(strain_av_cc)       [expr {$mat(strain_D1) + $mat(strain_D3)}]

        # GPA (is unitless anyway ...)
        set mat(elastic_C11)        223.0
        set mat(elastic_C12)        115.0
        set mat(elastic_C13)        92.0
        set mat(elastic_C33)        224.0
        set mat(elastic_C44)        48.0

        set mat(lattice_a_T300)     3.545
        set mat(lattice_c_T300)     5.703

    }
}


#####################################################
# return AlN parameter [3]
#####################################################
proc param_AlN { variable } {

    global nitrides_use_vurgaftman_2007
    upvar $variable mat

    # warning, the values for vbo are arbitrary
    # fix it by yourself. i just needed effmass dispersions.
    puts "WARNING VBO IS BULLSHIT!!!!"
    puts "I REPEAT: WARNING VBO IS BULLSHIT!!!!"

    if { $nitrides_use_vurgaftman_2007 == 1 } {

        set mat(identifier)         "AlN"
        set mat(crystal)            "Wurtzite"

        # bandgaps
        set mat(bandgapT0)          6.1
        set mat(alpha)              26.3e-4
        set mat(beta)               2082.0
        # delta_cr
        set mat(spin_orbit_delta_1) -0.227
        # delta 2 = delta 3 = delta_so / 3
        set mat(spin_orbit_delta_2) 0.012
        set mat(spin_orbit_delta_3) 0.012

        #################################################
        # according to marco tomamichel,                #
        # I. Vurgaftman and J. Meyer, unpublished (2006)#
        # give 0.8 as cb offset for InN / GaN           #
        # so we use vbo = 0 for InN and 0.8 offset      #
        # at T = 300 K to get vbo of GaN                #
        #################################################
        set mat(vbo)                -0.746

        # effective masses
        set mat(mcb_c)               0.30
        set mat(mcb_a)               0.32
        set mat(A1)                 -3.86
        set mat(A2)                 -0.25
        set mat(A3)                  3.58
        set mat(A4)                 -1.32
        set mat(A5)                 -1.47
        set mat(A6)                 -1.64
        set mat(A7)                  0.0

        # strain potentials
        set mat(strain_D1)          -2.9
        set mat(strain_D2)           4.9
        set mat(strain_D3)           9.4
        set mat(strain_D4)          -4.0
        set mat(strain_D5)          -3.3
        set mat(strain_D6)          -2.7
        # a1 + D1 
        set a1 -3.4
        set mat(strain_ac_cc)       [expr {$a1 + $mat(strain_D1)}]
        # a2 + D2 
        set a2 -11.8
        set mat(strain_ac_aa)       [expr {$a2 + $mat(strain_D2)}]

        # approximate av_aa and av_cc
        set mat(strain_av_aa)       [expr {$mat(strain_D2) + $mat(strain_D4)}]
        set mat(strain_av_cc)       [expr {$mat(strain_D1) + $mat(strain_D3)}]

        # GPA (is unitless anyway ...)
        set mat(elastic_C11)        396.0
        set mat(elastic_C12)        137.0
        set mat(elastic_C13)        108.0
        set mat(elastic_C33)        373.0
        set mat(elastic_C44)        116.0

        set mat(lattice_a_T300)     3.112
        set mat(lattice_c_T300)     4.982

        # piezo coefficents in electron charge / nm^2
        # (C / m^2 * (1/(1.6e-19)) * 1 / (1.0e9^2 nm^2)) = 6.241457
        # taken from bernadrini in piprek, nitride semiconductor devices
        set mat(piezo_e31)          -3.32669
        set mat(piezo_e33)           9.3933
        set mat(piezo_e15)          -2.19057
        set mat(polarization_Psp)   [expr { -0.09 / (1.60219e-19) / (1.0e9 * 1.0e9)}]

        # directional permittivity [4]
        set mat(el_permittivity_xx) 8.3
        set mat(el_permittivity_yy) 8.3
        # c axis
        set mat(el_permittivity_zz) 8.9
        # average el permittivity
        set mat(el_permittivity)    [expr {($mat(el_permittivity_xx) + $mat(el_permittivity_yy) + $mat(el_permittivity_zz)) / 3.0}]

    } else {
        ##############################################
        # vurgaftman 2003 [2]                        #
        ##############################################
        set mat(identifier)         "AlN"
        set mat(crystal)            "Wurtzite"

        # bandgaps
        set mat(bandgapT0)          6.25
        set mat(alpha)              1.799e-3
        set mat(beta)               1462.0
        # delta_cr
        set mat(spin_orbit_delta_1) -0.169
        # delta 2 = delta 3 = delta_so / 3
        set mat(spin_orbit_delta_2) 0.00633333
        set mat(spin_orbit_delta_3) 0.00633333

        #################################################
        # according to marco tomamichel,                #
        # I. Vurgaftman and J. Meyer, unpublished (2006)#
        # give 0.8 as cb offset for InN / GaN           #
        # so we use vbo = 0 for InN and 0.8 offset      #
        # at T = 300 K to get vbo of GaN                #
        #################################################
        set mat(vbo)                -0.746

        # effective masses
        set mat(mcb_c)              0.30
        set mat(mcb_a)              0.32
        set mat(A1)                 -3.86
        set mat(A2)                 -0.25
        set mat(A3)                  3.58
        set mat(A4)                 -1.32
        set mat(A5)                 -1.47
        set mat(A6)                 -1.64
        set mat(A7)                 0.0

        # strain potentials
        set mat(strain_D1)          -17.1
        set mat(strain_D2)          7.9
        set mat(strain_D3)          8.8
        set mat(strain_D4)          -3.9
        set mat(strain_D5)          -3.4
        set mat(strain_D6)          -3.4
        # a1 + D1 
        set a1 -3.4
        set mat(strain_ac_cc)       [expr {$a1 + $mat(strain_D1)}]
        # a2 + D2 
        set a2 -11.8
        set mat(strain_ac_aa)       [expr {$a2 + $mat(strain_D2)}]

        # approximate av_aa and av_cc
        set mat(strain_av_aa)       [expr {$mat(strain_D2) + $mat(strain_D4)}]
        set mat(strain_av_cc)       [expr {$mat(strain_D1) + $mat(strain_D3)}]

        # GPA (is unitless anyway ...)
        set mat(elastic_C11)        396.0
        set mat(elastic_C12)        137.0
        set mat(elastic_C13)        108.0
        set mat(elastic_C33)        373.0
        set mat(elastic_C44)        116.0

        set mat(lattice_a_T300)     3.112
        set mat(lattice_c_T300)     4.982
    }

}

#####################################################
# bowing parameters for ternary alloy such that     #
# In(x)Ga(1-x)N = x*In + (1-x) Ga - x(1-x)*bowing   #
#####################################################
proc bowing_InGaN { variable } {


    global nitrides_use_vurgaftman_2007
    upvar $variable bowing
    if { $nitrides_use_vurgaftman_2007 == 1 } {
        set bowing(bandgapT0)        1.4
        set bowing(bandgap)          1.4
        set bowing(polarization_Psp) [expr { -0.037 / (1.60219e-19) / (1.0e9 * 1.0e9)}]
    } else {
        set bowing(bandgapT0)       1.4
        set bowing(bandgap)         1.4
    }

}

#####################################################
# bowing parameters for ternary alloy such that     #
# Al(x)Ga(1-x)N = x*Al + (1-x) Ga - x(1-x)*bowing   #
#####################################################
proc bowing_AlGaN { variable } {

    global nitrides_use_vurgaftman_2007
    upvar $variable bowing
    if { $nitrides_use_vurgaftman_2007 == 1 } {
        set bowing(bandgapT0)       0.8
        set bowing(bandgap)         0.8
        set bowing(polarization_Psp) [expr { -0.021 / (1.60219e-19) / (1.0e9 * 1.0e9)}]
    } else {
        set bowing(bandgapT0)       0.7
        set bowing(bandgap)         0.7
    }

}

#####################################################
# return Gas parameter
#####################################################
proc param_Gas { variable } {

    upvar $variable mat

    set mat(identifier)         "Gas"
    set mat(crystal)            "Gas"

    # available params (to set)
    set mat(permittivity_e)         1.0

    # set to very (arbitrary) small value
    set mat(elastic_C11)        0.000396
    set mat(elastic_C12)        0.000137
    set mat(elastic_C13)        0.000108
    set mat(elastic_C33)        0.000373
    set mat(elastic_C44)        0.000116

}
