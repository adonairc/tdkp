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


if {![info exists snapshot] && $::argc > 0} {
    set snapshot $::argv
}

if {![info exists snapshot]} {
    puts "please set the snapshot variable!"
    exit
}


###########################################
# read config                             #
###########################################
if { ![file exists "${snapshot}_config.txt"] } {
    puts "missing config file ${snapshot}_config.txt"
    quit
}
set config_str [exec cat "${snapshot}_config.txt" | grep "config=" | sed -e "s/^config=//"]
set gridfile   [exec basename [exec cat "${snapshot}_config.txt" | grep "gridname=" | sed -e "s/gridname=//"]]

set config [InterfaceConfiguration]
$config unserialize [char2string $config_str]
$config dump
puts "-> gridfile ${gridfile}"

set parser [InputParser]
set material_database [MaterialDatabase]
set reader [DFISEGridReader -args $gridfile]
set geometry [$parser read_geometry $reader]
$reader -delete
unset reader

###########################################
# check for snapshot files                #
###########################################
if { [file exists "${snapshot}_InterfaceImplementation.txt"] } {
    ###########################################
    # copy snapshot file and remove grid name #
    ###########################################
    set snap_interface_data [PropertyContainerdouble]
    $snap_interface_data read_data_from_file "${snapshot}_InterfaceImplementation.txt"
}

###########################################
# copy material files                     #
###########################################
set material_files [split [exec find . -name "${snapshot}*.mat"]]
foreach material $material_files {
    set target [exec echo "$material" | sed -e "s/${snapshot}_//"]
    exec cp -v $material $target
}


###########################################
# build infodesk                          #
###########################################
if { [file exists "${snapshot}_slc.txt"] } {
    set slc 1
    set infodesk [InformationDeskImplementation -args $material_database "${snapshot}_slc.txt"]
    set slc_data [PropertyContainerdouble]
    $slc_data read_data_from_file "${snapshot}_slc.txt"
} else {
    set slc 0
    set infodesk [InformationDeskImplementation -args $material_database]
}

###########################################
# build interface                         #
###########################################
set dim [$geometry get_dimension]
if { $dim == 1 } {
    if { $slc == 1 } {
        set interface [InterfaceFactory_create_well_radial_slc [char2string "${snapshot}"] $config $infodesk $gridfile]
    } else {
        set interface [InterfaceFactory_create_well_radial [char2string "${snapshot}"] $config $infodesk $gridfile]
    }
} elseif { $dim == 2 } {
    if { $slc == 1 } {
        set interface [InterfaceFactory_create_wire_slc [char2string "${snapshot}"] $config $infodesk $gridfile]
    } else {
        set interface [InterfaceFactory_create_wire [char2string "${snapshot}"] $config $infodesk $gridfile]
    }
} else {
    if { $slc == 1 } {
        set interface [InterfaceFactory_create_dot_slc [char2string "${snapshot}"] $config $infodesk $gridfile]
    } else {
        set interface [InterfaceFactory_create_dot [char2string "${snapshot}"] $config $infodesk $gridfile]
    }
}

###########################################
# read potential if it exists             #
###########################################
if { [file exists "${snapshot}_potential.dat"] } {
    set potential [read_vector_from_file "${snapshot}_potential.dat"]
    $interface set_potential $potential
}

###########################################
# set interface properties                #
###########################################
$interface set_minimal_edges [$snap_interface_data get "minimal_cb_edge"] [$snap_interface_data get "minimal_vb_edge"]
if { [$snap_interface_data get "band_barriers_set"] == 1 } {
    $interface set_maximal_edges [$snap_interface_data get "maximal_cb_edge"] [$snap_interface_data get "maximal_vb_edge"]
}

# well ... not yet finished
#set cb_edge [get_empty_vector_double]
#set vb_edge [get_empty_vector_double]
#$interface calculate_bandedges $cb_edge $vb_edge

$interface calculate

##########################################
# do slc                                 #
##########################################
if { $slc == 1 } {
    $interface calculate_luminescence [$slc_data get "cb_fermi"] [$slc_data get "vb_fermi"] [$slc_data get "temperature"]
    set omin [expr { [$interface get_omega_min] * $constants_hbar}]
    set omax [expr { [$interface get_omega_max] * $constants_hbar}]
    write_lumi_vector_to_file $omin $omax 3 [$interface get_absorption] "${snapshot}_absorption.dat"
    write_lumi_vector_to_file $omin $omax 3 [$interface get_spont_emission] "${snapshot}_spont_emission.dat"
}

##########################################
# evaluate result                        #
##########################################
set cb_valid [$interface get_number_of_valid_cb_subbands]
set vb_valid [$interface get_number_of_valid_vb_subbands]
set cb_bound [$interface get_number_of_bound_cb_subbands]
set vb_bound [$interface get_number_of_bound_vb_subbands]

puts [format "cb valid: %d, cb bound: %d" $cb_valid $cb_bound]
puts [format "vb valid: %d, vb bound: %d" $vb_valid $vb_bound]

$interface write_bandstructure "${snapshot}_bandstructure.dat"

set num_vertices [$geometry get_num_nodes]
write_lumi_vector_to_file 0.0 $num_vertices 1 [$interface get_cb_base_probability 0] "${snapshot}_cb_probability_0.dat"
write_lumi_vector_to_file 0.0 $num_vertices 1 [$interface get_vb_base_probability 0] "${snapshot}_vb_probability_0.dat"







