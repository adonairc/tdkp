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
#  tow.tcl                                                    #
#                                                             #
#  theory of eWerything                                       #
#                                                             #
#  (c) ratko veprek, iis/ethz, jan. 2008                      #
#      veprek@iis.ee.ethz.ch                                  #
#                                                             #
#  purpose:                                                   #
#  - provide a simple user interface for all kinds of         #
#    calculations using tdkp                                  #
#  - offers options to calculate luminescence, gain, matrix   #
#    elements, fermi levels,                                  #
#  - various plotting of data                                 #
#  - loading and storing results                              #
#                                                             #
# default configuration and material properites database can  #
# be found in TDKPCONFPATH/../input                           #
#                                                             #
###############################################################

###############################
# check environment and files #
###############################
if { ![info exists env(TDKPCONFPATH)] } {
    if { [file exists "/usr/share/tdkp/conf/tow.tcl"] } {
        set env(TDKPCONFPATH)  "/usr/share/tdkp/conf"
    } else {
        puts "-> ** error **: please set TDKPCONFPATH"
        exit
    }
}
set confpath $env(TDKPCONFPATH)
if { ![file exists "$confpath/../input/tdkp_helpers.tcl"] } {
    puts "-> ** error **: can't find tdkp_helpers.tcl in $confpath/../input"
    puts "                check your path"
    exit
}
set sourcepath "$confpath/../input"
source "$confpath/../input/tdkp_helpers.tcl"


#######################
# show script options #
#######################
proc show_bands_tcl_options { err } {
    puts ""
    puts "-> input error: $err"
    puts ""
    puts "-> usage: tow.tcl <options>"
    puts ""
    puts "   theory of everything - a simple tcl script access to tdkp"
    puts ""
    puts "   simulation options:"
    puts "      -bands                        calculate bandstructure"
    puts "      -matrix_elements              calculate matrix elements"
    puts "      -slc          <N:P<,N:P>>     calculate luminescence for given densities (3d density / nm^3)"
    puts "                                    bands and matrix elements are calculated or read from result"
    puts "                                    directory if available"
    puts "      -slc          <filename>      calculate luminescence for the densities in the file"
    puts "      -slc_fermi    <filename>      calculate luminescence for the fermi levels in the file"
    puts "      -clc          <filename>      calculate luminescence including many body effects for the densities in the file"
    puts "      -clc_fermi    <filename>      calculate luminescence including many body effects for the fermi levels in the file"
    puts "      -fermi        <N:P<,N:P>>     calculate the fermi level for the given densities"
    puts "      -fermi        <filename>      calculate the fermi level for the densities in the file"
    puts "      -strain                       calculate (or read from cache if available) strains"
    puts "      -strain       <reference mat> calculate strain using the lattice constant of the reference material"
    puts "      -delta_strain                 calculate bandedge change due to strain field (requires -strain)"
    puts "      -coulomb      diag | all      calculate coulomb matrix elements (diag: n1 = n4, n2 = n3)"
    puts "      -coulomb      <combo<:combo>> calculate selected coulomb matrix elements (eg. combo = cb0-cb0-cb1-cb1)"
    puts "      -iso          cb1:cb2:vb1...  create iso-surface plots of the given states (interactively asking for iso-surf values)"
    puts "      -piezo                        calculate (or read from cache if available) piezo and spontaneous polarizations"
    puts ""
    puts "   input options:"
    puts "      -bulk         <material>      use if you like to calculate the bulk bandstructure"
    puts "      -grid         <gridfile>      the grid file"
    puts "      -deg          <polynomial-d>  degree of the finite elements"
    puts "      -bc           <tcl-file>      use user-defined boundary conditions (bc_strain, bc_poisson, bc_bands)"
    puts "      -model        <kp model>      effmass, kp4x4, kp6x6, kp8x8, kp6x6WZ, kp8x8WZ"
    puts "      -nk           <value gt 0>    number of k values in bandstructure>"
    puts "      -krange       <kmin:kmax>     k-range in 2pi / nm"
    puts "      -tdir         <x:y:z>         well/bulk transverse direction"
    puts "      -qdir         <x:y:z>         well quantized direction"
    puts "      -ncb                          num bands to solve for in conduction band (0 = don't calc)"
    puts "      -nvb                          num bands to solve for in valence band (0 = don't calc)"
    puts "      -energy_guess <float:float>   energy guess for k = 0 for cb:vb"
    puts "      -scale_coords <float>         scale node coordinates (change system size)"
    puts "      -temp         <temperature>   temperature in kelvin, default 300K "
    puts "      -tclscript    <tcl script>    source tcl script at the end of the calculations"
    puts "                                    for user defined post processing"
    puts ""
    puts "   output options:"
    puts "      -resdir       <dirname>       directory where to store the results (./ default)"
    puts "      -verbose                      be very verbose on what you are doing"
    puts "      -quiet                        be very quiet (only show warnings and errors)"
    puts ""
    puts "   cache options:"
    puts "      -load_bands                   load bands from binary file"
    puts "      -load_matrix_elements         load matrix elements from binary file"
    puts "      -load_coulomb                 load coulomb matrix elements from binary file"
    puts "      -nocache                      don't use previously calculated bands, matrix elements, etc."
    puts "                                    the policy is to load results if required and available but "
    puts "                                    no explicit command was given."
    puts ""
    puts "   program control options:"
    puts "      -conf         <conf file>"
    puts "      -log          <logfile>"
    puts "      -plot         <itm1<:itm2>    plot data (all, probabilities, bands, bandedges, dos,"
	puts "                                    transitions, fermi_levels, strain, hydrostatic_strain,"
    puts "                                    displacement, delta_bandedges, optical, coulomb, piezo)"
	puts "                                    you can define the output format by appending _med, _ascii"
	puts "                                    to the output option"
    puts "      -load_material                load material from files instead from tcl db"
    puts "      -dump_material                dump material into files placed in result_directory"
    puts "      -help                         print this text"
    puts ""
    puts "-> another piece of software brought to you by ratko veprek"
    puts "   (c) 2008, iis, eth zuerich"
    puts ""
    puts ""
    puts "-> input error: $err"
    puts ""
}

#####################################################
# set up existing plot informations                 #
#####################################################
# this list contains our 'known' plot options
set plot_valid_options [list \
	"strain" "hydrostatic_strain" "displacement" "delta_bandedges" "piezo" "potential" \
	"bandedges" "bands" "probabilities" "dos" "fermi_levels" "transitions" "matrix_elements" \
	"coulomb" "optical" ]

########################
# determine user input #
########################
puts "-> you started tow.tcl $argv"
if { $argc == 0 && ![info exists options] } {
    show_bands_tcl_options "TCLtow: not enough command line arguments!"
    exit
} else {
    # determine command line options and store it into array options
    determine_command_line_options options
}

#######################################
# create result directory if required #
#######################################
set result_directory "."
if { [info exists options(-resdir)] } {
    set result_directory $options(-resdir)
    if { ![file exists $result_directory] } {
        file mkdir $result_directory
    }
}

###########################################
# setup config classes and output options #
###########################################
set config [get_configuration_object];
set logmsg [get_logger_object];

if { [info exists options(-verbose)] } {
    $logmsg set_level $LOG_INFO_DEVEL2;
} elseif { [info exists options(-quiet)] } {
    $logmsg set_level $LOG_WARN;
} else {
    $logmsg set_level $LOG_INFO;
}

###########################################################
# load std config file and then override with user config #
###########################################################
source [format "%s/../input/default_config.tcl" $env(TDKPCONFPATH)]
if { [info exists options(-conf)] } {
    source $options(-conf)
}

##########################################
# standard log file                      #
##########################################
if { $enable_log_file == 1 || [info exists options(-log)] } {
    if { [info exists options(-log)] } {
        set logfile "$result_directory/$options(-log)"
    } else {
        set logfile "$result_directory/tdkp.log"
    }
    if { $logfile == "1" } {
        show_bands_tcl_options "no log file name defined"
        exit
    }
    puts "-> TCLtow: copying output to file $logfile"
    logger_output_to_file $logfile
    # write command line arguments to logfile (and thus reset it)
    set d [clock format [clock seconds] -format "%D %T"]
    logger_emit " ------------------------------------------------------\n $argv"
    logger_emit " started at $d"
    logger_emit " ------------------------------------------------------"
}

##################################################
# set up material database and i/o parser/writer #
##################################################
set material_database [MaterialDatabase];
set parser [InputParser];

###########################################
# overwrite config value on users request #
###########################################
if { [info exists options(-help)] } {
    show_bands_tcl_options "TCLtow: you asked for help ... here is your help"
    exit
}
if { [info exists options(-nk)] } {
    set num_k_values $options(-nk)
}
if { [info exists options(-krange)] } {
    regexp {([0-9\.e+-]+):([0-9\.e+-]+)} $options(-krange) matches kmin kmax
}
#################################
# set coulomb_q_max to 2 * kmax #
#################################
if { $coulomb_q_max == -1 } {
    set coulomb_q_max $kmax
}

if { [info exists options(-tdir)] } {
    regexp {([0-9\.e+-]+):([0-9\.e+-]+):([0-9\.e+-]+)} $options(-tdir) matches x y z
    set bulk_transverse_direction [Vector3D -args $x $y $z]
    set well_transverse_direction $bulk_transverse_direction
}
if { [info exists options(-qdir)] } {
    regexp {([0-9\.e+-]+):([0-9\.e+-]+):([0-9\.e+-]+)} $options(-qdir) matches x y z
    set well_quantized_direction [Vector3D -args $x $y $z]
}

if { [info exists options(-ncb)] } {
    set num_cb_sol $options(-ncb)
}
if { [info exists options(-nvb)] } {
    set num_vb_sol $options(-nvb)
}

if { [info exists options(-plot)] } {
	# new output handling
	process_command_line_plot_requests [split $options(-plot) ":"]
	set w "TCLplot: the following data will, if available, be plotted:"
	foreach req $plot_request {
		set w "${w} ${req}"
	}
	logger_emit $w $LOG_INFO_DEVEL2
}

if { [info exists options(-scale_coords)] } {
    set scale_coordinates $options(-scale_coords)
}
if { [info exists options(-load_material)] } {
    set material_load_from_file 1
}
if { [info exists options(-dump_material)] } {
    set material_file_dump 1
}
if { [info exists options(-temp)] } {
    set temperature $options(-temp)
}
if { [info exists options(-model)] } {
    if { $options(-model) != "kp8x8" } {
        set material_tcl_use_elliptic_params 0
        set material_tcl_use_optimal_Nm      0
    }
}

####################
# check user input #
####################
check_and_process_user_input options
# ncb or nvb must be greater 0
if { $num_cb_sol == 0 && $num_vb_sol == 0 } {
    show_bands_tcl_options "TCLtow: either ncb or nvb must at least be > 0 ..."
    exit
}
# need cb and vb for matrix elements and slc
if { ($num_cb_sol == 0 || $num_vb_sol == 0) } {
    if { [info exists options(-clc)] || [info exists options(-slc)] || [info exists options(-matrix_elements) ] } {
        show_bands_tcl_options "TCLtow: either ncb and nvb must be > 0 ..."
        exit
    }
}
# band structure loading requires bands to be loadable
if { [info exists options(-load_bands)] && !([file exists "$result_directory/$vb_bands_binary_file"] || [file exists "$result_directory/$cb_bands_binary_file"])} {
    show_bands_tcl_options "TCLtow: can not load bands because the binary files are missing ..."
    exit
}
# matrix element loading?
if { [info exists options(-load_matrix_elements)] && !([file exists "$result_directory/$binary_matrix_elements_file"])} {
    show_bands_tcl_options "TCLtow: can not load matrix elements because the binary file is missing ..."
    exit
}
# coulomb matrix elements loading?
if { [info exists options(-load_coulomb)] && !([file exists "$result_directory/$coulomb_binary_file"])} {
    show_bands_tcl_options "TCLtow: can not load coulomb matrix elements because the binary file is missing ..."
    exit
}
# strain loading?
if { [info exists options(-load_strain)] && !([file exists "$result_directory/$binary_strain_file"])} {
    show_bands_tcl_options "TCLtow: can not load strains because the binary file is missing ..."
    exit
}


#########################################################
# load geometry and check dimensionality of the problem #
#########################################################
source "$sourcepath/grid.tcl"


##############################
# get the strain field ready #
##############################
if { [info exists options(-strain)] || [info exists options(-load_strain)] } {
    source "$sourcepath/strain.tcl"
}

######################################
# calculate piezo charge / potential #
######################################
if { [info exists options(-piezo)] || [info exists options(-load_piezo)] } {
    source "$sourcepath/piezo.tcl"
}

##############################
# calculate bandedges        #
##############################
if { ([info exists options(-bands)] || [ plot_requested "bandedges" ] || [info exists options(-delta_strain)]) && $dimension > 0 } {
    source "$sourcepath/bandedges.tcl"
}

##############################
# source band calculation    #
##############################
source "$sourcepath/bands.tcl"


##########################################################################
# second section: tdkp post processing                                   #
##########################################################################


###############################
# calculate density of states #
###############################
if { $dimension < 3 && [plot_requested "dos"] } {
    source "$sourcepath/dos.tcl"
}


###################################################################################
# create momentum operator if required by matrix_elements or load_matrix_elements #
###################################################################################
if { [info exists options(-matrix_elements)] || [info exists options(-load_matrix_elements)] } {
    source "$sourcepath/matrix_elements.tcl"
}

############################################
# calculate coulomb matrix elements        #
############################################
if { [info exists options(-coulomb)] || [info exists options(-load_coulomb)] } {
    source "$sourcepath/coulomb.tcl"
}

##########################
# calculate fermi levels #
##########################
if { [info exists options(-slc_fermi)] || [info exists options(-fermi)] } {
    source "$sourcepath/fermi.tcl"
}

##########################
# calculate luminescence #
##########################
if { [info exists options(-slc)] } {
    source "$sourcepath/slc.tcl"
}
if { [info exists options(-clc)] } {
    source "$sourcepath/clc.tcl"
}

########################################
# create isosurface plots if requested #
########################################
if { [info exists options(-iso)] } {
    create_isosurface_plots $options(-iso) $geometry $cb_bands $vb_bands
}

#############################################
# source user defined tcl script on request #
#############################################
if { [info exists options(-tclscript)] } {
    if { $options(-tclscript) != "1" } {
        source $options(-tclscript)
    }
}

########################################
# quit if noquit not requested         #
########################################
if { ![info exists options(-noquit)] } {
    quit
}




