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
# package fermi.tcl                                     #
#                                                       #
# calculate fermi levels                                #
#                                                       #
#########################################################

# check quantized volume
if { $quantized_volume == 0.0 } {
    logger_emit "TCLfermi: quantized volume is ZERO! means: please write the names of the quantized regions into the file." $LOG_ERROR
    quit
}

############################################
# read fermi levels from file if requested #
############################################
if { [info exists options(-slc_fermi)] || [info exists options(-clc_fermi)] } {

    if { [info exists cb_bands_disp] } {
        set cb_fermi [Fermi -args $cb_bands_disp $electrons $quantized_volume]
    } else {
        set cb_fermi [Fermi -args $cb_bands $electrons $quantized_volume]
    }
    if { [info exists vb_bands_disp] } {
        set vb_fermi [Fermi -args $vb_bands_disp $holes     $quantized_volume]
    } else {
        set vb_fermi [Fermi -args $vb_bands $holes     $quantized_volume]
    }
    set f [open $options(-slc_fermi) "r"]
    while { [gets $f line] >= 0 } {
        set line [string trim $line " "]
        if { $line != "" } {
            if { ![regexp {([0-9eE+-\.]+)\s+([0-9eE+-\.]+)} $line match cf hf] } {
                logger_emit "TCLfermi: can not parse line in fermi level file: $line" $LOG_ERROR
                exit
            }
            lappend fermi_levels [list $cf $hf $temperature]
            set n [$cb_fermi calculate_density $cf $temperature]
            set p [$vb_fermi calculate_density $hf $temperature]
            lappend densities [list $n $p]
        }
    }
    close $f


} elseif { [info exists options(-fermi)] } {

    #########################################
    # read densities from file if requested #
    #########################################
    if { [info exists densities] } {
        unset densities
    }
    if { [file exists $options(-fermi)] } {
        logger_emit "TCLfermi: reading densities from file $options(-fermi)"
        set f [open $options(-fermi) "r"]

        while { [gets $f line] >= 0 } {
            set line [string trim $line " "]

            if { $line != "" } {
                if { ![regexp {([0-9eE+-\.]+)\s+([0-9eE+-\.]+)} $line match n p] } {
                    logger_emit "TCLfermi: can not parse line in density file: $line" $LOG_ERROR
                    exit
                }
                lappend densities [list $n $p]
            }

        }
        close $f
    } else {
        foreach f [split $options(-fermi) ","] {
            lappend densities [split $f ":"]
        }
    }
    if { [info exists cb_bands_disp] } {
        set cb_fermi [Fermi -args $cb_bands_disp $electrons $quantized_volume]
    } else {
        set cb_fermi [Fermi -args $cb_bands $electrons $quantized_volume]
    }
    if { [info exists vb_bands_disp] } {
        set vb_fermi [Fermi -args $vb_bands_disp $holes     $quantized_volume]
    } else {
        set vb_fermi [Fermi -args $vb_bands $holes     $quantized_volume]
    }

    if { [info exists fermi_levels] } {
        unset fermi_levels
    }

    ######################################################
    # loop over all densities and calculate fermi levels #
    ######################################################
    foreach dset $densities {
        if { [llength $dset] != 2 } {
            logger_emit "TCLfermi: can not parse density: $dset" $LOG_ERROR
            exit
        }
        set n   [lindex $dset 0]
        set p   [lindex $dset 1]

        set nEf [$cb_fermi calculate_fermi_level $n $temperature]
        set pEf [$vb_fermi calculate_fermi_level $p $temperature]

        lappend fermi_levels [list $nEf $pEf]

    }

}

##################################
# plot fermi levels if requested #
##################################
if { [plot_requested "fermi_levels"] && [info exists fermi_levels] && [info exists densities] } {

    logger_emit "TCLfermi: writing fermi levels to file"

    set cb_f [open "$result_directory/${plot_cb_fermi_file}.dat" "w"]
    set vb_f [open "$result_directory/${plot_vb_fermi_file}.dat" "w"]

    puts $cb_f "# ndens  cb_fermi"
    puts $vb_f "# pdens  vb_fermi"

    for {set ii 0} {$ii < [llength $fermi_levels]} {incr ii} {
        set dens [lindex $densities    $ii]
        set ferm [lindex $fermi_levels $ii]
        puts $cb_f [format "%s %s " [lindex $dens 0] [lindex $ferm 0] ]
        puts $vb_f [format "%s %s " [lindex $dens 1] [lindex $ferm 1] ]
    }
    close $cb_f
    close $vb_f

}
