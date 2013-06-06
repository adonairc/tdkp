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
# package scl.tcl                                       #
#                                                       #
# calculate simple luminescence                         #
#                                                       #
#########################################################

# check quantized volume
if { $quantized_volume == 0.0 } {
    logger_emit "TCLslc: quantized volume is ZERO! means: please write the names of the quantized regions into the file." $LOG_ERROR
    quit
}
# create slc object
if { [info exists vb_bands_disp] } {
    set matrix_elements_disp [$matrix_elements get_disp_matrix_elements [$cb_bands_disp get_domain]]
    set slc [SLC -args $cb_bands_disp $vb_bands_disp $matrix_elements_disp $slc_homogenous_broadening $slc_refractive_index $quantized_volume]
} elseif { [info exists cb_bands_disp] } {
    set slc [SLC -args $cb_bands_disp $vb_bands $matrix_elements $slc_homogenous_broadening $slc_refractive_index $quantized_volume]
} else {
    set slc [SLC -args $cb_bands $vb_bands $matrix_elements $slc_homogenous_broadening $slc_refractive_index $quantized_volume]
}
$slc set_interpolation_density $slc_delta_k
$slc set_omega_num $slc_omega_num
if { $slc_omega_range_mode == 1} {
    $slc set_omega_range $slc_omega_min $slc_omega_max $slc_omega_num
} elseif { $slc_omega_range_mode == 2 } {
    set cb0 [complex2real [$cb_bands get_energy 0 0]]
    set vb0 [complex2real [$vb_bands get_energy 0 0]]
    set my_omega_min [expr {($cb0 - $vb0 - 0.15) / $constants_hbar}]
}
# evaluate luminescence
for {set ii 0} {$ii < [llength $fermi_levels]} {incr ii} {
    set f [lindex $fermi_levels $ii]
    set d [lindex $densities $ii]
    if { $slc_omega_range_mode == 2 } {
        set cf [lindex $f 0]
        set vf [lindex $f 1]
        $slc set_omega_range $my_omega_min [expr {($cf - $vf + 0.7) / $constants_hbar}] $slc_omega_num
    }
    $slc calculate [lindex $f 0] [lindex $f 1] $temperature
	plot_to_file "optical" $slc [format "$result_directory/$plot_slc_optics_file" [lindex $d 0] [lindex $d 1]] "real"    

    #########################################################
    # evaluate B coefficent if requested (9.7.19) in chuang #
    #########################################################
    if { $slc_evaluate_B == 1 } {
        set Bnp 0.0
        set omega_d [expr {([$slc get_omega_max] - [$slc get_omega_min]) / ([$slc get_omega_num] - 1.0)}]
        for {set jj 1} {$jj < [$slc get_omega_num]} {incr jj} {
            set sp_0 [$slc get_spont_emission_average [expr {$jj - 1}]]
            set sp_1 [$slc get_spont_emission_average $jj]
            set Bnp [expr {($sp_0 + $sp_1) / 2.0 * $omega_d * $constants_hbar + $Bnp}]
        }
        set n [expr {[lindex $d 0] }]
        set p [expr {[lindex $d 1] }]
        set B_coeff [expr {$Bnp * $slc_B_scaling/ ($n * $p)}]
        logger_emit [format "TCLslc: spontaneous emission B coefficent = %g (cm^3/s)" $B_coeff]
    }
    #########################################################
    # emit some statistics                                  #
    #########################################################
    for {set jj 0} {$jj < 7} {incr jj} {
        set idxs($jj) 0
        set pval($jj) 0
    }
    for {set ww 1} {$ww < [$slc get_omega_num]} {incr ww} {
        for {set jj 0} {$jj < 7} {incr jj} {
            if {$jj < 3} {
                set valmax [$slc get_absorption $jj $idxs($jj)]
                set val    [$slc get_absorption $jj $ww]
            } elseif { $jj < 6 } {
                set valmax [$slc get_spont_emission [expr {$jj - 3}] $idxs($jj)]
                set val    [$slc get_spont_emission [expr {$jj - 3}] $ww]
            } else {
                set valmax [$slc get_spont_emission_average $idxs($jj)]
                set val    [$slc get_spont_emission_average $ww]
            }
            if {($jj >= 3 && $val > $valmax) || ($jj < 3 && $val < $valmax)} {
                set idxs($jj) $ww
                set valmax $val
            }
            set pval($jj) $valmax
        }
    }
    set omega_d [expr {([$slc get_omega_max] - [$slc get_omega_min]) / ([$slc get_omega_num] - 1.0)}]
    set omega_min [$slc get_omega_min]
    for {set jj 0} {$jj < 7} {incr jj} {
        set omx($jj) [expr {($omega_min + $omega_d * $idxs($jj)) * $constants_hbar}]
    }
    logger_emit [format "TCLslc: peak values for gain:\n  px = %g at %g eV\n  py = %g at %g eV\n  pz = %g at %g eV" $pval(0) $omx(0) $pval(1) $omx(1) $pval(2) $omx(2)]
    logger_emit [format "TCLslc: peak values for spont emission:\n  px = %g at %g eV\n  py = %g at %g eV\n  pz = %g at %g eV\n  avg = %g at %g eV" $pval(3) $omx(3) $pval(4) $omx(4) $pval(5) $omx(5) $pval(6) $omx(6)]
}

