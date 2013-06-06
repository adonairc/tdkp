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



if { $dimension == 3 } {
    logger_output "TCLclc: clc does not work for quantum dots!" $LOG_ERROR
    exit
}
###################################
# build clc ingredients container #
###################################
if { [info exists cb_bands_disp] } {
    if { [info exists vb_bands_disp] } {
        set matrix_elements_disp [$matrix_elements get_disp_matrix_elements [$cb_bands_disp get_domain]]
        set clc_ingredients [CLCIngredients1DR -args $cb_bands_disp $vb_bands_disp $matrix_elements_disp $coulomb_data $slc_homogenous_broadening $slc_refractive_index $clc_static_dielectric_constant $quantized_volume]
    } else {
        set clc_ingredients [CLCIngredients1DR -args $cb_bands_disp $vb_bands $matrix_elements $coulomb_data $slc_homogenous_broadening $slc_refractive_index $clc_static_dielectric_constant $quantized_volume]
    }
} else {
    set clc_ingredients [CLCIngredients1DR -args $cb_bands $vb_bands $matrix_elements $coulomb_data $slc_homogenous_broadening $slc_refractive_index $clc_static_dielectric_constant $quantized_volume]
}

###################################
# build solver                    #
###################################
set clc_shf [CLCSHF1DR -args $clc_ingredients]
$clc_shf set_interpolation_density $slc_delta_k

###################################
# prepare omega determination     #
###################################
if { $slc_omega_range_mode == 1} {
    $clc_shf set_omega_range $slc_omega_min $slc_omega_max $slc_omega_num
} elseif { $slc_omega_range_mode == 2 } {
    set cb0 [complex2real [$cb_bands get_energy 0 0]]
    set vb0 [complex2real [$vb_bands get_energy 0 0]]
    set my_omega_min [expr {($cb0 - $vb0 - 0.15) / $constants_hbar}]
} else {
    set cbii 0
    set vbii 0
    ##############################################
    # find tightest bound bands                  #
    ##############################################
    for {set ii 1} {$ii < [$cb_bands get_number_of_bands]} {incr ii} {
        set e1 [complex2real [$cb_bands get_energy 0 $ii]]
        set e2 [complex2real [$cb_bands get_energy 0 $cbii]]
        if { $e1 < $e2 } {
            set cbii $ii
        }
    }
    for {set ii 1} {$ii < [$vb_bands get_number_of_bands]} {incr ii} {
        set e1 [complex2real [$vb_bands get_energy 0 $ii]]
        set e2 [complex2real [$vb_bands get_energy 0 $vbii]]
        if { $e1 > $e2 } {
            set vbii $ii
        }
    }
    set cb0 [complex2real [$cb_bands get_energy 0 $cbii]]
    set vb0 [complex2real [$vb_bands get_energy 0 $vbii]]
    set my_omega_min [expr {($cb0 - $vb0 - 0.15) / $constants_hbar}]
    set cb0 [complex2real [$cb_bands get_energy [expr {[$cb_bands get_number_of_k_values] - 1}] $cbii]]
    set vb0 [complex2real [$vb_bands get_energy [expr {[$vb_bands get_number_of_k_values] - 1}] $vbii]]
    set my_omega_max [expr {($cb0 - $vb0) / $constants_hbar}]
    $clc_shf set_omega_range $my_omega_min $my_omega_max $slc_omega_num
}

####################################
# solve for available fermi levels #
####################################
for {set ii 0} {$ii < [llength $fermi_levels]} {incr ii} {
    set f [lindex $fermi_levels $ii]
    set d [lindex $densities $ii]
    if { $slc_omega_range_mode == 2 } {
        set cf [lindex $f 0]
        set vf [lindex $f 1]
        $clc_shf set_omega_range $my_omega_min [expr {($cf - $vf + 0.7) / $constants_hbar}] $slc_omega_num
    }
    set res [$clc_shf calculate [lindex $f 0] [lindex $f 1] $temperature]
	plot_to_file "optical" $res [format "$result_directory/$plot_clc_optics_file" [lindex $d 0] [lindex $d 1]] "real"    
    $res -delete
}

