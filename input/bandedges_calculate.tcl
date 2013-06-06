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

#####################################################
# package bandegdes_calculate.tcl                    #
#                                                    #
# core commands calculating the bandedges, called    #
# by bandedges.tcl with appropriately prepared       #
# objects                                            #
######################################################

set bandedges [MergedElementDatadouble]
if { $cb_effmass == 1 && $vb_effmass == 1} {
    $em_problem set_solution_type $holes
    $em_problem prepare
    set vb_bandedges [$em_problem get_bandedges]
    $bandedges add_data $vb_bandedges "vb"
    $em_problem set_solution_type $electrons
    $em_problem prepare
    set cb_bandedges [$em_problem get_bandedges]
    $bandedges add_data $cb_bandedges "cb"
    set vb_idx 0
    set cb_idx 1
} else {
    $vb_problem prepare
    set vb_bandedges [$vb_problem get_bandedges]
    #########################################
    # vb kp / cb effmass                    #
    #########################################
    if { $cb_effmass == 1 } {
        $bandedges add_data $vb_bandedges "vb"
        $em_problem set_solution_type $electrons
        $em_problem prepare
        set cb_bandedges [$em_problem get_bandedges]
        $bandedges add_data $cb_bandedges "cb"
        set vb_idx [expr {[$vb_bandedges get_num_data_per_element] - 1}]
        set cb_idx [expr {$vb_idx + 1}]
    } else {
    #########################################
    # kp 8x8                                #
    #########################################
        $bandedges -delete
        set bandedges $vb_bandedges
        set cb_idx 6
        set vb_idx 5
    }
}