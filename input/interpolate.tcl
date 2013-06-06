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
#  interpolate.tcl                                            #
#                                                             #
#  interpolate data from one grid to another                  #
#                                                             #
#  (c) ratko veprek, iis/ethz, september. 2008                #
#      veprek@iis.ee.ethz.ch                                  #
#                                                             #
#  purpose:                                                   #
#  - as said: interpolate data. currently used for extracting #
#    3D strain data for lower dimensional aqua slices         #
#                                                             #
###############################################################


if { $argc != 5 } {
    puts " usage: interpolate.tcl <source grid> <source data> <target_grid> <target data> <x:y:z (translation vector)>"
    exit
}

set logmsg [get_logger_object]
$logmsg set_level $LOG_INFO_DEVEL2


set fsourcegrid [lindex $::argv 0]
set fsourcedata [lindex $::argv 1]
set ftargetgrid [lindex $::argv 2]
set ftargetdata [lindex $::argv 3]
set translation [lindex $::argv 4]

set parser [InputParser]

##########################################
# tell me if we are interpolating strain #
##########################################
set its_strain 1

#######################################
# build rotation matrix               #
#######################################
set rmatrix [RotationMatrix -args 3 3]
$rmatrix set 2 0 1.0
$rmatrix set 1 1 1.0
$rmatrix set 0 2 1.0

#######################################
# build translation                   #
#######################################
set vl [split $translation ':']
if { [llength $vl] != 3 } {
    puts "-> can't parse translation vector $translation"
    exit
} else {
    set x [lindex $vl 0]
    set y [lindex $vl 1]
    set z [lindex $vl 2]
    set translation [Vector3D -args $x $y $z]
}

#######################################
# read grids                          #
#######################################
set sourcegrid [$parser read_geometry [DFISEGridReader -args $fsourcegrid] 2]
set targetgrid [$parser read_geometry [DFISEGridReader -args $ftargetgrid]]

#######################################
# custom scaling for grids            #
#######################################
set scale_source 1000.0
if { $scale_source != 1.0 } {
    $sourcegrid rescale_node_coordinates $scale_source
}

#######################################
# read data                           #
#######################################
if { $its_strain == 1 } {
    set sourcedata [StrainField -args $fsourcedata]
    set targetdata [StrainField -args $targetgrid]
} else {
    set sourcedata [StdNodeDatadouble]
    set targetdata [StdNodeDatadouble]
    $parser read_binary $sourcedata $fsourcedata
}

#######################################
# build interpolator                  #
#######################################
set ipol [GridInterpolator -args $sourcegrid $targetgrid $rmatrix $translation]
if { $its_strain == 1 } {
    $ipol interpolate $sourcedata $targetdata
} else {
    $ipol interpolate_double $sourcedata $targetdata
}

#######################################
# rotate strain (if it is strain)     #
#######################################
if { $its_strain == 1 } {

    set transpose [$rmatrix get_transpose]
    $targetdata rotate $transpose
    #######################################
    # write data                          #
    #######################################
    $targetdata write_binary $ftargetdata
    $parser write_to_dfise_file $targetgrid $targetdata "elem_strain.dat"
    $parser write_ascii_double $targetgrid $targetdata "elem_strain_ascii.dat"
} else {
    $parser write_binary $targetdata $ftargetdata
    $parser write_to_dfise_file $targetgrid $targetdata "elem_strain.dat"
    $parser write_ascii_double $targetgrid $targetdata "elem_strain_ascii.dat"
}




