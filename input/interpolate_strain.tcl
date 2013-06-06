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
#  interpolate_strain.tcl                                     #
#                                                             #
#  interpolate and rotate data between elements               #
#                                                             #
#  (c) ratko veprek, iis/ethz, july. 2008                     #
#      veprek@iis.ee.ethz.ch                                  #
#                                                             #
#  purpose:                                                   #
#  - interpolate element data from a higher dimensional grid  #
#    onto lower dimensional grid using element mappings       #
#    (a list of highm elem number -> low dim elem number      #
#  - used by aqua to interpolate strain from higher dim to    #
#    lower structures
#                                                             #
###############################################################

proc show_options { msg } {
	puts ""
	puts "-> $msg"
	puts ""
	puts "usage: interpolate.tcl <element map (ascii file)> <element data file from> <element data file to>"
	puts ""
	puts ""
}

proc read_element_map { mylist myfile elem_num } {
	upvar $mylist target
	upvar $elem_num num
	set num 0
	set f [open $myfile "r"]
	set expected_element 0
	while { [gets $f line] >= 0} {
		set it [split $line]
		if { [llength $it] != 2 } {
			puts "-> oopss ... could not split $line into a list of two numbers ..."
			quit
		}
		if {$expected_element != [lindex $it 0]} {
	     	puts "-> wooops, expected element $expected_element on line $num, but i got the string $line"
			quit
		}
		set target([lindex $it 0]) [lindex $it 1]
		incr num
		incr expected_element
	}
	close $f
}

if { $::argc != 3 } {
	show_options "not enough arguments supplied!"
	quit
} else {

	set mapfile    [lindex $::argv 0]
	set sourcefile [lindex $::argv 1]
	set targetfile [lindex $::argv 2]
	set parser     [InputParser]
    set scalar_data 0

	#############################################
	# initialize rotation matrix to rotate      #
	# strain onto kp grid!                      #
	# default: identity matrix                  #
	#
	# rotation definition:
	# the strain is given in the current coord  #
	# system x. R maps x to x': x' = Rx         #
	# x' will be the new coordinate system of   #
	# lower dimension. so you map the strain    #
	# there on (by e' = R e R^t)                #
	#############################################
	set rotation    [RotationMatrix -args 3 3]
	set no_rotation 0
	if { $no_rotation == 1 } {
		$rotation set 0 0 1
		$rotation set 1 1 1
		$rotation set 2 2 1
		puts "-> DONT FORGET TO SET THE ROTATION APPROPRIATELY"
	} else {
		puts "-> applying custom rotation"
		$rotation set 0 2 1
		$rotation set 1 1 1
		$rotation set 2 0 1
		$rotation print
	}

	#############################################
	# check if files exist                      #
	#############################################
	if { ![file exists $mapfile] } {
		puts "-> file $mapfile does not exist"
		quit
	}
	#############################################
	# read element map                          #
	#############################################
	read_element_map elem_map $mapfile num_elements

	#############################################
	# load element data                         #
	#############################################
    if { $scalar_data == 1 } {
        set sourcedata [StdElementDatadouble]
        set targetdata [StdElementDatadouble]
    } else {
        set sourcedata [StrainField]
        set targetdata [StrainField -args $num_elements]
    }
	$parser read_binary $sourcedata $sourcefile

	#############################################
	# map data                                  #
	#############################################
    $targetdata set_length $num_elements [$sourcedata get_num_data_per_element]
	for {set ii 0} {$ii < $num_elements} {incr ii} {
		if { $elem_map($ii) != -1 } {
			for {set jj 0} {$jj < [$sourcedata get_num_data_per_element]} {incr jj} {
				$targetdata set_element_value $ii $jj [eval_dref [$sourcedata get_element_value $elem_map($ii) $jj]]
			}
		}
	}

	#############################################
	# rotate strain                             #
	#############################################
    if { $scalar_data == 0 } {
    	$targetdata rotate $rotation
    }

	#############################################
	# store mapped data                         #
	#############################################
	$parser write_binary $targetdata $targetfile
	puts "-> stored interpolated data to $targetfile"
}







