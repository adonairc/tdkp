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
#  default_config.tcl                                         #
#                                                             #
#  (c) ratko veprek, iis/ethz, jan. 2008                      #
#      veprek@iis.ee.ethz.ch                                  #
#                                                             #
#  purpose:                                                   #
#  - default configuration for tdkp tcl scripts               #
###############################################################

######################
# simulation options #
######################
# number of k values in bandstructure
set num_k_values 16
set kmin         0.0e0
set kmax         [expr {0.25 * (2.0 * $constants_pi / 0.5)}];

# number of solutions we initially solve for in kp calculation
set num_cb_sol   12
set num_vb_sol   16

# standard transverse direction of a bulk bandstructure
set bulk_transverse_direction [Vector3D -args 0.0 0.0 1.0];

# standard transverse direction of a well bandstructure
set well_transverse_direction [Vector3D -args 0.0 1.0 1.0];
set well_quantized_direction  [Vector3D -args 1.0 0.0 0.0];

# standard directions for wire bandstructures
set wire_quantized_direction_x  $ex
set wire_quantized_direction_y  $ey
set wire_transverse_direction_z $ez

# wurzite intrinsic strain axes
# wurzite kp hamiltonians are [a,a,c] (so x axis is a per default)
# but for strain calculation that may be permuted
# values are also used in effective mass calculation for the hydrostatic
# strain offset
set wurzite_strain_axis_1  0
set wurzite_strain_axis_2  1
set wurzite_strain_axis_3  2

# standard direction for the dot
set dot_direction_x             $ex
set dot_direction_y             $ey
set dot_direction_z             $ez

# strained axes
# when calculating strains, we assume that the following access are
# intrinsically strained
set axis_x_strained             1
set axis_y_strained             1
set axis_z_strained             1

# scale system coordinates (1 == don't scale anything)
set scale_coordinates           1

# refractive index
set slc_refractive_index        3.65
# static dielectric constant
set clc_static_dielectric_constant 10.3
# broadening in meV => meV / hbar / 1000 =
set slc_homogenous_broadening   9e12

# optical matrix parameter for effective mass matrix elements (in [eV])
set effmass_optical_matrix_param 25.5

# mode 0: default (omega range is taken from bandstructure)
# mode 1: using user values
# mode 2: setting omega range according to bandstrcture and fermi level
set slc_omega_range_mode        0
set slc_omega_min               [expr { 0.4 / $constants_hbar } ]
set slc_omega_max               [expr { 2.0 / $constants_hbar } ]
set slc_omega_num               800
# interpolation density in k-space (if absorption / luminescence oscillate, use a lower value)
set slc_delta_k                 0.0025
# evaluate B coefficent
set slc_evaluate_B              1
# Bnp -> B = Bnp / np and Bnp is 1/(s nm^3) and B (nm^3/s)
set slc_B_scaling               1.0e-21

# default temperature in K
set temperature                 300.0

# coulomb matrix element evaluation params
# don't start at 0
set coulomb_q_min               0.0001
# warning, coulomb_q_max MUST BE 2*kmax for clc calculations.
# -1 means: please set coulomb_q_max to 2.0 * kmax
set coulomb_q_max               -1
set coulomb_q_num               26
set coulomb_progression         3.0

# band degeneracy threshold consider bands that start at same energy to be degenerate
set band_degeneracy_threshold   1.0e-5

# decross spin-degenerate bands (means: we resort and
# check that the same subband with same spin idx is
# orthognal to the one with the oder spin idx
set decross_bands 1

# use user defined settings for estimate of barrier
# bandedges (for determination of bound / unbound states)
# this is required e.g. if a quantized region touches the
# system boundary. any region touching the system boundary
# is treated as barrier material. in that case, you may define
# the barrier edges by yourself.
set use_user_defined_unbound_edges 0
set user_defined_unbound_cb_edge   3.0
set user_defined_unbound_vb_edge  -1.0

# disable regions with material "Gas" in band structure calculations
set disable_gas_for_bandstructure_calc 1

# set to true if metis ND should be used to reorder vertices
# when doing the bandstructure calculation
set graph_reordering_using_metis 0

##################################
# quantized volumes              #
##################################
# densities etc. are usually given in 3D densities.
# to calculate the fermi level, I therefore need the
# quantized volume, which can be calculated automatically.
# in some cases however, it is desirable to control this
# value by hand. so if quantized_volume_pref > 0, the given
# value is used
set quantized_volume_predef -1.0

###################
# pml options     #
###################
# if you wan't to use pml, simply set enable_pml 1 and
# name some of your regions PML_<id>_mx (for pml at left boundary, _px for
# pml at right boundary)
set enable_pml 1
set pml_print_bands_orderered_by_tau 1
set pml_alpha 1.0
set pml_beta  1.4
set pml_power 2
set pml_threshold_stable_states_max_imag_value 0.001
set pml_resort_according_to_lifetime 1

###################
# conf options    #
###################
set enable_log_file 1

# read quantized regions from file (if file exists)
# else we use the standard algorithm
set read_quantized_regions_from_file "quantized_regions.dat"

# load material from files instead from tcl database
set material_load_from_file     0

# write material properties to screen
set material_tcl_dump           1

# write material properties to files (only in conjuction with tcl material database)
set material_file_dump          0

# use optimal Nm for kp 8x8
set material_tcl_use_optimal_Nm 1

# use adjusted values for optmat and N for kp 8x8 (is disabled if model is not kp 8x8)
set material_tcl_use_elliptic_params 1
set material_tcl_elliptic_threshold  0.01

# nitrides: use vurgaftman 2007 parameteres
set nitrides_use_vurgaftman_2007 0


############################################
# output options NEW                       # 
############################################

# primary plot format
set plot_default_format       "med"
set plot_file_formats [list "binary" "sebise" "ascii" "med"]

# -----------------------------------------------
# default plot requests, required for result caching 
# within tdkp
# -----------------------------------------------
set plot_request [list "bandedges_binary" "optical_ascii" ]

# ----------------------------------------------
# strains
# ----------------------------------------------
set binary_strain_file          	"strain.bin"
set plot_strain_file            	"strain"
set plot_displacement_file      	"displacement"
set plot_hydrostatic_strain_file	"hydrostatic_strain"

# ----------------------------------------------
# bandedges
# ----------------------------------------------
# plot bandedges
set plot_bandedges_file       	"bandedges"
# the plot file is nodal
set plot_delta_bandedges_file   "bandedges_delta_nodal"
# binary file is stored element wise!
set binary_delta_bandedges_file "bandedges_delta.bin"


# -----------------------------------------------
# bands
# -----------------------------------------------
# keep only valid
set binary_keep_only_valid      1
# keep only bound AND valid states
set binary_keep_only_bound      1
# store binary cache file (i.e. all the bs is stored in here for later postprocessing)
set binary_write_bands          1
set vb_bands_binary_file        "vb_bands.bin"
set cb_bands_binary_file        "cb_bands.bin"
# default: plot bandstructure to ascii file
lappend plot_request "bands_ascii"
# filenames
set cb_bands_file         "cb_disp"
set vb_bands_file         "vb_disp"

# ----------------------------------------------
# probabilities
# ----------------------------------------------
# for which k points (index of it!) should we plot the values?
# possible options are:
#   integer n: plot probabilities given k-index point
#   string of two integers m:n gives probabilities of all k points inbetween
#   string all: plots all
set plot_probabilities_k_points 0
set cb_probability_file    "cb_probability_b%02d_k%s"
set vb_probability_file    "vb_probability_b%02d_k%s"


# ----------------------------------------------
# optics output file
# ----------------------------------------------
set plot_slc_optics_file  		"optics_n%g_p%g"
set plot_clc_optics_file 		"clc_optics_n%g_p%g"

# ----------------------------------------------
# coulomb output
# ----------------------------------------------
set plot_coulomb_file           "coulomb"
set coulomb_binary_file			"coulomb.bin"	

# ----------------------------------------------
# dos output file
# ----------------------------------------------
set plot_cb_dos_file           	"cb_dos"
set plot_vb_dos_file           	"vb_dos"

# ----------------------------------------------
# fermi levels (pure ascii)
# ----------------------------------------------
set plot_cb_fermi_file			"cb_fermi"
set plot_vb_fermi_file			"vb_fermi"

# ----------------------------------------------
# matrix elements and transitions
# ----------------------------------------------
# calculate matrix elements only for bound bands
set matrix_elements_bound_only   1
# matrix elements binary value caching file
set binary_write_matrix_elements 1
set binary_matrix_elements_file	 "matrix_elements.bin"
# plot files
set plot_matrix_elements_file    "matrix_elements"
set plot_transitions_file			"transitions.dat"
# transition strength at kmin (matrix cb <-> vb) for degenerate transitions
set plot_transition_strength_file 	"transition_matrix.dat"

# ----------------------------------------------
# polarization output
# ----------------------------------------------
set binary_pol_surf_charge_file  "pol_surface_charges.bin"
set binary_pol_vol_charge_file   "pol_volume_charges.bin"
set plot_all_charge_file         "polarization_charges"

# poisson output
set plot_potential_file  		"potential"
set binary_element_potential	"piezo_potential_element.bin"


