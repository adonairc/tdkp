
set nitrides_use_vurgaftman_2007 1

set wire_quantized_direction_x  $ey
set wire_quantized_direction_y  $ez
set wire_transverse_direction_z $ex

# wurzite intrinsic strain axes
# wurzite kp hamiltonians are [a,a,c] (so x axis is a per default)
# but for strain calculation that may be permuted
# values are also used in effective mass calculation for the hydrostatic
# strain offset
set wurzite_strain_axis_1  2
set wurzite_strain_axis_2  0
set wurzite_strain_axis_3  1

