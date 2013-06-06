
# set c-axis to y
set wurzite_strain_axis_1  0
set wurzite_strain_axis_2  2
set wurzite_strain_axis_3  1

# set directions 
set wire_quantized_direction_x  $ex
set wire_quantized_direction_y  $ez
set wire_transverse_direction_z [Vector3D -args 0.0 -1.0 0.0]

# grid is in microns
set scale_coordinates 1000.0

# use recent parameter set
set nitrides_use_vurgaftman_2007 1

# y-axis is unstrained, structure may relax
set axis_x_strained             1
set axis_y_strained             0
set axis_z_strained             1
