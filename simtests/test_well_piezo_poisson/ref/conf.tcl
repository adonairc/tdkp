
set nitrides_use_vurgaftman_2007 1

# standard transverse direction of a well bandstructure
set well_transverse_direction [Vector3D -args 1.0 1.0 0.0];
set well_quantized_direction  [Vector3D -args 0.0 0.0 1.0];

# wurzite intrinsic strain axes
# wurzite kp hamiltonians are [a,a,c] (so x axis is a per default)
# but for strain calculation that may be permuted
# values are also used in effective mass calculation for the hydrostatic
# strain offset
set wurzite_strain_axis_1  1
set wurzite_strain_axis_2  2
set wurzite_strain_axis_3  0

