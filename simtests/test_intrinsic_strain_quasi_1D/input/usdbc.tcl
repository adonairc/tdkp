
set bc_strain [BCDirichlet2DLines -args $geometry]
# generally no movement at bottom
$bc_strain add_dirichlet_line 0.0 0.0 100.0 0.0
# no movement in x direction at left boundary
$bc_strain add_dirichlet_line 0.0 0.0 0.0 250.0 0
# no movement in x direction at right boundary
$bc_strain add_dirichlet_line 100.0 0.0 100.0 250.0 0
$bc_strain prepare

set bc_poisson [BCDirichlet2DLines -args $geometry]
$bc_poisson add_dirichlet_line 0.0 0.0 100.0 0.0
$bc_poisson add_dirichlet_line 0.0 250.0 100.0 250.0
$bc_poisson prepare
