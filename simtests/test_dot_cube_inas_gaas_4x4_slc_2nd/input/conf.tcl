set slc_homogenous_broadening   5e12
set slc_omega_num 1200
# pardiso currently fails with symmetric matrices
$config set "assembly_build_nonsymmetric_matrices" 1
