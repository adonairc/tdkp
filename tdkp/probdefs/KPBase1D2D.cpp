// ------------------------------------------------------------
//
// This file is part of tdkp, a simulation tool for nanostrutctures
// of optoelectronics developed at ETH Zurich
//
// (C) 2005-2009 Ratko G. Veprek, ETH Zurich, veprek@iis.ee.ethz.ch
//
// 1) As of 18.6.2009 this code is property of ETH Zurich and must not be
// transferred, modified or used by third parties without appropriate
// licenses issued by authorized agents of ETH Zurich.
//
// 2) Violation of this will result in judicial action according to civil
// and penal law.
//
// 3) Any claim of authorship other than by the author himself is
// strictly forbidden.
//
// 4) The source code must retain the copyright notice, this list
// of conditions and the following disclaimer.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS
// BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
// IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------------------------------------


#include "tdkp/probdefs/KPBase1D2D.h"
#include "tdkp/main/FEMSolverGEVP.h"
#include "tdkp/main/FEMSolverGEVPRemote.h"
#include "tdkp/solvers/RemoteGEVPSolver.h"
#include <omp.h>

namespace tdkp {

KPBase1D2D::KPBase1D2D(const Geometry& geometry_, MaterialDatabase& material_database_)
: KPBase1D2DParent(geometry_, material_database_),
  k_idx_current(0) 
{	
	this->k_transversal    = 0.0;
	this->first_order_ignore_diagonal = false;
	this->upper_energy_limit_set      = false;	
	this->energy_of_last_solution     = 0.0;
	this->energy_offset      = 0.0;
	this->energy_barrier     = 0.0;
	this->energy_barrier_set = false;
}



KPBase1D2D::~KPBase1D2D() {
	this->delete_solutions();
}


void KPBase1D2D::calculate_element_matrices(const Element* elem, complex<double>* lhs, double *rhs, int* internal_idx, int &n) const {

	TDKP_ASSERT(this->kp_matrices.size() <= (unsigned)this->material_db.get_num_materials(),"this->properties.size() <= (unsigned)this->material_db.get_num_materials()");

	double lstiff[3][3][Element::max_num_nodes][Element::max_num_nodes]; /* stiff matrix */
	double lfirst[3][Element::max_num_nodes][Element::max_num_nodes];    /* first order matrix */
	double lmass[Element::max_num_nodes][Element::max_num_nodes]; 	     /* mass matrix */
	double lmass_nodal[Element::max_num_nodes][Element::max_num_nodes][Element::max_num_nodes];  /* nodal mass matrix */
	int    nnode;              /* number of nodes */
	int    neq;                /* number of kp equations */
	int    lsize; 			   /* size of lhs, rhs */
	int    nsparse;            /* number of sparse interaction matrix elements */
	cplx   i(0.0, 1.0);
	bool   local_first_order = this->first_order_terms_exist; // whether first order term exist

	KPMatrixBase* mat = this->kp_matrices[elem->get_region().get_material().get_id()];
	n     = 0;
	nnode = (signed)elem->get_num_nodes();
	neq   = this->get_num_equations_per_node();
	lsize = nnode * nnode * neq * neq;

	this->get_node_sparsity_pattern(nsparse);

	TDKP_BOUNDS_ASSERT(mat->ready(), "kp matrix not ready");
	TDKP_BOUNDS_ASSERT(nnode <= Element::max_num_nodes, "nnode <= Element::max_num_nodes (lhs rhs working arrays to small ...)");
	TDKP_ASSERT((signed)mat->get_second_order(0,0).size() == nsparse, "mat->second_order[0][0].size() == nsparse");
	TDKP_ASSERT((signed)mat->get_first_order(KPMatrixBase::op_left, 0).size() == nsparse, "mat->first_order[0].size() == nsparse");
	TDKP_ASSERT((signed)mat->get_zero_order().size() == nsparse, "mat->zero_order.size() == nsparse");
	// --------------------------------------------------
	// set element fields to kp matrix and recalculate
	// --------------------------------------------------
	// in the present way here we calculate the element matrices by integrating analytically.
	// varying constants which change linearly over the element can be replaced by their
	// average over the element.
	// only update if element is considered in calculation
	// in order to have this function to be const (so we can use openmp)
	// we need the kp matrices to be stored locally ...
	vector<cplx> second_order[3][3];
	vector<cplx> first_order[2][3];
	vector<cplx> zero_order[Element::max_num_nodes]; 
			
	// only update if element is considered in calculation
	if(nnode > 0) {		
		// ------------------------------------------------
		// first, calculate element wise constant stuff (i assume constant eff mass params)
		// ------------------------------------------------
		if(this->strain_field_set()) {		
			// take average strain	
			const StrainTensor& strain_tensor = get_strain_field().get(elem->get_index_global());			
			for(short ii = 0; ii < 3; ii++) {
				for(short jj = 0; jj < 3; jj++) {
					mat->build_second_order(strain_tensor, ii, jj, second_order[ii][jj]);
				}	
				mat->build_first_order(KPMatrixBase::op_left, strain_tensor, ii, first_order[KPMatrixBase::op_left][ii]);
				mat->build_first_order(KPMatrixBase::op_right, strain_tensor, ii, first_order[KPMatrixBase::op_right][ii]);
			}	
		} else {
			for(short ii = 0; ii < 3; ii++) {
				for(short jj = 0; jj < 3; jj++) {
					second_order[ii][jj] = mat->get_second_order(ii,jj);
				}	
				first_order[KPMatrixBase::op_left][ii]  = mat->get_first_order(KPMatrixBase::op_left, ii);
				first_order[KPMatrixBase::op_right][ii] = mat->get_first_order(KPMatrixBase::op_right, ii);
			}
		}
		// ------------------------------------------------
		// second, calculate nonconstant stuff (such as (2nd order) strain and potential)
		// ------------------------------------------------															
		for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
			double potential = 0.0;
			if(this->potential_energy_field_set()) {			
				potential = get_potential_energy_field().get_node_value(elem->get_node(mm).get_index_global());							
			}
			if(this->strain_field_set()) {							
				const StrainTensor& strain_tensor = get_strain_field().get_nodal_strain(elem->get_index_global(), mm);			
				mat->build_zero_order(strain_tensor, potential, zero_order[mm]);
			} else {
				zero_order[mm] = mat->get_zero_order();
				if(potential != 0.0) {
					mat->add_zero_order_energy(potential,zero_order[mm]);
				}			
			}
		}	
	}
		
	// ------------------------------------------------------
	// workspace used to store parametric dependence on k_t when calculating element matrices
	// ------------------------------------------------------
	/** workspace for first order terms standing LEFT of the differential operator */	
	vector<cplx> first_order_workspace_left [2];
	/** workspace for first order terms standing RIGHT of the differential operator (and thus will be partially integrated) */
	vector<cplx> first_order_workspace_right[2];
	vector<cplx> zero_order_workspace[Element::max_num_nodes];
		
	first_order_workspace_left[0].resize(zero_order[0].size());
	first_order_workspace_left[1].resize(zero_order[0].size());
	first_order_workspace_right[0].resize(zero_order[0].size());
	first_order_workspace_right[1].resize(zero_order[0].size());
	for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
		zero_order_workspace[mm].resize(zero_order[0].size());
	}
		

	// ------------------------------------------------------------------------------
	// perform transformation for k_transversal
	// ------------------------------------------------------------------------------
	if(this->geometry.get_dimension() == 2) {		 	
		// quantum wire: we assume quantization is in x,y plane, free direction is in
		// z direction. so, H(2)iz * k_t -> is added to H(1)i
		// and H(2)zz -> H0 and H(1)z -> H0
		// first H(2)iz -> H(1)i				
		cplx kz = - this->k_transversal / i;	
		for(int ss = 0; ss < nsparse; ss++) {
			for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
				// first order terms (assuming that they are standing left of the differential operator)
				// get the second order Hzi terms, minus sign is due to partial integration of second order terms
				first_order_workspace_left[aa][ss] = first_order[KPMatrixBase::op_left][aa][ss] - (second_order[D_DZ][aa][ss] * kz);
				// first order terms standing right of the differential operator get the
				// second order term Hiz terms (and we must ensure that the integral is partially integrated in 
				// order to get rid of the delta function that else would arise, therefore no minus sign
				// in the second order terms and a minus sign in the first order term)
				first_order_workspace_right[aa][ss] = - first_order[KPMatrixBase::op_right][aa][ss] + second_order[aa][D_DZ][ss] * kz;
			}
		}
		for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
			for(int ss = 0; ss < nsparse; ss++) {			
				// H(2)zz -> H(0) // the minus sign is because H(2) is the
				// partially integrated H(2), means - d^2/dx^2 -> d/dx d/dx
				zero_order_workspace[mm][ss] = zero_order[mm][ss] - second_order[D_DZ][D_DZ][ss] * kz * kz;
				// H(1)z -> H(0)
				zero_order_workspace[mm][ss] += (first_order[KPMatrixBase::op_left][D_DZ][ss] + first_order[KPMatrixBase::op_right][D_DZ][ss]) * kz;
			}			
		}
		local_first_order = true;
	} else if(this->geometry.get_dimension() == 1) {
		cplx kz = - this->k_transversal / i;
		for(int ss = 0; ss < nsparse; ss++) {
			for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {						
				// first order terms (assuming that they are standing left of the differential operator)
				// get the second order Hzi terms, minus sign is due to partial integration of second order terms
				first_order_workspace_left[aa][ss] =  first_order[KPMatrixBase::op_left][aa][ss] - (second_order[D_DZ][aa][ss] * kz);
				// first order terms standing right of the differential operator get the
				// second order term Hiz terms (and we must ensure that the integral is partially integrated in 
				// order to get rid of the delta function that else would arise, therefore no minus sign
				// in second order but a minus sign for first order terms!)
				first_order_workspace_right[aa][ss] = - first_order[KPMatrixBase::op_right][aa][ss] + second_order[aa][D_DZ][ss] * kz;
			}
		}
		for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
			for(int ss = 0; ss < nsparse; ss++) {
				// H(2)zz -> H(0) // the minus sign is because H(2) is the
				// partially integrated H(2), means - d^2/dx^2 -> d/dx d/dx
				zero_order_workspace[mm][ss] = zero_order[mm][ss] - second_order[D_DZ][D_DZ][ss] * kz * kz;
				// H(1)z -> H(0)
				zero_order_workspace[mm][ss] += (first_order[KPMatrixBase::op_left][D_DZ][ss] + first_order[KPMatrixBase::op_right][D_DZ][ss]) * kz;
			}			
		} 
		local_first_order = true;		
	} else {
		TDKP_GENERAL_EXCEPTION("invalid dimension");
	}

	// ------------------------------------------------------------------------------
	// build node_internal_indices (nonzero nodes, which we really calculate)
	// ------------------------------------------------------------------------------
	for(int ii = 0; ii < nnode; ii++) {
		if(elem->get_node(ii).get_index_internal() != -1) {
			internal_idx[n++] = ii;
		}			
	}
																
	// ------------------------------------------------------------------------------
	// set lhs and rhs to zero
	// ------------------------------------------------------------------------------
	for(int ii = 0; ii < lsize; ii++) {
		lhs[ii] = 0.0;
		rhs[ii] = 0.0;
	}		
	// ------------------------------------------------------------------------------
	// calculate local stiffness and mass matrix for nonzero indices
	// ------------------------------------------------------------------------------	
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {
			// build stiffness matrix for all derivatives
			for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
				for(unsigned int bb = 0; bb < this->geometry.get_dimension(); bb++) {					
					lstiff[aa][bb][ii][jj] = elem->get_element_integral_2nd_order(aa, internal_idx[ii], bb, internal_idx[jj]);
				}
			}
			lmass[ii][jj] = elem->get_element_integral_0th_order(internal_idx[ii],internal_idx[jj]);
			for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
				lmass_nodal[mm][ii][jj] = elem->get_element_integral_0th_order_nodal_data(mm, internal_idx[ii],internal_idx[jj]);				
			}
		}
	}	
	if(local_first_order) {		
		for(int ii = 0; ii < n; ii++) {
			for(int jj = 0; jj < n; jj++) {
				for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
					lfirst[aa][ii][jj] = elem->get_element_integral_1st_order(aa, internal_idx[ii],internal_idx[jj]);
				}
			}
		}
	}
	// ------------------------------------------------------------------------------
	// assemble ....
	// ------------------------------------------------------------------------------
	int offset = 0;	
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {					
			offset = (ii * n + jj) * nsparse;
			for(int ss = 0; ss < nsparse; ss++) {
				// assemble "stiffness" contributions
				for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
					for(unsigned int bb = 0; bb < this->geometry.get_dimension(); bb++) {
						lhs[offset + ss] += second_order[aa][bb][ss] *  lstiff[aa][bb][ii][jj];
					}
				}
				// assemble "mass" contributions
				for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
					lhs[offset + ss] += zero_order_workspace[mm][ss] * lmass_nodal[mm][ii][jj];
				}
				
				// assemble first order contributions
				if(local_first_order) { 					
					for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
						// for the left term we have v H dx u, where v is the test function
						// therefore, we need the integral  int Ni dx Nj, but the 
						// elements give dx Ni * Nj, therefore swap here jj and ii
						lhs[offset + ss] += first_order_workspace_left[aa][ss] * lfirst[aa][jj][ii];
						// [rv] no delta function at interface code
						// minus sign comes from partial integration (but is included in the first_order_workspace data)
						lhs[offset + ss] += first_order_workspace_right[aa][ss] * lfirst[aa][ii][jj]; // switched ii jj
					}
				}
			}
			rhs[ii * n + jj] = lmass[ii][jj];
		}
	}
	
}

void KPBase1D2D::add_solution(cplx solution_value, const cplx* solution_vector, int length) {

	// shift energy back
	double eshift   = this->get_energy_shift();
	solution_value -= eshift;
	
	// calculate length of solution vector
	int vector_length = this->geometry.get_num_nonzero_nodes() * this->get_num_equations_per_node();
	TDKP_ASSERT(length == vector_length, "invalid length of solution vector");
	int next = this->solution_cache.size();
	
	// inform user
	ostringstream sout;
	sout << "KPBase1D2D: added solution " << next << " with ev: " << solution_value;
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str()); 
	// write to cache
	this->solution_cache.push_back(vector<cplx>(vector_length + 1));	
	TDKP_ASSERT((signed)this->solution_cache[next].size() == vector_length + 1, "solution cache not long enough");
	this->solution_cache[next][0] = solution_value; // store the energy at first place
	// copy vector and normalize it
#pragma omp parallel for default(shared) schedule(static, 5000)	
	for(int ii = 0; ii < vector_length; ii++) {
		this->solution_cache[next][ii+1] = solution_vector[ii];
	}
	this->normalize_solution(&this->solution_cache[next][1]);

}

void KPBase1D2D::display_solution_info() const {
	ostringstream sout;
	switch(this->bandstructures.size()) {
		case 0:
			sout << "there is no solution available!";
			break;
		case 1:
			sout << "there is one solution available:\n";
			break;
		default:
			sout << "there are " << this->bandstructures.size() << " solutions available:\n";
			break;
	}
	Logger::get_instance()->emit(LOG_INFO, sout.str());
}


DomainMaster KPBase1D2D::build_kspace_domain(int num_subbands, double kmin, double kmax, int num_k_values) const {

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");
	}
	DomainMaster domain;	
	if(this->geometry.get_dimension() == 1) {
		// quantum well in radial approximation
		if(num_k_values == 1) {
			// just a point
			domain.add_node(new DomainNodeSingularPoint(true, 0.0, kmin));
			domain.update();
		} else {
			create_2D_domain_radial(domain, kmin, kmax, num_k_values);		
		}
	} else if(this->geometry.get_dimension() == 2) {		
		if(num_k_values == 1) {
			// just a point
			domain.add_node(new DomainNodeSingularPoint(false, 0.0, kmin));
			domain.update();			
		} else {		
			// quantum wire, approximate trapezoidal
			create_1D_domain_wire_bands(domain, kmin, kmax, num_k_values);	
		}		
	} else {
		TDKP_GENERAL_EXCEPTION("dimension of problem must be 1 or 2  ... ");	
	}
	return domain;
}
void KPBase1D2D::solve(int num_subbands, double kmin, double kmax, int num_k_values) {				
	this->solve(num_subbands, build_kspace_domain(num_subbands, kmin, kmax, num_k_values));
}

void KPBase1D2D::solving_preinform_user(int num_subbands, const DomainMaster& domain) const {
	// -----------------------------------------------
	// output to user
	// -----------------------------------------------
  	ostringstream sout;
  	switch(this->geometry.get_dimension()) {
    	case 2:
		  	sout << "KPBase1D2D: solving kp quantum wire equation:\n"
		  	     << "transversal direction: " << this->k_directions[D_DZ] << "\n"
				 << "k-space-range:         from " << domain.get_first_point().get_coord_abs() 
				 << " to " << domain.get_last_point().get_coord_abs() << " [1/nm]\n";
     		break;
    	case 1:
		  	sout << "KPBase1D2D: solving kp quantum well equation ";
		  	if(domain.radial()) {
		  		sout << "in the radial approximation:\n"
		  	         << "transversal direction: " << this->k_directions[D_DZ] << "\n"
					 << "k-space-range:         from " << domain.get_first_point().get_coord_abs() 
					 << " to " << domain.get_last_point().get_coord_abs() << " [1/nm]\n";		  	         
		  	} else {
		  		sout << "for the full plane. the domain axes are given by:\n"
		  		     << "x-axis:                " << this->k_directions[D_DZ] << "\n"
		  		     << "y-axis:                " << this->k_directions[D_DY] << "\n";
		  	}
      		break;
    	default:
      		TDKP_GENERAL_EXCEPTION("invalid type of geometry dimension for this class");
  	}
	sout << "number of subbands:    " << num_subbands << "\n"
		 << "type of subbands:      " << this->get_solution_type() << "\n"
		 << "strain field:          " << (strain_field_set() ? "included" : "not set") << "\n" 
		 << "potential energy:      " << (potential_energy_field_set() ? "included" : "not set");  
		 
	Logger::get_instance()->emit(LOG_INFO, sout.str());
}

void KPBase1D2D::solve(int num_subbands, const DomainMaster& domain) {

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");
	}

    this->solving_preinform_user(num_subbands, domain);  	 
  	
	// -------------------------------------------------------------
	// prepare bandstructure object
	// -------------------------------------------------------------
	BandstructureDomain<cplx>* band = new BandstructureDomain<cplx>(this->get_num_equations_per_node(), num_subbands, this->geometry.get_num_nodes(), domain);
  	this->bandstructures.push_back(band);
	
	// -------------------------------------------------------------
	// reset energy of last solution to the place where we expect the 
	// new solutions (bandedges is not a good idea when potential 
	// energy fields are present
	// -------------------------------------------------------------
	this->k_idx_current = 0; // nasty, must be 0 before we get the right energy shift (in subsequent calc., the k idx may still be set > 0)
	if(this->use_user_defined_energy_guess()) {
		TDKP_ASSERT(this->energy_guess.size() > 0, "this->energy_guess.size() > 0");
		this->energy_of_last_solution = this->energy_guess[0];
	} else {
		this->energy_of_last_solution = this->get_minimum_bandedges();
	}  

	// -------------------------------------------------------------
	// create fem assembler (parallel or non-parallel)
	// -------------------------------------------------------------
	bool remote_parallel = Configuration::get_instance()->get("desired_eigenvalue_solver") == 4.0;
	// vector to store parallel controllers 
	vector<RemoteGEVPSolverController*> remote_controllers;
	// vector to store threads performing the solve process of the controllers	
	vector<pthread_t> threads;
	// standard solver object
	FEMSolverGEVP<complex<double>, complex<double>, double>* solver        = 0;
	// specialiced remote solver
	FEMSolverGEVPRemote*                              remote_solver = 0;
	if(remote_parallel) {
		// check that user defined energy guesses are set
		if(Configuration::get_instance()->get("assembly_parallel_disable_automatic_energy_guess") == 0.0) {
			this->auto_prepare_energy_guess(domain);
		}
		TDKP_ASSERT(energy_guess_set.size() >= domain.get_number_of_points(), "parallel k-space calculation requires you to give energy guesses for all k-points!");
		for(unsigned int ii = 0; ii < domain.get_number_of_points(); ii++) {
			TDKP_ASSERT(energy_guess_set[ii], "energy guess for kidx " << ii << " is missing but required for parallel k-space calculation");	
		}		
		// create remote parallel fem solver
		remote_solver = new FEMSolverGEVPRemote(geometry, *this);
		solver = remote_solver;
		if(this->upper_energy_limit_set) {
			TDKP_LOGMSG(LOG_WARN, "you have requested upper limit termination while using parallel solvers where upper limit is not avialable. i will calculate full bandstructure.");	
		} 		
	} else {
		solver = new FEMSolverGEVP<complex<double>, complex<double>, double>(geometry, *this);	
	}
	
	
	// -------------------------------------------------------------
	// create matrix structure
	// -------------------------------------------------------------
	solver->create_matrix_structures();	
	if(this->get_solution_type() == electrons) {
		solver->set_ordering(ascending);
	} else {
		solver->set_ordering(descending);	
	}
	if(Logger::get_instance()->get_level() == LOG_INFO) {
		Logger::get_instance()->init_progress_bar("KPBase1D2D: calculating bandstructure in k space, using points: ", domain.get_number_of_points());
	}	
		
	this->prepare();
	
	// -------------------------------------------------------------
	// tell the matrices to shut up
	// -------------------------------------------------------------
	for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {	
		this->kp_matrices[ii]->set_output_surpression(true);		
	}
		
	int next    = 0;
	int current = 0;
	// -------------------------------------------------------------
	// loop over domain points
	// -------------------------------------------------------------
  	for(int kk = 0; kk < (signed)domain.get_number_of_points(); kk++) {
		adaptive_omp_threading();  		  		
  		this->k_idx_current = kk; // setting current index to store bandstructure at the right place
  		
  		// --------------------------------------------------------
  		// set the length value of the transversal direction
  		// --------------------------------------------------------
  		this->k_transversal = domain.get_point(kk).get_coord_abs();  		

		// --------------------------------------------------------
		// if we have a non-radial well, rotate the kp matrices for
		// the next k point
		// --------------------------------------------------------
		if(geometry.get_dimension() == 1 && !domain.radial()) {
			// ----------------------------------------------
			// kz is the free direction in the assembly
			// so, kz' = kz * cos(phi) - ky * sin(phi)
			// and this equals
			//         = kz * dir_x - ky * dir_y
			// where dir is the unit direction of the point 			
			// ----------------------------------------------
			double dir_x = domain.get_point(kk).get_coord(0) / this->k_transversal;
			double dir_y = domain.get_point(kk).get_coord(1) / this->k_transversal;
			Vector3D new_kz = dir_x * this->k_directions[D_DZ] - dir_y * this->k_directions[D_DY];
			Vector3D new_ky = dir_x * this->k_directions[D_DY] + dir_y * this->k_directions[D_DZ];
			RMatrix<double> tmp_rotation(3,3);
			tmp_rotation.set_row(0, this->k_directions[D_DX].get_all());
			tmp_rotation.set_row(1, new_ky.get_all());
			tmp_rotation.set_row(2, new_kz.get_all());
			// ----------------------------------------------
			// set rotation to kp matrices
			// ----------------------------------------------
			for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {	
				this->kp_matrices[ii]->set_rotation(tmp_rotation);				
			}
  			// --------------------------------------------------------
			// information for user
			// --------------------------------------------------------
			ostringstream sout;			
			sout << "KPBase1D2D: iteration " << kk + 1 << " for kx,ky = (" << domain.get_point(kk).get_coord(0)
			     << ", " << domain.get_point(kk).get_coord(1) <<") "
		     	 << " [1/nm]. expecting kp states at " << - this->get_energy_shift() << " [eV] ";
			if(this->use_user_defined_energy_guess()) {
				sout << "(user defined) ";	
			}		     
			sout << "for the next cycle.";
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());									   	
		} else {
  			// --------------------------------------------------------
			// information for user
			// --------------------------------------------------------			
			ostringstream sout;
			sout << "KPBase1D2D: iteration " << kk + 1 << " of " << domain.get_number_of_points() << " for k = " << this->k_transversal 
		     	 << " [1/nm]. expecting kp states at " << - this->get_energy_shift() << " [eV] ";
			if(this->use_user_defined_energy_guess()) {
				sout << "(user defined) ";	
			}		     
			sout << "for the next cycle.";
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());				
		}
		
		// --------------------------------------------------------		
    	// prepare kp matrices for the next iteration 
    	// (give them the new energy shift)
    	// --------------------------------------------------------		    	
		for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {						
			this->kp_matrices[ii]->set_energy_shift(this->get_energy_shift());		
			this->kp_matrices[ii]->calculate();
			TDKP_ASSERT(this->kp_matrices[ii]->surpress_output(), "somehow the kp matrices where deleted ... this is an inconsistency in prepare and it should happen ... ");			
		}  		
  		
		// --------------------------------------------------------		  		
  		// assemble system  		
  		// --------------------------------------------------------		
    	solver->assemble_system();
    	    	    	    
    	// --------------------------------------------------------
    	// solver serial or create parallel solvers
		// --------------------------------------------------------
		if(remote_parallel) {	
#ifndef NOREMOTESOLVER			
			TDKP_ASSERT(remote_solver != 0, "");			
			// create parallel solve controller
			remote_controllers.push_back(remote_solver->get_remote_solve_controller(num_subbands, this->get_problem_type()));
			// create thread executing solver
			pthread_t thread;			
        	int rc = pthread_create(&thread, NULL, parallel_controller_solve, remote_controllers.back());
			if(rc) {
				TDKP_GENERAL_EXCEPTION("thread creation failed");
			}
        	threads.push_back(thread);
#else
			TDKP_GENERAL_EXCEPTION("remote solver functionality not provided (tdkp was compiled with NOREMOTESOLVER)");
#endif	
				
		} else {    	
	    	solver->solve_system(num_subbands);
	    	
	    	// --------------------------------------------------------		
	    	// add solutions from cache to bs object    	
	    	// --------------------------------------------------------		
	    	this->add_cache_to_bandstructure_object(band, kk);
	    			  	
	   		// --------------------------------------------------------		
			// check if the bandstructure passed the upper limit
			// --------------------------------------------------------		
			if(this->upper_energy_limit_set) {
				bool do_a_break = false;
				if(this->get_solution_type() == electrons) {
					// energy of lowest band above the limit? 
					if(band->get_energy(kk,0).real() > this->upper_energy_limit) {
						do_a_break = true;	
					}	
				} else {				
					// energy of lowest band below the limit? 
					if(band->get_energy(kk,0).real() < this->upper_energy_limit) {
						do_a_break = true;	
					}		
				}
				if(do_a_break) {
					if(Logger::get_instance()->get_level() == LOG_INFO) {
	  					Logger::get_instance()->end_progress_bar();
	  				}	
	  				ostringstream sout;
	  				sout << "KPBase1D2D: the most bound band has now passed the energy "
	  				     << "limit of " << this->upper_energy_limit << " you provided. therefore "
	  				     << "i stop the calculation now. it is on you to handle the "
	  				     << "discontinuity in the bandstructure!";
					Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
					break;
				}					  	
			}
		}		  	
    	
    	if(current == next) {
    		if(Logger::get_instance()->get_level() == LOG_INFO) {
    			next = Logger::get_instance()->set_progress_bar(current, domain.get_number_of_points());
    		}
    	}
    	current++;
   		solver->reset_matrix_to_zero();   		   		
  	}
  	  	  	
  	if(Logger::get_instance()->get_level() == LOG_INFO) {
  		Logger::get_instance()->end_progress_bar();
  	}
  	
  	// -------------------------------------------------------------
  	// receive remote results
  	// -------------------------------------------------------------
  	if(remote_parallel) {
#ifndef NOREMOTESOLVER			
  		TDKP_ASSERT(domain.get_number_of_points() == threads.size(),"");
  		TDKP_LOGMSG(LOG_INFO_DEVEL2, "KPBase1D2D: waiting for " << threads.size() << " threads to finish");  		  		
  		for(int kk = 0; kk < (signed)domain.get_number_of_points(); kk++) {
			this->k_idx_current = kk;
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "KPBase1D2D: waiting for results from thread " << kk << " (executed on " << remote_controllers[kk]->get_remote_host() << ")");
			// wait for thread to finish
			pthread_join(threads[kk], 0);			
			// set results
			remote_controllers[kk]->set_results_to_problem_object(*this);
			// add cache to bandstructure object		
	    	this->add_cache_to_bandstructure_object(band, kk);			
			// kill object
			delete remote_controllers[kk]; remote_controllers[kk] = 0;							
  		}
  		TDKP_LOGMSG(LOG_INFO_DEVEL2, "KPBase1D2D: all remote processes terminated and data received.");
#else
		TDKP_GENERAL_EXCEPTION("remote solver functionality not provided (tdkp was compiled with NOREMOTESOLVER)");
#endif	
  	}
  	
  	// -------------------------------------------------------------
	// so, the matrices are allowed to talk again ...
	// -------------------------------------------------------------
	for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {	
		this->kp_matrices[ii]->set_output_surpression(false);
	}
	
	// delete solver object
	delete solver;
}

/** thread function: call solve on controller */
void* KPBase1D2D::parallel_controller_solve(void *arg) {
#ifndef NOREMOTESOLVER			
	RemoteGEVPSolverController* ctrl = (RemoteGEVPSolverController*) arg;
	int tries = 3;
	while(tries-- > 0) { 
		if(ctrl->solve()) {
			break;	
		}
		TDKP_LOGMSG(LOG_WARN, "KPBase1D2D: remote solver failed! trying again!");
	}
	if(tries == 0) {
		TDKP_GENERAL_EXCEPTION("remote solver failed again! giving up!");
	}
	pthread_exit(NULL);	
#else
		TDKP_GENERAL_EXCEPTION("remote solver functionality not provided (tdkp was compiled with NOREMOTESOLVER)");
#endif	
}

/** try to guess the dispersion energy for the whole domain 
 * 
 * does not override present values and expects kidx 0 to be
 * the zero point energy guess 
 */
void KPBase1D2D::auto_prepare_energy_guess(const DomainMaster& domain) {
	
	// ---------------------------------------------
	// check that kidx 0 is k = 0 and Eguess(0) is set
	// ---------------------------------------------
	if(tdkp_math::abs(domain.get_point(0).get_coord_abs()) > 1.0e-3) {
		TDKP_LOGMSG(LOG_WARN, "KPBase1D2D: auto prepare energy guess expects kidx 0 to be k = 0 point because it needs the 0 point as zero point energy guess");
	}
	TDKP_ASSERT(energy_guess_set.size() > 0 && energy_guess_set[0], "auto_prepare_energy_guess requires you to set the energy guess for k = 0"); 
		
	// ---------------------------------------------
	// get polynomial
	// ---------------------------------------------	
	double pc_c = energy_guess[0];
	double pc_b, pc_a;
	if(this->get_solution_type() == electrons) {
		pc_a = Configuration::get_instance()->get("assembly_parallel_energy_guess_cb_poly_coeff_a");
		pc_b = Configuration::get_instance()->get("assembly_parallel_energy_guess_cb_poly_coeff_b");	
	} else {
		pc_a = Configuration::get_instance()->get("assembly_parallel_energy_guess_vb_poly_coeff_a");
		pc_b = Configuration::get_instance()->get("assembly_parallel_energy_guess_vb_poly_coeff_b");			
	}

	// ---------------------------------------------
	// set energy guesses
	// ---------------------------------------------
	bool user_warned = false;	
	for(unsigned int ii = 1; ii < domain.get_number_of_points(); ii++) {
		if(!user_warned && energy_guess_set.size() > ii && energy_guess_set[ii]) {
			TDKP_LOGMSG(LOG_WARN, "KPBase1D2D: auto_prepare_energy_guess overrides old energy guesses");
			user_warned = true;	
		} 
		double guess = pc_a * domain.get_point(ii).get_coord_abs() * domain.get_point(ii).get_coord_abs()
		             + pc_b * domain.get_point(ii).get_coord_abs()
		             + pc_c;
		this->set_energy_guess(ii, guess);
	}		
		
}

/** write results into bandstructure object
 * 
 * the fem solver does solve the eigensystems and passes
 * the values via add_solution to the problem class
 * 
 * unfortunately i first have to renormalize the solution to its quantum
 * mechanical norm. therefore i have to first copy it, renormalize it
 * and then copy to the eigensolution object using global node indices.
 * finally the eigensolution object is added to the band structure object ...
 *  
 */
void KPBase1D2D::add_cache_to_bandstructure_object(Bandstructure<cplx>* band, int kidx) {
	
	int neq = this->get_num_equations_per_node();
		
	// prepare to store the energy closest to the bandedge (for the next iteration)
	bool ascending_order;
	const double& (*closest)(const double& a, const double& b);
	if(this->get_solution_type() == electrons) {		
		closest         = min;	
		ascending_order = true;
	} else {
		closest = max;
		ascending_order = false;	
	}
	
	//this->energy_of_last_solution = this->solution_cache[0][0].real();
	
	// determine correct ordering of our solutions as they may come
	// disordered from arpack due to the fact that we were looking
	// for eigenvalues next to a specific sigma
	vector<double>       values;
	vector<unsigned int> indexes;
	for(unsigned int ii = 0; ii < this->solution_cache.size(); ii++) {
		values.push_back(this->solution_cache[ii][0].real());
		indexes.push_back(ii);
	}

	tdkp_math::tracked_sort<unsigned int,double>(indexes.begin(), indexes.end(), values.begin(), values.end(), ascending_order);	
	
	// -------------------------------------------------------------------
	// if the energy barrier is set, shift states on the wrong side of the
	// barrier to the trash states at the end of the result array
	// -------------------------------------------------------------------
	if(this->energy_barrier_set) {
		unsigned int ii   = 0;
		double tmp_energy = 0;
		ostringstream eout;
		eout << "KPBase1D2D: barrier (" << this->energy_barrier << ") sorting kicked the following states to the end:\n";
		// swap states with energies on the other side of the
		// barrier to the end
		unsigned int last = indexes.size(); 
		while(ii < last) {
			tmp_energy = this->solution_cache[indexes[ii]][0].real();
			if(( ascending_order && tmp_energy < this->energy_barrier) || 
			   (!ascending_order && tmp_energy > this->energy_barrier)) {
			   	eout << "state " << setw(3) << indexes[ii] << " with energy " 
			   	     << setw(10) << tmp_energy << " [eV]\n"; 			   			   	
				indexes.push_back(indexes[ii]);				
				indexes.erase(indexes.begin() + ii);
				last--;		
			} else {
				ii++;	
			}
		}
		if(last != indexes.size()) {
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, eout.str());	
		}
	}

	int min_dist_idx = -1; 
	double min_dist_value = 1.0e100; // should be big enough ;-)			
	// for all solutions
	for(unsigned int jj = 0; jj < this->solution_cache.size(); jj++) {
		// use resorted bands
		int ii = indexes[jj];
		// find value closest to last solution
		if(!this->energy_barrier_set || ((ascending_order && this->solution_cache[ii][0].real() >= this->energy_barrier) || (!ascending_order && this->solution_cache[ii][0].real() <= this->energy_barrier))) {			         
    		double tmp_dist = fabs(this->solution_cache[ii][0].real() - this->energy_of_last_solution);
			if(min_dist_idx == -1 || tmp_dist < min_dist_value) {
				min_dist_value = tmp_dist;
				min_dist_idx   = ii;		
			}        														     
		} 	
		   
		EigenSolution<cplx>* tmp_sol = new EigenSolution<cplx>(this->geometry.get_num_nodes(), neq);
		tmp_sol->set_energy(this->solution_cache[ii][0]);	
		
		// regroup node values by index global
#pragma omp parallel for default(shared)
		for(int vv = 0; vv < (signed)this->geometry.get_num_nodes(); vv++) {
			const Node& node = this->geometry.get_node(vv);
			if(node.get_index_internal() == -1) {
				// for every equation
				for(int nn = 0; nn < neq; nn++) {
					tmp_sol->get_node_value(node.get_index_global(), nn) = 0.0;
				}
			} else {
				for(int nn = 0; nn < neq; nn++) {
					// the + 1 at the end is due to the fact that at pos 0 the energy is stored in the solution cache
					tmp_sol->get_node_value(node.get_index_global(), nn) = this->solution_cache[ii][node.get_index_internal() * neq + nn + 1];
				}
			}
		}
		band->add_eigensolution(kidx, jj, tmp_sol);
	}
	TDKP_ASSERT(min_dist_idx >= 0, "algorithm could not find a state close enough to previous result!");
	this->energy_of_last_solution = this->solution_cache[min_dist_idx][0].real();	
	this->solution_cache.clear();
	this->solution_cache.resize(0);

}

const int* KPBase1D2D::get_node_sparsity_pattern(int &num) const {
	num = this->sparsity_num;
	return this->sparsity_copy;
}


/** get bandedges 
 * 
 * particularly interesting when strain is present
 * @param filename name of file where to write
 */
 /*
StdElementData<double>* KPBase1D2D::get_bandedges() throw(Exception*) {

	if(this->kp_matrices.size() == 0) {
		TDKP_GENERAL_EXCEPTION("object is not prepared! call .prepare() before getting bandedges");
	}
	// -------------------------------------------------------------------
	// create dfDatasets
	// as we store element wise we have one dset per region
	// -------------------------------------------------------------------
	const int num_equations = this->get_num_equations_per_node();	
	StdElementData<double>* ret = new StdElementData<double>(this->get_num_equations_per_node(), this->geometry.get_num_elements());
	for(int bb = 0; bb < num_equations; bb++) {
		ostringstream sout;
		sout << "BandEdge" << bb;
		ret->set_identifier(bb, sout.str().c_str());		
  	}	
	// -------------------------------------------------------------------
	// calculate element matrices and extract bandedges
	// -------------------------------------------------------------------
	double shift                    = this->get_energy_shift(); // bandedges may have been shifted
	KPMatrixBase*         mat       = 0; 	
	RMatrix<cplx>*        bg_matrix = 0;
	vector<cplx> eigenvalues;
	vector<cplx> eigenvectors;	
	
	int next, current;
	next = current = 0;
	Logger::get_instance()->init_progress_bar("calculating the bandedge in all elements", this->geometry.get_num_elements());
			
	for(Geometry::element_const_iterator eit = this->geometry.elements_begin(); eit != this->geometry.elements_end(); eit++) {
		// get kp matrix and set strains and potential energy
		mat = this->kp_matrices[(*eit)->get_region().get_material().get_id()];	
		if(this->strain_field_set()) {
			mat->set_strains(get_strain_field().get((*eit)->get_index_global()));
		}
		if(this->potential_energy_field_set()) {
			mat->set_potential(get_potential_energy_field().get_element_value((*eit)->get_index_global(), 0));			
		}
		// recalculate and get zero order matrix
		mat->calculate();	
		bg_matrix = mat->get_zero_order_matrix();
		
		// get eigenvalues of that zero order matrix
		RMatrix<cplx>::get_eigensystem(*bg_matrix, eigenvalues, eigenvectors);		
		for(int bb = 0; bb < num_equations; bb++) {
			ret->set_element_value((*eit)->get_index_global(), bb, eigenvalues[bb].real() - shift);
		}
		delete bg_matrix; bg_matrix = 0;
		if(next == current) {
			next = Logger::get_instance()->set_progress_bar(current, this->geometry.get_num_elements());
		}
		current++;
	}  	
	Logger::get_instance()->end_progress_bar();
	return ret;
	
}
*/


/** set energy guess
 *
 * @param kidx   k index of the energy guess (at what energy do you want me to set it)
 * @param energy set energy where you expect the solutions
 */
void KPBase1D2D::set_energy_guess(unsigned int kidx, double energy) {
	
	TDKP_ASSERT(kidx < 10000, "sorry, but kidx must be < 10000 (some arbitrary value)");
	
	if(kidx >= this->energy_guess.size()) {
		this->energy_guess.resize(kidx + 5);
		this->energy_guess_set.resize(kidx + 5,false);	
	}

	this->energy_guess[kidx]     = energy;
	this->energy_guess_set[kidx] = true;
	
}

/** unset energy guess */
void KPBase1D2D::remove_energy_guess(unsigned int kidx) {
	if(kidx < this->energy_guess.size()) {
		this->energy_guess_set[kidx] = false;
	}	
}

void KPBase1D2D::remove_all_energy_guesses() {
	for(vector<bool>::iterator it = this->energy_guess_set.begin(); 
	    it != this->energy_guess_set.end(); it++) {
		*it = false;    	
    }    
}

/** return the energy range where we expect the states to be 
 * 
 * why do I shift the spectrum? because the eigensolver needs 
 * (A - sigma M)^-1 to solve for the desired eigenvalues
 * therefore i just set set sigma to zero and shift the A
 * to the appropriate place ...
 */
double KPBase1D2D::get_energy_shift() const {
	ostringstream sout;
	if(this->use_user_defined_energy_guess()) {		
		return - this->energy_guess[this->k_idx_current];
	} else { 
		return - (this->energy_of_last_solution + this->energy_offset);
	}
}

/** are we using a user defined energy guess at the moment? */
bool KPBase1D2D::use_user_defined_energy_guess() const {
	if(this->k_idx_current < this->energy_guess_set.size() && this->energy_guess_set[this->k_idx_current]) {
		return true;	
	} else {
		return false;	
	}	
}

/** basic init process which must be called by derived class constructors */ 
void KPBase1D2D::init_base1D2D() {
	KPMatrixBase* tmp_num          = this->get_matrix();
	this->num_equations_per_node = tmp_num->get_number_of_bands();
	this->first_order_terms_exist  = false;
	// copy sparsity pattern			
	const int* sparsity_pattern = tmp_num->get_sparsity_pattern(this->sparsity_num);
	this->sparsity_copy	= new int[this->sparsity_num * 2];
	for(int ii = 0; ii < (this->sparsity_num * 2); ii++) {
		this->sparsity_copy[ii] = sparsity_pattern[ii];
	}	
	delete tmp_num;
}

/** basic delete process which must be called by derived class destructors */
void KPBase1D2D::delete_base1D2D() {
	this->delete_solutions();	
	if(this->sparsity_copy) {
		delete[] this->sparsity_copy; this->sparsity_copy = 0;
	}
	for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
		delete this->kp_matrices[ii]; this->kp_matrices[ii] = 0;
	}
}


void KPBase1D2D::prepare() {
	
	
	if(this->geometry.get_dimension() != 2 && this->geometry.get_dimension() != 1) {
		TDKP_GENERAL_EXCEPTION("this problem class is only suited for 1D/2D problems");	
	}	
			
	// remove any existing kp matrices
	if(this->kp_matrices.size() > 0) {		
		for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
			delete this->kp_matrices[ii]; this->kp_matrices[ii] = 0;
		}		
		this->kp_matrices.clear();
	}
	
	// create new matrices
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {	
		KPMatrixBase* mat = this->get_matrix();
		mat->set_material(this->material_db.get_material(ii));
		mat->set_material_name(this->material_db.get_material_name(ii));
		mat->set_energy_shift(this->get_energy_shift());
		// attention, if you solve a quantum well in the nonradial approx,
		// the system is rotated for every new point
		mat->set_rotation(this->rotation_matrix);
		mat->calculate();
		// dump kp matrix on request
		if(Configuration::get_instance()->get("output_kpmatrix_dump_to_file") == 1) {
			ostringstream kout;
			kout << "kp" << this->num_equations_per_node << "x" 
			     << this->num_equations_per_node << "_matrix_" 
			     << this->material_db.get_material_name(ii)
			     << ".txt";
			mat->dump(kout.str().c_str()); 
		}		
		this->kp_matrices.push_back(mat);		
	}
			
	// allocate work space for kt dependence
	int nsparse = 0;
	this->get_node_sparsity_pattern(nsparse);

	this->ready = true;
}


/** determine the smallest bandedge (without strain) in present problem
 * 
 * needed for as a starting point where to look for states
 */
double KPBase1D2D::get_minimum_bandedges() const {
	
	const char* key;
	double asc;
	double bandedge;
	
	if(this->get_solution_type() == electrons) {
		key = "conduction_band_edge";
		asc = 1.0;
	} else if(this->get_solution_type() == holes) {
		key = "valence_band_edge";
		asc = -1.0;
	} else {
		TDKP_GENERAL_EXCEPTION("unknown solution type requested");	
	}
	
	if(this->material_db.get_num_materials() < 1) {
		TDKP_GENERAL_EXCEPTION("no materials available");	
	}	
	bandedge = this->material_db.get_material(0)->get(key);	
	for(int ii = 1; ii < this->material_db.get_num_materials(); ii++) {
		if(bandedge * asc > this->material_db.get_material(ii)->get(key) * asc) {
			bandedge = this->material_db.get_material(ii)->get(key);	
		}	
	}
	return bandedge;		
}

/** set energy offset when calculating bandstructure 
 * 
 * the problem is that for small bandgap materials in kp8x8,
 * holes may be nearer located to the first bound electron than
 * the other electrons, therefore the solver finds an electron
 * band and the rest will be hole bands.
 * 
 * therefore when we proceed along the k-space, we need to increase
 * the energy target level from the bandedge up to some reasonable 
 * target. 
 * 
 * the strategy now is to use the energy of the lowest bound state
 * as the initial target and add some offset value which is
 * by default 0.0 but may be changed to some appropriate value 
 */ 
void KPBase1D2D::set_target_energy_offset(double energy_offset_) {
	this->energy_offset = energy_offset_;
}

/** set the energy barrier 
 * 
 * when calculating the bandstructure, sometimes undesired
 * solutions may be found and the whole shifting strategy 
 * would go into the wrong direction. therefore one can set the 
 * energy barrier. solutions with energy of the wrong side of 
 * the barrier (below for electrons, above for holes) will be 
 * counted to the highest band (therefore the highest bands 
 * may just be noise)
 * 
 * but nothing is deleted. the user still may resort by himself
 */

void KPBase1D2D::set_energy_barrier(double energy) {
	this->energy_barrier_set = true;
	this->energy_barrier     = energy;
} 
void KPBase1D2D::remove_energy_barrier() {
	this->energy_barrier_set = false;	
}

/** set axes of the quantized system
 * 
 * the code is capable to calculate the bandstructure of any system
 * oriented in any direction of quantization. in the code itself, it 
 * is just assumed that the quantization of a wire is given in the
 * x/y plane and the transversal direction is the z direction. for a
 * quantum well, the assumption is that the quantization is along the
 * x direction and the transversal direction is the z direction.
 * 
 * the kp matrices can be rotated into any direction and therefore allowing
 * to calculate the bandstructure in any direction. 
 * 
 * there are two functions, one that takes three Vector3D objects and one 
 * only taking two. the former must be used for quantum wires and the latter
 * for quantum wells. 
 * 
 * the way to calculate the rotation axes is very easy. it's the users job,
 * as the rotation matrix is given by 
 *            R = [k_x^T; 
 *                 k_y^T; 
 *                 k_transversal^T]
 * 
 * e.g. x -> [010], y -> [001], transversal -> [100]
 * then R is given by
 *  
 *       [0 1 0]
 *   R = [0 0 1]
 *       [1 0 0]
 * 
 */
void KPBase1D2D::set_axes(const Vector3D& k_transversal, const Vector3D& k_x, const Vector3D& k_y) {
	
	double        dot_product;
	double        determinante;
	ostringstream sout;
	const string  named_directions[] = {string("k_x"), string("k_y"), string("k_transversal")};
	const double  tolerance = Configuration::get_instance()->get("assembly_rotation_matrix_orthogonality_threshold");
		
	if(this->geometry.get_dimension() != 2) {
		TDKP_GENERAL_EXCEPTION("you can not use the function for 2D problems for a non 2D structure");	
	}

	sout << "KPBase1D2D: using the following axes of quantization:\n"
	     << "quantized 1 (x) -> " << k_x << "\n"
	     << "quantized 2 (y) -> " << k_y << "\n"
	     << "transversal     -> " << k_transversal;	      
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str()); 	   
	sout.str("");  
	
	// -----------------------------------------------
	// store vectors and normalize them
	// -----------------------------------------------
	this->k_directions[D_DX] = k_x;
	this->k_directions[D_DY] = k_y;
	this->k_directions[D_DZ] = k_transversal;
	for(int ii = 0; ii < 3; ii++) {
		this->k_directions[ii].normalize();	
	}
	// ------------------------------------------------
	// test vectors (they must be orthogonal!)
	// -----------------------------------------------	
	for(int ii = 0; ii < 2; ii++) {
		dot_product = Vector3D::dot_product(this->k_directions[ii], this->k_directions[ii + 1]);
		if(fabs(dot_product) > tolerance) {
			sout << "KPBase1D2D: your directions " << named_directions[ii] << " and " 
			     << named_directions[ii + 1] << " are not orthogonal. "
			     << "their dot product is " << dot_product << " which is above the "
			     << "tolerance " << tolerance;
			TDKP_GENERAL_EXCEPTION(sout.str()); 
		}	
	}
		
	// ----------------------------------------------
	// and calculate the rotation matrix
	// ----------------------------------------------
	for(int ii = 0; ii < 3; ii++) {
		this->rotation_matrix.set_row(ii, this->k_directions[ii].get_all());
	}
		
	// ----------------------------------------------
	// finally test if the determinante is 1
	// ----------------------------------------------
	determinante = this->rotation_matrix(0,0) * this->rotation_matrix(1,1) * this->rotation_matrix(2,2)
	             + this->rotation_matrix(0,1) * this->rotation_matrix(1,2) * this->rotation_matrix(2,0)
	             + this->rotation_matrix(0,2) * this->rotation_matrix(1,0) * this->rotation_matrix(2,1)
	             - this->rotation_matrix(2,0) * this->rotation_matrix(1,1) * this->rotation_matrix(0,2)
	             - this->rotation_matrix(2,1) * this->rotation_matrix(1,2) * this->rotation_matrix(0,0)
	             - this->rotation_matrix(2,2) * this->rotation_matrix(1,0) * this->rotation_matrix(0,1);

	if(determinante < 0) {
		sout << "KPBase1D2D: the matrix of the resulting rotation given by your directions\n"
		     << "has a determinante of " << determinante << " which is definitely wrong.\n"
		     << "the matrix itself looks that way: " << this->rotation_matrix;
		TDKP_GENERAL_EXCEPTION(sout.str());
	} else if(fabs(determinante - 1.0) > tolerance) {
		sout << "KPBase1D2D: the rotation matrix given by your directions has a determinante of\n"
		     << determinante << " which deviates from 1.0 more than the tolerance\n" 
		     << tolerance << "\n" 
		     << "the matrix itself looks that way: " << this->rotation_matrix;
		TDKP_GENERAL_EXCEPTION(sout.str());		      	
	}	
	sout.str("");
	sout << "KPBase1D2D: the resulting rotation matrix is given by " << this->rotation_matrix;
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str()); 	

}

/** set the quantization axes of a quantum well 
 * 
 * the axis k_q is the direction of the quantized axis
 * k_t is the primary direction of the transversal direction.
 * so, the ktx direction
 * internall, the kty direction is set to kt ^ kq
 * 
 * if you solve for the radial bandstructure, you get the ktx direction,
 * if you solve for a planar domain, you get it in terms of ktx and kty
 */ 
void KPBase1D2D::set_axes(const Vector3D& k_transversal, const Vector3D& k_q) {
	
	double        dot_product;
	ostringstream sout;
	const double  tolerance = Configuration::get_instance()->get("assembly_rotation_matrix_orthogonality_threshold");	
	
	if(this->geometry.get_dimension() != 1) {
		TDKP_GENERAL_EXCEPTION("you can not use the function for 1D problems for a non 1D structure (is " << this->geometry.get_dimension() << ") ");	
	}	
	this->k_directions[D_DX] = k_q;	
	this->k_directions[D_DZ] = k_transversal;
	this->k_directions[D_DY] = Vector3D::cross_product(k_transversal,k_q); //Vector3D(0.0, 0.0, 0.0);
	
	this->k_directions[D_DX].normalize();
	this->k_directions[D_DY].normalize();
	this->k_directions[D_DZ].normalize();	
	
	dot_product = Vector3D::dot_product(this->k_directions[D_DX], this->k_directions[D_DZ]);
	if(fabs(dot_product) > tolerance) {
		sout << "KPBase1D2D: your directions are not orthogonal! their dot product is "
		     << dot_product << " which is more than the tolerance " << tolerance;
		TDKP_GENERAL_EXCEPTION(sout.str());
	}
	
	// ----------------------------------------------
	// and calculate the rotation matrix
	// ----------------------------------------------	
	for(int ii = 0; ii < 3; ii++) {
		this->rotation_matrix.set_row(ii, this->k_directions[ii].get_all());
	}
		
	sout.str("");
	sout << "KPBase1D2D: the resulting rotation matrix is given by " << this->rotation_matrix;
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str()); 	
	

}

/** refresh kp matrices (if you e.g. changed some material properties in the material database 
 *//*
void KPBase1D2D::recalculate_kp_matrices() {
	for(vector<KPMatrixBase*>::iterator it = this->kp_matrices.begin(); 
	    it != this->kp_matrices.end(); it++) {
		(*it)->enforce_recalculation();    	
	}	
}*/
/*
void KPBase1D2D::copy_kp_matrices(vector<KPMatrixBase*>& copy_into) const {
	copy_into.clear();
	for(vector<KPMatrixBase*>::const_iterator it = this->kp_matrices.begin(); 
	    it != this->kp_matrices.end(); it++) {
	    KPMatrixBase* mat = this->get_matrix();
	    *mat = **it;
		copy_into.push_back(mat);	
	}	
}*/

/** limit the bandstructure calculation to a certain energy range.
 * 
 *  if we are interested only in bound states up to a certain energy,
 *  this limit can be set. the bandstructure calculation terminates
 *  when all subbands passed the limit. 
 * 
 */
void KPBase1D2D::set_upper_energy_limit(const double& upper_limit) {	 
	upper_energy_limit = upper_limit;
	upper_energy_limit_set = true;
}

/** return eigenproblem type (needed for eigensolver) (hole == negative definite, elec = positive definite) */
EigenProblemType KPBase1D2D::get_problem_type() const {
	return negative_definite;	
}

/** delete calculated results */ 
void KPBase1D2D::delete_solutions() {	
	for(unsigned int ii = 0; ii < this->bandstructures.size(); ii++) {	
		if(this->bandstructures[ii] != 0) {	
			delete this->bandstructures[ii];
			this->bandstructures[ii] = 0;
		}	
	}
	this->bandstructures.resize(0);		
}

/** return pointer to single dispersion bandstructure object */
const BandstructureDomain<cplx>& KPBase1D2D::get_bandstructure(int idx) const {
	TDKP_ASSERT(idx < (signed)this->bandstructures.size(), "bandstructure index out of range " << idx << " max: " << this->bandstructures.size());
	if(idx < 0) {
		idx = this->bandstructures.size() - 1;
	}  
	return *this->bandstructures[idx]; 
}

} // end of namespace
