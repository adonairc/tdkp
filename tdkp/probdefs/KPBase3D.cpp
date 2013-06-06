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

#include "tdkp/probdefs/KPBase3D.h"
#include "tdkp/main/FEMSolverGEVP.h"

namespace tdkp {

KPBase3D::KPBase3D(const Geometry& geometry_, MaterialDatabase& material_database_)
: KPBase3DParent(geometry_, material_database_),
  first_order_terms_exist(false),
  energy_guess_set(false),
  energy_guess(0.0)
{
}

KPBase3D::~KPBase3D() {
	this->delete_3D();
}

void KPBase3D::init_3D() {	
	KPMatrixBase* tmp_mat          = this->get_matrix(); 
	this->num_equations_per_node = tmp_mat->get_number_of_bands();			
	const int* sparsity_pattern    = tmp_mat->get_sparsity_pattern(this->sparsity_num);
	this->sparsity_copy	           = new int[this->sparsity_num * 2];
	for(int ii = 0; ii < (this->sparsity_num * 2); ii++) {
		this->sparsity_copy[ii] = sparsity_pattern[ii];
	}
	delete tmp_mat;	
}

void KPBase3D::delete_3D() {	
	this->delete_solutions();			
}	

/** prepare the problem class for calculation */
void KPBase3D::prepare() {
						
	if(this->geometry.get_dimension() != 3) {
		TDKP_GENERAL_EXCEPTION("this problem class is only suited for 3D problems");	
	}	
	
	if(this->kp_matrices.size() > 0) {
		for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
			delete this->kp_matrices[ii]; this->kp_matrices[ii] = 0;
		}		
		this->kp_matrices.clear();
	}
	
	double eshift;					
	// the eigensolvers search eigenvalues at 0 (e.g. arpack uses shift-innode method to search
	// for the eigenvalue auf A^(-1), then eigenvalues at 0 are the LM eigenvalues and well 
	// separated. this leads to good convergence ...
	eshift = this->get_energy_shift();		
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {	
		KPMatrixBase* mat = this->get_matrix();
		mat->set_material(this->material_db.get_material(ii));
		mat->set_material_name(this->material_db.get_material_name(ii));
		mat->set_rotation(this->rotation_matrix);
		mat->set_energy_shift(eshift);
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
			
	this->ready = true;
	
}

void KPBase3D::calculate_element_matrices(const Element* elem, complex<double>* lhs, double *rhs, int* internal_idx, int &n) const {

	TDKP_ASSERT(this->kp_matrices.size() <= (unsigned)this->material_db.get_num_materials(),"this->properties.size() <= (unsigned)this->material_db.get_num_materials()");
	
	double lstiff[3][3][Element::max_num_nodes][Element::max_num_nodes]; /* stiff matrix */
	double lfirst[3][Element::max_num_nodes][Element::max_num_nodes];    /* first order matrix */	
	double lmass[Element::max_num_nodes][Element::max_num_nodes]; 	   /* mass matrix */
	double lmass_nodal[Element::max_num_nodes][Element::max_num_nodes][Element::max_num_nodes];  /* nodal mass matrix */	
	int    nnode;              /* number of nodes */
	int    neq;                /* number of kp equations */
	int    lsize; 			   /* size of lhs, rhs */
	int    nsparse;            /* number of sparse interaction matrix elements */
	                            	
	const KPMatrixBase* mat = this->kp_matrices[elem->get_region().get_material().get_id()];				
	n     = 0;
	nnode = (signed)elem->get_num_nodes();
	neq   = this->get_num_equations_per_node();
	lsize = nnode * nnode * neq * neq;
	
	this->get_node_sparsity_pattern(nsparse);
	
	TDKP_ASSERT(mat->ready(), "kp matrix not ready");
	TDKP_ASSERT(nnode <= Element::max_num_nodes, "nnode <= Element::max_num_nodes (working arrays to small ...)");
	TDKP_ASSERT((signed)mat->get_second_order(0,0).size() == nsparse, "mat->second_order[0][0].size() == nsparse");
	TDKP_ASSERT((signed)mat->get_first_order(KPMatrixBase::op_left, 0).size() == nsparse, "mat->first_order[0].size() == nsparse");
	TDKP_ASSERT((signed)mat->get_zero_order().size() == nsparse, "mat->zero_order.size() == nsparse");
	
	// --------------------------------------------------
	// set element fields to kp matrix and recalculate
	// --------------------------------------------------
	// in the present way here we calculate the element matrices by integrating analytically.  
	// varying constants which change linearly over the element can be replaced by their 
	// average over the element.
	
	// in order to have this function to be const (so we can use openmp)
	// we need the kp matrices to be stored locally ...
	vector<cplx> second_order[3][3];
	vector<cplx> first_order[2][3];
	vector<cplx> zero_order[Element::max_num_nodes]; 
	
	// only update if element is considered in calculation
	if(nnode > 0) {
		
		// ------------------------------------------------
		// first, calculate stuff thats constant over the element
		// ------------------------------------------------		
		if(this->strain_field_set()) {
			const StrainTensor& strain_tensor = get_strain_field().get(elem->get_index_global());			
			for(short ii = 0; ii < 3; ii++) {
				for(short jj = 0; jj < 3; jj++) {
					mat->build_second_order(strain_tensor, ii, jj, second_order[ii][jj]);
				}	
				mat->build_first_order(KPMatrixBase::op_left,  strain_tensor, ii, first_order[KPMatrixBase::op_left][ii]);
				mat->build_first_order(KPMatrixBase::op_right, strain_tensor, ii, first_order[KPMatrixBase::op_right][ii]);
			}			
		} else {
			for(short ii = 0; ii < 3; ii++) {
				for(short jj = 0; jj < 3; jj++) {
					second_order[ii][jj] = mat->get_second_order(ii,jj);
				}	
				first_order[KPMatrixBase::op_left][ii]  = mat->get_first_order(KPMatrixBase::op_left,  ii);
				first_order[KPMatrixBase::op_right][ii] = mat->get_first_order(KPMatrixBase::op_right, ii);
			}
		}
		// ----------------------------------------------
		// second include nodal zero order data
		// ----------------------------------------------
		for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
			double potential = 0.0;
			if(this->potential_energy_field_set()) {
				potential = get_potential_energy_field().get_node_value(elem->get_node(mm).get_index_global(), 0);			
			}
			if(this->strain_field_set()) {
				const StrainTensor& strain_tensor = get_strain_field().get(elem->get_index_global());			
				mat->build_zero_order(strain_tensor, potential, zero_order[mm]);
			} else {
				zero_order[mm] = mat->get_zero_order();
				if(potential != 0.0) {
					mat->add_zero_order_energy(potential,zero_order[mm]);
				}			
			}
		}			
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
			for(int aa = 0; aa < 3; aa++) {
				for(int bb = 0; bb < 3; bb++) {
					lstiff[aa][bb][ii][jj] = elem->get_element_integral_2nd_order(aa, internal_idx[ii], bb, internal_idx[jj]);								
				}	
			}
			lmass[ii][jj] = elem->get_element_integral_0th_order(internal_idx[ii],internal_idx[jj]);
			for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
				lmass_nodal[mm][ii][jj] = elem->get_element_integral_0th_order_nodal_data(mm, internal_idx[ii],internal_idx[jj]);				
			}
		}	
	}	
	if(this->first_order_terms_exist) {
		for(int ii = 0; ii < n; ii++) {
			for(int jj = 0; jj < n; jj++) {
				for(int aa = 0; aa < 3; aa++) {
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
				for(int aa = 0; aa < 3; aa++) {
					for(int bb = 0; bb < 3; bb++) {
						lhs[offset + ss] += second_order[aa][bb][ss] *  lstiff[aa][bb][ii][jj];	
					}
				}
				// assemble "mass" contributions
				for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
					lhs[offset + ss] += zero_order[mm][ss] * lmass_nodal[mm][ii][jj];
				}				
				// assemble first order contributions
				if(this->first_order_terms_exist) {					
					for(int aa = 0; aa < 3; aa++) {
						lhs[offset + ss] += first_order[KPMatrixBase::op_left][aa][ss]  * lfirst[aa][jj][ii];
						// the minus sign is due to partial integration
						lhs[offset + ss] += - first_order[KPMatrixBase::op_right][aa][ss] * lfirst[aa][ii][jj]; // switched ii jj
					}
				}
			}
			rhs[ii * n + jj] = lmass[ii][jj];
		}
	}										
}

void KPBase3D::solve(int num_solutions) {

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");	
	}

	this->delete_solutions(); // delete any old solutions
	this->prepare();
	
	FEMSolverGEVP<complex<double>, complex<double>, double> solver(geometry, *this);
	if(this->get_solution_type() == electrons) {
		solver.set_ordering(ascending);
	} else {
		solver.set_ordering(descending);	
	}	
	solver.create_matrix_structures();
	solver.assemble_system();
	solver.solve_system(num_solutions);	
		
	this->update_bandstructure_container();
			
}

void KPBase3D::add_solution(cplx solution_value, const cplx* solution_vector, int length) {
	// ----------------------------
	// apply back shift to energies
	// ----------------------------
	solution_value -= this->get_energy_shift();
	KPBase<SchroedingerProblem<EigenProblem3D<complex<double>, complex<double> > > >::add_solution(solution_value, solution_vector, length);	
}

const int* KPBase3D::get_node_sparsity_pattern(int &num) const {
	num = this->sparsity_num;
	return this->sparsity_copy;		
}


/** set energy guess
 * 
 * @param energy set energy where you expect the solutions
 */
void KPBase3D::set_energy_guess(double energy) {
	this->energy_guess = energy;
	this->energy_guess_set = true;	
}

/** unset energy guess */
void KPBase3D::drop_energy_guess() {
	this->energy_guess_set = false;	
}

/** determine the energy shift 
 * 
 * depending on the type of solution, the energy of the system needs to be shifted
 * in order to have the solutions of interest at energy 0.0
 * depending on solution_type, we return here either the max. valence band edge or minimum 
 * conduction band edge
 * 
 * @return depending on solution_type, we return here either the max. valence band edge or minimum conduction band edge (in [eV])
 */
double KPBase3D::get_energy_shift() const {
	
	const char* key;
	double asc;
	double shift;
	
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
	if(this->energy_guess_set) {
		ostringstream sout;
		sout << "will expect states to be around " << this->energy_guess << " [eV]";
		Logger::get_instance()->emit(LOG_INFO, sout.str());
		return - this->energy_guess;
	} else {
		shift = this->material_db.get_material(0)->get(key);
		for(int ii = 1; ii < this->material_db.get_num_materials(); ii++) {
			if(shift * asc > this->material_db.get_material(ii)->get(key) * asc) {
				shift = this->material_db.get_material(ii)->get(key);	
			}
		}
		return (- shift);	
	}
}

/** set crystal axes
 *
 * set into which crystal axis correspond to the system axis  
 * @param k_x	crystal direction of x axis
 * @param k_y   crystal direction of y axis
 * @param k_z   crystal direction of z axis
 */
void KPBase3D::set_axes(const Vector3D& k_x, const Vector3D& k_y, const Vector3D& k_z) {
	
	double        dot_product;
	double        determinante;
	ostringstream sout;
	const string  named_directions[] = {string("k_x"), string("k_y"), string("k_z")};
	const double  tolerance = 1.0e-15;

	sout << "using the following axes of quantization:\n"
	     << "quantized 1 (x) -> " << k_x << "\n"
	     << "quantized 2 (y) -> " << k_y << "\n"
	     << "quantized 3 (z) -> " << k_z;	      
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str()); 	   
	sout.str("");  
	
	// -----------------------------------------------
	// store vectors and normalize them
	// -----------------------------------------------
	this->k_directions[D_DX] = k_x;
	this->k_directions[D_DY] = k_y;
	this->k_directions[D_DZ] = k_z;
	for(int ii = 0; ii < 3; ii++) {
		this->k_directions[ii].normalize();	
	}
	
	// ------------------------------------------------
	// test vectors (they must be orthogonal!)
	// -----------------------------------------------	
	for(int ii = 0; ii < 2; ii++) {
		dot_product = Vector3D::dot_product(this->k_directions[ii], this->k_directions[ii + 1]);
		if(fabs(dot_product) > tolerance) {
			sout << "your directions " << named_directions[ii] << " and " 
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
		sout << "the matrix of the resulting rotation given by your directions\n"
		     << "has a determinante of " << determinante << " which is definitely wrong.\n"
		     << "the matrix itself looks that way: " << this->rotation_matrix;
		TDKP_GENERAL_EXCEPTION(sout.str());
	} else if(fabs(determinante - 1.0) > tolerance) {
		sout << "the rotation matrix given by your directions has a determinante of\n"
		     << determinante << " which deviates from 1.0 more than the tolerance\n" 
		     << tolerance << "\n" 
		     << "the matrix itself looks that way: " << this->rotation_matrix;
		TDKP_GENERAL_EXCEPTION(sout.str());		      	
	}	
	sout.str("");
	sout << "the resulting rotation matrix is given by " << this->rotation_matrix;
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str()); 	

}


/** return eigenproblem type (needed for eigensolver) (hole == negative definite, elec = positive definite) */
EigenProblemType KPBase3D::get_problem_type() const {
	return negative_definite;	
}

} // end namespace
