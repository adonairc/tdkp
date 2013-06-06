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

#include "tdkp/probdefs/PoissonEquation.h"


namespace tdkp {

PoissonEquation::PoissonEquation(Geometry& geometry_, MaterialDatabase& material_database_)
: LinearProblem<double>(geometry_, material_database_),
  sparsity_pattern(2,0),
  element_charge_density(0),
  node_charge_density(0),
  surface_charge_density(0),
  dirichlet_bnd_values(0),
  solver(0),
  vacuum_permittivity(constants::vacuum_permittivity)
{
	geometry_.prepare_boundaries();
	this->num_equations_per_node = 1;
	this->solver = new FEMSolverLE<double>(geometry_, *this);	
	this->solver->create_matrix_structures();
	this->solver->assemble_system();
}

/** protected constructor, used by derived class for poisson equation solving using fem in aqua */
PoissonEquation::PoissonEquation(
	Geometry& geometry_, 
	MaterialDatabase& material_database_, 
	NoSolver<double>* no_solver,
	const double& vacuum_permittivity_
)
: LinearProblem<double>(geometry_, material_database_),
  sparsity_pattern(2,0),
  element_charge_density(0),
  node_charge_density(0),
  surface_charge_density(0),
  dirichlet_bnd_values(0),
  solver(0),
  vacuum_permittivity(vacuum_permittivity_)
{
	geometry_.prepare_boundaries();
	this->num_equations_per_node = 1;
	this->solver = new FEMSolverLE<double>(geometry_, *this, no_solver);		
	this->solver->create_matrix_structures();
	this->solver->assemble_system();	
}

PoissonEquation::~PoissonEquation() {
	if(this->solver != 0) {
		delete this->solver; 
		this->solver = 0;	
	}
}

/** overwrite tdkp standard vacuum permittivity */
void PoissonEquation::set_vacuum_permittivity(const double& vacuum_permittivity_) {
	this->vacuum_permittivity = vacuum_permittivity_;	
}

const int* PoissonEquation::get_node_sparsity_pattern(int &num) const {
	num = 1;
	return &this->sparsity_pattern[0];
}

void PoissonEquation::prepare() {
			
	if(element_charge_density) {
		TDKP_ASSERT(element_charge_density->get_length() == (signed)geometry.get_num_elements(), "");
		TDKP_ASSERT(element_charge_density->get_num_data_per_element() == 1, "");	
	}
	if(node_charge_density) {
		TDKP_ASSERT(node_charge_density->get_length() == (signed)geometry.get_num_nodes(), "");
		TDKP_ASSERT(node_charge_density->get_num_data_per_node() == 1, "");	
	}	

	// ------------------------------------------
	// read from database if not set by user
	// ------------------------------------------
	if(permittivities.size() == 0) {
		for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {		
			permittivities.push_back(
				Permittivity(this->material_db.get_material(ii)->get("electrostatic_permittivity"))
			);									
		}
	}
	// ------------------------------------------
	// determine maximum material name length
	// ------------------------------------------
	unsigned int max_length = 0;
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
		if(max_length < this->material_db.get_material_name(ii).size()) {
			max_length = this->material_db.get_material_name(ii).size();	
		}
	}
	// ------------------------------------------
	// list and check values
	// ------------------------------------------
	
	ostringstream sout;
	sout << "PoissonEquation: i will use the following permittivities:";
	bool unset   = false;
	bool nonpos = false;
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
		sout << "\n  " << setw(max_length) << this->material_db.get_material_name(ii);
		if(ii < (signed)permittivities.size()) {
			sout << " - eps_xx: " << setw(6) << permittivities[ii].diag[0] 
			     << ", eps_yy: "  << setw(6) << permittivities[ii].diag[1]
			     << ", eps_zz: "  << setw(6) << permittivities[ii].diag[2];
			for(unsigned int jj = 0; jj < 3; jj++) {
				if(permittivities[ii].diag[jj] <= 0.0) {
					nonpos = true;
					break;	
				}
			}
			if(nonpos) {
				sout << " - INVALID!";	
			}
		} else {
			sout << " - NOT SET!";
			unset = true;	
		}		
	}		
	if(unset || nonpos) {
		TDKP_GENERAL_EXCEPTION(sout.str());	
	} else {
		TDKP_LOGMSG(LOG_INFO_DEVEL1, sout.str());
	}
	sout.str("");
	sout << "PoissonEquation: charge sources are defined on:";
	char comma = ' ';
	if(element_charge_density) {
		sout << " elements";
		comma = ','; 	
	}
	if(node_charge_density) {
		sout << comma << " nodes";
		comma = ','; 	
	}
	if(surface_charge_density) {
		sout << comma << " element boundaries";
		comma = ',';	
	}
	
	TDKP_LOGMSG(LOG_INFO_DEVEL1, sout.str()); 	   		
	this->ready = true;
		
}

/** set values of diagonal electrostatic permittivity tensor
 * 
 * you may define the values of the electrostatic permittivity tensor 
 * by yourself (e.g. for wurtzite crystals). if you call this function
 * once, you need to set all values by yourself for all materials
 */
void PoissonEquation::set_permittivity(unsigned int material_index, const double& epsilon_xx, const double& epsilon_yy, const double& epsilon_zz) {
	if(material_index >= permittivities.size()) {
		permittivities.resize(material_index + 1);	
	}
	permittivities[material_index].diag[0] = epsilon_xx;
	permittivities[material_index].diag[1] = epsilon_yy;
	permittivities[material_index].diag[2] = epsilon_zz;
	this->ready = false;
}

void PoissonEquation::set_element_charge_density(const ElementData<double>* element_charge_density_) {
	this->element_charge_density = element_charge_density_;	
	TDKP_ASSERT(element_charge_density->get_length() == (signed)geometry.get_num_elements(), "");
	TDKP_ASSERT(element_charge_density->get_num_data_per_element() == 1, "");	
}
const ElementData<double>* PoissonEquation::get_element_charge_density() const {
	return this->element_charge_density;	
}

/** set dirichlet values (index node -> value), 0 else */
void PoissonEquation::set_dirichlet_bnd_values(const NodeData<double>* dirichlet_values_) {
	this->dirichlet_bnd_values = dirichlet_values_;
}

void PoissonEquation::set_node_charge_density(const NodeData<double>* node_charge_density_) {
	this->node_charge_density = node_charge_density_;
	TDKP_ASSERT(node_charge_density->get_length() == (signed)geometry.get_num_nodes(), "");
	TDKP_ASSERT(node_charge_density->get_num_data_per_node() == 1, "");		
}
const NodeData<double>* PoissonEquation::get_node_charge_density() const {
	return this->node_charge_density;	
}

void PoissonEquation::set_surface_charge_density(const NodeData<double>* surface_charge_density_) {
	this->surface_charge_density = surface_charge_density_;
	TDKP_ASSERT(this->surface_charge_density->get_length() == (signed)geometry.get_num_boundaries(), "");
}
const NodeData<double>* PoissonEquation::get_surface_charge_density() const {
	return this->surface_charge_density;	
}

/** calculate the local element matrices (stiff + mass) and overlap 
 * 
 * stiff and mass is stored in lhs, overlap in rhs
 * lhs and rhs must be min. of length (elem->get_num_nodes()^2 * nonzeros in interaction matrix)
 * in node_internal_indices, the nonzero node indices of the ELEMENT are stored (not globals)
 * means: if there are 4 nodes ([0 ..3]) and node 2 is on the boundary (internal index == -1)
 * then node_internal_indicies will be [0 1 3].
 *
 * @param elem the element where the matrices should be calculated
 * @param lhs array where the element matrix will be stored: lhs[(ii * kpn + jj) * num_sparse + sparse_idx]
 * @param rhs array where the element overlap matrix will be stored
 * @param internal_idx the nonzero indices -> lhs[ ii * kpn + jj 
 */ 
void PoissonEquation::calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* internal_idx, int &n) const {
	
	TDKP_ASSERT(this->permittivities.size() <= (unsigned)this->material_db.get_num_materials(),"");
	
	double tmp; 
	bool   rhs_only = this->get_build_rhs_only();	
	int ii, jj;
	n = 0;
	
	// get relative permittivity
	const Permittivity& permittivity = this->permittivities[elem->get_region().get_material().get_id()];
	 				
	// build node_internal_indices (nonzero nodes, which we really calculate)
	for(ii = 0; ii < (signed)elem->get_num_nodes(); ii++) {
		if(elem->get_node(ii).get_index_internal() != -1) {			
			internal_idx[n++] = ii;	
		}
	}
				
	// its only one equation, therefore skipping the local interaction matrices stuff
	for(ii = 0; ii < n; ii++) {
		if(!rhs_only) {		
			for(jj = 0; jj < n; jj++) {
				// build stiffness matrix
				tmp = 0.0;
				for(unsigned int kk = 0; kk < this->geometry.get_dimension(); kk++) {
					tmp += permittivity.diag[kk] 
					     * elem->get_element_integral_2nd_order(kk, internal_idx[ii], kk, internal_idx[jj]);								
				}			
				lhs[ii * n + jj] = tmp * vacuum_permittivity;
			}
		}
		rhs[ii] = 0.0;
		if(element_charge_density) {			 		
			rhs[ii] = elem->get_single_integral_0th_order(ii) 
			        * this->element_charge_density->get_element_value(elem->get_index_global());
		}
		if(node_charge_density) {			
			for(jj = 0; jj < n; jj++) {
				rhs[ii] += this->node_charge_density->get_node_value(elem->get_node(internal_idx[jj]).get_index_global())
				         * elem->get_element_integral_0th_order(internal_idx[ii], internal_idx[jj]);	
			}	
		}
	}	
	
	// add surface charges if existing
	if(surface_charge_density) {
		// -------------------------------------------------------
		// surface charges are: C * delta(x) and delta is the delta
		// function on the boundary. 
		// boundary contribution must only be calculated once, but each 
		// boundary belongs to two elements. therefore we only add surface
		// charges to rhs if this element is element 0 in the boundary
		// -------------------------------------------------------
		int elemidx2bndidx[Element::max_num_nodes];
		for(unsigned int bb = 0; bb < elem->get_num_boundaries(); bb++) {
			const ElementBoundary& bnd = elem->get_boundary(bb);
			if(bnd.get_num_elements() == 2 && bnd.get_element(0).get_index_global() == elem->get_index_global()) {
				// --------------------------------------------------
				// build map between boundary nodes and element nodes
				// --------------------------------------------------
				for(unsigned int nn = 0; nn < elem->get_num_nodes(); nn++) {
					elemidx2bndidx[nn] = -1;		
					for(unsigned int mm = 0; mm < bnd.get_num_nodes(); mm++) {
						if(bnd.get_node(mm).get_index_global() == elem->get_node(nn).get_index_global()) {
							elemidx2bndidx[nn] = mm;
							break;
						}		
					}
				}
				// --------------------------------------------------
				// add charge to rhs
				// int_bnd Ni * (sum_j Nj cj)
				// where cj is the nodal charge and Nj is to interpolate
				// the nodal charge (so charge c(x) = sum_j Nj cj)
				// --------------------------------------------------
				for(ii = 0; ii < n; ii++) {					
					if(elemidx2bndidx[internal_idx[ii]] != -1) {						
						double tmp2 = 0.0;
						for(unsigned int nn = 0; nn < bnd.get_num_nodes(); nn++) {
							tmp2 += surface_charge_density->get_node_value(bnd.get_index_global(), nn)
							         * bnd.get_boundary_integral_0th_order(elemidx2bndidx[internal_idx[ii]], nn);
						}
						rhs[ii] += tmp2;
					}
				} 	
			}				 	
		}
	}
	
	// ----------------------------------------------
	// add terms resulting from transformation of
	// inhomogeneous equation to homogeneous equation
	// load -= integral u0_j * eps * grad Nj . grad Ni 
	// ----------------------------------------------
	if(dirichlet_bnd_values) {		
		for(ii = 0; ii < n; ii++) {
			tmp = 0.0;
			// ----------------------------------------
			// sum over ALL element nodes (this goes to
			// rhs and is NOT discared in matrix equation
			// during enforcement of DC)
			// ---------------------------------------- 
			for(jj = 0; jj < (signed)elem->get_num_nodes(); jj++) {
				const double& dc_val = dirichlet_bnd_values->get_node_value(elem->get_node(jj).get_index_global(), 0);
				for(unsigned int kk = 0; kk < this->geometry.get_dimension(); kk++) {
					tmp += dc_val * permittivity.diag[kk] *  elem->get_element_integral_2nd_order(kk, internal_idx[ii], kk, jj);
				}
			}
			rhs[ii] -= tmp * vacuum_permittivity;
		}
	}	
}

StdNodeData<double>* PoissonEquation::get_solution() const throw(Exception* ) {

	// let parent object create my solution
	StdNodeData<double>* ret = LinearProblem<double>::get_solution();
	
	// go and add dirichlet values
	if(dirichlet_bnd_values) {
		TDKP_ASSERT(ret->get_length() == dirichlet_bnd_values->get_length(), "");	
		for(int ii = 0; ii < ret->get_length(); ii++) {
			ret->get_node_value(ii, 0) += dirichlet_bnd_values->get_node_value(ii, 0);		
		}
	}
	return ret;
	
}

string PoissonEquation::get_equation_label(int idx) const throw(Exception*) {
	TDKP_ASSERT(idx == 0, "");
	return string("potential");	
}

void PoissonEquation::solve(int num_solutions) {

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");	
	}
	num_solutions = 1; // surpress compiler warnings on unused stuff
	
	// ---------------------------------
	// evaluate netto charge
	// ---------------------------------
	solver->assemble_rhs_only();
	ostringstream sout;
	sout << "PoissonEquation: netto charge in geometry is " << solver->sum_rhs();
	TDKP_LOGMSG(LOG_INFO_DEVEL1, sout.str());
	solver->solve_system(num_solutions);	

}

/** assemble rhs and return pointer to it */
const double* PoissonEquation::assemble_and_return_rhs() const {
	solver->assemble_rhs_only();
	return solver->get_rhs();	
}

} // end of namespace
