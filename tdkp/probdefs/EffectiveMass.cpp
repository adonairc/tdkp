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

#include "tdkp/probdefs/EffectiveMass.h"
#include "tdkp/main/FEMSolverGEVP.h"

namespace tdkp
{

/** sparsity pattern of the interaction matrix
 * 
 * in the kp theory, we solve interacting schroedinger equations. this leads
 * to a local kp matrix at every node. therefore, if we have an element with
 * say 4 nodes and we use e.g. 6x6 kp method, the element matrix will have 
 * a size of 24 * 24 elements. there is at each node - combination a 6x6
 * block matrix. 
 * 
 * as the interaction matrices are not really full, we only assemble the 
 * nonzero components then. 
 * 
 * the sparsity pattern array is therefore an array of integer pairs. 
 * the two integers denote the row/col index which is nonzero in the interaction 
 * matrix
 * 
 * example: for a 3x3 band case where only the 1 and 3 band would interact, the 
 * kp matrix would have nonzero entries one the diagonal and on (0,2) and (2,0)
 * therefore the sparsity pattern would be 
 * {0,0, 0,2, 1,1, 2,0, 2,2}
 * 
 * 
 * buf for eff mass calculation, there is only one band, so one equation and
 * therfore only nonzero at 0,0 and one equation per node
 */ 
const int EffectiveMass::sparsity_pattern[2] = {0, 0};

EffectiveMass::EffectiveMass(const Geometry& geometry_, MaterialDatabase& material_database_) 
: SchroedingerProblem<EigenProblem3D<complex<double>, double> >(geometry_, material_database_),
  solution_type(electrons),
  energy_guess_set(false),
  energy_guess(0.0)
{
	this->num_equations_per_node = 1;
	this->sparsity_copy			   = new int[2];
	for(int ii = 0; ii < 2; ii++) {
		this->sparsity_copy[ii] = sparsity_pattern[ii];
	}
	for(int ii = 0; ii < 3; ii++) {
		hydro_strain_potential_field_name_cb[ii] = string("strain_potential_ac");
		hydro_strain_potential_field_name_vb[ii] = string("strain_potential_av"); 
	}
}
			
EffectiveMass::~EffectiveMass() {
	delete[] this->sparsity_copy; this->sparsity_copy = 0;
	this->delete_solutions();
}

/** set hydrostatic strain potential field name 
 * 
 * wurtzite has different strain potentials for a-axis and c-axis strain. therefore
 * you may set that in the calculation here accordingly
 */
void EffectiveMass::set_hydro_strain_potential_field_name_cb(unsigned int axis, const char* field_name) {
	TDKP_ASSERT(axis < 3, "");
	hydro_strain_potential_field_name_cb[axis] = string(field_name);
	this->ready = false;	
}

void EffectiveMass::set_hydro_strain_potential_field_name_vb(unsigned int axis, const char* field_name) {
	TDKP_ASSERT(axis < 3, "");
	hydro_strain_potential_field_name_vb[axis] = string(field_name);
	this->ready = false;	
}

const int* EffectiveMass::get_node_sparsity_pattern(int &num) const {
	num = 1;
	return this->sparsity_copy;
}

void EffectiveMass::prepare() {
				
	if(!this->ready) {				
		const char*   effmass_field;
		const char*   bandedge_field;
		const string* strain_potential_field;
		double        soltype_sign;
			
		// -----------------------------------------------
		// some user information
		// -----------------------------------------------
		ostringstream sout;
		sout << "EffectiveMass: solving for " << this->solution_type;
		if(this->strain_field_set() && this->potential_energy_field_set()) {
			sout << ", including strain and a potential\n";	
		} else if(this->strain_field_set()) {
			sout << ", including strain";	
		} else if(this->potential_energy_field_set()) {
			sout << ", including a potential";	
		}
		if(Configuration::get_instance()->get("kpmatrix_ignore_second_order_strain_dependence") == 1.0) {
			sout << ", but IGNORING 2nd order strain dependence";	
		}
		sout << " using following params:";		
			
		if(this->solution_type == electrons) {		
			effmass_field          = "electron_effective_mass";
			bandedge_field         = "conduction_band_edge";
			strain_potential_field = hydro_strain_potential_field_name_cb; 
			soltype_sign   = 1.0;		
		} else {
			effmass_field          = "hole_effective_mass";
			bandedge_field         = "valence_band_edge";
			strain_potential_field = hydro_strain_potential_field_name_vb;
			soltype_sign           = -1.0;
		}
								
		this->properties.clear();
		RegionProperties prop;		
		for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
			sout << "\n" << this->material_db.get_material_name(ii) << ": meff = "
			     << this->material_db.get_material(ii)->get(effmass_field);		     
			// change sign for solution				
			prop.inv_effective_mass(0,0) = prop.inv_effective_mass(1,1) = prop.inv_effective_mass(2,2) 
			  = soltype_sign /  this->material_db.get_material(ii)->get(effmass_field);
			prop.particle_potential_energy    = this->material_db.get_material(ii)->get(bandedge_field);
			// get strain potentials
			for(unsigned int jj = 0; jj < 3; jj++) {
				// ---------------------------------
				// catch incorrect material key
				// ---------------------------------
				if(this->strain_field_set()) {
					if(!this->material_db.get_material(ii)->valid_key(strain_potential_field[jj].c_str())) {
						TDKP_GENERAL_EXCEPTION("invalid material property " << strain_potential_field[jj].c_str() << ". probably you are trying to use effmass for wurtzite, without setting the right field names!");	
					}
					if(this->material_db.get_material(ii)->is_set(strain_potential_field[jj].c_str())) {						
						prop.hydrostatic_strain_potential[jj] = this->material_db.get_material(ii)->get(strain_potential_field[jj].c_str());
						sout << ", a" << (this->solution_type == electrons ? "c" :"v") << jj << " = " << prop.hydrostatic_strain_potential[jj]; 
					} else {			
						ostringstream sout;
						sout << "hydrostatic strain potential " << jj << " in material " << this->material_db.get_material_name(ii) << " is missing. strain dependence will not be available.";
						Logger::get_instance()->emit(LOG_WARN, sout.str());			
						prop.hydrostatic_strain_potential[jj] = 0.0;
					}
				} 
			}
			this->properties.push_back(prop);
		}
		TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());
	}					
	this->ready = true;

}

/** return eigenproblem type (needed for eigensolver) (hole == negative definite, elec = positive definite) */
EigenProblemType EffectiveMass::get_problem_type() const {
	if(this->get_solution_type() == holes) {
		return negative_definite;	
	} else {
		return positive_definite;	
	}
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
void EffectiveMass::calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* internal_idx, int &n) const {
	
	TDKP_ASSERT(this->properties.size() <= (unsigned)this->material_db.get_num_materials(),"this->properties.size() <= (unsigned)this->material_db.get_num_materials()");
	
	double particle_potential_energy[Element::max_num_nodes], tmp, tmp_massmat;
	  	
	n = 0;
							
	// -----------------------------------------------
	// get properties and add energy shift
	// -----------------------------------------------
	RMatrix<double> inv_effective_mass = this->properties[elem->get_region().get_material().get_id()].inv_effective_mass;
	const double& tp    = this->properties[elem->get_region().get_material().get_id()].particle_potential_energy;
	double energy_shift = get_energy_shift();
	for(unsigned int nn = 0; nn < elem->get_num_nodes(); nn++) {	
		particle_potential_energy[nn] = tp + energy_shift;
	}
	
	// -----------------------------------------------
	// check for additional potential
	// -----------------------------------------------
	if(this->potential_energy_field_set()) {
		for(unsigned int nn = 0; nn < elem->get_num_nodes(); nn++) {				
			particle_potential_energy[nn] += this->get_potential_energy_field().get_node_value(elem->get_node(nn).get_index_global(), 0);
		}			
	}
	
	// -----------------------------------------------
	// check for hydrostatic strain potential
	// -----------------------------------------------
	if(this->strain_field_set()) {
		const double* strain_potential = this->properties[elem->get_region().get_material().get_id()].hydrostatic_strain_potential;		
		for(unsigned int nn = 0; nn < elem->get_num_nodes(); nn++) {							
			const StrainTensor& strain_tensor = this->get_strain_field().get_nodal_strain(elem->get_index_global(), nn);
			// strain shift
			for(unsigned int ii = 0; ii < 3; ii++) {
				particle_potential_energy[nn] += strain_potential[ii] * strain_tensor.get(ii,ii);  
			}
		}		
		// -----------------------------------------------------
		// calculate (1-e) (use element average strain)
		// -----------------------------------------------------
		if(Configuration::get_instance()->get("kpmatrix_ignore_second_order_strain_dependence") == 0.0) {
			const StrainTensor& strain_tensor = this->get_strain_field().get(elem->get_index_global());
			RMatrix<double> strain_transformation(3,3);
			for(short ee = 0; ee < 3; ee++) {
				strain_transformation(ee,ee) = 1.0e0;
				for(short ff = 0; ff < 3; ff++) {
					strain_transformation(ee,ff) -= strain_tensor.get(ee,ff);	
				}	
			}		
			// -----------------------------------------------------
			// calculate transformed effective mass tensor
			// given by (1-e)^T M (1-e)
			// -----------------------------------------------------
			inv_effective_mass = strain_transformation.get_transpose() * inv_effective_mass * strain_transformation;
		} 	
	}
			
	// build node_internal_indices (nonzero nodes, which we really calculate)
	int ii, jj;
	for(ii = 0; ii < (signed)elem->get_num_nodes(); ii++) {
		if(elem->get_node(ii).get_index_internal() != -1) {
			internal_idx[n++] = ii;
		}	
	}
	
	// its only one equation, therefore skipping the local interaction matrices stuff	
	for(ii = 0; ii < n; ii++) {
		for(jj = 0; jj < n; jj++) {
			// build stiffness matrix
			tmp = 0.0;
			for(unsigned int kk = 0; kk < geometry.get_dimension(); kk++) {
				for(unsigned int ll = 0; ll < geometry.get_dimension(); ll++) {
					tmp += inv_effective_mass(kk,ll) * elem->get_element_integral_2nd_order(kk, internal_idx[ii], ll, internal_idx[jj]);
				}								
			}			
			tmp_massmat = elem->get_element_integral_0th_order(internal_idx[ii],internal_idx[jj]);
			lhs[ii * n + jj] = (0.5 * constants::hbar_square / (constants::m0)) * tmp; 			                 
			rhs[ii * n + jj] = tmp_massmat;
			for(unsigned int mm = 0; mm < elem->get_num_nodes(); mm++) {
				tmp_massmat = elem->get_element_integral_0th_order_nodal_data(mm, internal_idx[ii],internal_idx[jj]);
				lhs[ii * n + jj] += tmp_massmat * particle_potential_energy[mm];
			}						                
		}	
	}	
}

void EffectiveMass::solve(int num_solutions) {

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");	
	}
	
	delete_solutions();
	
	FEMSolverGEVP<complex<double>, double, double> solver(geometry, *this);

	solver.create_matrix_structures();
	solver.assemble_system();
	
	if(this->solution_type == electrons) {
		solver.set_ordering(ascending);	
	} else {
		solver.set_ordering(descending);			
	}
	
	solver.solve_system(num_solutions);	
		
	this->update_bandstructure_container();
	
}

void EffectiveMass::add_solution(cplx solution_value, const cplx* solution_vector, int length) {
	// ----------------------------
	// apply back shift to energies
	// ----------------------------
	solution_value -= this->get_energy_shift();
	SchroedingerProblem<EigenProblem3D<complex<double>, double> >::add_solution(solution_value, solution_vector, length);	
}

/** set energy guess
 * 
 * @param energy set energy where you expect the solutions
 */
void EffectiveMass::set_energy_guess(double energy) {
	this->energy_guess = energy;
	this->energy_guess_set = true;	
}

/** unset energy guess */
void EffectiveMass::drop_energy_guess() {
	this->energy_guess_set = false;	
}

double EffectiveMass::get_energy_shift() const {
	if(this->energy_guess_set) {
		return - this->energy_guess;	
	} else {
		return 0.0;	
	}
}

StdElementData<double>* EffectiveMass::get_bandedges() throw(Exception*) {
	

	// ------------------------------------------
	// init if necessary
	// ------------------------------------------
	if(!this->ready) {
		TDKP_GENERAL_EXCEPTION("object is not prepare! call prepare before calculating bandedges");
	}
 	 		
	// ------------------------------------------
    // init new element data object
    // ------------------------------------------
    StdElementData<double>* bandedge = new StdElementData<double>(1, this->geometry.get_num_elements());
    bandedge->set_identifier(0, "BandEdge0");

	// ------------------------------------------
	// loop over all elements 
	// ------------------------------------------
	for(unsigned int ii = 0; ii < this->geometry.get_num_elements(); ii++) {
		const Element& elem = this->geometry.get_element(ii);
		
		if(elem.enabled()) {
		
			// -----------------------------------------------
			// get properties
			// -----------------------------------------------
			double particle_potential_energy = this->properties[elem.get_region().get_material().get_id()].particle_potential_energy;
	
			// -----------------------------------------------
			// check for additional potential
			// -----------------------------------------------
			if(this->potential_energy_field_set()) {
				double tmp_potential = 0.0;
				for(unsigned mm = 0; mm < elem.get_num_nodes(); mm++) {
					tmp_potential += get_potential_energy_field().get_node_value(elem.get_node(mm).get_index_global(), 0);
				} 	
				particle_potential_energy += (tmp_potential / static_cast<double>(elem.get_num_nodes()));			
			}
			// -----------------------------------------------
			// check for hydrostatic strain potential
			// -----------------------------------------------
			if(this->strain_field_set()) {				
				const double* strain_potential = this->properties[elem.get_region().get_material().get_id()].hydrostatic_strain_potential;		
				const StrainTensor& strain_tensor = this->get_strain_field().get(elem.get_index_global());
				// strain shift
				for(unsigned int jj = 0; jj < 3; jj++) {
					particle_potential_energy += strain_potential[jj] * strain_tensor.get(jj,jj);  
				}				
			}
			bandedge->set_element_value(ii, 0, particle_potential_energy);
		} 			
	}			
	return bandedge;
}


BandstructureDomain<complex<double> >* 
create_effmass_dispersion(const double& transverse_effmass,
						  const DomainMaster& domain,
                          const Bandstructure<complex<double> >& bands) {

	TDKP_ASSERT(transverse_effmass != 0.0, "sorry, but the effective mass you passed is exactly 0. i have therefore to assume that this is an error.");
	TDKP_ASSERT(bands.get_number_of_k_values() == 1, "bands.get_number_of_k_values() == 1 (does only work for eff mass)");
	// ------------------------------------------------
	// create new bandstructure object
	// ------------------------------------------------ 
	BandstructureDomain<complex<double> >* ret = new BandstructureDomain<complex<double> >(
		bands.get_basis_size(),
		bands.get_number_of_bands(),
		bands.get_solution_length(),
		domain
	); 
	// ------------------------------------------------
	// create bands (but no wavefunction!)
	// ------------------------------------------------
	for(int ii = 0; ii < bands.get_number_of_bands(); ii++) {
		for(unsigned int kk = 0; kk < domain.get_number_of_points(); kk++) {
			double kval = domain.get_point(kk).get_coord_abs();
			complex<double> tmp = bands.get_energy(0,ii) 
			                    + (constants::hbar * constants::hbar) * kval * kval 
			                    / (2.0 * constants::m0 * transverse_effmass);
			ret->get_energy(kk,ii) = tmp;
		}	
	}
	return ret;	
}

void EffectiveMass::get_minmax_edges(double& cb_min, double& cb_max, double& vb_min, double& vb_max) const {
	TDKP_GENERAL_EXCEPTION("implement me");	
}


} // end namespace
