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

#include "tdkp/probdefs/NEGFEffectiveMass.h"
#include "tdkp/main/FEMSolverGEVP.h"

namespace tdkp {

/** basic constructor for the NEGF effective mass hamiltonian */
NEGFEffectiveMass::NEGFEffectiveMass(const Geometry& geometry_, MaterialDatabase& material_database_)
: SchroedingerProblem<EigenProblem3D<complex<double>, complex<double> > >(geometry_, material_database_),
  cb_effmass(geometry_, material_database_),
  vb_effmass(geometry_, material_database_),
  energy_shift_set(false),
  energy_shift(0.0e0),
  k_transversal(0.0e0)
{
	num_equations_per_node   = 2;
	node_sparsity_pattern[0] = node_sparsity_pattern[1] = 0;
	node_sparsity_pattern[2] = node_sparsity_pattern[3] = 1;
	cb_effmass.set_solution_type(electrons);
	vb_effmass.set_solution_type(holes); 
}

NEGFEffectiveMass::~NEGFEffectiveMass()
{
}


/** solve hamiltonian (function serves soley for interface testing) */
void NEGFEffectiveMass::solve(int num_solutions) {

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");	
	}
	
	delete_solutions();	
	FEMSolverGEVP<complex<double>, complex<double>, double> solver(geometry, *this);
	solver.create_matrix_structures();
	solver.assemble_system();
	solver.set_ordering(ascending);	
	solver.solve_system(num_solutions);
		
	// ------------------------------------------------
	// apply back shift of energies 
	// ------------------------------------------------	
	double eshift = this->get_energy_shift();
	int counter = 0;
	ostringstream sout;	
	for(std::vector<cplx>::iterator it = this->solution_value.begin(); it != this->solution_value.end(); it++) {		
		(*it) -= eshift;
		sout << "added solution " << counter ++ << " with ev: " << (*it);
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
		sout.str("");
	}	
	
	this->update_bandstructure_container();
	
}

const int* NEGFEffectiveMass::get_node_sparsity_pattern(int &num) const {
	num = 2; 
	return node_sparsity_pattern;	
}   	

/** combine cb and vb effective mass element matrices into common 2x2 block matrix */
void NEGFEffectiveMass::calculate_element_matrices(
	const Element* elem, complex<double>* lhs, double *rhs, int* node_internal_indices, int &n
) const {
	
	// ---------------------------------------------
	// create temporary space 
	// ---------------------------------------------
	vector<double> cb_lhs(Element::max_num_nodes * Element::max_num_nodes);
	vector<double> vb_lhs(Element::max_num_nodes * Element::max_num_nodes);
	
	// ---------------------------------------------
	// calculate transversal energy
	// ---------------------------------------------
	double cb_transversal = constants::hbar_square * k_transversal * k_transversal
	                      / (2.0 * constants::m0 * elem->get_material().get("electron_effective_mass_transverse"));
	double vb_transversal = - constants::hbar_square * k_transversal * k_transversal
	                      / (2.0 * constants::m0 * elem->get_material().get("hole_effective_mass_transverse"));	                      
		
	// ---------------------------------------------
	// build cb and vb element matrix
	// --------------------------------------------- 
	cb_effmass.calculate_element_matrices(elem, &cb_lhs[0], rhs, node_internal_indices, n);
	vb_effmass.calculate_element_matrices(elem, &vb_lhs[0], rhs, node_internal_indices, n);
		
	// --------------------------------------------
	// combine cb and vb element matrices
	// --------------------------------------------
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {
			// cb contrib
			lhs[(n * ii + jj) * 2] = cb_lhs[(n * ii + jj)] 
			                       + cb_transversal * elem->get_element_integral_0th_order(node_internal_indices[ii], node_internal_indices[jj]);			                            
			// vb contrib
			lhs[(n * ii + jj) * 2 + 1] = vb_lhs[(n * ii + jj)]
									   + vb_transversal * elem->get_element_integral_0th_order(node_internal_indices[ii], node_internal_indices[jj]);			 	
		}
	} 		
}
   	  
void NEGFEffectiveMass::prepare() {

	// ------------------------------------------
	// set potential energy field, complain on strains
	// ------------------------------------------
	TDKP_ASSERT(!strain_field_set(), "NEGF class can not yet handle strains. implement it (set strains to decorated objects)");
	cb_effmass.set_field(&get_potential_energy_field());
	vb_effmass.set_field(&get_potential_energy_field());
	
	// ------------------------------------------
	// set energy shifts 
	// ------------------------------------------
	cb_effmass.set_energy_guess(- get_energy_shift());
	vb_effmass.set_energy_guess(- get_energy_shift());
	
	// ------------------------------------------
	// prepare objects ...
	// ------------------------------------------	 
	cb_effmass.prepare();
	vb_effmass.prepare();
	
	this->ready = true;
		
}


/** set energy shift
 */
void NEGFEffectiveMass::set_energy_shift(double energy) {
	this->energy_shift = energy;
	this->energy_shift_set = true;	
}

/** unset energy shift */
void NEGFEffectiveMass::drop_energy_shift() {
	this->energy_shift_set = false;	
}

double NEGFEffectiveMass::get_energy_shift() const {
	if(this->energy_shift_set) {
		return this->energy_shift;	
	} else {
		return 0.0;	
	}
}

}
