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

#ifndef POISSONEQUATION_H_
#define POISSONEQUATION_H_

#include "tdkp/probdefs/LinearProblem.h"
#include "tdkp/main/FEMSolverLE.h"

namespace tdkp {


class PoissonEquation : public tdkp::LinearProblem<double> {
public:
	PoissonEquation(Geometry& geometry_, MaterialDatabase& material_database_);
	PoissonEquation(Geometry& geometry_, MaterialDatabase& material_database_, NoSolver<double>* no_solver, const double& vacuum_permittivity);
		
	virtual ~PoissonEquation();
			
	virtual StdNodeData<double>* get_solution() const throw(Exception* );		
		
	// -----------------------------------------------------------
	// charge field settings
	// -----------------------------------------------------------
	/** set charge density defined on elements
	 */
	void set_element_charge_density(const ElementData<double>* element_charge_density_);
	/** return current element charge density
	 * 
	 * returns pointer to current element charge density or null pointer
	 * of no element charge density is set
	 */
	const ElementData<double>* get_element_charge_density() const;
	/** set charge density defined on nodes	 
	 */
	void set_node_charge_density(const NodeData<double>* node_charge_density_);
	/** return current node charge density
	 * 
	 * returns pointer to current node charge density or null pointer
	 * of no node charge density is set
	 */
	const NodeData<double>* get_node_charge_density() const;
	/** set surface charge density defined on nodes for every element boundary */
	void set_surface_charge_density(const NodeData<double>* surface_charge_density_);
	/** return current surface charge density
	 * 
	 * returns pointer to current surface charge density or null pointer
	 * of no surface charge density is set
	 */
	const NodeData<double>* get_surface_charge_density() const;	
	
	// -----------------------------------------------------------
	// non-isotropic permittivity settings (but only on diagonal
	// terms)
	// -----------------------------------------------------------
	void set_permittivity(unsigned int material_index, const double& epsilon_xx, const double& epsilon_yy, const double& epsilon_zz); 			
		
	// -----------------------------------------------------------
	// set dirichlet values (values of potential at dirichlet bnd)
	// -----------------------------------------------------------
	void set_dirichlet_bnd_values(const NodeData<double>* dirichlet_values);
	
	// -----------------------------------------------------------
	// set vacuum permittivity
	// -----------------------------------------------------------
	void set_vacuum_permittivity(const double& vacuum_permittivity_);
								
	// -----------------------------------------------------------
	// for communication between FEMSolver and problem class
	// -----------------------------------------------------------
	virtual const int* get_node_sparsity_pattern(int &num) const;	
	virtual void calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* node_internal_indices, int &n) const;   	
	virtual void prepare();
	virtual void solve(int num_solutions); 
	void display_solution_info() const {};
	virtual string get_unique_identifier() const { return string("PoissonEquation"); }
		
	const double* assemble_and_return_rhs() const;
		
protected:	
	virtual string get_equation_label(int idx) const throw(Exception*);

	class Permittivity {
	public:
		Permittivity() { diag[0] = diag[1] = diag[2] = 0.0e0; }
		Permittivity(const double& value) { diag[0] = diag[1] = diag[2] = value; }	
		double diag[3];	
	};
	vector<int> sparsity_pattern;	
	vector<Permittivity> permittivities;	
	const ElementData<double>* element_charge_density;
	const NodeData<double>*    node_charge_density;
	const NodeData<double>*    surface_charge_density;
	const NodeData<double>*    dirichlet_bnd_values;

	FEMSolverLE<double>* solver;
	double vacuum_permittivity;

};



}

#endif /*POISSONEQUATION_H_*/
