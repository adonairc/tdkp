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

#ifndef EIGENPROBLEM_H_
#define EIGENPROBLEM_H_

#include "tdkp/probdefs/ProblemDefinition.h"

namespace tdkp {


template<class T, class TMAT, class TRHS> 
class EigenProblem : public ProblemDefinition<T, TMAT, TRHS> {
public:
	EigenProblem(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~EigenProblem() {}
	
	virtual EigenProblemType get_problem_type() const = 0;
	
	// ----------------------------------------------------
	// intercommunication with fem solver
	// ----------------------------------------------------
	virtual void add_solution(T solution_value, const T* solution_vector, int length) = 0;
	
protected: 
	void normalize_solution(T* solution); 
};

/** base container class for handling potential energy and strain fields */ 
template<class T>
class SchroedingerProblem : public T {
public:
	SchroedingerProblem(const Geometry& geometry_, MaterialDatabase& material_database_);	
	void set_field(const PotentialEnergyField* field) { potential_energy_field = field; }
	void set_field(const StrainField* field) { strain_field = field; this->ready = false; }
	bool potential_energy_field_set() const { return potential_energy_field != 0; }
	bool strain_field_set() const { return strain_field != 0; }
	void remove_potential_energy_field() { potential_energy_field = 0; }
	void remove_strain_field() { strain_field = 0; }
	const StrainField& get_strain_field() const;  
	const PotentialEnergyField& get_potential_energy_field() const;
	
	virtual void set_solution_type(KPSolutionType type);
	KPSolutionType get_solution_type() const;
	 
private:
	KPSolutionType kp_solution_type;
	const StrainField* strain_field;
	const PotentialEnergyField* potential_energy_field;
	
		
}; 


template<class T, class TMAT, class TRHS>
EigenProblem<T, TMAT, TRHS>::EigenProblem(const Geometry& geometry_, MaterialDatabase& material_database_)
: ProblemDefinition<T, TMAT, TRHS>(geometry_, material_database_)
{
}

/** normalize solution ATTENTION: ONLY NONZERO VALUES 
 * 
 * although the resulting eigenvectors are M (rhs of kp equation system)
 * orthogonal and therefore almost normalized, we still have to normalize 
 * the solutions again. the reason is that the M orthogonal eigenvector
 * x has unit length for all integrals over nonzero nodes instead
 * over all nodes ...  
 */
template<class T, class TMAT, class TRHS>	
void EigenProblem<T, TMAT, TRHS>::normalize_solution(T* solution) {	
	typename Geometry::element_const_iterator it;
	double* values = new double[Geometry::max_nodes_per_element];
	double result  = 0.0;	
	int num_eq     = this->get_num_equations_per_node();
	int length     = num_eq * this->geometry.get_num_nonzero_nodes();	
	double tmp;
	// integrate norm of vector
	// for all elements
	for(it = this->geometry.elements_begin(); it != this->geometry.elements_end(); it++) {
		// only for enabled elements (non-contacts)
		if((*it)->enabled()) {
			// for all bands in the solution
			for(int nn = 0; nn < num_eq; nn++) {		
				// for all nodes in 
				for(unsigned int ii = 0; ii < (*it)->get_num_nodes(); ii++) {			
					// only non dirichlet nodes
					if((*it)->get_node(ii).get_index_internal() != -1) {
						tmp = tdkp_math::abs(solution[(*it)->get_node(ii).get_index_internal() * num_eq + nn]);
						values[ii] = tmp * tmp;
					} else {
						values[ii] = 0.0;	
					}
				}
				result += (*it)->integrate_solution(values);			
			}
		}				
	}
	TDKP_ASSERT(result > 0, "result > 0");
	result = sqrt(result);	
	for(int ii = 0; ii < length; ii++) {
		solution[ii] /= result;	
	}
	delete[] values;
}

template<class T>
SchroedingerProblem<T>::SchroedingerProblem(const Geometry& geometry_, MaterialDatabase& material_database_)
: T(geometry_, material_database_),
  kp_solution_type(holes),
  strain_field(0),
  potential_energy_field(0)
{	
}

/** return strain field */
template<class T>
const StrainField& SchroedingerProblem<T>::get_strain_field() const {
	TDKP_ASSERT(strain_field_set(), "strain_field_set()");	
	return *strain_field;
}

/** return potential energy field */
template<class T>
const PotentialEnergyField& SchroedingerProblem<T>::get_potential_energy_field() const {
	TDKP_ASSERT(potential_energy_field_set(), "potential_energy_field_set()");	
	return *potential_energy_field;	
}
	 
/** return type of solution (electrons / holes ... for conduction and valence band) */
template<class T>
KPSolutionType SchroedingerProblem<T>::get_solution_type() const {
	return this->kp_solution_type;	
}

/** set solution type */
template<class T>
void SchroedingerProblem<T>::set_solution_type(KPSolutionType type) {
	this->kp_solution_type = type;
}

} // end of namespace

#endif /*EIGENPROBLEM_H_*/
