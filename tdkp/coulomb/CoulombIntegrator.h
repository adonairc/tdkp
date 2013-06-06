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

#ifndef COULOMBINTEGRATOR_H_
#define COULOMBINTEGRATOR_H_

#include "tdkp/common/all.h"
#include "tdkp/solvers/LinearSolver.h"
#include "tdkp/main/FEMSolverME.h"
#include "tdkp/main/CSRMatrix.h"
#include "tdkp/coulomb/CoulombFunction.h"

namespace tdkp {

/** we evaluate the integrator at the mid points of elements ... */
class ElementMiddlePoints {
public:
	ElementMiddlePoints(const Geometry& geometry);
	const double* get_coords(unsigned int elem_idx) const { return &coordinates[elem_idx * dimension]; }
	const vector<double>& get_coordinates() const { return coordinates; } 
private:
	vector<double> coordinates;
	unsigned int   dimension;
};


/** object to evaluate dot products of two form factors for special matrices 
 */ 
template<class T>
class CoulombIntegrator : public LinearProblem<T> {
public:
	CoulombIntegrator(
		const Geometry& geometry_, 
		MaterialDatabase& material_database_,
		const ElementMiddlePoints& element_middle_points_
	);
	virtual ~CoulombIntegrator() {}
	void build_matrix(const CoulombFunction& function_object);
	void build_matrix(const vector<T>& values_at_midpoints);
	
	template<class V>
	void multiply_with_lhs(const vector<V>& in, vector<V>& out) const;
	unsigned int get_matrix_size() const;
			
	// ----------------------------------------------
	// intercommunication with FEMSolverLE
	// ----------------------------------------------
	virtual void calculate_element_matrices(const Element* elem, T* lhs, T *rhs, int* node_internal_indices, int &n) const;
	virtual const int* get_node_sparsity_pattern(int& n) const { n = 1; return sparsity_pattern; }
	virtual string get_unique_identifier() const { return string("CoulombIntegrator"); }  	
						
private:
	const ElementMiddlePoints& element_middle_points;
	FEMSolverME<T> assembler;
	static const int sparsity_pattern[2];
	vector<T> element_values;

};

template<>
const int CoulombIntegrator<cplx>::sparsity_pattern[2];
template<>
const int CoulombIntegrator<double>::sparsity_pattern[2];
 
template<class T>
CoulombIntegrator<T>::CoulombIntegrator(
	const Geometry& geometry_, MaterialDatabase& material_database_, 
	const ElementMiddlePoints& element_middle_points_)	
: LinearProblem<T>(geometry_, material_database_),
  element_middle_points(element_middle_points_),
  assembler(
  	geometry_, *this,
  	symmetric_matrix
  ),
  element_values(geometry_.get_num_elements())
{
	this->assembler.create_matrix_structures();
}

/** multiply vector with assembled matrix */
template<class T> template<class V>
void CoulombIntegrator<T>::multiply_with_lhs(const vector<V>& in, vector<V>& out) const {
	assembler.multiply_with_lhs(in,out);	
}

/** return size of matrix */
template<class T>
unsigned int CoulombIntegrator<T>::get_matrix_size() const {
	return assembler.get_matrix_size();
}

/** build new matrix for a function given by its values at the midpoints of all elements*/
template<class T>
void CoulombIntegrator<T>::build_matrix(const vector<T>& node_values) {

	TDKP_ASSERT(node_values.size() == this->geometry.get_num_nonzero_nodes(), "node_values.size() == geometry.get_num_nonzero_nodes()");
	
	// ---------------------------------------
	// calculate element values
	// ---------------------------------------
	element_values.assign(this->geometry.get_num_elements(), 0.0);	
	unsigned int nnode = 0;
	for(unsigned int ii = 0; ii < this->geometry.get_num_elements(); ii++) {
		const Element& elem = this->geometry.get_element(ii);
		nnode = 0;
		for(unsigned int jj = 0; jj < elem.get_num_nodes(); jj++) {
			if(elem.get_node(jj).get_index_internal() != -1) {
				element_values[elem.get_index_global()] += node_values[elem.get_node(jj).get_index_internal()];
				nnode++;
			}					
		}
		element_values[elem.get_index_global()] /= static_cast<double>(nnode);						
	}
	
	
	// --------------------------------------- 
	// assemble system
	// ---------------------------------------
	assembler.reset_matrix_to_zero();
	assembler.assemble_system();	
	 	
}

/** build new matrix for new q and z' */
template<class T>
void CoulombIntegrator<T>::build_matrix(const CoulombFunction& function_object) {

	TDKP_ASSERT(function_object.get_dimension() == this->geometry.get_dimension(), "function_object.get_dimension() == geometry.get_dimension()");

	// ---------------------------------------
	// evaluate for each element
	// ---------------------------------------
	vector<double> tmp(element_values.size());
	function_object.evaluate(element_middle_points.get_coordinates(), tmp);
	for(unsigned int ii = 0; ii < tmp.size(); ii++) {
		element_values[ii] = tmp[ii];	
	}
	
	// --------------------------------------- 
	// assemble system
	// ---------------------------------------
	assembler.reset_matrix_to_zero();
	assembler.assemble_system();
	
}	

template<class T>
void CoulombIntegrator<T>::calculate_element_matrices(const Element* elem, T* lhs, 
	T *rhs, int* internal_idx, int &n) const {
		
	double lmass[Element::max_num_nodes][Element::max_num_nodes]; 	   /* mass matrix */
	int    nnode;              /* number of nodes */
	int    lsize; 			   /* size of lhs, rhs */

	n     = 0;
	nnode = (signed)elem->get_num_nodes();
	lsize = nnode * nnode;
	
	TDKP_ASSERT(nnode <= Element::max_num_nodes, "nnode <= Element::max_num_nodes (working arrays to small ...)");
	
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
			lmass[ii][jj] = elem->get_element_integral_0th_order(internal_idx[ii],internal_idx[jj]);
		}
	}

	// ------------------------------------------------------------------------------
	// assemble ....
	// ------------------------------------------------------------------------------
	// for all nonzero nodes in element
	const T& element_value = element_values[elem->get_index_global()];
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {			
			lhs[ii * n + jj] += element_value * lmass[ii][jj];			
		}
	}	
}

}

#endif /*COULOMBINTEGRATOR_H_*/
