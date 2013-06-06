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

#ifndef LINEARPROBLEM_H_
#define LINEARPROBLEM_H_

#include "tdkp/probdefs/ProblemDefinition.h"

namespace tdkp {

/** problem class for pde's resulting in Ax = b equations */
template<class T> 
class LinearProblem : public tdkp::ProblemDefinition<T, T, T> {
public:
	LinearProblem(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~LinearProblem();
	
	virtual void solve(int num_solutions);
	virtual StdNodeData<T>* get_solution() const throw(Exception* );
	virtual void display_solution_info() const;
	
	
	// ----------------------------------------------------
	// intercommunication with fem solver
	// ----------------------------------------------------
	virtual const int*  get_node_sparsity_pattern(int &num) const = 0;
	virtual void 		calculate_element_matrices(const Element* elem, T* lhs, T* rhs, int* node_internal_indices, int &n) const = 0;
	virtual void add_solution(const T* solution_vector, int length);		
	void set_build_rhs_only(bool rhs_only_) { rhs_only_ = rhs_only; }
	bool get_build_rhs_only() const { return rhs_only; } 

protected:	
	vector<T> solution;
	
	virtual string get_equation_label(int idx) const throw(Exception*);
private:
	bool rhs_only;
	
};

template<class T>
LinearProblem<T>::LinearProblem(const Geometry& geometry_, MaterialDatabase& material_database_)
: ProblemDefinition<T,T,T>(geometry_, material_database_), 
  rhs_only(false)
{	
}

template<class T> 
LinearProblem<T>::~LinearProblem() {
} 

/** return StdNodeData class containing solution
 */
template<class T> 
StdNodeData<T>* LinearProblem<T>::get_solution() const throw(Exception* ) {
	
	TDKP_ASSERT(this->solution.size() > 0, "no solution available");
	TDKP_ASSERT(this->solution.size() == this->get_num_equations_per_node() * this->geometry.get_num_nonzero_nodes(), "wrong solution length");
	 	
	StdNodeData<T>* ret = new StdNodeData<T>(this->get_num_equations_per_node(), this->geometry.get_num_nodes());
	int neq = this->get_num_equations_per_node();	
	for(int nn = 0; nn < neq; nn++) {
		ret->set_identifier(nn, this->get_equation_label(nn).c_str());
	}
	for(Geometry::node_const_iterator it = this->geometry.nodes_begin(); it != this->geometry.nodes_end(); it++) {
		if((*it)->get_index_internal() != -1) {					
			for(int ii = 0; ii < this->get_num_equations_per_node(); ii++) {
				ret->set_node_value((*it)->get_index_global(), ii, this->solution[neq * (*it)->get_index_internal() + ii] );
			}
		}
    }    
	return ret;
		
}

/** internal function 
 * 
 * used by FEMSolverLE to add obtained solution to problem class
 */
template<class T>
void LinearProblem<T>::add_solution(const T* solution_vector, int length) {
	solution.resize(length);
	for(int ii = 0; ii < length; ii++) {
		this->solution[ii] = solution_vector[ii];	
	}
}

/** returns a string with a possibly meaningful label for the equation with index idx */
template<class T> 
string LinearProblem<T>::get_equation_label(int idx) const throw(Exception*) {
	if(idx < 0 || idx >= this->get_num_equations_per_node()) {
		TDKP_GENERAL_EXCEPTION("label index does not exist");	
	}
	ostringstream sout;
	sout << "solution_" << idx;
	return sout.str();
}

/** dummy solution info */
template<class T>
void LinearProblem<T>::display_solution_info() const {
	if(this->solution.size() > 0) {
		Logger::get_instance()->emit(LOG_INFO, "solution available");	
	} else {
		Logger::get_instance()->emit(LOG_INFO, "no solution available");
	}
}

} // end ofnamespace

#endif /*LINEARPROBLEM_H_*/
