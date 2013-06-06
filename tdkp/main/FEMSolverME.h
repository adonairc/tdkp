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

#ifndef FEMSOLVERME_H_
#define FEMSOLVERME_H_

#include "tdkp/main/FEMSolverLE.h"

namespace tdkp {

/** fem assembler for matrix element and coulomb matrix element solving */
template<class T>
class FEMSolverME : public FEMSolverLE<T> {
public:	
	FEMSolverME(const Geometry& geometry_, LinearProblem<T>& problem_, SparseMatrixProperties matrix_type);
	virtual ~FEMSolverME();
	template<class V> 
	void multiply_with_lhs(const vector<V>& in, vector<V>& out) const;
	unsigned int get_matrix_size() const { return nosolver->get_matrix().get_size(); }
	const CSRMatrix<T>& get_csr_matrix() const { return nosolver->get_csr_matrix(); }	
private:
	NoSolver<T>* nosolver;
	
};


/** constructor for matrix element matrix assembly class */
template<class T>
FEMSolverME<T>::FEMSolverME(const Geometry& geometry_, LinearProblem<T>& problem_, SparseMatrixProperties matrix_type)
: FEMSolverLE<T>(geometry_, problem_),
  nosolver(0)
{
	// ------------------------------------------
	// initialize no solver and get matrix
	// ------------------------------------------
	int matrix_size = this->problem.get_num_equations_per_node()
	                * this->geometry.get_num_nonzero_nodes();
	nosolver     = new NoSolver<T>(matrix_size, matrix_type);
	this->solver = nosolver;
	this->stiff  = &nosolver->get_matrix();
	this->set_be_quiet(true);
	                
}

/** destructor deletes nosolver object */
template<class T>
FEMSolverME<T>::~FEMSolverME() {	
	nosolver = 0;	
}

template<class T> template<class V> 
void FEMSolverME<T>::multiply_with_lhs(const vector<V>& in, vector<V>& out) const {
	TDKP_ASSERT(in.size() == this->stiff->get_size(), "input vector size is NOT correct");
	TDKP_ASSERT(out.size() == this->stiff->get_size(), "output vector size is NOT correct");
	this->nosolver->get_csr_matrix().mult_vec(&in[0], &out[0]);	
}


} // end of namespace tdkp

#endif /*FEMSOLVERLE_H_*/
