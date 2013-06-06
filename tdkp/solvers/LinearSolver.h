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

#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_


#include "tdkp/common/all.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/main/CSRMatrix.h"

// ------------------------------------------
// include or remove eigensolvers
// LAPACK is always included
// ------------------------------------------
//#define LINSOLV_INCLUDE_PARDISO
//#define LINSOLV_INCLUDE_UMFPACK
//#define LINSOLV_INCLUDE_ILS
//#define LINSOLV_INCLUDE_AZTECOO
//#define LINSOLV_INCLUDE_SUPERLU

namespace tdkp {

enum LinearSolverFlags {
	direct      = 1,
	iterative   = 2,
	distributed	= 4
};

// factory implementation is in EigenSolver.cpp

template<class T> 
class LinearSolver {

public:	
	LinearSolver() { }

	virtual ~LinearSolver() throw(Exception*) { };		
	virtual void prepare() throw(Exception*) = 0;
	virtual void release() throw(Exception*) {}
	virtual void solve_equation(T* res, T* rhs, int num = 1) throw(Exception*) = 0;		
	static LinearSolver<T>* factory(unsigned int size_); 
	virtual SparseMatrixInterface<T>& get_matrix() = 0;
			
};

template<>
LinearSolver<cplx>* LinearSolver<cplx>::factory(unsigned int size_); 
template<>
LinearSolver<double>* LinearSolver<double>::factory(unsigned int size_);

/** the no solver class
 *  
 * if only matrix assembly is needed, the no-solver class provides CSR matrices
 * and access to it. used for matrix element and coulomb potential calculation
 */
template<class T>
class NoSolver : public LinearSolver<T> {
public:	
	NoSolver(unsigned int size, SparseMatrixProperties matrix_type);
	virtual ~NoSolver() throw(Exception*) {}		
	virtual void prepare() throw(Exception*) {}
	virtual void release() throw(Exception*) {}
	virtual void solve_equation(T* res, T* rhs, int num = 1) throw(Exception*) {
		TDKP_GENERAL_EXCEPTION("NoSolver class does not solve problems ...");	
	}			 
	virtual SparseMatrixInterface<T>& get_matrix() { return matrix; }
	virtual CSRMatrix<T>& get_csr_matrix() { return matrix; }
private:
	CSRMatrix<T> matrix;			
};

template<class T>
NoSolver<T>::NoSolver(unsigned int size, SparseMatrixProperties matrix_type) 
 : matrix(size, matrix_type) {
	matrix.set_block_size(1);
}

} // end namespace

#endif /*LINEARSOLVER_H_*/

