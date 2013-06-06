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

#ifndef EIGENSOLVER_H_
#define EIGENSOLVER_H_

#include "tdkp/common/all.h"
#include "tdkp/main/CSRMatrix.h"

// ------------------------------------------
// include or remove eigensolvers
// LAPACK is always included
// ------------------------------------------
//#define EIGNSOLV_INCLUDE_ARPACK
//#define EIGNSOLV_INCLUDE_JDQZ
//#define EIGNSOLV_INCLUDE_JACOBI_DAVIDSON

namespace tdkp {


/** base class for eigensolvers. T denote the matrix type (cplx/double) and V denotes the resulting vectors
 * 
 * e.g. solvers may take real matrices, but return complex results
 * eigensolvers are created using the factory pattern. by using the configuration
 * object, the factories determine the requested type of eigensolver
 * 
 * the eigensolver again uses the factory pattern to get an appropriate 
 * linear solver. 
 * 
 * we do not care about the real matrix format but only use the
 * sparse matrix interface. thus i try to make a flexible integration
 * of trilinos software routines ...  
 */
template<class A, class M, class V>
class EigenSolver {
public:


	EigenSolver(unsigned int size);
	virtual ~EigenSolver();
	virtual void assign(int nev, EigenProblemType type);
	virtual void set_ordering(Ordering order_) { this->order = order_; }
	virtual Ordering get_ordering() const { return this->order; };

	virtual bool find_eigenvectors() throw(Exception*) = 0;
	virtual int converged_eigenvalues() const;
	virtual cplx eigenvalue(int nn) const throw(Exception*) = 0;
	virtual const V& eigenvector(int nn, int vv) const throw(Exception*) = 0;
	 
	// -------------------------------------------------------
	// factory and abstraction
	// -------------------------------------------------------
	static EigenSolver<A,M,V>* factory(unsigned int size_, unsigned int block_size); 
	virtual SparseMatrixInterface<A>& get_stiff() = 0;
	virtual SparseMatrixInterface<M>& get_overlap() = 0;
		
	EigenProblemType get_current_problem_type() const { return this->problem_type; }
				
protected:

			
	void sort_solutions(int num_calculated, V* &myvecs, V* &myvals);
	int  nev;
	int  msize;
	int  converged_ev;

private:
				
	Ordering         order;
	EigenProblemType problem_type;
						
};

/** base class constructor */
template<class A, class M, class V>
EigenSolver<A, M, V>::EigenSolver(unsigned int size)
: nev(0),
  msize(size),
  converged_ev(0),
  order(descending),
  problem_type(indefinite)
{	
	
}



template<class A, class M, class V>
EigenSolver<A, M, V>::~EigenSolver() {

}

/** assign generalized eigenvalue problem
 * 
 * to calculate the eigenvalue of Ax = lambda Bx
 * 
 * @param nev_   number of eigenvalues to calculate
 * @param stiff_ the A in the eigenvalue calculation
 * @param overlap_ the B in the eigenvalue calculation
 */
template<class A, class M, class V>
void EigenSolver<A, M, V>::assign(int nev_, EigenProblemType type) {		
	this->nev          = nev_;	
	this->msize        = this->get_stiff().get_size();
	this->problem_type = type;
}

/** returns the number of converged eigenvalues
 */
template<class A, class M, class V>
int EigenSolver<A, M, V>::converged_eigenvalues() const {
	return this->converged_ev;	
}


template<class A, class M, class V> 
void EigenSolver<A, M, V>::sort_solutions(int num_calculated, V* &myvecs, V* &myvals) {
	
	if(!myvecs || !myvals) {
		TDKP_GENERAL_EXCEPTION("sort_solutions called on empty evec/eval");
	}		
	
	TDKP_ASSERT(num_calculated >= this->nev, "num_calculated >= this->nev");
			
	int*   idx_order         = new int[num_calculated];
   	int    pos;
   	V*     new_eigenvalues  = new cplx[num_calculated];
   	V*     new_eigenvectors = new V[num_calculated * this->get_stiff().get_size()];
   	double ascend            = (this->order == ascending ? 1.0 : -1.0);
   	int    size              = (signed)this->get_stiff().get_size();
   
   	// build idx order
   	for(int ii = 0; ii < num_calculated; ii++) {
      	pos = 0;
      	for(int jj = 0; jj < num_calculated; jj++) {
         	if(ii != jj) {
            	if(tdkp_math::only_real(myvals[ii]) * ascend > tdkp_math::only_real(myvals[jj]) * ascend 
               		|| (tdkp_math::only_real(myvals[ii]) == tdkp_math::only_real(myvals[jj]) && ii > jj)) {
               		pos++;
            	}
         	}
      	}
      	idx_order[ii] = pos;
   	}
	
   	for(int ii = 0; ii < num_calculated; ii++) {   		        		 
      	new_eigenvalues[idx_order[ii]] = myvals[ii];
#pragma omp parallel for default(shared)      	
      	for(int jj = 0; jj < size; jj++) {
         	new_eigenvectors[idx_order[ii] * size + jj] = myvecs[ii * size + jj];
      	}
   	}

   	delete[] idx_order;
   	delete[] myvals;
   	myvals = new_eigenvalues;
   	delete[] myvecs;
   	myvecs = new_eigenvectors;
			
}



} // end namespace

#endif /*EIGENSOLVER_H_*/
