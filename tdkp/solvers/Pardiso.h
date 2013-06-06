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

#ifndef PARDISO_H_
#define PARDISO_H_



#include "tdkp/common/all.h"
#include "tdkp/main/CSRMatrix.h"
#include "tdkp/solvers/LinearSolver.h"

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING PARDISO
// -----------------------------------------------
#ifdef LINSOLV_INCLUDE_PARDISO

namespace tdkp {

// ---------------------------------------------------
// extern fortran function
// ---------------------------------------------------
extern "C" {

// ---------------------------------------------------	
// deprecated pardiso 3.0 interface 	
// ---------------------------------------------------
#ifdef PARDISO30

	int pardiso_(void* pt[64],     /* pardiso internal memory pointer's */
				  int* maxfct ,    /* max. number of factorized matrices having the SAME structure that should be kept in mem*/
				  int* mnum,       /* actual matrix for the solution phase */
				  int* mtype,      /* matrix type */
				  int* phase,      /* control of the solver (fill in reduction, factorization, fwd/backwrd subst, term) */
				  int* n,          /* size of A */
				  void* A,         /* the matrix ... in sparse format */
				  int* prow,       /* sparse pointer to rows */
				  int* icol,       /* colum indices of sparse entries in A */
				  int* perm,       /* user fill in reducing permutation ... */
				  int* nrhs,       /* number of right hand sides that need to be solved for */
				  int  iparam[64], /* pardiso control array */
				  int* msglvl, 	   /* message level information */
				  void* b,         /* b(n,nrhs) rhs vector/matrix */
				  void* x,         /* x(n,nrhs) */
				  int* error);     /* error */

	int pardisoinit_(
		void* pt[64],
		int*  mtype,
		int iparam[64]
	);
#else
	int pardiso_(void* pt[64],     /* pardiso internal memory pointer's */
				  int* maxfct ,    /* max. number of factorized matrices having the SAME structure that should be kept in mem*/
				  int* mnum,       /* actual matrix for the solution phase */
				  int* mtype,      /* matrix type */
				  int* phase,      /* control of the solver (fill in reduction, factorization, fwd/backwrd subst, term) */
				  int* n,          /* size of A */
				  void* A,         /* the matrix ... in sparse format */
				  int* prow,       /* sparse pointer to rows */
				  int* icol,       /* colum indices of sparse entries in A */
				  int* perm,       /* user fill in reducing permutation ... */
				  int* nrhs,       /* number of right hand sides that need to be solved for */
				  int  iparam[64], /* pardiso control array */
				  int* msglvl, 	   /* message level information */
				  void* b,         /* b(n,nrhs) rhs vector/matrix */
				  void* x,         /* x(n,nrhs) */
				  int* error,      /* error */
				  double dparm[64]); /* numerical parameters */ 

	int pardisoinit_(
		void* pt[64],
		int*  mtype,
		int*  solver,
		int iparam[64],
		double dparm[64],
		int* error		
	);
#endif
}


/** Pardiso wrapper class for CSR matrices
 *
 * may be used for matrices of type
 * - real and structurally symmetric
 * - real and symmetric indefinite
 * - complex and hermitian indefinite
 * - complex and structurally symmetric
 * */
template<class T>
class Pardiso : public LinearSolver<T> {

public:
	Pardiso(unsigned int size) throw(Exception*);
	virtual ~Pardiso() throw(Exception*);
	virtual void prepare() throw(Exception*);
	virtual void solve_equation(T* res, T* rhs, int num = 1) throw(Exception*);

	virtual SparseMatrixInterface<T>& get_matrix() { return *this->matrix; }

private:

	CSRMatrix<T>* matrix;

	void*  handle[64];
	int    maxfct;
	int    mnum;
	int    mtype;
	int    phase;
	int    n;
	int*   perm;
	int    nrhs;
	int    iparam[64];
	double dparam[64];
	int    msglvl;
	int    error;

	bool  factorized;

	void  handle_error();
	int   determine_mtype() const;

	int*  icol;
	int*  prow;
	vector<int> icol_space;
	vector<int> prow_space;

};

/** constructor for symmetric pardiso direct solver
 *
 * Pardiso is a high end direct linear equation solver for solving Ax = b where A is sparse
 * in our case it's used for the eigensolvers which from time to time need to solve an
 * equation system using a matrix several times.
 *
 * @mat complex lhs matrix
 */
template<class T>
Pardiso<T>::Pardiso(unsigned int size_) throw(Exception*) {

	const Configuration* controls = Configuration::get_instance();

	// init iparam and handle to zero
	for(int ii = 0; ii < 64; ii++) {
		this->handle[ii] = 0;
		this->iparam[ii] = 0;
	}

	if(Configuration::get_instance()->get("assembly_build_nonsymmetric_matrices") == 0.0) {
		this->matrix = new CSRMatrix<T>(size_, symmetric_matrix);
	} else {
		Logger::get_instance()->emit(LOG_INFO, "building nonsymmetric matrices");
		this->matrix = new CSRMatrix<T>(size_, nonsymmetric_matrix);
	}

	this->maxfct = 1; // keep only single matrix
	this->mnum   = 1;
	this->mtype  = this->determine_mtype(); // determine type of matrix

#ifdef PARDISO30
	pardisoinit_(this->handle, &mtype, iparam);
#else
	int solver = 0;
	pardisoinit_(this->handle, &mtype, &solver, iparam, dparam, &error);
	TDKP_ASSERT(error == 0, "");
#endif	

	this->n      = matrix->get_size();
	this->perm   = new int[this->n];
	this->nrhs   = 1;
	this->msglvl = (int)controls->get("pardiso_message_level");// 1; // yeah, talk to me, pardiso
	this->error  = 0; // just init ...

	// init permutations array to zero
	for(int ii = 0; ii < this->n; ii++) {
		this->perm[ii] = 0;
	}

	this->iparam[0]  = 1;  // not user supplied values
	this->iparam[1]  = 2;  // use metis
	this->iparam[2]  = 1;  // number of threads (default 1 (setting it to omp num threads below))
	this->iparam[3]  = 0;  // don't do preconditioned cgs
	this->iparam[4]  = 0;  // no user permutation
	this->iparam[5]  = 0;  // write solution to X
	this->iparam[7]  = 1;  // max numbers of iterative refinement steps
	this->iparam[9]  = 8;
	this->iparam[10] = 0;
	//this->iparam[11] = 0;
	//this->iparam[12] = 0;
	this->iparam[17] = -1;
	this->iparam[18] = -1;
	this->iparam[20] = 1;
	this->iparam[23] = 1;
	this->iparam[24] = 1;

	char* tmp;
	if((tmp = getenv("OMP_NUM_THREADS")) != 0) {
		this->iparam[2] = atoi(getenv("OMP_NUM_THREADS"));
		ostringstream sout;
		sout << "Pardiso: setting num threads to " << this->iparam[2];
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	}

	this->factorized = false;

	icol = 0;
	prow = 0;

}

/** deallocates all memory (including the memory generated via fortran */
template<class T>
Pardiso<T>::~Pardiso() throw(Exception*) {
	this->phase = -1;
	T* res = new T[this->n];
	T* rhs = new T[this->n];

	pardiso_(
		this->handle, &this->maxfct, &this->mnum, &this->mtype, &this->phase, &this->n,
		(void*)matrix->get_nonzeros(), prow, icol, this->perm, &this->nrhs,
		this->iparam, &this->msglvl, (void*)rhs, (void*)res, &this->error,
#ifndef PARDISO30
		this->dparam
#endif			 
	);

	delete[] res;
	delete[] rhs;
	delete[] this->perm;
	delete this->matrix;
	if(this->error != 0) this->handle_error();
}

/** factor matrix into LU componentes
 */
template<class T>
void Pardiso<T>::prepare() throw(Exception*) {
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, "Pardiso: starting with matrix factorization");
	this->phase = 12; // analysis + numerical factorization
	T* res = new T[this->n];
	T* rhs = new T[this->n];

	// check if matrix is built using fortran ordering
	if(matrix->get_prow()[0] == 0) {
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, "Pardiso: mapping to fortran indices");
		// its c ordering -> create some arrays ...
		int* matrix_icol =  matrix->get_icol();
		int* matrix_prow =  matrix->get_prow();
		icol_space.assign(matrix->get_num_nonzeros(), 0);
		prow_space.assign(matrix->get_size() + 1, 0);
		icol = &icol_space[0];
		prow = &prow_space[0];
		for(int ii = 0; ii < n + 1; ii++) {
			prow[ii] = matrix_prow[ii] + 1;
		}
		int nnz = matrix->get_num_nonzeros();
		for(int ii = 0; ii < nnz; ii++) {
			icol[ii] = matrix_icol[ii] + 1;
		}

	} else {
		prow = matrix->get_prow();
		icol = matrix->get_icol();

	}
	TDKP_ASSERT(prow[n] == (int)matrix->get_num_nonzeros() + 1, "prow[n] " << prow[n] << " == matrix->get_num_nonzeros() (" << matrix->get_num_nonzeros() << ") + 1");
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "Pardiso: starting factorization");
	pardiso_(
		this->handle, &this->maxfct, &this->mnum, &this->mtype, &this->phase, &this->n,
		(void*)matrix->get_nonzeros(), prow, icol, this->perm, &this->nrhs,
		this->iparam, &this->msglvl, (void*)rhs, (void*)res, &this->error,
#ifndef PARDISO30
		this->dparam
#endif					 
	);

	if(Configuration::get_instance()->get("output_eigensolver_statistics") == 1) {
		std::ostringstream sout;
		sout << "pardiso statistics:\n"
		     << "error: " << this->error << "\n"
		     << "number iterative refin steps: " << this->iparam[7] << "\n"
			 << "number pertubed pivots:       " << this->iparam[13] << "\n"
			 << "peak memory:                  " << this->iparam[14] << "\n"
			 << "permanent memory:             " << this->iparam[15] << "\n"
			 << "fact memory:                  " << this->iparam[16] << "\n"
			 << "nonzeros:                     " << this->iparam[17] << "\n";
		Logger::get_instance()->emit(LOG_INFO_DEVEL1, sout.str());
	}
	delete[] res;
	delete[] rhs;
	res = 0;
	if(this->error != 0) this->handle_error();
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, "pardiso finished factoring matrix");
	this->factorized = true;

}

/** using the factorized matrix, determine solutions to linear equation systems
 *
 * via fwd./bkwd substitution
 * allows to determine multiple solutions at same time
 *
 * @param res the x in Ax = b of length N * num while N is the size of matrix A
 * @param arg_rhs the b in Ax = b of length N * num while N is the size of matrix A
 * @param num number of solutions to determine at same time
 * @return on success, the result(s) of the linear equation system stored in res
 */
template<class T>
void Pardiso<T>::solve_equation(T* res, T* arg_rhs, int num) throw(Exception*) {

	// copy arguments
	this->phase = 33;
	pardiso_(this->handle, &this->maxfct, &this->mnum, &this->mtype, &this->phase, &this->n,
		(void*)matrix->get_nonzeros(), prow, icol, this->perm, &num,
		this->iparam, &this->msglvl, (void*)arg_rhs, (void*)res, &this->error,
#ifndef PARDISO30
		this->dparam
#endif					 
	);
	if(this->error != 0) this->handle_error();

}



/** throws exceptions containing the error messages returned from pardiso
 */
template<class T>
void Pardiso<T>::handle_error() {
	if(this->error != 0) {
		string err("pardiso says: ");
		switch(this->error) {
			case -1:
				err.append("input incosistent");
				break;
			case -2:
				err.append("not enough memory");
				break;
			case -3:
				err.append("reordering problem");
				break;
			case -4:
				err.append("zero pivot numerical fact. or iterative refinement problem");
				break;
			case -5:
				err.append("unclassified (internal) error");
				break;
			case -6:
				err.append("preordering failed");
				break;
			case -7:
				err.append("diagonal matrix problem");
				break;
			default:
				err.append("unknown error code");
		}
		TDKP_GENERAL_EXCEPTION(err);
	}
}



} // end of namespace

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING PARDISO
// -----------------------------------------------
#endif

#endif /*PARDISO_H_*/
