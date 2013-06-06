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

#ifndef PETSCSOLVER_H_
#define PETSCSOLVER_H_

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING UMFPACK
// -----------------------------------------------
#ifdef LINSOLV_INCLUDE_PETSC
#ifndef ENABLE_MPI
#error "if you enable PETSC, you have to include MPI via -DENABLE_MPI as well!"
#endif

#include "petscksp.h"
#include "tdkp/main/CSRMatrix.h"
#include "tdkp/main/CSRSparseMatrixInterface.h"
#include "tdkp/solvers/LinearSolver.h"

namespace tdkp {

/** petsc client, handling everything from matrix distribution to solving */
class PetscSolverClient {
public:
	PetscSolverClient();
	PetscSolverClient(unsigned int size);
	~PetscSolverClient();
	void set_matrix(CSRMatrix<double>& matrix);
	void solve_equation(double* res, const double* rhs) throw(Exception*);
private:
	void init();
	void send_matrix_slave();
	void solve_equation_mpi();
	void set_local_rows(const int* local_prow, const int* local_icol, const double* local_nonzeros);


	// -------------------------------------
	// petsc stuff
	// -------------------------------------
	Vec 			x,b;   		// solution, rhs 
	PC  			pc;   		// preconditioner context   
	Mat 			A;			// linear system matrix
	KSP 			ksp;		// linear solver context
	PetscErrorCode 	ierr;		// error code
	
	// -------------------------------------
	// other stuff
	// -------------------------------------
	bool preallocated;          // true if matrix has already been preallocated
	int  rank;	
	int  size;					// matrix size
	int  row_start;				// local row start
	int  row_end;				// local row end
	
	vector<int>    local_row_idxs; // local row indices [row_start,row_end[ (means, < row_end), required to collect data
	vector<double> res_cache;      // used in solve_equation_mpi to send result to master thread
	
	// -------------------------------------
	// master process stuff
	// -------------------------------------
	void send_task_to_slaves(int task) const;
	void recv_task_from_master(int* task) const;
	vector<int> slave_row_start;			// starting row index for slave rows
	vector<int> slave_row_end;         		// ending row index for slave rows
	vector<int> dummy_vector_insert_idx; 	// vector for inserting continuous values int petsc vectors
	vector<int> partition;
		
};

/** petsc linear solver class */
class PetscSolver : public LinearSolver<cplx> { 
public:
	PetscSolver(unsigned int size);
	virtual ~PetscSolver() throw(Exception*);
	virtual void prepare() throw(Exception*);
	virtual void solve_equation(cplx* res, cplx* rhs, int num = 1) throw(Exception*);
	virtual SparseMatrixInterface<cplx>& get_matrix() { return matrix_interface; }
private:	
	CSRMatrix<double> matrix;
	C2LocalDCSRSparseMatrixInterface matrix_interface;	
	PetscSolverClient* petsc_solver_client;
	vector<double> rhs_real;
	vector<double> res_real;
	// nnz_chache is used to check whether the matrix changed 
	// (i assume that if the sparsity pattern changes, the number of nnz is changed to)
	// at the moment, changing the sparsity pattern is not permitted
	int nnz_cache;  
		
};

}

#endif /* LINSOLV_INLCUDE_PETSC */

#endif /*PETSCSOLVER_H_*/
