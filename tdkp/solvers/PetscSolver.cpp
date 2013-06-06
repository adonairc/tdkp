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

#include "tdkp/solvers/PetscSolver.h"
#include "tdkp/common/MPIJobManager.h"

#ifdef LINSOLV_INCLUDE_PETSC

extern "C" {
#include <metis.h>
}

namespace tdkp {


PetscSolver::PetscSolver(unsigned int size)
: matrix(size * 2, nonsymmetric_matrix),
  matrix_interface(matrix),
  petsc_solver_client(0),
  rhs_real(size * 2),
  res_real(size * 2),
  nnz_cache(-1) 
{
	Logger::get_instance()->emit(LOG_INFO, "PetscSolver: using real valued sparse iterative petsc solver.");
}

PetscSolver::~PetscSolver() throw(Exception*) {
	if(petsc_solver_client != 0) {
		delete petsc_solver_client; petsc_solver_client = 0;
	} 
}

void PetscSolver::prepare() throw(Exception*) {
	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "PetscSolver: starting prepare");
	// -----------------------------------
	// create petsc solver client and let him handle everything
	// -----------------------------------
	if(petsc_solver_client == 0) {
		petsc_solver_client = new PetscSolverClient(matrix.get_size());
	}

	// -----------------------------------
	// ensure that sparsity pattern didnt change
	// -----------------------------------
	if(nnz_cache == -1) {
		nnz_cache = matrix.get_num_nonzeros();	
	} else {
		// sparsity pattern is not permitted to change, so quit here the number of nonzeros changed
		TDKP_ASSERT(nnz_cache == (signed)matrix.get_num_nonzeros(), "");
	}	
	// -----------------------------------
	// send matrix
	// -----------------------------------
	petsc_solver_client->set_matrix(matrix);
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "PetscSolver: matrix set and distributed");	

}

void PetscSolver::solve_equation(cplx* res, cplx* rhs, int num) throw(Exception*) {
				
	TDKP_ASSERT(num == 1, "sorry, but i can only solve one equation at time");
	TDKP_ASSERT(petsc_solver_client != 0, "petsc solver client has not yet been initialized. call prepare first");
		
	int size  = this->matrix.get_size() / 2;
	#pragma omp parallel for default(shared)
	for(int ii = 0; ii < size; ii++) {
		this->res_real[ii * 2] = this->res_real[ii * 2 + 1] = 0.0;
		this->rhs_real[ii * 2] = rhs[ii].real();
		this->rhs_real[ii * 2 + 1] = rhs[ii].imag();			
	}
	petsc_solver_client->solve_equation(&res_real[0], &rhs_real[0]);
	
	// copy back
	#pragma omp parallel for default(shared)
	for(int ii = 0; ii < size; ii++) {
		res[ii] = cplx(this->res_real[2 * ii],this->res_real[2*ii+1]); 	
	}
	TimeMeasurements::get_instance().track_memory_usage();		
							
}



// ------------------------------------------------------
// client implementation
// ------------------------------------------------------

/** Petsc solver client constructor without params, to be called by slave threads */
PetscSolverClient::PetscSolverClient()
: preallocated(false) 
{

	// -----------------------------------------
	// initialize petsc (this is only done once)
	// -----------------------------------------
	MPIJobManager::get_instance().start_petsc();	

	// -----------------------------------------
	// init
	// -----------------------------------------
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	this->init();
		
	// -----------------------------------------
	// solver loop
	// -----------------------------------------
	bool loop = true;
	while(loop) {
		// receive instruction from master thread
		int task = -1;
		recv_task_from_master(&task);		
		switch(task) {
			// 0 = quit,	
			case 0:
				loop = false;
				break;		 
			// 1 = solve
			case 1:
				this->solve_equation_mpi();
				break;
			// 2 = new matrix
			case 2:
				this->send_matrix_slave();
				break;
			default:
				TDKP_GENERAL_EXCEPTION("unexpected task " << task << " received from master thread");
		}
	}

}
	
/** Petsc solver client constructor, to be called by master process */
PetscSolverClient::PetscSolverClient(unsigned int size_)
: preallocated(false),
  size(size_) 
{

	// -----------------------------------------
	// get mpi info
	// -----------------------------------------
	int num_proc = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	TDKP_ASSERT(rank == 0,"");
	
	// -----------------------------------------
	// establish connection to slaves, start
	// petsc clients 
	// -----------------------------------------
	int mpi_job_manager_task = 1; // leads MPIJobManager to start the petsc solver client on rank > 0 threads
	send_task_to_slaves(mpi_job_manager_task);
		
	// -----------------------------------------
	// initialize petsc
	// -----------------------------------------
	MPIJobManager::get_instance().start_petsc();		
			
	// -----------------------------------------
	// create dummy vector for vector value insertion
	// -----------------------------------------
	dummy_vector_insert_idx.resize(size);
	for(int ii = 0; ii < size; ii++) {
		dummy_vector_insert_idx[ii] = ii;		
	}			
			
	// -----------------------------------------
	// init vectors, matrices and solvers
	// -----------------------------------------
	this->init();
		
}

/** master threads argument to initiate slaves for different action 
 * 
 * the non-blocking sends and receives together with usleep on the
 * slaves highly reduces the load on waiting threads, allowing a 
 * better mix of the openmp parallel stuff on the master thread
 * and the mpi slave nodes
 */
void PetscSolverClient::send_task_to_slaves(int task) const {
	
	int num_proc;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	TDKP_ASSERT(rank == 0, "");
	
	vector<MPI_Request> requests(num_proc);
	vector<MPI_Status>  statuses(num_proc);
	for(int ii = 1; ii < num_proc; ii++) {
     	MPI_Isend(&task, 1, MPI_INT, ii, 0, MPI_COMM_WORLD, &requests[ii]);
	}
	MPI_Waitall(num_proc - 1, &requests[1], &statuses[1]);
		
}

void PetscSolverClient::recv_task_from_master(int* task) const {
	
	TDKP_ASSERT(rank > 0, "");
	MPI_Request request;
	MPI_Status  status;
	// wait for request loop (sleeping between tests -> reduces load on sleeping nodes from 100% to irrelevant low values)			
	MPI_Irecv(task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
	int flag = 0;
	MPI_Test(&request, &flag, &status);
	while(flag == 0) {
		usleep(100);
		MPI_Test(&request, &flag, &status);
	}
				
}

void PetscSolverClient::init() {

	// -----------------------------------------
	// send/recv matrix size
	// -----------------------------------------
	MPI_Bcast(&size, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	TDKP_ASSERT(size > 0, "");

	// -----------------------------------------
	// init vectors, matrices and contexts
	// -----------------------------------------
	ierr = VecCreate(PETSC_COMM_WORLD,&x); 		TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE, size); 	TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;
	ierr = VecSetFromOptions(x);                TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = VecDuplicate(x,&b); 					TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;
			
	// -----------------------------------------
	// create matrices
	// -----------------------------------------
	ierr = MatCreate(PETSC_COMM_WORLD, &A); 						TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, size, size); 	TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = MatSetFromOptions(A);									TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);	
	ierr = MatSetOption(A, MAT_SYMMETRIC); 							TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = MatSetOption(A, MAT_COLUMNS_SORTED); 					TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);	
	
	// -----------------------------------------
	// init solvers
	// -----------------------------------------
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); 	TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = KSPSetType(ksp, KSPGMRES);           TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = KSPGetPC(ksp, &pc); 					TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = PCSetType(pc, PCASM);                TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	// hypre segfaults with > 24 processors, using ASM instead (which performs quite well)
	//ierr = PCHYPRESetType(this->pc, "pilut");   TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	// "pilut" (default), "parasails", "boomeramg", "euclid" 
	//ierr = PetscOptionsSetValue("-pc_hypre_pilut_tol","1.e-6"); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);		
	ierr = KSPSetTolerances(
				ksp, 
				Configuration::get_instance()->get("petsc_ksp_tolerance"), 
				PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT
		   );  TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;
	ierr = KSPSetFromOptions(ksp); 				TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	 		
	// -----------------------------------------
	// get local rows
	// -----------------------------------------	
	ierr = VecGetOwnershipRange(this->x, &row_start, &row_end);	TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;		 
	local_row_idxs.resize(row_end - row_start);
	res_cache.resize(row_end - row_start);
	for(int ii = row_start; ii < row_end; ii++) {
		local_row_idxs[ii - row_start] = ii;
	}
	 
	// TDKP_TRACE("juhu"); 		
}

PetscSolverClient::~PetscSolverClient() {
		
	// --------------------------------------
	// quit petsc
	// --------------------------------------
	if(rank == 0) {
		int stop = 0;
		send_task_to_slaves(stop);
	}
	
	// -----------------------------------------------
	// delete petsc stuff
	// -----------------------------------------------
	ierr = VecDestroy(x); 		TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = VecDestroy(b);  		TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = MatDestroy(A);  		TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	ierr = KSPDestroy(ksp); 	TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
			
}
	
	
void PetscSolverClient::set_matrix(CSRMatrix<double>& matrix) {
	
	int num_proc = 0;
	MPI_Comm_size(PETSC_COMM_WORLD, &num_proc);	
		
	// --------------------------------------------------
	// let slave petsc solvers enter send_matrix_slave
	// --------------------------------------------------
	int task = 2; // (slaves enter matrix distribution mode)
	send_task_to_slaves(task);
		
	slave_row_start.assign(num_proc, 0);
	slave_row_end.assign(num_proc, 0);	
	
	// ---------------------------------------------------
	// distribute matrix
	// ---------------------------------------------------
	MPI_Status status;
	int* prow        = matrix.get_prow();
	int* icol        = matrix.get_icol();
	double* nonzeros = matrix.get_nonzeros();
	
	int frame[2] = {0,0};
	for(int ii = 1; ii < num_proc; ii++) {
		// get rows from slave thread
		MPI_Recv(frame, 2, MPI_INT, ii, 0, PETSC_COMM_WORLD, &status);
		slave_row_start[ii] = frame[0];
		slave_row_end[ii]   = frame[1];
		// calculate length of array we send
		int sending_num = prow[frame[1]] - prow[frame[0]];
		// send length of arrays
		MPI_Send(&sending_num, 1, MPI_INT, ii, 0, PETSC_COMM_WORLD);
		// send prow array
		MPI_Send(&prow[frame[0]], frame[1] - frame[0] + 1, MPI_INT, ii, 0, PETSC_COMM_WORLD);
		// send icol array
		MPI_Send(&icol[prow[frame[0]]], sending_num, MPI_INT, ii, 0, PETSC_COMM_WORLD);
		// send nonzeros array
		MPI_Send(&nonzeros[prow[frame[0]]], sending_num, MPI_DOUBLE, ii, 0, PETSC_COMM_WORLD);						
	}
	// TDKP_TRACE("juhu");
	// --------------------------------------------------
	// init local rows (the master processes local rows)
	// --------------------------------------------------
	set_local_rows(&prow[row_start], &icol[prow[row_start]], &nonzeros[prow[row_start]]);
	// TDKP_TRACE("juhu");

}
		

void PetscSolverClient::send_matrix_slave() {

	// --------------------------------------------
	// send master thread my rows
	// --------------------------------------------
	int frame[2]; frame[0] = row_start; frame[1] = row_end;	
	MPI_Send(frame, 2, MPI_INT, 0, 0, PETSC_COMM_WORLD); 
	
	// --------------------------------------------
	// receive matrix
	// --------------------------------------------
	MPI_Status     status;
	vector<int>    prow;
	vector<int>    icol;
	vector<double> nonzeros;
	int local_nnz;
	// recv length of arrays
	MPI_Recv(&local_nnz, 1, MPI_INT, 0, 0, PETSC_COMM_WORLD, &status);
	// init arrays
	prow.resize(row_end - row_start + 1);
	icol.resize(local_nnz);
	nonzeros.resize(local_nnz);
	// recv prow array
	MPI_Recv(&prow[0], row_end - row_start + 1, MPI_INT, 0, 0, PETSC_COMM_WORLD, &status);
	// recv icol array
	MPI_Recv(&icol[0], local_nnz, MPI_INT, 0, 0, PETSC_COMM_WORLD, &status);
	// recv nonzeros array
	MPI_Recv(&nonzeros[0], local_nnz, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD, &status);		
	// --------------------------------------------------
	// init local rows (slave process local rows)
	// --------------------------------------------------	
	set_local_rows(&prow[0], &icol[0], &nonzeros[0]);
		
}

/** set local matrix values */
void PetscSolverClient::set_local_rows(
	const int* local_prow, const int* local_icol, const double* local_nonzeros
) {
	// ----------------------------------------------
	// set matrix preallocation
	// ----------------------------------------------
	// for every row
	// TDKP_TRACE("juhu");
	if(!preallocated) {
		// TDKP_TRACE("juhu");
		vector<int> d_nnz(row_end - row_start, 0);
		vector<int> o_nnz(row_end - row_start, 0);
		for(int ii = 0; ii < row_end - row_start; ii++) {
			// calculate diagonal box values
			int diagonal    = 0;
			int offdiagonal = 0;
			// for all entries in row
			for(int jj = local_prow[ii]; jj < local_prow[ii + 1]; jj++) {
				// check if this is inside the diagonal submatrix
				if(local_icol[jj - local_prow[0]] >= row_start && local_icol[jj - local_prow[0]] < row_end) {
					diagonal++;	
				} else {
					offdiagonal++;	
				}			
			}	
			// store
			d_nnz[ii] = diagonal;
			o_nnz[ii] = offdiagonal;
		}
		// TDKP_TRACE("juhu");
		// announce in matrix
		ierr = MatMPIAIJSetPreallocation(A, 0, &d_nnz[0], 0, &o_nnz[0]); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;
		preallocated = true;
	}
	// -------------------------------------------------
	// insert values
	// -------------------------------------------------
	// for every local row	
	int num_proc = 0;
	MPI_Comm_size(PETSC_COMM_WORLD, &num_proc);
	int col_entries = 0;	
	for(int ii = 0; ii < row_end - row_start; ii++) {				
		col_entries = local_prow[ii + 1] - local_prow[ii];			
		int II = ii + row_start;	 				
		ierr = MatSetValues(A, 1, &II, col_entries,	&local_icol[local_prow[ii] - local_prow[0]], 						
				&local_nonzeros[local_prow[ii]  - local_prow[0]], INSERT_VALUES); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);    					
	}		
	// TDKP_TRACE("juhu");

	// ---------------------------------------------------
	// assemble matrix
	// ---------------------------------------------------
	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;
	ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;
	// TDKP_TRACE("juhu");
	// ---------------------------------------------------
	// set operator to solvers
	// ---------------------------------------------------
	ierr = KSPSetOperators(ksp,A,A, SAME_NONZERO_PATTERN);	
	// TDKP_TRACE("juhu");
}

void PetscSolverClient::solve_equation_mpi() {
	
	// -----------------------------------------
	// assemble rhs (b has been set by master thread)
	// -----------------------------------------
	// TDKP_TRACE("juhu");
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
	// TDKP_TRACE("juhu");

	// -----------------------------------------
	// solve ...
	// -----------------------------------------
	ierr = KSPSolve(ksp, b, x); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;
	
	// TDKP_TRACE("juhu");
	// ------------------------------------------
	// extract values
	// send values to master 
	// ------------------------------------------
	if(rank > 0) {		
		ierr = VecGetValues(x, local_row_idxs.size(), &local_row_idxs[0], &res_cache[0]); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);		
		MPI_Send(&res_cache[0], res_cache.size(), MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD);
	}
	// TDKP_TRACE("juhu");	
}
			
/** master process solve_equation */
void PetscSolverClient::solve_equation(double* res, const double* rhs) throw(Exception*) {

	// -------------------------------------
	// set rhs. values 
	// -------------------------------------
	double t1,t2,t3,t4,t5,t6;
	t1 = TimeMeasurements::tic();
	ierr = VecSetValues(b, size, &dummy_vector_insert_idx[0], rhs, INSERT_VALUES);
	t2 = TimeMeasurements::tic(); 	
	TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
		// TDKP_TRACE("juhu");
	// -------------------------------------
	// tell slaves to enter solve_equation_mpi
	// -------------------------------------
	int task = 1;
	send_task_to_slaves(task);
	t3 = TimeMeasurements::tic();
	// TDKP_TRACE("juhu");				
	// -------------------------------------
	// solve equation
	// -------------------------------------
	solve_equation_mpi();
	t4 = TimeMeasurements::tic();
	// TDKP_TRACE("juhu");
	// -------------------------------------
	// assemble solution
	// -------------------------------------
	// put master solution into the array
	ierr = VecGetValues(x, local_row_idxs.size(), &local_row_idxs[0], &res[row_start]); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);
	t5 = TimeMeasurements::tic();
	// receive slaves
	MPI_Status status;
	int num_proc = 0;
	MPI_Comm_size(PETSC_COMM_WORLD, &num_proc);	
	for(int pp = 1; pp < num_proc; pp++) {
		// recv slave data
		MPI_Recv(&res[slave_row_start[pp]], slave_row_end[pp] - slave_row_start[pp], MPI_DOUBLE, pp, 0, PETSC_COMM_WORLD, &status);	
	}	
	t6 = TimeMeasurements::tic();
	TDKP_LOGMSG(LOG_INFO, "\nPetscSolver: timings\n 1-6 = " << t6 - t1 << ",\n 1-2 = " << t2 - t1 << ",\n 2-3 = " << t3 - t2 << ",\n 3-4 = " << t4 - t3 << ",\n 4-5 = " << t5 - t4 << "\n 5-6 = " << t6 - t5);
	// TDKP_TRACE("juhu");	
}	
	

} /** namespace tdkp */

#endif /* LINSOLV_INCLUDE_PETSC */
