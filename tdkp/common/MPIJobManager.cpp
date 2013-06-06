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

#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/MPIJobManager.h"
#include "tdkp/solvers/PetscSolver.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

namespace tdkp {

const char* MPIJobManager::help = "master mpi job manager started petsc ...";
MPIJobManager* MPIJobManager::mpijobmanager_singleton = 0;

#ifdef ENABLE_MPI

MPIJobManager::MPIJobManager(int* pargc_, char*** pargv_) 
: petsc_started(false),
  pargc(pargc_),
  pargv(pargv_) 
{
	
	int my_rank;
	int num_proc;
	
	// -------------------------------------
	// init mpi
	// -------------------------------------
	MPI_Init(pargc,pargv);
		
	// -------------------------------------
	// get my rank and num proc
	// -------------------------------------
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	// -------------------------------------
	// return to main thread for rank 0
	// -------------------------------------
	if(my_rank == 0) {
		if(num_proc > 1) {
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "MPIJobManager: tdkp has been started via mpi using " << num_proc << " parallel processes");			
		}	
		return;
	} 
}


void MPIJobManager::mpi_loop() {
	
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Request request;
	MPI_Status  status;	
	if(my_rank > 0) {
		single_mode_operation = false;
		// --------------------------------------------------
		// infinite job loop: receive calculation requests
		// --------------------------------------------------
		int  task = 0;
		bool loop = true;
		PetscSolverClient* petsc_client = 0;		
		while(loop) {
			// wait for request loop (sleeping between tests -> reduces load on sleeping nodes from 100% to irrelevant low values)			
     		MPI_Irecv(&task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
			int flag = 0;
			MPI_Test(&request, &flag, &status);
			while(flag == 0) {
				usleep(10000);
				MPI_Test(&request, &flag, &status);
			}			
			switch(task) {				
#ifdef LINSOLV_INCLUDE_PETSC				
				case 1:
					petsc_client = new PetscSolverClient();
					delete petsc_client; petsc_client = 0;
					break;
#endif					
				case 0:
					loop = false;
					break;
				default:
					TDKP_GENERAL_EXCEPTION("unknown task " << task << " received. don't know what to do so i'll just quit");
			}
		}
		if(petsc_started) {
			PetscErrorCode 	ierr;		// error code
			ierr = PetscFinalize(); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;	
		}
		MPI_Finalize();
		exit(0);
						
	}
	
}

void MPIJobManager::start_petsc() {
	if(!petsc_started) {
		PetscInitialize(pargc,pargv,(char*)0,help);
		petsc_started = true;		
	}
}	

MPIJobManager& MPIJobManager::get_instance() {

	TDKP_ASSERT(mpijobmanager_singleton != 0, "MPIJobManager::get_instance() can only be called after the job manager has been properly created");
	return *mpijobmanager_singleton; 	
}

void MPIJobManager::create_instance(int* argc, char*** argv) {	
	
	TDKP_ASSERT(mpijobmanager_singleton == 0, "create instances may only be called once!");	
	mpijobmanager_singleton = new MPIJobManager(argc,argv);	
}

MPIJobManager::~MPIJobManager() {
	
}

void MPIJobManager::finalize() {

	int my_rank;
	int num_proc;

	// -------------------------------------------
	// thread 0 finishes here
	// other threads never reach the point here
	// -------------------------------------------
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	
	if(my_rank == 0) {
		// --------------------------------------------
		// send quit signal to other processes
		// --------------------------------------------
		int stop_signal = 0;
		vector<MPI_Request> requests(num_proc);
		vector<MPI_Status>  statuses(num_proc);
		for(int ii = 1; ii < num_proc; ii++) {
     		MPI_Isend(&stop_signal, 1, MPI_INT, ii, 0, MPI_COMM_WORLD, &requests[ii]);
		}
		MPI_Waitall(num_proc - 1, &requests[1], &statuses[1]);				
		if(petsc_started) {
			PetscErrorCode 	ierr;		// error code
			ierr = PetscFinalize(); TDKP_ASSERT(ierr == 0, "petsc function returned nonzero " << ierr);;	
		}		
		MPI_Finalize();	
	} else {
		TDKP_GENERAL_EXCEPTION("MPIJobManager: rank " << my_rank << " != 0 should never reach finalize!");	
	}
}

#else

// --------------------------------------------------
// no mpi, don't do anything
// --------------------------------------------------
MPIJobManager::MPIJobManager(int* argc, char*** argv) {}
MPIJobManager::~MPIJobManager() {}
void MPIJobManager::finalize() {}

void MPIJobManager::create_instance(int* argc, char*** argv) {
	mpijobmanager_singleton = new MPIJobManager(argc,argv);	
}
MPIJobManager& MPIJobManager::get_instance() {
	TDKP_ASSERT(mpijobmanager_singleton != 0, "MPIJobManager::get_instance() can only be called after the job manager has been properly created");
	return *mpijobmanager_singleton;	
} 

void MPIJobManager::mpi_loop() {
	
}


#endif /* ENABLE_MPI */

}
