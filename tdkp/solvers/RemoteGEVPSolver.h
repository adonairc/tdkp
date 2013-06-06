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

#ifndef REMOTEGEVPSOLVER_H_
#define REMOTEGEVPSOLVER_H_

#include "tdkp/common/all.h"
#include "tdkp/solvers/EigenSolver.h"
#include "tdkp/probdefs/EigenProblem.h"

namespace tdkp {

class RemoteGEVPSolver;

#ifndef NOREMOTESOLVER

/** thread safe control object for non-thread safe eigenvalue solver 
 * 
 * arpack is not thread safe but i need to parallelize it properly.
 * therefore i decided to assemble the matrices locally and create
 * some remote solvers in their own variable space communicating
 * using sockets
 */
class RemoteGEVPSolverController {
public:
	RemoteGEVPSolverController( 
		const CSRMatrix<cplx>&   A, 
		const CSRMatrix<double>& M, 
		Ordering order, 
		unsigned int nev,
		unsigned int max_parallel
	);
	virtual ~RemoteGEVPSolverController();
	
	bool  solve();
	void  set_results_to_problem_object(EigenProblem<complex<double>, complex<double>, double>& problem);
	const string get_remote_host() const { return remote_host; }
	
private:
	Ordering order;
	vector<cplx> eigenvalues;
	vector<cplx> eigenvectors;

	CSRMatrix<cplx>   A_copy;
	CSRMatrix<double> M_copy;
	
	unsigned int max_parallel;
	unsigned int my_number;
	
	// may only be modified using mutex_solve_queue
	static unsigned int parallel_last_started;
	static unsigned int parallel_next_number; 
	static unsigned int parallel_running;
		
	string remote_host;		
};

class RemoteGEVPSolver : public EigenSolver<cplx,double,cplx> {
public:	
	RemoteGEVPSolver(unsigned int size, unsigned int block_size);
	virtual ~RemoteGEVPSolver();

	// -------------------------------------------------------
	// inherited but forbidden functions
	// -------------------------------------------------------
	virtual bool find_eigenvectors() throw(Exception*) { TDKP_GENERAL_EXCEPTION("remote solver must be controlled via control objects"); }
	virtual int converged_eigenvalues() const { TDKP_GENERAL_EXCEPTION("remote solver must be controlled via control objects"); }
	virtual cplx eigenvalue(int nn) const throw(Exception*) { TDKP_GENERAL_EXCEPTION("remote solver must be controlled via control objects"); }
	virtual const cplx& eigenvector(int nn, int vv) const throw(Exception*) { TDKP_GENERAL_EXCEPTION("remote solver must be controlled via control objects"); }
	 
	// -------------------------------------------------------
	// matrix access
	// ------------------------------------------------------- 
	virtual SparseMatrixInterface<cplx>&   get_stiff() { return *A; }
	virtual SparseMatrixInterface<double>& get_overlap() { return *M; }
	
	// -------------------------------------------------------
	// controller generation
	// -------------------------------------------------------
	RemoteGEVPSolverController* get_controller(Ordering order, unsigned int nev);
		
	
private:
	CSRMatrix<cplx>*   A;
	CSRMatrix<double>* M;
				
};

class RemoteGEVPSolverConnectionHandler {
public:

	static RemoteGEVPSolverConnectionHandler& get_instance();

	// -------------------------------------------------------
	// controller / base class communication
	// ------------------------------------------------------- 
	int get_and_lock_available_solver(string& remote_host_name);
	void release_solver(int client_fd); 
	unsigned int get_max_parallel_solvers() const { return max_parallel_solves; }

private:

	static RemoteGEVPSolverConnectionHandler* singleton;

	RemoteGEVPSolverConnectionHandler();
	~RemoteGEVPSolverConnectionHandler();

	const string&  get_remotehost(unsigned int current_index) const;
	bool           read_remotehosts_list(const string& filename);
	vector<string> remotehosts;
	unsigned int   max_parallel_solves;
	
	void start_parallel_solvers();
	void disconnect_parallel_solvers();
			
	bool bind_and_listen(int& sockfd, int& port) const;
	vector<int>  parallel_solver_sockets;
	vector<int>  parallel_solver_client_fd;
	vector<bool> parallel_solver_in_use;
	vector<int>  parallel_solver_statistics;

	
};

#else
// -----------------------------------------------------
// dummy implementation required by sebastian for negf on
// manno palu (ibm cray) as cray dislikes the use of 
// pthreads, openmp and sockets
// -----------------------------------------------------
class RemoteGEVPSolverController {
public:
	RemoteGEVPSolverController( 
		const CSRMatrix<cplx>&   A, 
		const CSRMatrix<double>& M, 
		Ordering order, 
		unsigned int nev,
		RemoteGEVPSolver& starter,
		unsigned int max_parallel
	) {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual ~RemoteGEVPSolverController() {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}	
	bool solve() {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	void set_results_to_problem_object(EigenProblem<complex<double>, complex<double>, double>& problem) {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	const string& get_remote_host() const {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");		
	}

private:

		
};

class RemoteGEVPSolver : public EigenSolver<cplx,double,cplx> {
public:	
	RemoteGEVPSolver(unsigned int size_, unsigned int block_size) : EigenSolver<cplx,double,cplx>(size_) {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual ~RemoteGEVPSolver() {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual bool find_eigenvectors() throw(Exception*) {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual int converged_eigenvalues() const {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual cplx eigenvalue(int nn) const throw(Exception*) {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual const cplx& eigenvector(int nn, int vv) const throw(Exception*) {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual SparseMatrixInterface<cplx>&   get_stiff() {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	virtual SparseMatrixInterface<double>& get_overlap() {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
	RemoteGEVPSolverController* get_controller(Ordering order, unsigned int nev) {
		TDKP_GENERAL_EXCEPTION("RemoteGEVPSolverController: dummy object disabled at compile time");	
	}
private:
			
};

#endif

}

#endif /*REMOTEGEVPSOLVER_H_*/
