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

#include "tdkp/solvers/RemoteGEVPSolver.h"
#include "tdkp/utilities/SocketIO.h"
#include <errno.h>
#include <unistd.h>
#include <omp.h>

#ifndef NOREMOTESOLVER

namespace tdkp {

// ---------------------------------------------
// static mutex initializer
// i need a waiting queue for the parallel solvers
// but there may be some parallel assemblies running
// so i have to do this globally
// ---------------------------------------------
pthread_mutex_t mutex_solve_queue = PTHREAD_MUTEX_INITIALIZER;;

// ---------------------------------------------
// the next mutex is required when creating sockets,
// because socket functions are not thread safe
// ---------------------------------------------
pthread_mutex_t remote_gevp_socket_mutex = PTHREAD_MUTEX_INITIALIZER;

RemoteGEVPSolver::RemoteGEVPSolver(unsigned int size_, unsigned int block_size_)
: EigenSolver<cplx,double,cplx>(size_),
  A(0),
  M(0)
{
	if(Configuration::get_instance()->get("assembly_build_nonsymmetric_matrices") == 0.0 ||
	   Configuration::get_instance()->get("remote_solver_umfpack_upper_size_limit") < size_) {
		this->A = new CSRMatrix<cplx>(size_, symmetric_matrix);
		this->M = new CSRMatrix<double>(size_ / block_size_, symmetric_matrix);
	} else {
		Logger::get_instance()->emit(LOG_INFO, "building nonsymmetric matrices");
		this->A = new CSRMatrix<cplx>(size_, nonsymmetric_matrix);
		this->M = new CSRMatrix<double>(size_ / block_size_, nonsymmetric_matrix);		
	}
	// make sure the connections are set up 
	RemoteGEVPSolverConnectionHandler::get_instance();	
}
RemoteGEVPSolver::~RemoteGEVPSolver() {

	if(this->A != 0) {
		delete A; A = 0;
	}
	if(this->M != 0) {
		delete M; M = 0;	
	}			
}
	
RemoteGEVPSolverController* RemoteGEVPSolver::get_controller(Ordering order_, unsigned int nev_) {	
	return new RemoteGEVPSolverController(
		*A, *M, order_, nev_, RemoteGEVPSolverConnectionHandler::get_instance().get_max_parallel_solvers()
	);	
}

unsigned int RemoteGEVPSolverController::parallel_running      = 0;
unsigned int RemoteGEVPSolverController::parallel_last_started = 0;
unsigned int RemoteGEVPSolverController::parallel_next_number  = 1; 

RemoteGEVPSolverController::RemoteGEVPSolverController( 
	const CSRMatrix<cplx>&   A, 
	const CSRMatrix<double>& M, 
	Ordering order_, 
	unsigned int nev,
	unsigned int max_parallel_
) : order(order_),
    A_copy(A),
    M_copy(M),
    max_parallel(max_parallel_),
    my_number(0) // 0 means not started yet
{
	TDKP_ASSERT(max_parallel > 0, "");
	TDKP_ASSERT(nev > 0,"");	
	// allocate result storage
	eigenvalues.assign(nev, 0.0);
	eigenvectors.assign(nev * A.get_size(), 0.0);
}

RemoteGEVPSolverController::~RemoteGEVPSolverController() {}
	
bool RemoteGEVPSolverController::solve() {

	bool wait_in_queue = true;
	int  clnt_fd       = -1;
	while(wait_in_queue) { 
		// lock mutex (so im the only one)
		pthread_mutex_lock(&mutex_solve_queue);
		// ----------------------------------------
		// get my number in the queue
		// ----------------------------------------
		if(my_number == 0) {
			my_number = parallel_next_number;
			parallel_next_number++;	
		}		
		// ----------------------------------------
		// go if there is some spare cpu power
		// ----------------------------------------
		if(parallel_running < max_parallel && my_number == parallel_last_started + 1) {
			++parallel_running;
			parallel_last_started = my_number;
			wait_in_queue = false;
			// -------------------------------------------
			// get socket connection
			// -------------------------------------------
			clnt_fd = RemoteGEVPSolverConnectionHandler::get_instance().get_and_lock_available_solver(remote_host);
			TDKP_ASSERT(clnt_fd != -1, "expected to get valid connection to atze remote client from my connection keeper, but received " << clnt_fd << ". sorry, this is unexpected and i have to terminate."); 			
			pthread_mutex_unlock(&mutex_solve_queue);				
		} else {
			pthread_mutex_unlock(&mutex_solve_queue);
			usleep(10000); // wait 0.01 sec	
		}
	}
								 
	try {
		
		// send a 1 as hello to atze remote gevp (0 terminates atze)
		int keepalive = 1;
		send_int_to_socket(clnt_fd, keepalive);
				
		// send if matrices are symmetric or nonsymmetric
		int symmetric = A_copy.symmetric() ? 1 : 0;
		send_int_to_socket(clnt_fd, symmetric);
			
		// send matrix sizes
		send_int_to_socket(clnt_fd, A_copy.get_size());
		send_int_to_socket(clnt_fd, M_copy.get_size()); 
			
		// send matrix A and M
		send_matrix_to_socket(clnt_fd, &A_copy);
		send_matrix_to_socket(clnt_fd, &M_copy);
	
		// send number of eigenvalues
		send_int_to_socket(clnt_fd, static_cast<int>(eigenvalues.size()));
		
		// send ordering
		if(order == ascending) {
			send_int_to_socket(clnt_fd, 0);
		} else {
			send_int_to_socket(clnt_fd, 1);
		} 
		
		// recv number of obtained eigenvalues
		int nev_obtained = read_int_from_socket(clnt_fd);
		
		if((unsigned)nev_obtained != eigenvalues.size()) {
			TDKP_LOGMSG(LOG_WARN, "only " << nev_obtained << " out of " << eigenvalues.size() << " eigenvalues converged!");
	    	pthread_mutex_lock(&mutex_solve_queue);
	    	--parallel_running;
	    	pthread_mutex_unlock(&mutex_solve_queue);						
			return false;
		}
			
		// recv eigenvalues
		read_array_from_socket(clnt_fd, &eigenvalues[0], eigenvalues.size());
		
		// recv eigenvectors
		read_array_from_socket(clnt_fd, &eigenvectors[0], eigenvectors.size());
			
		// deregister connection
		pthread_mutex_lock(&mutex_solve_queue);
		--parallel_running;
		RemoteGEVPSolverConnectionHandler::get_instance().release_solver(clnt_fd);		
		pthread_mutex_unlock(&mutex_solve_queue);
		
	} catch (Exception *e) {
		TDKP_LOGMSG(LOG_WARN, e->get_reason());
    	pthread_mutex_lock(&mutex_solve_queue);
    	--parallel_running;
    	pthread_mutex_unlock(&mutex_solve_queue);		
		return false;	
	}
	
	return true;
}

void RemoteGEVPSolverController::set_results_to_problem_object(EigenProblem<complex<double>, complex<double>, double>& problem) {
	// add solutions to problem
	int size = eigenvectors.size() / eigenvalues.size();
	for(unsigned int ii = 0; ii < eigenvalues.size(); ii++) {
		problem.add_solution(eigenvalues[ii], &eigenvectors[ii * size], size);
	}	
}

RemoteGEVPSolverConnectionHandler& RemoteGEVPSolverConnectionHandler::get_instance() {
	if(singleton == 0) {
		singleton = new RemoteGEVPSolverConnectionHandler();	
	}	
	return *singleton;
}

RemoteGEVPSolverConnectionHandler::RemoteGEVPSolverConnectionHandler() 
: max_parallel_solves(1)
{
	// ----------------------------------
	// test for atze_remote_gevp
	// ----------------------------------
#ifdef TCSH
	TDKP_ASSERT(system("which atze_remote_gevp >& /dev/null") == 0,"could not find atze_remote_gevp executable in path! (which atze_remote_gevp failed)");
#else
	TDKP_ASSERT(system("which atze_remote_gevp >/dev/null 2&>1") == 0,"could not find atze_remote_gevp executable in path! (which atze_remote_gevp failed)");
#endif
	// ----------------------------------
	// get environment data
	// ----------------------------------
	int config_max_solves = static_cast<int>(Configuration::get_instance()->get("remote_solver_max_number_of_parallel_solves"));
	if(config_max_solves == -1) {
		max_parallel_solves = omp_get_max_threads();
	} else {
		// use config value
		max_parallel_solves = config_max_solves;
	} 	
	TDKP_ASSERT(config_max_solves != 0, "moron, remote_solver_max_number_of_parallel_solves is set to 0");
	string remote_hosts("remote_hosts.dat");
	if(getenv("TDKP_REMOTE_HOSTS_FILE") != 0) {
		remote_hosts = getenv("TDKP_REMOTE_HOSTS_FILE");		 		
	}
		
	// ----------------------------------
	// load remote hostname list
	// ----------------------------------	
	if(Configuration::get_instance()->get("remote_solver_read_hostname_list") == 1.0) {
	   	TDKP_LOGMSG(LOG_INFO_DEVEL2, "RemoteGEVPSolver: using " << remote_hosts << " as remote host file");
		if(!this->read_remotehosts_list(remote_hosts)) {
			TDKP_LOGMSG(LOG_WARN, "you requested me to read remote hosts from remote_hosts.dat but the file does not exist or is mal-formatted. using localhost instead.");
			remotehosts.push_back(string("localhost"));	 
		} else {
			if(config_max_solves == -1) {
				// use size of list		
				max_parallel_solves = remotehosts.size();				
			} 
		}
	} else {
		remotehosts.push_back(string("localhost"));		
	}
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "RemoteGEVPSolver: using " << max_parallel_solves << " parallel eigenvalue solvers");

	// -------------------------------------------------------
	// start parallel eigensolvers
	// -------------------------------------------------------	
	this->start_parallel_solvers();
					 	
	
	
}
RemoteGEVPSolverConnectionHandler::~RemoteGEVPSolverConnectionHandler() {
	//cout << "DISCONNECTING\n";
	this->disconnect_parallel_solvers();	
}

/** bind and listen to socket */
bool RemoteGEVPSolverConnectionHandler::bind_and_listen(int& sockfd, int& port) const {

    pthread_mutex_lock(&remote_gevp_socket_mutex);
	int yes = 1;
	int max_port = static_cast<int>(Configuration::get_instance()->get("remote_solver_port_max"));

    // create socket
    if ((sockfd = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
        TDKP_GENERAL_EXCEPTION("cannot create socket");
    }
    
    // set reusable 
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1) {
        TDKP_GENERAL_EXCEPTION("cannot set setsockopt");
    }

    struct sockaddr_in my_addr;
    // init my address
    my_addr.sin_family = AF_INET;			// host byte order    
    my_addr.sin_addr.s_addr = INADDR_ANY;	// automatically fill with my IP
    
    int bind_ret = 0;
    
    do {
    	port++;
    	// init my address
    	my_addr.sin_family = AF_INET;			// host byte order    
    	my_addr.sin_addr.s_addr = INADDR_ANY;	// automatically fill with my IP    	
		my_addr.sin_port = htons(port);			// short, network byte order    
		memset(&(my_addr.sin_zero), '\0', 8);	// zero the rest of the struct    	
		if(port > max_port) {
			pthread_mutex_unlock(&remote_gevp_socket_mutex);
			return false;
		}
		bind_ret = bind(sockfd, (struct sockaddr *)&my_addr, sizeof(struct sockaddr));
		if(bind_ret == -1 && errno != EADDRINUSE) {
			if(errno == EADDRINUSE) {
				TDKP_LOGMSG(LOG_WARN, "bind says that the socket");
			} else {
				TDKP_LOGMSG(LOG_WARN, "bind returned errno " << errno << "!");
			}	
		}
    } while (bind_ret == -1);
        
    // start to listen 
    if (listen(sockfd, 1) == -1) {
    	if(errno == EADDRINUSE) {
    		TDKP_LOGMSG(LOG_WARN, "listen failed because another socket is already listening ");		
    	} else {
    		TDKP_LOGMSG(LOG_WARN, "listen failed with errno " << errno);	
    	}   
    	pthread_mutex_unlock(&remote_gevp_socket_mutex); 	    	
    	return false;
    }
    pthread_mutex_unlock(&remote_gevp_socket_mutex);
    return true;
	
}

RemoteGEVPSolverConnectionHandler* RemoteGEVPSolverConnectionHandler::singleton = 0;

void RemoteGEVPSolverConnectionHandler::start_parallel_solvers() {

	TDKP_ASSERT(max_parallel_solves > 0, "");
	
	// ---------------------------------------
	// init arrays
	// ---------------------------------------
	parallel_solver_in_use.assign(max_parallel_solves, false);
	parallel_solver_client_fd.assign(max_parallel_solves, -1);
	parallel_solver_sockets.assign(max_parallel_solves, -1);
	parallel_solver_statistics.assign(max_parallel_solves, 0);	
	const int port_min = static_cast<int>(Configuration::get_instance()->get("remote_solver_port_min"));
	
	// ---------------------------------------
	// for every remote process
	// --------------------------------------- 
	for(unsigned int ii = 0; ii < max_parallel_solves; ii++) {
		string remotehost = get_remotehost(ii);
		// ---------------------------------------
		// bind to socket
		// ---------------------------------------
	    int tries = 0;
		int port  = port_min; 
	    const int max_tries = 10;
	    while(!bind_and_listen(parallel_solver_sockets[ii],port) && tries < max_tries) {
	    	TDKP_LOGMSG(LOG_WARN, "could not listen to port " << port << ". try " << ++tries);
	    	port = port_min;
	    	close(parallel_solver_sockets[ii]);
	    	sleep(500);	
	    } 		
	    if(tries == max_tries) {
	    	close(parallel_solver_sockets[ii]);
	    	TDKP_GENERAL_EXCEPTION("could not establish socket connection after several attempts. terminating!");	    		
	    }	    
	
		// ---------------------------------------
		// start remote atzes
		// ---------------------------------------
	    // output file
	    ostringstream remote_log_file;
	    if(Configuration::get_instance()->get("remote_solver_create_logfiles") == 1.0) {
	    	remote_log_file << "remote." << port << ".out";
	    } else {
	    	remote_log_file << "/dev/null ";
	    }
			
		// create remote solver
		ostringstream rcmd;	
		if(remotehost == string("localhost")) {
#ifdef TCSH
			rcmd << "atze_remote_gevp localhost " << port << " >& " << remote_log_file.str() << " &";
#else
			rcmd << "atze_remote_gevp localhost " << port << " >" << remote_log_file.str() << " 2>&1 &";
#endif
		} else {
			// get my hostname
			char buf[1024];
			TDKP_ASSERT(gethostname(buf,1024) == 0, "could not determine my hostname"); 
#ifdef TCSH
			rcmd << "rsh " << remotehost << " \"atze_remote_gevp " << buf << " " << port << " >& " << remote_log_file.str() << " &\" >& /dev/null &";
#else
			rcmd << "rsh " << remotehost << " \"atze_remote_gevp " << buf << " " << port << " >" << remote_log_file.str() << " 2>&1 &\" >/dev/null 2>&1 &";				
#endif	
		}
		int ret = system(rcmd.str().c_str());
		TDKP_ASSERT(ret == 0, "starting atze via '" << rcmd.str() << "' failed");	
	}	
	
	// ---------------------------------------
	// accept connections
	// ---------------------------------------		
    // wait for clients
    struct sockaddr_in clnt_addr;
    socklen_t sin_size;    
    sin_size = sizeof(struct sockaddr_in);
    for(unsigned int ii = 0; ii < max_parallel_solves; ii++) {
	    if((parallel_solver_client_fd[ii] = accept(parallel_solver_sockets[ii], (struct sockaddr *)&clnt_addr, &sin_size)) == -1) {
	    	close(parallel_solver_sockets[ii]);	    	
			TDKP_GENERAL_EXCEPTION("accept connection " << ii << " of atze remote gevp failed. this is unexpected, so i terminate.");	
	    }
    }	
	
	// ---------------------------------------
	// ready to rumble
	// ---------------------------------------
	TDKP_LOGMSG(LOG_INFO, "RemoteGEVPSolver: successfully started " << max_parallel_solves << " remote atzes");
		
}
void RemoteGEVPSolverConnectionHandler::disconnect_parallel_solvers() {
	
	int term = 0;
	ostringstream sout;
	sout << "RemoteGEVPCH: remote solver statistics:";	     
	for(unsigned int ii = 0; ii < max_parallel_solves; ii++) {
		// terminate if client is still in use (which should not be)
		TDKP_ASSERT(!parallel_solver_in_use[ii], "remote atze " << ii << " is still in use while i tried to kill it.");		
		if(parallel_solver_client_fd[ii] > 0) {
			send_int_to_socket(parallel_solver_client_fd[ii], term);
			sout << "\n" << "  slot " << setw(3) << ii << " solved " <<  setw(7) << parallel_solver_statistics[ii] << " on " << get_remotehost(ii);  				
		}
	}
	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());
	
	usleep(200000); // wait 0.20 sec
	
	
	for(unsigned int ii = 0; ii < max_parallel_solves; ii++) {
		// terminate client
		if(parallel_solver_client_fd[ii] > 0) {						
			close(parallel_solver_client_fd[ii]);
			parallel_solver_client_fd[ii] = -1;
		}
		// close socket
		if(parallel_solver_sockets[ii] > 0) {
			close(parallel_solver_sockets[ii]);
			parallel_solver_sockets[ii] = -1;	
		}
	}			
}


int RemoteGEVPSolverConnectionHandler::get_and_lock_available_solver(string& remote_host_name) {
	for(unsigned int ii = 0; ii < max_parallel_solves; ii++) {
		if(!parallel_solver_in_use[ii]) {
			parallel_solver_in_use[ii] = true;
			remote_host_name = get_remotehost(ii);
			parallel_solver_statistics[ii]++;
			return parallel_solver_client_fd[ii];	
		}	
	}
	return -1;
}

void RemoteGEVPSolverConnectionHandler::release_solver(int client_fd) {
	// find corresponding idx
	for(unsigned int ii = 0; ii < max_parallel_solves; ii++) {
		if(parallel_solver_client_fd[ii] == client_fd) {
			TDKP_ASSERT(parallel_solver_in_use[ii], "tried to release parallel solver " << ii << ". but the solver isn't marked as in use.");
			parallel_solver_in_use[ii] = false;
			break;				
		}	
	}		
} 


const string& RemoteGEVPSolverConnectionHandler::get_remotehost(unsigned int current_index) const {
	return remotehosts[current_index % remotehosts.size()];
}

bool RemoteGEVPSolverConnectionHandler::read_remotehosts_list(const string& filename) {
	ifstream fin(filename.c_str());
	if(fin) {
		const int max_length = 8000;
		int added = 0;
		char cline[max_length];
		// read entries by line
		while(fin.getline(cline, max_length)) {
			string line(cline);			
			// skip lines starting with #
			if(line[0] != '#') {
				// find semicolon
				int ii = 0;
				while(ii < (signed)line.size() && line[ii] != ':' && line[ii] != '\n' && line[ii] != '\0') {
					ii++;	
				}
				// find newline or end
				int nn = 0;
				while(nn < (signed)line.size() && line[nn] != '\n' && line[nn] != '\0') {
					nn++;
				}
				string host;
				int ntimes = 1;
				// no semicolon						
				if(ii == nn) {
					host = line.substr(0, nn);
				} else {
				// with semicolon
					host = line.substr(0, ii);
					ntimes = lexical_cast<int>(line.substr(ii + 1, nn - ii - 1));
				}
				if(host.size() > 0) {
					for(int cc = 0; cc < ntimes; cc++) {
						remotehosts.push_back(host);
						added++;	
					}
				}
			} 	
		}
		if(added == 0) {
			TDKP_LOGMSG(LOG_WARN, "no valid host in remote host file.");
			return false;	
		} else {
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "using " << remotehosts.size() << " remote host entries");
			return true;	
		}				
	} else {
		return false;	
	}	
}



}

#endif // NOREMOTESOLVER

