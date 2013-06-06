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


#include <vector>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <netdb.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <linux/kernel.h>
#include <sys/sysinfo.h>
#include <sys/time.h>


#include "tdkp/common/all.h"
#include "tdkp/solvers/ArpackSIPardisoSolver.h"
#include "tdkp/main/CSRMatrix.h"
#include "tdkp/utilities/SocketIO.h"

using namespace std;
using namespace tdkp;

#define stop_right_now(msg) { std::cerr << msg << endl;	abort(); }	


int socket_connect(const char* hostname_string, int dest_port) {

	int    clnt_fd; 
    struct sockaddr_in their_addr;    // my address information
    struct hostent *he;         	
    clnt_fd = 0;
    
    string hostname(hostname_string);
    
  	// --------------------------------------------------------
  	// open socket for client server communication
	// --------------------------------------------------------
	ostringstream sout;  	
	cout << " - - - - - - - - - - - - - - - - - - - - - - - - - \n"
	     << " hello, my name is atze remote gevp, the remote arpack solver\n"
	     << " -> connecting to " << hostname << ":" << dest_port << "\n";
	if ((clnt_fd = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
    	stop_right_now("could not create socket");
	}


    if ((he=gethostbyname(hostname.c_str())) == NULL) {  // get the host info 
		sout << "could not get ip for " << hostname;
		stop_right_now(sout.str().c_str());
	} 
    their_addr.sin_family = AF_INET;    // host byte order 
	their_addr.sin_port = htons(dest_port);  // short, network byte order 
    their_addr.sin_addr = *((struct in_addr *)he->h_addr);
	memset(&(their_addr.sin_zero), '\0', 8);  // zero the rest of the struct 
	if (connect(clnt_fd, (struct sockaddr *)&their_addr,sizeof(struct sockaddr)) == -1) {    
		sout.str("");
		sout << " ******** error ******* - there is no tdkp process on " << hostname << " waiting!";
		stop_right_now(sout.str().c_str());			  
	}
	sout.str("");
	cout << "successfully connected!\n";
	return clnt_fd;
}
  	    	

int main(int argc, char* argv[]) {
	
	if(argc != 3) {
		cout << "usage: atze_remote_gevp.cpp <host> <port>\n";
		return 0; 
	}
	
	setenv("OMP_NUM_THREADS", "1", 1);
	single_mode_operation = false;
	
	Logger::get_instance()->set_level(LOG_INFO_DEVEL2);
	Logger::get_instance()->add_listener(&std::cout);
	
	TDKP_LOGMSG(LOG_INFO, "ATZE REMOTE STARTED");
		
	// open port and connect to tdkp
	int clnt_fd = socket_connect(argv[1], atoi(argv[2]));
	
	TDKP_LOGMSG(LOG_INFO, "ATZE REMOTE CONNECTED TO " << argv[1] << " on port " << atoi(argv[2]));;
		
	// -------------------------------------------------------
	// loop for ever
	// -------------------------------------------------------
	while(true) {		
		
		// expect 1 if we should keep the connection alive and 0 if 
		// the connection should be terminated
		int keepalive = read_int_from_socket(clnt_fd);
		if(keepalive != 1) {
			break;	
			TDKP_LOGMSG(LOG_INFO, "ATZE RECEIVED TERMINATE.");
		}		
		
		// receive matrix properties
		int sym    = read_int_from_socket(clnt_fd);
		int A_size = read_int_from_socket(clnt_fd);
		int M_size = read_int_from_socket(clnt_fd);
		
		TDKP_LOGMSG(LOG_INFO, "ATZE RECEIVED MATRIX WITH SIZES " << A_size << " and " << M_size);;
		
		// set fixed values for linear and eigensolver
		Configuration::get_instance()->set("desired_eigenvalue_solver", 2);		
		if(sym == 1) {
			// pardiso
			Configuration::get_instance()->set("desired_linear_equation_solver", 8);
			Configuration::get_instance()->set("assembly_build_nonsymmetric_matrices", 0.0);
		} else {
			// umfpack
			Configuration::get_instance()->set("desired_linear_equation_solver", 2);
			Configuration::get_instance()->set("assembly_build_nonsymmetric_matrices", 1.0);	
		}	
		
		TDKP_LOGMSG(LOG_INFO, "ATZE REMOTE CREATES SOLVER");;
		
		// create solver object
		ArpackSIPardisoSolver solver(A_size, A_size / M_size);
		
		TDKP_LOGMSG(LOG_INFO, "ATZE REMOTE READS MATRICES");
			
		// recieve matrices
		read_matrix_from_socket(clnt_fd, &(solver.get_stiff()));
		read_matrix_from_socket(clnt_fd, &(solver.get_overlap()));
		
		TDKP_LOGMSG(LOG_INFO,  "ATZE REMOTE READS NEV");		
				
		// send number of requested eigenvalues
		int nev = read_int_from_socket(clnt_fd);
	
		TDKP_LOGMSG(LOG_INFO,  "ATZE REMOTE READS ORDER");
	
		// read ordering
		Ordering order = read_int_from_socket(clnt_fd) == 0 ? ascending : descending;	
	
		TDKP_LOGMSG(LOG_INFO, "ATZE REMOTE WILL SOLVER FOR NEV: " << nev);
			
		// solve eigenproblem
		solver.assign(nev, indefinite);
		solver.set_ordering(order);
		solver.find_eigenvectors();
		int converged = solver.converged_eigenvalues();
		
		TDKP_LOGMSG(LOG_INFO, "ATZE REMOTE CONVERGED " << converged << " AND WILL SEND ");
			
		// send results
		send_int_to_socket(clnt_fd, converged);	
		send_array_to_socket(clnt_fd, solver.get_evals_data(), nev);
		send_array_to_socket(clnt_fd, solver.get_evecs_data(), nev * A_size); 
			
		TDKP_LOGMSG(LOG_INFO, "ATZE REMOTE SUCCSESSFULLY FINISHED");
				
	}		
	// and quit
	close(clnt_fd);
	
}
