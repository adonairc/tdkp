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

#ifndef SOCKETIO_H_
#define SOCKETIO_H_

#ifndef NOREMOTESOLVER

#include <vector>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <errno.h>

#include "tdkp/common/all.h"
#include "tdkp/main/CSRMatrix.h"

namespace tdkp {

static const int packet_size = 512;

int  read_int_from_socket(int sockfd);
void send_int_to_socket(int sockfd, int num);

template<class T>
void send_array_to_socket(int sockfd, const T* array, int size) {
    // send array in predefined chungs
    int length;
    int sent = 0;
    int current_chunk;
    int len = sizeof(T);
    int max_size_per_send = packet_size / len;
    int problems = 0;

    while (sent < size) {
        current_chunk = size - sent;
        if (current_chunk > max_size_per_send) current_chunk = max_size_per_send;
        length = send(sockfd, (void*)&(array[sent]), len * current_chunk, 0);
        if (length != len * current_chunk) {
        	if(length < 0) {
        		TDKP_GENERAL_EXCEPTION("send_array_to_socket failed. tried to send " << len * current_chunk << " while return value was " << length << ", errno is " << errno);
        	} else {
        		problems++;
        		sent += length;			
        	}
        } else {
        	sent += current_chunk;	
        }        
    }
    if(problems > 0) {
    	TDKP_LOGMSG(LOG_WARN, "send_array_to_socket had some problems " << problems << " times sending all data at once");	
    }
}

template<class T> 
void read_array_from_socket(int sockfd, T* array, int size) {
    // read array in predefined chunks
    int length;
    int received = 0;
    int expect;
    int len = sizeof(T);
    int max_size_per_recv = packet_size / len;

    while (received < size) {
        expect = size - received;
        if (expect > max_size_per_recv) expect = max_size_per_recv;

        length = recv(sockfd, (void*)&(array[received]), len * expect, MSG_WAITALL);
        if (length != len * expect) TDKP_GENERAL_EXCEPTION("read_array_from_socket failed: recv " << length);

        received += expect;
    }
}

template<class T> 
void send_matrix_to_socket(int sockfd, CSRMatrix<T> *mat, int row_start, int row_end) {
	T *nonzeros = mat->get_nonzeros();
	int *icol   = mat->get_icol();
	int *prow   = mat->get_prow();
	int fidx    = mat->get_fidx();

	// send numnz
	int numnz = prow[row_end + 1] - prow[row_start];
	send_int_to_socket(sockfd, numnz);

	// send fidx
	send_int_to_socket(sockfd, fidx);

	// send part of the nonzeros array
	send_array_to_socket(sockfd, nonzeros + prow[row_start] - fidx, numnz);

	// send icol array
	send_array_to_socket(sockfd, icol + prow[row_start] - fidx, numnz);

	// send prow array
	send_array_to_socket(sockfd, prow + row_start, row_end - row_start + 2);
}

template<class T>
void send_matrix_to_socket(int sockfd, const CSRMatrix<T>* mat) {

	const T *nonzeros = mat->get_nonzeros();
	const int *icol   = mat->get_icol();
	const int *prow   = mat->get_prow();
	int fidx    = mat->get_fidx();
	int sym     = mat->symmetric() ? 1 : 0;
	int numnz   = mat->get_num_nonzeros();
	int size    = mat->get_size();
	
	// send symmetric or not
	send_int_to_socket(sockfd, sym);
		
	// send size
	send_int_to_socket(sockfd, size);		
		
	// send numnz
	send_int_to_socket(sockfd, numnz);

	// send fidx
	send_int_to_socket(sockfd, fidx);

	// send nonzeros array
	send_array_to_socket(sockfd, nonzeros, numnz);

	// send icol array
	send_array_to_socket(sockfd, icol, numnz);

	// send prow array
	send_array_to_socket(sockfd, prow, size + 1);
			
}

template<class T>
void read_matrix_from_socket(int sockfd, SparseMatrixInterface<T>* mat) {

	TDKP_ASSERT(mat->init_from_csr_available(), "matrix type does not support initialization from csr data");
	
	// get symmetric or not
	int sym = read_int_from_socket(sockfd);
		
	// get size
	int size = read_int_from_socket(sockfd);		
		
	// get numnz
	int numnz = read_int_from_socket(sockfd);

	// read fidx
	int fidx = read_int_from_socket(sockfd);
	fidx++; // get rid of warning
	
	vector<int> icol(numnz);
	vector<int> prow(size + 1);
	vector<T>   nonzeros(numnz);

	// read nonzeros array
	read_array_from_socket(sockfd, &nonzeros[0], numnz);

	// read icol array
	read_array_from_socket(sockfd, &icol[0], numnz);

	// read prow array
	read_array_from_socket(sockfd, &prow[0], size + 1);

	// init matrix
	mat->init_from_csr(sym == 1, size, numnz, &prow[0], &icol[0], &nonzeros[0]);
					
}

	
}

#endif // NOREMOTESOLVER

#endif /*SOCKETIO_H_*/
