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


#include "tdkp/solvers/JDQZSolver.h"

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING JDQZ
// -----------------------------------------------
#ifdef EIGNSOLV_INCLUDE_JDQZ

namespace tdkp {

JDQZSolver<cplx>*   jdqz_global_cplx   = 0;
JDQZSolver<double>* jdqz_global_double = 0;
vector<cplx>        jdqz_tmp_storage;

/** check that only one solvers exists at the time */
void jdqz_check_only_one() {
	TDKP_ASSERT(tdkp::jdqz_global_cplx == 0 || tdkp::jdqz_global_double == 0, "there may only be one jdqz solver living at the time");
	TDKP_ASSERT(tdkp::jdqz_global_cplx != 0 || tdkp::jdqz_global_double != 0, "no jdqz object defined globally!");	
}


/** aquire global variable jdqz_global_cplx*/
template<>
void JDQZSolver<cplx>::aquire() {
	#pragma omp critical 
	{ 
		jdqz_global_cplx = this;
		jdqz_check_only_one();
		jdqz_tmp_storage.assign(get_stiff().get_size(), 0.0);
	}
}
/** release global variable jdqz_global_cplx*/
template<>
void JDQZSolver<cplx>::release() {
	TDKP_ASSERT(jdqz_global_cplx == this, "");
	jdqz_global_cplx = 0;
	jdqz_tmp_storage.clear();		
}
/** aquire global variable jdqz_global_double*/
template<>
void JDQZSolver<double>::aquire() {
	#pragma omp critical 
	{ 	
		jdqz_global_double = this;
		jdqz_check_only_one();
		jdqz_tmp_storage.assign(get_stiff().get_size(), 0.0);
	}	
}
/** release global variable jdqz_global_double*/
template<>
void JDQZSolver<double>::release() {
	TDKP_ASSERT(jdqz_global_double == this, "");
	jdqz_global_double = 0;
	jdqz_tmp_storage.clear();	
}


} // end of namespace

typedef complex<double> cplx;

extern "C" {
	void F77_Name(amul)(int* n, cplx* q, cplx* r);
	void F77_Name(bmul)(int* n, cplx* q, cplx* r);
	void F77_Name(precon)(int* n, cplx* q);
}




/** multiply q with A and store in r
 * 
 * external routine, used by jdqz 
 */
void F77_Name(amul)(int* n, cplx* q, cplx* r) {
	tdkp::jdqz_check_only_one();
	tdkp::SparseMatrixInterface<cplx>* stiff = 0;
	if(tdkp::jdqz_global_cplx != 0) {
		stiff = &tdkp::jdqz_global_cplx->get_stiff();
	} else {
		stiff = &tdkp::jdqz_global_double->get_stiff();
	}
	stiff->mult_vec(q,r);		
};

/** multiply q with B and store in r *
 *
 * external routine, used by jdqz 
 */
void F77_Name(bmul)(int* n, cplx* q, cplx* r) {
	
	tdkp::jdqz_check_only_one();	
	if(tdkp::jdqz_global_cplx != 0) {
		tdkp::jdqz_global_cplx->get_overlap_csr().mult_vec_multiple(q, r, tdkp::jdqz_global_cplx->get_block_size());
	} else {
		tdkp::jdqz_global_double->get_overlap_csr().mult_vec_multiple(q, r, tdkp::jdqz_global_double->get_block_size());
	}	
};

/** precondition q (multiply inv(A)q -> q)
 *
 * external routine, used by jdqz 
 */
void F77_Name(precon)(int* n, cplx* q) {
	tdkp::jdqz_check_only_one();
	
	if(*n != (signed)(tdkp::jdqz_tmp_storage.size())) {
		throw new tdkp::Exception("*n != (signed)(tdkp::jdqz_tmp_storage.size()", __LINE__, __FILE__, __DATE__, __TIME__,__func__);
	}	
	// --------------------------------------
	// copy q into tmp storage
	// --------------------------------------
	memcpy(&tdkp::jdqz_tmp_storage[0], q, sizeof(cplx) * (*n));
	// --------------------------------------
	// solve
	// --------------------------------------
	if(tdkp::jdqz_global_cplx != 0) {
		tdkp::jdqz_global_cplx->get_solver().solve_equation(q, &tdkp::jdqz_tmp_storage[0], 1);
	} else {
		tdkp::jdqz_global_double->get_solver().solve_equation(q, &tdkp::jdqz_tmp_storage[0], 1);
	}	
} 

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING JDQZ
// -----------------------------------------------
#endif

