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

#ifndef UMFPACKSOLVER_H_
#define UMFPACKSOLVER_H_



#include "tdkp/common/all.h"
#include "tdkp/solvers/LinearSolver.h"
#include "tdkp/main/CSRSparseMatrixInterface.h"

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING UMFPACK
// -----------------------------------------------
#ifdef LINSOLV_INCLUDE_UMFPACK


namespace tdkp {

class UmfpackSolverComplex : public LinearSolver<cplx> {
public:
	UmfpackSolverComplex(unsigned int size);
	virtual ~UmfpackSolverComplex()  throw(Exception*);
	virtual void prepare() throw(Exception*);
	virtual void solve_equation(cplx* res, cplx* rhs, int num = 1) throw(Exception*);
	virtual SparseMatrixInterface<cplx>& get_matrix() { return matrix; }

private:

	TwinComplexCSR matrix;
	void* UmfpackNumeric;
	void* UmfpackSymbolic; // symbolic factorization (reused in subsequent factorizations)	
	vector<double> rhs_real;
	vector<double> rhs_imag;
	vector<double> res_real;
	vector<double> res_imag;

	// TODO: add CSC matrix storage
	// temporary hack: umfpack needs CSC, and not CSR
	void csr2csc();	
	vector<long> 	csc_irow;
	vector<long>	csc_pcol;
	vector<double>	csc_real;
	vector<double>  csc_imag;
	int             csr_structure;
	

};


} // end of namespace tdkp

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING UMFPACK
// -----------------------------------------------
#endif

#endif /*UMFPACKSOLVER_H_*/
