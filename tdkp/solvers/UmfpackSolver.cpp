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


#include "tdkp/solvers/UmfpackSolver.h"

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING UMFPACK
// -----------------------------------------------
#ifdef LINSOLV_INCLUDE_UMFPACK

#include "umfpack.h"

namespace tdkp {

UmfpackSolverComplex::UmfpackSolverComplex(unsigned int size)
: matrix(size), // no fortran ordering
  UmfpackNumeric(0),
  UmfpackSymbolic(0),
  rhs_real(size),
  rhs_imag(size),
  res_real(size),
  res_imag(size),
  csr_structure(0)  
{

}

UmfpackSolverComplex::~UmfpackSolverComplex()  throw(Exception*) {
	if(UmfpackNumeric != 0) {
		umfpack_zl_free_numeric(&UmfpackNumeric);
		UmfpackNumeric = 0;
	}
	if(UmfpackSymbolic != 0) {
		umfpack_zl_free_symbolic(&UmfpackSymbolic);
	}
}

void UmfpackSolverComplex::csr2csc() {

	const unsigned int size = matrix.get_size();
	const unsigned int nnz  = matrix.get_number_of_nonzeros();

	const int*    icol     = matrix.get_icol();
	const int*    prow     = matrix.get_prow();
	const double* row_real = matrix.get_nz_real();
	const double* row_imag = matrix.get_nz_imag();
	ostringstream sout;
	sout << "UmfpackSolverComplex: reordering csr to csc";
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());

	// -------------------------------------------------
	// structure principally changed?
	// -------------------------------------------------
	bool structure_changed = this->csr_structure != matrix.get_current_structure();
	if(this->csr_structure != matrix.get_current_structure()) {

		// -----------------------------------------------
		// resort csr to csc
		// -----------------------------------------------
		// calculate elements in each colum
		vector<int> col_count;
		col_count.assign(size, 0);
		for(unsigned int ii = 0; ii < nnz; ii++) {
			col_count[icol[ii]]++;
		}
		// init new space
		csc_irow.assign(nnz, 0);
		csc_pcol.assign(size + 1, 0);
		csc_real.assign(nnz, 0.0);
		csc_imag.assign(nnz, 0.0); 
						
		// build pcol from col_count (sum up col counts, gives me number of elements in each col)
		for(unsigned int ii = 0; ii < size; ii++) {
			csc_pcol[ii + 1] = csc_pcol[ii] + col_count[ii];	
		}	
		this->csr_structure = matrix.get_current_structure();
		
	}
	
	TimeMeasurements::get_instance().start("csr2csc");
	vector<int> offsets(size);	
	// reorder arrays
	int from, until, csc_index;		 
	for(int rr = 0; rr < (signed)size; rr++) {
		from = prow[rr]; until = prow[rr + 1];
		for(int ii = from; ii < until; ii++) {					
			// so,  (rr,icol[ii]) is the element i have now
			// to save this in csc, i do the following:
			// 1. get starting index of colum icol[ii] in csc			
			csc_index = csc_pcol[icol[ii]] + offsets[icol[ii]];			
			// 2. store values at starting index + offset
			csc_real[csc_index] = row_real[ii];
			csc_imag[csc_index] = row_imag[ii];
			// 3. store current row in irow at index + offset
			csc_irow[csc_index] = rr;
			// 4. increase offset
			offsets[icol[ii]]++;		
		} 
	}
	TimeMeasurements::get_instance().stop("csr2csc");

	// --------------------------------------------------
	// perform symbolic factorization
	// --------------------------------------------------
	if(structure_changed) {
		if(UmfpackSymbolic != 0) {
			umfpack_zl_free_symbolic(&UmfpackSymbolic);
		}		
		double* null = (double*) NULL;
		int     size = matrix.get_size();		
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, "UmfpackSolverComplex: new matrix structure, symbolic factorization phase needed (reordering)");
		umfpack_zl_symbolic(size,size, &csc_pcol[0], &csc_irow[0], &csc_real[0], &csc_imag[0], &UmfpackSymbolic, null, null);
	}	
	
}

void UmfpackSolverComplex::prepare()   throw(Exception*){
	
	this->csr2csc();
	
	double* null = (double*) NULL;
	TDKP_ASSERT(this->csc_real.size() == matrix.get_number_of_nonzeros(), "csc not properly initialized");
	TDKP_ASSERT(this->csc_irow.size() == matrix.get_number_of_nonzeros(), "csc not properly initialized");
	TDKP_ASSERT(this->csc_pcol.size() == matrix.get_size() + 1,  "csc not properly initialized");	

	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "UmfpackSolverComplex: numeric factorization phase");
	umfpack_zl_numeric(&csc_pcol[0], &csc_irow[0], &csc_real[0], &csc_imag[0], UmfpackSymbolic, &UmfpackNumeric, null, null);
	
}

void UmfpackSolverComplex::solve_equation(cplx* res, cplx* rhs, int num)  throw(Exception*) {
	TDKP_ASSERT(num == 1, "sorry, but i can only solve one equation at time");
	double* null = (double*) NULL;
	// --------------------------------------
	// copy complex values into double arrays
	// --------------------------------------
	const int size = matrix.get_size();
#pragma omp parallel for default(shared)
	for(int ii = 0; ii < size; ii++) {
		rhs_real[ii] = rhs[ii].real();
		rhs_imag[ii] = rhs[ii].imag();
	}
	res_real.assign(size, 0.0e0);
	res_imag.assign(size, 0.0e0);

	umfpack_zl_solve(
		UMFPACK_A,
		&csc_pcol[0], &csc_irow[0], &csc_real[0], &csc_imag[0],
		&res_real[0], &res_imag[0],
		&rhs_real[0], &rhs_imag[0],
		UmfpackNumeric,
		null,
		null
	);	
	// ----------------------------------------
	// copy result back
	// ----------------------------------------
#pragma omp parallel for default(shared)
	for(int ii = 0; ii < size; ii++) {
		res[ii] = cplx(res_real[ii],res_imag[ii]);
	}

}

}


#endif
