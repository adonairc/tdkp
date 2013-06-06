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

#include "tdkp/solvers/LapackBandSolver.h"

namespace tdkp {

// ---------------------------------------------------
// there are some convergence problems in the eigenvalue
// calculation Ax = lambda M x of lapack if pmls are included. 
// interestingly, i had no problems to get the evs using
// matlab (which is based on lapack). 
// as a result, i lowered the convergence criterion (1e-15)
// in the respective lapack routines
// ---------------------------------------------------
#define USE_TDKP_ZGGEV

extern "C" {
	// -----------------------------------------------
	// blas matrix multiplication routine
	// for the complex nonsymmetric case
	// -----------------------------------------------
	void F77_Name(zgbmv)(
		char* TRANS,
		int* M,
		int* N,
		int* KL,
		int* KU,
		cplx* ALPHA,
		cplx* A,
		int* LDA,
		const cplx* X,
		int* INCX,
		cplx* BETA,
		cplx* Y,
		int* INCY
	);
	// -----------------------------------------------
	// blas matrix multiplication routine
	// for the complex HERMITIAN case
	// -----------------------------------------------
	void F77_Name(ztbmv)(
		char* UPLO, 
		char* TRANS, 
		char* DIAG, 
		int*  N, 
		int*  K, 
		cplx* A, 
		int*  LDA,
		cplx* X, 
		int*  INCX
	);
	// ------------------------------------------------
	// lapack nonsymmetric complex LU factorization expert driver
	// yeah man, i am an expert!
	// ------------------------------------------------
	void F77_Name(zgbsvx)(
		char* FACT, 
		char* TRANS,
		int * N,
		int*  KL,
		int*  KU,
		int*  NRHS,
		cplx* AB,
		int*  LDAB,
		cplx* AFB,
		int*  LDAFB,
		int*  IPIV,
		char* EQUED, 
		double* R,
		double* C,
		cplx*	B,
		int*    LDB,
		cplx*	X, 
		int* 	LDX,
		double*	RCOND,
		double* FERR,
		double* BERR,
		cplx*   WORK, 
		double* RWORK,
		int*    INFO 
	);	

#ifndef NOACML
	// ------------------------------------------
	// acml version of lapack's linear equation solver for full matrices
	// ------------------------------------------
	void zhesvx(char fact, char uplo, int n, int nrhs, cplx *a, int lda, 
			    cplx *af, int ldaf, int *ipiv, cplx *b, int ldb, cplx *x, 
			    int ldx, double *rcond, double *ferr, double *berr, int *info);
	// ------------------------------------------
	// acml lapack's eigenvalue solver for hermitian banded matrices
	// using the acml version
	// ------------------------------------------ 			    
	void zhbgvx(char jobz, char range, char uplo, int n, int ka, int kb, 
		cplx *ab, int ldab, cplx *bb, int ldbb, cplx *q, 
		int ldq, double vl, double vu, int il, int iu, double abstol, int *m, 
		double *w, cplx *z, int ldz, int *ifail, int *info);

	// ------------------------------------------
	// lapack's generalized eigenvalue solver for
	// the 'any kind of matrix'
	// so A and M are both complex -> used for PML
	// ------------------------------------------	
	extern void zggev(char jobvl, char jobvr, int n, cplx *a, int lda, 
		cplx *b, int ldb, cplx *alpha, cplx *beta, 
		cplx *vl, int ldvl, cplx *vr, int ldvr, int *info);
#else
	// ------------------------------------------
	// lapack's linear equation solver for full matrices
	// ------------------------------------------
	extern void F77_Name(zhesvx)(
		char* fact, char* uplo, int* n, int* nrhs, cplx *a, int* lda, 
	    cplx *af, int* ldaf, int* ipiv, cplx* b, int* ldb, cplx* x, 
	    int* ldx, double* rcond, double* ferr, double* berr, 
		cplx* work, int* lwork, double* rwork, int* info
	);
	// ------------------------------------------
	// lapack's eigenvalue solver for hermitian banded matrices
	// ------------------------------------------ 			    
	extern void F77_Name(zhbgvx)(
		char* jobz, char* range, char* uplo, int* n, int* ka, int* kb, 
		cplx* ab, int* ldab, cplx* bb, int* ldbb, cplx* q, 
		int* ldq, double* vl, double* vu, int* il, int* iu, double* abstol, 
		int* m, double* w, cplx* z, int* ldz, cplx* work, double* rwork, int* iwork,
		int* ifail, int* info
	);


#endif

	// direct lapack (MY VERSION) -> MODIFIED QZ TOLERANCE TO ACHIEVE CONVERGENCE
	extern void F77_Name(zggevr)(char* jobvl, char* jobvr, int* n, cplx *a, int* lda, 
		cplx *b, int* ldb, cplx *alpha, cplx *beta, 
		cplx *vl, int* ldvl, cplx *vr, int* ldvr,
		cplx* WORK, int* LWORK, double* RWORK, 
		int *info);
	// direct lapack 
	extern void F77_Name(zggev)(char* jobvl, char* jobvr, int* n, cplx *a, int* lda, 
		cplx *b, int* ldb, cplx *alpha, cplx *beta, 
		cplx *vl, int* ldvl, cplx *vr, int* ldvr,
		cplx* WORK, int* LWORK, double* RWORK, 
		int *info);

	// ------------------------------------------
	// expert driver for lapack's generalized eigenvalue solver for
	// the 'any kind of matrix'
	// so A and M are both complex -> used for PML
	// ------------------------------------------
	// this is via acml	
#ifndef NOACML		
	extern void zggevx(char balanc, char jobvl, char jobvr, char sense, int n, cplx *a, int lda, 
		cplx *b, int ldb, cplx *alpha, cplx *beta, 
		cplx *vl, int ldvl, cplx *vr, int ldvr, 
		int *ilo, int *ihi, double *lscale, double *rscale, 
		double *abnrm, double *bbnrm, double *rconde, 
		double *rcondv, int *info
	);
#endif
	// this is directly fortran (STANDARD VERSION)
	extern void F77_Name(zggevx)(char* balanc, char* jobvl, char* jobvr, char* sense, int* n, cplx *a, int* lda, 
		cplx *b, int* ldb, cplx *alpha, cplx *beta, 
		cplx *vl, int* ldvl, cplx *vr, int* ldvr, 
		int *ilo, int *ihi, double *lscale, double *rscale, 
		double *abnrm, double *bbnrm, double *rconde, 
		double *rcondv,
		cplx* work, int* LWORK, double* RWORK, int* IWORK, int* BWORK,  
		int *info
	);					
	// this is directly fortran (MY MODIFIED VERSION) -> MODIFIED QZ TOLERANCE TO ACHIEVE CONVERGENCE
	extern void F77_Name(zggevv)(char* balanc, char* jobvl, char* jobvr, char* sense, int* n, cplx *a, int* lda, 
		cplx *b, int* ldb, cplx *alpha, cplx *beta, 
		cplx *vl, int* ldvl, cplx *vr, int* ldvr, 
		int *ilo, int *ihi, double *lscale, double *rscale, 
		double *abnrm, double *bbnrm, double *rconde, 
		double *rcondv,
		cplx* work, int* LWORK, double* RWORK, int* IWORK, int* BWORK,  
		int *info
	);				



	extern void F77_Name(ilaver)(int* major, int* minor, int* patch);
	extern double F77_Name(dlamch)(char* ch);
								
}



// -----------------------------------------------------
// begin of implementation of BandMatrix
// -----------------------------------------------------
template<>
void BandMatrix<cplx>::mult_vec(const cplx* in, cplx* out) const {
	if(symmetric) {
		#pragma omp parallel for default(shared) schedule(static, 5000) 
		for(int ii = 0; ii < (signed)size; ii++) {
			out[ii] = in[ii];
		}		
		char uplo  = 'U';
		char trans = 'N';
		char diag  = 'N';
		int  N     = size;
		int  K     = band_width;
		int lda    = band_width + 1;
		int inc    = 1; 
		F77_Name(ztbmv)(&uplo, &trans, &diag, &N, &K, data, &lda, out, &inc);
		
	} else {
#pragma omp parallel for default(shared) schedule(static, 5000) 
		for(int ii = 0; ii < (signed)size; ii++) {
			out[ii] = 0.0;
		}
		// TODO: replace that with own, probably parallel code
		// use ZGBMV for matrix multiplication
		// set BETA to 1 (then out is preserved)
		char trans = 'N';
		int N  = size;
		int KL = band_width;
		cplx scalar(1.0, 0.0);
		int LDA = KL * 2 + 1;
		int INC = 1;
		F77_Name(zgbmv)(&trans, &N, &N, &KL, &KL, &scalar, data, &LDA, in, &INC, &scalar, out, &INC);
		
	}	
}

template<>
void BandMatrix<double>::mult_vec(const double* in, double* out) const {
	TDKP_GENERAL_EXCEPTION("NOT YET IMPLEMENTED");	
}
// -----------------------------------------------------
// end of implementation of BandMatrix
// -----------------------------------------------------





// -----------------------------------------------------
// begin of implementation of linear equation solver 
// LapackBandSolverComplex for full matrices
// -----------------------------------------------------
LapackBandSolverComplex::LapackBandSolverComplex(unsigned int size) 
: matrix(size),
  factorized(false),
  // -----------------------------------------
  // lapack zgbsvx variables
  // -----------------------------------------
  lp_lu_A(0),
  lp_ldafb(0),
  lp_afb(0),
  lp_ipiv(0),
  lp_equed('\0'), 
  lp_r(0),
  lp_c(0),
  lp_work(0),
  lp_rwork(0)  
{	
	matrix.set_property(nonsymmetric_matrix);
}

void LapackBandSolverComplex::release_lapack_space() {
  if(lp_lu_A  != 0) { delete[] lp_lu_A;  lp_lu_A = 0; }    
  if(lp_afb   != 0) { delete[] lp_afb;   lp_afb = 0;  } 
  if(lp_ipiv  != 0) { delete[] lp_ipiv;  lp_ipiv = 0;  }    
  if(lp_r     != 0) { delete[] lp_r;     lp_r = 0;  } 
  if(lp_c     != 0) { delete[] lp_c;     lp_c = 0;  } 
  if(lp_work  != 0) { delete[] lp_work;  lp_work = 0;  } 
  if(lp_rwork != 0) { delete[] lp_rwork; lp_rwork = 0;  } 	
}

LapackBandSolverComplex::~LapackBandSolverComplex() throw(Exception*) {
	this->release_lapack_space();	
}


void LapackBandSolverComplex::prepare() throw(Exception*) {
	this->factorized = false;
}

void LapackBandSolverComplex::factorize(cplx* res, cplx* rhs) {
	
	this->release_lapack_space();
	
	int data_length;
	cplx* matrix_data = matrix.get_data(data_length);
	
	if(matrix.property_is_set(symmetric_matrix)) {
		TDKP_GENERAL_EXCEPTION("sorry, but there is NO lapack routine for hermitian band matrices!");	
	} else {
		char   fact  = 'E'; // equilibrate (if necessary) and factorize
		char   trans = 'N'; // no transpose or conjugate
		int    n     = matrix.get_size();
		int    kb    = matrix.get_band_width(); // band width
		int    nrhs  = 1; 	// number of right hand sides (we make dummy factorization now ...)
		int    ldab  = kb + kb + 1;
		int    ldb   = n;
		int    ldx   = n;
		double rcond;
		double ferr;
		double berr;
		int    info;
		// -------------------------------
		// global
		// -------------------------------
		lp_lu_A = new cplx[data_length];
		for(int ii = 0; ii < data_length; ii++) {
			lp_lu_A[ii] = matrix_data[ii];	
		}
		lp_ldafb = 3*kb + 1;
		lp_afb   = new cplx[lp_ldafb * n];
		lp_ipiv  = new int[n];
		lp_r     = new double[n];
		lp_c     = new double[n];
		lp_work  = new cplx[2*n];
		lp_rwork = new double[n];
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, "LapackBandSolverComplex: factorizing matrix using zgbsvx_");
		F77_Name(zgbsvx)(
			&fact, &trans, &n, &kb,	&kb,
			&nrhs, lp_lu_A, &ldab, lp_afb,
			&lp_ldafb, lp_ipiv, &lp_equed,	lp_r,
			lp_c, rhs, &ldb, res, &ldx,
			&rcond,	&ferr, &berr, lp_work,
			lp_rwork, &info
		);
		
		if(info == 0) {
			factorized = true;	
		} else {
			ostringstream sout;
			sout << "there was an error during the LU factorization using "
			     << "zgbsvx. INFO returned: "  << info;
			TDKP_GENERAL_EXCEPTION(sout.str()); 	
		}		 
	}	
}

void LapackBandSolverComplex::solve_equation(cplx* res, cplx* rhs, int num) throw(Exception*) {

	TDKP_ASSERT(num == 1, "sorry, only implemented to solve one equation after one");

	if(matrix.property_is_set(symmetric_matrix)) {
		TDKP_GENERAL_EXCEPTION("IMPLEMENT ME");	
	} else {
		// nonsymmetric via zgbsvx
		if(!factorized) {
			factorize(res,rhs); // solves the first iteration, but keeps the results	
		} else {
			// use factorization results of last round						
			char   fact  = 'F'; // equilibrate (if necessary) and factorize
			char   trans = 'N'; // no transpose or conjugate
			int    n     = matrix.get_size();
			int    kb    = matrix.get_band_width(); // band width
			int    nrhs  = 1; 	// number of right hand sides (we make dummy factorization now ...)
			int    ldab  = kb + kb + 1;
			int    ldb   = n;
			int    ldx   = n;
			double rcond;
			double ferr;
			double berr;
			int    info;
			F77_Name(zgbsvx)(
				&fact, &trans, &n, &kb,	&kb,
				&nrhs, lp_lu_A, &ldab, lp_afb,
				&lp_ldafb, lp_ipiv, &lp_equed,	lp_r,
				lp_c, rhs, &ldb, res, &ldx,
				&rcond,	&ferr, &berr, lp_work,
				lp_rwork, &info
			);			
		}
	}
}
// -----------------------------------------------------
// end of implementation of linear equation solver 
// LapackBandSolverComplex for full matrices
// -----------------------------------------------------	







// -----------------------------------------------------
// begin of implementation of linear equation solver 
// LapackFullSolver for full matrices
// -----------------------------------------------------	

LapackFullSolver::LapackFullSolver(unsigned int size_) 
: matrix(size_),
  factorized(false)
{
	matrix.set_property(symmetric_matrix);
}
LapackFullSolver::~LapackFullSolver() throw(Exception*) {}
void LapackFullSolver::prepare() throw(Exception*) {
	factorized  = false;
	const int n = this->matrix.get_size();
	// aquire space		
	matrix_data = this->matrix.get_data_vector();
	af.resize(n * n);
	ipiv.resize(n);			
}
void LapackFullSolver::solve_equation(cplx* res, cplx* rhs, int num) throw(Exception*) {
	TDKP_ASSERT(num == 1, "does only work for one equation at the time");
	if(!factorized) {
		this->factorize_and_solve(res,rhs);	
	} else {
		this->solve_only(res,rhs);
	}	
}			
void LapackFullSolver::factorize_and_solve(cplx* res, cplx* rhs) {
	
	char fact = 'N'; // matrix is not factorized
	char uplo = 'U'; // use the upper halve of the matrix
	int  n    = this->matrix.get_size();
	int  nrhs = 1;
	int  lda = n;	
	int  ldaf = n;
	int  ldb = n;
	int  ldx = n;
	
	double rcond;
	double ferr; // vector of length 1
	double berr; // same
 	int    info;
#ifdef NOACML
	int lwork = 8 * n;
	vector<cplx> work(lwork, 0.0);
	vector<double> rwork(n, 0.0);
	F77_Name(zhesvx)(
		&fact, &uplo, &n, &nrhs, &matrix_data[0], &lda, &af[0],
		&ldaf, &ipiv[0], rhs, &ldb, res, &ldx, &rcond, &ferr, &berr, 
		&work[0], &lwork, &rwork[0], &info
	);
#else	
	zhesvx(fact, uplo, n, nrhs, &matrix_data[0], lda, &af[0],
		   ldaf, &ipiv[0], rhs, ldb, res, ldx, &rcond, &ferr, &berr, &info);
#endif	
	TDKP_ASSERT(info == 0, "zhesvx return nonzero info");		   
	this->factorized = true;
}

void LapackFullSolver::solve_only(cplx* res, cplx* rhs) {
	char fact = 'F'; // matrix is factorized
	char uplo = 'U'; // use the upper halve of the matrix
	int  n    = this->matrix.get_size();
	int  nrhs = 1;
	int  lda = n;	
	int  ldaf = n;
	int  ldb = n;
	int  ldx = n;
	
	double rcond;
	double ferr; // vector of length 1
	double berr; // same
 	int    info;
#ifdef NOACML
	int lwork = 8 * n;
	vector<cplx> work(lwork, 0.0);
	vector<double> rwork(n, 0.0);
	F77_Name(zhesvx)(
		&fact, &uplo, &n, &nrhs, &matrix_data[0], &lda, &af[0],
		&ldaf, &ipiv[0], rhs, &ldb, res, &ldx, &rcond, &ferr, &berr, 
		&work[0], &lwork, &rwork[0], &info
	);
#else	
	zhesvx(fact, uplo, n, nrhs, &matrix_data[0], lda, &af[0],
		   ldaf, &ipiv[0], rhs, ldb, res, ldx, &rcond, &ferr, &berr, &info);
#endif	
	TDKP_ASSERT(info == 0, "zhesvx return nonzero info");		   
}
	
// -----------------------------------------------------
// end of implementation of linear equation solver 
// LapackFullSolver for full matrices
// -----------------------------------------------------	


// ------------------------------------------------------
// helper function for full spectrum solvers
// ------------------------------------------------------

/** take evals and determine ordering we would get from an iterative method */
template<class T>
void determine_ordering(const vector<T>& evals, vector<unsigned int>& ordering) {

	ordering.assign(evals.size(), 0);
	// ---------------------------------------
	// compute distance to zero
	// ---------------------------------------
	vector<double> distance(evals.size(), 0);
	for(unsigned int ii = 0; ii < evals.size(); ii++) {
		ordering[ii] = ii;
		distance[ii] = tdkp_math::abs(tdkp_math::only_real(evals[ii]));		
	}
	// ---------------------------------------
	// sort according to distance
	// ---------------------------------------
	tdkp_math::tracked_sort<unsigned int, double>(ordering.begin(), ordering.end(), distance.begin(), distance.end(), true);
			
}


	
// -----------------------------------------------------
// begin of implementation of LapackComplexBandEigenSolver
// -----------------------------------------------------	
LapackComplexBandEigenSolver::LapackComplexBandEigenSolver(unsigned int size_,unsigned int block_size_) 
: EigenSolver<cplx,double,cplx>(size_),
  matrix_stiff(size_),
  matrix_overlap(size_),
  eigenvalues(0),
  eigenvectors(0)   
{
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "LapackComplexBandEigenSolver: is used as the eigensolver (using zhbgvx)");
	matrix_stiff.set_property(symmetric_matrix);
	matrix_overlap.set_property(symmetric_matrix);	
}  

LapackComplexBandEigenSolver::~LapackComplexBandEigenSolver() {
	if(eigenvalues != 0) {
		delete[] eigenvalues;
		eigenvalues = 0;	
	}
	if(eigenvectors != 0) {
		delete[] eigenvectors;
		eigenvectors = 0;	
	}			
}

cplx LapackComplexBandEigenSolver::eigenvalue(int nn) const throw(Exception*) {
	return eigenvalues[nn]; 	
}
const cplx& LapackComplexBandEigenSolver::eigenvector(int nn, int vv) const throw(Exception*) {	
	return eigenvectors[nn * msize + vv];	  	
}


bool LapackComplexBandEigenSolver::find_eigenvectors() throw(Exception*) {
		
	// ---------------------------------------------		
	// create data copies of matrices
	// ---------------------------------------------
	vector<cplx> A_copy;
	vector<cplx> B_copy;
	int data_length;

	cplx* a_data = matrix_stiff.get_data(data_length);	
	A_copy.resize(data_length);
	for(int ii = 0; ii < data_length; ii++) {
		A_copy[ii] = a_data[ii];	
	}
	double* b_data = matrix_overlap.get_data(data_length);
	B_copy.resize(data_length);
	for(int ii = 0; ii < data_length; ii++) {
		B_copy[ii] = b_data[ii];	
	}
	vector<cplx> Q(msize * msize);
	vector<int>  ifail(msize);
	
	vector<double> evals(msize);
	vector<cplx>   evecs(msize * msize);
	
	// ---------------------------------------------
	// initialize required values
	// ---------------------------------------------
	char jobz  = 'V'; // compute eigenvalues and vectors
	char range = 'V'; // eigenvalues in the half open interval (VL, VU]
	char uplo  = 'U'; // storing upper part of matrix
	int info = 0;
	
	// ---------------------------------------------
	// check if we seek for the full spectrum
	// ---------------------------------------------
	if(Configuration::get_instance()->get("lapack_zhbgvx_calculate_full_spectrum") == 1.0) {
		range = 'A';	
	}
		
	double lower_bound, upper_bound;
	
	// ---------------------------------------------
	// get bounds
	// ---------------------------------------------
	if(get_ordering() == ascending) {
		lower_bound = Configuration::get_instance()->get("lapack_zhbgvx_cb_calc_lower_bound");
		upper_bound = Configuration::get_instance()->get("lapack_zhbgvx_cb_calc_upper_bound");
		TDKP_ASSERT(lower_bound < upper_bound, "lower_bound < upper_bound failed. please check lapack_zhbgvx_cb_calc_{lower|upper}_bound"); 
	} else {
		lower_bound = Configuration::get_instance()->get("lapack_zhbgvx_vb_calc_lower_bound");
		upper_bound = Configuration::get_instance()->get("lapack_zhbgvx_vb_calc_upper_bound");
		TDKP_ASSERT(lower_bound < upper_bound, "lower_bound < upper_bound failed. please check lapack_zhbgvx_vb_calc_{lower|upper}_bound");
	}
					
	// ---------------------------------------------
	// ask lapack to calculate
	// ---------------------------------------------
#ifdef NOACML
	double abstol = 1.0e-14;
	int    illu   = 0;
	int    bwidth = matrix_stiff.get_band_width();
	int    ldab   = bwidth + 1; 
	vector<cplx>   work(msize, 0.0);
	vector<double> rwork(msize * 7, 0.0);
	vector<int>    iwork(msize * 5, 0);
	F77_Name(zhbgvx)(&jobz, &range, &uplo, &msize, &bwidth, 
		&bwidth, &A_copy[0], &ldab,	&B_copy[0], &ldab, &Q[0], &msize,
		&lower_bound, &upper_bound, &illu, &illu, &abstol, &converged_ev, 
		&evals[0], &evecs[0], &msize, &work[0], &rwork[0], &iwork[0],
		&ifail[0], &info);
#else							
	zhbgvx(jobz, range, uplo, msize, matrix_stiff.get_band_width(), 
		matrix_overlap.get_band_width(), &A_copy[0], matrix_stiff.get_band_width() + 1,
		&B_copy[0], matrix_overlap.get_band_width() + 1, &Q[0], msize,
		lower_bound, upper_bound, 0, 0, 1.0e-14, &converged_ev, 
		&evals[0], &evecs[0], msize, &ifail[0],
		&info);
#endif		
		
	// ---------------------------------------------
	// quit on error or if we didn't find enough eigenvalues
	// ---------------------------------------------			
	if(info != 0) {
		TDKP_GENERAL_EXCEPTION("LapackComplexBandEigenSolver: lapacks zhbgvx failed with " << info);	
	}
	if(converged_ev < nev) {
		TDKP_GENERAL_EXCEPTION("LapackComplexBandEigenSolver: lapacks zhbgvx could only find " << converged_ev << " eigenvalues out of " << nev << " requested eigenvalues. i guess the bounds were too small. so please increase lapack_zhbgvx_{cb|vb}_calc_{lower|upper}_bound or modify the target value.");  	
	}	
	if(eigenvalues != 0) {
		delete[] eigenvalues; eigenvalues = 0;	
	}
	if(eigenvectors != 0) {
		delete[] eigenvectors; eigenvectors = 0; 
	}
			
	// ---------------------------------------------
	// determine the nev states with real part close
	// to 0 (as we would get it from an iterative method
	// ---------------------------------------------
	vector<unsigned int> ordering;
	determine_ordering(evals, ordering);
		
	eigenvalues  = new cplx[nev];
	eigenvectors = new cplx[nev * msize];
	for(int ii = 0; ii < nev; ii++) {
		unsigned int nn = ordering[ii];
		eigenvalues[ii] = evals[nn];
		for(int vv = 0; vv < msize; vv++) {
			eigenvectors[msize * ii + vv] = evecs[msize * nn + vv];	
		}	
	}
	converged_ev = nev; 
	this->sort_solutions(nev, eigenvectors, eigenvalues);
	return true;
}

// -----------------------------------------------------
// end of implementation of LapackComplexBandEigenSolver
// -----------------------------------------------------	
	
// -----------------------------------------------------
// begin of implementation of LapackZGGEVEigensolver
// -----------------------------------------------------	
LapackZGGEVEigenSolver::LapackZGGEVEigenSolver(unsigned int size_, unsigned int block_size_)
: EigenSolver<cplx,cplx,cplx>(size_),
  matrix_stiff(size_),
  matrix_overlap(size_),
  eigenvalues(0),
  eigenvectors(0)   
{
	matrix_stiff.set_property(nonsymmetric_matrix);
	matrix_overlap.set_property(nonsymmetric_matrix);
#ifdef NOACML	
	int major, minor, patch;
	F77_Name(ilaver)(&major, &minor, &patch);
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "LapackZGGEVEigenSolver: LAPACK VERSION " << major << "." << minor << "." << patch);
#endif			
	 
}  

LapackZGGEVEigenSolver::~LapackZGGEVEigenSolver() {
	if(eigenvalues != 0) {
		delete[] eigenvalues; eigenvalues = 0;	
	}
	if(eigenvectors != 0) {
		delete[] eigenvectors; eigenvectors = 0; 
	}
}
bool LapackZGGEVEigenSolver::find_eigenvectors() throw(Exception*) {
	
	// then create data copies of matrices
	vector<cplx> A_copy(matrix_stiff.get_data_vector());
	vector<cplx> B_copy(matrix_overlap.get_data_vector());
				
	vector<cplx> eval_alpha(msize, cplx(0.0e0,0.0e0));
	vector<cplx> eval_beta(msize, cplx(0.0e0,0.0e0));
	vector<cplx> evecs(msize * msize, cplx(0.0e0,0.0e0));
	vector<cplx> evecs_vl(msize * msize, cplx(0.0e0,0.0e0));
	
	// initialize required values
	char jobvl = 'N'; // don't compute left eigenvectors
	char jobvr = 'V'; // compute right generalized eigenvectors
	
	int info = 0;

#ifdef USE_TDKP_ZGGEV
	int lwork = 33 * msize;
	vector<cplx>   work(lwork > 1 ? lwork : 100, cplx(0.0e0,0.0e0));
	vector<double> rwork(8 * msize, 0.0e0);	

	TDKP_LOGMSG(LOG_INFO_DEVEL2, "LapackZGGEVEigenSolver: is used as the eigensolver (using MY zggev)");		
	F77_Name(zggevr)(&jobvl, &jobvr, &msize, &A_copy[0], &msize, 
		&B_copy[0], &msize, &eval_alpha[0], &eval_beta[0], 
		&evecs_vl[0], &msize, &evecs[0], &msize,
		&work[0], &lwork, &rwork[0],
		&info);
#else
	#ifdef NOACML
		int lwork = 33 * msize;
		vector<cplx>   work(lwork > 1 ? lwork : 100, cplx(0.0e0,0.0e0));
		vector<double> rwork(8 * msize, 0.0e0);			
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LapackZGGEVEigenSolver: is used as the eigensolver");
		F77_Name(zggev)(&jobvl, &jobvr, &msize, &A_copy[0], &msize, 
			&B_copy[0], &msize, &eval_alpha[0], &eval_beta[0], 
			&evecs_vl[0], &msize, &evecs[0], &msize,
			&work[0], &lwork, &rwork[0],
			&info);
	#else			
		zggev(jobvl, jobvr, msize, &A_copy[0], msize, 
			&B_copy[0], msize, &eval_alpha[0], &eval_beta[0], 
			&evecs_vl[0], msize, &evecs[0], msize, &info);						
	#endif
#endif	


	if(info != 0) {
		TDKP_GENERAL_EXCEPTION("zggev failed with " << info);	
	}	
	if(eigenvalues != 0) {
		delete[] eigenvalues; eigenvalues = 0;	
	}
	if(eigenvectors != 0) {
		delete[] eigenvectors; eigenvectors = 0; 
	}	
	eigenvalues  = new cplx[nev];
	eigenvectors = new cplx[nev * msize];
		
	for(int nn = 0; nn < msize; nn++) {
		// calculate eigenvalue given by alpha / beta
		if(abs(eval_beta[nn]) > 0) {
			eval_alpha[nn] = eval_alpha[nn] / eval_beta[nn];
		} else {
			if(abs(eval_alpha[nn]) != 0.0) {
				TDKP_GENERAL_EXCEPTION("LapackZGGEVEigenSolver: i have a problem. eigenvalue pair nr " << nn << " gives me a right eigenvalue of 0 but " << eval_alpha[nn] << " to the left, so the eigenvalue is not defined!");
			}
		}	
	}
	
	// ---------------------------------------------
	// determine the nev states with real part close
	// to 0 (as we would get it from an iterative method
	// ---------------------------------------------
	vector<unsigned int> ordering;
	determine_ordering(eval_alpha, ordering);	
	eigenvalues  = new cplx[nev];
	eigenvectors = new cplx[nev * msize];
	vector<cplx> Mx(msize);
	for(int ii = 0; ii < nev; ii++) {
		unsigned int nn = ordering[ii];
		eigenvalues[ii] = eval_alpha[nn];
		
		for(int vv = 0; vv < msize; vv++) {
			eigenvectors[msize * ii + vv] = evecs[msize * nn + vv];	
		}
		// ------------------------------
		// normalize eigenvectors
		// so that conj(zT)Mz == 1
		// ------------------------------
		matrix_overlap.mult_vec(&eigenvectors[msize*ii], &Mx[0]);
		cplx tst(0.0);
		for(int mm = 0; mm < msize; mm++) {
			tst += conj(eigenvectors[msize*ii+mm])*Mx[mm];	
		}	
		tst = sqrt(tst);
		TDKP_ASSERT(abs(tst) > 0.0, "norm of eigenvector: conj(zT)Mz equals zero! i didn't expect that");	
		for(int mm = 0; mm < msize; mm++) {
			eigenvectors[msize*ii+mm] /= tst;
		}
		
		// ----------------------------------
		// test eigenvalue
		// ----------------------------------		
	} 
	converged_ev = nev;
	this->sort_solutions(nev, eigenvectors, eigenvalues);
	return true;
	
}
cplx LapackZGGEVEigenSolver::eigenvalue(int nn) const throw(Exception*) {
	TDKP_BOUNDS_ASSERT(eigenvalues != 0, "");
	return eigenvalues[nn];
}

const cplx& LapackZGGEVEigenSolver::eigenvector(int nn, int vv) const throw(Exception*) {
	TDKP_BOUNDS_ASSERT(eigenvectors != 0, "");
	return eigenvectors[nn * msize + vv];
}
// -----------------------------------------------------
// end of implementation of LapackZGGEVEigensolver
// -----------------------------------------------------	





// -----------------------------------------------------
// begin of implementation of LapackZGGEVXEigensolver
// -----------------------------------------------------	
LapackZGGEVXEigenSolver::LapackZGGEVXEigenSolver(unsigned int size_, unsigned int block_size_)
: LapackZGGEVEigenSolver(size_, block_size_)
{		
}  


bool LapackZGGEVXEigenSolver::find_eigenvectors() throw(Exception*) {
	
	
	// create data copies of matrices
	vector<cplx> A_copy(matrix_stiff.get_data_vector());
	vector<cplx> B_copy(matrix_overlap.get_data_vector());
			
	vector<cplx> eval_alpha(msize);
	vector<cplx> eval_beta(msize);
	vector<cplx> evecs(msize * msize);
	vector<cplx> evecs_vl(msize * msize);
	
	// initialize required values
	char balanc = 'B'; // scale and permute to improve convergence
	char jobvl  = 'N'; // don't compute left eigenvectors
	char jobvr  = 'V'; // compute right generalized eigenvectors
	char sense  = 'N'; // no reciprocal condition number
	
	// output values
	vector<double> lscale(msize, 0);
	vector<double> rscale(msize, 0);
	vector<double> rconde(msize, 0);
	vector<double> rcondv(msize, 0);
	int ilo  = 0;
	int ihi  = 0;
	int info = 0;
	double abnrm = 0;
	double bbnrm = 0;
		
#ifdef USE_TDKP_ZGGEV
	int lwork = 3*msize*msize + 3*msize;
	vector<cplx> work(lwork);
	vector<double> rwork(lwork);
	vector<int>	   iwork(msize + 5);
	vector<int>    bwork(msize + 5);
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "LapackZGGEVXEigenSolver: is used as the eigensolver (using MY zggevx)");
	F77_Name(zggevv)(&balanc, &jobvl, &jobvr, &sense, &msize, &A_copy[0], &msize, 
		&B_copy[0], &msize, &eval_alpha[0], &eval_beta[0],
		&evecs_vl[0], &msize, &evecs[0], &msize, 
		&ilo, &ihi, &lscale[0], &rscale[0], 
		&abnrm, &bbnrm, &rconde[0], 
		&rcondv[0], 
		&work[0], &lwork, &rwork[0], &iwork[0], &bwork[0],   
		&info
	);
#else
	#ifdef NOACML
		int lwork = 3*msize*msize + 3*msize;
		vector<cplx> work(lwork);
		vector<double> rwork(lwork);
		vector<int>	   iwork(msize + 5);
		vector<int>    bwork(msize + 5);	 
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LapackZGGEVXEigenSolver: is used as the eigensolver (using zggevx)");		
		F77_Name(zggevx)(&balanc, &jobvl, &jobvr, &sense, &msize, &A_copy[0], &msize, 
			&B_copy[0], &msize, &eval_alpha[0], &eval_beta[0],
			&evecs_vl[0], &msize, &evecs[0], &msize, 
			&ilo, &ihi, &lscale[0], &rscale[0], 
			&abnrm, &bbnrm, &rconde[0], 
			&rcondv[0], 
			&work[0], &lwork, &rwork[0], &iwork[0], &bwork[0],   
			&info
		);
	#else
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LapackZGGEVXEigenSolver: is used as the eigensolver (acml version)");
		// not using acml ... looks like a deep error down in acml			
		zggevx(balanc, jobvl, jobvr, sense, msize, &A_copy[0], msize, 
			&B_copy[0], msize, &eval_alpha[0], &eval_beta[0],
			&evecs_vl[0], msize, &evecs[0], msize, 
			&ilo, &ihi, &lscale[0], &rscale[0], 
			&abnrm, &bbnrm, &rconde[0], 
			&rcondv[0], &info
		);						
	#endif
#endif						


	if(info != 0) {
		TDKP_GENERAL_EXCEPTION("zggevx failed with " << info);	
	}	
	if(eigenvalues != 0) {
		delete[] eigenvalues; eigenvalues = 0;	
	}
	if(eigenvectors != 0) {
		delete[] eigenvectors; eigenvectors = 0; 
	}	
	eigenvalues  = new cplx[nev];
	eigenvectors = new cplx[nev * msize];
		
	for(int nn = 0; nn < msize; nn++) {
		// calculate eigenvalue given by alpha / beta
		if(abs(eval_beta[nn]) > 0) {
			eval_alpha[nn] = eval_alpha[nn] / eval_beta[nn];
		} else {
			if(abs(eval_alpha[nn]) != 0.0) {
				TDKP_GENERAL_EXCEPTION("LapackZGGEVEigenSolver: i have a problem. eigenvalue pair nr " << nn << " gives me a right eigenvalue of 0 but " << eval_alpha[nn] << " to the left, so the eigenvalue is not defined!");
			}
		}	
	}
	
	// ---------------------------------------------
	// determine the nev states with real part close
	// to 0 (as we would get it from an iterative method
	// ---------------------------------------------
	vector<unsigned int> ordering;
	determine_ordering(eval_alpha, ordering);	
	eigenvalues  = new cplx[nev];
	eigenvectors = new cplx[nev * msize];
	vector<cplx> Mx(msize);
	for(int ii = 0; ii < nev; ii++) {
		unsigned int nn = ordering[ii];
		eigenvalues[ii] = eval_alpha[nn];
		
		for(int vv = 0; vv < msize; vv++) {
			eigenvectors[msize * ii + vv] = evecs[msize * nn + vv];	
		}
		// ------------------------------
		// normalize eigenvectors
		// so that conj(zT)Mz == 1
		// ------------------------------
		matrix_overlap.mult_vec(&eigenvectors[msize*ii], &Mx[0]);
		cplx tst(0.0);
		for(int mm = 0; mm < msize; mm++) {
			tst += conj(eigenvectors[msize*ii+mm])*Mx[mm];	
		}	
		tst = sqrt(tst);
		TDKP_ASSERT(abs(tst) > 0.0, "norm of eigenvector: conj(zT)Mz equals zero! i didn't expect that");	
		for(int mm = 0; mm < msize; mm++) {
			eigenvectors[msize*ii+mm] /= tst;
		}
		
		// ----------------------------------
		// test eigenvalue
		// ----------------------------------		
	} 
	converged_ev = nev;
	this->sort_solutions(nev, eigenvectors, eigenvalues);
	return true;
	
}



// -----------------------------------------------------
// end of implementation of LapackZGGEVXEigensolver
// -----------------------------------------------------
	
template<>
void FullMatrix<cplx>::save_to_file(const char* filename) const {

	const int n = this->get_size();
	ofstream fout(filename);
	if(fout) {	
		for(int ii = 0; ii < n; ii++) {
			for(int jj = 0; jj < n; jj++) {
				cplx val;
				if(symmetric && jj < ii) {
					val = get(jj,ii);
				} else {
					val = get(ii,jj);	
				}				
				fout << val.real() << " " << val.imag() << "   ";
			}
			fout << "\n";	
		}
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to file " << filename);	
	}	
}
	

template<>
void FullMatrix<cplx>::set(int ii, int jj, const cplx& val) {
	TDKP_BOUNDS_ASSERT(!((ii > jj) && symmetric), "!((ii(=" << ii << ") > jj(="<<jj<<") && symmetric(" << symmetric << "))");
	TDKP_BOUNDS_ASSERT(this->mdata.size() > jj * this->get_size() + ii, "");
	this->mdata[jj * this->get_size() + ii] = val;
	if(symmetric && ii != jj) {
		this->mdata[ii * this->get_size() + jj] = conj(this->mdata[jj * this->get_size() + ii]); 		
	}	
}

template<>
void FullMatrix<cplx>::add(int ii, int jj, const cplx& val) {
	TDKP_BOUNDS_ASSERT(!((ii > jj) && symmetric), "!((ii(=" << ii << ") > jj(="<<jj<<") && symmetric(" << symmetric << "))");
	TDKP_BOUNDS_ASSERT(this->mdata.size() > jj * this->get_size() + ii, "");
	this->mdata[jj * this->get_size() + ii] += val;
	if(symmetric && ii != jj) {
		this->mdata[ii * this->get_size() + jj] = conj(this->mdata[jj * this->get_size() + ii]); 		
	}
}
	
template<>
void FullMatrix<cplx>::perform_symmetry_analysis() {
	
	if(!this->symmetric) {
		Logger::get_instance()->init_progress_bar("FullMatrix: checking matrix symmetry",size);
		int next = 0;

		double tol = Configuration::get_instance()->get("assembly_check_matrix_for_symmetry_tolerance");

		double num_wrong   = 0;
		double num_correct = 0;
		double num_offdiag = 0;

		double avg_value_wrong   = 0.0;
		double avg_value_correct = 0.0;
		double avg_value_offdiag = 0.0;
		double avg_value_diag    = 0.0;
		double avg_error         = 0.0;
		cplx value_ii_kk, value_kk_ii;

		double abs_value_ii_kk, abs_value_kk_ii;
		
		for(unsigned int ii = 0; ii < this->size; ii++) {

			avg_value_diag    = (avg_value_diag * ii + tdkp_math::abs(get(ii,ii))) / (ii + 1);

			// loop over upper right triangle
			for(unsigned int kk = ii + 1; kk < size; kk++) {
				
				value_ii_kk = this->get(ii,kk);
				value_kk_ii = this->get(kk,ii);

				abs_value_ii_kk = tdkp_math::abs(value_ii_kk);
				abs_value_kk_ii = tdkp_math::abs(value_kk_ii);
				// check if correct
				if(tdkp_math::abs(conj(value_ii_kk) - value_kk_ii) < tol) {
					avg_value_correct = ((avg_value_correct * num_correct) + abs_value_ii_kk) / (num_correct + 1);
					num_correct += 1;
				} else {
					avg_value_wrong = ((avg_value_wrong * num_wrong) + abs_value_ii_kk) / (num_wrong + 1);
					avg_error       = ((avg_error       * num_wrong) + tdkp_math::abs(value_ii_kk-value_kk_ii)) / (num_wrong + 1);
					num_wrong += 1;
				}
				avg_value_offdiag = (avg_value_offdiag * num_offdiag + abs_value_ii_kk + abs_value_kk_ii) / (num_offdiag + 2);
				num_offdiag += 2;
			}
			if((signed)ii == next) {
				next = Logger::get_instance()->set_progress_bar(ii, size);
			}

		}
		double average_error_ratio = 0.0;
		if(avg_value_wrong != 0.0) {
			average_error_ratio = (avg_error / avg_value_wrong);
		}
		Logger::get_instance()->end_progress_bar();


		ostringstream sout;
		sout << "FullMatrix: the " << size << "x" << size << " matrix is " << (num_wrong > 0 ? " NOT " : "") << "symmetric (hermitian)\n"
		     << "diagonal entries:     " << setw(10) << size             << ", avg size: " << avg_value_diag << "\n"
		     << "offdiagonal entries:  " << setw(10) << int(num_offdiag) << ", avg size: " << avg_value_offdiag << "\n"
		     << "wrong pairs:          " << setw(10) << int(num_wrong)   << ", avg size: " << avg_value_wrong << "\n"
		     << "average error:        " << setw(10) << avg_error        << ", ratio:    " << average_error_ratio << "\n"
		     << "correct pairs:        " << setw(10) << int(num_correct) << ", avg size: " << avg_value_correct;

		Logger::get_instance()->emit(LOG_INFO, sout.str());

	} else {
		Logger::get_instance()->emit(LOG_INFO, "symmetry analysis: matrix is already defined symmetric. therefore not analysing anything.");
	}	
}

template<>
void FullMatrix<double>::perform_symmetry_analysis() {
	TDKP_LOGMSG(LOG_WARN, "FullMatrix: testing for symmetry not implemented for double matrix");	
}	
	
	
	
} // end of namespace
