
#define TDKP_UNITTEST

// ------------------------------------------
// standard includes
// ------------------------------------------
#include <cxxtest/TestSuite.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
// ------------------------------------------
// standard tdkp includes and logger def
// ------------------------------------------
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"


// ------------------------------------------
// class includes
// ------------------------------------------
#include <math.h>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <boost/lexical_cast.hpp>
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/solvers/LinearSolver.h"
#include "tdkp/solvers/UmfpackSolver.h"
#include "tdkp/solvers/LapackBandSolver.h"
#include "MatrixCreation.h"
#ifdef LINSOLV_INCLUDE_ILS
	#include "tdkp/solvers/IlsSolver.h"
#endif
#ifdef LINSOLV_INCLUDE_SUPERLU
	#include "tdkp/solvers/SuperLUMTSolver.h"
#endif
#ifdef LINSOLV_INCLUDE_AZTECOO
	#include "tdkp/solvers/AztecOOSolver.h"
#endif
#include "tdkp/common/Vector3D.h"

using namespace tdkp;

const bool be_chatty = false;

extern "C" {
	extern void zgesv_(int* n, int* nrhs, complex<double> *a, int* lda, int *ipiv, complex<double> *b, int* ldb, int *info);
	extern void dgesv_(int* n, int* nrhs, double *a, int* lda, int *ipiv, double *b, int* ldb, int *info);
}

class LinearSolversTest : public CxxTest::TestSuite {
public:	
			
	LinearSolversTest() { 
		fout.open("13_LinearSolvers_output.log", ios::app | ios::out);		
		Logger::get_instance()->add_listener(&fout);
		if(!be_chatty) {
			Logger::get_instance()->del_listener(&std::cout);
		}
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);	
	}
	
	~LinearSolversTest() {
		Logger::get_instance()->del_listener(&fout);
		fout.close();	
	}
	
		
	std::fstream fout;	
	void setUp() {		

	}
		
	void tearDown() {

	}		

	void notest_superlu_random();
	void test_superlu_simple();

	void test_full_simple();
	void test_full_random();

	void test_rmatrix();
	void test_lapack();
	void test_aztec_simple();
	void test_aztec_random();

		
	void test_ils_simple();
	void test_ils_random();
	void test_umfpack_simple();
	void test_umfpack_random();
	void test_band_simple();

	void test_ils_simple_real();		
		
	template<class T>		
	void linear_test_random(const char* msg);
	template<class T>		
	void linear_test_simple(const char* msg);
	template<class T>		
	void linear_test_simple_real(const char* msg);	
	
															
	void solve_using_lapack(const complex<double>* matrix, const complex<double>* rhs, complex<double>* res, unsigned int size);		
	void solve_using_lapack(const double* matrix, const double* rhs, double* res, unsigned int size);
	void check_linear_solver(const char* msg, LinearSolver<cplx>& solver, complex<double>* matrix, complex<double>* rhs);
			
	void check_linear_solver(const char* msg, LinearSolver<double>& solver, double* matrix, double* rhs);
		
				
};

void LinearSolversTest::solve_using_lapack(const complex<double>* matrix, const complex<double>* rhs, complex<double>* res, unsigned int size) {

	complex<double>* copy_matrix = new complex<double>[size * size];
	
	for(unsigned int ii = 0; ii < size * size; ii++) {
		copy_matrix[ii] = matrix[ii];	
	}
	for(unsigned int ii = 0; ii < size; ii++) {
		res[ii] = rhs[ii]; 
	}

	int n    = size;
	int nrhs = 1;
	int lda  = n;
	vector<int> ipiv(n);
	int ldb  = n;
	int info;
	
	zgesv_(&n, &nrhs, copy_matrix, &lda, &ipiv[0], res, &ldb, &info);
	
	if(info != 0) {
		TDKP_GENERAL_EXCEPTION("zgesv_ returned " << info);	
	} 	
	
	delete[] copy_matrix;
				
}

void LinearSolversTest::solve_using_lapack(const double* matrix, const double* rhs, double* res, unsigned int size) {

	double* copy_matrix = new double[size * size];
	
	for(unsigned int ii = 0; ii < size * size; ii++) {
		copy_matrix[ii] = matrix[ii];	
	}
	for(unsigned int ii = 0; ii < size; ii++) {
		res[ii] = rhs[ii]; 
	}

	int n    = size;
	int nrhs = 1;
	int lda  = n;
	vector<int> ipiv(n);
	int ldb  = n;
	int info;
	
	dgesv_(&n, &nrhs, copy_matrix, &lda, &ipiv[0], res, &ldb, &info);
	
	if(info != 0) {
		TDKP_GENERAL_EXCEPTION("dgesv_ returned " << info);	
	} 	
	
	delete[] copy_matrix;
	
			
}


void LinearSolversTest::test_rmatrix() {
	
	
	// init random number generator
	timeval tp;
	gettimeofday(&tp, NULL);
	srand48(tp.tv_usec);
	const unsigned int size = 199;	
	RMatrix<cplx> matrix(size, size);
	vector<cplx>  rhs(size);
	vector<cplx>  xref(size);
	vector<cplx>  x(size, 0.0);
	// build matrix and random rhs
	for(unsigned int ii = 0; ii < size; ii++) {
		xref[ii] = 3.0 * drand48();		
		for(unsigned int jj = 0; jj < size; jj++) {
			matrix(ii,jj) = complex<double>((drand48() - 0.5), drand48() - 0.5);
		}			
	}
	// calculate Ax -> d
	for(unsigned int ii = 0; ii < size; ii++) {			
		for(unsigned int jj = 0; jj < size; jj++) {	
			rhs[ii] += matrix(ii,jj) * xref[jj]; 
		}
	}
	// solve Ax = rhs
	solve(matrix, rhs, x);
	for(unsigned int ii = 0; ii < size; ii++) {
		TS_ASSERT_DELTA(x[ii].real(), xref[ii].real(), 1.0e-3);
		TS_ASSERT_DELTA(x[ii].imag(), xref[ii].imag(), 1.0e-3);
	}

}

void LinearSolversTest::check_linear_solver(const char* msg, LinearSolver<cplx>& solver, complex<double>* matrix, complex<double>* rhs) {
		 	
	const int size = solver.get_matrix().get_size();
	cplx* res_lapack = get_random_vector<cplx>(size);
	cplx* res_solver = get_random_vector<cplx>(size);
	solver.get_matrix().reset();
		
	// set structure to matrix	
	for(int ii = 0; ii < size; ii++) {		
		int start_jj = 0;
		if(solver.get_matrix().property_is_set(symmetric_matrix)) {
			start_jj = ii;	
		}
		for(int jj = start_jj; jj < size; jj++) {
			if(abs(matrix[ii * size + jj]) > 0.0) {
				solver.get_matrix().announce(ii,jj);
			}		
		}			
	}	
	solver.get_matrix().set_structure();
	// set values to matrix
	for(int ii = 0; ii < size; ii++) {		
		int start_jj = 0;
		if(solver.get_matrix().property_is_set(symmetric_matrix)) {			
			start_jj = ii;	
		}
		for(int jj = start_jj; jj < size; jj++) {
			if(abs(matrix[ii * size + jj]) > 0.0) {
				solver.get_matrix().set(ii,jj, matrix[ii * size + jj]);
			}
		}
	}

	// compare matrix vector product for linear solver and brute force
	vector<cplx> solver_mult_vec(size, 0.0);
	vector<cplx> byhand_mult_vec(size, 0.0);
	cplx* in_mult_vec = get_random_vector<cplx>(size);
	for(int ii = 0; ii < size; ii++) {
		for(int jj = 0; jj < size; jj++) {
			byhand_mult_vec[ii] += matrix[ii * size + jj] * in_mult_vec[jj];	
		}
	}
	solver.get_matrix().mult_vec(&in_mult_vec[0], &solver_mult_vec[0]);
	// compare mult vecs
	for(int ii = 0; ii < size; ii++) {
		TSM_ASSERT_DELTA(msg, solver_mult_vec[ii].real(), byhand_mult_vec[ii].real(), 1.0e-7);
		TSM_ASSERT_DELTA(msg, solver_mult_vec[ii].imag(), byhand_mult_vec[ii].imag(), 1.0e-7);	
	}
	
	// solve with solver
	solver.prepare();		
	solver.solve_equation(res_solver, rhs);
	
				
	// solve with lapack
	reorder_for_fortran(matrix, size);			
	solve_using_lapack(matrix, rhs, res_lapack, size);
		
	for(int ii = 0; ii < size; ii++) {
		TSM_ASSERT_DELTA(msg, res_lapack[ii].real(), res_solver[ii].real(), 1.0e-7);
		TSM_ASSERT_DELTA(msg, res_lapack[ii].imag(), res_solver[ii].imag(), 1.0e-7);
	}		
		
	delete[] res_lapack;
	delete[] res_solver;		
		
}
		

void LinearSolversTest::check_linear_solver(const char* msg, LinearSolver<double>& solver, double* matrix, double* rhs) {
		 	
	const int size = solver.get_matrix().get_size();
	double* res_lapack = get_random_vector<double>(size);
	double* res_solver = get_random_vector<double>(size);
	solver.get_matrix().reset();

	// set structure to matrix	
	for(int ii = 0; ii < size; ii++) {		
		int start_jj = 0;
		if(solver.get_matrix().property_is_set(symmetric_matrix)) {
			start_jj = ii;	
		}
		for(int jj = start_jj; jj < size; jj++) {
			if(abs(matrix[ii * size + jj]) > 0.0) {
				solver.get_matrix().announce(ii,jj);
			}		
		}			
	}	
	solver.get_matrix().set_structure();
	// set values to matrix
	for(int ii = 0; ii < size; ii++) {		
		int start_jj = 0;
		if(solver.get_matrix().property_is_set(symmetric_matrix)) {
			start_jj = ii;	
		}
		for(int jj = start_jj; jj < size; jj++) {
			if(abs(matrix[ii * size + jj]) > 0.0) {
				solver.get_matrix().set(ii,jj, matrix[ii * size + jj]);
			}
		}
	}

	// compare matrix vector product for linear solver and brute force
	vector<double> solver_mult_vec(size, 0.0);
	vector<double> byhand_mult_vec(size, 0.0);
	double* in_mult_vec = get_random_vector<double>(size);
	for(int ii = 0; ii < size; ii++) {
		for(int jj = 0; jj < size; jj++) {
			byhand_mult_vec[ii] += matrix[ii * size + jj] * in_mult_vec[jj];	
		}
	}
	solver.get_matrix().mult_vec(&in_mult_vec[0], &solver_mult_vec[0]);
	// compare mult vecs
	for(int ii = 0; ii < size; ii++) {
		TSM_ASSERT_DELTA(msg, solver_mult_vec[ii], byhand_mult_vec[ii], 1.0e-7);
	}
	
	// solve with solver
	solver.prepare();		
	solver.solve_equation(res_solver, rhs);
					
	// solve with lapack
	reorder_for_fortran(matrix, size);			
	solve_using_lapack(matrix, rhs, res_lapack, size);
		
	for(int ii = 0; ii < size; ii++) {
		TSM_ASSERT_DELTA(msg, res_lapack[ii], res_solver[ii], 1.0e-7);	
	}		
		
	delete[] res_lapack;
	delete[] res_solver;		
		
}


		
 void LinearSolversTest::test_lapack() {
		
	complex<double> expect[] = {
		0.857142857142857,
		cplx(0,0.714285714285714),			
		-0.571428571428572,
		cplx(0,- 0.428571428571428), 			
		0.285714285714286,
		cplx(0,0.142857142857143)
	};

	unsigned int n = 6;
	cplx* matrix = create_laplace_hermitian(n);
	vector<complex<double> > rhs(n);
	vector<complex<double> > res(n);
	rhs[0] = 1;		
	reorder_for_fortran(matrix, n);						
	solve_using_lapack(matrix, &rhs[0], &res[0], n);
	
	for(unsigned int ii = 0; ii < n; ii++) {
		TS_ASSERT_DELTA(res[ii].real(), expect[ii].real(), 1.0e-10);
		TS_ASSERT_DELTA(res[ii].imag(), expect[ii].imag(), 1.0e-10);
	} 
	
	delete[] matrix;
			 
}

template<class MySolver> 
void LinearSolversTest::linear_test_simple(const char* msg) {
	//cerr << "linear_test_simple " << msg << "\n";
	int sizes[] = {100, 333, 500, 888}; 
	for(int ii = 0; ii < 4; ii++ ) {
		int n = sizes[ii];
		cplx* matrix = create_laplace_hermitian(n);
		MySolver solver(n);		
		for(int jj = 0; jj < 3; jj++) {
			cplx* rhs    = get_random_vector<cplx>(n);					
			check_linear_solver(msg, solver, matrix, rhs);
			delete[] rhs;
		}		
		delete[] matrix;		
	}	
}

template<class MySolver> 
void LinearSolversTest::linear_test_simple_real(const char* msg) {
	//cerr << "linear_test_simple_real " << msg << "\n";
	int sizes[] = {100, 333, 500, 888}; 
	for(int ii = 0; ii < 4; ii++ ) {
		int n = sizes[ii];
		double* matrix = create_laplace(n);
		MySolver solver(n);		
		for(int jj = 0; jj < 3; jj++) {
			double* rhs    = get_random_vector<double>(n);					
			check_linear_solver(msg, solver, matrix, rhs);
			delete[] rhs;
		}		
		delete[] matrix;		
	}	
}

template<class MySolver> 
void LinearSolversTest::linear_test_random(const char* msg) {
	//cerr << "linear_test_random " << msg << "\n";
	int sizes[] = {113, 310, 520, 897}; 
	for(int ii = 0; ii < 4; ii++ ) {		
		int n = sizes[ii];		
		cplx* matrix = create_random_hermitian(n);
		MySolver solver(n);		
		for(int jj = 0; jj < 3; jj++) {
			cplx* rhs    = get_random_vector<cplx>(n);					
			check_linear_solver(msg, solver, matrix, rhs);
			delete[] rhs;
		}		
		delete[] matrix;		
	}	
}		
		
void LinearSolversTest::test_ils_simple() {
#ifdef LINSOLV_INCLUDE_ILS
	linear_test_simple<IlsSolver<cplx> >("test_ils_simple");	
#endif
}
void LinearSolversTest::test_ils_random() {
#ifdef LINSOLV_INCLUDE_ILS
	linear_test_random<IlsSolver<cplx> >("test_ils_random");
#endif
}		
		
void LinearSolversTest::test_umfpack_simple() {
	linear_test_simple<UmfpackSolverComplex >("test_umfpack_simple");
}

void LinearSolversTest::test_umfpack_random() {
	linear_test_random<UmfpackSolverComplex >("test_umfpack_random");
}

void LinearSolversTest::test_band_simple() {
	try {
		linear_test_simple<LapackBandSolverComplex>("test_LapackBandSolverComplex");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}	
}

void LinearSolversTest::test_full_simple() {
	try {
		linear_test_simple<LapackFullSolver>("test_LapackFullSolver");		
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}	
}
void LinearSolversTest::test_full_random() {
	try {
		linear_test_random<LapackFullSolver>("test_LapackFullSolver");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}	
}

void LinearSolversTest::test_ils_simple_real() {
	try {
		linear_test_simple<LapackFullSolver>("test_LapackFullSolver");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}
}

void LinearSolversTest::test_superlu_simple() {
#ifdef LINSOLV_INCLUDE_SUPERLU
	try {
		linear_test_simple<SuperLUMTSolver>("test_SuperLUSolver");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}	
#endif	
}

/** super lu had problems calculating random matrices and actually its slow anyway ... */
void LinearSolversTest::notest_superlu_random() {	
#ifdef LINSOLV_INCLUDE_SUPERLU
	try {
		linear_test_random<SuperLUMTSolver>("test_SuperLUSolver");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}		
#endif
}

void LinearSolversTest::test_aztec_simple() {
#ifdef LINSOLV_INCLUDE_AZTECOO
	try {
		linear_test_simple<AztecOOSolver<cplx> >("test_AztecOOSolver");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}	
#endif
}

void LinearSolversTest::test_aztec_random() {
#ifdef LINSOLV_INCLUDE_AZTECOO
	return;
	try {
		linear_test_random<AztecOOSolver<cplx> >("test_AztecOOSolver");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}	
#endif
}

