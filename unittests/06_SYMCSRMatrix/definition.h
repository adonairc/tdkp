
#define TDKP_UNITTEST

// ------------------------------------------
// standard includes
// ------------------------------------------
#include <cxxtest/TestSuite.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <map>
#include <string>
#include <boost/lexical_cast.hpp>

// ------------------------------------------
// standard tdkp includes and logger def
// ------------------------------------------
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/main/CSRSparseMatrixInterface.h"
#include "tdkp/solvers/LapackBandSolver.h"

// ------------------------------------------
// class includes
// ------------------------------------------
#include "tdkp/main/CSRMatrix.h"



using namespace tdkp;

class SYMCSRtest : public CxxTest::TestSuite {
public:
	std::fstream fout;
	void setUp() {
		fout.open("06_SYMCSRTest_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);
		Logger::get_instance()->del_listener(&std::cout);
	}
	void tearDown() {
		fout.close();
	}

	static double** create_random_symmetric_matrix(int size);
	static double** create_random_nonsymmetric_matrix(int size);
	static cplx**   create_random_hermitian_matrix(int size);
	static cplx**   create_random_nonsymmetric_cplx_matrix(int size);
	
	template<class F> static void destroy_random_matrix(F** mat, int size);
	template<class F> static void set_matrix(SparseMatrixInterface<F>* mat, F** cmp, bool symmetric);

	//void test_c2locald_announce();
	void test_full_matrix();
	void test_c2locald_interface();
	void test_read_and_save_structure();
	void test_SYMCSRMatrix_generation();
	void test_operation();
	void test_hermitian_operation();
	void test_nonsym_operation();
	void test_symmetric_detection();

	void test_copy();
	void test_copy_sym_nonsym();

	//void test_singular();
	void test_MatVec();
	void test_mat_vec_nonsymm();
	void test_hermitian_mat_vec();

	void test_mult_mat_vec();

	void test_read_write_to_file();



};
/* structure can not be singular as we force the diagonal to be set
void SYMCSRtest::test_singular() {
	int N = 99;
	SYMCSRMatrix<double> pcsc(N);
	for(int run = 0; run < 10; run++) {
		pcsc.clear_all();
		for(int mm = 0; mm < N; mm++) {
			for(int nn = mm; nn < N; nn++) {
				if((mm != 11 && mm != 67) && drand48() < 0.25) {
					pcsc.announce(mm, nn);
				}
			}
		}
		// must be singular!
		TS_ASSERT_THROWS_ANYTHING(pcsc.set_structure());
	}
}*/

// --------------------------------------------------------------
// generate matrices of size 0 (should throw)
// --------------------------------------------------------------
void SYMCSRtest::test_SYMCSRMatrix_generation() {
	// create empty matrix
	CSRMatrix<double>*  pdcsc = NULL;
	TS_ASSERT_THROWS_ANYTHING(pdcsc = new CSRMatrix<double>((unsigned)0, nonsymmetric_matrix));
	delete pdcsc; pdcsc = 0;
	TS_ASSERT_THROWS_ANYTHING(pdcsc = new CSRMatrix<double>((unsigned)0, symmetric_matrix));
	delete pdcsc;
	pdcsc = NULL;
}

// --------------------------------------------------------------
// test announce/set/get for symmetric matrices
// --------------------------------------------------------------
void SYMCSRtest::test_operation() {
	// ----------------------------------------------
	// create randomly populated matrices
	// ----------------------------------------------
	unsigned int sizes[] =  {5, 9, 18, 25, 33, 77, 145};
	int N;
	int num     = 6;
	CSRMatrix<double>* pcsc = 0;
	for(int ii = 0; ii < num; ii++) {
		// create full matrix
		N = sizes[ii];
		double** cmp = SYMCSRtest::create_random_symmetric_matrix(N);

		// create sparse matrix
		TS_ASSERT_THROWS_NOTHING(pcsc = new CSRMatrix<double>(N, symmetric_matrix));
		// announce
//		std::cout << "matrix: \n";
		for(int mm = 0; mm < N; mm++) {
			for(int nn = mm; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {
					TS_ASSERT_THROWS_NOTHING(pcsc->announce(mm,nn));
				}
//				std::cout << cmp[mm][nn] << "  ";
			}
//			std::cout << "\n";
		}
		// set structure
		try {
			pcsc->set_structure();
		} catch(Exception* e) {
			TS_FAIL((e->get_reason()));
			return;
		}
//		TS_ASSERT_THROWS_NOTHING(pcsc->set_structure());
//		pcsc->out();
		TS_ASSERT(pcsc->check_matrix());

		// set and get value
		for(int mm = 0; mm < N; mm++) {
			for(int nn = mm; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {
					try {
						pcsc->set(mm,nn,cmp[mm][nn]);
					} catch (Exception* e) {
						TS_FAIL(e->get_reason());
						break;
					}
				}
				try {
					ETS_ASSERT_EQUALS(pcsc->get(mm,nn), cmp[mm][nn]);
				} catch (Exception* e) {
					TS_FAIL(e->get_reason());
					break;
				}
			}
		}

		// create full matrix and compare
		double* full = pcsc->create_full_matrix();
		for(unsigned int mm = 0; mm < pcsc->get_size(); mm++) {
			for(unsigned int nn = 0; nn < pcsc->get_size(); nn++) {
				TS_ASSERT_EQUALS(cmp[mm][nn], full[mm * pcsc->get_size() + nn]);
			}
		}
		delete[] full;

		SYMCSRtest::destroy_random_matrix(cmp, N);

		delete pcsc;

	}
}

// --------------------------------------------------------------
// test announce/get/set for nonsymmetric matrices
// --------------------------------------------------------------
void SYMCSRtest::test_nonsym_operation() {

	// ----------------------------------------------
	// create randomly populated matrices
	// ----------------------------------------------
	unsigned int sizes[] =  {5, 9, 18, 25, 33, 77, 145};
	int N;
	int num     = 6;
	CSRMatrix<double>* pcsc = 0;
	for(int ii = 0; ii < num; ii++) {
		// create full matrix
		N = sizes[ii];
		double** cmp = SYMCSRtest::create_random_nonsymmetric_matrix(N);

		// create sparse matrix
		TS_ASSERT_THROWS_NOTHING(pcsc = new CSRMatrix<double>(N, nonsymmetric_matrix));
		// announce
//		std::cout << "matrix: \n";
		for(int mm = 0; mm < N; mm++) {
			for(int nn = 0; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {
					TS_ASSERT_THROWS_NOTHING(pcsc->announce(mm,nn));
				}
//				std::cout << cmp[mm][nn] << "  ";
			}
//			std::cout << "\n";
		}
		// set structure
		try {
			pcsc->set_structure();
		} catch(Exception* e) {
			TS_FAIL((e->get_reason()));
			return;
		}
//		TS_ASSERT_THROWS_NOTHING(pcsc->set_structure());
//		pcsc->out();
		TS_ASSERT(pcsc->check_matrix());

		// set and get value
		for(int mm = 0; mm < N; mm++) {
			for(int nn = 0; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {
					try {
						pcsc->set(mm,nn,cmp[mm][nn]);
					} catch (Exception* e) {
						TS_FAIL(e->get_reason().c_str());
						break;
					}
				}
				try {
					ETS_ASSERT_EQUALS(pcsc->get(mm,nn), cmp[mm][nn]);
				} catch (Exception* e) {
					TS_FAIL(e->get_reason());
					break;
				}
			}
		}

		// create full matrix and compare
		double* full = pcsc->create_full_matrix();
		for(unsigned int mm = 0; mm < pcsc->get_size(); mm++) {
			for(unsigned int nn = 0; nn < pcsc->get_size(); nn++) {
				TS_ASSERT_EQUALS(cmp[mm][nn], full[mm * pcsc->get_size() + nn]);
			}
		}
		delete[] full;

		SYMCSRtest::destroy_random_matrix(cmp, N);

		delete pcsc;
	}
}

// --------------------------------------------------------------
// test detection whether matrix checks if its symmetric
// --------------------------------------------------------------
void SYMCSRtest::test_symmetric_detection() {

	CSRMatrix<double>* pcsr = 0;

	int sizes[] = {5, 9, 18, 25, 33, 77, 145};
	int N;
	int num = 6;
	double** cmp;
	// for different sizes
	for(int ii = 0; ii < num; ii++) {
		N   = sizes[ii];
		// first, create symmetric matrix from symmetric data
		cmp  = SYMCSRtest::create_random_symmetric_matrix(N);
		pcsr = new CSRMatrix<double>(N, symmetric_matrix);
		SYMCSRtest::set_matrix(pcsr, cmp, true);
		TS_ASSERT(pcsr->symmetric());
		TS_ASSERT(pcsr->symmetric_by_value());
		delete pcsr;
		// next, create nonsymmetric matrix from symmetric data
		pcsr = new CSRMatrix<double>(N, nonsymmetric_matrix);
		SYMCSRtest::set_matrix(pcsr, cmp, false);
		TS_ASSERT(!pcsr->symmetric());
		TS_ASSERT(pcsr->symmetric_by_value());
		delete pcsr;
		SYMCSRtest::destroy_random_matrix(cmp,N);
	}
}


void SYMCSRtest::test_copy() {

	// ----------------------------------------------
	// create randomly populated matrices
	// ----------------------------------------------
	unsigned int sizes[] = {5, 9, 18, 25, 33, 77, 145};
	int N;
	int num     = 6;
	CSRMatrix<double>* pcsc = 0;
	CSRMatrix<double>* copy = 0;
	for(int ii = 0; ii < num; ii++) {
		// create full matrix
		N = sizes[ii];
		double** cmp = SYMCSRtest::create_random_symmetric_matrix(N);
		// create sparse matrix
		TS_ASSERT_THROWS_NOTHING(pcsc = new CSRMatrix<double>(N, symmetric_matrix));
		// announce
//		std::cout << "matrix: \n";
		for(int mm = 0; mm < N; mm++) {
			for(int nn = mm; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {
					TS_ASSERT_THROWS_NOTHING(pcsc->announce(mm,nn));
				}
//				std::cout << cmp[mm][nn] << "  ";
			}
//			std::cout << "\n";
		}
		// set structure
		try {
			pcsc->set_structure();
		} catch(Exception* e) {
			TS_FAIL((e->get_reason()));
			return;
		}
		TS_ASSERT(pcsc->check_matrix());
		// set and get value
		for(int mm = 0; mm < N; mm++) {
			for(int nn = mm; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {
					try {
						pcsc->set(mm,nn,cmp[mm][nn]);
					} catch (Exception* e) {
						TS_FAIL(e->get_reason().c_str());
						break;
					}
				}
			}
		}
		// copy matrix
		TS_ASSERT_THROWS_NOTHING(copy = new CSRMatrix<double>((*pcsc)));
		// compare
		for(int mm = 0; mm < N; mm++) {
			for(int nn = 0; nn < N; nn++) {
				TS_ASSERT_EQUALS(copy->get(mm,nn), cmp[mm][nn]);
			}
		}
		SYMCSRtest::destroy_random_matrix(cmp,N);
		delete copy;
		delete pcsc;
	}
}

double** SYMCSRtest::create_random_symmetric_matrix(int N) {
	// init rnd
	timeval* tp = new timeval();
	gettimeofday(tp,NULL);
	srand48(tp->tv_usec);
	delete tp; tp = NULL;

	double** cmp = new double*[N];
	for(int mm = 0; mm < N; mm++) {
		cmp[mm] = new double[N];
	}
	for(int mm = 0; mm < N; mm++) {
		cmp[mm][mm] = (drand48() - 0.5) * 25.0;
		for(int nn = mm + 1; nn < N; nn++) {
			if(drand48() < 0.10) {
				cmp[nn][mm] = cmp[mm][nn] = (drand48() - 0.5) * 25.0;
			} else {
				cmp[nn][mm] = cmp[mm][nn] = 0.0;
			}
		}
	}
	return cmp;
}

cplx** SYMCSRtest::create_random_hermitian_matrix(int N) {
	// init rnd
	timeval* tp = new timeval();
	gettimeofday(tp,NULL);
	srand48(tp->tv_usec);
	delete tp; tp = NULL;

	cplx** cmp = new cplx*[N];

	for(int mm = 0; mm < N; mm++) {
		cmp[mm] = new cplx[N];
	}
	for(int mm = 0; mm < N; mm++) {
		cmp[mm][mm] = (drand48() - 0.5) * 25.0;
		for(int nn = mm + 1; nn < N; nn++) {
			if(drand48() < 0.05) {
				cmp[nn][mm] = (drand48() - 0.5) * 25.0 + cplx(0.0, 1.0) * (drand48() - 0.5) * 25.0;
				cmp[mm][nn] = conj(cmp[nn][mm]);
			} else {
				cmp[nn][mm] = cmp[mm][nn] = 0.0;
			}
		}
	}
	return cmp;
}

cplx** SYMCSRtest::create_random_nonsymmetric_cplx_matrix(int N) {
	// init rnd
	timeval* tp = new timeval();
	gettimeofday(tp,NULL);
	srand48(tp->tv_usec);
	delete tp; tp = NULL;

	cplx** cmp = new cplx*[N];

	for(int mm = 0; mm < N; mm++) {
		cmp[mm] = new cplx[N];
	}
	for(int mm = 0; mm < N; mm++) {
		cmp[mm][mm] = (drand48() - 0.5) * 25.0;
		for(int nn = 0; nn < N; nn++) {
			if(drand48() < 0.05 || mm == nn) {
				cmp[nn][mm] = (drand48() - 0.5) * 25.0 + cplx(0.0, 1.0) * (drand48() - 0.5) * 25.0;
				cmp[mm][nn] = conj(cmp[nn][mm]);
			} else {
				cmp[nn][mm] = cmp[mm][nn] = 0.0;
			}
		}
	}
	return cmp;

}

double** SYMCSRtest::create_random_nonsymmetric_matrix(int N) {
	// init rnd
	timeval* tp = new timeval();
	gettimeofday(tp,NULL);
	srand48(tp->tv_usec);
	delete tp; tp = NULL;

	double** cmp = new double*[N];
	for(int mm = 0; mm < N; mm++) {
		cmp[mm] = new double[N];
	}
	for(int mm = 0; mm < N; mm++) {
		cmp[mm][mm] = (drand48() - 0.5) * 25.0;
		for(int nn = 0; nn < N; nn++) {
			if(drand48() < 0.10) {
				cmp[nn][mm] = (drand48() - 0.5) * 25.0;
			} else {
				cmp[nn][mm] = 0.0;
			}
		}
	}
	return cmp;
}

template<class F>
void SYMCSRtest::destroy_random_matrix(F** mat, int size) {
	for(int ii = 0; ii < size; ii++) {
		delete[] mat[ii];
	}
	delete[] mat;
}

// -------------------------------------------------------
// set dense matrix to sparse matrix
// -------------------------------------------------------
template<class F>
void SYMCSRtest::set_matrix(SparseMatrixInterface<F>* mat, F** cmp, bool symmetric) {
	int N = mat->get_size();
	int next = 0;
	// announce
	for(int ii = 0; ii < N; ii++) {
		next = (symmetric ? ii:0);
		for(int jj = next; jj < N; jj++) {
			if(cmp[ii][jj] != 0.0) {
				mat->announce(ii,jj);
			}
		}
	}
	// set structure
	mat->set_structure();
	// set values
	for(int ii = 0; ii < N; ii++) {
		next = (symmetric ? ii:0);
		for(int jj = next; jj < N; jj++) {
			if(cmp[ii][jj] != 0.0) {
				mat->set(ii,jj, cmp[ii][jj]);
			}
		}
	}
}


void SYMCSRtest::test_hermitian_mat_vec() {
	int sizes[] = {10, 100, 333, 677};
	cplx** cmp;
	cplx*  in;
	cplx*  out_class;
	cplx*  out_cmp;
	int N;
	for(int nn = 0; nn < 4; nn++) {
		try {
			N   = sizes[nn];
			cmp = SYMCSRtest::create_random_hermitian_matrix(N);

			CSRMatrix<cplx> mat(N, symmetric_matrix);
			in        = new cplx[N];
			out_class = new cplx[N];
			out_cmp   = new cplx[N];

			for(int ii = 0; ii < N; ii++) {
				in[ii] = (drand48() - 0.5) * 2.222 + (drand48() - 0.5) * 1.0 * cplx(0.0,1.0);
				for(int jj = ii; jj < N; jj++) {
					if(cmp[ii][jj] != 0.0) {
						mat.announce(ii,jj);
					}
				}
			}
			mat.set_structure();
			for(int ii = 0; ii < N; ii++) {
				out_cmp[ii] = 0.0;
				for(int jj = ii; jj < N; jj++) {
					// set values
					if(cmp[ii][jj] != 0.0) {
						mat.set(ii,jj, cmp[ii][jj]);
					}
				}
				for(int jj = 0; jj < N; jj++) {
					// calculate cmp
					out_cmp[ii] += cmp[ii][jj] * in[jj];
				}
			}
			/*
			std::cout << "=========================================\n";
			std::cout.precision(6);
			for(int ii = 0; ii < N; ii++) {
				for(int jj = 0; jj < N; jj++) {
					std::cout << std::setw(9) << cmp[ii][jj] << " "
					          << std::setw(9) << mat.get(ii,jj) << " | ";
				}
				std::cout << "\n";
			}*/

			mat.mult_vec(in, out_class);
			// compare
			for(int ii = 0; ii < N; ii++) {
				TS_ASSERT_EQUALS(out_class[ii], out_cmp[ii]);
			}

			delete[] in;
			delete[] out_class;
			delete[] out_cmp;

			SYMCSRtest::destroy_random_matrix(cmp, N);

		} catch(Exception* e) {
			TS_FAIL(e->get_reason());
		}
	}
}


void SYMCSRtest::test_hermitian_operation() {
	int sizes[] = {10, 100, 333, 677};
	cplx** cmp;
	int N;
	for(int nn = 0; nn < 4; nn++) {
		try {
			N   = sizes[nn];
			cmp = SYMCSRtest::create_random_hermitian_matrix(N);
			CSRMatrix<cplx> mat(N, symmetric_matrix);
			set_matrix(&mat, cmp, true);
			int wrong = 0;
			for(int ii = 0; ii < N; ii++) {
				for(int jj = 0; jj < N; jj++) {
					if(conj(mat.get(jj,ii)) != mat.get(ii,jj)) {
						wrong++;
					}
				}
			}
			if(wrong > 0) {
				TS_FAIL("hermitian symmetric mat.get(jj,ii) != mat.get(ii,jj)");
			}
			wrong = 0;
			for(int ii = 0; ii < N; ii++) {
				for(int jj = 0; jj < N; jj++) {
					if(cmp[ii][jj] != mat.get(ii,jj)) {
						wrong++;
					}
				}
			}
			if(wrong > 0) {
				TS_FAIL("hermitian symmetric cmp[ii][jj] != mat.get(ii,jj)");
			}

			// test symmetric by value
			CSRMatrix<cplx> nonmat(N, nonsymmetric_matrix);
			set_matrix(&nonmat, cmp, false);

			TS_ASSERT(nonmat.symmetric_by_value());

			destroy_random_matrix(cmp,N);


		} catch (Exception*e) {
			TS_FAIL(e->get_reason());
		}
	}

}
void SYMCSRtest::test_MatVec() {
	int sizes[] = {10, 100, 333, 677};
	double** cmp;
	double*  in;
	double*  out_class;
	double*  out_cmp;
	int N;
	for(int nn = 0; nn < 4; nn++) {
		try {
			N   = sizes[nn];
			cmp = SYMCSRtest::create_random_symmetric_matrix(N);

			CSRMatrix<double> mat(N, symmetric_matrix);
			in        = new double[N];
			out_class = new double[N];
			out_cmp   = new double[N];

			for(int ii = 0; ii < N; ii++) {
				in[ii] = (drand48() - 0.5) * 2.222;
				for(int jj = ii; jj < N; jj++) {
					if(cmp[ii][jj] != 0.0) {
						mat.announce(ii,jj);
					}
				}
			}
			mat.set_structure();
			for(int ii = 0; ii < N; ii++) {
				out_cmp[ii] = 0.0;
				for(int jj = ii; jj < N; jj++) {
					// set values
					if(cmp[ii][jj] != 0.0) {
						mat.set(ii,jj, cmp[ii][jj]);
					}
				}
				for(int jj = 0; jj < N; jj++) {
					// calculate cmp
					out_cmp[ii] += cmp[ii][jj] * in[jj];
				}
			}
			/*
			std::cout << "=========================================\n";
			std::cout.precision(6);
			for(int ii = 0; ii < N; ii++) {
				for(int jj = 0; jj < N; jj++) {
					std::cout << std::setw(9) << cmp[ii][jj] << " "
					          << std::setw(9) << mat.get(ii,jj) << " | ";
				}
				std::cout << "\n";
			}*/

			mat.mult_vec(in, out_class);
			// compare
			for(int ii = 0; ii < N; ii++) {
				TS_ASSERT_EQUALS(out_class[ii], out_cmp[ii]);
			}

			delete[] in;
			delete[] out_class;
			delete[] out_cmp;

			SYMCSRtest::destroy_random_matrix(cmp, N);

		} catch(Exception* e) {
			TS_FAIL(e->get_reason());
		}
	}
}

// ---------------------------------------------------------
// tests nonsymmetric mat vec product
// ---------------------------------------------------------
void SYMCSRtest::test_mat_vec_nonsymm() {

	int sizes[] = {10, 100, 333, 677};
	double** cmp;
	double*  in;
	double*  out_class;
	double*  out_cmp;
	int N;
	for(int nn = 0; nn < 4; nn++) {
		try {
			N   = sizes[nn];
			cmp = SYMCSRtest::create_random_nonsymmetric_matrix(N);
			CSRMatrix<double> mat(N, nonsymmetric_matrix);
			in        = new double[N];
			out_class = new double[N];
			out_cmp   = new double[N];

			SYMCSRtest::set_matrix(&mat, cmp, false);

			for(int ii = 0; ii < N; ii++) {
				in[ii] = (drand48() - 0.5) * 2.222;
			}

			for(int ii = 0; ii < N; ii++) {
				out_cmp[ii] = 0.0;
				for(int jj = 0; jj < N; jj++) {
					// calculate cmp
					out_cmp[ii] += cmp[ii][jj] * in[jj];
				}
			}

			mat.mult_vec(in, out_class);
			// compare
			for(int ii = 0; ii < N; ii++) {
				TS_ASSERT_EQUALS(out_class[ii], out_cmp[ii]);
			}

			delete[] in;
			delete[] out_class;
			delete[] out_cmp;

			SYMCSRtest::destroy_random_matrix(cmp, N);

		} catch(Exception* e) {
			TS_FAIL(e->get_reason());
		}
	}
}

// ----------------------------------------------------------
// test 'enhanced' copy constructors (which allow to alter
// matrix type (symmetric/nonsymmetric)
// ----------------------------------------------------------
void SYMCSRtest::test_copy_sym_nonsym() {

	int sizes[] = {10, 100, 333, 677};
	int num     = 4;
	int N;
	double** cmp;
	for(int ii = 0; ii < num; ii++) {
		N = sizes[ii];
		try {
			// create symmetric matrix
			cmp = SYMCSRtest::create_random_symmetric_matrix(N);

			CSRMatrix<double> mat(N, symmetric_matrix);

			SYMCSRtest::set_matrix(&mat, cmp, true);
			// create nonsymmetric out of symmetric matrix

			CSRMatrix<double> symnonsym(mat, nonsymmetric_matrix);
			// compare values
			for(int aa = 0; aa < N; aa++) {
				for(int bb = 0; bb < N; bb++) {
					TS_ASSERT_EQUALS(mat.get(aa,bb), symnonsym.get(aa,bb));
				}
			}
			SYMCSRtest::destroy_random_matrix(cmp, N);
			// create nonsymmetric dense matrix and copy it to nonsymmetric matrix
			cmp = SYMCSRtest::create_random_nonsymmetric_matrix(N);

			CSRMatrix<double> nonsym(N, nonsymmetric_matrix);

			SYMCSRtest::set_matrix(&nonsym, cmp, false);
			SYMCSRtest::destroy_random_matrix(cmp, N);
			// create nonsymmetric sparse matrix via copy

			CSRMatrix<double> nonsym_copy(nonsym, nonsymmetric_matrix);
			// compare values
			for(int aa = 0; aa < N; aa++) {
				for(int bb = 0; bb < N; bb++) {
					TS_ASSERT_EQUALS(mat.get(aa,bb), symnonsym.get(aa,bb));
				}
			}
			// check that we get an exception when creating symmetric matrix out of a nonsymmetric one
			TS_ASSERT_THROWS_ANYTHING(CSRMatrix<double> tmp(nonsym, symmetric_matrix));
		} catch(Exception* e) {
			TS_FAIL(e->get_reason());
			return;
		}
	}
}

// -------------------------------------------------------------
// test read write of matrix to file
// -------------------------------------------------------------
void SYMCSRtest::test_read_write_to_file() {

	int sizes[] = {10, 100, 333, 677};
	int num     = 4;
	int N;
	double** cmp;
	for(int ii = 0; ii < num; ii++) {
		N = sizes[ii];
		try {
			// create non symmetric matrix
			cmp = SYMCSRtest::create_random_nonsymmetric_matrix(N);
			CSRMatrix<double> mat(N, symmetric_matrix);
			SYMCSRtest::set_matrix(&mat, cmp, true);
			// store to file
			mat.save_to_file("06unittest.dat");
			// init from file
			CSRMatrix<double> read("06unittest.dat");
			remove("06unittest.dat");
			// compare
			for(int aa = 0; aa < N; aa++) {
				for(int bb = 0; bb < N; bb++) {
					TS_ASSERT_EQUALS(mat.get(aa,bb),read.get(aa,bb));
				}
			}
			SYMCSRtest::destroy_random_matrix(cmp, N);
		} catch(Exception* e) {
			TS_FAIL(e->get_reason());
			return;
		}
	}
}

void SYMCSRtest::test_read_and_save_structure() {

	try {
		int N = 1000;
		// create non symmetric matrix
		double** cmp = SYMCSRtest::create_random_nonsymmetric_matrix(N);
		CSRMatrix<double> mat(N, nonsymmetric_matrix);
		CSRMatrix<double> mat2(N, nonsymmetric_matrix);
		SYMCSRtest::set_matrix(&mat, cmp, true);

		mat.save_structure("test_structure");
		mat2.load_structure("test_structure");

		for(unsigned int ii = 0; ii < mat.get_size() + 1; ii++) {
			TS_ASSERT(mat.get_prow()[ii] == mat2.get_prow()[ii]);
		}
		for(unsigned int ii = 0; ii < mat.get_num_nonzeros(); ii++) {
			TS_ASSERT(mat.get_icol()[ii] == mat2.get_icol()[ii]);
		}
		SYMCSRtest::destroy_random_matrix(cmp, N);

	} catch(Exception* e) {
		TS_FAIL(e->get_reason());
		return;
	}

}

void SYMCSRtest::test_mult_mat_vec() {
	int sizes[] = {10, 100, 333, 677, 6000};

	double** cmp;
	double*  in_single;
	double*  out_single;
	double*  in_mult;
	double*  out_mult;
	std::ostringstream sout;

	int N;
	for(int nn = 0; nn < 4; nn++) {
		try {
			N   = sizes[nn];

			in_single  = new double[N];
			out_single = new double[N];

			// ----------------------------------------------------
			// test nonsymmetric
			// ----------------------------------------------------
			cmp = SYMCSRtest::create_random_nonsymmetric_matrix(N);
			CSRMatrix<double> nonsymmat(N, nonsymmetric_matrix);
			SYMCSRtest::set_matrix(&nonsymmat, cmp, false);
			SYMCSRtest::destroy_random_matrix(cmp, N);
			cmp = SYMCSRtest::create_random_symmetric_matrix(N);
			CSRMatrix<double> symmat(N, symmetric_matrix);
			SYMCSRtest::set_matrix(&symmat, cmp, true);
			SYMCSRtest::destroy_random_matrix(cmp, N);

			// try for various multiple sizes
			for(int kpn = 2; kpn < 5; kpn++) {
				in_mult  = new double[N * kpn];
				out_mult = new double[N * kpn];
				// populate mult in
				for(int ii = 0; ii < kpn * N; ii++) {
					in_mult[ii] = (drand48() - 0.5) * 2.0;
					out_mult[ii] = 0.0;
				}
				// multiply symmetric
				symmat.mult_vec_multiple(in_mult, out_mult, kpn);
				// compare to single results
				for(int aa = 0; aa < kpn; aa++) {
					// copy in
					for(int ii = 0; ii < N; ii++) {
						in_single[ii] = in_mult[ii * kpn + aa];
					}
					// multiply single
					symmat.mult_vec(in_single, out_single);
					// compare to multi result
					int wrong = 0;
					for(int ii = 0; ii < N; ii++) {
						if(out_single[ii] != out_mult[ii * kpn + aa]) {
							wrong++;
						}
					}
					if(wrong > 0) {
						sout.str(""); sout << "symmetric mult_vec_multiple failed for N: " << N << " / kpn: " << kpn;
						TS_FAIL(sout.str().c_str());
					}

				}
				// multiply nonsymmetric
				nonsymmat.mult_vec_multiple(in_mult, out_mult, kpn);
				// compare to single results
				for(int aa = 0; aa < kpn; aa++) {
					// copy in
					for(int ii = 0; ii < N; ii++) {
						in_single[ii] = in_mult[ii * kpn + aa];
					}
					// multiply single
					nonsymmat.mult_vec(in_single, out_single);
					// compare to multi result
					int wrong = 0;
					for(int ii = 0; ii < N; ii++) {
						if(out_single[ii] != out_mult[ii * kpn + aa]) {
							wrong++;
						}
					}
					if(wrong > 0) {
						sout.str(""); sout << "nonsymmetric mult_vec_multiple failed for N: " << N << " / kpn: " << kpn;
						TS_FAIL(sout.str().c_str());
					}

				}
				delete[] in_mult;
				delete[] out_mult;
			}

			delete[] in_single;
			delete[] out_single;

		} catch(Exception* e) {
			TS_FAIL(e->get_reason());
			return;
		}
	}
}

void SYMCSRtest::test_c2locald_interface() {
	
	try {
		int N = 50;
		cplx** cmp = SYMCSRtest::create_random_hermitian_matrix(N);
		CSRMatrix<double> mat(N * 2, nonsymmetric_matrix);
		C2LocalDCSRSparseMatrixInterface interface(mat);
		vector<cplx> in(N);
		vector<cplx> out_cmp(N);
		vector<cplx> out(N);
		// assure that cmp has nonzeros on diagonal
		for(int ii = 0; ii < N; ii++) {
			cmp[ii][ii] = (drand48() - 0.5) * 2.222;
			in[ii] = cplx((drand48() - 0.5) * 2.222, (drand48() - 0.5));
		}		
		SYMCSRtest::set_matrix(&interface, cmp, false);
		// matrix vector product
		for(int ii = 0; ii < N; ii++) {
			out_cmp[ii] = 0.0;
			for(int jj = 0; jj < N; jj++) {
				out_cmp[ii] += cmp[ii][jj] * in[jj];
			}
		}
		interface.mult_vec(&in[0], &out[0]);
		for(int ii = 0; ii < N; ii++) {
			TS_ASSERT_DELTA(tdkp_math::abs(out_cmp[ii] - out[ii]), 0.0, 1.0e-12);	
		}
		SYMCSRtest::destroy_random_matrix(cmp, N);
		
	} catch(Exception* e) {
		TS_FAIL(e->get_reason());
		return;
	}
	
}

/** test FullMatrix object */
void SYMCSRtest::test_full_matrix() {
	try { 
		cplx i(0.0, 1.0);
		cplx wurst[5][5] = {
			{1.0, i, 0, 0, 0},
			{i, 2.0, i, 0, 0},
			{0,   i, 3.0, -i, 0},
			{0,  0,    i, 4.0, -i},
			{0, 0, 0, i, 5.0}			 	
		};
		cplx tv[] = {1.0, 1.0/3.0, -i, i/4.0, -3.0};
		
		FullMatrix<cplx> M(5);
		M.set_property(nonsymmetric_matrix);
		M.set_structure();
		for(unsigned int ii = 0; ii < 5; ii++) {
			for(unsigned int jj = 0; jj < 5; jj++) {
				M.set(ii,jj,wurst[ii][jj]);
			}
		} 		
		for(unsigned int ii = 0; ii < 5; ii++) {
			for(unsigned int jj = 0; jj < 5; jj++) {
				TS_ASSERT_DELTA(M.get(ii,jj).real(), wurst[ii][jj].real(), 1.0e-12);
				TS_ASSERT_DELTA(M.get(ii,jj).imag(), wurst[ii][jj].imag(), 1.0e-12);
			} 
		}
		cplx ref[5];
		cplx mv[5];
		for(unsigned int ii = 0; ii < 5; ii++) {
			ref[ii] = 0.0;
			for(unsigned int jj = 0; jj < 5; jj++) {
				ref[ii] += wurst[ii][jj] * tv[jj];	
			}
		}		
		M.mult_vec(tv,mv);	
		for(unsigned int ii = 0; ii < 5; ii++) {
			TS_ASSERT_DELTA(mv[ii].imag(),ref[ii].imag(),1.0e-12);			
		}
		
	} catch(Exception* e) {
		TS_FAIL(e->get_reason());
		return;
	}		
}
