
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
#include "tdkp/common/Configuration.h"
#include "tdkp/main/CSRMatrix.h"
#include "tdkp/solvers/EigenSolver.h"
#include "tdkp/solvers/ArpackSIPardisoSolver.h"
//#include "tdkp/solvers/JacobiDavidsonSolver.h"
#include "tdkp/solvers/LapackBandSolver.h"
#include "tdkp/solvers/JDQZSolver.h"
#include "tdkp/main/EigenSolution.h"
#include "tdkp/common/Vector3D.h"

using namespace tdkp;

const double matlab_values_AI[50] = {
   0.77700306501495,
   0.57979516721223,
   0.23543188324653,
   0.06674650102938,
   0.01461347063888,
   0.00261186830891,
   0.00039483184429,
   0.00005174398226,
   0.00000598709009,
   0.00000062027023,
   0.00000005818406,
   0.00000000498694,
   0.00000000039351,
   0.00000000002877,
   0.00000000000196,
   0.00000000000012,
   0.00000000000001,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
   0.00000000000000,
  -0.00000000000000,
  -0.00000000000000	
};

const double matlab_values_AB[50] = {
   0.52911324245681,
   0.39663626758081,
   0.16191355244630,
   0.04616673145654,
   0.01016853268454,
   0.00182872141762,
   0.00027820525945,
   0.00003669619726,
   0.00000427391748,
   0.00000044573232,
   0.00000004209284,
   0.00000000363223,
   0.00000000028857,
   0.00000000002124,
   0.00000000000146,
   0.00000000000009,
   0.00000000000001,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000,
   0.00000000000000   
};

class EigenSolverTest : public CxxTest::TestSuite {
public:	

	
	static const bool write_dat = false;
	const int length;
		
	fstream fout;		
		
	EigenSolverTest() : length(50) { 
		fout.open("07_EigenSolvers_output.log", ios::out);		
		Logger::get_instance()->add_listener(&fout);
		Logger::get_instance()->del_listener(&std::cout);
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);	
	}		
	
	~EigenSolverTest() {
		fout.close();	
	}
			
	template<class A, class M>
	void init_AB(EigenSolver<A,M,cplx>& solver) {
		SparseMatrixInterface<A>*   A = &solver.get_stiff();
		SparseMatrixInterface<M>* B = &solver.get_overlap();	
		bool symmetric = A->property_is_set(symmetric_matrix);
		A->reset();
		B->reset();
		for(int ii = 0; ii < length - 1; ii++) {
			A->announce(ii,ii+1);
			B->announce(ii,ii+1);
			if(!symmetric) {				
				A->announce(ii+1,ii);				
				B->announce(ii+1,ii);
			}				
		}
		for(int ii = 0; ii < length; ii++) {
			A->announce(ii,ii);
			B->announce(ii,ii);	
		}
		A->set_structure();
		B->set_structure();		
		cplx i(0,1.0);
		for(int ii = 0; ii < length - 1; ii++) {
			A->set(ii,ii+1, i);
			B->set(ii,ii+1,1.0);
			if(!symmetric) {
				A->set(ii+1,ii,-i);
				B->set(ii+1,ii,1.0);
			}
							
		}		
		for(int ii = 0; ii < length; ii++) {
			A->set(ii,ii,(double)ii + 1.0);
			B->set(ii,ii, 2.0);	
		}		
		
	} 
	 
	template<class A, class M>		
	void init_Aeye(EigenSolver<A,M,cplx>& solver) {
				
		SparseMatrixInterface<A>* A = &solver.get_stiff();
		SparseMatrixInterface<M>* B = &solver.get_overlap();	
		bool symmetric = A->property_is_set(symmetric_matrix);
		A->reset();
		B->reset();
						
		for(int ii = 0; ii < length - 1; ii++) {
			A->announce(ii,ii+1);
			if(!symmetric) {
				A->announce(ii+1,ii);
			}					
		}
		for(int ii = 0; ii < length; ii++) {
			A->announce(ii,ii);
			B->announce(ii,ii);	
		}
		A->set_structure();
		B->set_structure();
		cplx i(0,1.0);
		for(int ii = 0; ii < length - 1; ii++) {
			A->set(ii,ii+1, i);
			if(!symmetric) {
				A->set(ii+1,ii,-i);
			}						
		}		
		for(int ii = 0; ii < length; ii++) {
			A->set(ii,ii,(double)ii + 1.0);
			B->set(ii,ii, 1.0);	
		}		
		
	}
	void tearDown() {

	}

	void test_jdqz_double() {
		JDQZSolver<double> solver(length, 1);
		solver.set_ordering(ascending);		
		this->evaluate_Aeye(solver);
		this->evaluate_AB(solver);							
	}

	void test_jdqz_cplx() {
		JDQZSolver<cplx> solver(length, 1);
		solver.set_ordering(ascending);		
		this->evaluate_Aeye(solver);
		this->evaluate_AB(solver);							
	}	

	void test_lapack_zggevx() {
		LapackZGGEVXEigenSolver solver(length, 1);
		solver.set_ordering(ascending);		
		this->evaluate_Aeye(solver);
		this->evaluate_AB(solver);
					
	}	
	
	void test_RMatrix_complex();
	void test_lapack_band() {
		LapackComplexBandEigenSolver solver(length, 1);
		Configuration::get_instance()->set("lapack_zhbgvx_calculate_full_spectrum", 1.0);
		solver.set_ordering(ascending);		
		this->evaluate_Aeye(solver);	
		this->evaluate_AB(solver);			
	}


	void test_lapack_zggev() {
		LapackZGGEVEigenSolver solver(length, 1);
		solver.set_ordering(ascending);		
		this->evaluate_Aeye(solver);
		this->evaluate_AB(solver);			
	}	
			
	void test_arpack() {
		ArpackSIPardisoSolver solver(length, 1);
		solver.set_ordering(ascending);		
		this->evaluate_Aeye(solver);	
		this->evaluate_AB(solver);			
	}
	
	void no_test_jacobi_davidson() {
		//JacobiDavidsonSolver solver(length, 1);
		//solver.set_ordering(ascending);		
		//this->evaluate_Aeye(solver);	
		//this->evaluate_AB(solver);
	}
	
	/** evaluate with A = [1 idx 1] B = I  */
	template<class A, class M>
	void evaluate_Aeye(EigenSolver<A,M,cplx>& solver) {
		try {
			init_Aeye(solver);
			solver.assign(5, indefinite);
			solver.find_eigenvectors();
			TS_ASSERT(solver.converged_eigenvalues() == 5);			
			// compare
			// attention: the calculated eigenvectors are complex and therefore
			// are equal up to a phase factor. therefore we compare here absolute values
			if(write_dat) {				
				// this one is to write the solution data
				double* tmp = new double[solver.get_stiff().get_size()];
				for(unsigned int ii = 0; ii < solver.get_stiff().get_size(); ii++) {			
					tmp[ii] = abs(solver.eigenvector(0,ii));
				}
				EigenSolution<double> compare(solver.eigenvalue(0).real(), tmp, solver.get_stiff().get_size());
				fstream fout("07_EigenSolvers/data/SolutionAI_0.dat", ios::out);
				compare.write_binary(fout);
				fout.close();							
			} else {
				EigenSolution<double> compare;
				fstream fin("07_EigenSolvers/data/SolutionAI_0.dat", ios::in);
				compare.read_binary(fin); fin.close();
				TS_ASSERT_DELTA(solver.eigenvalue(0).real(), compare.get_energy(), 1.0e-8);
				// matlab result
				TS_ASSERT_DELTA(solver.eigenvalue(0).real(), 0.25380581709664, 1.0e-8);									
				for(unsigned int ii = 0; ii < solver.get_stiff().get_size(); ii++) {   					
					TS_ASSERT_DELTA(abs(solver.eigenvector(0,ii)), compare(ii), 1.0e-8);	
					TS_ASSERT_DELTA(abs(solver.eigenvector(0,ii)),matlab_values_AI[ii], 1.0e-8);
				}		
			}	
		} catch(Exception* e) {
			TS_FAIL(e->get_reason());	
		}
	}
	
	/** evaluate with B = [1 2 1] and A = [1 idx 1] */
	template<class A, class M>
	void evaluate_AB(EigenSolver<A,M,cplx>& solver) {
		try {
			init_AB(solver);
			solver.assign(5, positive_definite);
			solver.find_eigenvectors();
			TS_ASSERT(solver.converged_eigenvalues() == 5);			
			// compare
			// attention: the calculated eigenvectors are complex and therefore
			// are equal up to a phase factor. therefore we compare here absolute values
			if(write_dat) {				
				// this one is to write the solution data
				double* tmp = new double[solver.get_stiff().get_size()];
				for(unsigned int ii = 0; ii < solver.get_stiff().get_size(); ii++) {			
					tmp[ii] = abs(solver.eigenvector(0,ii));
				}
				EigenSolution<double> compare(solver.eigenvalue(0).real(), tmp, solver.get_stiff().get_size());
				fstream fin("07_EigenSolvers/data/SolutionAB_0.dat", ios::out);
				compare.write_binary(fin); fin.close();							
			} else {
				EigenSolution<double> compare;
				fstream fin("07_EigenSolvers/data/SolutionAB_0.dat", ios::in);
				compare.read_binary(fin); fin.close();			
				for(unsigned int ii = 0; ii < solver.get_stiff().get_size(); ii++) {			
					TS_ASSERT_DELTA(abs(solver.eigenvector(0,ii)), compare(ii), 1.0e-8);					
					TS_ASSERT_DELTA(abs(solver.eigenvector(0,ii)),matlab_values_AB[ii], 1.0e-8);					
				}
				TS_ASSERT_DELTA(solver.eigenvalue(0).real(), compare.get_energy(), 1.0e-8);
				TS_ASSERT_DELTA(solver.eigenvalue(0).real(), 0.12239091129195, 1.0e-8);
			}	
		} catch(Exception* e) {
			TS_FAIL(e->get_reason());	
		}
	}
	
		
};

void  EigenSolverTest::test_RMatrix_complex() {
	
	unsigned int size = 10;
	
	RMatrix<cplx> matrix(size,size);
	
	for(unsigned int ii = 0; ii < size; ii++) {
		matrix(ii,ii) = cplx(ii+1,0.0);	
	}
		
	vector<cplx> eigenvalues;
	vector<cplx> eigenvectors;
	
	RMatrix<cplx>::get_eigensystem(matrix, eigenvalues, eigenvectors);
	
	for(unsigned int ii = 0; ii < 10; ii++) {
		TS_ASSERT_DELTA(eigenvalues[ii].real(), double(ii + 1), 1.0e-12);	
	}
			
		
}
		
