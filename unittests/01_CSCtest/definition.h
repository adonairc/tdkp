
#define TDKP_UNITTEST

// ------------------------------------------
// standard includes
// ------------------------------------------
#include <cxxtest/TestSuite.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include "arlnsmat.h"

// ------------------------------------------
// standard tdkp includes and logger def
// ------------------------------------------
#include "all.h"
#include "Logger.h"
#include "Exception.h"
#define EXTERN_LOGMSG
namespace tdkp {
	Logger* logmsg = new Logger(LOG_INFO_DEVEL2);	
}
// ------------------------------------------
// class includes
// ------------------------------------------
#include <CSCMatrix.h>




using namespace tdkp;

class CSCtest : public CxxTest::TestSuite {
public:	
	std::fstream fout;	
	void setUp() {
		fout.open("01_CSCTest_output.log", ios::app | ios::out);
		logmsg->add_listener(&fout);	
	}
	void tearDown() {
		fout.close();	
	}
	void test_CSCMatrix_generation();
	void test_operation();
	void test_copy();
	void test_singular();
};

void CSCtest::test_singular() {
	int N = 99;
	CSCMatrix<double> pcsc(N);		
	for(int run = 0; run < 10; run++) {
		pcsc.clear_all();
		for(int mm = 0; mm < N; mm++) {
			for(int nn = 0; nn < N; nn++) {
				if((mm != 11 && mm != 67) && drand48() < 0.25) {
					pcsc.announce(mm, nn);
				} 
			}
		}				
		// must be singular!
		pcsc.set_structure();
		TS_ASSERT(pcsc.check_singular_structure());
	}	

}

void CSCtest::test_CSCMatrix_generation() {
	// create empty matrix
	CSCMatrix<double>*  pdcsc = NULL;
	TS_ASSERT_THROWS_ANYTHING(pdcsc = new CSCMatrix<double>(0));
	delete pdcsc;
	pdcsc = NULL;
	

}

void CSCtest::test_operation() {
	// init rnd
	timeval* tp = new timeval();
	gettimeofday(tp,NULL);
	srand48(tp->tv_usec);
	delete tp; tp = NULL;
	// ----------------------------------------------
	// create randomly populated matrices
	// ----------------------------------------------
	unsigned int sizes[] =  {5, 9, 18, 25, 33, 77, 145, 1342};
	int N;
	int num     = 6;
	CSCMatrix<double>* pcsc;
	for(int ii = 0; ii < num; ii++) {
		// create full matrix
		N = sizes[ii];
		double** cmp = new double*[N];
		for(int mm = 0; mm < N; mm++) {
			cmp[mm] = new double[N];
			for(int nn = 0; nn < N; nn++) {
				if(drand48() < 0.25) {
					cmp[mm][nn] = (drand48() - 0.5) * 25.0;	
				} else {
					cmp[mm][nn] = 0.0;	
				}	
			}
		}		
		// create sparse matrix
		TS_ASSERT_THROWS_NOTHING(pcsc = new CSCMatrix<double>(N));	
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
		TS_ASSERT_THROWS_NOTHING(pcsc->set_structure());
//		pcsc->out();
		TS_ASSERT(pcsc->check_matrix());
		
		// set and get value
		for(int mm = 0; mm < N; mm++) {
			for(int nn = 0; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {				
					try {
						pcsc->set(mm,nn,cmp[mm][nn]);
					} catch (Exception* e) {
						TS_FAIL(e->get_reason()->c_str());
						break;
					}
				}
				try {
					ETS_ASSERT_EQUALS(pcsc->get(mm,nn), cmp[mm][nn]);
				} catch (Exception* e) {
					TS_FAIL(e->get_reason()->c_str());
					break;	
				}
			}
		}				
		
		// create full matrix and compare
		double* full = pcsc->create_full_matrix();
		for(unsigned int mm = 0; mm < pcsc->get_size(); mm++) {
			for(unsigned int nn = 0; nn < pcsc->get_size(); nn++) {
				TS_ASSERT_EQUALS(pcsc->get(mm,nn), full[mm * pcsc->get_size() + nn]);	
			}	
		}
		delete[] full;
							
		for(int mm = 0; mm < N; mm++) {
			delete[] cmp[mm];
		}						
		delete pcsc;								
	}			

}

void CSCtest::test_copy() {
	// init rnd
	timeval* tp = new timeval();
	gettimeofday(tp,NULL);
	srand48(tp->tv_usec);
	delete tp; tp = NULL;

	// ----------------------------------------------
	// create randomly populated matrices
	// ----------------------------------------------
	unsigned int sizes[] = {5, 9, 18, 25, 33, 77, 145, 1342};
	int N;
	int num     = 6;
	CSCMatrix<double>* pcsc;
	CSCMatrix<double>* copy;
	for(int ii = 0; ii < num; ii++) {
		// create full matrix
		N = sizes[ii];
		double** cmp = new double*[N];
		for(int mm = 0; mm < N; mm++) {
			cmp[mm] = new double[N];
			for(int nn = 0; nn < N; nn++) {
				if(drand48() < 0.25) {
					cmp[mm][nn] = (drand48() - 0.5) * 25.0;	
				} else {
					cmp[mm][nn] = 0.0;	
				}	
			}
		}		
		// create sparse matrix
		TS_ASSERT_THROWS_NOTHING(pcsc = new CSCMatrix<double>(N));	
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
		TS_ASSERT_THROWS_NOTHING(pcsc->set_structure());
		TS_ASSERT(pcsc->check_matrix());
		// set and get value
		for(int mm = 0; mm < N; mm++) {
			for(int nn = 0; nn < N; nn++) {
				if(cmp[mm][nn] != 0.0) {
					try {	
						pcsc->set(mm,nn,cmp[mm][nn]);
					} catch (Exception* e) {
						TS_FAIL(e->get_reason()->c_str());
						break;
					}		
				}					
			}
		}				
		// copy matrix
		TS_ASSERT_THROWS_NOTHING(copy = new CSCMatrix<double>((*pcsc)));
		// compare 
		for(int mm = 0; mm < N; mm++) {
			for(int nn = 0; nn < N; nn++) {				
				TS_ASSERT_EQUALS(copy->get(mm,nn), cmp[mm][nn]);
			}
			delete[] cmp[mm];
		}						
		delete copy;
		delete pcsc;								
	}			

}
