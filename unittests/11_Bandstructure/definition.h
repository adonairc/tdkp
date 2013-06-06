#ifndef DEFINITION_H_
#define DEFINITION_H_

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
#include <time.h>

// ------------------------------------------
// standard tdkp includes and logger def
// ------------------------------------------
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"


// ------------------------------------------
// class includes
// ------------------------------------------
#include "tdkp/common/DataTypes.h"
#include "tdkp/main/EigenSolution.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/main/Bandstructure.h"

using namespace tdkp;

class BandstructureTest : public CxxTest::TestSuite {
public:	

	BandstructureTest() {
		timeval time;
		gettimeofday(&time, NULL);
		srand48(time.tv_usec);
	}

	std::fstream fout;	
	void setUp() {
		fout.open("11_BandstructureTest_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);	
		Logger::get_instance()->del_listener(&std::cout);		
	}
	void tearDown() {
		fout.close();	
	}
	
	EigenSolution<double>* get_random_solution(int length, int num_per_node);

	
	void test_eigensolution_basics();
	void test_bandstructure_basics();
	void test_bandstructure_copy();
	void test_bandstructure_domain();

	
	
};


void BandstructureTest::test_eigensolution_basics() {
	
	int sizes[] = {10, 100, 645, 777};
	
	for(int ii = 0; ii < 4; ii++) {
		for(int ee = 1; ee < 4; ee++) {
			try {
				// test eigensolution
				EigenSolution<double>* test = get_random_solution(sizes[ii], ee);
				test->set_energy(1.0);
				TS_ASSERT_EQUALS(test->get_energy(), 1.0);
				TS_ASSERT_EQUALS(test->get_length(), sizes[ii]);				
				// test copy constructur and comparison
				{ 
					EigenSolution<double> two(*test);
					TS_ASSERT(two.compare(*test));
					test->set_energy(2.0);
					TS_ASSERT(!two.compare(*test));
					test->set_energy(two.get_energy());
					TS_ASSERT(two.compare(*test));
					test->set_node_value(0,0, -12.0);
					TS_ASSERT(!two.compare(*test));
				}
				// test assignment constructor
				{
					EigenSolution<double> two;
					two = *test;
					TS_ASSERT(two.compare(*test));
					test->set_energy(2.0);
					TS_ASSERT(!two.compare(*test));
					test->set_energy(two.get_energy());
					TS_ASSERT(two.compare(*test));					
					test->set_node_value(0,0, 12345.0);
					TS_ASSERT_EQUALS(test->get_node_value(0,0), 12345.0);					
					TS_ASSERT(!two.compare(*test));
				}
				delete test;
			} catch (Exception*e) {
				TS_FAIL(e->get_reason());	
			}
		}
	}				
}

EigenSolution<double>* BandstructureTest::get_random_solution(int length, int num_per_node) {

	vector<double> data;
	data.resize(num_per_node * length);
	for(int ii = 0; ii < (signed)data.size(); ii++) {
		data[ii] = drand48(); 	
	}
	return new EigenSolution<double>(drand48(), &data[0], length, num_per_node);

}

void BandstructureTest::test_bandstructure_domain() {

	DomainMaster domain(new DomainNodeLine(0.0, 2.0));
	domain.refine(); domain.refine();
	domain.update();
	cout << "domain: " << domain.get_number_of_points() << "\n";
	int num_bands = 3;
	int sol_length = 100;
	
	try {
	
		BandstructureDomain<double> bands(num_bands,5, sol_length, domain);	
		TS_ASSERT((signed)domain.get_number_of_points() == bands.get_number_of_k_values());
		
	} catch(Exception* e) {
		TS_FAIL(e->get_reason());	
	}
	
	
}

void BandstructureTest::test_bandstructure_basics() {
	
	try {

		Vector3D dir(1.0, 0.0, 0.0);
		int basis_size      = 6;
		int num_bands       = 8;
		int solution_length = 100; 
		double kmin         = 0.0;
		double kmax         = 2.0 * constants::pi / 0.5;
		int num_k_values    = 10; 
		
		DomainMaster domain;
		create_2D_domain_radial(domain, kmin, kmax, num_k_values);
		BandstructureDomain<double> band(basis_size, num_bands, solution_length, domain); 		
		for(int bb = 0; bb < num_bands; bb++) {
			for(int kk = 0; kk < num_k_values; kk++) {
				TS_ASSERT_THROWS_NOTHING(band.add_eigensolution(kk, bb, get_random_solution(solution_length, basis_size)));	
			}	
		}
	
		// compare with input
		TS_ASSERT(band.get_number_of_bands() == num_bands);
		TS_ASSERT(band.get_basis_size()      == basis_size); 
		TS_ASSERT(band.get_number_of_k_values() == num_k_values); 
		TS_ASSERT(band.get_solution_length() == solution_length);
		TS_ASSERT(band.get_solution_length() == band.get_eigensolution(0,0).get_length()); 
		 
		TS_ASSERT_THROWS_NOTHING(band.write_binary("bandstructure.dat")); 
		BandstructureDomain<double> band_read("bandstructure.dat");
		
		TS_ASSERT(band_read.get_number_of_bands() == num_bands);
		TS_ASSERT(band_read.get_basis_size()      == basis_size); 
		TS_ASSERT(band_read.get_number_of_k_values() == num_k_values); 
		TS_ASSERT(band_read.get_solution_length() == solution_length);
		TS_ASSERT(band_read.get_solution_length() == band.get_eigensolution(0,0).get_length()); 
		 	
			
		// check dispersion relations (get_energy functions)
		for(int bb = 0; bb < num_bands; bb++) {
			for(int kk = 0; kk < num_k_values; kk++) {
				TS_ASSERT(band_read.get_energy(kk,bb) == band.get_energy(kk,bb)); 
			}
		}		
		
		remove("bandstructure.dat");
		
	} catch(Exception* e) {
		TS_FAIL(e->get_reason()); 	
	}
}

void BandstructureTest::test_bandstructure_copy() {
	
	try {
		Vector3D dir(1.0, 0.0, 0.0);
		int basis_size      = 6;
		int num_bands       = 8;
		int solution_length = 100; 
		double kmin         = 0.0;
		double kmax         = 2.0 * constants::pi / 0.5;
		int num_k_values    = 10; 
		
		DomainMaster domain;
		create_2D_domain_radial(domain, kmin, kmax, num_k_values);		
		BandstructureDomain<double> band(basis_size, num_bands, solution_length, domain); 		
		for(int bb = 0; bb < num_bands; bb++) {
			for(int kk = 0; kk < num_k_values; kk++) {
				TS_ASSERT_THROWS_NOTHING(band.add_eigensolution(kk, bb, get_random_solution(solution_length, basis_size)));	
			}	
		}
	
		BandstructureDomain<double> band2(band);
		
		TS_ASSERT(band == band2);

	} catch(Exception* e) {
		TS_FAIL(e->get_reason()); 	
	}		
}


#endif /*DEFINITION_H_*/
