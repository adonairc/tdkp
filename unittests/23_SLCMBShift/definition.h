#ifndef DEFINITION_H_
#define DEFINITION_H_


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
#include "tdkp/utilities/SLCMBShift.h"

using namespace tdkp;


class SLCMBShiftTest : public CxxTest::TestSuite {
public:
	SLCMBShiftTest();
	void test_basic_operation();
	void test_external_shift_file();
private:
	void create_data_file_simple(double nmin, double nmax, unsigned int num, const char* filename);
	ofstream fout;
	BandstructureDomain<cplx> cb_bands;
	BandstructureDomain<cplx> vb_bands;

};

SLCMBShiftTest::SLCMBShiftTest() {
	Logger::get_instance()->set_level(LOG_INFO_DEVEL2);
	Logger::get_instance()->del_listener(&std::cout);
	fout.open("logs/23_SLCMBShiftTest.log");
	Logger::get_instance()->add_listener(&fout);
	
	try {
	
		// --------------------------------------
		// create dummy band structure
		// --------------------------------------
		DomainMaster domain;
		create_2D_domain_radial(domain, 0.0, 2.0, 40);	
		cb_bands.reinit(1,2,1, domain);
		vb_bands.reinit(1,4,1, domain);
		const double h2m0 = constants::hbar * constants::hbar / 2.0 / constants::m0;
		for(unsigned int kk = 0; kk < domain.get_number_of_points(); kk++) {
			cb_bands.set_energy(kk, 0, 1.0 + domain.get_point(kk).get_coord_abs() * h2m0 / 0.1);
			cb_bands.set_energy(kk, 1, 1.2 + domain.get_point(kk).get_coord_abs() * h2m0 / 0.1); 	
			vb_bands.set_energy(kk, 0, 0.0 - domain.get_point(kk).get_coord_abs() * h2m0 / 0.25);
			vb_bands.set_energy(kk, 1, -0.025 - domain.get_point(kk).get_coord_abs() * h2m0 / 0.25);
			vb_bands.set_energy(kk, 2, -0.05 - domain.get_point(kk).get_coord_abs() * h2m0 / 0.25);
			vb_bands.set_energy(kk, 3, -0.1 - domain.get_point(kk).get_coord_abs() * h2m0 / 0.25);		
		}

	} catch (Exception* e) {
		cerr << e->get_reason();
		exit(1);	
	} catch (string s) {
		cerr << s;
		exit(1);
	}
		
}

void SLCMBShiftTest::create_data_file_simple(double nmin, double nmax, unsigned int num, const char* filename) {

	ofstream fout(filename);
	
	double nminlog = log10(nmin);
	double nmaxlog = log10(nmax);
	
	for(unsigned int ii = 0; ii < num; ii++) {
		for(unsigned int jj = 0; jj < num; jj++) {
			double nval = (nminlog) + (nmaxlog - nminlog) / (static_cast<double>(num) - 1.0) * static_cast<double>(ii);
			double pval = (nminlog) + (nmaxlog - nminlog) / (static_cast<double>(num) - 1.0) * static_cast<double>(jj);
			fout << pow(10.0,nval) << "   " << pow(10.0,pval) << "  " << 1.0 * (static_cast<double>(ii * jj) / static_cast<double>((num - 1) * (num - 1))) << "  1.0\n";				
		}	
	}
	fout.close();
	
}

void SLCMBShiftTest::test_basic_operation() {

	try {
		this->create_data_file_simple(1.0e-5,1.5, 20, "test_slc_shifts.dat");
		SLCMBShift mbshift("test_slc_shifts.dat");
		
		// --------------------------------------------------
		// test via density
		// --------------------------------------------------
		for(unsigned int ii = 0; ii < 10; ii++) {
			double ndens = 1.5 / 8.0 * ii;
			double shift_peak, ratio_peak;
			mbshift.calculate_mb_shift_density(shift_peak, ratio_peak, ndens, ndens);
			TS_ASSERT(shift_peak >= 0);
			TS_ASSERT(shift_peak <= 1);			  	
		}		
		
		mbshift.update_bandstructure(cb_bands,vb_bands);

		for(unsigned int ii = 0; ii < 10; ii++) {
			for(unsigned int jj = 0; jj < 10; jj++) {
				double nfermi = 1.0 + 0.4 / 10.0 * ii;
				double pfermi = 0.0 - 0.08 / 10.0 * jj;
				double shift_peak, ratio_peak;
				mbshift.calculate_mb_shift_fermi(shift_peak, ratio_peak, nfermi, pfermi, 300.0);			
				TS_ASSERT(shift_peak >= 0);
				TS_ASSERT(shift_peak <= 1);				 
			}			  	
		}
				
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	} catch (string s) {
		TS_FAIL(s);	
	}
}

void SLCMBShiftTest::test_external_shift_file() {

	try {
		SLCMBShift mbshift("23_SLCMBShift/slc_mb_shifts.dat");
		
		mbshift.update_bandstructure(cb_bands,vb_bands);
		unsigned int num = 50;	
		for(unsigned int ii = 0; ii < num; ii++) {
			for(unsigned int jj = 0; jj < num; jj++) {
				double nfermi = 1.0 + 0.4  / (num - 5) * ii;
				double pfermi = 0.0 - 0.08 / (num - 5) * jj;
				double shift_peak, ratio_peak;
				mbshift.calculate_mb_shift_fermi(shift_peak, ratio_peak, nfermi, pfermi, 300.0);
				TS_ASSERT(shift_peak > -0.2);
				TS_ASSERT(shift_peak <= 0.0); 				
				//cout << "n/p fermi = " << nfermi << "/" << pfermi << " shift peak = " << shift_peak << "  ratio = " <<ratio_peak << "\n";
			}			  	
		}
		/*
		ofstream fout("wurst.dat");		
		for(unsigned int ii = 0; ii < 100; ii++) {
			double nfermi = 1.0 + 0.4  / 100.0 * ii;
			double pfermi = 0.0 - 0.08 / 100.0 * ii;
			double shift_peak, ratio_peak;
			mbshift.calculate_mb_shift_fermi(shift_peak, ratio_peak, 0.8, pfermi, 300.0);												
			fout << pfermi << "   " << shift_peak << "\n";
		}
		fout.close();
		*/				
		
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	} catch (string s) {
		TS_FAIL(s);	
	}	
}

#endif /*DEFINITION_H_*/

