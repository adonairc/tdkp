
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
#include "tdkp/utilities/DensityOfStates.h"
#include "tdkp/io/InputParser.h"
#include "ValueChecker.h"

using namespace tdkp;


class DensityOfStatesTest : public CxxTest::TestSuite {

public:			
	DensityOfStatesTest() { 
		fout.open("15_DensityOfStates_output.log", ios::app | ios::out);		
		Logger::get_instance()->add_listener(&fout);
		//Logger::get_instance()->del_listener(&std::cout);
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);	
	}
	
	~DensityOfStatesTest() {
		Logger::get_instance()->del_listener(&fout);
		fout.close();	
	}
			
	std::fstream fout;
	InputParser parser;
	
	void setUp() {		

	}
		
	void tearDown() {

	}		
	
	void test_linear_curves() {
		vector<double> x;
		vector<double> y;
		double dx = 0.02;
		
		for(unsigned int ii = 0; ii < 100; ii++) {
			x.push_back(double(ii) * dx);
			y.push_back(x.back() * x.back());	
		}
		
		// ------------------------------------
		// general evaluation
		// ------------------------------------
		LinearCurve a(x,y);
		vector<double> res;
		a.EvalAt(x,res);
		ValueChecker<double>::compare_vectors(y, res, 1.0e-15);
		// ------------------------------------
		// derivative
		// ------------------------------------
		a.EvalDerivativeAt(x,res);
		for(unsigned int ii = 1; ii < res.size() - 1; ii++) {
			TS_ASSERT_DELTA(res[ii], 2.0 * x[ii], 1.0e-10); 	
		}
		// ------------------------------------
		// multiplication
		// ------------------------------------
		LinearCurve two = a * 2.0;
		two.EvalAt(x,res);
		for(unsigned int ii = 0; ii < res.size(); ii++) {
			TS_ASSERT_DELTA(res[ii], 2.0 * y[ii], 1.0e-10);	
		}
		// ------------------------------------
		// division
		// ------------------------------------
		LinearCurve three = a / 2.0;
		two.EvalAt(x,res);
		for(unsigned int ii = 0; ii < res.size(); ii++) {
			TS_ASSERT_DELTA(res[ii], 2.0 * y[ii], 1.0e-10);
			res[ii] = 2.0; // init for next	
		}		
		// -----------------------------------
		// adding two curves
		// -----------------------------------
		LinearCurve one(x,res); // inits one to 2
		LinearCurve added = one * a;
		added.EvalAt(x,res);
		for(unsigned int ii = 1; ii < res.size(); ii++) {
			TS_ASSERT_DELTA(res[ii], y[ii] * 2.0, 1.0e-10);
			double xi = (x[ii] + x[ii - 1]) / 2.0;
			TS_ASSERT_DELTA(added.EvalAt(xi), a.EvalAt(xi) * 2.0, 1.0e-10);
		}
		// -----------------------------------
		// testing sqrt
		// -----------------------------------
		LinearCurve square_root = a.sqrt();
		square_root.EvalAt(x,res);
		for(unsigned int ii = 0; ii < res.size(); ii++) {
			TS_ASSERT_DELTA(res[ii], x[ii], 1.0e-10);	
		}
		
	}	
	
	void test_spline_curves() {
		vector<double> x;
		vector<double> y;
		double dx = 0.02;
		
		for(unsigned int ii = 0; ii < 500; ii++) {
			x.push_back(double(ii) * dx);
			y.push_back(x.back() * x.back());	
		}
		
		// ------------------------------------
		// general evaluation
		// ------------------------------------
		LinearCurve a(x,y);
		vector<double> res;
		a.EvalAt(x,res);
		ValueChecker<double>::compare_vectors(y, res, 1.0e-15);
		// ------------------------------------
		// derivative
		// ------------------------------------
		a.EvalDerivativeAt(x,res);
		for(unsigned int ii = 1; ii < res.size() - 1; ii++) {
			TS_ASSERT_DELTA(res[ii], 2.0 * x[ii], 1.0e-7); 	
		}											
	}		
	
	LinearCurve build_3D_analytical_dos(const vector<double> energy_x, const double& effmass);

	void test_density_of_states3D_effmass_singleband() {
		DensityOfStates dos(3);
		vector<double> k_values;
		vector<double> bands;
		build_effmass_bandstructure(0.068, 0.0, k_values, bands, 1);				
		try {
			dos.set_bandstructure(k_values.front(), k_values.back(), k_values.size(), 1, bands);
			dos.calculate();		
		//	parser.write_ascii(&dos, "dos_single_effmass_3D.dat");
						
			// ------------------------------------------------
			// single band should be equal
			// ------------------------------------------------
			try {
				ValueChecker<double>::compare_ICurves(dos.get_dos(0), dos.get_full_dos(), 2.0e-3);
			} catch(ValueCheckerError& e) {				
				cerr << e.get_message();
				TS_FAIL("this one probably failed because linear interpolation is used");					
			}			
			// analytical 3D formula
			vector<double> energy_x;
			tdkp_math::linear_space(energy_x, dos.get_full_dos().GetGridReference().front(), dos.get_full_dos().GetGridReference().back(), dos.get_full_dos().GetGridReference().size());						
			LinearCurve analytical_3d_dos = build_3D_analytical_dos(energy_x, 0.068);
			
			for(unsigned int ii = 1; ii < energy_x.size() - 2; ii++) {
				TS_ASSERT_DELTA(analytical_3d_dos.EvalAt(energy_x[ii]), dos.get_full_dos().EvalAt(energy_x[ii]), 1.0e-3);	
			}				 
									
		} catch(Exception *e) {
			TS_FAIL(e->get_reason());	
		}		
	}

	void test_density_of_states3D_effmass_single_offzero();	
	void build_effmass_bandstructure(double effmass, double deltaE, vector<double>& kvalues, vector<double>& bands, unsigned int num_subbands) const;
};

void DensityOfStatesTest::test_density_of_states3D_effmass_single_offzero() {
		DensityOfStates dos(3);
		vector<double> k_values;
		vector<double> bands;
		build_effmass_bandstructure(0.068, 0.0, k_values, bands, 1);	
		for(unsigned int ii = 0; ii < k_values.size(); ii++) {
			bands[ii] += 1.0;	
		}			
		try {
			dos.set_bandstructure(k_values.front(), k_values.back(), k_values.size(), 1, bands);
			dos.calculate();		
		//	parser.write_ascii(&dos, "dos_single_offset_effmass_3D.dat");
			
		//	parser.write_ascii(&(dos.get_dos(0)), "dos_single_band_offset_effmass_3D.dat");			
						
			// ------------------------------------------------
			// single band should be equal
			// ------------------------------------------------
			try {
				ValueChecker<double>::compare_ICurves(dos.get_dos(0), dos.get_full_dos(), 2.0e-3);
			} catch(ValueCheckerError& e) {
				cerr << e.get_message();
				TS_FAIL("this one probably failed because linear interpolation is used");				
			}			
			// analytical 3D formula
			vector<double> energy_x;
			tdkp_math::linear_space(energy_x, dos.get_full_dos().GetGridReference().front(), dos.get_full_dos().GetGridReference().back(), dos.get_full_dos().GetGridReference().size());						
			LinearCurve analytical_3d_dos = build_3D_analytical_dos(energy_x, 0.068);
														
			for(unsigned int ii = 2; ii < energy_x.size() - 3; ii++) {
				TS_ASSERT_DELTA(analytical_3d_dos.EvalAt(energy_x[ii]), dos.get_full_dos().EvalAt(energy_x[ii]), 2.0e-3);	
			}				 
									
		} catch(Exception *e) {
			TS_FAIL(e->get_reason());	
		}		
}


void  DensityOfStatesTest::build_effmass_bandstructure(double effmass, double deltaE, vector<double>& kvalues, vector<double>& bands, unsigned int num_subbands) const {
	double kmin    = 0.0;
	double kmax    = 2.0 * constants::pi / 0.5 * 0.2;		
	unsigned int numk = 64;
	double dk      = (kmax - kmin) / double(numk - 1);
	
	kvalues.resize(numk);
	bands.resize(numk * num_subbands);
	for(unsigned int ii = 0; ii < numk; ii++) {
		kvalues[ii] = kmin + double(ii) * dk;
		for(unsigned int ss = 0; ss < num_subbands; ss++) {
			bands[ss * numk + ii] = constants::hbar * constants::hbar * kvalues[ii] * kvalues[ii] / (2.0 * constants::m0 * effmass)
			                      + double(ss) * deltaE; // subband														
		}							
	} 
} 

LinearCurve DensityOfStatesTest::build_3D_analytical_dos(const vector<double> energy_x, const double& effmass) {
	
	vector<double> energy_offset = energy_x;
	for(unsigned int ii = 0; ii < energy_offset.size(); ii++) {
		energy_offset[ii] -= energy_x.front();	
	}
	
	LinearCurve tmp(energy_x, energy_offset);						
	LinearCurve analytical_3d_dos = tmp.sqrt() * 
	                              (
	                              	(constants::pi * 4.0)
	                              	/ tdkp_math::pow((2.0 * constants::pi * constants::hbar), 3.0)
	                              	* sqrt(2) * tdkp_math::pow(constants::m0 * effmass, 3.0/2.0)			                              
	                              );
	return analytical_3d_dos;	
}
