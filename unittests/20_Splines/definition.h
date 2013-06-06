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
#include <math.h>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <boost/lexical_cast.hpp>
#include "tdkp/clc/Splines.h"
#include "tdkp/common/ICurves.h"

using namespace tdkp;

namespace tdkp {


}

class SplineTest : public CxxTest::TestSuite {

public:			
	SplineTest() {
		fout.open("20_SplineTest_output.log", ios::app | ios::out);		
		Logger::get_instance()->add_listener(&fout);
		//Logger::get_instance()->del_listener(&std::cout);
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);
        timeval time;
        gettimeofday(&time, NULL);
        srand48(time.tv_usec);				 
	}
	
	~SplineTest() {			
		Logger::get_instance()->del_listener(&fout);
		fout.close();		
	}
	std::fstream fout;
	
			
	void setUp() {}				
	void tearDown() {}		
	
	void test_clclinearcurve();
	void test_spline_1D();
	void test_nonuniform_spline_1D();
	void test_timing_1D_spline();
	
};

void SplineTest::test_spline_1D() {
	
	// test 1D spline class
	vector<double> x_lean;
	vector<double> x_dense;
	vector<double> y_lean;
	
	
	tdkp_math::linear_space(x_lean,  0.0, 6.0, 35);
	tdkp_math::linear_space(x_dense, 0.0, 6.0, 50);
	
	for(unsigned int ii = 0; ii < x_lean.size(); ii++) {
		y_lean.push_back(cos(x_lean[ii]));	
	}
	
	Spline1D spline(x_lean, y_lean, 0.0);
	
	for(unsigned int ii = 0; ii < x_dense.size(); ii++) {
		TS_ASSERT_DELTA(spline(x_dense[ii]), cos(x_dense[ii]), 0.02);	
	}
	
			
}

void SplineTest::test_nonuniform_spline_1D() {

	// test 1D spline class
	vector<double> x_lean;
	vector<double> x_dense;
	vector<double> y_lean;
	
	
	tdkp_math::linear_space(x_lean,  0.0, 6.0, 12);
	tdkp_math::linear_space(x_dense, 0.0, 6.0, 50);
			
	for(unsigned int ii = 0; ii < x_lean.size(); ii++) {
		// offset uniform points by some small random number
		if(ii > 0 && ii < x_lean.size() - 1) {
			double d = (drand48() - 0.5) * 0.2;
			x_lean[ii] += d; 	
		}
		y_lean.push_back(cos(x_lean[ii]));		
	}
		
	Spline1D spline(x_lean, y_lean, 0.0);
	
	for(unsigned int ii = 0; ii < x_dense.size(); ii++) {
		TS_ASSERT_DELTA(spline(x_dense[ii]), cos(x_dense[ii]), 0.02);
		
	}
	
}

/** thats only timing test to check wheter implementation is fast enough ... */
void SplineTest::test_timing_1D_spline() {

	// test 1D spline class
	vector<double> x;
	vector<double> y;
		
	vector<unsigned int> sizes;
	sizes.push_back(6);
	sizes.push_back(15);
	sizes.push_back(25);
	sizes.push_back(50);
	sizes.push_back(100);		

	double t1,t2,t3;
	double tspline[] = {0.0, 0.0};
	double tcurve[]  = {0.0, 0.0}; 
	
	for(unsigned int ii = 0; ii < sizes.size(); ii++) {			
		tdkp_math::linear_space(x,  0.0, 6.0, sizes[ii]);
		y.resize(sizes[ii]);		
		for(unsigned int jj = 0; jj < x.size(); jj++) {
			y[jj] = cos(x[jj]);
		}
		t1 = TimeMeasurements::tic();
		Spline1D spline(x,y,0.0);
		t2 = TimeMeasurements::tic(); 
		SplineCurve curve(x,y);
		t3 = TimeMeasurements::tic();
		
		tspline[0] = t2 - t1;
		tcurve[0]  = t3 - t2;
		cout << "size " << ii << ", spline = " << t2 - t1 << " [s], curve = " << t3 - t2 << ", rel " << (t2 - t1) / (t3 - t2) << "\n";
		unsigned int tries = 10000;
		t1 = TimeMeasurements::tic();
		double d = 0.0;
		for(unsigned int jj = 0; jj < tries; jj++) {
			 d += curve.EvalAt(6.0 * drand48());							
		}
		t2 = TimeMeasurements::tic();
		double d1 = 0.0;
		for(unsigned int jj = 0; jj < tries; jj++) {
			 d1 += spline(6.0 * drand48());							
		}
		t3 = TimeMeasurements::tic();						
		tspline[0] = t3 - t2;
		tcurve[0]  = t2 - t1;
		cout << tries << " evals, spline = " << t3 - t2 << " [s], curve = " << t2 - t1 << ", rel " << (t3 - t2) / (t2 - t1) << "\n";		
	}
}

void SplineTest::test_clclinearcurve() {

	vector<double> x;
	vector<complex<double> > yc;
	
	x.push_back(0);
	cplx a(0.2, 0.4);
	double dar = 0.02;
	double dai = 0.05;
	yc.push_back(a);
	double dx = 0.1;
	while(x.back() < 2.0) {
		x.push_back(x.back() + dx);
		yc.push_back(a + cplx(dar * x.back() * x.back(), dai * x.back() * x.back()));			
	}
	CLCLinearCurve curve(x,yc);
	for(unsigned int ii = 0; ii < x.size(); ii++) {
		TS_ASSERT_DELTA(abs(yc[ii]), abs(curve(x[ii])), 1.0e-4);		 	
	}
	for(unsigned int ii = 1; ii < x.size(); ii++) {
		cplx t = curve((x[ii] + x[ii - 1]) / 2.0);
		double minr = min(abs(yc[ii]), abs(yc[ii - 1]));
		double maxr = max(abs(yc[ii]), abs(yc[ii - 1]));
		TS_ASSERT(abs(t) >= minr);
		TS_ASSERT(abs(t) <= maxr);	
	}
	

}

#endif /*DEFINITION_H_*/
