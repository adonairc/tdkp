#ifndef DEFINITION_H_
#define DEFINITION_H_

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
#include "tdkp/coulomb/CoulombIntegrator.h"
#include "tdkp/coulomb/CoulombFunction.h"
#include "tdkp/coulomb/CoulombMatrixElement.h"
#include "tdkp/clc/CLCCoulombLambda.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/io/InputParser.h"
#include "tdkp/io/AsciiGridReader.h"

using namespace tdkp;

const bool be_chatty = false;

namespace tdkp {

/** coulomb function MOC returning 1 */
class CoulombFunctionMOC : public CoulombFunction {
public:
	CoulombFunctionMOC(double value_) : CoulombFunction(1), value(value_) {}
	virtual void evaluate(const vector<double>& z, vector<double>& result) const { result.assign(z.size(), value); }
	void set_value(double value_) { value = value_; }
private:
	double value;
};

}

class CoulombTest : public CxxTest::TestSuite {

public:
	CoulombTest() {
		fout.open("18_CoulombMatrixTest_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);
		if(be_chatty) {
			Logger::get_instance()->set_level(LOG_INFO_DEVEL2);
		} else {
			Logger::get_instance()->del_listener(&std::cout);
		}
		// --------------------------------------
		// init material database and geometry
		// --------------------------------------
		material_database.add_search_path("18_CoulombMatrix/");
		InputParser parser;
		geometry = parser.read_geometry(AsciiGridReader("18_CoulombMatrix/well384.asc.gz"));
		geometry->set_materials(material_database);
		geometry->set_boundary_conditions(new BCIncludeAll(*geometry));

		// -----------------------------------------------
		// construct four wavefunctions
		// [cos(z), 0]
		// [cos(z), sin(z)]
		// [sin(z), cos(z)]
		// [0,      sin(z)]
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < 4; ii++) {
			dummy_solutions.push_back(
				new EigenSolution<cplx>(
					geometry->get_num_nonzero_nodes(), 2
				)
			);
		}
		double fact = constants::pi / 7.5;
		for(unsigned int ii = 0; ii < geometry->get_num_nodes(); ii++) {
			const Node& vert = geometry->get_node(ii);
			if(vert.get_index_internal() != -1) {
				// sol 1
				dummy_solutions[0]->set_node_value(
					vert.get_index_internal(),
					0,
					cos(fact * vert.get_coord(0))
				);
				// sol 2
				dummy_solutions[1]->set_node_value(
					vert.get_index_internal(),
					0,
					cos(fact * vert.get_coord(0))
				);
				dummy_solutions[1]->set_node_value(
					vert.get_index_internal(),
					1,
					sin(fact * vert.get_coord(0))
				);
				// sol 3
				dummy_solutions[2]->set_node_value(
					vert.get_index_internal(),
					0,
					sin(fact * vert.get_coord(0))
				);
				dummy_solutions[2]->set_node_value(
					vert.get_index_internal(),
					1,
					cos(fact * vert.get_coord(0))
				);
				// sol 4
				dummy_solutions[3]->set_node_value(
					vert.get_index_internal(),
					1,
					sin(fact * vert.get_coord(0))
				);

			}
		}



	}


	~CoulombTest() {
		Logger::get_instance()->del_listener(&fout);
		fout.close();
		delete geometry;
	}

	Geometry* geometry;
	MaterialDatabase material_database;

	vector<EigenSolution<cplx>*> dummy_solutions;

	std::fstream fout;
	void setUp() {}
	void tearDown() {}
	void test_coulomb_matrix_element();
	void test_element_middle_points();
	void test_coulomb_function_well();
	void test_coulomb_function_wire();
	void test_coulomb_integrator();
	void test_coulomb_lambda();

};


void CoulombTest::test_element_middle_points() {

	try {
		ElementMiddlePoints midpoints(*geometry);

		// test coordinates
		double coord[3];
		for(unsigned int ii = 0; ii < geometry->get_num_elements(); ii++) {
			for(unsigned int dd = 0; dd < geometry->get_dimension(); dd++) {
				coord[dd] = 0.0;
			}
			const Element& elem = geometry->get_element(ii);
			for(unsigned int vv = 0; vv < elem.get_num_nodes(); vv++) {
				for(unsigned int dd = 0; dd < geometry->get_dimension(); dd++) {
					coord[dd] += elem.get_node(vv).get_coord(dd);
				}
			}
			const double* ref_coords = midpoints.get_coords(ii);
			for(unsigned int dd = 0; dd < geometry->get_dimension(); dd++) {
				coord[dd] /= elem.get_num_nodes();
				TS_ASSERT_DELTA(coord[dd], ref_coords[dd], 1.0e-10);
			}

		}
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}
}

void CoulombTest::test_coulomb_function_well() {
	try {
		// ----------------------------------------
		// for testing we simply integrate for q = 2
		// and z' such that we have 2pi exp(-2x) x in [0, inf]
		// because that must be pi
		// ----------------------------------------
		int num = 100000;
		vector<double> z;
		vector<double> r(num, 0.0e0);
		double dz = 0.005;
		for(int ii = 0; ii < num; ii++) {
			z.push_back(static_cast<double>(ii) * dz);
		}
		CoulombFunctionWell cf;
		cf.set(2.0, 0.0);
		cf.evaluate(z,r);
		double test = 0.0;
		for(int ii = 0; ii < num-1; ii++) {
			test += dz * (r[ii] + r[ii+1]) / 2;
		}
		TS_ASSERT_DELTA(test, constants::pi, 1.0e-2);
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}
}

void CoulombTest::test_coulomb_function_wire() {
	try {
		// ----------------------------------------
		// for testing we simply integrate from 0 to
		// infinity for q = 1 (must be pi) because
		// bessel K0 integrated gives pi / 2 and
		// we have 2*K0
		// ----------------------------------------
		int num = 100000;
		vector<double> z;
		vector<double> r(num, 0.0e0);
		double dz = 0.005;
		for(int ii = 0; ii < num; ii++) {
			z.push_back(static_cast<double>(ii) * dz);
			z.push_back(0.0);
		}
		CoulombFunctionWire cf;
		double z_prime[2] = {0.0, 0.0};
		cf.set(1.0, z_prime);
		cf.evaluate(z,r);
		double test = 0.0;
		for(int ii = 0; ii < num; ii++) {
			test += dz * r[ii];
		}
		TS_ASSERT_DELTA(test, constants::pi, 0.05);
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}
}

void CoulombTest::test_coulomb_integrator() {

	try {
		ElementMiddlePoints midpoints(*geometry);
		CoulombIntegrator<double> integrator(*geometry, material_database, midpoints);
		CoulombFunctionMOC cfmoc(1.0);
		vector<cplx> ones(geometry->get_num_nodes(), 1.0);
		vector<cplx> out(geometry->get_num_nodes());
		double dn = 0.5;
		double vl = 0.5;
		while(vl < 2) {
			cfmoc.set_value(vl);
			integrator.build_matrix(cfmoc);
			integrator.multiply_with_lhs(ones, out);
			cplx res = 0.0;
			for(unsigned int ii = 0; ii < out.size(); ii++) {
				res += out[ii];
			}
			TS_ASSERT_DELTA(abs(res), vl * 15.0, 1.0e-5);
			vl += dn;
		}
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}
}





void CoulombTest::test_coulomb_matrix_element() {

	try {
		CoulombMatrixElementQuantized matelem(*geometry, material_database, 2);
		matelem.set_q_range(1.0e-6, 1.0, 5);
		for(unsigned int ii = 0; ii < dummy_solutions.size(); ii++) {
			matelem.add_wavefunction(dummy_solutions[ii]);
		}
		matelem.request_matrix_element(0,0,0,0);
		matelem.request_matrix_element(1,0,0,1);
		matelem.request_matrix_element(0,2,2,0);
		matelem.request_matrix_element(0,2,0,2);
		matelem.calculate();
		// -----------------------------------------------
		// construct four wavefunctions
		// [cos(z), 0]
		// [cos(z), sin(z)]
		// [sin(z), cos(z)]
		// [0,      sin(z)]
		// -----------------------------------------------
		const vector<cplx>& nzzzz = matelem.get_matrix_element(0,0,0,0);
		const vector<cplx>& nezze = matelem.get_matrix_element(1,0,0,1);
		const vector<cplx>& nzttz = matelem.get_matrix_element(0,2,2,0);
		const vector<cplx>& nztzt = matelem.get_matrix_element(0,2,0,2);
		// mathematica results ...
		double tol = 1.0e-1;

		// no eval at q = 0 TS_ASSERT_DELTA(1.0, nzzzz[0].real() / (2.0 * constants::pi / 0.01 * 56.25),   tol);
		TS_ASSERT_DELTA(1.0, nzzzz[1].real() / (2.0 * constants::pi * 22.0849), tol); // q = 0.25
		TS_ASSERT_DELTA(1.0, nzzzz[2].real() / (2.0 * constants::pi * 13.7831), tol); // q = 0.5
		TS_ASSERT_DELTA(1.0, nzzzz[3].real() / (2.0 * constants::pi * 10.3688), tol); // q = 0.75
		TS_ASSERT_DELTA(1.0, nzzzz[4].real() / (2.0 * constants::pi * 8.44326), tol); // q = 1.0

		// no eval at q = 0 TS_ASSERT_DELTA(1.0, nezze[0].real() / (2.0 * constants::pi / 0.01 * 112.5),   tol);
		TS_ASSERT_DELTA(1.0, nezze[1].real() / (2.0 * constants::pi * 43.0987), tol); // q = 0.25
		TS_ASSERT_DELTA(1.0, nezze[2].real() / (2.0 * constants::pi * 24.9522), tol); // q = 0.5
		TS_ASSERT_DELTA(1.0, nezze[3].real() / (2.0 * constants::pi * 17.4313), tol); // q = 0.75
		TS_ASSERT_DELTA(1.0, nezze[4].real() / (2.0 * constants::pi * 13.4124), tol); // q = 1.0

		// no eval at q = 0 TS_ASSERT_DELTA(1.0, nzttz[0].real() / (2.0 * constants::pi / 0.01 * 112.5),   tol);
		TS_ASSERT_DELTA(1.0, nzttz[1].real() / (2.0 * constants::pi * 43.0987), tol);
		TS_ASSERT_DELTA(1.0, nzttz[2].real() / (2.0 * constants::pi * 24.9522), tol);
		TS_ASSERT_DELTA(1.0, nzttz[3].real() / (2.0 * constants::pi * 17.4313), tol);
		TS_ASSERT_DELTA(1.0, nzttz[4].real() / (2.0 * constants::pi * 13.4124), tol);

		// no eval at q = 0 TS_ASSERT_DELTA(0.0, nztzt[0].real(),   tol);
		TS_ASSERT_DELTA(1.0, nztzt[1].real() / (2.0 * constants::pi * 1.81309), tol);
		TS_ASSERT_DELTA(1.0, nztzt[2].real() / (2.0 * constants::pi * 2.35699), tol);
		TS_ASSERT_DELTA(1.0, nztzt[3].real() / (2.0 * constants::pi * 2.444), tol);
		TS_ASSERT_DELTA(1.0, nztzt[4].real() / (2.0 * constants::pi * 2.32466), tol);

		// -------------------------------------------
		// write binary data
		// -------------------------------------------
		const CoulombMatrixElementData& writing = matelem.get_data_object();
		writing.write_binary("coulomb.bin");

		// -------------------------------------------
		// read binary data and compare
		// -------------------------------------------
		CoulombMatrixElementData reading("coulomb.bin");

		// compare q values
		TS_ASSERT_EQUALS(reading.get_q_values().size(), writing.get_q_values().size());
		for(unsigned int ii = 0; ii < reading.get_q_values().size() && ii < writing.get_q_values().size(); ii++) {
			TS_ASSERT_EQUALS(reading.get_q_values()[ii], writing.get_q_values()[ii]);
		}
		// compare matrix elements (last set)
		const vector<cplx>& reading_data = reading.get_matrix_element(0,2,2,0);
		const vector<cplx>& writing_data = writing.get_matrix_element(0,2,2,0);
		if(reading_data.size() == writing_data.size()) {
			for(unsigned int ii = 0; ii < reading_data.size(); ii++) {
				TS_ASSERT_EQUALS(reading_data[ii], writing_data[ii]);
			}
		} else {
			TS_FAIL("reading_data.size() == writing_data.size() failed!");
		}

	} catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}

}

void CoulombTest::test_coulomb_lambda() {
	try {
		// ----------------------------------
		// test wire
		// ----------------------------------
		// build q-spaces and 1/q^2 values <- thats what we feed to mimic bessel function
		vector<double> q_values;
		vector<cplx>   values;
		tdkp_math::linear_space(q_values, 1.0e-3, 2.4, 50);

		for(unsigned int ii = 0; ii < q_values.size(); ii++) {
			values.push_back(1.0 / (q_values[ii] * q_values[ii]));

		}
		CLCCoulombLambdaWire lambda(q_values, values);

		// evaluate ... (test if we get the bessel function out)
		// note, the coulomb lambda function gives F(|k-k'|) + F(|k+k')
		// so we should expect values * 2		
		for(unsigned int ii = 3; ii < q_values.size(); ii++) {
			// check at data points
			TS_ASSERT_DELTA(2.0 * values[ii].real(), lambda.get_coulomb_lambda(q_values[ii], 0.0), 1.0e-5);
			// and verify in between
			double mid = lambda.get_coulomb_lambda(0.5 * (q_values[ii] + q_values[ii - 1]) , 0.0);
			TS_ASSERT(mid <= 2 * values[ii - 1].real());
			TS_ASSERT(mid >= 2 * values[ii].real());
			double qtmp = 0.5 * (q_values[ii] + q_values[ii - 1]);
			// divison by 2 -> the coulomb lambda returns us 2 times the values we gave
			TS_ASSERT_DELTA(1.0, (1.0 / (qtmp * qtmp)) / (mid / 2.0), 0.03);
		}

	} catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}
	// ---------------------------------------------
	// test quantum well
	// ---------------------------------------------
	try {
		vector<double> q_values;
		vector<cplx>   values;
		// we create 1/q
		tdkp_math::linear_space(q_values, 1.0e-2, 2.4, 50);
		for(unsigned int ii = 0; ii < q_values.size(); ii++) {
			values.push_back(1.0/ q_values[ii]);
		}
		// and now integrate 1/|k - k'| over the sphere
		CLCRadialCoulombLambdaWell lambda(q_values, values, 450);
		// mathematica results
		TS_ASSERT_DELTA(lambda.get_coulomb_lambda(1.0, 0.5) / 6.743, 1.0,   1.0e-2);
		TS_ASSERT_DELTA(lambda.get_coulomb_lambda(1.0, 0.8) / 7.98121, 1.0, 1.0e-2);


	}  catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}
	// ---------------------------------------------
	// test bulk
	// ---------------------------------------------
	try {
		{ // constant integration of 1 (means result is 4pi (as 1/q^2 was moved out of coulomb function)
			vector<double> q_values;
			vector<cplx>   values;
			// we create 1/q
			tdkp_math::linear_space(q_values, 1.0e-2, 2.4, 50);
			for(unsigned int ii = 0; ii < q_values.size(); ii++) {
				values.push_back(1.0);
			}
			// and now integrate 1/|k - k'| over the sphere
			CLCRadialCoulombLambdaBulk lambda(q_values, values, 450);
			TS_ASSERT_DELTA(lambda.get_coulomb_lambda(1.0, 0.0), 4.0 * constants::pi, 0.03);
		}
		{ // some more complicated integration, using SPLINES
			vector<double> q_values;
			vector<cplx>   values;
			// we create 1/q
			tdkp_math::linear_space(q_values, 1.0e-2, 2.4, 50);
			for(unsigned int ii = 0; ii < q_values.size(); ii++) {
				values.push_back(1.0 / (q_values[ii] * q_values[ii]));
			}
			// and now integrate 1/|k - k'| over the sphere
			CLCRadialCoulombLambdaBulk lambda(q_values, values, 900);
			TS_ASSERT_DELTA(lambda.get_coulomb_lambda(1.0, 0.0), 4.0 * constants::pi, 0.03);

			// mathematica results
			TS_ASSERT_DELTA(lambda.get_coulomb_lambda(1.0, 0.5) / 13.8056, 1.0, 1.0e-2);
			TS_ASSERT_DELTA(lambda.get_coulomb_lambda(1.0, 0.8) / 17.257,  1.0, 1.0e-2);
		}

	}  catch (Exception* e) {
		TS_FAIL(e->get_reason());
	}
}

#endif /*DEFINITION_H_*/
