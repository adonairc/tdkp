#ifndef DEFINITION_H_
#define DEFINITION_H_

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
#include "tdkp/probdefs/BulkBandstructureSolver.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/kpmatrices/KPMatrixBase.h"
#include "tdkp/kpmatrices/KPMatrix4x4EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix6x6EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix8x8EndersForeman.h"
#include "ValueChecker.h"

using namespace tdkp;


class BulkBandstructureTest : public CxxTest::TestSuite, public ResultTestBase {


public:	
		
	BulkBandstructureTest();
	
	
	void setUp() {
		fout.open("12_BulkBandstructureTest_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);	
		Logger::get_instance()->del_listener(&std::cout);		
	}
	void tearDown() {
		fout.close();	
	}
	
	
	
	void test_kp4x4_simple();	
	void test_kp6x6_simple();
	void test_kp8x8_simple();
	void test_kp4x4_rotated();
	void test_kp6x6_rotated();
	void test_kp8x8_rotated();	
	void test_kp4x4_rotated_vs_simple();
	void test_kp6x6_rotated_vs_simple();
	void test_kp8x8_rotated_vs_simple();
	
	void check_bands(KPMatrixBase& kp_matrix, const char* kptype);
	void check_bands_rotated(const string& message, KPMatrixBase& kp_matrix, KPMatrixBase& kp_copy);
	
protected:
	struct CalcSetup {
		string material_name;
		string direction;
	};
	vector<CalcSetup>    calculations;
	map<string,Vector3D> direction_map;		
	RMatrix<double>      rotation;
	
	
};

// -------------------------------------------------------
// BulkBandstructureTest implementation
// -------------------------------------------------------
BulkBandstructureTest::BulkBandstructureTest()
: ResultTestBase("./12_BulkBandstructure/"),
  rotation(3,3) 
{
	this->write_mode = false;
	
	direction_map[string("100")] = one_zero_zero;
	direction_map[string("110")] = one_one_zero;
	direction_map[string("111")] = one_one_one;
	
	CalcSetup s;
	s.material_name = "InAs";
	s.direction     = "100";
	calculations.push_back(s);
	s.direction     = "110";
	calculations.push_back(s);
	s.direction     = "111";
	calculations.push_back(s);
	s.material_name = "GaAs";
	s.direction     = "100";
	calculations.push_back(s);
	s.direction     = "110";
	calculations.push_back(s);
	s.direction     = "111";
	calculations.push_back(s);
	
	
	double  angle = 2.0 * constants::pi * 45.0 / 360.0;	
	rotation(0,0) = rotation(1,1) = cos(angle);
	rotation(0,1) = - sin(angle);
	rotation(1,0) = - rotation(0,1);
	rotation(2,2) = 1.0;
	
}

void BulkBandstructureTest::test_kp4x4_simple() {
	KPMatrix4x4EndersForeman kp_matrix;
	check_bands(kp_matrix, "kp4x4_simple");	
}

void BulkBandstructureTest::test_kp6x6_simple() {
	KPMatrix6x6EndersForeman kp_matrix;
	check_bands(kp_matrix, "kp6x6_simple");	
}

void BulkBandstructureTest::test_kp8x8_simple() {
	KPMatrix8x8EndersForeman kp_matrix;
	check_bands(kp_matrix, "kp8x8_simple");	
}

void BulkBandstructureTest::test_kp4x4_rotated() {
	KPMatrix4x4EndersForeman kp_matrix;
	kp_matrix.set_rotation(rotation);
	check_bands(kp_matrix, "kp4x4_rotated");	
}

void BulkBandstructureTest::test_kp6x6_rotated() {
	KPMatrix6x6EndersForeman kp_matrix;	
	kp_matrix.set_rotation(rotation);	
	check_bands(kp_matrix, "kp6x6_rotated");	
}

void BulkBandstructureTest::test_kp8x8_rotated() {
	KPMatrix8x8EndersForeman kp_matrix;
	kp_matrix.set_rotation(rotation);	
	check_bands(kp_matrix, "kp8x8_rotated");
}

void BulkBandstructureTest::test_kp4x4_rotated_vs_simple() {
	KPMatrix4x4EndersForeman kp_matrix;
	KPMatrix4x4EndersForeman kp_copy;	
	check_bands_rotated("kp4x4_rotated_vs_simple", kp_matrix, kp_copy);	
}

void BulkBandstructureTest::test_kp6x6_rotated_vs_simple() {
	KPMatrix6x6EndersForeman kp_matrix;
	KPMatrix6x6EndersForeman kp_copy;	
	check_bands_rotated("kp6x6_rotated_vs_simple", kp_matrix, kp_copy);	
}

void BulkBandstructureTest::test_kp8x8_rotated_vs_simple() {
	KPMatrix8x8EndersForeman kp_matrix;
	KPMatrix8x8EndersForeman kp_copy;	
	check_bands_rotated("kp8x8_rotated_vs_simple", kp_matrix, kp_copy);	
}

void BulkBandstructureTest::check_bands(KPMatrixBase& kp_matrix, const char* kptype) {
	
	try {					
		for(unsigned int ii = 0; ii < calculations.size(); ii++) {		
			kp_matrix.set_material(matdb.get_material(calculations[ii].material_name.c_str()));
			BulkBandstructureSolver  solver(&kp_matrix);

			solver.solve(direction_map[calculations[ii].direction], kmin, kmax, knum);
			// compare
			const BandstructureDomain<complex<double> >& bs = solver.get_bandstructure();
			
			ostringstream sout;
			sout << kptype << "_" << calculations[ii].material_name << "_" 
			     << calculations[ii].direction << ".dat";
			if(!compare_bandstructures(kptype, sout.str().c_str(), bs)) {
				ostringstream smsg;
				smsg << kptype << " bandstructure comparison failed for " << sout.str();
				TS_FAIL(smsg.str());	
			}
		}	
	} catch (Exception* e) {		
		TS_FAIL(e->get_reason()); 	
	} 
}

void BulkBandstructureTest::check_bands_rotated(const string& message, KPMatrixBase& kp_matrix, KPMatrixBase& kp_copy) {
	
	
	Vector3D random_orientation = ValueChecker<double>::get_random_dir();
	random_orientation.normalize();
	double angle = 2.0 * constants::pi * ((0.5 - drand48()) * 2.0);
	
	RMatrix<double> my_rotation = RMatrix<double>::get_rotation_matrix(random_orientation, cos(angle));	
	

	Vector3D rotated_dir;

	// rotate rotated system
	kp_matrix.set_rotation(my_rotation);
		
	cout << "rotation: " << my_rotation << endl;
			
	ostringstream fout1;
				
				
	try {					
		for(unsigned int ii = 0; ii < calculations.size(); ii++) {
			
			/*
			TDKP_TRACE("passed here with ii " << ii);			
			// output
			cout << "unrotated direction:\n"
			     << direction_map[calculations[ii].direction] << "\n"
			     << "rotated direction:\n"
			     << rotation * direction_map[calculations[ii].direction]			     
			     << endl; 
			*/
			// set materials
			kp_matrix.set_material(matdb.get_material(calculations[ii].material_name.c_str()));
			/*
			fout1.str(""); 
			fout1 << "kpmat_rotated_" << message << "_"
			      << calculations[ii].material_name << "_"
			      << calculations[ii].direction << ".dat";
			kp_matrix.dump(fout1.str().c_str());
			*/						
			kp_copy.set_material(matdb.get_material(calculations[ii].material_name.c_str()));
			/*
			fout1.str(""); 			
			fout1 << "kpmat_straight_" << message << "_"
			      << calculations[ii].material_name << "_"
			      << calculations[ii].direction << ".dat";
			kp_matrix.dump(fout1.str().c_str());
			*/
			
			
			// evaluate rotated system					
			BulkBandstructureSolver  solver_rotated(&kp_matrix);
			rotated_dir = my_rotation * direction_map[calculations[ii].direction];			
			solver_rotated.solve(rotated_dir, kmin, kmax, knum);
						
			const BandstructureDomain<complex<double> >& bs_rotated = solver_rotated.get_bandstructure();
			// evaulate straight system
			BulkBandstructureSolver solver_straight(&kp_copy);
			solver_straight.solve(direction_map[calculations[ii].direction], kmin, kmax, knum);
			const BandstructureDomain<complex<double> >&bs_straight = solver_straight.get_bandstructure();
			// compare
			if(!compare_bandstructures(
					message.c_str(), 
					bs_rotated,
					bs_straight
			  )) {
				ostringstream smsg;
				smsg << "bandstructure comparison failed for " << smsg.str();
				TS_FAIL(smsg.str());	
			}
		}	
	} catch (Exception* e) {		
		TS_FAIL(e->get_reason()); 	
	} 
}

#endif /*DEFINITION_H_*/
