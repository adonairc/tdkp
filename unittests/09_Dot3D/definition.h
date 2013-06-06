
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
#include <map>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <stdio.h>
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/geometry/Region.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/geometry/Element3DTR.h"
#include "tdkp/geometry/Element2DRect.h"
#include "tdkp/geometry/Element2DTriangle.h"
#include "tdkp/main/Fields.h"
#include "tdkp/io/InputParser.h"
#include "tdkp/io/AsciiGridReader.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/main/EigenSolution.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/kpmatrices/KPMatrix4x4EndersForeman.h"
#include "tdkp/probdefs/EffectiveMass.h"
#include "tdkp/probdefs/KP4x43D.h"
#include "tdkp/probdefs/KP8x83D.h"


using namespace tdkp;

class Dot3DTest : public CxxTest::TestSuite {
public:	

	fstream fout;
		
	Dot3DTest() { 				
 		Configuration::get_instance()->set("output_eigensolver_statistics", 		1);
		Configuration::get_instance()->set("assembly_build_nonsymmetric_matrices", 	1);
		Configuration::get_instance()->set("assembly_check_matrix_for_symmetry",   	1);
		Configuration::get_instance()->set("assembly_save_matrices_to_file",		0);					
		fout.open("09_Dot3DTest_output.log", ios::app | ios::out);

		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);
		Logger::get_instance()->add_listener(&fout);
		Logger::get_instance()->del_listener(&std::cout);	
		if(!write_mode) { TS_FAIL("WARNING, I AM IN WRITE MODE"); }		
	}	
	
	~Dot3DTest() {
		fout.close();
	}	
	
	static const bool write_mode = false;
	
	void test_effmass_3D() {
		int num = 8;		
		InputParser 		parser;
		MaterialDatabase 	mat_db;
		Geometry*           geometry = 0;		
		mat_db.add_search_path("./09_Dot3D/");
		geometry = parser.read_geometry(AsciiGridReader("./09_Dot3D/dot1400.asc.gz"));
		geometry->set_materials(mat_db);
		EffectiveMass problem1(*geometry, mat_db);
		try {
			problem1.solve(num);
			problem1.display_solution_info();
			string filename = "./09_Dot3D/effmass_compare_data.dat";
			if(write_mode) {
				EigenSolutionSet<double>* set = problem1.get_probability_objects();
				set->write(filename.c_str());
				delete set;
			} else {
				EigenSolutionSet<double> set;
				set.read(filename.c_str());
				for(int ii = 0; ii < num; ii++) {				
					EigenSolution<double>* tmp = problem1.get_probability_object(ii);
					TS_ASSERT(set.get(ii).compare(*tmp, 1.0e-8, true));
					delete tmp;	
				}	
			}
		} catch(Exception*e) {
			TS_FAIL(e->get_reason());	
		}
	}
	
	void test_effmass_hole_3D() {
		int num = 8;		
		InputParser 		parser;
		MaterialDatabase 	mat_db;
		Geometry*           geometry = 0;
		mat_db.add_search_path("./09_Dot3D/");				
		geometry = parser.read_geometry(AsciiGridReader("./09_Dot3D/dot1400.asc.gz"));
		geometry->set_materials(mat_db);
		EffectiveMass problem1(*geometry, mat_db);		
		problem1.set_solution_type(holes);
		
		try {
			problem1.solve(num);
			problem1.display_solution_info();
			string filename = "./09_Dot3D/effmasshole_compare_data.dat";
			if(write_mode) {
				EigenSolutionSet<double>* set = problem1.get_probability_objects();
				set->write(filename.c_str());
				delete set;
			} else {
				EigenSolutionSet<double> set;
				set.read(filename.c_str());
				for(int ii = 0; ii < num; ii++) {				
					EigenSolution<double>* tmp = problem1.get_probability_object(ii); 
					TS_ASSERT(set.get(ii).compare(*tmp, 1.0e-8));
					delete tmp;	
				}	
			}
		} catch(Exception*e) {
			TS_FAIL(e->get_reason());	
		}			
	}
	void test_kp4x4_3D() {
		int num = 8;		
		InputParser 		parser;
		MaterialDatabase 	mat_db;
		Geometry*           geometry = 0;
		mat_db.add_search_path("./09_Dot3D/");
		geometry = parser.read_geometry(AsciiGridReader("./09_Dot3D/dot1400.asc.gz"));
		geometry->set_materials(mat_db);		
		KP4x43D             problem1(*geometry, mat_db);

		try {	
			problem1.set_axes(
				Vector3D(1.0,0.0,0.0),
				Vector3D(0.0,1.0,0.0),
				Vector3D(0.0,0.0,1.0)
			);			 		
			problem1.solve(num);
			problem1.display_solution_info();
			string filename = "./09_Dot3D/kp4x4_compare_data.dat";		
			if(write_mode) {
				EigenSolutionSet<double>* set = problem1.get_probability_objects();
				set->write(filename.c_str());
				delete set;
			} else {
				EigenSolutionSet<double> set;
				set.read(filename.c_str());
				for(int ii = 0; ii < num; ii++) {				
					EigenSolution<double>* tmp = problem1.get_probability_object(ii); 
					TS_ASSERT(set.get(ii).compare(*tmp));
					delete tmp;	
				}	
			}
		} catch(Exception*e) {
			TS_FAIL(e->get_reason());	
		}		
	}
	
	void test_kp8x8e_3D() {
		
		int num = 8;		
		InputParser 		parser;
		MaterialDatabase 	mat_db;
		Geometry*           geometry = 0;
		mat_db.add_search_path("./09_Dot3D/");
		geometry = parser.read_geometry(AsciiGridReader("./09_Dot3D/dot1400.asc.gz"));
		geometry->set_materials(mat_db);		
		KP8x83D             problem1(*geometry, mat_db);		
		problem1.set_solution_type(electrons);
		TS_ASSERT_THROWS_NOTHING(
			problem1.set_axes(
				Vector3D(1.0,0.0,0.0),
				Vector3D(0.0,1.0,0.0),
				Vector3D(0.0,0.0,1.0)
			)
		); 		
		TS_ASSERT_THROWS_NOTHING(problem1.solve(num));
		TS_ASSERT_THROWS_NOTHING(problem1.display_solution_info());
		
		string filename = "./09_Dot3D/kp8x8e_compare_data.dat";		
		if(write_mode) {
			EigenSolutionSet<double>* set = problem1.get_probability_objects();
			set->write(filename.c_str());
			delete set;
		} else {
			EigenSolutionSet<double> set;
			set.read(filename.c_str());
			for(int ii = 0; ii < num; ii++) {				
				EigenSolution<double>* tmp = problem1.get_probability_object(ii); 
				TS_ASSERT(set.get(ii).compare(*tmp));
				delete tmp;	
			}	
		}
	}	
	
	void test_kp8x8_holes_3D() {
		int num = 8;		
		InputParser 		parser;
		MaterialDatabase 	mat_db;
		Geometry*           geometry = 0;
		mat_db.add_search_path("./09_Dot3D/");
		geometry = parser.read_geometry(AsciiGridReader("./09_Dot3D/dot1400.asc.gz"));
		geometry->set_materials(mat_db);		
		KP8x83D             problem1(*geometry, mat_db);
		problem1.set_solution_type(holes);
	
		TS_ASSERT_THROWS_NOTHING(
			problem1.set_axes(
				Vector3D(1.0,0.0,0.0),
				Vector3D(0.0,1.0,0.0),
				Vector3D(0.0,0.0,1.0)
			)
		); 		
		problem1.solve(num);
		problem1.display_solution_info();
		string filename = "./09_Dot3D/kp8x8h_compare_data.dat";
		if(write_mode) {
			EigenSolutionSet<double>* set = problem1.get_probability_objects();
			set->write(filename.c_str());
			delete set;
		} else {
			EigenSolutionSet<double> set;
			set.read(filename.c_str());
			for(int ii = 0; ii < num; ii++) {				
				EigenSolution<double>* tmp = problem1.get_probability_object(ii); 
				TS_ASSERT(set.get(ii).compare(*tmp, 1.0e-8));
				delete tmp;	
			}	
		}
	}				
};

