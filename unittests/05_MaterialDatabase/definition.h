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
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <map>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/main/MaterialDatabase.h"

using namespace tdkp;

class MaterialTest : public CxxTest::TestSuite {
public:	
	std::fstream fout;	
	void setUp() {
		fout.open("05_MaterialDatabase_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);
		Logger::get_instance()->del_listener(&std::cout);	
		system("cp 05_MaterialDatabase/testfiles/valid.cnf .");
		system("cp 05_MaterialDatabase/testfiles/material_a.mat .");
		system("cp 05_MaterialDatabase/testfiles/material_b.mat .");
		system("cp 05_MaterialDatabase/testfiles/material_invalid.mat .");		
	}
	void tearDown() {
		fout.close();	
		system("rm valid.cnf material_a.mat material_b.mat material_invalid.mat");
	}	
	
	void test_general() {
		MaterialDatabase mat;
		string tmp;
		
		TS_ASSERT_THROWS_NOTHING(mat.load_material("material_a"));
		//TS_ASSERT_THROWS_NOTHING(mat.load_material(tmp = "material_b"));		
		// material is invalid (wurstsalat is missing, but it should load
		TS_ASSERT_THROWS_NOTHING(mat.load_material("material_invalid"));
		
		mat.add_search_path("05_MaterialDatabase/testfiles/");
		TS_ASSERT_THROWS_NOTHING(mat.load_material("material_otherdir"));
		
		TS_ASSERT(mat.material_exists("material_a"));
	//	TS_ASSERT(mat.material_exists(tmp = "material_a"));
		
		TS_ASSERT(!mat.material_is_valid("material_invalid"));
		//TS_ASSERT(!mat.material_is_valid(tmp = "material_invalid"));
		
		TS_ASSERT(mat.material_is_valid("material_a"));
		//TS_ASSERT(mat.material_is_valid(tmp = "material_a"));		
		
		TS_ASSERT(mat.property_is_set("material_a", "wurstsalat"));
		TS_ASSERT_THROWS_ANYTHING(!mat.property_is_set("material_a", "notset"));
		TS_ASSERT(!mat.property_is_set("material_a", "bohnensalat"));
		
		TS_ASSERT_THROWS_ANYTHING(mat.property_is_set("nonexistent", "nonexistent"));
				
		TS_ASSERT_EQUALS(mat.get("material_a", "wurstsalat"), 4.0);
		TS_ASSERT_EQUALS(mat.get("material_a", "krautsalat"), 1.5);
						 				
	}
	
};

