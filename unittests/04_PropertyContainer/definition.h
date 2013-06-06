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
#include <boost/lexical_cast.hpp>
#include "tdkp/main/PropertyContainer.h"

using namespace tdkp;

class CSCtest : public CxxTest::TestSuite {
public:	
	std::fstream fout;	
	void setUp() {
		fout.open("04_PropertyContainer_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);	
		Logger::get_instance()->del_listener(&std::cout);
	}
	void tearDown() {
		fout.close();	
	}	
	
	// --------------------------------------------------------
	// file valid.cnf
	// # this is a funny comment
	// wurstsalat; mandatory; 0.0; 5.0; 0.1; 4.5;
	// bohnensalat; optional; 0.0; 10.0; 1.0; 2.0;
	// krautsalat; optional; 1.0; 1111.123; 20.0; 1000; 34.0;
	// --------------------------------------------------------		
	void test_valid() {
		string tmp;
				
		try {
			PropertyContainer<double> cont2("04_PropertyContainer/testfiles/valid.cnf");
		} catch(Exception* e) {
			TS_FAIL(e->get_reason().c_str());
			TS_FAIL("stopping");
			return;	
		}
		PropertyContainer<double> cont("04_PropertyContainer/testfiles/valid.cnf");
		
		// test valid keys
		TS_ASSERT(cont.valid_key("wurstsalat"));
		TS_ASSERT(cont.valid_key("bohnensalat"));
		TS_ASSERT(cont.valid_key("krautsalat"));
		TS_ASSERT(!cont.valid_key("xxx_wurstsalat"));
		TS_ASSERT(!cont.valid_key("xxx_bohnensalat"));
		TS_ASSERT(!cont.valid_key("xxx_krautsalat"));

		TS_ASSERT(cont.valid_key(tmp = "bohnensalat"));
		TS_ASSERT(!cont.valid_key(tmp = "aaabohnensalat"));
		
		// test if valid value throws exception if illegal key
		TS_ASSERT_THROWS_ANYTHING(cont.valid_value("wurstsalattttttttt", 2.5));		
		TS_ASSERT_THROWS_ANYTHING(cont.valid_value(tmp = "wurstsalattttttttt", 2.5));
		
		// test if valid value returns true on correct values and false on wrong
		TS_ASSERT(cont.valid_value("wurstsalat", 2.5));
		TS_ASSERT(cont.valid_value(tmp = "wurstsalat", 2.5));
		TS_ASSERT(!cont.valid_value("wurstsalat", -2.0));
		TS_ASSERT(!cont.valid_value(tmp = "wurstsalat", -2.0));		
		
		// should only warn but return true
		TS_ASSERT(cont.valid_value("wurstsalat", 0.05));
		TS_ASSERT(cont.valid_value(tmp = "wurstsalat", 0.05));	
		
		// should evaluate as invalid as wurstsalat not set yet
		TS_ASSERT(!cont.valid());
				
		// try to set value to invalid key
		TS_ASSERT_THROWS_ANYTHING(cont.set("____invalid", 2.0));
		TS_ASSERT_THROWS_ANYTHING(cont.set(tmp = "____invalid", 2.0));		
		
		// try to set invalid value to valid key
		TS_ASSERT_THROWS_ANYTHING(cont.set("wurstsalat", -2.0));
		TS_ASSERT_THROWS_ANYTHING(cont.set(tmp = "wurstsalat", -2.0));		
		
		// try to set valid value to valid key
		TS_ASSERT_THROWS_NOTHING(cont.set("wurstsalat", 2.0));
		TS_ASSERT_THROWS_NOTHING(cont.set(tmp = "wurstsalat", 2.0));		
		
		// now as the only mandatory field is set, it should evaluate as valid
		TS_ASSERT(cont.valid());
		
		// try to get value from invalid key
		TS_ASSERT_THROWS_ANYTHING(cont.get("__invalid__"));		
		TS_ASSERT_THROWS_ANYTHING(cont.get(tmp = "__invalid__"));
						
		// try to get valid value from valid key
		TS_ASSERT_EQUALS(cont.get("wurstsalat"), 2.0);
		TS_ASSERT_EQUALS(cont.get(tmp = "wurstsalat"), 2.0);
		
		// try to get value from valid but unset key
		TS_ASSERT_THROWS_ANYTHING(cont.get("bohnensalat"));
		TS_ASSERT_THROWS_ANYTHING(cont.get(tmp = "bohnensalat"));
		
		// check if there is a default value for krautsalat
		TS_ASSERT_EQUALS(cont.get("krautsalat"), 34.0);
		TS_ASSERT_EQUALS(cont.get(tmp = "krautsalat"), 34.0);		
						
	}
	
			
	void test_invalid_cf() {
		string files[] = {"invalid_wrong_bounds.cnf",
						 "invalid_wrong_bounds2.cnf",
						 "invalid_wrong_bounds3.cnf",
						 "mandatory_with_default.cnf"};
		std::string path = "04_PropertyContainer/testfiles/";
		std::string file;
		for(int ii = 0; ii < 4; ii++) { 
			file = path;
			file.append(files[ii].c_str());
			TS_ASSERT_THROWS_ANYTHING(PropertyContainer<double> cond(file.c_str()));
		}									 			
	}			
	
	void test_read_data_from_valid_file() {
		std::string path = "04_PropertyContainer/testfiles/";
		std::string file;						
		// read in definition
		PropertyContainer<double> cont("04_PropertyContainer/testfiles/valid.cnf");
		
		cont.clear();
		
		// try to read from nonexsitent file
		TS_ASSERT_THROWS_ANYTHING(cont.read_data_from_file("asldjflasdjflkasfd"));
		cont.clear();
		
		// init from file and read valid file	
		TS_ASSERT_THROWS_NOTHING(cont.read_data_from_file((file = path + "test1.dat").c_str()));
		
		// test all entries
		TS_ASSERT_EQUALS(cont.get("wurstsalat"), 2.0);
		TS_ASSERT_EQUALS(cont.get("bohnensalat"), 1.5);								
	}
	
	void test_read_invalid_files() {
		PropertyContainer<double> cont("04_PropertyContainer/testfiles/valid.cnf");		
		TS_ASSERT_THROWS_ANYTHING(cont.read_data_from_file("04_PropertyContainer/testfiles/test2i.dat"));
		cont.clear();
		TS_ASSERT_THROWS_NOTHING(cont.read_data_from_file("04_PropertyContainer/testfiles/test3i.dat"));
		TS_ASSERT_EQUALS(cont.get("wurstsalat"), 4.9);
		TS_ASSERT_EQUALS(cont.get("bohnensalat"), 1.5);
		TS_ASSERT_EQUALS(cont.get("krautsalat"), 34.0);
	}
	
	
	void test_read_data_from_files_without_prior_init() {
		PropertyContainer<double> cont;
		PropertyContainer<double> cont2;
		// test valid
		try {
			cont.read_data_from_file("04_PropertyContainer/testfiles/test_init.dat");
		} catch (Exception* e) {
			TS_FAIL(e->get_reason());
		}
						
		TS_ASSERT_EQUALS(cont.get("wurstsalat"), 2.0);
		TS_ASSERT_EQUALS(cont.get("bohnensalat"), 1.5);								
		// test invalid 
		TS_ASSERT_THROWS_ANYTHING(cont2.read_data_from_file("04_PropertyContainer/testfiles/test_init_invalid.dat"));						
	}
	
	
	void test_with_conf_path() {
		char* old = getenv("TDKPCONFPATH");
		setenv("TDKPCONFPATH", "04_PropertyContainer/testfiles/", 1);
		TS_ASSERT_THROWS_NOTHING(PropertyContainer<double> cont("valid.cnf"));
		if(old != NULL) {
			setenv("TDKPCONFPATH", 	old, 1);
		}
			
	}
	
	void test_double_key() {

		// read doublette cnf file
		TS_ASSERT_THROWS_ANYTHING(PropertyContainer<double> cont("04_PropertyContainer/testfiles/doublette.cnf"));
		
		// read doublette dat file
		PropertyContainer<double> cont;
		TS_ASSERT_THROWS_ANYTHING(cont.read_data_from_file("04_PropertyContainer/testfiles/test_init_doublette.dat"));						
	}
				
};

