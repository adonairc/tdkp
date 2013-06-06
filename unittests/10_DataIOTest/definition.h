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
#include <list>
#include <map>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/main/EigenSolution.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/main/Fields.h"
//#include "dfdat_c.h"
#include "tdkp/io/InputParser.h"

using namespace tdkp;

class DataIOTest : public CxxTest::TestSuite {
public:	

	DataIOTest() {
		timeval time;
		gettimeofday(&time, NULL);
		srand48(time.tv_usec);
	}

	std::fstream fout;	
	void setUp() {
		fout.open("10_DataIOTest_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);
		Logger::get_instance()->del_listener(&std::cout);				
	}
	void tearDown() {
		fout.close();	
	}
	
	void test_std_node_data();
	void test_std_element_data();	
	void test_eigensolution();
	void test_eigensolution_set();
	
	
	EigenSolution<double>* random_solution(int length, int num_per_node);
	
};


EigenSolution<double>* DataIOTest::random_solution(int length, int num_per_node) {

	vector<double> data;
	data.resize(num_per_node * length);
	for(int ii = 0; ii < (signed)data.size(); ii++) {
		data[ii] = drand48(); 	
	}
	return new EigenSolution<double>(drand48(), &data[0], length, num_per_node);

}

void DataIOTest::test_std_node_data() {
	
	int num_data_per_node = 5;
	InputParser parser;

	const int lengths[] = {10, 444, 784, 5345};
	int num = 4;
	for(int nn = 0; nn < num; nn++) {								
		int length = lengths[nn];
		
		// test simple stuff	
		StdNodeData<double> data2;
		TS_ASSERT_THROWS_NOTHING(data2.set_length(length, num_data_per_node));
		TS_ASSERT(data2.get_length() == length);
		TS_ASSERT(data2.get_num_data_per_node() == num_data_per_node);
		for(int qq = 0; qq < num_data_per_node; qq++) {
			ostringstream sout;
			sout << "blabla" << qq;
			TS_ASSERT_THROWS_NOTHING(data2.set_identifier(qq, sout.str().c_str()));
			TS_ASSERT(sout.str().compare(data2.get_identifier(qq)) == 0);
		}	
		// test constructor
		StdNodeData<double> data3(num_data_per_node, length);
		TS_ASSERT(data3.get_length() == length);
		TS_ASSERT(data3.get_num_data_per_node() == num_data_per_node);
		for(int qq = 0; qq < num_data_per_node; qq++) {
			ostringstream sout;
			sout << "blabla" << qq;
			TS_ASSERT_THROWS_NOTHING(data3.set_identifier(qq, sout.str().c_str()));
			TS_ASSERT(sout.str().compare(data3.get_identifier(qq)) == 0);
		}	
		
		// create data 
		StdNodeData<double> data(num_data_per_node, length);
		TS_ASSERT(data.get_length() == length);
		TS_ASSERT(data.get_num_data_per_node() == 5);
		double tmp;
		for(int ii = 0; ii < length; ii++) {			
			for(int qq = 0; qq < num_data_per_node; qq++) {
				tmp = drand48();
				TS_ASSERT_THROWS_NOTHING(data.set_node_value(ii, qq, tmp));
				TS_ASSERT(data.get_node_value(ii,qq) == tmp);
			}	
		}
		// store data to file		
		ofstream fout("10_StdNodeData.dat", ios::binary);
		TS_ASSERT_THROWS_NOTHING(parser.write_binary(data, fout));
		fout.close();
		// now test reading
		// first, try to get exception when reading data of different type
	 	ifstream fin("10_StdNodeData.dat", ios::binary);
		StdNodeData<cplx> data_cplx;
		TS_ASSERT_THROWS_ANYTHING(parser.read_binary(data_cplx, fin));		
		// seek back and read
		StdNodeData<double> data_compare;
		fin.seekg(0);
		try {
			parser.read_binary(data_compare, fin);
		} catch(Exception* e) {
			TS_FAIL(e->get_reason());	
			return;
		}
		fin.close();
		TS_ASSERT(data_compare.get_length() == length);
		TS_ASSERT(data_compare.get_num_data_per_node() == 5);		
		for(int ii = 0; ii < length; ii++) {			
			for(int qq = 0; qq < num_data_per_node; qq++) {
				TS_ASSERT(data.get_node_value(ii,qq) == data_compare.get_node_value(ii,qq));
			}	
		}						
		remove("10_StdNodeData.dat");
	}			
}



void DataIOTest::test_std_element_data() {
	
	int num_data_per_element = 5;
	InputParser parser;

	const int lengths[] = {10, 444, 784, 5345};
	int num = 4;
	for(int nn = 0; nn < num; nn++) {								
		int length = lengths[nn];
		
		// test simple stuff	
		StdElementData<double> data2;
		TS_ASSERT_THROWS_NOTHING(data2.set_length(length, num_data_per_element));
		TS_ASSERT(data2.get_length() == length);
		TS_ASSERT(data2.get_num_data_per_element() == num_data_per_element);
		for(int qq = 0; qq < num_data_per_element; qq++) {
			ostringstream sout;
			sout << "blabla" << qq;
			TS_ASSERT_THROWS_NOTHING(data2.set_identifier(qq, sout.str().c_str()));
			TS_ASSERT(sout.str().compare(data2.get_identifier(qq)) == 0);
		}	
		// test constructor
		StdElementData<double> data3(num_data_per_element, length);
		TS_ASSERT(data3.get_length() == length);
		TS_ASSERT(data3.get_num_data_per_element() == num_data_per_element);
		for(int qq = 0; qq < num_data_per_element; qq++) {
			ostringstream sout;
			sout << "blabla" << qq;
			TS_ASSERT_THROWS_NOTHING(data3.set_identifier(qq, sout.str().c_str()));
			TS_ASSERT(sout.str().compare(data3.get_identifier(qq)) == 0);
		}	
		
		// create data 
		StdElementData<double> data(num_data_per_element, length);
		TS_ASSERT(data.get_length() == length);
		TS_ASSERT(data.get_num_data_per_element() == 5);
		double tmp;
		for(int ii = 0; ii < length; ii++) {			
			for(int qq = 0; qq < num_data_per_element; qq++) {
				tmp = drand48();
				TS_ASSERT_THROWS_NOTHING(data.set_element_value(ii, qq, tmp));
				TS_ASSERT(data.get_element_value(ii,qq) == tmp);
			}	
		}
		// store data to file		
		ofstream fout("10_StdElementData.dat", ios::binary);
		TS_ASSERT_THROWS_NOTHING(parser.write_binary(data, fout));
		fout.close();
		// now test reading
		// first, try to get exception when reading data of different type
	 	ifstream fin("10_StdElementData.dat", ios::binary);
		StdElementData<cplx> data_cplx;
		TS_ASSERT_THROWS_ANYTHING(parser.read_binary(data_cplx, fin));		
		// seek back and read
		StdElementData<double> data_compare;
		fin.seekg(0);
		try {
			parser.read_binary(data_compare, fin);
		} catch(Exception* e) {
			TS_FAIL(e->get_reason());	
			return;
		}
		fin.close();
		TS_ASSERT(data_compare.get_length() == length);
		TS_ASSERT(data_compare.get_num_data_per_element() == 5);		
		for(int ii = 0; ii < length; ii++) {			
			for(int qq = 0; qq < num_data_per_element; qq++) {
				TS_ASSERT(data.get_element_value(ii,qq) == data_compare.get_element_value(ii,qq));
			}	
		}			
		remove("10_StdElementData.dat");
	}			
}

 
void DataIOTest::test_eigensolution() {
	
	vector<double> data;
	
	int lengths[] = {100, 250, 666, 1000, 5512};
	int num = 2;

	for(int ll = 0; ll < num; ll++) {
		try {
			int num_per_node = ll + 1;	
			int length = lengths[ll];	
			data.resize(num_per_node * length);
			for(int ii = 0; ii < num_per_node * length; ii++) {
				data[ii] = drand48();
			}
			EigenSolution<double> solution(ll, &data[0], length, num_per_node);
			// first, test function compare
			EigenSolution<double> compare(ll, &data[0], length, num_per_node);			
			TS_ASSERT(solution.compare(compare));
			// test copy constructor
			EigenSolution<double> copy(solution);
			TS_ASSERT(solution.compare(copy));
			compare.set_node_value(0,0, 12232.1213);
			TS_ASSERT(!solution.compare(compare));
			TS_ASSERT(solution.size() == data.size()); 
			TS_ASSERT(solution.get_energy() == ll); 
			TS_ASSERT(solution.get_length() == length);
			TS_ASSERT_THROWS_NOTHING(solution.write("10_DataIOTest.dat"));
			EigenSolution<double> readin;
			TS_ASSERT_THROWS_NOTHING(readin.read("10_DataIOTest.dat"));
			TS_ASSERT(readin.size() == data.size()); 
			TS_ASSERT(readin.get_energy() == ll); 
			TS_ASSERT(readin.get_length() == length);			
			TS_ASSERT(readin.compare(solution));			
		} catch(Exception* e) {
			TS_FAIL(e->get_reason());	
		}		
	}
}

void DataIOTest::test_eigensolution_set() {
	int num = 6;
	int num_per_node = 4;
	int length = 1000;
	vector<EigenSolution<double>* > solutions;
	vector<EigenSolution<double>* > copies;

	EigenSolutionSet<double> set;

	// build eigensolutionms
	for(int nn = 0; nn < num; nn++) {
		solutions.push_back(this->random_solution(length, num_per_node));
		copies.push_back(new EigenSolution<double>(*(solutions[solutions.size() - 1])));
		// add to set
		TS_ASSERT_THROWS_NOTHING(set.add(solutions[solutions.size() - 1]));
	}
	TS_ASSERT(set.get_length() == solutions[0]->get_length());
	TS_ASSERT(set.get_num_solutions() == num);
	// compare
	for(int nn = 0; nn < num; nn++) {
		TS_ASSERT(set.get(nn).compare(*copies[nn]));	
	}
	// store to file
	TS_ASSERT_THROWS_NOTHING(set.write("eigensolutionset.dat"));
	// read new set from file
	EigenSolutionSet<double> set_read;
	TS_ASSERT_THROWS_NOTHING(set_read.read("eigensolutionset.dat"));
	// compare stored against existing
	for(int nn = 0; nn < num; nn++) {
		TS_ASSERT(set.get(nn).compare(set_read.get(nn)));	
	}	
	// try to read with wrong type
	EigenSolutionSet<cplx> set_wrong;
	TS_ASSERT_THROWS_ANYTHING(set_wrong.read("eigensolutionset.dat"));
	
	remove("eingesolutionset.dat");
	// that should be enough
			
}



#endif /*DEFINITION_H_*/
