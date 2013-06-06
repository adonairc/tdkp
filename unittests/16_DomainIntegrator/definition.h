
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
#include "tdkp/common/Domain.h"

using namespace tdkp;


class DomainIntegratorTest : public CxxTest::TestSuite {

public:			
	DomainIntegratorTest() { 
		fout.open("16_DomainIntegratorTest_output.log", ios::app | ios::out);		
		Logger::get_instance()->add_listener(&fout);
		//Logger::get_instance()->del_listener(&std::cout);
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);	
	}
	
	~DomainIntegratorTest() {
		Logger::get_instance()->del_listener(&fout);
		fout.close();	
	}
			
	std::fstream fout;
		
	void setUp() {		

	}
		
	void tearDown() {

	}		
	
	void test_rectangle_refine();
	void test_rectangle_domain();
	void test_radial2D_domain();
	void test_cuboid_domain();
	void test_line_domain();
	void test_domain_master();
	void test_domain_copy();
	void test_domain_master_serialize();
	
	void test_create_radial();
	void test_radial3D_domain_factory();
	void test_radial2D_domain_factory();
	void test_1D_domain_factory();  
		
};

void DomainIntegratorTest::test_rectangle_refine() {
	DomainNodeRectangle domain(0.0, 0.0, 3.0, 3.0);
	
	TS_ASSERT_DELTA(domain.get_weight(), 9.0, 1.0e-12);
	domain.refine();
	TS_ASSERT_DELTA(domain.get_weight(), 9.0, 1.0e-12);
	

	
	
}

void DomainIntegratorTest::test_rectangle_domain() {
	DomainNodeRectangle domain(-1.0, -1.0, 1.0, 1.0);
	
	// standard
	TS_ASSERT_DELTA(domain.get_weight(),4.0, 1.0e-12);
	TS_ASSERT(domain.leaf());
	TS_ASSERT_THROWS_NOTHING(domain.get_point());
	TS_ASSERT_DELTA(domain.get_weight(), domain.get_point().get_weight(), 1.0e-12);
	TS_ASSERT(domain.get_number_of_children() == 0);
	
	// refine
	domain.refine();
	TS_ASSERT(domain.get_number_of_children() == 9);
	TS_ASSERT_DELTA(domain.get_weight(),4.0, 1.0e-12);
	
	
	double weight = 0;
	for(unsigned ii = 0; ii < domain.get_number_of_children(); ii++) {
		TS_ASSERT(domain.get_child(ii).leaf());
		TS_ASSERT_THROWS_NOTHING(domain.get_child(ii).get_point());
		TS_ASSERT_DELTA(domain.get_child(ii).get_weight(), domain.get_child(ii).get_point().get_weight(), 1.0e-12);
		weight += domain.get_child(ii).get_weight();	
	}
	
	TS_ASSERT_DELTA(domain.get_weight(), weight, 1.0e-12);

}

void DomainIntegratorTest::test_1D_domain_factory() {

	DomainMaster domain;
	create_1D_domain_trapezoidal(domain, -22.0, 1.342, 100);
	double length = 0.0;
	for(unsigned int ii = 0; ii < domain.get_number_of_points(); ii++) {
		length += domain.get_point(ii).get_weight();	
	}
	TS_ASSERT_DELTA(length, 22.0 + 1.342, 1.0e-6);
}
	

void DomainIntegratorTest::test_radial3D_domain_factory() {

	DomainMaster domain;
	Vector3D direction(drand48(), drand48(), drand48());
	create_3D_domain_radial(domain, direction, 0.0, 5, 100);
	double length = 0.0;
	for(unsigned int ii = 0; ii < domain.get_number_of_points(); ii++) {
		length += domain.get_point(ii).get_weight();	
	}

	TS_ASSERT_DELTA(length, 4.0 / 3.0 * constants::pi * 5.0 * 5.0 * 5.0, 1.0e-4);
}

void DomainIntegratorTest::test_radial2D_domain_factory() {

	DomainMaster domain;
	create_2D_domain_radial(domain, 0.0, 5, 100);
	double length = 0.0;
	for(unsigned int ii = 0; ii < domain.get_number_of_points(); ii++) {
		length += domain.get_point(ii).get_weight();	
	}
	TS_ASSERT_DELTA(length,  constants::pi * 5.0 * 5.0, 1.0e-4);
	
}
	

void DomainIntegratorTest::test_radial2D_domain() {
	
	DomainNodeRadialPlane domain(1.0, 1.0, 0.0, 2.0);
	
	// area of circle is pi*r^2
	TS_ASSERT_DELTA(domain.get_weight(), constants::pi * 4.0, 1.0e-12);
	TS_ASSERT(domain.leaf());
	TS_ASSERT_THROWS_NOTHING(domain.get_point());
	TS_ASSERT_DELTA(domain.get_weight(), domain.get_point().get_weight(), 1.0e-12);
	TS_ASSERT(domain.get_number_of_children() == 0);
	
	// refine
	domain.refine();
	TS_ASSERT(domain.get_number_of_children() == 3);
	TS_ASSERT_DELTA(domain.get_weight(), constants::pi * 4.0, 1.0e-12);
		
	double weight = 0;
	for(unsigned ii = 0; ii < domain.get_number_of_children(); ii++) {
		TS_ASSERT(domain.get_child(ii).leaf());
		TS_ASSERT_THROWS_NOTHING(domain.get_child(ii).get_point());
		TS_ASSERT_DELTA(domain.get_child(ii).get_weight(), domain.get_child(ii).get_point().get_weight(), 1.0e-12);
		weight += domain.get_child(ii).get_weight();	
	}
	
	TS_ASSERT_DELTA(domain.get_weight(), weight, 1.0e-12);
 
}

void DomainIntegratorTest::test_line_domain() {
	DomainNodeLine domain(-1.0, 1.0);
	
	// standard
	TS_ASSERT_DELTA(domain.get_weight(),2.0, 1.0e-12);
	TS_ASSERT(domain.leaf());
	TS_ASSERT_THROWS_NOTHING(domain.get_point());
	TS_ASSERT_DELTA(domain.get_weight(), domain.get_point().get_weight(), 1.0e-12);
	TS_ASSERT(domain.get_number_of_children() == 0);
	
	// refine
	domain.refine();
	TS_ASSERT(domain.get_number_of_children() == 3);
	TS_ASSERT_DELTA(domain.get_weight(),2.0, 1.0e-12);
	
	
	double weight = 0;
	for(unsigned ii = 0; ii < domain.get_number_of_children(); ii++) {
		TS_ASSERT(domain.get_child(ii).leaf());
		TS_ASSERT_THROWS_NOTHING(domain.get_child(ii).get_point());
		TS_ASSERT_DELTA(domain.get_child(ii).get_weight(), domain.get_child(ii).get_point().get_weight(), 1.0e-12);
		weight += domain.get_child(ii).get_weight();	
	}
	
	TS_ASSERT_DELTA(domain.get_weight(), weight, 1.0e-12);

}

void DomainIntegratorTest::test_cuboid_domain() {
	DomainNodeCuboid domain(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
	
	// standard
	TS_ASSERT_DELTA(domain.get_weight(), 8.0, 1.0e-12);
	TS_ASSERT(domain.leaf());
	TS_ASSERT_THROWS_NOTHING(domain.get_point());
	TS_ASSERT_DELTA(domain.get_weight(), domain.get_point().get_weight(), 1.0e-12);
	TS_ASSERT(domain.get_number_of_children() == 0);
	
	// refine
	domain.refine();
	TS_ASSERT(domain.get_number_of_children() == 27);
	TS_ASSERT_DELTA(domain.get_weight(), 8.0, 1.0e-12);
		
	double weight = 0;
	for(unsigned ii = 0; ii < domain.get_number_of_children(); ii++) {
		TS_ASSERT(domain.get_child(ii).leaf());
		TS_ASSERT_THROWS_NOTHING(domain.get_child(ii).get_point());
		TS_ASSERT_DELTA(domain.get_child(ii).get_weight(), domain.get_child(ii).get_point().get_weight(), 1.0e-12);
		weight += domain.get_child(ii).get_weight();	
	}	
	TS_ASSERT_DELTA(domain.get_weight(), weight, 1.0e-12);
	TS_ASSERT_DELTA(domain.get_weight(), 8.0, 1.0e-12);		

}


void DomainIntegratorTest::test_domain_copy() {

	DomainMaster master(new DomainNodeRectangle(-1.0, -1.0, 1.0, 1.0));

	// ---------------------------------------------------
	// refine, make copy and then compare pointers!
	// ---------------------------------------------------
	for(unsigned int ii = 0; ii < 3; ii++) {

		master.refine();
		master.collapse();
		DomainMaster copy(master);			
		TS_ASSERT(copy.get_number_of_root_nodes() == master.get_number_of_root_nodes());
		TS_ASSERT(copy.get_number_of_points() == master.get_number_of_points());
		// checking pointers to points and root nodes
		for(unsigned int ii = 0; ii < master.get_number_of_points(); ii++) {
			for(unsigned int jj = 0; jj < copy.get_number_of_points(); jj++) {
				TS_ASSERT(&master.get_point(ii) != &copy.get_point(jj));	
			}	
		}
		// checking pointers to root nodes
		for(unsigned int ii = 0; ii < master.get_number_of_root_nodes(); ii++) {
			for(unsigned int jj = 0; jj < copy.get_number_of_root_nodes(); jj++) {
				TS_ASSERT(&master.get_root_node(ii) != &copy.get_root_node(jj));	
			}	
		}		
		
	} 
		
		
	
}

void DomainIntegratorTest::test_domain_master_serialize() {
	DomainMaster master(new DomainNodeRectangle(-1.0, -1.0, 1.0, 1.0));
	master.update();
	master.refine();
	master.refine();
	master.freeze();
		
	ofstream fout("domain_master_serialized.bin", ios::binary);
	master.write_binary(fout);
	fout.close();
	
	ifstream fin("domain_master_serialized.bin", ios::binary);
	DomainMaster slave(fin);
	
	// compare
	TS_ASSERT_DELTA(slave.get_total_weight(), master.get_total_weight(), 1.0e-12);
	TS_ASSERT(slave.get_dimension() == master.get_dimension());
	TS_ASSERT(slave.get_number_of_root_nodes() == master.get_number_of_root_nodes());
	TS_ASSERT(slave.get_number_of_points() == master.get_number_of_points());
	
	// compare points
	for(unsigned int ii = 0; ii < slave.get_number_of_points(); ii++) {
		TS_ASSERT_DELTA(slave.get_point(ii).get_weight(), master.get_point(ii).get_weight(), 1.0e-12);
		TS_ASSERT_DELTA(slave.get_point(ii).get_coord(0), master.get_point(ii).get_coord(0), 1.0e-12);
		TS_ASSERT_DELTA(slave.get_point(ii).get_coord(1), master.get_point(ii).get_coord(1), 1.0e-12);
		TS_ASSERT_DELTA(slave.get_point(ii).get_coord(2), master.get_point(ii).get_coord(2), 1.0e-12);
		TS_ASSERT(slave.get_point(ii).get_index() == master.get_point(ii).get_index());				  		
	}
		
}
	
void DomainIntegratorTest::test_domain_master() {

	DomainMaster master(new DomainNodeRectangle(-1.0, -1.0, 1.0, 1.0));
	
	TS_ASSERT_DELTA(master.get_total_weight(), 4.0, 1.0e-12);
	TS_ASSERT(master.get_number_of_root_nodes() == 1);
	TS_ASSERT(master.get_number_of_points() == 1);

	master.update();
	TS_ASSERT(master.get_number_of_root_nodes() == 1);
	TS_ASSERT(master.get_number_of_points() == 1);		
	master.refine();
	TS_ASSERT(master.get_number_of_root_nodes() == 1);
	TS_ASSERT(master.get_number_of_points() == 9);
	TS_ASSERT_DELTA(master.get_total_weight(), 4.0, 1.0e-12);
	master.collapse();
	TS_ASSERT(master.get_number_of_root_nodes() == 9);
	TS_ASSERT(master.get_number_of_points() == 9);
	TS_ASSERT_DELTA(master.get_total_weight(), 4.0, 1.0e-12);		
	
	
	double weight = 0.0;
	for(unsigned int ii = 0; ii < master.get_number_of_points(); ii++) {
		weight += master.get_point(ii).get_weight();
	}
	TS_ASSERT_DELTA(master.get_total_weight(), weight, 1.0e-10);	
	
	
}

void DomainIntegratorTest::test_create_radial() {
	DomainMaster master(new DomainNodeRadialPlane(1.0, 0.0, 0.0, 3.0));
	master.refine(); master.refine(); master.refine();
	master.update();
	
	DomainMaster empty;
	create_2D_domain_radial(empty, 0.0, 3.0, 10);
	TS_ASSERT(empty.get_number_of_points() == 10);
	TS_ASSERT_DELTA(empty.get_total_weight(), master.get_total_weight(), 1.0e-12);
		
}
