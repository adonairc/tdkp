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
#include "tdkp/geometry/ElementBoundary.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/io/InputParser.h"

using namespace tdkp;

#include "21_ElementBoundary/moc_grid_reader.h"



class ElementBoundaryTest : public CxxTest::TestSuite {
public:
	ElementBoundaryTest();
	void test_1d_1st_line_boundary();
	void test_1d_2nd_line_boundary();
	void test_2d_1st_triangle();
	void test_2d_2nd_triangle();
	void test_2d_1st_rectangle();
	void test_2d_2nd_rectangle();
	void test_3d_1st_tetrahedron();
	void test_3d_2nd_tetrahedron();

private:
	ofstream fout;
	InputParser parser;	
};

ElementBoundaryTest::ElementBoundaryTest() {
	Logger::get_instance()->set_level(LOG_INFO_DEVEL2);
	Logger::get_instance()->del_listener(&std::cout);
	fout.open("21_ElementBoundary.log");
	Logger::get_instance()->add_listener(&fout);	
}

/** test boundary integral of 1d line object */
void ElementBoundaryTest::test_1d_1st_line_boundary() {
	
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader1D(), 1);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)3);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)3);
		// prepare "integrals"
		geometry->get_element_boundary(0).prepare_integrals();
		geometry->get_element_boundary(1).prepare_integrals();
		geometry->get_element_boundary(2).prepare_integrals();
		
		// get center boundary
		ElementBoundary& bnd = geometry->get_element_boundary(1);
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)1);
		TS_ASSERT_EQUALS(bnd.get_boundary_integral_0th_order(0,0), 1.0);
						
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}
		
}

/** test boundary integral of 1d line object */
void ElementBoundaryTest::test_1d_2nd_line_boundary() {
	
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader1D(), 2);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)5);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)3);
		// prepare "integrals"
		geometry->get_element_boundary(0).prepare_integrals();
		geometry->get_element_boundary(1).prepare_integrals();
		geometry->get_element_boundary(2).prepare_integrals();
		
		// get center boundary
		ElementBoundary& bnd = geometry->get_element_boundary(1);
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)1);
		TS_ASSERT_EQUALS(bnd.get_boundary_integral_0th_order(0,0), 1.0);
		TS_ASSERT_EQUALS(bnd.get_node(0).get_index_vertex(), geometry->get_node(1).get_index_vertex());
																	
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}
		
}

/** test boundary integral of 2d triangles */
void ElementBoundaryTest::test_2d_1st_triangle() {
	
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader2DTriangle(), 1);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)4);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)5);
		// prepare "integrals"
		int inner_boundaries = 0;
		unsigned int bidx = 0;
		for(unsigned int ii = 0; ii < geometry->get_num_boundaries(); ii++) {		
			if(geometry->get_element_boundary(ii).get_num_elements() == 2) {
				inner_boundaries++;
				bidx = ii;	
			}
			geometry->get_element_boundary(ii).prepare_integrals();
		}
		TS_ASSERT_EQUALS(inner_boundaries, 1);
		ElementBoundary& bnd = geometry->get_element_boundary(bidx);
		
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)2);
					
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,0), 1.0 / 6.0, 1.0e-9);
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,1), 1.0 / 6.0, 1.0e-9);					
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,1), 1.0 / 12.0, 1.0e-9);
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,0), 1.0 / 12.0, 1.0e-9);		
					
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}		
}


/** test boundary integral of 2d triangles */
void ElementBoundaryTest::test_2d_2nd_triangle() {
	
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader2DTriangle(), 2);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)9);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)5);
		// prepare "integrals"
		int inner_boundaries = 0;
		unsigned int bidx = 0;
		for(unsigned int ii = 0; ii < geometry->get_num_boundaries(); ii++) {		
			if(geometry->get_element_boundary(ii).get_num_elements() == 2) {
				inner_boundaries++;
				bidx = ii;	
			}
			geometry->get_element_boundary(ii).prepare_integrals();
		}
		TS_ASSERT_EQUALS(inner_boundaries, 1);
		ElementBoundary& bnd = geometry->get_element_boundary(bidx);
		
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)3);
		
		for(unsigned int ii = 0; ii < 3; ii++) {
			for(unsigned int jj = ii + 1; jj < 3; jj++) {
				TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(ii,jj), bnd.get_boundary_integral_0th_order(jj,ii), 1.0e-9);	
			}	
		}	
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,0), 1.0 / 15.0, 1.0e-9);			
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,1), 1.0 / 30.0, 1.0e-9);
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,2), - 1.0 / 60.0, 1.0e-9);
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,1), 4.0 / 15.0, 1.0e-9);
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,2), 1.0 / 30.0, 1.0e-9);
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(2,2), 1.0 / 15.0, 1.0e-9);
							
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}		
}


/** test boundary integral of 2d rectangels */
void ElementBoundaryTest::test_2d_1st_rectangle() {
	
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader2DRectangle(), 1);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)6);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)7);
		// prepare "integrals"
		int inner_boundaries = 0;
		unsigned int bidx = 0;
		for(unsigned int ii = 0; ii < geometry->get_num_boundaries(); ii++) {		
			if(geometry->get_element_boundary(ii).get_num_elements() == 2) {
				inner_boundaries++;
				bidx = ii;	
			}
			geometry->get_element_boundary(ii).prepare_integrals();
		}
		TS_ASSERT_EQUALS(inner_boundaries, 1);
		ElementBoundary& bnd = geometry->get_element_boundary(bidx);
		
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)2);
		
		for(unsigned int ii = 0; ii < 2; ii++) {
			for(unsigned int jj = ii + 1; jj < 2; jj++) {
				TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(ii,jj), bnd.get_boundary_integral_0th_order(jj,ii), 1.0e-9);	
			}	
		}	
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,0), 1.0 / 5.0, 1.0e-9);			
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,1), 1.0 / 10.0, 1.0e-9);
		TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,1), 1.0 / 5.0, 1.0e-9);
							
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}		
}


/** test boundary integral of 2d rectangels */
void ElementBoundaryTest::test_2d_2nd_rectangle() {
	
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader2DRectangle(), 2);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)13);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)7);
		// prepare "integrals"
		int inner_boundaries = 0;
		unsigned int bidx = 0;
		for(unsigned int ii = 0; ii < geometry->get_num_boundaries(); ii++) {		
			if(geometry->get_element_boundary(ii).get_num_elements() == 2) {
				inner_boundaries++;
				bidx = ii;	
			}
			geometry->get_element_boundary(ii).prepare_integrals();
		}
		TS_ASSERT_EQUALS(inner_boundaries, 1);
		ElementBoundary& bnd = geometry->get_element_boundary(bidx);
		
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)3);
		
		for(unsigned int ii = 0; ii < 3; ii++) {
			for(unsigned int jj = ii + 1; jj < 3; jj++) {
				TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(ii,jj), bnd.get_boundary_integral_0th_order(jj,ii), 1.0e-9);	
			}	
		}	
		
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,0), 2.0 / 25.0, 1.0e-9);          
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,1), 1.0 / 25.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,2), - 1.0 / 50.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,1), 8.0 / 25.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,2), 1.0 / 25.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(2,2), 2.0 / 25.0, 1.0e-9);
							
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}		
}

/** test boundary integral of 3d tetrahedrons */
void ElementBoundaryTest::test_3d_1st_tetrahedron() {
	
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader3DTetrahedron(), 1);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)5);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)7);
		// prepare "integrals"
		int inner_boundaries = 0;
		unsigned int bidx = 0;
		for(unsigned int ii = 0; ii < geometry->get_num_boundaries(); ii++) {		
			if(geometry->get_element_boundary(ii).get_num_elements() == 2) {
				inner_boundaries++;
				bidx = ii;	
			}
			geometry->get_element_boundary(ii).prepare_integrals();
		}
		TS_ASSERT_EQUALS(inner_boundaries, 1);
		ElementBoundary& bnd = geometry->get_element_boundary(bidx);
		
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)3);
		
		for(unsigned int ii = 0; ii < 3; ii++) {
			for(unsigned int jj = ii + 1; jj < 3; jj++) {
				TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(ii,jj), bnd.get_boundary_integral_0th_order(jj,ii), 1.0e-9);	
			}	
		}	
		
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,0), 1.0 / 10.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,1), 1.0 / 20.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,2), 1.0 / 20.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,1), 1.0 / 10.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,2), 1.0 / 20.0, 1.0e-9);        
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(2,2), 1.0 / 10.0, 1.0e-9);          

							
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}		
}

void ElementBoundaryTest::test_3d_2nd_tetrahedron() {
	try {
		// setup mini grid
		Geometry* geometry = parser.read_geometry(MocGridReader3DTetrahedron(), 2);
		TS_ASSERT_EQUALS(geometry->get_num_nodes(), (unsigned int)14);
		TS_ASSERT_EQUALS(geometry->get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(geometry->get_num_boundaries(), (unsigned int)7);
		// prepare "integrals"
		int inner_boundaries = 0;
		unsigned int bidx = 0;
		for(unsigned int ii = 0; ii < geometry->get_num_boundaries(); ii++) {		
			if(geometry->get_element_boundary(ii).get_num_elements() == 2) {
				inner_boundaries++;
				bidx = ii;	
			}
			geometry->get_element_boundary(ii).prepare_integrals();
		}
		TS_ASSERT_EQUALS(inner_boundaries, 1);
		ElementBoundary& bnd = geometry->get_element_boundary(bidx);
		
		TS_ASSERT_EQUALS(bnd.get_num_elements(), (unsigned int)2);
		TS_ASSERT_EQUALS(bnd.get_num_nodes(), (unsigned int)6);
		
		for(unsigned int ii = 0; ii < 6; ii++) {
			//cout << "NODE " << ii << bnd.get_node(ii);
			for(unsigned int jj = ii + 1; jj < 6; jj++) {
				TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(ii,jj), bnd.get_boundary_integral_0th_order(jj,ii), 1.0e-9);	
			}	
		}	
		
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,0), 1.0 / 50.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,1), 0.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,2), - 1.0 / 300.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,3), - 1.0 / 75.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,4), - 1.0 / 300.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(0,5), 0.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,1), 8.0 / 75.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,2), 0.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,3), 4.0/75.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,4), -1.0/75.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(1,5), 4.0 / 75.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(2,2), 1.0 / 50.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(2,3), 0.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(2,4), - 1.0 / 300.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(2,5), -1.0 / 75.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(3,3), 8.0 / 75.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(3,4), 0.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(3,5), 4.0 / 75.0, 1.0e-9);                  
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(4,4), 1.0 / 50.0, 1.0e-9);
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(4,5), 0.00, 1.0e-9);               
        TS_ASSERT_DELTA(bnd.get_boundary_integral_0th_order(5,5), 8.0 / 75.0, 1.0e-9);
							
		delete geometry;
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}		
}
	
#endif /*DEFINITION_H_*/
