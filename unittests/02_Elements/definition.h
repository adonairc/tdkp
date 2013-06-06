#define TDKP_UNITTEST


// skip 3d tests (testing every element integral is numerically very expensive ....)
const bool skip_3d_tests = false;

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
#include <complex>
#include <map>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "tdkp/geometry/Node.h"
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/geometry/Region.h"
#include "tdkp/main/Fields.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/geometry/Element3DTR.h"
#include "tdkp/geometry/Element3DTetrahedron2nd.h"
#include "tdkp/geometry/Element2DRect.h"
#include "tdkp/geometry/Element2DRect2nd.h"
#include "tdkp/geometry/Element2DTriangle.h"
#include "tdkp/geometry/Element2DTriangle2nd.h"
#include "tdkp/geometry/Element2DTriangleHermiteCubic.h"
#include "tdkp/geometry/Element1DLine.h"
#include "tdkp/geometry/Element1DLine2nd.h"
#include "tdkp/geometry/Element1DHermiteCubic.h"
#include "tdkp/geometry/ElementBoundary.h"

using namespace tdkp;

#define ERROR_THRESHOLD 0.05

#include "validator.h"

class ElementTest : public CxxTest::TestSuite {
public:	
	std::fstream fout;	
	Region region;	
	ElementTest() : region("dummy") {
		fout.open("02_Elements_output.log", ios::app | ios::out);
		Logger::get_instance()->add_listener(&fout);
		Logger::get_instance()->del_listener(&std::cout);			
	}
	~ElementTest() {
		fout.close();	
	}

	static double test_integration(const double* coords) { return 1.0; }
	static double test_integration_nonsimple(const double* coords) {
		return coords[0] + 0.5 * coords[1] + 2.0 * coords[0] * coords[1]; 
	}		
	
	void set_additional_nodes(Element& elem);

		
	
	void test_2d_triangle_all_numerical();
	void test_1d_cubic_hermite_numerical();			
	void test_2d_cubic_hermite_triangle_numerical();			
			
	void test_1d_boundary();
	void test_2d_boundary();
	void test_3d_boundary();

	void test_2d_rectangle_element();	
	void test_2d_rectangle_element_nonsimple();
	void test_2d_triangle_element();	
	void test_2d_triangle_element_nonsimple();
	void test_basic_operation_3d_tetraeder();	
	void test_stretched_element_3d_tetraeder();
	void test_3d_2nd_order_tetraeder_numerical();				 
	void test_1d_line_element_simple();
	void test_1d_line_element_simple_but_other_place();
	void test_1d_line_element_stretched();
	

	void test_1d_line_all_numerical();
	void test_1d_2nd_line_all_numerical();
	
		
	void test_2d_2nd_order_global2local();
	
	void test_2d_rectangle_all_numerical();	
	void test_2d_2nd_order_rectangle_all_numerical();
	void test_2d_2nd_order_triangle_all_numerical();	
		
	
	void test_3d_tetraeder_all_numerical();
	void test_3d_tetraeder_all_numerical_nonsimple();
	
		
};

template<class T>
void test_if_basis_is_nodal(T& elem, const string& msg) {		
	for(unsigned int ii = 0; ii < elem.get_num_nodes(); ii++) {
		const Node& node = elem.get_node(ii);
		for(unsigned int jj = 0; jj < 6; jj++) {
			double val = elem.evaluate_form_function_global(jj, node.get_coord(0), node.get_coord(1), node.get_coord(2)); 
			if(ii == jj) {
				TSM_ASSERT_DELTA(msg, val, 1.0, 1.0e-8);
			} else {
				TSM_ASSERT_DELTA(msg, val, 0.0, 1.0e-8);
			}			
		}	
	}	
}


void ElementTest::test_2d_2nd_order_global2local() {
	try 
	{
		// create basic nodes
		Node vert0((unsigned)0, 0.0, 0.0);
		Node vert1((unsigned)1, 2.0, 0.0);
		Node vert2((unsigned)2, 0.0, 2.0);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		
		
		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		ElementBoundary bnd1(ElementBoundary::edge, 1, 2);
		bnd1.set_location('i');
		ElementBoundary bnd2(ElementBoundary::edge, 2, 2);
		bnd2.set_location('i');
		
		bnd0.set_corner_node(0, &vert0);
		bnd0.set_corner_node(1, &vert1);
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert2);
		bnd2.set_corner_node(0, &vert2);
		bnd2.set_corner_node(1, &vert0);
						
		// add nodes to element
		Element2DTriangle2nd elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_element_boundary(0, &bnd0);
		elem.set_element_boundary(1, &bnd1);
		elem.set_element_boundary(2, &bnd2);
		
		TS_ASSERT(elem.get_num_additional_nodes() == 3);		
		vector<unsigned int> locator_vertices;
		Element::AdditionalNodeLocation loc;
		vector<double> coords;
		unsigned int tag = 0;
		
		for(unsigned int ii = 0; ii < 3; ii++) {
			elem.get_additional_node_locator(ii, loc, locator_vertices, coords, tag);
			elem.set_additional_node(ii, new Node(3 + ii, coords[0], coords[1]));
		}		
		elem.prepare();
		
		double x = 0.0;
		double y = 0.0;
		while(x < 2.0) {
			y = 0.0;
			while(y < 2.0 - x) {
				double ag[2];
				double al[2];
				ag[0] = x;
				ag[1] = y;
				elem.global2local(ag,al);
				TS_ASSERT_DELTA(al[0], ag[0] / 2.0, 1.0e-8);
				TS_ASSERT_DELTA(al[1], ag[1] / 2.0, 1.0e-8);
				y += 0.1;
			}
			x += 0.1;
		}
		
		// --------------------------------------------------
		// check if form functions are correct labelled
		// --------------------------------------------------
		TS_ASSERT(elem.get_num_nodes() == 6);
		test_if_basis_is_nodal(elem, "test_2d_2nd_order_global2local"); 		
				
	} catch (Exception *e) {
		TS_FAIL(e->get_reason());	
	}	
}

void ElementTest::test_2d_2nd_order_rectangle_all_numerical() {
	
	try 
	{				
		// create basic nodes
		Node vert0((unsigned)0,  0.0,  0.0);
		Node vert1((unsigned)1,  1.0,  0.0);
		Node vert2((unsigned)2,  1.0,  1.0);
		Node vert3((unsigned)3,  0.0,  1.0);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);
				
		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		ElementBoundary bnd1(ElementBoundary::edge, 1, 2);
		bnd1.set_location('i');
		ElementBoundary bnd2(ElementBoundary::edge, 2, 2);
		bnd2.set_location('i');
		ElementBoundary bnd3(ElementBoundary::edge, 3, 2);
		bnd3.set_location('i');		
		
		bnd0.set_corner_node(0, &vert0);
		bnd0.set_corner_node(1, &vert1);
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert2);
		bnd2.set_corner_node(0, &vert2);
		bnd2.set_corner_node(1, &vert3);
		bnd3.set_corner_node(0, &vert3);
		bnd3.set_corner_node(1, &vert0);
								
		// add nodes to element
		Element2DRect2nd elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		elem.set_element_boundary(0, &bnd0);
		elem.set_element_boundary(1, &bnd1);
		elem.set_element_boundary(2, &bnd2);
		elem.set_element_boundary(3, &bnd3);
		
		TS_ASSERT(elem.get_num_additional_nodes() == 4);		
		vector<unsigned int> locator_vertices;
		Element::AdditionalNodeLocation loc;
		vector<double> coords;
		unsigned int tag = 0;
		
		for(unsigned int ii = 0; ii < 4; ii++) {
			elem.get_additional_node_locator(ii, loc, locator_vertices, coords, tag);
			TS_ASSERT_EQUALS(tag, (unsigned int)0);
			TS_ASSERT_EQUALS(loc, Element::edge_node);
			TS_ASSERT_EQUALS(locator_vertices.size(), (unsigned int)2);
			TS_ASSERT_EQUALS(locator_vertices[0], ii);
			TS_ASSERT_EQUALS(locator_vertices[1], (ii + 1) % 4);											
			if(ii == 0) {
				TS_ASSERT_EQUALS(coords[0], 0.5);
				TS_ASSERT_EQUALS(coords[1], 0.0);	
			} else if(ii == 1) {
				TS_ASSERT_EQUALS(coords[0], 1.0);
				TS_ASSERT_EQUALS(coords[1], 0.5);							
			} else if(ii == 2) {
				TS_ASSERT_EQUALS(coords[0], 0.5);
				TS_ASSERT_EQUALS(coords[1], 1.0);				
			} else if(ii == 3) {
				TS_ASSERT_EQUALS(coords[0], 0.0);
				TS_ASSERT_EQUALS(coords[1], 0.5);				
			}
			elem.set_additional_node(ii, new Node(4 + ii, coords[0], coords[1]));
		}
		
		elem.prepare();
		TS_ASSERT_EQUALS(elem.get_volume(), 1.0);
		
		// ------------------------------------
		// test nodal basis
		// ------------------------------------
		for(unsigned int ii = 0; ii < 8; ii++) {
			double gx = elem.get_node(ii).get_coord(0);
			double gy = elem.get_node(ii).get_coord(1);
			for(unsigned int jj = 0; jj < 8; jj++) {
				double tmp = elem.evaluate_form_function_global(jj, gx, gy, 0.0);
				if(ii == jj) {
					TS_ASSERT_DELTA(tmp, 1.0, 1.0e-8);
				} else {
					TS_ASSERT_DELTA(tmp, 0.0, 1.0e-8);
				}	
			}	
		}
				
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_2nd_order_rectangle_all_numerical 1");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());		
	}
	try 
	{
   
		// create nontrivial
		Node vert0((unsigned)0,  1.0,  1.0);
		Node vert1((unsigned)1,  1.809016994374947,  1.587785252292473);
		Node vert2((unsigned)2,  1.221231742082474,  2.396802246667421);
		Node vert3((unsigned)3,  0.412214747707527,  1.809016994374947);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);
				
		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		ElementBoundary bnd1(ElementBoundary::edge, 1, 2);
		bnd1.set_location('i');
		ElementBoundary bnd2(ElementBoundary::edge, 2, 2);
		bnd2.set_location('i');
		ElementBoundary bnd3(ElementBoundary::edge, 3, 2);
		bnd3.set_location('i');		
		
		bnd0.set_corner_node(0, &vert0);
		bnd0.set_corner_node(1, &vert1);
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert2);
		bnd2.set_corner_node(0, &vert2);
		bnd2.set_corner_node(1, &vert3);
		bnd3.set_corner_node(0, &vert3);
		bnd3.set_corner_node(1, &vert0);
								
		// add nodes to element
		Element2DRect2nd elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		elem.set_element_boundary(0, &bnd0);
		elem.set_element_boundary(1, &bnd1);
		elem.set_element_boundary(2, &bnd2);
		elem.set_element_boundary(3, &bnd3);
		
		TS_ASSERT(elem.get_num_additional_nodes() == 4);		
		vector<unsigned int> locator_vertices;
		Element::AdditionalNodeLocation loc;
		vector<double> coords;
		unsigned int tag = 0;
		
		for(unsigned int ii = 0; ii < 4; ii++) {
			elem.get_additional_node_locator(ii, loc, locator_vertices, coords, tag);
			TS_ASSERT_EQUALS(tag, (unsigned int)0);
			TS_ASSERT_EQUALS(loc, Element::edge_node);			
			TS_ASSERT_EQUALS(locator_vertices.size(), (unsigned int)2);
			TS_ASSERT_EQUALS(locator_vertices[0], ii);
			TS_ASSERT_EQUALS(locator_vertices[1], (ii + 1) % 4);
			elem.set_additional_node(ii, new Node(4 + ii, coords[0], coords[1]));
		}
		
		elem.prepare();
		TS_ASSERT_DELTA(elem.get_volume(), 1.0, 1.0e-4);
		
		// ------------------------------------
		// test nodal basis
		// ------------------------------------
		for(unsigned int ii = 0; ii < 8; ii++) {
			double gx = elem.get_node(ii).get_coord(0);
			double gy = elem.get_node(ii).get_coord(1);
			for(unsigned int jj = 0; jj < 8; jj++) {
				double tmp = elem.evaluate_form_function_global(jj, gx, gy, 0.0);
				if(ii == jj) {
					TS_ASSERT_DELTA(tmp, 1.0, 1.0e-8);
				} else {
					TS_ASSERT_DELTA(tmp, 0.0, 1.0e-8);
				}	
			}	
		}
				
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_2nd_order_rectangle_all_numerical 1");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());		
	}		
}

void ElementTest::test_2d_2nd_order_triangle_all_numerical() {
	try 
	{
		// create basic nodes
		Node vert0((unsigned)0, 0.0, 0.0);
		Node vert1((unsigned)1, 0.0, 1.0);
		Node vert2((unsigned)2, 1.0, 0.0);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		
		
		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		ElementBoundary bnd1(ElementBoundary::edge, 1, 2);
		bnd1.set_location('i');
		ElementBoundary bnd2(ElementBoundary::edge, 2, 2);
		bnd2.set_location('i');
		
		bnd0.set_corner_node(0, &vert0);
		bnd0.set_corner_node(1, &vert1);
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert2);
		bnd2.set_corner_node(0, &vert2);
		bnd2.set_corner_node(1, &vert0);
						
		// add nodes to element
		Element2DTriangle2nd elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_element_boundary(0, &bnd0);
		elem.set_element_boundary(1, &bnd1);
		elem.set_element_boundary(2, &bnd2);
		
		TS_ASSERT(elem.get_num_additional_nodes() == 3);		
		vector<unsigned int> locator_vertices;
		Element::AdditionalNodeLocation loc;
		vector<double> coords;
		unsigned int tag = 0;
		
		for(unsigned int ii = 0; ii < 3; ii++) {
			elem.get_additional_node_locator(ii, loc, locator_vertices, coords, tag);
			TS_ASSERT_EQUALS(tag, (unsigned int)0);
			TS_ASSERT_EQUALS(loc, Element::edge_node);		
			TS_ASSERT_EQUALS(locator_vertices.size(), (unsigned int)2);
			TS_ASSERT_EQUALS(locator_vertices[0], ii);
			TS_ASSERT_EQUALS(locator_vertices[1], (ii + 1) % 3);
			/*		
			if(ii == 0) {
				TS_ASSERT_EQUALS(coords[0], 0.5);
				TS_ASSERT_EQUALS(coords[1], 0.0);	
			} else if(ii == 1) {
				TS_ASSERT_EQUALS(coords[0], 0.5);
				TS_ASSERT_EQUALS(coords[1], 0.5);				
			}  else if(ii == 2) {
				TS_ASSERT_EQUALS(coords[0], 0.0);
				TS_ASSERT_EQUALS(coords[1], 0.5);				
			}*/
			elem.set_additional_node(ii, new Node(3 + ii, coords[0], coords[1]));
		}
		
		elem.prepare();
		elem.print();
		TS_ASSERT_EQUALS(elem.get_volume(), 1.0 / 2.0);
		test_if_basis_is_nodal(elem, string("test_2d_2nd_order_triangle_all_numerical 1"));				
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_2nd_order_triangle_all_numerical 1");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());		
	}
	try 
	{
		// create more complicated nodes
		Node vert0((unsigned)37, 3.125, -3.125);
		Node vert1((unsigned)24,  3.0, -4.0);
		Node vert2((unsigned)51,  2.34375, -3.28125);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);

		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		ElementBoundary bnd1(ElementBoundary::edge, 1, 2);
		bnd1.set_location('i');
		ElementBoundary bnd2(ElementBoundary::edge, 2, 2);
		bnd2.set_location('i');
		
		bnd0.set_corner_node(0, &vert0);
		bnd0.set_corner_node(1, &vert1);
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert2);
		bnd2.set_corner_node(0, &vert2);
		bnd2.set_corner_node(1, &vert0);
						
		// add nodes to element
		Element2DTriangle2nd elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_element_boundary(0, &bnd0);
		elem.set_element_boundary(1, &bnd1);
		elem.set_element_boundary(2, &bnd2);
		
		TS_ASSERT(elem.get_num_additional_nodes() == 3);		
		vector<unsigned int> locator_vertices;
		Element::AdditionalNodeLocation loc;
		vector<double> coords;
		unsigned int tag = 0;
		
		for(unsigned int ii = 0; ii < 3; ii++) {
			elem.get_additional_node_locator(ii, loc, locator_vertices, coords, tag);
			TS_ASSERT_EQUALS(tag, (unsigned int)0);
			TS_ASSERT_EQUALS(loc, Element::edge_node);
			TS_ASSERT_EQUALS(locator_vertices.size(), (unsigned int)2);
			TS_ASSERT_EQUALS(locator_vertices[0], ii);
			TS_ASSERT_EQUALS(locator_vertices[1], (ii + 1) % 3);
			elem.set_additional_node(ii, new Node(3 + ii, coords[0], coords[1]));
		}
		
		elem.prepare();
						
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_2nd_order_triangle_all_numerical 2");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());		
	}			
	try 
	{
		// create more complicated nodes
		Node vert0((unsigned)0, 0.6, 0.2);
		Node vert1((unsigned)1, 0.5, 0.6);
		Node vert2((unsigned)2, 0.1, 0.3);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);

		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		ElementBoundary bnd1(ElementBoundary::edge, 1, 2);
		bnd1.set_location('i');
		ElementBoundary bnd2(ElementBoundary::edge, 2, 2);
		bnd2.set_location('i');
		
		bnd0.set_corner_node(0, &vert0);
		bnd0.set_corner_node(1, &vert1);
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert2);
		bnd2.set_corner_node(0, &vert2);
		bnd2.set_corner_node(1, &vert0);
						
		// add nodes to element
		Element2DTriangle2nd elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_element_boundary(0, &bnd0);
		elem.set_element_boundary(1, &bnd1);
		elem.set_element_boundary(2, &bnd2);
		
		TS_ASSERT(elem.get_num_additional_nodes() == 3);		
		vector<unsigned int> locator_vertices;
		Element::AdditionalNodeLocation loc;
		vector<double> coords;
		unsigned int tag = 0;
		
		for(unsigned int ii = 0; ii < 3; ii++) {
			elem.get_additional_node_locator(ii, loc, locator_vertices, coords, tag);
			TS_ASSERT_EQUALS(tag, (unsigned int)0);
			TS_ASSERT_EQUALS(loc, Element::edge_node);		
			if(ii == 0) {
				TS_ASSERT_DELTA(coords[0], 0.55, 1.0e-8);
				TS_ASSERT_DELTA(coords[1], 0.4, 1.0e-8);	
			} else if(ii == 1) {
				TS_ASSERT_DELTA(coords[0], 0.3, 1.0e-8);
				TS_ASSERT_DELTA(coords[1], 0.45, 1.0e-8);				
			}  else if(ii == 2) {
				TS_ASSERT_DELTA(coords[0], 0.35, 1.0e-8);
				TS_ASSERT_DELTA(coords[1], 0.25, 1.0e-8);				
			}
			elem.set_additional_node(ii, new Node(3 + ii, coords[0], coords[1]));
		}
		
		elem.prepare();
						
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_2nd_order_triangle_all_numerical 3");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());		
	}		

	
}

void ElementTest::test_2d_cubic_hermite_triangle_numerical() {
		
	// simple test
	try {
		// create basic nodes
		Node vert0((unsigned)0, -0.1, 0.1023);
		Node vert1((unsigned)1, 1.54, 0.5);
		Node vert2((unsigned)2, 0.4, 1.23);
				
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);		
		
		// add nodes to element
		Element2DTriangleHermiteCubic elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		set_additional_nodes(elem);
		elem.prepare();
		
		ElementNumericalValidator validator(elem);
		validator.plot_shape_functions("2d_cubic_hermite_triangle.dat");
		validator.test_evaluate_integral("test_2d_cubic_hermite_order_triangle_all_numerical 1");
		
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}
	try {
		// same element, but orient nodes other way around!
		Node vert0((unsigned)0, -0.1, 0.1023);
		Node vert1((unsigned)1, 1.54, 0.5);
		Node vert2((unsigned)2, 0.4, 1.23);
				
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);		
		
		// add nodes to element
		Element2DTriangleHermiteCubic elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert2);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert0);
		set_additional_nodes(elem);
		elem.prepare();
		
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_cubic_hermite_order_triangle_all_numerical 2");		
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());	
	}		
}

void ElementTest::test_2d_triangle_all_numerical() {
	
	try 
	{
		// create basic nodes
		Node vert0((unsigned)0, 0.0, 0.0);
		Node vert1((unsigned)1, 1.0, 0.0);
		Node vert2((unsigned)2, 0.0, 1.0);
				
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		
		// add nodes to element
		Element2DTriangle elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.prepare();
		
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_triangle_all_numerical 1");
	} catch(Exception*e) {
		TS_FAIL(e->get_reason());		
	}		
	{
		// create basic nodes
		Node vert0((unsigned)0,  1.1, 1.05);
		Node vert1((unsigned)1, -0.1, 3.0);
		Node vert2((unsigned)2, -1.2,-1.1);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		
		// add nodes to element
		Element2DTriangle elem(0);
		elem.set_region(&region,0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.prepare();
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_triangle_all_numerical 2");		
	}
	
}

void ElementTest::test_2d_rectangle_all_numerical() {
	
	{
		// create basic nodes
		Node vert0((unsigned)0, -1.0, -1.0);
		Node vert1((unsigned)1,  1.0, -1.0);		
		Node vert2((unsigned)2,  1.0,  1.0);
		Node vert3((unsigned)3, -1.0,  1.0);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);		
								
		// add nodes to element
		Element2DRect elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		try {
			elem.prepare();
		} catch(Exception *e) {
			cout << e->get_reason();	
		}
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_rectangle_all_numerical 1");		
	}
	{		
		// create basic nodes	
		Node vert0((unsigned)0, 0.0, 0.0);
		Node vert1((unsigned)1, 0.8660, 0.5);
		Node vert2((unsigned)2, 0.3660, 1.3660);
		Node vert3((unsigned)3, -0.5, 0.8660);
						
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);	
		
		// add nodes to element
		Element2DRect elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		try {
			elem.prepare();
		} catch(Exception *e) {
			cout << e->get_reason();	
		}
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_rectangle_all_numerical 2");
	
	}
	{
		
		// create basic nodes		
		Node vert0((unsigned)0, 5.00000000000000, -5.00000000000000);
		Node vert1((unsigned)1, 4.00000000000000, -3.26794919243112);
		Node vert2((unsigned)2, 1.40192378864668, -4.76794919243112);
		Node vert3((unsigned)3, 2.40192378864668, -6.50000000000000);
			
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);					
		
		// add nodes to element
		Element2DRect elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		try {
			elem.prepare();
		} catch(Exception *e) {
			cout << e->get_reason();	
		}
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_2d_rectangle_all_numerical 3");
	
	}	
}

// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------
void ElementTest::test_1d_line_element_simple() {
	
	Node a((unsigned int)0, 0.0);
	Node b((unsigned int)1, 1.0);
	
	a.set_index_vertex(0);
	b.set_index_vertex(1);
	
	Element1DLine elem(0);
	elem.set_region(&region,0);
	elem.set_corner_node(0,&a);
	elem.set_corner_node(1,&b);
	elem.prepare(); 
	
	TS_ASSERT_DELTA(elem.get_volume(), 1.0, 1.0e-10);	
	
	const double expect_2nd_order[2][2] = {{1.0, -1.0}, {-1.0,  1.0}};
	const double expect_1st_order[2][2] = {{-0.5, -0.5}, { 0.5, 0.5}};
    const double expect_0th_order[2][2] = {{1.0/3.0, 1.0/6.0}, {1.0/6.0, 1.0/3.0}};
	
	for(short ii = 0; ii < 2; ii++) {
		for(short jj = 0; jj < 2; jj++) {
			TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 0, jj), expect_2nd_order[ii][jj], 1.0e-12);
			TS_ASSERT_DELTA(elem.get_element_integral_1st_order(0, ii, jj),     expect_1st_order[ii][jj], 1.0e-12);
			TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii, jj), expect_0th_order[ii][jj], 1.0e-12);
		}	
	}				
}	

// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------
void ElementTest::test_1d_line_element_simple_but_other_place() {
	Node a((unsigned int)0, -2.0);
	Node b((unsigned int)1, -1.0);
	
	a.set_index_vertex(0);
	b.set_index_vertex(1);
	
	Element1DLine elem(0);
	elem.set_region(&region,0);
	elem.set_corner_node(0,&a);
	elem.set_corner_node(1,&b);
	elem.prepare(); 
	
	TS_ASSERT_DELTA(elem.get_volume(), 1.0, 1.0e-10);	
	
	const double expect_2nd_order[2][2] = {{1.0, -1.0}, {-1.0,  1.0}};
	const double expect_1st_order[2][2] = {{-0.5, -0.5}, {0.5, 0.5}};
    const double expect_0th_order[2][2] = {{1.0/3.0, 1.0/6.0}, {1.0/6.0, 1.0/3.0}};
	
	for(short ii = 0; ii < 2; ii++) {
		for(short jj = 0; jj < 2; jj++) {
			TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 0, jj), expect_2nd_order[ii][jj], 1.0e-12);
			TS_ASSERT_DELTA(elem.get_element_integral_1st_order(0, ii, jj),     expect_1st_order[ii][jj], 1.0e-12);
			TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii, jj), expect_0th_order[ii][jj], 1.0e-12);
		}	
	}	
}

// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------
void ElementTest::test_1d_line_element_stretched() {
	Node a((unsigned int)0, -1.0);
	Node b((unsigned int)1,  1.5);
	
	double h = 2.5;
	
	a.set_index_vertex(0);
	b.set_index_vertex(1);
	
	Element1DLine elem(0);
	elem.set_region(&region,0);
	elem.set_corner_node(0,&a);
	elem.set_corner_node(1,&b);
	elem.prepare(); 
	
	TS_ASSERT_DELTA(elem.get_volume(), h, 1.0e-10);	
	
	const double expect_2nd_order[2][2] = {{1.0, -1.0}, {-1.0,  1.0}};
	const double expect_1st_order[2][2] = {{-0.5, -0.5}, { 0.5, 0.5}};
    const double expect_0th_order[2][2] = {{1.0/3.0, 1.0/6.0}, {1.0/6.0, 1.0/3.0}};
	
	for(short ii = 0; ii < 2; ii++) {
		for(short jj = 0; jj < 2; jj++) {
			TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 0, jj), (1.0 / h) * expect_2nd_order[ii][jj], 1.0e-12);
			TS_ASSERT_DELTA(elem.get_element_integral_1st_order(0, ii, jj),     expect_1st_order[ii][jj], 1.0e-12);
			TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii, jj), h * expect_0th_order[ii][jj], 1.0e-12);
		}	
	}	
}

// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------	
	void ElementTest::test_stretched_element_3d_tetraeder() {
		double x[] = {0.0, 0.0, 1.5,  0.3, -0.1};
		double y[] = {0.0, 0.0, 0.1,  0.1, -2.0};
		double z[] = {0.0, 0.0, -0.1, 1.33, -0.1};
						
		// add nodes to element
		Element3DTR elem(0);
		elem.set_region(&region,0);

		Node* pvert;
		for(int ii = 0; ii < 4; ii++) {
			pvert = new Node(ii, x[ii + 1], y[ii + 1], z[ii + 1]);
			pvert->set_index_vertex(ii);
			elem.set_corner_node(ii, pvert);
		}
		elem.prepare();
		
		// check volume
		double volume = - x[3] * y[2] * z[1] + x[4] * y[2] * z[1] 
		                + x[2] * y[3] * z[1] - x[4] * y[3] * z[1] 
		                - x[2] * y[4] * z[1] + x[3] * y[4] * z[1] 
		                + x[3] * y[1] * z[2] - x[4] * y[1] * z[2] 
		                - x[1] * y[3] * z[2] + x[4] * y[3] * z[2] 
		                + x[1] * y[4] * z[2] - x[3] * y[4] * z[2] 
		                - x[2] * y[1] * z[3] + x[4] * y[1] * z[3] 
		                + x[1] * y[2] * z[3] - x[4] * y[2] * z[3] 
		                - x[1] * y[4] * z[3] + x[2] * y[4] * z[3] 
		                + x[2] * y[1] * z[4] - x[3] * y[1] * z[4] 
		                - x[1] * y[2] * z[4] + x[3] * y[2] * z[4] 
		                + x[1] * y[3] * z[4] - x[2] * y[3] * z[4];
		                
		volume = 1.0 / 6.0 * fabs(volume);
		TS_ASSERT_DELTA(elem.get_volume(), volume, 1.0e-10);		                
		                
		
		TS_ASSERT_DELTA(elem.get_jacobi_det(), 4.0237, 1.0e-8);
		
		double expect_0th_ord[4][4] = {
			{0.0670617, 0.0335308, 0.0335308, 0.0335308},
			{0.0335308, 0.0670617, 0.0335308, 0.0335308},
			{0.0335308, 0.0335308, 0.0670617, 0.0335308},
			{0.0335308, 0.0335308, 0.0335308, 0.0670617}
		};
		
		double expect_2nd_ord_dx[4][4] = {
			{0.373537, -0.329628, -0.0261215, -0.0177875},
			{-0.329628, 0.290881, 0.0230509, 0.0156966},
			{-0.0261215, 0.0230509, 0.00182668, 0.00124388},
			{-0.0177875, 0.0156966, 0.00124388, 0.000847023}			
		};
		
		double expect_2nd_ord_dy[4][4] = {
			{0.216838, -0.0097615, -0.0151635, -0.191913},
			{-0.0097615, 0.000439438, 0.000682622, 0.00863944},
			{-0.0151635, 0.000682622, 0.00106038, 0.0134205},
			{-0.191913, 0.00863944, 0.0134205, 0.169853}
		};
				
		double expect_2nd_ord_dz[4][4] = {
			{0.263041, 0.0615851, -0.312101, -0.0125258},
			{0.0615851, 0.0144187, -0.0730712, -0.00293262},
			{-0.312101, -0.0730712, 0.37031, 0.0148619},
			{-0.0125258, -0.00293262, 0.0148619, 0.000596466}
		};
		
		double expect_2nd_ord_dxdy[4][4] = {
			{-0.2846, 0.012812, 0.0199021, 0.251886},
			{0.251145, -0.0113059, -0.0175626, -0.222277}, 
			{0.0199021, -0.000895942, -0.00139175, -0.0176144}, 
			{0.0135524, -0.000610094, -0.000947718, -0.0119946}
		};

		double expect_2nd_ord_dxdz[4][4] = {
			{0.313458, 0.0733889, -0.37192, -0.0149266}, 
			{-0.276611, -0.0647621, 0.328201, 0.013172}, 
			{-0.0219201, -0.00513209, 0.0260084, 0.00104382},
			{-0.0149266, -0.00349471, 0.0177105, 0.000710789}
		};
		
		double expect_2nd_ord_dydz[4][4] = {
			{-0.238825, -0.0559154, 0.283368, 0.0113726}, 
			{0.0107513, 0.00251717, -0.0127565, -0.000511967}, 
			{0.016701, 0.00391017, -0.0198159, -0.000795288}, 
			{0.211373, 0.049488, -0.250795, -0.0100654}
		};
		
		double expect_1st_ord_dx[4][4] = {
			{-0.125125, -0.125125,-0.125125,-0.125125},
			{ 0.110417,  0.110417, 0.110417, 0.110417},
			{ 0.00875,   0.00875,  0.00875,  0.00875},
			{ 0.00595833,0.00595833, 0.00595833, 0.00595833}
		};
		
		double expect_1st_ord_dy[4][4] = {
			{0.0953333, 0.0953333, 0.0953333, 0.0953333},	
			{-0.00429167, -0.00429167, -0.00429167, -0.00429167},
			{-0.00666667, -0.00666667, -0.00666667, -0.00666667},
			{-0.084375, -0.084375, -0.084375, -0.084375}						
		};

		double expect_1st_ord_dz[4][4] = {
			{-0.105, -0.105, -0.105, -0.105},
			{ -0.0245833, -0.0245833,  -0.0245833, -0.0245833},
			{0.124583, 0.124583, 0.124583, 0.124583},						
			{0.005, 0.005, 0.005, 0.005},			
		};		
		
		for(short ii = 0; ii < 4; ii++) {
			for(short jj = 0; jj < 4; jj++) {
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 0, jj), expect_2nd_ord_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(1, ii, 1, jj), expect_2nd_ord_dy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(2, ii, 2, jj), expect_2nd_ord_dz[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii,jj), expect_0th_ord[ii][jj], 1.0e-6);			
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(0, ii, jj), expect_1st_ord_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(1, ii, jj), expect_1st_ord_dy[ii][jj], 1.0e-6);				
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(2, ii, jj), expect_1st_ord_dz[ii][jj], 1.0e-6);				
				
				
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 1, jj), expect_2nd_ord_dxdy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 2, jj), expect_2nd_ord_dxdz[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(1, ii, 2, jj), expect_2nd_ord_dydz[ii][jj], 1.0e-6);				
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 1, jj), elem.get_element_integral_2nd_order(1, jj, 0, ii), 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(0, ii, 2, jj), elem.get_element_integral_2nd_order(2, jj, 0, ii), 1.0e-6);				
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(2, ii, 1, jj), elem.get_element_integral_2nd_order(1, jj, 2, ii), 1.0e-6);				
			}
		}
		
		
		
		// 0.167654 P0 + 0.167654 P1 + 0.167654 P2 + 0.167654 P3						
		double P[4] = {12.0, 43.3, -4.23124, 0.000001};
		TS_ASSERT_DELTA(elem.integrate_solution(P) , 0.167654 * P[0] + 0.167654 * P[1] + 0.167654 * P[2] + 0.167654 * P[3], 1.0e-4);
		complex<double> D[4] = {complex<double>(12.0, 0), complex<double>(43.3, 0.0), complex<double>(-4.23124, 0.0), complex<double>(0.000001, 0.0)};
		TS_ASSERT_DELTA(elem.integrate_solution(D).real() ,0.167654 * P[0] + 0.167654 * P[1] + 0.167654 * P[2] + 0.167654 * P[3], 1.0e-4);		
								
								
		// ------------------------------------------------------------------
		// single integral
		// ------------------------------------------------------------------
		double expect_1st_single[3][4] = {
			{-0.5005, 0.441667, 0.035, 0.0238333},
			{0.381333, -0.0171667, -0.0266667, -0.3375},
			{-0.42, -0.0983333,  0.498333, 0.02}
		};								
		
		for(int dd = 0; dd < 3; dd++) {
			for(int ii = 0; ii < 4; ii++) {
				TS_ASSERT_DELTA(elem.get_single_integral_1st_order(dd,ii),expect_1st_single[dd][ii], 1.0e-5);	
			}	
		}
	}	
// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------	
	void ElementTest::test_2d_rectangle_element() {
		// create basic nodes
		Node vert0((unsigned)0, 0.0, 0.0);
		Node vert1((unsigned)1, 1.0, 0.0);
		Node vert2((unsigned)2, 1.0, 1.0);
		Node vert3((unsigned)3, 0.0, 1.0);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);					
		
		// add nodes to element
		Element2DRect elem(0);
		elem.set_region(&region,0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		elem.prepare();
		
		// check element volume (V == 1)
		TS_ASSERT_DELTA(elem.get_volume(), 1.0, 1.0e-9);
		
		// check jacobi (reference element is [-1,1])
		TS_ASSERT_DELTA(elem.get_jacobi_det(), 1.0/4.0, 1.0e-9);
		
		// zero order
		for(short ii = 0; ii < 4; ii++) {
			for(short jj = 0; jj < 4; jj++) {
				short dd = (short)abs(ii - jj);
				double val = 0.0;
				switch(dd) {
					case 0:
						val = 0.1111111;
						break;
					case 1: case 3:
						val = 0.0555556;
						break;
					case 2:
						val = 0.0277778;
				}
				TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii,jj),val,1.0e-6); 					
			}	
		}	
		
		const double second_order_dx_dx[4][4] = {
			{0.333333, -0.333333, -0.166667, 0.166667},
			{-0.333333, 0.333333, 0.166667, -0.166667},
			{-0.166667, 0.166667, 0.333333, -0.333333},
			{0.166667, -0.166667, -0.333333, 0.333333} 			
		};
		const double second_order_dy_dy[4][4] = {
			{0.333333, 0.166667, -0.166667, -0.333333},
			{0.166667, 0.333333, -0.333333, -0.166667},
			{-0.166667, -0.333333, 0.333333, 0.166667},
			{-0.333333, -0.166667, 0.166667, 0.333333} 			
		};
		const double second_order_dx_dy[4][4] = {
			{0.25, 0.25, -0.25, -0.25},
			{-0.25, -0.25, 0.25, 0.25},
			{-0.25, -0.25, 0.25, 0.25},
			{0.25, 0.25, -0.25, -0.25}
		};
		const double second_order_dy_dx[4][4] = {
			{0.25, -0.25, -0.25, 0.25},
			{0.25, -0.25, -0.25, 0.25},
			{-0.25, 0.25, 0.25, -0.25},
			{-0.25, 0.25, 0.25, -0.25}
		};
		
		const double first_order_dx[4][4] = {
			{-1.0/6.0, -1.0/6.0, -1.0/12.0, -1.0/12.0},
			{1.0/6.0, 1.0/6.0, 1.0/12.0, 1.0/12.0},
			{1.0/12.0, 1.0/12.0, 1.0/6.0, 1.0/6.0},
			{-1.0/12.0, -1.0/12.0, -1.0/6.0, -1.0/6.0}
		};
		
		const double first_order_dy[4][4] = {
			{-1.0/6.0, -1.0/12.0, -1.0/12.0, -1.0/6.0},
			{-1.0/12.0,-1.0/6.0,-1.0/6.0,-1.0/12.0},
			{1.0/12.0, 1.0/6.0, 1.0/6.0, 1.0/12.0},
			{1.0/6.0, 1.0/12.0, 1.0/12.0, 1.0/6.0}
		};		
		
		for(short ii = 0; ii < 4; ii++) {
			for(short jj = 0; jj < 4; jj++) {
				short dd = (short)abs(ii - jj);
				double val = 0.0;
				switch(dd) {
					case 0:
						val = 0.1111111;
						break;
					case 1: case 3:
						val = 0.0555556;
						break;
					case 2:
						val = 0.0277778;
				}
				TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii,jj),val,1.0e-6); 					
			}	
		}			
		
		for(short ii = 0; ii < 4; ii++) {
			for(short jj = 0; jj < 4; jj++) {
	//			cout << " TESTING " << ii << " " << jj << endl;
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DX,ii,jj), first_order_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DY,ii,jj), first_order_dy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DX,jj), second_order_dx_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DY,jj), second_order_dx_dy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DX,jj), second_order_dy_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DY,jj), second_order_dy_dy[ii][jj], 1.0e-6);												
			}
		}				
	}
	
	
// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------	
	void ElementTest::test_2d_rectangle_element_nonsimple() {
		Node vert0((unsigned)0, 5.00000000000000, -5.00000000000000);
		Node vert1((unsigned)1, 4.00000000000000, -3.26794919243112);
		Node vert2((unsigned)2, 1.40192378864668, -4.76794919243112);
		Node vert3((unsigned)3, 2.40192378864668, -6.50000000000000);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);					
		
		// add nodes to element
		Element2DRect elem(0);
		elem.set_region(&region,0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		elem.prepare();		

		// check element volume (V == 1)
		TS_ASSERT_DELTA(elem.get_volume(), 6.0, 1.0e-9);
		
		// check jacobi (reference element is [-1,1])
		TS_ASSERT_DELTA(elem.get_jacobi_det(), 1.5, 1.0e-9);

		// mathematica nintegrate data
		const double zero_order[4][4] = {
			{2.0/3.0, 1.0/3.0, 1.0/6.0, 1.0/3.0},
			{1.0/3.0, 2.0/3.0, 1.0/3.0, 1.0/6.0},
			{1.0/6.0, 1.0/3.0, 2.0/3.0, 1.0/3.0},
			{1.0/3.0, 1.0/6.0, 1.0/3.0, 2.0/3.0},									
		};
		
		const double second_order_dx_dx[4][4] = {
			{0.508173, -0.0416667, -0.36234, -0.104167},
			{-0.0416667, 0.0751603, -0.104167, 0.070673},
			{-0.36234, -0.104167, 0.508173, -0.0416667},
			{-0.104167, 0.070673, -0.0416667, 0.0751603} 
		};
		const double second_order_dx_dy[4][4] = {
			{-0.245281, 0.514619, 0.185141, -0.454478},
			{0.0146189, 0.00471869, 0.0455218,-0.0648593},
			{0.185141, -0.454478, -0.245281, 0.514619},  
			{0.0455218, -0.0648593, 0.0146189, 0.00471869}
		};
		const double second_order_dy_dx[4][4] = {
			{-0.245281, 0.0146189, 0.185141, 0.0455218},
			{0.514619,  0.00471869,-0.454478,-0.0648593},
			{0.185141, 0.0455218,  -0.245281, 0.0146189},
			{-0.454478, -0.0648593, 0.514619, 0.00471869}	
		};
		const double second_order_dy_dy[4][4] = {
			{0.214049, -0.347222, 0.00122857, 0.131944},
			{-0.347222, 0.647062, 0.131944, -0.431784},
			{0.00122857,0.131944, 0.214049, -0.347222},
			{0.131944, -0.431784, -0.347222, 0.647062}     	
		};
		const double first_order_dx[4][4] = {
			{0.538675,0.394338, 0.269338, 0.413675},
			{-0.105662,0.0386751,0.163675,0.0193376},
			{-0.269338,-0.413675, -0.538675, -0.394338},
			{-0.163675, -0.0193376, 0.105662, -0.0386751}     
		};	
		const double first_order_dy[4][4] = {
			{-0.266346, -0.349679, -0.133173, -0.0498397},
			{0.516346, 0.599679,0.383173,0.29984},
			{0.133173, 0.0498397, 0.266346, 0.349679},
			{-0.383173,-0.29984,-0.516346,-0.599679}			 
		};	
		
		for(short ii = 0; ii < 4; ii++) {
			for(short jj = 0; jj < 4; jj++) {
//				cout << " TESTING " << ii << " " << jj << endl;
				TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii,jj),	zero_order[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DX,ii,jj), first_order_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DY,ii,jj), first_order_dy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DX,jj), second_order_dx_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DY,jj), second_order_dx_dy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DX,jj), second_order_dy_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DY,jj), second_order_dy_dy[ii][jj], 1.0e-6);															}
		}
	}
// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------
	void ElementTest::test_2d_triangle_element() {
		// create basic nodes
		Node vert0((unsigned)0, 0.0, 0.0);
		Node vert1((unsigned)1, 1.0, 0.0);
		Node vert2((unsigned)2, 0.0, 1.0);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		
		// add nodes to element
		Element2DTriangle elem(0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_region(&region,0);
		elem.prepare();
		
		// check element volume (V == 1)
		TS_ASSERT_DELTA(elem.get_volume(), 0.5, 1.0e-9);
		
		// check jacobi (reference element is [-1,1])
		TS_ASSERT_DELTA(elem.get_jacobi_det(), 1.0, 1.0e-9);
		
		double val;
		for(short ii = 0; ii < 3; ii++) {
			for(short jj = 0; jj < 3; jj++) {
				if(ii == jj) {
					val = 1.0/12.0;	
				} else {
					val = 1.0/24.0;	
				}
				TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii,jj),val,1.0e-12);	
			}	
		}
		
		
		const double second_order_dx_dx[3][3] = {
			{0.5, -0.5, 0},
			{-0.5, 0.5, 0},
			{0, 0, 0}
		};		
		const double second_order_dx_dy[3][3] = {
			{0.5, 0, -0.5},
			{-0.5, 0, 0.5},
			{0, 0, 0}	
		};
		const double second_order_dy_dx[3][3] = {
			{0.5, -0.5, 0},
			{0.0, 0.0, 0.0},
			{-0.5, 0.5, 0.0}	
		};
		const double second_order_dy_dy[3][3] = {
			{0.5, 0, -0.5}, 
			{0, 0, 0},
			{-0.5, 0, 0.5}
		};
		const double first_order_dx[3][3] = {
			{-1.0/6.0, -1.0/6.0, -1.0/6.0},
			{1.0/6.0, 1.0/6.0, 1.0/6.0},				
			{0, 0, 0}
		};
		const double first_order_dy[3][3] = {
			{-1.0/6.0, -1.0/6.0, -1.0/6.0},			
			{0, 0, 0},
			{1.0/6.0, 1.0/6.0, 1.0/6.0}
		};		

		for(short ii = 0; ii < 3; ii++) {
			for(short jj = 0; jj < 3; jj++) {
	//			cout << " TESTING " << ii << " " << jj << endl;
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DX,ii,jj), first_order_dx[ii][jj], 1.0e-8);
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DY,ii,jj), first_order_dy[ii][jj], 1.0e-8);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DX,jj), second_order_dx_dx[ii][jj], 1.0e-8);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DY,jj), second_order_dx_dy[ii][jj], 1.0e-8);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DX,jj), second_order_dy_dx[ii][jj], 1.0e-8);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DY,jj), second_order_dy_dy[ii][jj], 1.0e-8);												
			}
		}									
	}	
// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------	
	void ElementTest::test_2d_triangle_element_nonsimple() {
		// create basic nodes
		Node vert0((unsigned)0,  1.1, 1.05);
		Node vert1((unsigned)1, -0.1, 3.0);
		Node vert2((unsigned)2, -1.2,-1.1);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);		
		// add nodes to element
		Element2DTriangle elem(0);
		elem.set_region(&region,0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.prepare();
		
		// check element volume (V == 1)
		TS_ASSERT_DELTA(elem.get_volume(), 3.5325, 1.0e-6);
		
		// check jacobi (reference element is [-1,1])
		TS_ASSERT_DELTA(elem.get_jacobi_det(), 7.065, 1.0e-6);
		
		const double zero_order[3][3] = {
			{0.58875, 0.294375, 0.294375},
			{0.294375, 0.58875, 0.294375},
			{0.294375, 0.294375, 0.58875}			
		};
				
		const double second_order_dx_dx[3][3] = {
			{1.189667, -0.62385, -0.565817},
			{-0.62385, 0.327141, 0.296709},
			{-0.565817, 0.296709, 0.269108}
		};		
		const double second_order_dx_dy[3][3] = {
			{-0.319179, 0.667374, -0.348195},
			{0.167374, -0.349965, 0.18259},
			{0.151805, -0.31741, 0.165605}	
		};
		const double second_order_dy_dx[3][3] = {
			{-0.319179, 0.167374, 0.151805},
			{0.667374, -0.349965, -0.31741},
			{-0.348195, 0.18259, 0.165605}	
		};
		const double second_order_dy_dy[3][3] = {
			{0.0856334, -0.179052, 0.0934183}, 
			{-0.179052, 0.374381, -0.195329},
			{0.0934183, -0.195329, 0.101911}
		};
		const double first_order_dx[3][3] = {
			{0.683333, 0.683333, 0.683333},
			{-0.358333, -0.358333, -0.358333},				
			{-0.325, -0.325, -0.325}
		};
		const double first_order_dy[3][3] = {
			{-0.183333, -0.183333, -0.183333},			
			{0.383333, 0.383333,0.383333},
			{-0.2, -0.2, -0.2}
		};		

		for(short ii = 0; ii < 3; ii++) {
			for(short jj = 0; jj < 3; jj++) {
		//		cout << " TESTING " << ii << " " << jj << endl;
				TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii,jj),	zero_order[ii][jj], 1.0e-6);	
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DX,ii,jj), first_order_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_1st_order(D_DY,ii,jj), first_order_dy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DX,jj), second_order_dx_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DX,ii,D_DY,jj), second_order_dx_dy[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DX,jj), second_order_dy_dx[ii][jj], 1.0e-6);
				TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(D_DY,ii,D_DY,jj), second_order_dy_dy[ii][jj], 1.0e-6);												
			}
		}			
						
	}
		



// ----------------------------------------------------------
// implementation 		
// ----------------------------------------------------------		
	void ElementTest::test_basic_operation_3d_tetraeder() {
		// create basic nodes
		Node vert0(0, 0.0, 0.0, 0.0);
		Node vert1(1, 1.0, 0.0, 0.0);
		Node vert2(2, 0.0, 1.0, 0.0);
		Node vert3(3, 0.0, 0.0, 1.0);				
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(-1);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);		
		
		// add nodes to element
		Element3DTR elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		elem.prepare();
		
		// check element volume (V == 1/6)
		TS_ASSERT_DELTA(elem.get_volume(), 1.0/6.0, 1.0e-9);
		
		// check jacobi
		TS_ASSERT_DELTA(elem.get_jacobi_det(), 1.0, 1.0e-9);
		

		double cmp = 0;
		for(short ii = 0; ii < 4; ii++) {
//			std::cout << "\n";
			for(short jj = 0; jj < 4; jj++) {
//				std::cout << elem.get_element_integral_2nd_order(1, ii, 0, jj) << " ";
				// check zero order functions
				if(ii == jj) {
					cmp = 1.0 / 60.0;	
				} else {
					cmp = 1.0 / 120.0;							
				}
				TS_ASSERT_DELTA(elem.get_element_integral_0th_order(ii,jj),cmp, 1.0e-9);	
				
				// check second order functions (d/dx d/dx, d/dy d/dy, d/dz d/dz)
				for(short diffop = 0; diffop < 3; diffop++) {
					if(ii == jj) {
						if(ii == 0 || ii == diffop + 1) {
							cmp =  1.0 / 6.0;	
						} else  {
							cmp = 0.0;
						}	
					} else {
						if((ii == 0 && jj == diffop	+ 1) || (jj == 0 && ii == diffop + 1)) {
							cmp = -1.0 / 6.0;	
						} else {
							cmp = 0.0;	
						}
					}
					TS_ASSERT_DELTA(elem.get_element_integral_2nd_order(diffop, ii, diffop, jj), cmp, 1.0e-9);
				}	
				// check first order functions
				for(short diffop = 0; diffop < 3; diffop++) {
					if(ii == 0) {
						cmp = - 1.0 / 24.0;	
					} else if(ii == diffop + 1) {
						cmp = 1.0 / 24.0;	
					} else {
						cmp = 0.0;	
					}
					TS_ASSERT_DELTA(elem.get_element_integral_1st_order(diffop, ii, jj), cmp, 1.0e-9);
				}
			}
		}																												
	}

	
void ElementTest::test_3d_2nd_order_tetraeder_numerical() {

	try {
		
		// create basic nodes
		vector<Node*> nodes;
		nodes.push_back(new Node((unsigned)0,  0.0,  0.0, -1.0));
		nodes.push_back(new Node((unsigned)1,  -1.0, 0.5, 0.0));		
		nodes.push_back(new Node((unsigned)2,  0.1,  1.0, 0.05));
		nodes.push_back(new Node((unsigned)3,  0.03,  0.06, 0.1));
																								
		// add nodes to element
		Element3DTetrahedron2nd elem(0);
		elem.set_region(&region, 0);
		for(unsigned int ii = 0; ii < 4; ii++) {
			nodes[ii]->set_index_internal(ii);
			nodes[ii]->set_index_vertex(ii);			
			elem.set_corner_node(ii, nodes[ii]);
		}
		// add additional nodes
		set_additional_nodes(elem);		
		elem.prepare();
		test_if_basis_is_nodal(elem, "test_3d_2nd_order_tetraeder_all_numerical");
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_3d_2nd_order_tetraeder_all_numerical");		
			
		for(unsigned int ii = 0; ii < nodes.size(); ii++) {
			delete nodes[ii];	
		}			
			
	} catch (Exception* e) {
		TS_FAIL(e->get_reason());	
	}	
	
}
	
void ElementTest::test_3d_tetraeder_all_numerical() {
	{
		if(skip_3d_tests) {
			TS_WARN("SKIPPING 3D TETRAEDER ALL NUMERICAL CHECK!");
			return;
		}
				
		// create basic nodes
		Node vert0((unsigned)0,  0.0,  0.0, 0.0);
		Node vert1((unsigned)1,  1.0,  0.0, 0.0);		
		Node vert2((unsigned)2,  0.0,  1.0, 0.0);
		Node vert3((unsigned)3,  0.0,  0.0, 1.0);
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);			
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(2);
				
								
		// add nodes to element
		Element3DTR elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		try {
			elem.prepare();
		} catch(Exception *e) {
			cout << e->get_reason();	
		}
		TDKP_TRACE("building element numerical validator");
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_3d_tetraeder_all_numerical");		
	}
}

void ElementTest::test_3d_tetraeder_all_numerical_nonsimple() {
	{
		
		if(skip_3d_tests) {
			TS_WARN("SKIPPING 3D TETRAEDER 2ND ORDER ALL NUMERICAL CHECK!");
			return;
		}
		
		// create basic nodes
		Node vert0((unsigned)0,  1.2,   1.3,-0.3);
		Node vert1((unsigned)1, -0.1,   1.1, 0.1);		
		Node vert2((unsigned)2, -0.12,  0.1, -0.2);
		Node vert3((unsigned)3,  0.5,  0.5, 0.8);
		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);

		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);
								
		// add nodes to element
		Element3DTR elem(0);
		elem.set_region(&region, 0);
		elem.set_corner_node(0, &vert0);
		elem.set_corner_node(1, &vert1);		
		elem.set_corner_node(2, &vert2);
		elem.set_corner_node(3, &vert3);
		try {
			elem.prepare();
		} catch(Exception *e) {
			cout << e->get_reason();	
		}		
		
		// ----------------------------------------------
		// evaluate edges / faces / corners
		// ----------------------------------------------
		TS_ASSERT_EQUALS(elem.get_shape(), Element::tetrahedron);
		TS_ASSERT_EQUALS(elem.get_num_corners(), (unsigned int)4);
		TS_ASSERT_EQUALS(elem.get_num_edges(),   (unsigned int)6);
		TS_ASSERT_EQUALS(elem.get_num_faces(),   (unsigned int)4);
		vector<unsigned int> vertex_indices;
		vector<unsigned int> vertex_indices_other;
		for(unsigned int ii = 0; ii < 4; ii++) {
			elem.get_face(ii, vertex_indices);
			// not twice the same index			
			TS_ASSERT_EQUALS(vertex_indices.size(), (unsigned int)3);
			for(unsigned int vv = 0; vv < 3; vv++) {
				for(unsigned int ww = vv + 1; ww < 3; ww++) {
					TS_ASSERT(vertex_indices[vv] != vertex_indices[ww]);	
				}	
			}
			// at least one vertex must differ compared to other faces
			for(unsigned int jj = ii + 1; jj < 4; jj++) {
				elem.get_face(jj, vertex_indices_other);
				TS_ASSERT_EQUALS(vertex_indices_other.size(), (unsigned int)3);
				bool all_same = true;
				for(unsigned int vv = 0; vv < 3; vv++) {
					bool found_it = false;
					for(unsigned int ww = 0; ww < 3; ww++) {
						if(vertex_indices[vv] == vertex_indices_other[ww]) {
							found_it = true;
							break;
						} 							
					}	
					// if i didnt found the vertex, the faces are not equal
					if(!found_it) {
						all_same = false;							
					}
				}
				TS_ASSERT(!all_same);								
			}
		}
		// edge check
		for(unsigned int ii = 0; ii < 6; ii++) {
			elem.get_edge(ii, vertex_indices);
			TS_ASSERT_EQUALS(vertex_indices.size(), (unsigned int)2);
			TS_ASSERT(vertex_indices[0] != vertex_indices[1]);
			// compare with other edges
			for(unsigned int jj = ii + 1; jj < 6; jj++) {
				elem.get_edge(jj, vertex_indices_other);
				TS_ASSERT_EQUALS(vertex_indices_other.size(), (unsigned int)2);
				TS_ASSERT(vertex_indices_other[0] != vertex_indices_other[1]);										
				bool all_same = true;
				for(unsigned int vv = 0; vv < 3; vv++) {
					bool found_it = false;
					for(unsigned int ww = 0; ww < 3; ww++) {
						if(vertex_indices[vv] == vertex_indices_other[ww]) {
							found_it = true;
							break;
						} 							
					}	
					// if i didnt found the vertex, the faces are not equal
					if(!found_it) {
						all_same = false;							
					}
				}
				TS_ASSERT(!all_same);
			}											
		}				
		
		
		ElementNumericalValidator validator(elem);
		validator.test_evaluate_integral("test_3d_tetraeder_all_numerical_nonsimple");
		
	}
}

void ElementTest::test_1d_boundary() {
	try {
		Node vert0((unsigned)0,  -1.0);
		Node vert1((unsigned)1,   0.0);		
		Node vert2((unsigned)2,   1.0);
			
			
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);	
		vert2.set_index_vertex(2);				
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);

		// build boundaries
		ElementBoundary bnd0(ElementBoundary::vertex, 0, 1);
		ElementBoundary bnd1(ElementBoundary::vertex, 0, 1);
		bnd0.set_location('i');
		bnd0.set_corner_node(0, &vert1);
		// swapped 
		bnd1.set_location('i');
		bnd1.set_corner_node(0, &vert1);

									
		// add nodes to elements
		Element1DLine elem0(0);
		Element1DLine elem1(0);
		elem0.set_region(&region, 0);
		elem1.set_region(&region, 0);
		elem0.set_corner_node(0, &vert0);
		elem0.set_corner_node(1, &vert1);
		elem1.set_corner_node(0, &vert1);		
		elem1.set_corner_node(1, &vert2);									
		bnd0.add_element(&elem0);
		bnd0.add_element(&elem1);
		// other way around in bnd1
		bnd1.add_element(&elem1);
		bnd1.add_element(&elem0);
		bnd0.prepare();
		bnd1.prepare();
		// check normal
		TS_ASSERT_DELTA(bnd0.get_normal()(0), 1.0, 1.0e-8);
		TS_ASSERT_DELTA(bnd1.get_normal()(0), -1.0, 1.0e-8);
		
	} catch (Exception *e) {
		TS_FAIL(e->get_reason());	
	}
	try {
		// -----------------------------------------------
		// the same for the 1d element with additional
		// nodes
		// -----------------------------------------------
		Node vert0((unsigned)0,  -1.0);
		Node vert0a((unsigned)0, -0.5);
		Node vert1((unsigned)1,   0.0);
		Node vert1a((unsigned)1,  0.5);				
		Node vert2((unsigned)2,   1.0);

		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);	
		vert2.set_index_vertex(2);
				
		vert0.set_index_internal(0);
		vert0a.set_index_internal(1);
		vert1.set_index_internal(2);
		vert1a.set_index_internal(3);		
		vert2.set_index_internal(4);

		// build boundaries
		ElementBoundary bnd0(ElementBoundary::vertex, 0, 1);
		ElementBoundary bnd1(ElementBoundary::vertex, 0, 1);
		bnd0.set_location('i');
		bnd0.set_corner_node(0, &vert1);
		// swapped 
		bnd1.set_location('i');
		bnd1.set_corner_node(0, &vert1);

									
		// add nodes to elements
		Element1DLine2nd elem0(0);
		Element1DLine2nd elem1(0);
		elem0.set_region(&region, 0);
		elem1.set_region(&region, 0);
		elem0.set_corner_node(0, &vert0);
		elem0.set_corner_node(1, &vert1);
		elem0.set_additional_node(0, &vert0a);
		elem1.set_corner_node(0, &vert1);		
		elem1.set_corner_node(1, &vert2);		
		elem1.set_additional_node(0, &vert1a);							
		bnd0.add_element(&elem0);
		bnd0.add_element(&elem1);
		// other way around in bnd1
		bnd1.add_element(&elem1);
		bnd1.add_element(&elem0);
		bnd0.prepare();
		bnd1.prepare();		
		// check normal
		TS_ASSERT_DELTA(bnd0.get_normal()(0), 1.0, 1.0e-8);
		TS_ASSERT_DELTA(bnd1.get_normal()(0), -1.0, 1.0e-8);
		
	} catch (Exception *e) {
		TS_FAIL(e->get_reason());	
	}			
}

void ElementTest::test_2d_boundary() {
	try {
		
		Node vert0((unsigned)0, 0.0, 0.0);
		Node vert1((unsigned)1, 1.0, 0.0);		
		Node vert2((unsigned)2, 1.0, 1.0);
		Node vert3((unsigned)3, 0.0, 1.0);
						
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);	
		vert2.set_index_vertex(2);				
		vert3.set_index_vertex(3);		
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);

		// build boundaries
		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		ElementBoundary bnd1(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		bnd0.set_corner_node(0, &vert1);
		bnd0.set_corner_node(1, &vert3);
		// swapped 
		bnd1.set_location('i');
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert3);

									
		// add nodes to elements
		Element2DTriangle elem0(0);
		Element2DTriangle elem1(0);
		elem0.set_region(&region, 0);
		elem1.set_region(&region, 0);
		elem0.set_corner_node(0, &vert0);
		elem0.set_corner_node(1, &vert1);
		elem0.set_corner_node(2, &vert3);		
		elem1.set_corner_node(0, &vert1);		
		elem1.set_corner_node(1, &vert2);
		elem1.set_corner_node(2, &vert3);		
											
		bnd0.add_element(&elem0);
		bnd0.add_element(&elem1);
		// other way around in bnd1
		bnd1.add_element(&elem1);
		bnd1.add_element(&elem0);
		bnd0.prepare();
		bnd1.prepare();				
		// check normal
		TS_ASSERT_DELTA(bnd0.get_normal()(0),  1.0 / constants::sqrt2, 1.0e-8);
		TS_ASSERT_DELTA(bnd0.get_normal()(1),  1.0 / constants::sqrt2, 1.0e-8);
		TS_ASSERT_DELTA(bnd1.get_normal()(0), -1.0 / constants::sqrt2, 1.0e-8);		
		TS_ASSERT_DELTA(bnd1.get_normal()(1), -1.0 / constants::sqrt2, 1.0e-8);		
		
	} catch (Exception *e) {
		TS_FAIL(e->get_reason());	
	}
	try {
		// -----------------------------------------------
		// the same for the 2d element with additional
		// nodes
		// -----------------------------------------------
		Node vert0 ((unsigned)0, 0.0, 0.0);
		Node vert0a((unsigned)1, 0.5, 0.0);
		Node vert1 ((unsigned)2, 1.0, 0.0);
		Node vert1a((unsigned)3, 1.0, 0.5);		
		Node vert2 ((unsigned)4, 1.0, 1.0);
		Node vert2a((unsigned)5, 0.5, 1.0);
		Node vert3 ((unsigned)6, 0.0, 1.0);
		Node vert3a((unsigned)7, 0.0, 0.5);
		Node vert4a((unsigned)8, 0.5, 0.5);

		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);	
		vert2.set_index_vertex(2);
		vert3.set_index_vertex(3);
				
		vert0. set_index_internal(0);
		vert0a.set_index_internal(1);
		vert1. set_index_internal(2);
		vert1a.set_index_internal(3);		
		vert2. set_index_internal(4);
		vert2a.set_index_internal(5);
		vert3. set_index_internal(6);
		vert3a.set_index_internal(7);
		vert4a.set_index_internal(8);		

		// build boundaries
		ElementBoundary bnd0(ElementBoundary::edge, 0, 2);
		ElementBoundary bnd1(ElementBoundary::edge, 0, 2);
		bnd0.set_location('i');
		bnd0.set_corner_node(0, &vert1);
		bnd0.set_corner_node(1, &vert3);
		// swapped 
		bnd1.set_location('i');
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert3);
								
		// add nodes to elements
		Element2DTriangle2nd elem0(0);
		Element2DTriangle2nd elem1(0);
		elem0.set_region(&region, 0);
		elem1.set_region(&region, 0);
		elem0.set_corner_node(0, &vert0);
		elem0.set_corner_node(1, &vert1);
		elem0.set_corner_node(2, &vert3);
		elem0.set_additional_node(0, &vert0a);
		elem0.set_additional_node(1, &vert4a);
		elem0.set_additional_node(2, &vert3a);		
		elem1.set_corner_node(0, &vert1);		
		elem1.set_corner_node(1, &vert2);
		elem1.set_corner_node(2, &vert3);		
		elem1.set_additional_node(0, &vert1a);
		elem1.set_additional_node(1, &vert2a);
		elem1.set_additional_node(2, &vert4a);

		bnd0.add_element(&elem0);
		bnd0.add_element(&elem1);
		// other way around in bnd1
		bnd1.add_element(&elem1);
		bnd1.add_element(&elem0);
		bnd0.prepare();
		bnd1.prepare();		
		// check normal
		TS_ASSERT_DELTA(bnd0.get_normal()(0),  1.0 / constants::sqrt2, 1.0e-8);
		TS_ASSERT_DELTA(bnd0.get_normal()(1),  1.0 / constants::sqrt2, 1.0e-8);
		TS_ASSERT_DELTA(bnd1.get_normal()(0), -1.0 / constants::sqrt2, 1.0e-8);		
		TS_ASSERT_DELTA(bnd1.get_normal()(1), -1.0 / constants::sqrt2, 1.0e-8);	
		
	} catch (Exception *e) {
		TS_FAIL(e->get_reason());	
	}			
}

void ElementTest::test_3d_boundary() {
	try {
		
		Node vert0((unsigned)0, 0.0, 0.0, 0.0);
		Node vert1((unsigned)1, 1.0, 0.0, 0.0);		
		Node vert2((unsigned)2, 0.0, 1.0, 0.0);
		Node vert3((unsigned)3, 0.0, 0.0, 1.0);
		Node vert4((unsigned)3, 1.0, 1.0, 1.0);
							
		vert0.set_index_vertex(0);
		vert1.set_index_vertex(1);	
		vert2.set_index_vertex(2);				
		vert3.set_index_vertex(3);
		vert4.set_index_vertex(4);			
		vert0.set_index_internal(0);
		vert1.set_index_internal(1);
		vert2.set_index_internal(2);
		vert3.set_index_internal(3);
		vert4.set_index_internal(4);	

		// build boundaries
		ElementBoundary bnd0(ElementBoundary::face, 0, 3);
		ElementBoundary bnd1(ElementBoundary::face, 0, 3);
		bnd0.set_location('i');
		bnd0.set_corner_node(0, &vert1);
		bnd0.set_corner_node(1, &vert2);
		bnd0.set_corner_node(2, &vert3);
		// swapped 
		bnd1.set_location('i');
		bnd1.set_corner_node(0, &vert1);
		bnd1.set_corner_node(1, &vert2);
		bnd1.set_corner_node(2, &vert3);
									
		// add nodes to elements
		Element3DTR elem0(0);
		Element3DTR elem1(1);
		elem0.set_region(&region, 0);
		elem1.set_region(&region, 0);
		elem0.set_corner_node(0, &vert0);
		elem0.set_corner_node(1, &vert1);
		elem0.set_corner_node(2, &vert2);
		elem0.set_corner_node(3, &vert3);			
		elem1.set_corner_node(0, &vert1);		
		elem1.set_corner_node(1, &vert2);
		elem1.set_corner_node(2, &vert3);
		elem1.set_corner_node(3, &vert4);		
											
		bnd0.add_element(&elem0);
		bnd0.add_element(&elem1);
		// other way around in bnd1
		bnd1.add_element(&elem1);
		bnd1.add_element(&elem0);
		bnd0.prepare();
		bnd1.prepare();				
		// check normal
		TS_ASSERT_DELTA(bnd0.get_normal()(0),  1.0 / constants::sqrt3, 1.0e-8);
		TS_ASSERT_DELTA(bnd0.get_normal()(1),  1.0 / constants::sqrt3, 1.0e-8);
		TS_ASSERT_DELTA(bnd0.get_normal()(2),  1.0 / constants::sqrt3, 1.0e-8);
		TS_ASSERT_DELTA(bnd1.get_normal()(0), -1.0 / constants::sqrt3, 1.0e-8);		
		TS_ASSERT_DELTA(bnd1.get_normal()(1), -1.0 / constants::sqrt3, 1.0e-8);		
		TS_ASSERT_DELTA(bnd1.get_normal()(2), -1.0 / constants::sqrt3, 1.0e-8);
		
	} catch (Exception *e) {
		TS_FAIL(e->get_reason());	
	}
}

void ElementTest::test_1d_line_all_numerical() {
	
	Node a((unsigned int)0,  -1.0);
	Node b((unsigned int)1,  2.3);
	
	a.set_index_vertex(0);
	b.set_index_vertex(1);
	
	Element1DLine elem(0);
	elem.set_region(&region,0);
	elem.set_corner_node(0,&a);
	elem.set_corner_node(1,&b);
	elem.prepare(); 
							
	ElementNumericalValidator validator(elem);
	validator.test_evaluate_integral("test_1d_line_all_numerical");
	
}

void ElementTest::test_1d_2nd_line_all_numerical() {

	Node a((unsigned int)0, -1.0);
	Node b((unsigned int)1,  1.5);
	Node c((unsigned int)2,  0.25);

	a.set_index_vertex(0);
	b.set_index_vertex(1);
		
	Element1DLine2nd elem(0);
	elem.set_region(&region,0);
	elem.set_corner_node(0,&a);
	elem.set_corner_node(1,&b);
	elem.set_additional_node(0, &c);
	elem.prepare(); 
							
	ElementNumericalValidator validator(elem);
	validator.test_evaluate_integral("test_1d_2nd_line_all_numerical");
	
}

void ElementTest::test_1d_cubic_hermite_numerical() {
	
	Node a((unsigned int)0, -1.121);
	Node b((unsigned int)1,  1.5);
	
	a.set_index_vertex(0);
	b.set_index_vertex(1);
		
	Element1DHermiteCubic elem(0);
	elem.set_region(&region,0);
	elem.set_corner_node(0,&a);
	elem.set_corner_node(1,&b);
	set_additional_nodes(elem);
	elem.prepare(); 
							
	ElementNumericalValidator validator(elem);	
	validator.test_evaluate_integral("test_1d_cubic_hermite_numerical");
	validator.plot_shape_functions("1d_cubic_hermite.dat");
	
}

void ElementTest::set_additional_nodes(Element& elem) {
	Element::AdditionalNodeLocation location_type;
	vector<unsigned int> involved_vertices;
	vector<double> coords;
	unsigned int tag;
	
	for(unsigned int ii = 0; ii < elem.get_num_additional_nodes(); ii++) {
		elem.get_additional_node_locator(ii, location_type, involved_vertices, coords, tag);
		Node* node = 0;
		switch(elem.get_dimension()) {
			case 1:
				node = new Node(5 + ii, coords[0]);
				break;
			case 2:
				node = new Node(5 + ii, coords[0], coords[1]);
				break;
			case 3:
				node = new Node(5 + ii, coords[0], coords[1], coords[2]);
				break;									
		}
		elem.set_additional_node(ii, node);
	}	
}


