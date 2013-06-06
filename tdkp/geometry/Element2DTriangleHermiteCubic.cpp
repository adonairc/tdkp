// ------------------------------------------------------------
//
// This file is part of tdkp, a simulation tool for nanostrutctures
// of optoelectronics developed at ETH Zurich
//
// (C) 2005-2009 Ratko G. Veprek, ETH Zurich, veprek@iis.ee.ethz.ch
//
// 1) As of 18.6.2009 this code is property of ETH Zurich and must not be
// transferred, modified or used by third parties without appropriate
// licenses issued by authorized agents of ETH Zurich.
//
// 2) Violation of this will result in judicial action according to civil
// and penal law.
//
// 3) Any claim of authorship other than by the author himself is
// strictly forbidden.
//
// 4) The source code must retain the copyright notice, this list
// of conditions and the following disclaimer.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS
// BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
// IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------------------------------------

#include "tdkp/geometry/Element2DTriangleHermiteCubic.h"

namespace tdkp {

Element2DTriangleHermiteCubic::Element2DTriangleHermiteCubic(unsigned int index_global_) 
: Element2DTriangleBase(2,10,3) 
{
	this->index_global = index_global_;	
}
	
Element2DTriangleHermiteCubic::~Element2DTriangleHermiteCubic() {
	
}

void Element2DTriangleHermiteCubic::set_corner_node(unsigned short corner_id, Node* node) {
	unsigned int node_idx = corner_id * 3;
	TDKP_BOUNDS_ASSERT(node_idx < nodes.size(), "");
	TDKP_BOUNDS_ASSERT(nodes[node_idx] == 0, "");
	nodes[node_idx] = node;	
}

const Node& Element2DTriangleHermiteCubic::get_corner_node(unsigned short corner_idx) const {
	unsigned int node_idx = corner_idx * 3;
	TDKP_BOUNDS_ASSERT(node_idx < nodes.size(), "");
	return *nodes[node_idx];		
}			
 
void Element2DTriangleHermiteCubic::get_additional_node_locator(
	unsigned int additional_node_idx, 
	AdditionalNodeLocation& location_type, 
	vector<unsigned int>& involved_vertices, 
	vector<double>& coords, 
	unsigned int& tag
) const {
	coords.assign(2,0.0);
	// -----------------------------------------
	// mid point node is different
	// -----------------------------------------
	if(additional_node_idx == 6) {
		location_type = Element::inner_node;
		involved_vertices.resize(0);		
		for(unsigned int ii = 0; ii < 3; ii++) {
			coords[0] += get_corner_node(ii).get_coord(0) / 3.0;
			coords[1] += get_corner_node(ii).get_coord(1) / 3.0;			
		}
		tag = 0;
	} else if(additional_node_idx < 6) {
		// ----------------------------------------
		// first derivatives at vertices
		// ----------------------------------------
		unsigned int corner_idx = (additional_node_idx - (additional_node_idx % 2)) / 2;
		location_type = Element::vertex_node;				
		involved_vertices.resize(1);
		// 0,1 are at corner 0, 2,3 are at corner 1, and 4,5 are at corner 2
		involved_vertices[0] = get_corner_node(corner_idx).get_index_global();
		coords[0] = get_corner_node(corner_idx).get_coord(0);
		coords[1] = get_corner_node(corner_idx).get_coord(1);
		tag = additional_node_idx % 2;
	} else {
		TDKP_GENERAL_EXCEPTION("invalid additional node index");	
	}
}

void Element2DTriangleHermiteCubic::set_additional_node(
	unsigned int additional_node_idx, Node* node
) {
	// -----------------------------------------
	// mid point node
	// -----------------------------------------
	unsigned int node_idx = 0;
	Node::NodeValueType value_type = Node::FunctionValue;
	if(additional_node_idx == 6) {
		// N9
		node_idx = 9; 
	} else if(additional_node_idx < 6) {
		// ----------------------------------------
		// first derivatives at vertices
		// ----------------------------------------		
		if(additional_node_idx % 2 == 0) {
			// A0,A2,A4, -> N1, N4, N7 
			value_type = Node::Derivative_X;
			// node idx  = corner * 3 + 1	
			node_idx = (additional_node_idx / 2) * 3 + 1; 
		} else {
			value_type = Node::Derivative_Y;
			// node idx  = corner * 3 + 2	
			node_idx = ((additional_node_idx - 1) / 2) * 3 + 2;			
		}		
	} else {
		TDKP_GENERAL_EXCEPTION("invalid additional node index");	
	}	
	TDKP_BOUNDS_ASSERT(nodes[node_idx] == 0, "");
	nodes[node_idx] = node;
	if(value_type != Node::FunctionValue) {
		node->set_value_type(value_type);	
	}
}

void Element2DTriangleHermiteCubic::get_node_local_coords(
	unsigned short lid, 
	vector<double>& local_coords
) const {
	local_coords.assign(2, 0.0);
	if(lid < 3) {
		// 0,0
	} else if(lid < 6) {
		local_coords[0] = 1.0; // (1,0)
	} else if(lid < 9) {
		local_coords[1] = 1.0; // (0,1)
	} else if(lid == 9) {
		local_coords[0] = local_coords[1] = 1.0 / 3.0; // (1/3,1/3)
	} else {
		TDKP_GENERAL_EXCEPTION("invalid node index");	
	}
}

double Element2DTriangleHermiteCubic::get_element_integral_2nd_order(
	short diffop_1, 
	short elem_shape_func_1, 
	short diffop_2, 
	short elem_shape_func_2
) const {
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1
	                       + 10 * 10 * diffop_2
	                       + 10 * 10 * 2 * diffop_1;
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
        case 0: // da1db1Nii1Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.990000e+02) * (jm[2] * jm[2])) + ((-3.490000e+02) * jm[2] * jm[3]) + ((1.990000e+02) * (jm[3] * jm[3])));
        case 10: // da1db1Nii2Njj1
        case 1: // da1db1Nii1Njj2
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.500000e+01) * (jm[2] * jm[2])) + ((-3.100000e+01) * jm[2] * jm[3]) + ((3.000000e+00) * (jm[3] * jm[3])));
        case 20: // da1db1Nii3Njj1
        case 2: // da1db1Nii1Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.000000e+00) * (jm[2] * jm[2])) + ((-3.100000e+01) * jm[2] * jm[3]) + ((3.500000e+01) * (jm[3] * jm[3])));
        case 30: // da1db1Nii4Njj1
        case 3: // da1db1Nii1Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+01) * (jm[2] * jm[2])) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[3] * jm[3]));
        case 40: // da1db1Nii5Njj1
        case 4: // da1db1Nii1Njj5
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((1.000000e+01) * (jm[2] * jm[2])) + ((-4.000000e+00) * jm[2] * jm[3]) + (jm[3] * jm[3]));
        case 50: // da1db1Nii6Njj1
        case 5: // da1db1Nii1Njj6
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (jm[2] + ((-1.000000e+00) * jm[3])) * (((2.000000e+00) * jm[2]) + ((5.000000e+00) * jm[3]));
        case 60: // da1db1Nii7Njj1
        case 6: // da1db1Nii1Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * ((jm[2] * jm[2]) + ((-1.000000e+00) * jm[2] * jm[3]) + ((7.000000e+01) * (jm[3] * jm[3])));
        case 70: // da1db1Nii8Njj1
        case 7: // da1db1Nii1Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (jm[2] + ((-1.000000e+00) * jm[3])) * (((5.000000e+00) * jm[2]) + ((2.000000e+00) * jm[3]));
        case 80: // da1db1Nii9Njj1
        case 8: // da1db1Nii1Njj9
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * ((jm[2] * jm[2]) + ((-4.000000e+00) * jm[2] * jm[3]) + ((1.000000e+01) * (jm[3] * jm[3])));
        case 90: // da1db1Nii10Njj1
        case 9: // da1db1Nii1Njj10
            return ((-1.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((3.000000e+01) * (jm[2] * jm[2])) + ((-3.900000e+01) * jm[2] * jm[3]) + ((3.000000e+01) * (jm[3] * jm[3])));
        case 11: // da1db1Nii2Njj2
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[2] * jm[2])) + ((-3.000000e+00) * jm[2] * jm[3]) + ((3.000000e+00) * (jm[3] * jm[3])));
        case 21: // da1db1Nii3Njj2
        case 12: // da1db1Nii2Njj3
            return ((-1.000000e+00) / (3.600000e+01)) * (1.0 / jacobi_det) * jm[2] * jm[3];
        case 31: // da1db1Nii4Njj2
        case 13: // da1db1Nii2Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[2]) + ((-3.000000e+00) * jm[3])) * (jm[2] + jm[3]);
        case 41: // da1db1Nii5Njj2
        case 14: // da1db1Nii2Njj5
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[2] * (((4.000000e+00) * jm[2]) + jm[3]);
        case 51: // da1db1Nii6Njj2
        case 15: // da1db1Nii2Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[2] * (((-1.000000e+00) * jm[2]) + jm[3]);
        case 61: // da1db1Nii7Njj2
        case 16: // da1db1Nii2Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[2] * (((5.000000e+00) * jm[2]) + ((-7.000000e+00) * jm[3]));
        case 71: // da1db1Nii8Njj2
        case 17: // da1db1Nii2Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[2] * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 81: // da1db1Nii9Njj2
        case 18: // da1db1Nii2Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[2] * (jm[2] + ((-2.000000e+00) * jm[3]));
        case 91: // da1db1Nii10Njj2
        case 19: // da1db1Nii2Njj10
            return ((3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * jm[2] * (((-2.000000e+00) * jm[2]) + jm[3]);
        case 22: // da1db1Nii3Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.000000e+00) * (jm[2] * jm[2])) + ((-3.000000e+00) * jm[2] * jm[3]) + ((7.000000e+00) * (jm[3] * jm[3])));
        case 32: // da1db1Nii4Njj3
        case 23: // da1db1Nii3Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[3] * (((-7.000000e+00) * jm[2]) + ((5.000000e+00) * jm[3]));
        case 42: // da1db1Nii5Njj3
        case 24: // da1db1Nii3Njj5
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[3] * (((-2.000000e+00) * jm[2]) + jm[3]);
        case 52: // da1db1Nii6Njj3
        case 25: // da1db1Nii3Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[3] * (((-1.000000e+00) * jm[2]) + jm[3]);
        case 62: // da1db1Nii7Njj3
        case 26: // da1db1Nii3Njj7
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.000000e+00) * jm[2]) + ((-1.400000e+01) * jm[3])) * (jm[2] + jm[3]);
        case 72: // da1db1Nii8Njj3
        case 27: // da1db1Nii3Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (jm[2] + ((-1.000000e+00) * jm[3])) * jm[3];
        case 82: // da1db1Nii9Njj3
        case 28: // da1db1Nii3Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[3] * (jm[2] + ((4.000000e+00) * jm[3]));
        case 92: // da1db1Nii10Njj3
        case 29: // da1db1Nii3Njj10
            return ((3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (jm[2] + ((-2.000000e+00) * jm[3])) * jm[3];
        case 33: // da1db1Nii4Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.900000e+01) * (jm[2] * jm[2])) + ((-4.900000e+01) * jm[2] * jm[3]) + ((1.990000e+02) * (jm[3] * jm[3])));
        case 43: // da1db1Nii5Njj4
        case 34: // da1db1Nii4Njj5
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[2] * jm[2])) + ((-7.000000e+00) * jm[2] * jm[3]) + ((1.900000e+01) * (jm[3] * jm[3])));
        case 53: // da1db1Nii6Njj4
        case 35: // da1db1Nii4Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[2] * jm[2])) + ((-3.900000e+01) * jm[2] * jm[3]) + ((3.500000e+01) * (jm[3] * jm[3])));
        case 63: // da1db1Nii7Njj4
        case 36: // da1db1Nii4Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+01) * (jm[2] * jm[2])) + ((-1.390000e+02) * jm[2] * jm[3]) + ((7.000000e+01) * (jm[3] * jm[3])));
        case 73: // da1db1Nii8Njj4
        case 37: // da1db1Nii4Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[2]) + ((-1.100000e+01) * jm[3])) * (jm[2] + ((-2.000000e+00) * jm[3]));
        case 83: // da1db1Nii9Njj4
        case 38: // da1db1Nii4Njj9
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[2] * jm[2])) + ((-1.600000e+01) * jm[2] * jm[3]) + ((1.000000e+01) * (jm[3] * jm[3])));
        case 93: // da1db1Nii10Njj4
        case 39: // da1db1Nii4Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[2] * jm[2])) + ((-7.000000e+00) * jm[2] * jm[3]) + ((1.000000e+01) * (jm[3] * jm[3])));
        case 44: // da1db1Nii5Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((2.000000e+00) * (jm[2] * jm[2])) + ((-2.000000e+00) * jm[2] * jm[3]) + ((5.000000e+00) * (jm[3] * jm[3])));
        case 54: // da1db1Nii6Njj5
        case 45: // da1db1Nii5Njj6
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * (jm[2] * jm[2])) + ((-6.000000e+00) * jm[2] * jm[3]) + ((7.000000e+00) * (jm[3] * jm[3])));
        case 64: // da1db1Nii7Njj5
        case 46: // da1db1Nii5Njj7
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((1.000000e+01) * (jm[2] * jm[2])) + ((-1.600000e+01) * jm[2] * jm[3]) + ((7.000000e+00) * (jm[3] * jm[3])));
        case 74: // da1db1Nii8Njj5
        case 47: // da1db1Nii5Njj8
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[2]) + ((-5.000000e+00) * jm[3])) * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 84: // da1db1Nii9Njj5
        case 48: // da1db1Nii5Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * (jm[2] * jm[2])) + ((-7.000000e+00) * jm[2] * jm[3]) + ((4.000000e+00) * (jm[3] * jm[3])));
        case 89: // da1db1Nii9Njj10
        case 94: // da1db1Nii10Njj5
        case 98: // da1db1Nii10Njj9
        case 49: // da1db1Nii5Njj10
            return ((3.000000e+00) / (1.000000e+01)) * (1.0 / jacobi_det) * ((jm[2] * jm[2]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[3] * jm[3]));
        case 77: // da1db1Nii8Njj8
        case 55: // da1db1Nii6Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[2] * jm[2])) + ((-1.100000e+01) * jm[2] * jm[3]) + ((7.000000e+00) * (jm[3] * jm[3])));
        case 65: // da1db1Nii7Njj6
        case 56: // da1db1Nii6Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.200000e+01) * (jm[2] * jm[2])) + ((-3.900000e+01) * jm[2] * jm[3]) + ((1.400000e+01) * (jm[3] * jm[3])));
        case 75: // da1db1Nii8Njj6
        case 57: // da1db1Nii6Njj8
            return ((1.000000e+00) / (3.600000e+01)) * (1.0 / jacobi_det) * ((jm[2] + ((-1.000000e+00) * jm[3])) * (jm[2] + ((-1.000000e+00) * jm[3])));
        case 85: // da1db1Nii9Njj6
        case 58: // da1db1Nii6Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((5.000000e+00) * jm[2]) + ((-4.000000e+00) * jm[3])) * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 95: // da1db1Nii10Njj6
        case 59: // da1db1Nii6Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (jm[2] + ((-2.000000e+00) * jm[3])) * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 66: // da1db1Nii7Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.990000e+02) * (jm[2] * jm[2])) + ((-4.900000e+01) * jm[2] * jm[3]) + ((4.900000e+01) * (jm[3] * jm[3])));
        case 76: // da1db1Nii8Njj7
        case 67: // da1db1Nii7Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.500000e+01) * (jm[2] * jm[2])) + ((-3.900000e+01) * jm[2] * jm[3]) + ((7.000000e+00) * (jm[3] * jm[3])));
        case 86: // da1db1Nii9Njj7
        case 68: // da1db1Nii7Njj9
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((1.900000e+01) * (jm[2] * jm[2])) + ((-7.000000e+00) * jm[2] * jm[3]) + ((7.000000e+00) * (jm[3] * jm[3])));
        case 96: // da1db1Nii10Njj7
        case 69: // da1db1Nii7Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((1.000000e+01) * (jm[2] * jm[2])) + ((-7.000000e+00) * jm[2] * jm[3]) + ((7.000000e+00) * (jm[3] * jm[3])));
        case 87: // da1db1Nii9Njj8
        case 78: // da1db1Nii8Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[2] * jm[2])) + ((-6.000000e+00) * jm[2] * jm[3]) + ((2.000000e+00) * (jm[3] * jm[3])));
        case 97: // da1db1Nii10Njj8
        case 79: // da1db1Nii8Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((2.000000e+00) * (jm[2] * jm[2])) + ((-3.000000e+00) * jm[2] * jm[3]) + (jm[3] * jm[3]));
        case 88: // da1db1Nii9Njj9
            return ((1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((5.000000e+00) * (jm[2] * jm[2])) + ((-2.000000e+00) * jm[2] * jm[3]) + ((2.000000e+00) * (jm[3] * jm[3])));
        case 99: // da1db1Nii10Njj10
            return ((8.100000e+01) / (2.000000e+01)) * (1.0 / jacobi_det) * ((jm[2] * jm[2]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[3] * jm[3]));
        case 200: // da2db1Nii1Njj1
        case 100: // da1db2Nii1Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-3.980000e+02) * jm[0] * jm[2]) + ((3.490000e+02) * jm[1] * jm[2]) + ((3.490000e+02) * jm[0] * jm[3]) + ((-3.980000e+02) * jm[1] * jm[3]));
        case 210: // da2db1Nii2Njj1
        case 101: // da1db2Nii1Njj2
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-7.000000e+01) * jm[0] * jm[2]) + ((1.300000e+01) * jm[1] * jm[2]) + ((4.900000e+01) * jm[0] * jm[3]) + ((-6.000000e+00) * jm[1] * jm[3]));
        case 220: // da2db1Nii3Njj1
        case 102: // da1db2Nii1Njj3
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-6.000000e+00) * jm[0] * jm[2]) + ((4.900000e+01) * jm[1] * jm[2]) + ((1.300000e+01) * jm[0] * jm[3]) + ((-7.000000e+01) * jm[1] * jm[3]));
        case 230: // da2db1Nii4Njj1
        case 103: // da1db2Nii1Njj4
            return ((-1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((1.400000e+02) * jm[0] * jm[2]) + ((8.900000e+01) * jm[1] * jm[2]) + ((-9.100000e+01) * jm[0] * jm[3]) + ((2.000000e+00) * jm[1] * jm[3]));
        case 240: // da2db1Nii5Njj1
        case 104: // da1db2Nii1Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+01) * jm[0] * jm[2]) + ((5.000000e+00) * jm[1] * jm[2]) + ((-1.300000e+01) * jm[0] * jm[3]) + ((2.000000e+00) * jm[1] * jm[3]));
        case 150: // da1db2Nii6Njj1
        case 205: // da2db1Nii1Njj6
        case 250: // da2db1Nii6Njj1
        case 105: // da1db2Nii1Njj6
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0] * jm[2]) + ((3.000000e+00) * jm[1] * jm[2]) + ((3.000000e+00) * jm[0] * jm[3]) + ((-1.000000e+01) * jm[1] * jm[3]));
        case 260: // da2db1Nii7Njj1
        case 106: // da1db2Nii1Njj7
            return ((-1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[0] * jm[2]) + ((-9.100000e+01) * jm[1] * jm[2]) + ((8.900000e+01) * jm[0] * jm[3]) + ((1.400000e+02) * jm[1] * jm[3]));
        case 170: // da1db2Nii8Njj1
        case 207: // da2db1Nii1Njj8
        case 270: // da2db1Nii8Njj1
        case 107: // da1db2Nii1Njj8
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.000000e+01) * jm[0] * jm[2]) + ((3.000000e+00) * jm[1] * jm[2]) + ((3.000000e+00) * jm[0] * jm[3]) + ((4.000000e+00) * jm[1] * jm[3]));
        case 280: // da2db1Nii9Njj1
        case 108: // da1db2Nii1Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[0] * jm[2]) + ((-1.300000e+01) * jm[1] * jm[2]) + ((5.000000e+00) * jm[0] * jm[3]) + ((2.000000e+01) * jm[1] * jm[3]));
        case 190: // da1db2Nii10Njj1
        case 209: // da2db1Nii1Njj10
        case 290: // da2db1Nii10Njj1
        case 109: // da1db2Nii1Njj10
            return ((1.000000e+00) / (4.000000e+01)) * (1.0 / jacobi_det) * (((6.000000e+01) * jm[0] * jm[2]) + ((-3.900000e+01) * jm[1] * jm[2]) + ((-3.900000e+01) * jm[0] * jm[3]) + ((6.000000e+01) * jm[1] * jm[3]));
        case 201: // da2db1Nii1Njj2
        case 110: // da1db2Nii2Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-7.000000e+01) * jm[0] * jm[2]) + ((4.900000e+01) * jm[1] * jm[2]) + ((1.300000e+01) * jm[0] * jm[3]) + ((-6.000000e+00) * jm[1] * jm[3]));
        case 211: // da2db1Nii2Njj2
        case 111: // da1db2Nii2Njj2
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.400000e+01) * jm[0] * jm[2]) + ((3.000000e+00) * jm[1] * jm[2]) + ((3.000000e+00) * jm[0] * jm[3]) + ((-6.000000e+00) * jm[1] * jm[3]));
        case 121: // da1db2Nii3Njj2
        case 212: // da2db1Nii2Njj3
        case 221: // da2db1Nii3Njj2
        case 112: // da1db2Nii2Njj3
            return ((1.000000e+00) / (7.200000e+01)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * jm[3]));
        case 231: // da2db1Nii4Njj2
        case 113: // da1db2Nii2Njj4
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * jm[0] * (((-4.000000e+00) * jm[2]) + jm[3])) + (jm[1] * (((-2.900000e+01) * jm[2]) + ((6.000000e+00) * jm[3]))));
        case 241: // da2db1Nii5Njj2
        case 114: // da1db2Nii2Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0] * jm[2]) + ((2.000000e+00) * jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]));
        case 151: // da1db2Nii6Njj2
        case 215: // da2db1Nii2Njj6
        case 251: // da2db1Nii6Njj2
        case 115: // da1db2Nii2Njj6
            return ((-1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 161: // da1db2Nii7Njj2
        case 216: // da2db1Nii2Njj7
        case 261: // da2db1Nii7Njj2
        case 116: // da1db2Nii2Njj7
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.000000e+01) * jm[0] * jm[2]) + ((7.000000e+00) * jm[1] * jm[2]) + ((7.000000e+00) * jm[0] * jm[3]));
        case 171: // da1db2Nii8Njj2
        case 217: // da2db1Nii2Njj8
        case 271: // da2db1Nii8Njj2
        case 117: // da1db2Nii2Njj8
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 181: // da1db2Nii9Njj2
        case 218: // da2db1Nii2Njj9
        case 281: // da2db1Nii9Njj2
        case 118: // da1db2Nii2Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * (((-1.000000e+00) * jm[2]) + jm[3])));
        case 191: // da1db2Nii10Njj2
        case 219: // da2db1Nii2Njj10
        case 291: // da2db1Nii10Njj2
        case 119: // da1db2Nii2Njj10
            return ((3.000000e+00) / (4.000000e+01)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0] * jm[2]) + ((-1.000000e+00) * jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]));
        case 202: // da2db1Nii1Njj3
        case 120: // da1db2Nii3Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-6.000000e+00) * jm[0] * jm[2]) + ((1.300000e+01) * jm[1] * jm[2]) + ((4.900000e+01) * jm[0] * jm[3]) + ((-7.000000e+01) * jm[1] * jm[3]));
        case 222: // da2db1Nii3Njj3
        case 122: // da1db2Nii3Njj3
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-6.000000e+00) * jm[0] * jm[2]) + ((3.000000e+00) * jm[1] * jm[2]) + ((3.000000e+00) * jm[0] * jm[3]) + ((-1.400000e+01) * jm[1] * jm[3]));
        case 132: // da1db2Nii4Njj3
        case 223: // da2db1Nii3Njj4
        case 232: // da2db1Nii4Njj3
        case 123: // da1db2Nii3Njj4
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * jm[1] * jm[2]) + ((7.000000e+00) * jm[0] * jm[3]) + ((-1.000000e+01) * jm[1] * jm[3]));
        case 142: // da1db2Nii5Njj3
        case 224: // da2db1Nii3Njj5
        case 242: // da2db1Nii5Njj3
        case 124: // da1db2Nii3Njj5
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * jm[3]) + ((-1.000000e+00) * jm[1] * jm[3]));
        case 152: // da1db2Nii6Njj3
        case 225: // da2db1Nii3Njj6
        case 252: // da2db1Nii6Njj3
        case 125: // da1db2Nii3Njj6
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * jm[3]) + ((-2.000000e+00) * jm[1] * jm[3]));
        case 262: // da2db1Nii7Njj3
        case 126: // da1db2Nii3Njj7
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((6.000000e+00) * jm[0] * jm[2]) + ((7.000000e+00) * jm[1] * jm[2]) + ((-2.900000e+01) * jm[0] * jm[3]) + ((-2.800000e+01) * jm[1] * jm[3]));
        case 172: // da1db2Nii8Njj3
        case 227: // da2db1Nii3Njj8
        case 272: // da2db1Nii8Njj3
        case 127: // da1db2Nii3Njj8
            return ((-1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * jm[3]) + ((-2.000000e+00) * jm[1] * jm[3]));
        case 282: // da2db1Nii9Njj3
        case 128: // da1db2Nii3Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((-1.000000e+00) * jm[1] * jm[2]) + ((2.000000e+00) * jm[0] * jm[3]) + ((4.000000e+00) * jm[1] * jm[3]));
        case 192: // da1db2Nii10Njj3
        case 229: // da2db1Nii3Njj10
        case 292: // da2db1Nii10Njj3
        case 129: // da1db2Nii3Njj10
            return ((-3.000000e+00) / (4.000000e+01)) * (1.0 / jacobi_det) * ((jm[1] * (jm[2] + ((-4.000000e+00) * jm[3]))) + (jm[0] * jm[3]));
        case 203: // da2db1Nii1Njj4
        case 130: // da1db2Nii4Njj1
            return ((-1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((1.400000e+02) * jm[0] * jm[2]) + ((-9.100000e+01) * jm[1] * jm[2]) + ((8.900000e+01) * jm[0] * jm[3]) + ((2.000000e+00) * jm[1] * jm[3]));
        case 213: // da2db1Nii2Njj4
        case 131: // da1db2Nii4Njj2
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-2.800000e+01) * jm[0] * jm[2]) + ((7.000000e+00) * jm[1] * jm[2]) + ((-2.900000e+01) * jm[0] * jm[3]) + ((6.000000e+00) * jm[1] * jm[3]));
        case 233: // da2db1Nii4Njj4
        case 133: // da1db2Nii4Njj4
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * ((jm[1] * (((4.900000e+01) * jm[2]) + ((-3.980000e+02) * jm[3]))) + ((4.900000e+01) * jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 143: // da1db2Nii5Njj4
        case 234: // da2db1Nii4Njj5
        case 243: // da2db1Nii5Njj4
        case 134: // da1db2Nii4Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[0] * jm[2]) + ((-7.000000e+00) * jm[1] * jm[2]) + ((-7.000000e+00) * jm[0] * jm[3]) + ((3.800000e+01) * jm[1] * jm[3]));
        case 253: // da2db1Nii6Njj4
        case 135: // da1db2Nii4Njj6
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.400000e+01) * jm[0] * jm[2]) + ((2.100000e+01) * jm[1] * jm[2]) + ((5.700000e+01) * jm[0] * jm[3]) + ((-7.000000e+01) * jm[1] * jm[3]));
        case 263: // da2db1Nii7Njj4
        case 136: // da1db2Nii4Njj7
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.400000e+02) * jm[0] * jm[2]) + ((4.900000e+01) * jm[1] * jm[2]) + ((2.290000e+02) * jm[0] * jm[3]) + ((-1.400000e+02) * jm[1] * jm[3]));
        case 273: // da2db1Nii8Njj4
        case 137: // da1db2Nii4Njj8
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-2.800000e+01) * jm[0] * jm[2]) + ((2.100000e+01) * jm[1] * jm[2]) + ((5.700000e+01) * jm[0] * jm[3]) + ((-4.400000e+01) * jm[1] * jm[3]));
        case 283: // da2db1Nii9Njj4
        case 138: // da1db2Nii4Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[0] * jm[2]) + ((-7.000000e+00) * jm[1] * jm[2]) + ((-2.500000e+01) * jm[0] * jm[3]) + ((2.000000e+01) * jm[1] * jm[3]));
        case 193: // da1db2Nii10Njj4
        case 239: // da2db1Nii4Njj10
        case 293: // da2db1Nii10Njj4
        case 139: // da1db2Nii4Njj10
            return ((3.000000e+00) / (4.000000e+01)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[0] * jm[2]) + ((-7.000000e+00) * jm[1] * jm[2]) + ((-7.000000e+00) * jm[0] * jm[3]) + ((2.000000e+01) * jm[1] * jm[3]));
        case 204: // da2db1Nii1Njj5
        case 140: // da1db2Nii5Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+01) * jm[0] * jm[2]) + ((-1.300000e+01) * jm[1] * jm[2]) + ((5.000000e+00) * jm[0] * jm[3]) + ((2.000000e+00) * jm[1] * jm[3]));
        case 214: // da2db1Nii2Njj5
        case 141: // da1db2Nii5Njj2
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0] * jm[2]) + ((-1.000000e+00) * jm[1] * jm[2]) + ((2.000000e+00) * jm[0] * jm[3]));
        case 244: // da2db1Nii5Njj5
        case 144: // da1db2Nii5Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * ((jm[1] * (jm[2] + ((-5.000000e+00) * jm[3]))) + (jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 154: // da1db2Nii6Njj5
        case 245: // da2db1Nii5Njj6
        case 254: // da2db1Nii6Njj5
        case 145: // da1db2Nii5Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[0] * jm[2]) + ((-3.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3]) + ((7.000000e+00) * jm[1] * jm[3]));
        case 264: // da2db1Nii7Njj5
        case 146: // da1db2Nii5Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+01) * jm[0] * jm[2]) + ((-7.000000e+00) * jm[1] * jm[2]) + ((-2.500000e+01) * jm[0] * jm[3]) + ((1.400000e+01) * jm[1] * jm[3]));
        case 274: // da2db1Nii8Njj5
        case 147: // da1db2Nii5Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0] * jm[2]) + ((-3.000000e+00) * jm[1] * jm[2]) + ((-6.000000e+00) * jm[0] * jm[3]) + ((5.000000e+00) * jm[1] * jm[3]));
        case 284: // da2db1Nii9Njj5
        case 148: // da1db2Nii5Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((-4.000000e+00) * jm[0] * jm[2]) + ((2.000000e+00) * jm[1] * jm[2]) + ((5.000000e+00) * jm[0] * jm[3]) + ((-4.000000e+00) * jm[1] * jm[3]));
        case 189: // da1db2Nii9Njj10
        case 194: // da1db2Nii10Njj5
        case 198: // da1db2Nii10Njj9
        case 249: // da2db1Nii5Njj10
        case 289: // da2db1Nii9Njj10
        case 294: // da2db1Nii10Njj5
        case 298: // da2db1Nii10Njj9
        case 149: // da1db2Nii5Njj10
            return ((3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * ((jm[1] * (jm[2] + ((-2.000000e+00) * jm[3]))) + (jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 235: // da2db1Nii4Njj6
        case 153: // da1db2Nii6Njj4
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.400000e+01) * jm[0] * jm[2]) + ((5.700000e+01) * jm[1] * jm[2]) + ((2.100000e+01) * jm[0] * jm[3]) + ((-7.000000e+01) * jm[1] * jm[3]));
        case 177: // da1db2Nii8Njj8
        case 255: // da2db1Nii6Njj6
        case 277: // da2db1Nii8Njj8
        case 155: // da1db2Nii6Njj6
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.400000e+01) * jm[0] * jm[2]) + ((1.100000e+01) * jm[1] * jm[2]) + ((1.100000e+01) * jm[0] * jm[3]) + ((-1.400000e+01) * jm[1] * jm[3]));
        case 265: // da2db1Nii7Njj6
        case 156: // da1db2Nii6Njj7
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-4.400000e+01) * jm[0] * jm[2]) + ((2.100000e+01) * jm[1] * jm[2]) + ((5.700000e+01) * jm[0] * jm[3]) + ((-2.800000e+01) * jm[1] * jm[3]));
        case 275: // da2db1Nii8Njj6
        case 157: // da1db2Nii6Njj8
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.000000e+01) * jm[0] * jm[2]) + ((7.000000e+00) * jm[1] * jm[2]) + ((1.300000e+01) * jm[0] * jm[3]) + ((-1.000000e+01) * jm[1] * jm[3]));
        case 285: // da2db1Nii9Njj6
        case 158: // da1db2Nii6Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((5.000000e+00) * jm[0] * jm[2]) + ((-3.000000e+00) * jm[1] * jm[2]) + ((-6.000000e+00) * jm[0] * jm[3]) + ((4.000000e+00) * jm[1] * jm[3]));
        case 195: // da1db2Nii10Njj6
        case 259: // da2db1Nii6Njj10
        case 295: // da2db1Nii10Njj6
        case 159: // da1db2Nii6Njj10
            return ((3.000000e+00) / (4.000000e+01)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[0] * jm[2]) + ((-3.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3]) + ((4.000000e+00) * jm[1] * jm[3]));
        case 206: // da2db1Nii1Njj7
        case 160: // da1db2Nii7Njj1
            return ((-1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[0] * jm[2]) + ((8.900000e+01) * jm[1] * jm[2]) + ((-9.100000e+01) * jm[0] * jm[3]) + ((1.400000e+02) * jm[1] * jm[3]));
        case 226: // da2db1Nii3Njj7
        case 162: // da1db2Nii7Njj3
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((6.000000e+00) * jm[0] * jm[2]) + ((-2.900000e+01) * jm[1] * jm[2]) + ((7.000000e+00) * jm[0] * jm[3]) + ((-2.800000e+01) * jm[1] * jm[3]));
        case 236: // da2db1Nii4Njj7
        case 163: // da1db2Nii7Njj4
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.400000e+02) * jm[0] * jm[2]) + ((2.290000e+02) * jm[1] * jm[2]) + ((4.900000e+01) * jm[0] * jm[3]) + ((-1.400000e+02) * jm[1] * jm[3]));
        case 246: // da2db1Nii5Njj7
        case 164: // da1db2Nii7Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+01) * jm[0] * jm[2]) + ((-2.500000e+01) * jm[1] * jm[2]) + ((-7.000000e+00) * jm[0] * jm[3]) + ((1.400000e+01) * jm[1] * jm[3]));
        case 256: // da2db1Nii6Njj7
        case 165: // da1db2Nii7Njj6
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-4.400000e+01) * jm[0] * jm[2]) + ((5.700000e+01) * jm[1] * jm[2]) + ((2.100000e+01) * jm[0] * jm[3]) + ((-2.800000e+01) * jm[1] * jm[3]));
        case 266: // da2db1Nii7Njj7
        case 166: // da1db2Nii7Njj7
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((4.900000e+01) * jm[1] * (jm[2] + ((-2.000000e+00) * jm[3]))) + (jm[0] * (((-3.980000e+02) * jm[2]) + ((4.900000e+01) * jm[3]))));
        case 276: // da2db1Nii8Njj7
        case 167: // da1db2Nii7Njj8
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-7.000000e+01) * jm[0] * jm[2]) + ((5.700000e+01) * jm[1] * jm[2]) + ((2.100000e+01) * jm[0] * jm[3]) + ((-1.400000e+01) * jm[1] * jm[3]));
        case 186: // da1db2Nii9Njj7
        case 268: // da2db1Nii7Njj9
        case 286: // da2db1Nii9Njj7
        case 168: // da1db2Nii7Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.800000e+01) * jm[0] * jm[2]) + ((-7.000000e+00) * jm[1] * jm[2]) + ((-7.000000e+00) * jm[0] * jm[3]) + ((1.400000e+01) * jm[1] * jm[3]));
        case 196: // da1db2Nii10Njj7
        case 269: // da2db1Nii7Njj10
        case 296: // da2db1Nii10Njj7
        case 169: // da1db2Nii7Njj10
            return ((3.000000e+00) / (4.000000e+01)) * (1.0 / jacobi_det) * (((2.000000e+01) * jm[0] * jm[2]) + ((-7.000000e+00) * jm[1] * jm[2]) + ((-7.000000e+00) * jm[0] * jm[3]) + ((1.400000e+01) * jm[1] * jm[3]));
        case 237: // da2db1Nii4Njj8
        case 173: // da1db2Nii8Njj4
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-2.800000e+01) * jm[0] * jm[2]) + ((5.700000e+01) * jm[1] * jm[2]) + ((2.100000e+01) * jm[0] * jm[3]) + ((-4.400000e+01) * jm[1] * jm[3]));
        case 247: // da2db1Nii5Njj8
        case 174: // da1db2Nii8Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0] * jm[2]) + ((-6.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3]) + ((5.000000e+00) * jm[1] * jm[3]));
        case 257: // da2db1Nii6Njj8
        case 175: // da1db2Nii8Njj6
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-1.000000e+01) * jm[0] * jm[2]) + ((1.300000e+01) * jm[1] * jm[2]) + ((7.000000e+00) * jm[0] * jm[3]) + ((-1.000000e+01) * jm[1] * jm[3]));
        case 267: // da2db1Nii7Njj8
        case 176: // da1db2Nii8Njj7
            return ((1.000000e+00) / (3.600000e+02)) * (1.0 / jacobi_det) * (((-7.000000e+01) * jm[0] * jm[2]) + ((2.100000e+01) * jm[1] * jm[2]) + ((5.700000e+01) * jm[0] * jm[3]) + ((-1.400000e+01) * jm[1] * jm[3]));
        case 187: // da1db2Nii9Njj8
        case 278: // da2db1Nii8Njj9
        case 287: // da2db1Nii9Njj8
        case 178: // da1db2Nii8Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * jm[0] * jm[2]) + ((-3.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3]) + ((2.000000e+00) * jm[1] * jm[3]));
        case 197: // da1db2Nii10Njj8
        case 279: // da2db1Nii8Njj10
        case 297: // da2db1Nii10Njj8
        case 179: // da1db2Nii8Njj10
            return ((3.000000e+00) / (4.000000e+01)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0] * jm[2]) + ((-3.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3]) + ((2.000000e+00) * jm[1] * jm[3]));
        case 208: // da2db1Nii1Njj9
        case 180: // da1db2Nii9Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[0] * jm[2]) + ((5.000000e+00) * jm[1] * jm[2]) + ((-1.300000e+01) * jm[0] * jm[3]) + ((2.000000e+01) * jm[1] * jm[3]));
        case 228: // da2db1Nii3Njj9
        case 182: // da1db2Nii9Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]) + ((4.000000e+00) * jm[1] * jm[3]));
        case 238: // da2db1Nii4Njj9
        case 183: // da1db2Nii9Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[0] * jm[2]) + ((-2.500000e+01) * jm[1] * jm[2]) + ((-7.000000e+00) * jm[0] * jm[3]) + ((2.000000e+01) * jm[1] * jm[3]));
        case 248: // da2db1Nii5Njj9
        case 184: // da1db2Nii9Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((-4.000000e+00) * jm[0] * jm[2]) + ((5.000000e+00) * jm[1] * jm[2]) + ((2.000000e+00) * jm[0] * jm[3]) + ((-4.000000e+00) * jm[1] * jm[3]));
        case 258: // da2db1Nii6Njj9
        case 185: // da1db2Nii9Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((5.000000e+00) * jm[0] * jm[2]) + ((-6.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3]) + ((4.000000e+00) * jm[1] * jm[3]));
        case 288: // da2db1Nii9Njj9
        case 188: // da1db2Nii9Njj9
            return ((1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * ((jm[1] * (jm[2] + ((-2.000000e+00) * jm[3]))) + (jm[0] * (((-5.000000e+00) * jm[2]) + jm[3])));
        case 299: // da2db1Nii10Njj10
        case 199: // da1db2Nii10Njj10
            return ((8.100000e+01) / (4.000000e+01)) * (1.0 / jacobi_det) * ((jm[1] * (jm[2] + ((-2.000000e+00) * jm[3]))) + (jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 300: // da2db2Nii1Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.990000e+02) * (jm[0] * jm[0])) + ((-3.490000e+02) * jm[0] * jm[1]) + ((1.990000e+02) * (jm[1] * jm[1])));
        case 310: // da2db2Nii2Njj1
        case 301: // da2db2Nii1Njj2
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.500000e+01) * (jm[0] * jm[0])) + ((-3.100000e+01) * jm[0] * jm[1]) + ((3.000000e+00) * (jm[1] * jm[1])));
        case 320: // da2db2Nii3Njj1
        case 302: // da2db2Nii1Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.000000e+00) * (jm[0] * jm[0])) + ((-3.100000e+01) * jm[0] * jm[1]) + ((3.500000e+01) * (jm[1] * jm[1])));
        case 330: // da2db2Nii4Njj1
        case 303: // da2db2Nii1Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+01) * (jm[0] * jm[0])) + ((-1.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1]));
        case 340: // da2db2Nii5Njj1
        case 304: // da2db2Nii1Njj5
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((1.000000e+01) * (jm[0] * jm[0])) + ((-4.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1]));
        case 350: // da2db2Nii6Njj1
        case 305: // da2db2Nii1Njj6
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * (((2.000000e+00) * jm[0]) + ((5.000000e+00) * jm[1]));
        case 360: // da2db2Nii7Njj1
        case 306: // da2db2Nii1Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * ((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[1]) + ((7.000000e+01) * (jm[1] * jm[1])));
        case 370: // da2db2Nii8Njj1
        case 307: // da2db2Nii1Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * (((5.000000e+00) * jm[0]) + ((2.000000e+00) * jm[1]));
        case 380: // da2db2Nii9Njj1
        case 308: // da2db2Nii1Njj9
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * ((jm[0] * jm[0]) + ((-4.000000e+00) * jm[0] * jm[1]) + ((1.000000e+01) * (jm[1] * jm[1])));
        case 390: // da2db2Nii10Njj1
        case 309: // da2db2Nii1Njj10
            return ((-1.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((3.000000e+01) * (jm[0] * jm[0])) + ((-3.900000e+01) * jm[0] * jm[1]) + ((3.000000e+01) * (jm[1] * jm[1])));
        case 311: // da2db2Nii2Njj2
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[0] * jm[0])) + ((-3.000000e+00) * jm[0] * jm[1]) + ((3.000000e+00) * (jm[1] * jm[1])));
        case 321: // da2db2Nii3Njj2
        case 312: // da2db2Nii2Njj3
            return ((-1.000000e+00) / (3.600000e+01)) * (1.0 / jacobi_det) * jm[0] * jm[1];
        case 331: // da2db2Nii4Njj2
        case 313: // da2db2Nii2Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[0]) + ((-3.000000e+00) * jm[1])) * (jm[0] + jm[1]);
        case 341: // da2db2Nii5Njj2
        case 314: // da2db2Nii2Njj5
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[0] * (((4.000000e+00) * jm[0]) + jm[1]);
        case 351: // da2db2Nii6Njj2
        case 315: // da2db2Nii2Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[0] * (((-1.000000e+00) * jm[0]) + jm[1]);
        case 361: // da2db2Nii7Njj2
        case 316: // da2db2Nii2Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[0] * (((5.000000e+00) * jm[0]) + ((-7.000000e+00) * jm[1]));
        case 371: // da2db2Nii8Njj2
        case 317: // da2db2Nii2Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[0] * (jm[0] + ((-1.000000e+00) * jm[1]));
        case 381: // da2db2Nii9Njj2
        case 318: // da2db2Nii2Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[0] * (jm[0] + ((-2.000000e+00) * jm[1]));
        case 391: // da2db2Nii10Njj2
        case 319: // da2db2Nii2Njj10
            return ((3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * jm[0] * (((-2.000000e+00) * jm[0]) + jm[1]);
        case 322: // da2db2Nii3Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.000000e+00) * (jm[0] * jm[0])) + ((-3.000000e+00) * jm[0] * jm[1]) + ((7.000000e+00) * (jm[1] * jm[1])));
        case 332: // da2db2Nii4Njj3
        case 323: // da2db2Nii3Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[1] * (((-7.000000e+00) * jm[0]) + ((5.000000e+00) * jm[1]));
        case 342: // da2db2Nii5Njj3
        case 324: // da2db2Nii3Njj5
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[1] * (((-2.000000e+00) * jm[0]) + jm[1]);
        case 352: // da2db2Nii6Njj3
        case 325: // da2db2Nii3Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[1] * (((-1.000000e+00) * jm[0]) + jm[1]);
        case 362: // da2db2Nii7Njj3
        case 326: // da2db2Nii3Njj7
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.000000e+00) * jm[0]) + ((-1.400000e+01) * jm[1])) * (jm[0] + jm[1]);
        case 372: // da2db2Nii8Njj3
        case 327: // da2db2Nii3Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[1];
        case 382: // da2db2Nii9Njj3
        case 328: // da2db2Nii3Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * jm[1] * (jm[0] + ((4.000000e+00) * jm[1]));
        case 392: // da2db2Nii10Njj3
        case 329: // da2db2Nii3Njj10
            return ((3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (jm[0] + ((-2.000000e+00) * jm[1])) * jm[1];
        case 333: // da2db2Nii4Njj4
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.900000e+01) * (jm[0] * jm[0])) + ((-4.900000e+01) * jm[0] * jm[1]) + ((1.990000e+02) * (jm[1] * jm[1])));
        case 343: // da2db2Nii5Njj4
        case 334: // da2db2Nii4Njj5
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[0] * jm[0])) + ((-7.000000e+00) * jm[0] * jm[1]) + ((1.900000e+01) * (jm[1] * jm[1])));
        case 353: // da2db2Nii6Njj4
        case 335: // da2db2Nii4Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[0] * jm[0])) + ((-3.900000e+01) * jm[0] * jm[1]) + ((3.500000e+01) * (jm[1] * jm[1])));
        case 363: // da2db2Nii7Njj4
        case 336: // da2db2Nii4Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+01) * (jm[0] * jm[0])) + ((-1.390000e+02) * jm[0] * jm[1]) + ((7.000000e+01) * (jm[1] * jm[1])));
        case 373: // da2db2Nii8Njj4
        case 337: // da2db2Nii4Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.400000e+01) * jm[0]) + ((-1.100000e+01) * jm[1])) * (jm[0] + ((-2.000000e+00) * jm[1]));
        case 383: // da2db2Nii9Njj4
        case 338: // da2db2Nii4Njj9
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[0] * jm[0])) + ((-1.600000e+01) * jm[0] * jm[1]) + ((1.000000e+01) * (jm[1] * jm[1])));
        case 393: // da2db2Nii10Njj4
        case 339: // da2db2Nii4Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[0] * jm[0])) + ((-7.000000e+00) * jm[0] * jm[1]) + ((1.000000e+01) * (jm[1] * jm[1])));
        case 344: // da2db2Nii5Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((2.000000e+00) * (jm[0] * jm[0])) + ((-2.000000e+00) * jm[0] * jm[1]) + ((5.000000e+00) * (jm[1] * jm[1])));
        case 354: // da2db2Nii6Njj5
        case 345: // da2db2Nii5Njj6
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.000000e+00) * (jm[0] * jm[0])) + ((-6.000000e+00) * jm[0] * jm[1]) + ((7.000000e+00) * (jm[1] * jm[1])));
        case 364: // da2db2Nii7Njj5
        case 346: // da2db2Nii5Njj7
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((1.000000e+01) * (jm[0] * jm[0])) + ((-1.600000e+01) * jm[0] * jm[1]) + ((7.000000e+00) * (jm[1] * jm[1])));
        case 374: // da2db2Nii8Njj5
        case 347: // da2db2Nii5Njj8
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * jm[0]) + ((-5.000000e+00) * jm[1])) * (jm[0] + ((-1.000000e+00) * jm[1]));
        case 384: // da2db2Nii9Njj5
        case 348: // da2db2Nii5Njj9
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((4.000000e+00) * (jm[0] * jm[0])) + ((-7.000000e+00) * jm[0] * jm[1]) + ((4.000000e+00) * (jm[1] * jm[1])));
        case 389: // da2db2Nii9Njj10
        case 394: // da2db2Nii10Njj5
        case 398: // da2db2Nii10Njj9
        case 349: // da2db2Nii5Njj10
            return ((3.000000e+00) / (1.000000e+01)) * (1.0 / jacobi_det) * ((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1]));
        case 377: // da2db2Nii8Njj8
        case 355: // da2db2Nii6Njj6
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[0] * jm[0])) + ((-1.100000e+01) * jm[0] * jm[1]) + ((7.000000e+00) * (jm[1] * jm[1])));
        case 365: // da2db2Nii7Njj6
        case 356: // da2db2Nii6Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((2.200000e+01) * (jm[0] * jm[0])) + ((-3.900000e+01) * jm[0] * jm[1]) + ((1.400000e+01) * (jm[1] * jm[1])));
        case 375: // da2db2Nii8Njj6
        case 357: // da2db2Nii6Njj8
            return ((1.000000e+00) / (3.600000e+01)) * (1.0 / jacobi_det) * ((jm[0] + ((-1.000000e+00) * jm[1])) * (jm[0] + ((-1.000000e+00) * jm[1])));
        case 385: // da2db2Nii9Njj6
        case 358: // da2db2Nii6Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((5.000000e+00) * jm[0]) + ((-4.000000e+00) * jm[1])) * (jm[0] + ((-1.000000e+00) * jm[1]));
        case 395: // da2db2Nii10Njj6
        case 359: // da2db2Nii6Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (jm[0] + ((-2.000000e+00) * jm[1])) * (jm[0] + ((-1.000000e+00) * jm[1]));
        case 366: // da2db2Nii7Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((1.990000e+02) * (jm[0] * jm[0])) + ((-4.900000e+01) * jm[0] * jm[1]) + ((4.900000e+01) * (jm[1] * jm[1])));
        case 376: // da2db2Nii8Njj7
        case 367: // da2db2Nii7Njj8
            return ((1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((3.500000e+01) * (jm[0] * jm[0])) + ((-3.900000e+01) * jm[0] * jm[1]) + ((7.000000e+00) * (jm[1] * jm[1])));
        case 386: // da2db2Nii9Njj7
        case 368: // da2db2Nii7Njj9
            return ((-1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((1.900000e+01) * (jm[0] * jm[0])) + ((-7.000000e+00) * jm[0] * jm[1]) + ((7.000000e+00) * (jm[1] * jm[1])));
        case 396: // da2db2Nii10Njj7
        case 369: // da2db2Nii7Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((1.000000e+01) * (jm[0] * jm[0])) + ((-7.000000e+00) * jm[0] * jm[1]) + ((7.000000e+00) * (jm[1] * jm[1])));
        case 387: // da2db2Nii9Njj8
        case 378: // da2db2Nii8Njj9
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / jacobi_det) * (((7.000000e+00) * (jm[0] * jm[0])) + ((-6.000000e+00) * jm[0] * jm[1]) + ((2.000000e+00) * (jm[1] * jm[1])));
        case 397: // da2db2Nii10Njj8
        case 379: // da2db2Nii8Njj10
            return ((-3.000000e+00) / (2.000000e+01)) * (1.0 / jacobi_det) * (((2.000000e+00) * (jm[0] * jm[0])) + ((-3.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1]));
        case 388: // da2db2Nii9Njj9
            return ((1.000000e+00) / (9.000000e+01)) * (1.0 / jacobi_det) * (((5.000000e+00) * (jm[0] * jm[0])) + ((-2.000000e+00) * jm[0] * jm[1]) + ((2.000000e+00) * (jm[1] * jm[1])));
        case 399: // da2db2Nii10Njj10
            return ((8.100000e+01) / (2.000000e+01)) * (1.0 / jacobi_det) * ((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1]));
		
		default:
			TDKP_GENERAL_EXCEPTION("invalid object index");	
	}	
}

double Element2DTriangleHermiteCubic::get_element_integral_1st_order(
	short diffop, 
	short elem_shape_func_1, 
	short elem_shape_func_2
) const {
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1
	                       + 10 * 10 * diffop;
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
        case 0: // da1Nii1Njj1
            return ((1.300000e+01) / (7.000000e+01)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 1: // da1Nii1Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.800000e+01) * jm[2]) + ((-3.700000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 2: // da1Nii1Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.700000e+01) * jm[2]) + ((-5.800000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 3: // da1Nii1Njj4
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((1.300000e+01) * jm[2]) + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 4: // da1Nii1Njj5
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (((6.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 5: // da1Nii1Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.500000e+01) * jm[2]) + ((-2.200000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 6: // da1Nii1Njj7
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (jm[2] + ((1.300000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 7: // da1Nii1Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[2]) + ((-1.500000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 8: // da1Nii1Njj9
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (jm[2] + ((6.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 9: // da1Nii1Njj10
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 10: // da1Nii2Njj1
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[2]) + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 11: // da1Nii2Njj2
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 12: // da1Nii2Njj3
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (jm[2] + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 13: // da1Nii2Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+01) * jm[2]) + ((-1.700000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 14: // da1Nii2Njj5
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.000000e+01) * jm[2]) + ((-3.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 27: // da1Nii3Njj8
        case 15: // da1Nii2Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 16: // da1Nii2Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[2]) + ((-7.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 17: // da1Nii2Njj8
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((-2.000000e+00) * jm[2]) + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 18: // da1Nii2Njj9
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+00) * jm[2]) + ((-2.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 19: // da1Nii2Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((-4.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 20: // da1Nii3Njj1
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (jm[2] + ((2.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 21: // da1Nii3Njj2
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (jm[2] + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 22: // da1Nii3Njj3
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * jm[3] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 23: // da1Nii3Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[2]) + ((-2.200000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 24: // da1Nii3Njj5
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[2]) + ((-5.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 25: // da1Nii3Njj6
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (jm[2] + ((-2.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 26: // da1Nii3Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.700000e+01) * jm[2]) + ((-5.000000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 28: // da1Nii3Njj9
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.000000e+00) * jm[2]) + ((-1.000000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 29: // da1Nii3Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (jm[2] + ((-4.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 30: // da1Nii4Njj1
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((-1.400000e+01) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 31: // da1Nii4Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.800000e+01) * jm[2]) + ((1.700000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 32: // da1Nii4Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[2]) + ((-2.200000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 33: // da1Nii4Njj4
            return ((1.300000e+01) / (7.000000e+01)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 34: // da1Nii4Njj5
            return ((-1.900000e+01) / (5.040000e+02)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 35: // da1Nii4Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.100000e+01) * jm[2]) + ((-5.800000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 36: // da1Nii4Njj7
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((1.400000e+01) * jm[2]) + ((-1.300000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 37: // da1Nii4Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.800000e+01) * jm[2]) + ((-4.500000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 38: // da1Nii4Njj9
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (((7.000000e+00) * jm[2]) + ((-6.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 39: // da1Nii4Njj10
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 40: // da1Nii5Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.800000e+01) * jm[2]) + ((-5.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 41: // da1Nii5Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((8.000000e+00) * jm[2]) + ((3.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 42: // da1Nii5Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[2]) + ((-5.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 43: // da1Nii5Njj4
            return ((-3.700000e+01) / (2.520000e+03)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 44: // da1Nii5Njj5
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 45: // da1Nii5Njj6
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((2.000000e+00) * jm[2]) + ((-3.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 46: // da1Nii5Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.800000e+01) * jm[2]) + ((-4.300000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 47: // da1Nii5Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((8.000000e+00) * jm[2]) + ((-1.100000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 48: // da1Nii5Njj9
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((5.000000e+00) * jm[2]) + ((-4.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 49: // da1Nii5Njj10
            return ((3.000000e+00) / (5.600000e+01)) * jacobi_det * jm[3] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 50: // da1Nii6Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.500000e+01) * jm[2]) + ((-2.200000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 72: // da1Nii8Njj3
        case 51: // da1Nii6Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((-1.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 52: // da1Nii6Njj3
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (jm[2] + ((-2.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 53: // da1Nii6Njj4
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (((3.000000e+00) * jm[2]) + ((-2.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 54: // da1Nii6Njj5
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((6.000000e+00) * jm[2]) + ((-5.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 77: // da1Nii8Njj8
        case 55: // da1Nii6Njj6
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 56: // da1Nii6Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.300000e+01) * jm[2]) + ((-5.000000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 57: // da1Nii6Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[2]) + ((-1.100000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 58: // da1Nii6Njj9
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[2]) + ((-1.000000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 59: // da1Nii6Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((3.000000e+00) * jm[2]) + ((-4.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 60: // da1Nii7Njj1
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (jm[2] + ((-1.400000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 61: // da1Nii7Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[2]) + ((-7.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 62: // da1Nii7Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.700000e+01) * jm[2]) + ((2.800000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 63: // da1Nii7Njj4
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((1.300000e+01) * jm[2]) + ((-1.400000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 64: // da1Nii7Njj5
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (((6.000000e+00) * jm[2]) + ((-7.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 65: // da1Nii7Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.500000e+01) * jm[2]) + ((-2.800000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 66: // da1Nii7Njj7
            return (1.300000e+01) * jacobi_det * jm[2] * (1.0 / (((7.000000e+01) * jm[1] * jm[2]) + ((-7.000000e+01) * jm[0] * jm[3])));
        case 67: // da1Nii7Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.800000e+01) * jm[2]) + ((-2.100000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 68: // da1Nii7Njj9
            return ((1.900000e+01) / (5.040000e+02)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 69: // da1Nii7Njj10
            return (9.000000e+00) * jacobi_det * jm[2] * (1.0 / (((3.500000e+01) * jm[1] * jm[2]) + ((-3.500000e+01) * jm[0] * jm[3])));
        case 70: // da1Nii8Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[2]) + ((-1.500000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 71: // da1Nii8Njj2
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((-2.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 73: // da1Nii8Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+01) * jm[2]) + ((-3.300000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 74: // da1Nii8Njj5
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.000000e+01) * jm[2]) + ((-7.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 75: // da1Nii8Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.100000e+01) * jm[2]) + ((-7.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 76: // da1Nii8Njj7
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[2]) + ((-3.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 78: // da1Nii8Njj9
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((5.000000e+00) * jm[2]) + ((-6.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 79: // da1Nii8Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((4.000000e+00) * jm[2]) + ((-3.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 80: // da1Nii9Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+00) * jm[2]) + ((-4.800000e+01) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 81: // da1Nii9Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+00) * jm[2]) + ((-2.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 82: // da1Nii9Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.000000e+00) * jm[2]) + ((8.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 83: // da1Nii9Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.300000e+01) * jm[2]) + ((-4.800000e+01) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 84: // da1Nii9Njj5
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((4.000000e+00) * jm[2]) + ((-5.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 85: // da1Nii9Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.100000e+01) * jm[2]) + ((-8.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 86: // da1Nii9Njj7
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 87: // da1Nii9Njj8
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((3.000000e+00) * jm[2]) + ((-2.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 88: // da1Nii9Njj9
            return jacobi_det * jm[2] * (1.0 / (((2.100000e+02) * jm[1] * jm[2]) + ((-2.100000e+02) * jm[0] * jm[3])));
        case 89: // da1Nii9Njj10
            return ((3.000000e+00) / (5.600000e+01)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 90: // da1Nii10Njj1
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 91: // da1Nii10Njj2
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((-4.000000e+00) * jm[2]) + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 92: // da1Nii10Njj3
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (jm[2] + ((-4.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 93: // da1Nii10Njj4
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * jm[3] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 94: // da1Nii10Njj5
            return ((3.000000e+00) / (5.600000e+01)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 95: // da1Nii10Njj6
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((3.000000e+00) * jm[2]) + ((-4.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 96: // da1Nii10Njj7
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 97: // da1Nii10Njj8
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((4.000000e+00) * jm[2]) + ((-3.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 98: // da1Nii10Njj9
            return (3.000000e+00) * jacobi_det * jm[2] * (1.0 / (((5.600000e+01) * jm[1] * jm[2]) + ((-5.600000e+01) * jm[0] * jm[3])));
        case 100: // da2Nii1Njj1
            return ((1.300000e+01) / (7.000000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 101: // da2Nii1Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.800000e+01) * jm[0]) + ((-3.700000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 102: // da2Nii1Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.700000e+01) * jm[0]) + ((-5.800000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 103: // da2Nii1Njj4
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((1.300000e+01) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 104: // da2Nii1Njj5
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (((6.000000e+00) * jm[0]) + jm[1]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 105: // da2Nii1Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.500000e+01) * jm[0]) + ((-2.200000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 106: // da2Nii1Njj7
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (jm[0] + ((1.300000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 107: // da2Nii1Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[0]) + ((-1.500000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 108: // da2Nii1Njj9
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (jm[0] + ((6.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 109: // da2Nii1Njj10
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 110: // da2Nii2Njj1
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 111: // da2Nii2Njj2
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * jm[0] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 112: // da2Nii2Njj3
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (jm[0] + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 113: // da2Nii2Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+01) * jm[0]) + ((-1.700000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 114: // da2Nii2Njj5
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.000000e+01) * jm[0]) + ((-3.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 127: // da2Nii3Njj8
        case 115: // da2Nii2Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 116: // da2Nii2Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[0]) + ((-7.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 117: // da2Nii2Njj8
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((-2.000000e+00) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 118: // da2Nii2Njj9
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+00) * jm[0]) + ((-2.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 119: // da2Nii2Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((-4.000000e+00) * jm[0]) + jm[1]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 120: // da2Nii3Njj1
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (jm[0] + ((2.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 121: // da2Nii3Njj2
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (jm[0] + jm[1]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 122: // da2Nii3Njj3
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 123: // da2Nii3Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[0]) + ((-2.200000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 124: // da2Nii3Njj5
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[0]) + ((-5.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 125: // da2Nii3Njj6
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (jm[0] + ((-2.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 126: // da2Nii3Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.700000e+01) * jm[0]) + ((-5.000000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 128: // da2Nii3Njj9
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.000000e+00) * jm[0]) + ((-1.000000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 129: // da2Nii3Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (jm[0] + ((-4.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 130: // da2Nii4Njj1
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((-1.400000e+01) * jm[0]) + jm[1]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 131: // da2Nii4Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.800000e+01) * jm[0]) + ((1.700000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 132: // da2Nii4Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[0]) + ((-2.200000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 133: // da2Nii4Njj4
            return (1.300000e+01) * jacobi_det * jm[1] * (1.0 / (((7.000000e+01) * jm[1] * jm[2]) + ((-7.000000e+01) * jm[0] * jm[3])));
        case 134: // da2Nii4Njj5
            return ((1.900000e+01) / (5.040000e+02)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 135: // da2Nii4Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.100000e+01) * jm[0]) + ((-5.800000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 136: // da2Nii4Njj7
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((1.400000e+01) * jm[0]) + ((-1.300000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 137: // da2Nii4Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.800000e+01) * jm[0]) + ((-4.500000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 138: // da2Nii4Njj9
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (((7.000000e+00) * jm[0]) + ((-6.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 139: // da2Nii4Njj10
            return (9.000000e+00) * jacobi_det * jm[1] * (1.0 / (((3.500000e+01) * jm[1] * jm[2]) + ((-3.500000e+01) * jm[0] * jm[3])));
        case 140: // da2Nii5Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.800000e+01) * jm[0]) + ((-5.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 141: // da2Nii5Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((8.000000e+00) * jm[0]) + ((3.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 142: // da2Nii5Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[0]) + ((-5.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 143: // da2Nii5Njj4
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 144: // da2Nii5Njj5
            return jacobi_det * jm[1] * (1.0 / (((2.100000e+02) * jm[1] * jm[2]) + ((-2.100000e+02) * jm[0] * jm[3])));
        case 145: // da2Nii5Njj6
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((2.000000e+00) * jm[0]) + ((-3.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 146: // da2Nii5Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.800000e+01) * jm[0]) + ((-4.300000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 147: // da2Nii5Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((8.000000e+00) * jm[0]) + ((-1.100000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 148: // da2Nii5Njj9
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((5.000000e+00) * jm[0]) + ((-4.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 149: // da2Nii5Njj10
            return ((3.000000e+00) / (5.600000e+01)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 150: // da2Nii6Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.500000e+01) * jm[0]) + ((-2.200000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 172: // da2Nii8Njj3
        case 151: // da2Nii6Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 152: // da2Nii6Njj3
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (jm[0] + ((-2.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 153: // da2Nii6Njj4
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (((3.000000e+00) * jm[0]) + ((-2.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 154: // da2Nii6Njj5
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((6.000000e+00) * jm[0]) + ((-5.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 177: // da2Nii8Njj8
        case 155: // da2Nii6Njj6
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 156: // da2Nii6Njj7
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.300000e+01) * jm[0]) + ((-5.000000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 157: // da2Nii6Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[0]) + ((-1.100000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 158: // da2Nii6Njj9
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((7.000000e+00) * jm[0]) + ((-1.000000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 159: // da2Nii6Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((3.000000e+00) * jm[0]) + ((-4.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 160: // da2Nii7Njj1
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (jm[0] + ((-1.400000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 161: // da2Nii7Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[0]) + ((-7.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 162: // da2Nii7Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.700000e+01) * jm[0]) + ((2.800000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 163: // da2Nii7Njj4
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * (((1.300000e+01) * jm[0]) + ((-1.400000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 164: // da2Nii7Njj5
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det * (((6.000000e+00) * jm[0]) + ((-7.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 165: // da2Nii7Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.500000e+01) * jm[0]) + ((-2.800000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 166: // da2Nii7Njj7
            return ((1.300000e+01) / (7.000000e+01)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 167: // da2Nii7Njj8
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.800000e+01) * jm[0]) + ((-2.100000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 168: // da2Nii7Njj9
            return ((-1.900000e+01) / (5.040000e+02)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 169: // da2Nii7Njj10
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 170: // da2Nii8Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((2.200000e+01) * jm[0]) + ((-1.500000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 171: // da2Nii8Njj2
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((-2.000000e+00) * jm[0]) + jm[1]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 173: // da2Nii8Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+01) * jm[0]) + ((-3.300000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 174: // da2Nii8Njj5
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.000000e+01) * jm[0]) + ((-7.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 175: // da2Nii8Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.100000e+01) * jm[0]) + ((-7.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 176: // da2Nii8Njj7
            return ((3.700000e+01) / (2.520000e+03)) * jacobi_det * (((2.000000e+00) * jm[0]) + ((-3.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 178: // da2Nii8Njj9
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((5.000000e+00) * jm[0]) + ((-6.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 179: // da2Nii8Njj10
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((4.000000e+00) * jm[0]) + ((-3.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 180: // da2Nii9Njj1
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+00) * jm[0]) + ((-4.800000e+01) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 181: // da2Nii9Njj2
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((5.000000e+00) * jm[0]) + ((-2.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 182: // da2Nii9Njj3
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((3.000000e+00) * jm[0]) + ((8.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 183: // da2Nii9Njj4
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((4.300000e+01) * jm[0]) + ((-4.800000e+01) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 184: // da2Nii9Njj5
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det * (((4.000000e+00) * jm[0]) + ((-5.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 185: // da2Nii9Njj6
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det * (((1.100000e+01) * jm[0]) + ((-8.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 186: // da2Nii9Njj7
            return ((-3.700000e+01) / (2.520000e+03)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 187: // da2Nii9Njj8
            return ((1.000000e+00) / (8.400000e+02)) * jacobi_det * (((3.000000e+00) * jm[0]) + ((-2.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 188: // da2Nii9Njj9
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 189: // da2Nii9Njj10
            return ((3.000000e+00) / (5.600000e+01)) * jacobi_det * jm[0] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 190: // da2Nii10Njj1
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 191: // da2Nii10Njj2
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((-4.000000e+00) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 192: // da2Nii10Njj3
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (jm[0] + ((-4.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 193: // da2Nii10Njj4
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 194: // da2Nii10Njj5
            return (3.000000e+00) * jacobi_det * jm[1] * (1.0 / (((5.600000e+01) * jm[1] * jm[2]) + ((-5.600000e+01) * jm[0] * jm[3])));
        case 195: // da2Nii10Njj6
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((3.000000e+00) * jm[0]) + ((-4.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 196: // da2Nii10Njj7
            return ((9.000000e+00) / (3.500000e+01)) * jacobi_det * jm[0] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 197: // da2Nii10Njj8
            return ((3.000000e+00) / (2.800000e+02)) * jacobi_det * (((4.000000e+00) * jm[0]) + ((-3.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 198: // da2Nii10Njj9
            return ((3.000000e+00) / (5.600000e+01)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
		
		default:
			return 0.0;	
	}	
}

double Element2DTriangleHermiteCubic::get_element_integral_0th_order(
	short elem_shape_func_1, 
	short elem_shape_func_2
) const {
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1;
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
        case 33: // Nii4Njj4
        case 66: // Nii7Njj7
        case 0: // Nii1Njj1
            return ((3.130000e+02) / (5.040000e+03)) * jacobi_det;
        case 2: // Nii1Njj3
        case 10: // Nii2Njj1
        case 20: // Nii3Njj1
        case 35: // Nii4Njj6
        case 53: // Nii6Njj4
        case 67: // Nii7Njj8
        case 76: // Nii8Njj7
        case 1: // Nii1Njj2
            return ((5.300000e+01) / (1.008000e+04)) * jacobi_det;
        case 6: // Nii1Njj7
        case 30: // Nii4Njj1
        case 36: // Nii4Njj7
        case 60: // Nii7Njj1
        case 63: // Nii7Njj4
        case 3: // Nii1Njj4
            return ((1.000000e+00) / (7.200000e+02)) * jacobi_det;
        case 8: // Nii1Njj9
        case 14: // Nii2Njj5
        case 28: // Nii3Njj9
        case 38: // Nii4Njj9
        case 40: // Nii5Njj1
        case 41: // Nii5Njj2
        case 46: // Nii5Njj7
        case 47: // Nii5Njj8
        case 58: // Nii6Njj9
        case 64: // Nii7Njj5
        case 74: // Nii8Njj5
        case 80: // Nii9Njj1
        case 82: // Nii9Njj3
        case 83: // Nii9Njj4
        case 85: // Nii9Njj6
        case 4: // Nii1Njj5
            return ((-1.000000e+00) / (2.520000e+03)) * jacobi_det;
        case 7: // Nii1Njj8
        case 16: // Nii2Njj7
        case 23: // Nii3Njj4
        case 32: // Nii4Njj3
        case 50: // Nii6Njj1
        case 61: // Nii7Njj2
        case 70: // Nii8Njj1
        case 5: // Nii1Njj6
            return ((-1.300000e+01) / (1.008000e+04)) * jacobi_det;
        case 39: // Nii4Njj10
        case 69: // Nii7Njj10
        case 90: // Nii10Njj1
        case 93: // Nii10Njj4
        case 96: // Nii10Njj7
        case 9: // Nii1Njj10
            return ((3.000000e+00) / (1.120000e+02)) * jacobi_det;
        case 22: // Nii3Njj3
        case 55: // Nii6Njj6
        case 77: // Nii8Njj8
        case 11: // Nii2Njj2
            return ((1.000000e+00) / (1.260000e+03)) * jacobi_det;
        case 21: // Nii3Njj2
        case 12: // Nii2Njj3
            return ((1.000000e+00) / (5.040000e+03)) * jacobi_det;
        case 26: // Nii3Njj7
        case 31: // Nii4Njj2
        case 37: // Nii4Njj8
        case 56: // Nii6Njj7
        case 62: // Nii7Njj3
        case 65: // Nii7Njj6
        case 73: // Nii8Njj4
        case 13: // Nii2Njj4
            return ((1.700000e+01) / (1.008000e+04)) * jacobi_det;
        case 27: // Nii3Njj8
        case 51: // Nii6Njj2
        case 72: // Nii8Njj3
        case 15: // Nii2Njj6
            return ((-1.000000e+00) / (1.008000e+04)) * jacobi_det;
        case 25: // Nii3Njj6
        case 52: // Nii6Njj3
        case 71: // Nii8Njj2
        case 17: // Nii2Njj8
            return ((-1.000000e+00) / (5.040000e+03)) * jacobi_det;
        case 24: // Nii3Njj5
        case 42: // Nii5Njj3
        case 81: // Nii9Njj2
        case 18: // Nii2Njj9
            return ((1.000000e+00) / (3.360000e+03)) * jacobi_det;
        case 29: // Nii3Njj10
        case 59: // Nii6Njj10
        case 79: // Nii8Njj10
        case 91: // Nii10Njj2
        case 92: // Nii10Njj3
        case 95: // Nii10Njj6
        case 97: // Nii10Njj8
        case 19: // Nii2Njj10
            return ((3.000000e+00) / (1.120000e+03)) * jacobi_det;
        case 43: // Nii5Njj4
        case 68: // Nii7Njj9
        case 86: // Nii9Njj7
        case 34: // Nii4Njj5
            return ((-5.300000e+01) / (5.040000e+03)) * jacobi_det;
        case 88: // Nii9Njj9
        case 44: // Nii5Njj5
            return ((1.000000e+00) / (5.040000e+02)) * jacobi_det;
        case 54: // Nii6Njj5
        case 78: // Nii8Njj9
        case 87: // Nii9Njj8
        case 45: // Nii5Njj6
            return ((-1.000000e+00) / (1.008000e+03)) * jacobi_det;
        case 84: // Nii9Njj5
        case 48: // Nii5Njj9
            return ((1.000000e+00) / (1.008000e+04)) * jacobi_det;
        case 89: // Nii9Njj10
        case 94: // Nii10Njj5
        case 98: // Nii10Njj9
        case 49: // Nii5Njj10
            return ((-3.000000e+00) / (5.600000e+02)) * jacobi_det;
        case 75: // Nii8Njj6
        case 57: // Nii6Njj8
            return ((1.000000e+00) / (2.016000e+03)) * jacobi_det;
        case 99: // Nii10Njj10
            return ((8.100000e+01) / (5.600000e+02)) * jacobi_det;
		
		default:
			TDKP_GENERAL_EXCEPTION("invalid object index");	
	}	
}
										
double Element2DTriangleHermiteCubic::get_single_integral_1st_order(
	short diffop, 
	short elem_shape_func_1
) const {
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_1 
	                       + 10 * diffop;
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
        case 0: // da1Nii1
            return jacobi_det * (((-1.000000e+00) * jm[2]) + jm[3]) * (1.0 / (((2.000000e+00) * jm[1] * jm[2]) + ((-2.000000e+00) * jm[0] * jm[3])));
        case 8: // da1Nii9
        case 1: // da1Nii2
            return ((1.000000e+00) / (1.200000e+01)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 4: // da1Nii5
        case 2: // da1Nii3
            return jacobi_det * jm[3] * (1.0 / (((1.200000e+01) * jm[1] * jm[2]) + ((-1.200000e+01) * jm[0] * jm[3])));
        case 3: // da1Nii4
            return jacobi_det * jm[3] * (1.0 / (((-2.000000e+00) * jm[1] * jm[2]) + ((2.000000e+00) * jm[0] * jm[3])));
        case 7: // da1Nii8
        case 5: // da1Nii6
            return ((1.000000e+00) / (1.200000e+01)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 6: // da1Nii7
            return jacobi_det * jm[2] * (1.0 / (((2.000000e+00) * jm[1] * jm[2]) + ((-2.000000e+00) * jm[0] * jm[3])));
        case 10: // da2Nii1
            return jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / (((2.000000e+00) * jm[1] * jm[2]) + ((-2.000000e+00) * jm[0] * jm[3])));
        case 18: // da2Nii9
        case 11: // da2Nii2
            return jacobi_det * jm[0] * (1.0 / (((1.200000e+01) * jm[1] * jm[2]) + ((-1.200000e+01) * jm[0] * jm[3])));
        case 14: // da2Nii5
        case 12: // da2Nii3
            return ((1.000000e+00) / (1.200000e+01)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 13: // da2Nii4
            return jacobi_det * jm[1] * (1.0 / (((2.000000e+00) * jm[1] * jm[2]) + ((-2.000000e+00) * jm[0] * jm[3])));
        case 15: // da2Nii6
            return ((1.000000e+00) / (2.000000e+00)) * jacobi_det * (((-1.000000e+00) * jm[0]) + jm[1]) * (1.0 / (((6.000000e+00) * jm[1] * jm[2]) + ((-6.000000e+00) * jm[0] * jm[3])));
        case 16: // da2Nii7
            return jacobi_det * jm[0] * (1.0 / (((-2.000000e+00) * jm[1] * jm[2]) + ((2.000000e+00) * jm[0] * jm[3])));
        case 17: // da2Nii8
            return ((1.000000e+00) / (2.000000e+00)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / (((-6.000000e+00) * jm[1] * jm[2]) + ((6.000000e+00) * jm[0] * jm[3])));		
		default:
			return 0.0;	
	}	
}

double Element2DTriangleHermiteCubic::get_single_integral_0th_order(
	short elem_shape_func_1
) const {
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_1;
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
        case 3: // Nii4
        case 6: // Nii7
        case 0: // Nii1
            return ((1.100000e+01) / (1.200000e+02)) * jacobi_det;
        case 2: // Nii3
        case 5: // Nii6
        case 7: // Nii8
        case 1: // Nii2
            return ((1.000000e+00) / (1.200000e+02)) * jacobi_det;
        case 8: // Nii9
        case 4: // Nii5
            return ((-1.000000e+00) / (6.000000e+01)) * jacobi_det;
        case 9: // Nii10
            return ((9.000000e+00) / (4.000000e+01)) * jacobi_det;		
		default:
			TDKP_GENERAL_EXCEPTION("invalid object index");	
	}	
}
	
double Element2DTriangleHermiteCubic::get_element_integral_0th_order_nodal_data(
	short nodal_data_point, 
	short elem_shape_func_1, 
	short elem_shape_func_2
) const {
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1
	                       + 10 * 10 * nodal_data_point;
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
        case 333: // da4Nii4Njj4
        case 666: // da7Nii7Njj7
        case 0: // da1Nii1Njj1
            return ((2.360900e+04) / (5.544000e+05)) * jacobi_det;
        case 2: // da1Nii1Njj3
        case 10: // da1Nii2Njj1
        case 20: // da1Nii3Njj1
        case 100: // da2Nii1Njj1
        case 200: // da3Nii1Njj1
        case 335: // da4Nii4Njj6
        case 353: // da4Nii6Njj4
        case 533: // da6Nii4Njj4
        case 667: // da7Nii7Njj8
        case 676: // da7Nii8Njj7
        case 766: // da8Nii7Njj7
        case 1: // da1Nii1Njj2
            return ((1.073000e+03) / (3.326400e+05)) * jacobi_det;
        case 6: // da1Nii1Njj7
        case 30: // da1Nii4Njj1
        case 33: // da1Nii4Njj4
        case 60: // da1Nii7Njj1
        case 66: // da1Nii7Njj7
        case 300: // da4Nii1Njj1
        case 303: // da4Nii1Njj4
        case 330: // da4Nii4Njj1
        case 336: // da4Nii4Njj7
        case 363: // da4Nii7Njj4
        case 366: // da4Nii7Njj7
        case 600: // da7Nii1Njj1
        case 606: // da7Nii1Njj7
        case 633: // da7Nii4Njj4
        case 636: // da7Nii4Njj7
        case 660: // da7Nii7Njj1
        case 663: // da7Nii7Njj4
        case 3: // da1Nii1Njj4
            return ((4.700000e+01) / (2.640000e+04)) * jacobi_det;
        case 8: // da1Nii1Njj9
        case 40: // da1Nii5Njj1
        case 80: // da1Nii9Njj1
        case 338: // da4Nii4Njj9
        case 383: // da4Nii9Njj4
        case 400: // da5Nii1Njj1
        case 466: // da5Nii7Njj7
        case 646: // da7Nii5Njj7
        case 664: // da7Nii7Njj5
        case 800: // da9Nii1Njj1
        case 833: // da9Nii4Njj4
        case 4: // da1Nii1Njj5
            return ((-1.970000e+02) / (4.158000e+05)) * jacobi_det;
        case 7: // da1Nii1Njj8
        case 50: // da1Nii6Njj1
        case 70: // da1Nii8Njj1
        case 166: // da2Nii7Njj7
        case 233: // da3Nii4Njj4
        case 323: // da4Nii3Njj4
        case 332: // da4Nii4Njj3
        case 500: // da6Nii1Njj1
        case 616: // da7Nii2Njj7
        case 661: // da7Nii7Njj2
        case 700: // da8Nii1Njj1
        case 5: // da1Nii1Njj6
            return ((-1.090000e+02) / (3.326400e+05)) * jacobi_det;
        case 90: // da1Nii10Njj1
        case 339: // da4Nii4Njj10
        case 393: // da4Nii10Njj4
        case 669: // da7Nii7Njj10
        case 696: // da7Nii10Njj7
        case 900: // da10Nii1Njj1
        case 933: // da10Nii4Njj4
        case 966: // da10Nii7Njj7
        case 9: // da1Nii1Njj10
            return ((9.830000e+02) / (6.160000e+04)) * jacobi_det;
        case 22: // da1Nii3Njj3
        case 101: // da2Nii1Njj2
        case 110: // da2Nii2Njj1
        case 202: // da3Nii1Njj3
        case 220: // da3Nii3Njj1
        case 355: // da4Nii6Njj6
        case 535: // da6Nii4Njj6
        case 553: // da6Nii6Njj4
        case 677: // da7Nii8Njj8
        case 767: // da8Nii7Njj8
        case 776: // da8Nii8Njj7
        case 11: // da1Nii2Njj2
            return ((1.010000e+02) / (2.376000e+05)) * jacobi_det;
        case 21: // da1Nii3Njj2
        case 102: // da2Nii1Njj3
        case 120: // da2Nii3Njj1
        case 201: // da3Nii1Njj2
        case 210: // da3Nii2Njj1
        case 12: // da1Nii2Njj3
            return ((2.300000e+01) / (1.848000e+05)) * jacobi_det;
        case 26: // da1Nii3Njj7
        case 31: // da1Nii4Njj2
        case 62: // da1Nii7Njj3
        case 103: // da2Nii1Njj4
        case 130: // da2Nii4Njj1
        case 206: // da3Nii1Njj7
        case 260: // da3Nii7Njj1
        case 301: // da4Nii1Njj2
        case 310: // da4Nii2Njj1
        case 356: // da4Nii6Njj7
        case 365: // da4Nii7Njj6
        case 367: // da4Nii7Njj8
        case 376: // da4Nii8Njj7
        case 536: // da6Nii4Njj7
        case 563: // da6Nii7Njj4
        case 602: // da7Nii1Njj3
        case 620: // da7Nii3Njj1
        case 635: // da7Nii4Njj6
        case 637: // da7Nii4Njj8
        case 653: // da7Nii6Njj4
        case 673: // da7Nii8Njj4
        case 736: // da8Nii4Njj7
        case 763: // da8Nii7Njj4
        case 13: // da1Nii2Njj4
            return ((5.330000e+02) / (8.316000e+05)) * jacobi_det;
        case 28: // da1Nii3Njj9
        case 41: // da1Nii5Njj2
        case 82: // da1Nii9Njj3
        case 104: // da2Nii1Njj5
        case 140: // da2Nii5Njj1
        case 208: // da3Nii1Njj9
        case 280: // da3Nii9Njj1
        case 358: // da4Nii6Njj9
        case 385: // da4Nii9Njj6
        case 401: // da5Nii1Njj2
        case 410: // da5Nii2Njj1
        case 467: // da5Nii7Njj8
        case 476: // da5Nii8Njj7
        case 538: // da6Nii4Njj9
        case 583: // da6Nii9Njj4
        case 647: // da7Nii5Njj8
        case 674: // da7Nii8Njj5
        case 746: // da8Nii5Njj7
        case 764: // da8Nii7Njj5
        case 802: // da9Nii1Njj3
        case 820: // da9Nii3Njj1
        case 835: // da9Nii4Njj6
        case 853: // da9Nii6Njj4
        case 14: // da1Nii2Njj5
            return ((-2.770000e+02) / (1.663200e+06)) * jacobi_det;
        case 27: // da1Nii3Njj8
        case 51: // da1Nii6Njj2
        case 72: // da1Nii8Njj3
        case 105: // da2Nii1Njj6
        case 150: // da2Nii6Njj1
        case 207: // da3Nii1Njj8
        case 270: // da3Nii8Njj1
        case 501: // da6Nii1Njj2
        case 510: // da6Nii2Njj1
        case 702: // da8Nii1Njj3
        case 720: // da8Nii3Njj1
        case 15: // da1Nii2Njj6
            return ((-1.700000e+01) / (8.316000e+05)) * jacobi_det;
        case 23: // da1Nii3Njj4
        case 32: // da1Nii4Njj3
        case 35: // da1Nii4Njj6
        case 53: // da1Nii6Njj4
        case 61: // da1Nii7Njj2
        case 67: // da1Nii7Njj8
        case 76: // da1Nii8Njj7
        case 106: // da2Nii1Njj7
        case 160: // da2Nii7Njj1
        case 203: // da3Nii1Njj4
        case 230: // da3Nii4Njj1
        case 302: // da4Nii1Njj3
        case 305: // da4Nii1Njj6
        case 320: // da4Nii3Njj1
        case 350: // da4Nii6Njj1
        case 503: // da6Nii1Njj4
        case 530: // da6Nii4Njj1
        case 601: // da7Nii1Njj2
        case 607: // da7Nii1Njj8
        case 610: // da7Nii2Njj1
        case 670: // da7Nii8Njj1
        case 706: // da8Nii1Njj7
        case 760: // da8Nii7Njj1
        case 16: // da1Nii2Njj7
            return ((-6.100000e+01) / (3.326400e+05)) * jacobi_det;
        case 25: // da1Nii3Njj6
        case 52: // da1Nii6Njj3
        case 71: // da1Nii8Njj2
        case 107: // da2Nii1Njj8
        case 167: // da2Nii7Njj8
        case 170: // da2Nii8Njj1
        case 176: // da2Nii8Njj7
        case 205: // da3Nii1Njj6
        case 235: // da3Nii4Njj6
        case 250: // da3Nii6Njj1
        case 253: // da3Nii6Njj4
        case 325: // da4Nii3Njj6
        case 352: // da4Nii6Njj3
        case 502: // da6Nii1Njj3
        case 520: // da6Nii3Njj1
        case 523: // da6Nii3Njj4
        case 532: // da6Nii4Njj3
        case 617: // da7Nii2Njj8
        case 671: // da7Nii8Njj2
        case 701: // da8Nii1Njj2
        case 710: // da8Nii2Njj1
        case 716: // da8Nii2Njj7
        case 761: // da8Nii7Njj2
        case 17: // da1Nii2Njj8
            return ((-5.900000e+01) / (1.663200e+06)) * jacobi_det;
        case 24: // da1Nii3Njj5
        case 42: // da1Nii5Njj3
        case 81: // da1Nii9Njj2
        case 108: // da2Nii1Njj9
        case 180: // da2Nii9Njj1
        case 204: // da3Nii1Njj5
        case 240: // da3Nii5Njj1
        case 402: // da5Nii1Njj3
        case 420: // da5Nii3Njj1
        case 801: // da9Nii1Njj2
        case 810: // da9Nii2Njj1
        case 18: // da1Nii2Njj9
            return ((1.000000e+00) / (1.980000e+04)) * jacobi_det;
        case 29: // da1Nii3Njj10
        case 91: // da1Nii10Njj2
        case 92: // da1Nii10Njj3
        case 109: // da2Nii1Njj10
        case 190: // da2Nii10Njj1
        case 209: // da3Nii1Njj10
        case 290: // da3Nii10Njj1
        case 359: // da4Nii6Njj10
        case 395: // da4Nii10Njj6
        case 539: // da6Nii4Njj10
        case 593: // da6Nii10Njj4
        case 679: // da7Nii8Njj10
        case 697: // da7Nii10Njj8
        case 769: // da8Nii7Njj10
        case 796: // da8Nii10Njj7
        case 901: // da10Nii1Njj2
        case 902: // da10Nii1Njj3
        case 910: // da10Nii2Njj1
        case 920: // da10Nii3Njj1
        case 935: // da10Nii4Njj6
        case 953: // da10Nii6Njj4
        case 967: // da10Nii7Njj8
        case 976: // da10Nii8Njj7
        case 19: // da1Nii2Njj10
            return ((9.700000e+01) / (6.160000e+04)) * jacobi_det;
        case 43: // da1Nii5Njj4
        case 68: // da1Nii7Njj9
        case 86: // da1Nii9Njj7
        case 304: // da4Nii1Njj5
        case 340: // da4Nii5Njj1
        case 346: // da4Nii5Njj7
        case 364: // da4Nii7Njj5
        case 368: // da4Nii7Njj9
        case 386: // da4Nii9Njj7
        case 403: // da5Nii1Njj4
        case 430: // da5Nii4Njj1
        case 436: // da5Nii4Njj7
        case 463: // da5Nii7Njj4
        case 608: // da7Nii1Njj9
        case 634: // da7Nii4Njj5
        case 638: // da7Nii4Njj9
        case 643: // da7Nii5Njj4
        case 680: // da7Nii9Njj1
        case 683: // da7Nii9Njj4
        case 806: // da9Nii1Njj7
        case 836: // da9Nii4Njj7
        case 860: // da9Nii7Njj1
        case 863: // da9Nii7Njj4
        case 34: // da1Nii4Njj5
            return ((-7.610000e+02) / (1.663200e+06)) * jacobi_det;
        case 63: // da1Nii7Njj4
        case 306: // da4Nii1Njj7
        case 360: // da4Nii7Njj1
        case 603: // da7Nii1Njj4
        case 630: // da7Nii4Njj1
        case 36: // da1Nii4Njj7
            return ((-3.100000e+01) / (2.520000e+04)) * jacobi_det;
        case 56: // da1Nii6Njj7
        case 65: // da1Nii7Njj6
        case 73: // da1Nii8Njj4
        case 136: // da2Nii4Njj7
        case 163: // da2Nii7Njj4
        case 236: // da3Nii4Njj7
        case 263: // da3Nii7Njj4
        case 307: // da4Nii1Njj8
        case 316: // da4Nii2Njj7
        case 326: // da4Nii3Njj7
        case 361: // da4Nii7Njj2
        case 362: // da4Nii7Njj3
        case 370: // da4Nii8Njj1
        case 506: // da6Nii1Njj7
        case 560: // da6Nii7Njj1
        case 605: // da7Nii1Njj6
        case 613: // da7Nii2Njj4
        case 623: // da7Nii3Njj4
        case 631: // da7Nii4Njj2
        case 632: // da7Nii4Njj3
        case 650: // da7Nii6Njj1
        case 703: // da8Nii1Njj4
        case 730: // da8Nii4Njj1
        case 37: // da1Nii4Njj8
            return ((-2.690000e+02) / (1.663200e+06)) * jacobi_det;
        case 46: // da1Nii5Njj7
        case 64: // da1Nii7Njj5
        case 83: // da1Nii9Njj4
        case 308: // da4Nii1Njj9
        case 380: // da4Nii9Njj1
        case 406: // da5Nii1Njj7
        case 460: // da5Nii7Njj1
        case 604: // da7Nii1Njj5
        case 640: // da7Nii5Njj1
        case 803: // da9Nii1Njj4
        case 830: // da9Nii4Njj1
        case 38: // da1Nii4Njj9
            return ((2.690000e+02) / (8.316000e+05)) * jacobi_det;
        case 69: // da1Nii7Njj10
        case 93: // da1Nii10Njj4
        case 96: // da1Nii10Njj7
        case 309: // da4Nii1Njj10
        case 369: // da4Nii7Njj10
        case 390: // da4Nii10Njj1
        case 396: // da4Nii10Njj7
        case 609: // da7Nii1Njj10
        case 639: // da7Nii4Njj10
        case 690: // da7Nii10Njj1
        case 693: // da7Nii10Njj4
        case 903: // da10Nii1Njj4
        case 906: // da10Nii1Njj7
        case 930: // da10Nii4Njj1
        case 936: // da10Nii4Njj7
        case 960: // da10Nii7Njj1
        case 963: // da10Nii7Njj4
        case 39: // da1Nii4Njj10
            return ((-2.900000e+01) / (3.080000e+04)) * jacobi_det;
        case 88: // da1Nii9Njj9
        case 388: // da4Nii9Njj9
        case 404: // da5Nii1Njj5
        case 440: // da5Nii5Njj1
        case 446: // da5Nii5Njj7
        case 464: // da5Nii7Njj5
        case 644: // da7Nii5Njj5
        case 808: // da9Nii1Njj9
        case 838: // da9Nii4Njj9
        case 880: // da9Nii9Njj1
        case 883: // da9Nii9Njj4
        case 44: // da1Nii5Njj5
            return ((7.000000e+00) / (5.940000e+04)) * jacobi_det;
        case 54: // da1Nii6Njj5
        case 78: // da1Nii8Njj9
        case 87: // da1Nii9Njj8
        case 405: // da5Nii1Njj6
        case 450: // da5Nii6Njj1
        case 504: // da6Nii1Njj5
        case 540: // da6Nii5Njj1
        case 708: // da8Nii1Njj9
        case 780: // da8Nii9Njj1
        case 807: // da9Nii1Njj8
        case 870: // da9Nii8Njj1
        case 45: // da1Nii5Njj6
            return ((6.100000e+01) / (1.663200e+06)) * jacobi_det;
        case 58: // da1Nii6Njj9
        case 74: // da1Nii8Njj5
        case 85: // da1Nii9Njj6
        case 146: // da2Nii5Njj7
        case 164: // da2Nii7Njj5
        case 238: // da3Nii4Njj9
        case 283: // da3Nii9Njj4
        case 328: // da4Nii3Njj9
        case 382: // da4Nii9Njj3
        case 407: // da5Nii1Njj8
        case 416: // da5Nii2Njj7
        case 458: // da5Nii6Njj9
        case 461: // da5Nii7Njj2
        case 470: // da5Nii8Njj1
        case 478: // da5Nii8Njj9
        case 485: // da5Nii9Njj6
        case 487: // da5Nii9Njj8
        case 508: // da6Nii1Njj9
        case 548: // da6Nii5Njj9
        case 580: // da6Nii9Njj1
        case 584: // da6Nii9Njj5
        case 614: // da7Nii2Njj5
        case 641: // da7Nii5Njj2
        case 704: // da8Nii1Njj5
        case 740: // da8Nii5Njj1
        case 748: // da8Nii5Njj9
        case 784: // da8Nii9Njj5
        case 805: // da9Nii1Njj6
        case 823: // da9Nii3Njj4
        case 832: // da9Nii4Njj3
        case 845: // da9Nii5Njj6
        case 847: // da9Nii5Njj8
        case 850: // da9Nii6Njj1
        case 854: // da9Nii6Njj5
        case 874: // da9Nii8Njj5
        case 47: // da1Nii5Njj8
            return ((1.300000e+01) / (3.326400e+05)) * jacobi_det;
        case 84: // da1Nii9Njj5
        case 408: // da5Nii1Njj9
        case 480: // da5Nii9Njj1
        case 804: // da9Nii1Njj5
        case 840: // da9Nii5Njj1
        case 48: // da1Nii5Njj9
            return ((-1.000000e+00) / (1.188000e+04)) * jacobi_det;
        case 89: // da1Nii9Njj10
        case 94: // da1Nii10Njj5
        case 98: // da1Nii10Njj9
        case 389: // da4Nii9Njj10
        case 398: // da4Nii10Njj9
        case 409: // da5Nii1Njj10
        case 469: // da5Nii7Njj10
        case 490: // da5Nii10Njj1
        case 496: // da5Nii10Njj7
        case 649: // da7Nii5Njj10
        case 694: // da7Nii10Njj5
        case 809: // da9Nii1Njj10
        case 839: // da9Nii4Njj10
        case 890: // da9Nii10Njj1
        case 893: // da9Nii10Njj4
        case 904: // da10Nii1Njj5
        case 908: // da10Nii1Njj9
        case 938: // da10Nii4Njj9
        case 940: // da10Nii5Njj1
        case 946: // da10Nii5Njj7
        case 964: // da10Nii7Njj5
        case 980: // da10Nii9Njj1
        case 983: // da10Nii9Njj4
        case 49: // da1Nii5Njj10
            return ((1.300000e+01) / (6.160000e+04)) * jacobi_det;
        case 77: // da1Nii8Njj8
        case 116: // da2Nii2Njj7
        case 161: // da2Nii7Njj2
        case 223: // da3Nii3Njj4
        case 232: // da3Nii4Njj3
        case 322: // da4Nii3Njj3
        case 505: // da6Nii1Njj6
        case 550: // da6Nii6Njj1
        case 611: // da7Nii2Njj2
        case 707: // da8Nii1Njj8
        case 770: // da8Nii8Njj1
        case 55: // da1Nii6Njj6
            return ((-1.000000e+00) / (3.024000e+04)) * jacobi_det;
        case 75: // da1Nii8Njj6
        case 507: // da6Nii1Njj8
        case 570: // da6Nii8Njj1
        case 705: // da8Nii1Njj6
        case 750: // da8Nii6Njj1
        case 57: // da1Nii6Njj8
            return ((-1.000000e+00) / (4.158000e+04)) * jacobi_det;
        case 79: // da1Nii8Njj10
        case 95: // da1Nii10Njj6
        case 97: // da1Nii10Njj8
        case 169: // da2Nii7Njj10
        case 196: // da2Nii10Njj7
        case 239: // da3Nii4Njj10
        case 293: // da3Nii10Njj4
        case 329: // da4Nii3Njj10
        case 392: // da4Nii10Njj3
        case 509: // da6Nii1Njj10
        case 590: // da6Nii10Njj1
        case 619: // da7Nii2Njj10
        case 691: // da7Nii10Njj2
        case 709: // da8Nii1Njj10
        case 790: // da8Nii10Njj1
        case 905: // da10Nii1Njj6
        case 907: // da10Nii1Njj8
        case 916: // da10Nii2Njj7
        case 923: // da10Nii3Njj4
        case 932: // da10Nii4Njj3
        case 950: // da10Nii6Njj1
        case 961: // da10Nii7Njj2
        case 970: // da10Nii8Njj1
        case 59: // da1Nii6Njj10
            return ((-1.900000e+01) / (3.080000e+04)) * jacobi_det;
        case 399: // da4Nii10Njj10
        case 699: // da7Nii10Njj10
        case 909: // da10Nii1Njj10
        case 939: // da10Nii4Njj10
        case 969: // da10Nii7Njj10
        case 990: // da10Nii10Njj1
        case 993: // da10Nii10Njj4
        case 996: // da10Nii10Njj7
        case 99: // da1Nii10Njj10
            return ((7.830000e+02) / (6.160000e+04)) * jacobi_det;
        case 222: // da3Nii3Njj3
        case 555: // da6Nii6Njj6
        case 777: // da8Nii8Njj8
        case 111: // da2Nii2Njj2
            return ((1.300000e+01) / (1.848000e+05)) * jacobi_det;
        case 121: // da2Nii3Njj2
        case 122: // da2Nii3Njj3
        case 211: // da3Nii2Njj2
        case 212: // da3Nii2Njj3
        case 221: // da3Nii3Njj2
        case 112: // da2Nii2Njj3
            return ((1.700000e+01) / (1.663200e+06)) * jacobi_det;
        case 131: // da2Nii4Njj2
        case 226: // da3Nii3Njj7
        case 262: // da3Nii7Njj3
        case 311: // da4Nii2Njj2
        case 377: // da4Nii8Njj8
        case 556: // da6Nii6Njj7
        case 565: // da6Nii7Njj6
        case 622: // da7Nii3Njj3
        case 655: // da7Nii6Njj6
        case 737: // da8Nii4Njj8
        case 773: // da8Nii8Njj4
        case 113: // da2Nii2Njj4
            return ((2.630000e+02) / (1.663200e+06)) * jacobi_det;
        case 141: // da2Nii5Njj2
        case 228: // da3Nii3Njj9
        case 282: // da3Nii9Njj3
        case 411: // da5Nii2Njj2
        case 477: // da5Nii8Njj8
        case 558: // da6Nii6Njj9
        case 585: // da6Nii9Njj6
        case 747: // da8Nii5Njj8
        case 774: // da8Nii8Njj5
        case 822: // da9Nii3Njj3
        case 855: // da9Nii6Njj6
        case 114: // da2Nii2Njj5
            return ((-1.000000e+00) / (2.520000e+04)) * jacobi_det;
        case 151: // da2Nii6Njj2
        case 227: // da3Nii3Njj8
        case 272: // da3Nii8Njj3
        case 511: // da6Nii2Njj2
        case 722: // da8Nii3Njj3
        case 115: // da2Nii2Njj6
            return ((-1.000000e+00) / (5.544000e+05)) * jacobi_det;
        case 171: // da2Nii8Njj2
        case 177: // da2Nii8Njj8
        case 225: // da3Nii3Njj6
        case 252: // da3Nii6Njj3
        case 255: // da3Nii6Njj6
        case 522: // da6Nii3Njj3
        case 525: // da6Nii3Njj6
        case 552: // da6Nii6Njj3
        case 711: // da8Nii2Njj2
        case 717: // da8Nii2Njj8
        case 771: // da8Nii8Njj2
        case 117: // da2Nii2Njj8
            return ((-1.000000e+00) / (1.848000e+05)) * jacobi_det;
        case 181: // da2Nii9Njj2
        case 224: // da3Nii3Njj5
        case 242: // da3Nii5Njj3
        case 422: // da5Nii3Njj3
        case 811: // da9Nii2Njj2
        case 118: // da2Nii2Njj9
            return ((1.000000e+00) / (1.039500e+05)) * jacobi_det;
        case 191: // da2Nii10Njj2
        case 229: // da3Nii3Njj10
        case 292: // da3Nii10Njj3
        case 559: // da6Nii6Njj10
        case 595: // da6Nii10Njj6
        case 779: // da8Nii8Njj10
        case 797: // da8Nii10Njj8
        case 911: // da10Nii2Njj2
        case 922: // da10Nii3Njj3
        case 955: // da10Nii6Njj6
        case 977: // da10Nii8Njj8
        case 119: // da2Nii2Njj10
            return ((3.000000e+00) / (1.232000e+04)) * jacobi_det;
        case 126: // da2Nii3Njj7
        case 132: // da2Nii4Njj3
        case 157: // da2Nii6Njj8
        case 162: // da2Nii7Njj3
        case 175: // da2Nii8Njj6
        case 213: // da3Nii2Njj4
        case 216: // da3Nii2Njj7
        case 231: // da3Nii4Njj2
        case 257: // da3Nii6Njj8
        case 261: // da3Nii7Njj2
        case 275: // da3Nii8Njj6
        case 312: // da4Nii2Njj3
        case 321: // da4Nii3Njj2
        case 517: // da6Nii2Njj8
        case 527: // da6Nii3Njj8
        case 571: // da6Nii8Njj2
        case 572: // da6Nii8Njj3
        case 612: // da7Nii2Njj3
        case 621: // da7Nii3Njj2
        case 715: // da8Nii2Njj6
        case 725: // da8Nii3Njj6
        case 751: // da8Nii6Njj2
        case 752: // da8Nii6Njj3
        case 123: // da2Nii3Njj4
            return ((-1.000000e+00) / (2.772000e+05)) * jacobi_det;
        case 128: // da2Nii3Njj9
        case 142: // da2Nii5Njj3
        case 182: // da2Nii9Njj3
        case 214: // da3Nii2Njj5
        case 218: // da3Nii2Njj9
        case 241: // da3Nii5Njj2
        case 281: // da3Nii9Njj2
        case 412: // da5Nii2Njj3
        case 421: // da5Nii3Njj2
        case 812: // da9Nii2Njj3
        case 821: // da9Nii3Njj2
        case 124: // da2Nii3Njj5
            return ((1.000000e+00) / (1.663200e+06)) * jacobi_det;
        case 127: // da2Nii3Njj8
        case 152: // da2Nii6Njj3
        case 172: // da2Nii8Njj3
        case 215: // da3Nii2Njj6
        case 217: // da3Nii2Njj8
        case 251: // da3Nii6Njj2
        case 271: // da3Nii8Njj2
        case 512: // da6Nii2Njj3
        case 521: // da6Nii3Njj2
        case 712: // da8Nii2Njj3
        case 721: // da8Nii3Njj2
        case 125: // da2Nii3Njj6
            return ((-1.000000e+00) / (8.316000e+05)) * jacobi_det;
        case 192: // da2Nii10Njj3
        case 219: // da3Nii2Njj10
        case 291: // da3Nii10Njj2
        case 912: // da10Nii2Njj3
        case 921: // da10Nii3Njj2
        case 129: // da2Nii3Njj10
            return ((1.000000e+00) / (1.232000e+04)) * jacobi_det;
        case 266: // da3Nii7Njj7
        case 313: // da4Nii2Njj4
        case 331: // da4Nii4Njj2
        case 337: // da4Nii4Njj8
        case 373: // da4Nii8Njj4
        case 566: // da6Nii7Njj7
        case 626: // da7Nii3Njj7
        case 656: // da7Nii6Njj7
        case 662: // da7Nii7Njj3
        case 665: // da7Nii7Njj6
        case 733: // da8Nii4Njj4
        case 133: // da2Nii4Njj4
            return ((1.333000e+03) / (1.663200e+06)) * jacobi_det;
        case 143: // da2Nii5Njj4
        case 268: // da3Nii7Njj9
        case 286: // da3Nii9Njj7
        case 314: // da4Nii2Njj5
        case 341: // da4Nii5Njj2
        case 347: // da4Nii5Njj8
        case 374: // da4Nii8Njj5
        case 413: // da5Nii2Njj4
        case 431: // da5Nii4Njj2
        case 437: // da5Nii4Njj8
        case 473: // da5Nii8Njj4
        case 568: // da6Nii7Njj9
        case 586: // da6Nii9Njj7
        case 628: // da7Nii3Njj9
        case 658: // da7Nii6Njj9
        case 682: // da7Nii9Njj3
        case 685: // da7Nii9Njj6
        case 734: // da8Nii4Njj5
        case 743: // da8Nii5Njj4
        case 826: // da9Nii3Njj7
        case 856: // da9Nii6Njj7
        case 862: // da9Nii7Njj3
        case 865: // da9Nii7Njj6
        case 134: // da2Nii4Njj5
            return ((-1.300000e+01) / (7.560000e+04)) * jacobi_det;
        case 153: // da2Nii6Njj4
        case 156: // da2Nii6Njj7
        case 165: // da2Nii7Njj6
        case 237: // da3Nii4Njj8
        case 267: // da3Nii7Njj8
        case 273: // da3Nii8Njj4
        case 276: // da3Nii8Njj7
        case 315: // da4Nii2Njj6
        case 327: // da4Nii3Njj8
        case 351: // da4Nii6Njj2
        case 372: // da4Nii8Njj3
        case 513: // da6Nii2Njj4
        case 516: // da6Nii2Njj7
        case 531: // da6Nii4Njj2
        case 561: // da6Nii7Njj2
        case 615: // da7Nii2Njj6
        case 627: // da7Nii3Njj8
        case 651: // da7Nii6Njj2
        case 672: // da7Nii8Njj3
        case 723: // da8Nii3Njj4
        case 726: // da8Nii3Njj7
        case 732: // da8Nii4Njj3
        case 762: // da8Nii7Njj3
        case 135: // da2Nii4Njj6
            return ((-1.000000e+00) / (6.652800e+04)) * jacobi_det;
        case 173: // da2Nii8Njj4
        case 256: // da3Nii6Njj7
        case 265: // da3Nii7Njj6
        case 317: // da4Nii2Njj8
        case 371: // da4Nii8Njj2
        case 526: // da6Nii3Njj7
        case 562: // da6Nii7Njj3
        case 625: // da7Nii3Njj6
        case 652: // da7Nii6Njj3
        case 713: // da8Nii2Njj4
        case 731: // da8Nii4Njj2
        case 137: // da2Nii4Njj8
            return ((-1.000000e+00) / (3.326400e+04)) * jacobi_det;
        case 183: // da2Nii9Njj4
        case 246: // da3Nii5Njj7
        case 264: // da3Nii7Njj5
        case 318: // da4Nii2Njj9
        case 381: // da4Nii9Njj2
        case 426: // da5Nii3Njj7
        case 462: // da5Nii7Njj3
        case 624: // da7Nii3Njj5
        case 642: // da7Nii5Njj3
        case 813: // da9Nii2Njj4
        case 831: // da9Nii4Njj2
        case 138: // da2Nii4Njj9
            return ((1.000000e+00) / (2.217600e+04)) * jacobi_det;
        case 193: // da2Nii10Njj4
        case 269: // da3Nii7Njj10
        case 296: // da3Nii10Njj7
        case 319: // da4Nii2Njj10
        case 379: // da4Nii8Njj10
        case 391: // da4Nii10Njj2
        case 397: // da4Nii10Njj8
        case 569: // da6Nii7Njj10
        case 596: // da6Nii10Njj7
        case 629: // da7Nii3Njj10
        case 659: // da7Nii6Njj10
        case 692: // da7Nii10Njj3
        case 695: // da7Nii10Njj6
        case 739: // da8Nii4Njj10
        case 793: // da8Nii10Njj4
        case 913: // da10Nii2Njj4
        case 926: // da10Nii3Njj7
        case 931: // da10Nii4Njj2
        case 937: // da10Nii4Njj8
        case 956: // da10Nii6Njj7
        case 962: // da10Nii7Njj3
        case 965: // da10Nii7Njj6
        case 973: // da10Nii8Njj4
        case 139: // da2Nii4Njj10
            return ((1.000000e+00) / (2.464000e+03)) * jacobi_det;
        case 288: // da3Nii9Njj9
        case 414: // da5Nii2Njj5
        case 441: // da5Nii5Njj2
        case 447: // da5Nii5Njj8
        case 474: // da5Nii8Njj5
        case 588: // da6Nii9Njj9
        case 744: // da8Nii5Njj5
        case 828: // da9Nii3Njj9
        case 858: // da9Nii6Njj9
        case 882: // da9Nii9Njj3
        case 885: // da9Nii9Njj6
        case 144: // da2Nii5Njj5
            return ((2.000000e+00) / (5.197500e+04)) * jacobi_det;
        case 154: // da2Nii6Njj5
        case 158: // da2Nii6Njj9
        case 185: // da2Nii9Njj6
        case 247: // da3Nii5Njj8
        case 274: // da3Nii8Njj5
        case 278: // da3Nii8Njj9
        case 287: // da3Nii9Njj8
        case 415: // da5Nii2Njj6
        case 427: // da5Nii3Njj8
        case 451: // da5Nii6Njj2
        case 472: // da5Nii8Njj3
        case 514: // da6Nii2Njj5
        case 518: // da6Nii2Njj9
        case 541: // da6Nii5Njj2
        case 581: // da6Nii9Njj2
        case 724: // da8Nii3Njj5
        case 728: // da8Nii3Njj9
        case 742: // da8Nii5Njj3
        case 782: // da8Nii9Njj3
        case 815: // da9Nii2Njj6
        case 827: // da9Nii3Njj8
        case 851: // da9Nii6Njj2
        case 872: // da9Nii8Njj3
        case 145: // da2Nii5Njj6
            return ((1.000000e+00) / (2.772000e+05)) * jacobi_det;
        case 174: // da2Nii8Njj5
        case 258: // da3Nii6Njj9
        case 285: // da3Nii9Njj6
        case 417: // da5Nii2Njj8
        case 471: // da5Nii8Njj2
        case 528: // da6Nii3Njj9
        case 582: // da6Nii9Njj3
        case 714: // da8Nii2Njj5
        case 741: // da8Nii5Njj2
        case 825: // da9Nii3Njj6
        case 852: // da9Nii6Njj3
        case 147: // da2Nii5Njj8
            return ((1.000000e+00) / (1.386000e+05)) * jacobi_det;
        case 184: // da2Nii9Njj5
        case 248: // da3Nii5Njj9
        case 284: // da3Nii9Njj5
        case 418: // da5Nii2Njj9
        case 428: // da5Nii3Njj9
        case 481: // da5Nii9Njj2
        case 482: // da5Nii9Njj3
        case 814: // da9Nii2Njj5
        case 824: // da9Nii3Njj5
        case 841: // da9Nii5Njj2
        case 842: // da9Nii5Njj3
        case 148: // da2Nii5Njj9
            return ((-1.000000e+00) / (9.240000e+04)) * jacobi_det;
        case 179: // da2Nii8Njj10
        case 194: // da2Nii10Njj5
        case 197: // da2Nii10Njj8
        case 259: // da3Nii6Njj10
        case 289: // da3Nii9Njj10
        case 295: // da3Nii10Njj6
        case 298: // da3Nii10Njj9
        case 419: // da5Nii2Njj10
        case 479: // da5Nii8Njj10
        case 491: // da5Nii10Njj2
        case 497: // da5Nii10Njj8
        case 529: // da6Nii3Njj10
        case 589: // da6Nii9Njj10
        case 592: // da6Nii10Njj3
        case 598: // da6Nii10Njj9
        case 719: // da8Nii2Njj10
        case 749: // da8Nii5Njj10
        case 791: // da8Nii10Njj2
        case 794: // da8Nii10Njj5
        case 829: // da9Nii3Njj10
        case 859: // da9Nii6Njj10
        case 892: // da9Nii10Njj3
        case 895: // da9Nii10Njj6
        case 914: // da10Nii2Njj5
        case 917: // da10Nii2Njj8
        case 925: // da10Nii3Njj6
        case 928: // da10Nii3Njj9
        case 941: // da10Nii5Njj2
        case 947: // da10Nii5Njj8
        case 952: // da10Nii6Njj3
        case 958: // da10Nii6Njj9
        case 971: // da10Nii8Njj2
        case 974: // da10Nii8Njj5
        case 982: // da10Nii9Njj3
        case 985: // da10Nii9Njj6
        case 149: // da2Nii5Njj10
            return ((-3.000000e+00) / (3.080000e+04)) * jacobi_det;
        case 277: // da3Nii8Njj8
        case 515: // da6Nii2Njj6
        case 551: // da6Nii6Njj2
        case 727: // da8Nii3Njj8
        case 772: // da8Nii8Njj3
        case 155: // da2Nii6Njj6
            return ((-1.000000e+00) / (2.376000e+05)) * jacobi_det;
        case 195: // da2Nii10Njj6
        case 279: // da3Nii8Njj10
        case 297: // da3Nii10Njj8
        case 489: // da5Nii9Njj10
        case 498: // da5Nii10Njj9
        case 519: // da6Nii2Njj10
        case 591: // da6Nii10Njj2
        case 729: // da8Nii3Njj10
        case 792: // da8Nii10Njj3
        case 849: // da9Nii5Njj10
        case 894: // da9Nii10Njj5
        case 915: // da10Nii2Njj6
        case 927: // da10Nii3Njj8
        case 948: // da10Nii5Njj9
        case 951: // da10Nii6Njj2
        case 972: // da10Nii8Njj3
        case 984: // da10Nii9Njj5
        case 159: // da2Nii6Njj10
            return ((-3.000000e+00) / (6.160000e+04)) * jacobi_det;
        case 186: // da2Nii9Njj7
        case 234: // da3Nii4Njj5
        case 243: // da3Nii5Njj4
        case 324: // da4Nii3Njj5
        case 342: // da4Nii5Njj3
        case 423: // da5Nii3Njj4
        case 432: // da5Nii4Njj3
        case 618: // da7Nii2Njj9
        case 681: // da7Nii9Njj2
        case 816: // da9Nii2Njj7
        case 861: // da9Nii7Njj2
        case 168: // da2Nii7Njj9
            return ((3.100000e+01) / (5.544000e+05)) * jacobi_det;
        case 187: // da2Nii9Njj8
        case 245: // da3Nii5Njj6
        case 254: // da3Nii6Njj5
        case 425: // da5Nii3Njj6
        case 452: // da5Nii6Njj3
        case 524: // da6Nii3Njj5
        case 542: // da6Nii5Njj3
        case 718: // da8Nii2Njj9
        case 781: // da8Nii9Njj2
        case 817: // da9Nii2Njj8
        case 871: // da9Nii8Njj2
        case 178: // da2Nii8Njj9
            return ((1.000000e+00) / (1.512000e+05)) * jacobi_det;
        case 244: // da3Nii5Njj5
        case 424: // da5Nii3Njj5
        case 442: // da5Nii5Njj3
        case 818: // da9Nii2Njj9
        case 881: // da9Nii9Njj2
        case 188: // da2Nii9Njj9
            return ((-1.000000e+00) / (1.039500e+05)) * jacobi_det;
        case 198: // da2Nii10Njj9
        case 249: // da3Nii5Njj10
        case 294: // da3Nii10Njj5
        case 429: // da5Nii3Njj10
        case 492: // da5Nii10Njj3
        case 579: // da6Nii8Njj10
        case 597: // da6Nii10Njj8
        case 759: // da8Nii6Njj10
        case 795: // da8Nii10Njj6
        case 819: // da9Nii2Njj10
        case 891: // da9Nii10Njj2
        case 918: // da10Nii2Njj9
        case 924: // da10Nii3Njj5
        case 942: // da10Nii5Njj3
        case 957: // da10Nii6Njj8
        case 975: // da10Nii8Njj6
        case 981: // da10Nii9Njj2
        case 189: // da2Nii9Njj10
            return ((9.000000e+00) / (6.160000e+04)) * jacobi_det;
        case 299: // da3Nii10Njj10
        case 599: // da6Nii10Njj10
        case 799: // da8Nii10Njj10
        case 919: // da10Nii2Njj10
        case 929: // da10Nii3Njj10
        case 959: // da10Nii6Njj10
        case 979: // da10Nii8Njj10
        case 991: // da10Nii10Njj2
        case 992: // da10Nii10Njj3
        case 995: // da10Nii10Njj6
        case 997: // da10Nii10Njj8
        case 199: // da2Nii10Njj10
            return ((8.100000e+01) / (6.160000e+04)) * jacobi_det;
        case 343: // da4Nii5Njj4
        case 433: // da5Nii4Njj4
        case 668: // da7Nii7Njj9
        case 686: // da7Nii9Njj7
        case 866: // da9Nii7Njj7
        case 334: // da4Nii4Njj5
            return ((-1.073000e+03) / (1.663200e+05)) * jacobi_det;
        case 434: // da5Nii4Njj5
        case 443: // da5Nii5Njj4
        case 688: // da7Nii9Njj9
        case 868: // da9Nii7Njj9
        case 886: // da9Nii9Njj7
        case 344: // da4Nii5Njj5
            return ((4.570000e+02) / (4.158000e+05)) * jacobi_det;
        case 354: // da4Nii6Njj5
        case 435: // da5Nii4Njj6
        case 453: // da5Nii6Njj4
        case 534: // da6Nii4Njj5
        case 543: // da6Nii5Njj4
        case 678: // da7Nii8Njj9
        case 687: // da7Nii9Njj8
        case 768: // da8Nii7Njj9
        case 786: // da8Nii9Njj7
        case 867: // da9Nii7Njj8
        case 876: // da9Nii8Njj7
        case 345: // da4Nii5Njj6
            return ((-4.570000e+02) / (8.316000e+05)) * jacobi_det;
        case 384: // da4Nii9Njj5
        case 438: // da5Nii4Njj9
        case 468: // da5Nii7Njj9
        case 483: // da5Nii9Njj4
        case 486: // da5Nii9Njj7
        case 648: // da7Nii5Njj9
        case 684: // da7Nii9Njj5
        case 834: // da9Nii4Njj5
        case 843: // da9Nii5Njj4
        case 846: // da9Nii5Njj7
        case 864: // da9Nii7Njj5
        case 348: // da4Nii5Njj9
            return ((1.930000e+02) / (1.663200e+06)) * jacobi_det;
        case 394: // da4Nii10Njj5
        case 439: // da5Nii4Njj10
        case 493: // da5Nii10Njj4
        case 689: // da7Nii9Njj10
        case 698: // da7Nii10Njj9
        case 869: // da9Nii7Njj10
        case 896: // da9Nii10Njj7
        case 934: // da10Nii4Njj5
        case 943: // da10Nii5Njj4
        case 968: // da10Nii7Njj9
        case 986: // da10Nii9Njj7
        case 349: // da4Nii5Njj10
            return ((-9.700000e+01) / (3.080000e+04)) * jacobi_det;
        case 375: // da4Nii8Njj6
        case 537: // da6Nii4Njj8
        case 567: // da6Nii7Njj8
        case 573: // da6Nii8Njj4
        case 576: // da6Nii8Njj7
        case 657: // da7Nii6Njj8
        case 675: // da7Nii8Njj6
        case 735: // da8Nii4Njj6
        case 753: // da8Nii6Njj4
        case 756: // da8Nii6Njj7
        case 765: // da8Nii7Njj6
        case 357: // da4Nii6Njj8
            return ((3.110000e+02) / (1.663200e+06)) * jacobi_det;
        case 387: // da4Nii9Njj8
        case 456: // da5Nii6Njj7
        case 465: // da5Nii7Njj6
        case 546: // da6Nii5Njj7
        case 564: // da6Nii7Njj5
        case 645: // da7Nii5Njj6
        case 654: // da7Nii6Njj5
        case 738: // da8Nii4Njj9
        case 783: // da8Nii9Njj4
        case 837: // da9Nii4Njj8
        case 873: // da9Nii8Njj4
        case 378: // da4Nii8Njj9
            return ((-2.570000e+02) / (1.663200e+06)) * jacobi_det;
        case 888: // da9Nii9Njj9
        case 444: // da5Nii5Njj5
            return ((-1.000000e+00) / (4.950000e+03)) * jacobi_det;
        case 454: // da5Nii6Njj5
        case 544: // da6Nii5Njj5
        case 788: // da8Nii9Njj9
        case 878: // da9Nii8Njj9
        case 887: // da9Nii9Njj8
        case 445: // da5Nii5Njj6
            return ((1.000000e+00) / (9.900000e+03)) * jacobi_det;
        case 484: // da5Nii9Njj5
        case 488: // da5Nii9Njj9
        case 844: // da9Nii5Njj5
        case 848: // da9Nii5Njj9
        case 884: // da9Nii9Njj5
        case 448: // da5Nii5Njj9
            return ((-1.000000e+00) / (3.465000e+04)) * jacobi_det;
        case 494: // da5Nii10Njj5
        case 889: // da9Nii9Njj10
        case 898: // da9Nii10Njj9
        case 944: // da10Nii5Njj5
        case 988: // da10Nii9Njj9
        case 449: // da5Nii5Njj10
            return ((1.000000e+00) / (1.540000e+03)) * jacobi_det;
        case 545: // da6Nii5Njj6
        case 554: // da6Nii6Njj5
        case 778: // da8Nii8Njj9
        case 787: // da8Nii9Njj8
        case 877: // da9Nii8Njj8
        case 455: // da5Nii6Njj6
            return ((-6.700000e+01) / (8.316000e+05)) * jacobi_det;
        case 475: // da5Nii8Njj6
        case 547: // da6Nii5Njj8
        case 574: // da6Nii8Njj5
        case 578: // da6Nii8Njj9
        case 587: // da6Nii9Njj8
        case 745: // da8Nii5Njj6
        case 754: // da8Nii6Njj5
        case 758: // da8Nii6Njj9
        case 785: // da8Nii9Njj6
        case 857: // da9Nii6Njj8
        case 875: // da9Nii8Njj6
        case 457: // da5Nii6Njj8
            return ((-1.000000e+00) / (2.376000e+04)) * jacobi_det;
        case 495: // da5Nii10Njj6
        case 549: // da6Nii5Njj10
        case 594: // da6Nii10Njj5
        case 789: // da8Nii9Njj10
        case 798: // da8Nii10Njj9
        case 879: // da9Nii8Njj10
        case 897: // da9Nii10Njj8
        case 945: // da10Nii5Njj6
        case 954: // da10Nii6Njj5
        case 978: // da10Nii8Njj9
        case 987: // da10Nii9Njj8
        case 459: // da5Nii6Njj10
            return ((-1.000000e+00) / (3.080000e+03)) * jacobi_det;
        case 899: // da9Nii10Njj10
        case 949: // da10Nii5Njj10
        case 989: // da10Nii9Njj10
        case 994: // da10Nii10Njj5
        case 998: // da10Nii10Njj9
        case 499: // da5Nii10Njj10
            return ((-8.100000e+01) / (3.080000e+04)) * jacobi_det;
        case 575: // da6Nii8Njj6
        case 577: // da6Nii8Njj8
        case 755: // da8Nii6Njj6
        case 757: // da8Nii6Njj8
        case 775: // da8Nii8Njj6
        case 557: // da6Nii6Njj8
            return ((2.300000e+01) / (5.544000e+05)) * jacobi_det;
        case 999: // da10Nii10Njj10
            return ((6.561000e+03) / (6.160000e+04)) * jacobi_det;
		
		default:
			TDKP_GENERAL_EXCEPTION("invalid object index");	
	}
}
			
double Element2DTriangleHermiteCubic::evaluate_form_function(
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	
	
	const double& x = local_reference_element_coords[0];
	const double& y = local_reference_element_coords[1];	
	
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func;
	
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
		case 0: // F(x0)
			return 2.0*x*x*x + 2.0*y*y*y + 13.0*x*x*y + 13.0*x*y*y - 3.0*x*x - 3.0*y*y - 13.0*x*y + 1.0;
		case 1: // d/dx F(x0)
			return 1.0*x*x*x + 3.0*x*x*y + 2.0*x*y*y - 2.0*x*x - 3.0*x*y + 1.0*x;
		case 2: // d/dy F(x0)
			return 1.0*y*y*y + 2.0*x*x*y + 3.0*x*y*y - 2.0*y*y - 3.0*x*y + 1.0*y;
		case 3: // Fjm[0]
			return -2.0*x*x*x + 7.0*x*x*y + 7.0*x*y*y + 3.0*x*x - 7.0*x*y;
		case 4: // d/dx Fjm[0]
			return 1.0*x*x*x - 2.0*x*x*y - 2.0*x*y*y - 1.0*x*x + 2.0*x*y;
		case 5: // d/dy Fjm[0]
			return 2.0*x*x*y + 1.0*x*y*y - 1.0*x*y;
		case 6: // Fjm[1]
			return -2.0*y*y*y + 7.0*x*x*y + 7.0*x*y*y + 3.0*y*y - 7.0*x*y;
		case 7: // d/dx Fjm[1]
			return 1.0*x*x*y + 2.0*x*y*y - 1.0*x*y;
		case 8: // d/dy Fjm[1]
			return 1.0*y*y*y - 2.0*x*x*y - 2.0*x*y*y - 1.0*y*y + 2.0*x*y;
		case 9: // bubble function
			return -27.0*x*x*y - 27.0*x*y*y + 27.0*x*y;		
		default:
			TDKP_GENERAL_EXCEPTION("invalid object index");	
	}	
}

double Element2DTriangleHermiteCubic::evaluate_form_function_derivative(
	short diffop, short elem_shape_func, const double* local_reference_element_coords
) const {
	
	const double& x = local_reference_element_coords[0];
	const double& y = local_reference_element_coords[1];
	
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func 
	                       + 10 * diffop;
	// ---------------------------------------------
	// return appropriate integral
	// ---------------------------------------------		
	switch(obj_index) {
		
        case 0: // da1Nii1
            return (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]))) * ((((-1.000000e+00) + (y)) * (y) * (((6.000000e+00) * jm[2]) + ((-1.300000e+01) * jm[3]))) + (((x) * (x)) * (((1.300000e+01) * jm[2]) + ((-6.000000e+00) * jm[3]))) + ((x) * (((1.300000e+01) * ((-1.000000e+00) + ((2.000000e+00) * (y))) * jm[2]) + ((6.000000e+00) * jm[3]) + ((-2.600000e+01) * (y) * jm[3]))));
        case 1: // da1Nii2
            return (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * (((-3.000000e+00) * ((x) * (x)) * (jm[2] + ((-1.000000e+00) * jm[3]))) + (((1.000000e+00) + ((-3.000000e+00) * (y)) + ((2.000000e+00) * ((y) * (y)))) * jm[3]) + ((x) * (((3.000000e+00) * jm[2]) + ((-4.000000e+00) * (y) * jm[2]) + ((-4.000000e+00) * jm[3]) + ((6.000000e+00) * (y) * jm[3]))));
        case 2: // da1Nii3
            return (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]))) * ((((1.000000e+00) + ((2.000000e+00) * ((x) * (x))) + ((-4.000000e+00) * (y)) + ((3.000000e+00) * ((y) * (y))) + ((x) * ((-3.000000e+00) + ((6.000000e+00) * (y))))) * jm[2]) + (((3.000000e+00) + ((-4.000000e+00) * (x)) + ((-3.000000e+00) * (y))) * (y) * jm[3]));
        case 3: // da1Nii4
            return (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * (((7.000000e+00) * ((-1.000000e+00) + (y)) * (y) * jm[3]) + ((-1.000000e+00) * ((x) * (x)) * (((7.000000e+00) * jm[2]) + ((6.000000e+00) * jm[3]))) + ((x) * ((((7.000000e+00) + ((-1.400000e+01) * (y))) * jm[2]) + ((2.000000e+00) * ((3.000000e+00) + ((7.000000e+00) * (y))) * jm[3]))));
        case 4: // da1Nii5
            return (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]))) * (((2.000000e+00) * ((-1.000000e+00) + (y)) * (y) * jm[3]) + ((-1.000000e+00) * ((x) * (x)) * (((2.000000e+00) * jm[2]) + ((3.000000e+00) * jm[3]))) + ((2.000000e+00) * (x) * (jm[2] + ((-2.000000e+00) * (y) * jm[2]) + jm[3] + ((2.000000e+00) * (y) * jm[3]))));
        case 5: // da1Nii6
            return (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * (((-2.000000e+00) * ((x) * (x)) * jm[2]) + (((-1.000000e+00) + (y)) * (y) * jm[3]) + ((x) * (jm[2] + ((-2.000000e+00) * (y) * jm[2]) + ((4.000000e+00) * (y) * jm[3]))));
        case 6: // da1Nii7
            return (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * (((-7.000000e+00) * ((x) * (x)) * jm[2]) + (((-1.000000e+00) + (y)) * (y) * (((6.000000e+00) * jm[2]) + ((7.000000e+00) * jm[3]))) + ((7.000000e+00) * (x) * (jm[2] + ((-2.000000e+00) * (y) * jm[2]) + ((2.000000e+00) * (y) * jm[3]))));
        case 7: // da1Nii8
            return (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * (((-1.000000e+00) * ((x) * (x)) * jm[2]) + ((y) * ((-1.000000e+00) + ((2.000000e+00) * (y))) * jm[3]) + ((x) * (jm[2] + ((-4.000000e+00) * (y) * jm[2]) + ((2.000000e+00) * (y) * jm[3]))));
        case 8: // da1Nii9
            return (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]))) * (((-2.000000e+00) * ((x) * (x)) * jm[2]) + ((y) * ((((-2.000000e+00) + ((3.000000e+00) * (y))) * jm[2]) + ((2.000000e+00) * ((-1.000000e+00) + (y)) * jm[3]))) + ((x) * ((((2.000000e+00) + ((-4.000000e+00) * (y))) * jm[2]) + ((4.000000e+00) * (y) * jm[3]))));
        case 9: // da1Nii10
            return (2.700000e+01) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * ((((x) * (x)) * jm[2]) + ((-1.000000e+00) * ((-1.000000e+00) + (y)) * (y) * jm[3]) + ((x) * ((((-1.000000e+00) + ((2.000000e+00) * (y))) * jm[2]) + ((-2.000000e+00) * (y) * jm[3]))));
        case 10: // da2Nii1
            return ((((x) * (x)) * (((1.300000e+01) * jm[0]) + ((-6.000000e+00) * jm[1]))) + ((((6.000000e+00) * jm[0]) + ((-1.300000e+01) * jm[1])) * ((-1.000000e+00) + (y)) * (y)) + ((x) * (((6.000000e+00) * jm[1]) + ((-2.600000e+01) * jm[1] * (y)) + ((1.300000e+01) * jm[0] * ((-1.000000e+00) + ((2.000000e+00) * (y))))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 11: // da2Nii2
            return (((3.000000e+00) * ((x) * (x)) * (jm[0] + ((-1.000000e+00) * jm[1]))) + (jm[1] * ((-1.000000e+00) + ((3.000000e+00) * (y)) + ((-2.000000e+00) * ((y) * (y))))) + ((x) * (((4.000000e+00) * jm[1]) + ((-6.000000e+00) * jm[1] * (y)) + (jm[0] * ((-3.000000e+00) + ((4.000000e+00) * (y))))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 12: // da2Nii3
            return ((jm[1] * ((3.000000e+00) + ((-4.000000e+00) * (x)) + ((-3.000000e+00) * (y))) * (y)) + (jm[0] * ((1.000000e+00) + ((2.000000e+00) * ((x) * (x))) + ((-4.000000e+00) * (y)) + ((3.000000e+00) * ((y) * (y))) + ((x) * ((-3.000000e+00) + ((6.000000e+00) * (y))))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 13: // da2Nii4
            return ((((x) * (x)) * (((7.000000e+00) * jm[0]) + ((6.000000e+00) * jm[1]))) + ((-7.000000e+00) * jm[1] * ((-1.000000e+00) + (y)) * (y)) + ((x) * (((7.000000e+00) * jm[0] * ((-1.000000e+00) + ((2.000000e+00) * (y)))) + ((-2.000000e+00) * jm[1] * ((3.000000e+00) + ((7.000000e+00) * (y))))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 14: // da2Nii5
            return (((-1.000000e+00) * ((x) * (x)) * (((2.000000e+00) * jm[0]) + ((3.000000e+00) * jm[1]))) + ((2.000000e+00) * jm[1] * ((-1.000000e+00) + (y)) * (y)) + ((2.000000e+00) * (x) * (jm[0] + jm[1] + ((-2.000000e+00) * jm[0] * (y)) + ((2.000000e+00) * jm[1] * (y))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 15: // da2Nii6
            return (-1.000000e+00) * (((-2.000000e+00) * ((x) * (x)) * jm[0]) + (jm[1] * ((-1.000000e+00) + (y)) * (y)) + ((x) * (jm[0] + ((-2.000000e+00) * jm[0] * (y)) + ((4.000000e+00) * jm[1] * (y))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 16: // da2Nii7
            return (((7.000000e+00) * ((x) * (x)) * jm[0]) + ((-1.000000e+00) * (((6.000000e+00) * jm[0]) + ((7.000000e+00) * jm[1])) * ((-1.000000e+00) + (y)) * (y)) + ((-7.000000e+00) * (x) * (jm[0] + ((-2.000000e+00) * jm[0] * (y)) + ((2.000000e+00) * jm[1] * (y))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 17: // da2Nii8
            return ((((x) * (x)) * jm[0]) + (jm[1] * ((1.000000e+00) + ((-2.000000e+00) * (y))) * (y)) + ((-1.000000e+00) * (x) * (jm[0] + ((-4.000000e+00) * jm[0] * (y)) + ((2.000000e+00) * jm[1] * (y))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 18: // da2Nii9
            return (((-2.000000e+00) * ((x) * (x)) * jm[0]) + ((x) * ((jm[0] * ((2.000000e+00) + ((-4.000000e+00) * (y)))) + ((4.000000e+00) * jm[1] * (y)))) + ((y) * (((2.000000e+00) * jm[1] * ((-1.000000e+00) + (y))) + (jm[0] * ((-2.000000e+00) + ((3.000000e+00) * (y))))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 19: // da2Nii10
            return (-2.700000e+01) * ((((x) * (x)) * jm[0]) + ((-1.000000e+00) * jm[1] * ((-1.000000e+00) + (y)) * (y)) + ((-1.000000e+00) * (x) * (jm[0] + ((-2.000000e+00) * jm[0] * (y)) + ((2.000000e+00) * jm[1] * (y))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
		
		default:
			TDKP_GENERAL_EXCEPTION("invalid object index");	
	}	
}					



}
