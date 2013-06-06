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

#include "Element2DRect2nd.h"

namespace tdkp {

// ---------------------------------------------------
// all analytical integrals have been obtained via 
// mathematica (1344 integrals (sic!) and the mathematica output
// has been parsed into c++ code
// ---------------------------------------------------

Element2DRect2nd::Element2DRect2nd(unsigned int index_)
: Element2DRectBase(2,8,4)
{
	this->index_global    = index_;
	this->vertex_offset   = 2; // every other node is a vertex
}

void Element2DRect2nd::set_corner_node(unsigned short corner_id, Node* node) {
	TDKP_ASSERT(corner_id < 4, "");
	TDKP_ASSERT(nodes.size() > corner_id, "");
	TDKP_ASSERT(nodes[corner_id * vertex_offset] == 0, "");
	nodes[corner_id * vertex_offset] = node;	
}		

const Node& Element2DRect2nd::get_corner_node(unsigned short corner_idx) const {
	TDKP_ASSERT(corner_idx < 4, "");
	TDKP_ASSERT(nodes.size() > corner_idx, "");
	return *nodes[corner_idx * vertex_offset];	
}

void Element2DRect2nd::get_additional_node_locator(unsigned int additional_node_idx, AdditionalNodeLocation& location_type, vector<unsigned int>& vertex_indices, vector<double>& coords, unsigned int& tag) const {

	TDKP_ASSERT(additional_node_idx < 4, "");
	location_type = edge_node;
	tag           = 0;
	// ------------------------------------------
	// find boundary
	// ------------------------------------------
	int vidx0 = nodes[additional_node_idx * 2]->get_index_vertex();
	int vidx1 = nodes[(additional_node_idx * 2 + 2) % 8]->get_index_vertex();
	vertex_indices.resize(2);
	vertex_indices[0] = vidx0;
	vertex_indices[1] = vidx1;
	// -------------------------------------------
	// set coords
	// -------------------------------------------
	coords.resize(2);
	for(unsigned int ii = 0; ii < 2; ii++) {
		coords[ii] = 0.5 * (nodes[additional_node_idx * 2]->get_coord(ii) + nodes[(additional_node_idx * 2 + 2) % 8]->get_coord(ii));
	}	
}


void Element2DRect2nd::set_additional_node(unsigned int additional_node_idx, Node* node) {
	unsigned int lidx = additional_node_idx * 2 + 1;
	TDKP_BOUNDS_ASSERT(lidx < nodes.size(), "");
	TDKP_BOUNDS_ASSERT(nodes[lidx] == 0, "");
	TDKP_BOUNDS_ASSERT(test_additional_node(*this, *node, additional_node_idx), "wrong additional node supplied: " << additional_node_idx);
	nodes[lidx] = node;
	// -----------------------------------
	// set linear interpolation coefficents for additional node
	// -----------------------------------
	if(node->get_num_contributions() == 0) {
		// node sits in middle of edge, so halve the value of each vertex
		node->set_contribution(nodes[lidx - 1]->get_index_vertex(), 0.5);
		node->set_contribution(nodes[(lidx + 1) % 8]->get_index_vertex(), 0.5);	
	}	
}




double Element2DRect2nd::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {

	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 8 * elem_shape_func_1
	                       + 8 * 8 * diffop_2
	                       + 8 * 8 * 2 * diffop_1;
	
	// ----------------------------------------------
	// return analytical integral
	// ----------------------------------------------
	switch(obj_index) {
        case 36: // da1db1Nii5Njj5
        case 0: // da1db1Nii1Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (((5.200000e+01) * ((ay) * (ay))) + ((-8.500000e+01) * (ay) * (by)) + ((5.200000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 8: // da1db1Nii2Njj1
        case 37: // da1db1Nii5Njj6
        case 44: // da1db1Nii6Njj5
        case 1: // da1db1Nii1Njj2
            return ((1.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ay) * (ay))) + ((2.000000e+01) * (ay) * (by)) + ((-4.000000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 16: // da1db1Nii3Njj1
        case 38: // da1db1Nii5Njj7
        case 52: // da1db1Nii7Njj5
        case 2: // da1db1Nii1Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (((1.700000e+01) * ((ay) * (ay))) + ((2.800000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 24: // da1db1Nii4Njj1
        case 39: // da1db1Nii5Njj8
        case 60: // da1db1Nii8Njj5
        case 3: // da1db1Nii1Njj4
            return (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by)))) * (((2.000000e+01) * ((ay) * (ay))) + ((-1.000000e+01) * (ay) * (by)) + ((3.000000e+00) * ((by) * (by))));
        case 32: // da1db1Nii5Njj1
        case 4: // da1db1Nii1Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (((2.300000e+01) * ((ay) * (ay))) + ((-3.500000e+01) * (ay) * (by)) + ((2.300000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 12: // da1db1Nii2Njj5
        case 33: // da1db1Nii5Njj2
        case 40: // da1db1Nii6Njj1
        case 5: // da1db1Nii1Njj6
            return (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by)))) * (((3.000000e+00) * ((ay) * (ay))) + ((-1.000000e+01) * (ay) * (by)) + ((2.000000e+01) * ((by) * (by))));
        case 20: // da1db1Nii3Njj5
        case 34: // da1db1Nii5Njj3
        case 48: // da1db1Nii7Njj1
        case 6: // da1db1Nii1Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (((2.800000e+01) * ((ay) * (ay))) + ((1.700000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 28: // da1db1Nii4Njj5
        case 35: // da1db1Nii5Njj4
        case 56: // da1db1Nii8Njj1
        case 7: // da1db1Nii1Njj8
            return (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by)))) * (((4.000000e+01) * ((ay) * (ay))) + ((-2.000000e+01) * (ay) * (by)) + ((-3.000000e+00) * ((by) * (by))));
        case 45: // da1db1Nii6Njj6
        case 9: // da1db1Nii2Njj2
            return ((8.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ay) * (ay))) + ((1.000000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 17: // da1db1Nii3Njj2
        case 46: // da1db1Nii6Njj7
        case 53: // da1db1Nii7Njj6
        case 10: // da1db1Nii2Njj3
            return ((1.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ay) * (ay))) + ((-2.000000e+01) * (ay) * (by)) + ((-4.000000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 25: // da1db1Nii4Njj2
        case 47: // da1db1Nii6Njj8
        case 61: // da1db1Nii8Njj6
        case 11: // da1db1Nii2Njj4
            return (8.000000e+00) * (ay) * (by) * (1.0 / (((-9.000000e+00) * (ay) * (bx)) + ((9.000000e+00) * (ax) * (by))));
        case 41: // da1db1Nii6Njj2
        case 13: // da1db1Nii2Njj6
            return ((-8.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ay) * (ay))) + ((-5.000000e+00) * ((by) * (by)))) * (1.0 / (element_volume));
        case 21: // da1db1Nii3Njj6
        case 42: // da1db1Nii6Njj3
        case 49: // da1db1Nii7Njj2
        case 14: // da1db1Nii2Njj7
            return (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by)))) * (((3.000000e+00) * ((ay) * (ay))) + ((1.000000e+01) * (ay) * (by)) + ((2.000000e+01) * ((by) * (by))));
        case 29: // da1db1Nii4Njj6
        case 43: // da1db1Nii6Njj4
        case 57: // da1db1Nii8Njj2
        case 15: // da1db1Nii2Njj8
            return (8.000000e+00) * (ay) * (by) * (1.0 / (((9.000000e+00) * (ay) * (bx)) + ((-9.000000e+00) * (ax) * (by))));
        case 54: // da1db1Nii7Njj7
        case 18: // da1db1Nii3Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (((5.200000e+01) * ((ay) * (ay))) + ((8.500000e+01) * (ay) * (by)) + ((5.200000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 26: // da1db1Nii4Njj3
        case 55: // da1db1Nii7Njj8
        case 62: // da1db1Nii8Njj7
        case 19: // da1db1Nii3Njj4
            return (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by)))) * (((4.000000e+01) * ((ay) * (ay))) + ((2.000000e+01) * (ay) * (by)) + ((-3.000000e+00) * ((by) * (by))));
        case 50: // da1db1Nii7Njj3
        case 22: // da1db1Nii3Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (((2.300000e+01) * ((ay) * (ay))) + ((3.500000e+01) * (ay) * (by)) + ((2.300000e+01) * ((by) * (by)))) * (1.0 / (element_volume));
        case 30: // da1db1Nii4Njj7
        case 51: // da1db1Nii7Njj4
        case 58: // da1db1Nii8Njj3
        case 23: // da1db1Nii3Njj8
            return (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by)))) * (((2.000000e+01) * ((ay) * (ay))) + ((1.000000e+01) * (ay) * (by)) + ((3.000000e+00) * ((by) * (by))));
        case 63: // da1db1Nii8Njj8
        case 27: // da1db1Nii4Njj4
            return ((8.000000e+00) / (4.500000e+01)) * (((1.000000e+01) * ((ay) * (ay))) + ((3.000000e+00) * ((by) * (by)))) * (1.0 / (element_volume));
        case 59: // da1db1Nii8Njj4
        case 31: // da1db1Nii4Njj8
            return ((8.000000e+00) / (4.500000e+01)) * (((5.000000e+00) * ((ay) * (ay))) + ((-3.000000e+00) * ((by) * (by)))) * (1.0 / (element_volume));
        case 100: // da1db2Nii5Njj5
        case 128: // da2db1Nii1Njj1
        case 164: // da2db1Nii5Njj5
        case 64: // da1db2Nii1Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (((-1.040000e+02) * (ax) * (ay)) + ((8.500000e+01) * (ay) * (bx)) + ((8.500000e+01) * (ax) * (by)) + ((-1.040000e+02) * (bx) * (by))) * (1.0 / (element_volume));
        case 101: // da1db2Nii5Njj6
        case 136: // da2db1Nii2Njj1
        case 172: // da2db1Nii6Njj5
        case 65: // da1db2Nii1Njj2
            return ((1.000000e+00) / (4.500000e+01)) * (((-3.000000e+00) * (ax) * (ay)) + ((-2.500000e+01) * (ay) * (bx)) + ((5.000000e+00) * (ax) * (by)) + ((4.000000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 102: // da1db2Nii5Njj7
        case 144: // da2db1Nii3Njj1
        case 180: // da2db1Nii7Njj5
        case 66: // da1db2Nii1Njj3
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / (element_volume)) * (((3.400000e+01) * (ax) * (ay)) + ((5.600000e+01) * (bx) * (by)) + ((1.500000e+01) * (element_volume)));
        case 88: // da1db2Nii4Njj1
        case 103: // da1db2Nii5Njj8
        case 124: // da1db2Nii8Njj5
        case 131: // da2db1Nii1Njj4
        case 152: // da2db1Nii4Njj1
        case 167: // da2db1Nii5Njj8
        case 188: // da2db1Nii8Njj5
        case 67: // da1db2Nii1Njj4
            return ((1.000000e+00) / (4.500000e+01)) * (((2.000000e+01) * (ax) * (ay)) + ((-5.000000e+00) * (ay) * (bx)) + ((-5.000000e+00) * (ax) * (by)) + ((3.000000e+00) * (bx) * (by))) * (1.0 / (element_volume));
        case 96: // da1db2Nii5Njj1
        case 132: // da2db1Nii1Njj5
        case 160: // da2db1Nii5Njj1
        case 68: // da1db2Nii1Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (((-4.600000e+01) * (ax) * (ay)) + ((3.500000e+01) * (ay) * (bx)) + ((3.500000e+01) * (ax) * (by)) + ((-4.600000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 76: // da1db2Nii2Njj5
        case 97: // da1db2Nii5Njj2
        case 104: // da1db2Nii6Njj1
        case 133: // da2db1Nii1Njj6
        case 140: // da2db1Nii2Njj5
        case 161: // da2db1Nii5Njj2
        case 168: // da2db1Nii6Njj1
        case 69: // da1db2Nii1Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * (ax) * (ay)) + ((-5.000000e+00) * (ay) * (bx)) + ((-5.000000e+00) * (ax) * (by)) + ((2.000000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 98: // da1db2Nii5Njj3
        case 148: // da2db1Nii3Njj5
        case 176: // da2db1Nii7Njj1
        case 70: // da1db2Nii1Njj7
            return ((-1.000000e+00) / (1.800000e+02)) * (((5.600000e+01) * (ax) * (ay)) + ((3.400000e+01) * (bx) * (by)) + ((-1.500000e+01) * (element_volume))) * (1.0 / (element_volume));
        case 99: // da1db2Nii5Njj4
        case 156: // da2db1Nii4Njj5
        case 184: // da2db1Nii8Njj1
        case 71: // da1db2Nii1Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (((4.000000e+01) * (ax) * (ay)) + ((5.000000e+00) * (ay) * (bx)) + ((-2.500000e+01) * (ax) * (by)) + ((-3.000000e+00) * (bx) * (by))) * (1.0 / (element_volume));
        case 108: // da1db2Nii6Njj5
        case 129: // da2db1Nii1Njj2
        case 165: // da2db1Nii5Njj6
        case 72: // da1db2Nii2Njj1
            return ((1.000000e+00) / (4.500000e+01)) * (((5.000000e+00) * (bx) * ((ay) + ((8.000000e+00) * (by)))) + ((-1.000000e+00) * (ax) * (((3.000000e+00) * (ay)) + ((2.500000e+01) * (by))))) * (1.0 / (element_volume));
        case 109: // da1db2Nii6Njj6
        case 137: // da2db1Nii2Njj2
        case 173: // da2db1Nii6Njj6
        case 73: // da1db2Nii2Njj2
            return ((-8.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * (ax) * (ay)) + ((1.000000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 110: // da1db2Nii6Njj7
        case 145: // da2db1Nii3Njj2
        case 181: // da2db1Nii7Njj6
        case 74: // da1db2Nii2Njj3
            return ((1.000000e+00) / (4.500000e+01)) * (((-3.000000e+00) * (ax) * (ay)) + ((-5.000000e+00) * (ay) * (bx)) + ((2.500000e+01) * (ax) * (by)) + ((4.000000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 89: // da1db2Nii4Njj2
        case 111: // da1db2Nii6Njj8
        case 125: // da1db2Nii8Njj6
        case 139: // da2db1Nii2Njj4
        case 153: // da2db1Nii4Njj2
        case 175: // da2db1Nii6Njj8
        case 189: // da2db1Nii8Njj6
        case 75: // da1db2Nii2Njj4
            return ((-4.000000e+00) / (9.000000e+00)) * (((ay) * (bx)) + ((ax) * (by))) * (1.0 / (element_volume));
        case 105: // da1db2Nii6Njj2
        case 141: // da2db1Nii2Njj6
        case 169: // da2db1Nii6Njj2
        case 77: // da1db2Nii2Njj6
            return ((8.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * (ax) * (ay)) + ((-5.000000e+00) * (bx) * (by))) * (1.0 / (element_volume));
        case 85: // da1db2Nii3Njj6
        case 106: // da1db2Nii6Njj3
        case 113: // da1db2Nii7Njj2
        case 142: // da2db1Nii2Njj7
        case 149: // da2db1Nii3Njj6
        case 170: // da2db1Nii6Njj3
        case 177: // da2db1Nii7Njj2
        case 78: // da1db2Nii2Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * (ax) * (ay)) + ((5.000000e+00) * (ay) * (bx)) + ((5.000000e+00) * (ax) * (by)) + ((2.000000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 93: // da1db2Nii4Njj6
        case 107: // da1db2Nii6Njj4
        case 121: // da1db2Nii8Njj2
        case 143: // da2db1Nii2Njj8
        case 157: // da2db1Nii4Njj6
        case 171: // da2db1Nii6Njj4
        case 185: // da2db1Nii8Njj2
        case 79: // da1db2Nii2Njj8
            return (4.000000e+00) * (((ay) * (bx)) + ((ax) * (by))) * (1.0 / (((-9.000000e+00) * (ay) * (bx)) + ((9.000000e+00) * (ax) * (by))));
        case 116: // da1db2Nii7Njj5
        case 130: // da2db1Nii1Njj3
        case 166: // da2db1Nii5Njj7
        case 80: // da1db2Nii3Njj1
            return ((-1.000000e+00) / (1.800000e+02)) * (((3.400000e+01) * (ax) * (ay)) + ((5.600000e+01) * (bx) * (by)) + ((-1.500000e+01) * (element_volume))) * (1.0 / (element_volume));
        case 117: // da1db2Nii7Njj6
        case 138: // da2db1Nii2Njj3
        case 174: // da2db1Nii6Njj7
        case 81: // da1db2Nii3Njj2
            return ((1.000000e+00) / (4.500000e+01)) * (((-3.000000e+00) * (ax) * (ay)) + ((2.500000e+01) * (ay) * (bx)) + ((-5.000000e+00) * (ax) * (by)) + ((4.000000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 118: // da1db2Nii7Njj7
        case 146: // da2db1Nii3Njj3
        case 182: // da2db1Nii7Njj7
        case 82: // da1db2Nii3Njj3
            return ((-1.000000e+00) / (1.800000e+02)) * (((1.040000e+02) * (ax) * (ay)) + ((8.500000e+01) * (ay) * (bx)) + ((8.500000e+01) * (ax) * (by)) + ((1.040000e+02) * (bx) * (by))) * (1.0 / (element_volume));
        case 119: // da1db2Nii7Njj8
        case 154: // da2db1Nii4Njj3
        case 190: // da2db1Nii8Njj7
        case 83: // da1db2Nii3Njj4
            return ((1.000000e+00) / (4.500000e+01)) * (((4.000000e+01) * (ax) * (ay)) + ((-5.000000e+00) * (ay) * (bx)) + ((2.500000e+01) * (ax) * (by)) + ((-3.000000e+00) * (bx) * (by))) * (1.0 / (element_volume));
        case 112: // da1db2Nii7Njj1
        case 134: // da2db1Nii1Njj7
        case 162: // da2db1Nii5Njj3
        case 84: // da1db2Nii3Njj5
            return ((-1.000000e+00) / (1.800000e+02)) * (1.0 / (element_volume)) * (((5.600000e+01) * (ax) * (ay)) + ((3.400000e+01) * (bx) * (by)) + ((1.500000e+01) * (element_volume)));
        case 114: // da1db2Nii7Njj3
        case 150: // da2db1Nii3Njj7
        case 178: // da2db1Nii7Njj3
        case 86: // da1db2Nii3Njj7
            return ((-1.000000e+00) / (1.800000e+02)) * (((4.600000e+01) * (ax) * (ay)) + ((3.500000e+01) * (ay) * (bx)) + ((3.500000e+01) * (ax) * (by)) + ((4.600000e+01) * (bx) * (by))) * (1.0 / (element_volume));
        case 94: // da1db2Nii4Njj7
        case 115: // da1db2Nii7Njj4
        case 122: // da1db2Nii8Njj3
        case 151: // da2db1Nii3Njj8
        case 158: // da2db1Nii4Njj7
        case 179: // da2db1Nii7Njj4
        case 186: // da2db1Nii8Njj3
        case 87: // da1db2Nii3Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (((5.000000e+00) * (ay) * (bx)) + ((3.000000e+00) * (bx) * (by)) + ((5.000000e+00) * (ax) * (((4.000000e+00) * (ay)) + (by)))) * (1.0 / (element_volume));
        case 126: // da1db2Nii8Njj7
        case 147: // da2db1Nii3Njj4
        case 183: // da2db1Nii7Njj8
        case 90: // da1db2Nii4Njj3
            return ((1.000000e+00) / (4.500000e+01)) * (((4.000000e+01) * (ax) * (ay)) + ((2.500000e+01) * (ay) * (bx)) + ((-5.000000e+00) * (ax) * (by)) + ((-3.000000e+00) * (bx) * (by))) * (1.0 / (element_volume));
        case 127: // da1db2Nii8Njj8
        case 155: // da2db1Nii4Njj4
        case 191: // da2db1Nii8Njj8
        case 91: // da1db2Nii4Njj4
            return ((-8.000000e+00) / (4.500000e+01)) * (((1.000000e+01) * (ax) * (ay)) + ((3.000000e+00) * (bx) * (by))) * (1.0 / (element_volume));
        case 120: // da1db2Nii8Njj1
        case 135: // da2db1Nii1Njj8
        case 163: // da2db1Nii5Njj4
        case 92: // da1db2Nii4Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (((-2.500000e+01) * (ay) * (bx)) + ((-3.000000e+00) * (bx) * (by)) + ((5.000000e+00) * (ax) * (((8.000000e+00) * (ay)) + (by)))) * (1.0 / (element_volume));
        case 123: // da1db2Nii8Njj4
        case 159: // da2db1Nii4Njj8
        case 187: // da2db1Nii8Njj4
        case 95: // da1db2Nii4Njj8
            return ((-8.000000e+00) / (4.500000e+01)) * (((5.000000e+00) * (ax) * (ay)) + ((-3.000000e+00) * (bx) * (by))) * (1.0 / (element_volume));
        case 228: // da2db2Nii5Njj5
        case 192: // da2db2Nii1Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (((5.200000e+01) * ((ax) * (ax))) + ((-8.500000e+01) * (ax) * (bx)) + ((5.200000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 200: // da2db2Nii2Njj1
        case 229: // da2db2Nii5Njj6
        case 236: // da2db2Nii6Njj5
        case 193: // da2db2Nii1Njj2
            return ((1.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ax) * (ax))) + ((2.000000e+01) * (ax) * (bx)) + ((-4.000000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 208: // da2db2Nii3Njj1
        case 230: // da2db2Nii5Njj7
        case 244: // da2db2Nii7Njj5
        case 194: // da2db2Nii1Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (((1.700000e+01) * ((ax) * (ax))) + ((2.800000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 216: // da2db2Nii4Njj1
        case 231: // da2db2Nii5Njj8
        case 252: // da2db2Nii8Njj5
        case 195: // da2db2Nii1Njj4
            return (((2.000000e+01) * ((ax) * (ax))) + ((-1.000000e+01) * (ax) * (bx)) + ((3.000000e+00) * ((bx) * (bx)))) * (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by))));
        case 224: // da2db2Nii5Njj1
        case 196: // da2db2Nii1Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (((2.300000e+01) * ((ax) * (ax))) + ((-3.500000e+01) * (ax) * (bx)) + ((2.300000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 204: // da2db2Nii2Njj5
        case 225: // da2db2Nii5Njj2
        case 232: // da2db2Nii6Njj1
        case 197: // da2db2Nii1Njj6
            return (((3.000000e+00) * ((ax) * (ax))) + ((-1.000000e+01) * (ax) * (bx)) + ((2.000000e+01) * ((bx) * (bx)))) * (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by))));
        case 212: // da2db2Nii3Njj5
        case 226: // da2db2Nii5Njj3
        case 240: // da2db2Nii7Njj1
        case 198: // da2db2Nii1Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (((2.800000e+01) * ((ax) * (ax))) + ((1.700000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 220: // da2db2Nii4Njj5
        case 227: // da2db2Nii5Njj4
        case 248: // da2db2Nii8Njj1
        case 199: // da2db2Nii1Njj8
            return (((4.000000e+01) * ((ax) * (ax))) + ((-2.000000e+01) * (ax) * (bx)) + ((-3.000000e+00) * ((bx) * (bx)))) * (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by))));
        case 237: // da2db2Nii6Njj6
        case 201: // da2db2Nii2Njj2
            return ((8.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ax) * (ax))) + ((1.000000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 209: // da2db2Nii3Njj2
        case 238: // da2db2Nii6Njj7
        case 245: // da2db2Nii7Njj6
        case 202: // da2db2Nii2Njj3
            return ((1.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ax) * (ax))) + ((-2.000000e+01) * (ax) * (bx)) + ((-4.000000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 217: // da2db2Nii4Njj2
        case 239: // da2db2Nii6Njj8
        case 253: // da2db2Nii8Njj6
        case 203: // da2db2Nii2Njj4
            return (8.000000e+00) * (ax) * (bx) * (1.0 / (((-9.000000e+00) * (ay) * (bx)) + ((9.000000e+00) * (ax) * (by))));
        case 233: // da2db2Nii6Njj2
        case 205: // da2db2Nii2Njj6
            return ((-8.000000e+00) / (4.500000e+01)) * (((3.000000e+00) * ((ax) * (ax))) + ((-5.000000e+00) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 213: // da2db2Nii3Njj6
        case 234: // da2db2Nii6Njj3
        case 241: // da2db2Nii7Njj2
        case 206: // da2db2Nii2Njj7
            return (((3.000000e+00) * ((ax) * (ax))) + ((1.000000e+01) * (ax) * (bx)) + ((2.000000e+01) * ((bx) * (bx)))) * (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by))));
        case 221: // da2db2Nii4Njj6
        case 235: // da2db2Nii6Njj4
        case 249: // da2db2Nii8Njj2
        case 207: // da2db2Nii2Njj8
            return (8.000000e+00) * (ax) * (bx) * (1.0 / (((9.000000e+00) * (ay) * (bx)) + ((-9.000000e+00) * (ax) * (by))));
        case 246: // da2db2Nii7Njj7
        case 210: // da2db2Nii3Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (((5.200000e+01) * ((ax) * (ax))) + ((8.500000e+01) * (ax) * (bx)) + ((5.200000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 218: // da2db2Nii4Njj3
        case 247: // da2db2Nii7Njj8
        case 254: // da2db2Nii8Njj7
        case 211: // da2db2Nii3Njj4
            return (((4.000000e+01) * ((ax) * (ax))) + ((2.000000e+01) * (ax) * (bx)) + ((-3.000000e+00) * ((bx) * (bx)))) * (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by))));
        case 242: // da2db2Nii7Njj3
        case 214: // da2db2Nii3Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (((2.300000e+01) * ((ax) * (ax))) + ((3.500000e+01) * (ax) * (bx)) + ((2.300000e+01) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 222: // da2db2Nii4Njj7
        case 243: // da2db2Nii7Njj4
        case 250: // da2db2Nii8Njj3
        case 215: // da2db2Nii3Njj8
            return (((2.000000e+01) * ((ax) * (ax))) + ((1.000000e+01) * (ax) * (bx)) + ((3.000000e+00) * ((bx) * (bx)))) * (1.0 / (((4.500000e+01) * (ay) * (bx)) + ((-4.500000e+01) * (ax) * (by))));
        case 255: // da2db2Nii8Njj8
        case 219: // da2db2Nii4Njj4
            return ((8.000000e+00) / (4.500000e+01)) * (((1.000000e+01) * ((ax) * (ax))) + ((3.000000e+00) * ((bx) * (bx)))) * (1.0 / (element_volume));
        case 251: // da2db2Nii8Njj4
        case 223: // da2db2Nii4Njj8
            return ((8.000000e+00) / (4.500000e+01)) * (((5.000000e+00) * ((ax) * (ax))) + ((-3.000000e+00) * ((bx) * (bx)))) * (1.0 / (element_volume));
		
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");
	}	
	
}

double Element2DRect2nd::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {
	
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 8 * elem_shape_func_1
	                       + 8 * 8 * diffop;
	                       
	
	// ----------------------------------------------
	// return analytical integral
	// ----------------------------------------------
	switch(obj_index) {
        case 0: // da1Nii1Njj1
            return ((1.000000e+00) / (1.500000e+01)) * ((ay) + ((-1.000000e+00) * (by)));
        case 1: // da1Nii1Njj2
            return (((1.300000e+01) / (9.000000e+01)) * (ay)) + (((-1.000000e+00) / (9.000000e+00)) * (by));
        case 2: // da1Nii1Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (((-3.000000e+00) * (ay)) + ((8.000000e+00) * (by)));
        case 51: // da1Nii7Njj4
        case 58: // da1Nii8Njj3
        case 60: // da1Nii8Njj5
        case 3: // da1Nii1Njj4
            return ((7.000000e+00) / (9.000000e+01)) * (by);
        case 4: // da1Nii1Njj5
            return ((1.000000e+00) / (6.000000e+01)) * (((-1.000000e+00) * (ay)) + (by));
        case 12: // da1Nii2Njj5
        case 14: // da1Nii2Njj7
        case 21: // da1Nii3Njj6
        case 5: // da1Nii1Njj6
            return ((-7.000000e+00) / (9.000000e+01)) * (ay);
        case 6: // da1Nii1Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (((-8.000000e+00) * (ay)) + ((3.000000e+00) * (by)));
        case 7: // da1Nii1Njj8
            return (((1.000000e+00) / (9.000000e+00)) * (ay)) + (((-1.300000e+01) / (9.000000e+01)) * (by));
        case 8: // da1Nii2Njj1
            return (((-7.000000e+00) / (9.000000e+01)) * (ay)) + (((1.000000e+00) / (9.000000e+00)) * (by));
        case 13: // da1Nii2Njj6
        case 9: // da1Nii2Njj2
            return ((4.000000e+00) / (1.500000e+01)) * (ay);
        case 10: // da1Nii2Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (((-7.000000e+00) * (ay)) + ((-1.000000e+01) * (by)));
        case 61: // da1Nii8Njj6
        case 11: // da1Nii2Njj4
            return ((2.000000e+00) / (9.000000e+00)) * ((ay) + ((-1.000000e+00) * (by)));
        case 29: // da1Nii4Njj6
        case 15: // da1Nii2Njj8
            return ((2.000000e+00) / (9.000000e+00)) * ((ay) + (by));
        case 16: // da1Nii3Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (((-3.000000e+00) * (ay)) + ((-8.000000e+00) * (by)));
        case 17: // da1Nii3Njj2
            return (((1.300000e+01) / (9.000000e+01)) * (ay)) + (((1.000000e+00) / (9.000000e+00)) * (by));
        case 18: // da1Nii3Njj3
            return ((1.000000e+00) / (1.500000e+01)) * ((ay) + (by));
        case 19: // da1Nii3Njj4
            return (((1.000000e+00) / (9.000000e+00)) * (ay)) + (((1.300000e+01) / (9.000000e+01)) * (by));
        case 20: // da1Nii3Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (((-8.000000e+00) * (ay)) + ((-3.000000e+00) * (by)));
        case 22: // da1Nii3Njj7
            return ((1.000000e+00) / (6.000000e+01)) * (((-1.000000e+00) * (ay)) + ((-1.000000e+00) * (by)));
        case 24: // da1Nii4Njj1
        case 30: // da1Nii4Njj7
        case 39: // da1Nii5Njj8
        case 23: // da1Nii3Njj8
            return ((-7.000000e+00) / (9.000000e+01)) * (by);
        case 47: // da1Nii6Njj8
        case 25: // da1Nii4Njj2
            return ((-2.000000e+00) / (9.000000e+00)) * ((ay) + ((-1.000000e+00) * (by)));
        case 26: // da1Nii4Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (((-1.000000e+01) * (ay)) + ((-7.000000e+00) * (by)));
        case 31: // da1Nii4Njj8
        case 27: // da1Nii4Njj4
            return ((4.000000e+00) / (1.500000e+01)) * (by);
        case 28: // da1Nii4Njj5
            return (((1.000000e+00) / (9.000000e+00)) * (ay)) + (((-7.000000e+00) / (9.000000e+01)) * (by));
        case 32: // da1Nii5Njj1
            return ((1.000000e+00) / (6.000000e+01)) * ((ay) + ((-1.000000e+00) * (by)));
        case 40: // da1Nii6Njj1
        case 42: // da1Nii6Njj3
        case 49: // da1Nii7Njj2
        case 33: // da1Nii5Njj2
            return ((7.000000e+00) / (9.000000e+01)) * (ay);
        case 34: // da1Nii5Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (((8.000000e+00) * (ay)) + ((-3.000000e+00) * (by)));
        case 35: // da1Nii5Njj4
            return (((-1.000000e+00) / (9.000000e+00)) * (ay)) + (((1.300000e+01) / (9.000000e+01)) * (by));
        case 36: // da1Nii5Njj5
            return ((1.000000e+00) / (1.500000e+01)) * (((-1.000000e+00) * (ay)) + (by));
        case 37: // da1Nii5Njj6
            return (((-1.300000e+01) / (9.000000e+01)) * (ay)) + (((1.000000e+00) / (9.000000e+00)) * (by));
        case 38: // da1Nii5Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (((3.000000e+00) * (ay)) + ((-8.000000e+00) * (by)));
        case 45: // da1Nii6Njj6
        case 41: // da1Nii6Njj2
            return ((-4.000000e+00) / (1.500000e+01)) * (ay);
        case 57: // da1Nii8Njj2
        case 43: // da1Nii6Njj4
            return ((-2.000000e+00) / (9.000000e+00)) * ((ay) + (by));
        case 44: // da1Nii6Njj5
            return (((7.000000e+00) / (9.000000e+01)) * (ay)) + (((-1.000000e+00) / (9.000000e+00)) * (by));
        case 46: // da1Nii6Njj7
            return (((7.000000e+00) / (9.000000e+01)) * (ay)) + (((1.000000e+00) / (9.000000e+00)) * (by));
        case 48: // da1Nii7Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (((8.000000e+00) * (ay)) + ((3.000000e+00) * (by)));
        case 50: // da1Nii7Njj3
            return ((1.000000e+00) / (6.000000e+01)) * ((ay) + (by));
        case 52: // da1Nii7Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (((3.000000e+00) * (ay)) + ((8.000000e+00) * (by)));
        case 53: // da1Nii7Njj6
            return ((1.000000e+00) / (9.000000e+01)) * (((-1.300000e+01) * (ay)) + ((-1.000000e+01) * (by)));
        case 54: // da1Nii7Njj7
            return ((1.000000e+00) / (1.500000e+01)) * (((-1.000000e+00) * (ay)) + ((-1.000000e+00) * (by)));
        case 55: // da1Nii7Njj8
            return ((1.000000e+00) / (9.000000e+01)) * (((-1.000000e+01) * (ay)) + ((-1.300000e+01) * (by)));
        case 56: // da1Nii8Njj1
            return (((-1.000000e+00) / (9.000000e+00)) * (ay)) + (((7.000000e+00) / (9.000000e+01)) * (by));
        case 63: // da1Nii8Njj8
        case 59: // da1Nii8Njj4
            return ((-4.000000e+00) / (1.500000e+01)) * (by);
        case 62: // da1Nii8Njj7
            return (((1.000000e+00) / (9.000000e+00)) * (ay)) + (((7.000000e+00) / (9.000000e+01)) * (by));
        case 64: // da2Nii1Njj1
            return ((1.000000e+00) / (1.500000e+01)) * (((-1.000000e+00) * (ax)) + (bx));
        case 65: // da2Nii1Njj2
            return (((-1.300000e+01) / (9.000000e+01)) * (ax)) + (((1.000000e+00) / (9.000000e+00)) * (bx));
        case 66: // da2Nii1Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (((3.000000e+00) * (ax)) + ((-8.000000e+00) * (bx)));
        case 115: // da2Nii7Njj4
        case 122: // da2Nii8Njj3
        case 124: // da2Nii8Njj5
        case 67: // da2Nii1Njj4
            return ((-7.000000e+00) / (9.000000e+01)) * (bx);
        case 68: // da2Nii1Njj5
            return ((1.000000e+00) / (6.000000e+01)) * ((ax) + ((-1.000000e+00) * (bx)));
        case 76: // da2Nii2Njj5
        case 78: // da2Nii2Njj7
        case 85: // da2Nii3Njj6
        case 69: // da2Nii1Njj6
            return ((7.000000e+00) / (9.000000e+01)) * (ax);
        case 70: // da2Nii1Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (((8.000000e+00) * (ax)) + ((-3.000000e+00) * (bx)));
        case 71: // da2Nii1Njj8
            return (((-1.000000e+00) / (9.000000e+00)) * (ax)) + (((1.300000e+01) / (9.000000e+01)) * (bx));
        case 72: // da2Nii2Njj1
            return (((7.000000e+00) / (9.000000e+01)) * (ax)) + (((-1.000000e+00) / (9.000000e+00)) * (bx));
        case 77: // da2Nii2Njj6
        case 73: // da2Nii2Njj2
            return ((-4.000000e+00) / (1.500000e+01)) * (ax);
        case 74: // da2Nii2Njj3
            return (((7.000000e+00) / (9.000000e+01)) * (ax)) + (((1.000000e+00) / (9.000000e+00)) * (bx));
        case 125: // da2Nii8Njj6
        case 75: // da2Nii2Njj4
            return ((-2.000000e+00) / (9.000000e+00)) * ((ax) + ((-1.000000e+00) * (bx)));
        case 93: // da2Nii4Njj6
        case 79: // da2Nii2Njj8
            return ((-2.000000e+00) / (9.000000e+00)) * ((ax) + (bx));
        case 80: // da2Nii3Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (((3.000000e+00) * (ax)) + ((8.000000e+00) * (bx)));
        case 81: // da2Nii3Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (((-1.300000e+01) * (ax)) + ((-1.000000e+01) * (bx)));
        case 82: // da2Nii3Njj3
            return ((1.000000e+00) / (1.500000e+01)) * (((-1.000000e+00) * (ax)) + ((-1.000000e+00) * (bx)));
        case 83: // da2Nii3Njj4
            return ((1.000000e+00) / (9.000000e+01)) * (((-1.000000e+01) * (ax)) + ((-1.300000e+01) * (bx)));
        case 84: // da2Nii3Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (((8.000000e+00) * (ax)) + ((3.000000e+00) * (bx)));
        case 86: // da2Nii3Njj7
            return ((1.000000e+00) / (6.000000e+01)) * ((ax) + (bx));
        case 88: // da2Nii4Njj1
        case 94: // da2Nii4Njj7
        case 103: // da2Nii5Njj8
        case 87: // da2Nii3Njj8
            return ((7.000000e+00) / (9.000000e+01)) * (bx);
        case 111: // da2Nii6Njj8
        case 89: // da2Nii4Njj2
            return ((2.000000e+00) / (9.000000e+00)) * ((ax) + ((-1.000000e+00) * (bx)));
        case 90: // da2Nii4Njj3
            return (((1.000000e+00) / (9.000000e+00)) * (ax)) + (((7.000000e+00) / (9.000000e+01)) * (bx));
        case 95: // da2Nii4Njj8
        case 91: // da2Nii4Njj4
            return ((-4.000000e+00) / (1.500000e+01)) * (bx);
        case 92: // da2Nii4Njj5
            return (((-1.000000e+00) / (9.000000e+00)) * (ax)) + (((7.000000e+00) / (9.000000e+01)) * (bx));
        case 96: // da2Nii5Njj1
            return ((1.000000e+00) / (6.000000e+01)) * (((-1.000000e+00) * (ax)) + (bx));
        case 104: // da2Nii6Njj1
        case 106: // da2Nii6Njj3
        case 113: // da2Nii7Njj2
        case 97: // da2Nii5Njj2
            return ((-7.000000e+00) / (9.000000e+01)) * (ax);
        case 98: // da2Nii5Njj3
            return ((1.000000e+00) / (1.800000e+02)) * (((-8.000000e+00) * (ax)) + ((3.000000e+00) * (bx)));
        case 99: // da2Nii5Njj4
            return (((1.000000e+00) / (9.000000e+00)) * (ax)) + (((-1.300000e+01) / (9.000000e+01)) * (bx));
        case 100: // da2Nii5Njj5
            return ((1.000000e+00) / (1.500000e+01)) * ((ax) + ((-1.000000e+00) * (bx)));
        case 101: // da2Nii5Njj6
            return (((1.300000e+01) / (9.000000e+01)) * (ax)) + (((-1.000000e+00) / (9.000000e+00)) * (bx));
        case 102: // da2Nii5Njj7
            return ((1.000000e+00) / (1.800000e+02)) * (((-3.000000e+00) * (ax)) + ((8.000000e+00) * (bx)));
        case 109: // da2Nii6Njj6
        case 105: // da2Nii6Njj2
            return ((4.000000e+00) / (1.500000e+01)) * (ax);
        case 121: // da2Nii8Njj2
        case 107: // da2Nii6Njj4
            return ((2.000000e+00) / (9.000000e+00)) * ((ax) + (bx));
        case 108: // da2Nii6Njj5
            return (((-7.000000e+00) / (9.000000e+01)) * (ax)) + (((1.000000e+00) / (9.000000e+00)) * (bx));
        case 110: // da2Nii6Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (((-7.000000e+00) * (ax)) + ((-1.000000e+01) * (bx)));
        case 112: // da2Nii7Njj1
            return ((1.000000e+00) / (1.800000e+02)) * (((-8.000000e+00) * (ax)) + ((-3.000000e+00) * (bx)));
        case 114: // da2Nii7Njj3
            return ((1.000000e+00) / (6.000000e+01)) * (((-1.000000e+00) * (ax)) + ((-1.000000e+00) * (bx)));
        case 116: // da2Nii7Njj5
            return ((1.000000e+00) / (1.800000e+02)) * (((-3.000000e+00) * (ax)) + ((-8.000000e+00) * (bx)));
        case 117: // da2Nii7Njj6
            return (((1.300000e+01) / (9.000000e+01)) * (ax)) + (((1.000000e+00) / (9.000000e+00)) * (bx));
        case 118: // da2Nii7Njj7
            return ((1.000000e+00) / (1.500000e+01)) * ((ax) + (bx));
        case 119: // da2Nii7Njj8
            return (((1.000000e+00) / (9.000000e+00)) * (ax)) + (((1.300000e+01) / (9.000000e+01)) * (bx));
        case 120: // da2Nii8Njj1
            return (((1.000000e+00) / (9.000000e+00)) * (ax)) + (((-7.000000e+00) / (9.000000e+01)) * (bx));
        case 127: // da2Nii8Njj8
        case 123: // da2Nii8Njj4
            return ((4.000000e+00) / (1.500000e+01)) * (bx);
        case 126: // da2Nii8Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (((-1.000000e+01) * (ax)) + ((-7.000000e+00) * (bx)));
		
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");
	}	
	
}

double Element2DRect2nd::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 8 * elem_shape_func_1;
	                       
	
	// ----------------------------------------------
	// return analytical integral
	// ----------------------------------------------
	switch(obj_index) {
        case 18: // Nii3Njj3
        case 36: // Nii5Njj5
        case 54: // Nii7Njj7
        case 0: // Nii1Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (element_volume);
        case 7: // Nii1Njj8
        case 8: // Nii2Njj1
        case 10: // Nii2Njj3
        case 17: // Nii3Njj2
        case 19: // Nii3Njj4
        case 26: // Nii4Njj3
        case 28: // Nii4Njj5
        case 35: // Nii5Njj4
        case 37: // Nii5Njj6
        case 44: // Nii6Njj5
        case 46: // Nii6Njj7
        case 53: // Nii7Njj6
        case 55: // Nii7Njj8
        case 56: // Nii8Njj1
        case 62: // Nii8Njj7
        case 1: // Nii1Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (element_volume);
        case 6: // Nii1Njj7
        case 16: // Nii3Njj1
        case 20: // Nii3Njj5
        case 34: // Nii5Njj3
        case 38: // Nii5Njj7
        case 48: // Nii7Njj1
        case 52: // Nii7Njj5
        case 2: // Nii1Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (element_volume);
        case 5: // Nii1Njj6
        case 12: // Nii2Njj5
        case 14: // Nii2Njj7
        case 21: // Nii3Njj6
        case 23: // Nii3Njj8
        case 24: // Nii4Njj1
        case 30: // Nii4Njj7
        case 33: // Nii5Njj2
        case 39: // Nii5Njj8
        case 40: // Nii6Njj1
        case 42: // Nii6Njj3
        case 49: // Nii7Njj2
        case 51: // Nii7Njj4
        case 58: // Nii8Njj3
        case 60: // Nii8Njj5
        case 3: // Nii1Njj4
            return ((-2.000000e+00) / (4.500000e+01)) * (element_volume);
        case 22: // Nii3Njj7
        case 32: // Nii5Njj1
        case 50: // Nii7Njj3
        case 4: // Nii1Njj5
            return ((1.000000e+00) / (6.000000e+01)) * (element_volume);
        case 27: // Nii4Njj4
        case 45: // Nii6Njj6
        case 63: // Nii8Njj8
        case 9: // Nii2Njj2
            return ((8.000000e+00) / (4.500000e+01)) * (element_volume);
        case 15: // Nii2Njj8
        case 25: // Nii4Njj2
        case 29: // Nii4Njj6
        case 43: // Nii6Njj4
        case 47: // Nii6Njj8
        case 57: // Nii8Njj2
        case 61: // Nii8Njj6
        case 11: // Nii2Njj4
            return ((1.000000e+00) / (9.000000e+00)) * (element_volume);
        case 31: // Nii4Njj8
        case 41: // Nii6Njj2
        case 59: // Nii8Njj4
        case 13: // Nii2Njj6
            return ((4.000000e+00) / (4.500000e+01)) * (element_volume);

		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");
	}		
}


double Element2DRect2nd::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	
	unsigned int obj_index = elem_shape_func_2 
	                       + 8 * elem_shape_func_1
	                       + 64 * nodal_data_point;

	switch(obj_index) {
        case 146: // Ncc3Nii3Njj3
        case 292: // Ncc5Nii5Njj5
        case 438: // Ncc7Nii7Njj7
        case 0: // Ncc1Nii1Njj1
            return ((1.100000e+01) / (7.000000e+02)) * jacobi_det;
        case 7: // Ncc1Nii1Njj8
        case 8: // Ncc1Nii2Njj1
        case 56: // Ncc1Nii8Njj1
        case 64: // Ncc2Nii1Njj1
        case 82: // Ncc2Nii3Njj3
        case 138: // Ncc3Nii2Njj3
        case 145: // Ncc3Nii3Njj2
        case 147: // Ncc3Nii3Njj4
        case 154: // Ncc3Nii4Njj3
        case 210: // Ncc4Nii3Njj3
        case 228: // Ncc4Nii5Njj5
        case 284: // Ncc5Nii4Njj5
        case 291: // Ncc5Nii5Njj4
        case 293: // Ncc5Nii5Njj6
        case 300: // Ncc5Nii6Njj5
        case 356: // Ncc6Nii5Njj5
        case 374: // Ncc6Nii7Njj7
        case 430: // Ncc7Nii6Njj7
        case 437: // Ncc7Nii7Njj6
        case 439: // Ncc7Nii7Njj8
        case 446: // Ncc7Nii8Njj7
        case 448: // Ncc8Nii1Njj1
        case 502: // Ncc8Nii7Njj7
        case 1: // Ncc1Nii1Njj2
            return ((1.000000e+00) / (2.100000e+01)) * jacobi_det;
        case 6: // Ncc1Nii1Njj7
        case 16: // Ncc1Nii3Njj1
        case 18: // Ncc1Nii3Njj3
        case 48: // Ncc1Nii7Njj1
        case 54: // Ncc1Nii7Njj7
        case 128: // Ncc3Nii1Njj1
        case 130: // Ncc3Nii1Njj3
        case 144: // Ncc3Nii3Njj1
        case 148: // Ncc3Nii3Njj5
        case 162: // Ncc3Nii5Njj3
        case 164: // Ncc3Nii5Njj5
        case 274: // Ncc5Nii3Njj3
        case 276: // Ncc5Nii3Njj5
        case 290: // Ncc5Nii5Njj3
        case 294: // Ncc5Nii5Njj7
        case 308: // Ncc5Nii7Njj5
        case 310: // Ncc5Nii7Njj7
        case 384: // Ncc7Nii1Njj1
        case 390: // Ncc7Nii1Njj7
        case 420: // Ncc7Nii5Njj5
        case 422: // Ncc7Nii5Njj7
        case 432: // Ncc7Nii7Njj1
        case 436: // Ncc7Nii7Njj5
        case 2: // Ncc1Nii1Njj3
            return ((-1.900000e+01) / (1.260000e+03)) * jacobi_det;
        case 5: // Ncc1Nii1Njj6
        case 24: // Ncc1Nii4Njj1
        case 40: // Ncc1Nii6Njj1
        case 100: // Ncc2Nii5Njj5
        case 118: // Ncc2Nii7Njj7
        case 149: // Ncc3Nii3Njj6
        case 151: // Ncc3Nii3Njj8
        case 170: // Ncc3Nii6Njj3
        case 186: // Ncc3Nii8Njj3
        case 192: // Ncc4Nii1Njj1
        case 246: // Ncc4Nii7Njj7
        case 268: // Ncc5Nii2Njj5
        case 289: // Ncc5Nii5Njj2
        case 295: // Ncc5Nii5Njj8
        case 316: // Ncc5Nii8Njj5
        case 320: // Ncc6Nii1Njj1
        case 338: // Ncc6Nii3Njj3
        case 398: // Ncc7Nii2Njj7
        case 414: // Ncc7Nii4Njj7
        case 433: // Ncc7Nii7Njj2
        case 435: // Ncc7Nii7Njj4
        case 466: // Ncc8Nii3Njj3
        case 484: // Ncc8Nii5Njj5
        case 3: // Ncc1Nii1Njj4
            return ((5.300000e+01) / (1.575000e+03)) * jacobi_det;
        case 32: // Ncc1Nii5Njj1
        case 36: // Ncc1Nii5Njj5
        case 150: // Ncc3Nii3Njj7
        case 178: // Ncc3Nii7Njj3
        case 182: // Ncc3Nii7Njj7
        case 256: // Ncc5Nii1Njj1
        case 260: // Ncc5Nii1Njj5
        case 288: // Ncc5Nii5Njj1
        case 402: // Ncc7Nii3Njj3
        case 406: // Ncc7Nii3Njj7
        case 434: // Ncc7Nii7Njj3
        case 4: // Ncc1Nii1Njj5
            return ((-3.100000e+01) / (2.100000e+03)) * jacobi_det;
        case 63: // Ncc1Nii8Njj8
        case 65: // Ncc2Nii1Njj2
        case 72: // Ncc2Nii2Njj1
        case 74: // Ncc2Nii2Njj3
        case 81: // Ncc2Nii3Njj2
        case 137: // Ncc3Nii2Njj2
        case 155: // Ncc3Nii4Njj4
        case 211: // Ncc4Nii3Njj4
        case 218: // Ncc4Nii4Njj3
        case 220: // Ncc4Nii4Njj5
        case 227: // Ncc4Nii5Njj4
        case 283: // Ncc5Nii4Njj4
        case 301: // Ncc5Nii6Njj6
        case 357: // Ncc6Nii5Njj6
        case 364: // Ncc6Nii6Njj5
        case 366: // Ncc6Nii6Njj7
        case 373: // Ncc6Nii7Njj6
        case 429: // Ncc7Nii6Njj6
        case 447: // Ncc7Nii8Njj8
        case 455: // Ncc8Nii1Njj8
        case 503: // Ncc8Nii7Njj8
        case 504: // Ncc8Nii8Njj1
        case 510: // Ncc8Nii8Njj7
        case 9: // Ncc1Nii2Njj2
            return ((-1.200000e+01) / (1.750000e+02)) * jacobi_det;
        case 17: // Ncc1Nii3Njj2
        case 55: // Ncc1Nii7Njj8
        case 62: // Ncc1Nii8Njj7
        case 66: // Ncc2Nii1Njj3
        case 80: // Ncc2Nii3Njj1
        case 129: // Ncc3Nii1Njj2
        case 136: // Ncc3Nii2Njj1
        case 156: // Ncc3Nii4Njj5
        case 163: // Ncc3Nii5Njj4
        case 212: // Ncc4Nii3Njj5
        case 226: // Ncc4Nii5Njj3
        case 275: // Ncc5Nii3Njj4
        case 282: // Ncc5Nii4Njj3
        case 302: // Ncc5Nii6Njj7
        case 309: // Ncc5Nii7Njj6
        case 358: // Ncc6Nii5Njj7
        case 372: // Ncc6Nii7Njj5
        case 391: // Ncc7Nii1Njj8
        case 421: // Ncc7Nii5Njj6
        case 428: // Ncc7Nii6Njj5
        case 440: // Ncc7Nii8Njj1
        case 454: // Ncc8Nii1Njj7
        case 496: // Ncc8Nii7Njj1
        case 10: // Ncc1Nii2Njj3
            return ((2.600000e+01) / (1.575000e+03)) * jacobi_det;
        case 25: // Ncc1Nii4Njj2
        case 29: // Ncc1Nii4Njj6
        case 43: // Ncc1Nii6Njj4
        case 47: // Ncc1Nii6Njj8
        case 61: // Ncc1Nii8Njj6
        case 67: // Ncc2Nii1Njj4
        case 87: // Ncc2Nii3Njj8
        case 88: // Ncc2Nii4Njj1
        case 92: // Ncc2Nii4Njj5
        case 94: // Ncc2Nii4Njj7
        case 99: // Ncc2Nii5Njj4
        case 103: // Ncc2Nii5Njj8
        case 115: // Ncc2Nii7Njj4
        case 119: // Ncc2Nii7Njj8
        case 122: // Ncc2Nii8Njj3
        case 124: // Ncc2Nii8Njj5
        case 126: // Ncc2Nii8Njj7
        case 143: // Ncc3Nii2Njj8
        case 157: // Ncc3Nii4Njj6
        case 171: // Ncc3Nii6Njj4
        case 175: // Ncc3Nii6Njj8
        case 185: // Ncc3Nii8Njj2
        case 189: // Ncc3Nii8Njj6
        case 193: // Ncc4Nii1Njj2
        case 197: // Ncc4Nii1Njj6
        case 200: // Ncc4Nii2Njj1
        case 204: // Ncc4Nii2Njj5
        case 206: // Ncc4Nii2Njj7
        case 213: // Ncc4Nii3Njj6
        case 225: // Ncc4Nii5Njj2
        case 232: // Ncc4Nii6Njj1
        case 234: // Ncc4Nii6Njj3
        case 238: // Ncc4Nii6Njj7
        case 241: // Ncc4Nii7Njj2
        case 245: // Ncc4Nii7Njj6
        case 267: // Ncc5Nii2Njj4
        case 271: // Ncc5Nii2Njj8
        case 281: // Ncc5Nii4Njj2
        case 303: // Ncc5Nii6Njj8
        case 313: // Ncc5Nii8Njj2
        case 317: // Ncc5Nii8Njj6
        case 323: // Ncc6Nii1Njj4
        case 327: // Ncc6Nii1Njj8
        case 339: // Ncc6Nii3Njj4
        case 343: // Ncc6Nii3Njj8
        case 344: // Ncc6Nii4Njj1
        case 346: // Ncc6Nii4Njj3
        case 350: // Ncc6Nii4Njj7
        case 359: // Ncc6Nii5Njj8
        case 371: // Ncc6Nii7Njj4
        case 376: // Ncc6Nii8Njj1
        case 378: // Ncc6Nii8Njj3
        case 380: // Ncc6Nii8Njj5
        case 395: // Ncc7Nii2Njj4
        case 399: // Ncc7Nii2Njj8
        case 409: // Ncc7Nii4Njj2
        case 413: // Ncc7Nii4Njj6
        case 427: // Ncc7Nii6Njj4
        case 441: // Ncc7Nii8Njj2
        case 453: // Ncc8Nii1Njj6
        case 458: // Ncc8Nii2Njj3
        case 460: // Ncc8Nii2Njj5
        case 462: // Ncc8Nii2Njj7
        case 465: // Ncc8Nii3Njj2
        case 469: // Ncc8Nii3Njj6
        case 481: // Ncc8Nii5Njj2
        case 485: // Ncc8Nii5Njj6
        case 488: // Ncc8Nii6Njj1
        case 490: // Ncc8Nii6Njj3
        case 492: // Ncc8Nii6Njj5
        case 497: // Ncc8Nii7Njj2
        case 11: // Ncc1Nii2Njj4
            return ((-1.600000e+01) / (2.250000e+02)) * jacobi_det;
        case 21: // Ncc1Nii3Njj6
        case 28: // Ncc1Nii4Njj5
        case 30: // Ncc1Nii4Njj7
        case 33: // Ncc1Nii5Njj2
        case 35: // Ncc1Nii5Njj4
        case 37: // Ncc1Nii5Njj6
        case 39: // Ncc1Nii5Njj8
        case 42: // Ncc1Nii6Njj3
        case 44: // Ncc1Nii6Njj5
        case 51: // Ncc1Nii7Njj4
        case 60: // Ncc1Nii8Njj5
        case 68: // Ncc2Nii1Njj5
        case 86: // Ncc2Nii3Njj7
        case 96: // Ncc2Nii5Njj1
        case 102: // Ncc2Nii5Njj7
        case 114: // Ncc2Nii7Njj3
        case 116: // Ncc2Nii7Njj5
        case 133: // Ncc3Nii1Njj6
        case 142: // Ncc3Nii2Njj7
        case 158: // Ncc3Nii4Njj7
        case 167: // Ncc3Nii5Njj8
        case 168: // Ncc3Nii6Njj1
        case 174: // Ncc3Nii6Njj7
        case 177: // Ncc3Nii7Njj2
        case 179: // Ncc3Nii7Njj4
        case 181: // Ncc3Nii7Njj6
        case 183: // Ncc3Nii7Njj8
        case 188: // Ncc3Nii8Njj5
        case 190: // Ncc3Nii8Njj7
        case 196: // Ncc4Nii1Njj5
        case 198: // Ncc4Nii1Njj7
        case 214: // Ncc4Nii3Njj7
        case 224: // Ncc4Nii5Njj1
        case 240: // Ncc4Nii7Njj1
        case 242: // Ncc4Nii7Njj3
        case 257: // Ncc5Nii1Njj2
        case 259: // Ncc5Nii1Njj4
        case 261: // Ncc5Nii1Njj6
        case 263: // Ncc5Nii1Njj8
        case 264: // Ncc5Nii2Njj1
        case 270: // Ncc5Nii2Njj7
        case 279: // Ncc5Nii3Njj8
        case 280: // Ncc5Nii4Njj1
        case 296: // Ncc5Nii6Njj1
        case 305: // Ncc5Nii7Njj2
        case 312: // Ncc5Nii8Njj1
        case 314: // Ncc5Nii8Njj3
        case 322: // Ncc6Nii1Njj3
        case 324: // Ncc6Nii1Njj5
        case 336: // Ncc6Nii3Njj1
        case 342: // Ncc6Nii3Njj7
        case 352: // Ncc6Nii5Njj1
        case 370: // Ncc6Nii7Njj3
        case 387: // Ncc7Nii1Njj4
        case 394: // Ncc7Nii2Njj3
        case 396: // Ncc7Nii2Njj5
        case 401: // Ncc7Nii3Njj2
        case 403: // Ncc7Nii3Njj4
        case 405: // Ncc7Nii3Njj6
        case 407: // Ncc7Nii3Njj8
        case 408: // Ncc7Nii4Njj1
        case 410: // Ncc7Nii4Njj3
        case 417: // Ncc7Nii5Njj2
        case 426: // Ncc7Nii6Njj3
        case 442: // Ncc7Nii8Njj3
        case 452: // Ncc8Nii1Njj5
        case 468: // Ncc8Nii3Njj5
        case 470: // Ncc8Nii3Njj7
        case 480: // Ncc8Nii5Njj1
        case 482: // Ncc8Nii5Njj3
        case 498: // Ncc8Nii7Njj3
        case 12: // Ncc1Nii2Njj5
            return ((4.600000e+01) / (1.575000e+03)) * jacobi_det;
        case 31: // Ncc1Nii4Njj8
        case 41: // Ncc1Nii6Njj2
        case 59: // Ncc1Nii8Njj4
        case 69: // Ncc2Nii1Njj6
        case 85: // Ncc2Nii3Njj6
        case 101: // Ncc2Nii5Njj6
        case 104: // Ncc2Nii6Njj1
        case 106: // Ncc2Nii6Njj3
        case 108: // Ncc2Nii6Njj5
        case 110: // Ncc2Nii6Njj7
        case 117: // Ncc2Nii7Njj6
        case 141: // Ncc3Nii2Njj6
        case 159: // Ncc3Nii4Njj8
        case 169: // Ncc3Nii6Njj2
        case 187: // Ncc3Nii8Njj4
        case 199: // Ncc4Nii1Njj8
        case 215: // Ncc4Nii3Njj8
        case 231: // Ncc4Nii5Njj8
        case 247: // Ncc4Nii7Njj8
        case 248: // Ncc4Nii8Njj1
        case 250: // Ncc4Nii8Njj3
        case 252: // Ncc4Nii8Njj5
        case 254: // Ncc4Nii8Njj7
        case 269: // Ncc5Nii2Njj6
        case 287: // Ncc5Nii4Njj8
        case 297: // Ncc5Nii6Njj2
        case 315: // Ncc5Nii8Njj4
        case 321: // Ncc6Nii1Njj2
        case 328: // Ncc6Nii2Njj1
        case 330: // Ncc6Nii2Njj3
        case 332: // Ncc6Nii2Njj5
        case 334: // Ncc6Nii2Njj7
        case 337: // Ncc6Nii3Njj2
        case 353: // Ncc6Nii5Njj2
        case 369: // Ncc6Nii7Njj2
        case 397: // Ncc7Nii2Njj6
        case 415: // Ncc7Nii4Njj8
        case 425: // Ncc7Nii6Njj2
        case 443: // Ncc7Nii8Njj4
        case 451: // Ncc8Nii1Njj4
        case 467: // Ncc8Nii3Njj4
        case 472: // Ncc8Nii4Njj1
        case 474: // Ncc8Nii4Njj3
        case 476: // Ncc8Nii4Njj5
        case 478: // Ncc8Nii4Njj7
        case 483: // Ncc8Nii5Njj4
        case 499: // Ncc8Nii7Njj4
        case 13: // Ncc1Nii2Njj6
            return ((-9.200000e+01) / (1.575000e+03)) * jacobi_det;
        case 19: // Ncc1Nii3Njj4
        case 23: // Ncc1Nii3Njj8
        case 26: // Ncc1Nii4Njj3
        case 46: // Ncc1Nii6Njj7
        case 49: // Ncc1Nii7Njj2
        case 53: // Ncc1Nii7Njj6
        case 58: // Ncc1Nii8Njj3
        case 70: // Ncc2Nii1Njj7
        case 84: // Ncc2Nii3Njj5
        case 98: // Ncc2Nii5Njj3
        case 112: // Ncc2Nii7Njj1
        case 131: // Ncc3Nii1Njj4
        case 135: // Ncc3Nii1Njj8
        case 140: // Ncc3Nii2Njj5
        case 152: // Ncc3Nii4Njj1
        case 161: // Ncc3Nii5Njj2
        case 165: // Ncc3Nii5Njj6
        case 172: // Ncc3Nii6Njj5
        case 184: // Ncc3Nii8Njj1
        case 194: // Ncc4Nii1Njj3
        case 208: // Ncc4Nii3Njj1
        case 230: // Ncc4Nii5Njj7
        case 244: // Ncc4Nii7Njj5
        case 266: // Ncc5Nii2Njj3
        case 273: // Ncc5Nii3Njj2
        case 277: // Ncc5Nii3Njj6
        case 286: // Ncc5Nii4Njj7
        case 298: // Ncc5Nii6Njj3
        case 307: // Ncc5Nii7Njj4
        case 311: // Ncc5Nii7Njj8
        case 318: // Ncc5Nii8Njj7
        case 326: // Ncc6Nii1Njj7
        case 340: // Ncc6Nii3Njj5
        case 354: // Ncc6Nii5Njj3
        case 368: // Ncc6Nii7Njj1
        case 385: // Ncc7Nii1Njj2
        case 389: // Ncc7Nii1Njj6
        case 392: // Ncc7Nii2Njj1
        case 412: // Ncc7Nii4Njj5
        case 419: // Ncc7Nii5Njj4
        case 423: // Ncc7Nii5Njj8
        case 424: // Ncc7Nii6Njj1
        case 444: // Ncc7Nii8Njj5
        case 450: // Ncc8Nii1Njj3
        case 464: // Ncc8Nii3Njj1
        case 486: // Ncc8Nii5Njj7
        case 500: // Ncc8Nii7Njj5
        case 14: // Ncc1Nii2Njj7
            return ((1.300000e+01) / (5.250000e+02)) * jacobi_det;
        case 57: // Ncc1Nii8Njj2
        case 71: // Ncc2Nii1Njj8
        case 83: // Ncc2Nii3Njj4
        case 90: // Ncc2Nii4Njj3
        case 120: // Ncc2Nii8Njj1
        case 139: // Ncc3Nii2Njj4
        case 153: // Ncc3Nii4Njj2
        case 202: // Ncc4Nii2Njj3
        case 209: // Ncc4Nii3Njj2
        case 229: // Ncc4Nii5Njj6
        case 236: // Ncc4Nii6Njj5
        case 285: // Ncc5Nii4Njj6
        case 299: // Ncc5Nii6Njj4
        case 348: // Ncc6Nii4Njj5
        case 355: // Ncc6Nii5Njj4
        case 375: // Ncc6Nii7Njj8
        case 382: // Ncc6Nii8Njj7
        case 431: // Ncc7Nii6Njj8
        case 445: // Ncc7Nii8Njj6
        case 449: // Ncc8Nii1Njj2
        case 456: // Ncc8Nii2Njj1
        case 494: // Ncc8Nii6Njj7
        case 501: // Ncc8Nii7Njj6
        case 15: // Ncc1Nii2Njj8
            return ((-4.000000e+00) / (7.500000e+01)) * jacobi_det;
        case 22: // Ncc1Nii3Njj7
        case 34: // Ncc1Nii5Njj3
        case 38: // Ncc1Nii5Njj7
        case 50: // Ncc1Nii7Njj3
        case 52: // Ncc1Nii7Njj5
        case 132: // Ncc3Nii1Njj5
        case 134: // Ncc3Nii1Njj7
        case 160: // Ncc3Nii5Njj1
        case 166: // Ncc3Nii5Njj7
        case 176: // Ncc3Nii7Njj1
        case 180: // Ncc3Nii7Njj5
        case 258: // Ncc5Nii1Njj3
        case 262: // Ncc5Nii1Njj7
        case 272: // Ncc5Nii3Njj1
        case 278: // Ncc5Nii3Njj7
        case 304: // Ncc5Nii7Njj1
        case 306: // Ncc5Nii7Njj3
        case 386: // Ncc7Nii1Njj3
        case 388: // Ncc7Nii1Njj5
        case 400: // Ncc7Nii3Njj1
        case 404: // Ncc7Nii3Njj5
        case 416: // Ncc7Nii5Njj1
        case 418: // Ncc7Nii5Njj3
        case 20: // Ncc1Nii3Njj5
            return ((-1.300000e+01) / (1.260000e+03)) * jacobi_det;
        case 45: // Ncc1Nii6Njj6
        case 76: // Ncc2Nii2Njj5
        case 78: // Ncc2Nii2Njj7
        case 97: // Ncc2Nii5Njj2
        case 113: // Ncc2Nii7Njj2
        case 173: // Ncc3Nii6Njj6
        case 191: // Ncc3Nii8Njj8
        case 195: // Ncc4Nii1Njj4
        case 216: // Ncc4Nii4Njj1
        case 222: // Ncc4Nii4Njj7
        case 243: // Ncc4Nii7Njj4
        case 265: // Ncc5Nii2Njj2
        case 319: // Ncc5Nii8Njj8
        case 325: // Ncc6Nii1Njj6
        case 341: // Ncc6Nii3Njj6
        case 360: // Ncc6Nii6Njj1
        case 362: // Ncc6Nii6Njj3
        case 393: // Ncc7Nii2Njj2
        case 411: // Ncc7Nii4Njj4
        case 471: // Ncc8Nii3Njj8
        case 487: // Ncc8Nii5Njj8
        case 506: // Ncc8Nii8Njj3
        case 508: // Ncc8Nii8Njj5
        case 27: // Ncc1Nii4Njj4
            return ((-1.480000e+02) / (1.575000e+03)) * jacobi_det;
        case 219: // Ncc4Nii4Njj4
        case 365: // Ncc6Nii6Njj6
        case 511: // Ncc8Nii8Njj8
        case 73: // Ncc2Nii2Njj2
            return ((1.600000e+01) / (3.500000e+01)) * jacobi_det;
        case 79: // Ncc2Nii2Njj8
        case 89: // Ncc2Nii4Njj2
        case 91: // Ncc2Nii4Njj4
        case 121: // Ncc2Nii8Njj2
        case 127: // Ncc2Nii8Njj8
        case 201: // Ncc4Nii2Njj2
        case 203: // Ncc4Nii2Njj4
        case 217: // Ncc4Nii4Njj2
        case 221: // Ncc4Nii4Njj6
        case 235: // Ncc4Nii6Njj4
        case 237: // Ncc4Nii6Njj6
        case 347: // Ncc6Nii4Njj4
        case 349: // Ncc6Nii4Njj6
        case 363: // Ncc6Nii6Njj4
        case 367: // Ncc6Nii6Njj8
        case 381: // Ncc6Nii8Njj6
        case 383: // Ncc6Nii8Njj8
        case 457: // Ncc8Nii2Njj2
        case 463: // Ncc8Nii2Njj8
        case 493: // Ncc8Nii6Njj6
        case 495: // Ncc8Nii6Njj8
        case 505: // Ncc8Nii8Njj2
        case 509: // Ncc8Nii8Njj6
        case 75: // Ncc2Nii2Njj4
            return ((1.600000e+01) / (7.500000e+01)) * jacobi_det;
        case 105: // Ncc2Nii6Njj2
        case 109: // Ncc2Nii6Njj6
        case 223: // Ncc4Nii4Njj8
        case 251: // Ncc4Nii8Njj4
        case 255: // Ncc4Nii8Njj8
        case 329: // Ncc6Nii2Njj2
        case 333: // Ncc6Nii2Njj6
        case 361: // Ncc6Nii6Njj2
        case 475: // Ncc8Nii4Njj4
        case 479: // Ncc8Nii4Njj8
        case 507: // Ncc8Nii8Njj4
        case 77: // Ncc2Nii2Njj6
            return ((1.600000e+01) / (1.050000e+02)) * jacobi_det;
        case 95: // Ncc2Nii4Njj8
        case 107: // Ncc2Nii6Njj4
        case 111: // Ncc2Nii6Njj8
        case 123: // Ncc2Nii8Njj4
        case 125: // Ncc2Nii8Njj6
        case 205: // Ncc4Nii2Njj6
        case 207: // Ncc4Nii2Njj8
        case 233: // Ncc4Nii6Njj2
        case 239: // Ncc4Nii6Njj8
        case 249: // Ncc4Nii8Njj2
        case 253: // Ncc4Nii8Njj6
        case 331: // Ncc6Nii2Njj4
        case 335: // Ncc6Nii2Njj8
        case 345: // Ncc6Nii4Njj2
        case 351: // Ncc6Nii4Njj8
        case 377: // Ncc6Nii8Njj2
        case 379: // Ncc6Nii8Njj4
        case 459: // Ncc8Nii2Njj4
        case 461: // Ncc8Nii2Njj6
        case 473: // Ncc8Nii4Njj2
        case 477: // Ncc8Nii4Njj6
        case 489: // Ncc8Nii6Njj2
        case 491: // Ncc8Nii6Njj4
        case 93: // Ncc2Nii4Njj6
            return ((3.200000e+01) / (2.250000e+02)) * jacobi_det;
		default:
			TDKP_GENERAL_EXCEPTION("wrong element / nodal point index");		
	};

}

double Element2DRect2nd::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
	
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_1 
	                       + 8 * diffop;
	                       	
	// ----------------------------------------------
	// return analytical integral
	// ----------------------------------------------
	switch(obj_index) {
        case 0: // da1Nii1
            return ((1.000000e+00) / (6.000000e+00)) * ((ay) + ((-1.000000e+00) * (by)));
        case 1: // da1Nii2
            return ((2.000000e+00) / (3.000000e+00)) * (ay);
        case 2: // da1Nii3
            return ((1.000000e+00) / (6.000000e+00)) * ((ay) + (by));
        case 3: // da1Nii4
            return ((2.000000e+00) / (3.000000e+00)) * (by);
        case 4: // da1Nii5
            return ((1.000000e+00) / (6.000000e+00)) * (((-1.000000e+00) * (ay)) + (by));
        case 5: // da1Nii6
            return ((-2.000000e+00) / (3.000000e+00)) * (ay);
        case 6: // da1Nii7
            return ((1.000000e+00) / (6.000000e+00)) * (((-1.000000e+00) * (ay)) + ((-1.000000e+00) * (by)));
        case 7: // da1Nii8
            return ((-2.000000e+00) / (3.000000e+00)) * (by);
        case 8: // da2Nii1
            return ((1.000000e+00) / (6.000000e+00)) * (((-1.000000e+00) * (ax)) + (bx));
        case 9: // da2Nii2
            return ((-2.000000e+00) / (3.000000e+00)) * (ax);
        case 10: // da2Nii3
            return ((1.000000e+00) / (6.000000e+00)) * (((-1.000000e+00) * (ax)) + ((-1.000000e+00) * (bx)));
        case 11: // da2Nii4
            return ((-2.000000e+00) / (3.000000e+00)) * (bx);
        case 12: // da2Nii5
            return ((1.000000e+00) / (6.000000e+00)) * ((ax) + ((-1.000000e+00) * (bx)));
        case 13: // da2Nii6
            return ((2.000000e+00) / (3.000000e+00)) * (ax);
        case 14: // da2Nii7
            return ((1.000000e+00) / (6.000000e+00)) * ((ax) + (bx));
        case 15: // da2Nii8
            return ((2.000000e+00) / (3.000000e+00)) * (bx);		
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");
	}	
}

double Element2DRect2nd::get_single_integral_0th_order(short elem_shape_func_1) const {

	// ----------------------------------------------
	// return analytical integral
	// ----------------------------------------------
	switch(elem_shape_func_1) {
        case 2: // Nii3
        case 4: // Nii5
        case 6: // Nii7
        case 0: // Nii1
            return ((-1.000000e+00) / (1.200000e+01)) * (element_volume);
        case 3: // Nii4
        case 5: // Nii6
        case 7: // Nii8
        case 1: // Nii2
            return ((1.000000e+00) / (3.000000e+00)) * (element_volume);		
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");	
	}
	
}						

		
double Element2DRect2nd::evaluate_form_function_derivative(short diffop, short elem_shape_func, const double& x, const double& y) const {	
	
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func 
	                       + 8 * diffop;
	                       	
	// ----------------------------------------------
	// return analytical integral
	// ----------------------------------------------
	switch(obj_index) {
        case 0: // da1Nii1
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((-1.000000e+00) * (by) * ((-1.000000e+00) + (y)) * (((2.000000e+00) * (x)) + (y))) + ((ay) * ((-1.000000e+00) + (x)) * ((x) + ((2.000000e+00) * (y)))));
        case 1: // da1Nii2
            return (1.0 / (element_volume)) * ((ay) + ((-1.000000e+00) * (ay) * ((x) * (x))) + ((2.000000e+00) * (by) * (x) * ((-1.000000e+00) + (y))));
        case 2: // da1Nii3
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((ay) * ((1.000000e+00) + (x)) * ((x) + ((-2.000000e+00) * (y)))) + ((-1.000000e+00) * (by) * (((2.000000e+00) * (x)) + ((-1.000000e+00) * (y))) * ((-1.000000e+00) + (y))));
        case 3: // da1Nii4
            return (1.0 / (element_volume)) * ((by) + ((2.000000e+00) * (ay) * ((1.000000e+00) + (x)) * (y)) + ((-1.000000e+00) * (by) * ((y) * (y))));
        case 4: // da1Nii5
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((by) * ((1.000000e+00) + (y)) * (((2.000000e+00) * (x)) + (y))) + ((-1.000000e+00) * (ay) * ((1.000000e+00) + (x)) * ((x) + ((2.000000e+00) * (y)))));
        case 5: // da1Nii6
            return (1.0 / (element_volume)) * (((ay) * ((-1.000000e+00) + ((x) * (x)))) + ((-2.000000e+00) * (by) * (x) * ((1.000000e+00) + (y))));
        case 6: // da1Nii7
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((-1.000000e+00) * (ay) * ((-1.000000e+00) + (x)) * ((x) + ((-2.000000e+00) * (y)))) + ((by) * (((2.000000e+00) * (x)) + ((-1.000000e+00) * (y))) * ((1.000000e+00) + (y))));
        case 7: // da1Nii8
            return (1.0 / (element_volume)) * (((-2.000000e+00) * (ay) * ((-1.000000e+00) + (x)) * (y)) + ((by) * ((-1.000000e+00) + ((y) * (y)))));
        case 8: // da2Nii1
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((bx) * ((-1.000000e+00) + (y)) * (((2.000000e+00) * (x)) + (y))) + ((-1.000000e+00) * (ax) * ((-1.000000e+00) + (x)) * ((x) + ((2.000000e+00) * (y)))));
        case 9: // da2Nii2
            return (1.0 / (element_volume)) * (((ax) * ((-1.000000e+00) + ((x) * (x)))) + ((-2.000000e+00) * (bx) * (x) * ((-1.000000e+00) + (y))));
        case 10: // da2Nii3
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((-1.000000e+00) * (ax) * ((1.000000e+00) + (x)) * ((x) + ((-2.000000e+00) * (y)))) + ((bx) * (((2.000000e+00) * (x)) + ((-1.000000e+00) * (y))) * ((-1.000000e+00) + (y))));
        case 11: // da2Nii4
            return (1.0 / (element_volume)) * (((-2.000000e+00) * (ax) * ((1.000000e+00) + (x)) * (y)) + ((bx) * ((-1.000000e+00) + ((y) * (y)))));
        case 12: // da2Nii5
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((-1.000000e+00) * (bx) * ((1.000000e+00) + (y)) * (((2.000000e+00) * (x)) + (y))) + ((ax) * ((1.000000e+00) + (x)) * ((x) + ((2.000000e+00) * (y)))));
        case 13: // da2Nii6
            return (1.0 / (element_volume)) * ((ax) + ((-1.000000e+00) * (ax) * ((x) * (x))) + ((2.000000e+00) * (bx) * (x) * ((1.000000e+00) + (y))));
        case 14: // da2Nii7
            return (1.0 / (((-2.000000e+00) * (ay) * (bx)) + ((2.000000e+00) * (ax) * (by)))) * (((ax) * ((-1.000000e+00) + (x)) * ((x) + ((-2.000000e+00) * (y)))) + ((-1.000000e+00) * (bx) * (((2.000000e+00) * (x)) + ((-1.000000e+00) * (y))) * ((1.000000e+00) + (y))));
        case 15: // da2Nii8
            return (1.0 / (element_volume)) * ((bx) + ((2.000000e+00) * (ax) * ((-1.000000e+00) + (x)) * (y)) + ((-1.000000e+00) * (bx) * ((y) * (y))));
		
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");	
	}
	
}


	
double Element2DRect2nd::evaluate_form_function(
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {

	const double& x = local_reference_element_coords[0];
	const double& y = local_reference_element_coords[1];

	switch(elem_shape_func) {
		case 0:
			return -0.25 * x*x*y - 0.25 * x*y*y + 0.25 * x*x + 0.25 * y*y + 0.25 * x*y - 0.25;
		case 1:
			return 0.5 * x*x*y - 0.5 * x*x - 0.5 * y + 0.5; 
		case 2:
			return -0.25 * x*x*y + 0.25 * x*y*y + 0.25 * x*x + 0.25 * y*y - 0.25 * x*y - 0.25; 
		case 3:
			return - 0.5 * x*y*y - 0.5 * y*y + 0.5 * x + 0.5; 
		case 4: 
			return 0.25 * x*x*y + 0.25 * x*y*y + 0.25 * x*x + 0.25 * y*y + 0.25 * x*y - 0.25; 
		case 5: 
			return -0.5 * x*x*y - 0.5 * x*x + 0.5 * y + 0.5; 
		case 6:
			return 0.25 * x*x*y - 0.25 * x*y*y + 0.25 * x*x + 0.25 * y*y - 0.25 * x*y - 0.25; 
		case 7:
			return 0.5 * x*y*y - 0.5 * y*y - 0.5 * x + 0.5; 
		default:
			TDKP_GENERAL_EXCEPTION("unknown shape function");
	}	
}

double Element2DRect2nd::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	return evaluate_form_function_derivative(diffop, elem_shape_func, local_reference_element_coords[0], local_reference_element_coords[1]);
}				
					

/** evaluate form function in global coordinates */
double Element2DRect2nd::evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const {
	
	// -------------------------------------------
	// calculate local coordinates
	// -------------------------------------------	
	double global[] = {global_x, global_y};
	double local[2];
	global2local(global, local);
	// evaluate
	return evaluate_form_function(elem_shape_func, local);
	
}

void Element2DRect2nd::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	local_coords.resize(2);
	switch(lid) {
		case 0:
			local_coords[0] = -1.0;
			local_coords[1] = -1.0; 
			break;			
		case 1:
			local_coords[0] =  0.0;
			local_coords[1] = -1.0; 
			break;			
		case 2:
			local_coords[0] =  1.0;
			local_coords[1] = -1.0; 
			break;
		case 3:
			local_coords[0] =  1.0;
			local_coords[1] =  0.0; 
			break;
		case 4:
			local_coords[0] =  1.0;
			local_coords[1] =  1.0; 
			break;
		case 5:
			local_coords[0] =  0.0;
			local_coords[1] =  1.0; 
			break;			
		case 6:
			local_coords[0] = -1.0;
			local_coords[1] =  1.0; 
			break;
		case 7:
			local_coords[0] = -1.0;
			local_coords[1] =  0.0; 
			break;							
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
}



}
