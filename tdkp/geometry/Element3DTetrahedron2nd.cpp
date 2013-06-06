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

#include "tdkp/geometry/Element3DTetrahedron2nd.h"

namespace tdkp {
	
#define ADB jacobi_det
 		
Element3DTetrahedron2nd::Element3DTetrahedron2nd(unsigned int index_)
: Element3DTetrahedronBase(10)
{
	this->index_global    = index_;	
}

Element3DTetrahedron2nd::~Element3DTetrahedron2nd() {
	// naaaaaaaathhhhhhhiiiiiiinggggg
}

/** determine locator (global edge index) of node */		
void Element3DTetrahedron2nd::get_additional_node_locator(
	unsigned int additional_node_idx, 
	AdditionalNodeLocation& location_type, 
	 vector<unsigned int>& involved_vertices, 
	vector<double>& coords, 
	unsigned int& tag
) const {
	
	TDKP_BOUNDS_ASSERT(additional_node_idx < 6, "");
	// ---------------------------------
	// tag is not used
	// ---------------------------------
	tag = 0;
	location_type = edge_node;
	coords.resize(3);
	involved_vertices.resize(2);	
	// ---------------------------------
	// get coords
	// ---------------------------------
	unsigned int nidx0, nidx1;
	if(additional_node_idx < 3) {
		nidx0 = additional_node_idx;
		nidx1 = (additional_node_idx + 1) % 3;  	
	} else {
		nidx1 = 3;
		nidx0 = (additional_node_idx - 3);	
	}
	involved_vertices[0] = nodes[nidx0]->get_index_vertex();
	involved_vertices[1] = nodes[nidx1]->get_index_vertex();
	for(unsigned int ii = 0; ii < 3; ii++) {
		coords[ii] = 0.5 * (nodes[nidx0]->get_coord(ii) + nodes[nidx1]->get_coord(ii)); 	
	}		
}

void Element3DTetrahedron2nd::set_additional_node(unsigned int additional_node_idx, Node* node) {
	TDKP_ASSERT(additional_node_idx + 4 < nodes.size(), "");
	TDKP_ASSERT(nodes[additional_node_idx + 4] == 0, "");
	TDKP_BOUNDS_ASSERT(test_additional_node(*this, *node, additional_node_idx), "wrong additional node supplied: " << additional_node_idx);
	nodes[additional_node_idx + 4] = node;
	
	// -----------------------------------
	// set linear interpolation coefficents for additional node
	// -----------------------------------
	if(node->get_num_contributions() == 0) {
		unsigned int nidx0, nidx1;
		if(additional_node_idx < 3) {
			nidx0 = additional_node_idx;
			nidx1 = (additional_node_idx + 1) % 3;  	
		} else {
			nidx1 = 3;
			nidx0 = (additional_node_idx - 3);	
		}		
		
		// node sits in middle of edge, so halve the value of each vertex
		node->set_contribution(nodes[nidx0]->get_index_vertex(), 0.5);
		node->set_contribution(nodes[nidx1]->get_index_vertex(), 0.5);	
	}	
}


double  Element3DTetrahedron2nd::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {
	
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1
	                       + 10 * 10 * diffop_2
	                       + 10 * 10 * 3 * diffop_1;	
	switch(obj_index) {
        case 0: // da1db1Nii1Njj1
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8])));
        case 10: // da1db1Nii2Njj1
        case 1: // da1db1Nii1Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 20: // da1db1Nii3Njj1
        case 2: // da1db1Nii1Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 30: // da1db1Nii4Njj1
        case 3: // da1db1Nii1Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 40: // da1db1Nii5Njj1
        case 4: // da1db1Nii1Njj5
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((4.000000e+00) * jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-4.000000e+00) * jm[8]))) + (jm[3] * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 50: // da1db1Nii6Njj1
        case 5: // da1db1Nii1Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 60: // da1db1Nii7Njj1
        case 6: // da1db1Nii1Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-4.000000e+00) * jm[3] * jm[8]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 70: // da1db1Nii8Njj1
        case 7: // da1db1Nii1Njj8
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-4.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 80: // da1db1Nii9Njj1
        case 8: // da1db1Nii1Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 90: // da1db1Nii10Njj1
        case 9: // da1db1Nii1Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 11: // da1db1Nii2Njj2
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])));
        case 21: // da1db1Nii3Njj2
        case 12: // da1db1Nii2Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 31: // da1db1Nii4Njj2
        case 13: // da1db1Nii2Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 41: // da1db1Nii5Njj2
        case 14: // da1db1Nii2Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[4] * jm[6]) + ((-3.000000e+00) * jm[5] * jm[6]) + ((-3.000000e+00) * jm[3] * jm[7]) + ((4.000000e+00) * jm[5] * jm[7]) + ((3.000000e+00) * jm[3] * jm[8]) + ((-4.000000e+00) * jm[4] * jm[8])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 51: // da1db1Nii6Njj2
        case 15: // da1db1Nii2Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[5] * (((3.000000e+00) * jm[6]) + jm[7])) + ((-1.000000e+00) * (((3.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 61: // da1db1Nii7Njj2
        case 16: // da1db1Nii2Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])))) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 71: // da1db1Nii8Njj2
        case 17: // da1db1Nii2Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[3] + ((-1.000000e+00) * jm[4])) * jm[8])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 81: // da1db1Nii9Njj2
        case 18: // da1db1Nii2Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (((-1.000000e+00) * (((3.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (((3.000000e+00) * jm[6]) + jm[8])));
        case 91: // da1db1Nii10Njj2
        case 19: // da1db1Nii2Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 22: // da1db1Nii3Njj3
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])));
        case 32: // da1db1Nii4Njj3
        case 23: // da1db1Nii3Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 42: // da1db1Nii5Njj3
        case 24: // da1db1Nii3Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8])))) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 52: // da1db1Nii6Njj3
        case 25: // da1db1Nii3Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * ((jm[5] * (jm[6] + ((3.000000e+00) * jm[7]))) + ((-1.000000e+00) * (jm[3] + ((3.000000e+00) * jm[4])) * jm[8]));
        case 62: // da1db1Nii7Njj3
        case 26: // da1db1Nii3Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (((-4.000000e+00) * jm[5] * jm[6]) + ((-3.000000e+00) * jm[3] * jm[7]) + ((3.000000e+00) * jm[5] * jm[7]) + ((3.000000e+00) * jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((4.000000e+00) * jm[3] * jm[8]));
        case 72: // da1db1Nii8Njj3
        case 27: // da1db1Nii3Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 82: // da1db1Nii9Njj3
        case 28: // da1db1Nii3Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])))) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 92: // da1db1Nii10Njj3
        case 29: // da1db1Nii3Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (((-3.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((3.000000e+00) * jm[7]) + jm[8])));
        case 33: // da1db1Nii4Njj4
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])));
        case 43: // da1db1Nii5Njj4
        case 34: // da1db1Nii4Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 53: // da1db1Nii6Njj4
        case 35: // da1db1Nii4Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 63: // da1db1Nii7Njj4
        case 36: // da1db1Nii4Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 73: // da1db1Nii8Njj4
        case 37: // da1db1Nii4Njj8
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((4.000000e+00) * jm[4] * jm[6]) + ((-3.000000e+00) * jm[5] * jm[6]) + ((-4.000000e+00) * jm[3] * jm[7]) + ((3.000000e+00) * jm[5] * jm[7]) + ((3.000000e+00) * jm[3] * jm[8]) + ((-3.000000e+00) * jm[4] * jm[8]));
        case 83: // da1db1Nii9Njj4
        case 38: // da1db1Nii4Njj9
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((-1.000000e+00) * (jm[3] + ((3.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (jm[6] + ((3.000000e+00) * jm[8]))));
        case 93: // da1db1Nii10Njj4
        case 39: // da1db1Nii4Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((-1.000000e+00) * jm[4] * jm[6]) + ((-3.000000e+00) * jm[5] * jm[6]) + (jm[3] * (jm[7] + ((3.000000e+00) * jm[8]))));
        case 44: // da1db1Nii5Njj5
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[5] * jm[5]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + (jm[3] * jm[5] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-1.000000e+00) * jm[8]))) + ((jm[3] * jm[3]) * ((jm[7] + ((-1.000000e+00) * jm[8])) * (jm[7] + ((-1.000000e+00) * jm[8])))) + ((jm[4] * jm[4]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + (jm[4] * ((jm[3] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[8])) * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[5] * (((-2.000000e+00) * (jm[6] * jm[6])) + ((-2.000000e+00) * jm[7] * jm[8]) + (jm[6] * (jm[7] + jm[8])))))));
        case 54: // da1db1Nii6Njj5
        case 79: // da1db1Nii8Njj10
        case 97: // da1db1Nii10Njj8
        case 45: // da1db1Nii5Njj6
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[5] * jm[5]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * jm[5] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((-1.000000e+00) * (jm[4] * jm[4]) * jm[6] * jm[8]) + ((2.000000e+00) * (jm[3] * jm[3]) * jm[8] * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[4] * ((jm[3] * (((2.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[5] * jm[6] * (((-2.000000e+00) * jm[6]) + jm[7] + ((2.000000e+00) * jm[8]))))));
        case 64: // da1db1Nii7Njj5
        case 89: // da1db1Nii9Njj10
        case 98: // da1db1Nii10Njj9
        case 46: // da1db1Nii5Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[4] * jm[4]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * (jm[5] * jm[5]) * jm[6]) + ((jm[3] * jm[3]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[3] * jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))) + (jm[4] * ((jm[5] * jm[6] * (((-1.000000e+00) * jm[6]) + jm[7] + ((2.000000e+00) * jm[8]))) + (jm[3] * (((jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 59: // da1db1Nii6Njj10
        case 74: // da1db1Nii8Njj5
        case 95: // da1db1Nii10Njj6
        case 47: // da1db1Nii5Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[5] * jm[5]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((-2.000000e+00) * (jm[4] * jm[4]) * jm[6] * jm[8]) + ((jm[3] * jm[3]) * jm[8] * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[3] * jm[5] * ((jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[4] * ((jm[3] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[5] * jm[6] * (((-1.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])))));
        case 69: // da1db1Nii7Njj10
        case 84: // da1db1Nii9Njj5
        case 96: // da1db1Nii10Njj7
        case 48: // da1db1Nii5Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[4] * jm[4]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[7] * (((-1.000000e+00) * (jm[5] * jm[5]) * jm[6]) + ((2.000000e+00) * (jm[3] * jm[3]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[3] * jm[5] * (((2.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[4] * ((jm[3] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[8])) * (((-2.000000e+00) * jm[7]) + jm[8])) + (jm[5] * jm[6] * (((-2.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])))));
        case 94: // da1db1Nii10Njj5
        case 49: // da1db1Nii5Njj10
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8]))));
        case 55: // da1db1Nii6Njj6
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[5] * jm[5]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + (jm[5] * ((jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7]))) * jm[8]) + (((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[4]) + (jm[4] * jm[4])) * (jm[8] * jm[8])));
        case 65: // da1db1Nii7Njj6
        case 78: // da1db1Nii8Njj9
        case 87: // da1db1Nii9Njj8
        case 56: // da1db1Nii6Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[4] * jm[4]) * (jm[6] + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[4] * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[7])) * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[3] * jm[8] * (((-1.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))) + (jm[7] * (((2.000000e+00) * (jm[5] * jm[5]) * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((jm[3] * jm[3]) * jm[8]) + ((-1.000000e+00) * jm[3] * jm[5] * (jm[6] + ((-2.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))));
        case 75: // da1db1Nii8Njj6
        case 57: // da1db1Nii6Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8])));
        case 85: // da1db1Nii9Njj6
        case 58: // da1db1Nii6Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[4] * jm[4]) * jm[8] * (((-1.000000e+00) * jm[6]) + jm[8])) + (jm[7] * (((jm[5] * jm[5]) * (((-1.000000e+00) * jm[6]) + jm[7])) + ((-2.000000e+00) * (jm[3] * jm[3]) * jm[8]) + (jm[3] * jm[5] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[4] * ((jm[3] * (((2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[5] * (((-2.000000e+00) * (jm[6] * jm[6])) + ((-2.000000e+00) * jm[7] * jm[8]) + (jm[6] * (jm[7] + jm[8])))))));
        case 66: // da1db1Nii7Njj7
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[5] * jm[5]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + ((jm[4] * jm[4]) * ((jm[6] + ((-1.000000e+00) * jm[8])) * (jm[6] + ((-1.000000e+00) * jm[8])))) + ((-1.000000e+00) * jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])) * ((jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + ((jm[3] * jm[3]) * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[3] * jm[5] * ((jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[7]) + jm[8])))));
        case 76: // da1db1Nii8Njj7
        case 67: // da1db1Nii7Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[4] * jm[4]) * (jm[6] + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[7] * (((jm[5] * jm[5]) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * jm[5] * (((-2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8]))) + ((2.000000e+00) * (jm[3] * jm[3]) * jm[8]))) + (jm[4] * ((jm[3] * jm[8] * (((-2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7]) + jm[8])) + (jm[5] * (((2.000000e+00) * (jm[6] * jm[6])) + ((2.000000e+00) * jm[7] * jm[8]) + ((-1.000000e+00) * jm[6] * (jm[7] + jm[8])))))));
        case 86: // da1db1Nii9Njj7
        case 68: // da1db1Nii7Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8]))) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8]))));
        case 77: // da1db1Nii8Njj8
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[5] * jm[5]) * ((jm[6] + ((-1.000000e+00) * jm[7])) * (jm[6] + ((-1.000000e+00) * jm[7])))) + (jm[3] * jm[5] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((jm[4] * jm[4]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + ((jm[3] * jm[3]) * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[5] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[3] * (((jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 88: // da1db1Nii9Njj9
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[5]) + (jm[5] * jm[5])) * (jm[7] * jm[7])) + ((jm[4] * jm[4]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + (jm[4] * jm[7] * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[8])))));
        case 99: // da1db1Nii10Njj10
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[4] * jm[4]) * (jm[6] * jm[6])) + ((jm[5] * jm[5]) * (jm[6] * jm[6])) + (jm[3] * jm[5] * jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((jm[3] * jm[3]) * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[4] * jm[6] * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[7]) + jm[8])))));
        case 300: // da2db1Nii1Njj1
        case 100: // da1db2Nii1Njj1
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 310: // da2db1Nii2Njj1
        case 101: // da1db2Nii1Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 320: // da2db1Nii3Njj1
        case 102: // da1db2Nii1Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 330: // da2db1Nii4Njj1
        case 103: // da1db2Nii1Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 340: // da2db1Nii5Njj1
        case 104: // da1db2Nii1Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * jm[7]) + ((-4.000000e+00) * jm[2] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8]) + ((4.000000e+00) * jm[1] * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 350: // da2db1Nii6Njj1
        case 105: // da1db2Nii1Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 360: // da2db1Nii7Njj1
        case 106: // da1db2Nii1Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[2] * jm[6]) + (jm[0] * jm[7]) + ((-1.000000e+00) * jm[2] * jm[7]) + ((-4.000000e+00) * jm[0] * jm[8]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 370: // da2db1Nii8Njj1
        case 107: // da1db2Nii1Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + ((-4.000000e+00) * jm[0] * jm[7]) + (jm[2] * jm[7]) + (jm[0] * jm[8]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 380: // da2db1Nii9Njj1
        case 108: // da1db2Nii1Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 390: // da2db1Nii10Njj1
        case 109: // da1db2Nii1Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 301: // da2db1Nii1Njj2
        case 110: // da1db2Nii2Njj1
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 311: // da2db1Nii2Njj2
        case 111: // da1db2Nii2Njj2
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 321: // da2db1Nii3Njj2
        case 112: // da1db2Nii2Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 331: // da2db1Nii4Njj2
        case 113: // da1db2Nii2Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 341: // da2db1Nii5Njj2
        case 114: // da1db2Nii2Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[1] * jm[6]) + ((-3.000000e+00) * jm[2] * jm[6]) + ((-3.000000e+00) * jm[0] * jm[7]) + ((4.000000e+00) * jm[2] * jm[7]) + ((3.000000e+00) * jm[0] * jm[8]) + ((-4.000000e+00) * jm[1] * jm[8])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 351: // da2db1Nii6Njj2
        case 115: // da1db2Nii2Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((3.000000e+00) * jm[6]) + jm[7])) + ((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[1]) * jm[8])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 361: // da2db1Nii7Njj2
        case 116: // da1db2Nii2Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 371: // da2db1Nii8Njj2
        case 117: // da1db2Nii2Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 381: // da2db1Nii9Njj2
        case 118: // da1db2Nii2Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (((((3.000000e+00) * jm[0]) + jm[2]) * jm[7]) + ((-1.000000e+00) * jm[1] * (((3.000000e+00) * jm[6]) + jm[8])));
        case 391: // da2db1Nii10Njj2
        case 119: // da1db2Nii2Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 302: // da2db1Nii1Njj3
        case 120: // da1db2Nii3Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 312: // da2db1Nii2Njj3
        case 121: // da1db2Nii3Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 322: // da2db1Nii3Njj3
        case 122: // da1db2Nii3Njj3
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 332: // da2db1Nii4Njj3
        case 123: // da1db2Nii3Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 342: // da2db1Nii5Njj3
        case 124: // da1db2Nii3Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 352: // da2db1Nii6Njj3
        case 125: // da1db2Nii3Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((3.000000e+00) * jm[7]))) + ((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[1])) * jm[8])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 362: // da2db1Nii7Njj3
        case 126: // da1db2Nii3Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-3.000000e+00) * jm[1] * jm[6]) + ((4.000000e+00) * jm[2] * jm[6]) + ((3.000000e+00) * jm[0] * jm[7]) + ((-3.000000e+00) * jm[2] * jm[7]) + ((-4.000000e+00) * jm[0] * jm[8]) + ((3.000000e+00) * jm[1] * jm[8])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 372: // da2db1Nii8Njj3
        case 127: // da1db2Nii3Njj8
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 382: // da2db1Nii9Njj3
        case 128: // da1db2Nii3Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 392: // da2db1Nii10Njj3
        case 129: // da1db2Nii3Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (((3.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * (((3.000000e+00) * jm[7]) + jm[8])));
        case 303: // da2db1Nii1Njj4
        case 130: // da1db2Nii4Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 313: // da2db1Nii2Njj4
        case 131: // da1db2Nii4Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 323: // da2db1Nii3Njj4
        case 132: // da1db2Nii4Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 333: // da2db1Nii4Njj4
        case 133: // da1db2Nii4Njj4
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 343: // da2db1Nii5Njj4
        case 134: // da1db2Nii4Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 353: // da2db1Nii6Njj4
        case 135: // da1db2Nii4Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 363: // da2db1Nii7Njj4
        case 136: // da1db2Nii4Njj7
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 373: // da2db1Nii8Njj4
        case 137: // da1db2Nii4Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((4.000000e+00) * jm[1] * jm[6]) + ((-3.000000e+00) * jm[2] * jm[6]) + ((-4.000000e+00) * jm[0] * jm[7]) + ((3.000000e+00) * jm[2] * jm[7]) + ((3.000000e+00) * jm[0] * jm[8]) + ((-3.000000e+00) * jm[1] * jm[8]));
        case 383: // da2db1Nii9Njj4
        case 138: // da1db2Nii4Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (jm[6] + ((3.000000e+00) * jm[8]))));
        case 393: // da2db1Nii10Njj4
        case 139: // da1db2Nii4Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * ((jm[1] * jm[6]) + ((3.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * (jm[7] + ((3.000000e+00) * jm[8]))));
        case 304: // da2db1Nii1Njj5
        case 140: // da1db2Nii5Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-4.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((4.000000e+00) * jm[4] * jm[8])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 314: // da2db1Nii2Njj5
        case 141: // da1db2Nii5Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (((3.000000e+00) * jm[4] * jm[6]) + ((-3.000000e+00) * jm[5] * jm[6]) + ((-3.000000e+00) * jm[3] * jm[7]) + ((4.000000e+00) * jm[5] * jm[7]) + ((3.000000e+00) * jm[3] * jm[8]) + ((-4.000000e+00) * jm[4] * jm[8]));
        case 324: // da2db1Nii3Njj5
        case 142: // da1db2Nii5Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 334: // da2db1Nii4Njj5
        case 143: // da1db2Nii5Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 344: // da2db1Nii5Njj5
        case 144: // da1db2Nii5Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])) * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[2] * (((2.000000e+00) * jm[5] * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + (jm[3] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-2.000000e+00) * (jm[6] * jm[6])) + ((-2.000000e+00) * jm[7] * jm[8]) + (jm[6] * (jm[7] + jm[8])))))) + (jm[1] * ((jm[3] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[8])) * (((-1.000000e+00) * jm[7]) + jm[8])) + ((2.000000e+00) * jm[4] * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + (jm[5] * (((-2.000000e+00) * (jm[6] * jm[6])) + ((-2.000000e+00) * jm[7] * jm[8]) + (jm[6] * (jm[7] + jm[8])))))));
        case 354: // da2db1Nii6Njj5
        case 145: // da1db2Nii5Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((2.000000e+00) * jm[5] * jm[6] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * jm[6] * (((-2.000000e+00) * jm[6]) + jm[7] + jm[8])))) + (jm[8] * ((jm[1] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[0] * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))))));
        case 364: // da2db1Nii7Njj5
        case 146: // da1db2Nii5Njj7
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))) + (jm[2] * ((jm[4] * jm[6] * (jm[7] + jm[8])) + (jm[7] * (((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))))) + (jm[0] * (((-1.000000e+00) * jm[4] * ((jm[6] * jm[7]) + (jm[8] * jm[8]))) + (jm[7] * ((jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[5] * (jm[6] + jm[8])))))));
        case 374: // da2db1Nii8Njj5
        case 147: // da1db2Nii5Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))) + ((-1.000000e+00) * jm[1] * ((jm[5] * jm[6] * (jm[7] + jm[8])) + (jm[8] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))))) + ((-1.000000e+00) * jm[0] * (((-1.000000e+00) * jm[5] * ((jm[7] * jm[7]) + (jm[6] * jm[8]))) + (jm[8] * ((jm[4] * (jm[6] + jm[7])) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])))))));
        case 384: // da2db1Nii9Njj5
        case 148: // da1db2Nii5Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (((2.000000e+00) * jm[4] * jm[6] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[8])) * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[5] * jm[6] * (((-2.000000e+00) * jm[6]) + jm[7] + jm[8])))) + (jm[7] * ((jm[2] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))) + (jm[0] * (((2.000000e+00) * jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])))))));
        case 194: // da1db2Nii10Njj5
        case 349: // da2db1Nii5Njj10
        case 394: // da2db1Nii10Njj5
        case 149: // da1db2Nii5Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 305: // da2db1Nii1Njj6
        case 150: // da1db2Nii6Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 315: // da2db1Nii2Njj6
        case 151: // da1db2Nii6Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[5] * (((3.000000e+00) * jm[6]) + jm[7])) + ((-1.000000e+00) * (((3.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 325: // da2db1Nii3Njj6
        case 152: // da1db2Nii6Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[5] * (jm[6] + ((3.000000e+00) * jm[7]))) + ((-1.000000e+00) * (jm[3] + ((3.000000e+00) * jm[4])) * jm[8]));
        case 335: // da2db1Nii4Njj6
        case 153: // da1db2Nii6Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 345: // da2db1Nii5Njj6
        case 154: // da1db2Nii6Njj5
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])) * (((2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))) + (jm[2] * (((2.000000e+00) * jm[5] * jm[6] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-2.000000e+00) * jm[3] * jm[6]) + (jm[4] * jm[6]) + (jm[3] * jm[7])) * jm[8]))) + (jm[1] * ((jm[5] * jm[6] * (((-2.000000e+00) * jm[6]) + jm[7] + jm[8])) + ((-1.000000e+00) * jm[8] * (((-2.000000e+00) * jm[3] * jm[6]) + (jm[4] * jm[6]) + (jm[3] * jm[8]))))));
        case 355: // da2db1Nii6Njj6
        case 155: // da1db2Nii6Njj6
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((2.000000e+00) * jm[5] * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + ((((-2.000000e+00) * jm[3] * jm[6]) + (jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[7])) * jm[8]))) + (jm[8] * ((jm[0] * (((-2.000000e+00) * jm[5] * jm[6]) + (jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[1] * ((jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))))));
        case 365: // da2db1Nii7Njj6
        case 156: // da1db2Nii6Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])) * ((jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))) + (jm[2] * (((2.000000e+00) * jm[5] * jm[7] * (((-1.000000e+00) * jm[6]) + jm[7])) + (((jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[7])) * jm[8]))) + (jm[0] * ((jm[5] * jm[7] * (jm[6] + ((-2.000000e+00) * jm[7]) + jm[8])) + ((-1.000000e+00) * jm[8] * ((jm[3] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[7]) + (jm[4] * jm[8]))))));
        case 175: // da1db2Nii8Njj6
        case 357: // da2db1Nii6Njj8
        case 375: // da2db1Nii8Njj6
        case 157: // da1db2Nii6Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 385: // da2db1Nii9Njj6
        case 158: // da1db2Nii6Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[7] * ((jm[0] * (((-2.000000e+00) * jm[5] * jm[6]) + (jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[2] * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))))) + (jm[1] * ((jm[8] * ((jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[8])))) + (jm[5] * (((2.000000e+00) * (jm[6] * jm[6])) + (jm[7] * jm[8]) + ((-1.000000e+00) * jm[6] * (jm[7] + jm[8])))))));
        case 395: // da2db1Nii10Njj6
        case 159: // da1db2Nii6Njj10
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6] * (((-1.000000e+00) * jm[5] * jm[6]) + ((2.000000e+00) * jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[2] * jm[6] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]))) + (jm[0] * ((jm[8] * (((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[4] * jm[7]) + (jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[5] * ((jm[6] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 306: // da2db1Nii1Njj7
        case 160: // da1db2Nii7Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * (((4.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-4.000000e+00) * jm[3] * jm[8]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 316: // da2db1Nii2Njj7
        case 161: // da1db2Nii7Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 326: // da2db1Nii3Njj7
        case 162: // da1db2Nii7Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (((-3.000000e+00) * jm[4] * jm[6]) + ((4.000000e+00) * jm[5] * jm[6]) + ((3.000000e+00) * jm[3] * jm[7]) + ((-3.000000e+00) * jm[5] * jm[7]) + ((-4.000000e+00) * jm[3] * jm[8]) + ((3.000000e+00) * jm[4] * jm[8]));
        case 336: // da2db1Nii4Njj7
        case 163: // da1db2Nii7Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 346: // da2db1Nii5Njj7
        case 164: // da1db2Nii7Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])))) + (jm[1] * ((jm[4] * jm[6] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[5] * jm[6] * (jm[7] + jm[8])) + ((-1.000000e+00) * jm[3] * ((jm[6] * jm[7]) + (jm[8] * jm[8]))))) + (jm[2] * ((jm[4] * jm[6] * (((-1.000000e+00) * jm[6]) + jm[8])) + (jm[7] * (((-2.000000e+00) * jm[5] * jm[6]) + (jm[3] * (jm[6] + jm[8])))))));
        case 356: // da2db1Nii6Njj7
        case 165: // da1db2Nii7Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[8] * ((jm[0] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[1] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))))) + (jm[2] * (((-1.000000e+00) * jm[4] * (jm[6] + ((-2.000000e+00) * jm[7])) * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[7] * (((2.000000e+00) * jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[3] * (jm[6] + ((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 366: // da2db1Nii7Njj7
        case 166: // da1db2Nii7Njj7
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])) * (((-1.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * jm[8]))) + (jm[2] * (((2.000000e+00) * jm[5] * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + ((-1.000000e+00) * jm[4] * (jm[6] + ((-2.000000e+00) * jm[7])) * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * ((jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[7]) + jm[8])))))) + (jm[0] * ((jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])) * (((-2.000000e+00) * jm[7]) + jm[8])) + ((2.000000e+00) * jm[3] * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[5] * ((jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 376: // da2db1Nii8Njj7
        case 167: // da1db2Nii7Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7])) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))))) + ((-1.000000e+00) * jm[0] * ((jm[5] * jm[7] * (jm[6] + jm[8])) + (jm[8] * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))))) + (jm[1] * ((jm[5] * ((jm[6] * jm[6]) + (jm[7] * jm[8]))) + ((-1.000000e+00) * jm[8] * ((jm[3] * (jm[6] + jm[7])) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])))))));
        case 186: // da1db2Nii9Njj7
        case 368: // da2db1Nii7Njj9
        case 386: // da2db1Nii9Njj7
        case 168: // da1db2Nii7Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 396: // da2db1Nii10Njj7
        case 169: // da1db2Nii7Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6] * (((-1.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * jm[8]))) + (jm[2] * jm[6] * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])))) + (jm[0] * ((jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])) * (((-2.000000e+00) * jm[7]) + jm[8])) + (jm[7] * (((2.000000e+00) * jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[5] * (jm[6] + ((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 307: // da2db1Nii1Njj8
        case 170: // da1db2Nii8Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-4.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 317: // da2db1Nii2Njj8
        case 171: // da1db2Nii8Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[3] + ((-1.000000e+00) * jm[4])) * jm[8]));
        case 327: // da2db1Nii3Njj8
        case 172: // da1db2Nii8Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 337: // da2db1Nii4Njj8
        case 173: // da1db2Nii8Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((4.000000e+00) * jm[4] * jm[6]) + ((-3.000000e+00) * jm[5] * jm[6]) + ((-4.000000e+00) * jm[3] * jm[7]) + ((3.000000e+00) * jm[5] * jm[7]) + ((3.000000e+00) * jm[3] * jm[8]) + ((-3.000000e+00) * jm[4] * jm[8]));
        case 347: // da2db1Nii5Njj8
        case 174: // da1db2Nii8Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]))) + (jm[1] * ((jm[5] * jm[6] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((((-2.000000e+00) * jm[4] * jm[6]) + (jm[3] * (jm[6] + jm[7]))) * jm[8]))) + (jm[2] * ((jm[5] * jm[6] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[4] * jm[6] * (jm[7] + jm[8])) + ((-1.000000e+00) * jm[3] * ((jm[7] * jm[7]) + (jm[6] * jm[8]))))));
        case 367: // da2db1Nii7Njj8
        case 176: // da1db2Nii8Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]))) + ((-1.000000e+00) * jm[0] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7])) * jm[7]) + ((((-2.000000e+00) * jm[3] * jm[7]) + (jm[4] * (jm[6] + jm[7]))) * jm[8]))) + (jm[2] * ((jm[4] * ((jm[6] * jm[6]) + (jm[7] * jm[8]))) + ((-1.000000e+00) * jm[7] * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[3] * (jm[6] + jm[8])))))));
        case 377: // da2db1Nii8Njj8
        case 177: // da1db2Nii8Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7])) * (((-1.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[5] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[6] + ((-2.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[4] * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + (jm[3] * (((jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])))))) + (jm[0] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[3] * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[4] * (((jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 387: // da2db1Nii9Njj8
        case 178: // da1db2Nii8Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[7] * ((jm[2] * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))))) + (jm[1] * (((-1.000000e+00) * jm[5] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[8] * ((jm[3] * (jm[6] + jm[7] + ((-2.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])))))));
        case 397: // da2db1Nii10Njj8
        case 179: // da1db2Nii8Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6] * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[3] + ((-1.000000e+00) * jm[4])) * jm[8]))) + (jm[2] * jm[6] * (((-1.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[8] * ((jm[4] * (jm[6] + jm[7] + ((-2.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])))))));
        case 308: // da2db1Nii1Njj9
        case 180: // da1db2Nii9Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 318: // da2db1Nii2Njj9
        case 181: // da1db2Nii9Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (((((3.000000e+00) * jm[3]) + jm[5]) * jm[7]) + ((-1.000000e+00) * jm[4] * (((3.000000e+00) * jm[6]) + jm[8])));
        case 328: // da2db1Nii3Njj9
        case 182: // da1db2Nii9Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 338: // da2db1Nii4Njj9
        case 183: // da1db2Nii9Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((-1.000000e+00) * (jm[3] + ((3.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (jm[6] + ((3.000000e+00) * jm[8]))));
        case 348: // da2db1Nii5Njj9
        case 184: // da1db2Nii9Njj5
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])) * (((-2.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]))) + (jm[2] * (((-1.000000e+00) * jm[7] * (((-2.000000e+00) * jm[3] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]))) + (jm[4] * jm[6] * (((-2.000000e+00) * jm[6]) + jm[7] + jm[8])))) + (jm[1] * (((2.000000e+00) * jm[4] * jm[6] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[3] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[8]))))));
        case 358: // da2db1Nii6Njj9
        case 185: // da1db2Nii9Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[8] * ((jm[1] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]))))) + (jm[2] * ((jm[7] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7])))) + (jm[4] * (((2.000000e+00) * (jm[6] * jm[6])) + (jm[7] * jm[8]) + ((-1.000000e+00) * jm[6] * (jm[7] + jm[8])))))));
        case 378: // da2db1Nii8Njj9
        case 187: // da1db2Nii9Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[2] * (jm[6] + ((-1.000000e+00) * jm[7])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[5] * jm[7] * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((((-1.000000e+00) * jm[3] * jm[7]) + (jm[4] * (jm[6] + jm[7] + ((-2.000000e+00) * jm[8])))) * jm[8]))) + (jm[1] * ((jm[5] * jm[7] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[8] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]) + ((2.000000e+00) * jm[4] * jm[8]))))));
        case 388: // da2db1Nii9Njj9
        case 188: // da1db2Nii9Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[7] * ((jm[2] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]))))) + (jm[1] * ((jm[7] * (((-2.000000e+00) * jm[3] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[8]) + ((-2.000000e+00) * jm[5] * jm[8]))) + ((2.000000e+00) * jm[4] * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))))));
        case 398: // da2db1Nii10Njj9
        case 189: // da1db2Nii9Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[1] * jm[6] * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])))) + (jm[0] * ((jm[7] * (((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-2.000000e+00) * jm[5] * jm[8]))) + (jm[4] * ((jm[6] * jm[7]) + ((-1.000000e+00) * jm[6] * jm[8]) + ((-1.000000e+00) * jm[7] * jm[8]) + ((2.000000e+00) * (jm[8] * jm[8])))))));
        case 309: // da2db1Nii1Njj10
        case 190: // da1db2Nii10Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 319: // da2db1Nii2Njj10
        case 191: // da1db2Nii10Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 329: // da2db1Nii3Njj10
        case 192: // da1db2Nii10Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (((3.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * (((3.000000e+00) * jm[7]) + jm[8])));
        case 339: // da2db1Nii4Njj10
        case 193: // da1db2Nii10Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[4] * jm[6]) + ((3.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * (jm[7] + ((3.000000e+00) * jm[8]))));
        case 359: // da2db1Nii6Njj10
        case 195: // da1db2Nii10Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[8] * ((jm[0] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[1] * (((2.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))))) + (jm[2] * ((jm[4] * jm[6] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[5] * jm[6] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[3] * (((-1.000000e+00) * jm[6] * jm[7]) + ((2.000000e+00) * (jm[7] * jm[7])) + (jm[6] * jm[8]) + ((-1.000000e+00) * jm[7] * jm[8]))))));
        case 369: // da2db1Nii7Njj10
        case 196: // da1db2Nii10Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])) * (((2.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))) + (jm[0] * ((jm[4] * jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])) + (jm[7] * ((jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]))))) + (jm[2] * (((-1.000000e+00) * jm[4] * jm[6] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[7] * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (jm[6] + ((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 379: // da2db1Nii8Njj10
        case 197: // da1db2Nii10Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7])) * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]))) + (jm[1] * ((jm[5] * jm[6] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * (jm[6] + jm[7] + ((-2.000000e+00) * jm[8])))) * jm[8]))) + ((-1.000000e+00) * jm[0] * ((jm[5] * jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[8] * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]))))));
        case 389: // da2db1Nii9Njj10
        case 198: // da1db2Nii10Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[7] * ((jm[2] * (((-1.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]))) + (jm[0] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))))) + (jm[1] * ((jm[5] * jm[6] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[4] * jm[6] * (((-1.000000e+00) * jm[6]) + jm[8])) + (jm[3] * ((jm[6] * jm[7]) + ((-1.000000e+00) * jm[6] * jm[8]) + ((-1.000000e+00) * jm[7] * jm[8]) + ((2.000000e+00) * (jm[8] * jm[8])))))));
        case 399: // da2db1Nii10Njj10
        case 199: // da1db2Nii10Njj10
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6] * (((-1.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-2.000000e+00) * jm[8]))))) + (jm[1] * jm[6] * (((2.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[0] * ((jm[5] * jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[4] * jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])) + ((2.000000e+00) * jm[3] * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))))));
        case 600: // da3db1Nii1Njj1
        case 200: // da1db3Nii1Njj1
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 610: // da3db1Nii2Njj1
        case 201: // da1db3Nii1Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 620: // da3db1Nii3Njj1
        case 202: // da1db3Nii1Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 630: // da3db1Nii4Njj1
        case 203: // da1db3Nii1Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 640: // da3db1Nii5Njj1
        case 204: // da1db3Nii1Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + ((4.000000e+00) * jm[2] * jm[4]) + (jm[1] * (jm[3] + ((-4.000000e+00) * jm[5]))) + (jm[0] * jm[5])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 650: // da3db1Nii6Njj1
        case 205: // da1db3Nii1Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 660: // da3db1Nii7Njj1
        case 206: // da1db3Nii1Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-4.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + (jm[2] * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + ((4.000000e+00) * jm[0] * jm[5])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 670: // da3db1Nii8Njj1
        case 207: // da1db3Nii1Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + ((-4.000000e+00) * jm[0] * jm[4]) + (jm[2] * jm[4]) + (jm[0] * jm[5]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 680: // da3db1Nii9Njj1
        case 208: // da1db3Nii1Njj9
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 690: // da3db1Nii10Njj1
        case 209: // da1db3Nii1Njj10
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 601: // da3db1Nii1Njj2
        case 210: // da1db3Nii2Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 611: // da3db1Nii2Njj2
        case 211: // da1db3Nii2Njj2
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 621: // da3db1Nii3Njj2
        case 212: // da1db3Nii2Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]));
        case 631: // da3db1Nii4Njj2
        case 213: // da1db3Nii2Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 641: // da3db1Nii5Njj2
        case 214: // da1db3Nii2Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[1] * jm[3]) + ((-3.000000e+00) * jm[2] * jm[3]) + ((-3.000000e+00) * jm[0] * jm[4]) + ((4.000000e+00) * jm[2] * jm[4]) + ((3.000000e+00) * jm[0] * jm[5]) + ((-4.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 651: // da3db1Nii6Njj2
        case 215: // da1db3Nii2Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((3.000000e+00) * jm[3]) + jm[4])) + ((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 661: // da3db1Nii7Njj2
        case 216: // da1db3Nii2Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 671: // da3db1Nii8Njj2
        case 217: // da1db3Nii2Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 681: // da3db1Nii9Njj2
        case 218: // da1db3Nii2Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((3.000000e+00) * jm[0]) + jm[2]) * jm[4]) + ((-1.000000e+00) * jm[1] * (((3.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 691: // da3db1Nii10Njj2
        case 219: // da1db3Nii2Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]));
        case 602: // da3db1Nii1Njj3
        case 220: // da1db3Nii3Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 612: // da3db1Nii2Njj3
        case 221: // da1db3Nii3Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]));
        case 622: // da3db1Nii3Njj3
        case 222: // da1db3Nii3Njj3
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 632: // da3db1Nii4Njj3
        case 223: // da1db3Nii3Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 642: // da3db1Nii5Njj3
        case 224: // da1db3Nii3Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 652: // da3db1Nii6Njj3
        case 225: // da1db3Nii3Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((3.000000e+00) * jm[4]))) + ((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[1])) * jm[5])) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 662: // da3db1Nii7Njj3
        case 226: // da1db3Nii3Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-3.000000e+00) * jm[1] * jm[3]) + ((4.000000e+00) * jm[2] * jm[3]) + ((3.000000e+00) * jm[0] * jm[4]) + ((-3.000000e+00) * jm[2] * jm[4]) + ((-4.000000e+00) * jm[0] * jm[5]) + ((3.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 672: // da3db1Nii8Njj3
        case 227: // da1db3Nii3Njj8
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 682: // da3db1Nii9Njj3
        case 228: // da1db3Nii3Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 692: // da3db1Nii10Njj3
        case 229: // da1db3Nii3Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * (((3.000000e+00) * jm[4]) + jm[5]))) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8]));
        case 603: // da3db1Nii1Njj4
        case 230: // da1db3Nii4Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 613: // da3db1Nii2Njj4
        case 231: // da1db3Nii4Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]));
        case 623: // da3db1Nii3Njj4
        case 232: // da1db3Nii4Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 633: // da3db1Nii4Njj4
        case 233: // da1db3Nii4Njj4
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[0] * jm[4])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 643: // da3db1Nii5Njj4
        case 234: // da1db3Nii4Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])))) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 653: // da3db1Nii6Njj4
        case 235: // da1db3Nii4Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 663: // da3db1Nii7Njj4
        case 236: // da1db3Nii4Njj7
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]));
        case 673: // da3db1Nii8Njj4
        case 237: // da1db3Nii4Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[1] * jm[3]) + ((-3.000000e+00) * jm[2] * jm[3]) + ((-4.000000e+00) * jm[0] * jm[4]) + ((3.000000e+00) * jm[2] * jm[4]) + ((3.000000e+00) * jm[0] * jm[5]) + ((-3.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 683: // da3db1Nii9Njj4
        case 238: // da1db3Nii4Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (jm[3] + ((3.000000e+00) * jm[5])))) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 693: // da3db1Nii10Njj4
        case 239: // da1db3Nii4Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((3.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * (jm[4] + ((3.000000e+00) * jm[5])))) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]));
        case 604: // da3db1Nii1Njj5
        case 240: // da1db3Nii5Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-4.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((4.000000e+00) * jm[4] * jm[8]));
        case 614: // da3db1Nii2Njj5
        case 241: // da1db3Nii5Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((3.000000e+00) * jm[4] * jm[6]) + ((-3.000000e+00) * jm[5] * jm[6]) + ((-3.000000e+00) * jm[3] * jm[7]) + ((4.000000e+00) * jm[5] * jm[7]) + ((3.000000e+00) * jm[3] * jm[8]) + ((-4.000000e+00) * jm[4] * jm[8]));
        case 624: // da3db1Nii3Njj5
        case 242: // da1db3Nii5Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 634: // da3db1Nii4Njj5
        case 243: // da1db3Nii5Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 644: // da3db1Nii5Njj5
        case 244: // da1db3Nii5Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])) * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[2] * ((jm[3] * ((jm[5] * (((-2.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (((2.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8]))))) + ((-2.000000e+00) * (jm[3] * jm[3]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[4] * jm[8]))))) + (jm[1] * (((2.000000e+00) * (jm[3] * jm[3]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[5] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[3] * ((jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])) + (jm[5] * (((2.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 654: // da3db1Nii6Njj5
        case 245: // da1db3Nii5Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * ((jm[4] * (((-1.000000e+00) * jm[4]) + jm[5]) * jm[6]) + (jm[3] * ((jm[5] * (((-2.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (((2.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8]))))) + ((-2.000000e+00) * (jm[3] * jm[3]) * (jm[7] + ((-1.000000e+00) * jm[8]))))) + (jm[5] * ((jm[1] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))) + (jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))))));
        case 664: // da3db1Nii7Njj5
        case 246: // da1db3Nii5Njj7
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[2] * (((-1.000000e+00) * (jm[4] * jm[4]) * jm[6]) + (jm[3] * jm[5] * jm[7]) + (jm[4] * ((jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]))))) + (jm[0] * (((jm[4] * jm[4]) * jm[6]) + ((-1.000000e+00) * (jm[5] * jm[5]) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]) + (jm[5] * jm[8]))))));
        case 674: // da3db1Nii8Njj5
        case 247: // da1db3Nii5Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[5] * (((-1.000000e+00) * jm[4]) + jm[5]) * jm[6]) + ((-1.000000e+00) * (((-2.000000e+00) * jm[1] * jm[3]) + (jm[0] * (jm[3] + jm[4]))) * jm[5] * jm[7]) + (jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))))) + ((((-1.000000e+00) * jm[1] * jm[3] * (jm[4] + jm[5])) + (jm[0] * ((jm[4] * jm[4]) + (jm[3] * jm[5])))) * jm[8]));
        case 684: // da3db1Nii9Njj5
        case 248: // da1db3Nii5Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[4] * ((jm[2] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[0] * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))))) + (jm[1] * (((jm[4] + ((-1.000000e+00) * jm[5])) * jm[5] * jm[6]) + ((2.000000e+00) * (jm[3] * jm[3]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[3] * ((jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])) + (jm[5] * (((2.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 294: // da1db3Nii10Njj5
        case 649: // da3db1Nii5Njj10
        case 694: // da3db1Nii10Njj5
        case 249: // da1db3Nii5Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 605: // da3db1Nii1Njj6
        case 250: // da1db3Nii6Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 615: // da3db1Nii2Njj6
        case 251: // da1db3Nii6Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * jm[5] * (((3.000000e+00) * jm[6]) + jm[7])) + ((((3.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 625: // da3db1Nii3Njj6
        case 252: // da1db3Nii6Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-1.000000e+00) * jm[5] * (jm[6] + ((3.000000e+00) * jm[7]))) + ((jm[3] + ((3.000000e+00) * jm[4])) * jm[8]));
        case 635: // da3db1Nii4Njj6
        case 253: // da1db3Nii6Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 645: // da3db1Nii5Njj6
        case 254: // da1db3Nii6Njj5
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * ((jm[5] * ((jm[4] * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7])))) + ((2.000000e+00) * jm[3] * (jm[3] + ((-1.000000e+00) * jm[4])) * jm[8]))) + ((-1.000000e+00) * jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])) * (((2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[5] * ((jm[5] * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7])))) + (jm[3] * (((-2.000000e+00) * jm[3]) + jm[4] + jm[5]) * jm[8]))));
        case 655: // da3db1Nii6Njj6
        case 255: // da1db3Nii6Njj6
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[5] * ((jm[5] * (((2.000000e+00) * jm[0] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]) + ((2.000000e+00) * jm[1] * jm[7]))) + (jm[2] * ((jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7])))))) + ((((2.000000e+00) * jm[2] * ((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[4]) + (jm[4] * jm[4]))) + (((jm[1] * (jm[3] + ((-2.000000e+00) * jm[4]))) + (jm[0] * (((-2.000000e+00) * jm[3]) + jm[4]))) * jm[5])) * jm[8]));
        case 665: // da3db1Nii7Njj6
        case 256: // da1db3Nii6Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[5] * ((jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[3] * jm[7]))) + ((-1.000000e+00) * jm[0] * jm[5] * ((jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[5] * jm[7]))) + (jm[4] * (((2.000000e+00) * jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[0] * (jm[3] + ((-2.000000e+00) * jm[4]) + jm[5]))) * jm[8]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])) * ((jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))));
        case 275: // da1db3Nii8Njj6
        case 657: // da3db1Nii6Njj8
        case 675: // da3db1Nii8Njj6
        case 257: // da1db3Nii6Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[3] + ((-1.000000e+00) * jm[4])) * jm[8]));
        case 685: // da3db1Nii9Njj6
        case 258: // da1db3Nii6Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[5] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7])))) + (jm[1] * (((2.000000e+00) * (jm[3] * jm[3])) + (jm[4] * jm[5]) + ((-1.000000e+00) * jm[3] * (jm[4] + jm[5]))) * jm[8]) + (jm[2] * jm[4] * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[3] + ((-1.000000e+00) * jm[4])) * jm[8]))) + (jm[0] * jm[4] * (((2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))));
        case 695: // da3db1Nii10Njj6
        case 259: // da1db3Nii6Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3] * (((-1.000000e+00) * jm[5] * jm[6]) + ((2.000000e+00) * jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[2] * jm[3] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]))) + (jm[0] * (((2.000000e+00) * (jm[4] * jm[4]) * jm[8]) + (jm[4] * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8]))) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[5] * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[3] * jm[8]))))));
        case 606: // da3db1Nii1Njj7
        case 260: // da1db3Nii7Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-4.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((4.000000e+00) * jm[3] * jm[8]));
        case 616: // da3db1Nii2Njj7
        case 261: // da1db3Nii7Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 626: // da3db1Nii3Njj7
        case 262: // da1db3Nii7Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-4.000000e+00) * jm[5] * jm[6]) + ((-3.000000e+00) * jm[3] * jm[7]) + ((3.000000e+00) * jm[5] * jm[7]) + ((3.000000e+00) * jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((4.000000e+00) * jm[3] * jm[8]));
        case 636: // da3db1Nii4Njj7
        case 263: // da1db3Nii7Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 646: // da3db1Nii5Njj7
        case 264: // da1db3Nii7Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[2] * ((jm[3] * jm[4] * jm[6]) + (jm[4] * jm[5] * jm[6]) + ((-1.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[3] * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[4] * jm[8]))) + (jm[1] * (((-1.000000e+00) * (jm[5] * jm[5]) * jm[6]) + ((jm[3] * jm[3]) * jm[7]) + (jm[3] * (((-1.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]) + (jm[5] * jm[8]))))));
        case 656: // da3db1Nii6Njj7
        case 265: // da1db3Nii7Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[5] * ((jm[1] * (((2.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]))))) + (jm[2] * (((-1.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[3] * ((jm[5] * jm[7]) + (jm[4] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-2.000000e+00) * jm[8]))))) + (jm[4] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[4] * jm[8]))))));
        case 666: // da3db1Nii7Njj7
        case 266: // da1db3Nii7Njj7
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])) * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))) + (jm[2] * ((jm[3] * ((jm[5] * (((-2.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-2.000000e+00) * jm[8]))))) + ((-1.000000e+00) * (jm[3] * jm[3]) * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[4] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[4] * jm[8]))))) + (jm[0] * (((2.000000e+00) * (jm[4] * jm[4]) * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[5] * (((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]))) + (jm[4] * ((jm[3] * (((-2.000000e+00) * jm[7]) + jm[8])) + (jm[5] * (((-2.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])))))));
        case 676: // da3db1Nii8Njj7
        case 267: // da1db3Nii7Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((((-2.000000e+00) * jm[0] * jm[4]) + (jm[1] * (jm[3] + jm[4]))) * jm[5] * jm[6]) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5] * (((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))))) + ((-1.000000e+00) * (((-1.000000e+00) * jm[0] * jm[4] * (jm[3] + jm[5])) + (jm[1] * ((jm[3] * jm[3]) + (jm[4] * jm[5])))) * jm[8]));
        case 286: // da1db3Nii9Njj7
        case 668: // da3db1Nii7Njj9
        case 686: // da3db1Nii9Njj7
        case 268: // da1db3Nii7Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 696: // da3db1Nii10Njj7
        case 269: // da1db3Nii7Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3] * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))))) + (jm[1] * jm[3] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[5] * jm[7]) + ((2.000000e+00) * (jm[4] * jm[4]) * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[4] * ((jm[3] * (((-2.000000e+00) * jm[7]) + jm[8])) + (jm[5] * (((-2.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])))))));
        case 607: // da3db1Nii1Njj8
        case 270: // da1db3Nii8Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[5] * jm[6]) + ((4.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + (jm[4] * (((-4.000000e+00) * jm[6]) + jm[8])));
        case 617: // da3db1Nii2Njj8
        case 271: // da1db3Nii8Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 627: // da3db1Nii3Njj8
        case 272: // da1db3Nii8Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[3] + ((-1.000000e+00) * jm[4])) * jm[8]));
        case 637: // da3db1Nii4Njj8
        case 273: // da1db3Nii8Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-4.000000e+00) * jm[4] * jm[6]) + ((3.000000e+00) * jm[5] * jm[6]) + ((4.000000e+00) * jm[3] * jm[7]) + ((-3.000000e+00) * jm[5] * jm[7]) + ((-3.000000e+00) * jm[3] * jm[8]) + ((3.000000e+00) * jm[4] * jm[8]));
        case 647: // da3db1Nii5Njj8
        case 274: // da1db3Nii8Njj5
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[5] * (((-1.000000e+00) * (jm[3] + jm[4]) * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]))) + (jm[1] * jm[3] * (jm[3] + ((-1.000000e+00) * jm[4])) * jm[8]) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]))) + (jm[2] * (((jm[4] * jm[4]) * jm[6]) + (jm[3] * jm[4] * (((-1.000000e+00) * jm[7]) + jm[8])) + ((-1.000000e+00) * jm[3] * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[3] * jm[8]))))));
        case 667: // da3db1Nii7Njj8
        case 276: // da1db3Nii8Njj7
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[5] * (((-2.000000e+00) * jm[0] * jm[4] * jm[6]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[3] + jm[4]) * jm[7]))) + ((-1.000000e+00) * (jm[3] + ((-1.000000e+00) * jm[4])) * (((-1.000000e+00) * jm[0] * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * jm[8]) + (jm[2] * (((-1.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[3] * jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[4] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[4] * jm[8]))))));
        case 677: // da3db1Nii8Njj8
        case 277: // da1db3Nii8Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[3] * jm[8]))) + (jm[0] * (((2.000000e+00) * (jm[4] * jm[4]) * jm[6]) + ((-2.000000e+00) * jm[4] * jm[5] * jm[6]) + ((2.000000e+00) * (jm[5] * jm[5]) * jm[6]) + ((-2.000000e+00) * jm[3] * jm[4] * jm[7]) + (jm[3] * jm[5] * jm[7]) + (jm[4] * jm[5] * jm[7]) + ((-2.000000e+00) * (jm[5] * jm[5]) * jm[7]) + ((jm[3] + ((-1.000000e+00) * jm[4])) * (jm[4] + ((-2.000000e+00) * jm[5])) * jm[8]))) + (jm[1] * ((jm[5] * (((2.000000e+00) * jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[8]))))) + ((jm[3] * jm[3]) * (((2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8]))) + (jm[3] * ((jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])) + (jm[5] * (jm[6] + ((-2.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))))));
        case 687: // da3db1Nii9Njj8
        case 278: // da1db3Nii8Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[5] * ((jm[1] * (jm[3] + jm[4] + ((-2.000000e+00) * jm[5])) * jm[6]) + ((-2.000000e+00) * jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[0] * jm[4] * (((-1.000000e+00) * jm[6]) + jm[7])))) + ((-1.000000e+00) * (jm[3] + ((-1.000000e+00) * jm[4])) * (((-1.000000e+00) * jm[0] * jm[4]) + (jm[1] * (jm[3] + ((-2.000000e+00) * jm[5])))) * jm[8]) + (jm[2] * jm[4] * (((-1.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))));
        case 697: // da3db1Nii10Njj8
        case 279: // da1db3Nii8Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[5] * (((-2.000000e+00) * jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])) * jm[6]) + (jm[1] * jm[3] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[3] + jm[4] + ((-2.000000e+00) * jm[5])) * jm[7]))) + ((-1.000000e+00) * (jm[3] + ((-1.000000e+00) * jm[4])) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + ((2.000000e+00) * jm[0] * jm[5])) * jm[8]) + (jm[2] * jm[3] * (((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[3] * jm[8]))));
        case 608: // da3db1Nii1Njj9
        case 280: // da1db3Nii9Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 618: // da3db1Nii2Njj9
        case 281: // da1db3Nii9Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * (((3.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (((3.000000e+00) * jm[6]) + jm[8])));
        case 628: // da3db1Nii3Njj9
        case 282: // da1db3Nii9Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 638: // da3db1Nii4Njj9
        case 283: // da1db3Nii9Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((jm[3] + ((3.000000e+00) * jm[5])) * jm[7]) + ((-1.000000e+00) * jm[4] * (jm[6] + ((3.000000e+00) * jm[8]))));
        case 648: // da3db1Nii5Njj9
        case 284: // da1db3Nii9Njj5
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * (jm[4] * jm[4]) * jm[6]) + ((-2.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[3] * ((jm[5] * jm[7]) + (jm[4] * (((2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8]))))))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])) * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[1] * ((jm[4] * jm[5] * jm[6]) + ((2.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[3] * (((-2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]))))));
        case 658: // da3db1Nii6Njj9
        case 285: // da1db3Nii9Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((2.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[4] * ((jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))))) + ((-1.000000e+00) * jm[3] * ((jm[5] * jm[7]) + (jm[4] * (((2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8]))))))) + (jm[5] * ((jm[0] * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]))))));
        case 678: // da3db1Nii8Njj9
        case 287: // da1db3Nii9Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[5] * (((-1.000000e+00) * jm[0] * jm[4] * jm[6]) + (jm[0] * (jm[3] + jm[4] + ((-2.000000e+00) * jm[5])) * jm[7]) + (jm[1] * ((jm[4] * jm[6]) + ((2.000000e+00) * (((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]))))) + ((-1.000000e+00) * jm[4] * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[0] * jm[4]) + ((-2.000000e+00) * jm[0] * jm[5]) + ((2.000000e+00) * jm[1] * jm[5])) * jm[8]) + (jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]))));
        case 688: // da3db1Nii9Njj9
        case 288: // da1db3Nii9Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * (((2.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[5] * ((jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[3] * (((-2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8]))))) + (jm[4] * ((jm[0] * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[2] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[4] * jm[8]))))));
        case 698: // da3db1Nii10Njj9
        case 289: // da1db3Nii9Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3] * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))))) + (jm[2] * jm[3] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * ((jm[5] * (((-1.000000e+00) * jm[3]) + ((2.000000e+00) * jm[5])) * jm[7]) + (jm[4] * ((jm[3] * jm[7]) + (jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]) + ((-2.000000e+00) * jm[8]))))) + ((jm[4] * jm[4]) * (((-1.000000e+00) * jm[6]) + jm[8])))));
        case 609: // da3db1Nii1Njj10
        case 290: // da1db3Nii10Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 619: // da3db1Nii2Njj10
        case 291: // da1db3Nii10Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 629: // da3db1Nii3Njj10
        case 292: // da1db3Nii10Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-3.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((3.000000e+00) * jm[7]) + jm[8])));
        case 639: // da3db1Nii4Njj10
        case 293: // da1db3Nii10Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[4] * jm[6]) + ((-3.000000e+00) * jm[5] * jm[6]) + (jm[3] * (jm[7] + ((3.000000e+00) * jm[8]))));
        case 659: // da3db1Nii6Njj10
        case 295: // da1db3Nii10Njj6
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * ((jm[4] * (((-2.000000e+00) * jm[4]) + jm[5]) * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[4] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8]))))) + ((jm[3] * jm[3]) * (((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[5] * ((jm[0] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[1] * (((2.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))))));
        case 669: // da3db1Nii7Njj10
        case 296: // da1db3Nii10Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * ((jm[4] * (((-2.000000e+00) * jm[4]) + jm[5]) * jm[6]) + ((-1.000000e+00) * (jm[3] * jm[3]) * jm[7]) + (jm[3] * jm[4] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8]))))) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])) * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[0] * (((2.000000e+00) * (jm[4] * jm[4]) * jm[6]) + (jm[3] * jm[5] * jm[7]) + (jm[4] * (((-2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))))));
        case 679: // da3db1Nii8Njj10
        case 297: // da1db3Nii10Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[1] * (jm[3] + jm[4] + ((-2.000000e+00) * jm[5]))) + ((2.000000e+00) * jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])))) * jm[5] * jm[6]) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[3] * jm[5] * jm[7]) + (jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (((-1.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * (jm[7] + ((-2.000000e+00) * jm[8]))))) + (jm[3] * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + ((2.000000e+00) * jm[0] * jm[5]) + ((-2.000000e+00) * jm[1] * jm[5])) * jm[8]));
        case 689: // da3db1Nii9Njj10
        case 298: // da1db3Nii10Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * ((jm[5] * (((-1.000000e+00) * jm[4]) + ((2.000000e+00) * jm[5])) * jm[6]) + (jm[3] * ((jm[4] * jm[6]) + (jm[5] * (((-1.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8]))))) + ((jm[3] * jm[3]) * (((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[4] * ((jm[0] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[2] * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]))))));
        case 699: // da3db1Nii10Njj10
        case 299: // da1db3Nii10Njj10
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[2] * jm[3] * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]))) + (jm[0] * (((2.000000e+00) * (jm[4] * jm[4]) * jm[6]) + (jm[5] * (((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]))) + (jm[4] * (((-2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))))));
        case 400: // da2db2Nii1Njj1
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]) + (jm[2] * jm[7]) + (jm[0] * jm[8]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]) + (jm[2] * jm[7]) + (jm[0] * jm[8]) + ((-1.000000e+00) * jm[1] * jm[8])));
        case 410: // da2db2Nii2Njj1
        case 401: // da2db2Nii1Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 420: // da2db2Nii3Njj1
        case 402: // da2db2Nii1Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 430: // da2db2Nii4Njj1
        case 403: // da2db2Nii1Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 440: // da2db2Nii5Njj1
        case 404: // da2db2Nii1Njj5
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]) + ((4.000000e+00) * jm[2] * jm[7]) + (jm[1] * (jm[6] + ((-4.000000e+00) * jm[8]))) + (jm[0] * jm[8])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 450: // da2db2Nii6Njj1
        case 405: // da2db2Nii1Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 460: // da2db2Nii7Njj1
        case 406: // da2db2Nii1Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[2] * jm[6]) + (jm[0] * jm[7]) + ((-1.000000e+00) * jm[2] * jm[7]) + ((-4.000000e+00) * jm[0] * jm[8]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 470: // da2db2Nii8Njj1
        case 407: // da2db2Nii1Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + ((-4.000000e+00) * jm[0] * jm[7]) + (jm[2] * jm[7]) + (jm[0] * jm[8]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 480: // da2db2Nii9Njj1
        case 408: // da2db2Nii1Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 490: // da2db2Nii10Njj1
        case 409: // da2db2Nii1Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8]))) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 411: // da2db2Nii2Njj2
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])));
        case 421: // da2db2Nii3Njj2
        case 412: // da2db2Nii2Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 431: // da2db2Nii4Njj2
        case 413: // da2db2Nii2Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((-1.000000e+00) * jm[2] * jm[7]) + (jm[1] * jm[8]));
        case 441: // da2db2Nii5Njj2
        case 414: // da2db2Nii2Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[1] * jm[6]) + ((-3.000000e+00) * jm[2] * jm[6]) + ((-3.000000e+00) * jm[0] * jm[7]) + ((4.000000e+00) * jm[2] * jm[7]) + ((3.000000e+00) * jm[0] * jm[8]) + ((-4.000000e+00) * jm[1] * jm[8])) * (((-1.000000e+00) * jm[2] * jm[7]) + (jm[1] * jm[8]));
        case 451: // da2db2Nii6Njj2
        case 415: // da2db2Nii2Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * ((jm[2] * (((3.000000e+00) * jm[6]) + jm[7])) + ((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 461: // da2db2Nii7Njj2
        case 416: // da2db2Nii2Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])))) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 471: // da2db2Nii8Njj2
        case 417: // da2db2Nii2Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[0] + ((-1.000000e+00) * jm[1])) * jm[8])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 481: // da2db2Nii9Njj2
        case 418: // da2db2Nii2Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (((3.000000e+00) * jm[6]) + jm[8])));
        case 491: // da2db2Nii10Njj2
        case 419: // da2db2Nii2Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[2] * jm[7]) + (jm[1] * jm[8])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 422: // da2db2Nii3Njj3
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])));
        case 432: // da2db2Nii4Njj3
        case 423: // da2db2Nii3Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 442: // da2db2Nii5Njj3
        case 424: // da2db2Nii3Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])))) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 452: // da2db2Nii6Njj3
        case 425: // da2db2Nii3Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[2] * (jm[6] + ((3.000000e+00) * jm[7]))) + ((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[1])) * jm[8]));
        case 462: // da2db2Nii7Njj3
        case 426: // da2db2Nii3Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (((-4.000000e+00) * jm[2] * jm[6]) + ((-3.000000e+00) * jm[0] * jm[7]) + ((3.000000e+00) * jm[2] * jm[7]) + ((3.000000e+00) * jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((4.000000e+00) * jm[0] * jm[8]));
        case 472: // da2db2Nii8Njj3
        case 427: // da2db2Nii3Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 482: // da2db2Nii9Njj3
        case 428: // da2db2Nii3Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])))) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 492: // da2db2Nii10Njj3
        case 429: // da2db2Nii3Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (((-3.000000e+00) * jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((3.000000e+00) * jm[7]) + jm[8])));
        case 433: // da2db2Nii4Njj4
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])));
        case 443: // da2db2Nii5Njj4
        case 434: // da2db2Nii4Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 453: // da2db2Nii6Njj4
        case 435: // da2db2Nii4Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 463: // da2db2Nii7Njj4
        case 436: // da2db2Nii4Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 473: // da2db2Nii8Njj4
        case 437: // da2db2Nii4Njj8
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((4.000000e+00) * jm[1] * jm[6]) + ((-3.000000e+00) * jm[2] * jm[6]) + ((-4.000000e+00) * jm[0] * jm[7]) + ((3.000000e+00) * jm[2] * jm[7]) + ((3.000000e+00) * jm[0] * jm[8]) + ((-3.000000e+00) * jm[1] * jm[8]));
        case 483: // da2db2Nii9Njj4
        case 438: // da2db2Nii4Njj9
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (jm[6] + ((3.000000e+00) * jm[8]))));
        case 493: // da2db2Nii10Njj4
        case 439: // da2db2Nii4Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (((-1.000000e+00) * jm[1] * jm[6]) + ((-3.000000e+00) * jm[2] * jm[6]) + (jm[0] * (jm[7] + ((3.000000e+00) * jm[8]))));
        case 444: // da2db2Nii5Njj5
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-1.000000e+00) * jm[8]))) + ((jm[0] * jm[0]) * ((jm[7] + ((-1.000000e+00) * jm[8])) * (jm[7] + ((-1.000000e+00) * jm[8])))) + ((jm[1] * jm[1]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[8])) * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[2] * (((-2.000000e+00) * (jm[6] * jm[6])) + ((-2.000000e+00) * jm[7] * jm[8]) + (jm[6] * (jm[7] + jm[8])))))));
        case 454: // da2db2Nii6Njj5
        case 479: // da2db2Nii8Njj10
        case 497: // da2db2Nii10Njj8
        case 445: // da2db2Nii5Njj6
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[2] * jm[2]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((-1.000000e+00) * (jm[1] * jm[1]) * jm[6] * jm[8]) + ((2.000000e+00) * (jm[0] * jm[0]) * jm[8] * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[2] * jm[6] * (((-2.000000e+00) * jm[6]) + jm[7] + ((2.000000e+00) * jm[8]))))));
        case 464: // da2db2Nii7Njj5
        case 489: // da2db2Nii9Njj10
        case 498: // da2db2Nii10Njj9
        case 446: // da2db2Nii5Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * (jm[2] * jm[2]) * jm[6]) + ((jm[0] * jm[0]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[0] * jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))) + (jm[1] * ((jm[2] * jm[6] * (((-1.000000e+00) * jm[6]) + jm[7] + ((2.000000e+00) * jm[8]))) + (jm[0] * (((jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 459: // da2db2Nii6Njj10
        case 474: // da2db2Nii8Njj5
        case 495: // da2db2Nii10Njj6
        case 447: // da2db2Nii5Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((-2.000000e+00) * (jm[1] * jm[1]) * jm[6] * jm[8]) + ((jm[0] * jm[0]) * jm[8] * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[0] * jm[2] * ((jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[0] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[2] * jm[6] * (((-1.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])))));
        case 469: // da2db2Nii7Njj10
        case 484: // da2db2Nii9Njj5
        case 496: // da2db2Nii10Njj7
        case 448: // da2db2Nii5Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[1] * jm[1]) * jm[6] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[7] * (((-1.000000e+00) * (jm[2] * jm[2]) * jm[6]) + ((2.000000e+00) * (jm[0] * jm[0]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[8])) * (((-2.000000e+00) * jm[7]) + jm[8])) + (jm[2] * jm[6] * (((-2.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])))));
        case 494: // da2db2Nii10Njj5
        case 449: // da2db2Nii5Njj10
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8]))));
        case 455: // da2db2Nii6Njj6
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + (jm[2] * ((jm[1] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[0] * (((-2.000000e+00) * jm[6]) + jm[7]))) * jm[8]) + (((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1])) * (jm[8] * jm[8])));
        case 465: // da2db2Nii7Njj6
        case 478: // da2db2Nii8Njj9
        case 487: // da2db2Nii9Njj8
        case 456: // da2db2Nii6Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[1] * jm[1]) * (jm[6] + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[1] * ((jm[2] * (jm[6] + ((-2.000000e+00) * jm[7])) * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[0] * jm[8] * (((-1.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))) + (jm[7] * (((2.000000e+00) * (jm[2] * jm[2]) * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((jm[0] * jm[0]) * jm[8]) + ((-1.000000e+00) * jm[0] * jm[2] * (jm[6] + ((-2.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))));
        case 475: // da2db2Nii8Njj6
        case 457: // da2db2Nii6Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])));
        case 485: // da2db2Nii9Njj6
        case 458: // da2db2Nii6Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * jm[8] * (((-1.000000e+00) * jm[6]) + jm[8])) + (jm[7] * (((jm[2] * jm[2]) * (((-1.000000e+00) * jm[6]) + jm[7])) + ((-2.000000e+00) * (jm[0] * jm[0]) * jm[8]) + (jm[0] * jm[2] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[2] * (((-2.000000e+00) * (jm[6] * jm[6])) + ((-2.000000e+00) * jm[7] * jm[8]) + (jm[6] * (jm[7] + jm[8])))))));
        case 466: // da2db2Nii7Njj7
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[7]) + (jm[7] * jm[7]))) + ((jm[1] * jm[1]) * ((jm[6] + ((-1.000000e+00) * jm[8])) * (jm[6] + ((-1.000000e+00) * jm[8])))) + ((-1.000000e+00) * jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])) * ((jm[2] * jm[6]) + ((2.000000e+00) * jm[0] * jm[7]) + ((-2.000000e+00) * jm[2] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8]))) + ((jm[0] * jm[0]) * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[0] * jm[2] * ((jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[7] * (((-2.000000e+00) * jm[7]) + jm[8])))));
        case 476: // da2db2Nii8Njj7
        case 467: // da2db2Nii7Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (jm[6] + ((-1.000000e+00) * jm[8])) * jm[8]) + (jm[7] * (((jm[2] * jm[2]) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * jm[2] * (((-2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8]))) + ((2.000000e+00) * (jm[0] * jm[0]) * jm[8]))) + (jm[1] * ((jm[0] * jm[8] * (((-2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7]) + jm[8])) + (jm[2] * (((2.000000e+00) * (jm[6] * jm[6])) + ((2.000000e+00) * jm[7] * jm[8]) + ((-1.000000e+00) * jm[6] * (jm[7] + jm[8])))))));
        case 486: // da2db2Nii9Njj7
        case 468: // da2db2Nii7Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8]))));
        case 477: // da2db2Nii8Njj8
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[6] + ((-1.000000e+00) * jm[7])) * (jm[6] + ((-1.000000e+00) * jm[7])))) + (jm[0] * jm[2] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((jm[1] * jm[1]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + ((jm[0] * jm[0]) * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[2] * (jm[6] + ((-1.000000e+00) * jm[7])) * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[0] * (((jm[7] + ((-2.000000e+00) * jm[8])) * jm[8]) + (jm[6] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 488: // da2db2Nii9Njj9
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[2]) + (jm[2] * jm[2])) * (jm[7] * jm[7])) + ((jm[1] * jm[1]) * ((jm[6] * jm[6]) + ((-1.000000e+00) * jm[6] * jm[8]) + (jm[8] * jm[8]))) + (jm[1] * jm[7] * ((jm[2] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[0] * (((-2.000000e+00) * jm[6]) + jm[8])))));
        case 499: // da2db2Nii10Njj10
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (jm[6] * jm[6])) + ((jm[2] * jm[2]) * (jm[6] * jm[6])) + (jm[0] * jm[2] * jm[6] * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((jm[0] * jm[0]) * ((jm[7] * jm[7]) + ((-1.000000e+00) * jm[7] * jm[8]) + (jm[8] * jm[8]))) + (jm[1] * jm[6] * (((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-2.000000e+00) * jm[7]) + jm[8])))));
        case 700: // da3db2Nii1Njj1
        case 500: // da2db3Nii1Njj1
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 710: // da3db2Nii2Njj1
        case 501: // da2db3Nii1Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 720: // da3db2Nii3Njj1
        case 502: // da2db3Nii1Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 730: // da3db2Nii4Njj1
        case 503: // da2db3Nii1Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 740: // da3db2Nii5Njj1
        case 504: // da2db3Nii1Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + ((4.000000e+00) * jm[2] * jm[4]) + (jm[1] * (jm[3] + ((-4.000000e+00) * jm[5]))) + (jm[0] * jm[5])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 750: // da3db2Nii6Njj1
        case 505: // da2db3Nii1Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 760: // da3db2Nii7Njj1
        case 506: // da2db3Nii1Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[2] * jm[3]) + (jm[0] * jm[4]) + ((-1.000000e+00) * jm[2] * jm[4]) + ((-4.000000e+00) * jm[0] * jm[5]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 770: // da3db2Nii8Njj1
        case 507: // da2db3Nii1Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + ((-4.000000e+00) * jm[0] * jm[4]) + (jm[2] * jm[4]) + (jm[0] * jm[5]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 780: // da3db2Nii9Njj1
        case 508: // da2db3Nii1Njj9
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 790: // da3db2Nii10Njj1
        case 509: // da2db3Nii1Njj10
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 701: // da3db2Nii1Njj2
        case 510: // da2db3Nii2Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 711: // da3db2Nii2Njj2
        case 511: // da2db3Nii2Njj2
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 721: // da3db2Nii3Njj2
        case 512: // da2db3Nii2Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 731: // da3db2Nii4Njj2
        case 513: // da2db3Nii2Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[2] * jm[7]) + (jm[1] * jm[8]));
        case 741: // da3db2Nii5Njj2
        case 514: // da2db3Nii2Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[1] * jm[3]) + ((-3.000000e+00) * jm[2] * jm[3]) + ((-3.000000e+00) * jm[0] * jm[4]) + ((4.000000e+00) * jm[2] * jm[4]) + ((3.000000e+00) * jm[0] * jm[5]) + ((-4.000000e+00) * jm[1] * jm[5])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 751: // da3db2Nii6Njj2
        case 515: // da2db3Nii2Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((3.000000e+00) * jm[3]) + jm[4])) + ((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 761: // da3db2Nii7Njj2
        case 516: // da2db3Nii2Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 771: // da3db2Nii8Njj2
        case 517: // da2db3Nii2Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 781: // da3db2Nii9Njj2
        case 518: // da2db3Nii2Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (((3.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[2] * jm[7]) + (jm[1] * jm[8]));
        case 791: // da3db2Nii10Njj2
        case 519: // da2db3Nii2Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]));
        case 702: // da3db2Nii1Njj3
        case 520: // da2db3Nii3Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 712: // da3db2Nii2Njj3
        case 521: // da2db3Nii3Njj2
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 722: // da3db2Nii3Njj3
        case 522: // da2db3Nii3Njj3
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 732: // da3db2Nii4Njj3
        case 523: // da2db3Nii3Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * jm[8]));
        case 742: // da3db2Nii5Njj3
        case 524: // da2db3Nii3Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 752: // da3db2Nii6Njj3
        case 525: // da2db3Nii3Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((3.000000e+00) * jm[4]))) + ((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[1])) * jm[5])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 762: // da3db2Nii7Njj3
        case 526: // da2db3Nii3Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-3.000000e+00) * jm[1] * jm[3]) + ((4.000000e+00) * jm[2] * jm[3]) + ((3.000000e+00) * jm[0] * jm[4]) + ((-3.000000e+00) * jm[2] * jm[4]) + ((-4.000000e+00) * jm[0] * jm[5]) + ((3.000000e+00) * jm[1] * jm[5])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 772: // da3db2Nii8Njj3
        case 527: // da2db3Nii3Njj8
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 782: // da3db2Nii9Njj3
        case 528: // da2db3Nii3Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * (((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * jm[8]));
        case 792: // da3db2Nii10Njj3
        case 529: // da2db3Nii3Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * (((3.000000e+00) * jm[4]) + jm[5]))) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8]));
        case 703: // da3db2Nii1Njj4
        case 530: // da2db3Nii4Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]));
        case 713: // da3db2Nii2Njj4
        case 531: // da2db3Nii4Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]));
        case 723: // da3db2Nii3Njj4
        case 532: // da2db3Nii4Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[0] * jm[7]));
        case 733: // da3db2Nii4Njj4
        case 533: // da2db3Nii4Njj4
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[0] * jm[4])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]));
        case 743: // da3db2Nii5Njj4
        case 534: // da2db3Nii4Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]));
        case 753: // da3db2Nii6Njj4
        case 535: // da2db3Nii4Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[0] * jm[7]));
        case 763: // da3db2Nii7Njj4
        case 536: // da2db3Nii4Njj7
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[0] * jm[7]));
        case 773: // da3db2Nii8Njj4
        case 537: // da2db3Nii4Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[1] * jm[3]) + ((-3.000000e+00) * jm[2] * jm[3]) + ((-4.000000e+00) * jm[0] * jm[4]) + ((3.000000e+00) * jm[2] * jm[4]) + ((3.000000e+00) * jm[0] * jm[5]) + ((-3.000000e+00) * jm[1] * jm[5])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]));
        case 783: // da3db2Nii9Njj4
        case 538: // da2db3Nii4Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (jm[3] + ((3.000000e+00) * jm[5])))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]));
        case 793: // da3db2Nii10Njj4
        case 539: // da2db3Nii4Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((3.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * (jm[4] + ((3.000000e+00) * jm[5])))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]));
        case 704: // da3db2Nii1Njj5
        case 540: // da2db3Nii5Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * jm[7]) + ((-4.000000e+00) * jm[2] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8]) + ((4.000000e+00) * jm[1] * jm[8]));
        case 714: // da3db2Nii2Njj5
        case 541: // da2db3Nii5Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((3.000000e+00) * jm[1] * jm[6]) + ((-3.000000e+00) * jm[2] * jm[6]) + ((-3.000000e+00) * jm[0] * jm[7]) + ((4.000000e+00) * jm[2] * jm[7]) + ((3.000000e+00) * jm[0] * jm[8]) + ((-4.000000e+00) * jm[1] * jm[8]));
        case 724: // da3db2Nii3Njj5
        case 542: // da2db3Nii5Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 734: // da3db2Nii4Njj5
        case 543: // da2db3Nii5Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 744: // da3db2Nii5Njj5
        case 544: // da2db3Nii5Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[4] * jm[7]))) + ((2.000000e+00) * (jm[0] * jm[0]) * (jm[4] + ((-1.000000e+00) * jm[5])) * (jm[7] + ((-1.000000e+00) * jm[8]))) + ((jm[1] * jm[1]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))) + (jm[0] * jm[2] * ((jm[5] * (((-2.000000e+00) * jm[6]) + jm[7])) + ((2.000000e+00) * jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((2.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[0] * ((jm[5] * (((2.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])) + ((2.000000e+00) * jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[2] * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[3] * (((-4.000000e+00) * jm[6]) + jm[7] + jm[8])))))));
        case 754: // da3db2Nii6Njj5
        case 545: // da2db3Nii5Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[1] * jm[1]) * jm[5] * jm[6]) + ((jm[2] * jm[2]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]))) + (jm[0] * jm[2] * ((jm[5] * (((-2.000000e+00) * jm[6]) + jm[7])) + ((((2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[4])) * (jm[7] + ((-1.000000e+00) * jm[8]))))) + ((2.000000e+00) * (jm[0] * jm[0]) * jm[5] * (((-1.000000e+00) * jm[7]) + jm[8])) + (jm[1] * ((jm[0] * jm[5] * (((2.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[2] * (((jm[4] + jm[5]) * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[8])))))));
        case 764: // da3db2Nii7Njj5
        case 546: // da2db3Nii5Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (((-1.000000e+00) * jm[3]) + jm[5]) * jm[6]) + ((jm[2] * jm[2]) * ((jm[4] * jm[6]) + (jm[3] * jm[7]))) + ((jm[0] * jm[0]) * jm[4] * (((-1.000000e+00) * jm[7]) + jm[8])) + ((-1.000000e+00) * jm[0] * jm[2] * ((jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[2] * (((-1.000000e+00) * (jm[4] + jm[5]) * jm[6]) + (jm[3] * (jm[6] + ((-1.000000e+00) * jm[8]))))) + (jm[0] * ((jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))))));
        case 774: // da3db2Nii8Njj5
        case 547: // da2db3Nii5Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * (jm[3] + ((-1.000000e+00) * jm[4])) * jm[6]) + ((jm[0] * jm[0]) * jm[5] * (((-1.000000e+00) * jm[7]) + jm[8])) + ((-1.000000e+00) * (jm[1] * jm[1]) * ((jm[5] * jm[6]) + (jm[3] * jm[8]))) + (jm[0] * jm[2] * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))) + (jm[1] * ((jm[2] * (((jm[4] + jm[5]) * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[6]) + jm[7])))) + (jm[0] * ((jm[5] * (jm[6] + jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * jm[8]))))));
        case 784: // da3db2Nii9Njj5
        case 548: // da2db3Nii5Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[4] * (((-1.000000e+00) * (jm[2] * jm[2]) * jm[6]) + ((2.000000e+00) * (jm[0] * jm[0]) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[2] * (((jm[4] + jm[5]) * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7])))) + (jm[0] * (((-1.000000e+00) * (((2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[5])) * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])))))));
        case 594: // da2db3Nii10Njj5
        case 749: // da3db2Nii5Njj10
        case 794: // da3db2Nii10Njj5
        case 549: // da2db3Nii5Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 705: // da3db2Nii1Njj6
        case 550: // da2db3Nii6Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 715: // da3db2Nii2Njj6
        case 551: // da2db3Nii6Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * (((3.000000e+00) * jm[6]) + jm[7])) + ((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 725: // da3db2Nii3Njj6
        case 552: // da2db3Nii6Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * (jm[6] + ((3.000000e+00) * jm[7]))) + ((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[1])) * jm[8]));
        case 735: // da3db2Nii4Njj6
        case 553: // da2db3Nii6Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 745: // da3db2Nii5Njj6
        case 554: // da2db3Nii6Njj5
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]))) + ((-1.000000e+00) * (jm[1] * jm[1]) * jm[3] * jm[8]) + ((2.000000e+00) * (jm[0] * jm[0]) * (((-1.000000e+00) * jm[4]) + jm[5]) * jm[8]) + (jm[0] * jm[2] * (((-2.000000e+00) * jm[5] * jm[6]) + (jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * (((2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[3]) + jm[4] + ((-2.000000e+00) * jm[5])) * jm[8]) + (jm[2] * ((jm[5] * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7] + jm[8])))))));
        case 755: // da3db2Nii6Njj6
        case 555: // da2db3Nii6Njj6
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[4] * jm[7]))) + ((2.000000e+00) * ((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1])) * jm[5] * jm[8]) + (jm[2] * ((jm[1] * ((jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + (jm[3] * jm[8]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * (((-2.000000e+00) * jm[5] * jm[6]) + (jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))))));
        case 765: // da3db2Nii7Njj6
        case 556: // da2db3Nii6Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[2] * (jm[3] + ((-1.000000e+00) * jm[5])) * (jm[6] + ((-2.000000e+00) * jm[7]))) + ((-1.000000e+00) * (jm[2] * jm[2]) * ((jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[3] * jm[7]))) + ((-1.000000e+00) * (jm[0] * jm[0]) * jm[4] * jm[8]) + (jm[1] * ((jm[2] * (jm[3] + ((-2.000000e+00) * jm[4]))) + (jm[0] * (jm[3] + ((2.000000e+00) * jm[4]) + ((-2.000000e+00) * jm[5])))) * jm[8]) + ((-2.000000e+00) * (jm[1] * jm[1]) * (jm[3] + ((-1.000000e+00) * jm[5])) * jm[8]) + (jm[0] * jm[2] * ((jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]) + jm[8])))));
        case 575: // da2db3Nii8Njj6
        case 757: // da3db2Nii6Njj8
        case 775: // da3db2Nii8Njj6
        case 557: // da2db3Nii6Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 785: // da3db2Nii9Njj6
        case 558: // da2db3Nii6Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (jm[3] + ((-1.000000e+00) * jm[5])) * jm[8]) + (jm[4] * (((jm[2] * jm[2]) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * jm[2] * (((-2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8]))) + ((2.000000e+00) * (jm[0] * jm[0]) * jm[8]))) + (jm[1] * ((jm[0] * (((-2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[4]) + jm[5]) * jm[8]) + (jm[2] * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[4] * jm[8]))))));
        case 795: // da3db2Nii10Njj6
        case 559: // da2db3Nii6Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * jm[3] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((2.000000e+00) * (jm[1] * jm[1]) * jm[3] * jm[8]) + ((jm[0] * jm[0]) * (jm[4] + ((-1.000000e+00) * jm[5])) * jm[8]) + (jm[0] * jm[2] * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((2.000000e+00) * jm[4] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[3] * jm[8]))) + (jm[1] * ((jm[2] * jm[3] * (jm[6] + ((-2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[3]) + ((-2.000000e+00) * jm[4]) + jm[5]) * jm[8]))));
        case 706: // da3db2Nii1Njj7
        case 560: // da2db3Nii7Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((4.000000e+00) * jm[2] * jm[6]) + (jm[0] * jm[7]) + ((-1.000000e+00) * jm[2] * jm[7]) + ((-4.000000e+00) * jm[0] * jm[8]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 716: // da3db2Nii2Njj7
        case 561: // da2db3Nii7Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 726: // da3db2Nii3Njj7
        case 562: // da2db3Nii7Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-3.000000e+00) * jm[1] * jm[6]) + ((4.000000e+00) * jm[2] * jm[6]) + ((3.000000e+00) * jm[0] * jm[7]) + ((-3.000000e+00) * jm[2] * jm[7]) + ((-4.000000e+00) * jm[0] * jm[8]) + ((3.000000e+00) * jm[1] * jm[8]));
        case 736: // da3db2Nii4Njj7
        case 563: // da2db3Nii7Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 746: // da3db2Nii5Njj7
        case 564: // da2db3Nii7Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[0] * jm[0]) * (jm[4] + ((-1.000000e+00) * jm[5])) * jm[7]) + ((-1.000000e+00) * (jm[2] * jm[2]) * ((jm[4] * jm[6]) + (jm[3] * jm[7]))) + ((jm[1] * jm[1]) * jm[3] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * jm[2] * (((jm[3] + ((-1.000000e+00) * jm[4]) + jm[5]) * jm[7]) + (jm[4] * jm[8]))) + (jm[1] * ((jm[0] * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[8]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])))) + (jm[2] * ((jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[6]) + jm[7] + jm[8])))))));
        case 756: // da3db2Nii6Njj7
        case 565: // da2db3Nii7Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[0] * jm[0]) * jm[5] * jm[7]) + ((-1.000000e+00) * (jm[2] * jm[2]) * ((jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]))) + (jm[3] * jm[7]))) + ((2.000000e+00) * (jm[1] * jm[1]) * jm[5] * (((-1.000000e+00) * jm[6]) + jm[8])) + (jm[0] * jm[2] * (((jm[3] + ((-2.000000e+00) * jm[4]) + jm[5]) * jm[7]) + (jm[4] * jm[8]))) + (jm[1] * ((jm[0] * jm[5] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-2.000000e+00) * jm[8]))) + (jm[2] * (((2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]) + (jm[3] * (((-1.000000e+00) * jm[6]) + jm[8])))))));
        case 766: // da3db2Nii7Njj7
        case 566: // da2db3Nii7Njj7
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[4] * jm[7]))) + ((2.000000e+00) * (jm[1] * jm[1]) * (jm[3] + ((-1.000000e+00) * jm[5])) * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((jm[0] * jm[0]) * (((2.000000e+00) * jm[4] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))) + (jm[0] * jm[2] * ((jm[5] * (((-2.000000e+00) * jm[6]) + jm[7])) + (jm[3] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[4] * (jm[6] + ((-4.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((jm[3] + ((2.000000e+00) * jm[4]) + ((-2.000000e+00) * jm[5])) * jm[8]))) + (jm[2] * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[7]))) + ((2.000000e+00) * jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])))))));
        case 776: // da3db2Nii8Njj7
        case 567: // da2db3Nii7Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * (jm[3] + ((-1.000000e+00) * jm[4])) * jm[7]) + ((-1.000000e+00) * jm[0] * jm[2] * ((jm[4] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((jm[3] + jm[5]) * jm[7]))) + ((jm[1] * jm[1]) * jm[5] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((jm[0] * jm[0]) * ((jm[5] * jm[7]) + (jm[4] * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[0] * ((jm[5] * (jm[6] + jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[3] * jm[8]))) + (jm[2] * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))))));
        case 586: // da2db3Nii9Njj7
        case 768: // da3db2Nii7Njj9
        case 786: // da3db2Nii9Njj7
        case 568: // da2db3Nii7Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 796: // da3db2Nii10Njj7
        case 569: // da2db3Nii7Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[2] * jm[2]) * jm[3] * jm[7]) + (jm[0] * jm[2] * ((jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]))) + ((jm[3] + jm[5]) * jm[7]))) + ((2.000000e+00) * (jm[1] * jm[1]) * jm[3] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((jm[0] * jm[0]) * (((2.000000e+00) * jm[4] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[1] * ((jm[2] * jm[3] * (((-2.000000e+00) * jm[6]) + ((2.000000e+00) * jm[7]) + jm[8])) + (jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + ((jm[3] + ((2.000000e+00) * jm[4]) + ((-1.000000e+00) * jm[5])) * jm[8]))))));
        case 707: // da3db2Nii1Njj8
        case 570: // da2db3Nii8Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * jm[6]) + ((4.000000e+00) * jm[0] * jm[7]) + ((-1.000000e+00) * jm[2] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8]) + (jm[1] * (((-4.000000e+00) * jm[6]) + jm[8])));
        case 717: // da3db2Nii2Njj8
        case 571: // da2db3Nii8Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 727: // da3db2Nii3Njj8
        case 572: // da2db3Nii8Njj3
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8]));
        case 737: // da3db2Nii4Njj8
        case 573: // da2db3Nii8Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((4.000000e+00) * jm[1] * jm[6]) + ((-3.000000e+00) * jm[2] * jm[6]) + ((-4.000000e+00) * jm[0] * jm[7]) + ((3.000000e+00) * jm[2] * jm[7]) + ((3.000000e+00) * jm[0] * jm[8]) + ((-3.000000e+00) * jm[1] * jm[8]));
        case 747: // da3db2Nii5Njj8
        case 574: // da2db3Nii8Njj5
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * jm[3] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((jm[0] * jm[0]) * (((-1.000000e+00) * jm[4]) + jm[5]) * jm[8]) + (jm[0] * jm[2] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[4] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]))) + ((-1.000000e+00) * (jm[1] * jm[1]) * ((jm[5] * jm[6]) + (jm[3] * jm[8]))) + (jm[1] * ((jm[0] * ((jm[5] * (jm[7] + ((-1.000000e+00) * jm[8]))) + ((jm[3] + jm[4]) * jm[8]))) + (jm[2] * ((jm[4] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[6]) + jm[7] + jm[8])))))));
        case 767: // da3db2Nii7Njj8
        case 576: // da2db3Nii8Njj7
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * jm[4] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((jm[1] * jm[1]) * (((-1.000000e+00) * jm[3]) + jm[5]) * jm[8]) + ((-1.000000e+00) * (jm[0] * jm[0]) * ((jm[5] * jm[7]) + (jm[4] * jm[8]))) + (jm[0] * jm[2] * ((jm[3] * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[2] * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7])) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[0] * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((jm[3] + jm[4]) * jm[8]))))));
        case 777: // da3db2Nii8Njj8
        case 577: // da2db3Nii8Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[2] * jm[2]) * (jm[3] + ((-1.000000e+00) * jm[4])) * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((jm[1] * jm[1]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))) + ((jm[0] * jm[0]) * (((2.000000e+00) * jm[4] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))) + (jm[0] * jm[2] * (((-2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))) + (jm[1] * ((jm[0] * ((jm[5] * (jm[6] + jm[7] + ((-4.000000e+00) * jm[8]))) + (jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])) + (jm[3] * (((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[2] * (((2.000000e+00) * jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7] + ((2.000000e+00) * jm[8]))))))));
        case 787: // da3db2Nii9Njj8
        case 578: // da2db3Nii8Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[1] * jm[1]) * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[3] * jm[8]))) + (jm[4] * (((2.000000e+00) * (jm[2] * jm[2]) * (((-1.000000e+00) * jm[6]) + jm[7])) + ((-1.000000e+00) * (jm[0] * jm[0]) * jm[8]) + (jm[0] * jm[2] * (jm[6] + ((-2.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))) + (jm[1] * ((jm[2] * ((jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + (jm[3] * (((-1.000000e+00) * jm[6]) + jm[7])) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * ((jm[5] * (jm[7] + ((-2.000000e+00) * jm[8]))) + ((jm[3] + jm[4]) * jm[8]))))));
        case 797: // da3db2Nii10Njj8
        case 579: // da2db3Nii8Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[2] * jm[2]) * jm[3] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * jm[2] * (((jm[4] + ((-2.000000e+00) * jm[5])) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-2.000000e+00) * jm[8]))))) + ((-1.000000e+00) * (jm[1] * jm[1]) * jm[3] * jm[8]) + ((-1.000000e+00) * (jm[0] * jm[0]) * ((jm[5] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[4] * jm[8]))) + (jm[1] * ((jm[2] * jm[3] * (((-2.000000e+00) * jm[6]) + jm[7] + ((2.000000e+00) * jm[8]))) + (jm[0] * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[8]))) + ((jm[3] + jm[4]) * jm[8]))))));
        case 708: // da3db2Nii1Njj9
        case 580: // da2db3Nii9Njj1
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))));
        case 718: // da3db2Nii2Njj9
        case 581: // da2db3Nii9Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((((3.000000e+00) * jm[0]) + jm[2]) * jm[7]) + ((-1.000000e+00) * jm[1] * (((3.000000e+00) * jm[6]) + jm[8])));
        case 728: // da3db2Nii3Njj9
        case 582: // da2db3Nii9Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 738: // da3db2Nii4Njj9
        case 583: // da2db3Nii9Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (jm[6] + ((3.000000e+00) * jm[8]))));
        case 748: // da3db2Nii5Njj9
        case 584: // da2db3Nii9Njj5
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * (jm[2] * jm[2]) * jm[3]) + ((2.000000e+00) * (jm[0] * jm[0]) * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[3]) + ((-2.000000e+00) * jm[4]) + jm[5]))) * jm[7]) + ((jm[1] * jm[1]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]))) + (jm[1] * ((jm[0] * (((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * (((2.000000e+00) * jm[6]) + jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-2.000000e+00) * jm[6]) + jm[8])))) + (jm[2] * ((jm[4] * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[6]) + jm[7] + jm[8])))))));
        case 758: // da3db2Nii6Njj9
        case 585: // da2db3Nii9Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((((jm[2] * jm[2]) * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * jm[2] * (((-2.000000e+00) * jm[3]) + jm[4] + ((-1.000000e+00) * jm[5]))) + ((2.000000e+00) * (jm[0] * jm[0]) * jm[5])) * jm[7]) + ((jm[1] * jm[1]) * jm[5] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[1] * ((jm[0] * jm[5] * (((-2.000000e+00) * jm[6]) + ((-1.000000e+00) * jm[7]) + jm[8])) + (jm[2] * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8]))))));
        case 778: // da3db2Nii8Njj9
        case 587: // da2db3Nii9Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((((2.000000e+00) * (jm[2] * jm[2]) * (((-1.000000e+00) * jm[3]) + jm[4])) + ((-1.000000e+00) * (jm[0] * jm[0]) * jm[5]) + (jm[0] * jm[2] * (jm[3] + ((-2.000000e+00) * jm[4]) + ((2.000000e+00) * jm[5])))) * jm[7]) + ((-1.000000e+00) * (jm[1] * jm[1]) * ((jm[5] * (jm[6] + ((-2.000000e+00) * jm[8]))) + (jm[3] * jm[8]))) + (jm[1] * ((jm[0] * ((jm[5] * (jm[6] + jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[4] * jm[8]))) + (jm[2] * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]) + (jm[3] * (((-1.000000e+00) * jm[6]) + jm[7] + ((2.000000e+00) * jm[8]))))))));
        case 788: // da3db2Nii9Njj9
        case 588: // da2db3Nii9Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * ((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[2]) + (jm[2] * jm[2])) * jm[4] * jm[7]) + ((jm[1] * jm[1]) * (((2.000000e+00) * jm[3] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))) + (jm[1] * ((jm[2] * ((jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8]))) + (jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[4] * jm[8]))))));
        case 798: // da3db2Nii10Njj9
        case 589: // da2db3Nii9Njj10
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((((2.000000e+00) * (jm[2] * jm[2]) * jm[3]) + (jm[0] * jm[2] * (((-1.000000e+00) * jm[3]) + jm[4] + ((-2.000000e+00) * jm[5]))) + ((jm[0] * jm[0]) * (((-1.000000e+00) * jm[4]) + jm[5]))) * jm[7]) + ((jm[1] * jm[1]) * jm[3] * (((-1.000000e+00) * jm[6]) + jm[8])) + (jm[1] * ((jm[2] * jm[3] * (jm[6] + ((-1.000000e+00) * jm[7]) + ((-2.000000e+00) * jm[8]))) + (jm[0] * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))))));
        case 709: // da3db2Nii1Njj10
        case 590: // da2db3Nii10Njj1
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))));
        case 719: // da3db2Nii2Njj10
        case 591: // da2db3Nii10Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 729: // da3db2Nii3Njj10
        case 592: // da2db3Nii10Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((3.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * (((3.000000e+00) * jm[7]) + jm[8])));
        case 739: // da3db2Nii4Njj10
        case 593: // da2db3Nii10Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[1] * jm[6]) + ((3.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * (jm[7] + ((3.000000e+00) * jm[8]))));
        case 759: // da3db2Nii6Njj10
        case 595: // da2db3Nii10Njj6
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * (((-1.000000e+00) * jm[3]) + jm[4]) * jm[6]) + ((2.000000e+00) * (jm[1] * jm[1]) * jm[5] * jm[6]) + ((jm[0] * jm[0]) * jm[5] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[0] * jm[2] * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[4] * jm[7]) + (jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[1] * ((jm[2] * (jm[3] + ((-2.000000e+00) * jm[4]) + ((-1.000000e+00) * jm[5])) * jm[6]) + (jm[0] * jm[5] * (((-1.000000e+00) * jm[6]) + ((-2.000000e+00) * jm[7]) + jm[8])))));
        case 769: // da3db2Nii7Njj10
        case 596: // da2db3Nii10Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * (jm[2] * jm[2]) * jm[4] * jm[6]) + ((2.000000e+00) * (jm[1] * jm[1]) * (jm[3] + ((-1.000000e+00) * jm[5])) * jm[6]) + ((jm[0] * jm[0]) * (((2.000000e+00) * jm[4] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + (jm[0] * jm[2] * ((jm[3] * jm[7]) + (jm[4] * (jm[6] + ((-2.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[2] * (((-2.000000e+00) * jm[3]) + ((2.000000e+00) * jm[4]) + jm[5]) * jm[6]) + (jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * (jm[6] + ((2.000000e+00) * jm[7]) + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-2.000000e+00) * jm[7]) + jm[8])))))));
        case 779: // da3db2Nii8Njj10
        case 597: // da2db3Nii10Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[2] * jm[2]) * (jm[3] + ((-1.000000e+00) * jm[4])) * jm[6]) + ((-1.000000e+00) * (jm[1] * jm[1]) * jm[5] * jm[6]) + ((-1.000000e+00) * (jm[0] * jm[0]) * ((jm[5] * (jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[4] * jm[8]))) + (jm[0] * jm[2] * (((-2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[7]) + ((2.000000e+00) * jm[8]))))) + (jm[1] * ((jm[2] * (((-2.000000e+00) * jm[3]) + jm[4] + ((2.000000e+00) * jm[5])) * jm[6]) + (jm[0] * ((jm[5] * (jm[6] + jm[7] + ((-2.000000e+00) * jm[8]))) + (jm[3] * jm[8]))))));
        case 789: // da3db2Nii9Njj10
        case 598: // da2db3Nii10Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (((-1.000000e+00) * jm[3]) + jm[5]) * jm[6]) + (jm[4] * (((2.000000e+00) * (jm[2] * jm[2]) * jm[6]) + (jm[0] * jm[2] * (((-1.000000e+00) * jm[6]) + jm[7] + ((-2.000000e+00) * jm[8]))) + ((jm[0] * jm[0]) * (((-1.000000e+00) * jm[7]) + jm[8])))) + (jm[1] * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]) + ((-2.000000e+00) * jm[5])) * jm[6]) + (jm[0] * ((jm[4] * jm[6]) + (jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))))));
        case 799: // da3db2Nii10Njj10
        case 599: // da2db3Nii10Njj10
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[1] * jm[1]) * jm[3] * jm[6]) + ((2.000000e+00) * (jm[2] * jm[2]) * jm[3] * jm[6]) + (jm[0] * jm[2] * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]))) + ((jm[0] * jm[0]) * (((2.000000e+00) * jm[4] * jm[7]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]) + ((2.000000e+00) * jm[5] * jm[8]))) + (jm[1] * (((-2.000000e+00) * jm[2] * jm[3] * jm[6]) + (jm[0] * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[3] * jm[8]))))));
        case 800: // da3db3Nii1Njj1
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + (jm[2] * jm[4]) + (jm[0] * jm[5]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + (jm[2] * jm[4]) + (jm[0] * jm[5]) + ((-1.000000e+00) * jm[1] * jm[5])));
        case 810: // da3db3Nii2Njj1
        case 801: // da3db3Nii1Njj2
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5])));
        case 820: // da3db3Nii3Njj1
        case 802: // da3db3Nii1Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5])));
        case 830: // da3db3Nii4Njj1
        case 803: // da3db3Nii1Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5])));
        case 840: // da3db3Nii5Njj1
        case 804: // da3db3Nii1Njj5
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + ((4.000000e+00) * jm[2] * jm[4]) + (jm[1] * (jm[3] + ((-4.000000e+00) * jm[5]))) + (jm[0] * jm[5])) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5])));
        case 850: // da3db3Nii6Njj1
        case 805: // da3db3Nii1Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5])));
        case 860: // da3db3Nii7Njj1
        case 806: // da3db3Nii1Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[2] * jm[3]) + (jm[0] * jm[4]) + ((-1.000000e+00) * jm[2] * jm[4]) + ((-4.000000e+00) * jm[0] * jm[5]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5])));
        case 870: // da3db3Nii8Njj1
        case 807: // da3db3Nii1Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((4.000000e+00) * jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + ((-4.000000e+00) * jm[0] * jm[4]) + (jm[2] * jm[4]) + (jm[0] * jm[5]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5])));
        case 880: // da3db3Nii9Njj1
        case 808: // da3db3Nii1Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5])));
        case 890: // da3db3Nii10Njj1
        case 809: // da3db3Nii1Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5])));
        case 811: // da3db3Nii2Njj2
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])));
        case 821: // da3db3Nii3Njj2
        case 812: // da3db3Nii2Njj3
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5]));
        case 831: // da3db3Nii4Njj2
        case 813: // da3db3Nii2Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[2] * jm[4]) + (jm[1] * jm[5]));
        case 841: // da3db3Nii5Njj2
        case 814: // da3db3Nii2Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((3.000000e+00) * jm[1] * jm[3]) + ((-3.000000e+00) * jm[2] * jm[3]) + ((-3.000000e+00) * jm[0] * jm[4]) + ((4.000000e+00) * jm[2] * jm[4]) + ((3.000000e+00) * jm[0] * jm[5]) + ((-4.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * jm[2] * jm[4]) + (jm[1] * jm[5]));
        case 851: // da3db3Nii6Njj2
        case 815: // da3db3Nii2Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((jm[2] * (((3.000000e+00) * jm[3]) + jm[4])) + ((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[1]) * jm[5]));
        case 861: // da3db3Nii7Njj2
        case 816: // da3db3Nii2Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5]));
        case 871: // da3db3Nii8Njj2
        case 817: // da3db3Nii2Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + ((jm[0] + ((-1.000000e+00) * jm[1])) * jm[5])) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5]));
        case 881: // da3db3Nii9Njj2
        case 818: // da3db3Nii2Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (((-1.000000e+00) * (((3.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (((3.000000e+00) * jm[3]) + jm[5])));
        case 891: // da3db3Nii10Njj2
        case 819: // da3db3Nii2Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[2] * jm[4]) + (jm[1] * jm[5])) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5])));
        case 822: // da3db3Nii3Njj3
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])));
        case 832: // da3db3Nii4Njj3
        case 823: // da3db3Nii3Njj4
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5]));
        case 842: // da3db3Nii5Njj3
        case 824: // da3db3Nii3Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5])))) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5]));
        case 852: // da3db3Nii6Njj3
        case 825: // da3db3Nii3Njj6
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * (jm[3] + ((3.000000e+00) * jm[4]))) + ((-1.000000e+00) * (jm[0] + ((3.000000e+00) * jm[1])) * jm[5]));
        case 862: // da3db3Nii7Njj3
        case 826: // da3db3Nii3Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-4.000000e+00) * jm[2] * jm[3]) + ((-3.000000e+00) * jm[0] * jm[4]) + ((3.000000e+00) * jm[2] * jm[4]) + ((3.000000e+00) * jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + ((4.000000e+00) * jm[0] * jm[5]));
        case 872: // da3db3Nii8Njj3
        case 827: // da3db3Nii3Njj8
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5]));
        case 882: // da3db3Nii9Njj3
        case 828: // da3db3Nii3Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5]));
        case 892: // da3db3Nii10Njj3
        case 829: // da3db3Nii3Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (((-3.000000e+00) * jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((3.000000e+00) * jm[4]) + jm[5])));
        case 833: // da3db3Nii4Njj4
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])));
        case 843: // da3db3Nii5Njj4
        case 834: // da3db3Nii4Njj5
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5])));
        case 853: // da3db3Nii6Njj4
        case 835: // da3db3Nii4Njj6
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5]));
        case 863: // da3db3Nii7Njj4
        case 836: // da3db3Nii4Njj7
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))));
        case 873: // da3db3Nii8Njj4
        case 837: // da3db3Nii4Njj8
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((4.000000e+00) * jm[1] * jm[3]) + ((-3.000000e+00) * jm[2] * jm[3]) + ((-4.000000e+00) * jm[0] * jm[4]) + ((3.000000e+00) * jm[2] * jm[4]) + ((3.000000e+00) * jm[0] * jm[5]) + ((-3.000000e+00) * jm[1] * jm[5]));
        case 883: // da3db3Nii9Njj4
        case 838: // da3db3Nii4Njj9
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((jm[0] + ((3.000000e+00) * jm[2])) * jm[4]) + ((-1.000000e+00) * jm[1] * (jm[3] + ((3.000000e+00) * jm[5]))));
        case 893: // da3db3Nii10Njj4
        case 839: // da3db3Nii4Njj10
            return ((1.000000e+00) / (3.000000e+01)) * (1.0 / (ADB)) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (((-1.000000e+00) * jm[1] * jm[3]) + ((-3.000000e+00) * jm[2] * jm[3]) + (jm[0] * (jm[4] + ((3.000000e+00) * jm[5]))));
        case 844: // da3db3Nii5Njj5
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[4]) + (jm[4] * jm[4]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[4])) * (jm[4] + ((-1.000000e+00) * jm[5]))) + ((jm[0] * jm[0]) * ((jm[4] + ((-1.000000e+00) * jm[5])) * (jm[4] + ((-1.000000e+00) * jm[5])))) + ((jm[1] * jm[1]) * ((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[5]) + (jm[5] * jm[5]))) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[5])) * (((-1.000000e+00) * jm[4]) + jm[5])) + (jm[2] * (((-2.000000e+00) * (jm[3] * jm[3])) + ((-2.000000e+00) * jm[4] * jm[5]) + (jm[3] * (jm[4] + jm[5])))))));
        case 854: // da3db3Nii6Njj5
        case 879: // da3db3Nii8Njj10
        case 897: // da3db3Nii10Njj8
        case 845: // da3db3Nii5Njj6
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[2] * jm[2]) * jm[3] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[4])) * (jm[4] + ((-2.000000e+00) * jm[5]))) + ((-1.000000e+00) * (jm[1] * jm[1]) * jm[3] * jm[5]) + ((2.000000e+00) * (jm[0] * jm[0]) * jm[5] * (((-1.000000e+00) * jm[4]) + jm[5])) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[3]) + jm[4] + ((-2.000000e+00) * jm[5])) * jm[5]) + (jm[2] * jm[3] * (((-2.000000e+00) * jm[3]) + jm[4] + ((2.000000e+00) * jm[5]))))));
        case 864: // da3db3Nii7Njj5
        case 889: // da3db3Nii9Njj10
        case 898: // da3db3Nii10Njj9
        case 846: // da3db3Nii5Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * jm[3] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[4] * (((-2.000000e+00) * (jm[2] * jm[2]) * jm[3]) + ((jm[0] * jm[0]) * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[0] * jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]) + ((2.000000e+00) * jm[5]))))) + (jm[1] * ((jm[2] * jm[3] * (((-1.000000e+00) * jm[3]) + jm[4] + ((2.000000e+00) * jm[5]))) + (jm[0] * (((jm[4] + ((-2.000000e+00) * jm[5])) * jm[5]) + (jm[3] * (((-2.000000e+00) * jm[4]) + jm[5])))))));
        case 859: // da3db3Nii6Njj10
        case 874: // da3db3Nii8Njj5
        case 895: // da3db3Nii10Njj6
        case 847: // da3db3Nii5Njj8
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * jm[3] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((-2.000000e+00) * (jm[1] * jm[1]) * jm[3] * jm[5]) + ((jm[0] * jm[0]) * jm[5] * (((-1.000000e+00) * jm[4]) + jm[5])) + (jm[0] * jm[2] * ((jm[3] * (jm[4] + ((-2.000000e+00) * jm[5]))) + (jm[4] * (((-2.000000e+00) * jm[4]) + jm[5])))) + (jm[1] * ((jm[0] * (jm[3] + ((2.000000e+00) * jm[4]) + ((-1.000000e+00) * jm[5])) * jm[5]) + (jm[2] * jm[3] * (((-1.000000e+00) * jm[3]) + ((2.000000e+00) * jm[4]) + jm[5])))));
        case 869: // da3db3Nii7Njj10
        case 884: // da3db3Nii9Njj5
        case 896: // da3db3Nii10Njj7
        case 848: // da3db3Nii5Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[1] * jm[1]) * jm[3] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[4] * (((-1.000000e+00) * (jm[2] * jm[2]) * jm[3]) + ((2.000000e+00) * (jm[0] * jm[0]) * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[0] * jm[2] * (((2.000000e+00) * jm[3]) + ((-2.000000e+00) * jm[4]) + jm[5])))) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[5])) * (((-2.000000e+00) * jm[4]) + jm[5])) + (jm[2] * jm[3] * (((-2.000000e+00) * jm[3]) + ((2.000000e+00) * jm[4]) + jm[5])))));
        case 894: // da3db3Nii10Njj5
        case 849: // da3db3Nii5Njj10
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))));
        case 855: // da3db3Nii6Njj6
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[4]) + (jm[4] * jm[4]))) + (jm[2] * ((jm[1] * (jm[3] + ((-2.000000e+00) * jm[4]))) + (jm[0] * (((-2.000000e+00) * jm[3]) + jm[4]))) * jm[5]) + (((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1])) * (jm[5] * jm[5])));
        case 865: // da3db3Nii7Njj6
        case 878: // da3db3Nii8Njj9
        case 887: // da3db3Nii9Njj8
        case 856: // da3db3Nii6Njj7
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((2.000000e+00) * (jm[1] * jm[1]) * (jm[3] + ((-1.000000e+00) * jm[5])) * jm[5]) + (jm[1] * ((jm[2] * (jm[3] + ((-2.000000e+00) * jm[4])) * (jm[3] + ((-2.000000e+00) * jm[5]))) + (jm[0] * jm[5] * (((-1.000000e+00) * jm[3]) + ((-2.000000e+00) * jm[4]) + ((2.000000e+00) * jm[5]))))) + (jm[4] * (((2.000000e+00) * (jm[2] * jm[2]) * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((jm[0] * jm[0]) * jm[5]) + ((-1.000000e+00) * jm[0] * jm[2] * (jm[3] + ((-2.000000e+00) * jm[4]) + ((2.000000e+00) * jm[5]))))));
        case 875: // da3db3Nii8Njj6
        case 857: // da3db3Nii6Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])));
        case 885: // da3db3Nii9Njj6
        case 858: // da3db3Nii6Njj9
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * jm[5] * (((-1.000000e+00) * jm[3]) + jm[5])) + (jm[4] * (((jm[2] * jm[2]) * (((-1.000000e+00) * jm[3]) + jm[4])) + ((-2.000000e+00) * (jm[0] * jm[0]) * jm[5]) + (jm[0] * jm[2] * (((2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[4]) + jm[5])))) + (jm[1] * ((jm[0] * (((2.000000e+00) * jm[3]) + jm[4] + ((-1.000000e+00) * jm[5])) * jm[5]) + (jm[2] * (((-2.000000e+00) * (jm[3] * jm[3])) + ((-2.000000e+00) * jm[4] * jm[5]) + (jm[3] * (jm[4] + jm[5])))))));
        case 866: // da3db3Nii7Njj7
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[4]) + (jm[4] * jm[4]))) + ((jm[1] * jm[1]) * ((jm[3] + ((-1.000000e+00) * jm[5])) * (jm[3] + ((-1.000000e+00) * jm[5])))) + ((-1.000000e+00) * jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])) * ((jm[2] * jm[3]) + ((2.000000e+00) * jm[0] * jm[4]) + ((-2.000000e+00) * jm[2] * jm[4]) + ((-1.000000e+00) * jm[0] * jm[5]))) + ((jm[0] * jm[0]) * ((jm[4] * jm[4]) + ((-1.000000e+00) * jm[4] * jm[5]) + (jm[5] * jm[5]))) + (jm[0] * jm[2] * ((jm[3] * (jm[4] + ((-2.000000e+00) * jm[5]))) + (jm[4] * (((-2.000000e+00) * jm[4]) + jm[5])))));
        case 876: // da3db3Nii8Njj7
        case 867: // da3db3Nii7Njj8
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (jm[3] + ((-1.000000e+00) * jm[5])) * jm[5]) + (jm[4] * (((jm[2] * jm[2]) * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * jm[2] * (((-2.000000e+00) * jm[3]) + jm[4] + ((-1.000000e+00) * jm[5]))) + ((2.000000e+00) * (jm[0] * jm[0]) * jm[5]))) + (jm[1] * ((jm[0] * jm[5] * (((-2.000000e+00) * jm[3]) + ((-1.000000e+00) * jm[4]) + jm[5])) + (jm[2] * (((2.000000e+00) * (jm[3] * jm[3])) + ((2.000000e+00) * jm[4] * jm[5]) + ((-1.000000e+00) * jm[3] * (jm[4] + jm[5])))))));
        case 886: // da3db3Nii9Njj7
        case 868: // da3db3Nii7Njj9
            return ((-2.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))));
        case 877: // da3db3Nii8Njj8
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[2] * jm[2]) * ((jm[3] + ((-1.000000e+00) * jm[4])) * (jm[3] + ((-1.000000e+00) * jm[4])))) + (jm[0] * jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (jm[4] + ((-2.000000e+00) * jm[5]))) + ((jm[1] * jm[1]) * ((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[5]) + (jm[5] * jm[5]))) + ((jm[0] * jm[0]) * ((jm[4] * jm[4]) + ((-1.000000e+00) * jm[4] * jm[5]) + (jm[5] * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (jm[3] + ((-2.000000e+00) * jm[5]))) + (jm[0] * (((jm[4] + ((-2.000000e+00) * jm[5])) * jm[5]) + (jm[3] * (((-2.000000e+00) * jm[4]) + jm[5])))))));
        case 888: // da3db3Nii9Njj9
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * ((((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[2]) + (jm[2] * jm[2])) * (jm[4] * jm[4])) + ((jm[1] * jm[1]) * ((jm[3] * jm[3]) + ((-1.000000e+00) * jm[3] * jm[5]) + (jm[5] * jm[5]))) + (jm[1] * jm[4] * ((jm[2] * (jm[3] + ((-2.000000e+00) * jm[5]))) + (jm[0] * (((-2.000000e+00) * jm[3]) + jm[5])))));
        case 899: // da3db3Nii10Njj10
            return ((4.000000e+00) / (1.500000e+01)) * (1.0 / (ADB)) * (((jm[1] * jm[1]) * (jm[3] * jm[3])) + ((jm[2] * jm[2]) * (jm[3] * jm[3])) + (jm[0] * jm[2] * jm[3] * (jm[4] + ((-2.000000e+00) * jm[5]))) + ((jm[0] * jm[0]) * ((jm[4] * jm[4]) + ((-1.000000e+00) * jm[4] * jm[5]) + (jm[5] * jm[5]))) + (jm[1] * jm[3] * (((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-2.000000e+00) * jm[4]) + jm[5])))));
		
		default:
			TDKP_GENERAL_EXCEPTION("invalid index");
	}
}
double  Element3DTetrahedron2nd::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {

	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1
	                       + 10 * 10 * diffop;	
	switch(obj_index) {
        case 0: // da1Nii1Njj1
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 2: // da1Nii1Njj3
        case 3: // da1Nii1Njj4
        case 1: // da1Nii1Njj2
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 6: // da1Nii1Njj7
        case 7: // da1Nii1Njj8
        case 4: // da1Nii1Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 8: // da1Nii1Njj9
        case 9: // da1Nii1Njj10
        case 40: // da1Nii5Njj1
        case 60: // da1Nii7Njj1
        case 70: // da1Nii8Njj1
        case 5: // da1Nii1Njj6
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 12: // da1Nii2Njj3
        case 13: // da1Nii2Njj4
        case 10: // da1Nii2Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 11: // da1Nii2Njj2
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 15: // da1Nii2Njj6
        case 18: // da1Nii2Njj9
        case 14: // da1Nii2Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 17: // da1Nii2Njj8
        case 19: // da1Nii2Njj10
        case 41: // da1Nii5Njj2
        case 51: // da1Nii6Njj2
        case 81: // da1Nii9Njj2
        case 16: // da1Nii2Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[5] * jm[7]) + (jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 21: // da1Nii3Njj2
        case 23: // da1Nii3Njj4
        case 20: // da1Nii3Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 22: // da1Nii3Njj3
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 27: // da1Nii3Njj8
        case 28: // da1Nii3Njj9
        case 52: // da1Nii6Njj3
        case 62: // da1Nii7Njj3
        case 92: // da1Nii10Njj3
        case 24: // da1Nii3Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 25: // da1Nii3Njj6
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 29: // da1Nii3Njj10
        case 26: // da1Nii3Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 31: // da1Nii4Njj2
        case 32: // da1Nii4Njj3
        case 30: // da1Nii4Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 33: // da1Nii4Njj4
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 35: // da1Nii4Njj6
        case 36: // da1Nii4Njj7
        case 73: // da1Nii8Njj4
        case 83: // da1Nii9Njj4
        case 93: // da1Nii10Njj4
        case 34: // da1Nii4Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[3] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 38: // da1Nii4Njj9
        case 39: // da1Nii4Njj10
        case 37: // da1Nii4Njj8
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 43: // da1Nii5Njj4
        case 42: // da1Nii5Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 44: // da1Nii5Njj5
            return ((-2.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 48: // da1Nii5Njj9
        case 45: // da1Nii5Njj6
            return ((-1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8]) + ((-1.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 47: // da1Nii5Njj8
        case 46: // da1Nii5Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * jm[5] * (jm[6] + jm[7])) + (jm[4] * (jm[6] + jm[8])) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 49: // da1Nii5Njj10
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 53: // da1Nii6Njj4
        case 50: // da1Nii6Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 58: // da1Nii6Njj9
        case 54: // da1Nii6Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + (jm[4] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 55: // da1Nii6Njj6
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 59: // da1Nii6Njj10
        case 56: // da1Nii6Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[5] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 57: // da1Nii6Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 63: // da1Nii7Njj4
        case 61: // da1Nii7Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 67: // da1Nii7Njj8
        case 64: // da1Nii7Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((-1.000000e+00) * jm[3] * (jm[7] + jm[8])));
        case 69: // da1Nii7Njj10
        case 65: // da1Nii7Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-2.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + ((2.000000e+00) * jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 66: // da1Nii7Njj7
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 68: // da1Nii7Njj9
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 72: // da1Nii8Njj3
        case 71: // da1Nii8Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 76: // da1Nii8Njj7
        case 74: // da1Nii8Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[4] * (jm[6] + jm[8])) + ((-1.000000e+00) * jm[3] * (jm[7] + jm[8])));
        case 75: // da1Nii8Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 77: // da1Nii8Njj8
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 79: // da1Nii8Njj10
        case 78: // da1Nii8Njj9
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-1.000000e+00) * jm[4] * jm[6]) + ((2.000000e+00) * jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 82: // da1Nii9Njj3
        case 80: // da1Nii9Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 85: // da1Nii9Njj6
        case 84: // da1Nii9Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[4] * jm[6]) + ((-2.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 86: // da1Nii9Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 89: // da1Nii9Njj10
        case 87: // da1Nii9Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[5] * jm[7]) + ((-2.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 88: // da1Nii9Njj9
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 91: // da1Nii10Njj2
        case 90: // da1Nii10Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 94: // da1Nii10Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 96: // da1Nii10Njj7
        case 95: // da1Nii10Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((2.000000e+00) * jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-2.000000e+00) * jm[7]) + jm[8])));
        case 98: // da1Nii10Njj9
        case 97: // da1Nii10Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[4] * jm[6]) + ((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + ((2.000000e+00) * jm[3] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 99: // da1Nii10Njj10
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 100: // da2Nii1Njj1
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 102: // da2Nii1Njj3
        case 103: // da2Nii1Njj4
        case 101: // da2Nii1Njj2
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 106: // da2Nii1Njj7
        case 107: // da2Nii1Njj8
        case 104: // da2Nii1Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 108: // da2Nii1Njj9
        case 109: // da2Nii1Njj10
        case 140: // da2Nii5Njj1
        case 160: // da2Nii7Njj1
        case 170: // da2Nii8Njj1
        case 105: // da2Nii1Njj6
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 112: // da2Nii2Njj3
        case 113: // da2Nii2Njj4
        case 110: // da2Nii2Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 111: // da2Nii2Njj2
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[7]) + (jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 115: // da2Nii2Njj6
        case 118: // da2Nii2Njj9
        case 114: // da2Nii2Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[7]) + (jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 117: // da2Nii2Njj8
        case 119: // da2Nii2Njj10
        case 141: // da2Nii5Njj2
        case 151: // da2Nii6Njj2
        case 181: // da2Nii9Njj2
        case 116: // da2Nii2Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 121: // da2Nii3Njj2
        case 123: // da2Nii3Njj4
        case 120: // da2Nii3Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 122: // da2Nii3Njj3
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 127: // da2Nii3Njj8
        case 128: // da2Nii3Njj9
        case 152: // da2Nii6Njj3
        case 162: // da2Nii7Njj3
        case 192: // da2Nii10Njj3
        case 124: // da2Nii3Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 126: // da2Nii3Njj7
        case 129: // da2Nii3Njj10
        case 125: // da2Nii3Njj6
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 131: // da2Nii4Njj2
        case 132: // da2Nii4Njj3
        case 130: // da2Nii4Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 133: // da2Nii4Njj4
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[0] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 135: // da2Nii4Njj6
        case 136: // da2Nii4Njj7
        case 173: // da2Nii8Njj4
        case 183: // da2Nii9Njj4
        case 193: // da2Nii10Njj4
        case 134: // da2Nii4Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 138: // da2Nii4Njj9
        case 139: // da2Nii4Njj10
        case 137: // da2Nii4Njj8
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[0] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 143: // da2Nii5Njj4
        case 142: // da2Nii5Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 144: // da2Nii5Njj5
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 148: // da2Nii5Njj9
        case 145: // da2Nii5Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[1] * jm[6]) + ((-2.000000e+00) * jm[2] * jm[6]) + ((-2.000000e+00) * jm[0] * jm[7]) + (jm[2] * jm[7]) + ((2.000000e+00) * jm[0] * jm[8]) + ((-1.000000e+00) * jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 147: // da2Nii5Njj8
        case 146: // da2Nii5Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * jm[2] * (jm[6] + jm[7])) + (jm[1] * (jm[6] + jm[8])) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 149: // da2Nii5Njj10
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 153: // da2Nii6Njj4
        case 150: // da2Nii6Njj1
            return ((-1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 158: // da2Nii6Njj9
        case 154: // da2Nii6Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[7]) + ((-2.000000e+00) * jm[0] * jm[8]) + (jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 155: // da2Nii6Njj6
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 159: // da2Nii6Njj10
        case 156: // da2Nii6Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * jm[6]) + ((-2.000000e+00) * jm[2] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8]) + ((2.000000e+00) * jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 157: // da2Nii6Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 163: // da2Nii7Njj4
        case 161: // da2Nii7Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 167: // da2Nii7Njj8
        case 164: // da2Nii7Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (jm[6] + jm[7])) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((-1.000000e+00) * jm[0] * (jm[7] + jm[8])));
        case 169: // da2Nii7Njj10
        case 165: // da2Nii7Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[6]) + ((-2.000000e+00) * jm[0] * jm[7]) + ((2.000000e+00) * jm[2] * jm[7]) + ((2.000000e+00) * jm[1] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 166: // da2Nii7Njj7
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 168: // da2Nii7Njj9
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 172: // da2Nii8Njj3
        case 171: // da2Nii8Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 176: // da2Nii8Njj7
        case 174: // da2Nii8Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[1] * (jm[6] + jm[8])) + ((-1.000000e+00) * jm[0] * (jm[7] + jm[8])));
        case 175: // da2Nii8Njj6
            return ((-1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 177: // da2Nii8Njj8
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 179: // da2Nii8Njj10
        case 178: // da2Nii8Njj9
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-2.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]) + ((2.000000e+00) * jm[2] * jm[7]) + (jm[1] * (jm[6] + ((-2.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 182: // da2Nii9Njj3
        case 180: // da2Nii9Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 185: // da2Nii9Njj6
        case 184: // da2Nii9Njj5
            return ((-1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[1] * jm[6]) + ((-2.000000e+00) * jm[0] * jm[7]) + (jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 186: // da2Nii9Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 189: // da2Nii9Njj10
        case 187: // da2Nii9Njj8
            return ((-1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]) + ((2.000000e+00) * jm[2] * jm[7]) + ((-2.000000e+00) * jm[1] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 188: // da2Nii9Njj9
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 191: // da2Nii10Njj2
        case 190: // da2Nii10Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 194: // da2Nii10Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 196: // da2Nii10Njj7
        case 195: // da2Nii10Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-2.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + ((2.000000e+00) * jm[0] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 198: // da2Nii10Njj9
        case 197: // da2Nii10Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[6]) + ((-2.000000e+00) * jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7]) + ((2.000000e+00) * jm[0] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 199: // da2Nii10Njj10
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 200: // da3Nii1Njj1
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 202: // da3Nii1Njj3
        case 201: // da3Nii1Njj2
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 203: // da3Nii1Njj4
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 206: // da3Nii1Njj7
        case 207: // da3Nii1Njj8
        case 204: // da3Nii1Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 208: // da3Nii1Njj9
        case 260: // da3Nii7Njj1
        case 205: // da3Nii1Njj6
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * (((-1.000000e+00) * jm[3]) + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 240: // da3Nii5Njj1
        case 270: // da3Nii8Njj1
        case 209: // da3Nii1Njj10
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 212: // da3Nii2Njj3
        case 213: // da3Nii2Njj4
        case 210: // da3Nii2Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[4]) + (jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 211: // da3Nii2Njj2
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 215: // da3Nii2Njj6
        case 218: // da3Nii2Njj9
        case 214: // da3Nii2Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 217: // da3Nii2Njj8
        case 219: // da3Nii2Njj10
        case 241: // da3Nii5Njj2
        case 281: // da3Nii9Njj2
        case 216: // da3Nii2Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[4]) + (jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 221: // da3Nii3Njj2
        case 223: // da3Nii3Njj4
        case 220: // da3Nii3Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 222: // da3Nii3Njj3
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 227: // da3Nii3Njj8
        case 228: // da3Nii3Njj9
        case 252: // da3Nii6Njj3
        case 262: // da3Nii7Njj3
        case 292: // da3Nii10Njj3
        case 224: // da3Nii3Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 225: // da3Nii3Njj6
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 229: // da3Nii3Njj10
        case 226: // da3Nii3Njj7
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 231: // da3Nii4Njj2
        case 232: // da3Nii4Njj3
        case 230: // da3Nii4Njj1
            return ((1.000000e+00) / (3.600000e+02)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[0] * jm[4])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 233: // da3Nii4Njj4
            return ((1.000000e+00) / (1.200000e+02)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 235: // da3Nii4Njj6
        case 236: // da3Nii4Njj7
        case 283: // da3Nii9Njj4
        case 293: // da3Nii10Njj4
        case 234: // da3Nii4Njj5
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[0] * jm[4])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 238: // da3Nii4Njj9
        case 239: // da3Nii4Njj10
        case 237: // da3Nii4Njj8
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 243: // da3Nii5Njj4
        case 242: // da3Nii5Njj3
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 244: // da3Nii5Njj5
            return ((-2.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 248: // da3Nii5Njj9
        case 245: // da3Nii5Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[2] * jm[3]) + ((2.000000e+00) * jm[0] * jm[4]) + ((-1.000000e+00) * jm[2] * jm[4]) + ((-2.000000e+00) * jm[0] * jm[5]) + (jm[1] * (((-2.000000e+00) * jm[3]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 247: // da3Nii5Njj8
        case 246: // da3Nii5Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[3] + jm[4])) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + ((-1.000000e+00) * jm[1] * (jm[3] + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 249: // da3Nii5Njj10
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 253: // da3Nii6Njj4
        case 250: // da3Nii6Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 251: // da3Nii6Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 258: // da3Nii6Njj9
        case 254: // da3Nii6Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[4]) + ((-2.000000e+00) * jm[0] * jm[5]) + (jm[1] * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 255: // da3Nii6Njj6
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 259: // da3Nii6Njj10
        case 256: // da3Nii6Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * jm[3]) + ((-2.000000e+00) * jm[2] * jm[4]) + ((-1.000000e+00) * jm[0] * jm[5]) + ((2.000000e+00) * jm[1] * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 257: // da3Nii6Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 261: // da3Nii7Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 263: // da3Nii7Njj4
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 267: // da3Nii7Njj8
        case 264: // da3Nii7Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[3] + jm[4])) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5]))) + ((-1.000000e+00) * jm[0] * (jm[4] + jm[5]))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 269: // da3Nii7Njj10
        case 265: // da3Nii7Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-2.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + ((2.000000e+00) * jm[0] * jm[4]) + ((-2.000000e+00) * jm[2] * jm[4]) + ((-1.000000e+00) * jm[0] * jm[5]) + ((2.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 266: // da3Nii7Njj7
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 268: // da3Nii7Njj9
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 272: // da3Nii8Njj3
        case 271: // da3Nii8Njj2
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 273: // da3Nii8Njj4
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 276: // da3Nii8Njj7
        case 274: // da3Nii8Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[1] * (jm[3] + jm[5])) + ((-1.000000e+00) * jm[0] * (jm[4] + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 275: // da3Nii8Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 277: // da3Nii8Njj8
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 279: // da3Nii8Njj10
        case 278: // da3Nii8Njj9
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[3]) + ((2.000000e+00) * jm[2] * jm[3]) + (jm[0] * jm[4]) + ((-2.000000e+00) * jm[2] * jm[4]) + ((-2.000000e+00) * jm[0] * jm[5]) + ((2.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 282: // da3Nii9Njj3
        case 280: // da3Nii9Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 285: // da3Nii9Njj6
        case 284: // da3Nii9Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[1] * jm[3]) + ((-2.000000e+00) * jm[0] * jm[4]) + (jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 286: // da3Nii9Njj7
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 289: // da3Nii9Njj10
        case 287: // da3Nii9Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + ((2.000000e+00) * jm[2] * jm[4]) + ((-2.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 288: // da3Nii9Njj9
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 291: // da3Nii10Njj2
        case 290: // da3Nii10Njj1
            return ((1.000000e+00) / (9.000000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 294: // da3Nii10Njj5
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 296: // da3Nii10Njj7
        case 295: // da3Nii10Njj6
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * (((2.000000e+00) * jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-2.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 298: // da3Nii10Njj9
        case 297: // da3Nii10Njj8
            return ((1.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-2.000000e+00) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4]) + ((2.000000e+00) * jm[0] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 299: // da3Nii10Njj10
            return ((2.000000e+00) / (4.500000e+01)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
			
		default:
			TDKP_GENERAL_EXCEPTION("invalid index");
	}
	
}
double  Element3DTetrahedron2nd::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {

	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1;
	
	switch(obj_index) {
        case 11: // Nii2Njj2
        case 22: // Nii3Njj3
        case 33: // Nii4Njj4
        case 0: // Nii1Njj1
            return ((1.000000e+00) / (4.200000e+02)) * (ADB);
        case 2: // Nii1Njj3
        case 3: // Nii1Njj4
        case 10: // Nii2Njj1
        case 12: // Nii2Njj3
        case 13: // Nii2Njj4
        case 20: // Nii3Njj1
        case 21: // Nii3Njj2
        case 23: // Nii3Njj4
        case 30: // Nii4Njj1
        case 31: // Nii4Njj2
        case 32: // Nii4Njj3
        case 1: // Nii1Njj2
            return ((1.000000e+00) / (2.520000e+03)) * (ADB);
        case 6: // Nii1Njj7
        case 7: // Nii1Njj8
        case 14: // Nii2Njj5
        case 15: // Nii2Njj6
        case 18: // Nii2Njj9
        case 25: // Nii3Njj6
        case 26: // Nii3Njj7
        case 29: // Nii3Njj10
        case 37: // Nii4Njj8
        case 38: // Nii4Njj9
        case 39: // Nii4Njj10
        case 40: // Nii5Njj1
        case 41: // Nii5Njj2
        case 51: // Nii6Njj2
        case 52: // Nii6Njj3
        case 60: // Nii7Njj1
        case 62: // Nii7Njj3
        case 70: // Nii8Njj1
        case 73: // Nii8Njj4
        case 81: // Nii9Njj2
        case 83: // Nii9Njj4
        case 92: // Nii10Njj3
        case 93: // Nii10Njj4
        case 4: // Nii1Njj5
            return ((-1.000000e+00) / (6.300000e+02)) * (ADB);
        case 8: // Nii1Njj9
        case 9: // Nii1Njj10
        case 16: // Nii2Njj7
        case 17: // Nii2Njj8
        case 19: // Nii2Njj10
        case 24: // Nii3Njj5
        case 27: // Nii3Njj8
        case 28: // Nii3Njj9
        case 34: // Nii4Njj5
        case 35: // Nii4Njj6
        case 36: // Nii4Njj7
        case 42: // Nii5Njj3
        case 43: // Nii5Njj4
        case 50: // Nii6Njj1
        case 53: // Nii6Njj4
        case 61: // Nii7Njj2
        case 63: // Nii7Njj4
        case 71: // Nii8Njj2
        case 72: // Nii8Njj3
        case 80: // Nii9Njj1
        case 82: // Nii9Njj3
        case 90: // Nii10Njj1
        case 91: // Nii10Njj2
        case 5: // Nii1Njj6
            return ((-1.000000e+00) / (4.200000e+02)) * (ADB);
        case 55: // Nii6Njj6
        case 66: // Nii7Njj7
        case 77: // Nii8Njj8
        case 88: // Nii9Njj9
        case 99: // Nii10Njj10
        case 44: // Nii5Njj5
            return ((4.000000e+00) / (3.150000e+02)) * (ADB);
        case 46: // Nii5Njj7
        case 47: // Nii5Njj8
        case 48: // Nii5Njj9
        case 54: // Nii6Njj5
        case 56: // Nii6Njj7
        case 58: // Nii6Njj9
        case 59: // Nii6Njj10
        case 64: // Nii7Njj5
        case 65: // Nii7Njj6
        case 67: // Nii7Njj8
        case 69: // Nii7Njj10
        case 74: // Nii8Njj5
        case 76: // Nii8Njj7
        case 78: // Nii8Njj9
        case 79: // Nii8Njj10
        case 84: // Nii9Njj5
        case 85: // Nii9Njj6
        case 87: // Nii9Njj8
        case 89: // Nii9Njj10
        case 95: // Nii10Njj6
        case 96: // Nii10Njj7
        case 97: // Nii10Njj8
        case 98: // Nii10Njj9
        case 45: // Nii5Njj6
            return ((2.000000e+00) / (3.150000e+02)) * (ADB);
        case 57: // Nii6Njj8
        case 68: // Nii7Njj9
        case 75: // Nii8Njj6
        case 86: // Nii9Njj7
        case 94: // Nii10Njj5
        case 49: // Nii5Njj10
            return ((1.000000e+00) / (3.150000e+02)) * (ADB);
		
		default:
			TDKP_GENERAL_EXCEPTION("invalid index");
	}
	
}


double Element3DTetrahedron2nd::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	
	unsigned int obj_index = elem_shape_func_2 
	                       + 10 * elem_shape_func_1
	                       + 100 * nodal_data_point;
	                       
	switch(obj_index) {
        case 111: // Ncc2Nii2Njj2
        case 222: // Ncc3Nii3Njj3
        case 333: // Ncc4Nii4Njj4
        case 0: // Ncc1Nii1Njj1
            return ((1.000000e+00) / (2.520000e+03)) * (ADB);
        case 2: // Ncc1Nii1Njj3
        case 3: // Ncc1Nii1Njj4
        case 10: // Ncc1Nii2Njj1
        case 11: // Ncc1Nii2Njj2
        case 20: // Ncc1Nii3Njj1
        case 22: // Ncc1Nii3Njj3
        case 30: // Ncc1Nii4Njj1
        case 33: // Ncc1Nii4Njj4
        case 100: // Ncc2Nii1Njj1
        case 101: // Ncc2Nii1Njj2
        case 110: // Ncc2Nii2Njj1
        case 112: // Ncc2Nii2Njj3
        case 113: // Ncc2Nii2Njj4
        case 121: // Ncc2Nii3Njj2
        case 122: // Ncc2Nii3Njj3
        case 131: // Ncc2Nii4Njj2
        case 133: // Ncc2Nii4Njj4
        case 200: // Ncc3Nii1Njj1
        case 202: // Ncc3Nii1Njj3
        case 211: // Ncc3Nii2Njj2
        case 212: // Ncc3Nii2Njj3
        case 220: // Ncc3Nii3Njj1
        case 221: // Ncc3Nii3Njj2
        case 223: // Ncc3Nii3Njj4
        case 232: // Ncc3Nii4Njj3
        case 233: // Ncc3Nii4Njj4
        case 300: // Ncc4Nii1Njj1
        case 303: // Ncc4Nii1Njj4
        case 311: // Ncc4Nii2Njj2
        case 313: // Ncc4Nii2Njj4
        case 322: // Ncc4Nii3Njj3
        case 323: // Ncc4Nii3Njj4
        case 330: // Ncc4Nii4Njj1
        case 331: // Ncc4Nii4Njj2
        case 332: // Ncc4Nii4Njj3
        case 1: // Ncc1Nii1Njj2
            return ((-1.000000e+00) / (7.560000e+03)) * (ADB);
        case 6: // Ncc1Nii1Njj7
        case 7: // Ncc1Nii1Njj8
        case 40: // Ncc1Nii5Njj1
        case 60: // Ncc1Nii7Njj1
        case 70: // Ncc1Nii8Njj1
        case 114: // Ncc2Nii2Njj5
        case 115: // Ncc2Nii2Njj6
        case 118: // Ncc2Nii2Njj9
        case 141: // Ncc2Nii5Njj2
        case 151: // Ncc2Nii6Njj2
        case 181: // Ncc2Nii9Njj2
        case 225: // Ncc3Nii3Njj6
        case 226: // Ncc3Nii3Njj7
        case 229: // Ncc3Nii3Njj10
        case 252: // Ncc3Nii6Njj3
        case 262: // Ncc3Nii7Njj3
        case 292: // Ncc3Nii10Njj3
        case 337: // Ncc4Nii4Njj8
        case 338: // Ncc4Nii4Njj9
        case 339: // Ncc4Nii4Njj10
        case 373: // Ncc4Nii8Njj4
        case 383: // Ncc4Nii9Njj4
        case 393: // Ncc4Nii10Njj4
        case 400: // Ncc5Nii1Njj1
        case 411: // Ncc5Nii2Njj2
        case 511: // Ncc6Nii2Njj2
        case 522: // Ncc6Nii3Njj3
        case 600: // Ncc7Nii1Njj1
        case 622: // Ncc7Nii3Njj3
        case 700: // Ncc8Nii1Njj1
        case 733: // Ncc8Nii4Njj4
        case 811: // Ncc9Nii2Njj2
        case 833: // Ncc9Nii4Njj4
        case 922: // Ncc10Nii3Njj3
        case 933: // Ncc10Nii4Njj4
        case 4: // Ncc1Nii1Njj5
            return ((1.000000e+00) / (1.890000e+03)) * (ADB);
        case 8: // Ncc1Nii1Njj9
        case 9: // Ncc1Nii1Njj10
        case 50: // Ncc1Nii6Njj1
        case 80: // Ncc1Nii9Njj1
        case 90: // Ncc1Nii10Njj1
        case 116: // Ncc2Nii2Njj7
        case 117: // Ncc2Nii2Njj8
        case 119: // Ncc2Nii2Njj10
        case 161: // Ncc2Nii7Njj2
        case 171: // Ncc2Nii8Njj2
        case 191: // Ncc2Nii10Njj2
        case 224: // Ncc3Nii3Njj5
        case 227: // Ncc3Nii3Njj8
        case 228: // Ncc3Nii3Njj9
        case 242: // Ncc3Nii5Njj3
        case 272: // Ncc3Nii8Njj3
        case 282: // Ncc3Nii9Njj3
        case 334: // Ncc4Nii4Njj5
        case 335: // Ncc4Nii4Njj6
        case 336: // Ncc4Nii4Njj7
        case 343: // Ncc4Nii5Njj4
        case 353: // Ncc4Nii6Njj4
        case 363: // Ncc4Nii7Njj4
        case 422: // Ncc5Nii3Njj3
        case 433: // Ncc5Nii4Njj4
        case 500: // Ncc6Nii1Njj1
        case 533: // Ncc6Nii4Njj4
        case 611: // Ncc7Nii2Njj2
        case 633: // Ncc7Nii4Njj4
        case 711: // Ncc8Nii2Njj2
        case 722: // Ncc8Nii3Njj3
        case 800: // Ncc9Nii1Njj1
        case 822: // Ncc9Nii3Njj3
        case 900: // Ncc10Nii1Njj1
        case 911: // Ncc10Nii2Njj2
        case 5: // Ncc1Nii1Njj6
            return ((1.000000e+00) / (3.780000e+03)) * (ADB);
        case 13: // Ncc1Nii2Njj4
        case 21: // Ncc1Nii3Njj2
        case 23: // Ncc1Nii3Njj4
        case 31: // Ncc1Nii4Njj2
        case 32: // Ncc1Nii4Njj3
        case 102: // Ncc2Nii1Njj3
        case 103: // Ncc2Nii1Njj4
        case 120: // Ncc2Nii3Njj1
        case 123: // Ncc2Nii3Njj4
        case 130: // Ncc2Nii4Njj1
        case 132: // Ncc2Nii4Njj3
        case 201: // Ncc3Nii1Njj2
        case 203: // Ncc3Nii1Njj4
        case 210: // Ncc3Nii2Njj1
        case 213: // Ncc3Nii2Njj4
        case 230: // Ncc3Nii4Njj1
        case 231: // Ncc3Nii4Njj2
        case 301: // Ncc4Nii1Njj2
        case 302: // Ncc4Nii1Njj3
        case 310: // Ncc4Nii2Njj1
        case 312: // Ncc4Nii2Njj3
        case 320: // Ncc4Nii3Njj1
        case 321: // Ncc4Nii3Njj2
        case 12: // Ncc1Nii2Njj3
            return ((-1.000000e+00) / (4.536000e+04)) * (ADB);
        case 16: // Ncc1Nii2Njj7
        case 17: // Ncc1Nii2Njj8
        case 18: // Ncc1Nii2Njj9
        case 24: // Ncc1Nii3Njj5
        case 25: // Ncc1Nii3Njj6
        case 27: // Ncc1Nii3Njj8
        case 29: // Ncc1Nii3Njj10
        case 34: // Ncc1Nii4Njj5
        case 36: // Ncc1Nii4Njj7
        case 38: // Ncc1Nii4Njj9
        case 39: // Ncc1Nii4Njj10
        case 42: // Ncc1Nii5Njj3
        case 43: // Ncc1Nii5Njj4
        case 51: // Ncc1Nii6Njj2
        case 52: // Ncc1Nii6Njj3
        case 61: // Ncc1Nii7Njj2
        case 63: // Ncc1Nii7Njj4
        case 71: // Ncc1Nii8Njj2
        case 72: // Ncc1Nii8Njj3
        case 81: // Ncc1Nii9Njj2
        case 83: // Ncc1Nii9Njj4
        case 92: // Ncc1Nii10Njj3
        case 93: // Ncc1Nii10Njj4
        case 105: // Ncc2Nii1Njj6
        case 106: // Ncc2Nii1Njj7
        case 107: // Ncc2Nii1Njj8
        case 108: // Ncc2Nii1Njj9
        case 124: // Ncc2Nii3Njj5
        case 126: // Ncc2Nii3Njj7
        case 128: // Ncc2Nii3Njj9
        case 129: // Ncc2Nii3Njj10
        case 134: // Ncc2Nii4Njj5
        case 135: // Ncc2Nii4Njj6
        case 137: // Ncc2Nii4Njj8
        case 139: // Ncc2Nii4Njj10
        case 142: // Ncc2Nii5Njj3
        case 143: // Ncc2Nii5Njj4
        case 150: // Ncc2Nii6Njj1
        case 153: // Ncc2Nii6Njj4
        case 160: // Ncc2Nii7Njj1
        case 162: // Ncc2Nii7Njj3
        case 170: // Ncc2Nii8Njj1
        case 173: // Ncc2Nii8Njj4
        case 180: // Ncc2Nii9Njj1
        case 182: // Ncc2Nii9Njj3
        case 192: // Ncc2Nii10Njj3
        case 193: // Ncc2Nii10Njj4
        case 204: // Ncc3Nii1Njj5
        case 205: // Ncc3Nii1Njj6
        case 207: // Ncc3Nii1Njj8
        case 209: // Ncc3Nii1Njj10
        case 214: // Ncc3Nii2Njj5
        case 216: // Ncc3Nii2Njj7
        case 218: // Ncc3Nii2Njj9
        case 219: // Ncc3Nii2Njj10
        case 235: // Ncc3Nii4Njj6
        case 236: // Ncc3Nii4Njj7
        case 237: // Ncc3Nii4Njj8
        case 238: // Ncc3Nii4Njj9
        case 240: // Ncc3Nii5Njj1
        case 241: // Ncc3Nii5Njj2
        case 250: // Ncc3Nii6Njj1
        case 253: // Ncc3Nii6Njj4
        case 261: // Ncc3Nii7Njj2
        case 263: // Ncc3Nii7Njj4
        case 270: // Ncc3Nii8Njj1
        case 273: // Ncc3Nii8Njj4
        case 281: // Ncc3Nii9Njj2
        case 283: // Ncc3Nii9Njj4
        case 290: // Ncc3Nii10Njj1
        case 291: // Ncc3Nii10Njj2
        case 304: // Ncc4Nii1Njj5
        case 306: // Ncc4Nii1Njj7
        case 308: // Ncc4Nii1Njj9
        case 309: // Ncc4Nii1Njj10
        case 314: // Ncc4Nii2Njj5
        case 315: // Ncc4Nii2Njj6
        case 317: // Ncc4Nii2Njj8
        case 319: // Ncc4Nii2Njj10
        case 325: // Ncc4Nii3Njj6
        case 326: // Ncc4Nii3Njj7
        case 327: // Ncc4Nii3Njj8
        case 328: // Ncc4Nii3Njj9
        case 340: // Ncc4Nii5Njj1
        case 341: // Ncc4Nii5Njj2
        case 351: // Ncc4Nii6Njj2
        case 352: // Ncc4Nii6Njj3
        case 360: // Ncc4Nii7Njj1
        case 362: // Ncc4Nii7Njj3
        case 371: // Ncc4Nii8Njj2
        case 372: // Ncc4Nii8Njj3
        case 380: // Ncc4Nii9Njj1
        case 382: // Ncc4Nii9Njj3
        case 390: // Ncc4Nii10Njj1
        case 391: // Ncc4Nii10Njj2
        case 402: // Ncc5Nii1Njj3
        case 403: // Ncc5Nii1Njj4
        case 412: // Ncc5Nii2Njj3
        case 413: // Ncc5Nii2Njj4
        case 420: // Ncc5Nii3Njj1
        case 421: // Ncc5Nii3Njj2
        case 430: // Ncc5Nii4Njj1
        case 431: // Ncc5Nii4Njj2
        case 501: // Ncc6Nii1Njj2
        case 502: // Ncc6Nii1Njj3
        case 510: // Ncc6Nii2Njj1
        case 513: // Ncc6Nii2Njj4
        case 520: // Ncc6Nii3Njj1
        case 523: // Ncc6Nii3Njj4
        case 531: // Ncc6Nii4Njj2
        case 532: // Ncc6Nii4Njj3
        case 601: // Ncc7Nii1Njj2
        case 603: // Ncc7Nii1Njj4
        case 610: // Ncc7Nii2Njj1
        case 612: // Ncc7Nii2Njj3
        case 621: // Ncc7Nii3Njj2
        case 623: // Ncc7Nii3Njj4
        case 630: // Ncc7Nii4Njj1
        case 632: // Ncc7Nii4Njj3
        case 701: // Ncc8Nii1Njj2
        case 702: // Ncc8Nii1Njj3
        case 710: // Ncc8Nii2Njj1
        case 713: // Ncc8Nii2Njj4
        case 720: // Ncc8Nii3Njj1
        case 723: // Ncc8Nii3Njj4
        case 731: // Ncc8Nii4Njj2
        case 732: // Ncc8Nii4Njj3
        case 801: // Ncc9Nii1Njj2
        case 803: // Ncc9Nii1Njj4
        case 810: // Ncc9Nii2Njj1
        case 812: // Ncc9Nii2Njj3
        case 821: // Ncc9Nii3Njj2
        case 823: // Ncc9Nii3Njj4
        case 830: // Ncc9Nii4Njj1
        case 832: // Ncc9Nii4Njj3
        case 902: // Ncc10Nii1Njj3
        case 903: // Ncc10Nii1Njj4
        case 912: // Ncc10Nii2Njj3
        case 913: // Ncc10Nii2Njj4
        case 920: // Ncc10Nii3Njj1
        case 921: // Ncc10Nii3Njj2
        case 930: // Ncc10Nii4Njj1
        case 931: // Ncc10Nii4Njj2
        case 15: // Ncc1Nii2Njj6
            return ((1.000000e+00) / (7.560000e+03)) * (ADB);
        case 28: // Ncc1Nii3Njj9
        case 35: // Ncc1Nii4Njj6
        case 53: // Ncc1Nii6Njj4
        case 82: // Ncc1Nii9Njj3
        case 91: // Ncc1Nii10Njj2
        case 109: // Ncc2Nii1Njj10
        case 127: // Ncc2Nii3Njj8
        case 136: // Ncc2Nii4Njj7
        case 163: // Ncc2Nii7Njj4
        case 172: // Ncc2Nii8Njj3
        case 190: // Ncc2Nii10Njj1
        case 208: // Ncc3Nii1Njj9
        case 217: // Ncc3Nii2Njj8
        case 234: // Ncc3Nii4Njj5
        case 243: // Ncc3Nii5Njj4
        case 271: // Ncc3Nii8Njj2
        case 280: // Ncc3Nii9Njj1
        case 305: // Ncc4Nii1Njj6
        case 316: // Ncc4Nii2Njj7
        case 324: // Ncc4Nii3Njj5
        case 342: // Ncc4Nii5Njj3
        case 350: // Ncc4Nii6Njj1
        case 361: // Ncc4Nii7Njj2
        case 423: // Ncc5Nii3Njj4
        case 432: // Ncc5Nii4Njj3
        case 503: // Ncc6Nii1Njj4
        case 530: // Ncc6Nii4Njj1
        case 613: // Ncc7Nii2Njj4
        case 631: // Ncc7Nii4Njj2
        case 712: // Ncc8Nii2Njj3
        case 721: // Ncc8Nii3Njj2
        case 802: // Ncc9Nii1Njj3
        case 820: // Ncc9Nii3Njj1
        case 901: // Ncc10Nii1Njj2
        case 910: // Ncc10Nii2Njj1
        case 19: // Ncc1Nii2Njj10
            return ((1.000000e+00) / (5.670000e+03)) * (ADB);
        case 45: // Ncc1Nii5Njj6
        case 48: // Ncc1Nii5Njj9
        case 54: // Ncc1Nii6Njj5
        case 56: // Ncc1Nii6Njj7
        case 65: // Ncc1Nii7Njj6
        case 66: // Ncc1Nii7Njj7
        case 69: // Ncc1Nii7Njj10
        case 77: // Ncc1Nii8Njj8
        case 78: // Ncc1Nii8Njj9
        case 79: // Ncc1Nii8Njj10
        case 84: // Ncc1Nii9Njj5
        case 87: // Ncc1Nii9Njj8
        case 96: // Ncc1Nii10Njj7
        case 97: // Ncc1Nii10Njj8
        case 144: // Ncc2Nii5Njj5
        case 146: // Ncc2Nii5Njj7
        case 147: // Ncc2Nii5Njj8
        case 155: // Ncc2Nii6Njj6
        case 156: // Ncc2Nii6Njj7
        case 159: // Ncc2Nii6Njj10
        case 164: // Ncc2Nii7Njj5
        case 165: // Ncc2Nii7Njj6
        case 174: // Ncc2Nii8Njj5
        case 178: // Ncc2Nii8Njj9
        case 187: // Ncc2Nii9Njj8
        case 188: // Ncc2Nii9Njj9
        case 189: // Ncc2Nii9Njj10
        case 195: // Ncc2Nii10Njj6
        case 198: // Ncc2Nii10Njj9
        case 245: // Ncc3Nii5Njj6
        case 246: // Ncc3Nii5Njj7
        case 254: // Ncc3Nii6Njj5
        case 255: // Ncc3Nii6Njj6
        case 258: // Ncc3Nii6Njj9
        case 264: // Ncc3Nii7Njj5
        case 266: // Ncc3Nii7Njj7
        case 267: // Ncc3Nii7Njj8
        case 276: // Ncc3Nii8Njj7
        case 279: // Ncc3Nii8Njj10
        case 285: // Ncc3Nii9Njj6
        case 289: // Ncc3Nii9Njj10
        case 297: // Ncc3Nii10Njj8
        case 298: // Ncc3Nii10Njj9
        case 299: // Ncc3Nii10Njj10
        case 347: // Ncc4Nii5Njj8
        case 348: // Ncc4Nii5Njj9
        case 358: // Ncc4Nii6Njj9
        case 359: // Ncc4Nii6Njj10
        case 367: // Ncc4Nii7Njj8
        case 369: // Ncc4Nii7Njj10
        case 374: // Ncc4Nii8Njj5
        case 376: // Ncc4Nii8Njj7
        case 377: // Ncc4Nii8Njj8
        case 384: // Ncc4Nii9Njj5
        case 385: // Ncc4Nii9Njj6
        case 388: // Ncc4Nii9Njj9
        case 395: // Ncc4Nii10Njj6
        case 396: // Ncc4Nii10Njj7
        case 399: // Ncc4Nii10Njj10
        case 404: // Ncc5Nii1Njj5
        case 405: // Ncc5Nii1Njj6
        case 408: // Ncc5Nii1Njj9
        case 414: // Ncc5Nii2Njj5
        case 416: // Ncc5Nii2Njj7
        case 417: // Ncc5Nii2Njj8
        case 425: // Ncc5Nii3Njj6
        case 426: // Ncc5Nii3Njj7
        case 437: // Ncc5Nii4Njj8
        case 438: // Ncc5Nii4Njj9
        case 440: // Ncc5Nii5Njj1
        case 441: // Ncc5Nii5Njj2
        case 450: // Ncc5Nii6Njj1
        case 452: // Ncc5Nii6Njj3
        case 461: // Ncc5Nii7Njj2
        case 462: // Ncc5Nii7Njj3
        case 471: // Ncc5Nii8Njj2
        case 473: // Ncc5Nii8Njj4
        case 480: // Ncc5Nii9Njj1
        case 483: // Ncc5Nii9Njj4
        case 504: // Ncc6Nii1Njj5
        case 506: // Ncc6Nii1Njj7
        case 515: // Ncc6Nii2Njj6
        case 516: // Ncc6Nii2Njj7
        case 519: // Ncc6Nii2Njj10
        case 524: // Ncc6Nii3Njj5
        case 525: // Ncc6Nii3Njj6
        case 528: // Ncc6Nii3Njj9
        case 538: // Ncc6Nii4Njj9
        case 539: // Ncc6Nii4Njj10
        case 540: // Ncc6Nii5Njj1
        case 542: // Ncc6Nii5Njj3
        case 551: // Ncc6Nii6Njj2
        case 552: // Ncc6Nii6Njj3
        case 560: // Ncc6Nii7Njj1
        case 561: // Ncc6Nii7Njj2
        case 582: // Ncc6Nii9Njj3
        case 583: // Ncc6Nii9Njj4
        case 591: // Ncc6Nii10Njj2
        case 593: // Ncc6Nii10Njj4
        case 605: // Ncc7Nii1Njj6
        case 606: // Ncc7Nii1Njj7
        case 609: // Ncc7Nii1Njj10
        case 614: // Ncc7Nii2Njj5
        case 615: // Ncc7Nii2Njj6
        case 624: // Ncc7Nii3Njj5
        case 626: // Ncc7Nii3Njj7
        case 627: // Ncc7Nii3Njj8
        case 637: // Ncc7Nii4Njj8
        case 639: // Ncc7Nii4Njj10
        case 641: // Ncc7Nii5Njj2
        case 642: // Ncc7Nii5Njj3
        case 650: // Ncc7Nii6Njj1
        case 651: // Ncc7Nii6Njj2
        case 660: // Ncc7Nii7Njj1
        case 662: // Ncc7Nii7Njj3
        case 672: // Ncc7Nii8Njj3
        case 673: // Ncc7Nii8Njj4
        case 690: // Ncc7Nii10Njj1
        case 693: // Ncc7Nii10Njj4
        case 707: // Ncc8Nii1Njj8
        case 708: // Ncc8Nii1Njj9
        case 709: // Ncc8Nii1Njj10
        case 714: // Ncc8Nii2Njj5
        case 718: // Ncc8Nii2Njj9
        case 726: // Ncc8Nii3Njj7
        case 729: // Ncc8Nii3Njj10
        case 734: // Ncc8Nii4Njj5
        case 736: // Ncc8Nii4Njj7
        case 737: // Ncc8Nii4Njj8
        case 741: // Ncc8Nii5Njj2
        case 743: // Ncc8Nii5Njj4
        case 762: // Ncc8Nii7Njj3
        case 763: // Ncc8Nii7Njj4
        case 770: // Ncc8Nii8Njj1
        case 773: // Ncc8Nii8Njj4
        case 780: // Ncc8Nii9Njj1
        case 781: // Ncc8Nii9Njj2
        case 790: // Ncc8Nii10Njj1
        case 792: // Ncc8Nii10Njj3
        case 804: // Ncc9Nii1Njj5
        case 807: // Ncc9Nii1Njj8
        case 817: // Ncc9Nii2Njj8
        case 818: // Ncc9Nii2Njj9
        case 819: // Ncc9Nii2Njj10
        case 825: // Ncc9Nii3Njj6
        case 829: // Ncc9Nii3Njj10
        case 834: // Ncc9Nii4Njj5
        case 835: // Ncc9Nii4Njj6
        case 838: // Ncc9Nii4Njj9
        case 840: // Ncc9Nii5Njj1
        case 843: // Ncc9Nii5Njj4
        case 852: // Ncc9Nii6Njj3
        case 853: // Ncc9Nii6Njj4
        case 870: // Ncc9Nii8Njj1
        case 871: // Ncc9Nii8Njj2
        case 881: // Ncc9Nii9Njj2
        case 883: // Ncc9Nii9Njj4
        case 891: // Ncc9Nii10Njj2
        case 892: // Ncc9Nii10Njj3
        case 906: // Ncc10Nii1Njj7
        case 907: // Ncc10Nii1Njj8
        case 915: // Ncc10Nii2Njj6
        case 918: // Ncc10Nii2Njj9
        case 927: // Ncc10Nii3Njj8
        case 928: // Ncc10Nii3Njj9
        case 929: // Ncc10Nii3Njj10
        case 935: // Ncc10Nii4Njj6
        case 936: // Ncc10Nii4Njj7
        case 939: // Ncc10Nii4Njj10
        case 951: // Ncc10Nii6Njj2
        case 953: // Ncc10Nii6Njj4
        case 960: // Ncc10Nii7Njj1
        case 963: // Ncc10Nii7Njj4
        case 970: // Ncc10Nii8Njj1
        case 972: // Ncc10Nii8Njj3
        case 981: // Ncc10Nii9Njj2
        case 982: // Ncc10Nii9Njj3
        case 992: // Ncc10Nii10Njj3
        case 993: // Ncc10Nii10Njj4
        case 44: // Ncc1Nii5Njj5
            return ((-1.000000e+00) / (1.890000e+03)) * (ADB);
        case 47: // Ncc1Nii5Njj8
        case 49: // Ncc1Nii5Njj10
        case 57: // Ncc1Nii6Njj8
        case 64: // Ncc1Nii7Njj5
        case 67: // Ncc1Nii7Njj8
        case 68: // Ncc1Nii7Njj9
        case 74: // Ncc1Nii8Njj5
        case 75: // Ncc1Nii8Njj6
        case 76: // Ncc1Nii8Njj7
        case 86: // Ncc1Nii9Njj7
        case 94: // Ncc1Nii10Njj5
        case 145: // Ncc2Nii5Njj6
        case 148: // Ncc2Nii5Njj9
        case 149: // Ncc2Nii5Njj10
        case 154: // Ncc2Nii6Njj5
        case 157: // Ncc2Nii6Njj8
        case 158: // Ncc2Nii6Njj9
        case 168: // Ncc2Nii7Njj9
        case 175: // Ncc2Nii8Njj6
        case 184: // Ncc2Nii9Njj5
        case 185: // Ncc2Nii9Njj6
        case 186: // Ncc2Nii9Njj7
        case 194: // Ncc2Nii10Njj5
        case 249: // Ncc3Nii5Njj10
        case 256: // Ncc3Nii6Njj7
        case 257: // Ncc3Nii6Njj8
        case 259: // Ncc3Nii6Njj10
        case 265: // Ncc3Nii7Njj6
        case 268: // Ncc3Nii7Njj9
        case 269: // Ncc3Nii7Njj10
        case 275: // Ncc3Nii8Njj6
        case 286: // Ncc3Nii9Njj7
        case 294: // Ncc3Nii10Njj5
        case 295: // Ncc3Nii10Njj6
        case 296: // Ncc3Nii10Njj7
        case 349: // Ncc4Nii5Njj10
        case 357: // Ncc4Nii6Njj8
        case 368: // Ncc4Nii7Njj9
        case 375: // Ncc4Nii8Njj6
        case 378: // Ncc4Nii8Njj9
        case 379: // Ncc4Nii8Njj10
        case 386: // Ncc4Nii9Njj7
        case 387: // Ncc4Nii9Njj8
        case 389: // Ncc4Nii9Njj10
        case 394: // Ncc4Nii10Njj5
        case 397: // Ncc4Nii10Njj8
        case 398: // Ncc4Nii10Njj9
        case 406: // Ncc5Nii1Njj7
        case 407: // Ncc5Nii1Njj8
        case 409: // Ncc5Nii1Njj10
        case 415: // Ncc5Nii2Njj6
        case 418: // Ncc5Nii2Njj9
        case 419: // Ncc5Nii2Njj10
        case 429: // Ncc5Nii3Njj10
        case 439: // Ncc5Nii4Njj10
        case 451: // Ncc5Nii6Njj2
        case 460: // Ncc5Nii7Njj1
        case 470: // Ncc5Nii8Njj1
        case 481: // Ncc5Nii9Njj2
        case 490: // Ncc5Nii10Njj1
        case 491: // Ncc5Nii10Njj2
        case 492: // Ncc5Nii10Njj3
        case 493: // Ncc5Nii10Njj4
        case 507: // Ncc6Nii1Njj8
        case 514: // Ncc6Nii2Njj5
        case 517: // Ncc6Nii2Njj8
        case 518: // Ncc6Nii2Njj9
        case 526: // Ncc6Nii3Njj7
        case 527: // Ncc6Nii3Njj8
        case 529: // Ncc6Nii3Njj10
        case 537: // Ncc6Nii4Njj8
        case 541: // Ncc6Nii5Njj2
        case 562: // Ncc6Nii7Njj3
        case 570: // Ncc6Nii8Njj1
        case 571: // Ncc6Nii8Njj2
        case 572: // Ncc6Nii8Njj3
        case 573: // Ncc6Nii8Njj4
        case 581: // Ncc6Nii9Njj2
        case 592: // Ncc6Nii10Njj3
        case 604: // Ncc7Nii1Njj5
        case 607: // Ncc7Nii1Njj8
        case 608: // Ncc7Nii1Njj9
        case 618: // Ncc7Nii2Njj9
        case 625: // Ncc7Nii3Njj6
        case 628: // Ncc7Nii3Njj9
        case 629: // Ncc7Nii3Njj10
        case 638: // Ncc7Nii4Njj9
        case 640: // Ncc7Nii5Njj1
        case 652: // Ncc7Nii6Njj3
        case 670: // Ncc7Nii8Njj1
        case 680: // Ncc7Nii9Njj1
        case 681: // Ncc7Nii9Njj2
        case 682: // Ncc7Nii9Njj3
        case 683: // Ncc7Nii9Njj4
        case 692: // Ncc7Nii10Njj3
        case 704: // Ncc8Nii1Njj5
        case 705: // Ncc8Nii1Njj6
        case 706: // Ncc8Nii1Njj7
        case 715: // Ncc8Nii2Njj6
        case 725: // Ncc8Nii3Njj6
        case 735: // Ncc8Nii4Njj6
        case 738: // Ncc8Nii4Njj9
        case 739: // Ncc8Nii4Njj10
        case 740: // Ncc8Nii5Njj1
        case 750: // Ncc8Nii6Njj1
        case 751: // Ncc8Nii6Njj2
        case 752: // Ncc8Nii6Njj3
        case 753: // Ncc8Nii6Njj4
        case 760: // Ncc8Nii7Njj1
        case 783: // Ncc8Nii9Njj4
        case 793: // Ncc8Nii10Njj4
        case 806: // Ncc9Nii1Njj7
        case 814: // Ncc9Nii2Njj5
        case 815: // Ncc9Nii2Njj6
        case 816: // Ncc9Nii2Njj7
        case 826: // Ncc9Nii3Njj7
        case 836: // Ncc9Nii4Njj7
        case 837: // Ncc9Nii4Njj8
        case 839: // Ncc9Nii4Njj10
        case 841: // Ncc9Nii5Njj2
        case 851: // Ncc9Nii6Njj2
        case 860: // Ncc9Nii7Njj1
        case 861: // Ncc9Nii7Njj2
        case 862: // Ncc9Nii7Njj3
        case 863: // Ncc9Nii7Njj4
        case 873: // Ncc9Nii8Njj4
        case 893: // Ncc9Nii10Njj4
        case 904: // Ncc10Nii1Njj5
        case 914: // Ncc10Nii2Njj5
        case 924: // Ncc10Nii3Njj5
        case 925: // Ncc10Nii3Njj6
        case 926: // Ncc10Nii3Njj7
        case 934: // Ncc10Nii4Njj5
        case 937: // Ncc10Nii4Njj8
        case 938: // Ncc10Nii4Njj9
        case 940: // Ncc10Nii5Njj1
        case 941: // Ncc10Nii5Njj2
        case 942: // Ncc10Nii5Njj3
        case 943: // Ncc10Nii5Njj4
        case 952: // Ncc10Nii6Njj3
        case 962: // Ncc10Nii7Njj3
        case 973: // Ncc10Nii8Njj4
        case 983: // Ncc10Nii9Njj4
        case 46: // Ncc1Nii5Njj7
            return ((-1.000000e+00) / (3.780000e+03)) * (ADB);
        case 88: // Ncc1Nii9Njj9
        case 99: // Ncc1Nii10Njj10
        case 166: // Ncc2Nii7Njj7
        case 177: // Ncc2Nii8Njj8
        case 199: // Ncc2Nii10Njj10
        case 244: // Ncc3Nii5Njj5
        case 277: // Ncc3Nii8Njj8
        case 288: // Ncc3Nii9Njj9
        case 344: // Ncc4Nii5Njj5
        case 355: // Ncc4Nii6Njj6
        case 366: // Ncc4Nii7Njj7
        case 424: // Ncc5Nii3Njj5
        case 434: // Ncc5Nii4Njj5
        case 442: // Ncc5Nii5Njj3
        case 443: // Ncc5Nii5Njj4
        case 505: // Ncc6Nii1Njj6
        case 535: // Ncc6Nii4Njj6
        case 550: // Ncc6Nii6Njj1
        case 553: // Ncc6Nii6Njj4
        case 616: // Ncc7Nii2Njj7
        case 636: // Ncc7Nii4Njj7
        case 661: // Ncc7Nii7Njj2
        case 663: // Ncc7Nii7Njj4
        case 717: // Ncc8Nii2Njj8
        case 727: // Ncc8Nii3Njj8
        case 771: // Ncc8Nii8Njj2
        case 772: // Ncc8Nii8Njj3
        case 808: // Ncc9Nii1Njj9
        case 828: // Ncc9Nii3Njj9
        case 880: // Ncc9Nii9Njj1
        case 882: // Ncc9Nii9Njj3
        case 909: // Ncc10Nii1Njj10
        case 919: // Ncc10Nii2Njj10
        case 990: // Ncc10Nii10Njj1
        case 991: // Ncc10Nii10Njj2
        case 55: // Ncc1Nii6Njj6
            return ((-1.000000e+00) / (1.134000e+03)) * (ADB);
        case 59: // Ncc1Nii6Njj10
        case 85: // Ncc1Nii9Njj6
        case 89: // Ncc1Nii9Njj10
        case 95: // Ncc1Nii10Njj6
        case 98: // Ncc1Nii10Njj9
        case 167: // Ncc2Nii7Njj8
        case 169: // Ncc2Nii7Njj10
        case 176: // Ncc2Nii8Njj7
        case 179: // Ncc2Nii8Njj10
        case 196: // Ncc2Nii10Njj7
        case 197: // Ncc2Nii10Njj8
        case 247: // Ncc3Nii5Njj8
        case 248: // Ncc3Nii5Njj9
        case 274: // Ncc3Nii8Njj5
        case 278: // Ncc3Nii8Njj9
        case 284: // Ncc3Nii9Njj5
        case 287: // Ncc3Nii9Njj8
        case 345: // Ncc4Nii5Njj6
        case 346: // Ncc4Nii5Njj7
        case 354: // Ncc4Nii6Njj5
        case 356: // Ncc4Nii6Njj7
        case 364: // Ncc4Nii7Njj5
        case 365: // Ncc4Nii7Njj6
        case 427: // Ncc5Nii3Njj8
        case 428: // Ncc5Nii3Njj9
        case 435: // Ncc5Nii4Njj6
        case 436: // Ncc5Nii4Njj7
        case 453: // Ncc5Nii6Njj4
        case 463: // Ncc5Nii7Njj4
        case 472: // Ncc5Nii8Njj3
        case 482: // Ncc5Nii9Njj3
        case 508: // Ncc6Nii1Njj9
        case 509: // Ncc6Nii1Njj10
        case 534: // Ncc6Nii4Njj5
        case 536: // Ncc6Nii4Njj7
        case 543: // Ncc6Nii5Njj4
        case 563: // Ncc6Nii7Njj4
        case 580: // Ncc6Nii9Njj1
        case 590: // Ncc6Nii10Njj1
        case 617: // Ncc7Nii2Njj8
        case 619: // Ncc7Nii2Njj10
        case 634: // Ncc7Nii4Njj5
        case 635: // Ncc7Nii4Njj6
        case 643: // Ncc7Nii5Njj4
        case 653: // Ncc7Nii6Njj4
        case 671: // Ncc7Nii8Njj2
        case 691: // Ncc7Nii10Njj2
        case 716: // Ncc8Nii2Njj7
        case 719: // Ncc8Nii2Njj10
        case 724: // Ncc8Nii3Njj5
        case 728: // Ncc8Nii3Njj9
        case 742: // Ncc8Nii5Njj3
        case 761: // Ncc8Nii7Njj2
        case 782: // Ncc8Nii9Njj3
        case 791: // Ncc8Nii10Njj2
        case 805: // Ncc9Nii1Njj6
        case 809: // Ncc9Nii1Njj10
        case 824: // Ncc9Nii3Njj5
        case 827: // Ncc9Nii3Njj8
        case 842: // Ncc9Nii5Njj3
        case 850: // Ncc9Nii6Njj1
        case 872: // Ncc9Nii8Njj3
        case 890: // Ncc9Nii10Njj1
        case 905: // Ncc10Nii1Njj6
        case 908: // Ncc10Nii1Njj9
        case 916: // Ncc10Nii2Njj7
        case 917: // Ncc10Nii2Njj8
        case 950: // Ncc10Nii6Njj1
        case 961: // Ncc10Nii7Njj2
        case 971: // Ncc10Nii8Njj2
        case 980: // Ncc10Nii9Njj1
        case 58: // Ncc1Nii6Njj9
            return ((-1.000000e+00) / (2.268000e+03)) * (ADB);
        case 555: // Ncc6Nii6Njj6
        case 666: // Ncc7Nii7Njj7
        case 777: // Ncc8Nii8Njj8
        case 888: // Ncc9Nii9Njj9
        case 999: // Ncc10Nii10Njj10
        case 444: // Ncc5Nii5Njj5
            return ((2.000000e+00) / (3.150000e+02)) * (ADB);
        case 446: // Ncc5Nii5Njj7
        case 447: // Ncc5Nii5Njj8
        case 448: // Ncc5Nii5Njj9
        case 454: // Ncc5Nii6Njj5
        case 455: // Ncc5Nii6Njj6
        case 464: // Ncc5Nii7Njj5
        case 466: // Ncc5Nii7Njj7
        case 474: // Ncc5Nii8Njj5
        case 477: // Ncc5Nii8Njj8
        case 484: // Ncc5Nii9Njj5
        case 488: // Ncc5Nii9Njj9
        case 544: // Ncc6Nii5Njj5
        case 545: // Ncc6Nii5Njj6
        case 554: // Ncc6Nii6Njj5
        case 556: // Ncc6Nii6Njj7
        case 558: // Ncc6Nii6Njj9
        case 559: // Ncc6Nii6Njj10
        case 565: // Ncc6Nii7Njj6
        case 566: // Ncc6Nii7Njj7
        case 585: // Ncc6Nii9Njj6
        case 588: // Ncc6Nii9Njj9
        case 595: // Ncc6Nii10Njj6
        case 599: // Ncc6Nii10Njj10
        case 644: // Ncc7Nii5Njj5
        case 646: // Ncc7Nii5Njj7
        case 655: // Ncc7Nii6Njj6
        case 656: // Ncc7Nii6Njj7
        case 664: // Ncc7Nii7Njj5
        case 665: // Ncc7Nii7Njj6
        case 667: // Ncc7Nii7Njj8
        case 669: // Ncc7Nii7Njj10
        case 676: // Ncc7Nii8Njj7
        case 677: // Ncc7Nii8Njj8
        case 696: // Ncc7Nii10Njj7
        case 699: // Ncc7Nii10Njj10
        case 744: // Ncc8Nii5Njj5
        case 747: // Ncc8Nii5Njj8
        case 766: // Ncc8Nii7Njj7
        case 767: // Ncc8Nii7Njj8
        case 774: // Ncc8Nii8Njj5
        case 776: // Ncc8Nii8Njj7
        case 778: // Ncc8Nii8Njj9
        case 779: // Ncc8Nii8Njj10
        case 787: // Ncc8Nii9Njj8
        case 788: // Ncc8Nii9Njj9
        case 797: // Ncc8Nii10Njj8
        case 799: // Ncc8Nii10Njj10
        case 844: // Ncc9Nii5Njj5
        case 848: // Ncc9Nii5Njj9
        case 855: // Ncc9Nii6Njj6
        case 858: // Ncc9Nii6Njj9
        case 877: // Ncc9Nii8Njj8
        case 878: // Ncc9Nii8Njj9
        case 884: // Ncc9Nii9Njj5
        case 885: // Ncc9Nii9Njj6
        case 887: // Ncc9Nii9Njj8
        case 889: // Ncc9Nii9Njj10
        case 898: // Ncc9Nii10Njj9
        case 899: // Ncc9Nii10Njj10
        case 955: // Ncc10Nii6Njj6
        case 959: // Ncc10Nii6Njj10
        case 966: // Ncc10Nii7Njj7
        case 969: // Ncc10Nii7Njj10
        case 977: // Ncc10Nii8Njj8
        case 979: // Ncc10Nii8Njj10
        case 988: // Ncc10Nii9Njj9
        case 989: // Ncc10Nii9Njj10
        case 995: // Ncc10Nii10Njj6
        case 996: // Ncc10Nii10Njj7
        case 997: // Ncc10Nii10Njj8
        case 998: // Ncc10Nii10Njj9
        case 445: // Ncc5Nii5Njj6
            return ((2.000000e+00) / (9.450000e+02)) * (ADB);
        case 457: // Ncc5Nii6Njj8
        case 459: // Ncc5Nii6Njj10
        case 468: // Ncc5Nii7Njj9
        case 469: // Ncc5Nii7Njj10
        case 475: // Ncc5Nii8Njj6
        case 479: // Ncc5Nii8Njj10
        case 486: // Ncc5Nii9Njj7
        case 489: // Ncc5Nii9Njj10
        case 494: // Ncc5Nii10Njj5
        case 495: // Ncc5Nii10Njj6
        case 496: // Ncc5Nii10Njj7
        case 497: // Ncc5Nii10Njj8
        case 498: // Ncc5Nii10Njj9
        case 499: // Ncc5Nii10Njj10
        case 547: // Ncc6Nii5Njj8
        case 549: // Ncc6Nii5Njj10
        case 557: // Ncc6Nii6Njj8
        case 567: // Ncc6Nii7Njj8
        case 568: // Ncc6Nii7Njj9
        case 574: // Ncc6Nii8Njj5
        case 575: // Ncc6Nii8Njj6
        case 576: // Ncc6Nii8Njj7
        case 577: // Ncc6Nii8Njj8
        case 578: // Ncc6Nii8Njj9
        case 579: // Ncc6Nii8Njj10
        case 586: // Ncc6Nii9Njj7
        case 587: // Ncc6Nii9Njj8
        case 594: // Ncc6Nii10Njj5
        case 597: // Ncc6Nii10Njj8
        case 648: // Ncc7Nii5Njj9
        case 649: // Ncc7Nii5Njj10
        case 657: // Ncc7Nii6Njj8
        case 658: // Ncc7Nii6Njj9
        case 668: // Ncc7Nii7Njj9
        case 675: // Ncc7Nii8Njj6
        case 678: // Ncc7Nii8Njj9
        case 684: // Ncc7Nii9Njj5
        case 685: // Ncc7Nii9Njj6
        case 686: // Ncc7Nii9Njj7
        case 687: // Ncc7Nii9Njj8
        case 688: // Ncc7Nii9Njj9
        case 689: // Ncc7Nii9Njj10
        case 694: // Ncc7Nii10Njj5
        case 698: // Ncc7Nii10Njj9
        case 745: // Ncc8Nii5Njj6
        case 749: // Ncc8Nii5Njj10
        case 754: // Ncc8Nii6Njj5
        case 755: // Ncc8Nii6Njj6
        case 756: // Ncc8Nii6Njj7
        case 757: // Ncc8Nii6Njj8
        case 758: // Ncc8Nii6Njj9
        case 759: // Ncc8Nii6Njj10
        case 765: // Ncc8Nii7Njj6
        case 768: // Ncc8Nii7Njj9
        case 775: // Ncc8Nii8Njj6
        case 785: // Ncc8Nii9Njj6
        case 786: // Ncc8Nii9Njj7
        case 794: // Ncc8Nii10Njj5
        case 795: // Ncc8Nii10Njj6
        case 846: // Ncc9Nii5Njj7
        case 849: // Ncc9Nii5Njj10
        case 856: // Ncc9Nii6Njj7
        case 857: // Ncc9Nii6Njj8
        case 864: // Ncc9Nii7Njj5
        case 865: // Ncc9Nii7Njj6
        case 866: // Ncc9Nii7Njj7
        case 867: // Ncc9Nii7Njj8
        case 868: // Ncc9Nii7Njj9
        case 869: // Ncc9Nii7Njj10
        case 875: // Ncc9Nii8Njj6
        case 876: // Ncc9Nii8Njj7
        case 886: // Ncc9Nii9Njj7
        case 894: // Ncc9Nii10Njj5
        case 896: // Ncc9Nii10Njj7
        case 944: // Ncc10Nii5Njj5
        case 945: // Ncc10Nii5Njj6
        case 946: // Ncc10Nii5Njj7
        case 947: // Ncc10Nii5Njj8
        case 948: // Ncc10Nii5Njj9
        case 949: // Ncc10Nii5Njj10
        case 954: // Ncc10Nii6Njj5
        case 957: // Ncc10Nii6Njj8
        case 964: // Ncc10Nii7Njj5
        case 968: // Ncc10Nii7Njj9
        case 974: // Ncc10Nii8Njj5
        case 975: // Ncc10Nii8Njj6
        case 984: // Ncc10Nii9Njj5
        case 986: // Ncc10Nii9Njj7
        case 994: // Ncc10Nii10Njj5
        case 449: // Ncc5Nii5Njj10
            return ((2.000000e+00) / (2.835000e+03)) * (ADB);
        case 465: // Ncc5Nii7Njj6
        case 478: // Ncc5Nii8Njj9
        case 487: // Ncc5Nii9Njj8
        case 546: // Ncc6Nii5Njj7
        case 564: // Ncc6Nii7Njj5
        case 589: // Ncc6Nii9Njj10
        case 598: // Ncc6Nii10Njj9
        case 645: // Ncc7Nii5Njj6
        case 654: // Ncc7Nii6Njj5
        case 679: // Ncc7Nii8Njj10
        case 697: // Ncc7Nii10Njj8
        case 748: // Ncc8Nii5Njj9
        case 769: // Ncc8Nii7Njj10
        case 784: // Ncc8Nii9Njj5
        case 796: // Ncc8Nii10Njj7
        case 847: // Ncc9Nii5Njj8
        case 859: // Ncc9Nii6Njj10
        case 874: // Ncc9Nii8Njj5
        case 895: // Ncc9Nii10Njj6
        case 958: // Ncc10Nii6Njj9
        case 967: // Ncc10Nii7Njj8
        case 976: // Ncc10Nii8Njj7
        case 985: // Ncc10Nii9Njj6
        case 456: // Ncc5Nii6Njj7
            return ((4.000000e+00) / (2.835000e+03)) * (ADB);
        case 467: // Ncc5Nii7Njj8
        case 476: // Ncc5Nii8Njj7
        case 485: // Ncc5Nii9Njj6
        case 548: // Ncc6Nii5Njj9
        case 569: // Ncc6Nii7Njj10
        case 584: // Ncc6Nii9Njj5
        case 596: // Ncc6Nii10Njj7
        case 647: // Ncc7Nii5Njj8
        case 659: // Ncc7Nii6Njj10
        case 674: // Ncc7Nii8Njj5
        case 695: // Ncc7Nii10Njj6
        case 746: // Ncc8Nii5Njj7
        case 764: // Ncc8Nii7Njj5
        case 789: // Ncc8Nii9Njj10
        case 798: // Ncc8Nii10Njj9
        case 845: // Ncc9Nii5Njj6
        case 854: // Ncc9Nii6Njj5
        case 879: // Ncc9Nii8Njj10
        case 897: // Ncc9Nii10Njj8
        case 956: // Ncc10Nii6Njj7
        case 965: // Ncc10Nii7Njj6
        case 978: // Ncc10Nii8Njj9
        case 987: // Ncc10Nii9Njj8
        case 458: // Ncc5Nii6Njj9
            return ((1.000000e+00) / (9.450000e+02)) * (ADB);
		default:
			return 0.0;		
	}	                       
}

double  Element3DTetrahedron2nd::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {

	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_1
	                       + 10 * diffop;	
	switch(obj_index) {
        case 4: // da1Nii5
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 5: // da1Nii6
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 6: // da1Nii7
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[3] + ((-1.000000e+00) * jm[5])) * jm[7]) + (jm[4] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 7: // da1Nii8
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[3]) + jm[4]) * jm[8]));
        case 8: // da1Nii9
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (((((-1.000000e+00) * jm[3]) + jm[5]) * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 9: // da1Nii10
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[5] * jm[6]) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 14: // da2Nii5
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[6]) + (jm[0] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 15: // da2Nii6
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 16: // da2Nii7
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[7]) + (jm[1] * (jm[6] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 17: // da2Nii8
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 18: // da2Nii9
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[7]) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 19: // da2Nii10
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 24: // da3Nii5
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 25: // da3Nii6
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 26: // da3Nii7
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (((jm[0] + ((-1.000000e+00) * jm[2])) * jm[4]) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 27: // da3Nii8
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + ((((-1.000000e+00) * jm[0]) + jm[1]) * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 28: // da3Nii9
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * (((((-1.000000e+00) * jm[0]) + jm[2]) * jm[4]) + (jm[1] * (jm[3] + ((-1.000000e+00) * jm[5])))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 29: // da3Nii10
            return ((1.000000e+00) / (6.000000e+00)) * (ADB) * ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[0] * (((-1.000000e+00) * jm[4]) + jm[5]))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
		
		default:            		
			return 0.0;
	}
	
}
double Element3DTetrahedron2nd::get_single_integral_0th_order(short elem_shape_func_1) const {

	switch(elem_shape_func_1) {
        case 1: // Nii2
        case 2: // Nii3
        case 3: // Nii4
        case 0: // Nii1
            return ((-1.000000e+00) / (1.200000e+02)) * (ADB);
        case 5: // Nii6
        case 6: // Nii7
        case 7: // Nii8
        case 8: // Nii9
        case 9: // Nii10
        case 4: // Nii5
            return ((1.000000e+00) / (3.000000e+01)) * (ADB);		
		default:
			TDKP_GENERAL_EXCEPTION("invalid index");
	}
		
}
	
					
double Element3DTetrahedron2nd::evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const {

	double global[] = { global_x, global_y, global_z } ;
	double local[3];
	global2local(global, local);
	return evaluate_form_function(elem_shape_func, local);
}

	
double Element3DTetrahedron2nd::evaluate_form_function(
	short elem_shape_func, 
	const double* local
) const {
	
	const double& x = local[0];
	const double& y = local[1];
	const double& z = local[2];
	
	switch(elem_shape_func) {
		case 0:
			return 2.0*x*x + 2.0*y*y + 2.0*z*z + 4.0*x*y + 4.0*x*z + 4.0*y*z - 3.0*x - 3.0*y - 3.0*z + 1.0;
		case 1:
			return 2.0*x*x - 1.0*x;
		case 2:
			return 2.0*y*y - 1.0*y;
		case 3:
			return 2.0*z*z - 1.0*z;
		case 4:
			return -4.0*x*x - 4.0*x*y - 4.0*x*z + 4.0*x;
		case 5:
			return 4.0*x*y;
		case 6:
			return -4.0*y*y -4.0*x*y -4.0*y*z + 4.0*y;
		case 7:
			return -4.0*z*z -4.0*x*z -4.0*y*z + 4.0*z;
		case 8:
			return 4.0*x*z;
		case 9:
			return 4.0*y*z;
		default:
			TDKP_GENERAL_EXCEPTION("invalid shape function");						
	}	
	
}

double Element3DTetrahedron2nd::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local
) const {
	return evaluate_form_function_derivative(diffop, elem_shape_func, local[0], local[1], local[2]);  	
}

double Element3DTetrahedron2nd::evaluate_form_function_derivative(short diffop, short elem_shape_func, const double& x, const double& y, const double& z) const {
	
	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func
	                       + 10 * diffop;	
	switch(obj_index) {
        case 0: // da1Nii1
            return ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y)) + ((4.000000e+00) * (z))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 1: // da1Nii2
            return ((-1.000000e+00) + ((4.000000e+00) * (x))) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 2: // da1Nii3
            return ((-1.000000e+00) + ((4.000000e+00) * (y))) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 3: // da1Nii4
            return ((-1.000000e+00) + ((4.000000e+00) * (z))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 4: // da1Nii5
            return (4.000000e+00) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * ((-1.000000e+00) + (y) + (z)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + ((x) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))));
        case 5: // da1Nii6
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((x) * jm[5] * jm[6]) + ((-1.000000e+00) * (y) * jm[5] * jm[7]) + ((-1.000000e+00) * (x) * jm[3] * jm[8]) + ((y) * jm[4] * jm[8]));
        case 6: // da1Nii7
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * ((-1.000000e+00) + (x) + (z)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]))) + ((y) * (((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[3] * jm[8]))));
        case 7: // da1Nii8
            return (4.000000e+00) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (z) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * ((((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[7]) + ((-1.000000e+00) * (z) * jm[8]))) + (jm[4] * (((-1.000000e+00) * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[6]) + ((z) * jm[8]))));
        case 8: // da1Nii9
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * (x) * jm[4] * jm[6]) + ((x) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * (z) * jm[7]) + (jm[4] * (z) * jm[8]));
        case 9: // da1Nii10
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * (y) * jm[4] * jm[6]) + (jm[5] * (z) * jm[6]) + ((y) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * (z) * jm[8]));
        case 10: // da2Nii1
            return ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y)) + ((4.000000e+00) * (z))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 11: // da2Nii2
            return ((-1.000000e+00) + ((4.000000e+00) * (x))) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 12: // da2Nii3
            return ((-1.000000e+00) + ((4.000000e+00) * (y))) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 13: // da2Nii4
            return ((-1.000000e+00) + ((4.000000e+00) * (z))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 14: // da2Nii5
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * ((-1.000000e+00) + (y) + (z)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]))) + ((x) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * jm[7]) + ((-2.000000e+00) * jm[2] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8]) + ((2.000000e+00) * jm[1] * jm[8]))));
        case 15: // da2Nii6
            return (4.000000e+00) * (((-1.000000e+00) * (x) * jm[2] * jm[6]) + (jm[2] * (y) * jm[7]) + ((x) * jm[0] * jm[8]) + ((-1.000000e+00) * jm[1] * (y) * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 16: // da2Nii7
            return (4.000000e+00) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (((-1.000000e+00) * ((-1.000000e+00) + (x) + ((2.000000e+00) * (y)) + (z)) * jm[6]) + ((y) * jm[7]))) + (jm[1] * (y) * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * (y) * (jm[7] + ((-2.000000e+00) * jm[8]))) + (((-1.000000e+00) + (x) + (z)) * jm[8]))));
        case 17: // da2Nii8
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[2] * (z) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * ((((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[7]) + ((-1.000000e+00) * (z) * jm[8]))) + (jm[1] * (((-1.000000e+00) * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[6]) + ((z) * jm[8]))));
        case 18: // da2Nii9
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((x) * jm[1] * jm[6]) + ((-1.000000e+00) * (x) * jm[0] * jm[7]) + (jm[2] * (z) * jm[7]) + ((-1.000000e+00) * jm[1] * (z) * jm[8]));
        case 19: // da2Nii10
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[1] * (y) * jm[6]) + ((-1.000000e+00) * jm[2] * (z) * jm[6]) + ((-1.000000e+00) * jm[0] * (y) * jm[7]) + (jm[0] * (z) * jm[8]));
        case 20: // da3Nii1
            return ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y)) + ((4.000000e+00) * (z))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 21: // da3Nii2
            return ((-1.000000e+00) + ((4.000000e+00) * (x))) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 22: // da3Nii3
            return ((-1.000000e+00) + ((4.000000e+00) * (y))) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 23: // da3Nii4
            return ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((-1.000000e+00) + ((4.000000e+00) * (z))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 24: // da3Nii5
            return (4.000000e+00) * (((x) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + (jm[0] * jm[4]) + ((-2.000000e+00) * jm[2] * jm[4]) + ((-1.000000e+00) * jm[0] * jm[5]) + ((2.000000e+00) * jm[1] * jm[5]))) + ((-1.000000e+00) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((-1.000000e+00) + (y) + (z)))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 25: // da3Nii6
            return (4.000000e+00) * (((x) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[2] * (y) * jm[4]) + ((-1.000000e+00) * (x) * jm[0] * jm[5]) + (jm[1] * (y) * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 26: // da3Nii7
            return (4.000000e+00) * ((jm[1] * (y) * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * (y) * (jm[4] + ((-2.000000e+00) * jm[5]))) + (jm[5] * ((-1.000000e+00) + (x) + (z))))) + (jm[2] * (((y) * jm[4]) + ((-1.000000e+00) * jm[3] * ((-1.000000e+00) + (x) + ((2.000000e+00) * (y)) + (z)))))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 27: // da3Nii8
            return (-4.000000e+00) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (z)) + (jm[1] * ((jm[5] * (z)) + ((-1.000000e+00) * jm[3] * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z)))))) + (jm[0] * (((-1.000000e+00) * jm[5] * (z)) + (jm[4] * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))))))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 28: // da3Nii9
            return (4.000000e+00) * (((-1.000000e+00) * (x) * jm[1] * jm[3]) + ((x) * jm[0] * jm[4]) + ((-1.000000e+00) * jm[2] * jm[4] * (z)) + (jm[1] * jm[5] * (z))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 29: // da3Nii10
            return (4.000000e+00) * (((-1.000000e+00) * jm[1] * (y) * jm[3]) + (jm[0] * (y) * jm[4]) + (jm[2] * jm[3] * (z)) + ((-1.000000e+00) * jm[0] * jm[5] * (z))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
		
		/* OLD 
        case 0: // da1Nii1
            return ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y)) + ((4.000000e+00) * (z))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (((-1.000000e+00) * jm[6]) + jm[7])) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[3] * (((-1.000000e+00) * jm[7]) + jm[8])));
        case 1: // da1Nii2
            return ((-1.000000e+00) + ((4.000000e+00) * (x))) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 2: // da1Nii3
            return ((-1.000000e+00) + ((4.000000e+00) * (y))) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 3: // da1Nii4
            return ((-1.000000e+00) + ((4.000000e+00) * (z))) * ((jm[4] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 4: // da1Nii5
            return (4.000000e+00) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * ((-1.000000e+00) + (y) + (z)) * ((jm[5] * jm[7]) + ((-1.000000e+00) * jm[4] * jm[8]))) + ((x) * (((-1.000000e+00) * jm[4] * jm[6]) + (jm[5] * jm[6]) + (jm[3] * jm[7]) + ((-2.000000e+00) * jm[5] * jm[7]) + ((-1.000000e+00) * jm[3] * jm[8]) + ((2.000000e+00) * jm[4] * jm[8]))));
        case 5: // da1Nii6
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((x) * jm[5] * jm[6]) + ((-1.000000e+00) * (y) * jm[5] * jm[7]) + ((-1.000000e+00) * (x) * jm[3] * jm[8]) + ((y) * jm[4] * jm[8]));
        case 6: // da1Nii7
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * ((-1.000000e+00) + (x) + (z)) * ((jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[8]))) + ((y) * (((-2.000000e+00) * jm[5] * jm[6]) + ((-1.000000e+00) * jm[3] * jm[7]) + (jm[5] * jm[7]) + (jm[4] * (jm[6] + ((-1.000000e+00) * jm[8]))) + ((2.000000e+00) * jm[3] * jm[8]))));
        case 7: // da1Nii8
            return (4.000000e+00) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[5] * (z) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[3] * ((((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[7]) + ((-1.000000e+00) * (z) * jm[8]))) + (jm[4] * (((-1.000000e+00) * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[6]) + ((z) * jm[8]))));
        case 8: // da1Nii9
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * (x) * jm[4] * jm[6]) + ((x) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[5] * (z) * jm[7]) + (jm[4] * (z) * jm[8]));
        case 9: // da1Nii10
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * (y) * jm[4] * jm[6]) + (jm[5] * (z) * jm[6]) + ((y) * jm[3] * jm[7]) + ((-1.000000e+00) * jm[3] * (z) * jm[8]));
        case 10: // da2Nii1
            return ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y)) + ((4.000000e+00) * (z))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * (jm[7] + ((-1.000000e+00) * jm[8]))) + (jm[1] * (((-1.000000e+00) * jm[6]) + jm[8])));
        case 11: // da2Nii2
            return ((-1.000000e+00) + ((4.000000e+00) * (x))) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 12: // da2Nii3
            return ((-1.000000e+00) + ((4.000000e+00) * (y))) * ((jm[2] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[8])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 13: // da2Nii4
            return ((-1.000000e+00) + ((4.000000e+00) * (z))) * ((jm[1] * jm[6]) + ((-1.000000e+00) * jm[0] * jm[7])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 14: // da2Nii5
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((-1.000000e+00) * ((-1.000000e+00) + (y) + (z)) * ((jm[2] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[8]))) + ((x) * (((-1.000000e+00) * jm[1] * jm[6]) + (jm[2] * jm[6]) + (jm[0] * jm[7]) + ((-2.000000e+00) * jm[2] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[8]) + ((2.000000e+00) * jm[1] * jm[8]))));
        case 15: // da2Nii6
            return (4.000000e+00) * (((-1.000000e+00) * (x) * jm[2] * jm[6]) + (jm[2] * (y) * jm[7]) + ((x) * jm[0] * jm[8]) + ((-1.000000e+00) * jm[1] * (y) * jm[8])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 16: // da2Nii7
            return (4.000000e+00) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8]))) * ((jm[2] * (((-1.000000e+00) * ((-1.000000e+00) + (x) + ((2.000000e+00) * (y)) + (z)) * jm[6]) + ((y) * jm[7]))) + (jm[1] * (y) * (jm[6] + ((-1.000000e+00) * jm[8]))) + (jm[0] * (((-1.000000e+00) * (y) * (jm[7] + ((-2.000000e+00) * jm[8]))) + (((-1.000000e+00) + (x) + (z)) * jm[8]))));
        case 17: // da2Nii8
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[2] * (z) * (jm[6] + ((-1.000000e+00) * jm[7]))) + (jm[0] * ((((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[7]) + ((-1.000000e+00) * (z) * jm[8]))) + (jm[1] * (((-1.000000e+00) * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))) * jm[6]) + ((z) * jm[8]))));
        case 18: // da2Nii9
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * (((x) * jm[1] * jm[6]) + ((-1.000000e+00) * (x) * jm[0] * jm[7]) + (jm[2] * (z) * jm[7]) + ((-1.000000e+00) * jm[1] * (z) * jm[8]));
        case 19: // da2Nii10
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8]))) * ((jm[1] * (y) * jm[6]) + ((-1.000000e+00) * jm[2] * (z) * jm[6]) + ((-1.000000e+00) * jm[0] * (y) * jm[7]) + (jm[0] * (z) * jm[8]));
        case 20: // da3Nii1
            return ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4]))) + (jm[0] * (jm[4] + ((-1.000000e+00) * jm[5]))) + (jm[1] * (((-1.000000e+00) * jm[3]) + jm[5]))) * ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y)) + ((4.000000e+00) * (z))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 21: // da3Nii2
            return ((-1.000000e+00) + ((4.000000e+00) * (x))) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 22: // da3Nii3
            return ((-1.000000e+00) + ((4.000000e+00) * (y))) * ((jm[2] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 23: // da3Nii4
            return ((jm[1] * jm[3]) + ((-1.000000e+00) * jm[0] * jm[4])) * ((-1.000000e+00) + ((4.000000e+00) * (z))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 24: // da3Nii5
            return (4.000000e+00) * (((x) * (((-1.000000e+00) * jm[1] * jm[3]) + (jm[2] * jm[3]) + (jm[0] * jm[4]) + ((-2.000000e+00) * jm[2] * jm[4]) + ((-1.000000e+00) * jm[0] * jm[5]) + ((2.000000e+00) * jm[1] * jm[5]))) + ((-1.000000e+00) * ((jm[2] * jm[4]) + ((-1.000000e+00) * jm[1] * jm[5])) * ((-1.000000e+00) + (y) + (z)))) * (1.0 / ((jm[2] * jm[4] * jm[6]) + ((-1.000000e+00) * jm[1] * jm[5] * jm[6]) + ((-1.000000e+00) * jm[2] * jm[3] * jm[7]) + (jm[0] * jm[5] * jm[7]) + (jm[1] * jm[3] * jm[8]) + ((-1.000000e+00) * jm[0] * jm[4] * jm[8])));
        case 25: // da3Nii6
            return (4.000000e+00) * (((x) * jm[2] * jm[3]) + ((-1.000000e+00) * jm[2] * (y) * jm[4]) + ((-1.000000e+00) * (x) * jm[0] * jm[5]) + (jm[1] * (y) * jm[5])) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 26: // da3Nii7
            return (4.000000e+00) * ((jm[1] * (y) * (jm[3] + ((-1.000000e+00) * jm[5]))) + (jm[0] * (((-1.000000e+00) * (y) * (jm[4] + ((-2.000000e+00) * jm[5]))) + (jm[5] * ((-1.000000e+00) + (x) + (z))))) + (jm[2] * (((y) * jm[4]) + ((-1.000000e+00) * jm[3] * ((-1.000000e+00) + (x) + ((2.000000e+00) * (y)) + (z)))))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 27: // da3Nii8
            return (-4.000000e+00) * ((jm[2] * (jm[3] + ((-1.000000e+00) * jm[4])) * (z)) + (jm[1] * ((jm[5] * (z)) + ((-1.000000e+00) * jm[3] * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z)))))) + (jm[0] * (((-1.000000e+00) * jm[5] * (z)) + (jm[4] * ((-1.000000e+00) + (x) + (y) + ((2.000000e+00) * (z))))))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 28: // da3Nii9
            return (4.000000e+00) * (((-1.000000e+00) * (x) * jm[1] * jm[3]) + ((x) * jm[0] * jm[4]) + ((-1.000000e+00) * jm[2] * jm[4] * (z)) + (jm[1] * jm[5] * (z))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
        case 29: // da3Nii10
            return (4.000000e+00) * (((-1.000000e+00) * jm[1] * (y) * jm[3]) + (jm[0] * (y) * jm[4]) + (jm[2] * jm[3] * (z)) + ((-1.000000e+00) * jm[0] * jm[5] * (z))) * (1.0 / (((-1.000000e+00) * jm[2] * jm[4] * jm[6]) + (jm[1] * jm[5] * jm[6]) + (jm[2] * jm[3] * jm[7]) + ((-1.000000e+00) * jm[0] * jm[5] * jm[7]) + ((-1.000000e+00) * jm[1] * jm[3] * jm[8]) + (jm[0] * jm[4] * jm[8])));
		*/	
		default:
			TDKP_GENERAL_EXCEPTION("invalid shape function / diffop");            

	}
}

void Element3DTetrahedron2nd::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	
	local_coords.resize(3);
	switch(lid) {
		case 0:
			local_coords[0] = 0.0;
			local_coords[1] = 0.0;
			local_coords[2] = 0.0; 
			break;
		case 1:
			local_coords[0] = 1.0;
			local_coords[1] = 0.0;
			local_coords[2] = 0.0;
			break;						
		case 2:
			local_coords[0] = 0.0;
			local_coords[1] = 1.0;
			local_coords[2] = 0.0; 
			break;
		case 3:
			local_coords[0] = 0.0;
			local_coords[1] = 0.0;
			local_coords[2] = 1.0; 
			break;					
		case 4:
			local_coords[0] = 0.5;
			local_coords[1] = 0.0;
			local_coords[2] = 0.0; 
			break;
		case 5:
			local_coords[0] = 0.5;
			local_coords[1] = 0.5;
			local_coords[2] = 0.0; 
			break;
		case 6:
			local_coords[0] = 0.0;
			local_coords[1] = 0.5;
			local_coords[2] = 0.0; 
			break;
		case 7:
			local_coords[0] = 0.0;
			local_coords[1] = 0.0;
			local_coords[2] = 0.5; 
			break;									
		case 8:
			local_coords[0] = 0.5;
			local_coords[1] = 0.0;
			local_coords[2] = 0.5; 
			break;
		case 9:
			local_coords[0] = 0.0;
			local_coords[1] = 0.5;
			local_coords[2] = 0.5; 
			break;			
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
}




}
