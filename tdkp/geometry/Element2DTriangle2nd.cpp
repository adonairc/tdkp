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

#include "Element2DTriangle2nd.h"

namespace tdkp {

// ---------------------------------------------------
// analytical integrals were obtained via mathematica and
// some python based parsing using some python scripts
// ---------------------------------------------------



/** 2d second order triangle element constructor */
Element2DTriangle2nd::Element2DTriangle2nd(unsigned int index_)
: Element2DTriangleBase(2,6,3) 
{
	this->index_global = index_;
	for(int ii = 0; ii < 4; ii++) {
		this->inverted_jacobi_matrix[ii] = 0;			
	}			
}

Element2DTriangle2nd::~Element2DTriangle2nd() {
	this->region = 0;		
}

/** set corner vertices (nodes 0, 2 and 4) */		
void Element2DTriangle2nd::set_corner_node(unsigned short corner_id, Node* node) {
	TDKP_BOUNDS_ASSERT(corner_id * 2 < (unsigned short)get_num_nodes(), "");
	nodes[corner_id * 2] = node;	
}

const Node& Element2DTriangle2nd::get_corner_node(unsigned short corner_idx) const {
	TDKP_BOUNDS_ASSERT(corner_idx * 2 < (unsigned short)get_num_nodes(), "");
	return *nodes[corner_idx * 2];	
}

/** get locator for nodes 1, 3, 4. see comment for 1D 2nd order  */
void Element2DTriangle2nd::get_additional_node_locator(
	unsigned int additional_node_idx, 
	AdditionalNodeLocation& location_type, 
	vector<unsigned int>& involved_vertices, 
	vector<double>& coords, 
	unsigned int& tag
) const {

	TDKP_ASSERT(nodes[0] != 0 && nodes[2] != 0 && nodes[4] != 0, "");		
	TDKP_ASSERT(additional_node_idx < 3, "");
	
	tag = 0; // no special tag required for lagrange elements
	location_type = edge_node; // nodes sit on the boundary

	int nd0 = additional_node_idx * 2;
	int nd1 = (additional_node_idx * 2 + 2) % 6;

	involved_vertices.resize(2);
	involved_vertices[0] = nodes[nd0]->get_index_vertex();
	involved_vertices[1] = nodes[nd1]->get_index_vertex();
		
	// ------------------------------------------
	// set coordinates
	// ------------------------------------------		
	coords.assign(2, 0.0);
	for(unsigned int ii = 0; ii < 2; ii++) {
		coords[ii] = (nodes[nd0]->get_coord(ii) + nodes[nd1]->get_coord(ii)) / 2.0;  	
	} 
		
}

void Element2DTriangle2nd::set_additional_node(unsigned int additional_node_idx, Node* node) {
	
	TDKP_ASSERT(additional_node_idx < 3, "");	
	TDKP_BOUNDS_ASSERT(additional_node_idx * 2 + 1 < get_num_nodes(), "");
	TDKP_BOUNDS_ASSERT(nodes[additional_node_idx * 2 + 1] == 0, "");
	// ------------------------------------------
	// check node position
	// ------------------------------------------
	int nd0 = additional_node_idx * 2;
	int nd1 = (additional_node_idx * 2 + 2) % 6;
	for(unsigned int ii = 0; ii < 2; ii++) {		
		double t = 0.5 * (nodes[nd0]->get_coord(ii) + nodes[nd1]->get_coord(ii));
		TDKP_ASSERT(tdkp_math::abs(t - node->get_coord(ii)) < 1.0e-6," ii = " << ii << " t = " << t);
	}	 		
	nodes[additional_node_idx * 2 + 1] = node;

	// -----------------------------------
	// set linear interpolation coefficents for additional node
	// -----------------------------------
	if(node->get_num_contributions() == 0) {
		int lidx = additional_node_idx * 2 + 1;
		// node sits in middle of edge, so halve the value of each vertex
		node->set_contribution(nodes[lidx - 1]->get_index_vertex(), 0.5);
		node->set_contribution(nodes[(lidx + 1) % 6]->get_index_vertex(), 0.5);	
	}		
				
}		
		
/** returns analytical second order integrals */
double Element2DTriangle2nd::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {
	
	TDKP_BOUNDS_ASSERT(diffop_1 >= 0 && diffop_1 < 2 && diffop_2 >= 0 && diffop_2 < 2, "diffop_1 >= 0 && diffop_1 < 2 && diffop_2 >= 0 && diffop_2 < 2");
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 6 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 6, "");

	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 6 * elem_shape_func_1
	                       + 6 * 6 * diffop_2
	                       + 6 * 6 * 2 * diffop_1;
	
	switch(obj_index) {
		// ---------------------------------------------
	    //  integration of all second order terms 
    	// ---------------------------------------------
        case 0: // da1db1Nii1Njj1
            return ((1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * ((jm[2] + ((-1.000000e+00) * jm[3])) * (jm[2] + ((-1.000000e+00) * jm[3])));
        case 6: // da1db1Nii2Njj1
        case 1: // da1db1Nii1Njj2
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (jm[2] + ((-1.000000e+00) * jm[3])) * jm[3];
        case 12: // da1db1Nii3Njj1
        case 2: // da1db1Nii1Njj3
            return ((-1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * (jm[2] + ((-1.000000e+00) * jm[3])) * jm[3];
        case 24: // da1db1Nii5Njj1
        case 4: // da1db1Nii1Njj5
            return ((1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[2] * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 30: // da1db1Nii6Njj1
        case 5: // da1db1Nii1Njj6
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[2] * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 21: // da1db1Nii4Njj4
        case 35: // da1db1Nii6Njj6
        case 7: // da1db1Nii2Njj2
            return ((4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * ((jm[2] * jm[2]) + ((-1.000000e+00) * jm[2] * jm[3]) + (jm[3] * jm[3]));
        case 13: // da1db1Nii3Njj2
        case 8: // da1db1Nii2Njj3
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[3] * (((-1.000000e+00) * jm[2]) + jm[3]);
        case 19: // da1db1Nii4Njj2
        case 9: // da1db1Nii2Njj4
            return ((4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[2] * (((-1.000000e+00) * jm[2]) + jm[3]);
        case 31: // da1db1Nii6Njj2
        case 11: // da1db1Nii2Njj6
            return ((-4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[2] * jm[3];
        case 14: // da1db1Nii3Njj3
            return ((1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * (jm[3] * jm[3]);
        case 20: // da1db1Nii4Njj3
        case 22: // da1db1Nii4Njj5
        case 27: // da1db1Nii5Njj4
        case 15: // da1db1Nii3Njj4
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[2] * jm[3];
        case 26: // da1db1Nii5Njj3
        case 16: // da1db1Nii3Njj5
            return ((1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[2] * jm[3];
        case 33: // da1db1Nii6Njj4
        case 23: // da1db1Nii4Njj6
            return ((4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (jm[2] + ((-1.000000e+00) * jm[3])) * jm[3];
        case 28: // da1db1Nii5Njj5
            return ((1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * (jm[2] * jm[2]);
        case 34: // da1db1Nii6Njj5
        case 29: // da1db1Nii5Njj6
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[2] * (((-1.000000e+00) * jm[2]) + jm[3]);
        case 72: // da2db1Nii1Njj1
        case 36: // da1db2Nii1Njj1
            return ((-1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 78: // da2db1Nii2Njj1
        case 37: // da1db2Nii1Njj2
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[1] * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 84: // da2db1Nii3Njj1
        case 38: // da1db2Nii1Njj3
            return ((1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[1] * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 96: // da2db1Nii5Njj1
        case 40: // da1db2Nii1Njj5
            return ((-1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[0] * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 70: // da1db2Nii6Njj5
        case 101: // da2db1Nii5Njj6
        case 102: // da2db1Nii6Njj1
        case 41: // da1db2Nii1Njj6
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[0] * (jm[2] + ((-1.000000e+00) * jm[3]));
        case 73: // da2db1Nii1Njj2
        case 42: // da1db2Nii2Njj1
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[3];
        case 57: // da1db2Nii4Njj4
        case 71: // da1db2Nii6Njj6
        case 79: // da2db1Nii2Njj2
        case 93: // da2db1Nii4Njj4
        case 107: // da2db1Nii6Njj6
        case 43: // da1db2Nii2Njj2
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * ((jm[1] * (jm[2] + ((-2.000000e+00) * jm[3]))) + (jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 85: // da2db1Nii3Njj2
        case 44: // da1db2Nii2Njj3
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[1] * (((-1.000000e+00) * jm[2]) + jm[3]);
        case 91: // da2db1Nii4Njj2
        case 45: // da1db2Nii2Njj4
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * (((-2.000000e+00) * jm[2]) + jm[3])));
        case 67: // da1db2Nii6Njj2
        case 83: // da2db1Nii2Njj6
        case 103: // da2db1Nii6Njj2
        case 47: // da1db2Nii2Njj6
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * ((jm[1] * jm[2]) + (jm[0] * jm[3]));
        case 74: // da2db1Nii1Njj3
        case 48: // da1db2Nii3Njj1
            return ((1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[3];
        case 80: // da2db1Nii2Njj3
        case 49: // da1db2Nii3Njj2
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (((-1.000000e+00) * jm[0]) + jm[1]) * jm[3];
        case 86: // da2db1Nii3Njj3
        case 50: // da1db2Nii3Njj3
            return ((-1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * jm[1] * jm[3];
        case 58: // da1db2Nii4Njj5
        case 92: // da2db1Nii4Njj3
        case 99: // da2db1Nii5Njj4
        case 51: // da1db2Nii3Njj4
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[0] * jm[3];
        case 98: // da2db1Nii5Njj3
        case 52: // da1db2Nii3Njj5
            return ((-1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[0] * jm[3];
        case 81: // da2db1Nii2Njj4
        case 55: // da1db2Nii4Njj2
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (((2.000000e+00) * jm[0] * jm[2]) + ((-1.000000e+00) * jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]));
        case 63: // da1db2Nii5Njj4
        case 87: // da2db1Nii3Njj4
        case 94: // da2db1Nii4Njj5
        case 56: // da1db2Nii4Njj3
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[1] * jm[2];
        case 69: // da1db2Nii6Njj4
        case 95: // da2db1Nii4Njj6
        case 105: // da2db1Nii6Njj4
        case 59: // da1db2Nii4Njj6
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * ((jm[1] * (jm[2] + ((-2.000000e+00) * jm[3]))) + (jm[0] * jm[3]));
        case 76: // da2db1Nii1Njj5
        case 60: // da1db2Nii5Njj1
            return ((-1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[2];
        case 88: // da2db1Nii3Njj5
        case 62: // da1db2Nii5Njj3
            return ((-1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[1] * jm[2];
        case 100: // da2db1Nii5Njj5
        case 64: // da1db2Nii5Njj5
            return ((-1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * jm[0] * jm[2];
        case 66: // da1db2Nii6Njj1
        case 77: // da2db1Nii1Njj6
        case 106: // da2db1Nii6Njj5
        case 65: // da1db2Nii5Njj6
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[2];
        case 108: // da2db2Nii1Njj1
            return ((1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * ((jm[0] + ((-1.000000e+00) * jm[1])) * (jm[0] + ((-1.000000e+00) * jm[1])));
        case 114: // da2db2Nii2Njj1
        case 109: // da2db2Nii1Njj2
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[1];
        case 120: // da2db2Nii3Njj1
        case 110: // da2db2Nii1Njj3
            return ((-1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[1];
        case 132: // da2db2Nii5Njj1
        case 112: // da2db2Nii1Njj5
            return ((1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[0] * (jm[0] + ((-1.000000e+00) * jm[1]));
        case 138: // da2db2Nii6Njj1
        case 113: // da2db2Nii1Njj6
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[0] * (jm[0] + ((-1.000000e+00) * jm[1]));
        case 129: // da2db2Nii4Njj4
        case 143: // da2db2Nii6Njj6
        case 115: // da2db2Nii2Njj2
            return ((4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * ((jm[0] * jm[0]) + ((-1.000000e+00) * jm[0] * jm[1]) + (jm[1] * jm[1]));
        case 121: // da2db2Nii3Njj2
        case 116: // da2db2Nii2Njj3
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[1] * (((-1.000000e+00) * jm[0]) + jm[1]);
        case 127: // da2db2Nii4Njj2
        case 117: // da2db2Nii2Njj4
            return ((-4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[0] * (jm[0] + ((-1.000000e+00) * jm[1]));
        case 139: // da2db2Nii6Njj2
        case 119: // da2db2Nii2Njj6
            return ((-4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[0] * jm[1];
        case 122: // da2db2Nii3Njj3
            return ((1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * (jm[1] * jm[1]);
        case 128: // da2db2Nii4Njj3
        case 130: // da2db2Nii4Njj5
        case 135: // da2db2Nii5Njj4
        case 123: // da2db2Nii3Njj4
            return ((-2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[0] * jm[1];
        case 134: // da2db2Nii5Njj3
        case 124: // da2db2Nii3Njj5
            return ((1.000000e+00) / (6.000000e+00)) * (1.0 / jacobi_det) * jm[0] * jm[1];
        case 141: // da2db2Nii6Njj4
        case 131: // da2db2Nii4Njj6
            return ((4.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * (jm[0] + ((-1.000000e+00) * jm[1])) * jm[1];
        case 136: // da2db2Nii5Njj5
            return ((1.000000e+00) / (2.000000e+00)) * (1.0 / jacobi_det) * (jm[0] * jm[0]);
        case 142: // da2db2Nii6Njj5
        case 137: // da2db2Nii5Njj6
            return ((2.000000e+00) / (3.000000e+00)) * (1.0 / jacobi_det) * jm[0] * (((-1.000000e+00) * jm[0]) + jm[1]);

        default:
			return 0.0; 
	};

}

void Element2DTriangle2nd::print() const {
	Element::print();
	ostringstream sout;
	sout << "jm:\n"
	     << "  " << jm[0] << " " << jm[1] << "\n"
	     << "  " << jm[2] << " " << jm[3] << "\n";	
	TDKP_LOGMSG(LOG_INFO, sout.str());
}

double Element2DTriangle2nd::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {
	
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 6 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 6, "");
	TDKP_BOUNDS_ASSERT(diffop >= 0 && diffop < 2, "diffop >= 0 && diffop < 2");

	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2
	                       + 6 * elem_shape_func_1
	                       + 6 * 6 * diffop;

	switch(obj_index) {

        case 0: // da1Nii1Njj1
            return ((1.000000e+00) / (1.500000e+01)) * jacobi_det * (((-1.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 5: // da1Nii1Njj6
        case 1: // da1Nii1Njj2
            return ((1.000000e+00) / (1.000000e+01)) * jacobi_det * (((-1.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 3: // da1Nii1Njj4
        case 4: // da1Nii1Njj5
        case 2: // da1Nii1Njj3
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 6: // da1Nii2Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (jm[2] + ((-3.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 7: // da1Nii2Njj2
            return ((4.000000e+00) / (1.500000e+01)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 8: // da1Nii2Njj3
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((2.000000e+00) * jm[2]) + ((-3.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 9: // da1Nii2Njj4
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (((-2.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 10: // da1Nii2Njj5
            return jacobi_det * jm[2] * (1.0 / (((3.000000e+01) * jm[1] * jm[2]) + ((-3.000000e+01) * jm[0] * jm[3])));
        case 11: // da1Nii2Njj6
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[2] + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 16: // da1Nii3Njj5
        case 17: // da1Nii3Njj6
        case 12: // da1Nii3Njj1
            return jacobi_det * jm[3] * (1.0 / (((3.000000e+01) * jm[1] * jm[2]) + ((-3.000000e+01) * jm[0] * jm[3])));
        case 15: // da1Nii3Njj4
        case 13: // da1Nii3Njj2
            return ((1.000000e+00) / (1.000000e+01)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 14: // da1Nii3Njj3
            return ((1.000000e+00) / (1.500000e+01)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 18: // da1Nii4Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((-1.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 19: // da1Nii4Njj2
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (((-2.000000e+00) * jm[2]) + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 20: // da1Nii4Njj3
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((2.000000e+00) * jm[2]) + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 21: // da1Nii4Njj4
            return ((4.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 22: // da1Nii4Njj5
            return ((-1.000000e+00) / (3.000000e+01)) * jacobi_det * (jm[2] + ((2.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 23: // da1Nii4Njj6
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[2] + ((-2.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 25: // da1Nii5Njj2
        case 26: // da1Nii5Njj3
        case 24: // da1Nii5Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * jm[2] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 29: // da1Nii5Njj6
        case 27: // da1Nii5Njj4
            return jacobi_det * jm[2] * (1.0 / (((1.000000e+01) * jm[1] * jm[2]) + ((-1.000000e+01) * jm[0] * jm[3])));
        case 28: // da1Nii5Njj5
            return jacobi_det * jm[2] * (1.0 / (((1.500000e+01) * jm[1] * jm[2]) + ((-1.500000e+01) * jm[0] * jm[3])));
        case 30: // da1Nii6Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((-3.000000e+00) * jm[2]) + jm[3]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 31: // da1Nii6Njj2
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[2] + jm[3]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 32: // da1Nii6Njj3
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * jm[3] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 33: // da1Nii6Njj4
            return ((-2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[2] + ((-2.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 34: // da1Nii6Njj5
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((3.000000e+00) * jm[2]) + ((-2.000000e+00) * jm[3])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 35: // da1Nii6Njj6
            return (4.000000e+00) * jacobi_det * jm[3] * (1.0 / (((1.500000e+01) * jm[1] * jm[2]) + ((-1.500000e+01) * jm[0] * jm[3])));
        case 36: // da2Nii1Njj1
            return ((1.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 41: // da2Nii1Njj6
        case 37: // da2Nii1Njj2
            return ((1.000000e+00) / (1.000000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 38: // da2Nii1Njj3
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 40: // da2Nii1Njj5
        case 39: // da2Nii1Njj4
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((-1.000000e+00) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 42: // da2Nii2Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (jm[0] + ((-3.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 43: // da2Nii2Njj2
            return (4.000000e+00) * jacobi_det * jm[0] * (1.0 / (((1.500000e+01) * jm[1] * jm[2]) + ((-1.500000e+01) * jm[0] * jm[3])));
        case 44: // da2Nii2Njj3
            return (((2.000000e+00) * jacobi_det * jm[0]) + ((-3.000000e+00) * jacobi_det * jm[1])) * (1.0 / (((3.000000e+01) * jm[1] * jm[2]) + ((-3.000000e+01) * jm[0] * jm[3])));
        case 45: // da2Nii2Njj4
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (((-2.000000e+00) * jm[0]) + jm[1]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 46: // da2Nii2Njj5
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 47: // da2Nii2Njj6
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[0] + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 52: // da2Nii3Njj5
        case 53: // da2Nii3Njj6
        case 48: // da2Nii3Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 51: // da2Nii3Njj4
        case 49: // da2Nii3Njj2
            return jacobi_det * jm[1] * (1.0 / (((1.000000e+01) * jm[1] * jm[2]) + ((-1.000000e+01) * jm[0] * jm[3])));
        case 50: // da2Nii3Njj3
            return jacobi_det * jm[1] * (1.0 / (((1.500000e+01) * jm[1] * jm[2]) + ((-1.500000e+01) * jm[0] * jm[3])));
        case 54: // da2Nii4Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 55: // da2Nii4Njj2
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (((-2.000000e+00) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 56: // da2Nii4Njj3
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((2.000000e+00) * jm[0]) + jm[1]) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 57: // da2Nii4Njj4
            return ((-4.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 58: // da2Nii4Njj5
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (jm[0] + ((2.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 59: // da2Nii4Njj6
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[0] + ((-2.000000e+00) * jm[1])) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 61: // da2Nii5Njj2
        case 62: // da2Nii5Njj3
        case 60: // da2Nii5Njj1
            return jacobi_det * jm[0] * (1.0 / (((3.000000e+01) * jm[1] * jm[2]) + ((-3.000000e+01) * jm[0] * jm[3])));
        case 65: // da2Nii5Njj6
        case 63: // da2Nii5Njj4
            return ((1.000000e+00) / (1.000000e+01)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 64: // da2Nii5Njj5
            return ((1.000000e+00) / (1.500000e+01)) * jacobi_det * jm[0] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 66: // da2Nii6Njj1
            return ((1.000000e+00) / (3.000000e+01)) * jacobi_det * (((-3.000000e+00) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 67: // da2Nii6Njj2
            return ((-2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[0] + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 68: // da2Nii6Njj3
            return jacobi_det * jm[1] * (1.0 / (((3.000000e+01) * jm[1] * jm[2]) + ((-3.000000e+01) * jm[0] * jm[3])));
        case 69: // da2Nii6Njj4
            return ((2.000000e+00) / (1.500000e+01)) * jacobi_det * (jm[0] + ((-2.000000e+00) * jm[1])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 70: // da2Nii6Njj5
            return (((3.000000e+00) * jacobi_det * jm[0]) + ((-2.000000e+00) * jacobi_det * jm[1])) * (1.0 / (((3.000000e+01) * jm[1] * jm[2]) + ((-3.000000e+01) * jm[0] * jm[3])));
        case 71: // da2Nii6Njj6
            return ((4.000000e+00) / (1.500000e+01)) * jacobi_det * jm[1] * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));

		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point!");
	}

}

double Element2DTriangle2nd::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < (signed)get_num_nodes() && elem_shape_func_2 >= 0 && elem_shape_func_2 < (signed)get_num_nodes(),	
		"elem_shape_func1 >= 0 && elem_shape_func1 < nnode && elem_shape_func2 >= 0 && elem_shape_func2 < nnode");

	// ----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_2 
	                       + 6 * elem_shape_func_1;

	switch(obj_index) {

        case 14: // Nii3Njj3
        case 28: // Nii5Njj5
        case 0: // Nii1Njj1
            return ((1.000000e+00) / (6.000000e+01)) * jacobi_det;
        case 4: // Nii1Njj5
        case 12: // Nii3Njj1
        case 16: // Nii3Njj5
        case 24: // Nii5Njj1
        case 26: // Nii5Njj3
        case 2: // Nii1Njj3
            return ((-1.000000e+00) / (3.600000e+02)) * jacobi_det;
        case 10: // Nii2Njj5
        case 17: // Nii3Njj6
        case 18: // Nii4Njj1
        case 25: // Nii5Njj2
        case 32: // Nii6Njj3
        case 3: // Nii1Njj4
            return ((-1.000000e+00) / (9.000000e+01)) * jacobi_det;
        case 21: // Nii4Njj4
        case 35: // Nii6Njj6
        case 7: // Nii2Njj2
            return ((4.000000e+00) / (4.500000e+01)) * jacobi_det;
        case 11: // Nii2Njj6
        case 19: // Nii4Njj2
        case 23: // Nii4Njj6
        case 31: // Nii6Njj2
        case 33: // Nii6Njj4
        case 9: // Nii2Njj4
            return ((2.000000e+00) / (4.500000e+01)) * jacobi_det;
        default:
			return 0.0e0; // some (unhandled ones) are 0
	}

}


double Element2DTriangle2nd::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	
	unsigned int obj_index = elem_shape_func_2 
	                       + 6 * elem_shape_func_1
	                       + 36 * nodal_data_point;
	                       
	switch(obj_index) {
        case 86: // Ncc3Nii3Njj3
        case 172: // Ncc5Nii5Njj5
        case 0: // Ncc1Nii1Njj1
            return ((1.000000e+00) / (1.400000e+02)) * jacobi_det;
        case 5: // Ncc1Nii1Njj6
        case 6: // Ncc1Nii2Njj1
        case 30: // Ncc1Nii6Njj1
        case 36: // Ncc2Nii1Njj1
        case 50: // Ncc2Nii3Njj3
        case 80: // Ncc3Nii2Njj3
        case 85: // Ncc3Nii3Njj2
        case 87: // Ncc3Nii3Njj4
        case 92: // Ncc3Nii4Njj3
        case 122: // Ncc4Nii3Njj3
        case 136: // Ncc4Nii5Njj5
        case 166: // Ncc5Nii4Njj5
        case 171: // Ncc5Nii5Njj4
        case 173: // Ncc5Nii5Njj6
        case 178: // Ncc5Nii6Njj5
        case 180: // Ncc6Nii1Njj1
        case 208: // Ncc6Nii5Njj5
        case 1: // Ncc1Nii1Njj2
            return ((1.000000e+00) / (2.100000e+02)) * jacobi_det;
        case 4: // Ncc1Nii1Njj5
        case 12: // Ncc1Nii3Njj1
        case 14: // Ncc1Nii3Njj3
        case 24: // Ncc1Nii5Njj1
        case 28: // Ncc1Nii5Njj5
        case 72: // Ncc3Nii1Njj1
        case 74: // Ncc3Nii1Njj3
        case 84: // Ncc3Nii3Njj1
        case 88: // Ncc3Nii3Njj5
        case 98: // Ncc3Nii5Njj3
        case 100: // Ncc3Nii5Njj5
        case 144: // Ncc5Nii1Njj1
        case 148: // Ncc5Nii1Njj5
        case 158: // Ncc5Nii3Njj3
        case 160: // Ncc5Nii3Njj5
        case 168: // Ncc5Nii5Njj1
        case 170: // Ncc5Nii5Njj3
        case 2: // Ncc1Nii1Njj3
            return ((-1.000000e+00) / (1.260000e+03)) * jacobi_det;
        case 18: // Ncc1Nii4Njj1
        case 64: // Ncc2Nii5Njj5
        case 89: // Ncc3Nii3Njj6
        case 104: // Ncc3Nii6Njj3
        case 108: // Ncc4Nii1Njj1
        case 154: // Ncc5Nii2Njj5
        case 169: // Ncc5Nii5Njj2
        case 194: // Ncc6Nii3Njj3
        case 3: // Ncc1Nii1Njj4
            return ((1.000000e+00) / (6.300000e+02)) * jacobi_det;
        case 13: // Ncc1Nii3Njj2
        case 29: // Ncc1Nii5Njj6
        case 34: // Ncc1Nii6Njj5
        case 38: // Ncc2Nii1Njj3
        case 48: // Ncc2Nii3Njj1
        case 73: // Ncc3Nii1Njj2
        case 78: // Ncc3Nii2Njj1
        case 94: // Ncc3Nii4Njj5
        case 99: // Ncc3Nii5Njj4
        case 124: // Ncc4Nii3Njj5
        case 134: // Ncc4Nii5Njj3
        case 149: // Ncc5Nii1Njj6
        case 159: // Ncc5Nii3Njj4
        case 164: // Ncc5Nii4Njj3
        case 174: // Ncc5Nii6Njj1
        case 184: // Ncc6Nii1Njj5
        case 204: // Ncc6Nii5Njj1
        case 8: // Ncc1Nii2Njj3
            return ((-1.000000e+00) / (6.300000e+02)) * jacobi_det;
        case 19: // Ncc1Nii4Njj2
        case 23: // Ncc1Nii4Njj6
        case 33: // Ncc1Nii6Njj4
        case 39: // Ncc2Nii1Njj4
        case 53: // Ncc2Nii3Njj6
        case 54: // Ncc2Nii4Njj1
        case 58: // Ncc2Nii4Njj5
        case 63: // Ncc2Nii5Njj4
        case 65: // Ncc2Nii5Njj6
        case 68: // Ncc2Nii6Njj3
        case 70: // Ncc2Nii6Njj5
        case 83: // Ncc3Nii2Njj6
        case 95: // Ncc3Nii4Njj6
        case 103: // Ncc3Nii6Njj2
        case 105: // Ncc3Nii6Njj4
        case 109: // Ncc4Nii1Njj2
        case 113: // Ncc4Nii1Njj6
        case 114: // Ncc4Nii2Njj1
        case 118: // Ncc4Nii2Njj5
        case 125: // Ncc4Nii3Njj6
        case 133: // Ncc4Nii5Njj2
        case 138: // Ncc4Nii6Njj1
        case 140: // Ncc4Nii6Njj3
        case 153: // Ncc5Nii2Njj4
        case 155: // Ncc5Nii2Njj6
        case 163: // Ncc5Nii4Njj2
        case 175: // Ncc5Nii6Njj2
        case 183: // Ncc6Nii1Njj4
        case 188: // Ncc6Nii2Njj3
        case 190: // Ncc6Nii2Njj5
        case 193: // Ncc6Nii3Njj2
        case 195: // Ncc6Nii3Njj4
        case 198: // Ncc6Nii4Njj1
        case 200: // Ncc6Nii4Njj3
        case 205: // Ncc6Nii5Njj2
        case 9: // Ncc1Nii2Njj4
            return ((-1.000000e+00) / (3.150000e+02)) * jacobi_det;
        case 26: // Ncc1Nii5Njj3
        case 76: // Ncc3Nii1Njj5
        case 96: // Ncc3Nii5Njj1
        case 146: // Ncc5Nii1Njj3
        case 156: // Ncc5Nii3Njj1
        case 16: // Ncc1Nii3Njj5
            return ((1.000000e+00) / (2.520000e+03)) * jacobi_det;
        case 46: // Ncc2Nii2Njj5
        case 61: // Ncc2Nii5Njj2
        case 107: // Ncc3Nii6Njj6
        case 111: // Ncc4Nii1Njj4
        case 126: // Ncc4Nii4Njj1
        case 151: // Ncc5Nii2Njj2
        case 197: // Ncc6Nii3Njj6
        case 212: // Ncc6Nii6Njj3
        case 21: // Ncc1Nii4Njj4
            return ((-2.000000e+00) / (3.150000e+02)) * jacobi_det;
        case 129: // Ncc4Nii4Njj4
        case 215: // Ncc6Nii6Njj6
        case 43: // Ncc2Nii2Njj2
            return ((2.000000e+00) / (3.500000e+01)) * jacobi_det;
        case 47: // Ncc2Nii2Njj6
        case 55: // Ncc2Nii4Njj2
        case 57: // Ncc2Nii4Njj4
        case 67: // Ncc2Nii6Njj2
        case 71: // Ncc2Nii6Njj6
        case 115: // Ncc4Nii2Njj2
        case 117: // Ncc4Nii2Njj4
        case 127: // Ncc4Nii4Njj2
        case 131: // Ncc4Nii4Njj6
        case 141: // Ncc4Nii6Njj4
        case 143: // Ncc4Nii6Njj6
        case 187: // Ncc6Nii2Njj2
        case 191: // Ncc6Nii2Njj6
        case 201: // Ncc6Nii4Njj4
        case 203: // Ncc6Nii4Njj6
        case 211: // Ncc6Nii6Njj2
        case 213: // Ncc6Nii6Njj4
        case 45: // Ncc2Nii2Njj4
            return ((2.000000e+00) / (1.050000e+02)) * jacobi_det;
        case 69: // Ncc2Nii6Njj4
        case 119: // Ncc4Nii2Njj6
        case 139: // Ncc4Nii6Njj2
        case 189: // Ncc6Nii2Njj4
        case 199: // Ncc6Nii4Njj2
        case 59: // Ncc2Nii4Njj6
            return ((4.000000e+00) / (3.150000e+02)) * jacobi_det;
		default:
			return 0.0;		
		
	}
	
}	                       
					
					
double Element2DTriangle2nd::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 6, "");
	TDKP_BOUNDS_ASSERT(diffop >= 0 && diffop < 2, "diffop >= 0 && diffop < 2");

	// -----------------------------------------------
	// build object identifier
	// ----------------------------------------------
	unsigned int obj_index = elem_shape_func_1 
	                       + 6 * diffop;

	switch(obj_index) {

        case 0: // da1Nii1
            return jacobi_det * (((-1.000000e+00) * jm[2]) + jm[3]) * (1.0 / (((6.000000e+00) * jm[1] * jm[2]) + ((-6.000000e+00) * jm[0] * jm[3])));
        case 1: // da1Nii2
            return (2.000000e+00) * jacobi_det * jm[2] * (1.0 / (((-3.000000e+00) * jm[1] * jm[2]) + ((3.000000e+00) * jm[0] * jm[3])));
        case 2: // da1Nii3
            return jacobi_det * jm[3] * (1.0 / (((-6.000000e+00) * jm[1] * jm[2]) + ((6.000000e+00) * jm[0] * jm[3])));
        case 3: // da1Nii4
            return (2.000000e+00) * jacobi_det * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / (((3.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3])));
        case 4: // da1Nii5
            return jacobi_det * jm[2] * (1.0 / (((6.000000e+00) * jm[1] * jm[2]) + ((-6.000000e+00) * jm[0] * jm[3])));
        case 5: // da1Nii6
            return (2.000000e+00) * jacobi_det * jm[3] * (1.0 / (((3.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3])));
        case 6: // da2Nii1
            return jacobi_det * (jm[0] + ((-1.000000e+00) * jm[1])) * (1.0 / (((6.000000e+00) * jm[1] * jm[2]) + ((-6.000000e+00) * jm[0] * jm[3])));
        case 7: // da2Nii2
            return (2.000000e+00) * jacobi_det * jm[0] * (1.0 / (((3.000000e+00) * jm[1] * jm[2]) + ((-3.000000e+00) * jm[0] * jm[3])));
        case 8: // da2Nii3
            return jacobi_det * jm[1] * (1.0 / (((6.000000e+00) * jm[1] * jm[2]) + ((-6.000000e+00) * jm[0] * jm[3])));
        case 9: // da2Nii4
            return ((2.000000e+00) / (3.000000e+00)) * jacobi_det * (((-1.000000e+00) * jm[0]) + jm[1]) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 10: // da2Nii5
            return jacobi_det * jm[0] * (1.0 / (((-6.000000e+00) * jm[1] * jm[2]) + ((6.000000e+00) * jm[0] * jm[3])));
        case 11: // da2Nii6
            return (2.000000e+00) * jacobi_det * jm[1] * (1.0 / (((-3.000000e+00) * jm[1] * jm[2]) + ((3.000000e+00) * jm[0] * jm[3])));
	        
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point!");
	}
	
}

double Element2DTriangle2nd::get_single_integral_0th_order(short elem_shape_func_1) const {

	switch(elem_shape_func_1) {
        case 3: // Nii4
        case 5: // Nii6
        case 1: // Nii2
            return ((1.000000e+00) / (6.000000e+00)) * jacobi_det;
		case 0:
		case 2:
		case 4:
			return 0.0;            
	 	default:
	 		TDKP_GENERAL_EXCEPTION("must not reach that point!");
	}

}
	
double Element2DTriangle2nd::evaluate_form_function(
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	
	const double& x = local_reference_element_coords[0];
	const double& y = local_reference_element_coords[1];
	switch(elem_shape_func) {
		case 0: // N_0 = 2*x^2 + 2*y^2 + 4 * x * y - 3* x - 3*y + 1
			return 2.0 * x * x + 2.0 * y * y + 4.0 * x * y - 3.0 * x - 3.0 * y + 1.0;
		case 1: // N_1 = - 4* x^2 - 4 * x * y + 4 * x
			return - 4.0 * x * x - 4.0 * x * y + 4.0 * x;								
		case 2: // N_2 = 2 * x^2 - x
			return 2.0 * x * x - x;
		case 3: // N_3 = 4 * x * y
			return 4.0 * x * y;				
		case 4: // N_4 = 2 * y^2 - y
			return 2.0 * y * y - y;
		case 5: // N_5 = -4*y^2 - 4 * x * y + 4 * y
			return - 4.0 * y * y - 4.0 * x * y + 4.0 * y;
		default:
			TDKP_GENERAL_EXCEPTION("this element has only three elem shape functions!");				
	}	
}

double Element2DTriangle2nd::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local
) const {
	return evaluate_form_function_derivative(diffop, elem_shape_func, local[0], local[1]);
}

/** calculate derivative (global) at local coordinates */
double Element2DTriangle2nd::evaluate_form_function_derivative(short diffop, short elem_shape_func, const double& x, const double& y) const {

	unsigned int obj_idx = elem_shape_func + 6 * diffop;
	
	switch(obj_idx) {
		
        case 0: // da1Nii1
            return ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y))) * (jm[2] + ((-1.000000e+00) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 1: // da1Nii2
            return (4.000000e+00) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * (((x) * (jm[2] + ((-2.000000e+00) * jm[3]))) + jm[3] + ((-1.000000e+00) * (y) * jm[3]));
        case 2: // da1Nii3
            return (jm[3] + ((-4.000000e+00) * (x) * jm[3])) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 3: // da1Nii4
            return (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3]))) * (((4.000000e+00) * (x) * jm[2]) + ((-4.000000e+00) * (y) * jm[3]));
        case 4: // da1Nii5
            return ((-1.000000e+00) + ((4.000000e+00) * (y))) * jm[2] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 5: // da1Nii6
            return (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3]))) * (((4.000000e+00) * ((-1.000000e+00) + (x) + ((2.000000e+00) * (y))) * jm[2]) + ((-4.000000e+00) * (y) * jm[3]));
        case 6: // da2Nii1
            return (jm[0] + ((-1.000000e+00) * jm[1])) * ((-3.000000e+00) + ((4.000000e+00) * (x)) + ((4.000000e+00) * (y))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 7: // da2Nii2
            return (-4.000000e+00) * (((x) * (jm[0] + ((-2.000000e+00) * jm[1]))) + jm[1] + ((-1.000000e+00) * jm[1] * (y))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 8: // da2Nii3
            return ((-1.000000e+00) + ((4.000000e+00) * (x))) * jm[1] * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 9: // da2Nii4
            return (4.000000e+00) * (((x) * jm[0]) + ((-1.000000e+00) * jm[1] * (y))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
        case 10: // da2Nii5
            return (jm[0] + ((-4.000000e+00) * jm[0] * (y))) * (1.0 / ((jm[1] * jm[2]) + ((-1.000000e+00) * jm[0] * jm[3])));
        case 11: // da2Nii6
            return (4.000000e+00) * ((jm[1] * (y)) + ((-1.000000e+00) * jm[0] * ((-1.000000e+00) + (x) + ((2.000000e+00) * (y))))) * (1.0 / (((-1.000000e+00) * jm[1] * jm[2]) + (jm[0] * jm[3])));
		 
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");
	}
}

void Element2DTriangle2nd::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	local_coords.resize(2);
	switch(lid) {
		case 0:
			local_coords[0] = 0.0;
			local_coords[1] = 0.0; 
			break;
		case 1:
			local_coords[0] = 0.5;
			local_coords[1] = 0.0; 
			break;						
		case 2:
			local_coords[0] = 1.0;
			local_coords[1] = 0.0; 
			break;
		case 3:
			local_coords[0] = 0.5;
			local_coords[1] = 0.5; 
			break;					
		case 4:
			local_coords[0] = 0.0;
			local_coords[1] = 1.0; 
			break;
		case 5:
			local_coords[0] = 0.0;
			local_coords[1] = 0.5; 
			break;				
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
}




}
