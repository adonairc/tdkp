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

#include "Element1DHermiteCubic.h"

namespace tdkp {

Element1DHermiteCubic::Element1DHermiteCubic(unsigned int index_global_) 
: Element1DLineBase(1,4,2)
{
	this->index_global  = index_global_;
	this->ready		    = false;
}

Element1DHermiteCubic::~Element1DHermiteCubic() {
	
}
	 
void Element1DHermiteCubic::get_node_local_coords(
	unsigned short lid, 
	vector<double>& local_coords
) const {
	local_coords.resize(1);
	if(lid < 2) {
		local_coords[0] = 0.0;
	} else {
		local_coords[0] = 1.0;
	}
}
	
void Element1DHermiteCubic::set_corner_node(unsigned short corner_id, Node* node) {	
	if(corner_id == 0) {
		TDKP_BOUNDS_ASSERT(this->nodes[0] == 0, "");
		this->nodes[0] = node;	
	} else if (corner_id  == 1) {
		TDKP_BOUNDS_ASSERT(this->nodes[2] == 0, "");
		this->nodes[2] = node;		
	} else {
		TDKP_GENERAL_EXCEPTION("invalid corner");	
	}	
}

const Node& Element1DHermiteCubic::get_corner_node(unsigned short corner_idx) const {
	if(corner_idx == 0) {
		TDKP_BOUNDS_ASSERT(this->nodes[0] != 0, "");
		return *this->nodes[0];	
	} else if(corner_idx == 1) {
		TDKP_BOUNDS_ASSERT(this->nodes[2] != 0, "");
		return *this->nodes[2];	
	} else {
		TDKP_GENERAL_EXCEPTION("invalid corner");		
	}
}

void Element1DHermiteCubic::get_additional_node_locator(
	unsigned int additional_node_idx, AdditionalNodeLocation& location_type, 
	vector<unsigned int>& involved_vertices, vector<double>& coords, unsigned int& tag
) const {
	// additional nodes sit at vertices!
	if(additional_node_idx < 2) {
		coords.resize(1);
		coords[0]     = get_corner_node(additional_node_idx).get_coord(0);
		tag           = 0;
		location_type = vertex_node;
		involved_vertices.resize(1);
		involved_vertices[0] = get_corner_node(additional_node_idx).get_index_vertex();
	} else {
		TDKP_GENERAL_EXCEPTION("invalid additional node index");	
	}	
}

void Element1DHermiteCubic::set_additional_node(unsigned int additional_node_idx, Node* node) {
	if(additional_node_idx < 2) {
		unsigned int target_idx = additional_node_idx * 2 + 1;
		TDKP_BOUNDS_ASSERT(this->nodes[target_idx] == 0, "");	
		this->nodes[target_idx] = node;
		node->set_value_type(Node::Derivative_X);		
	} else {
		TDKP_GENERAL_EXCEPTION("invalid additional node index");
	}	
}

double Element1DHermiteCubic::get_element_integral_2nd_order(
	short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2
) const {
	unsigned int obj_identifier = elem_shape_func_2
	                            + 4 * elem_shape_func_1;
	
	switch(obj_identifier) {
        case 10: // da1db1Nii3Njj3
        case 0: // da1db1Nii1Njj1
            return ((6.000000e+00) / (5.000000e+00)) * (1.0 / this->element_volume);
        case 3: // da1db1Nii1Njj4
        case 4: // da1db1Nii2Njj1
        case 12: // da1db1Nii4Njj1
        case 1: // da1db1Nii1Njj2
            return ((1.000000e+00) / (1.000000e+01)) * (1.0 / this->element_volume);
        case 8: // da1db1Nii3Njj1
        case 2: // da1db1Nii1Njj3
            return ((-6.000000e+00) / (5.000000e+00)) * (1.0 / this->element_volume);
        case 15: // da1db1Nii4Njj4
        case 5: // da1db1Nii2Njj2
            return ((2.000000e+00) / (1.500000e+01)) * (1.0 / this->element_volume);
        case 9: // da1db1Nii3Njj2
        case 11: // da1db1Nii3Njj4
        case 14: // da1db1Nii4Njj3
        case 6: // da1db1Nii2Njj3
            return ((-1.000000e+00) / (1.000000e+01)) * (1.0 / this->element_volume);
        case 13: // da1db1Nii4Njj2
        case 7: // da1db1Nii2Njj4
            return ((-1.000000e+00) / (3.000000e+01)) * (1.0 / this->element_volume);			
		default:
			TDKP_GENERAL_EXCEPTION("invalid obj identifier");
	}	
}
double Element1DHermiteCubic::get_element_integral_1st_order(
	short diffop, short elem_shape_func_1, short elem_shape_func_2
) const {
	unsigned int obj_identifier = elem_shape_func_2
	                            + 4 * elem_shape_func_1
	                            + 16 * diffop;
	
	switch(obj_identifier) {

        case 2: // da1Nii1Njj3
        case 0: // da1Nii1Njj1
            return (-1.000000e+00) / (2.000000e+00);
        case 6: // da1Nii2Njj3
        case 11: // da1Nii3Njj4
        case 12: // da1Nii4Njj1
        case 1: // da1Nii1Njj2
            return (-1.000000e+00) / (1.000000e+01);
        case 4: // da1Nii2Njj1
        case 9: // da1Nii3Njj2
        case 14: // da1Nii4Njj3
        case 3: // da1Nii1Njj4
            return (1.000000e+00) / (1.000000e+01);
        case 7: // da1Nii2Njj4
            return (1.000000e+00) / (6.000000e+01);
        case 10: // da1Nii3Njj3
        case 8: // da1Nii3Njj1
            return (1.000000e+00) / (2.000000e+00);
        case 13: // da1Nii4Njj2
            return (-1.000000e+00) / (6.000000e+01);
				
		default:
			return 0.0;
			
	}	
}
double Element1DHermiteCubic::get_element_integral_0th_order(
	short elem_shape_func_1, short elem_shape_func_2
) const {
	
	
	unsigned int obj_identifier = elem_shape_func_2
	                            + 4 * elem_shape_func_1;
	
	switch(obj_identifier) {
			
        case 10: // Nii3Njj3
        case 0: // Nii1Njj1
            return ((1.300000e+01) / (3.500000e+01)) * this->element_volume;
        case 4: // Nii2Njj1
        case 1: // Nii1Njj2
            return ((1.100000e+01) / (2.100000e+02)) * this->element_volume;
        case 8: // Nii3Njj1
        case 2: // Nii1Njj3
            return ((9.000000e+00) / (7.000000e+01)) * this->element_volume;
        case 12: // Nii4Njj1
        case 3: // Nii1Njj4
            return ((-1.300000e+01) / (4.200000e+02)) * this->element_volume;
        case 15: // Nii4Njj4
        case 5: // Nii2Njj2
            return ((1.000000e+00) / (1.050000e+02)) * this->element_volume;
        case 9: // Nii3Njj2
        case 6: // Nii2Njj3
            return ((1.300000e+01) / (4.200000e+02)) * this->element_volume;
        case 13: // Nii4Njj2
        case 7: // Nii2Njj4
            return ((-1.000000e+00) / (1.400000e+02)) * this->element_volume;
        case 14: // Nii4Njj3
        case 11: // Nii3Njj4
            return ((-1.100000e+01) / (2.100000e+02)) * this->element_volume;
		default:
			TDKP_GENERAL_EXCEPTION("invalid obj identifier");
	}            
            
		            	
}
double Element1DHermiteCubic::get_single_integral_1st_order(
	short diffop, short elem_shape_func_1
) const {

	unsigned int obj_identifier = elem_shape_func_1
	                            + 4 * diffop;
	
	switch(obj_identifier) {
		
		case 0:
			return - 1.0;
		case 1: 
		case 3:
			return 0;
		case 2:
			return 1.0;				
		default:
			TDKP_GENERAL_EXCEPTION("invalid obj identifier");
	}
	
}


double Element1DHermiteCubic::get_single_integral_0th_order(short elem_shape_func_1) const {
	
	unsigned int obj_identifier = elem_shape_func_1;
	
	switch(obj_identifier) {
        case 2: // Nii3
        case 0: // Nii1
            return ((1.000000e+00) / (2.000000e+00)) * this->element_volume;
        case 1: // Nii2
            return ((1.000000e+00) / (1.200000e+01)) * this->element_volume;
        case 3: // Nii4
            return ((-1.000000e+00) / (1.200000e+01)) * this->element_volume;		
		default:
			TDKP_GENERAL_EXCEPTION("invalid obj identifier");
	}
}




double Element1DHermiteCubic::get_element_integral_0th_order_nodal_data(
	short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2
) const {
	
	unsigned int obj_identifier = elem_shape_func_2
	                            + 4 * elem_shape_func_1
	                            + 16 * nodal_data_point;
	
	switch(obj_identifier) {
	
        case 42: // da3Nii3Njj3
        case 0: // da1Nii1Njj1
            return ((4.300000e+01) / (1.400000e+02)) * this->element_volume;
        case 4: // da1Nii2Njj1
        case 16: // da2Nii1Njj1
        case 1: // da1Nii1Njj2
            return ((9.700000e+01) / (2.520000e+03)) * this->element_volume;
        case 8: // da1Nii3Njj1
        case 10: // da1Nii3Njj3
        case 32: // da3Nii1Njj1
        case 34: // da3Nii1Njj3
        case 40: // da3Nii3Njj1
        case 2: // da1Nii1Njj3
            return ((9.000000e+00) / (1.400000e+02)) * this->element_volume;
        case 12: // da1Nii4Njj1
        case 48: // da4Nii1Njj1
        case 3: // da1Nii1Njj4
            return ((-4.300000e+01) / (2.520000e+03)) * this->element_volume;
        case 17: // da2Nii1Njj2
        case 20: // da2Nii2Njj1
        case 47: // da3Nii4Njj4
        case 59: // da4Nii3Njj4
        case 62: // da4Nii4Njj3
        case 5: // da1Nii2Njj2
            return ((2.000000e+00) / (3.150000e+02)) * this->element_volume;
        case 9: // da1Nii3Njj2
        case 18: // da2Nii1Njj3
        case 24: // da2Nii3Njj1
        case 33: // da3Nii1Njj2
        case 36: // da3Nii2Njj1
        case 6: // da1Nii2Njj3
            return ((1.000000e+00) / (7.200000e+01)) * this->element_volume;
        case 13: // da1Nii4Njj2
        case 19: // da2Nii1Njj4
        case 27: // da2Nii3Njj4
        case 28: // da2Nii4Njj1
        case 30: // da2Nii4Njj3
        case 39: // da3Nii2Njj4
        case 45: // da3Nii4Njj2
        case 49: // da4Nii1Njj2
        case 52: // da4Nii2Njj1
        case 54: // da4Nii2Njj3
        case 57: // da4Nii3Njj2
        case 7: // da1Nii2Njj4
            return ((-1.000000e+00) / (2.800000e+02)) * this->element_volume;
        case 14: // da1Nii4Njj3
        case 35: // da3Nii1Njj4
        case 44: // da3Nii4Njj1
        case 50: // da4Nii1Njj3
        case 56: // da4Nii3Njj1
        case 11: // da1Nii3Njj4
            return ((-1.000000e+00) / (7.200000e+01)) * this->element_volume;
        case 22: // da2Nii2Njj3
        case 25: // da2Nii3Njj2
        case 37: // da3Nii2Njj2
        case 51: // da4Nii1Njj4
        case 60: // da4Nii4Njj1
        case 15: // da1Nii4Njj4
            return ((1.000000e+00) / (3.150000e+02)) * this->element_volume;
        case 21: // da2Nii2Njj2
            return ((1.000000e+00) / (8.400000e+02)) * this->element_volume;
        case 29: // da2Nii4Njj2
        case 53: // da4Nii2Njj2
        case 23: // da2Nii2Njj4
            return ((-1.000000e+00) / (1.260000e+03)) * this->element_volume;
        case 38: // da3Nii2Njj3
        case 41: // da3Nii3Njj2
        case 26: // da2Nii3Njj3
            return ((4.300000e+01) / (2.520000e+03)) * this->element_volume;
        case 55: // da4Nii2Njj4
        case 61: // da4Nii4Njj2
        case 31: // da2Nii4Njj4
            return ((1.000000e+00) / (1.260000e+03)) * this->element_volume;
        case 46: // da3Nii4Njj3
        case 58: // da4Nii3Njj3
        case 43: // da3Nii3Njj4
            return ((-9.700000e+01) / (2.520000e+03)) * this->element_volume;
        case 63: // da4Nii4Njj4
            return ((-1.000000e+00) / (8.400000e+02)) * this->element_volume;
            
		default:
			TDKP_GENERAL_EXCEPTION("invalid obj identifier");            
	}
}
		
double Element1DHermiteCubic::evaluate_form_function(
	short elem_shape_func, const double* local_reference_element_coords
) const {
	const double& x = local_reference_element_coords[0];
	switch(elem_shape_func) {
		case 0:
			return 2.0 * x * x * x - 3.0 * x *x + 1.0;
 		case 1: 
 			return 1.0 * x * x * x - 2.0 * x * x + 1.0 * x; 
		case 2: 			
 			return - 2.0 * x * x * x + 3.0 * x * x;
 		case 3:
 			return 1.0 * x * x * x - 1.0 * x *x;
 		default:
 			TDKP_GENERAL_EXCEPTION("invalid shape function requested");
	}
}

double Element1DHermiteCubic::evaluate_form_function_derivative(
	short diffop, short elem_shape_func, const double* local_reference_element_coords
) const {
	const double& x = local_reference_element_coords[0];
	TDKP_BOUNDS_ASSERT(diffop == 0, "");
	switch(elem_shape_func) {
		case 0:
			return (1.0 / this->element_volume) * (6.0 * (x - 1.0) * x);
		case 1:
			return (1.0 / this->element_volume) * (1.0 - 4.0 * x + 3.0 * x * x);
		case 2:
			return (1.0 / this->element_volume) * (- 6.0 * (x - 1.0) * x);
		case 3:
			return (1.0 / this->element_volume) * (x * (3.0 * x - 2.0));
 		default:
 			TDKP_GENERAL_EXCEPTION("invalid shape function requested");		
	}
}		
	


}
