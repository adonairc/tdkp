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

#include "Element1DLine2nd.h"

namespace tdkp {

/** analytically evalulated second order integrals */
const double Element1DLine2nd::reference_element_integral_second_order[3][3] = {
	{7.0  / 3.0, -8.0 / 3.0,  1.0 / 3.0},
	{-8.0 / 3.0, 16.0 / 3.0, -8.0 / 3.0},
	{1.0  / 3.0, -8.0 / 3.0,  7.0 / 3.0} 	
};
/** analytically evalulated first order integrals */
const double Element1DLine2nd::reference_element_integral_first_order[3][3] = {
	{-1.0 / 2.0, -2.0 / 3.0,  1.0 / 6.0},
	{2.0  / 3.0,        0.0, -2.0 / 3.0},
	{-1.0 / 6.0,  2.0 / 3.0,  1.0 / 2.0}
};
/** analytically evalulated zero order integrals */
const double Element1DLine2nd::reference_element_integral_zero_order[3][3] = {
	{2.0 / 15.0,  1.0 / 15.0, - 1.0 / 30.0},
	{1.0 / 15.0,  8.0 / 15.0,   1.0 / 15.0},	
	{-1.0 / 30.0, 1.0 / 15.0,   2.0 / 15.0}
};	
/** analytically evaluated zero order integrals multiplied by nodal element function */
const double Element1DLine2nd::reference_element_integral_zero_order_nodal_data[3][3][3] = {
	{
		{ 13.0 / 140.0,   1.0 /  21.0, - 1.0 / 140.0},
		{  1.0 /  21.0,   4.0 / 105.0, - 2.0 / 105.0},
		{- 1.0 / 140.0, - 2.0 / 105.0, - 1.0 / 140.0} 	
	},
	{
		{  1.0 /  21.0,   4.0 / 105.0, - 2.0 / 105.0},
		{  4.0 / 105.0,  16.0 /  35.0,   4.0 / 105.0},
		{- 2.0 / 105.0,   4.0 / 105.0,   1.0 /  21.0} 
	},
	{
		{- 1.0 / 140.0, - 2.0 / 105.0, - 1.0 / 140.0},
		{- 2.0 / 105.0,   4.0 / 105.0, 	 1.0 /  21.0},
		{- 1.0 / 140.0,   1.0 /  21.0,  13.0 / 140.0}
	}	
};
/** analytically evalulated single first order integrals */
const double Element1DLine2nd::reference_element_single_first_order[3] = {
	-1.0, 0.0, 1.0	
};
/** analytically evalulated single zero order integrals */
const double Element1DLine2nd::reference_element_single_zero_order[3] = {
	1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0		
};

Element1DLine2nd::Element1DLine2nd(unsigned int index_global_)
: Element1DLineBase(1,3,2)
{
	this->index_global  = index_global_;
	this->ready		    = false;
}
Element1DLine2nd::~Element1DLine2nd() {
	
}

bool Element1DLine2nd::verify() const {
	return Element::verify();	
}

void Element1DLine2nd::set_corner_node(unsigned short corner_id, Node* node) {
	TDKP_BOUNDS_ASSERT(corner_id < 2, "");
	if(corner_id == 0) {
		this->nodes[0] = node;	
	} else {
		this->nodes[2] = node;	
	}			
}
const Node& Element1DLine2nd::get_corner_node(unsigned short corner_idx) const {
	TDKP_BOUNDS_ASSERT(corner_idx < 2, "");
	if(corner_idx == 0) {
		return *this->nodes[0];	
	} else {
		return *this->nodes[2];	
	}		
}

/** get location information for additional nodes
 * 
 * higher order elements need additional nodes. sometimes, these nodes sit
 * on edges or planes and are shared between elements. in order to identify 
 * these nodes, a location scheme is applied.
 * 
 * there are three stages of the location. the first stage is the type
 * of location (inner / vertex / element boundary) and the indicies of the
 * vertices of vertex / edge / plane where the node is located.
 * 
 * the second step is the coordinate of the node that is passed. so if an
 * edge has e.g. two nodes, the nodes can be distinguished according to their
 * coordinates.
 * 
 * the third step is a tag. so, if there are multiple nodes on the same spot,
 * they can get an additional tag that can be used for their identification.
 * 
 * this is particularly necessary for hermite elements, where additional nodes
 * are located at the vertex location but correspond to the derivatives of
 * the form-function
 * 
 * @param additional_node_idx (input) index of additional node (its number 0 < num_additional_nodes)
 * @param location_type (output) type of the location where the node is sitting on (inner, vertex, boundary)
 * @param involved vertices (output) vertex indices of edge / face / point where the node is sitting
 * @param coords (output) coordinates of the node
 * @param tag (output) additional tag for the given node 
 */
void Element1DLine2nd::get_additional_node_locator(
	unsigned int additional_node_idx, 
	AdditionalNodeLocation& location_type, 
	vector<unsigned int>& involved_vertices, 
	vector<double>& coords, 
	unsigned int& tag
) const {

	TDKP_BOUNDS_ASSERT(additional_node_idx == 0,"");
	TDKP_BOUNDS_ASSERT(nodes[0] != 0,"");
	TDKP_BOUNDS_ASSERT(nodes[2] != 0,"");
	// inner node -> belongs only to that element
	location_type = inner_node;
	// involved vertices? none
	involved_vertices.resize(0);	
	// set coordinate for additional node (at element mid point)
	coords.resize(1);
	coords[0] = 0.5 * (nodes[0]->get_coord(0) + nodes[2]->get_coord(0));
	// tag is zero (not a special node)
	tag = 0; 			
}

void Element1DLine2nd::set_additional_node(unsigned int additional_node_idx, Node* node) {
	TDKP_BOUNDS_ASSERT(additional_node_idx == 0, "");
	TDKP_BOUNDS_ASSERT(test_additional_node(*this, *node, additional_node_idx), "wrong additional node supplied: " << additional_node_idx);
	nodes[1] = node;
	// -----------------------------------
	// set linear interpolation coefficents for additional node
	// -----------------------------------
	if(node->get_num_contributions() == 0) {
		// node sits in middle of element, so halve the value of each vertex
		node->set_contribution(nodes[0]->get_index_vertex(), 0.5);
		node->set_contribution(nodes[2]->get_index_vertex(), 0.5);	
	}
}

// ----------------------------------------------
// evaluation functions
// ----------------------------------------------
double  Element1DLine2nd::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 3 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 3, "");
	return reference_element_integral_second_order[elem_shape_func_1][elem_shape_func_2]
	       / this->element_volume;	
}

double  Element1DLine2nd::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 3 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 3, "element shape function index " << elem_shape_func_1 << "/" << elem_shape_func_2 << " is invalid");
 	return reference_element_integral_first_order[elem_shape_func_1][elem_shape_func_2];	
}

double  Element1DLine2nd::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 3 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 3, "element shape function index is invalid");	
	return this->element_volume * reference_element_integral_zero_order[elem_shape_func_1][elem_shape_func_2];	
}
									
double  Element1DLine2nd::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 3, "element shape function index is invalid");	
	TDKP_BOUNDS_ASSERT(this->nodes[1]->get_coord(0) - this->nodes[0]->get_coord(0) > 0.0, "nodes are not ordered properly!");	 
	return reference_element_single_first_order[elem_shape_func_1];		
}
double  Element1DLine2nd::get_single_integral_0th_order(short elem_shape_func_1) const {
	return reference_element_single_zero_order[elem_shape_func_1] * this->element_volume;	
}

double Element1DLine2nd::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(nodal_data_point >= 0 && nodal_data_point < 3 && elem_shape_func_1 >= 0 && elem_shape_func_1 < 3 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 3, "");
	return this->element_volume * reference_element_integral_zero_order_nodal_data[nodal_data_point][elem_shape_func_1][elem_shape_func_2];
}
													

/** function to evaluate the form function inside the reference element */	
double Element1DLine2nd::evaluate_form_function(
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	const double& x = local_reference_element_coords[0];
	switch(elem_shape_func) {
		case 0:
			return 2.0 * x * x - 3.0 * x + 1.0;
		case 1:
			return - 4.0 * x * x + 4.0 * x;
		case 2:
			return 2.0 * x * x - x;
		default:
			TDKP_GENERAL_EXCEPTION("invalid form function index");
	}	
}

double Element1DLine2nd::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	
	const double& x = local_reference_element_coords[0];
	// d/dx -> 1/h d/de
	switch(elem_shape_func) {
		case 0:
			return (1.0 / element_volume) * (4.0 * x - 3.0);
		case 1:
			return (1.0 / element_volume) * (-8.0 * x + 4.0);
		case 2:
			return (1.0 / element_volume) * (4.0 * x - 1.0);
		default:
			TDKP_GENERAL_EXCEPTION("invalid form function index");
	}	
}

void Element1DLine2nd::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	local_coords.resize(1);
	switch(lid) {
		case 0:
			local_coords[0] = 0.0; 
			break;			
		case 1:
			local_coords[0] = 0.5;
			break;
		case 2:
			local_coords[0] = 1.0;
			break;			
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
}


}
