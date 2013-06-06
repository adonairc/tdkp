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

#include "tdkp/geometry/Element1DLine.h"

namespace tdkp {

Element1DLineBase::Element1DLineBase(unsigned int dimension, unsigned int nnodes, unsigned int nboundaries)
: Element(dimension, nnodes, nboundaries)
{
}

void Element1DLineBase::get_bounding_box(Node& low, Node& high) const {
	if(get_corner_node(0).get_coord(0) < get_corner_node(1).get_coord(0)) {
		low  = get_corner_node(0);	
		high = get_corner_node(1);
	} else {
		low  = get_corner_node(1);	
		high = get_corner_node(0);		
	}
}	 	
bool Element1DLineBase::inside_element(const double& global_x, const double& global_y, const double& global_z) const {
	
	double tol = 1.0e-10;
	double x0 = get_corner_node(0).get_coord(0);
	double x1 = get_corner_node(1).get_coord(0);	
	if(x0 < x1 && global_x + tol >= x0 && global_x - tol <= x1) {
		return true;	
	} else if(x0 >= x1 && global_x + tol >= x1 && global_x - tol <= x0) {
		return true;
	}
	return false;
}
double Element1DLineBase::evaluate_form_function_global(
	short elem_shape_func, 
	const double& global_x, 
	const double& global_y, 
	const double& global_z
) const {
	
	double local_x = (global_x - get_corner_node(0).get_coord(0)) 
	               / (get_corner_node(1).get_coord(0) - get_corner_node(0).get_coord(0)); 
	return evaluate_form_function(elem_shape_func, &local_x);	
} 

/** prepare element for calculation
 * 
 * in that case, calculate element size (wuhuuu)
 */
void Element1DLineBase::prepare() {
	TDKP_ASSERT(this->region != 0, "");
	this->element_mid_point[0] = 0.5;
	this->element_volume = get_corner_node(1).get_coord(0) - get_corner_node(0).get_coord(0);
	this->jacobi_det = this->element_volume; 
	this->ready = true;	
	TDKP_ASSERT(this->element_volume > 0, "sorry, but i need a properly ordered 1D mesh (ascending node coordinate).");
} 

void Element1DLineBase::get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const {
	vertex_indices.resize(1);
	if(idx == 0) {
		vertex_indices[0] = get_corner_node(0).get_index_vertex();			
	} else if(idx == 1) {	
		vertex_indices[0] = get_corner_node(1).get_index_vertex();			
	} else {
		TDKP_GENERAL_EXCEPTION("invalid boundary index ");
	}	
}

/** return integration points for integration according to simpsons rule */ 
DomainMaster Element1DLineBase::get_numerical_integration_points(unsigned int npoints) const {
	DomainMaster domain;
	create_1D_domain_simpson(domain, 0, 1, npoints);
	domain.freeze();
	return domain;
}

const double Element1DLine::reference_element_integral_second_order[2][2] = {
	{1.0, -1.0}, {-1.0, 1.0}
};
const double Element1DLine::reference_element_integral_first_order[2][2] = {
	{-1.0/2.0, -1.0/2.0}, {1.0/2.0, 1.0/2.0}
};
const double Element1DLine::reference_element_integral_zero_order[2][2] = {
	{1.0/3.0, 1.0/6.0}, {1.0/6.0, 1.0/3.0}	
};
const double Element1DLine::reference_element_integral_zero_order_nodal[2][2][2] = {
	{
		{1.0 / 4.0,  1.0 / 12.0},
		{1.0 / 12.0, 1.0 / 12.0}
	},
	{
		{1.0 / 12.0, 1.0 / 12.0},
		{1.0 / 12.0, 1.0 /  4.0}
	}	
};

const double Element1DLine::reference_element_single_first_order[2] = {-1.0, 1.0};

Element1DLine::Element1DLine(unsigned int index_global_)
: Element1DLineBase(1,2,2) 
{
	this->index_global  = index_global_;
	this->ready		    = false;
}

bool Element1DLine::verify() const {
	TDKP_ASSERT(region != 0, "no region set to element");
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		TDKP_ASSERT(nodes[ii] != 0, "node in element is null");
		for(unsigned int jj = ii + 1; jj < get_num_nodes(); jj++) {
			if(nodes[ii] == nodes[jj]) {
				Logger::get_instance()->emit(LOG_ERROR, 1000, "nodes %d and %d in element are equal", ii, jj);
				return false;	
			}
		}	
	}
	TDKP_ASSERT(this->element_volume != 0, "element volume is not calculated or is zero");
	return true;
}


double  Element1DLine::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 2 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 2, "element shape function index is invalid");
	return reference_element_integral_second_order[elem_shape_func_1][elem_shape_func_2]
	       / this->element_volume;   	
}

double  Element1DLine::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 2 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 2, "element shape function index " << elem_shape_func_1 << "/" << elem_shape_func_2 << " is invalid");
 	return reference_element_integral_first_order[elem_shape_func_1][elem_shape_func_2];	
}
double  Element1DLine::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 2 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 2, "element shape function index is invalid");	
	return this->element_volume * reference_element_integral_zero_order[elem_shape_func_1][elem_shape_func_2];	
}
double Element1DLine::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(nodal_data_point >= 0 && nodal_data_point < 2 && elem_shape_func_1 >= 0 && elem_shape_func_1 < 2 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 2, "");
	return this->element_volume * reference_element_integral_zero_order_nodal[nodal_data_point][elem_shape_func_1][elem_shape_func_2];		
}			

double  Element1DLine::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 2, "element shape function index is invalid");	
	TDKP_BOUNDS_ASSERT(this->nodes[1]->get_coord(0) - this->nodes[0]->get_coord(0) > 0.0, "nodes are not ordered properly!");
	 
	return reference_element_single_first_order[elem_shape_func_1];	
}
double  Element1DLine::get_single_integral_0th_order(short elem_shape_func_1) const {
	return 0.5 * this->element_volume;	
}
														
	
double Element1DLine::evaluate_form_function(
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	const double& x = local_reference_element_coords[0];
	if(elem_shape_func == 0) {
		return 1.0 - x;
	} else if(elem_shape_func == 1) {
		return x;	
	} else {
		TDKP_GENERAL_EXCEPTION("invalid form function index");	
	}
}

double Element1DLine::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	// d/dx_real = 1/h d/dx_reference
	if(elem_shape_func == 0) {
		return -1.0 / element_volume;	
	} else if(elem_shape_func == 1) {
		return  1.0 / element_volume;	
	} else {
		TDKP_GENERAL_EXCEPTION("invalid form function index");	
	}
}

void Element1DLine::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	local_coords.resize(1);
	switch(lid) {
		case 0:
			local_coords[0] = 0.0; 
			break;			
		case 1:
			local_coords[0] = 1.0;
			break;
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
}

} // end of namespace
