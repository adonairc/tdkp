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


#include "tdkp/geometry/Element2DTriangle.h"

namespace tdkp {

// ----------------------------------------------------------------
// analytical solutions of reference element integrations
// (obtained via mathematica and using paper and a pen
// ----------------------------------------------------------------
const double Element2DTriangle::reference_element_integral_zero[3][3] = {
	{1.0 / 12.0,  1.0 / 24.0, 1.0 / 24.0},
	{1.0 / 24.0,  1.0 / 12.0, 1.0 / 24.0},
	{1.0 / 24.0,  1.0 / 24.0, 1.0 / 12.0},		
};

const double Element2DTriangle::reference_element_gradients[3][2] = {
	{-1.0, -1.0},
	{ 1.0,  0.0},
	{ 0.0,  1.0}	
};


/** 2d linear rectangle constructor */
Element2DTriangle::Element2DTriangle(unsigned int index_)
: Element2DTriangleBase(2,3,3) 
{
	this->index_global    = index_;
	for(int ii = 0; ii < 4; ii++) {
		this->inverted_jacobi_matrix[ii] = 0;			
	}			
}

Element2DTriangle::~Element2DTriangle() {
	this->region = 0;		
}
		



double Element2DTriangle::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {
	
	TDKP_BOUNDS_ASSERT(diffop_1 >= 0 && diffop_1 < 2 && diffop_2 >= 0 && diffop_2 < 2, "diffop_1 >= 0 && diffop_1 < 2 && diffop_2 >= 0 && diffop_2 < 2");
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 4 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 4, "elem_shape_func_1 >= 0 && elem_shape_func_1 < 4 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 4");

	double left, right;
	left = right = 0;
	for(int ii = 0; ii < 2; ii++) {
		left  += inverted_jacobi_matrix[diffop_1 + ii * 2] * reference_element_gradients[elem_shape_func_1][ii];
		right += inverted_jacobi_matrix[diffop_2 + ii * 2] * reference_element_gradients[elem_shape_func_2][ii];
	}
	return (left * right * element_volume);

}

double Element2DTriangle::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {
	
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 3 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 3, "elem_shape_func_1 >= 0 && elem_shape_func_1 < 3 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 3");
	TDKP_BOUNDS_ASSERT(diffop >= 0 && diffop < 2, "diffop >= 0 && diffop < 2");
	double left = 0.0;
	for(int ii = 0; ii < 2; ii++) {
		left += inverted_jacobi_matrix[diffop + ii * 2] * reference_element_gradients[elem_shape_func_1][ii];
	}
	return jacobi_det * left * (1.0/6.0); 		
}

double Element2DTriangle::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < (signed)get_num_nodes() && elem_shape_func_2 >= 0 && elem_shape_func_2 < (signed)get_num_nodes(),	
		"elem_shape_func1 >= 0 && elem_shape_func1 < nnode && elem_shape_func2 >= 0 && elem_shape_func2 < nnode");
	return reference_element_integral_zero[elem_shape_func_1][elem_shape_func_2] * jacobi_det;
}
				
double Element2DTriangle::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	
	unsigned int obj_index = elem_shape_func_2 
	                       + 3 * elem_shape_func_1
	                       + 9 * nodal_data_point;
	switch(obj_index) {		
        case 13: // Ncc2Nii2Njj2
        case 26: // Ncc3Nii3Njj3
        case 0: // Ncc1Nii1Njj1
            return ((1.000000e+00) / (2.000000e+01)) * jacobi_det;
        case 2: // Ncc1Nii1Njj3
        case 3: // Ncc1Nii2Njj1
        case 4: // Ncc1Nii2Njj2
        case 6: // Ncc1Nii3Njj1
        case 8: // Ncc1Nii3Njj3
        case 9: // Ncc2Nii1Njj1
        case 10: // Ncc2Nii1Njj2
        case 12: // Ncc2Nii2Njj1
        case 14: // Ncc2Nii2Njj3
        case 16: // Ncc2Nii3Njj2
        case 17: // Ncc2Nii3Njj3
        case 18: // Ncc3Nii1Njj1
        case 20: // Ncc3Nii1Njj3
        case 22: // Ncc3Nii2Njj2
        case 23: // Ncc3Nii2Njj3
        case 24: // Ncc3Nii3Njj1
        case 25: // Ncc3Nii3Njj2
        case 1: // Ncc1Nii1Njj2
            return ((1.000000e+00) / (6.000000e+01)) * jacobi_det;
        case 7: // Ncc1Nii3Njj2
        case 11: // Ncc2Nii1Njj3
        case 15: // Ncc2Nii3Njj1
        case 19: // Ncc3Nii1Njj2
        case 21: // Ncc3Nii2Njj1
        case 5: // Ncc1Nii2Njj3
            return ((1.000000e+00) / (1.200000e+02)) * jacobi_det;
		default:            
			TDKP_GENERAL_EXCEPTION("invalid node data point / element shape function");		
	}	                       

}	                       				
					
					
double Element2DTriangle::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 3, "elem_shape_func_1 >= 0 && elem_shape_func_1 < 3 ");
	TDKP_BOUNDS_ASSERT(diffop >= 0 && diffop < 2, "diffop >= 0 && diffop < 2");
	double left = 0.0;
	for(int ii = 0; ii < 2; ii++) {
		left += inverted_jacobi_matrix[diffop + ii * 2] * reference_element_gradients[elem_shape_func_1][ii];
	}	
	return jacobi_det * left * 1.0 / 2.0;
}

double Element2DTriangle::get_single_integral_0th_order(short elem_shape_func_1) const {
	return jacobi_det * 1.0 / 6.0;
}

double Element2DTriangle::evaluate_form_function(short elem_shape_func_1, const double& x, const double& y, const double& z) const {
	switch(elem_shape_func_1) {
		case 0:
			return 1.0 - x - y;
		case 1:
			return x;
		case 2:
			return y;	
		default:
			TDKP_GENERAL_EXCEPTION("unknown shape function");
	}
}

		


ostream& operator<<(ostream& out, const Element2DTriangle &elem) {

	out << "2d triangle " << elem.get_index_global() << " region: " << elem.region->get_name() << "\nnodes: ";
	for(unsigned int ii = 0; ii < elem.get_num_nodes(); ii++) {
		out << " " << elem.get_node(ii);
	}
	return out;
}



/** function to evaluate the form function inside the reference element */	
double Element2DTriangle::evaluate_form_function(
	short elem_shape_func, 
	const double* local
) const {
	return evaluate_form_function(elem_shape_func, local[0], local[1], 0.0);	
}

/** function to evaluate the form function derivative inside the reference element */
double Element2DTriangle::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	// constant over element
	double ret = 0.0;
	for(unsigned int ii = 0; ii < 2; ii++) {
		ret += inverted_jacobi_matrix[diffop + ii * 2] * reference_element_gradients[elem_shape_func][ii];
	}
	return ret;		
}					


/** base constructor forward to element */
Element2DTriangleBase::Element2DTriangleBase(unsigned int dimension_, unsigned int nnodes_, unsigned int nboundaries_)
: Element(dimension_, nnodes_, nboundaries_)
{
	for(int ii = 0; ii < 4; ii++) {
		this->inverted_jacobi_matrix[ii] = 0.0;
		this->jm[0] = 0.0;			
	}		
}

/** test whether global coordinates x,y,z are inside the given element
 */
bool Element2DTriangleBase::inside_element(const double& global_x, const double& global_y, const double& global_z) const {
	double global[2];
	global[0] = global_x;
	global[1] = global_y;
	double local[] = {0.0, 0,0};
	return global2local(global,local);		
}

/** deterimines local coordinates
 * 
 * @param global array of x,y global coordiantes 
 * @param local array to store final local coordinates
 * @return bool if point x,y is inside element 
 */
bool Element2DTriangleBase::global2local(const double global[2], double local[2]) const {
	
	const double tol = 1.0e-5;
	
	const double* P0 = nodes[0]->get_coords();	
	double tmp_coords[2];
	for(short ii = 0; ii < 2; ii++) {
		tmp_coords[ii] = global[ii] - P0[ii];
		local[ii] = 0.0; 
	}
	// calc coordinates
	for(short ii = 0; ii < 2; ii++) {
		for(short jj = 0; jj < 2; jj++) {
			local[ii] += inverted_jacobi_matrix[ii * 2 + jj] * tmp_coords[jj];	
		}	
	}
	// test if inside
	for(short ii = 0; ii < 2; ii++) {	
		if(local[ii] < (0.0 - tol) || local[ii] > (1.0 + tol)) {
			return false;	
		}
	}
	double test_inside = 1.0 - local[0] - local[1];
			
	if(test_inside >= (0.0 - tol) && test_inside <= (1.0 + tol)) {
		return true;	
	} else {
		return false;
	}
	
}

bool Element2DTriangleBase::verify() const {
	if(region == 0) {
		Logger::get_instance()->emit(LOG_ERROR, "no region set to element");
		return false;	
	}
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		if(nodes[ii] == 0) {
			Logger::get_instance()->emit(LOG_ERROR, 1000, "node %d in element is null", ii);			
			return false;	
		}
		if(this->nodes[ii]->get_dimension() != 2) {
			ostringstream sout;
			sout << "node " << ii << " has wrong dimension " << this->nodes[ii]->get_dimension() << " != 2";
			Logger::get_instance()->emit(LOG_ERROR, sout.str());
			return false;
		}	
		for(unsigned int jj = ii + 1; jj < get_num_nodes(); jj++) {
			if(nodes[ii] == nodes[jj]) {
				Logger::get_instance()->emit(LOG_ERROR, 1000, "nodes %d and %d in element are equal", ii, jj);
				return false;	
			}
		}			
	}	
	if(element_volume <= 0.0) {
		Logger::get_instance()->emit(LOG_ERROR, 1000, "volume of element is %5.5g and therefore incorrect", element_volume);
		return false;	
	}
	if(jacobi_det <= 0.0) {
		Logger::get_instance()->emit(LOG_ERROR, 1000, "jacobi deterimnante is %5.5g and therefore incorrect", jacobi_det);
		return false;	
	} 
	return true;
}
							
void Element2DTriangleBase::get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const {
	TDKP_BOUNDS_ASSERT(edge_idx < 3, "");
	vertex_indices.resize(2);
	vertex_indices[0] = get_corner_node(edge_idx).get_index_vertex();
	vertex_indices[1] = get_corner_node((edge_idx + 1) % 3).get_index_vertex();		
}

void Element2DTriangleBase::get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const {
	get_edge(idx, vertex_indices);	
}

void Element2DTriangle::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	local_coords.resize(2);
	switch(lid) {
		case 0:
			local_coords[0] = 0.0;
			local_coords[1] = 0.0; 
			break;			
		case 1:
			local_coords[0] = 1.0;
			local_coords[1] = 0.0; 
			break;	
		case 2:
			local_coords[0] = 0.0;
			local_coords[1] = 1.0; 
			break;
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
}



/** prepares the element for assembly
 * 
 * 
 */
void Element2DTriangleBase::prepare() {
	
	// ----------------------------------------------
	// check if nodes are set
	// ----------------------------------------------
	for(unsigned int ii = 0; ii < nodes.size(); ii++) {
		TDKP_ASSERT(nodes[ii] != 0, "node " << ii << " is not set"); 	
	}	
	
	TDKP_ASSERT(this->region != 0, "");
	
	this->element_mid_point[0] = 1.0 / 3.0; // local coordinates
	this->element_mid_point[1] = 1.0 / 3.0; // local coordinates
	
	// ---------------------------------------------------
	// create jacobi matrix	
	// ---------------------------------------------------

	const double* P0 = nodes[0]->get_coords();
	const double* PII;
	// (jacobi is: Pi - P0, i = 1..3)
	for(unsigned int jj = 0; jj < 2; jj++) {		
		PII = get_corner_node(jj + 1).get_coords();
		for(unsigned int ii = 0; ii < 2; ii++) {
			jm[ii * 2 + jj] =  PII[ii] - P0[ii];
		}								
	}	
	// ----------------------------------------------
	// calculate jacobian and element volume
	// ----------------------------------------------
	this->jacobi_det = jm[0]*jm[3] - jm[1]*jm[2];
	this->jacobi_det = (jacobi_det < 0.0 ? - jacobi_det : jacobi_det);
	this->element_volume = 0.5 * this->jacobi_det;
		           			
	// ----------------------------------------------
	// calculate inverse jacobian (ii * 2 + jj)
	// ----------------------------------------------
	double bc_m_ad = jm[1]*jm[2] - jm[0]*jm[3];
	this->inverted_jacobi_matrix[0] = - jm[3] / bc_m_ad;
	this->inverted_jacobi_matrix[1] =   jm[1] / bc_m_ad;
	this->inverted_jacobi_matrix[2] =   jm[2] / bc_m_ad;
	this->inverted_jacobi_matrix[3] = - jm[0] / bc_m_ad;		
		
	this->ready = true;
		
}			

/** evaluate global form function
 * 
 * user has to ensure that coordinates are inside element 
 */
double Element2DTriangleBase::evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const {
	double global[2];
	global[0] = global_x;
	global[1] = global_y;
	double local[] = {0.0, 0,0};
	global2local(global,local);
	return evaluate_form_function(elem_shape_func, local);
} 

																			
} // end of namespace

