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

#include "tdkp/geometry/Element3DTR.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/common/Logger.h"

using namespace tdkp;

// ----------------------------------------------
// implementation of base tetrahedron functions
// ----------------------------------------------
Element3DTetrahedronBase::Element3DTetrahedronBase(unsigned int nnodes)
: Element(3,nnodes,4),
  jm(jacobi_matrix)
{
	this->ready		   	  = false;	
	for(int ii = 0; ii < 9; ii++) {
		this->inverted_jacobi_matrix[ii] = 0;
		this->jacobi_matrix[ii] = 0;			
	}
}


void Element3DTetrahedronBase::prepare() {
		
	// ----------------------------------------------
	// check if nodes are set
	// ----------------------------------------------
	for(unsigned int ii = 0; ii < nodes.size(); ii++) {
		TDKP_ASSERT(nodes[ii] != 0, "node " << ii << " is not set"); 	
	}	
	
	TDKP_ASSERT(this->region != 0, "");
		
	this->element_mid_point[0] = 1.0 / 4.0; // local coordinates
	this->element_mid_point[1] = 1.0 / 4.0; // local coordinates
	this->element_mid_point[2] = 1.0 / 4.0; // local coordinates
				
	// ---------------------------------------------------
	// create jacobi matrix	
	// ---------------------------------------------------	
	const double* P0 = nodes[0]->get_coords();
	const double* PII;
	// (jacobi is: Pi - P0, i = 1..3)
	for(unsigned int jj = 0; jj < 3; jj++) {		
		PII = nodes[jj + 1]->get_coords();
		for(unsigned int ii = 0; ii < 3; ii++) {
			jm[ii * 3 + jj] =  PII[ii] - P0[ii];
		}								
	}	
	// ---------------------------------------------------
	// calculate deterimante
	// ---------------------------------------------------
	jacobi_det = jm[0] * jm[4] * jm[8] 
	           + jm[1] * jm[5] * jm[6]
	           + jm[2] * jm[3] * jm[7]
	           - jm[6] * jm[4] * jm[2]
	           - jm[7] * jm[5] * jm[0]
	           - jm[8] * jm[3] * jm[1];
	jacobi_det = (jacobi_det < 0.0 ? - jacobi_det : jacobi_det);	           
   
		           
	// ------------------------------------------------------------------
	// inverted jacobi matrix (this is really the inverse of a 3x3 matrix)
	// ------------------------------------------------------------------
   	double div;
    div = - jm[2] * jm[4] * jm[6]
          + jm[1] * jm[5] * jm[6]
          + jm[2] * jm[3] * jm[7]
          - jm[0] * jm[5] * jm[7]
          - jm[1] * jm[3] * jm[8]
          + jm[0] * jm[4] * jm[8];
	inverted_jacobi_matrix[0] = (- jm[5] * jm[7] + jm[4] * jm[8]) / div;
   	inverted_jacobi_matrix[1] = (  jm[2] * jm[7] - jm[1] * jm[8]) / div;
   	inverted_jacobi_matrix[2] = (- jm[2] * jm[4] + jm[1] * jm[5]) / div;
   	inverted_jacobi_matrix[3] = (  jm[5] * jm[6] - jm[3] * jm[8]) / div;
   	inverted_jacobi_matrix[4] = (- jm[2] * jm[6] + jm[0] * jm[8]) / div;
   	inverted_jacobi_matrix[5] = (  jm[2] * jm[3] - jm[0] * jm[5]) / div;
   	inverted_jacobi_matrix[6] = (- jm[4] * jm[6] + jm[3] * jm[7]) / div;
   	inverted_jacobi_matrix[7] = (  jm[1] * jm[6] - jm[0] * jm[7]) / div;
   	inverted_jacobi_matrix[8] = (- jm[1] * jm[3] + jm[0] * jm[4]) / div;
   	
   	// -----------------------------------------------------------------
   	// calculate element volume
   	// -----------------------------------------------------------------
   	element_volume = 1.0 / 6.0 * (jacobi_det > 0.0 ? jacobi_det : - jacobi_det);
	ready = true;		           
	
}


/** check whether element is properly initialized 
 * 
 * checking whether all nodes are set, no node two times set,  
 * if element is assigned to a region, element volume was calculated
 * and if jacobi det is set
 * */
bool Element3DTetrahedronBase::verify() const {
	if(region == 0) {
		Logger::get_instance()->emit(LOG_ERROR, "no region set to element");
		return false;	
	}
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		if(nodes[ii] == 0) {
			Logger::get_instance()->emit(LOG_ERROR, 1000, "node %d in element is null", ii);			
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

/** returns true of global coordinates are inside the element */
bool Element3DTetrahedronBase::inside_element(const double& global_x, const double& global_y, const double& global_z) const {
	double global[3];
	global[0] = global_x;
	global[1] = global_y;
	global[2] = global_z;
	double local[] = {0.0, 0.0, 0.0};
	return global2local(global,local);			
}

/** deterimines local coordinates
 * 
 * @param global array of x,y global coordiantes 
 * @param local array to store final local coordinates
 * @return bool if point x,y is inside element 
 */
bool Element3DTetrahedronBase::global2local(const double global[3], double local[3]) const {
	
	const double tol = 1.0e-10;
	
	const double* P0 = nodes[0]->get_coords();
	double tmp_coords[3];
	for(short ii = 0; ii < 3; ii++) {
		tmp_coords[ii] = global[ii] - P0[ii];
		local[ii] = 0.0; 
	}
	// calc coordinates
	for(short ii = 0; ii < 3; ii++) {
		for(short jj = 0; jj < 3; jj++) {
			local[ii] += inverted_jacobi_matrix[ii * 3 + jj] * tmp_coords[jj];	
		}	
	}
	/*
	cout << "GLOBAL: " << global[0] << ", " << global[1] << ", " << global[2] << " " 
	     << "LOCAL: " << local[0] << ", " << local[1] << ", " << local[2] << "\n";
	*/
	// test if inside
	for(short ii = 0; ii < 3; ii++) {	
		if(local[ii] < (0.0 - tol) || local[ii] > (1.0 + tol)) {
			return false;	
		}
	}
			
	double test_inside = 1.0 - local[0] - local[1] - local[2];		
	return test_inside >= (0.0 - tol) && test_inside <= (1.0 + tol);
	
}

/** return vertex indices for the requested edge */
void Element3DTetrahedronBase::get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const {
	TDKP_BOUNDS_ASSERT(edge_idx < 6, "");
	vertex_indices.resize(2);
	if(edge_idx < 3) {
		vertex_indices[1] = get_corner_node((edge_idx + 1) % 3).get_index_vertex();
		vertex_indices[0] = get_corner_node(edge_idx).get_index_vertex();		
	} else {
		vertex_indices[1] = get_corner_node(3).get_index_vertex();
		vertex_indices[0] = get_corner_node(edge_idx - 3).get_index_vertex();
	}	
}

void Element3DTetrahedronBase::get_face(unsigned int face_idx, vector<unsigned int>& vertex_indices) const {

 	vertex_indices.resize(3);
	switch(face_idx) {
		case 0:
			vertex_indices[0] = get_corner_node(0).get_index_vertex();
			vertex_indices[1] = get_corner_node(1).get_index_vertex();
			vertex_indices[2] = get_corner_node(2).get_index_vertex();
			break;
		case 1:
			vertex_indices[0] = get_corner_node(0).get_index_vertex();
			vertex_indices[1] = get_corner_node(1).get_index_vertex();
			vertex_indices[2] = get_corner_node(3).get_index_vertex();
			break;		
		case 2:
			vertex_indices[0] = get_corner_node(1).get_index_vertex();
			vertex_indices[1] = get_corner_node(2).get_index_vertex();
			vertex_indices[2] = get_corner_node(3).get_index_vertex();
			break;		
		case 3:
			vertex_indices[0] = get_corner_node(2).get_index_vertex();
			vertex_indices[1] = get_corner_node(0).get_index_vertex();
			vertex_indices[2] = get_corner_node(3).get_index_vertex();
			break;		
		default:
			TDKP_GENERAL_EXCEPTION("invalid face index");		
	}
	
}

void Element3DTetrahedronBase::get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const {
	get_face(idx, vertex_indices);	
}

// ---------------------------------------------------
// implementation of linear 3D element
// ---------------------------------------------------

/** element integrals \f$ \int_{\hat{K}} \hat{N}_{i}(\mathbf{\zeta}) \hat{N}_{j}(\mathbf{\zeta}) \node \bigg(\frac{\partial \mathbf{x}}{\partial{\mathbf{\zeta}}}\bigg) \node d\mathbf{\zeta} \f$ */
const double Element3DTR::reference_elem_integral_zero[4][4] = {
   {  1.0/60.0, 1.0/120.0, 1.0/120.0, 1.0/120.0 },
   { 1.0/120.0,  1.0/60.0, 1.0/120.0, 1.0/120.0 },
   { 1.0/120.0, 1.0/120.0,  1.0/60.0, 1.0/120.0 },
   { 1.0/120.0, 1.0/120.0, 1.0/120.0,  1.0/60.0 }
};
/** gradients of reference element shape functions */
const double Element3DTR::reference_elem_gradients[4][3] = {
	{ -1.0, -1.0, -1.0 }, 
	{  1.0,  0.0,  0.0 },
	{  0.0,  1.0,  0.0 }, 
	{  0.0,  0.0,  1.0 }
};

/** element constructor
 * 
 * @param index_ global element index (used to determine values from fields)
 */
Element3DTR::Element3DTR(unsigned int index_)
: Element3DTetrahedronBase(4) 
{		
	this->index_global    = index_;	
}

Element3DTR::~Element3DTR() {
	this->region = 0;	
}




/** caluclate first order element integrals (analytical solution)
 * 
 * differential operator is defined as macro-constant
 * d_dx = 0, d_dy = 1, d_dz = 2
 * 
 * the element integral evaulated is:
 * \f$ \int_{K} \frac{d}{dx_{a}} N_{i}(\mathbf{x}) \frac{d}{dx_{b}} N_{j}(\mathbf{x}) d\mathbf{x} \f$
 * therefore, the element is mapped to the reference element and then evaluated as
 * \f$ \int_{\hat{K}} \bigg(\frac{\partial \mathbf{\zeta}}{\partial x_{a}}\bigg)^{T} \nabla_{\zeta} \hat{N}_{i}(\mathbf{\zeta}) \bigg(\frac{\partial \mathbf{\zeta}}{\partial x_{b}}\bigg)^{T} \nabla_{\zeta} \hat{N}_{j}(\mathbf{\zeta}) \mathrm{det}\bigg(\frac{\partial \mathbf{x}}{\partial \mathbf{\zeta}}\bigg) d\mathbf{\zeta} \f$

 * @param diffop_1 differential operator for shape function one
 * @param elem_shape_func_1 the first shape function
 * @param diffop_2 differential operator for shape function two
 * @param elem_shape_func_2 the second shape function
 * */
double Element3DTR::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {
	
	TDKP_ASSERT(this->element_volume > 0, "element_volume > 0");
	TDKP_ASSERT(this->ready, "element is not ready");
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && (unsigned int)elem_shape_func_1 < get_num_nodes(), "elem_shape_func1 >= 0 && (unsigned int)elem_shape_func1 < get_num_nodes()");
	TDKP_BOUNDS_ASSERT(elem_shape_func_2 >= 0 && (unsigned int)elem_shape_func_2 < get_num_nodes(), "elem_shape_func2 >= 0 && (unsigned int)elem_shape_func2 < get_num_nodes()");
	TDKP_BOUNDS_ASSERT(diffop_1 >= D_DX && diffop_1 <= D_DZ, "diffop_1 >= D_DX && diffop_1 <= D_DZ");	
	TDKP_BOUNDS_ASSERT(diffop_2 >= D_DX && diffop_2 <= D_DZ, "diffop_2 >= D_DX && diffop_2 <= D_DZ");	
	
	double left, right;
	left = right = 0;
	for(int ii = 0; ii < 3; ii++) {
		left  += inverted_jacobi_matrix[diffop_1 + ii * 3] * reference_elem_gradients[elem_shape_func_1][ii];
		right += inverted_jacobi_matrix[diffop_2 + ii * 3] * reference_elem_gradients[elem_shape_func_2][ii];
	}
	return (left * right * element_volume);
	
}


/** return 1st order element integrals
 * 
 * in terms of the reference element, the integral may be written as
 * \f$ \int_{\hat{K}} \bigg(\frac{\partial \mathbf{zeta}}{\partial x_{a}}\bigg)^{T} \cdot \nabla_{\zeta} \hat{N}_{i} \hat{N}_{j}  \mathrm{det}\bigg(\frac{\partial \mathbf{x}}{\partial \mathbf{\zeta}}\bigg) d\mathbf{\zeta} \f$.
 * as the gradient of \f$\hat{N}_{i}\f$ and the jacobi matrix is constant for this element, we may evaluate everything analytically as 
 * \f$ \mathrm{det}\bigg(\frac{\partial \mathbf{x}}{\partial \mathbf{\zeta}}\bigg) \bigg(\frac{\partial \mathbf{zeta}}{\partial x_{a}}\bigg)^{T} \cdot \nabla_{\zeta} N_{i} \int_{\hat{K}} N_{j} d\mathbf{\zeta} \f$.
 * The integral gives \f$\frac{1}{24}\f$. 
 * 
 * 
 * returns \f$ \int_{K} \frac{\partial}{\partial x} N_{i} N_{j} d\mathbf{x}\f$
 */
double Element3DTR::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {

	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && (unsigned int)elem_shape_func_1 < get_num_nodes(), "elem_shape_func1 >= 0 && (unsigned int)elem_shape_func1 < get_num_nodes()");
	TDKP_BOUNDS_ASSERT(elem_shape_func_2 >= 0 && (unsigned int)elem_shape_func_2 < get_num_nodes(), "elem_shape_func2 >= 0 && (unsigned int)elem_shape_func2 < get_num_nodes()");

	double left = 0.0;
	for(int ii = 0; ii < 3; ii++) {
		left += inverted_jacobi_matrix[diffop + ii * 3] * reference_elem_gradients[elem_shape_func_1][ii];
	}
	return jacobi_det * left * (1.0/24.0); 		
		
}

/** get 0th order element integrals
 * 
 * returns \f$ \int_{K} N_{i} N_{j} d\mathbf{x} \f$
 * @param elem_shape_func_1 element shape function index i (= 0..3)
 * @param elem_shape_func_2 element shape function index j (= 0..3)
 * @return the integral
 */
double Element3DTR::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {

	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && (unsigned int)elem_shape_func_1 < get_num_nodes(), "elem_shape_func1 >= 0 && (unsigned int)elem_shape_func1 < get_num_nodes()");
	TDKP_BOUNDS_ASSERT(elem_shape_func_2 >= 0 && (unsigned int)elem_shape_func_2 < get_num_nodes(), "elem_shape_func2 >= 0 && (unsigned int)elem_shape_func2 < get_num_nodes()");
	
	return reference_elem_integral_zero[elem_shape_func_1][elem_shape_func_2] * jacobi_det;
		
}


double Element3DTR::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	
	unsigned int obj_index = elem_shape_func_2 
	                       + 4 * elem_shape_func_1
	                       + 16 * nodal_data_point;	                      
	switch(obj_index) {
        case 21: // Ncc2Nii2Njj2
        case 42: // Ncc3Nii3Njj3
        case 63: // Ncc4Nii4Njj4
        case 0: // Ncc1Nii1Njj1
            return ((1.000000e+00) / (1.200000e+02)) * jacobi_det;
        case 2: // Ncc1Nii1Njj3
        case 3: // Ncc1Nii1Njj4
        case 4: // Ncc1Nii2Njj1
        case 5: // Ncc1Nii2Njj2
        case 8: // Ncc1Nii3Njj1
        case 10: // Ncc1Nii3Njj3
        case 12: // Ncc1Nii4Njj1
        case 15: // Ncc1Nii4Njj4
        case 16: // Ncc2Nii1Njj1
        case 17: // Ncc2Nii1Njj2
        case 20: // Ncc2Nii2Njj1
        case 22: // Ncc2Nii2Njj3
        case 23: // Ncc2Nii2Njj4
        case 25: // Ncc2Nii3Njj2
        case 26: // Ncc2Nii3Njj3
        case 29: // Ncc2Nii4Njj2
        case 31: // Ncc2Nii4Njj4
        case 32: // Ncc3Nii1Njj1
        case 34: // Ncc3Nii1Njj3
        case 37: // Ncc3Nii2Njj2
        case 38: // Ncc3Nii2Njj3
        case 40: // Ncc3Nii3Njj1
        case 41: // Ncc3Nii3Njj2
        case 43: // Ncc3Nii3Njj4
        case 46: // Ncc3Nii4Njj3
        case 47: // Ncc3Nii4Njj4
        case 48: // Ncc4Nii1Njj1
        case 51: // Ncc4Nii1Njj4
        case 53: // Ncc4Nii2Njj2
        case 55: // Ncc4Nii2Njj4
        case 58: // Ncc4Nii3Njj3
        case 59: // Ncc4Nii3Njj4
        case 60: // Ncc4Nii4Njj1
        case 61: // Ncc4Nii4Njj2
        case 62: // Ncc4Nii4Njj3
        case 1: // Ncc1Nii1Njj2
            return ((1.000000e+00) / (3.600000e+02)) * jacobi_det;
        case 7: // Ncc1Nii2Njj4
        case 9: // Ncc1Nii3Njj2
        case 11: // Ncc1Nii3Njj4
        case 13: // Ncc1Nii4Njj2
        case 14: // Ncc1Nii4Njj3
        case 18: // Ncc2Nii1Njj3
        case 19: // Ncc2Nii1Njj4
        case 24: // Ncc2Nii3Njj1
        case 27: // Ncc2Nii3Njj4
        case 28: // Ncc2Nii4Njj1
        case 30: // Ncc2Nii4Njj3
        case 33: // Ncc3Nii1Njj2
        case 35: // Ncc3Nii1Njj4
        case 36: // Ncc3Nii2Njj1
        case 39: // Ncc3Nii2Njj4
        case 44: // Ncc3Nii4Njj1
        case 45: // Ncc3Nii4Njj2
        case 49: // Ncc4Nii1Njj2
        case 50: // Ncc4Nii1Njj3
        case 52: // Ncc4Nii2Njj1
        case 54: // Ncc4Nii2Njj3
        case 56: // Ncc4Nii3Njj1
        case 57: // Ncc4Nii3Njj2
        case 6: // Ncc1Nii2Njj3
            return ((1.000000e+00) / (7.200000e+02)) * jacobi_det;
		default:
			TDKP_GENERAL_EXCEPTION("invalid index");		
	}	                       
}

/** return 1st order single integral
 * 
 * mainly used for assembly of load vector. returns \f$ \int_{K} \partial_{\alpha} N_{i} d\mathbf{x} \f$
 */
double Element3DTR::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
	double left = 0.0;
	for(int ii = 0; ii < 3; ii++) {
		left += inverted_jacobi_matrix[diffop + ii * 3] * reference_elem_gradients[elem_shape_func_1][ii];
	}		
	return left * jacobi_det * 1.0 / 6.0;
}
double Element3DTR::get_single_integral_0th_order(short elem_shape_func_1) const {
	return 1.0 / 24.0 * jacobi_det; 	
}


	
/** deterimines interpolation coordinates for a point with pcoords
 * 
 * if point is in element, contrib will be filled with contributions of values at element node
 * to value at point. if point is not inside element, false is returned and nothing set to contrib
 * 
 * @param pcoords array of x,y,z coordiantes of point p
 * @param contrib array to store final contributions of nodes to given point
 * @return bool if point is inside element and contrib was set
 */
bool Element3DTR::get_contribution(const double pcoords[3], double contrib[4]) const {
	
	const double tol = 1.0e-12;
	
	const double* P0 = nodes[0]->get_coords();
	double tmp_coords[3];
	for(short ii = 0; ii < 3; ii++) {
		tmp_coords[ii] = pcoords[ii] - P0[ii];
		contrib[ii + 1] = 0.0; 
	}	
	for(short ii = 0; ii < 3; ii++) {
		for(short jj = 0; jj < 3; jj++) {
			contrib[ii + 1] += inverted_jacobi_matrix[ii * 3 + jj] * tmp_coords[jj];	
		}	
		if(contrib[ii + 1] < (0.0 - tol) || contrib[ii + 1] > (1.0 + tol)) {
			return false;	
		}
	}
	contrib[0] = 1.0 - contrib[1] - contrib[2] - contrib[3];		
	return contrib[0] >= (0.0 - tol) && contrib[0] <= (1.0 + tol);
	
}


/** function to evaluate the form function inside the reference element */	
double Element3DTR::evaluate_form_function(
	short elem_shape_func, 
	const double* local
) const {
	
	switch(elem_shape_func) {
		case 0: // N_0 = 1 - x - y
			return 1.0 - local[0] - local[1] - local[2];
			break;			
		case 1: // N_1 = x
			return local[0];
			break;			
		case 2: // N_2 = y
			return local[1];
			break;
		case 3: // N_3 = z
			return local[2];
			break;	
		default:
			TDKP_GENERAL_EXCEPTION("this element has only four elem shape functions!");				
	}	
}

/** function to evaluate the form function derivative inside the reference element */
double Element3DTR::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	
	TDKP_BOUNDS_ASSERT(diffop < 3 && diffop >= 0, ""); 		
	double ret = 0.0;
	for(int ii = 0; ii < 3; ii++) {
		ret += inverted_jacobi_matrix[diffop + ii * 3] * reference_elem_gradients[elem_shape_func][ii];
	}
	return ret;
}					




double Element3DTR::evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const {
	double global[3];
	global[0] = global_x;
	global[1] = global_y;
	global[2] = global_z;
	double local[] = {0.0, 0.0, 0.0};
	// it's not my problem if the coordinates are outside of the element ...
	global2local(global,local);
	return evaluate_form_function(elem_shape_func, local);
}


void Element3DTR::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	
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
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
	
}



