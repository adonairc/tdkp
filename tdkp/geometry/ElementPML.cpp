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

#include "tdkp/geometry/ElementPML.h"
#include "tdkp/probdefs/SchroedingerPML.h"

namespace tdkp {
	

ElementPML::ElementPML(const Element* base_element_)
: Element(
	base_element_->get_dimension(), 
	base_element_->get_num_nodes(),
	base_element_->get_num_boundaries()
  ),
  base_element(base_element_),
  imaginary_integration(false)
{
	stretches[0] = stretches[1] = stretches[2] = 0; 
	TDKP_ASSERT(base_element_->ready, "base_element_->ready failed. the base element must damn be prepared!");
	this->set_properties_from_base_element();
	this->initialize_numerical_integration_points(); // one time job	
}

ElementPML::ElementPML(const ElementPML& base_pml, const Element* new_element) 
: Element(
	base_pml.get_dimension(), 
	base_pml.get_num_nodes(),
	base_pml.get_num_boundaries()
  ),
  base_element(new_element),
  imaginary_integration(false),
  integration_points(base_pml.integration_points),
  shape_function_values(base_pml.shape_function_values),
  shape_function_derivatives(base_pml.shape_function_derivatives)
{
	// copy data
	stretches[0] = base_pml.stretches[0];  
	stretches[1] = base_pml.stretches[1];
	stretches[2] = base_pml.stretches[2];
	this->set_properties_from_base_element();
	this->ready = false;		
}
		
ElementPML::~ElementPML() {

}

void ElementPML::set_pml_stretch(unsigned int dir, const PMLStretch* stretch) {
	TDKP_ASSERT(dir < 3, "");
	stretches[dir] = stretch;	
}

void ElementPML::return_imaginary_part() {
	imaginary_integration = true;
}
void ElementPML::return_real_part() {
	imaginary_integration = false;
}

/** copy base element properties to current object */
void ElementPML::set_properties_from_base_element() {

	TDKP_BOUNDS_ASSERT(nodes.size() == base_element->nodes.size(), "");
	TDKP_BOUNDS_ASSERT(element_boundaries.size() == base_element->element_boundaries.size(), "");
	// ----------------------------------------------
	// copy data
	// ----------------------------------------------
	nodes.assign(
		base_element->nodes.begin(), 
		base_element->nodes.end()
	);
	element_boundaries.assign(
		base_element->element_boundaries.begin(), 
		base_element->element_boundaries.end()
	);
	region = base_element->region;
	index_global = base_element->index_global;
	index_region = base_element->index_region;
	jacobi_det   = base_element->jacobi_det;
	element_volume = base_element->element_volume;
	ready          = base_element->ready;	
	element_mid_point[0] = base_element->element_mid_point[0];
	element_mid_point[1] = base_element->element_mid_point[1];
	element_mid_point[2] = base_element->element_mid_point[2]; 
	TDKP_ASSERT(jacobi_det > 0, "");

}

void ElementPML::set_new_base_element(const Element* new_base_element) {
	TDKP_ASSERT(new_base_element->ready, "new_base_element->ready failed. the base element must damn be prepared!");
	TDKP_ASSERT(base_element->get_dimension() == new_base_element->get_dimension(), "");
	TDKP_ASSERT(base_element->get_shape() == new_base_element->get_shape(), "");
	TDKP_ASSERT(base_element->get_element_order() == new_base_element->get_element_order(), "");
	TDKP_ASSERT(base_element->get_num_nodes() == new_base_element->get_num_nodes(), "");	
	base_element = new_base_element;
	this->set_properties_from_base_element();
	this->ready = false;
	this->prepare();
}

void ElementPML::initialize_numerical_integration_points() {
	
	// ---------------------------------------
	// get integration points
	// ---------------------------------------
	integration_points = base_element->get_numerical_integration_points(20);
	
	// ---------------------------------------
	// init storage space
	// ---------------------------------------
	shape_function_values.resize(base_element->get_num_nodes());	
	for(unsigned int ii = 0; ii < get_dimension(); ii++) {
		shape_function_derivatives[ii].resize(base_element->get_num_nodes());
		pml_values[ii].resize(integration_points.get_number_of_points());			
	}
	for(unsigned int nn = 0; nn < get_num_nodes(); nn++) {
		shape_function_values[nn].resize(integration_points.get_number_of_points());
		for(unsigned int ii = 0; ii < get_dimension(); ii++) {
			shape_function_derivatives[ii][nn].resize(integration_points.get_number_of_points());
		}
	}	
	global_points.assign(integration_points.get_number_of_points(), Pt());
	pml_determinant.assign(integration_points.get_number_of_points(), 1.0);
	
	// --------------------------------------
	// evaluate form functions at integration points
	// --------------------------------------
	/*
	for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
		const DomainPoint& ipt = integration_points.get_point(ii);
		// evaluate form functions
		for(unsigned int nn = 0; nn < get_num_nodes(); nn++) {
			shape_function_values[nn][ii] = base_element->evaluate_form_function(nn, ipt.get_coords());
			// evaluate derivatives
			for(unsigned int dd = 0; dd < get_dimension(); dd++) {
				shape_function_derivatives[dd][nn][ii] = base_element->evaluate_form_function_derivative(dd, nn, ipt.get_coords());		
			}			
		}				
	}*/				
}

void ElementPML::prepare() {
	if(!ready) {	
		// --------------------------------------
		// evaluate form functions at integration points
		// --------------------------------------
		for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
			const DomainPoint& ipt = integration_points.get_point(ii);
			// evaluate form functions
			for(unsigned int nn = 0; nn < get_num_nodes(); nn++) {
				shape_function_values[nn][ii] = base_element->evaluate_form_function(nn, ipt.get_coords());
				// evaluate derivatives
				for(unsigned int dd = 0; dd < get_dimension(); dd++) {
					shape_function_derivatives[dd][nn][ii] = base_element->evaluate_form_function_derivative(dd, nn, ipt.get_coords());		
				}			
			}				
		}						
		// ------------------------------------------------
		// determine global positions of integration points
		// ------------------------------------------------
		// reset to zero
		global_points.assign(integration_points.get_number_of_points(), Pt());
		// calculate global positions
		for(unsigned int nn = 0; nn < get_num_nodes(); nn++) {
			//TDKP_TRACE("moep");
			// transform node to simple pt object 
			Pt nodal_point(get_node(nn).get_coords());					
			for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
				//TDKP_TRACE("whoop " << shape_function_values[nn][ii]);
				global_points[ii].multiply_and_add(shape_function_values[nn][ii], nodal_point);
			}
		}
		// ------------------------------------------------			
		// evaluate PML functions at integration points
		// ------------------------------------------------
		for(unsigned int dd = 0; dd < get_dimension(); dd++) {
			pml_values[dd].assign(integration_points.get_number_of_points(), 1.0);
			if(stretches[dd] != 0) {
				for(unsigned int ii = 0; ii < global_points.size(); ii++) {
				//	TDKP_TRACE("revaluate wurst " << global_points[ii].get(dd));					
					pml_values[dd][ii] = stretches[dd]->evaluate(global_points[ii].get(dd));	
				}	
			}
		}
		// ------------------------------------------------
		// evaluate determinante
		// ------------------------------------------------
		pml_determinant.assign(integration_points.get_number_of_points(), 1.0);
		for(unsigned int dd = 0; dd < get_dimension(); dd++) {
			for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {			
				pml_determinant[ii] *= pml_values[dd][ii];	
			}
		}

		
	}
	this->ready = true;
}


// ----------------------------------------------
// evaluation functions
// ---------------------------------------------
double ElementPML::get_element_integral_2nd_order(
	short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2
) const {
	cplx ret = 0.0;
	const vector<double>& shape_1 = shape_function_derivatives[diffop_1][elem_shape_func_1];
	const vector<double>& shape_2 = shape_function_derivatives[diffop_2][elem_shape_func_2];
	for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
		ret += shape_1[ii] * shape_2[ii] * jacobi_det 
		     * integration_points.get_point(ii).get_weight()
		     * 1.0 / pml_values[diffop_1][ii]
		     * 1.0 / pml_values[diffop_2][ii]
		     * pml_determinant[ii];		     
	}
	/*double cmp = base_element->get_element_integral_2nd_order(diffop_1, elem_shape_func_1, diffop_2, elem_shape_func_2);
	if(tdkp_math::abs(cmp - ret.real()) > 1.0e-2) {
		cout << "base element volume: " << base_element->get_volume() << ", jacobi det: " 
		     << jacobi_det << ", pml vals = " << 1.0 / pml_values[diffop_1][0] 
		     << ", " <<  1.0 / pml_values[diffop_2][0] 
		     << ", shape1: " << shape_1[0] << ", shape2 " << shape_2[0] << "\n";  
		cout << "2nd order: d" << diffop_1 << "N" << elem_shape_func_1
		     << "d" << diffop_2 << "N" << elem_shape_func_2
		     << " numerical = " << ret << ", analytical  = " << cmp << ", diff = "
		     << tdkp_math::abs(cmp - ret.real()) << "\n";
	}*/		
	if(imaginary_integration) {
		//return 0.0;
		return ret.imag();	
	} else {
		//return cmp;
		return ret.real();	
	}
}
double ElementPML::get_element_integral_1st_order(
	short diffop, short elem_shape_func_1, short elem_shape_func_2
) const {
	cplx ret = 0.0;
	const vector<double>& shape_1 = shape_function_derivatives[diffop][elem_shape_func_1];
	const vector<double>& shape_2 = shape_function_values[elem_shape_func_2];
	for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
		ret += shape_1[ii] * shape_2[ii] * jacobi_det 
		     * integration_points.get_point(ii).get_weight()
		     * 1.0 / pml_values[diffop][ii]
		     * pml_determinant[ii];		     
	}
	/*
	double cmp = base_element->get_element_integral_1st_order(diffop, elem_shape_func_1, elem_shape_func_2);
	
	if(tdkp_math::abs(cmp - ret.real()) > 1.0e-4) {
		cout << "1th order: d" << diffop << "N" << elem_shape_func_1
		     << "N" << elem_shape_func_2
		     << " numerical = " << ret << ", analytical  = " << cmp << ", diff = "
		     << tdkp_math::abs(cmp - ret.real()) << "\n";
	}*/		
		
	if(imaginary_integration) {
		return ret.imag();	
	} else {
		return ret.real();	
	}	
}
double ElementPML::get_element_integral_0th_order(
	short elem_shape_func_1, short elem_shape_func_2
) const {
	cplx ret = 0.0;
	const vector<double>& shape_1 = shape_function_values[elem_shape_func_1];
	const vector<double>& shape_2 = shape_function_values[elem_shape_func_2];
	for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
		ret += shape_1[ii] * shape_2[ii] * jacobi_det 
		     * integration_points.get_point(ii).get_weight()
		     * pml_determinant[ii];		     
	}
	/*
	double cmp = base_element->get_element_integral_0th_order(elem_shape_func_1, elem_shape_func_2);
	if(tdkp_math::abs(cmp - ret.real()) > 1.0e-2) {
		cout << "0th order: N" << elem_shape_func_1
		     << "N" << elem_shape_func_2
		     << " numerical = " << ret << ", analytical  = " << cmp << ", diff = "
		     << tdkp_math::abs(cmp - ret.real()) << "\n";
	}	*/	
	if(imaginary_integration) {
		return ret.imag();	
	} else {
		return ret.real();	
	}		
}
double ElementPML::get_element_integral_0th_order_nodal_data(
	short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2
) const {
	cplx ret = 0.0;
	const vector<double>& nodal   = shape_function_values[nodal_data_point];
	const vector<double>& shape_1 = shape_function_values[elem_shape_func_1];
	const vector<double>& shape_2 = shape_function_values[elem_shape_func_2];
	for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
		ret += nodal[ii] * shape_1[ii] * shape_2[ii] * jacobi_det 
		     * integration_points.get_point(ii).get_weight()
		     * pml_determinant[ii];		     
	}
	if(imaginary_integration) {
		return ret.imag();	
	} else {
		return ret.real();	
	}	
}

double ElementPML::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
	cplx ret = 0.0;
	const vector<double>& shape_1 = shape_function_derivatives[diffop][elem_shape_func_1];	
	for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
		ret += shape_1[ii] * jacobi_det
		     * 1.0 / pml_values[diffop][ii] 
		     * integration_points.get_point(ii).get_weight()
		     * pml_determinant[ii];		     
	}
	if(imaginary_integration) {
		return ret.imag();	
	} else {
		return ret.real();	
	}		
}
double ElementPML::get_single_integral_0th_order(short elem_shape_func_1) const {
	cplx ret = 0.0;
	const vector<double>& shape_1 = shape_function_values[elem_shape_func_1];
	for(unsigned int ii = 0; ii < integration_points.get_number_of_points(); ii++) {
		ret += shape_1[ii] * jacobi_det 
		     * integration_points.get_point(ii).get_weight()
		     * pml_determinant[ii];		     
	}
	if(imaginary_integration) {
		return ret.imag();	
	} else {
		return ret.real();	
	}	
}

// ---------------------------------------------
// collection of functions that have to be defined
// for a standard element but are more required for
// building and setting up structures
// ---------------------------------------------
void ElementPML::get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!"); 	
}
void ElementPML::get_face(unsigned int face_idx, vector<unsigned int>& vertex_indices) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");
}   		
bool ElementPML::inside_element(const double& global_x, const double& global_y, const double& global_z) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");
}   
double ElementPML::evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");
}   			
unsigned int ElementPML::get_num_additional_nodes() const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");
}   
void ElementPML::get_additional_node_locator(unsigned int additional_node_idx, AdditionalNodeLocation& location_type, vector<unsigned int>& involved_vertices, vector<double>& coords, unsigned int& tag) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");
}   
void ElementPML::set_additional_node(unsigned int additional_node_idx, Node* node) {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");
}   
void ElementPML::get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");
}   
void ElementPML::set_corner_node(unsigned short corner_id, Node* node) {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!"); 	
}
const Node&  ElementPML::get_corner_node(unsigned short corner_idx) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!"); 	
}
double ElementPML::evaluate_form_function(
	short elem_shape_func, const double* local_reference_element_coords
) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");		
}
double ElementPML::evaluate_form_function_derivative(
	short diffop, short elem_shape_func, const double* local_reference_element_coords
) const {
	TDKP_GENERAL_EXCEPTION("ElementPML: this function is disabled and calling it does not make sense on the masking PML element!");	
}

// ---------------------------------------------
// collection of functions simply returning
// base element properties
// ---------------------------------------------
Element::ElementShape ElementPML::get_shape() const {
	return base_element->get_shape();
}	
unsigned int ElementPML::get_num_corners() const {
	return base_element->get_num_corners();	
}
unsigned int ElementPML::get_num_edges() const {
	return base_element->get_num_edges();	
}
unsigned int ElementPML::get_num_faces() const {
	return base_element->get_num_faces();	
}			
void ElementPML::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	base_element->get_node_local_coords(lid, local_coords);	
}
const Material& ElementPML::get_material() const {
	return base_element->get_material();	
}		
bool ElementPML::enabled() const {
	return base_element->enabled();	
}
bool ElementPML::verify() const {
	return base_element->verify();	
}
unsigned int ElementPML::get_element_order() const {
	return base_element->get_element_order();	
}

}
