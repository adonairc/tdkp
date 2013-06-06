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

#ifndef ELEMENT3DTETRAHEDRON2ND_H_
#define ELEMENT3DTETRAHEDRON2ND_H_

#include "tdkp/common/all.h"
#include "tdkp/geometry/Element3DTR.h"

namespace tdkp {

/** second order tetrahedron object
 * 
 * the node number is as follows:
 * the vertices are the same as for the linear element: N0 - N3
 * N4: 0.5 0 0
 * N5: 0.5 0.5 0
 * N6: 0 0.5 0
 * N7: 0 0 0.5
 * N8: 0.5 0 0.5
 * N9: 0 0.5 0.5
 */ 
class Element3DTetrahedron2nd : public Element3DTetrahedronBase {
public:
	Element3DTetrahedron2nd(unsigned int index_);
	virtual ~Element3DTetrahedron2nd();

	// ----------------------------------------------
	// setup functions
	// ----------------------------------------------		
	virtual unsigned int get_num_additional_nodes() const { return 6; }
	virtual void get_additional_node_locator(unsigned int additional_node_idx, AdditionalNodeLocation& location_type, vector<unsigned int>& involved_vertices, vector<double>& coords, unsigned int& tag) const;
	virtual void set_additional_node(unsigned int additional_node_idx, Node* node);
	virtual void get_node_local_coords(unsigned short lid, vector<double>& local_coords) const;
	virtual unsigned int get_element_order() const { return 2; }
	virtual unsigned int get_element_unique_type_key() const { return 9; }
	// ----------------------------------------------
	// evaluation functions
	// ----------------------------------------------
	double  get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const;
	double  get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const;
	double  get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const;

	double  get_single_integral_1st_order(short diffop, short elem_shape_func_1) const;
	double  get_single_integral_0th_order(short elem_shape_func_1) const;
	
	double  get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const;
		
	/** function to evaluate the form function inside the reference element */	
	double evaluate_form_function(short elem_shape_func, const double* local_reference_element_coords) const;
	/** function to evaluate the form function derivative inside the reference element */
	double evaluate_form_function_derivative(short diffop, short elem_shape_func, const double* local_reference_element_coords) const;					
								
	virtual double evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const;
	
private:
	double evaluate_form_function_derivative(short diffop, short elem_shape_func, const double& lx, const double& ly, const double& lz) const;

	
	
};

}

#endif /*ELEMENT3DTETRAHEDRON2ND_H_*/
