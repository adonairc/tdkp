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

#ifndef ELEMENTCONTACT_H_
#define ELEMENTCONTACT_H_

#include "Element.h"

namespace tdkp {


/** contact element, regarded as dummy element without function */
class ElementContact : public Element {
public:
	ElementContact(unsigned int dimension_, unsigned int nnodes_, unsigned int nboundaries_, unsigned int index_global_);
	virtual ~ElementContact();
	virtual ElementShape get_shape() const;	
	virtual unsigned int get_num_corners() const;
	virtual unsigned int get_num_edges() const;
	virtual unsigned int get_num_faces() const;
	virtual void get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const;
	virtual void get_face(unsigned int face_idx, vector<unsigned int>& vertex_indices) const;
	virtual void get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const;
	virtual void get_node_local_coords(unsigned short lid, vector<double>& local_coords) const;
	bool verify() const;
	void prepare(); 
	virtual bool enabled() const { return false; }
	virtual unsigned int get_element_order() const { return 0; }
	virtual unsigned int get_element_unique_type_key() const { return 10; }
	// ----------------------------------------------
	// evaluation functions
	// ----------------------------------------------
	double  get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const;
	double  get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const;
	double  get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const;
	double  get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const;
										
	double  get_single_integral_1st_order(short diffop, short elem_shape_func_1) const;
	double  get_single_integral_0th_order(short elem_shape_func_1) const;
		
	double evaluate_form_function(short elem_shape_func, const double* local_reference_element_coords) const;
	double evaluate_form_function_derivative(short diffop, short elem_shape_func, const double* local_reference_element_coords) const;
														
	
};

}

#endif /*ELEMENTCONTACT_H_*/
