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

#ifndef ELEMENT1DHERMITECUBIC_H_
#define ELEMENT1DHERMITECUBIC_H_

#include "tdkp/geometry/Element1DLine.h"

namespace tdkp {

/** cubic hermite finite element
 * 
 * form functions
 * N1 =  2.0*x*x*x - 3.0*x*x + 1.00
 * N2 =  1.0*x*x*x - 2.0*x*x + 1.00*x 
 * N3 = -2.0*x*x*x + 3.0*x*x 
 * N4 =  1.0*x*x*x - 1.0*x*x  
 * 
 * so, N1 and N3 are vertex functions (1 at vertex, and N2/N4 are 
 * derivatives at the vertices. ATTENTION! DERIVATIVES MUST NOT BE
 * ZERO!
 * 
 */ 
class Element1DHermiteCubic : public Element1DLineBase {
public:
	Element1DHermiteCubic(unsigned int index_global_);
	virtual ~Element1DHermiteCubic();
		
	// ----------------------------------------------
	// setup functions		
	// ----------------------------------------------		
	virtual void get_node_local_coords(unsigned short lid, vector<double>& local_coords) const;
	virtual unsigned int get_element_unique_type_key() const { return 2; }
	// ----------------------------------------------
	// creation of additional nodes
	// ----------------------------------------------
	virtual void set_corner_node(unsigned short corner_id, Node* node);
	virtual const Node&  get_corner_node(unsigned short corner_idx) const;
	virtual unsigned int get_num_additional_nodes() const { return 2; }
	virtual void get_additional_node_locator(unsigned int additional_node_idx, AdditionalNodeLocation& location_type, vector<unsigned int>& involved_vertices, vector<double>& coords, unsigned int& tag) const;
	virtual void set_additional_node(unsigned int additional_node_idx, Node* node);

	virtual unsigned int get_element_order() const { return 3; }	

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
	
private:
	
};

}

#endif /*ELEMENT1DHERMITECUBIC_H_*/
