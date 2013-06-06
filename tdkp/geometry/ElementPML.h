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

#ifndef ELEMENTPML_H_
#define ELEMENTPML_H_

#include "tdkp/geometry/Element.h"

namespace tdkp {

/** forward declaration of pml stretch class (defined in SchroedingerPML.h) */
class PMLStretch;

/** pml element wrapper
 * 
 * this wrapper takes a standard element object and calculates the
 * element integrals via numerical integration. 
 *  
 * the problem is: element objects are defined to return real valued
 * values as a result of their integrals. that's deep down in the code.
 * templating is not a solution as this would require all analytical 
 * element integrals to be defined as templates.
 * 
 * but what i can do is the following:
 * i call calculate_element_matrices with an element returning only the
 * real part and then call calculate_element_matrices again with an element
 * returning the imaginary part of the element integral.
 * at the end, i combine them together and i'm done ...
 * 
 * now, the numerical integrals have something in common:
 * - we need to know the local integration points
 * - we need a constant jacobi determinant 
 * - we need to evaluate the global posititions of the local integration
 *   points in order to numerically integrate the global PML function
 * 
 * so the approach is to calculate for every element type the local integration
 * points and evaluate at every integration point every formfunction and every
 * derivative of the form function
 * 
 * then we can use that information for all elements of the same type. the global 
 * coordinates of the local points can easily be determined from 
 *   Xglobal(xlocal) = sum Pi Ni(xlocal) 
 */ 
class ElementPML : public Element {
public:
	ElementPML(const Element* base_element);
	ElementPML(const ElementPML& base_pml, const Element* new_element);	
	virtual ~ElementPML();
	
	// ----------------------------------------------
	// numerical integration controls
	// ----------------------------------------------
	bool returning_real_part() const;
	void return_imaginary_part();
	void return_real_part();
	void prepare(); 
	void set_new_base_element(const Element* base_element);
	void set_pml_stretch(unsigned int dir, const PMLStretch* stretch);
			
	// ----------------------------------------------
	// setup functions
	// ----------------------------------------------
	virtual void set_corner_node(unsigned short corner_id, Node* node);
	virtual const Node&  get_corner_node(unsigned short corner_idx) const;
	virtual ElementShape get_shape() const;	
	virtual unsigned int get_num_corners() const;
	virtual unsigned int get_num_edges() const;
	virtual unsigned int get_num_faces() const;
	virtual void get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const;
	virtual void get_face(unsigned int face_idx, vector<unsigned int>& vertex_indices) const;   
				
	virtual void get_node_local_coords(unsigned short lid, vector<double>& local_coords) const;
		
	virtual void get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const;
	
	virtual const Material& get_material() const;
			
	virtual bool enabled() const;
	virtual bool verify() const;

	virtual unsigned int get_element_order() const;
	// ----------------------------------------------
	// evaluation functions
	// ---------------------------------------------
	virtual double  get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const;
	virtual double  get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const;
	virtual double  get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const;
	virtual double  get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const;
	virtual double  get_single_integral_1st_order(short diffop, short elem_shape_func_1) const;
	virtual double  get_single_integral_0th_order(short elem_shape_func_1) const;	
	virtual double  evaluate_form_function(short elem_shape_func, const double* local_reference_element_coords) const;
	virtual double  evaluate_form_function_derivative(short diffop, short elem_shape_func, const double* local_reference_element_coords) const;		
			
	virtual bool   inside_element(const double& global_x, const double& global_y, const double& global_z) const;
	virtual double evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const;			
	
	// ----------------------------------------------
	// higher order elements -> additional nodes
	// ---------------------------------------------- 
	virtual unsigned int get_num_additional_nodes() const;
	virtual void get_additional_node_locator(unsigned int additional_node_idx, AdditionalNodeLocation& location_type, vector<unsigned int>& involved_vertices, vector<double>& coords, unsigned int& tag) const;
	virtual void set_additional_node(unsigned int additional_node_idx, Node* node);
	
private:
	void  initialize_numerical_integration_points();
	void  set_properties_from_base_element();
	const Element* base_element;
	bool  imaginary_integration;
	
	class Pt {
	public:
		Pt() { clear(); }
		explicit Pt(const double* cc) { set(cc); }
		void set(const double* cc) { X[0] = cc[0]; X[1] = cc[1]; X[2] = cc[2]; }		
		void clear() { X[0] = X[1] = X[2] = 0.0e0; }
		void add(const Pt& rhs) { X[0] += rhs.X[0]; X[1] += rhs.X[1]; X[2] += rhs.X[2]; }
		void multiply_and_add(const double& value, const Pt& rhs) { X[0] += value * rhs.X[0]; X[1] += value * rhs.X[1]; X[2] += value * rhs.X[2]; }
		const double& get(unsigned int ii) { return X[ii]; }		
	private:		 
		double X[3];
	};
	
	// ---------------------------------------
	// persistent data
	// ---------------------------------------
	DomainMaster integration_points;
	vector<vector<double> > shape_function_values;
	vector<vector<double> > shape_function_derivatives[3];
	
	// ---------------------------------------
	// per element data
	// ---------------------------------------
	vector<Pt>                global_points;
	vector<complex<double> >  pml_values[3]; // sx(X,Y,Z),sy(X,Y,Z),sz(X,Y,Z) 
	vector<complex<double> >  pml_determinant; // sx(X,Y,Z) * sy(X,Y,Z) * sz(X,Y,Z)
	
	const PMLStretch* stretches[3];
};

}

#endif /*ELEMENTPML_H_*/
