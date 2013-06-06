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

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Domain.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/geometry/Region.h"
#include "tdkp/geometry/ElementBoundary.h"

using namespace std;
namespace tdkp {

class Element {
	friend ostream& operator<<(ostream& out, const Element &elem);
	friend class ElementPML;
public:	
	/** possible element shapes of which tdkp knowns about
	 *
	 * note that vertex is used for dummy 0D contact elements,
	 * line is used for 1D elements, triangle and rectangle obviously for 2D
	 * and the only 3D element is the tetrahedron
	 *
	 * to add more element shapes, you have to add the corresponding
	 * new element class to the Element::factory and also check that
	 * all the io stuff knows about it (i.e. translation from MED
	 * shapes into tdkps internal shapes)
	 *
	 * btw: ElementShapes must be assigned to indexes [0,num_element_shapes]
	 *
	 */
	enum ElementShape {
		vertex = 0, line = 1, triangle = 2, rectangle = 3, tetrahedron = 4
	};
	/** number of ElementShape */
	static const int num_element_shapes = 5;
	/** maximum number of element shape functions. increase it if you add elements with more shape functions */
	static const int max_num_nodes = 10;
	/** maximum number of corners (== vertices) for all elements. increase if you add other functions. */
	static const int max_num_corners = 4;

	Element(unsigned int dimension, unsigned int nnodes, unsigned int nboundaries);	
	virtual ~Element();
	
	// ----------------------------------------------
	// setup functions
	// ----------------------------------------------
	virtual void set_corner_node(unsigned short corner_id, Node* node);
	virtual const Node&  get_corner_node(unsigned short corner_idx) const;
	virtual ElementShape get_shape() const = 0;
	virtual unsigned int get_element_unique_type_key() const { TDKP_GENERAL_EXCEPTION("element not properly implemented. every element MUST have a unique identifier"); }	
	virtual unsigned int get_num_corners() const = 0;
	virtual unsigned int get_num_edges() const = 0;
	virtual unsigned int get_num_faces() const = 0;
	virtual void get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const { TDKP_GENERAL_EXCEPTION("no edges avl. wrong call or not implemented in derived class."); }
	virtual void get_face(unsigned int face_idx, vector<unsigned int>& vertex_indices) const { TDKP_GENERAL_EXCEPTION("no faces avl. wrong call or not implemented in derived class."); }   
				
	const Node&  get_node(unsigned short lid) const;	
	Node&        get_node(unsigned short lid);
	virtual void get_node_local_coords(unsigned short lid, vector<double>& local_coords) const = 0;
	
	
	virtual void get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const = 0;
	void set_element_boundary(unsigned int idx, const ElementBoundary* element_boundary);	
	const ElementBoundary& get_boundary(unsigned int idx) const;	
	unsigned int get_num_boundaries() const { return element_boundaries.size(); }
	unsigned int get_local_boundary_index(const ElementBoundary& boundary) const;	
	
	void set_region(const Region* region, unsigned int index_region);
	const Region& get_region() const;
	virtual const Material& get_material() const;
			
	virtual bool enabled() const { return true; }
	virtual bool verify() const;
	virtual bool is_ready() const { return ready; }
	virtual void prepare() = 0;
	virtual unsigned int get_element_order() const = 0; 
	/** return index global inside geometry object */
	unsigned int get_index_global() const { return this->index_global; }
	/** return index inside region as defined in grid file */
	unsigned int get_index_region() const { return this->index_region; }
	unsigned int get_num_nodes() const { return this->nodes.size(); }
	unsigned int get_dimension() const { TDKP_ASSERT(this->dimension > 0, "this->dimension > 0"); return this->dimension; }
					
	bool compare(const Element& other) const;
	  
	// ----------------------------------------------
	// evaluation functions
	// ---------------------------------------------
	/** return d_m Ni d_n Nj */
	virtual double  get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const = 0;
	/** return d_m Ni Nj */
	virtual double  get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const = 0;
	/** return Ni Nj */
	virtual double  get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const = 0;
	/** return Nn * Ni * Nj
	 * 
	 * given a field (e.g. charge) given at nodes F(x) = sum_n f_n * Nn, 
	 * 0 order integrals can be built using this functions (preferrable over element wise constant fields)
	 */	 
	virtual double  get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const = 0;
 
	virtual double  get_single_integral_1st_order(short diffop, short elem_shape_func_1) const = 0;
	virtual double  get_single_integral_0th_order(short elem_shape_func_1) const = 0;	
	/** function to evaluate the form function (damn, form function is german, shape function would be correct) inside the reference element */		
	virtual double evaluate_form_function(short elem_shape_func, const double* local_reference_element_coords) const = 0;
	/** function to evaluate the form function derivative inside the reference element BUT FOR THE GLOBAL DIFFOP! */
	virtual double evaluate_form_function_derivative(short diffop, short elem_shape_func, const double* local_reference_element_coords) const = 0;		
	
	double          get_volume()     const { return element_volume; }
	double          get_jacobi_det() const { return jacobi_det; }
	
	/** integrate solution over element */
	template<class T>
	T  integrate_solution(const T *) const;
		
	/** return differential of formfunction with given values at given position */
	template<class T>	
	T differentiate_solution(short diffop, const T* values, const double* local_reference_element_coords) const;
		
	/** return average differential of formfunction with given values */
	template<class T>					
	T  differentiate_solution(short diffop, const T* ) const;

	// --------------------------------------------------
	// functions used for numerical revaliation of 
	// element integrals (numerically determination of element integrals)
	// used in unit test in order to evaluate analytical integration
	// --------------------------------------------------	
	void get_bounding_box(Node& low, Node& high) const;	 	
	virtual bool   inside_element(const double& global_x, const double& global_y, const double& global_z) const { TDKP_GENERAL_EXCEPTION("not implemented for base element class. probably the element you are using does not implement these"); }
	virtual double evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const { TDKP_GENERAL_EXCEPTION("not implemented for base element class. probably the element you are using does not implement these"); }				
	virtual DomainMaster get_numerical_integration_points(unsigned int npoints) const { TDKP_GENERAL_EXCEPTION("not implemented for base element class. probably the element you are using does not implement these"); }
	
	// ----------------------------------------------
	// higher order elements -> additional nodes
	// ----------------------------------------------
	enum AdditionalNodeLocation { inner_node, vertex_node, edge_node }; 
	virtual unsigned int get_num_additional_nodes() const { return 0; }
	virtual void get_additional_node_locator(unsigned int additional_node_idx, AdditionalNodeLocation& location_type, vector<unsigned int>& involved_vertices, vector<double>& coords, unsigned int& tag) const { TDKP_GENERAL_EXCEPTION("thats a bad element not deriving the class appropriately"); }
	virtual void set_additional_node(unsigned int additional_node_idx, Node* node) { TDKP_GENERAL_EXCEPTION("thats a bad element not deriving the class appropriately"); }
	
	virtual void print() const;
	
	static Element* factory(ElementShape shape, unsigned int element_index_global, unsigned int dimension, unsigned int polynom_degree);
	
	
protected:
	/** dimension of the element */
	const unsigned int dimension;
	/** array (length nnode) of pointers to nodes of */
	vector<Node*>  nodes;
	/** boundary objects of element */
	vector<const ElementBoundary*> element_boundaries;
	/** pointer to region where the element belongs to */
	const Region* region;
	/** global index number (corresponds to index in geometry) */
	unsigned int   index_global;
	/** index in the region where the element belongs to */
	unsigned int   index_region;
	/** jacobi determinante */
	double         jacobi_det;
	/** volume of element */
	double         element_volume;
	/** boolean to indicate whether the element is ready for assembly */
	bool           ready;
	/** mid point of element */
	double         element_mid_point[3];

};

 	
bool test_additional_node(const Element& elem, const Node& node, unsigned int additional_idx);

ostream& operator<<(ostream& out, const Element &elem);

inline const Node& Element::get_node(unsigned short lid) const {
	TDKP_BOUNDS_ASSERT(lid < get_num_nodes(), "lid (" << lid << ") < get_num_nodes() (" << get_num_nodes() << ") failed");
	TDKP_ASSERT(this->nodes[lid] != 0, "node " << lid << " has not yet been set!");
	return *this->nodes[lid];	
}
inline Node& Element::get_node(unsigned short lid) {
	return const_cast<Node&>(
		static_cast<const Element&>(*this).get_node(lid)
	);	
}

/** integrate solution over element */
template<class T>
T Element::integrate_solution(const T* values) const {
	T ret = 0.0;
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		ret += values[ii] * get_single_integral_0th_order(ii);	
	}
	return ret;
}
	
/** return differential of formfunction with given values at given position */
template<class T>	
T Element::differentiate_solution(short diffop, const T* values, const double* local_reference_element_coords) const {
	T ret = 0.0;
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		ret += values[ii] * evaluate_form_function_derivative(diffop, ii, local_reference_element_coords);	
	}
	return ret;	
}
	
/** return average differential of formfunction with given values */
template<class T>					
T Element::differentiate_solution(short diffop, const T* values ) const {
	return differentiate_solution(diffop, values, element_mid_point);
}

} // end of namespace tdkp

#endif /*ELEMENT_H_*/
