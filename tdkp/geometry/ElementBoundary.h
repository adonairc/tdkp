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

#ifndef ELEMENTBOUNDARY_H_
#define ELEMENTBOUNDARY_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Vector3D.h"


namespace tdkp {

class Element;
class Node;

class ElementBoundary {
	
public:
	enum BoundaryType { vertex, edge, face };
	
	ElementBoundary(BoundaryType type, unsigned int index_global_, unsigned int nvertex);
	virtual ~ElementBoundary();
	
	BoundaryType get_type() const { return boundary_type; }
	void set_corner_node(unsigned int idx_corner, const Node* node);
	void add_element(const Element* element);
	unsigned int get_num_elements() const { return elements[1] != 0 ? 2:(elements[0] != 0 ? 1:0); }
	
	unsigned int get_index_global() const { return index_global; }
	unsigned int get_num_vertices() const { return vertices.size(); }
	 	
	const Element& get_element(unsigned int ii) const;
	const Node& get_vertex(unsigned int ii) const;
	Node& get_vertex(unsigned int ii);
	const Vector3D& get_normal() const;
	
	char  get_location() const { return location; }
	void  set_location(char loc) { location = loc; }
	
	void prepare(); // prepare does not update the integrals
	void prepare_integrals();
	
	// ---------------------------------------------------
	// element boundary integrals
	// this is the boundary integral of all element shape function
	// which aren't zero on the boundary
	// i.e. all nodes that sit on the boundary for our
	// lagrangian elements
	// hermite elements would be much more complicated
	// ---------------------------------------------------	
	double get_boundary_integral_0th_order(short node_shape_func_1, short node_shape_func_2) const;
	unsigned int get_num_nodes() const { return nodes.size(); }
	const Node& get_node(unsigned int ii) const;
	const double& get_surface_volume() const;

		
private:
	void update_normal();	
				
	BoundaryType     boundary_type;
	unsigned int     index_global;
	char             location; // from dfise file
	/** normal of boundary */
	Vector3D         normal;
	vector<const Node*> vertices;
	/** every boundary element should obviously belong to two elements */
	const Element*   elements[2];
	/** cached boundary integral values, computed in prepare */
	vector<double> boundary_integrals_0th_order;
	/** the nodes with non-zero shapefunctions on the boundary */
	vector<const Node*> nodes;
	/** true if prepare_integrals was already executed */
	bool integrals_ready;
	/** surface volume of boundary (1 for point, l for edge with length l ...) */
	double surface_volume;
	
	// ---------------------------------------------------
	// functions used for evaluating single 0th order integrals
	// ---------------------------------------------------
	/** create temporary lower dimensional element */
	Element* create_temporary_element();	
	/** project coordinates onto boundary system, return false if coords are not on the boundary */ 
	bool project_to_boundary_coords(unsigned int dimension, const double* coords, double& x, double& y, const double& tolerance) const;
	
	
	
		
};

}

#endif /*ELEMENTBOUNDARY_H_*/
