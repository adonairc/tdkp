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

#include "tdkp/geometry/ElementBoundary.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/common/Configuration.h"

namespace tdkp {

ElementBoundary::ElementBoundary(BoundaryType type, unsigned int index_global_, unsigned int nvertices)
: boundary_type(type),
  index_global(index_global_),
  vertices(nvertices, 0),
  integrals_ready(false),
  surface_volume(0.0)
{
	elements[0] = 0;
	elements[1] = 0;
}
ElementBoundary::~ElementBoundary() {}

void ElementBoundary::set_corner_node(unsigned int idx_corner, const Node* node) {
	TDKP_BOUNDS_ASSERT(vertices.size() > idx_corner, "");
	vertices[idx_corner] = node;	
}
/** add elements (maximum of two ...) */
void ElementBoundary::add_element(const Element* element) {
	if(elements[0] == 0) {
		elements[0] = element;	
	} else if(elements[1] == 0) {
		elements[1] = element;
	} else {
		TDKP_GENERAL_EXCEPTION("tried to add more than 2 elements to element boundary!");
	}
}

/** calculate boundary normal pointing out of element 0 into element 1 
 *
 *  IF ELEMENT IS OUTER ELEMENT, NORMAL IS NOT CALCULATED! 
 */
void ElementBoundary::update_normal() {
	
	if(get_num_elements() != 2) {
		return; // don't calculate normal
	}
	
	// -----------------------------------------------
	// test if all vertices are set
	// -----------------------------------------------
	for(unsigned int ii = 0; ii < vertices.size(); ii++) {
		TDKP_ASSERT(vertices[ii] != 0, "vertices[ii] != 0 failed for ii = " << ii);	
	}
	
	// -----------------------------------------------
	// calculate oriented normal
	// -----------------------------------------------
	if(vertices.size() == 1 && elements[0]->get_dimension() == 1) {
		// -----------------------------------------
		// normal is clear, sign not, so check sign
		// -----------------------------------------		
		int e1nidx = -1;		
		for(unsigned int ii = 0; ii < elements[1]->get_num_nodes(); ii++) {
			if(elements[1]->get_node(ii).get_index_vertex() != -1 && 
			   elements[1]->get_node(ii).get_index_vertex() != vertices[0]->get_index_vertex()) {
			   	e1nidx = ii;
			   	break;
			}			
		}
		TDKP_ASSERT(e1nidx != -1, "");
		double P1 = elements[1]->get_node(e1nidx).get_coord(0);
		if(P1 > vertices[0]->get_coord(0)) {
			// boundary: left is element 0, right is element 1
			normal = Vector3D(1.0, 0.0, 0.0);
		} else {
			// boundary: left is element 1, right is element 0
			normal = Vector3D(-1.0, 0.0, 0.0);
		}
	} else if((vertices.size() == 2 && elements[0]->get_dimension() == 2) ||
	          (vertices.size() >= 3 && elements[0]->get_dimension() == 3)) {	          	
		if(elements[0]->get_dimension() == 2) {		          	
			// ------------------------------------------------
			// calculate edge normal 
			// ------------------------------------------------
			Vector3D edge(
				vertices[1]->get_coord(0) - vertices[0]->get_coord(0), 
				vertices[1]->get_coord(1) - vertices[0]->get_coord(1),
				0.0
			);
			// normal = edge ^ ez / |edge ^ ez|
			normal = Vector3D::cross_product(edge, Vector3D(0.0, 0.0, 1.0));
		} else {
			// ------------------------------------------------
			// calculate face normal 
			// ------------------------------------------------
			// normal = p1 ^ p2 / |p1 ^ p2|
			normal = Vector3D::cross_product(
				Vector3D(
					vertices[2]->get_coord(0) - vertices[0]->get_coord(0),
					vertices[2]->get_coord(1) - vertices[0]->get_coord(1),
					vertices[2]->get_coord(2) - vertices[0]->get_coord(2)
				),
				Vector3D(
					vertices[1]->get_coord(0) - vertices[0]->get_coord(0),
					vertices[1]->get_coord(1) - vertices[0]->get_coord(1),
					vertices[1]->get_coord(2) - vertices[0]->get_coord(2)
				)			
			);				
		}
		TDKP_ASSERT(normal.norm() > 0.0, "");
		normal.normalize();
		// ------------------------------------------------
		// check orientation
		// ------------------------------------------------
		// first, find node of element 1 thats not involved in the edge 
		int e1nidx = -1;
		for(unsigned int ii = 0; ii < elements[1]->get_num_nodes(); ii++) {
			if(elements[1]->get_node(ii).get_index_vertex() != -1) {
				bool candidate = true;
				for(unsigned int jj = 0; jj < vertices.size(); jj++) {
					if(elements[1]->get_node(ii).get_index_vertex() == vertices[jj]->get_index_vertex()) {
						candidate = false;
						break;													
					}	
				}
				// if node ii is not involved in the edge, stop here
				if(candidate) {
					e1nidx = ii;	
					break;
				}	
			}	
		}
		TDKP_ASSERT(e1nidx != -1, "");
		Vector3D vec_into_e1(
			elements[1]->get_node(e1nidx).get_coord(0) - vertices[0]->get_coord(0), 
			elements[1]->get_node(e1nidx).get_coord(1) - vertices[0]->get_coord(1),
			elements[1]->get_node(e1nidx).get_coord(2) - vertices[0]->get_coord(2)					
		);
		// if normal does not point into same direction, change sign of normal!
		if(Vector3D::dot_product(vec_into_e1, normal) < 0.0) {
			normal = -1.0 * normal;
		}		
	} else {
		TDKP_GENERAL_EXCEPTION("invalid number of vertices or dimension");	
	}

	
}

/** returns normal of element boundary 
 * 
 * the normal is oriented such that it points out of 
 * element 0 into element 1
 */
const Vector3D& ElementBoundary::get_normal() const {
	// normal points from 0 > 1. if only 0 is set, it points to nowhere ... 
	TDKP_ASSERT(elements[0] != 0 && elements[1] != 0, "normal is only available for interior boundaries");
	TDKP_BOUNDS_ASSERT(tdkp_math::abs(normal.norm() - 1.0) < 1.0e-4, "normal not calculated yet. you need to call update_normal()");
	return normal;
	
}
	

/** normal of boundary */

/** returns element of boundary 
 * 
 * @param ii element index < 2
 */
const Element& ElementBoundary::get_element(unsigned int ii) const {
	TDKP_BOUNDS_ASSERT(ii < 2, "");
	return *elements[ii];	
}
const Node& ElementBoundary::get_vertex(unsigned int ii) const {
	TDKP_BOUNDS_ASSERT(vertices.size() > ii,"");
	return *vertices[ii];	
}
Node& ElementBoundary::get_vertex(unsigned int ii) {
	return const_cast<Node&>(static_cast<const ElementBoundary&>(*this).get_vertex(ii));	
}

const Node& ElementBoundary::get_node(unsigned int ii) const {
	TDKP_BOUNDS_ASSERT(ii < nodes.size(), "");
	return *nodes[ii];	
}

/** calculate element single integrals, update normal */
void ElementBoundary::prepare() {	
	update_normal();	
}

/** returns procomputed boundary integral */
double ElementBoundary::get_boundary_integral_0th_order(
	short node_shape_func_1, short node_shape_func_2
) const {
	
	unsigned int obj_index = node_shape_func_1 * get_num_nodes()
	                       + node_shape_func_2;	
	TDKP_BOUNDS_ASSERT(obj_index < boundary_integrals_0th_order.size(), "");
	return boundary_integrals_0th_order[obj_index];	
}

void ElementBoundary::prepare_integrals() {
	
	if(get_num_elements() != 2 || integrals_ready) {
		return;	
	}
	
	TDKP_ASSERT(boundary_integrals_0th_order.size() == 0, "");
	TDKP_ASSERT(nodes.size() == 0, "");
 	integrals_ready = true;
 	const double projection_tolerance = Configuration::get_instance()->get("elementboundary_node_projection_tolerance");
 	
	// ----------------------------------------
	// 1D case: boundary is a vertex and our
	// basis is nodal, so the integral is simply
	// the form function at the vertex which is 
	// one. 
	// ----------------------------------------
	if(get_element(0).get_dimension() == 1) {
		TDKP_BOUNDS_ASSERT(vertices.size() == 1, "");
		boundary_integrals_0th_order.push_back(1.0);
		nodes.push_back(vertices[0]);
		// surface volume is just one ... because the charge is the charge ...
		surface_volume = 1.0;
		return;
	}
 	
	// ----------------------------------------
	// create temporary element
	// ----------------------------------------
	Element* tmp_elem = create_temporary_element();
	
	// ----------------------------------------
	// create clones of corner nodes and set to temporary element
	// ----------------------------------------
	vector<Node*> tmp_nodes;
	double x,y;
	for(unsigned int ii = 0; ii < get_num_vertices(); ii++) {
		// project vertex onto new boundary coordinates
		project_to_boundary_coords(
			get_element(0).get_dimension(), 
			get_vertex(ii).get_coords(),
			x, y, projection_tolerance
		);
		// create new node
		if(get_element(0).get_dimension() == 3) {
			tmp_nodes.push_back(
				new Node(static_cast<unsigned int>(get_vertex(ii).get_index_global()), x, y)
			);				
		} else {
			tmp_nodes.push_back(
				new Node(static_cast<unsigned int>(get_vertex(ii).get_index_global()), x)
			);			
		}
		tmp_nodes.back()->set_index_vertex(get_vertex(ii).get_index_vertex());		
		// set node
		tmp_elem->set_corner_node(ii, tmp_nodes.back());	
	} 
	
	// ----------------------------------------
	// - get additional node coordinates from low-dim element
	// - from node coordinates, determine appropriate
	//   node in element 0
	// - create clone of node and set to temp element 
	// ----------------------------------------
	vector<double>         tmp_node_coords;
	unsigned int           unused_tag;
	vector<unsigned int>   unused_involved;	
	Element::AdditionalNodeLocation unused_loc;	
	// ----------------------------------------
	// for every additional node in tmp element
	// ----------------------------------------
	for(unsigned int ii = 0; ii < tmp_elem->get_num_additional_nodes(); ii++) {
		
		// --------------------------------------
		// get node coordinates
		// --------------------------------------
		tmp_elem->get_additional_node_locator(ii, unused_loc, unused_involved, tmp_node_coords, unused_tag);
		
		// --------------------------------------
		// loop over all additional nodes in real 
		// element 0 and find matching node
		// --------------------------------------
		bool matched = false;
		bool good_projection = false;
		double cached_dist[Element::max_num_nodes]; 
		for(unsigned int jj = 0; jj < get_element(0).get_num_nodes(); jj++) {
			
			good_projection = project_to_boundary_coords(
				get_element(0).get_dimension(), 
				get_element(0).get_node(jj).get_coords(),
				x, y, projection_tolerance
			);		
						
			
			double dist = 0.0;
			if(get_element(0).get_dimension() == 2) {
				dist = tdkp_math::abs(x - tmp_node_coords[0]);
			} else {
				dist = sqrt((x - tmp_node_coords[0]) * (x - tmp_node_coords[0]) + 
				            (y - tmp_node_coords[1]) * (y - tmp_node_coords[1]));
			}
			cached_dist[jj] = dist;
		//	cout << "   projected coords: " << x << ", " << y << ", dist to add node: " << dist << " ";
			if(dist < projection_tolerance && good_projection) {
			//	cout << " TAKING IT\n";
				matched = true;
				// --------------------------------------
				// create new node
				// --------------------------------------
				switch(get_element(0).get_dimension()) {
					case 2:
						tmp_nodes.push_back(
							new Node(static_cast<unsigned int>(get_element(0).get_node(jj).get_index_global()), x)
						);
						break;			
					case 3:
						tmp_nodes.push_back(
							new Node(static_cast<unsigned int>(get_element(0).get_node(jj).get_index_global()), x, y)
							);
						break;
					default:
						TDKP_GENERAL_EXCEPTION("invalid element dimension"); 	
				}					
				// done here
				break;
			} else {
				//cout << " TOO FAR\n";	
			}	
		}
		if(!matched) {
			ostringstream sout;
			sout << "could not match additional node " << ii << " for element boundary " 
			      << get_index_global() << ". distance to element 0's nodes is:\n";			      
			for(unsigned int jj = 0; jj < get_element(0).get_num_nodes(); jj++) {				
				sout <<  jj << ": " << cached_dist[jj] << "\n";	
			}
			sout << "current projection tolerance is: " << projection_tolerance << ". maybe you should increase it,\n"
			     << "by increasing elementboundary_node_projection_tolerance.";
			TDKP_GENERAL_EXCEPTION(sout.str()); 
		}
		tmp_elem->set_additional_node(ii, tmp_nodes.back());		
	}
		
	// ----------------------------------------
	// build nodes array 
	// ----------------------------------------
	// for all nodes in temporary element
	for(unsigned int ii = 0; ii < tmp_elem->get_num_nodes(); ii++) {
		bool matched = false;
		// find corresponding node in element 0
		for(unsigned int jj = 0; jj < get_element(0).get_num_nodes(); jj++) {
			// on match, insert node into nodes array
			if(tmp_elem->get_node(ii).get_index_global() == get_element(0).get_node(jj).get_index_global()) {
				matched = true;
				nodes.push_back(&(get_element(0).get_node(jj)));	
			}	
		}	
		TDKP_ASSERT(matched, "could not match low dim node with element node!");
	}

	// ----------------------------------------
	// prepare element
	// ----------------------------------------
	tmp_elem->set_region(&(get_element(0).get_region()), 0);
	tmp_elem->prepare();
	surface_volume = tmp_elem->get_volume();
	
	// ----------------------------------------
	// evaluate element integral and store results
	// ----------------------------------------
	boundary_integrals_0th_order.assign(tmp_elem->get_num_nodes() * tmp_elem->get_num_nodes(), 0.0);
	for(unsigned int ii = 0; ii < tmp_elem->get_num_nodes(); ii++) {
		for(unsigned int jj = 0; jj < tmp_elem->get_num_nodes(); jj++) {
			boundary_integrals_0th_order[ii * tmp_elem->get_num_nodes() + jj] = tmp_elem->get_element_integral_0th_order(ii,jj);
		}	
	}
	// ----------------------------------------
	// clean up
	// ----------------------------------------
	delete tmp_elem;
	for(unsigned int ii = 0; ii < tmp_nodes.size(); ii++) {
		delete tmp_nodes[ii];	
	}
		
}

const double& ElementBoundary::get_surface_volume() const {
	TDKP_ASSERT(surface_volume > 0.0, "surface volume has not yet been calculated. call prepare_boundaries on geometry object to make it available.");
	return surface_volume;	
}

Element* ElementBoundary::create_temporary_element() {
	Element::ElementShape shape;
	if(boundary_type == edge) {
		shape = Element::line;
	} else if(boundary_type == face) {
		shape = Element::triangle;
	} else {
		TDKP_GENERAL_EXCEPTION("not implemented boundary type");	
	}
	return Element::factory(shape, 0, get_element(0).get_dimension() - 1, get_element(0).get_element_order());
}

/** project coordinates onto boundary system, return false if coords are not on the boundary */ 
bool ElementBoundary::project_to_boundary_coords(
	unsigned int dimension, const double* coords, double& x, double& y, const double& tolerance
) const {
	
	// ----------------------------------------------
	// we know two projections: transform a line down to
	// the x-axis (2D edges) and express face coordinates
	// as a plane
	// ----------------------------------------------
	if(dimension == 2) {
		TDKP_ASSERT(get_num_vertices() == 2, "");
		// V0 + (V1 - V0) / |V1 - V0| x = Pi
		Vector3D V1V0(
			get_vertex(1).get_coord(0) - get_vertex(0).get_coord(0),
			get_vertex(1).get_coord(1) - get_vertex(0).get_coord(1),
			0.0
		);
		TDKP_BOUNDS_ASSERT(V1V0.norm() > 0.0, "");
		Vector3D PiV0(
			coords[0] - get_vertex(0).get_coord(0),
			coords[1] - get_vertex(0).get_coord(1),
			0.0			
		);
		// test if point matches V0
		if(PiV0.norm() < tolerance) {
			x = 0.0; // V0 is zero point
			return true;	
		}
		// -----------------------------------------------
		// check if we may match the point
		// in other words: same vectors, different length
		// -----------------------------------------------
		Vector3D normalized_PiV0(PiV0); normalized_PiV0.normalize();
		Vector3D normalized_V1V0(V1V0); normalized_V1V0.normalize();
		if(tdkp_math::abs(Vector3D::dot_product(normalized_PiV0, normalized_V1V0) - 1.0) < tolerance) {
			x = PiV0.norm(); // length of vector
			return true;
		}
		return false;		
	} else if(dimension == 3) {
		TDKP_ASSERT(get_num_vertices() == 3, "3D -> 2D projection only implemented for triangular boundaries");
		// -----------------------------------------
		// triangle points
		// -----------------------------------------
		Vector3D P0(get_vertex(0).get_coords());
		Vector3D P1(get_vertex(1).get_coords());
		Vector3D P2(get_vertex(2).get_coords());
		Vector3D Pi(coords);
	/*	TDKP_TRACE("P0 " << P0);
		TDKP_TRACE("P1 " << P1);
		TDKP_TRACE("P2 " << P2);
		TDKP_TRACE("Pi " << Pi); */
		// -----------------------------------------
		// build spanning vectors
		// -----------------------------------------
		Vector3D avec = P1 - P0;
		Vector3D bvec = P2 - P0;
		TDKP_ASSERT(avec.norm() > 0.0, "");
		TDKP_ASSERT(bvec.norm() > 0.0, "");
		// -----------------------------------------
		// base vector n1 is avec normalized
		// -----------------------------------------
		Vector3D n1 = avec; n1.normalize();
		// -----------------------------------------
		// base vector n2 is (a x b) x a (normalized)
		// then n2 is orthogonal to n1 but in plane of a and b
		// -----------------------------------------
		Vector3D n2 = Vector3D::cross_product(Vector3D::cross_product(avec,bvec),avec);
		n2.normalize();
//		TDKP_TRACE("n1 " << n1);
//		TDKP_TRACE("n2 " << n2);
		// -----------------------------------------
		// desired position Pi
		// is solved via (n1,n2) (x,y)^T = Pi - P0
		// -----------------------------------------
		Pi = Pi - P0; // position relative to P0
		double nmat[4];
		unsigned int perm[][3] = {{0,1,2},{0,2,1},{1,2,0}};
		unsigned int pp = 0;
		unsigned int i0, i1, i2;		
		double det_nmat;
		do {
			i0 = perm[pp][0];
			i1 = perm[pp][1];
			i2 = perm[pp][2];
			nmat[0] = n1(i0); nmat[1] = n2(i0);
			nmat[2] = n1(i1); nmat[3] = n2(i1);
			det_nmat = tdkp_math::abs(tdkp_math::det_2x2_matrix(nmat));
			pp++;			
		} while(tdkp_math::abs(tdkp_math::det_2x2_matrix(nmat)) < 1.0e-8 && pp < 3); 
		TDKP_ASSERT(det_nmat > 1.0e-8, "n1 and n2 are parallel!!!"); 
		tdkp_math::invert_2x2_matrix(nmat);
		x = nmat[0] * Pi(i0) + nmat[1] * Pi(i1);
		y = nmat[2] * Pi(i0) + nmat[3] * Pi(i1);
		// test if remaining coordinate of Pi is matched 
		double z = n1(i2) * x + n2(i2) * y;
//		TDKP_TRACE("last coord check? " << z);
		if(tdkp_math::abs(z - Pi(i2)) < tolerance) {
			return true;	
		} else {
			return false;	
		}
	} else {
		TDKP_GENERAL_EXCEPTION("invalid dimension!");	
	}
}


} // end of namespace
