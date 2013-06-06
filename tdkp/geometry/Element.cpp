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


#include "tdkp/geometry/Element.h"
#include "tdkp/geometry/Element1DLine.h"
#include "tdkp/geometry/Element1DLine2nd.h"
#include "tdkp/geometry/Element1DHermiteCubic.h"
#include "tdkp/geometry/Element2DTriangle.h"
#include "tdkp/geometry/Element2DTriangle2nd.h"
#include "tdkp/geometry/Element2DTriangleHermiteCubic.h"
#include "tdkp/geometry/Element2DRect.h"
#include "tdkp/geometry/Element2DRect2nd.h"
#include "tdkp/geometry/Element3DTR.h"
#include "tdkp/geometry/Element3DTetrahedron2nd.h"
#include "tdkp/geometry/ElementContact.h"

namespace tdkp {


Element::Element(unsigned int dimension_, unsigned int nnodes, unsigned int nboundaries)
: dimension(dimension_),
  nodes(nnodes, 0),
  element_boundaries(nboundaries, 0),
  region(0),
  index_global(0),
  index_region(0), 
  jacobi_det(0),
  element_volume(0),  
  ready(false)
{
	for(short ii = 0; ii < 3; ii++) {
		element_mid_point[ii] = 0;
	}
}

Element::~Element() {
	
}

/** set vertex to element
 * 
 * in a usual gridfile, only corner vertices are defined while
 * internal nodes due to higher order elements have to be generated
 * internally. this function therefore takes the corner nodes.
 */	
void Element::set_corner_node(unsigned short lid, Node* node) {
	TDKP_ASSERT(lid < this->get_num_nodes(), "lid " << lid << " < get_num_nodes() " << get_num_nodes() << " failed for elem  " << index_global);
	this->nodes[lid] = node;
}


/** get vertex from element
 */
const Node& Element::get_corner_node(unsigned short corner_idx) const {
	TDKP_ASSERT(corner_idx < this->get_num_nodes(), "");
	TDKP_ASSERT(this->nodes[corner_idx]->get_index_vertex() != -1, "");
	return *this->nodes[corner_idx];			
}

/** set element boundary 
 * 
 * element boundaries are nodes for 1D, edges for 2D and faces for 3D
 */
void Element::set_element_boundary(unsigned int idx, const ElementBoundary* element_boundary) {
	TDKP_BOUNDS_ASSERT(element_boundaries.size() > idx, "");
	element_boundaries[idx] = element_boundary;
}

/** return boundary of element */
const ElementBoundary& Element::get_boundary(unsigned int idx) const {
	TDKP_BOUNDS_ASSERT(element_boundaries.size() > idx, "");
	TDKP_BOUNDS_ASSERT(element_boundaries[idx] != 0, "boundary " << idx << " does not exist!");
	return *element_boundaries[idx];	
}

/** return local index of element boundary
 * 
 * throws exception if boundary is not a part of the element 
 */
unsigned int Element::get_local_boundary_index(const ElementBoundary& boundary) const {
	for(unsigned int ii = 0; ii < get_num_boundaries(); ii++) {
		if(boundary.get_index_global() == get_boundary(ii).get_index_global()) {
			return ii;	
		}	
	}
	TDKP_GENERAL_EXCEPTION("could not find matching boundary object");		
}

/** set region and region element index
 * 
 * every element get's a pointer to its region.
 * values in dfDataset's are stored region wise. therefore we need the information, to
 * which element a dataset may belong. that's why we need the index in the region
 * 
 * @param region_ pointer to a region object
 * @param index_region_ index of element in region
 */
void Element::set_region(const Region* region_, unsigned int index_region_) {
	this->region = region_;
	this->index_region = index_region_;
}

const Region& Element::get_region() const {
	return *this->region;	
}


/** compare to elements if they are equal
 * 
 * two elements are equal when they have the same nodes
 */
bool Element::compare(const Element& other) const {
	unsigned int equal_nodes = 0;
	if(this->get_num_nodes() != other.get_num_nodes()) {
		return false;	
	}
	for(unsigned int ii = 0; ii < this->get_num_nodes(); ii++) {
		for(unsigned int jj = 0; jj < this->get_num_nodes(); jj++) {
			if(this->nodes[ii] == other.nodes[jj]) {
				equal_nodes++;	
				break;
			}
		}	
	}	
	return equal_nodes == this->get_num_nodes();
}


ostream& operator<<(ostream& out, const Element &elem) {
	out << "element " << elem.get_index_global() << " region: " << elem.region->get_name() << "\n nodes: ";	
	for(unsigned int ii = 0; ii < elem.get_num_nodes(); ii++) {
		out << " " << elem.get_node(ii);
	}
	return out;
}



/** return bounding box of element
 * 
 * @return low lowest left corner in low, highest in high
 */
void Element::get_bounding_box(Node& low, Node& high) const {
	TDKP_ASSERT(this->get_num_nodes() > 0, "element not properly initialized");
	// set both to first node
	for(unsigned short ii = 0; ii < this->dimension; ii++) {
		low.get_coord(ii)  = this->nodes[0]->get_coord(ii);
		high.get_coord(ii) = this->nodes[0]->get_coord(ii);	
	}
	// get bounding box for remaining nodes
	for(unsigned short ii = 1; ii < this->get_num_nodes(); ii++) {
		for(unsigned short jj = 0; jj < this->dimension; jj++) {
			low.get_coord(jj)  = min(low.get_coord(jj),  this->nodes[ii]->get_coord(jj));			
			high.get_coord(jj) = max(high.get_coord(jj), this->nodes[ii]->get_coord(jj));
		}			
	}
}

/** return elements material object */
const Material& Element::get_material() const {	
	// usually its the one of the region
	TDKP_ASSERT(enabled(), "must not call this function on disabled elements!"); 
	return this->get_region().get_material();		
}

/** standard element verifictation */
bool Element::verify() const {
	TDKP_ASSERT(region != 0, "no region set to element");
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		TDKP_ASSERT(nodes[ii] != 0, "node in element is null");
		for(unsigned int jj = ii + 1; jj < get_num_nodes(); jj++) {
			if(nodes[ii] == nodes[jj]) {
				Logger::get_instance()->emit(LOG_ERROR, 1000, "nodes %d and %d in element are equal", ii, jj);
				return false;	
			}
		}	
	}
	TDKP_ASSERT(this->element_volume != 0, "element volume is not calculated or is zero");
	return true;
}

void Element::print() const {
	ostringstream sout;
	sout << *this;
	Logger::get_instance()->emit(LOG_INFO, sout.str());	
}
/** element factory: create element of given dimension and polynom degree */
Element* Element::factory(ElementShape elem_shape, unsigned int element_index_global, unsigned int dimension, unsigned int polynom_degree) {

	Element* elem = 0;
	
	if((dimension == 1 &&  line == elem_shape) ||
	   (dimension == 2 && (elem_shape == triangle || elem_shape == rectangle)) ||
	   (dimension == 3 && elem_shape == tetrahedron)
	) {
		// -----------------------------------------------
		// create new element
		// -----------------------------------------------
		switch(elem_shape) {
			case line:
				TDKP_ASSERT(dimension == 1, "");
				if (polynom_degree == 3) {
					elem = new Element1DHermiteCubic(element_index_global);
				} else if(polynom_degree == 2) {
					elem = new Element1DLine2nd(element_index_global);
				} else {
					elem = new Element1DLine(element_index_global);	
				}
				break; 
			case triangle:
				TDKP_ASSERT(dimension == 2, "");
				if (polynom_degree == 3) {
					elem = new Element2DTriangleHermiteCubic(element_index_global);
				} else if(polynom_degree == 2) {
					elem = new Element2DTriangle2nd(element_index_global);
				} else {
					elem = new Element2DTriangle(element_index_global);
				}				
				break;
			case rectangle:
				TDKP_ASSERT(dimension == 2, "");
				if(polynom_degree == 2) {
					elem = new Element2DRect2nd(element_index_global);
				} else {
					elem = new Element2DRect(element_index_global);
				}
				break;
			case tetrahedron:
				TDKP_ASSERT(dimension == 3, "");						
				if(polynom_degree == 2) {
					elem = new Element3DTetrahedron2nd(element_index_global);
				} else {
					elem = new Element3DTR(element_index_global);
				}
				break;
			default:
				TDKP_GENERAL_EXCEPTION("unknown element type(" << elem_shape << ") for element (grid idx " << element_index_global << ")");			
		}		
	} else {
		unsigned int nnode = 0;
		// consistency check, thats what we expect
		if(dimension == 1) {
			TDKP_ASSERT(elem_shape == vertex, "");
			nnode = 1;		
		} else if(dimension == 2) {
			TDKP_ASSERT(elem_shape == line, "");
			nnode = 2;
		} else {
			TDKP_ASSERT(elem_shape == triangle || elem_shape == rectangle, "");
			if(elem_shape == triangle) {
				nnode = 3;	
			} else {
				nnode = 4;	
			}
		}   				
		// ----------------------------------------------------
		// assign a contact element for a low-dimensional element type
		// ----------------------------------------------------			 
		elem = new ElementContact(
			dimension - 1,
			nnode, // num nodes
			nnode, // num boundaries of the element ...
			element_index_global
		);
	}
	return elem;
}

/** check if supplied node has correct node coordinates */
bool test_additional_node(const Element& elem, const Node& node, unsigned int additional_idx) {
	Element::AdditionalNodeLocation location_type;
	vector<unsigned int> involved_vertices;
	vector<double> coords;
	unsigned int tag;
	elem.get_additional_node_locator(additional_idx, location_type, involved_vertices, coords, tag);
	double dist = 0.0;
	for(unsigned int ii = 0; ii < elem.get_dimension(); ii++) {
		double tmp = coords[ii] - node.get_coord(ii);
		dist += (tmp * tmp);	
	}		
	dist = sqrt(dist);
	if(dist < 1.0e-8) {
		return true;	
	} else {
		return false;	
	}
}

} // end of namespace
