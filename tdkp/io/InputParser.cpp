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

#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <stdio.h>

#include "tdkp/io/InputParser.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/common/Logger.h"
#include "tdkp/geometry/ElementContact.h"
#include "tdkp/io/BaseGridReader.h"

using namespace tdkp;
using namespace std;

const char* InputParser::NodeDataHeader  = "NodeDataBegin";
const char* InputParser::NodeDataFooter  = "NodeDataEnd";
const char* InputParser::ElementDataHeader = "ElementDataBegin";
const char* InputParser::ElementDataFooter = "ElementDataEnd";



namespace tdkp {


/** add vertex object to tree (if it does not exist) and add the element index */
void InputParser::build_vertex_objects_tree(
	unsigned int element_index, vector<VertexObject*>& objects, VertexObject& vobj
) {
	// check if it already exists
	if(objects[vobj.get_lowest_vertex_index()] != 0) {
		// check in tree
		VertexObject* located_obj = objects[vobj.get_lowest_vertex_index()]->locate(vobj);
		// create new if it does not exist yet
		if(located_obj == 0) {					
			located_obj = vobj.clone();						
			objects[vobj.get_lowest_vertex_index()]->insert(located_obj);																					
		}
		located_obj->add_element(element_index);
	} else {
		// add new tree
		objects[vobj.get_lowest_vertex_index()] = vobj.clone();
		objects[vobj.get_lowest_vertex_index()]->add_element(element_index);	
	}				
}

/** create global indices to objects */ 
unsigned int InputParser::set_global_indices_to_tree(vector<VertexObject*>& objects) {
	unsigned int index_global = 0;
	for(unsigned int ii = 0; ii < objects.size(); ii++) {
		VertexObject* vobj = objects[ii];
		while(vobj != 0) {
			vobj->set_index_global(index_global++);
			vobj = vobj->get_next();	
		}	
	}	
	return index_global;
}

// -------------------------------------------------------------
// new implementation of grid reading for general grid obj
// -------------------------------------------------------------
Geometry* InputParser::read_geometry(const BaseGridReader& grid_reader, unsigned int polynom_degree) {
	
	TimeMeasurements::get_instance().start("building of FE grid");
	
	// ------------------------------------------------
	// dump some information
	// ------------------------------------------------
	ostringstream sout;
	sout.setf(ios::left);
	sout  << "InputParser:\n"
	      << setw(22) << "  dimension: "         << grid_reader.get_dimension()    << "\n"
	      << setw(22) << "  vertices: "          << grid_reader.get_num_vertices() << "\n"
	      << setw(22) << "  elements: "          << grid_reader.get_num_elements() << "\n"
	      << setw(22) << "  regions: "           << grid_reader.get_num_regions()  << "\n"
	      << setw(22) << "  polynomial degree: " << polynom_degree << "\n";
	TDKP_LOGMSG(LOG_INFO_DEVEL1, sout.str());
	      
	// ------------------------------------------------
	// create basic geometry object
	// ------------------------------------------------
	Geometry* geometry = new Geometry(
		grid_reader.get_dimension(),
		grid_reader.get_num_elements(),
		grid_reader.get_num_vertices(),
		grid_reader.get_num_regions()
	);
	geometry->set_identifier(grid_reader.get_unique_identifier());
	const unsigned int dimension = grid_reader.get_dimension();
	bool faces_required        = false;
	bool edges_required        = false;
	bool bnd_vertices_required = false;
	
	// ------------------------------------------------
	// bnd vertices / faces / edges are required as boundary obj
	// and later use to create additional nodes
	// ------------------------------------------------
	if(dimension == 1) {
		bnd_vertices_required = true;	
	}
	if(dimension == 2) {
		bnd_vertices_required = true;
		edges_required = true;	
	}
	if(dimension == 3) {		
		bnd_vertices_required = true;
		faces_required = true;
		edges_required = true;	
	}	
		
	// ------------------------------------------------	
	// read vertices and create their node objects
	// ------------------------------------------------
	//TimeMeasurements::get_instance().start("vertex nodes");
	vector<Node*> vertices(grid_reader.get_num_vertices(), 0);
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: creating vertex nodes");
	#pragma omp parallel for
	for(int ii = 0; ii < (signed)grid_reader.get_num_vertices(); ii++) {
		Node* node = 0;
		switch(dimension) {
			case 1:				
				node = new Node(
					static_cast<unsigned int>(ii), 
					grid_reader.get_vertex_x_coord(ii)
				);	
				break;					
			case 2:
				node = new Node(
					static_cast<unsigned int>(ii), 
					grid_reader.get_vertex_x_coord(ii),
					grid_reader.get_vertex_y_coord(ii)
				);	
				break;
			case 3:
				node = new Node(
					static_cast<unsigned int>(ii), 
					grid_reader.get_vertex_x_coord(ii),
					grid_reader.get_vertex_y_coord(ii),
					grid_reader.get_vertex_z_coord(ii)
				);	
				break;
			default:
				TDKP_GENERAL_EXCEPTION("invalid dimension " << dimension);
		}
		// set vertex index
		node->set_index_vertex(ii);
		vertices[ii] = node;		
	}	
	// serial adding to geometry object
	for(unsigned int ii = 0; ii < vertices.size(); ii++) {
		geometry->add_node(vertices[ii]);	
	}
	//TimeMeasurements::get_instance().stop("vertex nodes"); 
	
	

	// ------------------------------------------------
	// build element objects
	// ------------------------------------------------
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: creating elements");
	//TimeMeasurements::get_instance().start("elements");
	vector<Element*> elements(grid_reader.get_num_elements(), 0);
	#pragma omp parallel for
	for(int ii = 0; ii < (signed)grid_reader.get_num_elements(); ii++) {				
		Element::ElementShape elem_shape = grid_reader.get_element_shape(ii);
		const unsigned int nvertex = grid_reader.get_num_vertices_in_element(ii);
		Element* elem  = Element::factory(elem_shape, ii, dimension, polynom_degree);
		// --------------------------------------------
		// set corner nodes/vertices
		// --------------------------------------------
		for(unsigned int jj = 0; jj < nvertex; jj++) {
			unsigned int vidx = grid_reader.get_vertex_in_element(ii,jj);
			TDKP_BOUNDS_ASSERT(vidx < grid_reader.get_num_vertices(), "");
			elem->set_corner_node(jj, vertices[vidx]);				
		}			
		elements[ii] = elem;   		       				
	}
	for(unsigned int ii = 0; ii < elements.size(); ii++) {
		geometry->add_element(elements[ii]);	
	}
	//TimeMeasurements::get_instance().stop("elements");
	
	// ------------------------------------------------
	// use elements to build edges and face list (if required)
	// ------------------------------------------------
	vector<VertexObject*> edges(vertices.size(), 0);
	vector<VertexObject*> faces(vertices.size(), 0);
	vector<VertexObject*> bnd_vertices(vertices.size(), 0);
	//TimeMeasurements::get_instance().start("boundary vertex edge face");
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: building boundary vertex, edges and faces list");
	// ------------------------------------------------
	// TODO: parallelize
	// ------------------------------------------------
	unsigned int nelements_enabled = 0;
	for(unsigned int ii = 0; ii < elements.size(); ii++) {
		// only enabled elements (non-contacts)
		if(elements[ii]->enabled()) {
			nelements_enabled++;
			vector<unsigned int> vertex_indices;
			// build face object
			VertexObject vobj;
			if(faces_required) {
				for(unsigned int ff = 0; ff < elements[ii]->get_num_faces(); ff++) {
					// build face tree
					elements[ii]->get_face(ff, vertex_indices);								
					vobj.init(vertex_indices);
					build_vertex_objects_tree(ii,faces,vobj);
				}
			}
			if(edges_required) {
				for(unsigned int ee = 0; ee < elements[ii]->get_num_edges(); ee++) {
					// build edge tree
					elements[ii]->get_edge(ee, vertex_indices);								
					vobj.init(vertex_indices);
					build_vertex_objects_tree(ii,edges,vobj);					
				}	
			}
			// 1D: boundary vertices, we do not enforce that elements must be numbered in a sequence,
			// so we play the same game for 1D grids and possibly also for 2D and 3D grids with
			// hermite functions 
			if(bnd_vertices_required) {
				for(unsigned int vv = 0; vv < elements[ii]->get_num_corners(); vv++) {
					vertex_indices.assign(1, elements[ii]->get_corner_node(vv).get_index_vertex());
					vobj.init(vertex_indices);
					build_vertex_objects_tree(ii,bnd_vertices,vobj);
				}			
			}
		}	
	}
	//TimeMeasurements::get_instance().stop("boundary vertex edge face");
	//TimeMeasurements::get_instance().start("global indicies");
	// ---------------------------------------------------
	// set global indices
	// ---------------------------------------------------
	unsigned int nfaces = 0;
	if(faces_required) {
		nfaces = set_global_indices_to_tree(faces);		
	}	
	unsigned int nedges = 0;
	if(edges_required) {
		nedges = set_global_indices_to_tree(edges);	
	}
	unsigned int nverts = 0;
	if(bnd_vertices_required) {
		nverts = set_global_indices_to_tree(bnd_vertices);
	}
	geometry->set_num_edges(nedges);
	geometry->set_num_faces(nfaces);
	//TimeMeasurements::get_instance().stop("global indicies");
	//TimeMeasurements::get_instance().start("boundary objects");
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: structure has " << nedges << " edges and " << nfaces << " faces.");
			
	// ------------------------------------------------
	// collect boundary vobjs and set location
	// ------------------------------------------------
	unsigned int nboundary = 0;
	vector<VertexObject*>* bnds;
	if(dimension == 1) {
		bnds = &bnd_vertices;	
	} else if(dimension == 2) {
		bnds = &edges;
	} else {
		bnds = &faces;	
	}
	vector<VertexObject*> vobj_boundaries;
	vobj_boundaries.reserve((*bnds).size());
	VertexObject* vobj = 0;
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: collecting boundary objects");	  
	for(unsigned int ii = 0; ii < bnds->size(); ii++) {
		vobj = (*bnds)[ii];								
		while(vobj != 0) {
			TDKP_BOUNDS_ASSERT(vobj->get_num_elements() > 0, "");
			if(vobj->get_num_elements() == 1) {
				vobj->set_location('e'); // external	
			} else {
				vobj->set_location('i'); // internal
			}
			vobj_boundaries.push_back(vobj);
			vobj = vobj->get_next();				
		}	
	}
	nboundary = vobj_boundaries.size();
	//TimeMeasurements::get_instance().stop("boundary objects");
	//TimeMeasurements::get_instance().start("vertex to face mapping edge location");
	// ------------------------------------------------
	// 3D elements need to know the location of the edge
	// but the location is bound to the face. 
	// ------------------------------------------------
	if(dimension == 3) {
		// ------------------------------------------
		// first, build vertex2face map as edges have
		// less vertices and therefore the 'fast' access
		// via lowest vertex index may not work
		// ------------------------------------------
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: building vertex to face mapping");
		vector<list<VertexObject*> > vertex2faces(vertices.size());
		// loop over all faces
		for(unsigned int ii = 0; ii < faces.size(); ii++) {
			VertexObject* fobj = faces[ii];
			while(fobj) {
				for(unsigned int nn = 0; nn	< fobj->get_num_vertices(); nn++) {
					vertex2faces[fobj->get_vertex_index(nn)].push_back(fobj);	
				}
				fobj = fobj->get_next();
			}	
		}
		
		
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: determing edge location");
		// for every edge
		#pragma omp parallel for
		for(int ii = 0; ii < (signed)edges.size(); ii++) {
			VertexObject* eobj = edges[ii];			
			list<VertexObject*>::iterator fobj_it;
			list<VertexObject*>::iterator fobj_end;
			if(eobj) {										
				// still, for every edge
				while(eobj) {
					// set initial location
					TDKP_BOUNDS_ASSERT(eobj->get_num_vertices() == 2, "");	
					eobj->set_location('i');	
					// find face
					fobj_it  = vertex2faces[eobj->get_vertex_index(0)].begin();
					fobj_end = vertex2faces[eobj->get_vertex_index(0)].end();
					TDKP_BOUNDS_ASSERT(fobj_it != fobj_end, "");
					bool found = false;
					int  checked = 0;
					// check all nearby faces
					while(fobj_it != fobj_end) {
						TDKP_BOUNDS_ASSERT((*fobj_it)->get_num_vertices() == 3, "");
						// check if vertices in eobj are contained in fobj
						if((*fobj_it)->all_rhs_vertices_match(*eobj)) {
							found = true;
							// if fobj is exterior, tell this the edge and break
							if((*fobj_it)->get_location() == 'e') {
								eobj->set_location('e');
								break;	
							}
						}
						checked++; 					
						fobj_it++;
					}
					TDKP_ASSERT(found, "could not find face for edge! checked " << checked << ". problems by edge " << eobj->get_index_global());
					TDKP_ASSERT(eobj->get_location() == 'i' || eobj->get_location() == 'e', "");				
					eobj = eobj->get_next();	
				}
			}		
		}
		
	} 	
	//TimeMeasurements::get_instance().stop("vertex to face mapping edge location");	
	// ------------------------------------------------
	// build element boundary objects
	// ------------------------------------------------
	TDKP_ASSERT(nboundary > nelements_enabled, "");	
	vector<ElementBoundary*> boundaries(nboundary, 0);
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: creating element boundaries");
	//TimeMeasurements::get_instance().start("element boundaries");		
	#pragma omp parallel 
	{	
		VertexObject* omp_vobj = 0;
		#pragma omp for
		for(int bb = 0; bb < (signed)nboundary; bb++) {	
			omp_vobj = vobj_boundaries[bb];
			// determine shape of boundary
			ElementBoundary::BoundaryType btype = ElementBoundary::vertex;
			if(dimension == 2) {
				btype = ElementBoundary::edge;
			} else if(dimension == 3) {
				btype = ElementBoundary::face;	
			}
			// create boundary and set corners
			ElementBoundary* boundary = new ElementBoundary(btype, bb, omp_vobj->get_num_vertices());		
			for(unsigned int ii = 0; ii < omp_vobj->get_num_vertices(); ii++) {			
				boundary->set_corner_node(ii, vertices[omp_vobj->get_vertex_index(ii)]);
			}
			// set adjacent elements
			for(unsigned int ii = 0; ii < omp_vobj->get_num_elements(); ii++) {
				boundary->add_element(elements[omp_vobj->get_element_index(ii)]);	
			}			
			// one element: thats an external boundary
			if(boundary->get_num_elements() == 1) {
				boundary->set_location('e');			
			} else {
				TDKP_BOUNDS_ASSERT(boundary->get_num_elements() == 2, "");
				boundary->set_location('i');	
			}
			boundaries[bb] = boundary;			
		}
	}			
	
	// ------------------------------------------------ 
	// add boundaries, set location to vertex
	// ------------------------------------------------
	for(unsigned int bb = 0; bb < nboundary; bb++) {	
		geometry->add_element_boundary(boundaries[bb]);
		ElementBoundary& ebnd = *boundaries[bb];
		char bnd_location = ebnd.get_location();
		TDKP_ASSERT(bnd_location != 'u', "")		
		for(unsigned int vv = 0; vv < ebnd.get_num_vertices(); vv++) {
			// external face -> external vertex
			Node& vert = ebnd.get_vertex(vv); 				
			if(bnd_location == 'e') {
				vert.set_location(location_exterior);
				if(bnd_vertices_required) {
					bnd_vertices[vert.get_index_vertex()]->set_location(bnd_location);			
				}
			// internal face -> internal vertex if not already exteranl ...			
			} else if(vert.get_location() != location_exterior) {
				vert.set_location(location_interior);
				if(bnd_vertices_required) {
					TDKP_ASSERT(bnd_vertices[vert.get_index_vertex()] != 0, "bnd_vertices[vert.get_index_vertex()] != 0 failed for vertex " << vert.get_index_vertex());					
					bnd_vertices[vert.get_index_vertex()]->set_location('i');			
				}										
			}
		}
	}
	//TimeMeasurements::get_instance().stop("element boundaries");
	//TimeMeasurements::get_instance().start("setting element boundaries");	
	
	// --------------------------------------------------
	// set boundaries to element and location to vertices
	// --------------------------------------------------
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: setting element boundaries to elements");
	for(unsigned int ee = 0; ee < elements.size(); ee++) {
		if(elements[ee]->enabled()) {
			vector<unsigned int> involved_vertices;
			VertexObject bvobj;		
			// for every element boundary
			for(unsigned int bb = 0; bb < elements[ee]->get_num_boundaries(); bb++) {
				elements[ee]->get_element_boundary_vertices(bb, involved_vertices);
				bvobj.init(involved_vertices);
				// locate boundary
				VertexObject* located_bvobj = (*bnds)[bvobj.get_lowest_vertex_index()]->locate(bvobj);
				TDKP_ASSERT(located_bvobj != 0, "");			
				elements[ee]->set_element_boundary(
					bb,
					boundaries[located_bvobj->get_index_global()]
				);
			}
		}
	}
	//TimeMeasurements::get_instance().stop("setting element boundaries");
	//TimeMeasurements::get_instance().start("additional nodes");
	
	// ------------------------------------------------
	// create additional nodes
	// ------------------------------------------------
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: creating additional nodes");
	int current_node_index = vertices.size(); // additional nodes start there
	for(unsigned int ee = 0; ee < elements.size(); ee++) {
		if(elements[ee]->enabled()) {						
			VertexObject tmp_vobj;
			VertexObject* located_vobj = 0;
			vector<unsigned int> involved_vertices;
			vector<double>       node_coords;
			unsigned int         tag = 0;
			Element::AdditionalNodeLocation location;
			// for every addiational node 
			for(unsigned int nn = 0; nn < elements[ee]->get_num_additional_nodes(); nn++) {
				// get locator
				elements[ee]->get_additional_node_locator(nn, location, involved_vertices, node_coords, tag);
				// find vertex/edge/face
				Node* node = 0;
				switch(location) {
					case Element::edge_node:
						// build vobj for search
						tmp_vobj.init(involved_vertices);
						located_vobj = edges[tmp_vobj.get_lowest_vertex_index()]->locate(tmp_vobj);
						// must not be zero
						TDKP_ASSERT(located_vobj != 0, "");					
						// use obj to create / get node
						elements[ee]->set_additional_node(
							nn,
							located_vobj->locate_node(*geometry, current_node_index, node_coords, tag, 1.0e-6)
						);					
						break;
					case Element::inner_node:
						TDKP_ASSERT(tag == 0, "tagged inner nodes are not handled");
						// only belongs to element. so create and add to geometry and element					
						if(dimension == 1) {							
							node = new Node(static_cast<unsigned int>(current_node_index++), node_coords[0]);
						} else if(dimension == 2) {											
							node = new Node(static_cast<unsigned int>(current_node_index++), node_coords[0], node_coords[1]);	
						} else {
							node = new Node(static_cast<unsigned int>(current_node_index++), node_coords[0], node_coords[1], node_coords[2]);	
						}
						// set vertex index
						node->set_location(location_interior);
						geometry->add_node(node);
						elements[ee]->set_additional_node(nn,node);
						break;					 					
					case Element::vertex_node:
						// build vobj for search
						TDKP_BOUNDS_ASSERT(involved_vertices.size() == 1, "involved_vertices.size() == 1, but is " << involved_vertices.size());
						TDKP_BOUNDS_ASSERT(involved_vertices[0] < bnd_vertices.size(), "");						
						located_vobj = bnd_vertices[involved_vertices[0]];
						// must not be zero
						TDKP_ASSERT(located_vobj != 0, "located_vobj != 0 failed for vertex " << involved_vertices[0]);					
						// use obj to create / get node
						elements[ee]->set_additional_node(
							nn,
							located_vobj->locate_node(*geometry, current_node_index, node_coords, tag, 1.0e-6)
						);								
						break;						
					default:
						TDKP_GENERAL_EXCEPTION("unhandled additional node location type. implement it");	
				}
			}
		}
	}
		
	// a second order mesh has nodes at vertices and edges, 
	if(polynom_degree == 2) {
		if(dimension > 1) {
			TDKP_ASSERT(nedges + vertices.size() == geometry->get_num_nodes(), "nedges + vertices.size() == geometry->get_num_nodes() failed. but 2nd order elements with nodes on edges should have that!");
		} else {
			TDKP_ASSERT(nelements_enabled + vertices.size() == geometry->get_num_nodes(), "nelements_enabled + vertices.size() == geometry->get_num_nodes() failed. but 1d 2nd order elements should have that!");				
		}	
	}		
	//TimeMeasurements::get_instance().stop("additional nodes");
	
		
	// ------------------------------------------------
	// build regions and set region information to elements	
	// ------------------------------------------------
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: setting region to elements");	
	for(unsigned int rr = 0; rr < grid_reader.get_num_regions(); rr++) {
		Region* region = new Region(grid_reader.get_region_name(rr)); 		
		region->set_index_global(rr);
		region->set_material_name(grid_reader.get_region_material_name(rr).c_str());
		if(string("Contact").compare(region->get_material_name()) == 0) {
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "disabled region " << region->get_name());				
			region->set_enabled(false);	
		}		
		// ---------------------------------------------
		// set region to elements
		// ---------------------------------------------
		for(unsigned int jj = 0; jj < grid_reader.get_num_elements_in_region(rr); jj++) {
			elements[grid_reader.get_element_index_global(rr,jj)]->set_region(region,jj);
		}
		region->set_num_elements(grid_reader.get_num_elements_in_region(rr));
		geometry->add_region(region);
	}
			
	sout.str("");
	sout << "InputParser: final grid has\n"
	     << "  " << setw(18) << "nodes: " << geometry->get_num_nodes() << "\n"
	     << "  " << setw(18) << "edges: " << geometry->get_num_edges() << "\n"
	     << "  " << setw(18) << "faces: " << geometry->get_num_faces();			
	TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());
					
	// ------------------------------------------------
	// prepare geometry and finish
	// ------------------------------------------------
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "InputParser: preparing geometry");
	geometry->prepare();
	
	TimeMeasurements::get_instance().stop("building of FE grid");
	
	// ------------------------------------------------
	// kill linked lists
	// ------------------------------------------------
	for(unsigned int ii = 0; ii < edges.size(); ii++) {
		if(edges[ii] != 0) {
			delete edges[ii];	
		}	
	}
	for(unsigned int ii = 0; ii < edges.size(); ii++) {
		if(faces[ii] != 0) {
			delete faces[ii];	
		}	
	}
	for(unsigned int ii = 0; ii < edges.size(); ii++) {
		if(bnd_vertices[ii] != 0) {
			delete bnd_vertices[ii];	
		}	
	}		
	
	return geometry;

} 

InputParser::VertexObject::VertexObject(const VertexObject& copy) { 
	TDKP_GENERAL_EXCEPTION("no copy!"); 
}

/** empty constructor for vertex object */
InputParser::VertexObject::VertexObject()
: nverts(0),
  element_indices(0),
  lowest_vertex_index(0),
  next(0),
  index_global(-1),
  nodes(0),
  node_tags(0),
  location('u') // undefined location
{
}

/** destructor, deletes tree BUT does not touch created nodes! (they belong to the grid) */
InputParser::VertexObject::~VertexObject() {
	if(next) {
		delete next;	
		next = 0;
	}	
}
  			
/** create vertex object from supplied indices */  			
InputParser::VertexObject::VertexObject(const vector<unsigned int>& vertex_indices_)
: nverts(0),
  element_indices(0),
  lowest_vertex_index(0),
  next(0),
  index_global(-1),
  nodes(0),
  node_tags(0),
  location('u') // undefined location
{
	this->init(vertex_indices_);	
}
/** create clone of vertex object. allowed only for objects not linked in lists or having nodes */
InputParser::VertexObject* InputParser::VertexObject::clone() const {
	TDKP_BOUNDS_ASSERT(next == 0 && nodes.size() == 0 && element_indices.size() == 0, "clone not allowed on objects having nodes / elements or involved in linked lists");
	TDKP_BOUNDS_ASSERT(index_global == -1, "index global must not be set on vertex objects that may be cloned");
	VertexObject* ret   = new VertexObject();	
	ret->nverts         = nverts;
	// copy vertex indices
	memcpy(ret->vertex_indices, vertex_indices, sizeof(unsigned int) * nverts);	
	ret->lowest_vertex_index = lowest_vertex_index; 
	return ret;			
}
/** copy vertex list and find lowest index */
void InputParser::VertexObject::init(const vector<unsigned int>& vertex_indices_) {
	TDKP_BOUNDS_ASSERT(vertex_indices_.size() > 0 && vertex_indices_.size() < 4, "");
	TDKP_BOUNDS_ASSERT(next == 0 && nodes.size() == 0 && element_indices.size() == 0, "init not allowed on objects having nodes / elements or involved in linked lists");
	nverts = vertex_indices_.size();
	memcpy(vertex_indices, &vertex_indices_[0], sizeof(unsigned int) * nverts);	
	// find lowest index
	lowest_vertex_index = vertex_indices[0];
	for(unsigned int ii = 1; ii < nverts; ii++) {
		if(lowest_vertex_index > vertex_indices[ii]) {
			lowest_vertex_index = vertex_indices[ii];
		} 
	}		
}
void InputParser::VertexObject::add_element(unsigned int element_index_global) {
#ifdef DEBUG	
	for(unsigned int ii = 0; ii < element_indices.size(); ii++) {
		if(element_indices[ii] == element_index_global) {
			TDKP_GENERAL_EXCEPTION("element has already been added!");	
		}			
	}
#endif
	element_indices.push_back(element_index_global);	
}

unsigned int InputParser::VertexObject::get_num_elements() const {
	return element_indices.size(); 	
}
unsigned int InputParser::VertexObject::get_element_index(unsigned int ii) const {
	TDKP_BOUNDS_ASSERT(ii < element_indices.size(), "");
	return element_indices[ii];	
}
unsigned int InputParser::VertexObject::get_num_vertices() const {
	return nverts; 	
}
unsigned int InputParser::VertexObject::get_vertex_index(unsigned int ii) const {
	TDKP_BOUNDS_ASSERT(ii < nverts, "");
	return vertex_indices[ii];	
}
/** returns true if objects have the same vertex indices */
bool InputParser::VertexObject::compare_vertices(const VertexObject& rhs) const {
	// false if not the same number of vertices
	if(nverts != rhs.nverts) {
		return false;	
	}
	// check for each pair
	bool good;
	for(unsigned int ii = 0; ii < nverts; ii++) {
		good = false;		
		for(unsigned int jj = 0; jj < nverts; jj++) {
			if(vertex_indices[ii] == rhs.vertex_indices[jj]) {
				good = true;
				break;	
			}	
		}
		// no match, return false
		if(!good) {
			return false;	
		}	
	}
	// all matched, return true
	return true;	
}

bool InputParser::VertexObject::all_rhs_vertices_match(const VertexObject& rhs) const {
	// exception if lhs has less nodes than rhs
	if(nverts < rhs.nverts) {
		TDKP_GENERAL_EXCEPTION("lhs has less nodes than rhs!");	
	}
	// check for each pair
	bool good;
	for(unsigned int ii = 0; ii < rhs.nverts; ii++) {
		good = false;		
		for(unsigned int jj = 0; jj < nverts; jj++) {
			if(vertex_indices[jj] == rhs.vertex_indices[ii]) {
				good = true;
				break;	
			}	
		}
		// no match, return false
		if(!good) {
			return false;	
		}	
	}
	// all matched, return true
	return true;
}


// linked list functions
InputParser::VertexObject* InputParser::VertexObject::locate(const VertexObject& obj) {
	if(compare_vertices(obj)) {
		// it's me
		return this;	
	} else if(next) {
		// ask neighbour
		return next->locate(obj);	
	} else {
		// i'm leaf, so not found
		return 0;	
	}
} 
/** insert to the right */
void InputParser::VertexObject::insert(VertexObject* obj) {
	if(next) {
		next->insert(obj);	
	} else {
		next = obj;
	}		
}

// node creation functions
Node* InputParser::VertexObject::locate_node(
	Geometry& geometry,
	int& current_node_idx, 
	const vector<double>& coords, 
	unsigned int tag, 
	const double& tolerance
) {
	
	double dist;
	// -------------------------------------------
	// find matching node
	// -------------------------------------------
	for(unsigned int ii = 0; ii < nodes.size(); ii++) {
		dist = 0.0;
		double tmp;
		for(int dd = 0; dd < nodes[ii]->get_dimension(); dd++) {
			tmp = coords[dd] - nodes[ii]->get_coord(dd);
			dist += tmp * tmp;	
		}
		if(sqrt(dist) < tolerance) {
			// check if its the right tag
			if(node_tags[ii] == tag) {	
				return nodes[ii];
			}	
		}		
	}
	
	// -------------------------------------------
	// create new node
	// -------------------------------------------
	Node* node = 0;
	switch(coords.size()) {
		case 1:
			node = new Node(static_cast<unsigned int>(current_node_idx), coords[0]);
			current_node_idx++;
			break;
		case 2:
			node = new Node(static_cast<unsigned int>(current_node_idx), coords[0], coords[1]);
			current_node_idx++;
			break;
		case 3:
			node = new Node(static_cast<unsigned int>(current_node_idx), coords[0], coords[1], coords[2]);
			current_node_idx++;
			break;
		default: 
			TDKP_GENERAL_EXCEPTION("the existing of that kind of space dimension has not yet been prooved");	
	}	
	switch(location) {
		case 'i':
			node->set_location(location_interior);
			break;
		case 'e':
			node->set_location(location_exterior);
			break;
		default:			
			TDKP_GENERAL_EXCEPTION("unknown location " << location);
	}
	nodes.push_back(node);
	node_tags.push_back(tag);
	geometry.add_node(node);
	return node;		
}		
		
	
} // end of namespace tdkp
