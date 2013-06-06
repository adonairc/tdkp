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

#include <vector>
#include <iostream>
#include "tdkp/geometry/Geometry.h"
#include "tdkp/geometry/Element3DTR.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/geometry/Region.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/geometry/BoundaryCondition.h"

using namespace tdkp;
using namespace std;


/** init empty */
Geometry::Geometry() {
	this->init(0,0,0,0);	
}

Geometry::Geometry(unsigned int dimension, unsigned int nelem, unsigned int nnode, unsigned int nreg) {
	this->init(dimension, nelem, nnode, nreg);	
}

Geometry::~Geometry() {	
	element_iterator eit;
	for(eit = elements.begin(); eit != elements.end(); eit++) {
		delete (*eit);
		(*eit) = 0;	
	}		
	node_iterator vit;
	for(vit = nodes.begin(); vit != nodes.end(); vit++) {
		delete (*vit);
		*vit = 0;
	}	
	region_iterator rit;
	for(rit = regions.begin(); rit != regions.end(); rit++) {
		delete (*rit);
		*rit = 0;	
	}
	boundary_iterator bit;
	for(bit = element_boundaries.begin(); bit != element_boundaries.end(); bit++) {
		delete (*bit);
		*bit = 0;	
	}	
	delete bc_handler; bc_handler = 0;
}
	
/** initialize geometry object and tell the number of expected elements and nodes
 * 
 * the number of elements and regions may differ as low-dimensional regions such 
 * as contacts are neglected
 *  
 * @param nelem_ number of dfise elements
 * @param nnode_ number of nodes
 * @param nreg_  number of dfise regions
 */
void Geometry::init(unsigned int dimension, unsigned int nelem_, unsigned int nnode_, unsigned int nreg_) {
	
	this->elements.reserve(nelem_);
	this->nodes.reserve(nnode_);		
	this->regions. reserve(nreg_);
	this->dimension = dimension;
	if(elements.capacity() < nelem_) {
		TDKP_GENERAL_EXCEPTION("not enough memory for element vector");			
	}
	if(nodes.capacity() < nnode_) {
		TDKP_GENERAL_EXCEPTION("not enough memory for node vector");
	}
	if(regions.capacity() < nreg_) {
		TDKP_GENERAL_EXCEPTION("not enough memory for region vector");	
	}
	this->nedges 		       = 0;
	this->nfaces               = 0;	
	this->next_elem            = 0;
	this->next_node            = 0;
	this->next_region          = 0;
	this->next_boundary        = 0;
	this->identifier           = "";
	this->bc_handler           = new BCDirichlet(*this);
	
}


/** prepare geometry for calculation
 * 
 * sets the internal indices of the nodes and calculates the jacobi matrices in the elements
 */
 void Geometry::prepare() {
 	this->set_internal_indices();

 	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "Geometry: preparing elements for calculation");
	#pragma omp parallel 
	{ 	
		#pragma omp for
	 	for(int ii = 0; ii < (signed)get_num_elements(); ii++) {
			get_element(ii).prepare();
			TDKP_BOUNDS_ASSERT(get_element(ii).is_ready(), "get_element(ii).is_ready() failed for element " << ii << " of type " << get_element(ii).get_element_unique_type_key());   		
	 	}
	 	#pragma omp for
	 	for(int ii = 0; ii < (signed)get_num_boundaries(); ii++) { 	
	 		get_element_boundary(ii).prepare();	 	
	 	}
	}
 }
 
void Geometry::prepare_boundaries() {
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "Geometry: preparing boundary integrals");
 	#pragma omp for
 	for(int ii = 0; ii < (signed)get_num_boundaries(); ii++) {
 		get_element_boundary(ii).prepare_integrals();	
 	} 	
} 
 
/** set internal indices to nodes 
 *
 * the internal index corresponds to the position of the node 
 * in the global matrix system. any node with index -1 is a boundary
 * node and therefore disregarded in the calculation
 */
void Geometry::set_internal_indices() {
	int current_index = 0;
	node_iterator it;
	for(it = nodes.begin(); it != nodes.end(); it++) {		
		if(bc_handler->ignore_node_fully(*(*it))) {
			(*it)->set_index_internal(-1);
		} else {
			(*it)->set_index_internal(current_index++);	
		} 	
	} 
}

/** verify parsed geometry
 * 
 * 1. check if 'announced' number of nodes, elements and regions corresponds to created
 *    size of nodes, elements and regions.
 * 2. check if there are no nodes with same index
 * 3. check if boundary conditions are set
 * 4. check if the global node index corresponds to the position it was added
 * 5. check if no element has two times the same node
 * 6. check if there are no elements with equal nodes
 */
bool Geometry::verify() const {
	// check sizes
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "Geometry: checking geometry");
	
	// -----------------------------------------------------
	// calculate number of boundary nodes and check whether all nodes have proper index
	// -----------------------------------------------------
	int current_index;  
	int bnode 		 = 0;
	int global_index = 0;
	node_const_iterator it, cit;
	for(it = nodes.begin(); it != nodes.end(); it++) {
		// check global index
		if(global_index != (*it)->get_index_global()) {
			Logger::get_instance()->emit(LOG_ERROR, 1000, "Geometry: inconsistent geometry - nodes were not added sequentially. node %d was added as %d", (*it)->get_index_global(), 	global_index);
			return false;
		}
		// check internal index
		current_index = (*it)->get_index_internal();
		if(current_index == -1) {
			bnode++;
		} else {			
			// check that no other node has my index
			for(cit = it + 1; cit != nodes.end(); cit++) {
				if(current_index == (*cit)->get_index_internal()) {
					Logger::get_instance()->emit(LOG_WARN, 1000, "Geometry: inconsistent geometry. node: %d has same internal index (%d) like node %d", (*it)->get_index_global(),(*it)->get_index_internal(), (*cit)->get_index_global());
					return false;					
				}
			}	
		}	
		global_index++;
	}
	ostringstream sout;
	sout << "Geometry: there are " <<  bnode << " nodes at boundary and " << get_num_nodes() - bnode << " nodes inside";
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	
	// -----------------------------------------------------
	// check elements
	// -----------------------------------------------------
	element_const_iterator eit, ecit;
	int elem_idx       = 0;
	int elem_idx_other = 0;
	for(eit = elements.begin(); eit != elements.end(); eit++) {
		// check if an element has two times the same node
		if(!(*eit)->verify()) {
			Logger::get_instance()->emit(LOG_ERROR, 1000, "Geometry: element %d claims to be invalid", elem_idx);
			return false;
		}
		elem_idx_other = elem_idx + 1;
		// check if there is another element with the same nodes
		for(ecit = eit + 1; ecit != elements.end(); ecit++) {
			if((*eit)->compare(*(*ecit))) {
				Logger::get_instance()->emit(LOG_ERROR, 1000, "Geometry: element %d has the same nodes as element %d", elem_idx, elem_idx_other);
				return false;
			}
			elem_idx_other++;			
		}
				
		
		elem_idx++;
	}
	// ---------------------------------------------------
	// check if there are any nonzero nodes
	// ---------------------------------------------------
	if(this->get_num_nonzero_nodes() == 0) {
		Logger::get_instance()->emit(LOG_ERROR, "Geometry: all nodes are on boundary and therefore we have nothing to calculate");
		return false; 	
	}
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, "Geometry: geometry looks good");
	return true;	
}

					
/** add node to geometry 
 * the passed node will be deleteted when the geometry object is deleted
 * nodes have to be added in the sequence of their global indices
 * @param node pointer to node object
 * @return index in vector nodes (== global index)
 */
unsigned int Geometry::add_node(Node* node) {
	nodes.push_back(node);
	return next_node++;		
}

/** add element to geometry
 * the passed element will be deleted when the geometry object is deleted
 * @param elem pointer to element object
 * @return index in vector elements (== global element index)
 */
unsigned int Geometry::add_element(Element* elem) {
	elements.push_back(elem);
	return next_elem++;	
}

/** add region to geometry
 * the passed region will be deleted when the geometry object is deleted
 * @param reg pointer to region object
 * @return region index
 */
unsigned int Geometry::add_region(Region* reg) {
	if(reg->get_num_elements() == 0) {
		ostringstream sout;
		sout << "region " << reg->get_name() << " has no elements assigned";
		TDKP_GENERAL_EXCEPTION(sout.str());
	}
	regions.push_back(reg);
	return next_region++;	
}
/** add element boundary to geometry
 * 
 * element boundaries can be nodes (for 1D elements), edges (for 2D elements)
 * or faces (for 3D elements)
 */
unsigned int Geometry::add_element_boundary(ElementBoundary* boundary) {
	this->element_boundaries.push_back(boundary);
	return next_boundary++;	
}
	
/** return region */
const Region& Geometry::get_region (unsigned int region_idx) const {
	TDKP_BOUNDS_ASSERT(region_idx >= 0 && region_idx < get_num_regions(), "");
	return *regions[region_idx];	
}

Region& Geometry::get_region (unsigned int region_idx) {
	return const_cast<Region&>(
		static_cast<const Geometry&>(*this).get_region(region_idx)
	);	
}

/** return element boundary */
const ElementBoundary& Geometry::get_element_boundary(unsigned int boundary_idx) const {
	TDKP_BOUNDS_ASSERT(boundary_idx < this->element_boundaries.size(), "");
	TDKP_BOUNDS_ASSERT(element_boundaries[boundary_idx] != 0, "boundary " << boundary_idx << " does not exist!");	
	return *element_boundaries[boundary_idx];
}
/** return element boundary */
ElementBoundary& Geometry::get_element_boundary(unsigned int boundary_idx) {
	return const_cast<ElementBoundary&>(
		static_cast<const Geometry&>(*this).get_element_boundary(boundary_idx)
	);	
}

/** return reference to node with global idx node_idx*/
Node& Geometry::get_node (unsigned int node_idx) {
	return const_cast<Node&>(
		static_cast<const Geometry&>(*this).get_node(node_idx)
	);	
}
/** return reference to element */
Element& Geometry::get_element(unsigned int elem_idx) {
	return const_cast<Element&>(
		static_cast<const Geometry&>(*this).get_element(elem_idx)
	);	
}

unsigned int Geometry::get_num_nonzero_nodes() const {
	int num = 0;
	node_const_iterator it;
	for(it = nodes.begin(); it != nodes.end(); it++) {
		if((*it)->get_index_internal() != -1) {			
			num++;	
		}
	}	
	return num;
}
 
void Geometry::set_num_faces(int nfaces_) {
	this->nfaces = nfaces_; 
	this->element_boundaries.reserve(max(this->nfaces, this->nedges));
} 

void Geometry::set_num_edges(int nedges_) {
	this->nedges = nedges_;
	this->element_boundaries.reserve(max(this->nfaces, this->nedges));	
}


int Geometry::find_next_node(double x, double y, double z) const {
	
	//cout << "lookinfg for node at " << x << ", " << y << ", " << z << "\n";
	
	double coords[]     = {x,y,z};
	double min_distance = 0.0;
	int    min_index;
	node_const_iterator it = this->nodes.begin(); 
	min_index = (*it)->get_index_global(); 
	double tmp;
	for(unsigned int ii = 0; ii < this->dimension; ii++) {
		tmp =  (*it)->get_coord(ii) - coords[ii];
		min_distance += tmp * tmp;
	}
	double distance;
	for(; it != this->nodes.end(); it++) {
		distance = 0.0; 
		const double* vcoords = (*it)->get_coords(); 
		for(unsigned int ii = 0; ii < this->dimension; ii++) {
			distance += (vcoords[ii] - coords[ii]) * (vcoords[ii] - coords[ii]);
		}
		if(min_distance > distance) {
			min_index    = (*it)->get_index_global(); 	
			min_distance = distance;
		}
	}
	cout << "next node (d:" << sqrt(min_distance) << "): " << *this->nodes[min_index] << "\n";
	
	return min_index;
}

void Geometry::list_elements(int node_index_global) const {
	for(element_const_iterator it = this->elements.begin(); it != this->elements.end(); it++) {
		for(unsigned int ii = 0; ii < (*it)->get_num_nodes(); ii++) {
			if((*it)->get_node(ii).get_index_global() == node_index_global) {
				cout << *(*it);
				break;
			}			
		}	
	}
}
void Geometry::rescale_node_coordinates(const double& rescale) {
	ostringstream sout;
	sout <<  "Geometry: rescaling node coordinates with a factor of " << rescale;
	Logger::get_instance()->emit(LOG_INFO, sout.str());
#pragma omp parallel for default(shared) schedule(static, 5000)	
	for(int vv = 0; vv < (signed)this->nodes.size(); vv++) {
		for(unsigned int ii = 0; ii < this->dimension; ii++) {
			this->nodes[vv]->get_coord(ii) *= rescale;	
		}
	}	
	this->prepare();
}


/** calculate the volume of the quantized region for kp problems 
 *
 * the quantized region is given by all regions with band edges below
 * the lowest conduction bandedge of a material that sits on the boundary  
 */
 /*
double Geometry::calculate_quantized_volume() const {

	TDKP_ASSERT(boundary_materials.size() > 0, "materials not set to object");
					
	// ------------------------------------------
	// get materials
	// ------------------------------------------
	vector<const Material*> materials(boundary_materials.size());
	for(unsigned int mm = 0; mm < get_num_regions(); mm++) {
		TDKP_ASSERT(get_region(mm).get_material().get_id() < materials.size(), "get_region(mm).get_material().get_id() " << get_region(mm).get_material().get_id() << " < materials.size()"  << materials.size());
		materials[get_region(mm).get_material().get_id()] = &(get_region(mm).get_material());  	
	}		 	
	
	// ------------------------------------------
	// find min bandedges of boundary material
	// ------------------------------------------
	double cb_max_edge      = 0.0;
	double vb_max_edge      = 0.0;
	bool   edges_set        = false;		
	for(unsigned int ii = 0; ii < boundary_materials.size(); ii++) {
		if(boundary_materials[ii]) {
			if(!edges_set || materials[ii]->get("conduction_band_edge") < cb_max_edge) {
				cb_max_edge = materials[ii]->get("conduction_band_edge"); 	
			}
			if(!edges_set || materials[ii]->get("valence_band_edge") > vb_max_edge) {
				vb_max_edge = materials[ii]->get("valence_band_edge"); 	
			}
			edges_set = true;
		}
	}
	
	// --------------------------------------------------
	// check which material belongs to quantized regions
	// --------------------------------------------------
	vector<bool> quantized_materials(boundary_materials.size(), false);
	for(unsigned int ii = 0; ii < boundary_materials.size(); ii++) {
		if(!boundary_materials[ii]) {
			if(materials[ii]->get("conduction_band_edge") < cb_max_edge  
			   && materials[ii]->get("valence_band_edge") > vb_max_edge) {
				quantized_materials[ii] = true;
			} else if(materials[ii]->get("conduction_band_edge") < cb_max_edge) {
				TDKP_LOGMSG(LOG_WARN, "material " << ii << " has cb quantized, but vb not. not counting it to quantized region."); 
			} else if(materials[ii]->get("valence_band_edge") > vb_max_edge) {
				TDKP_LOGMSG(LOG_WARN, "material " << ii << " has vb quantized, but cb not. not counting it to quantized region."); 				
			}					
		}
	}
											
	// ------------------------------------------
	// calculate volume
	// ------------------------------------------
	double volume = 0.0;
	for(element_const_iterator it = elements_begin();  
		it != elements_end(); it++) {
		if(!boundary_materials[(*it)->get_material().get_id()]) {		
			volume += (*it)->get_volume();	
		}				
	}	
	TDKP_ASSERT(volume > 0, "volume of quantized region is 0!"); 	
	return volume;	 
}*/

/** set materials to geometry */
void Geometry::set_materials(MaterialDatabase& material_database) {
	region_iterator it;
	for(it = this->regions_begin(); it != this->regions_end(); it++) {
		// ignore contact
		if((*it)->enabled()) {
			if(!material_database.material_exists((*it)->get_material_name().c_str())) {
				material_database.load_material((*it)->get_material_name().c_str());
			}
			// set material to region
			(*it)->set_material(material_database.get_material((*it)->get_material_name().c_str()));
		}
	}		
	update_boundary_materials(); 
}

/** update boolean array of boundary materials */
void Geometry::update_boundary_materials() {
					
	// ------------------------------------------
	// find num materials
	// ------------------------------------------
	unsigned int num_materials = 0;
	for(unsigned int mm = 0; mm < get_num_regions(); mm++) {
		if(get_region(mm).enabled()) {
			num_materials = max(num_materials, get_region(mm).get_material().get_id() + 1);
		}
	}	
		
	// ------------------------------------------
	// find boundary materials
	// ------------------------------------------
	boundary_materials.assign(num_materials, false);
	bool boundary = false;
	for(unsigned int ii = 0; ii < this->get_num_elements(); ii++) {
		boundary = false;
		const Element& element = get_element(ii);
		// only if we know thats this is not yet on the boundary
		if(element.enabled() && !boundary_materials[element.get_material().get_id()]) {
			for(unsigned int vv = 0; vv < element.get_num_nodes(); vv++) {			
				if(element.get_node(vv).get_index_internal() == -1) {
					boundary = true;	
				}
			}
			boundary_materials[element.get_material().get_id()] = boundary;
		}		
	}	
}
 
bool Geometry::boundary_material(unsigned int material_id) const {
	TDKP_BOUNDS_ASSERT(material_id < boundary_materials.size(), "material_id < boundary_conditions.size()");
	return boundary_materials[material_id];	
}	

/** set new boundary conditions to geometry */
void Geometry::set_boundary_conditions(BoundaryCondition* bc_handler_) {
	delete bc_handler;
	bc_handler = bc_handler_;	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "Geometry: setting new boundary conditions");
	this->update_boundary_materials();
	this->prepare();	
}

const BoundaryCondition& Geometry::get_boundary_conditions() const { 
	return *bc_handler; 
}

unsigned int Geometry::get_num_vertices() const {
	unsigned int nvertex = 0;
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		if(get_node(ii).get_index_vertex() != -1) {
			nvertex++;	
		}	
	}	
	return nvertex;
}

QuantizedVolumeCalculator::QuantizedVolumeCalculator(const Geometry& geometry_)
: geometry(geometry_) 
{
}

/** calculates the quantized volume by checking the minimal band edges on the outer boundary */	
double QuantizedVolumeCalculator::calculate_volume() const {
	
	const vector<bool>& boundary_materials = get_geometry().get_boundary_materials();	
	TDKP_ASSERT(boundary_materials.size() > 0, "materials not set to geometry object");
					
	// ------------------------------------------
	// get materials
	// ------------------------------------------
	vector<const Material*> materials(boundary_materials.size());
	for(unsigned int mm = 0; mm < get_geometry().get_num_regions(); mm++) {
		if(get_geometry().get_region(mm).enabled()) {
			TDKP_ASSERT(get_geometry().get_region(mm).get_material().get_id() < materials.size(), "get_geometry().get_region(mm).get_material().get_id() " << get_geometry().get_region(mm).get_material().get_id() << " < materials.size()"  << materials.size());
			materials[get_geometry().get_region(mm).get_material().get_id()] = &(get_geometry().get_region(mm).get_material());
		}  	
	}		 	
	
	// ------------------------------------------
	// find min bandedges of boundary material
	// ------------------------------------------
	double cb_max_edge      = 0.0;
	double vb_max_edge      = 0.0;
	bool   edges_set        = false;		
	for(unsigned int ii = 0; ii < boundary_materials.size(); ii++) {
		if(boundary_materials[ii] && materials[ii]->valid_key("conduction_band_edge") && materials[ii]->valid_key("valence_band_edge")) {
			if(!edges_set || materials[ii]->get("conduction_band_edge") < cb_max_edge) {
				cb_max_edge = materials[ii]->get("conduction_band_edge"); 	
			}
			if(!edges_set || materials[ii]->get("valence_band_edge") > vb_max_edge) {
				vb_max_edge = materials[ii]->get("valence_band_edge"); 	
			}
			edges_set = true;
		}
	}
	
	// --------------------------------------------------
	// check which material belongs to quantized regions
	// --------------------------------------------------
	vector<bool> quantized_materials(boundary_materials.size(), false);
	for(unsigned int ii = 0; ii < boundary_materials.size(); ii++) {
		if(!boundary_materials[ii]) {
			if(materials[ii]->valid_key("conduction_band_edge") && materials[ii]->valid_key("valence_band_edge")) {
				if(materials[ii]->get("conduction_band_edge") < cb_max_edge  
				   && materials[ii]->get("valence_band_edge") > vb_max_edge) {
					quantized_materials[ii] = true;
				} else if(materials[ii]->get("conduction_band_edge") < cb_max_edge) {
					TDKP_LOGMSG(LOG_WARN, "Geometry: material " << ii << " has cb quantized, but vb not. not counting it to quantized region."); 
				} else if(materials[ii]->get("valence_band_edge") > vb_max_edge) {
					TDKP_LOGMSG(LOG_WARN, "Geometry: material " << ii << " has vb quantized, but cb not. not counting it to quantized region."); 				
				}
			}					
		}
	}
											
	// ------------------------------------------
	// calculate volume
	// ------------------------------------------
	double volume = 0.0;
	for(Geometry::element_const_iterator it = get_geometry().elements_begin();  
		it != get_geometry().elements_end(); it++) {
		if((*it)->enabled() && !boundary_materials[(*it)->get_material().get_id()]) {		
			volume += (*it)->get_volume();	
		}				
	}
	if(volume <= 0.0) {
		TDKP_LOGMSG(LOG_WARN, "Geometry: volume of quantized region is " << volume << "!");	
	} 		 
	return volume;	 
}

UserDefinedQuantizedVolumeCalculator::UserDefinedQuantizedVolumeCalculator(const Geometry& geometry_)
: QuantizedVolumeCalculator(geometry_)
{
}

bool UserDefinedQuantizedVolumeCalculator::is_quantized(const string& region_name) const {
	if(find(quantized_region_names.begin(), 
			quantized_region_names.end(), 
			region_name) == quantized_region_names.end()) {
		return false;
	} 
	return true;				
}

/** calculate the volume of the quantized region defined by user input */
double UserDefinedQuantizedVolumeCalculator::calculate_volume() const {

	// ------------------------------------------
	// inform
	// ------------------------------------------

	TDKP_LOGMSG(LOG_INFO, "Geometry: calculating volume of quantized regions using user-defined region names.");
	ostringstream sout;
	sout << "Geometry: user-defined quantized region names:\n"; 		
	for(unsigned int ii = 0; ii < quantized_region_names.size(); ii++) {
		sout << "  " << quantized_region_names[ii] << "\n";
	}
	TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());
	// ------------------------------------------
	// build cheaper boolean array
	// ------------------------------------------
	vector<bool> region_quantized(get_geometry().get_num_regions());
	unsigned int num_quantized_regions = 0;
	for(unsigned int ii = 0; ii < get_geometry().get_num_regions(); ii++) {
		if(get_geometry().get_region(ii).enabled()) {
			region_quantized[ii] = is_quantized(get_geometry().get_region(ii).get_name());
			if(region_quantized[ii]) {
				num_quantized_regions++;	
			}
		}
	}
	TDKP_ASSERT(num_quantized_regions > 0, "no valid quantized region defined");
	// ------------------------------------------
	// loop over elements and calculate volume
	// ------------------------------------------
	double volume = 0.0;
	for(Geometry::element_const_iterator it = get_geometry().elements_begin();  
		it != get_geometry().elements_end(); it++) {
		if((*it)->enabled()) {
			if(region_quantized[(*it)->get_region().get_index_global()]) {		
				volume += (*it)->get_volume();	
			}
		}				
	}		
	TDKP_ASSERT(volume > 0, "volume of quantized region is 0!"); 	
	return volume;	
}
/** add new region that should be counted to the quantized regions */
void UserDefinedQuantizedVolumeCalculator::add(const char* quantized_region_name) {
	if(!is_quantized(string(quantized_region_name))) {
		quantized_region_names.push_back(string(quantized_region_name));			
	}	
}

