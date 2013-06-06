/*
 * MEDGridReader.cpp
 *
 *  Created on: Jul 22, 2010
 *      Author: vepi
 */

#include "tdkp/common/all.h"
#include "tdkp/io/MEDGridReader.h"
#include "MEDMEM_Med.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Field.hxx"

namespace tdkp {

MEDGridReader::MEDGridReader(const char* cfilename) :
	filename(cfilename), my_med(0), my_mesh(0)
{
	try {

		// -------------------------------------------------
		// open med file
		// -------------------------------------------------
		my_med = new MED(MED_DRIVER, filename);
		TDKP_ASSERT(my_med->getNumberOfMeshes() == 1, "MEDGridReader: your .med file has " << my_med->getNumberOfMeshes() << " meshes. Currently only files containing one mesh are supported");

		// -------------------------------------------------
		// get grid
		// -------------------------------------------------
		my_med->read(0);
		my_mesh = my_med->getMesh(my_med->getMeshNames().front());

		// -------------------------------------------------
		// get the already existing information
		// -------------------------------------------------
		num_elements 		= my_mesh->getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS);
		connectivity 		= my_mesh->getConnectivity(MED_FULL_INTERLACE, MED_NODAL, MED_CELL, MED_ALL_ELEMENTS);
		connectivity_index  = my_mesh->getConnectivityIndex(MED_NODAL, MED_CELL);
		coordinates 		= my_mesh->getCoordinates(MED_FULL_INTERLACE);
		space_dim 			= my_mesh->getSpaceDimension();
		num_vertices        = my_mesh->getNumberOfNodes();

		// -------------------------------------------------
		// scan for invalid elements
		// -------------------------------------------------
		for(unsigned int ii = 0; ii < num_elements; ii++) {
			this->get_element_shape(ii); // throws exception if shape is invalid!
		}

		// -------------------------------------------------
		// check that every element belongs to exactly one group
		// and build region to elements map
		// -------------------------------------------------
		vector<int>	elem_to_region_map(num_elements, -1);
		const vector<GROUP*> groups = my_mesh->getGroups(MED_CELL);
		region_names.resize(groups.size(), "unnamed");
		region_materials.resize(groups.size(), "");
		regions.resize(groups.size());
		for(unsigned int ii = 0; ii < groups.size(); ii++)
		{
			const GROUP* group = groups[ii];
			// set name
			region_names[ii] = group->getName();
			// get num elements in region
			unsigned int num_region_elements = group->getNumberOfElements(MED_ALL_ELEMENTS);
			// allocate space for region to element map
			regions[ii].resize(num_region_elements);
			// check if group covers all elements
			if(group->isOnAllElements()) {
				if(ii > 0) {
					TDKP_GENERAL_EXCEPTION("MEDGridReader: group no. " << ii << " '" << group->getName() << "' is defined on all elements but there is already also another group existing. Groups may not intersect!");
				}
				for(unsigned int ee = 0; ee < num_region_elements; ee++) {
					regions[ii][ee] = ee;
					elem_to_region_map[ee] = ii;
				}
			} else {
				// get elements
				const int* elements = group->getNumber(MED_ALL_ELEMENTS);
				// store
				for(unsigned int ee = 0; ee < num_region_elements; ee++) {
					// ensure that there is no second assignment
					TDKP_ASSERT( \
						elem_to_region_map[elements[ee] - 1] == -1, \
						"MEDGridReader: Element " << elements[ee] - 1 << " is already assigned to region/group " \
						<< elem_to_region_map[elements[ee] - 1] \
						<< " but there is also a second assignment to region/group " << ii \
						<< ". Check your grid! Each element must belong to exactly one group!" \
					);
					// note: element numbering in salome starts at 1!
					elem_to_region_map[elements[ee] - 1] = ii;
					regions[ii][ee] = elements[ee] - 1;
				}
			}
		}
	} catch(MEDMEM::MEDEXCEPTION e) {
		TDKP_GENERAL_EXCEPTION("MEDGridReader: MED threw exception: " << e.what());
	}

}

MEDGridReader::~MEDGridReader() {
	delete_object(my_med);
	my_mesh = 0; // mesh belongs to my_med, so he will kill it
}

/** set the material name for the given region */
void MEDGridReader::set_region_material_name(unsigned int region_index_global, const string& material_name) {
	TDKP_ASSERT(region_index_global < region_names.size(), "");
	this->region_materials[region_index_global] = material_name;
}

/** set the material name for the given region */
void MEDGridReader::set_region_material_name(unsigned int region_index_global, const char* material_name) {
	TDKP_ASSERT(region_index_global < region_names.size(), "");
	this->region_materials[region_index_global] = material_name;
}

/** return dimension of grid */
unsigned int MEDGridReader::get_dimension() const {
	return this->space_dim;
}

/** return number of vertices in geometry */
unsigned int MEDGridReader::get_num_vertices() const {
	return this->num_vertices;
}

/** return number of elements in geometry */
unsigned int MEDGridReader::get_num_elements() const {
	return this->num_elements;
}

/** return number of elements per region in geometry */
unsigned int MEDGridReader::get_num_elements_in_region(unsigned int region_index_global) const {
	TDKP_BOUNDS_ASSERT(region_index_global < region_names.size(), "");
	return this->regions[region_index_global].size();
}

/** return number of regions in geometry */
unsigned int MEDGridReader::get_num_regions() const {
	return this->region_names.size();
}

/** return region name */
string MEDGridReader::get_region_name(unsigned int region_index_global) const {
	TDKP_BOUNDS_ASSERT(region_index_global < region_names.size(), "");
	return this->region_names[region_index_global];
}

/** return material name of region */
string MEDGridReader::get_region_material_name(unsigned int region_index_global) const {
	TDKP_BOUNDS_ASSERT(region_index_global < region_names.size(), "");
	TDKP_ASSERT(region_materials[region_index_global] != "", "MEDGridReader: There has no material been assigned to region nr. " << region_index_global << " (" << region_names[region_index_global] << ")");
	return region_materials[region_index_global];
}

/** return global element index of element in region */
unsigned int MEDGridReader::get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const {
	TDKP_BOUNDS_ASSERT(region_index_global < regions.size(), "");
	TDKP_BOUNDS_ASSERT(regions[region_index_global].size() > element_index_region, "");
	return regions[region_index_global][element_index_region];
}

/** return element shape */
Element::ElementShape MEDGridReader::get_element_shape(unsigned int element_index) const {
	TDKP_ASSERT(element_index < num_elements, "element_index < num_elements failed because " << element_index << " < " << num_elements);
	medGeometryElement type = my_mesh->getElementType(MED_CELL, element_index + 1);
	switch(type) {
	case MED_TRIA3:
		TDKP_ASSERT(space_dim == 2, "MEDGridReader: Element " << element_index << " is a 2D triangle, but the structure is in 3D. Please ensure that you only have 3D elements in your mesh");
		return Element::triangle;
	case MED_QUAD4:
		TDKP_ASSERT(space_dim == 2, "MEDGridReader: Element " << element_index << " is a 2D rectangle, but the structure is in 3D. Please ensure that you only have 3D elements in your mesh");
		return Element::rectangle;
	case MED_TETRA4:
		TDKP_ASSERT(space_dim == 3, "MEDGridReader: Element " << element_index << " is a 3D tetrahedron, but the structure is in 2D. Please ensure that you only have 3D elements in your mesh");
		return Element::tetrahedron;
	case MED_TRIA6:
	case MED_QUAD8:
		TDKP_GENERAL_EXCEPTION("MEDGridReader: Element " << element_index << " within mesh from file '" << filename << "' is a second order element. Higher order elements (and the required nodes) are generated internally by tdkp. Therefore, only meshes with vertices on the corners of the elements can be used!");
	default:
		TDKP_GENERAL_EXCEPTION("MEDGridReader: Element " << element_index << " has an unknown element shape (MED code: " << type << ")");
	}
}

/** return number of vertices in element */
unsigned int MEDGridReader::get_num_vertices_in_element(unsigned int element_index) const {
	TDKP_BOUNDS_ASSERT(element_index < num_elements, "element_index < num_elements failed because " << element_index << " < " << num_elements);
	return connectivity_index[element_index + 1] - connectivity_index[element_index];
}

/** return global vertex index of local element vertex */
unsigned int MEDGridReader::get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const {
	TDKP_BOUNDS_ASSERT(element_index < num_elements, "element_index < num_elements failed because " << element_index << " < " << num_elements);
	return connectivity[connectivity_index[element_index] + local_vertex_index - 1] - 1;
}

/** return vertex coordinate x */
double MEDGridReader::get_vertex_x_coord(unsigned int vertex_index) const {
	return this->coordinates[space_dim * vertex_index];
}
/** return vertex coordinate x */
double MEDGridReader::get_vertex_y_coord(unsigned int vertex_index) const {
	if(space_dim > 1) {
		return this->coordinates[space_dim * vertex_index + 1];
	}
	return 0.0;
}

/** return vertex coordinate x */
double MEDGridReader::get_vertex_z_coord(unsigned int vertex_index) const {
	if(space_dim > 2) {
		return this->coordinates[space_dim * vertex_index + 2];
	}
	return 0.0;
}

/** return unique identifier for grid file */
string MEDGridReader::get_unique_identifier() const {
	return filename;
}

}
