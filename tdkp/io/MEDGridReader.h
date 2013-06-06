/*
 * MEDGridReader.h
 *
 *  Created on: Jul 22, 2010
 *      Author: vepi
 */

#ifndef MEDGRIDREADER_H_
#define MEDGRIDREADER_H_

#include "BaseGridReader.h"

// --------------------------------------------
// fwd declarations
// --------------------------------------------
namespace MEDMEM {
	class MED;
	class MESH;
}

namespace tdkp {

/**
 * class for reading grid files stored in MED format
 *
 * the MED format seems to be salome's (salome-platform.org) main
 * data storage format. this class here provides the interface to read
 * 2D and 3D grid files created within salome.
 *
 * material types must be defined via associating them to the
 * corresponding region. this is not done via the .med files.
 * there are two reasons why not:
 * - med files are binary. changing the material type is therefore cumbersome
 *   and involves messing around in salome
 * - creating groups as the material name leads to the case that
 *   several groups with the same name may exist. thats not that
 *   nice ...
 */
class MEDGridReader: public BaseGridReader {
public:
	MEDGridReader(const char* filename);
	virtual ~MEDGridReader();

	// -------------------------------------------------
	// material setting
	// -------------------------------------------------
	/** set the material name for the given region */
	void set_region_material_name(unsigned int region_index_global, const string& material_name);
	/** set the material name for the given region */
	void set_region_material_name(unsigned int region_index_global, const char* material_name);

	// -------------------------------------------------
	// Base Grid Reader Interface
	// -------------------------------------------------
	/** return dimension of grid */
	virtual unsigned int get_dimension() const;
	/** return number of vertices in geometry */
	virtual unsigned int get_num_vertices() const;
	/** return number of elements in geometry */
	virtual unsigned int get_num_elements() const;
	/** return number of elements per region in geometry */
	virtual unsigned int get_num_elements_in_region(unsigned int region_index_global) const;
	/** return number of regions in geometry */
	virtual unsigned int get_num_regions() const;
	/** return region name */
	virtual string get_region_name(unsigned int region_index_global) const;
	/** return material name of region */
	virtual string get_region_material_name(unsigned int region_index_global) const;
	/** return global element index of element in region */
	virtual unsigned int get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const;
	/** return element shape */
	virtual Element::ElementShape get_element_shape(unsigned int element_index) const;
	/** return number of vertices in element */
	virtual unsigned int get_num_vertices_in_element(unsigned int element_index) const;
	/** return global vertex index of local element vertex */
	virtual unsigned int get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const;
	/** return vertex coordinate x */
	virtual double get_vertex_x_coord(unsigned int vertex_index) const;
	/** return vertex coordinate x */
	virtual double get_vertex_y_coord(unsigned int vertex_index) const;
	/** return vertex coordinate x */
	virtual double get_vertex_z_coord(unsigned int vertex_index) const;
	/** return unique identifier for grid file */
	virtual string get_unique_identifier() const;

private:
	string 			filename;				//!< the grids file name
	MEDMEM::MED* 	my_med;					//!< med object (i.e. the hdf5 file object containing the grid)
	MEDMEM::MESH*	my_mesh;				//!< the mesh in medmem format
	unsigned int    num_elements;           //!< number of elements
	const int* 		connectivity;			//!< the element to node connectivity
	const int* 		connectivity_index;		//!< where does the element connectivity starts in the connectivity array?
	const double* 	coordinates;			//!< node coordinates
	unsigned int    num_vertices;			//!< number of vertices
	unsigned int    space_dim;				//!< space dim (required for node coordinates)
	vector<string>	region_names;			//!< region (group) names
	vector<string>  region_materials;		//!< materials used for each region
	vector<vector<unsigned int> > regions;	//!< relation of elements to regions


};

}

#endif /* MEDGRIDREADER_H_ */
