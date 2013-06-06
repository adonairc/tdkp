/*
 * AsciiGrid.h
 *
 *  Created on: Mar 3, 2012
 *      Author: vepi
 */

#ifndef ASCIIGRID_H_
#define ASCIIGRID_H_

#include "tdkp/io/BaseGridReader.h"

namespace tdkp {

class AsciiGridImpl;

/** Very simple ASCII format reader for 1-3D grids */
class AsciiGridReader : public BaseGridReader {
public:
	AsciiGridReader(const char * filename);
	virtual ~AsciiGridReader();

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
	AsciiGridImpl* impl;


};

}

#endif /* ASCIIGRID_H_ */
