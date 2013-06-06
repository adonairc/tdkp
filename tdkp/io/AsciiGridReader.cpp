/*
 * AsciiGrid.cpp
 *
 *  Created on: Mar 3, 2012
 *      Author: vepi
 */

#include "tdkp/io/AsciiGridReader.h"
#include "tdkp/geometry/Element.h"
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace boost::iostreams;

namespace tdkp {

#define CHECK(recv, expect) { if(recv != expect) { ostringstream __expcstr; __expcstr << "Error in file " << filename << " on line " << line_count << ": Expected " << expect << ", found " << recv; TDKP_THROW_EXCEPTION(__expcstr.str()); }}
#define PARSE_ASSERT(cond) { if(!(cond)) { ostringstream __expcstr; __expcstr << "Error in file " << filename << " on line " << line_count << ": Condition " << #cond << " failed"; TDKP_THROW_EXCEPTION(__expcstr.str()); } }

struct AsciiGridImpl {

	struct LVertex {
		LVertex() {
			memset(coords, 0, 3*sizeof(double));
		}
		double coords[3];
	};
	struct LElement {
		Element::ElementShape shape;
		vector<unsigned int> vertices;
		unsigned int region_index;
		LElement(Element::ElementShape shape_, int num_corners, int region_index_) {
			shape = shape_;
			vertices.resize(num_corners);
			region_index = region_index_;
		}
	};
	struct LRegion {
		string name;
		string material;
		vector<unsigned int> elements_global_indexes;
		LRegion(string name_, string material_) {
			name = name_;
			material = material_;
		}
	};
	string filename;
	unsigned int line_count;
	unsigned int dimension;
	vector<LVertex> vertices;
	vector<LElement> elements;
	vector<LRegion> regions;

	AsciiGridImpl(const char* filename_)
	: filename(filename_),
	  line_count(1),
	  dimension(0)
	{
		readfile();
	}

	/** return dimension of grid */
	unsigned int get_dimension() const {
		return dimension;
	}
	/** return number of vertices in geometry */
	unsigned int get_num_vertices() const {
		return vertices.size();
	}
	/** return number of elements in geometry */
	unsigned int get_num_elements() const {
		return elements.size();
	}
	/** return number of elements per region in geometry */
	unsigned int get_num_elements_in_region(unsigned int region_index_global) const {
		TDKP_ASSERT(regions.size() > region_index_global,"");
		return regions[region_index_global].elements_global_indexes.size();
	}
	/** return number of regions in geometry */
	unsigned int get_num_regions() const {
		return regions.size();
	}
	/** return region name */
	string get_region_name(unsigned int region_index_global) const {
		TDKP_ASSERT(regions.size() > region_index_global,"");
		return regions[region_index_global].name;
	}
	/** return material name of region */
	string get_region_material_name(unsigned int region_index_global) const {
		TDKP_ASSERT(regions.size() > region_index_global,"");
		return regions[region_index_global].material;
	}
	/** return global element index of element in region */
	unsigned int get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const {
		TDKP_ASSERT(regions.size() > region_index_global,"");
		TDKP_ASSERT(regions[region_index_global].elements_global_indexes.size() > element_index_region,"");
		return regions[region_index_global].elements_global_indexes[element_index_region];
	}
	/** return element shape */
	Element::ElementShape get_element_shape(unsigned int element_index) const {
		TDKP_ASSERT(elements.size() > element_index,"");
		return elements[element_index].shape;
	}
	/** return number of vertices in element */
	unsigned int get_num_vertices_in_element(unsigned int element_index) const {
		TDKP_ASSERT(elements.size() > element_index,"");
		return elements[element_index].vertices.size();
	}
	/** return global vertex index of local element vertex */
	unsigned int get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const {
		TDKP_ASSERT(elements.size() > element_index,"");
		TDKP_ASSERT(elements[element_index].vertices.size() > local_vertex_index,"");
		return elements[element_index].vertices[local_vertex_index];
	}
	/** return vertex coordinate x */
	double get_vertex_x_coord(unsigned int vertex_index) const {
		TDKP_ASSERT(vertices.size() > vertex_index,"");
		return vertices[vertex_index].coords[0];
	}
	/** return vertex coordinate x */
	double get_vertex_y_coord(unsigned int vertex_index) const {
		TDKP_ASSERT(vertices.size() > vertex_index,"");
		return vertices[vertex_index].coords[1];
	}
	/** return vertex coordinate x */
	double get_vertex_z_coord(unsigned int vertex_index) const {
		TDKP_ASSERT(vertices.size() > vertex_index,"");
		return vertices[vertex_index].coords[2];
	}
	/** return unique identifier for grid file */
	string get_unique_identifier() const {
		return filename;
	}

	template<class T>
	void read_stream(T& in) {

		string line;
		istringstream sstream;
		int stage = 0;
		int nelem = -1;
		int nvert = -1;
		bool skip = false;

		while(std::getline(in, line)) {
			if(line.size() > 0 && line[0] != '#') {
				if(stage < 2) {
					TDKP_LOGMSG(LOG_INFO_DEVEL2, "Stage is " << stage << " reading: " << line);
				}
				skip = false;
				sstream.clear(); sstream.str(line);
				if(stage == 0 && read_dimension(sstream)) {
					stage++;
					skip = true;
				}
				if(stage == 1 && !skip && read_regions(sstream)) {
					stage++;
					sstream.clear(); sstream.str(line);
				}
				if(stage == 2 && !skip && read_vertices(sstream, nvert)) {
					stage++;
					skip = true;
				}
				if(stage == 3 && !skip && read_elements(sstream, nelem)) {
					break;
				}
			}
			line_count++;
		}
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "Found " << elements.size() << " elements");
		TDKP_ASSERT(nelem == (int)elements.size(), "");
	}

	bool read_dimension(istream& in) {
		string cmd;
		in >> cmd;
		CHECK(cmd, "DIMENSION");
		in >> dimension;
		PARSE_ASSERT(dimension > 0 && dimension <= 3);
		return true;
	}

	bool read_regions(istream& in) {
		string cmd, region_name, region_material;
		in >> cmd;
		if(cmd != "REGION") {
			return true;
		}
		in >> region_name;
		PARSE_ASSERT(region_name.size() > 0);
		in >> region_material;
		PARSE_ASSERT(region_material.size() > 0);
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "Found region '" << region_name << "' consisting of '" << region_material << "'");
		regions.push_back(LRegion(region_name, region_material));
		return false;
	}

	bool read_vertices(istream& in, int& nvert) {
		if(nvert < 0) {
			string cmd;
			in >> cmd;
			CHECK(cmd, "VERTICES");
			in >> nvert;
			PARSE_ASSERT(nvert > 0);
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "Grid has " << nvert << " vertices");
			return false;
		} else {
			LVertex vert;
			for(unsigned int dd = 0; dd < dimension; dd++) {
				in >> vert.coords[dd];
			}
			vertices.push_back(vert);
			return (int)vertices.size() == nvert;
		}
	}

	bool read_elements(istream& in, int& nelem) {
		if(nelem < 0) {
			string cmd;
			in >> cmd;
			CHECK(cmd, "ELEMENTS");
			in >> nelem;
			PARSE_ASSERT(nelem > 0);
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "Grid has " << nelem << " elements");
			return false;
		} else {
			int region_index;
			Element::ElementShape shape;
			int num_vertices, vertex_index, tmp;
			in >> region_index;
			PARSE_ASSERT(region_index >= 0 && region_index < (int)regions.size());
			in >> tmp; shape = (Element::ElementShape)tmp;
			PARSE_ASSERT(shape < Element::num_element_shapes);
			in >> num_vertices;
			PARSE_ASSERT(num_vertices > 0);
			LElement elem(shape, num_vertices, region_index);
			for(int vv = 0; vv < num_vertices; vv++) {
				in >> vertex_index;
				PARSE_ASSERT(vertex_index >= 0 && vertex_index < (int)vertices.size());
				elem.vertices[vv] = vertex_index;
			}
			regions[region_index].elements_global_indexes.push_back(elements.size());
			elements.push_back(elem);
			return (int)elements.size() == nelem;
		}
	}

	void readfile() {
		if(filename.substr(filename.size() - 2) == "gz") {
			ifstream file(filename.c_str(), ios_base::in | ios_base::binary);
			boost::iostreams::filtering_istream in;
			in.push(gzip_decompressor());
			in.push(file);
			read_stream(in);
		} else {
			ifstream fin(filename.c_str());
			read_stream(fin);
		}
	}
};

AsciiGridReader::AsciiGridReader(const char * filename) {
	impl = new AsciiGridImpl(filename);
}

AsciiGridReader::~AsciiGridReader() {
	delete impl;
}

/** return dimension of grid */
unsigned int AsciiGridReader::get_dimension() const {
	return impl->get_dimension();
}
/** return number of vertices in geometry */
unsigned int AsciiGridReader::get_num_vertices() const {
	return impl->get_num_vertices();
}
/** return number of elements in geometry */
unsigned int AsciiGridReader::get_num_elements() const {
	return impl->get_num_elements();
}
/** return number of elements per region in geometry */
unsigned int AsciiGridReader::get_num_elements_in_region(unsigned int region_index_global) const {
	return impl->get_num_elements_in_region(region_index_global);
}
/** return number of regions in geometry */
unsigned int AsciiGridReader::get_num_regions() const {
	return impl->get_num_regions();
}
/** return region name */
string AsciiGridReader::get_region_name(unsigned int region_index_global) const {
	return impl->get_region_name(region_index_global);
}
/** return material name of region */
string AsciiGridReader::get_region_material_name(unsigned int region_index_global) const {
	return impl->get_region_material_name(region_index_global);
}
/** return global element index of element in region */
unsigned int AsciiGridReader::get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const {
	return impl->get_element_index_global(region_index_global, element_index_region);
}
/** return element shape */
Element::ElementShape AsciiGridReader::get_element_shape(unsigned int element_index) const {
	return impl->get_element_shape(element_index);
}
/** return number of vertices in element */
unsigned int AsciiGridReader::get_num_vertices_in_element(unsigned int element_index) const {
	return impl->get_num_vertices_in_element(element_index);
}
/** return global vertex index of local element vertex */
unsigned int AsciiGridReader::get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const {
	return impl->get_vertex_in_element(element_index, local_vertex_index);
}
/** return vertex coordinate x */
double AsciiGridReader::get_vertex_x_coord(unsigned int vertex_index) const {
	return impl->get_vertex_x_coord(vertex_index);
}
/** return vertex coordinate x */
double AsciiGridReader::get_vertex_y_coord(unsigned int vertex_index) const {
	return impl->get_vertex_y_coord(vertex_index);
}
/** return vertex coordinate x */
double AsciiGridReader::get_vertex_z_coord(unsigned int vertex_index) const {
	return impl->get_vertex_z_coord(vertex_index);
}
/** return unique identifier for grid file */
string AsciiGridReader::get_unique_identifier() const {
	return impl->get_unique_identifier();
}


}
