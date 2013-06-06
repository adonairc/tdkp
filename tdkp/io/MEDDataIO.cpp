/*
 * MEDDataIO.cpp
 *
 *  Created on: Jul 25, 2010
 *      Author: vepi
 */

#include "MEDDataIO.h"
#include "MEDMEM_Med.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Field.hxx"
#include "MEDMEM_Meshing.hxx"

namespace tdkp {

struct MEDDataIOPimpl {

	const Geometry& geometry;				//!< the geometry -> should correspond to the one in the med file
	MED*  			my_med;					//!< basic med file representation
	MESH*		 	my_mesh;				//!< the mesh (owned by my_med)
	int   			my_driver;				//!< idx for writing medmem files
	vector<int> 	idx_global_2_idx_med; 	//!< maps to geometry element id to the med mesh id

	/** write med data based on a geometry object
	 *
	 * this function actually creates a new file and stores the geometry
	 * into it. so creating such an object allows to convert any 2D/3D grid
	 * object into a .med style grid
	 */
	MEDDataIOPimpl(const Geometry& geometry_, const char* outputfile) :
		geometry(geometry_), my_mesh(0)
	{
		TDKP_ASSERT(geometry.get_dimension() > 1, "MEDDataIO: med output is only possible for 2D / 3D");

		try {
			my_med    = new MED();
			this->create_med_mesh_from_geometry();
			my_driver = my_med->addDriver(MED_DRIVER, outputfile);
			my_med->write(this->my_driver);
		} catch(MEDMEM::MEDEXCEPTION e) {
			TDKP_GENERAL_EXCEPTION("MEDGridReader: MED threw exception: " << e.what());
		}
	}

	~MEDDataIOPimpl() {
		my_med->write(this->my_driver);
		delete_object(my_med);
	}

	/** strips the directory from the varname (which is actually a file name) */
	string strip_directory(const string& varname) const;

	double adjust_value(const complex<double> & a) const {
		return abs(a)*abs(a);
	}
	double adjust_value(const double& b) const {
		return b;
	}

	/** creates new mesh from geometry object */
	void create_med_mesh_from_geometry();
	/** translates the tdkp element type into the med element type */
	medGeometryElement map_tdkp_to_med_elem_type(unsigned int shape) const throw(Exception*);
	/** builds the element idx global to idx med map */
	void create_elem_global_to_med_map();

	template<class T>
	void write_element(const ElementData<T>& data, const char* varname);
	template<class T>
	void write_nodal(const NodeData<T>& data, const char* varname);

};


/** strips the directory from the varname (which is actually a file name) */
string MEDDataIOPimpl::strip_directory(const string& varname) const {
	string::size_type pos = varname.find_last_of('/');
	if(pos != string::npos) {
		return varname.substr(pos + 1);
	} else {
		return varname;
	}
}

/** translates the tdkp element type into the med element type */
medGeometryElement MEDDataIOPimpl::map_tdkp_to_med_elem_type(unsigned int shape) const throw(Exception*) {
	switch(shape) {
	case Element::tetrahedron:
		return MED_TETRA4;
	case Element::triangle:
		return  MED_TRIA3;
	case Element::rectangle:
		return MED_QUAD4;
	case Element::line:
		return MED_SEG2;
	case Element::vertex:
		return MED_POINT1; break;
	default:
		TDKP_GENERAL_EXCEPTION("unknown element type. please fix me!");
	}
}

/** builds the element idx global to idx med map */
void MEDDataIOPimpl::create_elem_global_to_med_map() {

	// ---------------------------------------------
	// first problem: we only store elements of dimension space_dim.
	// so we have to build a map which maps the index_global into the index_med
	// ---------------------------------------------
	unsigned int space_dim = geometry.get_dimension();
	idx_global_2_idx_med.assign(geometry.get_num_elements(), -1);
	int midx = 1; // med index starts at 1
	for(int tt = 0; tt < Element::num_element_shapes; tt++) {
		for(unsigned int ee = 0; ee < geometry.get_num_elements(); ee++) {
			const Element& elem = geometry.get_element(ee);
			if(elem.get_dimension() == space_dim && elem.get_shape() == tt) {
				idx_global_2_idx_med[elem.get_index_global()] = midx++;
			}
		}
	}

}

/** creates new mesh from geometry object */
void MEDDataIOPimpl::create_med_mesh_from_geometry() {

	TDKP_BOUNDS_ASSERT(my_mesh == 0, "");
	MESHING* meshing = new MESHING();
	string identifier = geometry.get_identifier();
	if(identifier.size() > MED_TAILLE_NOM) {
		identifier = identifier.substr(0, MED_TAILLE_NOM);
	}
	meshing->setName(identifier);

	// ---------------------------------------------
	// 1. set coordinates
	// for this, we have to store them all in a huge array
	// ---------------------------------------------
	int space_dim = geometry.get_dimension();
	int num_verts = geometry.get_num_vertices();
	vector<double> coords(space_dim * num_verts, 0.0);
	meshing->setMeshDimension(space_dim);
	// loop over all nodes and find our vertices
	for(unsigned int nn = 0; nn < geometry.get_num_nodes(); nn++) {
		const Node& nd = geometry.get_node(nn);
		if(nd.get_index_vertex() != -1) {
			// store coords
			for(int ss = 0; ss < space_dim; ss++) {
				coords[space_dim * nd.get_index_vertex() + ss] = nd.get_coord(ss);
			}
		}
	}
	meshing->setCoordinates(space_dim, num_verts, &coords[0], "CARTESIAN", MED_FULL_INTERLACE);
	string cnames[3] = { "X","Y","Z" };
	meshing->setCoordinatesNames(cnames);
	string units[3] = { "nm","nm","nm" };
	meshing->setCoordinatesUnits(units);

	// ----------------------------------------------
	// count cell shapes (but only include elements
	// which have the correct dimensionality (discard
	// contacts) and build connectivity
	// ----------------------------------------------
	vector<unsigned int> cell_counts(Element::num_element_shapes, 0);
	vector<vector<int> > cell_connectivity(Element::num_element_shapes);
	unsigned int num_elements = 0;
	for(unsigned int ee = 0; ee < geometry.get_num_elements(); ee++) {
		const Element& elem = geometry.get_element(ee);
		if((signed)elem.get_dimension() == space_dim) {
			// store number of elements
			TDKP_BOUNDS_ASSERT(elem.get_shape() < (signed)cell_counts.size() && elem.get_shape() >= 0, "");
			cell_counts[elem.get_shape()]++;
			num_elements++;
			// store connectivity of the element (note, salome starts counting vertices at 1!)
			for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
				if(elem.get_node(nn).get_index_vertex() != -1) {
					// append vertex index
					cell_connectivity[elem.get_shape()].push_back(
						elem.get_node(nn).get_index_vertex() + 1
					);
				}
			}
		}
	}

	// ----------------------------------------------
	// build element type info
	// ----------------------------------------------
	vector<medGeometryElement> types;
	vector<int> reduced_counts; // thats cell_counts corresponding to types
	for(unsigned int tt = 0; tt < cell_counts.size(); tt++) {
		if(cell_counts[tt] > 0) {
			reduced_counts.push_back(cell_counts[tt]);
			types.push_back(this->map_tdkp_to_med_elem_type(tt));
		}
	}

	// ----------------------------------------------
	// write element type info to meshing object
	// ----------------------------------------------
	meshing->setNumberOfTypes(types.size(), MED_CELL);
	meshing->setTypes(&types[0], MED_CELL);
	meshing->setNumberOfElements(&reduced_counts[0], MED_CELL);

	// ----------------------------------------------
	// set connectivity
	// ----------------------------------------------
	unsigned int reduced_index = 0;
	for(unsigned int cc = 0; cc < cell_counts.size(); cc++) {
		if(cell_counts[cc] > 0) {
			meshing->setConnectivity(&cell_connectivity[cc][0], MED_CELL, types[reduced_index]);
			reduced_index++;
		}
	}

	// ---------------------------------------------
	// create map as our element numbering is not the
	// same as the one in the med file. med file sorts the
	// elements according to their shape
	// ---------------------------------------------
	this->create_elem_global_to_med_map();

	// ---------------------------------------------
	// create element groups for every region
	// ---------------------------------------------
	// ---------------------------------------------
	// the region information is actually not stored in
	// the geometry block directly, so we have to rebuild it from the
	// region assignment of each element
	// ---------------------------------------------
	// allocate storage
	vector<vector<vector<int> > > reg2shape2idxs;	// [region_idx][cell type][element_region_idx]
	reg2shape2idxs.resize(geometry.get_num_regions());
	for(unsigned int rr = 0; rr < geometry.get_num_regions(); rr++) {
		reg2shape2idxs[rr].resize(Element::num_element_shapes);
	}
	// build map (i.e. association between elements and regions, including the element type
	for(unsigned int ee = 0; ee < geometry.get_num_elements(); ee++) {
		const Element& elem = geometry.get_element(ee);
		if((signed)elem.get_dimension() == space_dim) {
			// store the med index
			reg2shape2idxs[elem.get_region().get_index_global()][elem.get_shape()].push_back(
				idx_global_2_idx_med[elem.get_index_global()]
			);
		}
	}
	// ----------------------------------------------
	// for every region
	// ----------------------------------------------
	vector<int> group_connectivity;
	vector<int> group_connectivity_index;
	for(unsigned int rr = 0; rr < geometry.get_num_regions(); rr++) {
		// create group and set basic information
		GROUP region;
		region.setName(geometry.get_region(rr).get_name());
		region.setMesh(meshing);
		region.setEntity(MED_CELL);

		// clean up some variables
		types.clear();
		reduced_counts.clear();
		group_connectivity.clear();
		group_connectivity_index.clear();
		// build med-style connectivity array
		int sidx = 1;
		for(unsigned int tt = 0; tt < reg2shape2idxs[rr].size(); tt++) {
			if(reg2shape2idxs[rr][tt].size() > 0) {
				// sort!
				sort(reg2shape2idxs[rr][tt].begin(), reg2shape2idxs[rr][tt].end());
				// store element type and number of elements of this type
				types.push_back(this->map_tdkp_to_med_elem_type(tt));
				reduced_counts.push_back(reg2shape2idxs[rr][tt].size());
				// store where this cell type begins (fortran ordering!)
				group_connectivity_index.push_back(sidx);
				// append elements to the connectivity array
				group_connectivity.insert(
					group_connectivity.end(),
					reg2shape2idxs[rr][tt].begin(),
					reg2shape2idxs[rr][tt].end()
				);
				sidx += reg2shape2idxs[rr][tt].size();
 			}
		}
		// and finish the connectivity index array
		group_connectivity_index.push_back(sidx);

		// set cell type information
		region.setNumberOfGeometricType(types.size());
		region.setGeometricType(&types[0]);
		// set connnectivity and the connectivity index
		region.setNumberOfElements(&reduced_counts[0]);
		region.setNumber(
			&group_connectivity_index[0],
			&group_connectivity[0]
		);

		meshing->addGroup(region);
		//cout << region;

	}
	//cout << *meshing;
	my_mesh = meshing;
	my_med->addMesh(meshing);

}

/** append element data to file */
template<class T>
void MEDDataIOPimpl::write_element(const ElementData<T>& data, const char* varname) {

	// get varname without result directory
	string short_varname = strip_directory(varname);

	// test
	TDKP_ASSERT(geometry.get_num_elements() == (unsigned)data.get_length() || (geometry.get_num_nodes() == (unsigned)data.get_length()), "MEDDataIO: the passed data does not fit to the geometry");

	// -------------------------------------------------
	// ensure that we don't overwrite a field!
	// -------------------------------------------------
	vector<string> field_names(my_med->getNumberOfFields(), "");
	my_med->getFieldNames(&field_names[0]);
	TDKP_ASSERT(find(field_names.begin(), field_names.end(), short_varname) == field_names.end(), "MEDDataIO: the output field '" << short_varname << "' already exists in the med data file");

	// -------------------------------------------------
	// setup med file:
	// - get support on all elements
	// - create field and set varnames
	// -------------------------------------------------
	SUPPORT* my_support     = new SUPPORT(my_mesh, "All elements", MED_CELL);
	FIELD<double>* my_field = new FIELD<double>(my_support, data.get_num_data_per_element());
	my_field->setName(short_varname.c_str());
	for(int ii = 0; ii < data.get_num_data_per_element(); ii++) {
		my_field->setComponentName(ii + 1, data.get_identifier(ii).c_str());
	}

	// -------------------------------------------------
	// loop over all available data
	// -------------------------------------------------
	vector<double> tmp_values(geometry.get_num_elements(), 0.0);
	for(int ii = 0; ii < data.get_num_data_per_element(); ii++) {
		// -------------------------------------------------
		// get values only on the vertices
		// -------------------------------------------------
		for(unsigned int ee = 0; ee < geometry.get_num_elements(); ee++) {
			tmp_values[this->idx_global_2_idx_med[geometry.get_element(ee).get_index_global()] - 1] = adjust_value(data.get_element_value(ee, ii));
		}
		my_field->setColumn(ii + 1, &tmp_values[0]);
	}
	my_med->addField(my_field);
	my_med->write(my_driver);
}
/** append node data to file */
template<class T>
void MEDDataIOPimpl::write_nodal(const NodeData<T>& data, const char* varname) {

	// test
	TDKP_ASSERT(geometry.get_num_vertices() == (unsigned)data.get_length() || (geometry.get_num_nodes() == (unsigned)data.get_length()), "MEDDataIO: the passed data does not fit to the geometry");

	// get varname without result directory
	string short_varname = strip_directory(varname);

	// -------------------------------------------------
	// ensure that we don't overwrite a field!
	// -------------------------------------------------
	vector<string> field_names(my_med->getNumberOfFields(), "");
	my_med->getFieldNames(&field_names[0]);
	TDKP_ASSERT(find(field_names.begin(), field_names.end(), short_varname) == field_names.end(), "MEDDataIO: the output field '" << short_varname << "' already exists in the med data file");

	// -------------------------------------------------
	// setup med file:
	// - get support (i.e. all nodes = all vertices, as we don't allow salome to create higher order elements)
	// - create field and set varnames
	// -------------------------------------------------
	SUPPORT* my_support     = new SUPPORT(my_mesh, "All nodes", MED_NODE);
	FIELD<double>* my_field = new FIELD<double>(my_support, data.get_num_data_per_node());
	my_field->setName(short_varname.c_str());
	for(int ii = 0; ii < data.get_num_data_per_node(); ii++) {
		my_field->setComponentName(ii + 1, data.get_identifier(ii).c_str());
	}

	// -------------------------------------------------
	// loop over all available data
	// -------------------------------------------------
	vector<double> tmp_values(geometry.get_num_vertices(), 0.0);
	for(int ii = 0; ii < data.get_num_data_per_node(); ii++) {
		// -------------------------------------------------
		// get values only on the vertices
		// -------------------------------------------------
		for(unsigned int nn = 0; nn < geometry.get_num_nodes(); nn++) {
			if(geometry.get_node(nn).get_index_vertex() != -1) {
				tmp_values[geometry.get_node(nn).get_index_vertex()] = adjust_value(data.get_node_value(nn, ii));
			}
		}
		my_field->setColumn(ii + 1, &tmp_values[0]);
	}
	my_med->addField(my_field);
	my_med->write(my_driver);
}

/** write med data based on a geometry object */
MEDDataIO::MEDDataIO(const Geometry& geometry, const char* outputfile) {
	impl = new MEDDataIOPimpl(geometry, outputfile);
}
/** finalizes med file (calls write) and deletes med object */
MEDDataIO::~MEDDataIO() {
	delete impl;
}

void MEDDataIO::write_element_cplx(const ElementData<cplx>& data, const char* varname) {
	impl->write_element(data, varname);
}
void MEDDataIO::write_nodal_cplx(const NodeData<cplx>& data, const char* varname) {
	impl->write_nodal(data, varname);
}
void MEDDataIO::write_element_dbl(const ElementData<double>& data, const char* varname) {
	impl->write_element(data, varname);
}
void MEDDataIO::write_nodal_dbl(const NodeData<double>& data, const char* varname) {
	impl->write_nodal(data, varname);
}


}
