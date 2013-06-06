
#define TDKP_UNITTEST

bool be_chatty = true;

// ------------------------------------------
// standard includes
// ------------------------------------------
#include <cxxtest/TestSuite.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
// ------------------------------------------
// standard tdkp includes and logger def
// ------------------------------------------
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"

// ------------------------------------------
// class includes
// ------------------------------------------
#include <vector>
#include <complex>
#include <map>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/geometry/Region.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/main/Fields.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/geometry/Element3DTR.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/io/InputParser.h"
#include "tdkp/io/AsciiGridReader.h"
#include "tdkp/geometry/GeometryHelperTemplates.h"
#include "tdkp/geometry/BoundaryCondition.h"


using namespace tdkp;
using namespace std;

//   cube_700_vertices.grd  cube_no_interior.grd  cube_only_interior.grd  default.grd

struct GridTestFile {
	int nvert, nelem, nregion;
	int dimension;
	string filename;
	double bounding_box[2][3];			
};

class GeometryTest : public CxxTest::TestSuite {
public:	
	
	std::fstream fout;		
	vector<GridTestFile> gridtests;	
		
	GeometryTest();
	~GeometryTest() { fout.close(); }		
	
	void test_mesh_with_dat_mf();
	void test_boundary_conditions_plane();	
	void test_input_parser_geometry();
	void test_helper_templates();

	
};


GeometryTest::GeometryTest() {
	
	fout.open("03_Geometry_output.log", ios::app | ios::out);
	Logger::get_instance()->add_listener(&fout);
	Logger::get_instance()->set_level(LOG_INFO_DEVEL2);
	if(!be_chatty) {	
		Logger::get_instance()->del_listener(&std::cout);
	}

	// -------------------------------------------------------		
	// setup grid tests
	// -------------------------------------------------------		
	GridTestFile cube67;
	cube67.nvert     = 68;
	cube67.nelem     = 234;
	cube67.nregion   = 1;
	cube67.dimension = 3;
	cube67.filename  = "03_Geometry/valid_grids/cube_67_points.asc.gz";
	gridtests.push_back(cube67);
	GridTestFile cube700;
	cube700.nvert     = 670;
	cube700.nelem     = 3655;
	cube700.nregion   = 1;
	cube700.dimension = 3;
	cube700.filename  = "03_Geometry/valid_grids/cube_700_vertices.asc.gz";
	gridtests.push_back(cube700);
	GridTestFile  dot;
	dot.nvert     = 1955;
	dot.nelem     = 12567;
	dot.nregion   = 2;
	dot.dimension = 3;
	dot.filename  = "03_Geometry/valid_grids/dot1955.asc.gz";
	gridtests.push_back(dot);
	GridTestFile wire;
	wire.nvert     = 596;
	wire.nelem     = 689;
	wire.nregion   = 2;
	wire.dimension = 2;
	wire.filename  = "03_Geometry/valid_grids/wire.asc.gz";
	gridtests.push_back(wire);
	GridTestFile well;
	well.nvert     = 596;
	well.nelem     = 689;
	well.nregion   = 2;
	well.dimension = 2;
	well.filename  = "03_Geometry/valid_grids/wire.asc.gz";
	gridtests.push_back(well);
	
	
 
}

void GeometryTest::test_input_parser_geometry() {
									
	InputParser parser;
	Geometry* geo;
		
	for(vector<GridTestFile>::const_iterator it = gridtests.begin(); it != gridtests.end(); it++) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "checking gridfile " << (*it).filename);
		try { 
			geo = parser.read_geometry(AsciiGridReader((*it).filename.c_str()));
									
		} catch(Exception* e) {
			std::cout << (e->get_reason()) << std::endl;
			return;	
		}
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "doing verify");
		TS_ASSERT(geo->verify());
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "geometry verified");
		TS_ASSERT_EQUALS((signed)geo->get_num_nodes(),(*it).nvert);
		TS_ASSERT_EQUALS((signed)geo->get_num_elements(),(*it).nelem);
		TS_ASSERT_EQUALS((signed)geo->get_num_regions(), (*it).nregion);
		TS_ASSERT_EQUALS((signed)geo->get_dimension(), (*it).dimension);	
		delete geo;			
	}	
	cerr << "test is over\n";	
}


void GeometryTest::test_boundary_conditions_plane() {
	try {
		InputParser parser;
		MaterialDatabase matdb;
		MaterialPropertyContainer* material = MaterialPropertyContainer::create_zincblende_material();
		matdb.add_material("GaAs", material);
		Geometry* geometry = parser.read_geometry(AsciiGridReader(("03_Geometry/valid_grids/cube_700_vertices.asc.gz")));
		geometry->set_materials(matdb);		
		BCDirichlet3DPlanes* pbc = new BCDirichlet3DPlanes(*geometry);
		BCDirichlet3DPlanes& bc = *pbc;
				
		// -------------------------------------------
		// set two planes
		// -------------------------------------------
		// first plane: x/y at bottom
		Vector3D a0(1.0, 0.0, 0.0);
		Vector3D b0(0.0, 1.0, 0.0);
		Vector3D r0(0.0, 0.0, -2.0);
		bc.add_dirichlet_plane(r0, a0, b0);
		// second plane: y/z at left side
		Vector3D a1(0.0, 1.0, 0.0);
		Vector3D b1(0.0, 0.0, 1.0);
		Vector3D r1(-2.0, 0.0, 0.0);
		bc.add_dirichlet_plane(r1, a1, b1);		
		bc.prepare();
		geometry->set_boundary_conditions(&bc);
		// test
		for(unsigned int vv = 0; vv < geometry->get_num_nodes(); vv++) {
			const Node& vert = geometry->get_node(vv);
			bool dv = false;
			if(tdkp_math::abs(vert.get_coord(2) + 2.0) < 1.0e-6) {
				dv = true;
			} 
			if(tdkp_math::abs(vert.get_coord(0) + 2.0) < 1.0e-6) {
				dv = true;	
			}	
			TS_ASSERT(dv == (vert.get_index_internal() == -1));
		}
		delete geometry;
			
	} catch(Exception  *e) {
		TS_FAIL(e->get_reason());	
	}	
}

void GeometryTest::test_mesh_with_dat_mf() {

	try {
		InputParser parser;
		Geometry* geometry = parser.read_geometry(AsciiGridReader("03_Geometry/valid_grids/quasi1D1SQW_msh.asc.gz"));
		TS_ASSERT(geometry->verify());
		delete geometry;				
	} catch(Exception  *e) {
		TS_FAIL(e->get_reason());	
	}
	exit(0);		
}

void GeometryTest::test_helper_templates() {
	
	InputParser parser;
	Geometry* geo;
	// loop over different grids
	  for(vector<GridTestFile>::const_iterator it = gridtests.begin(); it != gridtests.end(); it++) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "checking gridfile: " << (*it).filename);			
		try { 
			geo = parser.read_geometry(AsciiGridReader((*it).filename.c_str()));
		} catch(Exception* e) {
			std::cout << (e->get_reason()) << std::endl;
			return;	
		}
		// loop over different num equations
		for(int ee = 1; ee < 5; ee++) {
				
			StdNodeData<double> ref_internal(ee, geo->get_num_nonzero_nodes());
			StdNodeData<double> ref_external(ee, geo->get_num_nodes());
			// initialize with values
			for(int qq = 0; qq < ee; qq++) {
				for(int ii = 0; ii < (signed)geo->get_num_nonzero_nodes(); ii++) {
					ref_internal.set_node_value(ii,qq, double(ii * ee + qq)); 							
				}					
				for(int ii = 0; ii < (signed)geo->get_num_nodes(); ii++) {
					// zero nodes must be zero ;-)
					if(geo->get_node(ii).get_index_internal() != -1) {
						ref_external.set_node_value(ii,qq, double(ii * ee + qq));
					} else {
						ref_external.set_node_value(ii,qq, 0.0);
					} 							
				} 	
			}
			
			StdNodeData<double> to, back;
			// test map index global to index internal and back
			try {
				GeometryHelperTemplates::map_index_global_to_index_internal(*geo, ref_external, to);
			} catch(Exception*e) {
				TS_FAIL(e->get_reason());
			}				
			TS_ASSERT(to.get_num_data_per_node() == ref_internal.get_num_data_per_node());
			TS_ASSERT(to.get_length() == ref_internal.get_length());
			try {				
				GeometryHelperTemplates::map_index_internal_to_index_global(*geo, to, back);
			} catch(Exception*e) {
				TS_FAIL(e->get_reason());
			}
			TS_ASSERT(back.get_num_data_per_node() == ref_external.get_num_data_per_node());
			TS_ASSERT(back.get_length() == ref_external.get_length());				
			// now back must be equal with ref external				
			for(int ii = 0; ii < (signed)geo->get_num_nodes(); ii++) {
				for(int qq = 0; qq < ee; qq++) {
					TS_ASSERT_EQUALS(back.get_node_value(ii,qq), ref_external.get_node_value(ii,qq));
				}	
			}
			// other way around
			try {
				GeometryHelperTemplates::map_index_internal_to_index_global(*geo, ref_internal, to);
			} catch(Exception*e) {
				TS_FAIL(e->get_reason());
			}				
			TS_ASSERT(to.get_num_data_per_node() == ref_external.get_num_data_per_node());
			TS_ASSERT(to.get_length() == ref_external.get_length());
			try {												
				GeometryHelperTemplates::map_index_global_to_index_internal(*geo, to, back);
			} catch(Exception*e) {
				TS_FAIL(e->get_reason());
			}				
			TS_ASSERT(back.get_num_data_per_node() == ref_internal.get_num_data_per_node());
			TS_ASSERT(back.get_length() == ref_internal.get_length());
			
			// now back must be equal with ref external				
			for(int ii = 0; ii < (signed)geo->get_num_nonzero_nodes(); ii++) {
				for(int qq = 0; qq < ee; qq++) {
					TS_ASSERT_EQUALS(back.get_node_value(ii,qq), ref_internal.get_node_value(ii,qq));
				}	
			}				
			
		}				
		delete geo;			
	}		
	
}


