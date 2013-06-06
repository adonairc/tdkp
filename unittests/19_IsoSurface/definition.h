#ifndef DEFINITION_H_
#define DEFINITION_H_

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
#include <math.h>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <boost/lexical_cast.hpp>
#include "tdkp/povray/IsoSurfaceGenerator.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/io/InputParser.h"
#include "tdkp/io/AsciiGridReader.h"

using namespace tdkp;

namespace tdkp {


}

class IsoSurfaceTest : public CxxTest::TestSuite {

public:			
	IsoSurfaceTest() {
		fout.open("19_IsoSurface_output.log", ios::app | ios::out);		
		Logger::get_instance()->add_listener(&fout);
		//Logger::get_instance()->del_listener(&std::cout);
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);		 
	}
	
	~IsoSurfaceTest() {			
		Logger::get_instance()->del_listener(&fout);
		fout.close();		
	}
	std::fstream fout;
	StdNodeData<double>* build_spherical_data(const Geometry& geometry);
			
	void setUp() {}				
	void tearDown() {}		
	
	void test_isosurface();
};

void IsoSurfaceTest::test_isosurface() {

	InputParser parser;
	Geometry* geometry           = parser.read_geometry(AsciiGridReader("19_IsoSurface/dot_24105_verts.asc.gz"));
	StdNodeData<double>* data  = build_spherical_data(*geometry);
	//parser.write_to_dfise_file(*geometry, *data, "test.dat");
	IsoSurfaceGenerator iso;	
	// nice: ;-)
	iso.set_geometry(geometry);
	iso.set_data(data);	
	iso.create_lists();
	iso.create_surface_plot("test75.pov", 0.5);
	

}

StdNodeData<double>* IsoSurfaceTest::build_spherical_data(const Geometry& geometry) {
	StdNodeData<double>* data = new StdNodeData<double>(1, geometry.get_num_nodes());
	Node center(0.0, 0.0, 0.0);
	double value;
	for(unsigned int ii = 0; ii < geometry.get_num_nodes(); ii++) {
		const Node& vert = geometry.get_node(ii);
		value = Node::distance(vert, center);
		data->set_node_value(ii,0,value); 			
	}  
	return data;	
}

#endif /*DEFINITION_H_*/
