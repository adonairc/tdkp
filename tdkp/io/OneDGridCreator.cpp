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

#include "tdkp/io/OneDGridCreator.h"

namespace tdkp {

/** 1D tdkp grid constructors */
OneDGridCreator::OneDGridCreator(const char* filename) {

	parse_region_data(filename);
	build_grid();
	
	// store unique identifer
	ostringstream ss;
	ss << "grid_" << filename;
	unique_identifier = ss.str();	
	
}

void OneDGridCreator::parse_region_data(const char* filename) {

	// -------------------------------------------
	// open file
	// -------------------------------------------
	ifstream fin(filename);
	if(!fin) {
		TDKP_GENERAL_EXCEPTION("can not read given gridfile " << filename);	
	}
	
	// ------------------------------------------
	// read first line and check if the header is o.k.
	// ------------------------------------------ 
	const string expected_header("TDKP1DGRID");
	string sin;
	getline(fin,sin);
	TDKP_ASSERT(sin.size() >= expected_header.size(), "first line in grid file starts with " << sin << ", i expected " << expected_header); 	 
	bool good = true;	
	for(unsigned int ii = 0; ii < expected_header.size(); ii++) {
		if(expected_header[ii] != sin[ii]) {
			good = false;	
		}	
	}
	TDKP_ASSERT(good, "first line in grid file starts with " << sin << ", i expected " << expected_header);

	// ------------------------------------------------
	// parse specifications
	// ------------------------------------------------
	while(getline(fin,sin)) {
		if(sin.size() > 0) {
			// check number of entries (must be four)
			int  entry_idx = -1;
			bool entry     = false;
			vector<string> entries;
			for(unsigned int ii = 0; ii < sin.size(); ii++) {
				if(sin[ii] != ' ' && sin[ii] != '\t') {
					// new expression
					if(!entry) {
						entry = true;
						entry_idx++;
						entries.push_back(string(""));	
					}
					// add character
					entries[entry_idx].push_back(sin[ii]);	
				} else {
					entry = false;	
				}
			}	
			// check number entries in file (must be four per line)
			TDKP_ASSERT(entries.size() == 4, "can not parse line " << sin << ". expected four entries, got " << entries.size());
			// collect strings name and material
			RegionData regdata;
			regdata.name          = entries[0];
			regdata.material_name = entries[1];
			// parse width via istringstreams
			istringstream sparse(entries[2]); 
			sparse >> regdata.width;
			TDKP_ASSERT(regdata.width > 0.0, "width of region " << entries[0] << " is " << regdata.width << " which i consider to be bullshit");
			istringstream sparse2(entries[3]);			
			sparse2 >> regdata.element_width;
			TDKP_ASSERT(regdata.width > 0.0, "element width of region " << entries[0] << " is " << regdata.element_width << " which i consider to be bullshit");
			// add region
			region_data.push_back(regdata);
		}	
	}	
	TDKP_ASSERT(region_data.size() > 0, "no regions found in file " << filename);
	
}

void OneDGridCreator::build_grid() {
	
	// -----------------------------------------
	// build grid from region data
	// -----------------------------------------
	double startx = 0.0;
	double endx   = 0.0;
	double avg    = startx;
	// -----------------------------------------
	// create left vertex
	// -----------------------------------------
	vertices.push_back(startx);
	ostringstream sout;
	for(unsigned int ii = 0; ii < region_data.size(); ii++) {		
		// -----------------------------------------------
		// end point of region
		// -----------------------------------------------
		endx = startx + region_data[ii].width;		
		// -----------------------------------------------
		// calculate real number of elements
		// -----------------------------------------------	
		unsigned int num_elements = static_cast<unsigned int>(ceil(region_data[ii].width / region_data[ii].element_width));
		TDKP_ASSERT(num_elements < 100000, "num_elements < 100000 expected, but found " << num_elements << " which is highly doubious");
		// -----------------------------------------------
		// create vertices and elements
		// -----------------------------------------------
		unsigned int num_created = 0;
		for(unsigned int ee = 1; ee <= num_elements; ee++) {
			double x = static_cast<double>(ee) * (endx - startx) / static_cast<double>(num_elements) + startx;
			avg += x;
			vertices.push_back(x);
			element_region_idx.push_back(ii);
			num_created++; 	
		} 	
		// -----------------------------------------------
		// assemble for output
		// -----------------------------------------------
		sout << "\nregion " << region_data[ii].name << ", idx = " << ii << ", material = " 
		     << region_data[ii].material_name << ", width = " << region_data[ii].width << ", nelem = "
		     << num_created;
		     			
		// ----------------------------------------------
		// shift startx
		// ----------------------------------------------
		startx = endx;	
	}
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "OneDGridCreator: created the following regions: " << sout.str());	
	
	// ---------------------------------------------
	// center structure around 0
	// ---------------------------------------------
	avg = avg / static_cast<double>(vertices.size());
	for(unsigned int ii = 0; ii < vertices.size(); ii++) {
		vertices[ii] -= avg;	
	}



}

OneDGridCreator::~OneDGridCreator() {} 	
/** return number of vertices in geometry */
unsigned int OneDGridCreator::get_num_vertices() const {
	return vertices.size(); 	
}
/** return number of elements in geometry */
unsigned int OneDGridCreator::get_num_elements() const {
	return element_region_idx.size();	
}
/** return number of elements per region in geometry */
unsigned int OneDGridCreator::get_num_elements_in_region(unsigned int region_index_global) const {
	TDKP_ASSERT(region_index_global < region_data.size(), "");
	unsigned int num = 0; 
	for(unsigned int ii = 0; ii < element_region_idx.size(); ii++) {
		if(element_region_idx[ii] == (signed)region_index_global) {
			num++;	
		}
	} 
	return num;		
}
/** return number of regions in geometry */
unsigned int OneDGridCreator::get_num_regions()  const {
	return region_data.size();	
}
/** return region name */
string OneDGridCreator::get_region_name(unsigned int region_index_global) const {
	TDKP_ASSERT(region_index_global < region_data.size(), "");
	return region_data[region_index_global].name;
}
/** return material name of region */
string OneDGridCreator::get_region_material_name(unsigned int region_index_global) const {
	TDKP_ASSERT(region_index_global < region_data.size(), "");
	return region_data[region_index_global].material_name;
	
}
/** return global element index of element in region */
unsigned int OneDGridCreator::get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const {
	TDKP_ASSERT(region_index_global < region_data.size(), "");
	unsigned int index_region = 0; 
	// scan trough all elements
	for(unsigned int ii = 0; ii < element_region_idx.size(); ii++) {
		// if element is in region
		if(element_region_idx[ii] == (signed)region_index_global) {
			// if index region matches
			if(index_region == element_index_region) {
				return ii;	
			}
			index_region++;	
		}
	} 
	TDKP_GENERAL_EXCEPTION("no element with index region " << element_index_region << " found in region " << region_index_global); 		
}
/** return number of vertices in element */
unsigned int OneDGridCreator::get_num_vertices_in_element(unsigned int element_index) const {
	return 2;	
}
/** return global vertex index of local element vertex */
unsigned int OneDGridCreator::get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const {
	if(local_vertex_index == 0) {
		return element_index;	
	} else if(local_vertex_index == 1) {
		return element_index + 1;	
	} else {
		TDKP_GENERAL_EXCEPTION("local_vertex_index " << local_vertex_index << " requested which is not existing in a 1D line segment with 2 vertices!");	
	}
}
/** return vertex coordinate x */
double OneDGridCreator::get_vertex_x_coord(unsigned int vertex_index) const {
	TDKP_ASSERT(vertices.size() > vertex_index, "");
	return vertices[vertex_index];	
}
/** return unique identifier for grid file */
string OneDGridCreator::get_unique_identifier() const {
	return unique_identifier;
}


} /* end of namespace */
