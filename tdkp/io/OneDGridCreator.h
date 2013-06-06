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

#ifndef ONEDGRIDCREATOR_H_
#define ONEDGRIDCREATOR_H_

#include "tdkp/io/BaseGridReader.h"

namespace tdkp
{

class OneDGridCreator : public BaseGridReader
{
public:
	OneDGridCreator(const char* filename);
	virtual ~OneDGridCreator();	
	/** return dimension of grid */
	virtual unsigned int get_dimension() const { return 1; };
	/** return number of vertices in geometry */
	virtual unsigned int get_num_vertices() const;
	/** return number of elements in geometry */
	virtual unsigned int get_num_elements() const;
	/** return number of elements per region in geometry */
	virtual unsigned int get_num_elements_in_region(unsigned int region_index_global) const;
	/** return number of regions in geometry */
	virtual unsigned int get_num_regions()  const;
	/** return region name */
	virtual string get_region_name(unsigned int region_index_global) const;
	/** return material name of region */
	virtual string get_region_material_name(unsigned int region_index_global) const;
	/** return global element index of element in region */
	virtual unsigned int get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const;
	/** return element shape */
	virtual Element::ElementShape get_element_shape(unsigned int element_index) const { return Element::line; }
	/** return number of vertices in element */
	virtual unsigned int get_num_vertices_in_element(unsigned int element_index) const;
	/** return global vertex index of local element vertex */
	virtual unsigned int get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const;
	/** return vertex coordinate x */
	virtual double get_vertex_x_coord(unsigned int vertex_index) const;
	/** return vertex coordinate x */
	virtual double get_vertex_y_coord(unsigned int vertex_index) const { return 0.0; };
	/** return vertex coordinate x */
	virtual double get_vertex_z_coord(unsigned int vertex_index) const { return 0.0; };
	/** return unique identifier for grid file */
	virtual string get_unique_identifier() const;	

private:
	void parse_region_data(const char* filename);
	void build_grid();

	string unique_identifier;

	struct RegionData {
		string name;
		string material_name;
		double width;
		double element_width;
	};
	vector<RegionData> region_data;

	vector<int>    element_region_idx;
	vector<double> vertices;	

	
	
	
	
};

}

#endif /*ONEDGRIDCREATOR_H_*/
