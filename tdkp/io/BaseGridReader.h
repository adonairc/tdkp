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

#ifndef BASEGRIDREADER_H_
#define BASEGRIDREADER_H_

#include "tdkp/geometry/Element.h"

namespace tdkp {

class BaseGridReader {
public:
	
	virtual ~BaseGridReader() {}	
	/** return dimension of grid */
	virtual unsigned int get_dimension() const = 0;
	/** return number of vertices in geometry */
	virtual unsigned int get_num_vertices() const = 0;
	/** return number of elements in geometry */
	virtual unsigned int get_num_elements() const = 0;
	/** return number of elements per region in geometry */
	virtual unsigned int get_num_elements_in_region(unsigned int region_index_global) const = 0;
	/** return number of regions in geometry */
	virtual unsigned int get_num_regions()  const = 0;
	/** return region name */
	virtual string get_region_name(unsigned int region_index_global) const = 0;
	/** return material name of region */
	virtual string get_region_material_name(unsigned int region_index_global) const = 0;
	/** return global element index of element in region */
	virtual unsigned int get_element_index_global(unsigned int region_index_global, unsigned element_index_region) const = 0;
	/** return element shape */
	virtual Element::ElementShape get_element_shape(unsigned int element_index) const = 0;
	/** return number of vertices in element */
	virtual unsigned int get_num_vertices_in_element(unsigned int element_index) const = 0;
	/** return global vertex index of local element vertex */
	virtual unsigned int get_vertex_in_element(unsigned int element_index, unsigned int local_vertex_index) const = 0;
	/** return vertex coordinate x */
	virtual double get_vertex_x_coord(unsigned int vertex_index) const = 0;
	/** return vertex coordinate x */
	virtual double get_vertex_y_coord(unsigned int vertex_index) const = 0;
	/** return vertex coordinate x */
	virtual double get_vertex_z_coord(unsigned int vertex_index) const = 0;
	/** return unique identifier for grid file */
	virtual string get_unique_identifier() const = 0;	
};

}

#endif /*BASEGRIDREADER_H_*/
