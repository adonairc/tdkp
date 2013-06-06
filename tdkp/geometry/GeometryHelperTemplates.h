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

#ifndef GEOMETRYHELPERTEMPLATES_H_
#define GEOMETRYHELPERTEMPLATES_H_


#include "tdkp/geometry/Geometry.h"
#include "tdkp/common/DataTypes.h"


namespace tdkp {

class GeometryHelperTemplates {
public:

template<class T> 
static void map_index_global_to_index_internal(const Geometry& geom, const NodeData<T>& in, NodeData<T>& out);

template<class T> 
static void map_index_internal_to_index_global(const Geometry& geom, const NodeData<T>& in, NodeData<T>& out);

template<class T>
static void map_index_global_to_index_vertex(const Geometry& geom, const NodeData<T>& in, vector<T>& out);

};


/** map values from global to local indices */
template<class T> 	
void GeometryHelperTemplates::map_index_global_to_index_internal(const Geometry& gref, const NodeData<T>& in, NodeData<T>& out) {
	
	
	const int num_equations = in.get_num_data_per_node();
	int index_internal, index_global;
	
	// just assure that we do not connode something that does not match ...
	TDKP_ASSERT((unsigned int)in.get_length() == gref.get_num_nodes(), "NodeData has not the correct length!");
	// reinit StdNodeData
	out.set_length(gref.get_num_nonzero_nodes(), in.get_num_data_per_node());
	
	// map node data
#pragma omp parallel for default(none) private(index_internal, index_global)	
	for(int ii = 0; ii < (signed)gref.get_num_nodes(); ii++) {
		const Node& node = gref.get_node(ii);
		index_internal = node.get_index_internal();
		if(index_internal != -1) {
			index_global   = node.get_index_global();
			for(int ee = 0; ee < num_equations; ee++) {
				out.set_node_value(index_internal,ee,in.get_node_value(index_global,ee));
			}
		}    	
	}  	
}

template<class T> 	
void GeometryHelperTemplates::map_index_internal_to_index_global(const Geometry& gref, const NodeData<T>& in, NodeData<T>& out) {
		
	const unsigned int num_equations = in.get_num_data_per_node();
	int index_internal, index_global;
		
		// just assure that we do not connode something that does not match ...
	TDKP_ASSERT(in.get_length() == (signed)gref.get_num_nonzero_nodes(), "NodeData has not the correct length!");
	// reinit StdNodeData
	out.set_length(gref.get_num_nodes(), in.get_num_data_per_node());
	
	// map node data
#pragma omp parallel for default(none) private(index_internal, index_global)  	
	for(int ii = 0; ii < (signed)gref.get_num_nodes(); ii++) {
		const Node& node = gref.get_node(ii);
		index_internal = node.get_index_internal();
		index_global   = node.get_index_global();
		if(index_internal == -1) {
			for(unsigned int ee = 0; ee < num_equations; ee++) {
				out.set_node_value(index_global,ee, (T)0.0);
			}				
		} else {			
			for(unsigned int ee = 0; ee < num_equations; ee++) {
				out.set_node_value(index_global,ee,in.get_node_value(index_internal,ee));
			}
		}    	
	}
}

template<class T> 	
void GeometryHelperTemplates::map_index_global_to_index_vertex(const Geometry& geometry, const NodeData<T>& in, vector<T>& out) {

	// -----------------------------------------------
	// calculate num vertices
	// -----------------------------------------------
	unsigned int nvertex = 0;
	for(unsigned int ii = 0; ii < geometry.get_num_nodes(); ii++) {
		if(geometry.get_node(ii).get_index_vertex() != -1) {
			nvertex++;	
		}	
	}	
	out.resize(nvertex * in.get_num_data_per_node());
	
	// -----------------------------------------------
	// map data
	// -----------------------------------------------
	for(unsigned int ii = 0; ii < geometry.get_num_nodes(); ii++) {
		if(geometry.get_node(ii).get_index_vertex() != -1) {
			for(int jj = 0; jj < in.get_num_data_per_node(); jj++) {
				out[geometry.get_node(ii).get_index_vertex() * in.get_num_data_per_node() + jj] = in.get_node_value(ii,jj);
			}		
		}
	} 
		
}


} // end of namespace 

#endif /*GEOMETRYHELPERTEMPLATES_H_*/
