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

#ifndef _NODE3D_H
#define _NODE3D_H
#include "tdkp/common/all.h"

namespace tdkp {

	/** node location */	
	typedef char NodeLocation;
	const NodeLocation location_undefined = 'u';
	// interior and exterior are exclusive
	const NodeLocation location_interior  = 'i';
	const NodeLocation location_exterior  = 'e';	
		
	/** node class for 1D,2D and 3D problems */
	class Node {	
	public:
		enum NodeValueType {
			FunctionValue, Derivative_X, Derivative_Y, Derivative_Z
		};
		
		Node(const Node& copy);
		Node(double x, double y, double z);
		Node(unsigned int index_global, double x, double y, double z);
		Node(double x, double y);
		Node(unsigned int index_global, double x, double y);
		Node(double x);
		Node(unsigned int index_global, double x);		
		~Node();
				
		int  get_index_internal() const { return index_internal; }
		void set_index_internal(int index_internal);
		
		int  get_index_global() const { return index_global; }
		void set_index_global(unsigned int index_global);
		
		int  get_index_vertex() const { return index_vertex; }
		void set_index_vertex(int index_vertex);
		
		void           set_location(NodeLocation loc);
		NodeLocation   get_location() const;

		NodeValueType  get_value_type() const { return node_value_type; }
		void           set_value_type(NodeValueType node_value_type_);
									
		const double*  get_coords() const { return this->coord; }
		double&		   get_coord(unsigned short ii) { return this->coord[ii]; }
		const double&  get_coord(unsigned short ii) const { return this->coord[ii]; }
		short          get_dimension() const { return this->dimension; }
		
		static double  distance(const Node &a, const Node &b);
					
		void print() const;
		
		// ------------------------------------
		// vertex contribution
		// allows automatic (linear) interpolation of vertex values 
		// to nodal values
		// information is set by elements during grid building
		// if node is a non vertex. if it is a vertex, its 
		// clear a-priori ...
		// this does not work for hermite elements ...
		// ------------------------------------
		unsigned int get_num_contributions() const;
		void get_contribution(unsigned int ii, int& index_vertex, double& contrib) const;
		void set_contribution(const int index_vertex, const double& contrib);					
																
	private:
		/** coordinates */
		double coord[3];
		/** internal index in the global equation matrix. if set to -1, 
		 * node has dirichlet bc (= 0) and is therefore discarded */
		int index_internal;
		/** global index of node in geometry */
		unsigned int index_global;		
		/** index of node in grid file, so if it is a vertex, then index != -1 */
		int index_vertex;
		/** dimension */
		short dimension;
		/** location */
		NodeLocation location; 
		/** node value type */
		NodeValueType node_value_type;
		
		struct Contribution {
			int index_vertex;
			double contrib;
		};
		/** linear interpolation coefficents from other vertices */
		vector<Contribution> vertex_contributions;
		
				 
	};	
	
	ostream& operator<<(ostream& stream, const Node & node);
	
}


#endif
