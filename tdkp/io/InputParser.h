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

#ifndef _INPUTPARSER_H
#define _INPUTPARSER_H


#include <fstream>
#include <string>
#include "tdkp/common/all.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/main/Fields.h"
#include "tdkp/io/BaseGridReader.h"

namespace tdkp {
	
using namespace std;



/** altough the name suggests something different, this is the base class for reading and writing input and output files,
 *  depending on actual objects for parsing grid files.
 */				
class InputParser{
public:	
	InputParser() {}	
	Geometry*    read_geometry(const BaseGridReader& grid_reader, unsigned int polynom_degree = 1);
	~InputParser() {};
		
	template<class T> 
	void write_ascii(const XYData<T>& data, const char* filename, int precision = 12) const throw(Exception*);
	
	template<class T>
	void write_ascii(const Geometry& geom, const ElementData<T>& data, const char* filename) const;

	template<class T>
	void write_ascii(const Geometry& geom, const NodeData<T>& data, const char* filename) const;
	 
	template<class T>
	void write_ascii_real(const XYData<T>& data, const char* filename, int precision = 12) const throw(Exception*);
		
	template<class T> 
	void read_binary(NodeData<T>&  data, const char* filename) const throw(Exception*);
	template<class T> 		
	void write_binary(const NodeData<T>& data, const char* filename) const throw(Exception*);		
	template<class T> 
	void read_binary(ElementData<T>&  data, const char* filename) const throw(Exception*);
	template<class T> 		
	void write_binary(const ElementData<T>& data, const char* filename) const throw(Exception*);			
	template<class T> 
	void read_binary(NodeData<T>&  data, istream& fin) const throw(Exception*);
	template<class T> 		
	void write_binary(const NodeData<T>& data, ostream& fout) const throw(Exception*);		
	template<class T> 
	void read_binary(ElementData<T>&  data, istream& fin) const throw(Exception*);
	template<class T> 		
	void write_binary(const ElementData<T>& data, ostream& fout) const throw(Exception*);	

			
protected:
		
	/** binary data file header and footer definition, attention functions reading
	 *  have only a buffer of length 50 to read headers/footer. adjust if you change definitions here */
	static const char* NodeDataHeader;
	static const char* NodeDataFooter;
	static const char* ElementDataHeader;
	static const char* ElementDataFooter;

	void write_ascii_value(ofstream& stream, const double& value) const {
		stream << value;
	}

	void write_ascii_value(ofstream& stream, const cplx& value) const {
		stream << value.real() << " \t" << value.imag();
	}
			
private:
	
	/** base class for temporary building of edges / faces */
	class VertexObject {
	public:
		VertexObject();			
		~VertexObject();
		VertexObject(const vector<unsigned int>& vertex_indices_);
		VertexObject* clone() const;
		void init(const vector<unsigned int>& vertex_indices_);
		void add_element(unsigned int element_index_global);
		
		unsigned int get_num_elements() const;
		unsigned int get_element_index(unsigned int ii) const;
		unsigned int get_num_vertices() const;
		unsigned int get_vertex_index(unsigned int ii) const;
		unsigned int get_lowest_vertex_index() const { return lowest_vertex_index; }
		int get_index_global() const { return index_global; }
		void set_index_global(int index_global_) { index_global = index_global_; }
		/** check if all vertices match */
		bool compare_vertices(const VertexObject& rhs) const;
		/** check if rhs vertices are contained in object */
		bool all_rhs_vertices_match(const VertexObject& rhs) const;
		// linked list functions
		VertexObject* get_next() const { return next; }		
		VertexObject* locate(const VertexObject& obj);
		void insert(VertexObject* obj);
		
		// node creation functions
		void  set_location(char loc) { location = loc; }
		char  get_location() const { return location; }
		Node* locate_node(Geometry& geometry, int& current_node_idx, const vector<double>& coords, unsigned int tag, const double& tolerance);		
	private:
		VertexObject(const VertexObject& copy);
		// hard coded, face has maximum 4 vertices.
		unsigned int         vertex_indices[4];
		unsigned int         nverts;
		vector<unsigned int> element_indices;
		unsigned int lowest_vertex_index;
		VertexObject* next;
		/** vobj index: index of edge, face or vertex that we stand for */ 
		int index_global;
		/** nodes vector for building and locating additional nodes */
		vector<Node*>        nodes;
		vector<unsigned int> node_tags;				
		char location;
	};

	static void build_vertex_objects_tree(unsigned int element_index, vector<VertexObject*>& objects, VertexObject& vobj);
	static unsigned int set_global_indices_to_tree(vector<VertexObject*>& objects);
			
};	
			
} // end namespace tdkp

#include "tdkp/io/InputParser.tcc"

#endif


