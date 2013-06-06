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

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <vector>
#include "tdkp/common/Vector3D.h"
#include "tdkp/geometry/Region.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/geometry/ElementBoundary.h"
#include "tdkp/main/MaterialDatabase.h"


namespace tdkp {

	class BoundaryCondition;
					
	class Geometry{
	public:
				
		static const int max_nodes_per_element = 10;
	
		typedef std::vector<Element*>::iterator 	   element_iterator;
		typedef std::vector<Element*>::const_iterator  element_const_iterator;
		typedef std::vector<Node*>::iterator    	   node_iterator;
		typedef std::vector<Node*>::const_iterator     node_const_iterator;	
		typedef std::vector<Region*>::iterator         region_iterator;	
		typedef std::vector<Region*>::const_iterator   region_const_iterator;
		typedef std::vector<ElementBoundary*>::const_iterator boundary_const_iterator;
		typedef std::vector<ElementBoundary*>::iterator boundary_iterator;		
			
		Geometry();
		Geometry(unsigned int dimension, unsigned int nelem, unsigned int nnode, unsigned int nreg);
		~Geometry();
		
		void         init(unsigned int dimension, unsigned int nelem, unsigned int nnode, unsigned int nreg);
		bool         verify() const;
		void         prepare();
		void         prepare_boundaries();		
		void         set_internal_indices();
		void         rescale_node_coordinates(const double& rescale);	
				
		void         set_boundary_conditions(BoundaryCondition* bc_handler_);
		const BoundaryCondition& get_boundary_conditions() const;
		void         set_num_edges(int nedges);
		void         set_num_faces(int faces);	
									
		unsigned int add_node(Node*);
		unsigned int add_element(Element*);
		unsigned int add_region(Region*);
		unsigned int add_element_boundary(ElementBoundary*);
		
		const Node&    get_node (unsigned int node_idx) const;
		const Element& get_element(unsigned int elem_idx) const; 
		const Region&  get_region (unsigned int region_idx) const;
		const ElementBoundary& get_element_boundary(unsigned int boundary_idx) const;
		
		Node&          get_node (unsigned int node_idx);
		Element&       get_element(unsigned int elem_idx); 
		Region&        get_region (unsigned int region_idx);
		ElementBoundary& get_element_boundary(unsigned int boundary_idx);
					
		unsigned int get_num_nonzero_nodes() const;
		unsigned int get_num_elements() const { return elements.size(); }
		unsigned int get_num_nodes() const { return nodes.size(); }
		unsigned int get_num_vertices() const;
		unsigned int get_num_regions() const  { return regions.size(); }
		unsigned int get_num_boundaries() const { return element_boundaries.size(); }				
		
		unsigned int get_num_faces() const { return nfaces; }
		unsigned int get_num_edges() const { return nedges; }		
		unsigned int get_dimension() const { return this->dimension; }		
		int          find_next_node(double x, double y = 0.0, double z = 0.0) const;		
		void 		 list_elements(int node_index_global) const;
		
		element_iterator       elements_begin() { return elements.begin(); }
		element_const_iterator elements_begin() const { return elements.begin(); }
		element_const_iterator elements_end() const   { return elements.end(); }
				
		node_iterator        nodes_begin() { return nodes.begin(); }
		node_const_iterator  nodes_begin() const { return nodes.begin(); }
		node_const_iterator  nodes_end() const { return nodes.end(); }
				
		region_iterator  	   regions_begin() { return regions.begin(); }
		region_const_iterator  regions_begin() const { return regions.begin(); } 
		region_const_iterator  regions_end() const {  return regions.end(); }		

		boundary_iterator  	     element_boundaries_begin() { return element_boundaries.begin(); }
		boundary_const_iterator  element_boundaries_begin() const { return element_boundaries.begin(); } 
		boundary_const_iterator  element_boundaries_end() const {  return element_boundaries.end(); }
		 									
		void set_identifier(const string& identifier_) { identifier = identifier_; }
		const string& get_identifier() const { return identifier; }																	
		
		void   set_materials(MaterialDatabase& material_database);
		
		//double calculate_quantized_volume() const;						
		bool boundary_material(unsigned int material_id) const;
		const vector<bool>& get_boundary_materials() const { return boundary_materials; }
																			
																				
	private:	
	
		void update_boundary_materials();
				
		std::vector<Element*>         elements;
		std::vector<Node*>            nodes;
		std::vector<Region*>          regions;
		std::vector<ElementBoundary*> element_boundaries;		
							 		
		unsigned int nfaces;
		unsigned int nedges;
		
		unsigned int next_elem;	
		unsigned int next_node;
		unsigned int next_region;
		unsigned int next_boundary;	
					
		unsigned int dimension;	

		string identifier;
		
		BoundaryCondition* bc_handler;
		vector<bool>       boundary_materials;
			
	};


	/** standard quantized volume determination 
	 * 
	 * calculates minimum band edges on outer boundary. everything with 
	 * bandedge below these minima is counted to the quantized region
	 */
	class QuantizedVolumeCalculator {
	public:
		QuantizedVolumeCalculator(const Geometry& geometry);
		virtual ~QuantizedVolumeCalculator() {}
		virtual double calculate_volume() const;
		const Geometry& get_geometry() const { return geometry; }		
	private:
		const Geometry& geometry;		
	};
	
	/** quantized volume calculator based on user defined quantized region names */
	class UserDefinedQuantizedVolumeCalculator : public QuantizedVolumeCalculator {
	public:
		UserDefinedQuantizedVolumeCalculator(const Geometry& geometry);
		virtual ~UserDefinedQuantizedVolumeCalculator() {}
		virtual double calculate_volume() const;
		void add(const char* quantized_region_name);
	private:
		bool is_quantized(const string& region_name) const;
		vector<string> quantized_region_names;		
	};
	
	/** return node with global idx node_idx*/
	inline const Node& Geometry::get_node (unsigned int node_idx) const {
		TDKP_BOUNDS_ASSERT(node_idx >= 0 && node_idx < get_num_nodes(), "");
		return *nodes[node_idx];	
	}
	/** return element */
	inline const Element& Geometry::get_element(unsigned int elem_idx) const {
		TDKP_BOUNDS_ASSERT(elem_idx >= 0 && elem_idx < get_num_elements(), "");
		return *elements[elem_idx];	
	}
	
	
	
} // end of namespace
#endif

