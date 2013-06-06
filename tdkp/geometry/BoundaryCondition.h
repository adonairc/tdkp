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

#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include "tdkp/common/all.h"
#include "tdkp/geometry/Geometry.h"

namespace tdkp {
		
	class Geometry;		
		
	class BoundaryCondition	{
	public:
		BoundaryCondition(const Geometry& geometry_) : geometry(geometry_) {}		
		virtual ~BoundaryCondition() {}
		/** returns true if node should be discarded (e.g. dirichlet bc) */
		virtual bool ignore_node_fully(const Node& node) const = 0;
		/** returns true if only a specific dof should be disregarded (e.g. dirichlet in strain for only some directions) */
		virtual bool ignore_node_partially(const Node& node, unsigned int dof) const { return ignore_node_fully(node); }
		/** returns true if at least one node may be partially ignored */
		virtual bool have_ignore_partially() const { return false; } 		
		const Geometry& get_geometry() const { return geometry; }			
	private:
		const Geometry& geometry;
									
	};	

	/** standard dirichlet boundary condition enforcing zero at boundary */
	class BCDirichlet : public BoundaryCondition {
	public:
		BCDirichlet(const Geometry& geometry_) : BoundaryCondition(geometry_) {}		
		virtual ~BCDirichlet() {}
		virtual bool ignore_node_fully(const Node& node) const {
			return node.get_value_type() == Node::FunctionValue 
			    && node.get_location() == location_exterior;
		}
	};
	
	/** include all nodes (used for matrix elements) */
	class BCIncludeAll : public BoundaryCondition {
	public:
		BCIncludeAll(const Geometry& geometry_) : BoundaryCondition(geometry_) {}
		virtual ~BCIncludeAll() {}
		virtual bool ignore_node_fully(const Node& node) const { return false; }
	};
	
	/** 1D boundary condition that lets outer node float 
	 * 
	 * to be used for 1D biaxial strained quantum well strain calculation
	 */
	class BCWellStrainFloatRight : public BoundaryCondition {
	public:
		BCWellStrainFloatRight(const Geometry& geometry_);		
		virtual bool ignore_node_fully(const Node& node) const;
	private:
		int dirichlet_node;		
	};
	
	/** 2d boundary conditions for defining dirichlet u = 0 boundary conditions */
	class BCDirichlet2DLines : public BoundaryCondition {
	public:
		BCDirichlet2DLines(const Geometry& geometry_);
		virtual bool ignore_node_fully(const Node& node) const;
		virtual bool ignore_node_partially(const Node& node, unsigned int dof) const;
		virtual bool have_ignore_partially() const { return have_partial_dc; } 				
		void add_dirichlet_line(const double& start_x, const double& start_y, const double& end_x, const double& end_y, int dof = -1);
		void prepare();
		void set_tolerance(const double& tolerance_);
		
	private:
		struct LineSegment {
			double start[2];
			double end[2];
			int    dof;
		};
		vector<LineSegment>   segments;
		vector<bool>          dirichlet_nodes;
		vector<vector<bool> > partial_dirichlet_nodes;
		bool         have_partial_dc;
		double       tolerance;
	};
	
	/** user defined planes as dirichlet boundary conditions */
	class BCDirichlet3DPlanes : public BoundaryCondition {
	public:
		BCDirichlet3DPlanes(const Geometry& geometry_);
		void add_dirichlet_plane(const Vector3D& r0, const Vector3D& plane_x, const Vector3D& plane_y, int dof = -1);
		void set_tolerance(const double& tolerance);
		bool ignore_node_fully(const Node& node) const;
		bool ignore_node_partially(const Node& node, unsigned int dof) const;
		bool have_ignore_partially() const { return have_partial_dc; } 							
		void prepare();
				
	private:				
		struct Plane {
			Vector3D r0,a,b;
			int dof;	
		};		
		vector<Plane> planes;
		vector<bool>  dirichlet_node;
		vector<vector<bool> > partial_dirichlet_nodes;		
		double tolerance;
		bool have_partial_dc;
		
	};
	
	/** boundary condition which additionally ignores any dof belonging to gas */
	class BCIgnoreGas : public BoundaryCondition {
	public:
		BCIgnoreGas(BoundaryCondition* basic_boundary_conditions, const MaterialDatabase& material_database);
		~BCIgnoreGas();				
		virtual bool ignore_node_fully(const Node& node) const;		
		virtual bool ignore_node_partially(const Node& node, unsigned int dof) const;		
		virtual bool have_ignore_partially() const { return basic_bc->have_ignore_partially(); } 		
				
	private:
		BoundaryCondition* basic_bc;	
		vector<bool> gas_node;	
	};


}

#endif /*BOUNDARYCONDITION_H_*/
