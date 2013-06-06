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

#include "tdkp/geometry/BoundaryCondition.h"

using namespace tdkp;


BCWellStrainFloatRight::BCWellStrainFloatRight(const Geometry& geometry_) 
: BoundaryCondition(geometry_), dirichlet_node(0) {

	// -----------------------------------------
	// find left node
	// -----------------------------------------
	TDKP_ASSERT(get_geometry().get_dimension() == 1, "get_geometry.get_dimension() == 1");
	double coord_left = get_geometry().get_node(0).get_coord(0);
	dirichlet_node  = get_geometry().get_node(0).get_index_global();
	for(unsigned int ii = 1; ii < get_geometry().get_num_nodes(); ii++) {
		// but only funciton value (hermite have derivatives)
		if(get_geometry().get_node(ii).get_coord(0) <= coord_left 
		&& get_geometry().get_node(ii).get_value_type() == Node::FunctionValue) {
			coord_left      = get_geometry().get_node(ii).get_coord(0);
			dirichlet_node  = get_geometry().get_node(ii).get_index_global();				
		}				
	}

}

bool BCWellStrainFloatRight::ignore_node_fully(const Node& node) const { 
	if(dirichlet_node == node.get_index_global()) {
		return true;	
	}
	return false;
}

/** constructor for user defined 3D dirichlet planes */
BCDirichlet3DPlanes::BCDirichlet3DPlanes(const Geometry& geometry_)
: BoundaryCondition(geometry_), 
  tolerance(1.0e-6),
  have_partial_dc(false) 
{
}

/** add plane vectors defining plane of dirichlet nodes 
 * 
 * @param r0      vector to a point of the plane 
 * @param plane_x plane vector 1
 * @param plane_y plane vector 2
 * @param dof     enforce dirichlet only for degree of freedom dof; -1 means: dirichlet for all dofs
 */
void BCDirichlet3DPlanes::add_dirichlet_plane(const Vector3D& r0, const Vector3D& plane_x, const Vector3D& plane_y, int dof) {
	Plane p;
	p.r0 = r0; 
	p.a  = plane_x;
	p.b  = plane_y;
	p.a.normalize();
	p.b.normalize();
	p.dof = dof; // degree of freedom	
	planes.push_back(p);
	if(dirichlet_node.size() > 0) {
		dirichlet_node.clear();
	}
	// grow partial_dirichlet_nodes
	if(dof >= static_cast<int>(partial_dirichlet_nodes.size())) {
		partial_dirichlet_nodes.resize(dof + 1);
		have_partial_dc = true;	
	}	
	
}

/** set tolerance of how close points must be to be included into the plane */
void BCDirichlet3DPlanes::set_tolerance(const double& tolerance_) {
	this->tolerance = tolerance_;	
} 
/** virtual function return true if node belongs to user defined plane */
bool BCDirichlet3DPlanes::ignore_node_fully(const Node& node) const {
	TDKP_ASSERT(planes.size() > 0, "you didn't define any plane for BCDirichlet3DPlanes. so no boundary conditions set");
	TDKP_ASSERT(dirichlet_node.size() > 0, "you added planes to BCDirichlet3DPlanes but you didn't call prepare before using the class");
	TDKP_ASSERT((signed)dirichlet_node.size() > node.get_index_global(),"dirichlet_node.size() > node.get_index_global()"); 
	return dirichlet_node[node.get_index_global()]; 		
}

bool BCDirichlet3DPlanes::ignore_node_partially(const Node& node, unsigned int dof) const {
	if(dof < partial_dirichlet_nodes.size()) {		
		return partial_dirichlet_nodes[dof][node.get_index_global()];
	} else {
		// dof is free
		return false;	
	}
}

/** calculate dirichlet nodes for defined planes */ 
void BCDirichlet3DPlanes::prepare() {
		
	// --------------------------------------
	// reinit parital dirichlet nodes to false
	// --------------------------------------
	for(unsigned int ii = 0; ii < partial_dirichlet_nodes.size(); ii++) {
		partial_dirichlet_nodes[ii].assign(get_geometry().get_num_nodes(), false);	
	}		
	
	// --------------------------------------
	// set all nodes to neumann
	// --------------------------------------		
	dirichlet_node.assign(this->get_geometry().get_num_nodes(), false);

	// --------------------------------------
	// for all planes
	// --------------------------------------
	unsigned int ncount = 0;
	for(unsigned int ii = 0; ii < planes.size(); ii++) {
		// ---------------------------------------------
		// get plane vectors
		// ---------------------------------------------
		const Vector3D& p1 = planes[ii].a;
		const Vector3D& p2 = planes[ii].b;
		const Vector3D& r0 = planes[ii].r0;
		Vector3D p_ortho   = Vector3D::cross_product(p1,p2);
		p_ortho.normalize();
		
		double m[9];
		double coeff[3];
		for(unsigned int jj = 0; jj < 3; jj++) {
			m[jj * 3] = p1(jj);
			m[jj * 3 + 1] = p2(jj);
			m[jj * 3 + 2] = p_ortho(jj); 	
		}
		
		// -----------------------------------------------------
		// set correct vector where we write into
		// -----------------------------------------------------
		vector<bool>* tmp_seg_vec = 0;
		if(planes[ii].dof == -1) {
			tmp_seg_vec = &dirichlet_node;	
		} else {
			TDKP_BOUNDS_ASSERT(static_cast<int>(partial_dirichlet_nodes.size()) > planes[ii].dof, "");
			tmp_seg_vec = &(partial_dirichlet_nodes[planes[ii].dof]);	
		}
		vector<bool>& r_dc_nodes = *tmp_seg_vec;		
			
		// ---------------------------------------------
		// check all nodes
		// ---------------------------------------------
		for(unsigned int vv = 0; vv < this->get_geometry().get_num_nodes(); vv++) {
			if(!r_dc_nodes[vv]) {				
				Vector3D(this->get_geometry().get_node(vv).get_coords());
				Vector3D tp = Vector3D(this->get_geometry().get_node(vv).get_coords()) - r0;
				tdkp_math::solve_3x3_system(m, tp.get_all(), coeff);
				if(fabs(coeff[2]) < tolerance) {					
					r_dc_nodes[vv] = true;
					ncount++;	
				}
			}	
		}				
	}
		
	if(ncount == 0) {
		TDKP_GENERAL_EXCEPTION("BCDirichlet3DPlanes: your dirichlet boundary conditions do not contain any node.");	
	}

	// ------------------------------------------------
	// emit full dirichlet nodes
	// ------------------------------------------------
	ostringstream sout;	
	unsigned num_dc = 0;	
	for(unsigned int vv = 0; vv < this->get_geometry().get_num_nodes(); vv++) {
		if(dirichlet_node[vv]) {
			if(num_dc > 0) {
				sout << ", ";	
			}
			num_dc++;
			sout << this->get_geometry().get_node(vv).get_index_global();
		}
	} 
	if(num_dc > 0) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "BCDirichlet3DPlanes: on the following " << num_dc << " nodes, all dofs will be enforced to 0: " << sout.str()); 	
	} else {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "BCDirichlet3DPlanes: no nodes with all dofs enforced to 0");		
	}
	
	// ----------------------------------------------
	// emit partial dirichlet nodes
	// ----------------------------------------------
	for(unsigned int dd = 0; dd < partial_dirichlet_nodes.size(); dd++) {
		num_dc = 0;
		sout.str("");
		for(unsigned int vv = 0; vv < this->get_geometry().get_num_nodes(); vv++) {
			if(partial_dirichlet_nodes[dd][vv]) {
				if(num_dc > 0) {
					sout << ", ";	
				}
				num_dc++;
				sout << this->get_geometry().get_node(vv).get_index_global();
			}
		}
		if(num_dc > 0) {
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "BCDirichlet3DPlanes: on the following " << num_dc << " nodes, dof " << dd << " will be enforced to 0: " << sout.str()); 	
		} else {
			// quit if dirichlet condition was defined but never used
			for(unsigned int pp = 0; pp < planes.size(); pp++) {
				if(planes[pp].dof == (signed)dd) {										
					TDKP_GENERAL_EXCEPTION("BCDirichlet3DPlanes: no nodes will have dof " << dd <<" enforced to 0. i guess that this shouldn't be (as it wouldn't make sense to specify a boundary condition which does not affect anything ...)");
				}
			}		
		}						
	} 	
}

BCDirichlet2DLines::BCDirichlet2DLines(const Geometry& geometry_)
: BoundaryCondition(geometry_),
  have_partial_dc(false),
  tolerance(1.0e-6)  
{
}

/** set all points on the infinite line defined by the two points to dirichlet bc 
 * 
 * @param start_x  first point
 * @param start_y  first point
 * @param end_x    second point
 * @param end_y    second point
 * @param dof      degree of freedom (-1 means all)
 */
void BCDirichlet2DLines::add_dirichlet_line(
	const double& start_x, const double& start_y, const double& end_x, const double& end_y,
	int dof
) {
	// store segment
	LineSegment seg;
	seg.start[0] = start_x;
	seg.start[1] = start_y;
	seg.end[0]   = end_x;
	seg.end[1]   = end_y;
	seg.dof      = dof; // degree of freedom	
	segments.push_back(seg);
	// delete nodes
	if(dirichlet_nodes.size() > 0) {
		dirichlet_nodes.clear();	
	}
	// grow partial_dirichlet_nodes
	if(dof >= static_cast<int>(partial_dirichlet_nodes.size())) {
		partial_dirichlet_nodes.resize(dof + 1);
		have_partial_dc = true;	
	}
} 

void BCDirichlet2DLines::set_tolerance(const double& tolerance_) {
	tolerance = tolerance_;	
}

void BCDirichlet2DLines::prepare() {
	
	// reinit parital dirichlet nodes to false
	for(unsigned int ii = 0; ii < partial_dirichlet_nodes.size(); ii++) {
		partial_dirichlet_nodes[ii].assign(get_geometry().get_num_nodes(), false);	
	}
		
	// init dirichlet nodes to neuman
	dirichlet_nodes.assign(get_geometry().get_num_nodes(), false);
	unsigned int ncount = 0;
	
	// ----------------------------------------------------
	// for every segment
	// ----------------------------------------------------
	for(unsigned int ss = 0; ss < segments.size(); ss++) {
		const LineSegment& seg = segments[ss];
		Vector3D r0(seg.start[0], seg.start[1], 0.0);
		Vector3D a(seg.end[0] - seg.start[0], seg.end[1] - seg.start[1], 0.0);
		Vector3D an = a; an.normalize();

		// -----------------------------------------------------
		// set correct vector where we write into
		// -----------------------------------------------------
		vector<bool>* tmp_seg_vec = 0;
		if(seg.dof == -1) {
			tmp_seg_vec = &dirichlet_nodes;	
		} else {
			TDKP_BOUNDS_ASSERT(static_cast<int>(partial_dirichlet_nodes.size()) > seg.dof, "");
			tmp_seg_vec = &(partial_dirichlet_nodes[seg.dof]);	
		}
		vector<bool>& r_dc_nodes = *tmp_seg_vec;

		// ----------------------------------------------------
		// test every node
		// ----------------------------------------------------
		for(unsigned int ii = 0; ii < get_geometry().get_num_nodes(); ii++) {
			const Node& node = get_geometry().get_node(ii);
			if(!r_dc_nodes[node.get_index_global()]) {
				Vector3D nvec(node.get_coord(0) - seg.start[0], node.get_coord(1) - seg.start[1], 0.0);
				Vector3D nvecnormal = nvec;
				// if nvec is r0, we accept it as dirichlet node
				if(nvec.norm() < tolerance) {
					r_dc_nodes[node.get_index_global()] = true;
				} else {
					nvecnormal.normalize();
					// node position is aligned
					if(tdkp_math::abs(Vector3D::dot_product(an, nvecnormal) - 1.0) < tolerance && 
					   nvec.norm() - a.norm() < tolerance) { 					  					
						r_dc_nodes[node.get_index_global()] = true;
						ncount++;
					} 
				}
			} 
		}			 
	}
	if(ncount == 0) {
		TDKP_GENERAL_EXCEPTION("BCDirichlet2DLines: your dirichlet boundary conditions do not contain any node.");
	}
}
bool BCDirichlet2DLines::ignore_node_fully(const Node& node) const {
	TDKP_BOUNDS_ASSERT(node.get_index_global() < (signed)dirichlet_nodes.size(), "did you call prepare on the dirichlet 2D boundary object?");
	return dirichlet_nodes[node.get_index_global()];	
}

bool BCDirichlet2DLines::ignore_node_partially(const Node& node, unsigned int dof) const {
	if(dof < partial_dirichlet_nodes.size()) {
		return partial_dirichlet_nodes[dof][node.get_index_global()];
	} else {
		// dof is free
		return false;	
	}
}

// -------------------------------------------------------
// BCIgnoreGas implementation
// -------------------------------------------------------
BCIgnoreGas::BCIgnoreGas(BoundaryCondition* basic_boundary_conditions, const MaterialDatabase& material_database)
: BoundaryCondition(basic_boundary_conditions->get_geometry()),
  basic_bc(basic_boundary_conditions),
  gas_node(get_geometry().get_num_nodes(), false)
{
	// -----------------------------------------------
	// get index of material Gas
	// -----------------------------------------------
	unsigned int id_gas = 0;
	bool found_gas = false;	
	for(int ii = 0; ii < material_database.get_num_materials(); ii++) {
		if(material_database.get_material_name(ii) == string("Gas")) {
			TDKP_ASSERT(!found_gas, "material Gas is loaded twice in material database! first idx is " << id_gas << " and second idx is " << ii);
			found_gas = true;
			id_gas = ii;   			
		}	
	}
		
	// -----------------------------------------------
	// determine nodes belonging to gas
	// -----------------------------------------------
	if(found_gas) {
		
		for(unsigned int ii = 0; ii < get_geometry().get_num_elements(); ii++) {
			const Element& elem = get_geometry().get_element(ii);
			if(elem.enabled() && elem.get_material().get_id() == id_gas) {
				for(unsigned nn = 0; nn < elem.get_num_nodes(); nn++) {
					gas_node[elem.get_node(nn).get_index_global()] = true;
				}	
			}
		} 	
		
		// statistics
		int count_gas_nodes = 0;
		for(unsigned int ii = 0; ii < gas_node.size(); ii++) {
			if(gas_node[ii]) {
				++count_gas_nodes;	
			}						
		}
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "BCIgnoreGas: material " << id_gas << " is Gas and " << count_gas_nodes << " nodes therefrom will be ignored");
	}
	
}
BCIgnoreGas::~BCIgnoreGas() {
	if(basic_bc != 0) {
		delete basic_bc; basic_bc = 0;	
	}	
}
				
bool BCIgnoreGas::ignore_node_fully(const Node& node) const {
	return basic_bc->ignore_node_fully(node) || gas_node[node.get_index_global()];	
}
		
bool BCIgnoreGas::ignore_node_partially(const Node& node, unsigned int dof) const {
	return basic_bc->ignore_node_partially(node,dof) || gas_node[node.get_index_global()];
}

