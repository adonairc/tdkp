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

#include "tdkp/geometry/Element2DRect.h"
#include "tdkp/common/Vector3D.h"

namespace tdkp {
	
const short dxN0dxN0 = 0;
const short dxN0dxN1 = 1;
const short dxN0dxN2 = 2;
const short dxN0dxN3 = 3;

const short dxN1dxN0 = 4;
const short dxN1dxN1 = 5;
const short dxN1dxN2 = 6;
const short dxN1dxN3 = 7;

const short dxN2dxN0 = 8;
const short dxN2dxN1 = 9;
const short dxN2dxN2 = 10;
const short dxN2dxN3 = 11;

const short dxN3dxN0 = 12;
const short dxN3dxN1 = 13;
const short dxN3dxN2 = 14;
const short dxN3dxN3 = 15;

const short dxN0dyN0 = 16;
const short dxN0dyN1 = 17;
const short dxN0dyN2 = 18;
const short dxN0dyN3 = 19;

const short dxN1dyN0 = 20;
const short dxN1dyN1 = 21;
const short dxN1dyN2 = 22;
const short dxN1dyN3 = 23;
 
const short dxN2dyN0 = 24;
const short dxN2dyN1 = 25;
const short dxN2dyN2 = 26;
const short dxN2dyN3 = 27;
 
const short dxN3dyN0 = 28;
const short dxN3dyN1 = 29;
const short dxN3dyN2 = 30;
const short dxN3dyN3 = 31;

const short dyN0dxN0 = 32;
const short dyN0dxN1 = 33;
const short dyN0dxN2 = 34;
const short dyN0dxN3 = 35;

const short dyN1dxN0 = 36;
const short dyN1dxN1 = 37;
const short dyN1dxN2 = 38;
const short dyN1dxN3 = 39;

const short dyN2dxN0 = 40;
const short dyN2dxN1 = 41;
const short dyN2dxN2 = 42;
const short dyN2dxN3 = 43;

const short dyN3dxN0 = 44;
const short dyN3dxN1 = 45;
const short dyN3dxN2 = 46;
const short dyN3dxN3 = 47;

const short dyN0dyN0 = 48;
const short dyN0dyN1 = 49;
const short dyN0dyN2 = 50;
const short dyN0dyN3 = 51;

const short dyN1dyN0 = 52;
const short dyN1dyN1 = 53;
const short dyN1dyN2 = 54;
const short dyN1dyN3 = 55;

const short dyN2dyN0 = 56;
const short dyN2dyN1 = 57;
const short dyN2dyN2 = 58;
const short dyN2dyN3 = 59;

const short dyN3dyN0 = 60;
const short dyN3dyN1 = 61;
const short dyN3dyN2 = 62;
const short dyN3dyN3 = 63;


const short dxN0N0 = 0;
const short dxN0N1 = 1;
const short dxN0N2 = 2;
const short dxN0N3 = 3;

const short dxN1N0 = 4;
const short dxN1N1 = 5;
const short dxN1N2 = 6;
const short dxN1N3 = 7;

const short dxN2N0 = 8;
const short dxN2N1 = 9;
const short dxN2N2 = 10;
const short dxN2N3 = 11;

const short dxN3N0 = 12;
const short dxN3N1 = 13;
const short dxN3N2 = 14;
const short dxN3N3 = 15;

const short dyN0N0 = 16;
const short dyN0N1 = 17;
const short dyN0N2 = 18;
const short dyN0N3 = 19;

const short dyN1N0 = 20;
const short dyN1N1 = 21;
const short dyN1N2 = 22;
const short dyN1N3 = 23;

const short dyN2N0 = 24;
const short dyN2N1 = 25;
const short dyN2N2 = 26;
const short dyN2N3 = 27;
 
const short dyN3N0 = 28;
const short dyN3N1 = 29;
const short dyN3N2 = 30;
const short dyN3N3 = 31;


// ----------------------------------------------------------------
// ELEMNENT2dRectBASE
// ----------------------------------------------------------------

/** 2d basic rectangle constructor */
Element2DRectBase::Element2DRectBase(unsigned int dimension_, unsigned int nnodes_, unsigned int nboundaries_)
: Element(dimension_, nnodes_, nboundaries_),
  ax(0.0),
  ay(0.0),
  bx(0.0),
  by(0.0),
  axby_m_aybx(0.0),
  vertex_offset(1) 
{
	a[0] = a[1] = b[0] = b[1] = 0.0e0;
	this->ready		   	  = false;
	for(int ii = 0; ii < 4; ii++) {
		this->inverted_jacobi_matrix[ii] = 0;			
	}			
}

bool Element2DRectBase::verify() const {
	
	bool good = true;
	
	if(region == 0) {
		Logger::get_instance()->emit(LOG_ERROR, "no region set to element");
		good = false;	
	}
	for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
		if(nodes[ii] == 0) {
			Logger::get_instance()->emit(LOG_ERROR, 1000, "node %d in element is null", ii);			
			return false;	
		}
		if(this->nodes[ii]->get_dimension() != 2) {
			ostringstream sout;
			sout << "node " << ii << " has wrong dimension " << this->nodes[ii]->get_dimension() << " != 2";
			Logger::get_instance()->emit(LOG_ERROR, sout.str());
			return false;
		}	
		for(unsigned int jj = ii + 1; jj < get_num_nodes(); jj++) {
			if(nodes[ii] == nodes[jj]) {
				Logger::get_instance()->emit(LOG_ERROR, 1000, "nodes %d and %d in element are equal", ii, jj);
				return false;	
			}
		}			
	}	
	if(element_volume <= 0.0) {
		Logger::get_instance()->emit(LOG_ERROR, 1000, "volume of element is %5.5g and therefore incorrect", element_volume);
		return false;	
	}
	if(jacobi_det <= 0.0) {
		Logger::get_instance()->emit(LOG_ERROR, 1000, "jacobi deterimnante is %5.5g and therefore incorrect", jacobi_det);
		return false;	
	} 
	// --------------------------------------------------
	// check if corners are 90 deg
	// --------------------------------------------------
	Vector3D edges[4];
	int next;
	for(unsigned int ii = 0; ii < get_num_nodes(); ii += vertex_offset) {
		next = (ii + vertex_offset)	% this->get_num_nodes();
		edges[ii / vertex_offset] = Vector3D(this->nodes[next]->get_coords()) - Vector3D(this->nodes[ii]->get_coords());
		/* 
		edges[ii][0] = this->nodes[next]->get_coord(0) - this->nodes[ii]->get_coord(0);
		edges[ii][1] = this->nodes[next]->get_coord(1) - this->nodes[ii]->get_coord(1);
		*/		
	}
	double dot;
	static int bad_elements = 10;
	for(unsigned int ii = 0; ii < get_num_nodes(); ii += vertex_offset) {
		next = (ii + vertex_offset)	% this->get_num_nodes();		
		dot  = Vector3D::dot_product(edges[next / vertex_offset], edges[ii / vertex_offset]); 
		if(fabs(dot) > 1.0e-10) {
			ostringstream sout;
			sout << "edge " << ii / vertex_offset << " and " << next / vertex_offset << " are not orthogonal";
			Logger::get_instance()->emit(LOG_ERROR, sout.str());
			good = false;	
		}
		Vector3D cross = Vector3D::cross_product(edges[ii / vertex_offset], edges[next / vertex_offset]);
		TDKP_ASSERT(cross(0) == cross(1) && cross(0) == 0.0e0, "2D elements must not have a 3D component ...");
		
		if(cross(2) <= 0.0) {
			if(bad_elements >= 0) {
				ostringstream sout;
				sout << "bad rectangular element number " << get_index_global() <<":\n"
				     << "  nodes:\n";
				for(int vv = 0; vv < (signed)this->get_num_nodes(); vv++) {
					sout << "  " << vv << ": " << this->get_node(vv) << "\n"; 	
				}
				sout << "  edges:\n";
				for(int vv = 0; vv < (signed)this->get_num_nodes(); vv++) {
					sout << "  " << vv << " -> " << ((vv + 1) % this->get_num_nodes()) << ": " << edges[vv] << "\n"; 	
				}				     
				Logger::get_instance()->emit(LOG_ERROR, sout.str());							 				
				bad_elements--;
				if(bad_elements < 0) {
					Logger::get_instance()->emit(LOG_WARN, "stopping to show further bad elements!");	
				}
				return false;
			}
		}
	}
	return good;
}

/** return vertices of edge */
void Element2DRectBase::get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const {
	TDKP_BOUNDS_ASSERT(edge_idx < 4, "");
	vertex_indices.resize(2);
	vertex_indices[0] = get_corner_node(edge_idx).get_index_vertex();
	vertex_indices[1] = get_corner_node((edge_idx + 1) % 4).get_index_vertex();
}


/** return vertex indices of boundary */
void Element2DRectBase::get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const {
	get_edge(idx, vertex_indices);
}

/** prepares the element for assembly
 * 
 * we calculate the parameters of the element, the jacobian and 
 * prepare some terms which are necessary to evaluate the element 
 * integrals analytically.
 * 
 * in order to calculate the integrals analytically, one must express
 * the rectangle in terms of the two vectors a,b
 * so we say that
 * P0 (lower left corner)
 * P1 = P0 + a   
 * P2 = P0 + a + b
 * P3 = P0 + b
 * 
 * here Pi are vertex corners. the edge functions points are at mid-edge.
 *
 * so the affine transformation between the global and reference coordinate system
 * is simply given by 
 *    x = sum Pi Ni
 * as for the triangles ... 
 * 
 * using that transformation and the corresponding jacobian, the element integrals
 * can be evaluated analytically (using mathematica ...)
 * 
 */
void Element2DRectBase::prepare() {

	// ----------------------------------------------
	// check if nodes are set
	// ----------------------------------------------
	for(unsigned int ii = 0; ii < nodes.size(); ii++) {
		TDKP_ASSERT(nodes[ii] != 0, "node " << ii << " is not set");		 
	}
	
	TDKP_ASSERT(this->region != 0, "");
	
	this->element_mid_point[0] = 0.0; // local coordinates
	this->element_mid_point[1] = 0.0; // local coordinates

	// ----------------------------------------------
	// check if nodes are defined subsequently
	// ----------------------------------------------
	Vector3D edges[4];		
	// build edges
	for(unsigned int ii = 0; ii < get_num_nodes(); ii+= vertex_offset) {
		int next = (ii + vertex_offset)	% this->get_num_nodes();
		edges[ii / vertex_offset] = Vector3D(this->nodes[next]->get_coords()) - Vector3D(this->nodes[ii]->get_coords());
	}
	bool counterclockwise = true;								
	// calculate dot product between edges		
	for(unsigned int ii = 0; ii < get_num_nodes(); ii+= vertex_offset) {
		int next = (ii + vertex_offset)	% this->get_num_nodes();
		double dot = Vector3D::dot_product(edges[next / vertex_offset], edges[ii / vertex_offset]);
		if(fabs(dot) > 1.0e-10) {
			ostringstream sout;
			sout << "found a bad rectangular 2D element number " << this->get_index_global() << ":\n"
			     << "the vertices are not ordered subsequently counterclockwise!\n"
			     << "  the vertices are (in that particular order):\n";
			for(unsigned int jj = 0; jj < this->get_num_nodes(); jj+= vertex_offset) {
				sout << *this->nodes[jj] << "\n";	
			}
			sout << "the problem was detected on vertex " << (ii + vertex_offset) % this->get_num_nodes() << ".\n"
			     << "error means: the dot product between subsequent edges was > 1.0e-10."; 			     
			TDKP_GENERAL_EXCEPTION(sout.str());						     	
		}			
		Vector3D cross = Vector3D::cross_product(edges[next / vertex_offset], - edges[ii / vertex_offset]);

		TDKP_ASSERT(cross(0) == cross(1) && cross(0) == 0.0e0, "2D elements must not have a 3D component ...");		
		if(cross(2) <= 0.0) {
			counterclockwise = false;
		}
		TDKP_ASSERT(counterclockwise || cross(2) < 0.0, "either its counter clock wise or not, but not both!");
	}
	// ----------------------------------------------
	// if the nodes are ordered clockwise, we need to swap them!
	// ----------------------------------------------	
	if(!counterclockwise) {
		static int complain_clockwise = 10;
		if(complain_clockwise > 0) {
			complain_clockwise--;
			ostringstream scomplain;
			scomplain << "nodes in element " << this->get_index_global() << " are ordered clockwise. swapping it!";
			Logger::get_instance()->emit(LOG_WARN, scomplain.str());
			if(complain_clockwise == 0) {
				Logger::get_instance()->emit(LOG_WARN, "will not show further problems on wrong ordered elements!");				
			}			
		}
		// -----------------------------------
		// swap
		// -----------------------------------
		vector<Node*> tmp_nodes(get_num_nodes());
		for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
			tmp_nodes[ii] = nodes[ii];	
		}
		for(unsigned int ii = 0; ii < get_num_nodes(); ii++) {
			// that we, node 0 stays node 0, important for higher order elements as
			// at nodes ii * vertex_offset we expect vertices
			nodes[ii] = tmp_nodes[(get_num_nodes() - ii) % get_num_nodes()];	
		}
		/*
		Node* tmp;
		tmp = this->nodes[0]; 
		this->nodes[0] = this->nodes[3];
		this->nodes[3] = tmp;
		tmp = this->nodes[1];
		this->nodes[1] = this->nodes[2];
		this->nodes[2] = tmp;
		*/		
	}
	
	// ----------------------------------------------
	// check if vertices are where we expect them
	// ----------------------------------------------
	for(unsigned int ii = 0; ii < nodes.size(); ii++) {
		if(ii % vertex_offset == 0) {
			TDKP_ASSERT(nodes[ii]->get_index_vertex() != -1, "");
		} else {
			TDKP_ASSERT(nodes[ii]->get_index_vertex() == -1, "");	
		}
	}
	// ----------------------------------------------	
	// calculate a and b (a = P1 - P0, b = P3 - P0)
	// and just test that P0 + a + b gives P2 up to 1.0e-10 
	// ----------------------------------------------
	double P2[2];
	for(short ii = 0; ii < 2; ii++) {
		this->a[ii] = this->nodes[1 * vertex_offset]->get_coord(ii) - this->nodes[0]->get_coord(ii);
		this->b[ii] = this->nodes[3 * vertex_offset]->get_coord(ii) - this->nodes[0]->get_coord(ii);
		P2[ii] = this->nodes[0]->get_coord(ii) + this->a[ii] + this->b[ii];		  			
	}
	if(fabs(P2[0] - this->nodes[2 * vertex_offset]->get_coord(0)) > 1.0e-9 || fabs(P2[1] - this->nodes[2 * vertex_offset]->get_coord(1)) > 1.0e-9) {
		ostringstream sout;
		sout << "found a bad rectangular 2D element number " << this->get_index_global() << ":\n"
		     << "the sides aren't parallel!\n"
		     << "  the vertices are (in that particular order):\n";
		for(unsigned int jj = 0; jj < this->get_num_nodes(); jj += vertex_offset) {
			sout << *this->nodes[jj] << "\n";	
		}
		TDKP_GENERAL_EXCEPTION(sout.str());
	}
		
	this->ax = a[0]; this->ay = a[1];
	this->bx = b[0]; this->by = b[1];		

	this->axby_m_aybx = ax*by - ay*bx;
	TDKP_ASSERT(this->axby_m_aybx > 0, "check orientation of nodes in elmement. Det[B] must be positive (i assumed that during element integral derivation) or the element integrals may be wrong.");		
//	cout << "axby_m_aybx: " << this->axby_m_aybx << "\n";
				
	// ----------------------------------------------
	// calculate jacobian and element volume
	// ----------------------------------------------
	this->element_volume = tdkp_math::abs(axby_m_aybx);
	this->jacobi_det     = tdkp_math::abs(axby_m_aybx) / 4.0;
		
	// ----------------------------------------------
	// calculate inverse jacobian (ii * 2 + jj)
	// ----------------------------------------------
	this->inverted_jacobi_matrix[0] =   2.0 * by / axby_m_aybx;
	this->inverted_jacobi_matrix[1] = - 2.0 * bx / axby_m_aybx;
	this->inverted_jacobi_matrix[2] = - 2.0 * ay / axby_m_aybx;
	this->inverted_jacobi_matrix[3] =   2.0 * ax / axby_m_aybx;
	
	this->ready = true;
	
	this->verify();	
//	TDKP_ASSERT(this->verify(), "bad element!");
						
}			

/** map global to local variables
 * 
 * so, reference rectangle is at (+/- 1, +/- 1) and corners
 * are counterclockwise starting with N0 at (1,1) (= P0)
 * 
 * to calculate the local coordinates of a global point we simply use:
 * 
 * x = 2.0 * dot_product(P3 - P2, PX - P2) - 1
 * y = 2.0 * dot_product(P1 - P2, PX - P2) - 1 
 */
bool Element2DRectBase::global2local(const double global[2], double local[2]) const {
	
	double tol = 1.0e-12;
	
	Vector3D p0(this->get_node(0).get_coord(0), this->get_node(0).get_coord(1), 0.0);
	Vector3D p1(this->get_node(1 * vertex_offset).get_coord(0), this->get_node(1 * vertex_offset).get_coord(1), 0.0);
	Vector3D p3(this->get_node(3 * vertex_offset).get_coord(0), this->get_node(3 * vertex_offset).get_coord(1), 0.0);
	
	Vector3D px(
		global[0],
		global[1],
		0.0
	);

	Vector3D va = p1 - p0;
	Vector3D vb = p3 - p0;
	Vector3D vx = px - p0;
	
	 	
	local[0] = 2.0 * Vector3D::dot_product(va, vx) / (va.norm() * va.norm()) - 1.0;
	local[1] = 2.0 * Vector3D::dot_product(vb, vx) / (vb.norm() * vb.norm()) - 1.0;
	

	for(short ii = 0; ii < 2; ii++) {
		if(local[ii] > 1.0 + tol || local[ii] < -1.0 - tol) {
			return false;
		}
	} 
	return true;
}

bool   Element2DRectBase::inside_element(const double& global_x, const double& global_y, const double& global_z) const {
	double global[] = {global_x, global_y};
	double local[2];
	return global2local(global,local);
}
		

// ----------------------------------------------------------------
// N E W   I M P L E M E N T A T I O N
// ----------------------------------------------------------------

double Element2DRect::reference_element_integral_zero[4][4] = {
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0, 1.0 / 18.0},
	{1.0 / 18.0, 1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0},
	{1.0 / 36.0, 1.0 / 18.0, 1.0 / 9.0,  1.0 / 18.0},
	{1.0 / 18.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 9.0}
};

double Element2DRect::reference_single_integral_first[2][4][2] = {
	{{1.0,  -1.0}, {1.0,  1.0}, {-1.0, 1.0}, {-1.0, -1.0}},
	{{-1.0,  1.0}, {-1.0,-1.0}, { 1.0,-1.0}, { 1.0,  1.0}}
};

/** 2d linear rectangle constructor */
Element2DRect::Element2DRect(unsigned int index_)
: Element2DRectBase(2,4,4) 
{
	a[0] = a[1] = b[0] = b[1] = 0.0e0;
	this->ready		   	  = false;
	this->index_global    = index_;
	for(int ii = 0; ii < 4; ii++) {
		this->inverted_jacobi_matrix[ii] = 0;			
	}			
}

Element2DRect::~Element2DRect() {
	this->region = 0; 		
}
		



double Element2DRect::get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const {
	
	TDKP_BOUNDS_ASSERT(diffop_1 >= 0 && diffop_1 < 2 && diffop_2 >= 0 && diffop_2 < 2, "diffop_1 >= 0 && diffop_1 < 2 && diffop_2 >= 0 && diffop_2 < 2");
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 4 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 4, "elem_shape_func_1 >= 0 && elem_shape_func_1 < 4 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 4");

	// ----------------------------------------------------------
	// this one is complicated: the 2nd order element integrals
	// can be obtained analytically, but you get complex formulas
	// in terms of phi, a, b etc.	
	// but some of them are equal, 
	//    e.g. d/dx Ni d/dy Nj == d/dy Nj d/dx Ni
	//    which means that call with 0 i 1 j == call with 1 j 0 i
	// altough it is a slight performance increase compared to
	// numerical evaluation of the element integrals, it was a lot
	// of work to implement it ... 
	// ----------------------------------------------------------	
	short select_integral = elem_shape_func_2      // four element shape func 2 
	                      + 4 * elem_shape_func_1  // four element shape func 1 
	                      + 16 * diffop_2         // two diff operators
	                      + 32 * diffop_1;
	                      
	double ret;
	switch(select_integral) {
		case dxN0dxN0:
		case dxN2dxN2:
		    ret = 1.0/6.0 * (2.0 * ay*ay - 3.0 * ay * by + 2.0 * by*by);
		    break;
		
		case dxN0dxN1:
		case dxN3dxN2:
		case dxN1dxN0:
		case dxN2dxN3:
		    ret = 1.0/6.0 * (ay*ay - 2.0 * by*by);
		    break;
		
		case dxN0dxN2:
		case dxN2dxN0:
		    ret = 1.0/6.0 * (-ay*ay + 3.0 * ay * by - by*by);
		    break;
		
		case dxN0dxN3:
		case dxN3dxN0:
		case dxN1dxN2:
		case dxN2dxN1:
		    ret = 1.0/6.0 * (-2.0 * ay*ay + by*by);
		    break;
		
		case dxN0dyN0:
		case dyN2dxN2:
		case dyN0dxN0:
		case dxN2dyN2:
		    ret = 1.0/12.0 * (-4.0 * ax * ay + 3.0 * ay * bx + 3.0 * ax * by - 4.0 * bx * by);
		    break;
		
		case dxN0dyN1:
		case dyN3dxN2:
		case dyN1dxN0:
		case dxN2dyN3:
		    ret = 1.0/12.0 * (-2.0 * ax * ay - 3.0 * ay * bx + 3.0 * ax * by+ 4.0 * bx * by);
		    break;
		
		case dxN0dyN2:
		case dyN0dxN2:
		case dxN2dyN0:
		case dyN2dxN0:
		    ret = 1.0/12.0 * (2.0 * ax * ay - 3.0 * ay * bx - 3.0 * ax * by + 2.0 * bx * by);
		    break;
		
		case dxN0dyN3:
		case dyN1dxN2:
		case dxN2dyN1:
		case dyN3dxN0:
		    ret = 1.0/12.0 * (4.0 * ax * ay + 3.0 * ay * bx - 3.0 * ax * by - 2.0 * bx * by);
		    break;
		
		case dxN1dxN1:
		case dxN3dxN3:
		    ret = 1.0/6.0 * (2.0 * ay*ay + 3.0 * ay * by + 2.0 * by*by);
		    break;
		
		case dxN1dxN3:
		case dxN3dxN1:
		    ret = 1.0/6.0 * (-ay * ay - 3.0 * ay * by - by * by);
		    break;
		
		case dxN1dyN0:
		case dyN2dxN3:
		case dyN0dxN1:
		case dxN3dyN2:
		    ret = 1.0/12.0 * (-2.0 * ax * ay + 3.0 * ay * bx - 3.0 * ax * by + 4 * bx * by);
		    break;
		
		case dxN1dyN1:
		case dyN1dxN1:
		case dxN3dyN3:
		case dyN3dxN3:
		    ret = 1.0/12.0 * (-4.0 * ax * ay - 3.0 * ay * bx - 3.0 * ax * by - 4.0 * bx * by);
		    break;
		
		case dxN1dyN2:
		case dyN2dxN1:
		case dxN3dyN0:
		case dyN0dxN3:
		    ret = 1.0/12.0 * (4.0 * ax * ay-3.0 * ay * bx + 3.0 * ax * by - 2.0 * bx * by);
		    break;
		
		case dxN1dyN3:
		case dxN3dyN1:
		case dyN3dxN1:
		case dyN1dxN3:
		    ret = 1.0/12.0 * (2.0 * ax * ay + 3.0 * ay * bx + 3.0 * ax * by + 2.0 * bx * by);
		    break;
		
		case dyN0dyN0:
		case dyN2dyN2:
		    ret = 1.0/6.0 * (2.0 * ax*ax - 3.0 * ax * bx + 2.0 * bx * bx);
		    break;
		
		case dyN0dyN1:
		case dyN1dyN0:
		case dyN3dyN2:
		case dyN2dyN3:
		    ret = 1.0/6.0 * (ax*ax - 2.0 * bx*bx);
		    break;
		
		case dyN0dyN2:
		case dyN2dyN0:
		    ret = 1.0/6.0 * (-ax*ax + 3.0 * ax * bx - bx*bx);
		    break;
		
		case dyN0dyN3:
		case dyN3dyN0:
		case dyN1dyN2:
		case dyN2dyN1:
		    ret = 1.0/6.0 * (-2.0 * ax*ax + bx*bx);
		    break;
		
		case dyN1dyN1:
		case dyN3dyN3:
		    ret = 1.0/6.0 * (2.0 * ax*ax + 3.0 * ax * bx + 2.0 * bx*bx);
		    break;
		
		case dyN1dyN3:
		case dyN3dyN1:
		    ret = 1.0/6.0 * (-ax*ax - 3.0 * ax * bx - bx*bx);
		    break;
		default:
			TDKP_GENERAL_EXCEPTION("must not reach this point");
	}                      
	return ret / axby_m_aybx;
	TDKP_GENERAL_EXCEPTION("must not reach this point");
}

double Element2DRect::get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const {
	
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 4 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 4, "elem_shape_func_1 >= 0 && elem_shape_func_1 < 4 && elem_shape_func_2 >= 0 && elem_shape_func_2 < 4");
	TDKP_BOUNDS_ASSERT(diffop >= 0 && diffop < 2, "diffop >= 0 && diffop < 2");
	short select = elem_shape_func_2 
	             + 4  * elem_shape_func_1 
	             + 16 * diffop;
	
	double ret;
	switch(select) {	
	
		case dxN0N0:
		    ret = (ay-by) / 6.0;
		    break;
		
		case dxN0N1:
		    ret = 1.0/12.0 *  (ay-2.0 * by);
		    break;
		
		case dxN0N2:
		    ret = (ay-by) / 12.0;
		    break;
		
		case dxN0N3:
		    ret = 1.0/12.0 *  ( 2.0 * ay-by);
		    break;
		
		case dxN1N0:
		    ret = 1.0/12.0 *  (ay+ 2.0 * by);
		    break;
		
		case dxN1N1:
		    ret = (ay+by) / 6.0;
		    break;
		
		case dxN1N2:
		    ret = 1.0/12.0 *  ( 2.0 * ay+by);
		    break;
		
		case dxN1N3:
		    ret = (ay+by) / 12.0;
		    break;
		
		case dxN2N0:
		    ret = 1.0/12.0 *  (-ay+by);
		    break;
		
		case dxN2N1:
		    ret = 1.0/12.0 *  (- 2.0 * ay+by);
		    break;
		
		case dxN2N2:
		    ret = 1.0 / 6.0 * (-ay+by);
		    break;
		
		case dxN2N3:
		    ret = 1.0/12.0 *  (-ay+ 2.0 * by);
		    break;
		
		case dxN3N0:
		    ret = 1.0/12.0 *  (- 2.0 * ay-by);
		    break;
		
		case dxN3N1:
		    ret = 1.0/12.0 *  (-ay-by);
		    break;
		
		case dxN3N2:
		    ret = 1.0/12.0 *  (-ay- 2.0 * by);
		    break;
		
		case dxN3N3:
		    ret = 1.0 / 6.0 * (-ay-by);
		    break;
		
		case dyN0N0:
		    ret = 1.0 / 6.0 * (-ax+bx);
		    break;
		
		case dyN0N1:
		    ret = 1.0/12.0 *  (-ax+ 2.0 * bx);
		    break;
		
		case dyN0N2:
		    ret = 1.0/12.0 *  (-ax+bx);
		    break;
		
		case dyN0N3:
		    ret = 1.0/12.0 *  (- 2.0 * ax+bx);
		    break;
		
		case dyN1N0:
		    ret = 1.0/12.0 *  (-ax- 2.0 * bx);
		    break;
		
		case dyN1N1:
		    ret = 1.0 / 6.0 * (-ax-bx);
		    break;
		
		case dyN1N2:
		    ret = 1.0/12.0 *  (- 2.0 * ax-bx);
		    break;
		
		case dyN1N3:
		    ret = 1.0/12.0 *  (-ax-bx);
		    break;
		
		case dyN2N0:
		    ret = (ax-bx) / 12.0;
		    break;
		
		case dyN2N1:
		    ret = 1.0/12.0 *  ( 2.0 * ax-bx);
		    break;
		
		case dyN2N2:
		    ret = (ax-bx) / 6.0;
		    break;
		
		case dyN2N3:
		    ret = 1.0/12.0 *  (ax- 2.0 * bx);
		    break;
		
		case dyN3N0:
		    ret = 1.0/12.0 *  ( 2.0 * ax+bx);
		    break;
		
		case dyN3N1:
		    ret = (ax+bx) / 12.0;
		    break;
		
		case dyN3N2:
		    ret = 1.0/12.0 *  (ax+ 2.0 * bx);
		    break;
		
		case dyN3N3:
		    ret = (ax+bx) / 6.0;
		    break;
		default:
			TDKP_GENERAL_EXCEPTION("must not reach this point");		    
	}
	return ret;
	TDKP_GENERAL_EXCEPTION("must not reach that point");
}



double Element2DRect::get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const {	
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < (signed)get_num_nodes() && elem_shape_func_2 >= 0 && elem_shape_func_2 < (signed)get_num_nodes(),	
		"elem_shape_func1 >= 0 && elem_shape_func1 < get_num_nodes() && elem_shape_func2 >= 0 && elem_shape_func2 < get_num_nodes()");
	return axby_m_aybx * reference_element_integral_zero[elem_shape_func_1][elem_shape_func_2];		
}
				
double Element2DRect::get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const {
	
	unsigned int obj_index = elem_shape_func_2 
	                       + 4 * elem_shape_func_1
	                       + 16 * nodal_data_point;

	switch(obj_index) {
        case 21: // Ncc2Nii2Njj2
        case 42: // Ncc3Nii3Njj3
        case 63: // Ncc4Nii4Njj4
        case 0: // Ncc1Nii1Njj1
            return ((1.000000e+00) / (4.000000e+00)) * jacobi_det;
        case 3: // Ncc1Nii1Njj4
        case 4: // Ncc1Nii2Njj1
        case 5: // Ncc1Nii2Njj2
        case 12: // Ncc1Nii4Njj1
        case 15: // Ncc1Nii4Njj4
        case 16: // Ncc2Nii1Njj1
        case 17: // Ncc2Nii1Njj2
        case 20: // Ncc2Nii2Njj1
        case 22: // Ncc2Nii2Njj3
        case 25: // Ncc2Nii3Njj2
        case 26: // Ncc2Nii3Njj3
        case 37: // Ncc3Nii2Njj2
        case 38: // Ncc3Nii2Njj3
        case 41: // Ncc3Nii3Njj2
        case 43: // Ncc3Nii3Njj4
        case 46: // Ncc3Nii4Njj3
        case 47: // Ncc3Nii4Njj4
        case 48: // Ncc4Nii1Njj1
        case 51: // Ncc4Nii1Njj4
        case 58: // Ncc4Nii3Njj3
        case 59: // Ncc4Nii3Njj4
        case 60: // Ncc4Nii4Njj1
        case 62: // Ncc4Nii4Njj3
        case 1: // Ncc1Nii1Njj2
            return ((1.000000e+00) / (1.200000e+01)) * jacobi_det;
        case 6: // Ncc1Nii2Njj3
        case 7: // Ncc1Nii2Njj4
        case 8: // Ncc1Nii3Njj1
        case 9: // Ncc1Nii3Njj2
        case 10: // Ncc1Nii3Njj3
        case 11: // Ncc1Nii3Njj4
        case 13: // Ncc1Nii4Njj2
        case 14: // Ncc1Nii4Njj3
        case 18: // Ncc2Nii1Njj3
        case 19: // Ncc2Nii1Njj4
        case 23: // Ncc2Nii2Njj4
        case 24: // Ncc2Nii3Njj1
        case 27: // Ncc2Nii3Njj4
        case 28: // Ncc2Nii4Njj1
        case 29: // Ncc2Nii4Njj2
        case 30: // Ncc2Nii4Njj3
        case 31: // Ncc2Nii4Njj4
        case 32: // Ncc3Nii1Njj1
        case 33: // Ncc3Nii1Njj2
        case 34: // Ncc3Nii1Njj3
        case 35: // Ncc3Nii1Njj4
        case 36: // Ncc3Nii2Njj1
        case 39: // Ncc3Nii2Njj4
        case 40: // Ncc3Nii3Njj1
        case 44: // Ncc3Nii4Njj1
        case 45: // Ncc3Nii4Njj2
        case 49: // Ncc4Nii1Njj2
        case 50: // Ncc4Nii1Njj3
        case 52: // Ncc4Nii2Njj1
        case 53: // Ncc4Nii2Njj2
        case 54: // Ncc4Nii2Njj3
        case 55: // Ncc4Nii2Njj4
        case 56: // Ncc4Nii3Njj1
        case 57: // Ncc4Nii3Njj2
        case 61: // Ncc4Nii4Njj2
        case 2: // Ncc1Nii1Njj3
            return ((1.000000e+00) / (3.600000e+01)) * jacobi_det;
		default:
			TDKP_GENERAL_EXCEPTION("invalid element / node data point index");		
	};	                       
	                       	
}
				
double Element2DRect::get_single_integral_1st_order(short diffop, short elem_shape_func_1) const {
		
	
	TDKP_BOUNDS_ASSERT(elem_shape_func_1 >= 0 && elem_shape_func_1 < 4, "elem_shape_func_1 >= 0 && elem_shape_func_1 < 4 ");	
	if(diffop == 0) {
		return 0.5 * (  ay * this->reference_single_integral_first[diffop][elem_shape_func_1][0] 
	    	          + by * this->reference_single_integral_first[diffop][elem_shape_func_1][1]); 
	} else {
		return 0.5 * (  ax * this->reference_single_integral_first[diffop][elem_shape_func_1][0] 
	    	          + bx * this->reference_single_integral_first[diffop][elem_shape_func_1][1]); 		
	}	    	          
		     		   
}
double Element2DRect::get_single_integral_0th_order(short elem_shape_func_1) const {
	return axby_m_aybx / 4.0;
}
					
 	
/** function to evaluate the form function inside the reference element */	
double Element2DRect::evaluate_form_function(
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	const double& x = local_reference_element_coords[0];
	const double& y = local_reference_element_coords[1];
	switch(elem_shape_func) {
		case 0:
			return 0.25 * (1.0 - x) * (1.0 - y);
		case 1:
			return 0.25 * (1.0 + x) * (1.0 - y);
		case 2:
			return 0.25 * (1.0 + x) * (1.0 + y);
		case 3:
			return 0.25 * (1.0 - x) * (1.0 + y);
		default:
			TDKP_GENERAL_EXCEPTION("unknown shape function index");			
	}	
}

double Element2DRect::evaluate_form_function_derivative(
	short diffop, 
	short elem_shape_func, 
	const double* local_reference_element_coords
) const {
	TDKP_BOUNDS_ASSERT(diffop >= 0 && diffop < 2, "diffop >= 0 && diffop < 2");
	const double& x = local_reference_element_coords[0];
	const double& y = local_reference_element_coords[1];
	if(diffop == 0) {
		switch(elem_shape_func) {		
			case 0:
				return (ay - ay * x + by * (-1.0 + y)) / (2.0 * axby_m_aybx); 
			case 1:
				return (ay + by + ay * x - by * y) / (2.0 * axby_m_aybx);
			case 2:
				return ((-ay) * (1.0 + x) + by * (1.0 + y)) / (2.0 * axby_m_aybx);	
			case 3:
				return (ay + by - ay * x + by * y) / (- 2.0 * axby_m_aybx);
			default:
				TDKP_GENERAL_EXCEPTION("unknown shape function index");			
		}		
	} else {		
		switch(elem_shape_func) {		
			case 0:
				return (bx + ax * (-1.0 + x) - bx * y)/  (2.0 * axby_m_aybx);
			case 1:
				return (ax + bx + ax * x - bx * y) / (- 2.0 * axby_m_aybx);	
			case 2:
				return (ax * (1.0 + x) - bx * (1.0 + y))/ (2.0 * axby_m_aybx);		
			case 3:
				return (ax + bx - ax * x + bx * y) / (2.0 * axby_m_aybx);
			default:
				TDKP_GENERAL_EXCEPTION("unknown shape function index");			
		}
	}			
}

 
double Element2DRect::evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const {
	double global[] = {global_x, global_y};
	double local[2];
	if(global2local(global,local)) {
		return evaluate_form_function(elem_shape_func, local);
	} else {
		TDKP_GENERAL_EXCEPTION("global coordinates are not inside element");
	}
}			

void Element2DRect::get_node_local_coords(unsigned short lid, vector<double>& local_coords) const {
	local_coords.resize(2);
	switch(lid) {
		case 0:
			local_coords[0] = -1.0;
			local_coords[1] = -1.0; 
			break;			
		case 1:
			local_coords[0] =  1.0;
			local_coords[1] = -1.0; 
			break;	
		case 2:
			local_coords[0] =  1.0;
			local_coords[1] =  1.0; 
			break;
		case 3:
			local_coords[0] = -1.0;
			local_coords[1] =  1.0; 
			break;				
		default:
			TDKP_GENERAL_EXCEPTION("invalid node index");		
	}
}
																
} // end of namespace
