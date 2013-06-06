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

#ifndef _ELEMENT_3DTR_H
#define _ELEMENT_3DTR_H


#include "tdkp/geometry/Region.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/common/Exception.h"
#include "tdkp/main/Fields.h"


using namespace std;

namespace tdkp {

/** base tetrahedron class
 * 
 * vertex, edge and face number is as follows:
 * N0 = 0 0 0
 * N1 = 1 0 0
 * N2 = 0 1 0
 * N3 = 0 0 1
 * 
 * E0 = N0 N1
 * E1 = N1 N2
 * E2 = N2 N0
 * E3 = N0 N3
 * E4 = N1 N3
 * E5 = N2 N3
 * 
 * F0 = N0 N1 N2
 * F1 = N0 N1 N3
 * F2 = N1 N2 N3
 * F3 = N2 N0 N3
 */
class Element3DTetrahedronBase : public Element {
public:	
	// ----------------------------------------------
	// setup functions		
	// ----------------------------------------------
	virtual ElementShape get_shape() const { return tetrahedron; }	
	virtual unsigned int get_num_corners() const { return 4; }
	virtual unsigned int get_num_edges() const { return 6; }
	virtual unsigned int get_num_faces() const { return 4; }
	virtual void get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const;
	virtual void get_face(unsigned int face_idx, vector<unsigned int>& vertex_indices) const;
	virtual void get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const;	
	bool 	     verify() const;
	void         prepare(); 					
	bool 		 global2local(const double global[3], double local[3]) const;	
	virtual bool inside_element(const double& global_x, const double& global_y, const double& global_z) const ;			
						
protected:
	Element3DTetrahedronBase(unsigned int nnodes);
	double inverted_jacobi_matrix[9];	// inverted jacobian col 1 holds d (u,v,w) / dx etc. 
	double jacobi_matrix[9];
	double* const jm; // shortcut to jacobi_matrix
};
	
/** element 3DTR - 3 dimensional tetraeder, region based
 * 
 * three dimensional tetraeder 
 *  */
class Element3DTR : public Element3DTetrahedronBase {
							
public:				
	Element3DTR(unsigned int index);
	~Element3DTR();
	
	virtual unsigned int get_element_order() const { return 1; }	
				
	// ----------------------------------------------
	// evaluation functions
	// ----------------------------------------------
	double  get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const;
	double  get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shape_func_2) const;
	double  get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const;

	double  get_single_integral_1st_order(short diffop, short elem_shape_func_1) const;
	double  get_single_integral_0th_order(short elem_shape_func_1) const;
	
	double  get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const;
							
	double evaluate_form_function(short elem_shape_func, const double* local_reference_element_coords) const;
	double evaluate_form_function_derivative(short diffop, short elem_shape_func, const double* local_reference_element_coords) const;					
	
	bool   get_contribution(const double pcoords[3], double contrib[4]) const;
	virtual void get_node_local_coords(unsigned short lid, vector<double>& local_coords) const;				
	virtual double evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const;
	virtual unsigned int get_element_unique_type_key() const { return 8; }																		
private:		
					
	static const double reference_elem_integral_zero[4][4]; // zero order element integrals
	static const double reference_elem_gradients[4][3];     // derivatives of linear element fcts are constants ...
																
};
	
} // end namespace 




#endif
