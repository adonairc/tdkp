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

#ifndef ELEMENT2DRECT_H_
#define ELEMENT2DRECT_H_

#include "tdkp/geometry/Element.h"

namespace tdkp {
	
/** basic class for rectangular 2d elements */	
class Element2DRectBase : public Element {

public:	
	virtual ~Element2DRectBase() {};	

	// ----------------------------------------------
	// setup functions		
	// ----------------------------------------------		
	bool verify() const;
	void prepare(); 
	virtual ElementShape get_shape() const { return rectangle; }	
	virtual unsigned int get_num_corners() const { return 4; }
	virtual unsigned int get_num_edges() const { return 4; }
	virtual unsigned int get_num_faces() const { return 0; }
	virtual void get_edge(unsigned int edge_idx, vector<unsigned int>& vertex_indices) const;
	virtual void get_element_boundary_vertices(unsigned int idx, vector<unsigned int>& vertex_indices) const;
				
	// ---------------------------------------------
	// global numerical testing of formfunctions
	// ---------------------------------------------
	bool 		   global2local(const double global[2], double local[2]) const;	
	virtual bool   inside_element(const double& global_x, const double& global_y, const double& global_z) const;
	
protected:
	Element2DRectBase(unsigned int dimension, unsigned int nnodes, unsigned int nboundaries);
	double inverted_jacobi_matrix[4];	// inverted jacobian col 1 holds d (u,v,w) / dx etc.	
	/** parameters of the rectangle */
	double a[2];
	double b[2];
	/** parameters for the integrals (i just copy expressions from mathematica, 
	 *  where i named a[0] = ax etc .. )*/
	double ax,ay,bx,by;
	double axby_m_aybx;	
	/** vertex offset -> every vertex_offset node in nodes is a vertex  (1 for linear) */ 
	short vertex_offset;

};	

class Element2DRect : public Element2DRectBase {
public:

	Element2DRect(unsigned int index);
	virtual ~Element2DRect();

	virtual void get_node_local_coords(unsigned short lid, vector<double>& local_coords) const;
	virtual unsigned int get_element_order() const { return 1; }
	virtual unsigned int get_element_unique_type_key() const { return 3; }			
	// ----------------------------------------------
	// evaluation functions
	// ----------------------------------------------
	double  get_element_integral_2nd_order(short diffop_1, short elem_shape_func_1, short diffop_2, short elem_shape_func_2) const;
	double  get_element_integral_1st_order(short diffop, short elem_shape_func_1, short elem_shalpe_func_2) const;
	double  get_element_integral_0th_order(short elem_shape_func_1, short elem_shape_func_2) const;

	double  get_single_integral_1st_order(short diffop, short elem_shape_func_1) const;
	double  get_single_integral_0th_order(short elem_shape_func_1) const;
	
	double  get_element_integral_0th_order_nodal_data(short nodal_data_point, short elem_shape_func_1, short elem_shape_func_2) const;						
									
	/** function to evaluate the form function inside the reference element */	
	double evaluate_form_function(short elem_shape_func, const double* local_reference_element_coords) const;
	/** function to evaluate the form function derivative inside the reference element */
	double evaluate_form_function_derivative(short diffop, short elem_shape_func, const double* local_reference_element_coords) const;				
				
	// ---------------------------------------------
	// global numerical testing of formfunctions
	// --------------------------------------------- 
	virtual double evaluate_form_function_global(short elem_shape_func, const double& global_x, const double& global_y, const double& global_z) const;					
															
private:
	static double reference_element_integral_zero[4][4];
	static double reference_single_integral_first[2][4][2];
	
					
};


}

#endif /*ELEMENT2DRECT_H_*/
