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

#ifndef _PROBLEMDEFINITION_H_
#define _PROBLEMDEFINITION_H_

#include <math.h>
#include <map>

#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/main/EigenSolution.h"
#include "tdkp/kpmatrices/KPMatrixBase.h"

namespace tdkp {

/** template problem definition -> physics go here. 
 * 
 * @param T    result value type coming from eigensolver
 * @param TMat value type for matrices
 * @param TRHS value type for right hand side
 */  
template<class T, class TMAT, class TRHS>
class ProblemDefinition {
public:

	ProblemDefinition(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~ProblemDefinition();

	// -------------------------------------------------------
	// init & setup & actions
	// -------------------------------------------------------	
	const Geometry& get_geometry() const { return this->geometry; }
	const MaterialDatabase& get_material_database() const { return material_db; }
	MaterialDatabase& get_material_database() { return material_db; }
	
	virtual void solve(int num_solutions) = 0;
		
	// ------------------------------------------------------
	// result access 
	// -------------------------------------------------------	
	virtual  void   display_solution_info() const = 0;
			
	// -------------------------------------------------------
	// (intercommunication with FEMSolver during solver)
	// -------------------------------------------------------
	virtual void 		prepare();
	bool                is_ready() const;
	int 		 		get_num_equations_per_node() const { return num_equations_per_node; }	
	virtual const int*  get_node_sparsity_pattern(int &num) const = 0;
	virtual void 		calculate_element_matrices(const Element* elem, TMAT* lhs, TRHS *rhs, int* node_internal_indices, int &n) const = 0;   	
	virtual string	    get_unique_identifier() const = 0;

protected:	
	/** geometry object */
	const Geometry& geometry;
	/** material database from where material properties are read */
	MaterialDatabase& material_db;
	/** number of pde's defined */
	int	num_equations_per_node;
	/** boolean to mark whether object is ready for calculation */
	bool ready;
	
};

template<class T, class TMAT, class TRHS>
ProblemDefinition<T,TMAT,TRHS>::ProblemDefinition(const Geometry& geometry_, MaterialDatabase& material_database_) 
: geometry(geometry_),
  material_db(material_database_),
  num_equations_per_node(1),
  ready(false)
{
	// -------------------------------------
	// check if geometry knows all materials
	// -------------------------------------
	typename Geometry::region_const_iterator it;
	for(it = this->geometry.regions_begin(); it != this->geometry.regions_end(); it++) {
		if((*it)->enabled()) {
			TDKP_ASSERT((*it)->material_set(), "material is not set to geometry object!");
		}		
	}
}

template<class T, class TMAT, class TRHS>
ProblemDefinition<T,TMAT,TRHS>::~ProblemDefinition() {
	
}


/** prepare is called by tdkp_mail once before assembly
 */ 
template<class T, class TMAT, class TRHS>
void ProblemDefinition<T,TMAT,TRHS>::prepare() {
	
}

template<class T, class TMAT, class TRHS>
bool ProblemDefinition<T,TMAT,TRHS>::is_ready() const {
	return this->ready;	
}

}
#endif /*PROBLEMDEFINITION_H_*/
