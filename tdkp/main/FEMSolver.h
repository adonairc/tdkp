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

#ifndef FEMSOLVER_H_
#define FEMSOLVER_H_

#include "tdkp/common/all.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/probdefs/ProblemDefinition.h"


namespace tdkp {
	
/** base fem solver class 
 * 
 * is called from problem class definition. creates matrix structures depending on geometry and the problems
 * sparsity pattern, assembles the systems and solves the resulting equation. then returns the problem class
 * the obtained solution(s). 
 */
class FEMSolver {
public:
	FEMSolver(const Geometry& geometry_) : geometry(geometry_) {}
	virtual ~FEMSolver() {};

	// ---------------------------------------------------
	// simulation setup functions 
	// ---------------------------------------------------
	virtual bool verify_setup() const = 0;				
		
	// ---------------------------------------------------
	// assembling 
	// ---------------------------------------------------
	virtual void create_matrix_structures()      = 0;
	virtual void reset_matrix_to_zero()          = 0; 
	virtual void assemble_system()               = 0;
	virtual void solve_system(int num_solutions) = 0;
																			
protected:
	const Geometry& geometry;
	
};


} // end of namespace

#endif /*FEMSOLVER_H_*/
