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

#ifndef NEGFHAMILTONANASSEMBLER_H_
#define NEGFHAMILTONANASSEMBLER_H_

#include "tdkp/common/all.h"
#include "tdkp/probdefs/KPBase1D2D.h"
#include "tdkp/probdefs/NEGFEffectiveMass.h"
#include "tdkp/main/FEMSolverGEVP.h"


namespace tdkp {


/** NEGF hamiltonian assembler for sebis NEGF solver */
template<class KPProblem, class FEMSolver>
class NEGFHamiltonianAssembler : public KPProblem {	
	public:
		NEGFHamiltonianAssembler(const Geometry& geometry_, MaterialDatabase& material_database_);
		virtual ~NEGFHamiltonianAssembler();			
		void assemble(const double& k_transversal);
		virtual string get_unique_identifier() const { return string("NEGFHamiltonianAssembler"); }
		
		const typename FEMSolver::SMatrix& get_lhs_matrix() const;
		const typename FEMSolver::OMatrix& get_rhs_matrix() const;
					
	private:				
		FEMSolver solver;
		bool      solver_ready;

};


/** return assembled lhs matrix  */
template<class KPProblem, class FEMSolver>
const typename FEMSolver::SMatrix& NEGFHamiltonianAssembler<KPProblem, FEMSolver>::get_lhs_matrix() const {
	if(solver_ready) {
		return solver.get_lhs_matrix();		
	} else {
		TDKP_GENERAL_EXCEPTION("fem solver is not ready, so no lhs matrix is available");	
	}	
}

/** return assembled rhs matrix  */
template<class KPProblem, class FEMSolver>
const typename FEMSolver::OMatrix& NEGFHamiltonianAssembler<KPProblem, FEMSolver>::get_rhs_matrix() const {
	if(solver_ready) {
		return solver.get_rhs_matrix();		
	} else {
		TDKP_GENERAL_EXCEPTION("fem solver is not ready, so no rhs matrix is available");	
	}	
}
	
/** assembler construction */
template<class KPProblem, class FEMSolver>
NEGFHamiltonianAssembler<KPProblem, FEMSolver>::NEGFHamiltonianAssembler(const Geometry& geometry_, MaterialDatabase& material_database_)
: KPProblem(geometry_, material_database_),
  solver(geometry_, *this),
  solver_ready(false)
{
		   
}


/** assembler destruction */
template<class KPProblem, class FEMSolver>
NEGFHamiltonianAssembler<KPProblem, FEMSolver>::~NEGFHamiltonianAssembler() {

}


template<>
void NEGFHamiltonianAssembler<NEGFEffectiveMass, FEMSolverGEVP<cplx, complex<double>, double> >::assemble(const double& k_transversal_);

/** assemble matrix */
template<class KPProblem, class FEMSolver>
void NEGFHamiltonianAssembler<KPProblem, FEMSolver>::assemble(const double& k_transversal_) {

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");
	}
	
	this->k_idx_current = 0; // nasty, must be 0 before we get the right energy shift (in subsequent calc., the k idx may still be set > 0)
	
	if(!this->solver_ready) {
		solver.create_matrix_structures();
		this->solver_ready = true;
	} else {
		solver.reset_matrix_to_zero();	
	}
	this->prepare();
	// -------------------------------------------------------------
	// tell the matrices to shut up
	// -------------------------------------------------------------
	for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {	
		this->kp_matrices[ii]->set_output_surpression(true);		
	}
	this->k_transversal = k_transversal_;
	solver.assemble_system();
			
}
	
	
} // end of namespace


#endif /*NEGFHAMILTONANASSEMBLER_H_*/
