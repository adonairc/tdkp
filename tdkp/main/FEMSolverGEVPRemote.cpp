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

#include "FEMSolverGEVPRemote.h"

namespace tdkp
{

FEMSolverGEVPRemote::FEMSolverGEVPRemote(const Geometry& geometry_, EigenProblem<cplx,cplx,double>& problem_)
: FEMSolverGEVP<cplx,cplx, double>(geometry_, problem_),
  remote_solver(0)
{
} 

/** remote solver variable killed via eigensolver variable */
FEMSolverGEVPRemote::~FEMSolverGEVPRemote() {} 
 		
void FEMSolverGEVPRemote::solve_system(int num_solutions) {
	TDKP_GENERAL_EXCEPTION("you have to create remote solver objects. this is just an assembler");	
}

/** creates a solver controller with the current stiff and overlap matrix */
RemoteGEVPSolverController* FEMSolverGEVPRemote::get_remote_solve_controller(int num_subbands, EigenProblemType prbl_type) const {
	TDKP_ASSERT(this->remote_solver != 0, "");
	return this->remote_solver->get_controller(this->order, num_subbands);		
}

/** creates remote solver object */
EigenSolver<cplx,double,cplx>* FEMSolverGEVPRemote::create_eigensolver_object(unsigned int matrix_size, unsigned int block_size) {
	this->remote_solver = new RemoteGEVPSolver(matrix_size, block_size);
	return this->remote_solver;	
}

}
