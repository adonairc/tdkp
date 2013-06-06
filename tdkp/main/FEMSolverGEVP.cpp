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

#include "tdkp/main/FEMSolverGEVP.h"
#include "tdkp/common/Configuration.h"


using namespace std;

namespace tdkp {



/** solve complex generalized eigenvalue problem 
 */
template<> void FEMSolverGEVP<cplx, cplx, double>::solve_system(int num_solutions) {
	
	TDKP_ASSERT(this->stiff->get_size() > 0, "GEVP does not seem to be assembled yet!");
	eigensolver->set_ordering(this->order);
	TDKP_ASSERT(num_solutions > 0,"number of solutions must be greater than 0");
	cplx* result = new cplx[this->stiff->get_size()];
	TDKP_ASSERT(result > 0, "pointer: result > 0");
	TimeMeasurements::get_instance().start("searching eigenvalues");	
	eigensolver->assign(num_solutions, this->problem.get_problem_type());
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverGEVP: searching eigenvalues");		
	eigensolver->find_eigenvectors();
	TimeMeasurements::get_instance().stop("searching eigenvalues");	
	TimeMeasurements::get_instance().start("result postprocessing");	
	for(int ii = 0; ii < eigensolver->converged_eigenvalues(); ii++) {
		for(int jj = 0; jj < (signed)this->stiff->get_size(); jj++) {
			result[jj] = eigensolver->eigenvector(ii,jj); 
		}				
    	problem.add_solution(eigensolver->eigenvalue(ii), result, this->stiff->get_size());		     
	}
	TimeMeasurements::get_instance().stop("result postprocessing");	
	delete[] result;
}

/** solve complex eigenvalue problem
 */
template<> void FEMSolverGEVP<cplx, cplx, cplx>::solve_system(int num_solutions) {
	
	TDKP_ASSERT(this->stiff->get_size() > 0, "GEVP does not seem to be assembled yet!");
	eigensolver->set_ordering(this->order);
	TDKP_ASSERT(num_solutions > 0,"number of solutions must be greater than 0");
	cplx* result = new cplx[this->stiff->get_size()];
	TDKP_ASSERT(result > 0, "pointer: result > 0");
	TimeMeasurements::get_instance().start("searching eigenvalues");	
	eigensolver->assign(num_solutions, this->problem.get_problem_type());
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverGEVP: searching eigenvalues");		
	eigensolver->find_eigenvectors();
	TimeMeasurements::get_instance().stop("searching eigenvalues");	
	TimeMeasurements::get_instance().start("result postprocessing");	
	for(int ii = 0; ii < eigensolver->converged_eigenvalues(); ii++) {
		for(int jj = 0; jj < (signed)this->stiff->get_size(); jj++) {
			result[jj] = eigensolver->eigenvector(ii,jj); 
		}				
    	problem.add_solution(eigensolver->eigenvalue(ii), result, this->stiff->get_size());		     
	}
	TimeMeasurements::get_instance().stop("result postprocessing");	
	delete[] result;
}

/** solve real symmetric generalized eigenvalue problem 
 */
template<> void FEMSolverGEVP<cplx, double, double>::solve_system(int num_solutions) {
	
	TDKP_ASSERT(num_solutions > 0,"number of solutions must be greater than 0");
	
	// ----------------------------------------------------------------------------
	// all solvers here are complex (and actually it's quite unlike that somebody
	// will use the eff mass stuff ... ), so create a complex problem out of it ..
	// ----------------------------------------------------------------------------
	eigensolver->set_ordering(this->order);	
	cplx* result = new cplx[this->stiff->get_size()];
	TimeMeasurements::get_instance().start("searching eigenvalues");	
	eigensolver->assign(num_solutions, this->problem.get_problem_type());
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverGEVP: searching eigenvalues");
	eigensolver->find_eigenvectors();
	TimeMeasurements::get_instance().stop("searching eigenvalues");
	TimeMeasurements::get_instance().start("result postprocessing");
	// store eigenvalues 
	for(int ii = 0; ii < eigensolver->converged_eigenvalues(); ii++) {   	
		for(int jj = 0; jj < (signed)this->stiff->get_size(); jj++) {
			result[jj] = eigensolver->eigenvector(ii,jj); 
		}		
    	problem.add_solution(eigensolver->eigenvalue(ii), result, this->stiff->get_size());		
	}
	TimeMeasurements::get_instance().stop("result postprocessing");		
	delete[] result;

}



} // end of namespace
