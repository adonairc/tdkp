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


#include "tdkp/solvers/EigenSolver.h"
#include "tdkp/solvers/LapackBandSolver.h"
#ifdef EIGNSOLV_INCLUDE_ARPACK
#include "tdkp/solvers/ArpackSIPardisoSolver.h"
#endif
#ifdef EIGNSOLV_INCLUDE_JDQZ
#include "tdkp/solvers/JDQZSolver.h"
#endif
#ifdef LINSOLV_INCLUDE_PARDISO
#include "tdkp/solvers/Pardiso.h"
#endif
#ifdef LINSOLV_INCLUDE_UMFPACK
#include "tdkp/solvers/UmfpackSolver.h"
#endif
#ifdef LINSOLV_INCLUDE_PETSC
#include "tdkp/solvers/PetscSolver.h"
#endif


namespace tdkp {


template<>
EigenSolver<cplx,double,cplx>* EigenSolver<cplx,double,cplx>::factory(unsigned int size_, unsigned int block_size_) {

	if(Configuration::get_instance()->get("desired_eigenvalue_solver") == 3.0) {
    	return new LapackComplexBandEigenSolver(size_, block_size_);
    } else if(Configuration::get_instance()->get("desired_eigenvalue_solver") == 4.0) {
#ifdef EIGNSOLV_INCLUDE_ARPACK    	
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "inappropriate request for eigensolver (remote mode requested in factory). falling back to standard arpack");
    	return new ArpackSIPardisoSolver(size_, block_size_);
#else
		TDKP_GENERAL_EXCEPTION("Eigensolver::factory: ARPACK Eigensolver was disabled at compile time");
#endif    	
    } else if(Configuration::get_instance()->get("desired_eigenvalue_solver") == 5.0) {
#ifdef EIGNSOLV_INCLUDE_JDQZ    			
    	return new JDQZSolver<double>(size_, block_size_);
#else
		TDKP_GENERAL_EXCEPTION("Eigensolver::factory: JDQZ Eigensolver was disabled at compile time");
#endif    	    	
    } else {
#ifdef EIGNSOLV_INCLUDE_ARPACK    	
	   	return new ArpackSIPardisoSolver(size_, block_size_);
#else
		TDKP_GENERAL_EXCEPTION("Eigensolver::factory: ARPACK Eigensolver was disabled at compile time");
#endif	   	
    }
}

// currently we treat the real problem as complex ... 
template<>
EigenSolver<double,double,cplx>* EigenSolver<double,double,cplx>::factory(unsigned int size_, unsigned int block_size_) {
#ifdef EIGNSOLV_INCLUDE_ARPACK		
	return new ArpackSIRealWrapper(size_, block_size_);
#else
		TDKP_GENERAL_EXCEPTION("Eigensolver::factory: ARPACK Eigensolver was disabled at compile time");
#endif	
}  

// ----------------------------------------------------
// can't template these functions as Pardiso/Ils are 
// a derived class of linear solver ...
// probably i should move the LinearSolverFactory to its own class and file 
// ----------------------------------------------------
template<>  
LinearSolver<cplx>* LinearSolver<cplx>::factory(unsigned int size_) {
	Configuration* conf = Configuration::get_instance();
	double solver_id = conf->get("desired_linear_equation_solver");
	TDKP_ASSERT(size_ > 0, "size_ > 0");
	if(solver_id == 2.0) {
#ifdef LINSOLV_INCLUDE_UMFPACK		
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory UmfpackSolverComplex is used as linear equation solver");
		return new UmfpackSolverComplex(size_);
#else
		TDKP_GENERAL_EXCEPTION("LinearSolver::factory: Umfpack linear solver was disabled at compile time");
#endif				
	} else if(solver_id == 3.0) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory LapackBandSolverComplex is used as linear equation solver");
		return new LapackBandSolverComplex(size_);
	} else if(solver_id == 4.0) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory DummyComplex is used as linear equation solver");
		return new DummyComplex(size_);
	} else if(solver_id == 5.0) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory LapackFullSolver is used as linear equation solver");
		return new LapackFullSolver(size_);
	} else if(solver_id == 8.0) {
#ifdef LINSOLV_INCLUDE_PARDISO		
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory: Pardiso is used as linear equation solver");
		return new Pardiso<cplx>(size_);
#else
		TDKP_GENERAL_EXCEPTION("LinearSolver::factory: Pardiso linear solver was disabled at compile time");
#endif		
	} else if(solver_id == 9.0) {
#ifdef LINSOLV_INCLUDE_PETSC		
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory: Petsc is used as linear equation solver");
		return new PetscSolver(size_);
#else
		TDKP_GENERAL_EXCEPTION("LinearSolver::factory: Petsc linear solver was disabled at compile time");
#endif									
	} else {	
		// -------------------------------------------------
		// normal approach release: we assume to have pardiso and umfpack ...
		// -------------------------------------------------
		unsigned int umfpack_limit = static_cast<unsigned int>(Configuration::get_instance()->get("autosolver_umfpack_upper_size_limit"));
		if(size_ < umfpack_limit) {
#ifdef LINSOLV_INCLUDE_UMFPACK
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory: UmfpackSolverComplex is used as linear equation solver");		
			return new UmfpackSolverComplex(size_);
#else
			TDKP_GENERAL_EXCEPTION("LinearSolver::factory: UMFPACK was disabled at compile time");
#endif
		} else {
#ifdef LINSOLV_INCLUDE_PARDISO
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "LinearSolver::factory: Pardiso is used as linear equation solver");
			return new Pardiso<cplx>(size_);
#else
			TDKP_GENERAL_EXCEPTION("LinearSolver::factory: Pardiso linear solver was disabled at compile time");
#endif
		}
	}		
}

template<>  
LinearSolver<double>* LinearSolver<double>::factory(unsigned int size_) {

	Configuration* conf = Configuration::get_instance();
	double solver_id = conf->get("desired_linear_equation_solver");

#ifdef LINSOLV_INCLUDE_PARDISO		
	return new Pardiso<double>(size_);
#else
	TDKP_GENERAL_EXCEPTION("LinearSolver::factory: Pardiso linear solver was disabled at compile time");
#endif

}

template<>
EigenSolver<cplx,cplx,cplx>* EigenSolver<cplx,cplx,cplx>::factory(unsigned int size_, unsigned int block_size_) {
	if(Configuration::get_instance()->get("desired_eigenvalue_solver_pml") == 2.0) {
		return new LapackZGGEVXEigenSolver(size_, block_size_);
    } else if(Configuration::get_instance()->get("desired_eigenvalue_solver_pml") == 3.0) {
#ifdef EIGNSOLV_INCLUDE_JDQZ 
    	return new JDQZSolver<cplx>(size_, block_size_);
#else
		TDKP_GENERAL_EXCEPTION("Eigensolver::factory: JDQZ Eigensolver was disabled at compile time");
#endif    	 		
	} else {		
		return new LapackZGGEVEigenSolver(size_, block_size_);
	}
}  


}
