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

#ifndef ARPACKSIPARDISOSOLVER_H_
#define ARPACKSIPARDISOSOLVER_H_


#include "tdkp/common/all.h"
#include "tdkp/solvers/EigenSolver.h"
#include "tdkp/solvers/LinearSolver.h"
#include "tdkp/main/CSRSparseMatrixInterface.h"

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING ARPACK
// -----------------------------------------------
#ifdef EIGNSOLV_INCLUDE_ARPACK

namespace tdkp {

class ArpackSIPardisoSolver : public EigenSolver<cplx, double, cplx> {
public:
	ArpackSIPardisoSolver(unsigned int size, unsigned int block_size);
	virtual ~ArpackSIPardisoSolver();		
	virtual bool find_eigenvectors() throw(Exception*);	
	virtual cplx eigenvalue(int nn) const throw(Exception*);
	virtual const cplx& eigenvector(int nn, int vv) const throw(Exception*);	
	virtual void orthogonalize_solutions();	
	
	virtual SparseMatrixInterface<cplx>&   get_stiff() { return solver->get_matrix(); }
	virtual SparseMatrixInterface<double>& get_overlap() { return *this->overlap; }	
	
	const cplx* get_evals_data() const { return evals; }
	const cplx* get_evecs_data() const { return evecs; }
	
protected:	

private:
	
	cplx* evecs;
	cplx* evals;
	
	LinearSolver<cplx>* solver;
	CSRMatrix<double>* overlap;

	
};

/** wrapper that uses the complex version to calculate real problems 
 *
 * decorator pattern 
 */
class ArpackSIRealWrapper : public EigenSolver<double, double, cplx> {
public:
	ArpackSIRealWrapper(unsigned int size, unsigned int block_size);
	~ArpackSIRealWrapper();
	virtual void assign(int nev, EigenProblemType type);
	virtual void set_ordering(Ordering order_);
	virtual Ordering get_ordering() const;
	virtual bool find_eigenvectors() throw(Exception*);
	virtual int converged_eigenvalues() const;
	virtual cplx eigenvalue(int nn) const throw(Exception*);
	virtual const cplx& eigenvector(int nn, int vv) const throw(Exception*);
		
	// -------------------------------------------------------
	// factory and abstraction
	// -------------------------------------------------------	
	virtual SparseMatrixInterface<double>& get_stiff() { return matrix_wrapper; }
	virtual SparseMatrixInterface<double>& get_overlap() { return solver.get_overlap(); }		
		
protected:	
	ArpackSIPardisoSolver       solver;
	D2CCSRSparseMatrixInterface matrix_wrapper;
	
};

} // end of namespace

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING ARPACK
// -----------------------------------------------
#endif

#endif /*ARPACKSIPARDISOSOLVER_H_*/
