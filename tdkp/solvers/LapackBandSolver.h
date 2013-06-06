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

#ifndef LAPACKBANDSOLVER_H_
#define LAPACKBANDSOLVER_H_

#include "tdkp/common/all.h"
#include "tdkp/main/SparseMatrixInterface.h"
#include "tdkp/solvers/LinearSolver.h"
#include "tdkp/solvers/EigenSolver.h"

namespace tdkp {


// munteres ratespiel lead to the following band storages format for banded lapack matrices
#define FBAND_STORAGE(bd,ii,jj)  (((bd + 1) + (ii - jj)) + (jj) * (2 * bd + 1))

/** banded matrix, usable for 1d calculations */	
template<class T>
class BandMatrix : public SparseMatrixInterface<T> {
public:
	BandMatrix(unsigned int size);
	~BandMatrix();
	virtual bool			property_is_set(SparseMatrixProperties flag) const;
	virtual void			announce(int ii, int jj);
	virtual void 			set_structure();
	virtual void 			clear_but_keep_structure();
	virtual void 			reset();
	virtual void 			set_block_size(unsigned int block_size_);
	virtual unsigned int	get_block_size() const;
	virtual void     		set(int ii, int jj, const T& val);
	virtual void 			set_row_and_column_to_zero(int ii); 
	virtual void     		add(int ii, int jj, const T& val);
	virtual T 				get(int ii, int jj) const;	
	virtual unsigned int 	get_size() const { return size; }				
	virtual void mult_vec(const T* in, T* out) const;			
	void    set_property(SparseMatrixProperties property);
	
	T*	get_data(int& data_length_) { data_length_ = data_length; return data; }
	unsigned int get_band_width() const { return band_width; }
	
protected:
	bool is_symmetric() const { return symmetric; }	
	
private:
	T*       data;
	unsigned int  size;
	unsigned int  band_width;
	unsigned int  data_length;
	unsigned int  block_size; // no effect	
	bool     structure_is_set;
	bool     symmetric;
	
	inline int get_idx(int ii,int jj) const;
	
};

/** full matrix (so non sparse), used for lapack routines requiring full matrices */
template<class T>
class FullMatrix : public SparseMatrixInterface<T> { 
public:
	FullMatrix(unsigned int size_);
	virtual ~FullMatrix() {};	
	virtual void 		 set_structure();
	virtual void     	 set(int ii, int jj, const T& val);
	virtual void 	     set_row_and_column_to_zero(int ii);
	virtual void     	 add(int ii, int jj, const T& val);
	virtual T 			 get(int ii, int jj) const;				
	const vector<T>&     get_data_vector() const { return mdata; }
	virtual void         mult_vec(const T* in, T* out) const;		
		

	virtual void		 announce(int ii, int jj) {}
	virtual void 		 clear_but_keep_structure();
	virtual void 		 reset();
	virtual void 		 set_block_size(unsigned int block_size_) { block_size = block_size_; }
	virtual unsigned int get_block_size() const { return block_size; }
	virtual unsigned int get_size() const { return size; }
	virtual bool         property_is_set(SparseMatrixProperties flag) const;									
	virtual void         set_property(SparseMatrixProperties property);
	virtual void         save_to_file(const char* filename) const;
	virtual void         perform_symmetry_analysis();		
		
private:
	bool          symmetric;
	vector<T>     mdata;
	unsigned int  size;
	unsigned int  block_size; // no effect		
			
};


/** dummy complex linear solver (inoperable)
 *
 * used on places where we just assemble the matrix but don't
 * really calculate anything. (e.g. calculating matrix elements etc.) 
 */
class DummyComplex : public LinearSolver<cplx> {
public:
	DummyComplex(unsigned int size) : matrix(size, nonsymmetric_matrix) { }
	virtual ~DummyComplex() throw(Exception*) {}
	virtual void prepare() throw(Exception*) {}
	virtual void solve_equation(cplx* res, cplx* rhs, int num = 1) throw(Exception*) { TDKP_GENERAL_EXCEPTION("sorry, im just a dummy!"); }			
	virtual SparseMatrixInterface<cplx>& get_matrix() { return matrix; }
private:
	CSRMatrix<cplx> matrix;
};	

/** lapack linear equation solver for full matrices */
class LapackFullSolver : public LinearSolver<cplx> {
public:
	LapackFullSolver(unsigned int size);
	virtual ~LapackFullSolver() throw(Exception*);
	virtual void prepare() throw(Exception*);
	virtual void solve_equation(cplx* res, cplx* rhs, int num = 1) throw(Exception*);
	virtual SparseMatrixInterface<cplx>& get_matrix() { return matrix; }		
private:
	void factorize_and_solve(cplx* res, cplx* rhs);
	void solve_only(cplx* res, cplx* rhs);

	FullMatrix<cplx> matrix;	
	bool factorized;
	
	vector<cplx>   matrix_data;
	vector<cplx>   af;
	vector<int>    ipiv;

				
};

/** lapack solver for banded hermitian matrices */
class LapackBandSolverComplex : public LinearSolver<cplx> {
public:
	LapackBandSolverComplex(unsigned int size);
	virtual ~LapackBandSolverComplex() throw(Exception*);		
	virtual void prepare() throw(Exception*);
	virtual void solve_equation(cplx* res, cplx* rhs, int num = 1) throw(Exception*);			
	virtual SparseMatrixInterface<cplx>& get_matrix() { return matrix; }
private:
	void release_lapack_space();
	void factorize(cplx* res, cplx* rhs);
	BandMatrix<cplx> matrix;
	bool       factorized;
	
	
	// ---------------------------------------------------
	// variables used by zgbsvx to store LU factorization
	// nonsymmetric case
	// ---------------------------------------------------
	cplx* 	lp_lu_A;
	int   	lp_ldafb;
	cplx* 	lp_afb;
	int* 	lp_ipiv;
	char 	lp_equed; // output variable used in subsequent solving steps
	double* lp_r;
	double* lp_c;
	cplx* 	lp_work;
	double* lp_rwork;
	
	
	
};

/** use lapacks zhbgvx for calculating the eigenvalues of a hermitian general EVP 
 *
 * zhbgvx means: A,M is banded, A is hermitian and M is positive definite
 */
class LapackComplexBandEigenSolver : public EigenSolver<cplx,double,cplx> {
public:	 
	LapackComplexBandEigenSolver(unsigned int size, unsigned int block_size);
	virtual ~LapackComplexBandEigenSolver();
	virtual bool find_eigenvectors() throw(Exception*);

	virtual cplx eigenvalue(int nn) const throw(Exception*);
	virtual const cplx& eigenvector(int nn, int vv) const throw(Exception*);
	
	virtual SparseMatrixInterface<cplx>&   get_stiff() { return matrix_stiff; }
	virtual SparseMatrixInterface<double>& get_overlap() { return matrix_overlap; }
		
private:
	BandMatrix<cplx>   matrix_stiff;
	BandMatrix<double> matrix_overlap;
	
	cplx*   eigenvalues;
	cplx*   eigenvectors;

};

/** lapacks general eigenvalue solver routine
 * 
 * no requirements for A and M ... so i won't expect that it's fast ...
 */ 
class LapackZGGEVEigenSolver : public EigenSolver<cplx,cplx,cplx> {
public:
	LapackZGGEVEigenSolver(unsigned int size, unsigned int block_size);
	virtual ~LapackZGGEVEigenSolver();
	virtual bool find_eigenvectors() throw(Exception*);

	virtual cplx eigenvalue(int nn) const throw(Exception*);
	virtual const cplx& eigenvector(int nn, int vv) const throw(Exception*);
	
	virtual SparseMatrixInterface<cplx>& get_stiff() { return matrix_stiff; }
	virtual SparseMatrixInterface<cplx>& get_overlap() { return matrix_overlap; }
		
protected:
	FullMatrix<cplx> matrix_stiff;
	FullMatrix<cplx> matrix_overlap;
	
	cplx*   eigenvalues;
	cplx*   eigenvectors;
		
};
/** lapacks general eigenvalue solver routine (expert version) */
class LapackZGGEVXEigenSolver : public LapackZGGEVEigenSolver {
public:	
	LapackZGGEVXEigenSolver(unsigned int size, unsigned int block_size);
	virtual ~LapackZGGEVXEigenSolver() {};
	virtual bool find_eigenvectors() throw(Exception*);
};

#include "tdkp/solvers/LapackBandSolver.tcc"

} //end of namespace

#endif /*LAPACKBANDSOLVER_H_*/
