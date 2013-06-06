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

#ifndef CSRMATRIX_H_
#define CSRMATRIX_H_

#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/main/SparseMatrixInterface.h"

#define MAX_INT_SIZE 2147483637

namespace tdkp {

/** my compressed sparse row matrix */
template<class T>
class CSRMatrix : public SparseMatrixInterface<T> {

	friend class CSRMatrix<cplx>;
	friend class CSRMatrix<double>;

public:
	CSRMatrix() {}
	CSRMatrix(unsigned int size);
	CSRMatrix(unsigned int size, SparseMatrixProperties matrix_type);
	CSRMatrix(unsigned int size, SparseMatrixProperties matrix_type, bool use_fortran_indices);
	CSRMatrix(const CSRMatrix<T> &copy);
	CSRMatrix(const CSRMatrix<T> &copy, SparseMatrixProperties matrix_type);
	CSRMatrix(const char* filename);
	CSRMatrix(unsigned int size, unsigned int num_nz, const int* icol, const int* prow, const T* nonzeros, SparseMatrixProperties type);
	~CSRMatrix();

	// ---------------------------------------------------------
	// structure
	// ---------------------------------------------------------
	bool init_from_csr_available() const { return true; }
	void init_from_csr(bool symmetric, int size, int num_nz, int* prow, int* icol, T* nonzeros);
	void clear_but_keep_structure();
	void reset();
	void announce(int ii, int jj);
	void block_announce(int ii, int jj, const int* sparse_pat, int sparse_num);
	void set_structure();
	void print_sparse_structure() const;
	bool check_singular_structure() const; // deprecated [rv] (diagonal set anyway)
	bool check_matrix() const;
	bool property_is_set(SparseMatrixProperties flag) const;
	bool symmetric() const;
	bool symmetric_by_value() const;
	void perform_symmetry_analysis(); // is not const as we may write to track_nonhermitian_rows
	bool structure_is_set() const { return this->structure_set; }
	void set_block_size(unsigned int block_size_) { this->block_size = block_size_; }
	unsigned int get_block_size() const { return this->block_size; }



	template<class P>
	bool   symmetry_value_check(P a, P b, double tol) const;
	double symmetry_value_error(T a, T b) const;


	// ---------------------------------------------------------
	// access
	// ---------------------------------------------------------
	inline void  set(int ii, int jj, const T& val);
	inline void  add(int ii, int jj, const T& val);
	virtual void set_row_and_column_to_zero(int ii);
	inline T     get(int ii, int jj) const;
	inline int   pos(int ii, int jj) const;
	unsigned int get_size() const  { return size; }
	unsigned int get_num_nonzeros() const { return this->num_nz; }
	int          get_fidx() const  { return this->fidx; }

	T*   get_nonzeros() { return this->nonzeros; }
	int* get_icol()     { return this->icol; }
	int* get_prow()     { return this->prow; }
	const T*   get_nonzeros() const { return this->nonzeros; }
	const int* get_icol()     const { return this->icol; }
	const int* get_prow()     const { return this->prow; }

	const vector<double>& get_nonhermitian_rows() const;
	void  delete_nonhermitian_rows_data();

	// ---------------------------------------------------------
	// I/O
	// ---------------------------------------------------------
	void print_stats() const;
	T*   create_full_matrix() const;
	void save_to_file(const char* filename) const;
	void read_from_file(const char* filename);

	// --------------------------------------------------
	// structure cache functions
	// --------------------------------------------------
	virtual bool caching_structure_possible() const;
	virtual bool structure_cache_available(const string& identifier) const;
	virtual void load_structure(const string& identifier);
	virtual void save_structure(const string& identifier);



	// ---------------------------------------------------------
	// operations
	// ---------------------------------------------------------
	void mult_vec(const T* in, T* out) const;
	template<class B>
	void mult_vec(const B* in, B* out) const;

	template<class B>
	void mult_vec_multiple(const B* in, B* out, int kpn) const;


protected:

	class announce_list {
		public:
			announce_list() { col = -1; next = 0; prev = 0; }
			int insert(const int &col_);
			int col;
			announce_list* next;
			announce_list* prev;
	};
	void clear_all();
	void clear_announced();
	void init_announced();
	void copy_raw(const CSRMatrix<T> &copy);
	void init();
	T    internal_conj(T z) const;

	unsigned int size;
	unsigned int block_size;

	int* prow;		// array of length size + 1, gives position of row ii in nonzeros
	int* icol;  	// array of length nonzeros, gives column index of element ii at nonzeros(ii)
	T*   nonzeros;	// nonzeros ...
	int  num_nz;    // size of nonzeros
	int  fidx;	// fortran index -> if = 1, prow/icol use fortran indices 

	vector<double> track_nonhermitian_rows; // here we store whether there is a nonherm term in a row

	bool structure_set;
	bool is_symmetric;
	bool is_factorized;

	announce_list** announced;

};

/** block CSR matrix (for C.Voemel/P.Arbenz (debugging of ML in trilinos)) */
template<class T>
class BlockCSRMatrix : public CSRMatrix<T> {
public:
	BlockCSRMatrix(unsigned int size_) 
	: CSRMatrix<T>(size_) {}
	
	BlockCSRMatrix(unsigned int size_, SparseMatrixProperties matrix_type_)
	: CSRMatrix<T>(size_, matrix_type_) {}
	void announce(int ii, int jj);
	void block_announce(int ii, int jj, const int* sparse_pat, int sparse_num);
	void set_block_size(unsigned int block_size_) {	
		CSRMatrix<T>::set_block_size(block_size_);	
	}
private:

};		

template<class T>
void BlockCSRMatrix<T>::block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num) {	
	if(this->is_symmetric && gii > gjj) {
		TDKP_GENERAL_EXCEPTION("symmetric matrix but gii > gjj!");	
	} else {
		if(this->is_symmetric && gii == gjj) {
			for(unsigned int dd = 0; dd < this->block_size; dd++) {
				for(unsigned int cc = dd; cc < this->block_size; cc++) {				
					CSRMatrix<T>::announce(gii * this->block_size + cc, gjj * this->block_size + dd);	
				}
			}  				
		} else {
			for(unsigned int cc = 0; cc < this->block_size; cc++) {
				for(unsigned int dd = 0; dd < this->block_size; dd++) {
					CSRMatrix<T>::announce(gii * this->block_size + cc, gjj * this->block_size + dd);	
				}
			}  	
		}
	}
}

template<class T>
void BlockCSRMatrix<T>::announce(int ii, int jj) {	
	TDKP_BOUNDS_ASSERT(this->block_size > 1, "");
	int bii = (ii - (ii % this->block_size));
	int bjj = (jj - (jj % this->block_size));	
	for(unsigned int cc = 0; cc < this->block_size; cc++) {
		for(unsigned int dd = 0; dd < this->block_size; dd++) {
			CSRMatrix<T>::announce(bii + cc, bjj + dd);	
		}
	}  	
}
		
		
/** create CSRMatrix from sparse arrays
 *
 */
template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int size_, unsigned int num_nz_, const int* icol_, const int* prow_, const T* nonzeros_, SparseMatrixProperties type) {
	
	// attention icol and prow is swapped!
	this->init_from_csr(type == symmetric_matrix, size_, num_nz_, prow_, icol_, nonzeros_);  
	
}
/** create matrix from crs structures */
template<class T>
void CSRMatrix<T>::init_from_csr(bool symmetric_, int size_, int num_nz_, int* prow_, int* icol_, T* nonzeros_) {
	
	this->size 			= size_;
	this->is_symmetric 	= symmetric_;
	this->num_nz 		= num_nz_;
	this->structure_set = true;
	this->is_factorized = false;
	this->fidx          = prow_[0];
	
	TDKP_ASSERT(prow_[0] == 0 || prow_[0] == 1,"");

	this->icol          = new int[this->num_nz];
	this->prow          = new int[this->size + 1];
	this->nonzeros      = new T[this->num_nz];

	TDKP_POINTER_ASSERT(this->icol && this->prow && this->nonzeros);

	for(unsigned int ii = 0; ii < this->size + 1; ii++) {
		this->prow[ii] = prow_[ii];
	}
	for(int ii = 0; ii < this->num_nz; ii++) {
		this->icol[ii]     = icol_[ii];
		this->nonzeros[ii] = nonzeros_[ii];
	}
	
}

/** compatibility to SparseMatrixInterface */
template<class T>
bool CSRMatrix<T>::property_is_set(SparseMatrixProperties flag) const {
	if(flag == symmetric_matrix && this->symmetric()) {
		return true;
	}
	if(flag == nonsymmetric_matrix && !this->symmetric()) {
		return true;
	}
	return false;
}

/** create SYMMETRIC compressed sparse row matrix
 *
 * to allow switching between CSCMatrix used by arpack++ and symmetric CSR matrices
 * used by JDQZ and ARPACK + Pardiso, we have to provide the class with a constructor
 * which only takes the size as input value. therefore we default here to a symmetric
 * matrix to allow simple switching via redefinition of type GMatrix.
 * @param size the size of the symmetric matrix
 */
template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int size_) {
	if(size_ == 0) {
		TDKP_GENERAL_EXCEPTION("invalid matrix size(0) set");
	}
	this->init(); // sets all to zero except size
	this->is_symmetric = true;
	this->size         = size_;
	this->init_announced();
	this->clear_all();
}

/** create compressed sparse row matrix
 *
 * @param size_     size of matrix
 * @param symmetric boolean to indicate whether matrix is symmetric or not
 */
template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int size_, SparseMatrixProperties matrix_type) : size(size_) {
	if(size_ == 0) {
		TDKP_GENERAL_EXCEPTION("invalid matrix size(0) set");
	}
	this->init(); // sets all to zero except size
	this->is_symmetric = (matrix_type == symmetric_matrix ? true : false);
	this->init_announced();
	this->clear_all();
}

/** create compressed sparse row matrix
 * @param use_fortran_indices tell me whether i should set fidx to 0
 */
template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int size_, SparseMatrixProperties matrix_type, bool use_fortran_indices)
: size(size_) {
	if(size_ == 0) {
		TDKP_GENERAL_EXCEPTION("invalid matrix size(0) set");
	}
	this->init(); // sets all to zero except size
	this->init_announced();
	this->is_symmetric = (matrix_type == symmetric_matrix ? true : false);
	this->clear_all();
	if(use_fortran_indices) {
		this->fidx = 1;
	} else {
		this->fidx = 0;
	}
}

/** create matrix and read data from file
 *
 * @param matrix file where matrix was saved with @see CSRMatrix<T>::save_to_file
 */
template<class T>
CSRMatrix<T>::CSRMatrix(const char* filename) {
	this->init();
	this->clear_all();
	this->read_from_file(filename);
}


/** copy constructor to create a deep copy of a CSRMatrix
 *
 * only matrices with initialized structures may be copied. if
 * copy has no structure set, we throw an exception
 *
 * calls copy_raw
 */
template<class T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix<T> &copy) {
	this->init();
	this->copy_raw(copy);
}

/** copy constructor to create a deep copy of a CSRMatrix, but allowing to create a nonsymmetric matrix from a symmetric one
 *
 * only matrices with initialized structures may be copied. if
 * copy has no structure set, we throw an exception
 * also if symmetric == true && copy.symmetric == false, an exception is thrown
 *
 * @param copy matrix to copy
 * @param symmetric boolean to indicate if new matrix should be reagarded as symmetric or not
 */
template<class T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix<T> &copy, SparseMatrixProperties matrix_type) {
	this->init();
	if(copy.icol && copy.prow && copy.nonzeros) {
		if(matrix_type == symmetric_matrix) {
			// symmetric to symmetric -> copy raw
			if(copy.is_symmetric) {
				this->copy_raw(copy);
			} else {
				TDKP_GENERAL_EXCEPTION("can not create symmetric matrix from unsymmetric matrix");
			}
		} else {
			// nonsymmetric to nonsymmetric -> copy raw
			if(!copy.is_symmetric) {
				this->copy_raw(copy);
			} else {
				this->prow 	 	    = 0;
				this->icol  	 	= 0;
				this->nonzeros 	    = 0;
				this->size          = copy.size;
				this->structure_set = false;
				this->is_factorized = false;
				this->is_symmetric  = false; // structure is not symmetric
				this->init_announced();

				this->clear_all();
				int until, from, colidx;
				T tmp;
				// announce
				for(int ii = 0; ii < (signed)copy.size; ii++) {
					until    = copy.prow[ii + 1] - copy.fidx;
					from     = copy.prow[ii] - copy.fidx;
					this->announce(ii,ii);
					for(int jj = from + 1; jj < until; jj++) {
						colidx = copy.icol[jj] - copy.fidx;
						this->announce(ii,colidx);
						this->announce(colidx,ii);
					}
				}
				// create matrix
				this->set_structure();
				// fill in
				for(int ii = 0; ii < (signed)copy.size; ii++) {
					until    = copy.prow[ii + 1] - copy.fidx;
					from     = copy.prow[ii] - copy.fidx;
					tmp      = copy.nonzeros[from];
					this->set(ii,ii, tmp);
					TDKP_ASSERT((copy.icol[from] - copy.fidx) == ii, "(copy.icol[from] - copy.fidx) == ii");
					for(int jj = from + 1; jj < until; jj++) {
						tmp = copy.nonzeros[jj];
						colidx = copy.icol[jj] - copy.fidx;
						this->set(ii,colidx, tmp);
						// set cplx conjugate
						this->set(colidx,ii, internal_conj(tmp));
					}
				}
			}
		}
	} else {
		TDKP_GENERAL_EXCEPTION("can not create CSRMatrix from uninitialized one");
	}
}


/** destructor
 *
 * calls clear_all and clear_announced to kill allocated data in heap
 */
template<class T>
CSRMatrix<T>::~CSRMatrix() {
	this->clear_all();
	this->clear_announced();
}

/** copy function
 *
 * performs a deep copy of sparse data
 * @param copy CSRMatrix to copy
 */
template<class T>
void CSRMatrix<T>::copy_raw(const CSRMatrix<T> &copy) {

	// abort if object already initialized
	if(this->announced || this->prow || this->icol || this->nonzeros) {
		TDKP_GENERAL_EXCEPTION("copy_raw may not be used on initialized objects");
	}

	this->size 	   		= copy.size;
	this->fidx          = copy.fidx;
	this->num_nz   		= copy.num_nz;
	this->prow 	   		= 0;
 	this->icol 	   		= 0;
	this->nonzeros 		= 0;
	this->announced 	= 0;
	this->structure_set = copy.structure_set;
	this->is_symmetric  = copy.is_symmetric;
	this->is_factorized = copy.is_factorized;

	// check if we copy the matrix itself
	if(copy.icol && copy.prow && copy.nonzeros) {
		// create copy of full csc
		icol     = new int[num_nz];
		prow     = new int[size + 1];
		nonzeros = new T[num_nz];
		if(!icol || !prow || ! nonzeros) {
			TDKP_GENERAL_EXCEPTION("can not allocate memory for matrix");
		}
		for(int ii = 0; ii < num_nz; ii++) {
			this->icol[ii]     = copy.icol[ii];
			this->nonzeros[ii] = copy.nonzeros[ii];
		}
		for(unsigned int ii = 0; ii <= size; ii++) {
			this->prow[ii] = copy.prow[ii];
		}
	} else {
		// structure not set, throw exception
		TDKP_GENERAL_EXCEPTION("can not copy matrix with empty structure");
	}
}
/*
template<class T>
CSRMatrix<cplx>* CSRMatrix<T>::create_complex_matrix() const {

	CSRMatrix<cplx>* ret = new CSRMatrix<cplx>(this->get_size());
	ret->clear_announced();
	ret->size 	   	   = this->size;
	ret->fidx          = this->fidx;
	ret->num_nz   	   = this->num_nz;
	ret->prow 	   	   = new int[this->size + 1];
 	ret->icol 	   	   = new int[this->num_nz];
	ret->nonzeros 	   = new cplx[this->num_nz];
	ret->announced 	   = 0;
	ret->structure_set = this->structure_set;
	ret->is_symmetric  = this->is_symmetric;
	ret->is_factorized = this->is_factorized;
	if(!ret->icol || !ret->prow || !ret->nonzeros) {
		TDKP_GENERAL_EXCEPTION("can not allocate memory for matrix");
	}
	for(int ii = 0; ii < this->num_nz; ii++) {
		ret->icol[ii]     = this->icol[ii];
		ret->nonzeros[ii] = this->nonzeros[ii];
	}
	for(unsigned int ii = 0; ii <= this->size; ii++) {
		ret->prow[ii] = this->prow[ii];
	}

	return ret;

}
*/

/** init all member variables except size to zero
 */
template<class T>
void CSRMatrix<T>::init() {
	this->fidx          = 0; // prow/icol have fortran indices (for paradiso)
	this->block_size    = 1;
	this->prow 	 	    = 0;
	this->icol  	 	= 0;
	this->nonzeros 	    = 0;
	this->structure_set = false;
	this->is_factorized = false;
	this->is_symmetric  = false;
	this->announced     = 0;
	this->num_nz        = 0;
}

/** set value at (ii,jj)
 *
 * if matrix was initialized as symmetric, only the upper triangle matrix may be accessed. if lower triangle is accessed, an exception will be thrown
 * also, setting of an nonannouced position results in an exception
 *
 * @param ii row index
 * @param jj col index
 */
template<class T>
inline void CSRMatrix<T>::set(int ii, int jj, const T& val) {
	TDKP_BOUNDS_ASSERT(ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size, "ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size");

	if(this->is_symmetric && ii > jj) {
		TDKP_GENERAL_EXCEPTION("setting of lower triangle matrix for a symmetric matrix is forbidden");
	}

	TDKP_BOUNDS_ASSERT(prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx, "prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx");
	TDKP_BOUNDS_ASSERT(prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx, "prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx");

	for(int rr = prow[ii]; rr < prow[ii + 1]; rr++) {
		if(icol[rr - fidx] == (jj + fidx)) {
			nonzeros[rr - fidx] = val;
			return;
		}
	}

	throw new Exception(1000, "tried to set value to not announced pos: %d/%d", ii, jj);
}

/** sets all values in row to zero (required for dirichlet b.c.) */
template<class T>
void CSRMatrix<T>::set_row_and_column_to_zero(int aa) {
	
	// set row to 0
	for(int rr = prow[aa]; rr < prow[aa + 1]; rr++) {
		nonzeros[rr - fidx] = 0.0;
	}
	// set column to 0
	// for every row
	#pragma omp parallel for schedule(dynamic,20000)
	for(int ii = 0; ii < (signed)this->size; ii++) {
		// scan through row
		for(int rr = prow[ii] - fidx; rr < prow[ii + 1] - fidx; rr++) {
			// find column aa
			if((this->icol[rr] - fidx) == aa) {
				nonzeros[rr] = 0.0;
				break; // coulmn aa already found	
			}
		}
	}
}


/** return position in nonzero vector where to value of (ii,jj) resides! */
template<class T>
inline int CSRMatrix<T>::pos(int ii, int jj) const {
	for(int rr = prow[ii]; rr < prow[ii + 1]; rr++) {
		if(icol[rr - fidx] == (jj + fidx)) {
			return rr - fidx;
		}
	}
	return -1;
}

/** add to value at (ii,jj)
 *
 * if matrix was initialized as symmetric, only the upper triangle matrix may be accessed. if lower triangle is accessed, an exception will be thrown
 * also, setting of an nonannouced position results in an exception
 *
 * @param ii row index
 * @param jj col index
 */
template<class T>
inline void CSRMatrix<T>::add(int ii, int jj, const T& val) {

	TDKP_BOUNDS_ASSERT(ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size, "ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size");

	if(this->is_symmetric && ii > jj) {
		TDKP_GENERAL_EXCEPTION("setting of lower triangle matrix for a symmetric matrix is forbidden");
	}
	TDKP_BOUNDS_ASSERT(prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx, "prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx");
	TDKP_BOUNDS_ASSERT(prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx, "prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx");

	for(int rr = prow[ii]; rr < prow[ii + 1]; rr++) {
		if(icol[rr - fidx] == jj + fidx) {
			nonzeros[rr - fidx] += val;
			return;
		}
	}
	throw new Exception(1000, "tried to set value to not announced pos: %d/%d", ii, jj);
}


/** returns value at (ii,jj)
 *
 * @param ii row index [0, size - 1]
 * @param jj col index [0, size - 1]
 */
template<class T>
inline T CSRMatrix<T>::get(int ii, int jj) const {
	int tmp;
	bool return_conjugate = false;
	TDKP_BOUNDS_ASSERT(ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size, "ii >= 0 && ii < (signed)size && jj >= 0 && jj < (signed)size");

	if(this->is_symmetric && ii > jj) {
		tmp = jj;
		jj  = ii;
		ii  = tmp;
		return_conjugate = true;
	}
	TDKP_BOUNDS_ASSERT(prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx, "prow[ii]   >= 0 + fidx && prow[jj]   <= num_nz + fidx");
	TDKP_BOUNDS_ASSERT(prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx, "prow[ii+1] >= 0 + fidx && prow[ii+1] <= num_nz + fidx");

	for(int rr = prow[ii]; rr < prow[ii + 1]; rr++) {
		if(icol[rr - fidx] == jj + fidx) {
			if(return_conjugate) {
				return internal_conj(nonzeros[rr - fidx]);
			} else {
				return nonzeros[rr - fidx];
			}
		}
	}
	return (T)0;
}

/** delete sparse matrix data and structure
 *
 * attention: used for destructors, does not reinit the matrix
 */
template<class T>
void CSRMatrix<T>::clear_all() {
	if(icol)     { delete[] icol; icol = 0; }
	if(prow)     { delete[] prow; prow = 0; }
	if(nonzeros) { delete[] nonzeros; nonzeros = 0;}
	num_nz 		  = 0;
	structure_set = false;
}

/** creates the announce structure on the diagonal
 */
template<class T>
void CSRMatrix<T>::init_announced() {
	this->announced     = new announce_list*[this->size];
	for(unsigned int ii = 0; ii < this->size; ii++) {
		this->announced[ii] = new announce_list();
		this->announced[ii]->col = ii; // set the diagonal
	}
}

template<class T>
void CSRMatrix<T>::reset() {
	this->clear_all();
	this->clear_announced();
	this->init_announced();
}

/** destroy the announced structure
 */
template<class T>
void CSRMatrix<T>::clear_announced() {
	announce_list* ptr;
	announce_list* next;
	announce_list* prev;

	// if announced array is initialized
	if(this->announced) {
		// loop over announced array
		for(unsigned int ii = 0; ii < size; ii++) {
			if(this->announced[ii]) {
				// if element is set, descend tree and delete remaining list entries
				// delete prev
				ptr  = 0;
				prev = this->announced[ii]->prev;
				while(prev) {
					ptr  = prev;
					prev = ptr->prev;
					delete ptr;
				}
				// delete next
				ptr  = 0;
				next = this->announced[ii];
				while(next) {
					ptr = next;
					next = ptr->next;
					delete ptr;
				}
			}
		}
		delete[] this->announced;
		this->announced = 0;
	}
}

/** insert col_ into list tree
 *
 * @param col some value (here: column index)
 * @return 1 if col was inserted, 0 if col is already in tree (therefore allowing to count number of nonzeros)
 */
template<class T>
int CSRMatrix<T>::announce_list::insert(const int &col_) {

	if(col_ == this->col) {
		return 0;
	} else if(col_ > this->col) {
		// pass to next?
		if(this->next && this->next->col <= col_) {
			return this->next->insert(col_);
		// no next avialble, add new at end
		} else if(!next) {
			this->next       = new announce_list();
			this->next->col  = col_;
			this->next->next = 0;
			this->next->prev = this;
			return 1;
		// value is between this and next -> insert new one
		} else {
			announce_list* ptr = new announce_list();
			ptr->col         = col_;
			ptr->next        = this->next;
			ptr->prev  		 = this; // me is prev of next
			this->next->prev = ptr;  // old next has now new next as prev instead of me
			this->next = ptr;
			return 1;
		}
	} else if(col_ < this->col) {
		// previous available, pass to prev?
		if(this->prev && this->prev->col >= col_) {
			return this->prev->insert(col_);
		// no prev available, add new at beginning
		} else if(!this->prev) {
			announce_list* ptr = new announce_list();
			ptr->col 		   = col_;
			ptr->next          = this;
			ptr->prev          = 0;
			this->prev         = ptr;
			return 1;
		// add in between prev and this
		} else {
			announce_list* ptr = new announce_list();
			ptr->col         = col_;
			ptr->next        = this;
			ptr->prev        = this->prev;
			this->prev->next = ptr;
			this->prev       = ptr;
			return 1;
		}
	} else {
		TDKP_GENERAL_EXCEPTION("should not reach this point");
	}
}

/** set values to zero
 */
template<class T>
void CSRMatrix<T>::clear_but_keep_structure() {
	if(this->nonzeros) {
		for(int ii = 0; ii < num_nz; ii++) {
			this->nonzeros[ii] = (T)0.0;
		}
	}
}

/** announce value at position (ii,jj)
 *
 * announce that there is going to be a nonzero value at position (ii,jj)
 * it is mandatory that all nonzero values are announced first, then
 * the matrix structure is build with set_structure
 * after that, the functions set, add, and get can be used
 *
 * @param ii row index [0, size - 1]
 * @param jj col index [0, size - 1]
 */
template<class T>
void CSRMatrix<T>::announce(int ii, int jj) {
	if(this->is_symmetric && ii > jj) {
		TDKP_GENERAL_EXCEPTION("announcing lower triangle matrix elements for a symmetric matrix is forbidden");
	}
	if(!this->announced) {
		TDKP_GENERAL_EXCEPTION("announced not initialized");
	}
	if(!this->announced[ii]) {
		this->announced[ii] = new announce_list();
		this->announced[ii]->col = jj;

#pragma omp atomic
		this->num_nz++;

	} else {
		int tt = this->announced[ii]->insert(jj);
#pragma omp atomic
		this->num_nz += tt;

		TDKP_ASSERT(this->num_nz < MAX_INT_SIZE, "this>num_nz < MAX_INT_SIZE");
	}
}

/** announce a whole block to the matrix
 *
 * instead of announcing a single entry, we can also announce whole kp blocks
 * which is way cheaper in a shared object world
 * further, also the incrementing of the shared number of nonzeros
 * is reduced and therefore the code should be faster when running parallel
 */
template<class T>
void CSRMatrix<T>::block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num) {

	int added_num_nz = 0;
	int ii,jj;

	if(this->is_symmetric && gii > gjj) {
		TDKP_GENERAL_EXCEPTION("announcing lower triangle matrix elements for a symmetric matrix is forbidden");
	} else if(this->is_symmetric && gii == gjj) {
		// symmetric diagonal, only announce upper triangle
		for(int ss = 0; ss < sparse_num; ss++) {
			ii = gii * block_size + sparse_pat[2*ss];
			jj = gjj * block_size + sparse_pat[2*ss + 1];
			if(!(ii > jj)) {
				TDKP_BOUNDS_ASSERT(this->announced[ii] != 0, "this->announced[ii] != 0");
				added_num_nz += this->announced[ii]->insert(jj);
			}
		}
	} else {
		// nonsymmetric part
		for(int ss = 0; ss < sparse_num; ss++) {
			ii = gii * block_size + sparse_pat[2*ss];
			jj = gjj * block_size + sparse_pat[2*ss + 1];
			TDKP_BOUNDS_ASSERT(this->announced[ii] != 0, "this->announced[ii] != 0");
			added_num_nz += this->announced[ii]->insert(jj);
		}
	}
#pragma omp atomic
	this->num_nz += added_num_nz;

	TDKP_ASSERT(this->num_nz < MAX_INT_SIZE, "this>num_nz < MAX_INT_SIZE");
}


/** sets announced structure to matrix and deletes announcement tree
 */
template<class T>
void CSRMatrix<T>::set_structure() {

	// we need to start at a fresh point
	if(icol || prow || nonzeros) {
		TDKP_GENERAL_EXCEPTION("set structure called twice.");
	}
	// we need an announcement tree
	if(!this->announced) {
		TDKP_GENERAL_EXCEPTION("announcment tree not build");
	}

	// structure must have diagonal elements (requested by pardiso),
	// therefore we have to add them to the counter
	this->num_nz += this->size;

	// we need at least size  nonzero elements
	if((signed)this->size > this->num_nz) {
		TDKP_GENERAL_EXCEPTION("there are less nonzero elements than matrix size ... matrix is therefore surely singular :-(");
	}

	ostringstream statout;
	double tmp_size = static_cast<double>(size);
	long memusage = (num_nz * sizeof(T) + (num_nz + size + 1) * sizeof(int)) / 1024 / 1024;
	statout << "CSRMatrix: building matrix (N = " << size << ") with " << num_nz << " nonzeros. fill in is "
			<< (100.0 * ((double)num_nz / (tmp_size * tmp_size))) << "%, memory usage is " << memusage << "M."; 
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, statout.str());

	// allocate
	nonzeros = new   T[num_nz];
	icol     = new int[num_nz];
	prow     = new int[size + 1];

	if(!nonzeros || !icol || !prow) {
		TDKP_GENERAL_EXCEPTION("could not allocate memory for matrix");
	}

	announce_list* ptr;
	announce_list* ptr_del;

	int ii = 0;
	// loop over rows and create structure
	for(unsigned int rr = 0; rr < size; rr++) {
		// stop if row was not announced
		if(!this->announced[rr]) {
			std::ostringstream sout;
			sout << "error in row " << rr << ". no element announced";
			TDKP_GENERAL_EXCEPTION(sout.str());
		} else {
			// set current pos to prow
			prow[rr] = ii + fidx; // include fortran indices
			ptr = this->announced[rr];
			// get first entry
			while(ptr->prev) {
				ptr = ptr->prev;
			}
			// loop over entries, get colums and delete linked list
			while(ptr) {
				this->nonzeros[ii] = (T)0.0;
			 	this->icol[ii]     = ptr->col + fidx;
			 	ii++;
			 	ptr_del = ptr;
				ptr     = ptr->next;
				delete ptr_del;
			}
			this->announced[rr] = 0;
		}
	}
	// set last row pointer
	prow[size] 		= ii + fidx;
	TDKP_ASSERT(prow[size] == num_nz + fidx, "prow[size] == num_nz + fidx");

	// delete announced
	delete[] this->announced;
	this->announced 	= 0;
	this->structure_set = true;
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "CSRMatrix: finished creating matrix structure");

}

/** creates full matrix
 *
 * @return pointer to newly created array of length size * size where element (ii,jj) is at ret[ii * size + jj]
 */
template<class T>
T* CSRMatrix<T>::create_full_matrix() const {
	T* full = new T[this->size * this->size];
	for(unsigned int ii = 0; ii < this->size; ii++) {
		for(unsigned int jj = 0; jj < this->size; jj++) {
			full[ii * this->size + jj] = this->get(ii,jj);
		}
	}
	return full;
}


/** check whether the matrix structure has errors
 */
template<class T>
bool CSRMatrix<T>::check_matrix() const {
  int i, j, k;

  // Checking if prow is in ascending order.
  i = 0;
  while ((i!=(signed)this->size)&&(prow[i]<=prow[i+1])) i++;
  if (i!=(signed)this->size) {
  	 Logger::get_instance()->emit(LOG_INFO_DEVEL2, 1000, "prow not in ascendig order at i = %d: %d next %d, size: %d, num_nz: %d", i, prow[i], prow[i+1],this->size, this->num_nz) ;
  	 return false;
  }

  // Checking if icol components are in order and within bounds.
  for (i=0; i!=(signed)this->size; i++) {
    j = prow[i]     - fidx;
    k = prow[i+1]-1 - fidx;
    if (j<=k) {
      if ((icol[j] < (0 + fidx)) || (icol[k] >=((signed)this->size + fidx))) {
      	Logger::get_instance()->emit(LOG_INFO_DEVEL2, 1000, "icol[k] %d >= this->n %d,  k: %d", icol[k], this->size, k);
      	return false;
      }
      while ((j!=k)&&(icol[j]<icol[j+1])) j++;
      if (j!=k) {
      	Logger::get_instance()->emit(LOG_INFO_DEVEL2, 1000, "i != k %d %d", i, k);
      	return false;
      }
    }
  }
  return true;
}

/** deprecated function. the matrix is never singular as the diagonal is set anyway
 */
template<class T>
bool CSRMatrix<T>::check_singular_structure() const {
	if(prow && icol) {
		Logger::get_instance()->init_progress_bar("checking whether matrix structure is singular", size);
		int  next = 0;
		bool have_col_ii;
		for(int ii = 0; ii < (signed)size; ii++) {
			have_col_ii = false;
			if(prow[ii] == prow[ii + 1]) {
				Logger::get_instance()->end_progress_bar();
				Logger::get_instance()->emit(LOG_ERROR, 1000, "matrix is singular for row: %d", ii);
				return true;
			}
			for(int jj = 0; jj < num_nz; jj++) {
				if(icol[jj] == ii + fidx) {
					have_col_ii = true;
					break;
				}
			}
			if(!have_col_ii) {
				Logger::get_instance()->end_progress_bar();
				Logger::get_instance()->emit(LOG_ERROR, 1000, "matrix is singular for col: %d", ii);
				return true;
			}
			if(ii == next) {
				next = Logger::get_instance()->set_progress_bar(ii, size);
			}
		}
		Logger::get_instance()->end_progress_bar();
	} else {
		Logger::get_instance()->emit(LOG_WARN,"check_singular_structure called while structure not set yet");
	}
	return false;
}

/** returns whether the matrix has symmetric structure
 */
template<class T>
bool CSRMatrix<T>::symmetric() const {
	return this->is_symmetric;
}

/** returns whether the matrix has symmetric strcture OR the values of a nonsymmetric structure are symmetric
 */
template<class T>
bool CSRMatrix<T>::symmetric_by_value() const {
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, "checking matrix symmetry");
	if(this->is_symmetric) return true;
	if(this->size < 10000) {
		for(int ii = 0; ii < (signed)this->size; ii++) {
			for(int jj = ii; jj < (signed)this->size; jj++) {
				if(tdkp_math::abs(this->get(ii,jj) - internal_conj(this->get(jj,ii))) > 1.0e-9) {
					return false;
				}
			}
		}
		return true;
	} else {
		Logger::get_instance()->emit(LOG_WARN, "not checking whether matrix is symmetric");
		return true;
	}
}


template<class T>
void CSRMatrix<T>::perform_symmetry_analysis() {

	if(!this->is_symmetric) {
		Logger::get_instance()->init_progress_bar("checking matrix symmetry",size);
		int next = 0;
		// ------------------------------------------------
		// loop control variables
		// ------------------------------------------------
		int until,from,mid;
		int kk;
		int display_errors = 30;
		ostringstream sout;
		T value_ii_kk, value_kk_ii;
		double tol = Configuration::get_instance()->get("assembly_check_matrix_for_symmetry_tolerance");

		double num_wrong   = 0;
		double num_correct = 0;
		double num_offdiag = 0;

		double avg_value_wrong   = 0.0;
		double avg_value_correct = 0.0;
		double avg_value_offdiag = 0.0;
		double avg_value_diag    = 0.0;
		double avg_error         = 0.0;

		double abs_value_ii_kk, abs_value_kk_ii;

		// we may track and output to dfise file whether nodes are
		// nonhermitian or not ...
		bool track_nodes = Configuration::get_instance()->get("assembly_track_nonhermitian_nodes") != 0 ? true : false;
		if(track_nodes) {
			this->track_nonhermitian_rows.resize(this->size);
			for(vector<double>::iterator it = this->track_nonhermitian_rows.begin();
			    it != this->track_nonhermitian_rows.end(); it++) {
			    *it = 0.0;
			}
		}

		for(int ii = 0; ii < (signed)this->size; ii++) {
			until    = this->prow[ii + 1] - this->fidx;
			from     = this->prow[ii] - this->fidx;
			// advance to the middle element
			mid = -1;
			for(int jj = from; jj < until; jj++) {
				if((this->icol[jj] - this->fidx) == ii) {
					mid = jj;
					break;
				}
			}
			TDKP_ASSERT(mid != -1, "middle element != -1");

			avg_value_diag    = (avg_value_diag * ii + tdkp_math::abs(this->nonzeros[mid])) / (ii + 1);

			// loop over upper right triangle
			for(int jj = mid + 1; jj < until; jj++) {
				kk = this->icol[jj] - this->fidx;
				value_ii_kk = this->nonzeros[jj];
				value_kk_ii = this->get(kk,ii);

				abs_value_ii_kk = tdkp_math::abs(value_ii_kk);
				abs_value_kk_ii = tdkp_math::abs(value_kk_ii);
				// check if correct
				if(this->symmetry_value_check(value_ii_kk,value_kk_ii,tol)) {
					avg_value_correct = ((avg_value_correct * num_correct) + abs_value_ii_kk) / (num_correct + 1);
					num_correct += 1;
				} else {
					avg_value_wrong = ((avg_value_wrong * num_wrong) + abs_value_ii_kk) / (num_wrong + 1);
					avg_error       = ((avg_error       * num_wrong) + symmetry_value_error(value_ii_kk, value_kk_ii)) / (num_wrong + 1);
					num_wrong += 1;
					if(track_nodes) {
						this->track_nonhermitian_rows[ii] += fabs(symmetry_value_error(value_ii_kk, value_kk_ii));
						this->track_nonhermitian_rows[kk] += fabs(symmetry_value_error(value_ii_kk, value_kk_ii));
					}
					if(display_errors > 0) {
						sout.str("");
						sout << "nonhermitian in (" << ii;
						if(this->block_size > 1) {
							sout << " [" << ii % this->block_size << "]";
						}
						sout << "/" << kk;
						if(this->block_size > 1) {
							sout << " [" << ii % this->block_size << "]";
						}
						sout << "): "
						     << value_ii_kk << " vs. "
						     << value_kk_ii << " err: "
						     << symmetry_value_error(value_ii_kk, value_kk_ii);
						display_errors--;
						Logger::get_instance()->emit(LOG_WARN, sout.str());
					}
				}
				avg_value_offdiag = (avg_value_offdiag * num_offdiag + abs_value_ii_kk + abs_value_kk_ii) / (num_offdiag + 2);
				num_offdiag += 2;
			}
			if(ii == next) {
				next = Logger::get_instance()->set_progress_bar(ii, size);
			}

		}
		double average_error_ratio = 0.0;
		if(avg_value_wrong != 0.0) {
			average_error_ratio = (avg_error / avg_value_wrong);
		}
		Logger::get_instance()->end_progress_bar();

		// ------------------------------------------------------
		// calculate matrix waste
		// ------------------------------------------------------
		unsigned int entries_are_zero = 0;
		for(int ii = 0; ii < this->num_nz; ii++) {
			if(this->nonzeros[ii] == 0.0) {
				entries_are_zero++;
			}
		}


		sout.str("");
		sout << "the " << size << "x" << size << " matrix is " << (num_wrong > 0 ? " NOT " : "") << "symmetric (hermitian)\n"
		     << "wasted space (zeros): " << setw(10) << entries_are_zero << ", in prnt:  " << (100.0 * static_cast<double>(entries_are_zero) / static_cast<double>(num_nz)) << "\n"
		     << "diagonal entries:     " << setw(10) << size             << ", avg size: " << avg_value_diag << "\n"
		     << "offdiagonal entries:  " << setw(10) << int(num_offdiag) << ", avg size: " << avg_value_offdiag << "\n"
		     << "wrong pairs:          " << setw(10) << int(num_wrong)   << ", avg size: " << avg_value_wrong << "\n"
		     << "average error:        " << setw(10) << avg_error        << ", ratio:    " << average_error_ratio << "\n"
		     << "correct pairs:        " << setw(10) << int(num_correct) << ", avg size: " << avg_value_correct;

		Logger::get_instance()->emit(LOG_INFO, sout.str());

	} else {
		Logger::get_instance()->emit(LOG_INFO, "symmetry analysis: matrix is already defined symmetric. therefore not analysing anything.");
	}
}

/** return vector of matrix size containing the value of nonhermiticity in given row
 *
 * vector is only available if assembly_track_nonhermitian_nodes is set to 1
 */
template<class T>
const vector<double>& CSRMatrix<T>::get_nonhermitian_rows() const {
	return this->track_nonhermitian_rows;
}
/** purge nonhermitian row data */
template<class T>
void CSRMatrix<T>::delete_nonhermitian_rows_data() {
	this->track_nonhermitian_rows.resize(0);
}


template<class T>
void CSRMatrix<T>::print_sparse_structure() const {
	std::cout << "icol: ";
	for(int ii = 0; ii < num_nz; ii++) {
		std::cout << icol[ii] << "  ";
	}
	std::cout << "\nprow: ";
	for(int ii = 0; ii <= size; ii++) {
		std::cout << prow[ii] << "  ";
	}
	std::cout << " / size would be: " << num_nz << "\nnz: ";
	for(int ii = 0; ii < num_nz; ii++) {
		std::cout << nonzeros[ii] << "  ";
	}
	std::cout << "\n";
}

template<class T>
std::ostream& operator<<(std::ostream& stream,const CSRMatrix<T> &mat) {
	int width = 14;
	stream.precision(width - 5);
	stream.setf(std::ios::left);

	for(int ii = 0; ii < (signed)mat.get_size(); ii++) {
		for(int jj = 0; jj < (signed)mat.get_size(); jj++) {
			stream << std::setw(width) << mat.get(ii,jj) << " ";
		}
		stream << "\n";
	}
	return stream;
}

template<class T>
void CSRMatrix<T>::print_stats() const {
	T offdiag = 0;
	T diag    = 0;
	int num_offdiag = 0;
	for(int ii = 0; ii < (signed)this->get_size(); ii++) {
		diag += abs(this->get(ii,ii));
		for(int jj = ii + 1; jj < (signed)this->get_size(); jj++) {
			offdiag += abs(this->get(ii,jj));
			num_offdiag++;
		}
	}
	ostringstream sout;
	sout << "diagonal: num = " << this->get_size() << " val: " << diag / (double)this->get_size() << "  "
	     << "offdiag:  num = " << 2 * num_offdiag << " val: " << offdiag / (double)num_offdiag << "\n";
	Logger::get_instance()->emit(LOG_INFO, sout.str());
}

template<class T>
void CSRMatrix<T>::save_to_file(const char* filename) const {
	if(prow && icol	&& nonzeros) {
		std::fstream fout(filename, std::ios::out);
		if(fout) {
			fout << this->size << "  " << this->num_nz << "  "
				 << this->is_symmetric << "  " << this->is_factorized << "\n";
			for(unsigned int ii = 0; ii < size + 1; ii++) {
				fout << prow[ii] << "\n";
			}
			for(int ii = 0; ii < num_nz; ii++) {
				fout << icol[ii] << "\n";
			}
			fout.precision(20);
			for(int ii = 0; ii < num_nz; ii++) {
				fout << nonzeros[ii] << "\n";
			}
			fout << 101 << "\n";
			fout.close();
		} else {
			TDKP_FILE_EXCEPTION("can not open file", filename);
		}
	} else {
		TDKP_GENERAL_EXCEPTION("can not save noninitialized matrix");
	}
}

template<> void CSRMatrix<cplx>::mult_vec(const cplx* in, cplx* out) const;


/** multiply a vector with matrix Ain = out
 *
 * @param in  input vector
 * @param out output vector
 */
template<class T>
void CSRMatrix<T>::mult_vec(const T* in, T* out) const {

	if(this->prow && this->icol	&& this->nonzeros) {

		int until, from;
		// --------------------------------------------------
		// distinguish between symmetric and nonsymmetric
		// --------------------------------------------------
		if(this->is_symmetric) {
			// clean up out
			for(int ii = 0; ii < (signed)this->size; ii++) {
				out[ii] = (T)0.0;
			}
			// perform symmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - fidx;
				from     = this->prow[ii] - fidx;
				out[ii] += this->nonzeros[from] * in[ii]; // diagonal
		//		TDKP_BOUNDS_ASSERT(from >= 0 && from < this->num_nz && until >= 0 && until <= this->num_nz, "from >= 0 && from < this->size && until >= 0 && until <= this->size");
				for(int jj = from + 1; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - fidx)];
					out[(this->icol[jj] - fidx)] += this->nonzeros[jj] * in[ii];
				}
			}
		} else {
			// perform nonsymmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - this->fidx;
				from     = this->prow[ii] - this->fidx;
				out[ii]  = 0.0;
				for(int jj = from; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - this->fidx)];
				}
			}
		}
	} else {
		TDKP_GENERAL_EXCEPTION("can not multiply with noninitialized matrix");
	}
}

template<>template<>    
void CSRMatrix<double>::mult_vec<cplx>(const cplx* in, cplx* out) const;

template<class T> template<class B>
void CSRMatrix<T>::mult_vec(const B* in, B* out) const {
	TDKP_GENERAL_EXCEPTION("not yet implemented");	
}

/** multiply 'multiple' vectors with matrix
 *
 * in FEM formulation of multiple partial differential equations on a grid
 * as the kp stuff here, the rhs matrix depends on the mesh and leads to
 * a 'block' sparse matrix, where each block is a diagonal identity matrix.
 * this leads to an independence between the single solutions, therefore allowing
 * only to save this part and to perform the multiplication at the same time
 *
 * now, suppose that the matrix B is a (n * kpn) x (n * kpn) matrix, where
 * kpn denotes the number of partial differential equations solved at the same time
 *
 * @param in input vector of length n * kpn, values stored as [x1_0, x1_1, ..., x1_kpn, x2_0, x2_1, ...]
 * @param out output vector of length n * kpn
 * @param kpn number of partial differential equataions > 1
 */
template<class T> template<class B>
void CSRMatrix<T>::mult_vec_multiple(const B* in, B* out, int kpn) const {
	/*
	if(kpn <= 1) {
		TDKP_GENERAL_EXCEPTION("kpn must be greater than one");
	}
	*/

	int msize    = (signed)this->size;
	int vec_size = msize * kpn;

	if(this->prow && this->icol	&& this->nonzeros) {
		int until, from, off;
		// --------------------------------------------------
		// distinguish between symmetric and nonsymmetric
		// --------------------------------------------------
		if(this->is_symmetric) {
			// clean up out
			for(int ii = 0; ii < vec_size; ii++) {
				out[ii] = (T)0.0;
			}
			// perform symmetric matrix vector product
			for(int ii = 0; ii < msize; ii++) {
				until    = this->prow[ii + 1] - fidx;
				from     = this->prow[ii] - fidx;
				off      = ii * kpn;
				for(int nn = 0; nn < kpn; nn++) {
					out[off + nn] += this->nonzeros[from] * in[off + nn]; // diagonal
				}
				TDKP_BOUNDS_ASSERT(from >= 0 && from < this->num_nz && until >= 0 && until <= this->num_nz, "from >= 0 && from < this->size && until >= 0 && until <= this->size");
				for(int jj = from + 1; jj < until; jj++) {
					for(int nn = 0; nn < kpn; nn++) {
						out[off + nn] += this->nonzeros[jj] * in[(this->icol[jj] - fidx) * kpn + nn];
						out[(this->icol[jj] - fidx) * kpn + nn] += this->nonzeros[jj] * in[off + nn];
					}
				}
			}
		} else {
			B ss[MAX_NUM_EQUATIONS];
			TDKP_ASSERT(kpn <= MAX_NUM_EQUATIONS, "kpn <= MAX_NUM_EQUATIONS");

			// perform nonsymmetric matrix vector product
			int* copy_prow		= this->prow;
			int* copy_icol		= this->icol;
			T*   copy_nonzeros	= this->nonzeros;
#pragma omp parallel for default(none) private(until,from,ss) shared(copy_prow,copy_icol,copy_nonzeros,out,kpn,in) schedule(static,5000)
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = copy_prow[ii + 1] - this->fidx;
				from     = copy_prow[ii] - this->fidx;
				for(int nn = 0; nn < kpn; nn++) {
					ss[nn] = 0.0;
//					out[ii * kpn + nn] = 0.0;
				}
				for(int jj = from; jj < until; jj++) {
					for(int nn = 0; nn < kpn; nn++) {
						ss[nn] += copy_nonzeros[jj] * in[(copy_icol[jj] - this->fidx) * kpn + nn];
//						out[ii * kpn + nn] += this->nonzeros[jj] * in[(this->icol[jj] - this->fidx) * kpn + nn];
					}
				}
				for(int nn = 0; nn < kpn; nn++) {
					out[ii * kpn + nn]  = ss[nn];
				}
			}
		}
	} else {
		TDKP_GENERAL_EXCEPTION("can not multiply with noninitialized matrix");
	}

}



/** read matrix from file
 *
 * @param filename of matrix which was stored by save_to_file
 */
template<class T>
void CSRMatrix<T>::read_from_file(const char* filename) {
	if(this->prow || this->icol || this->nonzeros) {
		this->clear_all();
	}
	std::fstream fin(filename, std::ios::in);
	if(fin) {
		int test;
		fin >> this->size;
		fin >> this->num_nz;
		fin >> this->is_symmetric;
		fin >> this->is_factorized;
		if(this->size <= 0) {
			TDKP_GENERAL_EXCEPTION("invalid matrix size (<=0)");
		}
		if(this->num_nz <= 0) {
			TDKP_GENERAL_EXCEPTION("invalid number of nonzeros (<=0)");
		}
		this->prow     = new int[this->size + 1];
		this->icol     = new int[this->num_nz];
		this->nonzeros = new   T[this->num_nz];
		TDKP_ASSERT(this->prow && this->icol && this->nonzeros, "this->prow && this->icol && this->nonzeros");
		for(unsigned int ii = 0; ii < size + 1; ii++) {
			fin >> this->prow[ii];
		}
		for(int ii = 0; ii < num_nz; ii++) {
			fin >> this->icol[ii];
		}
		for(int ii = 0; ii < num_nz; ii++) {
			fin >> this->nonzeros[ii];
		}
		// should not be at end of file
		if(fin.eof()) {
			TDKP_GENERAL_EXCEPTION("file should end with 101 ... but does not");
		}
		fin >> test;
		if(test != 101) {
			TDKP_GENERAL_EXCEPTION("file should end with 101 ... but does not");
		}
		fin.close();
		this->structure_set = true;
		return;
	} else {
		TDKP_FILE_EXCEPTION("can not open file", filename);
	}
}

/** test if corresponding structure was cached to file */
template<class T>
bool CSRMatrix<T>::structure_cache_available(const string& identifier) const {
	ostringstream sout;
	// build new identifier + size + symmetric
	sout << "cache_" << identifier << "_" << size << "_"
	     << symmetric() << ".dat";
	// replace /  by _
	string tmp = sout.str();
	for(unsigned int ii = 0; ii < tmp.size(); ii++) {
		if(tmp[ii] == '/') {
			tmp[ii] = '_';	
		}	
	}
	ifstream fin(tmp.c_str());
	if(fin) {
		fin.close();
		return true;
	} else {
		return false;
	}
}

/** load structure from file cache */
template<class T>
void CSRMatrix<T>::load_structure(const string& identifier) {
	if(this->prow || this->icol || this->nonzeros) {
		this->clear_all();
	}
	if(structure_cache_available(identifier)) {
		ostringstream sout;
		// build new identifier + size + symmetric
		sout << "cache_" << identifier << "_" << size << "_"
	    	 << symmetric() << ".dat";
		// replace /  by _
		string tmp = sout.str();
		for(unsigned int ii = 0; ii < tmp.size(); ii++) {
			if(tmp[ii] == '/') {
				tmp[ii] = '_';	
			}	
		}	    	 
		ifstream fin(tmp.c_str(), ios::binary);
		if(fin) {
			fin.read((char*)&num_nz, sizeof(int));
			TDKP_ASSERT(num_nz > (signed)size, "num_nz > size");
			TDKP_POINTER_ASSERT(prow     = new int[size + 1])
			TDKP_POINTER_ASSERT(icol     = new int[num_nz]);
			TDKP_POINTER_ASSERT(nonzeros = new T[num_nz]);
			fin.read((char*)prow, (size + 1) * sizeof(int));
			fin.read((char*)icol, (num_nz) * sizeof(int));
			int tmp;
			fin.read((char*)&tmp, sizeof(int));
			TDKP_ASSERT(tmp == 7041978, "tmp == 7041978");
			fin.close();
			for(int ii = 0; ii < num_nz; ii++) {
				nonzeros[ii] = 0;
			}
		} else {
			TDKP_GENERAL_EXCEPTION("can not read structure from " << sout.str());
		}

	} else {
		TDKP_GENERAL_EXCEPTION("structure cache is not available");
	}
}
template<class T>
bool CSRMatrix<T>::caching_structure_possible() const {
	return true;
}
template<class T>
void CSRMatrix<T>::save_structure(const string& identifier) {
	TDKP_ASSERT(prow != 0 && icol != 0, "prow != 0 && icol != 0");
	ostringstream sout;
	// build new identifier + size + symmetric
	sout << "cache_" << identifier << "_" << size << "_"
	     << symmetric() << ".dat";
	// replace /  by _
	string tmp = sout.str();
	for(unsigned int ii = 0; ii < tmp.size(); ii++) {
		if(tmp[ii] == '/') {
			tmp[ii] = '_';	
		}	
	}	     
	ofstream fout(tmp.c_str(), ios::binary);
	if(fout) {
		int tmp = 7041978;
		fout.write((char*)&num_nz, sizeof(int));
		fout.write((char*)prow, (size + 1) * sizeof(int));
		fout.write((char*)icol, (num_nz) * sizeof(int));
		fout.write((char*)&tmp, sizeof(int));
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write structure to " << sout.str());
	}
}




} // end namespace

#endif /*SYMCSRMATRIX_H_*/
