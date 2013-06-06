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

#ifndef SPARSEMATRIXINTERFACE_H_
#define SPARSEMATRIXINTERFACE_H_

#include <vector>
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"

namespace tdkp {

/** abstract matrix interface for the general sparse matrix
 *
 * so, in order to get an abstraction between the fem
 * solver and the actual implementation of the sparse matrix
 * we introduce this abstract interface which is then
 * implemented for possible matrix types
 *
 * which matrix should be used is defined by the requested solver
 * that we will generate via factory patterns
 *
 */
template<class T>
class SparseMatrixInterface {
public:
	SparseMatrixInterface() {};
	virtual ~SparseMatrixInterface() {};

	virtual bool			property_is_set(SparseMatrixProperties flag) const = 0;
	virtual void			announce(int ii, int jj) = 0;
	virtual void            block_announce(int ii, int jj, const int* sparse_pat, int sparse_num);
	virtual void 			set_structure() = 0;
	virtual void 			clear_but_keep_structure() = 0;
	virtual void 			reset() = 0;
	virtual void 			set_block_size(unsigned int block_size_) = 0;
	virtual unsigned int	get_block_size() const = 0;
	virtual void     		set(int ii, int jj, const T& val) = 0;
	virtual void     		add(int ii, int jj, const T& val) = 0;
	virtual T 				get(int ii, int jj) const = 0;
	virtual void 			set_row_and_column_to_zero(int ii) = 0; 
	virtual unsigned int 	get_size() const = 0;


	virtual void mult_vec(const T* in, T* out) const = 0;
	
	// --------------------------------------------------
	// init from crs data
	// --------------------------------------------------
	virtual bool init_from_csr_available() const;
	virtual void init_from_csr(bool symmetric, int size, int num_nz, int* prow, int* icol, T* nonzeros);

	// --------------------------------------------------
	// maintenance and analysis functions
	// --------------------------------------------------
	virtual void  perform_symmetry_analysis();
	virtual void  save_to_file(const char* filename) const ;
	virtual const vector<double>& get_nonhermitian_rows() const;
	virtual void  delete_nonhermitian_rows_data() {};

	// --------------------------------------------------
	// structure cache functions
	// --------------------------------------------------
	virtual bool caching_structure_possible() const { return false; }
	virtual bool structure_cache_available(const string& identifier) const { return false; }
	virtual void load_structure(const string& identifier) { TDKP_GENERAL_EXCEPTION("not implemented for base class"); }
	virtual void save_structure(const string& identifier) { TDKP_GENERAL_EXCEPTION("not implemented for base class"); }


private:
	vector<double> dummy_nonhermitian_nodes;
};

template<class T>
void SparseMatrixInterface<T>::save_to_file(const char* filename) const {
	Logger::get_instance()->emit(LOG_WARN, "saving matrix to file is not supported for the current matrix type!");
}

template<class T>
const vector<double>& SparseMatrixInterface<T>::get_nonhermitian_rows() const {
	Logger::get_instance()->emit(LOG_WARN, "tracking non-hermitian entries in the matrix is not supported for the current matrix type!");
	return dummy_nonhermitian_nodes;
}

template<class T>
void SparseMatrixInterface<T>::perform_symmetry_analysis() {
	Logger::get_instance()->emit(LOG_WARN, "performing symmetry analysis is not supported for the current matrix type!");
}

template<class T>
bool SparseMatrixInterface<T>::init_from_csr_available() const {
	return false;	
}
template<class T>
void SparseMatrixInterface<T>::init_from_csr(bool symmetric, int size, int num_nz, int* prow, int* icol, T* nonzeros) {
	TDKP_GENERAL_EXCEPTION("init from crs is not available for this matrix type	");
}

/** slow wrapper for block announce using standard announce */
template<class T>
void SparseMatrixInterface<T>::block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num) {
	int block_size = get_block_size();
	bool symmetric = property_is_set(symmetric_matrix);
	if(symmetric && gii == gjj) {
		int ii,jj;
		for(int ss = 0; ss < sparse_num; ss++) {
			ii = gii * block_size + sparse_pat[2*ss];
			jj = gjj * block_size + sparse_pat[2*ss + 1];
			if(!(ii < jj)) {
				this->announce(ii,jj);
			}
		}
	} else {
		for(int ss = 0; ss < sparse_num; ss++) {
			this->announce(gii * block_size + sparse_pat[2*ss], gjj * block_size + sparse_pat[2*ss + 1]);
		}
	}
}

} // end of namespace

#endif /*SPARSEMATRIXINTERFACE_H_*/
