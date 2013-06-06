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

#ifndef C2DCSRSPARSEMATRIXINTERFACE_H_
#define C2DCSRSPARSEMATRIXINTERFACE_H_

#include "tdkp/common/all.h"
#include "tdkp/main/SparseMatrixInterface.h"
#include "tdkp/main/CSRMatrix.h"

namespace tdkp {

/** real representation of complex matrix */
class C2DCSRSparseMatrixInterface : public SparseMatrixInterface<cplx> {
public:
	C2DCSRSparseMatrixInterface(CSRMatrix<double>& base_matrix);
	virtual ~C2DCSRSparseMatrixInterface();
	virtual bool			property_is_set(SparseMatrixProperties flag) const;
	virtual void			announce(int ii, int jj);
	virtual void 			block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num);
	virtual void 			set_structure();
	virtual void 			clear_but_keep_structure();
	virtual void 			reset();
	virtual void 			set_block_size(unsigned int block_size_);
	virtual unsigned int	get_block_size() const;
	virtual void     		set(int ii, int jj, const cplx& val);
	virtual void            set_row_and_column_to_zero(int ii);
	virtual void     		add(int ii, int jj, const cplx& val);
	virtual cplx  			get(int ii, int jj) const;
	virtual unsigned int 	get_size() const;
	virtual void perform_symmetry_analysis();
	virtual const vector<double>& get_nonhermitian_rows() const;
	virtual void  delete_nonhermitian_rows_data();

	virtual void mult_vec(const cplx* in, cplx* out) const;

	virtual bool caching_structure_possible() const { return matrix.caching_structure_possible(); }
	virtual bool structure_cache_available(const string& identifier) const;
	virtual void load_structure(const string& identifier);
	virtual void save_structure(const string& identifier);

protected:
	CSRMatrix<double>& matrix;
	int original_size;
	vector<double> nonhermitian_rows;

};

/** locally structured real representation of complex matrix
 *
 * instead or [R -I;I R] globally, the real and imaginary part of a number
 * is stored locally
 */
class C2LocalDCSRSparseMatrixInterface : public C2DCSRSparseMatrixInterface {
public:
	C2LocalDCSRSparseMatrixInterface(CSRMatrix<double>& base_matrix);
	virtual ~C2LocalDCSRSparseMatrixInterface() {}

	virtual void			announce(int ii, int jj);
	virtual void 			block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num);
	virtual void     		set(int ii, int jj, const cplx& val);
	virtual void            set_row_and_column_to_zero(int ii);
	virtual void     		add(int ii, int jj, const cplx& val);
	virtual cplx  			get(int ii, int jj) const;
	virtual void mult_vec(const cplx* in, cplx* out) const;
	virtual bool structure_cache_available(const string& identifier) const;
	virtual void load_structure(const string& identifier);
	virtual void save_structure(const string& identifier);
private:
		
};

/* real matrix stored as a complex one ... */
class D2CCSRSparseMatrixInterface : public SparseMatrixInterface<double> {

public:
	D2CCSRSparseMatrixInterface(SparseMatrixInterface<cplx>& base_matrix);
	virtual ~D2CCSRSparseMatrixInterface();
	virtual bool			property_is_set(SparseMatrixProperties flag) const;
	virtual void			announce(int ii, int jj);
	virtual void 			set_structure();
	virtual void 			clear_but_keep_structure();
	virtual void 			reset();
	virtual void 			set_block_size(unsigned int block_size_);
	virtual unsigned int	get_block_size() const;
	virtual void     		set(int ii, int jj, const double& val);
	virtual void            set_row_and_column_to_zero(int ii);
	virtual void     		add(int ii, int jj, const double& val);
	virtual double 			get(int ii, int jj) const;
	virtual unsigned int 	get_size() const;
	virtual void save_to_file(const char* filename) const;
	virtual void  perform_symmetry_analysis();

	virtual void mult_vec(const double* in, double* out) const;

protected:
	SparseMatrixInterface<cplx>& matrix;
};

/** complex matrix, storing real and imaginary part in separate double arrays
 */
class TwinComplexCSR : public SparseMatrixInterface<cplx> {
public:
	TwinComplexCSR(unsigned int size);
	virtual ~TwinComplexCSR();
	virtual bool init_from_csr_available() const { return true; }
	virtual void init_from_csr(bool symmetric, int size, int num_nz, int* prow, int* icol, cplx* nonzeros);
	virtual bool			property_is_set(SparseMatrixProperties flag) const;
	virtual void			announce(int ii, int jj);
	virtual void 			block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num);
	virtual void 			set_structure();
	virtual void 			clear_but_keep_structure();
	virtual void 			reset();
	virtual void 			set_block_size(unsigned int block_size_);
	virtual unsigned int	get_block_size() const;
	virtual void     		set(int ii, int jj, const cplx& val);
	virtual void     		add(int ii, int jj, const cplx& val);
	virtual cplx 			get(int ii, int jj) const;
	virtual void            set_row_and_column_to_zero(int ii);
	virtual unsigned int 	get_size() const;
	virtual void            perform_symmetry_analysis();
	

	virtual void mult_vec(const cplx* in, cplx* out) const;

	const int* 			get_prow() const { TDKP_ASSERT(prow != 0, "matrix not assembled yet!"); return prow; }
	const int* 			get_icol() const { TDKP_ASSERT(icol != 0, "matrix not assembled yet!"); return icol; }
	const double*       get_nz_real() const { TDKP_ASSERT(real_array != 0, "matrix not assembled yet"); return real_array; }
	const double*       get_nz_imag() const { TDKP_ASSERT(imag_array.size() >= 0, "matrix not assembled yet!"); return &imag_array[0]; }
	unsigned int        get_number_of_nonzeros() const { return real_csr.get_num_nonzeros(); }

	/** saving matrix to file */
	void                save_to_file(const char* filename) const; 

	// return the structure count (if structure is the same, the umfpack conversion from csr to csc is cheaper ...)
	int 				get_current_structure() { return this->current_structure; }
	
	

private:
	CSRMatrix<double>  real_csr;
	vector<double>     imag_array;
	double*            real_array;
	// prow and icol are non const in CSR
	const int* 		   prow;
	const int* 		   icol;
	int                current_structure;

};


} // end of namespace tdkp
#endif /*C2DCSRSPARSEMATRIXINTERFACE_H_*/
