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

#include "tdkp/main/CSRSparseMatrixInterface.h"

namespace tdkp {


/** real valued representation of a complex matrix
 *
 * so, if A = R + iI then the real valued representation is
 * [R  -I]
 * [I   R]
 *
 */
C2DCSRSparseMatrixInterface::C2DCSRSparseMatrixInterface(CSRMatrix<double>& base_matrix)
: matrix(base_matrix),
  original_size(base_matrix.get_size() / 2)
{
	TDKP_ASSERT((signed)base_matrix.get_size() == this->original_size * 2, "base matrix should be a real representation of a complex matrix ... ");
	TDKP_ASSERT(!base_matrix.symmetric(), "real representation of complex matrices to the present only accepts nonsymmetric matrices");
}

C2DCSRSparseMatrixInterface::~C2DCSRSparseMatrixInterface() {
	// nothing to do ...
}

bool C2DCSRSparseMatrixInterface::property_is_set(SparseMatrixProperties flag) const {
	if(flag == symmetric_matrix && matrix.symmetric()) {
		return true;
	}
	if(flag == nonsymmetric_matrix && !matrix.symmetric()) {
		return true;
	}
	return false;
}
/** announce that element at (ii,jj) will be nonzero */
void C2DCSRSparseMatrixInterface::announce(int ii, int jj) {

	TDKP_BOUNDS_ASSERT(ii >= 0 && ii < this->original_size && jj >= 0 && jj < this->original_size, "indexes out of range");

	// announce upper triangle any way
	// if the matrix is symmetric and ii > jj, then announce here will
	// throw an exception ...
	matrix.announce(ii,jj);
	matrix.announce(ii,jj + original_size);
	matrix.announce(ii + original_size, jj + original_size);
	if(!matrix.symmetric()) {
		matrix.announce(ii + original_size, jj);
	}
}

void C2DCSRSparseMatrixInterface::block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num) {
	// announce upper triangle any way
	// if the matrix is symmetric and ii > jj, then announce here will
	// throw an exception ...
	int global_original_size = original_size / get_block_size();
	matrix.block_announce(gii,gjj, sparse_pat, sparse_num);
	matrix.block_announce(gii,gjj + global_original_size, sparse_pat, sparse_num);
	matrix.block_announce(gii + global_original_size, gjj + global_original_size, sparse_pat, sparse_num);
	if(!matrix.symmetric()) {
		matrix.block_announce(gii + global_original_size, gjj, sparse_pat, sparse_num);
	}
}

void C2DCSRSparseMatrixInterface::set_structure() {
	matrix.set_structure();
}
void C2DCSRSparseMatrixInterface::clear_but_keep_structure() {
	matrix.clear_but_keep_structure();
}
void C2DCSRSparseMatrixInterface::reset() {
	matrix.reset();
}
void C2DCSRSparseMatrixInterface::set_block_size(unsigned int block_size_) {
	matrix.set_block_size(block_size_);
}
unsigned int C2DCSRSparseMatrixInterface::get_block_size() const {
	return matrix.get_block_size();
}
void C2DCSRSparseMatrixInterface::set(int ii, int jj, const cplx& val) {
	if(val.real() != 0.0) {
		matrix.set(ii,jj, val.real());
		matrix.set(ii + original_size, jj + original_size, val.real());
	}
	if(val.imag() != 0.0) {
		matrix.set(ii,jj + original_size, - val.imag());
		if(!matrix.symmetric()) {
			matrix.set(ii + original_size, jj, val.imag());
		}
	}
}

void C2DCSRSparseMatrixInterface::set_row_and_column_to_zero(int ii) {
	matrix.set_row_and_column_to_zero(ii);
	matrix.set_row_and_column_to_zero(original_size + ii);		
}

void C2DCSRSparseMatrixInterface::add(int ii, int jj, const cplx& val) {
	if(val.real() != 0.0) {
		matrix.add(ii,jj, val.real());
		matrix.add(ii + original_size, jj + original_size, val.real());
	}
	if(val.imag() != 0.0) {
		matrix.add(ii,jj + original_size, - val.imag());
		if(!matrix.symmetric()) {
			matrix.add(ii + original_size, jj, val.imag());
		}
	}
}
cplx C2DCSRSparseMatrixInterface::get(int ii, int jj) const {
	return cplx(matrix.get(ii,jj), - matrix.get(ii, jj + original_size));
}
unsigned int C2DCSRSparseMatrixInterface::get_size() const {
	return this->original_size;
}

void C2DCSRSparseMatrixInterface::mult_vec(const cplx* in, cplx* out) const {

	vector<double> temporary_in(this->matrix.get_size());
	vector<double> temporary_out(this->matrix.get_size());

	// expand to real representation ...
	for(int ii = 0; ii < this->original_size; ii++) {
		temporary_in[ii] = in[ii].real();
		temporary_in[ii + this->original_size] = in[ii].imag();
	}
	this->matrix.mult_vec(&temporary_in[0], &temporary_out[0]);
	// back to complex representation ...
	for(int ii = 0; ii < this->original_size; ii++) {
		out[ii] = cplx(temporary_out[ii],temporary_out[ii + this->original_size]);
	}
}

/** return a vector with the abs value of nonhermitian terms per row */
const vector<double>& C2DCSRSparseMatrixInterface::get_nonhermitian_rows() const {
	return nonhermitian_rows;
}

void C2DCSRSparseMatrixInterface::perform_symmetry_analysis() {
	matrix.perform_symmetry_analysis();
	// reduce from real to complex representation
	const vector<double>& tmp = matrix.get_nonhermitian_rows();
	if(tmp.size() == matrix.get_size()) {
		nonhermitian_rows.assign(this->get_size(), 0.0);
		TDKP_ASSERT(this->get_size() * 2 == matrix.get_size(), "this->get_size() * 2 == matrix.get_size()");
		for(unsigned int ii = 0; ii < this->get_size(); ii++) {
			nonhermitian_rows[ii] += tmp[ii] + tmp[ii + this->get_size()];
		}
	}
}

void C2DCSRSparseMatrixInterface::delete_nonhermitian_rows_data() {
	matrix.delete_nonhermitian_rows_data();
	nonhermitian_rows.resize(0);
}

bool C2DCSRSparseMatrixInterface::structure_cache_available(const string& identifier) const {
	ostringstream sout;
	sout << identifier << "_c2dcsr";
	return matrix.structure_cache_available(sout.str());	
}

void C2DCSRSparseMatrixInterface::load_structure(const string& identifier) {
	ostringstream sout;
	sout << identifier << "_c2dcsr";
	matrix.load_structure(sout.str());
}
void C2DCSRSparseMatrixInterface::save_structure(const string& identifier) {
	ostringstream sout;
	sout << identifier << "_c2dcsr";
	matrix.save_structure(sout.str());
}

// ---------------------------------------------------------------
// implementation of C2LocalDCSRSparseMatrixInterface
// ---------------------------------------------------------------
C2LocalDCSRSparseMatrixInterface::C2LocalDCSRSparseMatrixInterface(CSRMatrix<double>& base_matrix)
: C2DCSRSparseMatrixInterface(base_matrix)
{
}


/** announce that element at (ii,jj) will be nonzero */
void C2LocalDCSRSparseMatrixInterface::announce(int ii, int jj) {

	TDKP_BOUNDS_ASSERT(ii >= 0 && ii < this->original_size && jj >= 0 && jj < this->original_size, "indexes out of range");

	// announce upper triangle any way
	// if the matrix is symmetric and ii > jj, then announce here will
	// throw an exception ...
	matrix.announce(ii * 2,     jj * 2);
	matrix.announce(ii * 2,     jj * 2 + 1);
	matrix.announce(ii * 2 + 1, jj * 2 + 1);
	if(!matrix.symmetric() || ii < jj) {
		matrix.announce(ii * 2 + 1, jj * 2);
	}
}

void C2LocalDCSRSparseMatrixInterface::block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num) {
		
	// announce upper triangle any way
	// if the matrix is symmetric and ii > jj, then announce here will
	// throw an exception ...
	const int block_size = matrix.get_block_size();
	if(matrix.symmetric() && gii == gjj) {
		for(int ss = 0; ss < sparse_num; ss++) {
			if(gii + sparse_pat[2 * ss] <= gjj + sparse_pat[2 * ss + 1]) {
				this->announce(gii * block_size + sparse_pat[2 * ss], gjj * block_size + sparse_pat[2 * ss]);
			}
		}
	} else {
		for(int ss = 0; ss < sparse_num; ss++) {
			this->announce(gii * block_size + sparse_pat[2 * ss], gjj * block_size + sparse_pat[2 * ss + 1]);
		}
	}
}
void C2LocalDCSRSparseMatrixInterface::set(int ii, int jj, const cplx& val) {
	if(val.real() != 0.0) {
		matrix.set(ii * 2,jj * 2, val.real());
		matrix.set(ii * 2 + 1, jj * 2 + 1, val.real());
	}
	if(val.imag() != 0.0) {
		matrix.set(ii * 2,jj * 2 + 1, - val.imag());
		if(ii == jj) {
			if(!matrix.symmetric()) {
				matrix.set(ii * 2 + 1, jj * 2, val.imag());
			}
		} else {
			matrix.set(ii * 2 + 1, jj * 2, val.imag());
		}
	}
}

void C2LocalDCSRSparseMatrixInterface::set_row_and_column_to_zero(int ii) {
	matrix.set_row_and_column_to_zero(2*ii);
	matrix.set_row_and_column_to_zero(2*ii + 1);		
}

void C2LocalDCSRSparseMatrixInterface::add(int ii, int jj, const cplx& val) {
	if(val.real() != 0.0) {
		matrix.add(ii * 2,jj * 2, val.real());
		matrix.add(ii * 2 + 1, jj * 2 + 1, val.real());
	}
	if(val.imag() != 0.0) {
		matrix.add(ii * 2,jj * 2 + 1, - val.imag());
		if(ii == jj) {
			if(!matrix.symmetric()) {
				matrix.add(ii * 2 + 1, jj * 2, val.imag());
			}
		} else {
			matrix.add(ii * 2 + 1, jj * 2, val.imag());
		}
	}
}

cplx C2LocalDCSRSparseMatrixInterface::get(int ii, int jj) const {
	return cplx(matrix.get(ii * 2,jj * 2), - matrix.get(ii * 2, jj * 2 + 1));
}



void C2LocalDCSRSparseMatrixInterface::mult_vec(const cplx* in, cplx* out) const {

	TDKP_ASSERT((signed)this->matrix.get_size() == this->original_size * 2, "this->matrix.get_size() == this->original_size * 2");
	vector<double> temporary_in(matrix.get_size());
	vector<double> temporary_out(matrix.get_size());	
	// expand to real representation ...
	for(int ii = 0; ii < this->original_size; ii++) {
		temporary_in[ii * 2] = in[ii].real();
		temporary_in[ii * 2 + 1] = in[ii].imag();
	}
	this->matrix.mult_vec(&temporary_in[0], &temporary_out[0]);
	// back to complex representation ...
	for(int ii = 0; ii < this->original_size; ii++) {
		out[ii] = cplx(temporary_out[ii * 2],temporary_out[ii * 2 + 1]);
	}
}

bool C2LocalDCSRSparseMatrixInterface::structure_cache_available(const string& identifier) const {
	ostringstream sout;
	sout << identifier << "_c2localdcsr";
	return matrix.structure_cache_available(sout.str());	
}
void C2LocalDCSRSparseMatrixInterface::load_structure(const string& identifier) {
	ostringstream sout;
	sout << identifier << "_c2localdcsr";
	matrix.load_structure(sout.str());
}
void C2LocalDCSRSparseMatrixInterface::save_structure(const string& identifier) {
	ostringstream sout;
	sout << identifier << "_c2localdcsr";
	matrix.save_structure(sout.str());
}

// ----------------------------------------------------------------
// implementation of double 2 complex matrix
// ----------------------------------------------------------------


D2CCSRSparseMatrixInterface::D2CCSRSparseMatrixInterface(SparseMatrixInterface<cplx>& base_matrix)
: matrix(base_matrix)
{

}

D2CCSRSparseMatrixInterface::~D2CCSRSparseMatrixInterface() {
	// nothing to do
}

bool D2CCSRSparseMatrixInterface::property_is_set(SparseMatrixProperties flag) const {
	return matrix.property_is_set(flag);
}

void D2CCSRSparseMatrixInterface::announce(int ii, int jj) {
	matrix.announce(ii,jj);
}
void D2CCSRSparseMatrixInterface::set_structure() {
	matrix.set_structure();
}
void D2CCSRSparseMatrixInterface::clear_but_keep_structure() {
	matrix.clear_but_keep_structure();
}
void D2CCSRSparseMatrixInterface::reset() {
	matrix.reset();
}
void D2CCSRSparseMatrixInterface::set_block_size(unsigned int block_size_) {
	matrix.set_block_size(block_size_);
}
unsigned int D2CCSRSparseMatrixInterface::get_block_size() const {
	return matrix.get_block_size();
}
void D2CCSRSparseMatrixInterface::set(int ii, int jj, const double& val) {
	matrix.set(ii,jj,cplx(val,0.0));
}
void D2CCSRSparseMatrixInterface::set_row_and_column_to_zero(int ii) {
	matrix.set_row_and_column_to_zero(ii);		
}

void D2CCSRSparseMatrixInterface::add(int ii, int jj, const double& val) {
	matrix.add(ii,jj,cplx(val,0.0));
}
double D2CCSRSparseMatrixInterface::get(int ii, int jj) const {
	return matrix.get(ii,jj).real();
}
unsigned int D2CCSRSparseMatrixInterface::get_size() const {
	return matrix.get_size();
}

void D2CCSRSparseMatrixInterface::mult_vec(const double* in, double* out) const {

	vector<cplx> temporary_in(this->matrix.get_size());
	vector<cplx> temporary_out(this->matrix.get_size());

	unsigned int size = matrix.get_size();

	// connode to complex representation
	for(unsigned int ii = 0; ii < size; ii++) {
		temporary_in[ii] = in[ii];
	}
	this->matrix.mult_vec(&temporary_in[0], &temporary_out[0]);
	// back to complex representation ...
	for(unsigned int ii = 0; ii < size; ii++) {
		out[ii] = temporary_out[ii].real();
		TDKP_BOUNDS_ASSERT(temporary_out[ii].imag() == 0.0, "temporary_out[ii].imag() == 0.0");
	}
}
void D2CCSRSparseMatrixInterface::save_to_file(const char* filename) const {
	matrix.save_to_file(filename);	
}

void D2CCSRSparseMatrixInterface::perform_symmetry_analysis() {
	matrix.perform_symmetry_analysis();	
}


// ----------------------------------------------------------
// TwinComplexCSR
// ----------------------------------------------------------
TwinComplexCSR::TwinComplexCSR(unsigned int size)
: real_csr(size, nonsymmetric_matrix, false),
  imag_array(0),
  real_array(0),
  prow(0),
  icol(0),
  current_structure(0)
{

}
TwinComplexCSR::~TwinComplexCSR() {}
bool TwinComplexCSR::property_is_set(SparseMatrixProperties flag) const {
	return real_csr.property_is_set(flag);
}
void TwinComplexCSR::announce(int ii, int jj) {
	real_csr.announce(ii,jj);
}

void TwinComplexCSR::block_announce(int gii, int gjj, const int* sparse_pat, int sparse_num) {
	real_csr.block_announce(gii,gjj,sparse_pat,sparse_num);
}

void TwinComplexCSR::set_structure() {
	real_csr.set_structure();
	real_array = real_csr.get_nonzeros();
	imag_array.assign(real_csr.get_num_nonzeros(), 0.0e0);
	prow = real_csr.get_prow();
	icol = real_csr.get_icol();
	current_structure++;
}

	
void TwinComplexCSR::init_from_csr(bool symmetric, int size_, int num_nz_, int* prow_, int* icol_, cplx* nonzeros_) {
	TDKP_ASSERT(!symmetric, "TwinComplexCSR does not work with symmetric matrices");
	TDKP_ASSERT(prow_[0] == 0, "no fortran indices allowed");
	vector<double> real_tmp(num_nz_);
	imag_array.resize(num_nz_);
	for(unsigned int ii = 0; ii < real_tmp.size(); ii++) {
		real_tmp[ii]   = nonzeros_[ii].real();
		imag_array[ii] = nonzeros_[ii].imag();
	}
	real_csr.init_from_csr(false, size_, num_nz_, prow_, icol_, &real_tmp[0]);
	real_array = real_csr.get_nonzeros();
	prow       = real_csr.get_prow();
	icol       = real_csr.get_icol();
	current_structure++;		
}

void TwinComplexCSR::clear_but_keep_structure() {
	real_csr.clear_but_keep_structure();
	imag_array.assign(real_csr.get_num_nonzeros(), 0.0e0);
	prow = real_csr.get_prow();
	icol = real_csr.get_icol();
}

void TwinComplexCSR::reset() {
	real_csr.reset();
	imag_array.clear();
	real_array = 0;
	prow = 0;
	icol = 0;
}
void TwinComplexCSR::set_block_size(unsigned int block_size_) {
	real_csr.set_block_size(block_size_);
}
unsigned int TwinComplexCSR::get_block_size() const {
	return real_csr.get_block_size();
}
void TwinComplexCSR::set(int ii, int jj, const cplx& val) {
	int pos = real_csr.pos(ii,jj);
	TDKP_BOUNDS_ASSERT(pos >= 0, "(ii,jj) does not exist!");
	TDKP_BOUNDS_ASSERT(real_array != 0, "real array is not assigned!");
	real_array[pos] = val.real();
	imag_array[pos] = val.imag();
}

void TwinComplexCSR::set_row_and_column_to_zero(int ii) {
	TDKP_GENERAL_EXCEPTION("sorry, but TwinComplexCSR::set_row_and_column_to_zero is not implemented.");	
}

void TwinComplexCSR::add(int ii, int jj, const cplx& val) {
	int pos = real_csr.pos(ii,jj);
	TDKP_BOUNDS_ASSERT(pos >= 0, "(ii,jj) does not exist!");
	TDKP_BOUNDS_ASSERT(real_array != 0, "real array is not assigned!");
	real_array[pos] += val.real();
	imag_array[pos] += val.imag();
}
cplx TwinComplexCSR::get(int ii, int jj) const {
	int pos = real_csr.pos(ii,jj);
	if(pos == -1) {
		return cplx(0.0, 0.0);
	} else {
		return cplx(real_array[pos], imag_array[pos]);
	}
}
unsigned int TwinComplexCSR::get_size() const {
	return real_csr.get_size();
}

void TwinComplexCSR::mult_vec(const cplx* in, cplx* out) const {

	TDKP_ASSERT(!real_csr.symmetric(), "matrix must be nonsymmetric!");

	const int  fidx = real_csr.get_fidx();
	const int  size = real_csr.get_size();

	int from, until;

	// perform nonsymmetric matrix vector product
#pragma omp parallel for default(shared) private(until, from)
	for(int ii = 0; ii < size; ii++) {
		until    = prow[ii + 1] - fidx;
		from     = prow[ii] - fidx;
		out[ii]  = 0.0;
		for(int jj = from; jj < until; jj++) {
			out[ii] += cplx(real_array[jj],imag_array[jj]) * in[(icol[jj] - fidx)];
		}
	}
}

void TwinComplexCSR::perform_symmetry_analysis() {
	Logger::get_instance()->emit(LOG_WARN, "performing symmetry analysis is not supported for the current matrix type!");
}

void TwinComplexCSR::save_to_file(const char* filename) const {
	
	unsigned int size = this->get_size();
	bool symmetric    = this->property_is_set(symmetric_matrix);
	int  num_nz       = this->get_number_of_nonzeros();
	bool factorized   = false;
	
	const int* prow   = this->get_prow();
	const int* icol   = this->get_icol();
	const double* nz_real = get_nz_real();
	const double* nz_imag = get_nz_imag();
	
	std::fstream fout(filename, std::ios::out);
	if(fout) {		
		fout << size << "  " << num_nz << "  "
			 << symmetric << "  " << factorized << "\n";
		for(unsigned int ii = 0; ii < size + 1; ii++) {
			fout << prow[ii] << "\n";
		}
		for(int ii = 0; ii < num_nz; ii++) {
			fout << icol[ii] << "\n";
		}
		fout.precision(20);
		for(int ii = 0; ii < num_nz; ii++) {
			fout << cplx(nz_real[ii],nz_imag[ii]) << "\n";
		}
		fout << 101 << "\n";
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to file " << filename);	
	}
}



} // end of namespace
