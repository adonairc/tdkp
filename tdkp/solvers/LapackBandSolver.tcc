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

template<class T>
FullMatrix<T>::FullMatrix(unsigned int size_)
: symmetric(symmetric_matrix),
  size(size_),
  block_size(1)
{
}
		
template<class T>
void FullMatrix<T>::clear_but_keep_structure() {
	mdata.assign(mdata.size(), 0.0e0);	
}

template<class T>
void FullMatrix<T>::reset() {
	mdata.clear();	
}

template<class T>
bool FullMatrix<T>::property_is_set(SparseMatrixProperties flag) const {
	if(flag == symmetric_matrix && symmetric == true) {
		return true;	
	} else if(flag == nonsymmetric_matrix && symmetric == false) {
		return true;	
	} else {
		return false;	
	}	
}	
								
template<class T>
void FullMatrix<T>::set_property(SparseMatrixProperties flag) {
	if(flag == symmetric_matrix) {
		symmetric = true;
	} else if(flag == nonsymmetric_matrix) {
		symmetric = false;	
	}		
}

template<class T>
void FullMatrix<T>::set_structure() {
	long memusage = (this->get_size() * this->get_size() * sizeof(T)) / 1024 / 1024;	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "FullMatrix: matrix (N = " << this->get_size() << ") memory usage is " << memusage << "M."); 		
	this->mdata.assign(this->get_size() * this->get_size(), 0.0e0);	
}

template<>
void FullMatrix<cplx>::set(int ii, int jj, const cplx& val);

template<class T>
void FullMatrix<T>::set(int ii, int jj, const T& val) {
	TDKP_BOUNDS_ASSERT(!((ii > jj) && symmetric), "!((ii(=" << ii << ") > jj(="<<jj<<") && symmetric(" << symmetric << "))");
	TDKP_BOUNDS_ASSERT(this->mdata.size() > jj * this->get_size() + ii, "");
	this->mdata[jj * this->get_size() + ii] = val;
	if(symmetric && ii != jj) {
		this->mdata[ii * this->get_size() + jj] = this->mdata[jj * this->get_size() + ii]; 		
	}	
}

template<class T>
void FullMatrix<T>::set_row_and_column_to_zero(int ii) {
	if(symmetric) {
		// set row to zero
		for(int jj = ii; jj < static_cast<int>(this->get_size()); jj++) {
			this->set(ii,jj,0.0);	
		}
		// set column to zero
		for(int jj = 0; jj <= ii; jj++) {
			this->set(jj,ii,0.0);
		}
	} else {
		// set row and column to zero
		for(int jj = 0; jj < static_cast<int>(this->get_size()); jj++) {
			this->set(ii,jj,0.0);
			this->set(jj,ii,0.0);	
		}
	}
}

template<>
void FullMatrix<cplx>::add(int ii, int jj, const cplx& val);

template<class T>
void FullMatrix<T>::add(int ii, int jj, const T& val) {
	TDKP_BOUNDS_ASSERT(!((ii > jj) && symmetric), "!((ii(=" << ii << ") > jj(="<<jj<<") && symmetric(" << symmetric << "))");
	TDKP_BOUNDS_ASSERT(this->mdata.size() > jj * this->get_size() + ii, "");
	this->mdata[jj * this->get_size() + ii] += val;
	if(symmetric && ii != jj) {
		this->mdata[ii * this->get_size() + jj] = this->mdata[jj * this->get_size() + ii]; 		
	}
}

template<class T>
T FullMatrix<T>::get(int ii, int jj) const {
	TDKP_BOUNDS_ASSERT(!((ii > jj) && symmetric), "!((ii(=" << ii << ") > jj(="<<jj<<") && symmetric(" << symmetric << "))");
	TDKP_BOUNDS_ASSERT(this->mdata.size() > jj * this->get_size() + ii, "");
	return this->mdata[jj * this->get_size() + ii];
}	

template<class T>
void FullMatrix<T>::mult_vec(const T* in, T* out) const {
	
	const int n = this->get_size();

#pragma omp parallel for default(shared) schedule(static, 100)	
	for(int ii = 0; ii < n; ii++) {
		out[ii] = 0;
		for(int jj = 0; jj < n; jj++) {
			out[ii] += this->mdata[jj * n + ii] * in[jj];	
		}	
	}
}

template<>
void FullMatrix<cplx>::save_to_file(const char* filename) const;

template<class T>
void FullMatrix<T>::save_to_file(const char* filename) const {
	const int n = this->get_size();
	ofstream fout(filename);
	if(fout) {	
		for(int ii = 0; ii < n; ii++) {
			for(int jj = 0; jj < n; jj++) {
				T val;
				if(symmetric && jj < ii) {
					val = get(jj,ii);
				} else {
					val = get(ii,jj);	
				}
				fout << val << "  ";
			}
			fout << "\n";	
		}
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to file " << filename);	
	}	
}
template<>
void FullMatrix<cplx>::perform_symmetry_analysis();

template<>
void FullMatrix<double>::perform_symmetry_analysis();

template<class T>
BandMatrix<T>::BandMatrix(unsigned int size_) 
: data(0),
  size(size_),
  band_width(0),
  data_length(0),
  block_size(1),
  structure_is_set(false),
  symmetric(false)
{
	
}

template<class T>
int BandMatrix<T>::get_idx(int ii,int jj) const {
	TDKP_BOUNDS_ASSERT(abs(ii - jj) <= (signed)band_width, "index (ii/jj) out of range!");	
	if(symmetric) {
		return (band_width + 1) * jj + band_width + (ii - jj);
	} else {
		return (2 * band_width + 1) * jj + band_width + (ii - jj);
	}	
}

template<class T>
BandMatrix<T>::~BandMatrix() {
	if(data != 0) {
		delete[] data;
		data = 0;	
	}	
}

template<class T>
bool BandMatrix<T>::property_is_set(SparseMatrixProperties flag) const {
	if(flag == symmetric_matrix && symmetric) {
		return true;	
	}
	if(flag == nonsymmetric_matrix && !symmetric) {
		return true;	
	}
	return false;
}

template<class T>
void BandMatrix<T>::announce(int ii, int jj) {
	TDKP_ASSERT(!structure_is_set, "you can not call announce on an initialized matrix!");
#pragma omp critical
	{ 
		band_width = max(abs(ii - jj), int(band_width));
	}	
}

template<class T>
void BandMatrix<T>::set_structure() {
	if(data != 0) {
		delete[] data;
		data = 0;	
	}
	if(symmetric) {
		data_length = size * (band_width + 1);		
	} else {
		data_length = size * (band_width * 2 + 1);
	}	
	data = new T[data_length];
	TDKP_POINTER_ASSERT(data);
	for(unsigned int ii = 0; ii < data_length; ii++) {
		data[ii] = 0.0;
	}
	structure_is_set = true;
	ostringstream sout;
	sout << "building ";
	if(symmetric) {
		sout << " symmetric";
	} else {
		sout << " nonsymmetric";	
	}	
	sout << " lapack band matrix of size " << size << " with "
	     << "band width of " << band_width;
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
}

template<class T>
void BandMatrix<T>::clear_but_keep_structure() {
	for(unsigned int ii = 0; ii < data_length; ii++) {
		data[ii] = 0.0;
	}
}


template<class T>
void BandMatrix<T>::reset() {
	data_length = 0; 
	band_width = 0;
	structure_is_set = false;
	delete[] data; data = 0;
		
}
template<class T>
void BandMatrix<T>::set_block_size(unsigned int block_size_) {
	block_size = block_size_;		
}
template<class T>
unsigned int BandMatrix<T>::get_block_size() const {
	return block_size;	
}
template<class T>
void BandMatrix<T>::set(int ii, int jj, const T& val) {
	int idx = get_idx(ii,jj);
	if(idx < 0 || idx >= (signed)data_length) {
		cout << "index out of range: " << idx << "  data length: " << data_length 
		     << " at ii: " << ii << "  jj: " << jj << "\n";
	}	
	TDKP_ASSERT(idx > -1 && idx < (signed)data_length, "(ii,jj) index out of range!: ");
	data[idx] = val;
}


template<class T>
void BandMatrix<T>::set_row_and_column_to_zero(int ii) {
	if(symmetric) {
		// set row to zero
		for(int jj = ii; jj < static_cast<int>(this->get_size()) && (jj - ii) < static_cast<int>(band_width); jj++) {			
			this->set(ii,jj,0.0);	
		}
		// set column to zero
		for(int jj = ii - static_cast<int>(band_width); jj <= ii; jj++) {
			if(jj >= 0) {
				this->set(jj,ii,0.0);
			}
		}
	} else {
		// set row and column to zero
		for(int jj = ii - band_width; jj < ii + static_cast<int>(band_width); jj++) {
			if(jj >= 0 && jj < static_cast<int>(this->get_size())) {
				this->set(ii,jj,0.0);
				this->set(jj,ii,0.0);
			}	
		}
	}
}

template<class T>
void BandMatrix<T>::add(int ii, int jj, const T& val) {
	int idx = get_idx(ii,jj);
	TDKP_ASSERT(idx > -1 && idx < (signed)data_length, "(ii,jj) index out of range!: ");
	data[idx] += val;	
}
template<class T>
T 	 BandMatrix<T>::get(int ii, int jj) const {
	int idx = get_idx(ii,jj);
	if(idx > -1 && idx < (signed)data_length) {
		return data[idx];	
	} else {
		return 0.0;	
	}	
}	
			
	
template<class T>
void BandMatrix<T>::set_property(SparseMatrixProperties property) {
	TDKP_ASSERT(!structure_is_set, "you may not alter the properties of the matrix after it is initialized!");
	if(property == symmetric_matrix) {
		symmetric = true;
	} else if(property == nonsymmetric_matrix) {
		symmetric = false;
	} else {
		TDKP_GENERAL_EXCEPTION("sorry, unhandled matrix property set: " << property);	
	}
}
