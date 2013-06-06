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

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include "tdkp/common/all.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

namespace tdkp {

class Vector3D {
public:
	Vector3D();
	Vector3D(const Vector3D &copy);
    Vector3D(const double* init);
    Vector3D(double x, double y, double z);
    Vector3D(const double* a, const double* b, double x);
    double          get(unsigned int ii) const  { return v[ii]; }    
    void            set(unsigned int ii, double val) { v[ii] = val; }
	double&         operator()(unsigned int ii) { return v[ii]; }
	const double&   operator()(unsigned int ii) const { return v[ii]; }
	Vector3D        operator-(const Vector3D &rhs) const;
	Vector3D        operator+(const Vector3D &rhs) const;
	Vector3D        operator*(const double& rhs) const;
	Vector3D&       operator=(const Vector3D &rhs);		
	bool            operator==(const Vector3D &rhs) const;
	bool            operator!=(const Vector3D &rhs) const;
	double          norm() const;	
	void            normalize();
	const double*   get_all() const { return this->v; }	
	static double   dot_product(const Vector3D &x, const Vector3D &y);
	static Vector3D cross_product(const Vector3D &x, const Vector3D &y);
	
	double v[3];
};


Vector3D operator*(const double a, const Vector3D &rhs);
Vector3D operator-(const Vector3D& rhs);
std::ostream& operator<<(std::ostream& out, const Vector3D& vec);


/** a simple matrix class */
template<class T>
class RMatrix {
public:
	RMatrix();
	RMatrix(unsigned int row, unsigned int cols);
	RMatrix(unsigned int row, unsigned int cols, const T* const init_data);
	RMatrix(const vector<T>& values, const int* sparsity_pattern, int sparsity_num);
	RMatrix(const RMatrix&);
	
	unsigned int rows() const;
	unsigned int cols() const;
	
	T&       operator()(unsigned int ii, unsigned int jj);
	const T& operator()(unsigned int ii, unsigned int jj) const;
	
	
	T&       get(unsigned int ii, unsigned int jj);
	const T& get(unsigned int ii, unsigned int jj) const;
	const vector<T>& get_data() const { return data; }
	
	void     set(unsigned int ii, unsigned int jj, const T& val);
		
	const RMatrix& operator=(const RMatrix& rhs);
	RMatrix<T> operator*(const RMatrix<T>& rhs) const;
	RMatrix<T> operator-(const RMatrix<T>& rhs) const;
	RMatrix<T> mm(const RMatrix<T>& rhs) const { return *this * rhs; }

	void     set_col(unsigned int jj, const T* data);
	void     set_row(unsigned int ii, const T* data);
	
	Vector3D operator*(const Vector3D&) const;
	Vector3D mv(const Vector3D&) const; // for tcl	
	RMatrix<T> get_transpose() const;		
			
	void print() const;
	
	static RMatrix<T> get_rotation_matrix(const Vector3D& n, const double& cos_a);

	static void get_eigensystem(const RMatrix<T>& matrix, vector<T>& eigenvalues, vector<T>& eigenvectors);	
	static bool hermitian(const RMatrix<cplx>& matrix);
	 			
private:
	unsigned int nrows;
	unsigned int ncols;
	vector<T>    data;
};

template<class T>
bool RMatrix<T>::hermitian(const RMatrix<cplx>& matrix) {	
	TDKP_ASSERT(matrix.cols() == matrix.rows(), "only a square matrix can be hermitian");
	for(unsigned int ii = 0; ii < matrix.rows(); ii++) {		
		if(matrix(ii,ii).imag() != 0.0) {
			return false;	
		}
		for(unsigned int jj = ii + 1; jj < matrix.cols(); jj++) {
			if(tdkp_math::abs(matrix(ii,jj) - conj(matrix(jj,ii))) > 1.0e-12) {
				return false;
			} 
		}	
	}
	return true;
}

// solve x = A^{-1}d
void solve(const RMatrix<cplx>& A, const vector<cplx>& d, vector<cplx>& x);

extern const RMatrix<double> identity_matrix;
extern const Vector3D ex,ey,ez;

template<class T>
std::ostream& operator<<(std::ostream& out, const RMatrix<T>& rmat) {
	out << "(RMatrix " << rmat.rows() << " x " << rmat.cols() << ")\n";
	out.precision(6);
	for(unsigned int ii = 0; ii < rmat.rows(); ii++) {
		for(unsigned int jj = 0; jj < rmat.cols(); jj++) {
			out << setw(10) << rmat(ii,jj) << " ";	
		}
		out << "\n";
	} 
	return out;
}


template<class T>
RMatrix<T>::RMatrix()
: nrows(0),
  ncols(0),
  data(1)
{
	
}

template<class T>
RMatrix<T>::RMatrix(unsigned int row, unsigned int col)
: nrows(row),
  ncols(col),
  data(row * col, 0.0)
{	
}

template<class T> 
RMatrix<T>::RMatrix(unsigned int rows, unsigned int cols, const T* const init_data)
: nrows(rows),
  ncols(cols),
  data(rows*cols, 0.0) 
{
	for(unsigned int ii = 0; ii < rows; ii++) {
		for(unsigned int jj = 0; jj < cols; jj++) {
			(*this)(ii,jj) = init_data[ii * rows + jj];	
		} 
	}	
}

/** create matrix from sparse matrix data 
 * 
 * in order to save a little space, the kp blocks are also reduced
 * to their sparsity pattern. this constructor here is used to create
 * an RMatrix from such sparse data storage
 */ 
template<class T>
RMatrix<T>::RMatrix(const vector<T>& values, const int* sparsity_pattern, int sparsity_num) 
: nrows(0),
  ncols(0),
  data(0)
{
	TDKP_ASSERT(sparsity_num == (signed)values.size(), "sparse_num == values.size()");
	// find nrows and ncols
	for(int ii = 0; ii < sparsity_num; ii++) {
		nrows = max(int(nrows), int(sparsity_pattern[2*ii]));
		ncols =	max(int(ncols), int(sparsity_pattern[2*ii+1]));
	}		
	nrows++; ncols++; // sparsity is index, here i need #
	data.assign(nrows * ncols, 0.0);
	for(int ii = 0; ii < sparsity_num; ii++) {
		(*this)(sparsity_pattern[2*ii], sparsity_pattern[2*ii+1]) = values[ii];	
	}
}	

template<class T>
RMatrix<T>::RMatrix(const RMatrix<T>& copy)
: nrows(copy.nrows),
  ncols(copy.ncols),
  data(copy.data.begin(), copy.data.end()) 
{		
}

template<class T>
unsigned int RMatrix<T>::rows() const { return this->nrows; }

template<class T>
unsigned int RMatrix<T>::cols() const { return this->ncols; }

template<class T>
T& RMatrix<T>::operator()(unsigned int ii, unsigned int jj) {
	return const_cast<T&>(static_cast<const RMatrix<T>&>(*this)(ii,jj)); 	
}

template<class T>
const T& RMatrix<T>::operator()(unsigned int ii, unsigned int jj) const {
	//TDKP_ASSERT(ii < this->nrows && jj < this->ncols, "ii < nrows && jj < ncols");
	if(!(ii < this->nrows && jj < this->ncols)) {
		cout << "ii: " << ii << " nrows: " << nrows << " jj: " << jj << " ncols: " << ncols << "\n";	
	}
	return this->data[ii * this->ncols + jj];
}

template<class T> 
T& RMatrix<T>::get(unsigned int ii, unsigned int jj) {
	return (*this)(ii,jj);
}

template<class T>
const T& RMatrix<T>::get(unsigned int ii, unsigned int jj) const {
	return (*this)(ii,jj);
}		

template<class T> 
void RMatrix<T>::set(unsigned int ii, unsigned int jj, const T& val) {
	(*this)(ii,jj) = val;	
}

template<class T>
const RMatrix<T>& RMatrix<T>::operator=(const RMatrix& rhs) {
	this->nrows = rhs.nrows;
	this->ncols = rhs.ncols;
	this->data.resize(nrows * ncols);
	this->data.insert(this->data.begin(), rhs.data.begin(), rhs.data.end());
	return *this;	
}

template<class T>
void RMatrix<T>::set_col(unsigned int jj, const T* cdata) {
	for(unsigned int ii = 0; ii < this->nrows; ii++) {
		(*this)(ii,jj) = cdata[ii];		
	} 
}

template<class T>
void RMatrix<T>::set_row(unsigned int ii, const T* rdata) {
	for(unsigned int jj = 0; jj < this->ncols; jj++) {
		(*this)(ii,jj) = rdata[jj];		
	}
}

template<class T> 
Vector3D RMatrix<T>::operator*(const Vector3D& rhs) const {
	TDKP_ASSERT(this->nrows == 3 && this->ncols == 3, "matrix is not 3x3!");	
	Vector3D res;
	for(unsigned int ii = 0; ii < 3; ii++) {
		for(unsigned int jj = 0; jj < 3; jj++) {
			res(ii) += (*this)(ii,jj) * rhs(jj);
		}	
	}	
	return res;
}

/** tcl does not like operators */
template<class T> 
Vector3D RMatrix<T>::mv(const Vector3D& rhs) const {
	return (*this) * rhs;	
} 

/** create rotation matrix around normal n using cos(a) where a is the angle
 */
template<class T> 
RMatrix<T> RMatrix<T>::get_rotation_matrix(const Vector3D& n, const double& cos_a) {
		
	TDKP_BOUNDS_ASSERT(abs(n.norm() - 1) < 1.0e-9, "abs(n.norm() - 1) < 1.0e-9");
	RMatrix<T> rotation_matrix(3,3);	
	double sin_a = sqrt(1.0 - (cos_a * cos_a));
	rotation_matrix(0,0) = cos_a + n(0) * n(0) * (1.0 - cos_a);
	rotation_matrix(0,1) = n(0) * n(1) * (1.0 - cos_a) - n(2) * sin_a;
	rotation_matrix(0,2) = n(0) * n(2) * (1.0 - cos_a) + n(1) * sin_a;
	rotation_matrix(1,0) = n(1) * n(0) * (1.0 - cos_a) + n(2) * sin_a;
	rotation_matrix(1,1) = cos_a + n(1) * n(1) * (1.0 - cos_a);
	rotation_matrix(1,2) = n(1) * n(2) * (1.0 - cos_a) - n(0) * sin_a;
	rotation_matrix(2,0) = n(2) * n(0) * (1.0 - cos_a) - n(1) * sin_a;
	rotation_matrix(2,1) = n(2) * n(1) * (1.0 - cos_a) + n(0) * sin_a;
	rotation_matrix(2,2) = cos_a + n(2) * n(2) * (1.0 - cos_a);	
 
 	return rotation_matrix;
 	
}
 
/** get transpose of a matrix ... */
template<class T> 
RMatrix<T> RMatrix<T>::get_transpose() const {
	RMatrix<T> trans(this->ncols, this->nrows);	
	for(unsigned int ii = 0; ii < this->nrows; ii++) {
		for(unsigned int jj = 0; jj < this->ncols; jj++) {
			trans(jj,ii) = (*this)(ii,jj);				
		}			
	}
	return trans;
}
template<class T> 
RMatrix<T> RMatrix<T>::operator*(const RMatrix<T>& rhs) const {
	TDKP_ASSERT(this->ncols == rhs.nrows, "this->ncols == rhs.nrows");
	TDKP_ASSERT(rhs.ncols == this->nrows, "rhs.ncols == this->nrows");
	RMatrix<T> tmp(this->nrows, rhs.nrows);

#ifndef DEADRAT	
#pragma omp parallel for default(shared)
#endif
	for(int ii = 0; ii < (signed)tmp.nrows; ii++) {
		for(unsigned int jj = 0; jj < tmp.ncols; jj++) {
			for(unsigned int kk = 0; kk < this->ncols; kk++) {
				tmp(ii,jj) += (*this)(ii,kk) * rhs(kk,jj);
			} 	
		}	
	}
	
	return tmp;
		
}

template<class T> 
RMatrix<T> RMatrix<T>::operator-(const RMatrix<T>& rhs) const {
	TDKP_ASSERT(this->ncols == rhs.ncols, "this->ncols == rhs.ncols");
	TDKP_ASSERT(rhs.nrows == this->nrows, "rhs.nrows == this->nrows");
	RMatrix<T> tmp(this->nrows, this->ncols);	

	for(unsigned int ii = 0; ii < tmp.nrows; ii++) {
		for(unsigned int jj = 0; jj < tmp.ncols; jj++) {
			tmp(ii,jj) = (*this)(ii,jj) - rhs(ii,jj);
		}	
	}
	
	return tmp;
		
}

template<class T>
void RMatrix<T>::print() const {
	cout << *this << endl;	
}
 
} // end namespace

#endif /*VECTOR3D_H_*/
