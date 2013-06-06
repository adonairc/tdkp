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

#include <math.h>

#include "tdkp/common/Vector3D.h"
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"

using namespace std;

extern "C" {
	extern void F77_Name(zgesv)(int* n, int* nrhs, complex<double> *a, int* lda, int *ipiv, complex<double> *b, int* ldb, int *info);
}


namespace tdkp {

// defining for acml acml takes doublecomplex as a struct ..

extern "C" {
	
#ifdef NOACML
	// ---------------------------------------------------------
	// interface definition of fortran lapack functions 
	// ---------------------------------------------------------
	void zheev_(char* jobz, char* uplo, int* n, cplx* a, int* lda, double* w, 
	            cplx* work, int* lwork, double* rwork, int* info);
	void zheevd_(char* jobz, char* uplo, int* n, cplx* a, int* lda, double* w, 
	            cplx* work, int* lwork, double* rwork, int* lrwork, int* iwork, 
	            int* liwork, int* info);	             
#else
	// ---------------------------------------------------------	
	// interface definition for acml cstyle call
	// ---------------------------------------------------------
	void zheev(char jobz, char uplo, int n, cplx *a, int lda, double *w, int *info);
	void zheevd(char jobz, char uplo, int n, cplx *a, int lda, double *w, int *info);
#endif	
}

	
const double init_identity[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};  	
	
const RMatrix<double> identity_matrix(3,3,init_identity);
const Vector3D ex(1.0, 0.0, 0.0);
const Vector3D ey(0.0, 1.0, 0.0);
const Vector3D ez(0.0, 0.0, 1.0);
	
Vector3D operator*(const double a, const Vector3D &rhs) {
  return rhs * a;
}

Vector3D operator-(const Vector3D& rhs) {
  Vector3D ret(rhs);
  return -1.0 * ret;    	  	
}

Vector3D::Vector3D() {
  v[0] = v[1] = v[2] = 0.0;
}

Vector3D::Vector3D(const Vector3D &copy) {
  v[0] = copy.v[0];
  v[1] = copy.v[1];
  v[2] = copy.v[2];
}

Vector3D::Vector3D(const double* init) {
  v[0] = init[0]; v[1] = init[1]; v[2] = init[2];
}

Vector3D::Vector3D(double x, double y, double z) {
  v[0] = x; v[1] = y; v[2] = z;
}

/** init vector as a * x + (1 - x) * b
 */
Vector3D::Vector3D(const double* a, const double* b, double x) {
	for(int ii = 0; ii < 3; ii++) {
		v[ii] = a[ii] * x + (1.0 - x) * b[ii];		
	}
}

Vector3D& Vector3D::operator=(const Vector3D &rhs) {
	if(this == &rhs) {
		return (*this);	
	}
	v[0] = rhs.v[0];
	v[1] = rhs.v[1];
	v[2] = rhs.v[2];
	return (*this);			
}

Vector3D Vector3D::operator-(const Vector3D &rhs) const {
  Vector3D ret;
  ret.v[0] = v[0] - rhs.v[0];
  ret.v[1] = v[1] - rhs.v[1];
  ret.v[2] = v[2] - rhs.v[2];
  return ret;
}

Vector3D Vector3D::operator+(const Vector3D &rhs) const {
  Vector3D ret;
  ret.v[0] = v[0] + rhs.v[0];
  ret.v[1] = v[1] + rhs.v[1];
  ret.v[2] = v[2] + rhs.v[2];
  return ret;
}

Vector3D Vector3D::operator*(const double& rhs) const {
  Vector3D ret((*this));
  ret.v[0] *= rhs;
  ret.v[1] *= rhs;
  ret.v[2] *= rhs;
  return ret;
}

double Vector3D::norm() const {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void Vector3D::normalize() {
   double nr = this->norm();
   TDKP_ASSERT(nr > 0.0, "norm of vector that should be normalized is zero!");
   v[0] /= nr;
   v[1] /= nr;
   v[2] /= nr;
}

double Vector3D::dot_product(const Vector3D &x, const Vector3D &y) {
	return x.v[0] * y.v[0] + x.v[1] * y.v[1] + x.v[2] * y.v[2];	
}
Vector3D Vector3D::cross_product(const Vector3D &x, const Vector3D &y) {
	Vector3D res;
	res.v[0] = x.v[1] * y.v[2] - x.v[2] * y.v[1];
	res.v[1] = x.v[2] * y.v[0] - x.v[0] * y.v[2];
	res.v[2] = x.v[0] * y.v[1] - x.v[1] * y.v[0];		
	return res;
}

ostream& operator<<(ostream& out, const Vector3D& vec) {
  out << "<" << vec.v[0] << ", " << vec.v[1] << ", " << vec.v[2] << ">";
  return out;
}

bool Vector3D::operator==(const Vector3D& rhs) const {
	if(this->norm() > 0.0) {
		return ((*this - rhs).norm() / this->norm()) < 1.0e-12;
	} else if(rhs.norm() < 1.0e-12) {
		return true; 	
	} else {
		return false;	
	}
		
}
bool Vector3D::operator!=(const Vector3D& rhs) const {
	return !(*this == rhs);	
}

template<>
void RMatrix<cplx>::get_eigensystem(const RMatrix<cplx>& matrix, vector<cplx>& eigenvalues, vector<cplx>& eigenvectors) {

	TDKP_ASSERT(matrix.rows() == matrix.cols(), "matrix is not a square matrix!");	
	eigenvectors.resize(matrix.rows() * matrix.rows());
	bool hermitian = true;
	for(unsigned int ii = 0; ii < matrix.rows(); ii++) {
		for(unsigned int jj = 0; jj < matrix.rows(); jj++) {
			if(ii <= jj && tdkp_math::abs(conj(matrix(ii,jj)) - matrix(jj,ii)) > 1.0e-13) {				 
				hermitian = false;	
			}			
			eigenvectors[jj * matrix.rows() + ii] = matrix.data[ii * matrix.rows() + jj];
		}			
	}	
	eigenvalues.resize(matrix.rows());
	if(hermitian) {
		int info = 0;
		vector<double> eigenvalues_double(eigenvalues.size());
#ifdef NOACML
		char jobz = 'V';
		char uplo = 'U';
		int  N    = matrix.rows();
		int  lwork  = 2 * (2*N + N*N);
		int  lrwork = 2 * (1 + 5 * N + 2 * N * N);
		int  liwork = 6 + 10 * N;
		vector<cplx>   work(lwork, 0.0);  
		vector<double> rwork(lrwork, 0.0);
		vector<int>    iwork(liwork, 0);						
		zheevd_(
			&jobz, &uplo, &N, &eigenvectors[0], &N, &eigenvalues_double[0], 
			&work[0], &lwork, &rwork[0], &lrwork, &iwork[0], &liwork, &info
		);
#else			
		zheevd('V', 'U', matrix.rows(), &eigenvectors[0], matrix.rows(), &eigenvalues_double[0], &info);
#endif		
		if(info != 0) {
			Logger::get_instance()->emit(LOG_WARN, "zheevd failed determining eigensystems. falling back to zheev");
			// try using other zheev routine
			for(unsigned int ii = 0; ii < matrix.rows(); ii++) {
				for(unsigned int jj = 0; jj < matrix.rows(); jj++) {
					eigenvectors[jj * matrix.rows() + ii] = matrix.data[ii * matrix.rows() + jj];
				}			
			}		
			info = 0;	
#ifdef NOACML
			char jobz = 'V';
			char uplo = 'U';
			int  N    = matrix.rows();
			int  lwork  = 2 * (2*N + N*N);
			int  lrwork = 2 * (3*N - 2);
			vector<cplx>   work(lwork, 0.0);  
			vector<double> rwork(lrwork, 0.0);
			zheev_(
				&jobz, &uplo, &N, &eigenvectors[0], &N, &eigenvalues_double[0],
				&work[0], &lwork, &rwork[0], &info
			);
			TDKP_GENERAL_EXCEPTION("sorry, but lapack full solvers uses ACML ...");
#else				
			zheev('V', 'U',matrix.rows(), &eigenvectors[0], matrix.rows(), &eigenvalues_double[0], &info);
#endif			
			if(info != 0) {							 		
				TDKP_GENERAL_EXCEPTION("zheev returned nonzero info: " << info);
			}
		}
		sort(eigenvalues_double.begin(), eigenvalues_double.end());
		for(unsigned int ii = 0; ii < eigenvalues_double.size(); ii++) {
			eigenvalues[ii] = eigenvalues_double[ii];		 
		}
#ifdef DEBUG		
		// rayleigh quotient must match eigenvalue!
		vector<cplx> a(matrix.rows());
		for(unsigned int ee = 0; ee < eigenvalues.size(); ee++) {
			cplx rq = 0.0;
			for(unsigned int ii = 0; ii < matrix.rows(); ii++) {
				for(unsigned int jj = 0; jj < matrix.rows(); jj++) {
					rq += conj(eigenvectors[ee * matrix.rows() + ii]) 	
					    * matrix(ii,jj)
					    * eigenvectors[ee * matrix.rows() + jj];
				}
			}
			TDKP_ASSERT(tdkp_math::abs(eigenvalues[ee] - rq) < 1.0e-10, "tdkp_math::abs(eigenvalues[ee] - rq) < 1.0e-10");				
		}
#endif		
	} else {
		TDKP_GENERAL_EXCEPTION("eigenvalues for non-hermitian matrices are not supported yet");	
	}			
	
}

template<>
void RMatrix<double>::get_eigensystem(const RMatrix<double>& matrix, vector<double>& eigenvalues, vector<double>& eigenvectors) {
	TDKP_GENERAL_EXCEPTION("implement me");	
}

// solve x = A^{-1}d 
void solve(const RMatrix<cplx>& A, const vector<cplx>& d, vector<cplx>& x) {
	
	unsigned int size = A.rows();
	TDKP_ASSERT(A.rows() == A.cols(), "");
	TDKP_ASSERT(d.size() == size, "");
	x.resize(size);
			
	complex<double>* copy_matrix = new complex<double>[size * size];	
	for(unsigned int ii = 0; ii < size; ii++) {
		for(unsigned int jj = 0; jj < size; jj++) {
			copy_matrix[jj * size + ii] = A(ii,jj);
		}
		x[ii] = d[ii];	
	}

	int n    = size;
	int nrhs = 1;
	int lda  = n;
	vector<int> ipiv(n);
	int ldb  = n;
	int info;
	
	F77_Name(zgesv)(&n, &nrhs, copy_matrix, &lda, &ipiv[0], &x[0], &ldb, &info);
	
	delete[] copy_matrix;
	if(info != 0) {
		TDKP_GENERAL_EXCEPTION("zgesv_ returned " << info);	
	}
	 		
}

} // end namespace
