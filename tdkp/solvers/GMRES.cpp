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

#include "tdkp/common/all.h"
#include "tdkp/solvers/GMRES.h"
#include "tdkp/common/Configuration.h"

#include <math.h> 

namespace tdkp {
// from netlib.org/templates
//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************

cplx dot(const GMRES::GMRESVector& x1, const GMRES::GMRESVector& x2) {
	cplx ret = 0.0;
	for(int ii = 0; ii < x1.size(); ii++) {
		ret += conj(x1(ii)) * x2(ii); 	
	}
	return ret;	
}

double norm(const GMRES::GMRESVector& x) {
	double ret = 0.0;
	for(int ii = 0; ii < x.size(); ii++) {
		ret += (x(ii).real() * x(ii).real() + x(ii).imag() * x(ii).imag()); 	
	}	
	return sqrt(ret);
}

void mat_vec(const RMatrix<cplx>& A, const GMRES::GMRESVector& x, GMRES::GMRESVector& w) {
	const vector<cplx>& data = A.get_data();
	for(int ii = 0; ii < x.size(); ii++) {
		w(ii) = 0.0;
		for(int jj = 0; jj < x.size(); jj++) {
			w(ii) += data[ii * x.size() + jj] * x(jj);				
		}		
	}	
}


template < class Matrix, class Vector >
void 
Update(Vector &x, int k, Matrix &h, const Vector &s, Vector v[])
{
  Vector y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y(i) /= h(i,i);
    for (int j = i - 1; j >= 0; j--) {      
      y(j) -= h(j,i) * y(i);
    }
  }

  for (int j = 0; j <= k; j++) { 
  	for(int n = 0; n < x.size(); n++) { 
    	x(n) += v[j](n) * y(j);
  	}
  } 
}


template < class Real >
Real 
abs(Real x)
{
  return (x > 0 ? x : -x);
}

template<class Real> 
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (abs(dy) > abs(dx)) {
    Real temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    Real temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}
template<class Real> 
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  Real temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

template < class Operator, class Vector, 
           class Matrix, class Real >
int 
NETLIB_GMRES(const Operator &A, Vector &x, const Vector &b,
      Matrix &H, int &m, int &max_iter,
      Real &tol)
{
  Real resid;
  int i, j = 1, k;
  Vector s(m+1), cs(m+1), sn(m+1), w(A.cols()), wtmp(A.cols());
  
  Real normb = norm(b);
  // r = b - Ax;  
  Vector r(x.size());
  mat_vec(A,x,r);
  r *= -1.0;
  r += b;
  Real beta  = norm(r);
  
  if (normb == 0.0) {
    normb = 1;
  }
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  Vector *v = new Vector[m+1];

  while (j <= max_iter) {
    v[0] = r;
    v[0] *= (1.0 / beta);    // ??? r / beta
    s.assign(0.0);
    s(0) = beta;
    
    for (i = 0; i < m && j <= max_iter; i++, j++) {
      // w = A * v[i];
      mat_vec(A,v[i],w);      
      for (k = 0; k <= i; k++) {
        H(k, i) = dot(w, v[k]);
        // w -= H(k, i) * v[k];
        wtmp = v[k];
        wtmp *= (-1.0) * H(k, i);
        w += wtmp;        
      }
      H(i+1, i) = norm(w);
      // w * (1.0 / H(i+1, i))      
      v[i+1] = w;
      v[i+1] *= (1.0 / H(i+1, i)); // ??? w / H(i+1, i)

      for (k = 0; k < i; k++) {
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
      }
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
      
      if ((resid = abs(s(i+1)) / normb) < tol) {
        Update(x, i, H, s, v);
        tol = resid;
        max_iter = j;
        delete[] v;
        return 0;
      }
    }
    Update(x, m - 1, H, s, v);
    // r = b - A * x;
    mat_vec(A,x,r);
  	r *= -1.0;
  	r += b;
  	
    beta = norm(r);
    if ((resid = beta / normb) < tol) {
      tol = resid;
      max_iter = j;
      delete [] v;
      return 0;
    }
  }
  
  tol = resid;
  delete [] v;
  return 1;
}

GMRES::GMRESVector::GMRESVector(unsigned int n_)
: n(n_), 
  my_array(true) 
{
	val = new cplx[n];
	memset(val, 0, sizeof(cplx) * n);
}

GMRES::GMRESVector::GMRESVector(vector<cplx>& target)
: n(target.size()),
  my_array(false)
{	
	val = &target[0];
}

GMRES::GMRESVector::GMRESVector(const vector<cplx>& source) 
: n(source.size()),
  my_array(true)
{	
	val = new cplx[n];
	memcpy(val, &source[0], sizeof(cplx) * source.size());		
}

GMRES::GMRESVector::GMRESVector() 
: n(0),
  my_array(false)
{
	
}

GMRES::GMRESVector::GMRESVector(const GMRESVector& copy)
: n(copy.n),
  my_array(true)
{
	val = new cplx[n];
	memcpy(val, copy.val, sizeof(cplx) * n);	
}

GMRES::GMRESVector::~GMRESVector() {
	if(my_array) {
		delete[] val; val = 0;	
	}	
}


void GMRES::GMRESVector::operator*=(const double& rhs) {
	for(int ii = 0; ii < n; ii++) {
		val[ii] *= rhs;	
	}	
}
void GMRES::GMRESVector::operator*=(const cplx& rhs) {
	for(int ii = 0; ii < n; ii++) {
		val[ii] *= rhs;	
	}	
}

void GMRES::GMRESVector::operator+=(const GMRESVector& rhs) {
	for(int ii = 0; ii < n; ii++) {
		val[ii] += rhs.val[ii];	
	}	
}

void GMRES::GMRESVector::assign(const double& x) {
	for(int ii = 0; ii < n; ii++) {
		val[ii] = x;	
	}	
}

void GMRES::GMRESVector::operator=(const GMRESVector& rhs) {
	if(n == 0) {
		n = rhs.n;
		val = new cplx[n];
		my_array = true;	
	}
	if(n != rhs.n) {
		TDKP_GENERAL_EXCEPTION("n != rhs.n");	
	}
	memcpy(val, rhs.val, sizeof(cplx) * n);	
}

void GMRES::solve(int nrestart, const RMatrix<cplx>& A, const vector<cplx>& d, vector<cplx>& x) {
	
	RMatrix<cplx> Householder(nrestart + 1, nrestart + 1);	
	GMRESVector Gx(x);
	GMRESVector Gb(d);
	double tol  = Configuration::get_instance()->get("clc_gmres_tolerance");
	int   maxit = 200;
	int ret = NETLIB_GMRES<RMatrix<cplx>, GMRESVector, RMatrix<cplx>, double>(
		A, Gx, Gb, Householder, nrestart, maxit, tol
	);	
	TDKP_ASSERT(ret == 0, "GMRES DID NOT CONVERGE");
	
}

} // end of namespace
