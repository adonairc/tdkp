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


#include <iostream>
#include <math.h>
#include <assert.h>

#include "tdkp/common/ICurves.h"


namespace tdkp {


int ICurve::get_x_length() const {
	// extremely inefficient ;-)
	GridPoints tmp;
	this->GetIntrinsicGrid(tmp);
	return tmp.size(); 	
}
int ICurve::get_num_y_sets() const {
	return 1;
} 
void ICurve::get_x(int xidx, vector<double> &x) const {
	TDKP_ASSERT(xidx == 0, "xidx == 0");
	this->GetIntrinsicGrid(x);		
}
void ICurve::get_y(int yidx, vector<double>& y) const {
	TDKP_ASSERT(yidx == 0, "yidx == 0");
	// extremely unefficient ;-)
	GridPoints tmp;
	this->GetIntrinsicGrid(tmp);
	this->EvalAt(tmp, y);	
} 
string ICurve::get_x_identifier(int xidx) const {
	return string("x");	
}
string ICurve::get_y_identifier(int yidx) const {
	return string("y");			
}	

// analytic cubic root finder
typedef vector<double> Roots;
void FindCubicRoots(double alpha, double beta, double gamma, Roots& x);
void FindQuadraticRoots(double alpha, double beta, Roots& x);

/** Free unused resources.
 */
void LinearCurve::Cleanup()
{
  grid.clear();
  alpha.clear();
  beta.clear();
  len = 0;
}

/** Create curve from another curve.
 *
 * Using intrinsic grid of ICurve object as our grid.
 *
 * @param curve  Source curve.
 */
void LinearCurve::CopyFrom(const ICurve& curve)
{
  Cleanup();

  curve.GetIntrinsicGrid(grid);
  len = grid.size();
  curve.EvalAt(grid, beta);
  CalcDerivative();
}

/** Create curve from grid and value vector.
 *
 * @param x  Grid to be used.
 * @param y  Curve values at grid points.
 */
void LinearCurve::CopyFrom(const GridPoints& x, const GridValues& y)
{
  assert(x.size() == y.size());
  assert(x.size() >= 2);

  Cleanup();

  grid = x;
  len = grid.size();
  beta = y;
  CalcDerivative();
}

/** Evaluate curve at position x using interpolation if necessary.
 *
 * @param x  Point where curve should be evaluated.
 */
double LinearCurve::EvalAt(const double& x) const
{
  assert(len >= 2);
  assert(alpha.size() == len-1);

  // find segment in which x lies
  int p;
  for (p = 0; p < (signed)len-1; p++) { if (x < grid[p+1]) break; }
  if (p == (signed)len-1) p--;
  // interpolate
  return beta[p] + (x - grid[p]) * alpha[p];
}

/** Evaluate curve at positions x using interpolation if necessary
 *
 * @param x  Points where curve should be evaluated.
 * @param y  Will store values at points given by x.
 */
void LinearCurve::EvalAt(const GridPoints& x, GridValues& y) const
{
  	y.resize(x.size(), 0.0);
  	// we couldn't do it more efficiently than this...
#pragma omp parallel for default(shared)  
  	for (int i = 0; i < (signed)x.size(); i++) {
	    y[i] = EvalAt(x[i]);
  	}
}

/** find inverse */
void LinearCurve::FindInverse(const double& y, GridPoints& x) const {
	x.resize(0);
	for(unsigned int ii = 0; ii < len - 1; ii++) {				
		if(beta[ii] <= y && beta[ii + 1] > y) {
			if(beta[ii] != beta[ii + 1]) {
				double t = (y - beta[ii]) / alpha[ii];
				TDKP_BOUNDS_ASSERT(t <= grid[ii+1] - grid[ii], "t <= grid[ii+1] - grid[ii]");
				x.push_back(grid[ii] + t);
			} else {
				x.push_back(grid[ii]);	
			}	
		} 
	}
	if(beta.back() == y) {
		x.push_back(grid.back());	
	}
#ifdef DEBUG	
	// --------------------------------------
	// test the inverse
	// --------------------------------------
	const double threshold = 1.0e-6;
	for(unsigned int ii = 0; ii < x.size(); ii++) {
		double yi = this->EvalAt(x[ii]);
		if(tdkp_math::abs(yi - y) > threshold) {
			double ratio = y != 0 ? y:1.0;
			ostringstream sout;
			sout << " bad linear inverse calculation: " << ii << " xi: " << x[ii] << " yi: " << yi << " y0_value: " << y;
			sout << " (" << tdkp_math::abs(yi - y) / ratio << ") ";
			Logger::get_instance()->emit(LOG_WARN, sout.str());	
		}	
	}	
#endif
}

void SplineCurve::FindInverse(const double& y0_value, GridPoints& xi_values) const {

	// ---------------------------------------
	// loop over all segments
	// ---------------------------------------
	vector<double> segments_xis;
	xi_values.clear();
	for(unsigned int ii = 0; ii < len - 1; ii++) {
		// check if y0 is in the segment's range
		if(this->segment_ranges[ii].min_value <= y0_value && this->segment_ranges[ii].max_value >= y0_value) {
			this->FindInverse(ii, y0_value, segments_xis);
			while(segments_xis.size() > 0) {
				xi_values.push_back(segments_xis.back());
				segments_xis.pop_back();	
			}
		}  	
	}	
#ifdef DEBUG	
	// --------------------------------------
	// test the inverse
	// --------------------------------------
	const double threshold = 1.0e-6;

	for(unsigned int ii = 0; ii < xi_values.size(); ii++) {
		double yi = this->EvalAt(xi_values[ii]);
		if(tdkp_math::abs(yi - y0_value) > threshold) {
			double ratio = y0_value;
			if(y0_value == 0.0) {
				ratio = 1.0;	
			}			
			ostringstream sout;
			sout << " bad spline interverse calculation: " << ii << " xi: " << xi_values[ii] << " yi: " << yi << " y0_value: " << y0_value;
			sout << " (" << tdkp_math::abs(yi - y0_value) / ratio << ") ";
			Logger::get_instance()->emit(LOG_WARN, sout.str());	
		}	
	}	
#endif	
}

/** prepare min and max of all segments (needed for inverse calculations) */
void SplineCurve::prepare_segments() {
	this->segment_ranges.resize(len - 1);	
	for(int ii = 0; ii < (signed)len - 1; ii++) {		
		this->GetMinMax(ii, this->segment_ranges[ii].min_value, this->segment_ranges[ii].max_value);
	}
}

double LinearCurve::EvalDerivativeAt(const double& x) const {
	
	if(x < grid.front()) {
		return alpha.front();	
	} else if(x > grid.back()) {
		return alpha.back();	
	} else {
		TDKP_BOUNDS_ASSERT(len >= 2, "len >= 2");
		TDKP_BOUNDS_ASSERT(alpha.size() == len-1, "alpha.size() == len-1");

		// find segment in which x lies
		int p;
		for (p = 0; p < (signed)len-1; p++) { if (x < grid[p+1]) break; }
		if (p == (signed)len-1) p--; 
		TDKP_BOUNDS_ASSERT( p < (signed)alpha.size(), " p < (signed)alpha.size()");
		if(p == 0) {
			return alpha[p];
		} else {	 
			return (alpha[p] + alpha[p - 1]) / 2.0;
		}		 		
	}
  
}

void LinearCurve::EvalDerivativeAt(const GridPoints& x, GridValues& y) const
{
  	y.resize(x.size(), 0.0);
  	// we couldn't do it more efficiently than this...
#pragma omp parallel for default(shared)  
  	for (int i = 0; i < (signed)x.size(); i++) {
	    y[i] = EvalDerivativeAt(x[i]);
  	}
}

/** Integrate over this curve.
 */
double LinearCurve::Integrate() const
{
  return Integrate(grid.front(), grid.back());
}

void LinearCurve::GetIntrinsicGridRange(double& xl, double& xr) const {
	xl = grid.front();
	xr = grid.back();	
}



/** Integrate over part of this curve.
 *
 * @param start  Starting point.
 * @param end  End point.
 */
double LinearCurve::Integrate(double start, double end) const
{
  assert(start >= grid.front());
  assert(end <= grid.back());
  assert(end >= start);

  double res = 0.0;
  for (unsigned int i = 0; i < len-1; i++) {
    if ((end > grid[i]) && (start < grid[i+1]))
      res += IntegrateSegment(i, start, end);
  }
  return res;
}

/** Integrate over part of this curve lying in segment seg.
 *
 * @param seg  Segment to be integrated over.
 * @param start  Starting point.
 * @param end  End point.
 */
double LinearCurve::IntegrateSegment(int seg, double start, double end) const
{
  // since curve is piecewise linear, this can be done analytically for each segment
  // f(x) = a (x - x0) + b   =>   int(f, x1, x2) = a/2 [x2^2 - x1^2 - 2 x0 (x2 - x1)] + b (x2 - x1)
  //                                             = dx (a/2 [x2 + x1 - 2 x0] + b)
  const double x1 = (grid[seg] > start) ? grid[seg] : start;
  const double x2 = (grid[seg+1] < end) ? grid[seg+1] : end;
  const double dx = x2 - x1;
  assert(dx > 0);

  return dx * (alpha[seg]/2.0*(x2 + x1 - 2*grid[seg]) + beta[seg]);
}





/** Calculate derivative used for linear interpolation in EvalAt method.
 */
void LinearCurve::CalcDerivative()
{
  assert(len >= 2);
  assert(alpha.size() == 0);

  alpha.resize(len-1, 0.0);
 
  for (int i = 0; i < (signed)len-1; i++) {
    assert(grid[i+1] > grid[i]);
    alpha[i] = (beta[i+1] - beta[i]) / (grid[i+1] - grid[i]);
  }
}




/** Stream a LinearCurve object.
 *
 * @param out  Stream object.
 * @param curve  Curve to be streamed.
 */
ostream& operator<<(ostream& out, const LinearCurve& curve) {
  out << "Curve [LinearCurve]: " << curve.len << " points." << endl;
  out << "y(x) = a (x - x0) + b     (x0 < x < x1)" << endl;
  out << endl;
  out << "  interval                       a              b" << endl;
  out << " --------------------------------------------------------" << endl;
  for (unsigned int i = 0; i < curve.len-1; i++) {
    char string[60];
    sprintf(string, "  %8.2E < x < %8.2E   %15.4E%15.4E", curve.grid[i], curve.grid[i+1],
        curve.alpha[i], curve.beta[i]);
    out << string << endl;
  }
  out << endl;
  return out;
}

/** element wise multiplication with a double */
LinearCurve LinearCurve::operator*(const double& mul) const {
	LinearCurve curve(*this);
	TDKP_ASSERT(alpha.size() == beta.size() - 1, "alpha.size() == beta.size()");
	
#pragma omp parallel for default(shared)	
	for(int ii = 0; ii < (signed)alpha.size(); ii++) {
		curve.alpha[ii] *= mul;
		curve.beta[ii]   *= mul;
	}
	curve.beta.back() *= mul;
	return curve;
}
/** element wise division with a double */
LinearCurve LinearCurve::operator/(const double& mul) const {
	return (*this) * (1.0 / mul);	
}
/** element by element multiplication of two curves 
 * 
 * the x ranges of both curves must be equal! 
 */
LinearCurve LinearCurve::operator*(const ICurve& rhs) const {
	
	TDKP_ASSERT(rhs.GetGridReference().front() == this->GetGridReference().front(), "rhs.GetGridReference().front() == this->GetGridReference().front()");
	TDKP_ASSERT(rhs.GetGridReference().back() == this->GetGridReference().back(), "rhs.GetGridReference().back() == this->GetGridReference().back()");
	
	vector<double> lhs_val;
	vector<double> rhs_val;
	
	rhs.EvalAt(this->GetGridReference(), rhs_val);
	this->EvalAt(this->GetGridReference(), lhs_val);
		
	vector<double> res_val(lhs_val.size());

#pragma omp parallel for default(shared)		
	for(int ii = 0; ii < (signed)lhs_val.size(); ii++) {
		res_val[ii] = lhs_val[ii] * rhs_val[ii];		
	}	
	
	return LinearCurve(this->GetGridReference(), res_val);
	
}

LinearCurve LinearCurve::sqrt() const {
	
	vector<double> vals(beta.size());

#pragma omp parallel for default(shared)
	for(int ii = 0; ii < (signed)vals.size(); ii++) {
		vals[ii] = std::sqrt(beta[ii]);	
	}	
	return LinearCurve(grid, vals);
	
		
}

// ---------------------------------------------------
// spline curves
// ---------------------------------------------------

/** Free unused resources.
 */
void SplineCurve::Cleanup()
{
  grid.clear();
  coefs.clear();
  len = 0;
}

/** Create curve from another curve.
 *
 * Using intrinsic grid of ICurve object as our grid.
 *
 * @param curve  Source curve.
 */
void SplineCurve::CopyFrom(const ICurve& curve)
{
  Cleanup();

  curve.GetIntrinsicGrid(grid);
  len = grid.size();
  GridValues val;
  curve.EvalAt(grid, val);
  CalcSplines(val);
  prepare_segments();
}


/** Create curve from grid and value vector.
 *
 * @param x  Grid to be used.
 * @param y  Curve values at grid points.
 */
void SplineCurve::CopyFrom(const GridPoints& x, const GridValues& y)
{
  assert(x.size() == y.size());
  assert(x.size() >= 2);

  Cleanup();

  grid = x;
  len = grid.size();
  CalcSplines(y);
  prepare_segments();
}

/** Evaluate curve at position x using interpolation if necessary.
 *
 * @param x  Point where curve should be evaluated.
 */
double SplineCurve::EvalAt(const double& x) const
{
  assert(len >= 2);
  // find grid interval in which x lies
  int p;
  for (p = 0; p < (signed)len-1; p++) { if (x < grid[p+1]) break; }
  if (p == (signed)len-1) p--;
  // interpolate
  const double dx = x - grid[p];
  return coefs[p].d + coefs[p].c*dx + coefs[p].b*dx*dx + coefs[p].a*dx*dx*dx;
}

/** Evaluate curve at positions x using interpolation if necessary
 *
 * @param x  Points where curve should be evaluated.
 * @param y  Will store values at points given by x.
 */
void SplineCurve::EvalAt(const GridPoints& x, GridValues& y) const {
	y.resize(x.size(), 0.0);
  	const int length = x.size();
#pragma omp parallel for default(shared)
  	for (int i = 0; i < length; i++) {
    	y[i] = EvalAt(x[i]);
  	}
}


/** Evaluate first derivative of curve at position x using interpolation if necessary.
 *
 * @param x  Point where first derivative should be evaluated.
 */
double SplineCurve::EvalDerivativeAt(const double& x) const
{
  assert(len >= 2);
  // find grid interval in which x lies
  int p;
  for (p = 0; p < (signed)len-1; p++) { if (x < grid[p+1]) break; }
  if (p == (signed)len-1) p--;
  // interpolate
  const double dx = x - grid[p];
  return coefs[p].c + 2*coefs[p].b*dx + 3*coefs[p].a*dx*dx;
}

/** Evaluate first derivative of curve at positions x using interpolation if necessary
 *
 * @param x  Points where first derivative should be evaluated.
 * @param y  Will store values at points given by x.
 */
void SplineCurve::EvalDerivativeAt(const GridPoints& x, GridValues& y) const {
	y.resize(x.size(), 0.0);
	const int length = x.size();
#pragma omp parallel for default(shared)
  	for (int i = 0; i < length; i++) {
    	y[i] = EvalDerivativeAt(x[i]);
  	}
}


double SplineCurve::Integrate() const
{
  assert(false);  // not implemented yet
  return 0.0;
}

double SplineCurve::Integrate(double start, double end) const
{
  assert(false);  // not implemented yet
  return 0.0;
}



/** Calculate spline coefficents for a given grid and curve values.
 *
 * @param val  Values of curve at grid points stored in grid.
 */
void SplineCurve::CalcSplines(const GridValues& val)
{
  assert(len >= 2);
  assert(coefs.size() == 0);
  assert((bc == bc_notaknot) || (bc == bc_natural));

  const int N = len - 1;
  // calculate distance and first derivative between grid points
  // h_i = grid_{i+1} - grid_i
  // der_i = (val_{i+1} - val{i}) / h_i
  GridPoints h = GridPoints(N);
  GridValues der = GridValues(N);
  for (int i = 0; i < N; i++) {
    h[i] = grid[i+1] - grid[i];
    assert(h[i] > 0);
    der[i] = (val[i+1] - val[i]) / h[i];
  }

  // N the number of coefficents / spline segments (N = len - 1)
  // A spline polynomial s_i(x), 0 <= i < N is valid in the interval
  //               grid_i < x < grid_{i+1}
  // and given by
  //               s_i(x) = a_i*(x - grid_i)^3 + b_i*(x - grid_i)^2 + c_i*(x - grid_i) + d_i
  // such that
  //         (1,2) s_i(grid_i) = val_i   and   s_i(grid_{i+1}) = val_{i+1}         (2*N conditions)
  //           (3) s_{i-1}'(grid_i) = s_i'(grid_i)                                 (N-1 conditions)
  //           (4) s_{i-1}''(grid_i) == s_i''(grid_i)				   (N-1 conditions)
  //
  // The remaining two conditions are given by the boundary conditions, see below.
  //
  // By introducing a new variable sigma, we can reduce the problem to the
  // solution of a trigonal matrix.
  // We define sigma_i = s''(grid_i) and can now formulate (1-4) in terms of sigma to get
  // a recursion formula for sigma:
  //         sigma_{i-1} * h_{i-1} + sigma_{i+1} * h_i + sigma_i * 2 * (h_i + h_{i-1}) = g_i
  // with
  //         g_i = 6 * (der_i - der_{i-1})
  //
  // The coefficients can then be expressed in terms of the sigma_i as follows:
  //         a_i = (sigma_{i+1} - sigma_i) / (6 * h_i)
  //         b_i = sigma_i / 2
  //         c_i = der_i - (h_i / 6) * (2 * sigma_i + sigma_{i+1})
  //         d_i = val_i

  TrigMatrix A = TrigMatrix(len);
  TrigVector b = TrigVector(len);
  TrigVector sigma;

  // (a) define recursive part
  for (int i = 1; i < N; i++) {
    b[i] = 6.0 * (der[i] - der[i-1]);
    A[i].low = h[i-1];
    A[i].up = h[i];
    A[i].diag = 2.0 * (h[i] + h[i-1]);
  }

  // (b) define boundary conditions
  switch (bc) {
  // not-a-knot boundary conditions
  // s_0'''(x_1) = s_1'''(x_1),  s_{N-1}'''(x_N) = s_N'''(x_N)
  case bc_notaknot:
    b[0] = b[1];
    A[0].low = 0.0;
    A[0].diag = h[0] - (h[1] * h[1]) / h[0];
    A[0].up = 2.0 * (h[1] + h[0]) + ((h[1] + h[0]) * h[1]) / h[0];
    b[N] = b[N-1];
    A[N].low = 2.0 * (h[N-1] + h[N-2]) + ((h[N-1] + h[N-2]) * h[N-2]) / h[N-1];
    A[N].diag = h[N-1] - (h[N-2] * h[N-2]) / h[N-1];
    A[N].up = 0.0;
    break;
  // natural boundary conditions
  // s_0''(x_0) = 0,  s_N''(x_{N+1}) = 0
  case bc_natural:
    b[0] = b[N] = 0;
    A[0].low = A[0].up = 0; A[0].diag = 1;
    A[N].low = A[N].up = 0; A[N].diag = 1;
    break;
  }

  // (c) solve
  SolveTrigonal(A, b, sigma);
//  for (int i = 0; i < len; i++) cout << A[i].low << '\t' << A[i].diag << '\t' << A[i].up << '\t' << b[i] << endl;
//  cout << endl;

  // (d) extract coefficients
  coefs.resize(N);
  for (int i = 0; i < N; i++) {
    coefs[i].a = (sigma[i+1] - sigma[i]) / (6.0 * h[i]);
    coefs[i].b = sigma[i] / 2.0;
    coefs[i].c = der[i] - (h[i] / 6.0) * (2 * sigma[i] + sigma[i+1]);
    coefs[i].d = val[i];
  }
}


/** Solve simple linear equation Ax = b with A a triagonal matrix.
 *
 * The trigonal Matrix A is given by
 *
 *               | A[0].d    A[0].up                      0         |
 *               | A[1].low  A[1].d    A[1].up                      |
 *           A = |           A[2].low     .      .                  |
 *               |                        .      .        A[N-2].up |
 *               |        0                   A[N-1].low  A[N-1].d  |
 *
 * @param A  Matrix of size N x 3, storing diagonal and
 *           first secondary diagonals of A.
 * @param b  Vector of size N.
 * @param x  Solution vector of size N set by this function.
 */
void SplineCurve::SolveTrigonal(TrigMatrix& A, TrigVector& b, TrigVector& x)
{
  assert(A.size() == b.size());
  assert((A[0].low == 0.0) && (A[A.size()-1].up == 0.0));

  // A matrix and b vector will be changed by Gauss reduction.
  // NOTE: The upper secondary diagonal of A (-> A[i].up) does
  // not change during this transformation.
  int N = A.size();
  x.resize(N);

  if (A[0].diag != 0.0) {
    // Gaussian reduction of trigonal matrix in A to get upper triangular matrix
    for (int i = 1; i < N; i++) {
      // matrix transformation (i)' = (i) - C * (i-1)
      // where C is chosen such that the lower secondary diagonal element vanishes
      A[i].diag -= A[i].low/A[i-1].diag * A[i-1].up;
      b[i] -= A[i].low/A[i-1].diag * b[i-1];
      //A[i].low = 0; // could be omitted, as it is never used
    }
    // With A in upper triangular form
    // it is now trivial to find x from bottom to top
    x[N-1] = b[N-1]/A[N-1].diag;
    for (int i = N - 2; i >= 0; i--) {
      x[i] = (b[i] - A[i].up*x[i+1]) / A[i].diag;
    }
  } else {
    // special treatment if diagonal is zero
    // keep second line till the end, it will be used to solve for x[0]...
    double low1 = A[1].low; double diag1 = A[1].diag; double up1 = A[1].up; double b1 = b[1];
    // exchange second line with first one
    A[1].diag = A[0].up; A[1].low = A[0].diag; A[1].up = 0.0;
    // starting from our new second line
    for (int i = 2; i < N; i++) {
      A[i].diag -= A[i].low/A[i-1].diag * A[i-1].up;
      b[i] -= A[i].low/A[i-1].diag * b[i-1];
    }
    x[N-1] = b[N-1]/A[N-1].diag;
    for (int i = N - 2; i >= 1; i--) {
      x[i] = (b[i] - A[i].up*x[i+1]) / A[i].diag;
    }
    // now at last, we can determine x[0]...
    x[0] = (b1 - low1 * x[1] - up1 * x[2]) / diag1;
  }
}


/** Find curve minimum and maximum inside an interval
 *
 * @param segment  Segment number between 0 and len-2.
 * @param min  To store minimum.
 * @param max  To store maximum.
 */
double SplineCurve::GetMinMax(int segment, double& min, double& max) const
{
				
  assert((segment >= 0) && (segment < (signed)len-1));

  const SplineSegmentCoefs& c = coefs[segment];
  const double dx = grid[segment+1] - grid[segment];
  double pos = 0.0;

  // find extrema in interval (x0, x1) of function given by coefs
  //   f   = a (x - x0)^3 + b (x - x0)^2 + c (x - x0) + d
  // df/dx = 3a (x - x0)^2 + 2b (x - x0) + c = 0
  //
  // normalize such that y = (x - x0)/dx in [0, 1] and solve
  // f(y) = y^2 + 2b/(3a dx) y + c/(3a dx^2) = 0
  //
  Roots roots;
  FindQuadraticRoots((2*c.b) / (3*c.a*dx), c.c/(3*c.a*dx*dx), roots);
  // boundaries are possible extrema
  roots.push_back(0.0);
  roots.push_back(1.0); 
  bool minmaxset = false;

  for (unsigned int i = 0; i < roots.size(); i++) {
    // check for roots inside this interval
    if ((roots[i] >= 0.0) && (roots[i] <= 1.0)) {
      double xr = roots[i]*dx;
      double yr = c.a * xr*xr*xr + c.b * xr*xr + c.c * xr + c.d;
      if(minmaxset) {
        if (yr < min) { min = yr; pos = xr + grid[segment]; }
        if (yr > max) max = yr;
      } else {
        min = max = yr;        
        pos = xr + grid[segment];
        minmaxset = true;
      }
    }
  }
  assert(minmaxset);
  assert(min <= max);
  return pos;
}


/** Find inverse of a curve in a given interval
 *
 * Find all values x for which holds f(x) = y, where f(x) is the spline function
 * in the given interval.
 *
 * @param segment  Segment number between 0 and len-2.
 * @param y  Value y, sucht that f(x) = y.
 * @param x  Stores x values (more then one if not bijective).
 */
void SplineCurve::FindInverse(int segment, double y, GridPoints& x) const
{
  	assert((segment >= 0) && (segment < (signed)len-1));
  	x.clear();

  	const SplineSegmentCoefs& c = coefs[segment];
  	const double dx = grid[segment+1] - grid[segment];

 	double low, up;
  	GetMinMax(segment, low, up);
  	// check if segment intersects y
  	if ((y >= low) && (y < up)) {
    	//
    	// Solve f(x) = a (x - x0)^3 + b (x - x0)^2 + c (x - x0) + d = y
    	// Normalize such that z = (x - x0)/dx in [0, 1] and solve
    	//
    	//   f(z) = alpha z^3 + beta z^2 + gamma z + delta = 0
    	//
    	//   alpha = a * dx^3, beta = b * dx^2, gamma = c * dx, delta = d - y
    	const double alpha = c.a*dx*dx*dx;
    	const double beta = c.b*dx*dx;
    	const double gamma = c.c*dx;
    	const double delta = c.d - y;

    	Roots roots;

    	// ignore leading term (cubic term) if it is considerably smaller than
    	// the largest of all other terms. Quadratic roots are numerically much
    	// more stable. Also, if constant term is very small, z = 0 will be a
    	// solution and we can solve a quadratic equation as well.
    	const double ignore_leading_term = 1e-6;
    	const double ignore_constant_term = 1e-16;

    	if (fabs(delta) < ignore_constant_term) {
      		FindQuadraticRoots(beta/alpha, gamma/alpha, roots);
      		roots.push_back(0.0);
    	} else if (fabs(alpha)/(fabs(beta) + fabs(gamma) + fabs(delta)) < ignore_leading_term) {
			FindQuadraticRoots(gamma/beta, delta/beta, roots);
    	} else {
      		FindCubicRoots(beta/alpha, gamma/alpha, delta/alpha, roots);
    	}

		const double threshold = 1.0e-12;

    	// check if roots lie inside this segment
    	for (unsigned int i = 0; i < roots.size(); i++) {
      		if ((roots[i] >= 0.0 - threshold) && (roots[i] <= 1.0 + threshold)) {
        		x.push_back(grid[segment] + dx*roots[i]);
      		}
    	}
    
    	// ------------------------------------
    	// make roots unique 
    	// ------------------------------------
    	vector<double> ux(x);
    	vector<double>::iterator uend = unique(ux.begin(), ux.end());
    	x.assign(ux.begin(), uend);

    	// ----------------------------------------------
    	// [rv] o.k., the cubic spline is extremly sensitive, e.g.
    	// if y is close to low or high we could miss a root
    	// if this happens, we can try to use a newton solver to locate
    	// the root via deltax = - f(x) / f'(x)
    	// ----------------------------------------------
    	if(x.size() == 0) {
    		double newton_root = CubicNewton(beta/alpha, gamma/alpha, delta/alpha, 0.5);
    		if(newton_root >= 0.0 && newton_root <= 1.0) {
    			ostringstream sout;
    			sout <<  "standard inverse determination failed. could find one root at "
    			     << newton_root << " in [0,1] using a newton approach. affected "
    			     << " segment is nr. " << segment  
		  			 << " with values (min < y < max): " << low << " < " << y << " < " << up;
				Logger::get_instance()->emit(LOG_WARN, sout.str());
    			x.push_back(grid[segment] + dx*newton_root);	
    		}
    	}

#ifndef NDEBUG
    	// [mt] debug trace code
    	// should never occur, something is wrong if it does anyways
	    if (x.size() == 0) {	
	    	
      		cout << endl << endl;
      		cout << "from: " << roots.size() << " to: " << x.size() << endl;
		  	cout.precision(20);
		  	cout << "segment: " << segment << endl;
		  	cout << "(min < y < max): " << low << " < " << y << " < " << up << endl;
		  	cout << "(y0, y1): " << c.d << ", " << c.a*dx*dx*dx+c.b*dx*dx+c.c*dx+c.d << endl;
		  	cout << roots.size() << " root(s): ";
		  	cout.precision(20);
	      	for (unsigned int i = 0; i < roots.size()-1; i++) cout << roots[i] << ", ";
	      	cout << roots[roots.size()-1] << endl;
	      	cout.precision(20);
	      	cout << "(alpha, beta, gamma, delta): " << alpha << ", "
	             << beta << ", " << gamma << ", " << delta << endl;
	      	cout.precision(6);
	      	cout << "x0, x1, dx: " << grid[segment] << ", " << grid[segment+1] << ", " << dx << endl;
	      	cout << "coefs: " << "a=" << c.a << " b=" << c.b << " c=" << c.c << " d=" << c.d << endl;
	      	FindQuadraticRoots(c.c/(c.b*dx), (c.d - y)/(c.b*dx*dx), roots);
	      	cout << "quad-roots (# " << roots.size() << "): " << endl;
	      	for (unsigned int i = 0; i < roots.size()-1; i++) {
	      		TDKP_BOUNDS_ASSERT(segment >= 0 && segment < (signed)grid.size(), "segment >= 0 && segment < grid.size()");
	      	  	TDKP_BOUNDS_ASSERT(i >= 0 && i < roots.size(), "i >= 0 && i < roots.size()");
			  	cout << roots[i] << ": " << EvalAt(grid[segment]+dx*roots[i]) - y << ", ";		 
	      	}
	      	cout << roots[roots.size()-1] << ": " << EvalAt(grid[segment]+dx*roots[roots.size()-1]) - y << endl;
	      	cout << endl << endl;
	      	TDKP_GENERAL_EXCEPTION("there seems to be a problem with a spline interpolation. it's probably a round off error, but they can get severe, so check it and fix it.");
    	}
#endif
    	// [mt] check results: we know that at least one root must lie insdide this segment
    	// still leads to errors in some rare cases
    	//assert(x.size() >= 1);
#ifndef NDEBUG
    	for (unsigned int i = 0; i < x.size(); i++) {
      // allowed relative error for y compared to f(f^-1(y))
      //const double relative_error = 1e-6;
//      double ycheck = EvalAt(x[i]);
//      assert(fabs(ycheck - y) / max(fabs(y),fabs(ycheck)) < relative_error);
	    }
#endif

  	}
}

/** try to solve cubic equation using a newton approach
 * 
 * sometimes, the cubic roots function does file due to numerical
 * reasons. in that cases, sometimes looking for roots using a newton
 * approach may help ....
 * 
 * solves f(x) = x^3 + alpha * x^2 + beta * x + gamma = 0
 * up 1.0e-14.
 * 
 * does not validate the return. quits after 200 iterations
 *
 */ 
double SplineCurve::CubicNewton(const double& alpha, const double& beta, const double& gamma, double x) const {

    const double threshold = 1.0e-14;
    
    double deltax;
    int iterations = 200;
    do {    	
        deltax = - (x*x*x + alpha*x*x + beta*x + gamma) / (3.0*x*x + 2.0*alpha*x + beta);
        x += deltax;
    } while(fabs(deltax) > threshold && iterations-- > 0);
    
    return x; 
}


/** Stream a SplineCurve object.
 *
 * @param out  Stream object.
 * @param curve  Curve to be streamed.
 */
ostream& operator<<(ostream& out, const SplineCurve& curve)
{
  out << "Curve [SplineCurve]: " << curve.len << " points." << endl;
  out << "y(x) = a (x - x0)^3 + b (x - x0) ^2 + c (x - x0) + d     (x0 < x < x1)" << endl;
  out << endl;
  out << "  interval                       a              b              c              d" << endl;
  out << " ------------------------------------------------------------------------------------------" << endl;
  for (unsigned int i = 0; i < curve.len-1; i++) {
    char string[100];
    sprintf(string, "  %8.2E < x < %8.2E   %15.4E%15.4E%15.4E%15.4E", curve.grid[i], curve.grid[i+1],
        curve.coefs[i].a, curve.coefs[i].b, curve.coefs[i].c, curve.coefs[i].d);
    out << string << endl;
  }
  out << endl;
  return out;
}


inline double Chebyshev(double t)
{
  assert(t >= -2.0);

  if (fabs(t) <= 2.0) {
    return 2.0 * cos(acos(t/2.0)/3.0);
  } else {
    return 2.0 * cosh(acosh(t/2.0)/3.0);
  }
}

inline double Chebyshevi(double t)
{
  return 2.0 * sinh(asinh(t/2.0)/3.0);
}

/** Find real roots of function f(x) = x^3 + alpha*x^2 + beta*x + gamma = 0
 *
 * @param alpha  Alpha parameter.
 * @param beta  Beta parameter.
 * @param gamma  Gamma parameter.
 * @param x  Vector storing results.
 */
void FindCubicRoots(double alpha, double beta, double gamma, Roots& x)
{
  // Chebyshev's method which avoids complex numbers, as found in
  // http://en.wikipedia.org/wiki/Cubic_equation
  x.clear();

  // f(x) = x^3 - 3 p x - q = 0
  const double p = (alpha*alpha - 3.0*beta)/9.0;
  const double q = -(2.0*alpha*alpha*alpha  - 9.0*alpha*beta + 27.0*gamma)/27.0;
  // [mt] cout << "(alpha, beta, gamma): " << alpha << ", " << beta << ", " << gamma << endl;
  // [mt] cout << "(p, q): " << p << ", " << q << endl;

  if (p == 0.0) {
    // solve x^3 = q -> only one real solution
    if (q > 0.0) {
      x.push_back(pow(q, 1.0/3.0));
    } else {
      x.push_back(-pow(-q, 1.0/3.0));
    }
  } else if (p > 0.0) {
    const double t = q / (p * sqrt(p));
    // [mt] cout << "t: " << t << endl;
    if (t > 2.0) {
      // only one real solution
      x.push_back(sqrt(p)*Chebyshev(t) - alpha/3.0);
    } else if (t < -2.0) {
      // only one real solution
      x.push_back(-sqrt(p)*Chebyshev(-t) - alpha/3.0);
    } else {
      // three real solutions
      const double x1 = sqrt(p)*Chebyshev(t) - alpha/3.0;
      const double x2 = -sqrt(p)*Chebyshev(-t) - alpha/3.0;
      x.push_back(x1); x.push_back(x2);
      x.push_back(-(x1 + x2 + alpha));
    }
  } else {
    // only one real solution
    const double t = - q / (p * sqrt(-p));
    // [mt] cout << "t: " << t << endl;
    x.push_back(sqrt(-p)*Chebyshevi(t) - alpha/3.0);
  }
}


/** Find roots of function f(x) = x^2 + alpha*x + beta = 0
 *
 * @param alpha  Alpha parameter.
 * @param beta  Beta parameter.
 * @param x  Vector storing results.
 */
void FindQuadraticRoots(double alpha, double beta, Roots& x)
{
  x.clear();

  const double d = alpha*alpha - 4.0*beta;
  // determine if roots are real
  if (d > 0) {
    x.push_back((-alpha + sqrt(d)) / 2.0);
    x.push_back((-alpha - sqrt(d)) / 2.0);
  }
}

void SplineCurve::GetIntrinsicGridRange(double& xl, double& xr) const {
	xl = grid.front();
	xr = grid.back();	
}



} // end of namespace


