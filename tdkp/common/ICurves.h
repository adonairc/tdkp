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

#ifndef ICURVES_H_
#define ICURVES_H_
#include <iostream>
#include <vector>
 
#include "tdkp/common/all.h"
#include "tdkp/common/DataTypes.h" 
 
namespace tdkp {

using namespace std;

// to store grid base point(s)
typedef vector<double> GridPoints;
// to store values of the curve at certain grid base point(s)
typedef vector<double> GridValues;

/** purely abstract class defining interface to a curve object 
 *  
 * marco tomamichels spline and curve class
 */
class ICurve : public XYData<double> {
	
public:
  	// actual construction of curve
  	virtual void CopyFrom(const ICurve& curve) = 0;
  	virtual void CopyFrom(const GridPoints& grid, const GridValues& val) = 0;
  
  	virtual ~ICurve() {};
  
  	// evaluate curve (interpolate if necessary)
  	virtual double EvalAt(const double& x) const = 0;
  	virtual void   EvalAt(const GridPoints& x, GridValues& y) const = 0;
  	virtual double EvalDerivativeAt(const double& x) const = 0;
  	virtual void   EvalDerivativeAt(const GridPoints& x, GridValues& y) const = 0;
  	virtual double Integrate() const = 0;
  	virtual double Integrate(double x1, double x2) const = 0;
  	
  	virtual void   FindInverse(const double& y, GridPoints& x) const = 0;
  	
  
  	// extract underlying base grid and values
  	// (= curve values that do not need to be interpolated)
  	virtual void   GetIntrinsicGrid(GridPoints& x) const = 0;
	virtual void   GetIntrinsicGridRange(double& xl, double& xr) const = 0; 
	virtual const double& GetIntrinsicGridMin() const = 0;
	virtual const double& GetIntrinsicGridMax() const = 0; 	
  	
  	virtual const GridPoints& GetGridReference() const = 0;
  
  	// --------------------------------------------------
  	// XYData interface
  	// --------------------------------------------------
	int    get_x_length()                const;
	int    get_num_y_sets()              const;
	int    get_num_x_sets()              const { return 1; } 
	void   get_x(int xidx, vector<double> &x) const;
	void   get_y(int yidx, vector<double>& y) const; 
	string get_x_identifier(int xidx)    const;
	string get_y_identifier(int yidx) 	 const;	  
  
};

class LinearCurve;

// simple curve object that stores values at certain grid base points
// and does linear interpolation in between

class LinearCurve : public ICurve {
			
public:

  // standard constructors
  LinearCurve() { len = 0; }
  LinearCurve(const ICurve& curve) { CopyFrom(curve); }
  LinearCurve(const LinearCurve& curve) { CopyFrom(curve); }
  LinearCurve(const GridPoints& x, const GridValues& y) { CopyFrom(x, y); }
  // standard destructor
  ~LinearCurve() { Cleanup(); }

  virtual void CopyFrom(const ICurve& curve);
  virtual void CopyFrom(const GridPoints& x, const GridValues& y);

  virtual double EvalAt(const double& x) const;
  virtual void   EvalAt(const GridPoints& x, GridValues& y) const;
  virtual double EvalDerivativeAt(const double& x) const;
  virtual void   EvalDerivativeAt(const GridPoints& x, GridValues& y) const;
  virtual double Integrate() const;
  virtual double Integrate(double x1, double x2) const;
  virtual void   FindInverse(const double& y, GridPoints& x) const;

  virtual void   GetIntrinsicGrid(GridPoints& x) const { x = grid; };  
  virtual void   GetIntrinsicGridRange(double& xl, double& xr) const;
  virtual const double& GetIntrinsicGridMin() const { return grid.front(); }
  virtual const double& GetIntrinsicGridMax() const { return grid.back(); }  
  const GridPoints& GetGridReference() const { return grid; }  

  // assignment operator
  const ICurve& operator=(const ICurve& curve) { CopyFrom(curve); return *this; }
  // stream output operator
  friend ostream& operator<<(ostream& out, const LinearCurve& curve);

  // -------------------------------------
  // some operators for easy manipulation
  // -------------------------------------
  // think before you write formulas!
  // e.g. (double * double * double * double / double) * (curve)
  // instead of double * (double * (double * curve))
  LinearCurve operator*(const double& mul) const;
  LinearCurve operator/(const double& mul) const;
  LinearCurve operator*(const ICurve& rhs) const; 
  LinearCurve sqrt() const; 
     
protected:

  void Cleanup();
  void CalcDerivative();
  
  virtual double IntegrateSegment(int seg, double x1, double x2) const;

  size_t len; 
  GridPoints grid;
  GridValues beta;  // values of curve at grid points
  // size of this array is len-1
  GridValues alpha;  // derivative in each segment

private:


};



// spline curve object that stores curve as a cubic spline and
// uses this spline to interpolate between base points

class SplineCurve : public ICurve {
public:

  typedef enum {
    bc_natural = 1,
    bc_notaknot = 2,
    bc_standard = bc_notaknot
  } BoundaryCondition;

  // standard constructors
  SplineCurve() { len = 0; bc = bc_standard; }
  SplineCurve(BoundaryCondition abc) { len = 0; bc = abc; }
  SplineCurve(const ICurve& curve) { bc = bc_standard; CopyFrom(curve); }
  SplineCurve(const GridPoints& x, const GridValues& y) { bc = bc_standard; CopyFrom(x, y); }
  // standard destructor
  ~SplineCurve() { Cleanup(); }

  virtual void CopyFrom(const ICurve& curve);
  virtual void CopyFrom(const GridPoints& x, const GridValues& y);

  virtual double EvalAt(const double& x) const;
  virtual void   EvalAt(const GridPoints& x, GridValues& y) const;
  virtual double EvalDerivativeAt(const double& x) const;
  virtual void   EvalDerivativeAt(const GridPoints& x, GridValues& y) const;
  virtual double Integrate() const;
  virtual double Integrate(double x1, double x2) const;
  virtual void   FindInverse(const double& y, GridPoints& x) const;

  virtual void GetIntrinsicGrid(GridPoints& x) const { x = grid; };
  const GridPoints& GetGridReference() const { return grid; }
  virtual void   GetIntrinsicGridRange(double& xl, double& xr) const;
  virtual const double& GetIntrinsicGridMin() const { return grid.front(); }
  virtual const double& GetIntrinsicGridMax() const { return grid.back(); }
  
  // assignment operator
  const ICurve& operator=(const ICurve& curve) { CopyFrom(curve); return *this; }
  // stream output operator
  friend ostream& operator<<(ostream& out, const SplineCurve& curve);

protected:

  void Cleanup();
  void CalcSplines(const GridValues& val);

  // to be used by DOSBroadeningTool
  double GetMinMax(int segment, double& min, double& max) const;
  void   FindInverse(int segment, double y, GridPoints& x) const;  
  double CubicNewton(const double& alpha, const double& beta, const double& gamma, double x) const;
  
  typedef struct {
    double a, b, c, d;
  } SplineSegmentCoefs;
  typedef vector<SplineSegmentCoefs> SplineCoefs;

  // trigonal matrix used to calculate spline coefficients
  // Ax = b
  typedef struct {
    double low, diag, up;
  } TrigMatrixLine;
  typedef vector<TrigMatrixLine> TrigMatrix;
  typedef vector<double> TrigVector;

  // during Gaussian transformations we will write into all matrices/vectors
  void SolveTrigonal(TrigMatrix& A, TrigVector& b, TrigVector& x);

  BoundaryCondition bc;
  size_t len;
  GridPoints grid;
  SplineCoefs coefs;

private:
	
	void prepare_segments();
	/** segments min max values needed for fast inverse calculations */	  		  		
	struct MinMaxValue {
		double min_value;
		double max_value;			
	};
	vector<MinMaxValue> segment_ranges;

};


ostream& operator<<(ostream& out, const LinearCurve& curve);
ostream& operator<<(ostream& out, const SplineCurve& curve);




}

#endif /*CURVES_H_*/
