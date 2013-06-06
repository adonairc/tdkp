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

#ifndef SPLINES_H_
#define SPLINES_H_

#include "tdkp/common/all.h"

namespace tdkp {

template<class T, class V>
class GenericCurve {
public:
	virtual T operator()(const V& x) const = 0;
					
};

/** spline class for clc
 * 
 * that's my second order (not cubic!) spline class for clc
 * which does not need any linear equation solving. the spline  
 * is given by
 * \f$	S_{i}(x) = a_{i} * (x - x_{i})^2 + b_{i} * (x - x_{i}) + c_{i} \f$.
 * 
 * the coefficents are obtained as follows:
 * 	c_{i} = f_{i} (the values provided)
 *  b_{i} = 2 * (c_{i} - c_{i - 1}) / (x_{i} - x_{i - 1}) - b_{i - 1}
 *  a_{i} = (b_{i+1} - b_{i}) / (2.0 * (x_{i + 1} - x_{i}))
 * 
 */
class Spline1D : public GenericCurve<double,double> {
public:
	Spline1D(const vector<double>& x, const vector<double>& f);
	Spline1D(const vector<double>& x, const vector<double>& f, const double& dSdx);
			
	double operator()(const double& x) const;
	void eval(const vector<double>& x, vector<double>& res) const;
	
private:	
	Spline1D(const Spline1D& copy) { TDKP_GENERAL_EXCEPTION("copy forbidden"); }
	void init(const vector<double>& x, const vector<double>& f, const double& dSdx);

	/** data container for the spline variables */
	struct Spline1DSegment {
		double xi, ai, bi, ci;
	};
	vector<Spline1DSegment> segments;
	
	/** segment locator for fast spline segment location using a pivot tree approach */
	class SegmentLocatorTree {
	public:
		typedef vector<Spline1DSegment>::const_iterator segment_iterator;
		SegmentLocatorTree();
		~SegmentLocatorTree();
		void create_tree(segment_iterator left, segment_iterator right);
		const Spline1DSegment& get(const double& x) const;
		void  test() const;
	private:
		double x_border;	
		SegmentLocatorTree* trees[2];
		const Spline1DSegment* segments[2];							
	};
	
	SegmentLocatorTree root;
	double xmin;
	double xmax;
			
};


inline const Spline1D::Spline1DSegment& Spline1D::SegmentLocatorTree::get(const double& x) const {
	int idx = 1;
	// left or right
	if(x < x_border) {
		idx = 0;
	}
	// if left is a leaf, return segment, else ask child
	if(segments[idx] != 0) {
		return *segments[idx];	
	} else {
		return trees[idx]->get(x);	
	}
}	
	
inline double Spline1D::operator()(const double& x) const {
	
	TDKP_BOUNDS_ASSERT(x >= xmin - 1.0e-6, "x (" << x << ") >= xmin (" << xmin << "), xmax = " << xmax << " this = (" << this << ") ");	
	TDKP_BOUNDS_ASSERT(x <= xmax + 1.0e-6, "x (" << x << ") <= xmax (" << xmax << "), xmin = " << xmin << " this = (" << this << ") ");	
	const Spline1DSegment& segment = this->root.get(x);
	const double dx = (x - segment.xi);
	return segment.ai * dx * dx + segment.bi * dx + segment.ci; 
			
}

/** linear curve object for linear interpolation between complex numbers 
 *
 * the radial and the angle phi in the complex plane are interpolated 
 */
class CLCLinearCurve {
public:	
	CLCLinearCurve(const vector<double>& x, const vector<complex<double> >& f);
	complex<double> operator()(const double& x) const;
	void eval(const vector<double>& x, vector<complex<double> >& res) const;
private:
	vector<double> x;
	vector<double> f_r;
	vector<double> f_phi;
				
};

}
#endif /*SPLINES_H_*/
