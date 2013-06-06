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

#include "tdkp/clc/Splines.h"
#include "tdkp/common/Logger.h"

#ifndef NOACML
#include "acml_mv.h"
#endif

namespace tdkp {
	
/** create 1D spline, use Finite difference first derivative as starting value */	
Spline1D::Spline1D(const vector<double>& x, const vector<double>& f) {
	TDKP_ASSERT(x.size() > 1 && f.size() > 1,"");
	TDKP_ASSERT(x[0] < x[1],"x[0](" << x[0] << ") < x[1](" << x[1] << ")");
	this->init(x,f,(f[1] - f[0]) / (x[1] - x[0]));
}

/** create 1D spline with user supplied first derivative */	
Spline1D::Spline1D(const vector<double>& x, const vector<double>& f, const double& dSdx) {
	this->init(x,f,dSdx);
}

/** actual init function
 * 
 * calculates spline and builds search tree
 */		
void Spline1D::init(const vector<double>& x, const vector<double>& f, const double& dSdx) {

	TDKP_ASSERT(x.size() == f.size(), "x and f vectors must be of same length!");
	TDKP_ASSERT(x.size() > 1, "");
	
	segments.resize(x.size() - 1);
	segments[0].bi = dSdx;
	
	// ---------------------------------------------
	// test ordering of x
	// ---------------------------------------------
	for(unsigned int ii = 1; ii < x.size(); ii++) {
		TDKP_ASSERT(x[ii - 1] < x[ii],"");	
		segments[ii - 1].xi = x[ii - 1];			
		segments[ii - 1].ci = f[ii - 1];
	}
	
	xmin = x.front();
	xmax = x.back();
		
	// ---------------------------------------------
	// calculate bi coefficent
	// ---------------------------------------------
	for(unsigned int ii = 1; ii < segments.size(); ii++) {
		segments[ii].bi = 2.0 * (f[ii] - f[ii - 1]) / (x[ii] - x[ii - 1]) - segments[ii - 1].bi;	
	}
	
	// ---------------------------------------------
	// calculate ai coefficent
	// ---------------------------------------------
	for(unsigned int ii = 0; ii < segments.size() - 1; ii++) {
		segments[ii].ai = (segments[ii + 1].bi - segments[ii].bi) / (2.0 * (x[ii + 1] - x[ii]));	
	}
	
	// ---------------------------------------------
	// last ai coeff (because we have b[ii + 1]
	// ---------------------------------------------
	double bnp = 2.0 * (f[f.size() - 1] - f[f.size() - 2]) / (x[x.size() - 1] - x[x.size() - 2]) - segments.back().bi;
	segments.back().ai = (bnp - segments.back().bi) / (2.0 * (x[x.size() - 1] - x[x.size() - 2]));   

	// ---------------------------------------------
	// create tree
	// ---------------------------------------------
	this->root.create_tree(segments.begin(), segments.end() - 1);
	this->root.test();
	
		
}	


Spline1D::SegmentLocatorTree::SegmentLocatorTree()
: x_border(0.0)
{
	trees[0] = trees[1] = 0;
	segments[0] = segments[1] = 0;
}
	
Spline1D::SegmentLocatorTree::~SegmentLocatorTree() {
	for(short ii = 0; ii < 2; ii++) {
		if(trees[ii] != 0) {
			delete trees[ii]; trees[ii] = 0;	
		}
	}
}

/** recursive function to create a tree */
void Spline1D::SegmentLocatorTree::create_tree(segment_iterator left, segment_iterator right) {

	unsigned int nseg = (right - left) + 1;
	TDKP_ASSERT(nseg > 1,"");
	
	// two segments left
	if(nseg == 2) {
		segments[0] = &(*left);
		segments[1] = &(*right);
		x_border = (*right).xi;
	// three segments left
	} else if(nseg == 3) {
		// left terminates, right continues one level
		segments[0] = &(*left);
		trees[1]    = new SegmentLocatorTree();
		trees[1]->create_tree(left + 1, right);
		x_border = (*(left+1)).xi;
	// more than two segments left		
	} else {
		trees[0] = new SegmentLocatorTree();
		trees[1] = new SegmentLocatorTree();
		segment_iterator mid = left + (nseg / 2);
		trees[0]->create_tree(left, mid - 1);
		trees[1]->create_tree(mid, right);
		x_border = (*mid).xi;
	}

}

/** test if segment tree is o.k (recurisvely) */	
void Spline1D::SegmentLocatorTree::test() const {
	TDKP_ASSERT(trees[0] || segments[0], "left tree and segment is zero");
	TDKP_ASSERT(trees[1] || segments[1], "right tree and segment is zero");
	TDKP_ASSERT(segments[1] == 0 || tdkp_math::abs(segments[1]->xi - x_border) < 1.0e-7, "border is not at its correct position");
	if(trees[0] != 0) {
		trees[0]->test();	
	}
	if(trees[1] != 0) {
		trees[1]->test();	
	}			
}	

/** vectorial evaluation of spline */	
void Spline1D::eval(const vector<double>& x_, vector<double>& res) const {
	res.resize(x_.size());
	for(unsigned int ii = 0; ii < x_.size(); ii++) {
		res[ii] = (*this)(x_[ii]);	
	}	
}	
	
/** linear complex curve constructur */	
CLCLinearCurve::CLCLinearCurve(const vector<double>& x_, const vector<complex<double> >& f_)
: x(x_),
  f_r(x_.size()),
  f_phi(x_.size())
{
	
	TDKP_ASSERT(x_.size() == f_.size(),"");
	TDKP_ASSERT(x_.size() > 1, "");
	// -------------------------------------------
	// calcualte r and phi
	// -------------------------------------------
	for(unsigned int ii = 0; ii < x.size(); ii++) {		
		f_r[ii]   = abs(f_[ii]);
		if(f_r[ii] > 1.0e-13) {
			f_phi[ii] = acos(f_[ii].real() / f_r[ii]);
			if(f_[ii].imag() < 0.0) {
				f_phi[ii] = 2.0 * constants::pi - f_phi[ii];	
			}
		} else {
			f_phi[ii] = 0.0;	
		}	
	} 		
		
}

/** vectorial evaluator */
void CLCLinearCurve::eval(const vector<double>& x_, vector<complex<double> >& res) const {
	res.resize(x_.size());
	for(unsigned int ii = 0; ii < x_.size(); ii++) {
		res[ii] = (*this)(x_[ii]);	
	}	
}

/** single evaluator */
complex<double> CLCLinearCurve::operator()(const double& x_) const {
	
	TDKP_ASSERT(this != 0, "");
	
	// -----------------------------------
	// check that x is within our range
	// -----------------------------------
	TDKP_BOUNDS_ASSERT(x_ >= x.front() - 1.0e-6, "x_ (" << x_ << ") >= x.front() (" << x.front() << ")");
	TDKP_BOUNDS_ASSERT(x_ <= x.back() + 1.0e-6, "x_ (" << x_ << ") <= x.back() (" << x.back() << ")");
	
	// -----------------------------------
	// find correct segment
	// -----------------------------------	
	unsigned int idx = 0;
	while(idx < x.size() - 1 && x[idx + 1] < x_) {		
		idx++;
	}	
	double r = 0.0e0, phi = 0.0e0;
	// -----------------------------------
	// linearly interpolate r and phi
	// -----------------------------------
	if(idx == 0 && x[0] >= x_) {
		r   = f_r[0];
		phi = f_phi[0];								
	} else if(idx == x.size() - 1) {
		r   = f_r.back();
		phi = f_phi.back();	
	} else {
		double t = (x_ - x[idx]) / (x[idx + 1] - x[idx]);
		r   = f_r[idx]   * (1.0 - t) + f_r[idx + 1]   * t;
		phi = f_phi[idx] * (1.0 - t) + f_phi[idx + 1] * t;
	}
	// -----------------------------------
	// map back to complex value
	// -----------------------------------
#ifdef NOACML
	return cplx(r * cos(phi), r * sin(phi));
#else	
	return cplx(r * fastcos(phi), r * fastsin(phi));
#endif	
			
}


	
	
}
