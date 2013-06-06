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

#ifndef CLCSCREENING_H_
#define CLCSCREENING_H_

#include "tdkp/common/DataTypes.h"
#include "tdkp/clc/CLCIngredients.h"
#include "tdkp/clc/Splines.h"

namespace tdkp {

/** calculate screening epsilon(q,0) using lindhard formula */
class CLCScreening : public XYData<double>, public GenericCurve<double,double> {
public:
	CLCScreening(
		const CLCIngredients1DR& ingredients,
		const DomainMaster& dense_k_space	
	);
	virtual ~CLCScreening() {}
	void update(
		unsigned int num_cb_bands, unsigned int num_vb_bands, 
		const double& cb_fermi, const double& vb_fermi, const double& temperature
	);
	
	/** returns e(q,0) (if calculated via update, else returns 1) */
	inline double get(const double& q) const;
	/** operator inherited from GenericFunction */
	double operator()(const double& q) const { return get(q); }
	
	// -----------------------------------	
	// interface XYData
	// -----------------------------------
	virtual int    get_x_length() const;
	virtual int    get_num_y_sets() const;
	virtual int    get_num_x_sets() const; 
	virtual void   get_x(int xidx, vector<double>& x) const;
	virtual void   get_y(int yidx, vector<double>& y) const; 
	virtual string get_x_identifier(int xidx) const;
	virtual string get_y_identifier(int yidx) const;	
			 
private:
	double evaluate_longitudinal_polarization(
		const double& qq, const Spline1D& band, const Spline1D& fermi_stats
	) const;

	const CLCIngredients1DR& ingredients;
	const DomainMaster&      k_space;
	vector<double>  		 epsilon_screening;
	Spline1D*                epsilon_spline;
	vector<double>			 cos_phi;
	vector<double>			 sin_theta;
	double                   cb_fermi;
	double                   vb_fermi;
	double                   temperature;	
	unsigned int num_angular_integration_points;
	double                   dphi;
	double                   dtheta;
	double                   qmin;
	

		
};


/** returns e(q,0)*(1|q|q^2) (if calculated via update, else returns 1/q/q^2) */
inline double CLCScreening::get(const double& q) const {
	TDKP_BOUNDS_ASSERT(q >= 0.0, "");
	if(epsilon_spline != 0) {
		if(q < qmin) {
			return (*epsilon_spline)(-qmin);
		} else { 
			return (*epsilon_spline)(-q);
		}	
	} else {
		if(k_space.get_dimension() == 2) {		
			return q;	
		} else if(k_space.get_dimension() == 3) {
			return q*q;	
		} else {
			return 1.0;	
		}
	}	
}

}

#endif /*CLCSCREENING_H_*/
