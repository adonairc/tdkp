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

#include "tdkp/coulomb/CoulombFunction.h"
#include "tdkp/common/Logger.h"

#include <fstream>
using namespace std;

extern "C" {
	// acml vector function
	void vrda_exp(int, double *, double *);	
}

//#define NOACML

namespace tdkp {


double bessel_params_zaehler[3][4] = {
     // first set for x < 0.3
     {
           7.769569290104101355609, // c
        2437.574569063437138538574, // * x
       94045.294570112280780449509, // * x^2
      570178.807248712750151753425  // * x^3
     },
     // second set for 0.3 <= x < 2
     {
        4.4187208281850303137616, // c
       29.6682463040275301580095, // x
      - 5.1845243980118205229246, // x^2
        0.1593405524022399322170  // x^3
     },
     // third set for x > 2
     {
        1.5871254215963990219506,
      - 0.3642657873668007639800,
        0.0290940683932740563888,
      - 0.0008071702656177604917
     }
};

double bessel_params_nenner[3][7] = {
    // first set
    {
           1.0, // c
         478.10420462491129001137, // x
       28107.47316441927614505402, // x^2
      322767.74011291493661701679, // x^3
      705253.38399367441888898611, // x^4
      260890.79729612809023819864, // x^5
      256752.48985054454533383250  // x^6
    },
    // second set
    {
       1.0,
      16.9684716904316736929558,
      32.2573444444241701489772,
      11.9427488391509530885059,
       6.2596072550498789155426,
       0.2830919765714991487293,
       0.3150933548470242206995
    },
    // third set
    {
        1.0,
        0.669936750334212960389,
        1.380339728881989191933,
      - 0.196919657508525136613,
        0.159833576834217117035,
      - 0.016005752348467430085,
        0.002697497621511030230
    }
};


/** evaluate modified bessel function of second kind
 *
 * the bessel function is split into 3 segments
 * x0 = 0.005 x1 = 0.3 x2 = 2 x3 = 8
 * and evaluated using mathematica's RationalInterpolation[BesselK[0,x],{3,6},{x,x0,x1}]
 */
double BesselK0(double x) {

    double zaehler = 0.0;
    double nenner  = 0.0;
    int    idx     = 0;
    if(x > 3.0) {
        idx = 2;
    } else if(x > 0.3) {
        idx = 1;
    }
    double w = 1.0;
    // first three powers
    for(short ii = 0; ii < 4; ii++) {
        zaehler += bessel_params_zaehler[idx][ii] * w;
        nenner  += bessel_params_nenner[idx][ii] * w;
        w *= x;
    }
    for(short ii = 4; ii < 7; ii++) {
        nenner  += bessel_params_nenner[idx][ii] * w;
        w *= x;
    }    
    return zaehler / nenner;
}	
		
CoulombFunction::CoulombFunction(unsigned int dimension_)
: dimension(dimension_)
{
}

CoulombFunctionWell::CoulombFunctionWell()
: CoulombFunction(1),
  z_prime(0.0),
  q(0.0)
{
}
	
/** 2pi exp(-qx) (x = |z - z'|) */
void CoulombFunctionWell::evaluate(const vector<double>& z, vector<double>& result) const {
		  	
	TDKP_BOUNDS_ASSERT(result.size() == z.size(), "result.size() == z.size()");   		  		
	for(unsigned int ii = 0; ii < z.size(); ii++) {
		result[ii] = - q * fabs(z[ii] - z_prime);
#ifdef NOACML
		result[ii] = 2.0 * constants::pi  * exp(result[ii]);		
#endif	
	}	
#ifndef NOACML
	vrda_exp(z.size(), &result[0], &result[0]);
	for(unsigned int ii = 0; ii < z.size(); ii++) {
		result[ii] *= 2.0 * constants::pi;
	}
#endif
		 	
}

void CoulombFunctionWell::set(const double& q_, const double z_prime_) {
	q       = q_;
	z_prime = z_prime_;		
}
	
CoulombFunctionWire::CoulombFunctionWire()
: CoulombFunction(2),
  q(0.0)
{
	z_prime[0] = z_prime[1] = 0.0e0;
}

  
void CoulombFunctionWire::evaluate(const vector<double>& z, vector<double>& result) const {
	
	unsigned int num = result.size();
	TDKP_BOUNDS_ASSERT(result.size() == z.size() / 2, "result.size() == z.size() / 2");
	double zd1, zd2;	
	for(unsigned int ii = 0; ii < num; ii++) {
		zd1 = z[ii * 2]     - z_prime[0];
		zd2 = z[ii * 2 + 1] - z_prime[1];
		result[ii] = 2.0 * BesselK0(q * sqrt(zd1 * zd1 + zd2 * zd2)); 		
	}	
	
}

void CoulombFunctionWire::set(const double& q_, const double z_prime_[2]) {
	q = q_;
	z_prime[0] = z_prime_[0];
	z_prime[1] = z_prime_[1];	
}
	
CoulombFunctionDot::CoulombFunctionDot()
: CoulombFunction(3)
{
	z_prime[0] = z_prime[1] = z_prime[2] = 0.0e0;
}
  
void CoulombFunctionDot::evaluate(const vector<double>& z, vector<double>& result) const {
	
	unsigned int num = result.size();
	TDKP_BOUNDS_ASSERT(result.size() == z.size() / 3, "result.size() == z.size() / 2");
	double zd1, zd2, zd3;
	static bool warned = false;	
	for(unsigned int ii = 0; ii < num; ii++) {
		zd1 = z[ii * 3]     - z_prime[0];
		zd2 = z[ii * 3 + 1] - z_prime[1];
		zd3 = z[ii * 3 + 2] - z_prime[2];
		zd1 = sqrt(zd1 * zd1 + zd2 * zd2 + zd3 * zd3);
		if(zd1 != 0.0) {
			result[ii] = 1.0 / zd1;	
		} else {
			if(!warned) {
				Logger::get_instance()->emit(LOG_WARN, "coulomb approx assumes that z is never z'!");	
			}
			result[ii] = 0.0;	
		}
	}		
}

void CoulombFunctionDot::set(const double& q_, const double* z_prime_) {	
	z_prime[0] = z_prime_[0];
	z_prime[1] = z_prime_[1];
	z_prime[2] = z_prime_[2];	
}
	
	
	
} // end of namespace
