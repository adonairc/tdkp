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

#include "tdkp/clc/CLCCoulombLambda.h"
#include "tdkp/common/Logger.h"


extern "C" {
	// acml vector function
	void vrda_cos(int, double *, double *);
	void vrda_sin(int, double *, double *);		
}

namespace tdkp
{

CLCRadialCoulombLambda::CLCRadialCoulombLambda(
	const vector<double>& q_values, 
	const vector<cplx>& values
) {
	
	TDKP_ASSERT(q_values.size() == values.size(), "");
	TDKP_ASSERT(q_values.size() > 2, "not enough q values for the spline interpolation ...");
	
	// ----------------------------------------
	// create spline BUT:
	// well clear, this guy diverges with q->0
	// but the spline needs a bc at left 
	// the splines goes from [xmin -> xmax].
	// but: lambda should converge to 0 for q->inf, so i simply
	// invert q -> - q ...
	// ----------------------------------------
	vector<double> inv_values(values.size());
	vector<double> inv_q_values(values.size());
	this->qmin = q_values.front();
	for(unsigned int ii = 0; ii < values.size(); ii++) {
		inv_values[ii] = values[values.size() - 1 - ii].real();
		TDKP_ASSERT(tdkp_math::abs(values[values.size() - 1 - ii].imag()) < 1.0e-8, "diagonal matrix elements with n1 = n4 and n2 = n3 should be real valued!");		
		inv_q_values[ii] = - q_values[values.size() - 1 - ii];		
	}	
	curve = new Spline1D(inv_q_values, inv_values, 0.0);
	
}



	
double CLCRadialCoulombLambda::get(const double& absq, const GenericCurve<double,double>* screening) const {
		
	TDKP_BOUNDS_ASSERT(absq >= 0.0, "");

	double qtmp = - absq;
	if(qmin > absq) {
		qtmp = - qmin;	
	} 

	if(screening != 0) {
		return (*curve)(qtmp) / (*screening)(absq);
	} else {
		return (*curve)(qtmp);
	}
}


CLCCoulombLambdaWire::CLCCoulombLambdaWire(
	const vector<double>& q_values, 
	const vector<cplx>& values
) : CLCRadialCoulombLambda(q_values,values) {}



/** return lambda_|k - k'| + lambda_|k + k'| */
double CLCCoulombLambdaWire::get_coulomb_lambda(const double& k, const double& kprime, const GenericCurve<double,double>* screening) const {
	return (this->get(tdkp_math::abs(k - kprime), screening) + this->get(tdkp_math::abs(k + kprime), screening));	
}

CLCRadialCoulombLambdaWell::CLCRadialCoulombLambdaWell(
	const vector<double>& q_values, 
	const vector<cplx>& values,
	unsigned int num_theta_points_
) : CLCRadialCoulombLambda(q_values, values), 
    num_theta_points(num_theta_points_),
    cos_theta(num_theta_points)    
{
	// ridicoulously small number, but anyway at least something
	TDKP_ASSERT(num_theta_points > 2, "");
	// precompute cos_theta	
#ifdef NOACML
	// no acml version
	for(unsigned int ii = 0; ii < num_theta_points; ii++) {
		cos_theta[ii] = cos((2.0 * constants::pi) / num_theta_points * ii);
	}
#else
	// acml version
	for(unsigned int ii = 0; ii < num_theta_points; ii++) {		
		cos_theta[ii] = (2.0 * constants::pi) / num_theta_points * ii;		
	}
	vrda_cos(num_theta_points, &cos_theta[0], &cos_theta[0]);
#endif
		 			
}


double CLCRadialCoulombLambdaWell::get_coulomb_lambda(const double& k, const double& kprime, const GenericCurve<double,double>* screening) const {
			
	// -----------------------------------------------
	// special handling for k = kprime (exclude q = 0)
	// cause: jellium model. if num_theta_points
	// is small enough, the integral should mathematically converge
	// -----------------------------------------------
	unsigned int theta_start = 0;
	if(tdkp_math::abs(k - kprime) < 1.0e-5) {
		theta_start = 1;		
	}
	// -----------------------------------------------
	// k and k prime are tooo small? ... integral over
	// k space would be zero (cause: k dk)
	// -----------------------------------------------
	if(tdkp_math::abs(k) < 1.0e-6 && tdkp_math::abs(kprime) < 1.0e-6) {
		return 0.0;	
	}	
	// --------------------------------------------
	// calculate q values
	// --------------------------------------------
	double k_square = k*k;
	double kprime_square = kprime * kprime;
	double ret = 0.0;
	double q = 0.0;
	double dtheta = 2.0 * constants::pi / (num_theta_points - 1);
	for(unsigned int ii = theta_start; ii < num_theta_points; ii++) {		
		q    = sqrt(k_square + kprime_square - 2.0 * k * kprime * cos_theta[ii]);
		ret += dtheta * this->get(q, screening);  
	}
	return ret;

}

CLCRadialCoulombLambdaBulk::CLCRadialCoulombLambdaBulk(
	const vector<double>& q_values, 
	const vector<cplx>& values,
	unsigned int num_points_	
) : CLCRadialCoulombLambda(q_values, values),
   	num_points(num_points_),
   	cos_phi(num_points_),
   	sin_theta(num_points_ / 2)
{
	// ridicoulously small number
	TDKP_ASSERT(num_points > 4, "");
		
	// -------------------------------------
	// prepare sin_theta and cos_phi 
	// -------------------------------------
#ifdef NOACML		
	for(unsigned int ii = 0; ii < cos_phi.size(); ii++) {
		cos_phi[ii] = cos((2.0 * constants::pi) / num_points * ii); 
	}
	for(unsigned int ii = 0; ii < sin_theta.size(); ii++) {
		sin_theta[ii] = sin(constants::pi / sin_theta.size() * ii); 
	}
#else
	for(unsigned int ii = 0; ii < cos_phi.size(); ii++) {
		cos_phi[ii] = (2.0 * constants::pi) / num_points  * ii; 
	}
	vrda_cos(cos_phi.size(), &cos_phi[0], &cos_phi[0]); 
	for(unsigned int ii = 0; ii < sin_theta.size(); ii++) {
		sin_theta[ii] = (constants::pi) / sin_theta.size() * ii; 
	}
	vrda_sin(sin_theta.size(), &sin_theta[0], &sin_theta[0]);
#endif

}

double CLCRadialCoulombLambdaBulk::get_coulomb_lambda(const double& k, const double& kprime, const GenericCurve<double,double>* screening) const {
	
	// -----------------------------------------------
	// k and k prime are tooo small? ... integral over
	// k space would be zero (cause: k^2 dk)
	// -----------------------------------------------
	if(tdkp_math::abs(k) < 1.0e-5 && tdkp_math::abs(kprime) < 1.0e-5) {
		return 0.0;	
	}
		
	const double k_square = k * k;
	const double k_prime_square = kprime * kprime;		
		
	// -----------------------------------------------
	// integrate over both angles
	// -----------------------------------------------
	double ret    = 0.0;
	const double dtheta = constants::pi / sin_theta.size();
	const double dphi   = 2.0 * constants::pi / cos_phi.size();
	unsigned int start_phi = 0;
	double q = 0.0;	  	
	
	for(unsigned int tt = 0; tt < sin_theta.size(); tt++) {
		// exclude q = 0!
		start_phi = sin_theta[tt] == 1.0 ? 1 : 0;			
		for(unsigned int pp = start_phi; pp < cos_phi.size(); pp++) {
			q = sqrt(k_square + k_prime_square - 2.0 * k * kprime * cos_phi[pp] * sin_theta[tt]);			  
			ret += sin_theta[tt] * dtheta * dphi * this->get(q, screening);	
		}
	}
	return ret;		
}

} // end of namespace
