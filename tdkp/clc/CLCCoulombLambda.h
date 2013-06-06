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

#ifndef CLCCOULOMBLAMBDA_H_
#define CLCCOULOMBLAMBDA_H_

#include "tdkp/common/all.h"
#include "tdkp/clc/Splines.h"

namespace tdkp {


/** coulomb lambda base class 
 * 
 * base container for the fourier transformed coulomb
 * interaction. the radial classes also contain the radial integration,
 * wire class contains return of N(ki,kj) + N(ki,-kj) 
 *
 * note, this classes do not correspond to Lambda in the documentation.
 * they correspond to N_dD(ki,kj)
 * 
 */
template<class T>
class CLCCoulombLambda {
public:
	/** return N_iD(k,k_prime) (q = k - k_prime, N_iD ~ Lambda_k-kprime,k,kprime) */
	virtual T get_coulomb_lambda(const double& k, const double& kprime, const GenericCurve<double,double>* screening = 0) const = 0;		
};


/** radial / diagonal (therefore real valued) fourier transformed coulomb potential base class
 * 
 * this class takes a diagonal set of q-dependent coulomb matrix elements
 * and generates a spline-interpolation routine for deriving subclasses
 * 
 * means: lambda_q_k_k' -> lambda_|k - k'|_0_0
 * 
 * this is the simplest approximation as
 * 1. the states should not change too much with increasing k
 * 2. and we only have lambda a long a specific direction but generalise it to
 *    the radial approximation
 * 
 * well / bulk / wire classes derive from this class and use the spline based
 * function for their possible integrations about the space-angles
 */ 
class CLCRadialCoulombLambda : public CLCCoulombLambda<double> {
public:
	CLCRadialCoulombLambda(const vector<double>& q_values, const vector<cplx>& values);
	
	virtual ~CLCRadialCoulombLambda() {
		delete curve;
	}	
	/** return raw and possibliy screened (not angularly integrated data) */	
	double get(const double& absq, const GenericCurve<double,double>* screening = 0) const;		
		
protected:
 	CLCRadialCoulombLambda() {};		
private:
	CLCRadialCoulombLambda(const CLCRadialCoulombLambda& copy) { TDKP_GENERAL_EXCEPTION("copy forbidden"); }
	Spline1D* curve;
	double    qmin;
};

/** coulomb class (diagonal elements) for quantum wires */ 
class CLCCoulombLambdaWire : public CLCRadialCoulombLambda {
public:
	CLCCoulombLambdaWire(const vector<double>& q_values, const vector<cplx>& values);
	double get_coulomb_lambda(const double& k, const double& kprime, const GenericCurve<double,double>* screening = 0) const;
private: 	
	CLCCoulombLambdaWire(const CLCCoulombLambdaWire& copy) { /* copy forbidden */ }
};



/** coulomb class (diagonal elements) for quantum wells 
 * 
 * don't return lambda but int_theta lambda_q(theta) dtheta
 * @param q_values q-values
 * @param values   values of Lambda_q
 * @param num_theta_points number of points to use in the integration
 */ 
class CLCRadialCoulombLambdaWell : public CLCRadialCoulombLambda {
public:
	CLCRadialCoulombLambdaWell(
		const vector<double>& q_values, 
		const vector<cplx>& values,
		unsigned int num_theta_points		
	);


	double get_coulomb_lambda(const double& k, const double& kprime, const GenericCurve<double,double>* screening = 0) const;
private:
	CLCRadialCoulombLambdaWell(const CLCRadialCoulombLambdaWell& copy) { /* copy forbidden */ } 	
	unsigned int num_theta_points;
	vector<double> cos_theta;
};


/** coulomb class (diagonal elements) for bulk 
 * 
 * don't return lambda but int_theta[0,pi] int_phi[0,2pi] lambda_q(theta,phi) sin_theta dtheta dphi
 * @param num_points we take num_points for the integration over phi and half of it for theta 
 */ 
class CLCRadialCoulombLambdaBulk : public CLCRadialCoulombLambda {
public:
	CLCRadialCoulombLambdaBulk(
		const vector<double>& q_values, 
		const vector<cplx>& values,
		unsigned int num_points	
	);	
	double get_coulomb_lambda(const double& k, const double& kprime, const GenericCurve<double,double>* screening = 0) const;
	
private:
	CLCRadialCoulombLambdaBulk(const CLCRadialCoulombLambdaBulk& copy) { /* copy forbidden */ }
	unsigned int num_points;
	vector<double> cos_phi;
	vector<double> sin_theta; 	
};


}

#endif /*CLCCOULOMBLAMBDA_H_*/
