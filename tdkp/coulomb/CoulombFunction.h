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

#ifndef COULOMBFUNCTION_
#define COULOMBFUNCTION_

#include "tdkp/common/all.h"

namespace tdkp {

/** basic function that will be integrated on the grid */
class CoulombFunction {
public:	
	virtual ~CoulombFunction() {}
	virtual void evaluate(const vector<double>& z, vector<double>& result) const = 0;
	unsigned int get_dimension() const { return dimension; }	
protected:
	CoulombFunction(unsigned int dimension);		
private:
	unsigned int dimension;	
};

/** form factor / fourier transform (multiplied by q) of coulomb potential for quantum wells
 * 
 * 
 * (the zeta * q in my documentation)
 * given by 2pi  * exp(-qx) (x = |z - z'|)
 * the * q is used to prevent divergence in the limit of q->0 
 */
class CoulombFunctionWell : public CoulombFunction {
public:
	CoulombFunctionWell();
	void evaluate(const vector<double>& z, vector<double>& result) const;
	void set(const double& q, const double z_prime_);
private:
	double z_prime;
	double q;	
};

/** form factor / fourier transform of coulomb potential for quantum wires
 *
 * (the zeta in my documentation) 
 * given by 2 * K0(|z - z'|q)
 */
class CoulombFunctionWire : public CoulombFunction {
public:
	CoulombFunctionWire();
	void evaluate(const vector<double>& z, vector<double>& result) const;
	void set(const double& q, const double z_prime_[2]);
private:
	double z_prime[2];
	double q;		
};

/** coulomb potential for quantum dots
 *
 * (the zeta in my documentation) 
 * given by 1/|z-z'|, the case z = z' is
 * currently ignored! (so set to 0) but that shouldn't be too band
 * cause i assume that 1/|z-z|' is constant over the element and
 * z is mid element while z' is node ...
 */
class CoulombFunctionDot : public CoulombFunction {
public:
	CoulombFunctionDot();
	void evaluate(const vector<double>& z, vector<double>& result) const;
	void set(const double& q, const double* z_prime_);
private:
	double z_prime[3];	
};
	
} // end of namespace

#endif /*COULOMBFUNCTION_*/
