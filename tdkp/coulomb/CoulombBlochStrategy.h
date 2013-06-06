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

#ifndef COULOMBBLOCHSTRATEGY_H_
#define COULOMBBLOCHSTRATEGY_H_

#include "tdkp/common/all.h"
#include "tdkp/coulomb/CoulombFunction.h"
#include "tdkp/coulomb/CoulombIntegrator.h"

namespace tdkp {

/** handle explicit form of coulomb integral
 * 
 */
class CoulombBlochStrategy {
public:
	CoulombBlochStrategy(unsigned int kp_basis_size_) : kp_basis_size(kp_basis_size_) {}
	virtual ~CoulombBlochStrategy() {}	
	template<class T>	
	cplx evaluate_M_product(const EigenSolution<cplx>* wave_left, 
		const EigenSolution<cplx>* wave_right, const CoulombIntegrator<T>& integrator) const;		
	virtual void prepare_coulomb_function(const double& q, const double* z_prime) = 0;		
	virtual const CoulombFunction& get_coulomb_function() const = 0;
	virtual CoulombBlochStrategy* clone_obj() const = 0;
private:
	unsigned int kp_basis_size;
};

/** basic strategy for wells */
class CoulombBlochStrategyWell : public CoulombBlochStrategy {
public:
	CoulombBlochStrategyWell(unsigned int kp_basis_size_) : CoulombBlochStrategy(kp_basis_size_) {}
	virtual void prepare_coulomb_function(const double& q, const double* z_prime);		
	virtual const CoulombFunction& get_coulomb_function() const { return coulomb_function; }
	virtual CoulombBlochStrategy* clone_obj() const { return new CoulombBlochStrategyWell(*this); }
private:
	CoulombFunctionWell coulomb_function;
};

/** basic strategy for wires */
class CoulombBlochStrategyWire : public CoulombBlochStrategy {
public:
	CoulombBlochStrategyWire(unsigned int kp_basis_size_) : CoulombBlochStrategy(kp_basis_size_) {}
	virtual void prepare_coulomb_function(const double& q, const double* z_prime);		
	virtual const CoulombFunction& get_coulomb_function() const { return coulomb_function; }
	virtual CoulombBlochStrategy* clone_obj() const { return new CoulombBlochStrategyWire(*this); }
private:
	CoulombFunctionWire coulomb_function;
};

/** basic strategy for dots */
class CoulombBlochStrategyDot : public CoulombBlochStrategy {
public:
	CoulombBlochStrategyDot(unsigned int kp_basis_size_) : CoulombBlochStrategy(kp_basis_size_){}
	virtual void prepare_coulomb_function(const double& q, const double* z_prime);		
	virtual const CoulombFunction& get_coulomb_function() const { return coulomb_function; }
	virtual CoulombBlochStrategy* clone_obj() const { return new CoulombBlochStrategyDot(*this); }
private:
	CoulombFunctionDot coulomb_function;
};

/** evaluate M product of two wave functions */
template<class T>
cplx CoulombBlochStrategy::evaluate_M_product(
	const EigenSolution<cplx>* wave_left, const EigenSolution<cplx>* wave_right, 
	const CoulombIntegrator<T>& integrator) const {
		
	// ----------------------------------------------
	// test compatibility
	// ----------------------------------------------
	const unsigned int vec_length = wave_left->get_length();
	unsigned int basis_size_left  = wave_left->get_basis_size();
	unsigned int basis_size_right = wave_right->get_basis_size();
	TDKP_ASSERT(basis_size_left == kp_basis_size || basis_size_left == 1, "basis_size_left == kp_basis_size || basis_size_left == 1");
	TDKP_ASSERT(basis_size_right == kp_basis_size || basis_size_right == 1, "basis_size_left == kp_basis_size || basis_size_left == 1");
	TDKP_ASSERT(wave_left->get_length() == wave_right->get_length(), "wave_left->get_length() == wave_right->get_length()");
	TDKP_ASSERT(wave_left->get_length() == (signed)integrator.get_matrix_size(), "wave_left->get_length() == integrator.get_matrix_size()"); 
	if(basis_size_left != basis_size_right) {
		Logger::get_instance()->emit(LOG_WARN, "dumb usage of coulomb matrix element calculation! matrix element is zero due to bloch lattice periodic part orthogonality!");
		return cplx(0.0,0.0);	
	}
	cplx res(0.0,0.0);
	// ----------------------------------------------
	// same basis size!
	// ----------------------------------------------
	const vector<cplx>& left_data  = wave_left->get_data();
	const vector<cplx>& right_data = wave_right->get_data();
	TDKP_ASSERT(basis_size_left *  vec_length == left_data.size(),  "basis_size_left (" << basis_size_left << ") * vec_length (" << vec_length << ") == left_data.size() (" << left_data.size() << ")");
	TDKP_ASSERT(basis_size_right *  vec_length == right_data.size(),  "basis_size_right (" << basis_size_right << ") * vec_length (" << vec_length << ") == right_data.size() (" << right_data.size() << ")");			
	 	
	vector<cplx> out(integrator.get_matrix_size());
	vector<cplx>  in(integrator.get_matrix_size());
	for(unsigned int bb = 0; bb < basis_size_left; bb++) {
		// extract rhs envelope		
		for(unsigned int ii = 0; ii < vec_length; ii++) {
			in[ii] = right_data[ii * basis_size_left + bb];
		}
		// calculate Mx		
		integrator.multiply_with_lhs(in, out);
		// dot product
		for(unsigned int ii = 0; ii < vec_length; ii++) {
			res += out[ii] * conj(left_data[ii * basis_size_left + bb]);		
		}
	}
	return res;
}

}

#endif /*COULOMBBLOCHSTRATEGY_H_*/
