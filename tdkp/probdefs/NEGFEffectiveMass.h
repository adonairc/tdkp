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

#ifndef NEGFEFFECTIVEMASS_H_
#define NEGFEFFECTIVEMASS_H_

#include "tdkp/probdefs/EigenProblem3D.h"
#include "tdkp/probdefs/EffectiveMass.h"

namespace tdkp {

/** effective mass decorator to create combined cb <-> vb hamiltonians for NEGF calculations 
 * 
 * the hamiltonian has a simple form of a 2x2 matrix with only diagonal entries, where the
 * upper entries correspond to the cb band and the lower denote the vb band 
 */
class NEGFEffectiveMass : public SchroedingerProblem<EigenProblem3D<complex<double>, complex<double> > > {
public:
	NEGFEffectiveMass(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~NEGFEffectiveMass();
	// -----------------------------------------------------------
	// for communication between FEMSolver and problem class
	// -----------------------------------------------------------
	virtual const int* get_node_sparsity_pattern(int &num) const;   	
	virtual void calculate_element_matrices(const Element* elem, complex<double>* lhs, double *rhs, int* node_internal_indices, int &n) const;   	
	virtual void prepare();
	virtual void solve(int num_solutions);
	virtual string get_unique_identifier() const { return string("NEGFEffectiveMass"); }
	virtual EigenProblemType get_problem_type() const { return indefinite; } 
	
	virtual void set_energy_shift(double energy);	
	virtual void drop_energy_shift();			

	/** dummy function, called by negf hamiltonian assembler (for kp problems) */
	void set_axes(const Vector3D& dummy_a, const Vector3D& dummy_b) const {}
	
	/** set transversal k for negf hamiltonian */
	void set_k_transversal(const double& k_transversal_) { k_transversal = k_transversal_; }
	
protected:
 	
	
private:
	EffectiveMass cb_effmass;
	EffectiveMass vb_effmass;	
	int node_sparsity_pattern[4];
	bool   energy_shift_set;
	double energy_shift;
	double get_energy_shift() const;
	double k_transversal;
	 		
};

}

#endif /*NEGFEFFECTIVEMASS_H_*/
