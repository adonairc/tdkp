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

#ifndef SCHROEDINGERPML_H_
#define SCHROEDINGERPML_H_

#include "tdkp/probdefs/EigenProblem.h"
#include "tdkp/geometry/ElementPML.h"

namespace tdkp {

class PMLStretch;

/** wrapper class for solving schroedinger equations including perfectly matched layers at domain boundaries
 * 
 * pmls can be used to mimic open boundaries. the result is that we get
 * an imaginary part in the energy of the state beeing proportional to the
 * lifetime of the state.
 * 
 * technically we create these pmls by stretching the coordinates into the complex plane,
 * i.e we transform in our integrals into x -> z = S(x). so the path in the complex plane
 * is parametrized using the real space x. the result is that
 * dz = dS(x)/dx dx = s(x)dx and d/dz -> 1/s(x) dx
 * 
 * in the present implementation, we offer s(x) = 1 + (alpha + i beta) * ((x - x0) / pml_length)^n
 *
 * the problem is: pmls lead to a very complex eigenequation, solving that is 
 * very cumbersome and costly. also matrix assembly time is more costly as 
 * we need to perform numerical intergration of the element integrals including
 * the complex coordinate stretch.
 * 
 * pml's must be oriented along a main coordinate axis.
 * 
 * so, here is now how you define a PML:
 * name your region PML_nr_px for a pml that extends from x0 to PLUS infinity
 * name your region PML_nr_mx for a pml that extends from MINUS infinity to x0
 * name your region PML_nr_pxmy and so on ...
 * nr is just a number helping you to give different region names 
 */ 
template<class PC>
class SchroedingerPML : public EigenProblem<complex<double>, complex<double>, complex<double> >{
public:
	SchroedingerPML(
		PC& problem
	);
	virtual ~SchroedingerPML();
	
	// ----------------------------------------------------
	// intercommunication with fem solver
	// ----------------------------------------------------
	virtual EigenProblemType get_problem_type() const;
	virtual void display_solution_info() const;	
	virtual void add_solution(cplx solution_value, const cplx* solution_vector, int length);
	virtual const int* get_node_sparsity_pattern(int &num) const;   	
	virtual void calculate_element_matrices(const Element* elem, cplx* lhs, cplx* rhs, int* node_internal_indices, int &n) const;   	
	virtual void prepare();
	virtual void solve(int num_solutions) = 0;
	virtual void solve(int num_subbands, double kmin, double kmax, int num_k_values) { TDKP_GENERAL_EXCEPTION("solve function not implemented in derived class. probably you should not use this one"); }	
	virtual void solve(int num_subbands, const DomainMaster& domain) { TDKP_GENERAL_EXCEPTION("solve function not implemented in derived class. probably you should not use this one"); }
	virtual string get_unique_identifier() const;
		
	// ----------------------------------------------------
	// pml controls
	// ----------------------------------------------------
	void set_stretch_parameters(const double& alpha, const double& beta, unsigned int pml_power);
	void dump_pml_functions(const char* filename) const;
	
protected:
	
	void build_coordinate_stretches();
	void delete_stretches();
	
	void build_element_wrappers();
	void delete_pml_element_wrappers();

	PC& base_problem;	
	vector<PMLStretch*> stretch[3]; // stretch for every region!
	double default_pml_alpha;
	double default_pml_beta;
	unsigned int default_pml_power;
	
	vector<ElementPML*> pml_element_wrappers;
	vector<bool>        pml_regions;
	
		
};

template<class PC>
class SchroedingerPML3D : public SchroedingerPML<PC> {
public:
	SchroedingerPML3D(
		PC& problem
	);
	virtual ~SchroedingerPML3D() {}
	virtual void solve(int num_solutions);
		
};

template<class PC>
class SchroedingerPML1D2D : public SchroedingerPML<PC> {
public:
	SchroedingerPML1D2D(
		PC& problem
	);
	virtual ~SchroedingerPML1D2D() {}
	virtual void solve(int num_solutions);
	virtual void solve(int num_subbands, double kmin, double kmax, int num_k_values);	
	virtual void solve(int num_subbands, const DomainMaster& domain);
	
};

/** pml complex coordinate stretch
 * 
 * we have complex coordinates z and express it in terms of 
 * the real coordinate x. 
 * this leads to dz = s(x)dx
 * and s(x) = (1 + (alpha + i beta)*(x-pml_start)/(pml_end - pml_start)^pml_power
 * 
 */ 
class PMLStretch {
public:
	PMLStretch(
		const double& pml_start, 
		const double& pml_end,
		const double& alpha, 
		const double& beta,
		unsigned int pml_power
	);
	complex<double> evaluate(const double& x) const;
private:
	double pml_start; 
	double pml_end;
	complex<double> fact;
	unsigned int pml_power;		
};

}

#include "tdkp/probdefs/SchroedingerPML.tcc"

#endif /*SCHROEDINGERPML_H_*/
