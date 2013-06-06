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

#ifndef EFFECTIVEMASS_H_
#define EFFECTIVEMASS_H_

#include "tdkp/probdefs/EigenProblem3D.h"

namespace tdkp {


/** using eff mass calculations and a k-space domain, we return a dispersive bandstructure (without envelopes!) */
BandstructureDomain<complex<double> >* 
create_effmass_dispersion(const double& transverse_effmass,
						  const DomainMaster& domain,
                          const Bandstructure<complex<double> >& bands);
						   

/** effective mass calculation of quantized regions 
 * 
 * this class works for wells, wires and dots
 */
class EffectiveMass : public SchroedingerProblem<EigenProblem3D<complex<double>, double> >  {
public:	
	EffectiveMass(const Geometry& geometry_, MaterialDatabase& material_database_);	
	virtual ~EffectiveMass();
	// -----------------------------------------------------------
	// for communication between FEMSolver and problem class
	// -----------------------------------------------------------
	virtual const int* get_node_sparsity_pattern(int &num) const;   	
	virtual void calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* node_internal_indices, int &n) const;   	
	virtual void prepare();
	virtual void solve(int num_solutions);
	virtual void add_solution(cplx solution_value, const cplx* solution_vector, int length);
	virtual string get_unique_identifier() const { return string("EffectiveMass"); }

	void           set_solution_type(KPSolutionType type_) { this->solution_type = type_; this->ready = false; }
	KPSolutionType get_solution_type() const { return this->solution_type; }
		
	virtual void set_energy_guess(double energy);	
	virtual void drop_energy_guess();			
	
	virtual EigenProblemType get_problem_type() const;	
	
	StdElementData<double>* get_bandedges() throw(Exception*);
	void get_minmax_edges(double& cb_min, double& cb_max, double& vb_min, double& vb_max) const;
	double get_energy_shift() const;

	void set_hydro_strain_potential_field_name_cb(unsigned int axis, const char* field_name);
	void set_hydro_strain_potential_field_name_vb(unsigned int axis, const char* field_name);
		
protected:
	class RegionProperties {
	public:
		RegionProperties() : inv_effective_mass(3,3), particle_potential_energy(0.0) 
		{ 
				hydrostatic_strain_potential[0] = hydrostatic_strain_potential[1] =  hydrostatic_strain_potential[2] = 0.0; 
		}	
		RMatrix<double> inv_effective_mass;
		double particle_potential_energy;
		double hydrostatic_strain_potential[3];	
	};

	static const int sparsity_pattern[2];
	int*			 sparsity_copy;
	
	vector<RegionProperties> properties;	
	KPSolutionType solution_type;
	
	bool   energy_guess_set;
	double energy_guess;	
	
	/** name of hydrostatic cb strain potential for given axis in material object */	
	string hydro_strain_potential_field_name_cb[3];
	/** name of hydrostatic vb strain potential for given axis in material object */
	string hydro_strain_potential_field_name_vb[3];		
			
};

}

#endif /*EffectiveMass3D_H_*/
