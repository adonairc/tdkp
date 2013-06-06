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

#ifndef SCHROEDINGERPOISSON_H_
#define SCHROEDINGERPOISSON_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Domain.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/probdefs/PoissonEquation.h"
#include "tdkp/main/Bandstructure.h"
#include "tdkp/utilities/GridInterpolator.h"

namespace tdkp {

template<class CBPB, class VBPB>
class SchroedingerPoisson {
public:
	SchroedingerPoisson(
		CBPB& cb_problem_,
		VBPB& vb_problem_,
		PoissonEquation& poisson_equation_,
		const DomainMaster& kspace
	);	
	virtual ~SchroedingerPoisson();
	
	bool calculate(
		unsigned int num_cb_bands, unsigned int num_vb_bands, 
		const double& electron_density, const double& hole_density, const double& temperature
	);
		
	const BandstructureDomain<complex<double> >& get_cb_bands() const;
	const BandstructureDomain<complex<double> >& get_vb_bands() const;
	
	const BandstructureDomain<complex<double> >& get_cb_dispersion() const;
	const BandstructureDomain<complex<double> >& get_vb_dispersion() const;
				
private:
		
	const DomainMaster& get_domain() const { return kspace_domain; }
	
	void calculate_cb_bandstructure();
	void calculate_vb_bandstructure();
	void post_cb_bandstructure();
	void post_vb_bandstructure();
	void init_relaxation();
			
	void calculate_fermi_levels();
	void calculate_carrier_distribution(bool predictor_corrector);
	void calculate_potential();	
	void update_kp_potential();
	bool check_convergence(bool predictor_corrector_mode = false);
	void do_predictor_corrector();
	void potential_update_with_slowdown(double& update_ratio);

	double get_cb_transverse_effective_mass() const;
	double get_vb_transverse_effective_mass() const;


	CBPB& 					cb_problem;
	VBPB& 					vb_problem;
	PoissonEquation&        poisson_equation;	
	const Geometry&         kp_geometry;
	const Geometry&         poisson_geometry;
	const DomainMaster      kspace_domain;
	
	PotentialEnergyField    potential_energy_field;
	GridInterpolator        interpolator_kp_to_poisson;
	GridInterpolator        interpolator_poisson_to_kp;
	const NodeData<double>* default_nodal_charges;
	StdNodeData<double>     kp_current_nodal_charges;
	StdNodeData<double>     poisson_current_nodal_charges;
	NodeData<double>*		current_poisson_solution;
	NodeData<double>*		last_poisson_solution;
	NodeData<double>*       predictor_corrector_reference_solution;
	
	const BandstructureDomain<complex<double> >* cb_bands;
	const BandstructureDomain<complex<double> >* vb_bands;	
	
	BandstructureDomain<complex<double> >* cb_effmass_disp;
	BandstructureDomain<complex<double> >* vb_effmass_disp;

	double electron_density;
	double hole_density;
	double temperature;	
	double cb_fermi;
	double vb_fermi;
	unsigned int number_of_valid_cb_subbands;
	unsigned int number_of_valid_vb_subbands;
	double cb_min_edge;
	double cb_max_edge;
	double vb_min_edge;
	double vb_max_edge;
	unsigned int num_cb_bands;
	unsigned int num_vb_bands;
	
	vector<double> cb_distribution;
	vector<double> vb_distribution;

	double poisson_solution_difference_norm;
	double poisson_solution_last_difference_norm;
	double poisson_solution_solution_norm;
	double sg_update_ratio;
	const double min_update_ratio;
	const double max_update_ratio;
	const double success_increase_ratio;
	const double failure_decrease_ratio;
	
	bool  include_lda_xc;
		
};


}

#include "tdkp/utilities/SchroedingerPoisson.tcc"

#endif /*SCHROEDINGERPOISSON_H_*/
