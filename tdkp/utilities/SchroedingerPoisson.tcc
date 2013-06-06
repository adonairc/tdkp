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


#include "tdkp/utilities/Fermi.h"
#include "tdkp/probdefs/KP8x83D.h"
#include "tdkp/probdefs/KP4x43D.h"
#include "tdkp/probdefs/KPBase1D2D.h"
#include "tdkp/probdefs/KP4x41D2D.h"
#include "tdkp/probdefs/KP6x61D2D.h"
#include "tdkp/probdefs/KP8x81D2D.h"
#include "tdkp/probdefs/EffectiveMass.h"

namespace tdkp {

template<class CBPB, class VBPB>
SchroedingerPoisson<CBPB,VBPB>::SchroedingerPoisson(
	CBPB& cb_problem_,
	VBPB& vb_problem_,
	PoissonEquation& poisson_equation_,
	const DomainMaster& kspace_
)
: cb_problem(cb_problem_),
  vb_problem(vb_problem_),
  poisson_equation(poisson_equation_),
  kp_geometry(vb_problem.get_geometry()),
  poisson_geometry(poisson_equation.get_geometry()),
  kspace_domain(kspace_),
  potential_energy_field(kp_geometry.get_num_nodes()),
  interpolator_kp_to_poisson(kp_geometry, poisson_geometry, identity_matrix, Vector3D(0.0, 0.0, 0.0), GridInterpolator::nodal_interpolations),
  interpolator_poisson_to_kp(poisson_geometry, kp_geometry, identity_matrix, Vector3D(0.0, 0.0, 0.0), GridInterpolator::nodal_interpolations),
  default_nodal_charges(poisson_equation.get_node_charge_density()),
  kp_current_nodal_charges(1, kp_geometry.get_num_nodes()),
  poisson_current_nodal_charges(1, poisson_geometry.get_num_nodes()),
  current_poisson_solution(0),
  last_poisson_solution(0),
  predictor_corrector_reference_solution(0),
  cb_bands(0),
  vb_bands(0),
  cb_effmass_disp(0),
  vb_effmass_disp(0),
  electron_density(0.0),
  hole_density(0.0),
  temperature(300.0),
  cb_fermi(0.0),
  vb_fermi(0.0),
  number_of_valid_cb_subbands(0),
  number_of_valid_vb_subbands(0),
  num_cb_bands(0),
  num_vb_bands(0),
  cb_distribution(kp_geometry.get_num_nodes()),
  vb_distribution(kp_geometry.get_num_nodes()),
  poisson_solution_difference_norm(0),
  poisson_solution_solution_norm(0),
  sg_update_ratio(Configuration::get_instance()->get("sp_potential_update_start_ratio")),
  min_update_ratio(Configuration::get_instance()->get("sp_potential_min_update_ratio")),
  max_update_ratio(Configuration::get_instance()->get("sp_potential_max_update_ratio")),
  success_increase_ratio(Configuration::get_instance()->get("sp_potential_successfull_update_ratio_increase")),
  failure_decrease_ratio(Configuration::get_instance()->get("sp_potential_bad_update_ratio_decrease")),
  include_lda_xc(Configuration::get_instance()->get("sp_include_lda_xc") == 1.0 ? true:false)
{	
	cb_problem.set_solution_type(electrons);
	vb_problem.set_solution_type(holes);
	cb_problem.prepare();
	vb_problem.prepare();
}

template<class CBPB, class VBPB>
SchroedingerPoisson<CBPB,VBPB>::~SchroedingerPoisson() {

	if(current_poisson_solution != 0) {
		delete current_poisson_solution; current_poisson_solution = 0;
	}
	if(last_poisson_solution != 0) {
		delete last_poisson_solution; last_poisson_solution = 0;
	}
	if(predictor_corrector_reference_solution != 0) {
		delete predictor_corrector_reference_solution; predictor_corrector_reference_solution = 0;
	}
	if(cb_effmass_disp != 0) {
		delete cb_effmass_disp; cb_effmass_disp = 0;
	}
	if(vb_effmass_disp != 0) {
		delete vb_effmass_disp; vb_effmass_disp = 0;
	}

}




/** schroedinger poisson iteration
 *
 * @param electron_density_ low dimensional electron density (e.g. sheet density for QW: 1/nm^2)
 * @param hole_density_ low dimensional electron density (e.g. sheet density for QW: 1/nm^2)
 * @param temperature temperature
 */
template<class CBPB, class VBPB>
bool SchroedingerPoisson<CBPB,VBPB>::calculate(
	unsigned int num_cb_bands_,
	unsigned int num_vb_bands_,
	const double& electron_density_,
	const double& hole_density_,
	const double& temperature_
) {

	const unsigned int max_iterations = 200;
	unsigned int iteration = 0;

	electron_density = electron_density_;
	hole_density = hole_density_;
	temperature = temperature_;
	num_cb_bands = num_cb_bands_;
	num_vb_bands = num_vb_bands_;

	TDKP_ASSERT(num_cb_bands > 0, "");
	TDKP_ASSERT(num_vb_bands > 0, "");
	TDKP_ASSERT(electron_density > 0, "");
	TDKP_ASSERT(hole_density > 0, "");

	// -------------------------------------------
	// precompute potential if there are default
	// charges
	// -------------------------------------------
	if(default_nodal_charges != 0 ||
	   	poisson_equation.get_surface_charge_density() ||
       	poisson_equation.get_element_charge_density()) {
    	this->calculate_potential();
	}

	bool predictor_corrector = Configuration::get_instance()->get("sp_use_predictor_corrector");
	sg_update_ratio = Configuration::get_instance()->get("sp_potential_update_start_ratio");
	poisson_solution_last_difference_norm = -1;
	poisson_solution_difference_norm      = -1;

	// -------------------------------------------
	// loop until convergence
	// -------------------------------------------
	do {
		this->update_kp_potential();
		this->calculate_cb_bandstructure();
		this->calculate_vb_bandstructure();
		this->calculate_fermi_levels();
		this->calculate_carrier_distribution(false);
		this->calculate_potential();
		if(predictor_corrector) {
			this->do_predictor_corrector();
		}
	} while(!this->check_convergence() && ++iteration < max_iterations);

	if(iteration < max_iterations) {
		TDKP_LOGMSG(LOG_INFO, "SchroedingerPoisson: converged after " << iteration << " iterations");
	} else {
		TDKP_LOGMSG(LOG_WARN, "SchroedingerPoisson: still not converged after " << iteration << " iterations. stopping");
	}

	return iteration < max_iterations;

}


template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::calculate_fermi_levels() {

	BandstructureDomain<cplx>* cb_valid = this->get_cb_dispersion().extract_bands(number_of_valid_cb_subbands);
	Fermi cb_fermi_calc(*cb_valid, electrons, 1.0);
	cb_fermi = cb_fermi_calc.calculate_fermi_level(electron_density, temperature);
	delete cb_valid;

	BandstructureDomain<cplx>* vb_valid = this->get_vb_dispersion().extract_bands(number_of_valid_vb_subbands);
	Fermi vb_fermi_calc(*vb_valid, holes, 1.0);
	vb_fermi = vb_fermi_calc.calculate_fermi_level(hole_density, temperature);
	delete vb_valid;

	TDKP_LOGMSG(LOG_INFO_DEVEL2, "SchroedingerPoisson: fermi levels are at: cb = " << cb_fermi << " (" << number_of_valid_cb_subbands << "bands), vb = " << vb_fermi << " (" << number_of_valid_vb_subbands << " bands)");

}

template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::calculate_carrier_distribution(bool predictor_corrector) {

	// --------------------------------------------
	// prepare some constants and required objects
	// --------------------------------------------
	const double kbT = constants::kb * temperature;

	// --------------------------------------------
	// reset carrier distribution to 0
	// --------------------------------------------
	kp_current_nodal_charges.get_data_vector().assign(kp_geometry.get_num_nodes(), 0.0);
	vector<double>& r_vec_charges = kp_current_nodal_charges.get_data_vector();
	cb_distribution.assign(kp_geometry.get_num_nodes(), 0.0);
	vb_distribution.assign(kp_geometry.get_num_nodes(), 0.0);


	TDKP_LOGMSG(LOG_INFO_DEVEL2, "SchroedingerPoisson: valid subbands cb: " << number_of_valid_cb_subbands << ", vb: " <<  number_of_valid_vb_subbands);
	TDKP_ASSERT(number_of_valid_cb_subbands > 0, "");
	TDKP_ASSERT(number_of_valid_vb_subbands > 0, "");

	// --------------------------------------------
	// precompute predictor corrector delta phi
	// --------------------------------------------
	StdNodeData<double> delta_phi_kp_grid;
	if(predictor_corrector) {
		// ----------------------------------------
		// calculate ec(phi - phi_k) (for trellakis, eq. 26)
		// ----------------------------------------
		interpolator_poisson_to_kp.interpolate(*current_poisson_solution, delta_phi_kp_grid);
		for(unsigned int ii = 0; ii < kp_geometry.get_num_nodes(); ii++) {
			delta_phi_kp_grid.set_node_value(ii,0,
				delta_phi_kp_grid.get_node_value(ii,0) + kp_current_nodal_charges.get_node_value(ii,0) / constants::ec
			);
		}
	} else {
		delta_phi_kp_grid.set_length(kp_geometry.get_num_nodes(), 1);
	}

	// --------------------------------------------
	// calculate electron charges of current window
	// --------------------------------------------
	if(this->get_cb_bands().get_basis_size() > 1) {
		// -------------------------------------------------
		// for all cb bands
		// -------------------------------------------------
		for(unsigned int ii = 0; ii < number_of_valid_cb_subbands; ii++) {
			// -------------------------------------------------
			// for all k values
			// -------------------------------------------------
			for(unsigned int kk = 0; kk < get_domain().get_number_of_points(); kk++) {
				// -------------------------------------------------
				// for all real space positions
				// -------------------------------------------------
				EigenSolution<double>* probability = this->get_cb_bands().get_probability(kk,ii);
				for(unsigned int xx = 0; xx < kp_geometry.get_num_nodes(); xx++) {
					double fermi_weight = 1.0 / (1.0 + exp((this->get_cb_bands().get_energy(kk,ii).real()
					                                       - cb_fermi - constants::ec * delta_phi_kp_grid.get_node_value(xx,0)) / kbT))
				                        * get_domain().get_point(kk).get_weight();
					cb_distribution[xx] += fermi_weight * probability->get_node_value(xx,0);
				}
				delete probability;
			}
		}
	} else {
		// -------------------------------------------------
		// for all cb bands
		// -------------------------------------------------
		for(unsigned int ii = 0; ii < number_of_valid_cb_subbands; ii++) {

			// -------------------------------------------------
			// for all real space positions
			// -------------------------------------------------
			EigenSolution<double>* probability = this->get_cb_bands().get_probability(0,ii);
			for(unsigned int xx = 0; xx < kp_geometry.get_num_nodes(); xx++) {
				// -------------------------------------------------
				// calculate occupation integral k f(E(k))
				// -------------------------------------------------
				double fermi_weight = 0.0;
				for(unsigned int kk = 0; kk < get_domain().get_number_of_points(); kk++) {
					fermi_weight += 1.0 / (1.0 + exp((this->get_cb_dispersion().get_energy(kk,ii).real() - cb_fermi - constants::ec * delta_phi_kp_grid.get_node_value(xx,0)) / kbT))
					              * get_domain().get_point(kk).get_weight();
				}
				cb_distribution[xx] += fermi_weight * probability->get_node_value(xx,0);
			}
			delete probability;
		}
	}

	// --------------------------------------------
	// add hole charges, start with non effmass
	// --------------------------------------------
	if(this->get_vb_bands().get_basis_size() > 1) {
		// -------------------------------------------------
		// for all vb bands
		// -------------------------------------------------
		for(unsigned int ii = 0; ii < number_of_valid_vb_subbands; ii++) {
			// -------------------------------------------------
			// for all k values
			// -------------------------------------------------
			for(unsigned int kk = 0; kk < get_domain().get_number_of_points(); kk++) {

				// -------------------------------------------------
				// for all real space positions
				// -------------------------------------------------
				EigenSolution<double>* probability = this->get_vb_bands().get_probability(kk,ii);
				for(unsigned int xx = 0; xx < kp_geometry.get_num_nodes(); xx++) {
					double fermi_weight = 1.0 / (1.0 + exp((vb_fermi - this->get_vb_bands().get_energy(kk,ii).real() + constants::ec * delta_phi_kp_grid.get_node_value(xx,0)) / kbT))
					                    * get_domain().get_point(kk).get_weight();
					vb_distribution[xx] += fermi_weight * probability->get_node_value(xx,0);
				}
				delete probability;
			}
		}
	} else {
		// -------------------------------------------------
		// for all vb bands
		// -------------------------------------------------
		for(unsigned int ii = 0; ii < number_of_valid_vb_subbands; ii++) {
			// -------------------------------------------------
			// for all real space positions
			// -------------------------------------------------
			EigenSolution<double>* probability = this->get_vb_bands().get_probability(0,ii);
			for(unsigned int xx = 0; xx < kp_geometry.get_num_nodes(); xx++) {
				// -------------------------------------------------
				// calculate occupation integral k f(E(k))
				// -------------------------------------------------
				double fermi_weight = 0.0;
				for(unsigned int kk = 0; kk < get_domain().get_number_of_points(); kk++) {
					fermi_weight += 1.0 / (1.0 + exp((vb_fermi - this->get_vb_dispersion().get_energy(kk,ii).real() + constants::ec * delta_phi_kp_grid.get_node_value(xx,0)) / kbT))
					              * get_domain().get_point(kk).get_weight();
				}
				vb_distribution[xx] += fermi_weight * probability->get_node_value(xx,0);
			}
			delete probability;
		}
	}

	// ------------------------------------------------
	// normalize distribution exactly via integration
	// on kp grid
	// ------------------------------------------------
	double cb_norm = 0.0;
	double vb_norm = 0.0;
	for(unsigned int ee = 0; ee < kp_geometry.get_num_elements(); ee++) {
		const Element& elem = kp_geometry.get_element(ee);
		if(elem.enabled()) {
			for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
				cb_norm += elem.get_single_integral_0th_order(nn)
				           * cb_distribution[elem.get_node(nn).get_index_global()];
				vb_norm += elem.get_single_integral_0th_order(nn)
				           * vb_distribution[elem.get_node(nn).get_index_global()];
			}
		}
	}

	// ------------------------------------------------
	// calculate underrelaxed shape of distribution
	// ------------------------------------------------
	for(unsigned int xx = 0; xx < kp_geometry.get_num_nodes(); xx++) {
		// renormalize current distribution to 1
		cb_distribution[xx] /= cb_norm;
		vb_distribution[xx] /= vb_norm;
	}

	// ------------------------------------------------
	// calculate density
	// ------------------------------------------------
	double total_charge = 0.0;
	for(unsigned int xx = 0; xx < kp_geometry.get_num_nodes(); xx++) {
		r_vec_charges[xx]  = (-1.0) * constants::ec * electron_density * cb_distribution[xx];
		r_vec_charges[xx] += ( 1.0) * constants::ec * hole_density * vb_distribution[xx];
	}
	double cb_charge = 0.0;
	for(unsigned int ee = 0; ee < kp_geometry.get_num_elements(); ee++) {
		const Element& elem = kp_geometry.get_element(ee);
		if(elem.enabled()) {
			for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
				total_charge += elem.get_single_integral_0th_order(nn)
			 				  * r_vec_charges[elem.get_node(nn).get_index_global()];
				cb_charge +=  elem.get_single_integral_0th_order(nn)
							* (-1.0) * constants::ec * electron_density * cb_distribution[elem.get_node(nn).get_index_global()];
			}
		}
	}
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "SchroedingerPoisson: total charge is " << total_charge << ", e: " << electron_density << ", h: " << hole_density);

	/*
	ostringstream fname;
	static int counter = 0;
	fname << "charge" << counter++ << ".dat";
	ofstream fout(fname.str().c_str());
	for(unsigned int ii = 0; ii < kp_geometry.get_num_nodes(); ii++) {
		fout << kp_geometry.get_node(ii).get_coord(0) << "  "
		     << r_vec_charges[ii] << " 	"
		     << (-1.0) * constants::ec * electron_density * cb_distribution[ii]
		     << "  "
		     << ( 1.0) * constants::ec * hole_density * vb_distribution[ii]
		     << "  "
		     << delta_phi_kp_grid.get_node_value(ii,0)
		     << "\n";
	}
	fout.close();
	*/
}

template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::do_predictor_corrector() {

	TDKP_LOGMSG(LOG_INFO_DEVEL2, "SchroedingerPoisson: entering PC");
	// --------------------------------------------------
	// prepare for pc iteration
	// --------------------------------------------------
	TDKP_ASSERT(current_poisson_solution != 0, "");
	TDKP_ASSERT(last_poisson_solution != 0, "");

	StdNodeData<double> last_solution_copy(*last_poisson_solution);
	const double cached_norm = poisson_solution_difference_norm; // this should be the norm of our last solution
	if(predictor_corrector_reference_solution == 0) {
		predictor_corrector_reference_solution = new StdNodeData<double>(*current_poisson_solution);
	} else {
		*predictor_corrector_reference_solution = *current_poisson_solution;
	}

	// ---------------------------------------------------
	// predictor corrector has its own update ratio
	// ---------------------------------------------------
	double pc_update_ratio = Configuration::get_instance()->get("sp_potential_update_start_ratio");

	// --------------------------------------------------
	// iterate potential and density
	// --------------------------------------------------
	int iteration = 0;
	int max_iteration = 100;
	do {
		if(iteration > 0) {
			this->potential_update_with_slowdown(pc_update_ratio);
		}
		this->calculate_carrier_distribution(true);
		this->calculate_potential();
	} while(!this->check_convergence(true) && iteration++ < max_iteration);

	// --------------------------------------------------
	// reset last solution (current last solution is the one of the
	// pc iteration)
	// --------------------------------------------------
	*last_poisson_solution = last_solution_copy;
	poisson_solution_difference_norm = cached_norm;

	ostringstream fname;
	static int counter = 0;
	fname << "pc_potential" << counter++ << ".dat";
	ofstream fout(fname.str().c_str());
	for(unsigned int ii = 0; ii < poisson_geometry.get_num_nodes(); ii++) {
		fout << poisson_geometry.get_node(ii).get_coord(0) << " "
		     << predictor_corrector_reference_solution->get_node_value(ii,0) << "  "
		     << current_poisson_solution->get_node_value(ii,0)
		     << "\n";
	}
	fout.close();

	TDKP_LOGMSG(LOG_INFO_DEVEL2, "SchroedingerPoisson: leaving PC");

}

template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::potential_update_with_slowdown(double& update_ratio) {

	// ---------------------------------------------
	// update current potential with new potential
	// ---------------------------------------------
	if(poisson_solution_last_difference_norm != -1) {
		// solution increased difference, reduce update
		if(poisson_solution_difference_norm > poisson_solution_last_difference_norm) {
			update_ratio *= failure_decrease_ratio;
		} else {
			update_ratio *= success_increase_ratio;
		}
		// ensure bounds
		if(update_ratio < min_update_ratio) {
			update_ratio = min_update_ratio;
		}
		if(update_ratio > max_update_ratio) {
			update_ratio = max_update_ratio;
		}
	}
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "SchroedingerPoisson: current update ratio is " << update_ratio << ",\nlast norm = " << poisson_solution_last_difference_norm << ", current norm = " << poisson_solution_difference_norm);

	if(last_poisson_solution != 0) {
		for(int ii = 0; ii < current_poisson_solution->get_length(); ii++) {
			current_poisson_solution->set_node_value(ii,0,
				update_ratio * current_poisson_solution->get_node_value(ii,0)
				+ (1.0 - update_ratio) * last_poisson_solution->get_node_value(ii,0)
			);
		}
	}

}

template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::calculate_potential() {

	// ---------------------------------------------
	// interpolate charge density to poisson grid
	// ---------------------------------------------
	interpolator_kp_to_poisson.interpolate(kp_current_nodal_charges, poisson_current_nodal_charges);

	// ---------------------------------------------
	// add default nodal charges
	// ---------------------------------------------
	if(default_nodal_charges != 0) {
		TDKP_ASSERT(default_nodal_charges->get_length() == poisson_current_nodal_charges.get_length(),"");
		vector<double>& rcharge_vec = poisson_current_nodal_charges.get_data_vector();
		for(int ii = 0; ii < poisson_current_nodal_charges.get_length(); ii++) {
			rcharge_vec[ii] += default_nodal_charges->get_node_value(ii,0);
		}
	}

	// ---------------------------------------------
	// set charges to poisson solver and solve
	// ---------------------------------------------
	poisson_equation.set_node_charge_density(&poisson_current_nodal_charges);
	poisson_equation.solve(0);

	// ---------------------------------------------
	// get solution
	// ---------------------------------------------
	if(last_poisson_solution != 0) {
		delete last_poisson_solution; last_poisson_solution = 0;
	}
	last_poisson_solution = current_poisson_solution;
	current_poisson_solution = poisson_equation.get_solution();

}

/** update current kp potential with new values */
template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::update_kp_potential() {

	// ---------------------------------------------
	// update current potential with new potential
	// ---------------------------------------------
	this->potential_update_with_slowdown(sg_update_ratio);

	// ---------------------------------------------
	// interpolate to new grid
	// ---------------------------------------------
	interpolator_poisson_to_kp.interpolate(*current_poisson_solution, potential_energy_field);

	// ---------------------------------------------
	// multiply potential by electron charge
	// ---------------------------------------------
	for(int ii = 0; ii < potential_energy_field.get_length(); ii++) {
		potential_energy_field.get_node_value(ii,0) *= constants::ec * (-1.0);		
	}

	static int counter = 0;
	ostringstream sout;
	sout << "potential" << counter++ << ".dat";
	ofstream fout(sout.str().c_str());
	for(int ii = 0; ii < potential_energy_field.get_length(); ii++) {
		fout << kp_geometry.get_node(ii).get_coord(0) << "  " << potential_energy_field.get_node_value(ii,0) << "\n";
	}
	fout.close();
}

template<class CBPB, class VBPB>
bool SchroedingerPoisson<CBPB,VBPB>::check_convergence(bool predictor_corrector_mode) {

	if(current_poisson_solution != 0 && last_poisson_solution != 0) {
		poisson_solution_last_difference_norm = poisson_solution_difference_norm;
		poisson_solution_solution_norm   = 0.0;
		poisson_solution_difference_norm = 0.0;
		double tmp = 0.0;
		const double tol = Configuration::get_instance()->get("sp_potential_relative_change_convergence_tolerance");
		for(int ii = 0; ii < current_poisson_solution->get_length(); ii++) {
			tmp = current_poisson_solution->get_node_value(ii,0);
			poisson_solution_solution_norm += tmp * tmp;
			tmp -= last_poisson_solution->get_node_value(ii,0);
			poisson_solution_difference_norm += tmp * tmp;
		}
		poisson_solution_solution_norm = sqrt(poisson_solution_solution_norm);
		poisson_solution_difference_norm = sqrt(poisson_solution_difference_norm);

		ostringstream sout;
		if(predictor_corrector_mode) {
			sout << "SchroedingerPoisson: Predictor corrector |new - old| = "
				 << poisson_solution_difference_norm << "\n|new| = "
				 << poisson_solution_solution_norm << ", tol = " << tol << ": ";
		} else {
			sout << "SchroedingerPoisson: checking field convergence |new - old| = "
				 << poisson_solution_difference_norm << "\n|new| = "
				 << poisson_solution_solution_norm << ", tol = " << tol << ": ";
		}
		bool converged = false;
		if(poisson_solution_solution_norm > 0.0 &&
		   poisson_solution_difference_norm / poisson_solution_solution_norm < tol) {
			converged = true;
		} else if(poisson_solution_difference_norm < tol) {
			converged = true;
		}
		if(converged) {
			sout << " converged!";
			TDKP_LOGMSG(LOG_INFO, sout.str());
		} else {
			sout << " not converged";
			TDKP_LOGMSG(LOG_INFO_DEVEL1, sout.str());
		}
		return converged;
	} else {
		TDKP_ASSERT(!predictor_corrector_mode, "");
		return false;
	}
}


template<class CBPB, class VBPB>
const BandstructureDomain<complex<double> >& SchroedingerPoisson<CBPB,VBPB>::get_cb_bands() const {
	TDKP_ASSERT(cb_bands != 0, "cb bands are not available yet!");
	return *cb_bands;
}

template<class CBPB, class VBPB>
const BandstructureDomain<complex<double> >& SchroedingerPoisson<CBPB,VBPB>::get_vb_bands() const {
	TDKP_ASSERT(vb_bands != 0, "vb bands are not available yet!");
	return *vb_bands;
}

template<class CBPB, class VBPB>
const BandstructureDomain<complex<double> >& SchroedingerPoisson<CBPB,VBPB>::get_cb_dispersion() const {
	if(cb_effmass_disp != 0) {
		return *cb_effmass_disp;
	} else {
		return get_cb_bands();
	}
}

template<class CBPB, class VBPB>
const BandstructureDomain<complex<double> >& SchroedingerPoisson<CBPB,VBPB>::get_vb_dispersion() const {
	if(vb_effmass_disp != 0) {
		return *vb_effmass_disp;
	} else {
		return get_vb_bands();
	}
}



/*
template<>
void SchroedingerPoisson<EffectiveMass,KP4x43D>::calculate_vb_bandstructure();
template<>
void SchroedingerPoisson<EffectiveMass,KP6x63D>::calculate_vb_bandstructure();
template<>
void SchroedingerPoisson<EffectiveMass,KP6x63DWZ>::calculate_vb_bandstructure();
*/

template<>
void SchroedingerPoisson<KP8x83D,KP8x83D>::calculate_vb_bandstructure();
template<>
void SchroedingerPoisson<KP8x81D2D,KP8x81D2D>::calculate_cb_bandstructure();
template<>
void SchroedingerPoisson<KP8x81D2DWZ,KP8x81D2DWZ>::calculate_cb_bandstructure();
template<>
void SchroedingerPoisson<EffectiveMass,EffectiveMass>::calculate_vb_bandstructure();

/** standard implementation for cb bandstructure calculation */
template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::calculate_cb_bandstructure() {

	// --------------------------------------
	// delete old data
	// --------------------------------------
	cb_bands = 0;
	if(cb_effmass_disp != 0) {
		delete cb_effmass_disp; cb_effmass_disp = 0;
	}

	// --------------------------------------
	// calculate bandstructure
	// --------------------------------------
	cb_problem.set_solution_type(electrons);
	cb_problem.set_field(&potential_energy_field);

	const double test_value = -111.0e0;
	double minmax[4];
	for(unsigned int ii = 0; ii < 4; ii++) {
		minmax[ii] = test_value;
	}

	cb_problem.get_minmax_edges(minmax[0], minmax[1], minmax[2], minmax[3]);
	TDKP_ASSERT(minmax[0] != -1 && minmax[1] != -1, "cb edges could not be determined");
	for(unsigned int ii = 0; ii < 4; ii++) {
		TDKP_ASSERT(minmax[ii] != test_value,"minmax[ii] == test_value for " << ii);
	}
	cb_min_edge = minmax[0];
	cb_max_edge = minmax[1];
	cb_problem.set_energy_guess(cb_min_edge + 0.3);
	cb_problem.solve(num_cb_bands);

	// -------------------------------------
	// extract data
	// -------------------------------------
	cb_bands = &cb_problem.get_bandstructure();

	// -------------------------------------
	// build effmass dispersion
	// -------------------------------------
	if(get_domain().get_number_of_points() > 1 && cb_problem.get_geometry().get_dimension() < 3) {
		cb_effmass_disp = create_effmass_dispersion(
			get_cb_transverse_effective_mass(),
			get_domain(),
			*cb_bands
		);
	}
	post_cb_bandstructure();

}

/** standard implementation for vb bandstructure calculation */
template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::calculate_vb_bandstructure() {

	// --------------------------------------
	// delete old data
	// --------------------------------------
	vb_bands = 0;
	vb_problem.delete_solutions();

	// --------------------------------------
	// calculate minmax edges
	// --------------------------------------
	const double test_value = -111.0e0;
	double minmax[4];
	for(unsigned int ii = 0; ii < 4; ii++) {
		minmax[ii] = test_value;
	}

	vb_problem.set_solution_type(holes);
	vb_problem.set_field(&potential_energy_field);
	vb_problem.get_minmax_edges(minmax[0], minmax[1], minmax[2], minmax[3]);

	if(minmax[0] != -1 && minmax[1] != -1) {
		cb_min_edge = minmax[0];
		cb_max_edge = minmax[1];
	}

	for(unsigned int ii = 0; ii < 4; ii++) {
		TDKP_ASSERT(minmax[ii] != test_value,"minmax[ii] == test_value for " << ii);
	}

	TDKP_ASSERT(minmax[2] != -1 && minmax[3] != -1, "");
	vb_min_edge = minmax[2];
	vb_max_edge = minmax[3];

	// --------------------------------------
	// set energy guess and solve
	// --------------------------------------
	vb_problem.set_energy_guess(0, vb_min_edge);
	if(minmax[0] != -1 && minmax[1] != -1) {
		if(cb_min_edge < vb_min_edge) {
			TDKP_LOGMSG(LOG_WARN, "wohaaa, this potential is too tilted and i maybe won't be able to distinguish between cb and vb subbands!");
		}
		vb_problem.set_energy_barrier(0.5 * (cb_min_edge + vb_min_edge));
	}
	vb_problem.solve(num_vb_bands, get_domain());

	// -------------------------------------
	// extract data
	// -------------------------------------
	vb_bands = &vb_problem.get_bandstructure();
	post_vb_bandstructure();

}

template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::post_cb_bandstructure() {
	// -----------------------------------------------------
    // determine number of valid bands due to bandgap energy barrier
    // -----------------------------------------------------
    number_of_valid_cb_subbands = this->get_cb_bands().get_number_of_bands();
    for(int kk = 0; kk < this->get_cb_bands().get_number_of_k_values(); kk++) {
        for(unsigned int bb = 1; bb < number_of_valid_cb_subbands; bb++) {
            if(get_cb_bands().get_energy(kk, bb - 1).real() > get_cb_bands().get_energy(kk, bb).real()) {
                number_of_valid_cb_subbands = bb;
                break;
            }
        }
    }
	TDKP_ASSERT(number_of_valid_cb_subbands > 0, "could not find any valid cb subband!");
}

template<class CBPB, class VBPB>
void SchroedingerPoisson<CBPB,VBPB>::post_vb_bandstructure() {

	// -----------------------------------------------------
    // determine number of valid bands due to bandgap energy barrier
    // -----------------------------------------------------
    number_of_valid_vb_subbands = get_vb_bands().get_number_of_bands();
    for(int kk = 0; kk < get_vb_bands().get_number_of_k_values(); kk++) {
        for(unsigned int bb = 1; bb < number_of_valid_vb_subbands; bb++) {
            if(get_vb_bands().get_energy(kk, bb - 1).real() < get_vb_bands().get_energy(kk, bb).real()) {
                number_of_valid_vb_subbands = bb;
                break;
            }
        }
    }
    TDKP_ASSERT(number_of_valid_vb_subbands > 0, "could not find any valid vb subband!");
}

/** returns the electron effective mass of the lowest bandedge material */
template<class CBPB, class VBPB>
double SchroedingerPoisson<CBPB,VBPB>::get_cb_transverse_effective_mass() const {
	double min_cb_edge = 0.0;
	int    min_region_idx = 0;
	for(unsigned int ii = 0; ii < this->kp_geometry.get_num_regions(); ii++) {
		double cb_edge = this->kp_geometry.get_region(ii).get_material().get("conduction_band_edge");
		if(ii == 0 || cb_edge < min_cb_edge) {
			cb_edge = min_cb_edge;
			min_region_idx = ii;
		}
	}
	return this->kp_geometry.get_region(min_region_idx).get_material().get("electron_effective_mass_transverse");
}
/** returns the electron effective mass of the lowest bandedge material */
template<class CBPB, class VBPB>
double SchroedingerPoisson<CBPB,VBPB>::get_vb_transverse_effective_mass() const {
	double min_vb_edge = 0.0;
	int    min_region_idx = 0;
	for(unsigned int ii = 0; ii < this->kp_geometry.get_num_regions(); ii++) {
		double vb_edge = this->kp_geometry.get_region(ii).get_material().get("valence_band_edge");
		if(ii == 0 || vb_edge > min_vb_edge) {
			vb_edge = min_vb_edge;
			min_region_idx = ii;
		}
	}
	return this->kp_geometry.get_region(min_region_idx).get_material().get("hole_effective_mass_transverse");
}


} // end of namespace

