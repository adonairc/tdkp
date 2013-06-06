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

#include "tdkp/clc/CLCScreening.h"
#include "tdkp/clc/CLCFermiStats.h"

namespace tdkp {
	
/** main constructor for lindhard screening object */
CLCScreening::CLCScreening(
	const CLCIngredients1DR& ingredients_,
	const DomainMaster& dense_k_space		
) 
: ingredients(ingredients_),
  k_space(dense_k_space),
  epsilon_spline(0),
  cb_fermi(0),
  vb_fermi(0),
  temperature(0),
  num_angular_integration_points(80),
  dphi(0.0),
  dtheta(0.0),
  qmin(0.0)
{
	
	// ------------------------------------------------
	// ensure we are radial and concistency
	// ------------------------------------------------
	TDKP_ASSERT(k_space.get_dimension() == 1 || k_space.radial(), "");
	TDKP_ASSERT(k_space.get_dimension() == ingredients.get_cb_bands().get_domain().get_dimension(), "");
		
	// ------------------------------------------------
	// deterimine num angular integration points
	// ------------------------------------------------
	switch(k_space.get_dimension()) {
		case 1:
			// no angular integartion
			break;
		case 2:
			num_angular_integration_points = static_cast<unsigned int>(
				Configuration::get_instance()->get("clc_radial_angular_integration_well_num_points")
			);
			break;		
		case 3:
			num_angular_integration_points = static_cast<unsigned int>(
				Configuration::get_instance()->get("clc_radial_angular_integration_bulk_num_points")
			);
			break;		
		default:
			TDKP_GENERAL_EXCEPTION("invalid dimension");		
	}
	TDKP_ASSERT(num_angular_integration_points > 1, "");
	
	// ------------------------------------------------
	// calculate sin phi (0,2pi) and cos phi (0,2pi)
	// and sin theta (0,pi)
	// ------------------------------------------------
	cos_phi.assign(num_angular_integration_points, 0.0);
	sin_theta.assign(num_angular_integration_points / 2, 0.0);	
	dphi   = (constants::pi * 2.0) / cos_phi.size();
	dtheta = (constants::pi) / sin_theta.size();
	double phi   = 0.0;
	double theta = 0.0;
	for(unsigned int ii = 0; ii < cos_phi.size(); ii++) {
		cos_phi[ii]   = cos(phi);
		phi           += dphi;
	}
	for(unsigned int ii = 0; ii < sin_theta.size(); ii++) {
		sin_theta[ii] = sin(theta);				
		theta         += dtheta;			
	}
	
	qmin = ingredients.get_coulomb_matrix_elements().get_q_values().front();

}

/** evaluate screening */
void CLCScreening::update(unsigned int num_cb_bands, unsigned int num_vb_bands, const double& cb_fermi_, const double& vb_fermi_, const double& temperature_) {
	
	cb_fermi    = cb_fermi_;
	vb_fermi    = vb_fermi_;
	temperature = temperature_;
			
	TDKP_ASSERT((signed)num_cb_bands <= ingredients.get_cb_bands().get_number_of_bands(), "");
	TDKP_ASSERT((signed)num_vb_bands <= ingredients.get_vb_bands().get_number_of_bands(), "");				
						
	// ------------------------------------------
	// check if we have some effmass bands
	// if so, double it (to include spin)
	// ------------------------------------------
	bool cb_effmass = ingredients.get_cb_bands().get_basis_size() == 1 ? true:false;
	bool vb_effmass = ingredients.get_vb_bands().get_basis_size() == 1 ? true:false;
	const int kspace_dimension = k_space.get_dimension();
	bool ideal_coulomb_potential = Configuration::get_instance()->get("clc_screening_use_ideal_coulomb_potential") == 1;
	
	if(kspace_dimension > 1 && ideal_coulomb_potential) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLCScreening: calculating lindhard screening with ideal coulomb potential");	
	} else {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLCScreening: calculating lindhard screening");
	}				

	// -------------------------------------------
	// create prefactor
	// -------------------------------------------	
	double cb_prefactor;
	double vb_prefactor;
	// 1. for the potential
	cb_prefactor = constants::ec * constants::ec 
	             / (4.0 * constants::pi * constants::vacuum_permittivity * ingredients.get_static_dielectric_constant());
	// 2. add factor from sum to integral conversion
	if(k_space.get_dimension() == 1) {
		cb_prefactor /= 2.0 * constants::pi;
	} else if(k_space.get_dimension() == 2) { // well
		// first is conversion to integral, 
		cb_prefactor /= (4.0 * constants::pi * constants::pi);   	
	} else if(k_space.get_dimension() == 3) { // bulk
		// first is conversion to integral
		cb_prefactor /= (8.0 * constants::pi * constants::pi * constants::pi);
	} else {
		TDKP_GENERAL_EXCEPTION("invalid dimension " << k_space.get_dimension()  << " of k space for screening"); 
	}
	vb_prefactor = cb_prefactor;
	// 3. for effective mass stuff
	if(cb_effmass) {
		cb_prefactor *= 2;	
	}
	if(vb_effmass) {
		vb_prefactor *= 2;	
	}	                     
	   	                       						
	// -------------------------------------------
	// kspace extension
	// -------------------------------------------
	unsigned int nkmax = static_cast<unsigned int>(Configuration::get_instance()->get("clc_bands_extrapolate_kmax_n_factor"));
	TDKP_ASSERT(nkmax >= 1, "clc_bands_extrapolate_kmax_n_factor >= 1");
		   	                       								   	                       						
	// -------------------------------------------
	// required data and objects 
	// -------------------------------------------
	CLCFermiStats     carrier_stats;
	vector<double>    tmp_energies(k_space.get_number_of_points() + nkmax - 1, 0);
	vector<double>    tmp_kvalues(k_space.get_number_of_points() + nkmax - 1, 0);
	vector<double>    tmp_stats(k_space.get_number_of_points() + nkmax - 1, 0);	
	
	const vector<double>& q_space = ingredients.get_coulomb_matrix_elements().get_q_values();
	vector<double>    tmp_coulomb_pot(q_space.size());
	
	
	// ---------------------------------------------
	// preinit q-space
	// ---------------------------------------------		
	if(k_space.get_dimension() == 1) {
		epsilon_screening.assign(q_space.size(), 1.0);
	} else if(k_space.get_dimension() == 2) { // well
		// epsilon screening is in fact q*epsilon_screening 
		epsilon_screening = q_space;   	
	} else if(k_space.get_dimension() == 3) { // bulk		
		// epsilon screening is in fact q*q*epsilon_screening 
		epsilon_screening.resize(q_space.size());
		for(unsigned int ii = 0; ii < q_space.size(); ii++) {
			epsilon_screening[ii] = q_space[ii] * q_space[ii];	
		}
	} else { 
		TDKP_GENERAL_EXCEPTION("invalid dimension " << k_space.get_dimension()  << " of k space for screening"); 
	}	
	 	
	// --------------------------------------------
	// copy k values into vector (for my spline ...)
	// --------------------------------------------
	for(unsigned int kk = 0; kk < k_space.get_number_of_points(); kk++) {
		tmp_kvalues[kk] = k_space.get_point(kk).get_coord_abs();	
	}
	// extend kspace
	for(unsigned int ii = 0; ii < nkmax - 1; ii++) {
		tmp_kvalues[k_space.get_number_of_points() + ii] = (tmp_kvalues[k_space.get_number_of_points() - 1]) * (ii + 2); 	
	}
	
	// -----------------------------------------
	// ensure that qspace is not too big
	// -----------------------------------------
	if(q_space.back() + k_space.get_last_point().get_coord_abs() > tmp_kvalues.back() + 1.0e-6) {
		ostringstream sout;
		sout << "CLCScreening: sorry, your q-space reaches " << q_space.back() 
		     << " and your k-space is max at " << k_space.get_last_point().get_coord_abs()
		     << ", but at the same time, your k-space extension (via clc_bands_extrapolate_kmax_n_factor) ends at "
		     << tmp_kvalues.back() << ", "
		     << "which is required for the screening evaluation is too small. please, "
		     << "either increase clc_bands_extrapolate_kmax_n_factor or kmax, or lower "
		     << "q_max. my condition is: kmax + qmax <= clc_bands_extrapolate_kmax_n_factor * kmax.";			     
		TDKP_GENERAL_EXCEPTION(sout.str());			        	
	}
	 		 	
	// -------------------------------------------
	// for all cb bands
	// -------------------------------------------
	for(unsigned int ii = 0; ii < num_cb_bands; ii++) {	
		const Spline1D& cb = ingredients.get_cb_band(ii);
		// -----------------------------------------------
		// evaluate bands for all k points in k space
		// -----------------------------------------------
		for(unsigned int kk = 0; kk < tmp_kvalues.size(); kk++) {
			tmp_energies[kk] = cb(tmp_kvalues[kk]);	
		}
		// -----------------------------------------------
		// evaluate carrier statistics
		// -----------------------------------------------
		carrier_stats.calculate(electrons, cb_fermi, temperature, tmp_energies, tmp_stats);		
		// -----------------------------------------------
		// create new spline
		// ----------------------------------------------- 			
		Spline1D fermi_spline(tmp_kvalues, tmp_stats);
		
		// -----------------------------------------------
		// preevaluate coulomb potential
		// -----------------------------------------------
		const CLCRadialCoulombLambda& cmat = ingredients.get_coulomb_matrix_element_cb_cb(ii); 
		for(unsigned int qq = 0; qq < q_space.size(); qq++) {
			if(ideal_coulomb_potential && kspace_dimension > 1) {
				if(kspace_dimension == 2) {
					tmp_coulomb_pot[qq] = 2.0 * constants::pi;
				} else if(kspace_dimension == 3) {		
					tmp_coulomb_pot[qq] = 4.0 * constants::pi;
				} else {
					TDKP_GENERAL_EXCEPTION("invalid kspace dimension");	
				}
			} else {   			
				tmp_coulomb_pot[qq] = cmat.get(q_space[qq]);
			}		
		}
		
		// -----------------------------------------------
		// evaluate for all q values
		// -----------------------------------------------
		#pragma omp parallel for			
		for(int qq = 0; qq < (signed)q_space.size(); qq++) {											
			epsilon_screening[qq] -=  
			 	cb_prefactor
			 	* tmp_coulomb_pot[qq]  
			 	* evaluate_longitudinal_polarization(
						q_space[qq], cb, fermi_spline
				);							
		}
			
	}

		
	// -------------------------------------------
	// for all vb bands
	// -------------------------------------------
	for(unsigned int ii = 0; ii < num_vb_bands; ii++) {
		const Spline1D& vb = ingredients.get_vb_band(ii);
		for(unsigned int kk = 0; kk < tmp_kvalues.size(); kk++) {
			tmp_energies[kk] = vb(tmp_kvalues[kk]);	
		}
		carrier_stats.calculate(holes, vb_fermi, temperature, tmp_energies, tmp_stats);			
		// -----------------------------------------------
		// create new spline
		// ----------------------------------------------- 			
		Spline1D fermi_spline(tmp_kvalues, tmp_stats);
		// -----------------------------------------------
		// preevaluate coulomb potential
		// -----------------------------------------------
		const CLCRadialCoulombLambda& cmat = ingredients.get_coulomb_matrix_element_vb_vb(ii);	
		for(unsigned int qq = 0; qq < q_space.size(); qq++) {					
			if(ideal_coulomb_potential && kspace_dimension > 1) {
				// ideal quantum well
				if(kspace_dimension == 2) {
					tmp_coulomb_pot[qq] = 2.0 * constants::pi;
				// ideal coulomb potential
				} else if(kspace_dimension == 3) {
					tmp_coulomb_pot[qq] = 4.0 * constants::pi;
				} else {
					TDKP_GENERAL_EXCEPTION("invalid kspace dimension");	
				}
			} else {   			
				tmp_coulomb_pot[qq] = cmat.get(q_space[qq]);
			}		
		}
							
		// -----------------------------------------------
		// evaluate for all q values
		// -----------------------------------------------
		#pragma omp parallel for			
		for(int qq = 0; qq < (signed)q_space.size(); qq++) {		
			epsilon_screening[qq] += // note, plus is here due to change in statistics 
			 	vb_prefactor
			 	* tmp_coulomb_pot[qq] 
			 	* evaluate_longitudinal_polarization(
					q_space[qq], vb, fermi_spline
				);							
		}	
	}	

	
	// ---------------------------------------------
	// create spline
	// ---------------------------------------------
	vector<double> inverse_q(q_space.size());
	vector<double> inverse_epsilon(q_space.size());
	for(unsigned int ii = 0; ii < q_space.size(); ii++) {
		inverse_q[ii] = - q_space[q_space.size() - ii - 1];
		inverse_epsilon[ii] = epsilon_screening[q_space.size() - ii - 1];	
	}
	if(epsilon_spline != 0) {
		delete epsilon_spline;	
	}

	epsilon_spline = new Spline1D(inverse_q, inverse_epsilon);
	
}	
  	
/** evaluate longitudinal polarization (note, has to be scaled by some prefactors, as done in the function update(..)) */  	
double CLCScreening::evaluate_longitudinal_polarization(
	const double& qq, const Spline1D& band, const Spline1D& fermi_stats
) const {
	
	// ------------------------------------------------
	// ensure we are radial
	// ------------------------------------------------
	TDKP_ASSERT(k_space.radial(), "");

	// ------------------------------------------------
	// integrate longitudinal polarization
	// ------------------------------------------------
	double res = 0.0;	
	if(k_space.get_dimension() == 1) {
		double up;
		double down;
		double kval;
		double kprime;
		double weight_correction = 0.5;
		// quantum wire is easy, simple integral, -k folded onto +k,
		// factor 2 in radial weight removed
		for(unsigned int kk = 0; kk < k_space.get_number_of_points(); kk++) {
			// first, +k
			kval   = k_space.get_point(kk).get_coord_abs();
			kprime = tdkp_math::abs(kval - qq);
			up = fermi_stats(kprime)
			   - fermi_stats(kval);
			down = band(kprime) - band(kval); 			   
			res += (k_space.get_point(kk).get_weight() * weight_correction) // factor 2.0, removing radial 2 
			     * up / down;
			// then, -k
			kprime = tdkp_math::abs(kval + qq);
			up = fermi_stats(kprime)
			   - fermi_stats(kval);
			down = band(kprime) - band(kval); 			   
			res += (k_space.get_point(kk).get_weight() * weight_correction) // factor 2.0, removing radial 2 
			     * up / down;						     	
		} 	
	} else if(k_space.get_dimension() == 2) {
		double up;
		double down;
		double kval;
		double kprime;
		double fermi_kval;
		double energy_kval;
		double weight_correction = 1.0 / (2.0 * constants::pi); // remove angular integration
		double pweight;
		// ---------------------------------------------------
		// quantum well longitudinal polarization in radial approximation needs
		// angular integration
		// --------------------------------------------------- 		
		for(unsigned int kk = 0; kk < k_space.get_number_of_points(); kk++) {					
			kval        = k_space.get_point(kk).get_coord_abs();			
			pweight     = k_space.get_point(kk).get_weight() * weight_correction;
			fermi_kval  = fermi_stats(kval);
			energy_kval = band(kval);			
			unsigned int ff_start = (qq == kval ? 1:0);
			// angular integration
			for(unsigned int ff = ff_start; ff < cos_phi.size(); ff++) {
				kprime = sqrt(kval*kval + qq*qq - 2.0 * kval * qq * cos_phi[ff]);				
				up   = fermi_stats(kprime) - fermi_kval;
				down = band(kprime) - energy_kval;				
				res += up / down * dphi * pweight;					
			}
		}		
	} else if(k_space.get_dimension() == 3) {
		double up;
		double down;
		double kval;
		double kprime;
		double fermi_kval;
		double energy_kval;
		double weight_correction = 1.0 / (4.0 * constants::pi); // remove angular integration
		double pweight;
		// ---------------------------------------------------
		// bulk longitudinal polarization in radial approximation needs
		// angular integration over phi and theta 
		// --------------------------------------------------- 		
		for(unsigned int kk = 0; kk < k_space.get_number_of_points(); kk++) {			
			kval        = k_space.get_point(kk).get_coord_abs();			
			pweight     = k_space.get_point(kk).get_weight() * weight_correction;
			fermi_kval  = fermi_stats(kval);
			energy_kval = band(kval);			
			// angular integration
			for(unsigned int gg = 0; gg < sin_theta.size(); gg++) {
				// exclude q = 0!
				int start_phi = sin_theta[gg] == 1.0 ? 1 : 0;	
				for(unsigned int ff = start_phi; ff < cos_phi.size(); ff++) {				
					kprime = sqrt(kval*kval + qq*qq - 2.0 * kval * qq * cos_phi[ff]*sin_theta[gg]);
					up   = fermi_stats(kprime) - fermi_kval;
					down = band(kprime) - energy_kval;
					res += up / down * sin_theta[gg] * dtheta * dphi * pweight;
				}					
			}
		}		
	} else {
		TDKP_GENERAL_EXCEPTION("invalid k-space dimension");	
	}	
	return res;
}

int CLCScreening::get_x_length() const {
	return ingredients.get_coulomb_matrix_elements().get_q_values().size();
}
int CLCScreening::get_num_y_sets() const {
	return 2;
}
int CLCScreening::get_num_x_sets() const {
	return 1;
} 
void CLCScreening::get_x(int xidx, vector<double> &x) const {	
	TDKP_BOUNDS_ASSERT(xidx == 0, "");
	x = ingredients.get_coulomb_matrix_elements().get_q_values();
}
void CLCScreening::get_y(int yidx, vector<double>& y) const {
	//TDKP_BOUNDS_ASSERT(yidx == 0, "");
	if(yidx == 0) {			
		y = epsilon_screening;
	}
} 
string CLCScreening::get_x_identifier(int xidx) const {
	return string("q");	
}
string CLCScreening::get_y_identifier(int yidx) const {
	return string("epsilon(q)"); 
}	

} // end of namespace
