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


#include <omp.h>
#include "tdkp/utilities/SLC.h"
#include "tdkp/common/Logger.h"
#include "tdkp/clc/Splines.h"

extern "C" {
	// acml vector function
	void vrda_exp(int, double *, double *);
}

namespace tdkp {




// --------------------------------------------------------
// TODO: for nonradial approximation:
//       - fix create dense domain
//       - implement DomainMapping
// --------------------------------------------------------

SLC::SLC(
	const BandstructureDomain<complex<double> >& cb_bands_,
	const BandstructureDomain<complex<double> >& vb_bands_,
	const MatrixElements& matrix_elements_,
	const double& homogeneous_broadening_,
	const double& refractive_index_,
	const double& geometry_output_division_factor_)
: dk(0.004),
  omega_min(0.0),
  omega_max(0.0),
  omega_num(100),
  homogeneous_broadening(homogeneous_broadening_),
  refractive_index(refractive_index_),
  geometry_output_division_factor(geometry_output_division_factor_),
  spont_emission(0),
  absorption(0),
  cb_bands(cb_bands_),
  vb_bands(vb_bands_),
  matrix_elements(matrix_elements_),
  num_cb_bands(0),
  num_vb_bands(0)
{
	// ----------------------------------------------------
	// check bandstructure and matrix elements compatibility
	// ----------------------------------------------------
	TDKP_ASSERT(cb_bands.get_domain().get_dimension() == vb_bands.get_domain().get_dimension(), "band dimensions do not match");
	TDKP_ASSERT(cb_bands.get_domain().get_dimension() == matrix_elements.get_domain().get_dimension(), "matrix element dimensions do not match! " << matrix_elements.get_domain().get_dimension());

	// ----------------------------------------------------
	// check that domain points match
	// ----------------------------------------------------
	TDKP_ASSERT(cb_bands.get_domain().compare_points(vb_bands.get_domain()), "cb_bands.get_domain().compare_points(vb_bands.get_domain())");
	TDKP_ASSERT(cb_bands.get_domain().compare_points(matrix_elements.get_domain()), "cb_bands.get_domain().compare_points(matrix_elements.get_domain())");

	// ----------------------------------------------------
	// check that number of transitions is not higher
	// than the number of bands
	// ----------------------------------------------------
	TDKP_ASSERT((signed)matrix_elements.get_num_cb_bands() <= cb_bands.get_number_of_bands(), "matrix_elements.get_num_cb_bands() (" << matrix_elements.get_num_cb_bands()  << ") <= cb_bands.get_number_of_bands() (" << cb_bands.get_number_of_bands() << ")");
	TDKP_ASSERT((signed)matrix_elements.get_num_vb_bands() <= vb_bands.get_number_of_bands(), "matrix_elements.get_num_vb_bands() <= vb_bands.get_number_of_bands()");

	// ----------------------------------------------------
	// initially we only work with radial approximations
	// ----------------------------------------------------
	if(cb_bands.get_domain().get_dimension() > 1) {
		TDKP_ASSERT(cb_bands.get_domain().radial(), "conduction band is not in radial approximation!");
		TDKP_ASSERT(vb_bands.get_domain().radial(), "valence band is not in radial approximation!");
		TDKP_ASSERT(matrix_elements.get_domain().radial(), "matrix elements are not in radial approximation!");
	}
	this->create_new_dense_domain();
	this->interpolate_values();

	// ----------------------------------------------------
	// set and check photon angular frequency
	// ----------------------------------------------------
	this->set_omega_range_from_bandstructure();
	TDKP_ASSERT(omega_min < omega_max, "omega_min < omega_max");

}

/** create new calculation domain */
void SLC::create_new_dense_domain() {

	// ----------------------------------------------------
	// determine kmax of our integration
	// ----------------------------------------------------
	double radial_kmax = cb_bands.get_domain().get_last_point().get_coord_abs();
	radial_kmax = min(radial_kmax, vb_bands.get_domain().get_last_point().get_coord_abs());
	radial_kmax = min(radial_kmax, matrix_elements.get_domain().get_last_point().get_coord_abs());
    double radial_kmin = cb_bands.get_domain().get_first_point().get_coord_abs();

	TDKP_ASSERT(cb_bands.get_domain().get_dimension() == 0 || radial_kmax > 0.0, "radial_kmax > 0.0");

	// ----------------------------------------------------
	// determine number of k points
	// ----------------------------------------------------
	unsigned int num_k_points = static_cast<unsigned int>(ceil(radial_kmax / dk));

 	ostringstream sout;
 	sout << "SLC: creating dense radial domain with " << num_k_points
 	     << " points within [0, " << radial_kmax << "]";
 	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());


	switch(cb_bands.get_domain().get_dimension()) {
		// --------------------------------------------------
		// quantum dot: just copy the cb grid ...
		// --------------------------------------------------
		case 0:
			dense_domain = cb_bands.get_domain();
			break;
		// --------------------------------------------------
		// quantum wire: create dense trapezoidal line
		// --------------------------------------------------
		case 1:
			create_1D_domain_wire_bands(dense_domain, radial_kmin, radial_kmax, num_k_points);
			break;
		// --------------------------------------------------
		// quantum well: create dense radial integration line
		// --------------------------------------------------
		case 2:
			create_2D_domain_radial(dense_domain, radial_kmin, radial_kmax, num_k_points);
			break;
		// --------------------------------------------------
		// bulk: create dense 3D radial integration line
		// --------------------------------------------------
		case 3:
			create_3D_domain_radial(dense_domain, Vector3D(1.0, 0.0, 0.0), radial_kmin, radial_kmax, num_k_points);
			break;
		default:
			TDKP_GENERAL_EXCEPTION("invalid dimension ... ");
	}


}

void SLC::map_data(const vector<double>& x, const vector<double>& y, const vector<double>& xd, double* yd) const {

	Spline1D sp(x,y,0.0);
	for(unsigned int ii = 0; ii < xd.size(); ii++) {
		yd[ii] = sp(xd[ii]);
	}
}

/** interpolate results onto new domain */
void SLC::interpolate_values() {

	// ------------------------------------------------
	// initialize interpolated arrays
	// ------------------------------------------------
	this->num_cb_bands = matrix_elements.get_num_cb_bands();
	this->num_vb_bands = matrix_elements.get_num_vb_bands();
	interpolated_cb_bands.assign(dense_domain.get_number_of_points() * num_cb_bands, 0.0);
	interpolated_vb_bands.assign(dense_domain.get_number_of_points() * num_vb_bands, 0.0);
	interpolated_matrix_elements.assign(3 * num_cb_bands * num_vb_bands * dense_domain.get_number_of_points(), 0.0);

	// ------------------------------------------------
	// interpolate values
	// ------------------------------------------------
	if(dense_domain.get_dimension() == 0) {
		// -----------------------------------------
		// dot does not need any interpolation, so just copy
		// -----------------------------------------
		TDKP_ASSERT(dense_domain.get_number_of_points() == 1, "");
		TDKP_ASSERT(dense_domain.get_point(0).get_weight() == 1.0, "");

		for(unsigned int cc = 0; cc < num_cb_bands; cc++) {
			interpolated_cb_bands[cc] = cb_bands.get_energy(0, cc).real();
			// for all vb bands (transitions)
			for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
				// for all polarization directions
				for(unsigned short pp = 0; pp < 3; pp++) {
					interpolated_matrix_elements[(cc * num_vb_bands + vv) * 3 + pp]
	  				 += matrix_elements.get_abs_square(cc,vv,pp,0);
				}
			}
		}
		// loop over all vb bands
		for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
			interpolated_vb_bands[vv] = vb_bands.get_energy(0, vv).real();
		}

	} else if(Configuration::get_instance()->get("slc_use_splines_for_interpolation") == 1.0) {
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, "interpolating bands and matrix elements on new dense grid using second order splines");
		// build k vector
		const unsigned int number_of_points = dense_domain.get_number_of_points();
		vector<double> kdense(dense_domain.get_number_of_points());
		vector<double> klean(cb_bands.get_domain().get_number_of_points());
		vector<double> values(cb_bands.get_domain().get_number_of_points());
		for(unsigned int ii = 0; ii < dense_domain.get_number_of_points(); ii++) {
			kdense[ii] = dense_domain.get_point(ii).get_coord_abs();
		}
		for(unsigned int ii = 0; ii < cb_bands.get_domain().get_number_of_points(); ii++) {
			klean[ii] = cb_bands.get_domain().get_point(ii).get_coord_abs();
		}

		// loop over all cb bands (and transitions)
		for(unsigned int cc = 0; cc < num_cb_bands; cc++) {
			// get value
			for(unsigned int ii = 0; ii < cb_bands.get_domain().get_number_of_points(); ii++) {
				values[ii] = cb_bands.get_energy(ii, cc).real();
			}
			map_data(klean,values,kdense, &interpolated_cb_bands[number_of_points * cc]);
			// for all vb bands (transitions)
			for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
				// for all polarization directions
				for(unsigned short pp = 0; pp < 3; pp++) {
					for(unsigned int ii = 0; ii < cb_bands.get_domain().get_number_of_points(); ii++) {
						values[ii] = matrix_elements.get_abs_square(cc,vv,pp,ii);
					}
					map_data(klean,values,kdense, &interpolated_matrix_elements[((cc * num_vb_bands + vv) * 3 + pp) * number_of_points]);
				}
			}
		}
		// loop over all vb bands
		for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
			for(unsigned int ii = 0; ii < vb_bands.get_domain().get_number_of_points(); ii++) {
				values[ii] = vb_bands.get_energy(ii, vv).real();
			}
			map_data(klean,values,kdense, &interpolated_vb_bands[number_of_points * vv]);
		}
	} else {
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, "interpolating bands and matrix elements on new dense grid using linear interpolation");
		const unsigned int number_of_points = dense_domain.get_number_of_points();
		DomainMap cb_map(cb_bands.get_domain(), dense_domain);
		// loop over all points in the new domain
		for(unsigned int ii = 0; ii < number_of_points; ii++) {
			// loop over all mapped points to current new domain point
			for(DomainMap::PointMapIterator it = cb_map.begin(ii); it != cb_map.end(ii); it++) {
				// loop over all cb bands (and transitions)
				for(unsigned int cc = 0; cc < num_cb_bands; cc++) {
					interpolated_cb_bands[number_of_points * cc + ii] += (*it).contribution * cb_bands.get_energy(it->point_idx, cc).real();
					// for all vb bands (transitions)
					for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
						// for all polarization directions
						for(unsigned short pp = 0; pp < 3; pp++) {
							interpolated_matrix_elements[((cc * num_vb_bands + vv) * 3 + pp) * number_of_points + ii]
	  						 += (*it).contribution * matrix_elements.get_abs_square(cc,vv,pp,it->point_idx);
						}
					}
				}
				// loop over all vb bands
				for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
					interpolated_vb_bands[number_of_points * vv + ii] += (*it).contribution * vb_bands.get_energy(it->point_idx, vv).real();
				}
			}
		}
	}
}


void SLC::set_omega_num(unsigned int omega_num_) {
	omega_num = omega_num_;
	spont_emission.resize(0);
	absorption.resize(0);
}

/** reset omega space (purges calculated data) */
void SLC::set_omega_range(const double& omega_min_, const double& omega_max_, unsigned int omega_num_) {
	TDKP_ASSERT(omega_min_ < omega_max_, "omega_min < omega_max");
	omega_min = omega_min_;
	omega_max = omega_max_;
	omega_num = omega_num_;
	spont_emission.resize(0);
	absorption.resize(0);
}


SLC::~SLC() {

}

/** reset bandstructure interpolation density */
void SLC::set_interpolation_density(const double& delta_k) {
	TDKP_ASSERT(delta_k > 0.0, "delta_k > 0.0");
	dk = delta_k;
	dense_domain.clean();
	this->create_new_dense_domain();
	this->interpolate_values();
}

/** determine omega range from provided bandstructure */
void SLC::set_omega_range_from_bandstructure() {

	// -------------------------------------------------
	// find smallest transition
	// -------------------------------------------------

	// roughly 0.3 eV below edge
	double min_trans = cb_bands.get_energy(0,0).real() - vb_bands.get_energy(0,0).real();
	for(unsigned int cc = 0; cc < num_cb_bands; cc++) {
		for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
			min_trans = min(min_trans, cb_bands.get_energy(0,cc).real() - vb_bands.get_energy(0,vv).real());
		}
	}
	omega_min = (min_trans - 0.2) / constants::hbar;

	// plus 2 eV
	omega_max = omega_min + 0.7 / constants::hbar;
	TDKP_ASSERT(omega_min < omega_max, "omega_min < omega_max failed - is the kmax vector of your bandstructure big enough?");

}



void SLC::calculate_fermi_distributions(const double& cb_fermi, const double& vb_fermi, const double& temperature, vector<double>& cb_fermi_dist, vector<double>& vb_fermi_dist ) const {

	TDKP_ASSERT(temperature >= 0.0, "smart ass! temperature must be >= 0 ... ");

	double kbT = constants::kb * temperature;
	TDKP_BOUNDS_ASSERT(cb_fermi_dist.size() == interpolated_cb_bands.size(), "cb_fermi_dist.size() == interpolated_cb_bands.size()");
	TDKP_BOUNDS_ASSERT(vb_fermi_dist.size() == interpolated_vb_bands.size(), "vb_fermi_dist.size() == interpolated_vb_bands.size()");

	// low temperature limit is treated separately
	if(temperature > 0.05) {
#ifdef NOACML
		for(unsigned int ii = 0; ii < interpolated_cb_bands.size(); ii++) {
			cb_fermi_dist[ii] = 1.0 / (1.0 + exp((interpolated_cb_bands[ii] - cb_fermi) / kbT));
		}
#else
		for(unsigned int ii = 0; ii < interpolated_cb_bands.size(); ii++) {
			cb_fermi_dist[ii] = (interpolated_cb_bands[ii] - cb_fermi) / kbT;
		}
		vrda_exp(interpolated_cb_bands.size(), &cb_fermi_dist[0], &cb_fermi_dist[0]);
		for(unsigned int ii = 0; ii < interpolated_cb_bands.size(); ii++) {
			cb_fermi_dist[ii] = 1.0 / (1.0 + cb_fermi_dist[ii]);
		}
#endif
	} else {
		for(unsigned int ii = 0; ii < interpolated_cb_bands.size(); ii++) {
			cb_fermi_dist[ii] = interpolated_cb_bands[ii] > cb_fermi ? 0.0 : 1.0;
		}
	}

	// low temperature limit is treated separately
	if(temperature > 0.05) {
#ifdef NOACML
		for(unsigned int ii = 0; ii < interpolated_vb_bands.size(); ii++) {
			vb_fermi_dist[ii] = 1.0 / (1.0 + exp((vb_fermi - interpolated_vb_bands[ii]) / kbT));
		}
#else
		for(unsigned int ii = 0; ii < interpolated_vb_bands.size(); ii++) {
			vb_fermi_dist[ii] = (vb_fermi - interpolated_vb_bands[ii]) / kbT;
		}
		vrda_exp(interpolated_vb_bands.size(), &vb_fermi_dist[0], &vb_fermi_dist[0]);
		for(unsigned int ii = 0; ii < interpolated_vb_bands.size(); ii++) {
			vb_fermi_dist[ii] = 1.0 / (1.0 + vb_fermi_dist[ii]);
		}
#endif
	} else {
		for(unsigned int ii = 0; ii < interpolated_vb_bands.size(); ii++) {
			vb_fermi_dist[ii] = interpolated_vb_bands[ii] < vb_fermi ? 0.0 : 1.0;
		}
	}
}

/** calculate luminescence / absorption */
void SLC::calculate(const double& cb_fermi, const double& vb_fermi, const double& temperature) {
	calculate(cb_fermi, vb_fermi, temperature, absorption, spont_emission);	
}

/** calculate luminescence / absorption for certain fermi levels on the fly
 *
 * (need that for sebastians parallel bulk evaluations)
 */
void SLC::calculate(const double& cb_fermi, const double& vb_fermi, const double& temperature, vector<double>& local_absorption, vector<double>& local_spont_emission) const {

	ostringstream sout;

	// -----------------------------------------------------
	// no spin degeneracy included here!
	// so, i don't take a factor of 2 for spin up / down for neither band
	// this basically means that either the bandstructure is spin
	// degenerate or the spin summation has been carried out in
	// the matrix element!
    // ------------------------------------------------------
	double prefactor = constants::ec * constants::ec
	                 / (constants::vacuum_permittivity
	                    * refractive_index
	                    * constants::c
	                    * constants::m0
	                    * constants::m0
	                    * constants::hbar);

	// -------------------------------------------------------
	// spont prefactor (see (9.7.18) chuang and replace
	// h by hbar (gets rid of 8pi) and E by hbar omega
	// -------------------------------------------------------
	double spont_prefactor = refractive_index * refractive_index
	                       / (constants::pi * constants::pi * constants::hbar
	                          * constants::c * constants::c);

	// -------------------------------------------------------
	// build dimension dependent prefactor
	// -------------------------------------------------------
	double dim_prefactor = 0.0;
	switch(dense_domain.get_dimension()) {
		case 0:
			// dot: sum still is a sum ...
			dim_prefactor = 1.0 / geometry_output_division_factor;
			break;
		// ----------------------------------------------------
		// -> the radial part is handled inside the DomainMaster,
		//    so this is only the volume stuff for the sum to integral!
		// ----------------------------------------------------
		case 1:
			// wire:
			dim_prefactor = 1.0 / (2.0 * constants::pi * geometry_output_division_factor);
			break;
		case 2:
			// well:
			dim_prefactor = 1.0 / (4.0 * constants::pi * constants::pi * geometry_output_division_factor);
			break;
		case 3:
			// bulk:
			dim_prefactor = 1.0 / (8.0 * constants::pi * constants::pi * constants::pi * geometry_output_division_factor);
			break;
		default:
			TDKP_GENERAL_EXCEPTION("invalid dimension");
	}

	local_absorption.assign(omega_num * 3, 0.0);
	local_spont_emission.assign(omega_num * 3, 0.0);
	sout << "SLC: evaluating spontaneous emission and absorption for "
	     << num_cb_bands << " cb bands and " << num_vb_bands << " vb bands "
	     << "for photons within omega range of "
	     << constants::hbar * omega_min << " -> "
	     << constants::hbar * omega_max << " [eV] and fermi levels of "
	     << cb_fermi << " (cb) and " << vb_fermi << " (vb) at "
	     << "a temperature of " << temperature << " Kelvin.";
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, sout.str());

	// -------------------------------------
	// calculate fermi distributions
	// -------------------------------------
	vector<double> cb_fermi_dist(interpolated_cb_bands.size());
	vector<double> vb_fermi_dist(interpolated_vb_bands.size());
	this->calculate_fermi_distributions(cb_fermi, vb_fermi, temperature, cb_fermi_dist, vb_fermi_dist);

	// -------------------------------------
	// possible restriction of luminescence 
	// evaluation to some bands
	// -------------------------------------
	const int cb_band_start = Configuration::get_instance()->get("slc_restrict_cb_range_start");
	const int cb_band_end   = Configuration::get_instance()->get("slc_restrict_cb_range_end");
	const int vb_band_start = Configuration::get_instance()->get("slc_restrict_vb_range_start");
	const int vb_band_end   = Configuration::get_instance()->get("slc_restrict_vb_range_end");

	// -------------------------------------
	// for every transition
	// -------------------------------------
	unsigned int number_of_points = dense_domain.get_number_of_points();
	unsigned int num_threads      = omp_get_max_threads();

	// -------------------------------------
	// parallel evaluation
	// -------------------------------------
	vector<vector<double> >  thread_absorptions(num_threads);
	vector<vector<double> >  thread_emissions(num_threads);
	for(unsigned int ii = 0; ii < num_threads; ii++) {
		thread_absorptions[ii].assign(omega_num * 3, 0.0);
		thread_emissions[ii].assign(omega_num * 3, 0.0);
	}

	#pragma omp parallel
	{
		vector<double> omega_k(number_of_points);
		vector<double> inversion(number_of_points);
		vector<double> spont_emission_stat(number_of_points); // bahhhh... don't know how to name that: fcb * fvb
		vector<double> work(number_of_points);
		vector<double> lineshape(number_of_points);
		unsigned int   thread_num         = omp_get_thread_num();
		vector<double>& thread_absorption = thread_absorptions[thread_num];
		vector<double>& thread_emission   = thread_emissions[thread_num];
		bool ignore_transition = false;
		#pragma omp for schedule(dynamic) nowait
		for(int vv = 0; vv < (signed)num_vb_bands; vv++) {
			
			// ----------------------------------------------
			// ignore vb transition on request
			// ----------------------------------------------
			ignore_transition = false;
			if(vb_band_start != -1 && vv < vb_band_start) {
				ignore_transition = true;	
			}
			if(vb_band_end != -1 && vv > vb_band_end) {
				ignore_transition = true;				 	
			} 			
			if(!ignore_transition) {
				for(unsigned int cc = 0; cc < num_cb_bands; cc++) {
					// ----------------------------------------------
					// ignore cb transition on request
					// ----------------------------------------------
					if(cb_band_start != -1 && static_cast<int>(cc) < cb_band_start) {
						ignore_transition = true;	
					}
					if(cb_band_end != -1 && static_cast<int>(cc) > cb_band_end) {
						ignore_transition = true;				 	
					} 			
					if(!ignore_transition) {					
						// calculate omega_k (Ecb(k) - Evb(k))
						calculate_omega_k(cc, vv, omega_k);
						// calculate inversion fck + fhk - 1
						calculate_fermi_stats(cc, vv, cb_fermi_dist, vb_fermi_dist, inversion, spont_emission_stat);
						// loop over all photon angular frequencies
						for(int ww = 0; ww < (signed)omega_num; ww++) {
							double omega = omega_min + (omega_max - omega_min) / (omega_num - 1.0)
							                           * double(ww);
							calculate_lineshape(omega, omega_k, lineshape, work);
							// integrate for all k points
							for(unsigned int kk = 0; kk < number_of_points; kk++) {
								const DomainPoint& point = dense_domain.get_point(kk);
								// for all polarization directions
								for(unsigned int pp = 0; pp < 3; pp++) {
									thread_absorption[ww * 3 + pp] +=
									       point.get_weight()
										 * interpolated_matrix_elements[((cc * num_vb_bands + vv) * 3 + pp) * number_of_points + kk]
										 / omega
										 * lineshape[kk]
										 * inversion[kk];
									thread_emission[ww * 3 + pp] +=
									       omega * omega
									     * point.get_weight() // weight contains the integratiion contribution!
										 * interpolated_matrix_elements[((cc * num_vb_bands + vv) * 3 + pp) * number_of_points + kk]
										 / omega
										 * lineshape[kk]
										 * spont_emission_stat[kk];
								}
							}
						}
					} else {
						TDKP_LOGMSG(LOG_INFO, "SLC: ignoring transition between cb " << cc << " and vb " << vv);		
					}
				}
			} else {
				TDKP_LOGMSG(LOG_INFO, "SLC: ignoring ALL transitions involving vb " << vv);	
			}
		}
	}
	// -------------------------------------------
	// master thread: collect
	// -------------------------------------------
	for(unsigned int tt = 0; tt < num_threads; tt++) {
		for(unsigned int ii = 0; ii < local_absorption.size(); ii++) {
			local_absorption[ii]     += thread_absorptions[tt][ii];
			local_spont_emission[ii] += thread_emissions[tt][ii];
		}
	}

	// -------------------------------------------
	// finally, multiply with prefactor
	// absorption gets -1 (i actually calculated gain ...)
	// -------------------------------------------
	for(unsigned int ii = 0; ii < local_absorption.size(); ii++) {
		local_absorption[ii]     *= (-1.0) * dim_prefactor * prefactor;
		local_spont_emission[ii] *= spont_prefactor * dim_prefactor * prefactor;
	}
		
}

/** calculate lineshape function (using cosh here)
 *
 * if you wanna do lorentzian, use 1/gamma Lorentizan
 */
void SLC::calculate_lineshape(const double& omega, const vector<double>& omega_k, vector<double>& lineshape, vector<double>& work) const {
	TDKP_BOUNDS_ASSERT(lineshape.size() == omega_k.size(), "lineshape.size() == omega_k.size()");

	if(Configuration::get_instance()->get("slc_selected_lineshape_function") == 2.0) {
		for(unsigned int ii = 0; ii < omega_k.size(); ii++) {
			lineshape[ii] = homogeneous_broadening / (homogeneous_broadening * homogeneous_broadening
			                                       + (omega - omega_k[ii]) * (omega - omega_k[ii]));
		}
	} else {
#ifdef NOACML
		for(unsigned int ii = 0; ii < omega_k.size(); ii++) {
			lineshape[ii] = 1.0 / (homogeneous_broadening * cosh((omega - omega_k[ii]) / homogeneous_broadening));
		}
#else
		double dx;
		for(unsigned int ii = 0; ii < omega_k.size(); ii++) {
			dx = (omega - omega_k[ii]) / homogeneous_broadening;
			lineshape[ii] = dx;
			work[ii]      = -dx;
		}
		// cosh = (exp(x) + exp(-x)) / 2
		vrda_exp(lineshape.size(), &lineshape[0], &lineshape[0]);
		vrda_exp(work.size(), &work[0], &work[0]);
		for(unsigned int ii = 0; ii < omega_k.size(); ii++) {
			lineshape[ii] = 2.0 / (homogeneous_broadening * (lineshape[ii] + work[ii]));
		}
#endif
	}

}

/** calculate inversion (fck + fhk - 1) and spont. emission factors (fck * fhk) */
void SLC::calculate_fermi_stats(unsigned int cc, unsigned int vv,  const vector<double>& cb_fermi_dist, const vector<double>& vb_fermi_dist, vector<double>& inversion, vector<double>& spont_emission_stat) const {
	unsigned int number_of_points = dense_domain.get_number_of_points();
	TDKP_ASSERT(cb_fermi_dist.size() >= (cc + 1) * number_of_points , "cb_fermi_dist.size() >= (cc + 1) * number_of_points");
	TDKP_ASSERT(vb_fermi_dist.size() >= (vv + 1) * number_of_points , "vb_fermi_dist.size() >= (vv + 1) * number_of_points");
	inversion.resize(number_of_points);
	spont_emission_stat.resize(number_of_points);
	for(unsigned int ii = 0; ii < number_of_points; ii++) {
		inversion[ii] = cb_fermi_dist[number_of_points * cc + ii] + vb_fermi_dist[number_of_points * vv + ii] - 1.0;
		spont_emission_stat[ii] = cb_fermi_dist[number_of_points * cc + ii] * vb_fermi_dist[number_of_points * vv + ii];
	}
}

/** calculate transition frequency (Ecb - Evb) / hbar */
void SLC::calculate_omega_k(unsigned int cc, unsigned int vv, vector<double>& omega_k) const {
	unsigned int number_of_points = dense_domain.get_number_of_points();
	TDKP_ASSERT(interpolated_cb_bands.size() >= (cc + 1) * number_of_points , "interpolated_cb_bands.size() >= (cc + 1) * number_of_points");
	TDKP_ASSERT(interpolated_vb_bands.size() >= (vv + 1) * number_of_points , "interpolated_vb_bands.size() >= (vv + 1) * number_of_points");
	omega_k.resize(number_of_points);
	for(unsigned int ii = 0; ii < number_of_points; ii++) {
		omega_k[ii] = (interpolated_cb_bands[number_of_points * cc + ii] - interpolated_vb_bands[number_of_points * vv + ii])
		            / constants::hbar;
	}
}

void SLC::dump(const char* filename) const {

	ofstream fout(filename);
	if(fout) {
		for(unsigned int ww = 0; ww < omega_num; ww++) {
			fout << setw(20) << constants::hbar * (omega_min + static_cast<double>(ww)
			                                       * (omega_max - omega_min)
			                                       / static_cast<double>(omega_num - 1))
				 << "  ";
			for(unsigned int pp = 0; pp < 3; pp++) {
				fout << setw(20) << absorption[ww * 3 + pp] << "  ";
			}
			for(unsigned int pp = 0; pp < 3; pp++) {
				fout << setw(20) << spont_emission[ww * 3 + pp] << "  ";
			}
			fout << "\n";
		}
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to file " << filename);
	}

}

int SLC::get_x_length() const {
	return this->get_omega_num();
}

int SLC::get_num_y_sets() const {
	return 7;
}
int SLC::get_num_x_sets() const {
	return 1;
}
void SLC::get_x(int xidx, vector<double> &x) const {
	TDKP_ASSERT(this->get_omega_num() > 1, "this->get_omega_num() > 1");
	TDKP_ASSERT(xidx == 0, "xidx == 0");
	double delta_omega = (this->get_omega_max() - this->get_omega_min()) / static_cast<double>(this->get_omega_num() - 1);
	double omega       = this->get_omega_min();
	x.assign(this->get_omega_num(), 0.0);
	for(unsigned int ii = 0; ii < this->get_omega_num(); ii++, omega += delta_omega) {
		x[ii] = omega * constants::hbar; // eV
	}
}
void SLC::get_y(int yidx, vector<double>& y) const {
	TDKP_ASSERT(yidx >= 0 && yidx < 7, "yidx >= 0 && yidx < 7");
	vector<double> tmp;
	y.assign(this->get_omega_num(), 0.0);
	int offset = 0;
	if(yidx == 6) {
		tmp = this->get_spont_emission();
		for(unsigned int ii = 0; ii < this->get_omega_num(); ii++) {
			y[ii] = (tmp[ii * 3 + 0] + tmp[ii * 3 + 1] + tmp[ii * 3 + 2]) / 3.0;
		}
		return;
	}
	if(yidx < 3) {
		tmp    = this->get_absorption();
		offset = yidx;
	} else {
		tmp = this->get_spont_emission();
		offset = yidx - 3;
	}
	for(unsigned int ii = 0; ii < this->get_omega_num(); ii++) {
		y[ii] = tmp[ii * 3 + offset];
	}
}
string SLC::get_x_identifier(int xidx) const {
	return string("hbar_omega");
}
string SLC::get_y_identifier(int yidx) const {
	switch(yidx) {
		case 0: return string("gx");
		case 1: return string("gy");
		case 2: return string("gz");
		case 3: return string("spx");
		case 4: return string("spy");
		case 5: return string("spz");
		case 6: return string("lumi");
		default: TDKP_GENERAL_EXCEPTION("invalid yidx");
	}
}

const double& SLC::get_absorption(unsigned int pp, unsigned int omega_idx) const {
	unsigned int offset = omega_idx * 3 + pp;
	TDKP_ASSERT(offset < absorption.size(), "offset < absorption.size()");
	return absorption[offset];
}
const double& SLC::get_spont_emission(unsigned int pp, unsigned int omega_idx) const {
	unsigned int offset = omega_idx * 3 +pp;
	TDKP_ASSERT(offset < spont_emission.size(), "offset < spont_emission.size()");
	return spont_emission[offset];
}

/** returns average of spont. emission over all directions */
double SLC::get_spont_emission_average(unsigned int omega_idx) const {
	return (get_spont_emission(0,omega_idx)
	      + get_spont_emission(1,omega_idx)
	      + get_spont_emission(2,omega_idx)) / 3.0;
}

/** shift spontaneous emission 
 *
 * to include some mb effects in luminescence spectras in aqua, we decided
 * to help us with a precomputed, density dependent shift.
 * 
 * Lnew(hbar*omega) = Lold(hbar*omega - shift_peak) * ratio_peak
 * 
 * @param shift_peak eV peak shift
 * @param ratio_peak ratio of peak   
 */ 
void SLC::shift_spont_emission(double shift_peak, double ratio_peak) {
	
	// -------------------------------------------------------------
	// calculate omega
	// -------------------------------------------------------------
	vector<double> tmp_omega(omega_num);
	TDKP_ASSERT(omega_num > 1, "");
	double domega = (omega_max - omega_min) / (omega_num - 1);
	for(unsigned int oo = 0; oo < omega_num; oo++) {
		tmp_omega[oo] = omega_min + domega * static_cast<double>(oo);	
	}

	// -------------------------------------------------------------
	// for all polarizations
	// -------------------------------------------------------------
	for(unsigned int pp = 0; pp < 3; pp++) {
		// -------------------------------------------------------------
		// copy luminescence and rescale
		// -------------------------------------------------------------
		vector<double> tmp_sp(omega_num);
		for(unsigned int oo = 0; oo < omega_num; oo++) {
			tmp_sp[oo] = get_spont_emission(pp,oo);				
		}
		// -------------------------------------------------------------
		// calculate luminesence spline
		// -------------------------------------------------------------
		Spline1D sp_spline(tmp_omega, tmp_sp);
		// -------------------------------------------------------------
		// back interpolation
		// -------------------------------------------------------------		
		for(unsigned int oo = 0; oo < omega_num; oo++) {
			if(tmp_omega[oo] - shift_peak / constants::hbar < get_omega_min()) {
				spont_emission[oo * 3 + pp] = ratio_peak * sp_spline(get_omega_min());
			} else if(tmp_omega[oo] - shift_peak / constants::hbar > get_omega_max()) {
				spont_emission[oo * 3 + pp] = ratio_peak * sp_spline(get_omega_max());
			} else {
				spont_emission[oo * 3 + pp] = ratio_peak * sp_spline(tmp_omega[oo] - shift_peak / constants::hbar);
			}
		}		 		
	} 	
}



} // end of namespace
