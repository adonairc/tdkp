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

#include "tdkp/common/Configuration.h"
#include "tdkp/utilities/Fermi.h"

extern "C" {
	// acml vector function
	void vrda_exp(int, double *, double *);
}

//#define NOACML

namespace tdkp {

/** create fermi solver class
 *
 * includes spin degeneracy for effective mass bands (basis size == 1)
 * @param geometry_scale_factor must be the size of the quantized system
 */

Fermi::Fermi(
	const BandstructureDomain<complex<double> >& bands_,
	KPSolutionType type_,
	const double& geometry_scale_factor_
)
: type(type_),
  num_bands(bands_.get_number_of_bands()),
  num_points(bands_.get_domain().get_number_of_points()),
  domain(bands_.get_domain()),
  start_value(0.0),
  geometry_scale_factor(geometry_scale_factor_),
  max_density(0.0)
{
	start_value = bands_.get_energy(0,0).real();

	// -----------------------------------------
	// calculate prefactor
	// -----------------------------------------
	switch(domain.get_dimension()) { // thats the dimension of the k-space
		case 0:
			prefactor = 1.0 / (geometry_scale_factor);
			break;
		case 1:
			prefactor = 1.0 / (geometry_scale_factor * (2.0 * constants::pi));
			break;
		case 2:
			prefactor = 1.0 / (geometry_scale_factor * (2.0 * constants::pi) * (2.0 * constants::pi));
			break;
		case 3:
			prefactor = 1.0 / ((2.0 * constants::pi) * (2.0 * constants::pi) * (2.0 * constants::pi));
			break;
		default:
			TDKP_GENERAL_EXCEPTION("unhandled dimension");
	}

	// -----------------------------------------
	// eff mass bands -> * 2 for spin degeneracy
	// -----------------------------------------
	if(bands_.get_basis_size() == 1) {
		prefactor *= 2.0;
	}

	TDKP_ASSERT(bands_.get_number_of_bands(), "a bandstructure without bands has no density ...");
	if(domain.get_dimension() > 0 && bands_.get_domain().get_number_of_points() < 3) {
		TDKP_GENERAL_EXCEPTION("the supplied bandstructure has less than 3 k-points. thats not suitable to calculate fermi levels");
	}
	// -----------------------------------------
	// copy dispersion
	// -----------------------------------------
	bands.resize(num_points * num_bands);
	for(unsigned int bb = 0; bb < num_bands; bb++) {
		for(unsigned int kk = 0; kk < num_points; kk++) {
			bands[bb * num_points + kk] = bands_.get_energy(kk,bb).real();
		}
	}

	// -------------------------------------------------
	// calculate max density (where its +/- reasonalbe)
	// -------------------------------------------------
    if(domain.get_dimension() == 0) {
       max_density = calculate_density(bands[num_bands - 1], 1.0);
    } else {
	   max_density = calculate_density(bands[num_points - 1], 1.0);
	}
}

/** newton solver to calculate fermi level */
double Fermi::calculate_fermi_level(const double& target_density, const double& temperature) const {

	//double delta_fermi  = 1.0e-12;
	double fermi_level  = start_value;
	if(type == electrons) {
		fermi_level += 0.1;
	} else {
		fermi_level -= 0.1;
	}
	const double density_precision = Configuration::get_instance()->get("fermi_density_solver_tolerance");
	double density      = calculate_density(fermi_level, temperature);
	double last_density = density;
	double update       = 0.3;
	int    loops        = 200;
	bool   spam         = Configuration::get_instance()->get("output_fermi_level_newton") == 1.0;
	//double density_plus, density_minus;
	double d_density_d_fermi;
	//double analytical;
	double fermi_update = 0.0;
	ostringstream sout;
	sout.precision(4);

	if(spam) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "-------- fermi level solver --------\n" << "max density = " << max_density << ", start value = " << start_value << ", first point = " << bands[0]);
	}

	// ---------------------------------------------------
	// emit waring if we are close to the maximum density
	// ---------------------------------------------------
	if(max_density * 0.9 < target_density) {
		TDKP_LOGMSG(LOG_WARN, "the supplied carrier density (" << target_density << ") is close to or above the maximum (" << max_density << ") that fits in the dispersion!");
	}

	do {
		//density_plus      = calculate_density(fermi_level + delta_fermi, temperature);
		//density_minus     = calculate_density(fermi_level - delta_fermi, temperature);
		//d_density_d_fermi = (density_plus - density_minus) / (2.0 * delta_fermi);
		d_density_d_fermi = calculate_density_derivative(fermi_level, temperature);


		if(spam) {
			sout.str("");
			sout << "N(t) = " << setw(5) << target_density << ", "
			     << "N(ef) = " << setw(5) << density << ", "
			     << "ef = " << setw(5) << fermi_level << ", "
			     << "dN/def = " << setw(5) << d_density_d_fermi << ", "
			     << "upd. = " << setw(5) << update << ", ";
			TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());
		}
		TDKP_ASSERT(d_density_d_fermi != 0.0, "d_density_d_fermi != 0.0");
		fermi_update      = update * (target_density - density) / d_density_d_fermi;
		if(abs(fermi_update) > 0.15) {
			fermi_update *= 0.15 / abs(fermi_update);
		}
		fermi_level += fermi_update;
		density = calculate_density(fermi_level, temperature);
		if(fabs(target_density - last_density) < fabs(target_density - density)) {
			update /= 2.0;
		} else {
			update *= 1.2;
		}
		if(update > 1.0) {
			update = 1.0;
		}
		last_density = density;
	} while(fabs(target_density - density) > density_precision && loops-- > 0);

	TDKP_ASSERT(loops > 0, "fermi solver did not converge: target = " << target_density << ", current = " << density << ", fermi_level = " << fermi_level);

	return fermi_level;

}

/** calculate density from fermi level */
double Fermi::calculate_density(const double& fermi_level, const double& temperature) const {

	double kbT = constants::kb * temperature;
	vector<double> fermi_dist(bands.size());

	// low temperature limit is treated separately
	if(temperature > 0.05) {
		if(type == electrons) {
#ifdef NOACML
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = 1.0 / (1.0 + exp((bands[ii] - fermi_level) / kbT));
			}
#else
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = (bands[ii] - fermi_level) / kbT;
			}
			vrda_exp(bands.size(), &fermi_dist[0], &fermi_dist[0]);
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = 1.0 / (1.0 + fermi_dist[ii]);
			}
#endif
		} else {
#ifdef NOACML
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = 1.0 / (1.0 + exp((fermi_level - bands[ii]) / kbT));
			}
#else
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = (fermi_level - bands[ii]) / kbT;
			}
			vrda_exp(bands.size(), &fermi_dist[0], &fermi_dist[0]);
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = 1.0 / (1.0 + fermi_dist[ii]);
			}
#endif
		}
	} else {
		// low temperature limit
		if(type == electrons) {
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = bands[ii] > fermi_level ? 0.0 : 1.0;
			}
		} else {
			for(unsigned int ii = 0; ii < bands.size(); ii++) {
				fermi_dist[ii] = bands[ii] < fermi_level ? 0.0 : 1.0;
			}
		}
	}

	double density = 0.0;
	double weight;
	for(unsigned int kk = 0; kk < num_points; kk++) {
		weight = domain.get_point(kk).get_weight();
		for(unsigned int bb = 0; bb < num_bands; bb++) {
			density += weight * fermi_dist[bb * num_points + kk];
		}
	}
	return density * prefactor;

}

/** put in fermi level in eV in here to get derivative of density */
double Fermi::calculate_density_derivative(const double& fermi_level, const double& temperature) const {

	// ------------------------------------------------------
	// derivative:
	// instead of integrating f, we integrate df/dEf
	// and df/dEf = 1 / (kbT * (2 + t + 1/t) ) and t = exp((E-Ef) / kbT) for electrons and
	//            = - 1 / (kbT * (2 + t + 1/t) ) and t = exp((Ef - E) / kbT) for holes
	// ------------------------------------------------------

	const double kbT = constants::kb * temperature;
	double t;
	vector<double> derivative_fermi_dist(bands.size());

	// low temperature limit is treated separately
	if(type == electrons) {
		for(unsigned int ii = 0; ii < bands.size(); ii++) {
			t = exp((bands[ii] - fermi_level) / kbT);
			derivative_fermi_dist[ii] = 1.0 / (kbT * (2.0 + t + 1.0 / t));
		}
	} else {
		for(unsigned int ii = 0; ii < bands.size(); ii++) {
			t = exp((fermi_level - bands[ii]) / kbT);
			derivative_fermi_dist[ii] = - 1.0 / (kbT * (2.0 + t + 1.0 / t));
		}
	}

	double density = 0.0;
	double weight;
	for(unsigned int kk = 0; kk < num_points; kk++) {
		weight = domain.get_point(kk).get_weight();
		for(unsigned int bb = 0; bb < num_bands; bb++) {
			density += weight * derivative_fermi_dist[bb * num_points + kk];
		}
	}
	return density * prefactor;

}



}
