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

#include "tdkp/clc/CLCOpticalResults.h"

namespace tdkp {

	
CLCOpticalResults::CLCOpticalResults(
	const CLCIngredients& ingredients_,
	const double& omega_min_,
	const double& omega_max_,
	const vector<cplx>& suszeptibility_,
	const vector<double>& spont_emission_,
	const double& cb_fermi_,
	const double& vb_fermi_,
	const double& temperature_	
) : absorption(suszeptibility_.size()),
    spont_emission(spont_emission_),
    delta_n(suszeptibility_.size()),
    suszeptibility(suszeptibility_),
    omega_min(omega_min_),
    omega_max(omega_max_),
    ingredients(ingredients_),
    cb_fermi(cb_fermi_),
    vb_fermi(vb_fermi_),
    temperature(temperature_)
{
	// ---------------------------------------------
	// calculate absorption via 
	// A = k0 Im(chi) k0 = omega * nb / c	
	// ---------------------------------------------	
	for(unsigned int ww = 0; ww < get_omega_num(); ww++) {
		const double omega = omega_min + (omega_max - omega_min) / static_cast<double>(get_omega_num() - 1) * ww;
		const double k0 = omega * ingredients.get_optical_refractive_index() / constants::c; 
		for(unsigned int pp = 0; pp < 3; pp++) {				
			absorption[ww * 3 + pp] = k0 * suszeptibility[ww * 3 + pp].imag();
			delta_n[ww * 3 + pp]    = suszeptibility[ww * 3 + pp].real() / 2.0;   	
		}
	}
	this->evaluate_spont_emission_kms();	
}
	
/** return absorption vector (intensity absorption, amplitude is half of it!) 
 * 
 * vector has length omega_num * 3 and stores values as
 * data[ww * 3 + pp] where ww is omega_idx and pp is polarization idx
 */		
const vector<double>& CLCOpticalResults::get_absorption() const {
	return absorption;	
}
/** return spontaneous emission
 * 
 * vector has length omega_num * 3 and stores values as
 * data[ww * 3 + pp] where ww is omega_idx and pp is polarization idx
 */		 
const vector<double>& CLCOpticalResults::get_spont_emission() const {
	return spont_emission;	
}
/** return spontaneous emission, evaluated from absorption via KMS
 * 
 * vector has length omega_num * 3 and stores values as
 * data[ww * 3 + pp] where ww is omega_idx and pp is polarization idx
 */		 
const vector<double>& CLCOpticalResults::get_spont_emission_kms() const {
	return spont_emission_kms;	
}
/** return refractive index change in fraction (delta_n / n)
 * 
 * vector has length omega_num * 3 and stores values as
 * data[ww * 3 + pp] where ww is omega_idx and pp is polarization idx
 */		
const vector<double>& CLCOpticalResults::get_delta_n() const {
	return delta_n;	
}
/** return optical suszeptibility
 * 
 * vector has length omega_num * 3 and stores values as
 * data[ww * 3 + pp] where ww is omega_idx and pp is polarization idx
 */		
const vector<cplx>&   CLCOpticalResults::get_suszeptibility() const {
	return suszeptibility;	
}
	
const double& CLCOpticalResults::get_omega_min() const {
	return omega_min;	
}
const double& CLCOpticalResults::get_omega_max() const {
	return omega_max;	
}
unsigned int  CLCOpticalResults::get_omega_num() const {
	return absorption.size() / 3;	
}	
	
	
int CLCOpticalResults::get_x_length() const {
	return get_omega_num();	
}
int CLCOpticalResults::get_num_y_sets() const {
	return 11;
}
int CLCOpticalResults::get_num_x_sets() const {
	return 1;
}
void CLCOpticalResults::get_x(int xidx, std::vector<double> &x) const {
	TDKP_ASSERT(xidx == 0, "");
	x.resize(get_omega_num());
	for(unsigned int ww = 0; ww < get_omega_num(); ww++) {
		x[ww] = constants::hbar * (omega_min + (omega_max - omega_min) / static_cast<double>(get_omega_num() - 1) * ww);
	}
}
void CLCOpticalResults::get_y(int yidx, std::vector<double>& y) const {
	y.resize(get_omega_num());
	// absorption
	if(yidx >= 0 && yidx < 3) {
		for(unsigned int ii = 0; ii < get_omega_num(); ii++) {
			y[ii] = absorption[ii * 3 + yidx];	
		}
	// spont emission	
	} else if(yidx >= 3 && yidx < 6) {
		yidx -= 3;
		for(unsigned int ii = 0; ii < get_omega_num(); ii++) {
			y[ii] = spont_emission[ii * 3 + yidx];	
		}		
	// spont emission average		
	} else if(yidx == 6) {
		for(unsigned int ii = 0; ii < get_omega_num(); ii++) {
			y[ii] = ( spont_emission[ii * 3 + 0]
			        + spont_emission[ii * 3 + 1]
			        + spont_emission[ii * 3 + 2] ) / 3.0;
		}
	// delta n
	} else if(yidx < 10) {
		yidx -= 7;
		for(unsigned int ii = 0; ii < get_omega_num(); ii++) {
			y[ii] = delta_n[ii * 3 + yidx];	
		}
	} else if(yidx == 10) {
		// spont emission kms average
		for(unsigned int ii = 0; ii < get_omega_num(); ii++) {
			y[ii] = ( spont_emission_kms[ii * 3 + 0]
			        + spont_emission_kms[ii * 3 + 1]
			        + spont_emission_kms[ii * 3 + 2] ) / 3.0;
		}						
	} else {
		TDKP_GENERAL_EXCEPTION("boooom ... unknown yidx ...");
	}
} 

string CLCOpticalResults::get_x_identifier(int xidx) const {
	return string("omega");
}
string CLCOpticalResults::get_y_identifier(int yidx) const {
	switch(yidx) {
		case 0:
			return string("a_px");
		case 1:
			return string("a_py");			
		case 2:
			return string("a_pz");
		case 3:
			return string("sp_px");
		case 4:
			return string("sp_py");
		case 5:
			return string("sp_pz");
		case 6:
			return string("sp_avg");
		case 7:
			return string("dn_px");
		case 8:
			return string("dn_py");
		case 9:
			return string("dn_pz");
		case 10:
			return string("sp_kms_avg");			
		default:
			TDKP_GENERAL_EXCEPTION("boooom ... unknown yidx ...");				
	}
}	
	
/** calculate spontaneous emission from gain data via Kubo Martin Schwinger */	
void CLCOpticalResults::evaluate_spont_emission_kms() {

	
	const double prefactor = ingredients.get_optical_refractive_index()
	                       * ingredients.get_optical_refractive_index()
	                       / (constants::pi * constants::pi * constants::hbar
	                          * constants::c * constants::c);
	const double deltaFermi = cb_fermi - vb_fermi;
	const double kbT		= constants::kb * temperature;
	
	TDKP_ASSERT(temperature > 0, "");
	 
	double domega = 0.0e0;
	double omega  = get_omega_min();	
	if(get_omega_num() > 1) {
		domega = (get_omega_max() - get_omega_min()) / get_omega_num();
	}
	
	spont_emission_kms.assign(absorption.size(), 0);
	
	for(unsigned int ww = 0; ww < get_omega_num(); ww++) {		
		for(unsigned int pp = 0; pp < 3; pp++) {
			double expfact = (constants::hbar * omega - deltaFermi) / kbT;
			if(expfact != 0.0) {
				spont_emission_kms[ww * 3 + pp] = - prefactor
				                                * absorption[ww * 3 + pp]
				                                * omega * omega
				                                / (1.0 - exp(expfact));  				                                
			}	
		}
		omega += domega;					
	}	
	
}	
	
	
} // end of namespace tdkp
