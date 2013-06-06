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

#include "tdkp/main/Bandstructure.h"

namespace tdkp {

/** returns newly created eigensolution object using real values (for plotting) */
template<>
EigenSolution<double>* Bandstructure<double>::get_eigensolution_real(int kidx, int band_idx) const throw(Exception*) {

	TDKP_ASSERT(kidx >= 0 && kidx < this->num_k_values && band_idx >= 0 && band_idx < this->num_bands, "kidx or band idx out of range");		

	const EigenSolution<double>& orig = this->get_eigensolution(kidx, band_idx); 	
	
	return new EigenSolution<double>(orig); 
	
}

/** returns newly created eigensolution object using real values (for plotting) */
template<>
EigenSolution<double>* Bandstructure<cplx>::get_eigensolution_real(int kidx, int band_idx) const throw(Exception*) {
	
	TDKP_ASSERT(kidx >= 0 && kidx < this->num_k_values && band_idx >= 0 && band_idx < this->num_bands, "kidx or band idx out of range");		

	const EigenSolution<cplx>& orig = this->get_eigensolution(kidx, band_idx); 	
	int length                      = this->get_solution_length();	
	
	vector<double> data(length * (this->basis_size * 2 + 1));
	
	cplx tmp;
	// rewrite data in real/imag components
	for(int ii = 0; ii < length; ii++) {
		double prob = 0.0; 
		for(int jj = 0; jj < this->basis_size; jj++) {
			tmp = orig.get_node_value(ii, jj);
			prob += abs(tmp) * abs(tmp); 
			data[ii * (this->basis_size * 2 + 1) + jj * 2 + 1] = tmp.real(); 
			data[ii * (this->basis_size * 2 + 1) + jj * 2 + 2] = tmp.imag();
		}
		data[ii * (this->basis_size * 2 + 1) + 0] = prob;		
	}
	// create new eigensolution set
	EigenSolution<double>* result = new EigenSolution<double>(orig.get_energy().real(), &data[0], length, this->basis_size * 2 + 1);
	// set identifiers
	ostringstream sident;
	sident << "probability";
	result->set_identifier(0, sident.str());	
	for(int ii = 0; ii < this->basis_size; ii++) {
		sident.str("");
		sident << "band_" << ii << "_real";
		result->set_identifier(ii*2 + 1, sident.str());
		sident.str("");
		sident << "band_" << ii << "_imag";
		result->set_identifier(ii*2 + 2, sident.str());		
	}

	return result;

}


/** extracts stable bands
 * 
 * means, we sort at k = 0, check how many bands are stable and then continue by sorting
 * according to the smallest imaginary absolute value. that should give somehow a reasonable
 * estimate
 */
void extract_stable_bands( 
	const double& maximum_imag_abs, 
	const BandstructureDomain<complex<double> >& source, 
	BandstructureDomain<complex<double> >& target
) {

	// ----------------------------------------------
	// reinit target to source values
	// ----------------------------------------------
	target.reinit(
		source.get_basis_size(),
		source.get_number_of_bands(),
		source.get_solution_length(),
		source.get_domain()
	);
	
	// ----------------------------------------------
	// for every k value, resort according to imaginary value
	// ----------------------------------------------
	vector<double> imag_values;
	vector<int>    indices;
	for(int kk = 0; kk < source.get_number_of_k_values(); kk++) {
		indices.clear();
		imag_values.clear();
		for(int ii =  0; ii < source.get_number_of_bands(); ii++) {
			double val_imag = tdkp_math::abs(source.get_energy(0,ii).imag());
			imag_values.push_back(val_imag);
			indices.push_back(ii);
		}
		// ----------------------------------------------				
		// sort
		// ----------------------------------------------
		tdkp_math::tracked_sort<int,double>(
			indices.begin(),
			indices.end(),
			imag_values.begin(),
			imag_values.end(),
			true
		);	
		// ----------------------------------------------
		// copy values according to order
		// ---------------------------------------------- 
		for(int ii = 0; ii < source.get_number_of_bands(); ii++) {
			target.add_eigensolution(
				kk, ii, new EigenSolution<cplx>(source.get_eigensolution(kk,indices[ii]))
			);				
		}
	}
	
	// ----------------------------------------------------
	// now, check which bands can be considered as stable and extract them 
	// ----------------------------------------------------
	unsigned int stable = 0;
	for(int ii = 0; ii < source.get_number_of_bands(); ii++) {
		if(tdkp_math::abs(target.get_energy(0,ii).imag()) < maximum_imag_abs) {
			stable++;	
		}
	} 	
	BandstructureDomain<cplx>* tmp = target.extract_bands(stable, 0);
	target = *tmp;
	delete tmp;
	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "stable band filtering selected " << target.get_number_of_bands() << " from " << source.get_number_of_bands() << " bands.");
			
}




} // end of namespace
