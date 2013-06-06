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

#include "CLCFermiStats.h"

extern "C" {
	// acml vector function
	void vrda_exp(int, double *, double *);	
}

namespace tdkp
{

void CLCFermiStats::calculate(
	KPSolutionType carrier_type,
	const double& fermi_level,
	const double& temperature,
	const vector<double>& energies,		
	vector<double>& fermi_stats
) const {
	
	TDKP_ASSERT(temperature >= 0.0, "smart ass! temperature must be >= 0 ... "); 
			
	const double kbT       = constants::kb * temperature;
	const double eh_factor = carrier_type == electrons ? 1.0 : -1.0;
		
	fermi_stats.resize(energies.size());
	 			
	// low temperature limit is treated separately
	if(temperature > 0.05) {		
#ifdef NOACML
		for(unsigned int ii = 0; ii < energies.size(); ii++) {
			fermi_stats[ii] = 1.0 / (1.0 + exp(eh_factor * (energies[ii] - fermi_level) / kbT));
		}
#else
		for(unsigned int ii = 0; ii < energies.size(); ii++) {
			fermi_stats[ii] = eh_factor * (energies[ii] - fermi_level) / kbT;	
		}
		vrda_exp(energies.size(), &fermi_stats[0], &fermi_stats[0]);
		for(unsigned int ii = 0; ii < energies.size(); ii++) {
			fermi_stats[ii] = 1.0 / (1.0 + fermi_stats[ii]);
		} 
#endif		
	} else {
		for(unsigned int ii = 0; ii < energies.size(); ii++) {
			if(carrier_type == electrons) { 
				fermi_stats[ii] = energies[ii] > fermi_level ? 0.0 : 1.0;
			} else {
				fermi_stats[ii] = energies[ii] < fermi_level ? 0.0 : 1.0;
			}
		}
	}
}

} // end of namespace
