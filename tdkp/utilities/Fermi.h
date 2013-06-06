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

#ifndef FERMI_H_
#define FERMI_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Exception.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Domain.h"
#include "tdkp/main/Bandstructure.h"

namespace tdkp {

/** newton solver for fermi level
 * 
 * ATTENTIONE:
 * the solver needs a size number (like in the SLC case) 
 * (given as geometry_scale_factor)
 * - in a bulk crystal: = 1
 * - in a quantum well: = well width
 * - in a quantum wire: = wire area
 * - in a quantum dot:  = dot volume
 * 
 * BUT: if you try to do a non-radial calculation and you calculated
 *      e.g. 1/4 of the C4 symmetric bandstructure, you can fix the
 *      missing part here (or by extending the domain point weights
 *      appropriately)
 */
class Fermi {
public:
	Fermi(
		const BandstructureDomain<complex<double> >& bands_, 
		KPSolutionType type,
		const double& geometry_scale_factor_	
	);
	virtual ~Fermi() {}

	/** put density in  1/nm^3 into here to get fermi level in eV */ 
	double calculate_fermi_level(const double& density, const double& temperature) const;
	/** put in fermi leven in eV here to get density in 1/nm^3 */
	double calculate_density(const double& fermi_level, const double& temperature) const;
	/** put in fermi level in eV in here to get derivative of density */
	double calculate_density_derivative(const double& fermi_level, const double& temperature) const;

private:
	KPSolutionType type;
	unsigned int   num_bands;
	unsigned int   num_points;
	DomainMaster   domain;
	vector<double> bands;
	double         start_value;
	const double   geometry_scale_factor;
	double         prefactor;
	double         max_density;
		
	const double& get(unsigned int kk, unsigned int bb) const {
		TDKP_BOUNDS_ASSERT(bb * num_points + kk < bands.size(), "bb * num_points + kk < bands.size()");
		return bands[bb * num_points + kk];    	
	}
	
};

}

#endif /*FERMI_H_*/
