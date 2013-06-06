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

#include "tdkp/utilities/SchroedingerPoisson.h"

namespace tdkp {

template<>
void SchroedingerPoisson<KP8x81D2D,KP8x81D2D>::calculate_cb_bandstructure() {
	// --------------------------------------
	// delete old data
	// --------------------------------------
	cb_bands = 0;	
	cb_problem.delete_solutions();

	// --------------------------------------
	// calculate minmax edges
	// --------------------------------------
	double minmax[4];
	cb_problem.set_solution_type(electrons);
	cb_problem.set_field(&potential_energy_field);	
	cb_problem.get_minmax_edges(minmax[0], minmax[1], minmax[2], minmax[3]);
	
	if(minmax[0] != -1 && minmax[1] != -1) {
		cb_min_edge = minmax[0];
		cb_max_edge = minmax[1];	
	}
	TDKP_ASSERT(minmax[2] != -1 && minmax[3] != -1, "");
	vb_min_edge = minmax[2];
	vb_min_edge = minmax[3];

	// --------------------------------------
	// set energy guess and solve
	// --------------------------------------
	cb_problem.set_energy_guess(0, cb_min_edge + 0.3);
	if(minmax[0] != -1 && minmax[1] != -1) {
		if(cb_min_edge < vb_min_edge) {		
			TDKP_LOGMSG(LOG_WARN, "wohaaa, this potential is too tilted and i maybe won't be able to distinguish between cb and vb subbands!");	
		} 
		cb_problem.set_energy_barrier(0.5 * (cb_min_edge + vb_min_edge));		
	}	
	cb_problem.solve(num_cb_bands, get_domain());
	
	// -------------------------------------	
	// extract data
	// -------------------------------------
	cb_bands = &cb_problem.get_bandstructure();
	post_cb_bandstructure();	
		
}

template<>	
void SchroedingerPoisson<KP8x83D,KP8x83D>::calculate_vb_bandstructure() {
	
	// --------------------------------------
	// delete old data
	// --------------------------------------
	vb_bands = 0;
	if(vb_effmass_disp != 0) {
		delete vb_effmass_disp; vb_effmass_disp = 0;	
	}	
	
	// --------------------------------------
	// calculate bandstructure
	// --------------------------------------
	vb_problem.set_solution_type(holes);
	vb_problem.set_field(&potential_energy_field);
	double minmax[4];
	vb_problem.get_minmax_edges(minmax[0], minmax[1], minmax[2], minmax[3]);
	TDKP_ASSERT(minmax[0] != -1 && minmax[1] != -1, "vb edges could not be determined");
	vb_problem.set_energy_guess(vb_min_edge);		
	vb_problem.solve(num_vb_bands);
	
	// -------------------------------------	
	// extract data
	// -------------------------------------
	vb_bands = &vb_problem.get_bandstructure();	
	post_vb_bandstructure();	
		
}
	
/*	
template<>	
void SchroedingerPoisson<KP8x83D,KP8x83D>::calculate_cb_bandstructure() {	
}

template<>	
void SchroedingerPoisson<EffectiveMass,KP4x43D>::calculate_vb_bandstructure() {
	TDKP_GENERAL_EXCEPTION("not implemented yet");	
}
template<>	
void SchroedingerPoisson<EffectiveMass,KP6x63D>::calculate_vb_bandstructure() {
	TDKP_GENERAL_EXCEPTION("not implemented yet");
}
template<>	
void SchroedingerPoisson<EffectiveMass,KP6x63DWZ>::calculate_vb_bandstructure() {
	TDKP_GENERAL_EXCEPTION("not implemented yet");
}
*/
template<>
void SchroedingerPoisson<KP8x81D2DWZ,KP8x81D2DWZ>::calculate_cb_bandstructure() {
	TDKP_GENERAL_EXCEPTION("not implemented yet");
}
template<>
void SchroedingerPoisson<EffectiveMass,EffectiveMass>::calculate_vb_bandstructure() {
	TDKP_GENERAL_EXCEPTION("not implemented yet");
}

} // end of namespace
