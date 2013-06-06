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
#include "tdkp/utilities/DensityOfStates.h"

namespace tdkp {
	
// -------------------------------------------------------
// DOS spline curve for fast lookup of x values corresponding
// to certain y values
// -------------------------------------------------------
DensityOfStates::DOSSplineCurve::DOSSplineCurve(BoundaryCondition abc)
: SplineCurve(abc) 
{
}
DensityOfStates::DOSSplineCurve::DOSSplineCurve(const ICurve& curve)
: SplineCurve(curve) 
{
}
DensityOfStates::DOSSplineCurve::DOSSplineCurve(const GridPoints& x, const GridValues& y)
: SplineCurve(x,y) 
{
}

/** find all xi values that match y(xi) = y0
 */
void DensityOfStates::DOSSplineCurve::find_x_values(const double& y0_value, vector<double>& xi_values) const {
	
	TDKP_ASSERT(this->segment_ranges.size() == len - 1, "this->segment_ranges.size() == len");
	// ---------------------------------------
	// loop over all segments
	// ---------------------------------------
	vector<double> segments_xis;
	xi_values.clear();
	for(unsigned int ii = 0; ii < len - 1; ii++) {
		// check if y0 is in the segment's range
		if(this->segment_ranges[ii].min_value <= y0_value && this->segment_ranges[ii].max_value >= y0_value) {
			this->FindInverse(ii, y0_value, segments_xis);
			while(segments_xis.size() > 0) {
				xi_values.push_back(segments_xis.back());
				segments_xis.pop_back();	
			}
		}  	
	} 		
#ifdef DEBUG	
	// --------------------------------------
	// test the inverse
	// --------------------------------------
	const double threshold = 1.0e-6;

	for(unsigned int ii = 0; ii < xi_values.size(); ii++) {
		double yi = this->EvalAt(xi_values[ii]);
		if(tdkp_math::abs(yi - y0_value) > threshold) {
			double ratio = y0_value;
			if(y0_value == 0.0) {
				ratio = 1.0;	
			}			
			ostringstream sout;
			sout << " bad spline interverse calculation: " << ii << " xi: " << xi_values[ii] << " yi: " << yi << " y0_value: " << y0_value;
			sout << " (" << tdkp_math::abs(yi - y0_value) / ratio << ") ";
			Logger::get_instance()->emit(LOG_WARN, sout.str());	
		}	
	}	

#endif
}

void DensityOfStates::DOSSplineCurve::prepare_segments() {
	this->segment_ranges.resize(len - 1);	
	for(int ii = 0; ii < (signed)len - 1; ii++) {		
		this->GetMinMax(ii, this->segment_ranges[ii].min_value, this->segment_ranges[ii].max_value);
	}
}

void DensityOfStates::DOSSplineCurve::get_minmax_all(double& miny, double& maxy) {
	TDKP_ASSERT(this->segment_ranges.size() == len - 1, "this->segment_ranges.size() == len");
	TDKP_ASSERT(len - 1 > 0, "len  - 1 > 0");
	
	miny = this->segment_ranges[0].min_value;
	maxy = this->segment_ranges[0].max_value;
	
	for(unsigned int ii = 1; ii < len - 1; ii++) {
		miny = min(this->segment_ranges[ii].min_value, miny);
		maxy = max(this->segment_ranges[ii].max_value, maxy);			
	}
		
}

// --------------------------------------------------------
// DOS implementation
// --------------------------------------------------------	
	
DensityOfStates::DensityOfStates() 
: dimension(0),
  k_values(0),
  bandstructure(0),
  number_of_subbands(0),
  band_density_of_states(0)  
{
	
}
DensityOfStates::DensityOfStates(unsigned int dim) 
: dimension(dim),
  k_values(0),
  bandstructure(0),
  number_of_subbands(0),
  band_density_of_states(0)
{	
}
DensityOfStates::~DensityOfStates() {
	
}	
	
// --------------------------------------------
// initialization and assignment
// --------------------------------------------
void DensityOfStates::set_dimension(unsigned int dim) {
	dimension = dim;
}	


/** set the bandstructure
 *
 * @param kmin           k0
 * @param kmax           km
 * @param num_k_values  
 * @param number_of_subbands the number of subbands to take from the passed bandstructure 
 * @param bandstructure_ bandstructure, where the dispersion relation is stored 
 *                       as [cb0k0 ... cb0kn, cb1k0 ... cb1kn, .... ] etc.
 */
void DensityOfStates::set_bandstructure(const double& kmin, const double& kmax, unsigned int num_k_values, unsigned int number_of_subbands_, const vector<double>& bandstructure_) {
	TDKP_ASSERT(kmin < kmax, 		"kmin < kmax failed");
	TDKP_ASSERT(kmin >= 0, 			"kmin >= 0 failed");
	TDKP_ASSERT(num_k_values > 1, 	"num_k_values > 1");
	
	tdkp_math::linear_space(k_values, kmin, kmax, num_k_values);
	this->number_of_subbands = number_of_subbands_;
	this->bandstructure.resize(num_k_values * number_of_subbands);
	for(unsigned int ii = 0; ii < this->bandstructure.size(); ii++) {
		this->bandstructure[ii] = bandstructure_[ii];	
	}	
}


void DensityOfStates::set_bandstructure(const double& kmin, const double& kmax, unsigned int num_k_values, unsigned int number_of_subbands_, const vector<cplx>& bandstructure_) {
	TDKP_ASSERT(kmin < kmax, 		"kmin < kmax failed");
	TDKP_ASSERT(kmin >= 0, 			"kmin >= 0 failed");
	TDKP_ASSERT(num_k_values > 1, 	"num_k_values > 1");
	
	tdkp_math::linear_space(k_values, kmin, kmax, num_k_values);
	this->number_of_subbands = number_of_subbands_;
	this->bandstructure.resize(num_k_values * number_of_subbands);
	for(unsigned int ii = 0; ii < this->bandstructure.size(); ii++) {
		this->bandstructure[ii] = bandstructure_[ii].real();	
	}	
}
	
/** init yourself from a bandstructure object 
 * 
 * @param number_of_subbands subbands to take (first [0, num .. [)
 * 
 * */	
void DensityOfStates::set_bandstructure(unsigned int number_of_subbands_, const BandstructureDomain<cplx>& bands) {
	
	TDKP_ASSERT((signed)number_of_subbands_ <= bands.get_number_of_bands(),  "number_of_subbands <= bands.get_number_of_bands()");
	TDKP_ASSERT(bands.get_domain().radial() || dimension == 1, "sorry, but density of states currently does only work for radial approximation to the bandstructure!");

	tdkp_math::linear_space(k_values, bands.get_domain().get_first_point().get_coord_abs(), bands.get_domain().get_last_point().get_coord_abs(), bands.get_number_of_k_values());	
	this->number_of_subbands = number_of_subbands_;
	this->bandstructure.resize(bands.get_number_of_k_values() * number_of_subbands);
#pragma omp parallel for default(shared)	
	for(int bb = 0; bb < (signed)this->number_of_subbands; bb++) {	
		for(int kk = 0; kk < bands.get_number_of_k_values(); kk++) {		
			this->bandstructure[bb * bands.get_number_of_k_values() + kk] = bands.get_energy(kk,bb).real();
		}
	}		
		
}

/** calculate density of states
 */
void DensityOfStates::calculate() {
	ostringstream sout;
	sout << "calculating density of states of a " << dimension << "D carrier gas";
	Logger::get_instance()->emit(LOG_INFO, sout.str());
	sout.str("");
	
	TDKP_ASSERT(this->number_of_subbands > 0, "from a logical point of view, i would need you first to give me some bandstructure before i can calculate the dos!");
	TDKP_ASSERT(this->dimension > 0 && this->dimension < 4, "dimension must be 1, 2 or 3!");

	this->band_density_of_states.resize(0);
	
	// ----------------------------------------------
	// build dense k space where we will search for peaks
	// ----------------------------------------------
	const double& kmin   = k_values[0];
	const double& kmax   = k_values.back();
	TDKP_ASSERT(kmin < kmax, "kmin < kmax");	
	double numk        = (kmax - kmin) / (Configuration::get_instance()->get("dos_peak_search_k_space_spacing"));
	// to next int
	numk = floor(numk);
	if(numk < 64.0) {
		Logger::get_instance()->emit(LOG_WARN, "your dos peak search spacing is too big and gives to few points. will use fixed number of 64 points");
		numk = 64.0;
	}		
	TDKP_ASSERT(numk > 1, "floor(deltak) > 1");	
	vector<double> dense_k_space;
	dense_k_space.resize(int(numk));
	const double deltak = (kmax - kmin) / (numk - 1.0);	
	for(unsigned int ii = 0; ii < dense_k_space.size(); ii++) {
		dense_k_space[ii] = kmin + deltak * double(ii);	
	} 

	// ---------------------------------------------------
	// get involved type of bands
	// --------------------------------------------------- 
	bool have_cb_bands = false;
	bool have_vb_bands = false;	
	for(unsigned int ss = 0; ss < this->number_of_subbands; ss++) {
		if(this->bandstructure[this->k_values.size() * ss] < this->bandstructure[this->k_values.size() * (ss + 1) - 1]) {
			have_cb_bands = true;			
		} else {
			have_vb_bands = true;	
		}
	}
	
	// ----------------------------------------------------
	// for all subbands
	// ----------------------------------------------------
//#pragma omp parallel for default(shared)	
	for(unsigned int ss = 0; ss < this->number_of_subbands; ss++) {
		
		// ------------------------------------------------
		// calculate spline
		// ------------------------------------------------		
		vector<double> subband(this->k_values.size());		
		for(unsigned int aa = 0; aa < this->k_values.size(); aa++) {
			subband[aa] = this->bandstructure[this->k_values.size() * ss + aa];
		}		
		DOSSplineCurve energy_vs_k(this->k_values, subband);
		energy_vs_k.prepare_segments();
							
		// ------------------------------------------------
		// e-space refinement in the problematic sections
		// where |E'(k)| < 1.0e-4
		// there we will sample heavier
		// ------------------------------------------------
		// first, determine problematic region		
		const double       threshold      = Configuration::get_instance()->get("dos_peak_threshold");						
		vector<double>     prblm(dense_k_space.size(), 0);
		unsigned int       problems   = 0;
		vector<double>     derivatives_dense;
		vector<double>     subband_dense;
		
		energy_vs_k.EvalDerivativeAt(dense_k_space, derivatives_dense);
		energy_vs_k.EvalAt(dense_k_space, subband_dense);

		// -------------------------------------------
		// try to find peaky regions
		// -------------------------------------------	
		for(unsigned int ii = 0; ii < derivatives_dense.size(); ii++) {		
			if(derivatives_dense[ii] == 0.0) {
				prblm[ii] = -30; // we treat 0.0 as 1.0e-30 (as inc_points will only go to 1.0e-20)
				problems++;
			} else if(tdkp_math::abs(derivatives_dense[ii]) < threshold) {
				prblm[ii] = floor(tdkp_math::abs(log(tdkp_math::abs(derivatives_dense[ii])) / log(10.0)));	
				problems++;
			}							
		}
		
		// ----------------------------------------------
		// now, build refined dos e space (x-axis)
		// where we refine when we have to expect peaks
		// ----------------------------------------------
		// determine min. discretization
		unsigned int min_points = (unsigned) int(Configuration::get_instance()->get("dos_min_energy_axis_points"));
 		const unsigned int inc_points     = (unsigned) int(Configuration::get_instance()->get("dos_peak_increments"));
		const unsigned int max_inc_points = (unsigned) int(Configuration::get_instance()->get("dos_peak_max_increments"));		
		if(this->k_values.size() > min_points) {
			min_points = this->k_values.size();
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, "your k space is more dense than dos_min_energy_axis_points. so i'll use the size of k space for the dos e space");
		}		
		double start_e;
		double end_e;
		// ------------------------------------------
		// determine energy range of subband (which is the range of the DOS)
		// ------------------------------------------
		energy_vs_k.get_minmax_all(start_e, end_e);
				
		double dE = (end_e - start_e) / double(min_points - 1);
		TDKP_ASSERT(start_e < end_e, "start_e < end_e");
		// build energy space
		vector<double> dos_energy_space;
		// --------------------------------------------
		// how many energy points did we skip? 
		// -> if we have a nasty peak, we do not sample,
		// but we sample dense close to the peak
		// --------------------------------------------
		int skipped_energy_points = 0; 
		
		// --------------------------------------------
		// build refined energy space for dos
		// --------------------------------------------
		if(problems > 0) {
			for(unsigned int ii = 0; ii < min_points; ii++) {			
				double segment_start_energy = start_e + double(ii) * dE;
				double segment_end_energy   = segment_start_energy + dE;
				// ------------------------------------------
				// check if there is a peak in the current segment
				// ------------------------------------------			
				double peak_level = 0.0;
				for(unsigned int kk = 0; kk < dense_k_space.size(); kk++) {
					// peaky region
					if(segment_start_energy <= subband_dense[kk] &&
					   segment_end_energy   >= subband_dense[kk] &&
					   prblm[kk] > 0.0) {
						peak_level = max(peak_level,prblm[kk]);
					}		   										
				}
				// ------------------------------------------
				// if it was peak, we refine
				// ------------------------------------------
				if(peak_level > 0) {
					// ------------------------------------------
					// refinement depends on the power of the peak
					// ------------------------------------------					
					unsigned int additional_points = int(peak_level) * inc_points;
					if(additional_points > max_inc_points) {
						additional_points = max_inc_points;	
					}				
					TDKP_ASSERT(additional_points > 1, "additional_points > 1");
					const double local_energy_spacing = dE / additional_points;
					// -----------------------------------------------
					// insert points but not where we expect a peak!
					// -----------------------------------------------									
					for(unsigned int ee = 0; ee < additional_points; ee++) {												
						double suggested_point = segment_start_energy + double(ee) * local_energy_spacing;
						// accept peak?
						if(tdkp_math::abs(energy_vs_k.EvalDerivativeAt(suggested_point)) < threshold) {
							dos_energy_space.push_back(suggested_point);
						} else {
							skipped_energy_points++;
						}
					}
				} else {
					dos_energy_space.push_back(start_e + double(ii) * dE);	
				}	
			}		
		} else {
			// no problems at all
			for(unsigned int ii = 0; ii < min_points; ii++) {
				dos_energy_space.push_back(start_e + double(ii) * dE);	
			}			
		}		
		
		// -------------------------------------------------------		
		// remove front or end
		// the approximation of the first derivative in the last two
		// points is usually crappy, so we forget it ...
		// -------------------------------------------------------
		if(subband.front() < subband.back()) {
			dos_energy_space.pop_back();			
			dos_energy_space.pop_back();			
		} else {
			dos_energy_space.erase(dos_energy_space.begin());
			dos_energy_space.erase(dos_energy_space.begin());		
		}

		
#ifdef DEBUG
		// test if dos e space is strict monotonic
		for(unsigned int ii = 1; ii < dos_energy_space.size(); ii++) {
		//	cout << "ii: " << ii << "  " << dos_energy_space[ii] << "\n";
			TDKP_ASSERT(dos_energy_space[ii] > dos_energy_space[ii - 1], "dos_energy_space[ii] > dos_energy_space[ii - 1]"); 	
		}
#endif


		

		vector<double> bands_dos(dos_energy_space.size());
		vector<double> kis;
		vector<double> kis_energy_derivatives;
		int stepping = 0;
		double stepping_step = 0.0;
		// and now build the density of states curve
		for(unsigned int ee = 0; ee < dos_energy_space.size(); ee++) {
			// find E(ki) = E0
			energy_vs_k.find_x_values(dos_energy_space[ee], kis);						
			// calculate dE/dk at kis
			if(kis.size() > 0) {				
				energy_vs_k.EvalDerivativeAt(kis, kis_energy_derivatives);
				bool problematic = false;												
				for(unsigned int kk = 0; kk < kis_energy_derivatives.size(); kk++) {
					if(kis_energy_derivatives[kk] == 0.0) {
						kis_energy_derivatives[kk] = 1.0e-20;
					} else {
						kis_energy_derivatives[kk] = tdkp_math::abs(kis_energy_derivatives[kk]);	
					}
					if(kis_energy_derivatives[kk] < threshold) {
						ostringstream eout;												
						eout << "problems with a derivative at E = " << dos_energy_space[ee] << " "
						     << "for kis: " << kk << " (val: " << kis[kk] << ") out of " << kis.size()
						     << ": i got a derivative of " << kis_energy_derivatives[kk] << ". ";
						if(stepping == 0) {
							eout << " trying to resolve by offsetting the energy.";	
						} else {
							eout << " offsetting try: " << stepping << ".";	
						}
						Logger::get_instance()->emit(LOG_INFO_DEVEL2, eout.str());
						problematic = true;
					}					
				}
				
				if(stepping > 10 && problematic) {
					// no further stepping ... just accept
					problematic = false;
					ostringstream eout;
					eout <<  "i have serious problems with peaks in subband " << ss << " which couldn't be resolved by energy offstepping.";
					Logger::get_instance()->emit(LOG_WARN, eout.str());
				} 					
				// ---------------------------------------------
				// so, i will try to offset the given energy point
				// and try another one ...
				// ---------------------------------------------			
				if(problematic) {
					stepping++;
					if(stepping == 1) {
						int step_index = ee;
						if(ee == dos_energy_space.size() - 1) {
							step_index = ee - 1;	
						}
						stepping_step = (dos_energy_space[step_index + 1] - dos_energy_space[step_index]) / 20.0;						
					}
					dos_energy_space[ee] += stepping_step;
					ee--; 
				} else {
					
					stepping = 0;
					
					// --------------------------------------------------------
					// DOS: 
					//   Generell: dos(E) = 1/(delta_k)^n integral (dk)^n delta(E-E(k))
					//   3D: 1/2pi^2 * sum ki^2 / |E'(ki)|  (E(ki) = E)
					//   2D: 1/2pi   * sum ki   / |E'(ki)|  (E(ki) = E)
					//   1D: 1/2pi   * sum  1   / |E'(ki)|  (E(ki) = E)
					// ATTENTION: THIS IS THE SINGLE BAND DOS (SPIN IGNORED!)
					// --------------------------------------------------------
					for(unsigned int kk = 0; kk < kis_energy_derivatives.size(); kk++) {
						switch(dimension) {
							case 3:
								bands_dos[ee] +=  1.0 / (2.0*constants::pi * constants::pi) * kis[kk]*kis[kk] / kis_energy_derivatives[kk];  	
								break;
							case 2:
								bands_dos[ee] +=  1.0 / (2.0*constants::pi) * kis[kk] / kis_energy_derivatives[kk];
								break;
							case 1:
								bands_dos[ee] +=  1.0 / (2.0*constants::pi) / kis_energy_derivatives[kk];	
								break;
							default:
								TDKP_GENERAL_EXCEPTION("invalid dimension");	
						}
					}
				}
			}	
		}

		LinearCurve tmp_dos(dos_energy_space, bands_dos);		        
		band_density_of_states.push_back(tmp_dos);		
		const int width = 30;
		sout.str("");
		sout << "finshed building dos for subband " << ss << ". stats:\n"
		     << setw(width) << "initial k discretization: "  << this->k_values.size() << "\n" 
		     << setw(width) << "peak search space points: "  << dense_k_space.size() << "\n"
		     << setw(width) << "peaky intervals: "           << problems << "\n"		     
		     << setw(width) << "  skipped refinement points: " << skipped_energy_points << "\n"
		     << setw(width) << "subband dos energy points: " << dos_energy_space.size() << "\n"
		     << setw(width) << "dos energy range: "          << dos_energy_space.front() << " -> " << dos_energy_space.back() << "\n"
		     << setw(width) << "max dos density: "           << tmp_dos.Integrate() << "\n";
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	}
	
	// ----------------------------------------------
	// merge subband DOS to total DOS
	// ----------------------------------------------			
	double total_dos_energy_min;
	double total_dos_energy_max;
	// get overall energy range
	this->band_density_of_states[0].GetIntrinsicGridRange(total_dos_energy_min, total_dos_energy_max);	
	for(unsigned int ss = 1; ss < this->band_density_of_states.size(); ss++) {
		double tmp_min, tmp_max;		
		this->band_density_of_states[ss].GetIntrinsicGridRange(tmp_min, tmp_max);
		total_dos_energy_min = min(tmp_min, total_dos_energy_min);
		total_dos_energy_max = max(tmp_max, total_dos_energy_max);
	}		

	// -------------------------------------------
	// shrink energy range to a range where
	// all bands are defined
	// -------------------------------------------	
	for(unsigned int ss = 0; ss < this->band_density_of_states.size(); ss++) {
		bool vb_band = this->bandstructure[this->k_values.size() * ss] > this->bandstructure[this->k_values.size() * (ss + 1) - 1];		
		double emin, emax;
		this->band_density_of_states[ss].GetIntrinsicGridRange(emin, emax);
						
		// valence band				 
		if(vb_band) {
			if(total_dos_energy_min < emin) {
				total_dos_energy_min = emin;
			}
			if(total_dos_energy_max < emax && !have_cb_bands) {
				total_dos_energy_max = emax;	
			}
		} else {
			// conduction band
			if(total_dos_energy_max > emax) {
				total_dos_energy_max = emax;	
			} 
			if(total_dos_energy_min > emin && !have_vb_bands) {				
				total_dos_energy_min = emin;
			}
		}

	}
	
	vector<double> new_total_dos_energy_space;
	vector<double> new_total_dos;
				
	tdkp_math::linear_space(new_total_dos_energy_space, total_dos_energy_min, total_dos_energy_max, int(Configuration::get_instance()->get("dos_min_energy_axis_points")));
	new_total_dos.assign(new_total_dos_energy_space.size(), 0.0);
	// init to zero
	this->total_density_of_states = LinearCurve(new_total_dos_energy_space, new_total_dos);			
		 						
	sout.str("");
	sout << "merging " << this->band_density_of_states.size() << " single subband dos into total dos. "
	     << "total dos space will extend from  " << total_dos_energy_min << " to " << total_dos_energy_max;	  
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	 	
	for(unsigned int ii = 0; ii < this->band_density_of_states.size(); ii++) {				
		// build new dos energy space
		const vector<double>& total_dos_energy_space = this->total_density_of_states.GetGridReference();	
		const vector<double>& band_dos_energy_space  = this->band_density_of_states[ii].GetGridReference();
				
		DOSSplineCurve tmp_band(this->band_density_of_states[ii]);
		
		this->merge_sort_and_unique(total_dos_energy_space, band_dos_energy_space, new_total_dos_energy_space);
		// calculate dos at these points
		this->total_density_of_states.EvalAt(new_total_dos_energy_space, new_total_dos);
		// and merge new values into it
		double emin, emax;		
		this->band_density_of_states[ii].GetIntrinsicGridRange(emin, emax);				
		for(unsigned int jj = 0; jj < new_total_dos_energy_space.size(); jj++) {
			if(emin <= new_total_dos_energy_space[jj] && emax >= new_total_dos_energy_space[jj]) {
				new_total_dos[jj] += this->band_density_of_states[ii].EvalAt(new_total_dos_energy_space[jj]);					
			}
		}
		this->total_density_of_states = LinearCurve(new_total_dos_energy_space, new_total_dos);
	}
	
			
}  

/** merge vector merge_into into vector base without extending the boundaries of vector base */
void DensityOfStates::merge_sort_and_unique(const vector<double>& base, const vector<double>& merge_into, vector<double>& return_target) const {
	vector<double> target(base.begin(), base.end());
	const double max_points = Configuration::get_instance()->get("dos_total_dos_maximum_points");
	const double min_separation = tdkp_math::abs(base.back() - base.front()) / max_points; 
	for(unsigned int ii = 0; ii < merge_into.size(); ii++) {
		if(merge_into[ii] > base.front() && merge_into[ii] < base.back()) {
			target.push_back(merge_into[ii]);	
		}
	}
	 
	vector<double>::iterator eit = unique(target.begin(), target.end());
	sort(target.begin(), eit);
	// build target, but don't get too dense	
	return_target.clear();
	vector<double>::iterator it  = target.begin();
	return_target.push_back((*it));
	it++;
	for(; it != eit; it++) {
		if(tdkp_math::abs(return_target.back() - (*it)) >= min_separation
			|| (it + 1) == eit) {
			return_target.push_back((*it));	
		}
	}	
}

			
// ---------------------------------------------
// result access
// ---------------------------------------------
const ICurve& DensityOfStates::get_dos(unsigned int subband) const {
	TDKP_ASSERT(subband < this->band_density_of_states.size(), "subband < this->band_density_of_states.size()");
	return this->band_density_of_states[subband];
}
const ICurve& DensityOfStates::get_full_dos() const {
	TDKP_ASSERT(this->band_density_of_states.size() > 0, "dos not yet calculated");	
	return this->total_density_of_states;		
}			
			
// ---------------------------------------------
// XY data interface
// ---------------------------------------------
// write DOS at @k using the total density of states grid
// we list <totaldos> <dos1> <dos2> etc. 
// ---------------------------------------------						
int    DensityOfStates::get_x_length()                const {
	return this->total_density_of_states.get_x_length();	
}
int    DensityOfStates::get_num_y_sets()              const {
	return 1 + (signed)this->band_density_of_states.size();	
} 
int    DensityOfStates::get_num_x_sets()              const {
	return 1; // x is energy	
} 
void   DensityOfStates::get_x(int xidx, vector<double> &x) const {
	this->total_density_of_states.get_x(xidx, x);	
}
void   DensityOfStates::get_y(int yidx, vector<double>& y) const {
	if(yidx == 0) {
		this->total_density_of_states.get_y(0, y);	
	} else {
		TDKP_ASSERT(yidx - 1 >= 0 && yidx - 1 < (signed)this->band_density_of_states.size(), "yidx - 1 >= 0 && yidx - 1 < this->band_density_of_states.size()");  
		const vector<double>& tmp_x = this->total_density_of_states.GetGridReference(); 
		y.assign(tmp_x.size(), 0.0);
		// get range
		double emin, emax;		
		this->band_density_of_states[yidx - 1].GetIntrinsicGridRange(emin, emax);
		for(unsigned int ii = 0; ii < tmp_x.size(); ii++) {
			if(tmp_x[ii] >= emin && tmp_x[ii] <= emax) {
				y[ii] =	this->band_density_of_states[yidx - 1].EvalAt(tmp_x[ii]);
			}
		}
	}	
} 
string DensityOfStates::get_x_identifier(int xidx) const {
	return string("k axis");
}
string DensityOfStates::get_y_identifier(int yidx) 	  const {
	if(yidx == 0) {
		return "totaldos";	
	} else {
		ostringstream sout;
		sout << "subband" << yidx;
		return sout.str();
	}		
}	

}
