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

#include "tdkp/utilities/SLCMBShift.h"

namespace tdkp
{

void SLCMBShift::update_bandstructure(
	const BandstructureDomain<complex<double> >& cb_bands, 
	const BandstructureDomain<complex<double> >& vb_bands
) {
	if(cb_fermi_calc != 0) {
		delete cb_fermi_calc; cb_fermi_calc = 0;	
	}
	if(vb_fermi_calc != 0) {
		delete vb_fermi_calc; vb_fermi_calc = 0;	
	}
	cb_fermi_calc = new Fermi(cb_bands, electrons, 1.0);
	vb_fermi_calc = new Fermi(vb_bands, holes, 1.0);	
}
SLCMBShift::SLCMBShift(
	const char* filename
) : cb_fermi_calc(0),
	vb_fermi_calc(0)
{	
	// -------------------------------------------
	// read file with mb shifts
	// -------------------------------------------
	ifstream fin(filename);
	
	if(!fin) { 
		TDKP_GENERAL_EXCEPTION("file " << filename << " can not be accessed");	
	}
	
	// -------------------------------------------
	// count number of values
	// (need Npoints = sqrt(num_values / 4)
	// -------------------------------------------
	int counter = 0;
	double tmp;
	vector<double> values;
	string line;
	while(getline(fin,line)) {
		if(line.size() > 0) {
			istringstream sin(line);
			for(int jj = 0; jj < 4; jj++) {
				sin >> tmp;					
				values.push_back(tmp);						
				counter++;					
			}	
		}		
	}
	TDKP_ASSERT(counter % 4 == 0, "number of values (" << counter << ") in file can not be divded by four ... (but i expect that)");
	TDKP_ASSERT(sqrt(counter / 4) * sqrt(counter / 4) == counter / 4, "sqrt(counter / 4)^2 == counter / 4 failed");	
	fin.close();
	
	// -------------------------------------------
	// read lines:
	// ndens pdens shift_peak ratio_peak
	// where stuff is ordered via ndens pdens
	// -------------------------------------------
	unsigned int num_points = sqrt(counter / 4);
	TDKP_ASSERT(num_points > 0, "");
	// first, read ndens array
	dens_points.assign(num_points, 0);
	for(unsigned int ii = 0; ii < num_points; ii++) {
		TDKP_BOUNDS_ASSERT(values.size() > ii * 4 * num_points, "");
		dens_points[ii] = log10(values[ii * 4 * num_points + 0]);
		if(ii > 0) {
			TDKP_BOUNDS_ASSERT(dens_points[ii] > dens_points[ii - 1], "");	
		}				
	}
	// then read for every ndens point all pdens and shift values und create splines
	peak_shifts.resize(num_points);
	peak_ratios.resize(num_points);
	vector<double> tmp_shifts(num_points), tmp_ratios(num_points);
/*
	ostringstream fname;
	static int cnt = 0;
	fname << "mbdebug" << cnt++ << ".dat";	
	ofstream fdebug(fname.str().c_str());
*/			
	for(unsigned int ii = 0; ii < num_points; ii++) {		  
		for(unsigned int jj = 0; jj < num_points; jj++) {
			tmp_shifts[jj] = values[ii * 4 * num_points + jj * 4 + 2];
			tmp_ratios[jj] = values[ii * 4 * num_points + jj * 4 + 3];
/*			fdebug << pow(10.0, dens_points[ii]) << "  " << pow(10.0, dens_points[jj]) << "  "
			       << tmp_shifts[jj] << "   " << tmp_ratios[jj] << "\n"; */ 				
		}	
		peak_shifts[ii] = new Spline1D(dens_points, tmp_shifts);
		peak_ratios[ii] = new Spline1D(dens_points, tmp_ratios);
	}
	//fdebug.close();
		
}

SLCMBShift::~SLCMBShift() {
	for(unsigned int ii = 0; ii < peak_shifts.size(); ii++) {
		if(peak_shifts[ii] != 0) {
			delete peak_shifts[ii]; peak_shifts[ii] = 0;	
		}
		if(peak_ratios[ii] != 0) {
			delete peak_ratios[ii]; peak_ratios[ii] = 0;	
		}				
	}	
	if(cb_fermi_calc != 0) {
		delete cb_fermi_calc; cb_fermi_calc = 0;	
	}
}
	

void SLCMBShift::calculate_mb_shift_density(
	double& shift_peak, 
	double& ratio_peak, 
	const double& ndens,
	const double& pdens
) {

	// calculate density logarithms
	double ndenslog = log10(ndens);
	double pdenslog = log10(pdens);
		
	// determine intervals
	unsigned int high = 0, low = 0;	
	if(dens_points.front() > ndenslog) {
		high = low = 0;	
	} else if(dens_points.back() < ndenslog) {
		high = low = dens_points.size() - 1;	
	} else {
		for(unsigned int ii = 1; ii < dens_points.size(); ii++) {
			if(dens_points[ii] > ndenslog) {
				low  = ii - 1;
				high = ii; 	
				break;
			}
		}
	}
	// determine weight
	double t = 0;
	if(high == low) {
		t = 0.5;  	
		TDKP_LOGMSG(LOG_INFO, "SLCMBShift: ndens = " << ndens << " out of bounds " << pow(10.0, dens_points.front()) << "/" << pow(10.0, dens_points.back()) << ", setting to min (or max)");
	} else {
		t = (ndenslog - dens_points[low]) / (dens_points[high] - dens_points[low]);
		TDKP_BOUNDS_ASSERT(t >= 0.0 && t <= 1.0, "t >= 0.0 && t <= 1.0 failed with t = " << t << " for " << ndenslog << " min/max = " << dens_points[low] << "/" << dens_points[high]);	
	}
	// restrict pdenslog
	if(pdenslog < dens_points.front()) {
		pdenslog = dens_points.front();	
	} else if(pdenslog > dens_points.back()) {
		pdenslog = dens_points.back();	
	}  
	// evaluate	
	shift_peak = t * (*peak_shifts[low])(pdenslog) + (1.0 - t) * (*peak_shifts[high])(pdenslog);
	ratio_peak = t * (*peak_ratios[low])(pdenslog) + (1.0 - t) * (*peak_ratios[high])(pdenslog);
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "SLCMBShift: shift_peak = " << shift_peak << ", t = " << t << ", low/high = " << low << "/" << high << ", ndens/pdens " << ndens << "/" << pdens);  
	
}

void SLCMBShift::calculate_mb_shift_fermi(
	double& shift_peak, 
	double& ratio_peak, 
	const double& cb_fermi,
	const double& vb_fermi,
	const double& temperature
) { 
	// -----------------------------------------
	// check if fermi density calculators are available
	// -----------------------------------------
	TDKP_ASSERT(cb_fermi_calc != 0, "you are trying to calculate the density shift via passing fermi levels but i miss the bandstructure for this operation. please update the bandstructure first before using this function.");
	
	// -----------------------------------------
	// determine corresponding density
	// -----------------------------------------
	double cb_dens = cb_fermi_calc->calculate_density(cb_fermi, temperature);
	double vb_dens = vb_fermi_calc->calculate_density(vb_fermi, temperature);
	this->calculate_mb_shift_density(shift_peak, ratio_peak, cb_dens, vb_dens);	
}	


}
