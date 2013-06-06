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

#ifndef SLC_H_
#define SLC_H_

#include <list>
#include "tdkp/common/all.h"
#include "tdkp/main/Bandstructure.h"
#include "tdkp/utilities/MatrixElements.h"


namespace tdkp {


/** simple luminescence calculator 
 * 
 * simple luminescence calculator using free carrier theory
 * and Kubo Martin Schwinger (or einstein formalism) for the 
 * relation between absorption and luminescence.
 * 
 * its equivalent to eq. (9.7.17) and (9.7.18) in chuang
 * with the difference that it was derived for bulk, well,
 * wire and dot systems (the summation over the k-vector
 * goes over different dimensions) 
 * 
 * ATTENTIONE:
 * the output has to be divided by a factor:
 * (given as geometry_output_division_factor)
 * - in a bulk crystal: = 1
 * - in a quantum well: = well width
 * - in a quantum wire: = wire area
 * - in a quantum dot:  = dot volume 
 * 
 */
class SLC : public XYData<double> {
public:
	SLC(
		const BandstructureDomain<complex<double> >& cb_bands, 
		const BandstructureDomain<complex<double> >& vb_bands,
		const MatrixElements& matrix_elements,
		const double& homogeneous_broadening,		
		const double& refractive_index,
		const double& geometry_output_division_factor		
	);		
	virtual ~SLC();
			
	void calculate(const double& cb_fermi, const double& vb_fermi, const double& temperature, vector<double>& absorption, vector<double>& spont_emission) const;		
	void calculate(const double& cb_fermi, const double& vb_fermi, const double& temperature);
	void dump(const char* filename) const;
	
	// ------------------------------------------------
	// spont. emission manipulation
	// ------------------------------------------------
	void shift_spont_emission(double delta_photon, double ratio_peak);
	
	// ------------------------------------------------
	// result access
	// ------------------------------------------------
	const vector<double>& get_absorption() const { return absorption; }
	const double& get_absorption(unsigned int pp, unsigned int omega_idx) const;	
	const vector<double>& get_spont_emission() const { return spont_emission; }
	const double& get_spont_emission(unsigned int pp, unsigned int omega_idx) const;
	double		  get_spont_emission_average(unsigned int omega_idx) const;
	const double& get_omega_min() const { return omega_min; }
	const double& get_omega_max() const { return omega_max; }
	unsigned int  get_omega_num() const { return omega_num; }
				
	// ------------------------------------------------
	// configuration functions
	// ------------------------------------------------
	void set_omega_range(const double& omega_min, const double& omega_max, unsigned int omega_num);
	void set_omega_num(unsigned int omega_num);
	void set_interpolation_density(const double& delta_k);

	// ------------------------------------------------
	// XYData interface functions
	// ------------------------------------------------	
	int    get_x_length() const;
	int    get_num_y_sets() const;
	int    get_num_x_sets() const; 
	void   get_x(int xidx, vector<double> &x) const;
	void   get_y(int yidx, vector<double>& y) const; 
	string get_x_identifier(int xidx) const;
	string get_y_identifier(int yidx) const;	
	

private:	
	void create_new_dense_domain();
	void interpolate_values();
	void set_omega_range_from_bandstructure();
	void calculate_omega_k(unsigned int cc, unsigned int vv, vector<double>& omega_k) const;
	void calculate_fermi_stats(unsigned int cc, unsigned int vv, const vector<double>& cb_fermi_dist, const vector<double>& vb_fermi_dist, vector<double>& inversion, vector<double>& spont_emission_stat) const;	
	void calculate_lineshape(const double& omega, const vector<double>& omega_k, vector<double>& lineshape, vector<double>& work) const;	
	void calculate_fermi_distributions(const double& cb_fermi, const double& vb_fermi, const double& temperature, vector<double>& cb_fermi_dist, vector<double>& vb_fermi_dist) const;
	void map_data(const vector<double>& x, const vector<double>& y, const vector<double>& xd, double* yd) const;
	
	/** max difference between k points in interpolated structures */
	double dk;
	/** minimum photon wavenumber value (auto determined from bandstructure) */
	double omega_min;
	/** maximum photon wavenumber (auto determined from bandstructure) */
	double omega_max;
	/** preset number of omega values (default = 100) */
	unsigned int omega_num;
	/** homogenous broadening */
	double homogeneous_broadening;
	/** refractive index of host material */
	double refractive_index;
	/** output division factor taking care of the system size where symmetry of crystal is broken */
	double geometry_output_division_factor;		
	/** spon. emission cache */
	vector<double> spont_emission;	
	/** absorption cache */
	vector<double> absorption;
	/** conduction band bands */
	const BandstructureDomain<complex<double> >& cb_bands;
	/** valence band bands */
	const BandstructureDomain<complex<double> >& vb_bands;
	/** matrix elements */	
	const MatrixElements& matrix_elements;
	/** dense integration domain */
	DomainMaster dense_domain;
	/** interpolated cb bands (given at dense domain points) */
	vector<double> interpolated_cb_bands; // E_cc(kk) = [dense_domain.get_number_of_points() * cc + kk]
	unsigned int num_cb_bands;
	/** interpolated vb bands */
	vector<double> interpolated_vb_bands; // E_vv(kk) = [dense_domain.get_number_of_points() * vv + kk]
	unsigned int num_vb_bands;
	/** interpolated matrix elements */	
	vector<double> interpolated_matrix_elements; 
	// [((cc * num_vb_bands) * 3 + pp) * number_of_points + kk]


						
};



}

#endif /*SLC_H_*/
