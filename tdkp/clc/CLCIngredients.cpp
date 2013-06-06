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

#include "CLCIngredients.h"

namespace tdkp
{

CLCIngredients::CLCIngredients(
	const BandstructureDomain<complex<double> >& cb_bands_, 
	const BandstructureDomain<complex<double> >& vb_bands_,
	const MatrixElements& matrix_elements_,
	const CoulombMatrixElementData& coulomb_matrix_elements_,
	const double& homogeneous_broadening_,		
	const double& optical_refractive_index_,
	const double& static_dielectric_constant_,
	const double& geometry_output_division_factor_) 
: cb_bands(cb_bands_),
  vb_bands(vb_bands_),
  matrix_elements(matrix_elements_),
  coulomb_matrix_elements(coulomb_matrix_elements_),
  homogeneous_broadening(homogeneous_broadening_),
  optical_refractive_index(optical_refractive_index_),
  static_dielectric_constant(static_dielectric_constant_),
  geometry_output_division_factor(geometry_output_division_factor_)
{
	// --------------------------------------------
	// check domains
	// --------------------------------------------
	TDKP_ASSERT(cb_bands.get_domain() == vb_bands.get_domain(), ""); 
	TDKP_ASSERT(cb_bands.get_domain() == matrix_elements.get_domain(), "");
	TDKP_ASSERT(cb_bands.get_number_of_bands() == (signed) matrix_elements.get_num_cb_bands(), "");
	TDKP_ASSERT(vb_bands.get_number_of_bands() == (signed) matrix_elements.get_num_vb_bands(), "");
	TDKP_ASSERT(homogeneous_broadening > 0.0, "");
	TDKP_ASSERT(optical_refractive_index > 0.0, "");
	TDKP_ASSERT(static_dielectric_constant > 0.0, "");
	TDKP_ASSERT(geometry_output_division_factor > 0.0, "");	
}

/** return domain of data objects (same for bands and momentum matrix element) */
const DomainMaster& CLCIngredients::get_domain() const {
	return cb_bands.get_domain();	
}

/** return cb bandstructure container */
const BandstructureDomain<complex<double> >& CLCIngredients::get_cb_bands() const {
	return cb_bands;	
}

/** return vb bandstructure container */ 
const BandstructureDomain<complex<double> >& CLCIngredients::get_vb_bands() const {
	return vb_bands;	
}

/** return matrix elements container */
const MatrixElements& CLCIngredients::get_matrix_elements() const {
	return matrix_elements;	
}

/** return coulomb matrix elements */
const CoulombMatrixElementData& CLCIngredients::get_coulomb_matrix_elements() const {
	return coulomb_matrix_elements;	
}

/** return homogeneous broadening factor gamma */
const double& CLCIngredients::get_homogeneous_broadening() const {
	return homogeneous_broadening;	
}		

/** return optical refractive index */
const double& CLCIngredients::get_optical_refractive_index() const {
	return optical_refractive_index;	
}

/** return static dielectric constant */
const double& CLCIngredients::get_static_dielectric_constant() const {
	return static_dielectric_constant;	
}

/** return length of quantized system (thats L) */
const double& CLCIngredients::get_geometry_output_division_factor() const {
	return geometry_output_division_factor;	
}

CLCIngredients1DR::CLCIngredients1DR(
	const BandstructureDomain<complex<double> >& cb_bands_, 
	const BandstructureDomain<complex<double> >& vb_bands_,
	const MatrixElements& matrix_elements_,
	const CoulombMatrixElementData& coulomb_matrix_elements_,
	const double& homogeneous_broadening_,		
	const double& optical_refractive_index_,
	const double& static_dielectric_constant_,
	const double& geometry_output_division_factor_
) : CLCIngredients(
		cb_bands_, 
		vb_bands_, 
		matrix_elements_, 
		coulomb_matrix_elements_,
		homogeneous_broadening_,
		optical_refractive_index_,
		static_dielectric_constant_,
		geometry_output_division_factor_
	)	
{
		
	TDKP_ASSERT(this->get_cb_bands().get_number_of_bands() == (signed)this->get_matrix_elements().get_num_cb_bands(), "");
	TDKP_ASSERT(this->get_vb_bands().get_number_of_bands() == (signed)this->get_matrix_elements().get_num_vb_bands(), "");		
		
	// ---------------------------------------------------
	// kmin must be somewhere around 0
	// ---------------------------------------------------
	TDKP_ASSERT(this->get_cb_bands().get_domain().get_first_point().get_coord_abs() >= 0.0, "");
	TDKP_ASSERT(this->get_cb_bands().get_domain().get_first_point().get_coord_abs() <= 1.0e-2, "this->get_cb_bands().get_domain().get_first_point().get_coord_abs() <= 1.0e-2 failed! kmin MUST BE somewhere around 0.0");
	TDKP_ASSERT(this->get_domain().get_number_of_points() > 1, "");
	TDKP_ASSERT(this->get_domain().get_last_point().get_coord_abs() > 0.0, ""); 		
		
	// ---------------------------------------------------
	// create spline objects
	// ---------------------------------------------------
	// extend bandstructure parabolically to ntimes*kmax
	unsigned int nkmax = static_cast<unsigned int>(Configuration::get_instance()->get("clc_bands_extrapolate_kmax_n_factor"));
	TDKP_ASSERT(nkmax >= 1, "clc_bands_extrapolate_kmax_n_factor >= 1");
	const unsigned int ndomain_points = this->get_domain().get_number_of_points(); 	
	vector<double> xext(ndomain_points + nkmax - 1);
	vector<double> xorig(ndomain_points);
	vector<double> y(ndomain_points + nkmax - 1);
	for(unsigned int ii = 0; ii < ndomain_points; ii++) {
		xorig[ii] = xext[ii] = this->get_domain().get_point(ii).get_coord_abs();			
	}
	for(unsigned int ii = 0; ii < nkmax - 1; ii++) {
		xext[ndomain_points + ii] = (xext[ndomain_points - 1]) * (ii + 2); 	
	}
		
	// ---------------------------------------------------
	// for all cb bands
	// ---------------------------------------------------
	for(int ii = 0; ii < this->get_cb_bands().get_number_of_bands(); ii++) {
		for(unsigned int kk = 0; kk < ndomain_points; kk++) {
			y[kk] = this->get_cb_bands().get_energy(kk,ii).real();			
		}
		// extension of bands 
		double acurve = (y[ndomain_points - 1] - y[0]) / (xext[ndomain_points - 1]*xext[ndomain_points - 1]);
		TDKP_BOUNDS_ASSERT(acurve > 0.0, "");
		for(unsigned int nn = 0; nn < nkmax - 1; nn++) {
			y[ndomain_points + nn] = y[0] + acurve * (xext[ndomain_points + nn]*xext[ndomain_points + nn]); 	
		}   
		curve_cb_bands.push_back(new Spline1D(xext,y, 0.0));				 					
	}
	// ---------------------------------------------------
	// for all vb bands
	// ---------------------------------------------------
	for(int ii = 0; ii < this->get_vb_bands().get_number_of_bands(); ii++) {
		for(unsigned int kk = 0; kk < this->get_domain().get_number_of_points(); kk++) {
			y[kk] = this->get_vb_bands().get_energy(kk,ii).real();			
		}
		// extension of bands 
		double acurve = (y[ndomain_points - 1] - y[0]) / (xext[ndomain_points - 1]*xext[ndomain_points - 1]);
		TDKP_BOUNDS_ASSERT(acurve < 0.0, "");
		for(unsigned int nn = 0; nn < nkmax - 1; nn++) {
			y[ndomain_points + nn] = y[0] + acurve * (xext[ndomain_points + nn]*xext[ndomain_points + nn]); 	
		}  		
		curve_vb_bands.push_back(new Spline1D(xext,y, 0.0));				 					
	}
	// ---------------------------------------------------
	// for all matrix elements
	// ---------------------------------------------------		
	vector<complex<double> > yc(this->get_domain().get_number_of_points());
	// spin up + spin down data?
	unsigned int num_value_sets;
	if(this->get_cb_bands().get_basis_size() == 1 && this->get_vb_bands().get_basis_size() == 1) {
		num_value_sets = 1;
	} else if(this->get_cb_bands().get_basis_size() == 1) {
		num_value_sets = 2; // matrix elements have mu_k spin up and spin down (for spin up and spin down cb band)
	} else {
		num_value_sets = 1;	
	}
	// resize vector
	curve_matrix_elements.resize(
		this->get_matrix_elements().get_num_cb_bands()
		* this->get_matrix_elements().get_num_vb_bands()
		* 3
		* num_value_sets 	
	);
	for(unsigned int cc = 0; cc < this->get_matrix_elements().get_num_cb_bands(); cc++) {
		for(unsigned int vv = 0; vv < this->get_matrix_elements().get_num_vb_bands(); vv++) {
			for(unsigned int pp = 0; pp < 3; pp++) {
				for(unsigned int ss = 0; ss < num_value_sets; ss++) {
					for(unsigned int kk = 0; kk < this->get_domain().get_number_of_points(); kk++) {
						yc[kk] = this->get_matrix_elements().get(cc,vv,pp,kk,ss);			
					}
					TDKP_BOUNDS_ASSERT(curve_matrix_elements.size() > get_storage_idx(cc,vv,pp,ss), "");
					TDKP_BOUNDS_ASSERT(curve_matrix_elements[get_storage_idx(cc,vv,pp,ss)] == 0, "");
					curve_matrix_elements[get_storage_idx(cc,vv,pp,ss)] = new CLCLinearCurve(xorig,yc);
				}
			}
		}						 				
	}	
	// ---------------------------------------------------
	// create coulomb matrix element objects
	// ---------------------------------------------------	
	const CoulombMatrixElementData& cmed = get_coulomb_matrix_elements();
	cb_cb_coulomb_matrix_elements.resize(get_cb_bands().get_number_of_bands(), 0);	
	vb_vb_coulomb_matrix_elements.resize(get_vb_bands().get_number_of_bands(), 0);
	vb_cb_coulomb_matrix_elements.resize(get_cb_bands().get_number_of_bands() * get_vb_bands().get_number_of_bands(), 0);
	
	// -----------------------------------------------
	// get type of every wavefunction we covered
	// -----------------------------------------------
	vector<unsigned int> cb_wf;
	vector<unsigned int> vb_wf;	
	for(unsigned int nn = 0; nn < cmed.get_num_wavefunctions(); nn++) {
		// -----------------------------------------------
		// check type
		// -----------------------------------------------
		bool is_cb_wf = false;
		if(cmed.get_wavefunction_type(nn).substr(0,2) == string("cb")) {
			is_cb_wf = true;	
		} else if(cmed.get_wavefunction_type(nn).substr(0,2) != string("vb")) {
			TDKP_GENERAL_EXCEPTION("unknown wavefunction type " << cmed.get_wavefunction_type(nn)); 	
		}
		// -----------------------------------------------
		// get index
		// -----------------------------------------------		
		unsigned int subband_idx = lexical_cast<unsigned int>(cmed.get_wavefunction_type(nn).substr(2));
		TDKP_ASSERT(subband_idx < cmed.get_num_wavefunctions(), "");
		if(is_cb_wf) {
			TDKP_ASSERT(subband_idx == cb_wf.size(), "data in coulomb matrix elements is not properly ordered! cb subband index must be increasing, starting from 0");
			cb_wf.push_back(nn);	
		} else {
			TDKP_ASSERT(subband_idx == vb_wf.size(), "data in coulomb matrix elements is not properly ordered! vb subband index must be increasing, starting from 0");
			vb_wf.push_back(nn);			
		}		
	}
	
	// -----------------------------------------------
	// build all objects
	// -----------------------------------------------
	TDKP_ASSERT(cb_cb_coulomb_matrix_elements.size() == cb_wf.size(), "");
	TDKP_ASSERT(vb_vb_coulomb_matrix_elements.size() == vb_wf.size(), "");
	// for all cb bands
	for(unsigned int cc = 0; cc < cb_wf.size(); cc++) {
		cb_cb_coulomb_matrix_elements[cc] = create_clc_coulomb_lambda_object(
			cmed.get_q_values(),
			cmed.get_matrix_element(cb_wf[cc],cb_wf[cc],cb_wf[cc],cb_wf[cc]),
			this->get_domain().get_dimension()
		);
		// for all vb bands
		for(unsigned int vv = 0; vv < vb_wf.size(); vv++) {
			vb_cb_coulomb_matrix_elements[get_cb_bands().get_number_of_bands() * vv + cc] = create_clc_coulomb_lambda_object(
				cmed.get_q_values(),
				cmed.get_matrix_element(vb_wf[vv],cb_wf[cc],cb_wf[cc],vb_wf[vv]),
				this->get_domain().get_dimension()
			);
		}								
	}
	// for all vb bands
	for(unsigned int vv = 0; vv < vb_wf.size(); vv++) {	
		vb_vb_coulomb_matrix_elements[vv] = create_clc_coulomb_lambda_object(
			cmed.get_q_values(),
			cmed.get_matrix_element(vb_wf[vv],vb_wf[vv],vb_wf[vv],vb_wf[vv]),
			this->get_domain().get_dimension()
		);		
	}	
	
}

/** create correct clc coulomb lambda object (well/wire/bulk) */
CLCRadialCoulombLambda* CLCIngredients1DR::create_clc_coulomb_lambda_object(
	const vector<double>& q_values, 
	const vector<cplx>& data, 
	unsigned int kspace_dimension
) const {
	
	switch(kspace_dimension) {
		case 1:
			return new CLCCoulombLambdaWire(q_values, data);
		case 2:
			return new CLCRadialCoulombLambdaWell(
				q_values, 
				data, 
				Configuration::get_instance()->get("clc_radial_angular_integration_well_num_points")
			);
		case 3:
			return new CLCRadialCoulombLambdaBulk(
				q_values, 
				data,
				Configuration::get_instance()->get("clc_radial_angular_integration_bulk_num_points")
			);
		default:
			TDKP_GENERAL_EXCEPTION("invalid k space dimension");			
	}
}

/** return diagonal coulomb matrix element object with n1 = n4 = n2 = n3 = cb_dix */
const CLCRadialCoulombLambda& CLCIngredients1DR::get_coulomb_matrix_element_cb_cb(
	unsigned int cb_idx
) const {
	TDKP_ASSERT(cb_idx < cb_cb_coulomb_matrix_elements.size(), "");
	return *cb_cb_coulomb_matrix_elements[cb_idx];	
}

/** return diagonal coulomb matrix element object with n1 = n4 = vb_idx, n2 = n3 = cb_dix */	
const CLCRadialCoulombLambda& CLCIngredients1DR::get_coulomb_matrix_element_vb_cb(
	unsigned int vb_idx,
	unsigned int cb_idx
) const {
	TDKP_ASSERT(vb_idx * get_cb_bands().get_number_of_bands() + cb_idx < vb_cb_coulomb_matrix_elements.size(), "");
	return *vb_cb_coulomb_matrix_elements[vb_idx * get_cb_bands().get_number_of_bands() + cb_idx];	
}

/** return diagonal coulomb matrix element object with n1 = n4 = n2 = n3 = vb_dix */
const CLCRadialCoulombLambda& CLCIngredients1DR::get_coulomb_matrix_element_vb_vb(
	unsigned int vb_idx
) const {
	TDKP_ASSERT(vb_idx < vb_vb_coulomb_matrix_elements.size(), "");
	return *vb_vb_coulomb_matrix_elements[vb_idx];	
}

/** return idx where we stored the appropriate curve */
unsigned int CLCIngredients1DR::get_storage_idx(		
	unsigned int cb_idx, 
	unsigned int vb_idx, 
	unsigned int dir, 
	unsigned int value_index
) const {
	
	unsigned int num_value_sets;
	if(this->get_cb_bands().get_basis_size() == 1 && this->get_vb_bands().get_basis_size() == 1) {
		num_value_sets = 1;
	} else if(this->get_cb_bands().get_basis_size() == 1) {
		num_value_sets = 2; // matrix elements have mu_k spin up and spin down (for spin up and spin down cb band)
	} else {
		num_value_sets = 1;	
	}	
	
	return value_index 
	     + num_value_sets * dir
	     + num_value_sets * 3 * vb_idx
	     + num_value_sets * 3 * this->get_matrix_elements().get_num_vb_bands() * cb_idx;	     
}

/** delete curve objects */
CLCIngredients1DR::~CLCIngredients1DR() {
		
	for(unsigned int ii = 0; ii < curve_cb_bands.size(); ii++) {
		delete curve_cb_bands[ii]; curve_cb_bands[ii] = 0; 	
	}
	for(unsigned int ii = 0; ii < curve_vb_bands.size(); ii++) {
		delete curve_vb_bands[ii]; curve_vb_bands[ii] = 0; 	
	}
	for(unsigned int ii = 0; ii < curve_matrix_elements.size(); ii++) {
		delete curve_matrix_elements[ii]; curve_matrix_elements[ii] = 0; 	
	} 			
	for(unsigned int ii = 0; ii < cb_cb_coulomb_matrix_elements.size(); ii++) {
		delete cb_cb_coulomb_matrix_elements[ii];
		cb_cb_coulomb_matrix_elements[ii] = 0; 	
	}
	for(unsigned int ii = 0; ii < vb_cb_coulomb_matrix_elements.size(); ii++) {
		delete vb_cb_coulomb_matrix_elements[ii];
		vb_cb_coulomb_matrix_elements[ii] = 0; 	
	}
	for(unsigned int ii = 0; ii < vb_vb_coulomb_matrix_elements.size(); ii++) {
		delete vb_vb_coulomb_matrix_elements[ii];
		vb_vb_coulomb_matrix_elements[ii] = 0; 	
	}		
}

const Spline1D& CLCIngredients1DR::get_cb_band(unsigned int cb_idx) const {
	TDKP_BOUNDS_ASSERT(curve_cb_bands.size() > cb_idx, "");
	return *curve_cb_bands[cb_idx];	
}
const Spline1D& CLCIngredients1DR::get_vb_band(unsigned int vb_idx) const {
	TDKP_BOUNDS_ASSERT(curve_vb_bands.size() > vb_idx, "");
	return *curve_vb_bands[vb_idx];
}
	
const CLCLinearCurve& CLCIngredients1DR::get_momentum_matrix_element(
	unsigned int cb_idx, 
	unsigned int vb_idx, 
	unsigned int dir, 
	unsigned int value_index
) const {
	unsigned int idx = get_storage_idx(cb_idx, vb_idx, dir, value_index);
	TDKP_BOUNDS_ASSERT(idx < curve_matrix_elements.size(), "");
	TDKP_BOUNDS_ASSERT(curve_matrix_elements[idx] != 0, "idx = " << idx << ", size = " << curve_matrix_elements.size() << ", cb_idx = " << cb_idx << ", vb_idx = " << vb_idx << ", dir = " << dir << ", value_idx = " << value_index);
	return *curve_matrix_elements[idx]; 	
}	



}
