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

#include "tdkp/kpmatrices/KPMatrix6x6Wurtzite.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/common/Vector3D.h"

namespace tdkp {
	
const int KPMatrix6x6Wurtzite::sparsity_pattern[44]= {
	0,0,	0,1,	0,2,	
	1,0,	1,1,	1,2,					1,5,
	2,0,	2,1,	2,2,			2,4,
							3,3,	3,4,	3,5,
					4,2,	4,3,	4,4,	4,5,
			5,1,			5,3,	5,4,	5,5
};	
const int KPMatrix6x6Wurtzite::diagonal_pattern[6] = {
	0, 4, 9, 11, 16, 21	
};	

bool KPMatrix6x6Wurtzite::check_solution_type_available(KPSolutionType type) const {
	return type == electrons ? false : true;	
}

inline int KPMatrix6x6Wurtzite::get_idx(int ii, int jj) const {
	for(int ss = 0; ss < 22; ss++) {
		if(sparsity_pattern[2*ss] == ii && sparsity_pattern[2*ss+1] == jj) {
			return ss;	
		}	
	}	
	TDKP_GENERAL_EXCEPTION("invalid idx: " << ii << " and " << jj);
}
		
const int* KPMatrix6x6Wurtzite::get_sparsity_pattern(int& length) const {
	length = 22;
	return sparsity_pattern;	
}	
const int* KPMatrix6x6Wurtzite::get_diagonal_pattern(int& length) const {
	length = 6;
	return diagonal_pattern;
}	

KPMatrix6x6Wurtzite::KPMatrix6x6Wurtzite() {
	this->allocate_matrix_space();
	this->material = 0;
	this->have_first_order_terms = true;
}

KPMatrix6x6Wurtzite::KPMatrix6x6Wurtzite(const Material* mat) {
	this->set_material(mat);			
	// matrix space must be allocated from derived class
	// allocate space calls virtual function, which is defined in derived class
	// so it's not available in the default constructor of the base class
	this->allocate_matrix_space();
	this->have_first_order_terms = true;			
}

KPMatrix6x6Wurtzite::~KPMatrix6x6Wurtzite() {
	this->material = 0;
}

void KPMatrix6x6Wurtzite::init_base_and_strain_matrix() {
	TDKP_ASSERT(this->material != 0,"this->material != 0");
	
	// ---------------------------------------------------------
	// wurtzite kp matrix
	// taken from phys rev b 54, 1996, chuang and chang,
	// kp method for strained wurzite semiconductors
	// 
	// the operator ordering is a straight forward assumption/extension
	// of the foreman ordering of the zinc-blende matrix 
	// 
	// the operator ordering was derived using (26) in chuang, 
	// applying the operator ordering split
	// N1 kxky -> kx N1p ky + ky N1m kx
	// and the same for N2
	// 
	// as A5 = N1 and A6 = N2 (33) and using the basis transformation
	// between X> Y> Z> and chuangs basis, one obtaines the operator
	// ordering here
	//
	// we assume here that the right distribution between A5p and A5m is
	// given by the optimal value of nonellipticity (same for A6).
	// and guess what: this is optimal if A5m and A6m = 0 and A5p = A5,
	// A6p = A6 (its really elliptic)
	// ---------------------------------------------------------		
	const complex<double> i(0,1);
	const double hbar_2m0 = constants::hbar_square / (2.0 * constants::m0);		
	const short splitting_symmetric    = 0;
	const short splitting_user_defined = 1;
	const short splitting_auto         = 2;
	short operator_choice = splitting_auto;
	if(Configuration::get_instance()->get("kpmatrix_disable_foreman_enable_symmetrization") == 1.0) {
		operator_choice = splitting_symmetric;		
	} else if(this->material->is_set("hole_effective_mass_A5_minus") && this->material->is_set("hole_effective_mass_A6_minus")) {
		operator_choice = splitting_user_defined;	
	}
			
	// ------------------------------------------------------
	// get material paramters
	// ------------------------------------------------------
	const double vb_edge  = this->material->get("valence_band_edge");
	const double delta_1  = this->material->get("spin_orbit_split_delta_1");
	const double delta_2  = this->material->get("spin_orbit_split_delta_2");
	const double delta_3  = this->material->get("spin_orbit_split_delta_3");	
	const double A1 = hbar_2m0 * this->material->get("hole_effective_mass_A1");
	const double A2 = hbar_2m0 * this->material->get("hole_effective_mass_A2");
	const double A3 = hbar_2m0 * this->material->get("hole_effective_mass_A3");
	const double A4 = hbar_2m0 * this->material->get("hole_effective_mass_A4");
	const double A5 = hbar_2m0 * this->material->get("hole_effective_mass_A5");
	const double A6 = hbar_2m0 * this->material->get("hole_effective_mass_A6");
	double       A7 = 0.0;
	
	if(Configuration::get_instance()->get("kpmatrix_wurtzite_include_spin_orbit_interaction") == 1.0) {
		A7 = this->material->get("hole_spin_orbit_A7");	
	}
				
	// ------------------------------------------------------
	// determine parameter splitting
	// ------------------------------------------------------
	double A5p, A5m, A6p, A6m;
	determine_splitted_parameters(A5p, A5m, A6p, A6m);
	vector<double> eigenvalues;
	
	// ------------------------------------------------------
	// check non-ellipticity
	// ------------------------------------------------------
	WurtziteEffectiveMassParams params;
	params.A1 = A1; params.A2 = A2; params.A3 = A3; params.A4 = A4; params.A5 = A5;
	params.A6 = A6; params.A5m = A5m; params.A6m = A6m;	
	determine_ellipticity_eigenvalues(params, eigenvalues);
	const double nonellipticity = get_nonellipticity_ratio(eigenvalues);	
	
	A5p *= hbar_2m0;
	A6p *= hbar_2m0;
	A5m *= hbar_2m0;
	A6m *= hbar_2m0; 
	
	// ------------------------------------------------------
	// get strain parameters
	// ------------------------------------------------------
	double D1, D2, D3, D4, D5, D6;
		
	if(this->material->is_set("strain_potential_D1") &&
	   this->material->is_set("strain_potential_D2") &&
	   this->material->is_set("strain_potential_D3") &&
	   this->material->is_set("strain_potential_D4") &&
	   this->material->is_set("strain_potential_D5") &&
	   this->material->is_set("strain_potential_D6")) {
	   	
		D1 = this->material->get("strain_potential_D1");
		D2 = this->material->get("strain_potential_D2");
		D3 = this->material->get("strain_potential_D3");
		D4 = this->material->get("strain_potential_D4");
		D5 = this->material->get("strain_potential_D5");
		D6 = this->material->get("strain_potential_D6"); 	
		   	
	   	this->strain_dependence_available = true;
	} else {
		this->strain_dependence_available = false;
		Logger::get_instance()->emit(LOG_WARN, "at least one of the strain potentials is not set. strain dependence is therefore not available!"); 	
		D1 = D2 = D3 = D4 = D5 = D6 = 0.0;
	}
	
	// -------------------------------------------------------
	// determine band edge shifiting of vb
	// attention: cb_edge - vb_edge gives the bandgap
	// but in chuangs kp matrix, the upper vb is not at 0 but at
	// Ev + delta_1 + delta_2
	// vb_edge really referes to the top. so what do we do?
	// we shift Ev' = Ev - delta_1 - delta_2
	// so we have the top at 0
	// this is important for the effmass cb edge
	// see eq. (11) in chuang and fig. 2.
	// -------------------------------------------------------
	// much more complicated, we include all possible zone center energies and take the top
	double E1 = delta_1 + delta_2;
	double tt = sqrt(((delta_1 - delta_2) / 2.0) * ((delta_1 - delta_2) / 2.0) + 2.0 * delta_3 * delta_3);
	double E2 = (delta_1 - delta_2) / 2.0 + tt;
	double E3 = (delta_1 - delta_2) / 2.0 - tt;
	double vb_shift = 0.0;
	if(E1 > E2 && E1 > E3) {
		vb_shift = E1;	
	} else if(E2 > E3) {
		vb_shift = E2;	
	} else {
		vb_shift = E3;	
	}
	
	// -------------------------------------------------------
	// tell the user the paramters
	// -------------------------------------------------------

	ostringstream sout; 
	sout.setf(ios::left);
 	sout << "kp 6x6 wurtzite matrix for material " << this->get_material_name() << " uses following parameters:\n"
 		 << "hole effective mass A1   = " << A1 / hbar_2m0<< "\n"
 		 << "hole effective mass A2   = " << A2 / hbar_2m0 << "\n"
 		 << "hole effective mass A3   = " << A3 / hbar_2m0 << "\n"
 		 << "hole effective mass A4   = " << A4 / hbar_2m0 << "\n"
 		 << "hole effective mass A5   = " << A5 / hbar_2m0 << "\n"
 		 << "hole effective mass A6   = " << A6 / hbar_2m0 << "\n"
 		 << "hole spin orbit     A7   = " << setw(15) << A7 << " " << (Configuration::get_instance()->get("kpmatrix_wurtzite_include_spin_orbit_interaction") == 1.0 ? "(enabled)" : "(disabled)") << "\n"
 		 << "spin orbit split delta 1 = " << delta_1 << "\n"
 		 << "spin orbit split delta 2 = " << delta_2 << "\n"
 		 << "spin orbit split delta 3 = " << delta_3 << "\n"
 		 << "valence band edge        = " << vb_edge << "\n"
 		 << "vb shift (delta offset)  = " << vb_shift << "\n"
 		 << "operator ordering        = ";
	if(operator_choice == splitting_symmetric) {
		sout << "symmetric splitting\n";
	} else if(operator_choice == splitting_user_defined) {
		sout << "user defined\n";
	} else {
		sout << "automatically setting A- to 0\n";	
	}
	sout  << "    => A5+ " << A5p / hbar_2m0 << ", A5- => " << A5m / hbar_2m0 
	      << ", A6+ => " << A6p  / hbar_2m0 << ", A6m => " << A6m / hbar_2m0 << "\n"
	      << "nonellipticity           = " << nonellipticity << "\n";

	if(this->strain_dependence_available) {
		sout << "strain potential D1      = " << D1 << "\n"
		     << "strain potential D2      = " << D2 << "\n"
		     << "strain potential D3      = " << D3 << "\n"
		     << "strain potential D4      = " << D4 << "\n"
		     << "strain potential D5      = " << D5 << "\n"
		     << "strain potential D6      = " << D6 << "\n";		     
	}	      
	       
	if(!surpress_output()) {
 		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	}
 	this->initialized = true;
 	
 	// ------------------------------------------------------
 	// build D_DX D_DX
 	// ------------------------------------------------------
/*(0,0)*/ this->second_order_strainless[D_DX][D_DX][get_idx(0,0)] =   A2 + A4;
/*(0,1)*/ this->second_order_strainless[D_DX][D_DX][get_idx(0,1)] = - A5;
/*(1,0)*/ this->second_order_strainless[D_DX][D_DX][get_idx(1,0)] = - A5;
/*(1,1)*/ this->second_order_strainless[D_DX][D_DX][get_idx(1,1)] =   A2 + A4;
/*(2,2)*/ this->second_order_strainless[D_DX][D_DX][get_idx(2,2)] =   A2;
/*(3,3)*/ this->second_order_strainless[D_DX][D_DX][get_idx(3,3)] =   A2 + A4;
/*(3,4)*/ this->second_order_strainless[D_DX][D_DX][get_idx(3,4)] = - A5;
/*(4,3)*/ this->second_order_strainless[D_DX][D_DX][get_idx(4,3)] = - A5;
/*(4,4)*/ this->second_order_strainless[D_DX][D_DX][get_idx(4,4)] =   A2 + A4;
/*(5,5)*/ this->second_order_strainless[D_DX][D_DX][get_idx(5,5)] =   A2;

 	// ------------------------------------------------------
 	// build D_DY D_DY
 	// ------------------------------------------------------
/*(0,0)*/ this->second_order_strainless[D_DY][D_DY][get_idx(0,0)] =   A2 + A4;
/*(0,1)*/ this->second_order_strainless[D_DY][D_DY][get_idx(0,1)] =   A5;
/*(1,0)*/ this->second_order_strainless[D_DY][D_DY][get_idx(1,0)] =   A5;
/*(1,1)*/ this->second_order_strainless[D_DY][D_DY][get_idx(1,1)] =   A2 + A4;
/*(2,2)*/ this->second_order_strainless[D_DY][D_DY][get_idx(2,2)] =   A2;
/*(3,3)*/ this->second_order_strainless[D_DY][D_DY][get_idx(3,3)] =   A2 + A4;
/*(3,4)*/ this->second_order_strainless[D_DY][D_DY][get_idx(3,4)] =   A5;
/*(4,3)*/ this->second_order_strainless[D_DY][D_DY][get_idx(4,3)] =   A5;
/*(4,4)*/ this->second_order_strainless[D_DY][D_DY][get_idx(4,4)] =   A2 + A4;
/*(5,5)*/ this->second_order_strainless[D_DY][D_DY][get_idx(5,5)] =   A2;
 	
 	// ------------------------------------------------------
 	// build D_DZ D_DZ
 	// ------------------------------------------------------
/*(0,0)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(0,0)] =   A1 + A3;
/*(1,1)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(1,1)] =   A1 + A3;
/*(2,2)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(2,2)] =   A1;
/*(3,3)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(3,3)] =   A1 + A3;
/*(4,4)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(4,4)] =   A1 + A3;
/*(5,5)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(5,5)] =   A1; 	
			
 	// ------------------------------------------------------
 	// build D_DX D_DY
 	// ------------------------------------------------------
/*(0,0)*/ this->second_order_strainless[D_DX][D_DY][get_idx(0,0)] =   i * (A5m - A5p);
/*(0,1)*/ this->second_order_strainless[D_DX][D_DY][get_idx(0,1)] = - i * (A5m + A5p);
/*(1,0)*/ this->second_order_strainless[D_DX][D_DY][get_idx(1,0)] =   i * (A5m + A5p);
/*(1,1)*/ this->second_order_strainless[D_DX][D_DY][get_idx(1,1)] = - i * (A5m - A5p);
/*(3,3)*/ this->second_order_strainless[D_DX][D_DY][get_idx(3,3)] = - i * (A5m - A5p);
/*(3,4)*/ this->second_order_strainless[D_DX][D_DY][get_idx(3,4)] =   i * (A5m + A5p);
/*(4,3)*/ this->second_order_strainless[D_DX][D_DY][get_idx(4,3)] = - i * (A5m + A5p);
/*(4,4)*/ this->second_order_strainless[D_DX][D_DY][get_idx(4,4)] =   i * (A5m - A5p);

 	// ------------------------------------------------------
 	// build D_DY D_DX
 	// ------------------------------------------------------
/*(0,0)*/ this->second_order_strainless[D_DY][D_DX][get_idx(0,0)] = - i * (A5m - A5p);
/*(0,1)*/ this->second_order_strainless[D_DY][D_DX][get_idx(0,1)] = - i * (A5m + A5p);
/*(1,0)*/ this->second_order_strainless[D_DY][D_DX][get_idx(1,0)] =   i * (A5m + A5p);
/*(1,1)*/ this->second_order_strainless[D_DY][D_DX][get_idx(1,1)] =   i * (A5m - A5p);
/*(3,3)*/ this->second_order_strainless[D_DY][D_DX][get_idx(3,3)] =   i * (A5m - A5p);
/*(3,4)*/ this->second_order_strainless[D_DY][D_DX][get_idx(3,4)] =   i * (A5m + A5p);
/*(4,3)*/ this->second_order_strainless[D_DY][D_DX][get_idx(4,3)] = - i * (A5m + A5p);
/*(4,4)*/ this->second_order_strainless[D_DY][D_DX][get_idx(4,4)] = - i * (A5m - A5p);

 	// ------------------------------------------------------
 	// build D_DX D_DZ
 	// ------------------------------------------------------
/*(0,2)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(0,2)] = - A6p;
/*(1,2)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(1,2)] =   A6p;
/*(2,0)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(2,0)] = - A6m;
/*(2,1)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(2,1)] =   A6m;
/*(3,5)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(3,5)] =   A6p;
/*(4,5)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(4,5)] = - A6p;
/*(5,3)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(5,3)] =   A6m;
/*(5,4)*/ this->second_order_strainless[D_DX][D_DZ][get_idx(5,4)] = - A6m;

 	// ------------------------------------------------------
 	// build D_DZ D_DX
 	// ------------------------------------------------------
/*(0,2)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(0,2)] = - A6m;
/*(1,2)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(1,2)] =   A6m;
/*(2,0)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(2,0)] = - A6p;
/*(2,1)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(2,1)] =   A6p;
/*(3,5)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(3,5)] =   A6m;
/*(4,5)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(4,5)] = - A6m;
/*(5,3)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(5,3)] =   A6p;
/*(5,4)*/ this->second_order_strainless[D_DZ][D_DX][get_idx(5,4)] = - A6p;

 	// ------------------------------------------------------
 	// build D_DY D_DZ
 	// ------------------------------------------------------
/*(0,2)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(0,2)] = - i * A6p;
/*(1,2)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(1,2)] = - i * A6p;
/*(2,0)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(2,0)] =   i * A6m;
/*(2,1)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(2,1)] =   i * A6m;
/*(3,5)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(3,5)] = - i * A6p;
/*(4,5)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(4,5)] = - i * A6p;
/*(5,3)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(5,3)] =   i * A6m;
/*(5,4)*/ this->second_order_strainless[D_DY][D_DZ][get_idx(5,4)] =   i * A6m;

 	// ------------------------------------------------------
 	// build D_DZ D_DY
 	// ------------------------------------------------------
/*(0,2)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(0,2)] = - i * A6m;
/*(1,2)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(1,2)] = - i * A6m;
/*(2,0)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(2,0)] =   i * A6p;
/*(2,1)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(2,1)] =   i * A6p;
/*(3,5)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(3,5)] = - i * A6m;
/*(4,5)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(4,5)] = - i * A6m;
/*(5,3)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(5,3)] =   i * A6p;
/*(5,4)*/ this->second_order_strainless[D_DZ][D_DY][get_idx(5,4)] =   i * A6p; 	
 	
	// ------------------------------------------------------
	// build zero order strainless
	// ------------------------------------------------------
/*(0,0)*/ this->zero_order_strainless[get_idx(0,0)] = vb_edge + delta_1 + delta_2 - vb_shift;
/*(1,1)*/ this->zero_order_strainless[get_idx(1,1)] = vb_edge + delta_1 - delta_2 - vb_shift;
/*(1,5)*/ this->zero_order_strainless[get_idx(1,5)] = constants::sqrt2 * delta_3;
/*(2,2)*/ this->zero_order_strainless[get_idx(2,2)] = vb_edge - vb_shift;
/*(2,4)*/ this->zero_order_strainless[get_idx(2,4)] = constants::sqrt2 * delta_3;
/*(3,3)*/ this->zero_order_strainless[get_idx(3,3)] = vb_edge + delta_1 + delta_2 - vb_shift;
/*(4,2)*/ this->zero_order_strainless[get_idx(4,2)] = constants::sqrt2 * delta_3;
/*(4,4)*/ this->zero_order_strainless[get_idx(4,4)] = vb_edge + delta_1 - delta_2 - vb_shift;
/*(5,1)*/ this->zero_order_strainless[get_idx(5,1)] = constants::sqrt2 * delta_3;
/*(5,5)*/ this->zero_order_strainless[get_idx(5,5)] = vb_edge - vb_shift;
		
	// ------------------------------------------------------
	// spin orbit coupling a la ren et al (appl. phys. lett 74 (1999))
	// according to ren, they recently found a new parameter A7 ;-)
	// in fact, chuang neglected the spin orbit coupling but introduced
	// it phenomenologically via a energy split at k = 0
	//
	// if you don't neglect it, you get these first order terms here 
	// by taking the bible (pikus and bir) on page 330, enter that matrix
	// into mathematica, build the transformation matrix for the unitary
	// basis transformation into the my basis (or chuangs basis, but
	// be careful, chuangs resulting kp matrix is wrong in the H and K terms
	// (or at least, i get something different. if you follow his derivation, 
	// you won't get his results). but anyway, my stuff is consistent with
	// pikus and bir, the bible. therefore the veppo basis transformation is 
	// given by
	//      {0, 0, 0, 0, 0, -I},
	//      {0, I, 0, 0, 0, 0},
	//      {0, 0, 0, -1, 0, 0},
    //		{1, 0, 0, 0, 0, 0},
    //		{0, 0, 0, 0, -1, 0},
    //		{0, 0, I, 0, 0, 0}
	// using pikus and birs matrix, the A7 terms in my basis
	// are given by   U HPB conjTrans(U)	
	// --------------------------------------------------------		
	if(Configuration::get_instance()->get("kpmatrix_wurtzite_include_spin_orbit_interaction") == 1.0) {
			
	// ------------------------------------------------------
	// first order strainless D_DX
	// ------------------------------------------------------	
/*(0,2)*/ this->first_order_strainless[op_left][D_DX][get_idx(0,2)]  = - i * A7;
/*(1,2)*/ this->first_order_strainless[op_left][D_DX][get_idx(1,2)]  = - i * A7;
/*(2,0)*/ this->first_order_strainless[op_right][D_DX][get_idx(2,0)] = - i * A7;
/*(2,1)*/ this->first_order_strainless[op_right][D_DX][get_idx(2,1)] = - i * A7;
/*(3,5)*/ this->first_order_strainless[op_left][D_DX][get_idx(3,5)]  = - i * A7;
/*(4,5)*/ this->first_order_strainless[op_left][D_DX][get_idx(4,5)]  = - i * A7;
/*(5,3)*/ this->first_order_strainless[op_right][D_DX][get_idx(5,3)] = - i * A7;
/*(5,4)*/ this->first_order_strainless[op_right][D_DX][get_idx(5,4)] = - i * A7;	

	// ------------------------------------------------------
	// first order strainless D_DY
	// ------------------------------------------------------
/*(0,2)*/ this->first_order_strainless[op_left][D_DY][get_idx(0,2)]  =   A7;
/*(1,2)*/ this->first_order_strainless[op_left][D_DY][get_idx(1,2)]  = - A7;
/*(2,0)*/ this->first_order_strainless[op_right][D_DY][get_idx(2,0)] = - A7;
/*(2,1)*/ this->first_order_strainless[op_right][D_DY][get_idx(2,1)] =   A7;
/*(3,5)*/ this->first_order_strainless[op_left][D_DY][get_idx(3,5)]  = - A7;
/*(4,5)*/ this->first_order_strainless[op_left][D_DY][get_idx(4,5)]  =   A7;
/*(5,3)*/ this->first_order_strainless[op_right][D_DY][get_idx(5,3)] =   A7;
/*(5,4)*/ this->first_order_strainless[op_right][D_DY][get_idx(5,4)] = - A7;
		
	}
	
	
	// ------------------------------------------------------
	// strain dependence
	// ------------------------------------------------------
	if(this->strain_dependence_available) {		
		const short EX = 0; 
		const short EY = 1;
		const short EZ = 2;	
				

		// ---------------------------------------------------
		// EXX
		// ---------------------------------------------------		
/*(0,0)*/ this->zero_order_strain_dependent[EX][EX][get_idx(0,0)] =   D2 + D4;
/*(0,1)*/ this->zero_order_strain_dependent[EX][EX][get_idx(0,1)] = - D5;
/*(1,0)*/ this->zero_order_strain_dependent[EX][EX][get_idx(1,0)] = - D5;
/*(1,1)*/ this->zero_order_strain_dependent[EX][EX][get_idx(1,1)] =   D2 + D4;
/*(2,2)*/ this->zero_order_strain_dependent[EX][EX][get_idx(2,2)] =   D2;
/*(3,3)*/ this->zero_order_strain_dependent[EX][EX][get_idx(3,3)] =   D2 + D4;
/*(3,4)*/ this->zero_order_strain_dependent[EX][EX][get_idx(3,4)] = - D5;
/*(4,3)*/ this->zero_order_strain_dependent[EX][EX][get_idx(4,3)] = - D5;
/*(4,4)*/ this->zero_order_strain_dependent[EX][EX][get_idx(4,4)] =   D2 + D4;
/*(5,5)*/ this->zero_order_strain_dependent[EX][EX][get_idx(5,5)] =   D2;

		// -------------------------------------------------
		// EYY
		// ---------------------------------------------------		
/*(0,0)*/ this->zero_order_strain_dependent[EY][EY][get_idx(0,0)] =   D2 + D4;
/*(0,1)*/ this->zero_order_strain_dependent[EY][EY][get_idx(0,1)] =   D5;
/*(1,0)*/ this->zero_order_strain_dependent[EY][EY][get_idx(1,0)] =   D5;
/*(1,1)*/ this->zero_order_strain_dependent[EY][EY][get_idx(1,1)] =   D2 + D4;
/*(2,2)*/ this->zero_order_strain_dependent[EY][EY][get_idx(2,2)] =   D2;
/*(3,3)*/ this->zero_order_strain_dependent[EY][EY][get_idx(3,3)] =   D2 + D4;
/*(3,4)*/ this->zero_order_strain_dependent[EY][EY][get_idx(3,4)] =   D5;
/*(4,3)*/ this->zero_order_strain_dependent[EY][EY][get_idx(4,3)] =   D5;
/*(4,4)*/ this->zero_order_strain_dependent[EY][EY][get_idx(4,4)] =   D2 + D4;
/*(5,5)*/ this->zero_order_strain_dependent[EY][EY][get_idx(5,5)] =   D2;

		// ---------------------------------------------------
		// EZZ
		// ---------------------------------------------------		
/*(0,0)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(0,0)] = D1 + D3 ;
/*(1,1)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(1,1)] = D1 + D3;
/*(2,2)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(2,2)] = D1;
/*(3,3)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(3,3)] = D1 + D3;
/*(4,4)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(4,4)] = D1 + D3;
/*(5,5)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(5,5)] = D1;

		// ---------------------------------------------------
		// EXY
		// ---------------------------------------------------
/*(0,1)*/ this->zero_order_strain_dependent[EX][EY][get_idx(0,1)] = - i * 2.0 * D5;
/*(1,0)*/ this->zero_order_strain_dependent[EX][EY][get_idx(1,0)] =   i * 2.0 * D5;
/*(3,4)*/ this->zero_order_strain_dependent[EX][EY][get_idx(3,4)] =   i * 2.0 * D5;
/*(4,3)*/ this->zero_order_strain_dependent[EX][EY][get_idx(4,3)] = - i * 2.0 * D5;

		// ---------------------------------------------------
		// EXZ
		// ---------------------------------------------------
/*(0,2)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(0,2)] =   i * D6;
/*(1,2)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(1,2)] =   i * D6;
/*(2,0)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(2,0)] = - i * D6;
/*(2,1)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(2,1)] = - i * D6;
/*(3,5)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(3,5)] =   i * D6;
/*(4,5)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(4,5)] =   i * D6;
/*(5,3)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(5,3)] = - i * D6;
/*(5,4)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(5,4)] = - i * D6;		  

		// ---------------------------------------------------
		// EYZ
		// ---------------------------------------------------		
/*(0,2)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(0,2)] = - D6;
/*(1,2)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(1,2)] =   D6;
/*(2,0)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(2,0)] = - D6;
/*(2,1)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(2,1)] =   D6;
/*(3,5)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(3,5)] =   D6;
/*(4,5)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(4,5)] = - D6;
/*(5,3)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(5,3)] =   D6;
/*(5,4)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(5,4)] = - D6;
		
	}			
						
}


void KPMatrix6x6Wurtzite::determine_splitted_parameters(double &A5p, double &A5m, double &A6p, double &A6m) const {
	const double A5 = this->material->get("hole_effective_mass_A5");
	const double A6 = this->material->get("hole_effective_mass_A6");
	// --------------------------------------------------
	// symmetrize of somebody is that stupid ...
	// --------------------------------------------------
	if(Configuration::get_instance()->get("kpmatrix_disable_foreman_enable_symmetrization") == 1.0) {
		A5m = A5p = A5 / 2.0;
		A6m = A6p = A6 / 2.0;			
	} else {
		// -------------------------------------------------
		// read splitting from material database
		// -------------------------------------------------
		if(this->material->is_set("hole_effective_mass_A5_minus") && this->material->is_set("hole_effective_mass_A6_minus")) {
			A5m = this->material->get("hole_effective_mass_A5_minus");
			A6m = this->material->get("hole_effective_mass_A6_minus");
			A5p = A5 - A5m;
			A6p = A6 - A6m;
		} else {
			// -------------------------------------------------
			// optimal splitting 
			// -------------------------------------------------
			// although we have all the code to determine the optimal
			// splitting, we don't need it: the optimal splitting is
			// A?+ = A? and A?- = 0 ...
			// seems to be something symmetry caused ... 
			A5p = A5;
			A6p = A6;
			A5m = 0.0;
			A6m = 0.0;
		}
	}
	
}



} // end of namespace
