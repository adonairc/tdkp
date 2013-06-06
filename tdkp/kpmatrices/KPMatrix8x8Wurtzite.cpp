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

#include "tdkp/kpmatrices/KPMatrix8x8Wurtzite.h"
#include "tdkp/common/Configuration.h"

namespace tdkp {

/* sparsity pattern of the kp matrix
 * 
 * attention, i moved the conduction bands to the higher indices
 * so, using the notation in chuang, we have his basis in the 
 * sequence: u1, u2, u3, u4, u5, iS>up, iS>down
 */ 	
const int KPMatrix8x8Wurtzite::sparsity_pattern[72]= {
	0,0,	0,1,	0,2,							0,6,
	1,0,	1,1,	1,2,					1,5,	1,6,
	2,0,	2,1,	2,2,			2,4,			2,6,
							3,3,	3,4,	3,5,			3,7,
					4,2,	4,3,	4,4,	4,5,			4,7,
			5,1,			5,3,	5,4,	5,5,			5,7,
	6,0,	6,1,	6,2,							6,6,
							7,3,	7,4,	7,5,			7,7			
};	

const int KPMatrix8x8Wurtzite::diagonal_pattern[8] = {
	0, 5, 11, 14, 20, 26, 31, 35	
};	

bool KPMatrix8x8Wurtzite::check_solution_type_available(KPSolutionType type) const {
	return true; // we have both ...	
}

inline int KPMatrix8x8Wurtzite::get_idx(int ii, int jj) const {
	for(int ss = 0; ss < 36; ss++) {
		if(sparsity_pattern[2*ss] == ii && sparsity_pattern[2*ss+1] == jj) {
			return ss;	
		}	
	}	
	TDKP_GENERAL_EXCEPTION("invalid idx: " << ii << " and " << jj);
}
		
const int* KPMatrix8x8Wurtzite::get_sparsity_pattern(int& length) const {
	length = 36;
	return sparsity_pattern;	
}	
const int* KPMatrix8x8Wurtzite::get_diagonal_pattern(int& length) const {
	length = 8;
	return diagonal_pattern;
}	

KPMatrix8x8Wurtzite::KPMatrix8x8Wurtzite() {
	this->allocate_matrix_space();
	this->material = 0;
	this->have_first_order_terms = true;

}

KPMatrix8x8Wurtzite::KPMatrix8x8Wurtzite(const Material* mat) {
	this->set_material(mat);			
	// matrix space must be allocated from derived class
	// allocate space calls virtual function, which is defined in derived class
	// so it's not available in the default constructor of the base class
	this->allocate_matrix_space();
	this->have_first_order_terms = true;	
}

KPMatrix8x8Wurtzite::~KPMatrix8x8Wurtzite() {
	this->material = 0;
}

void KPMatrix8x8Wurtzite::init_base_and_strain_matrix() {
	TDKP_ASSERT(this->material != 0,"this->material != 0");
	
	// ---------------------------------------------------------
	// wurtzite kp matrix
	// taken from phys rev b 54, 1996, chuang and chang,
	// kp method for strained wurzite semiconductors
	// its similar to J. Hader et al. model in 
	// Interband Transitions in InGaN Quantum Wells,  
	// Nitride Semiconductor Devices: Principles and Simulation
	// Joachim Piprek (Editor)
	// ISBN: 978-3-527-40667-8 
	//
	// o.k., it's the extension of chuangs 6x6 model to an 8x8 model
	// hader does a cubic approxmation, but we take here the full matrix
	//
	// the plan
	//  1. the kane matrix in chuang (5) gives the cb <-> vb stuff
	//     but warning, my definition here is different. see point 9.
	//
	//  2. the matrix elements P2 and P1 are calculated from cb c-plane
	//     and a plane effective masses (eq. (18)) in chuang
	//     this matrix elements should now give the right band bending 
	//     in the cb band
	//
	//  3. the valence band part is the same as in the 6x6 model
	//
	//  4. the hole effective mass parameters A1-A6 (A7) inserted 
	//     into the vb part need to be renormalized as they include
	//     the effect of the now explicitly treated conduction band
	//
	//  5. this renormalization is performed as in the case of 
	//     the zincblende 8x8 case
	//
	//  6. the renormalization is performed using this formula:
	//      Hij(valenceband)= Hij + 0.5 * (1/(Hii - Ec) + 1/(Hjj - Ec)) 
	//                                  * (Hi1*H1j + Hi2*H2j)
	//     cleary, the bands 1 and 2 here denote the twofold degenerate
	//     conduction band
	//
	//  7. now using the renormalized hamiltonian, the disperson should
	//     be equal as for the initial 6x6 hamiltonian. so this gives
	//     then equations for the new parameterrs
	//
	//  8. problem 1: you can not take Hii, as it differs for each band
	//     (due to the spin orbit splitting) and you end up with 
	//     different effective mass parameters for each basis function
	//     solution 1: if you say Hii - Ec = Ev +/- delta_so - Ec
	//     and as delta_so << Ev - Ec i set delta_so = 0
	//
	//  9. problem 2: i started for the renormalization with the kane
	//     hamiltonian in chuang, eq. (5). it turned out that the 
	//     renormalization needs the effective mass params to get complex
	//     as this shouldn't be i tried to get a lucky shot and took
	//     the transpose of the kane cb-vb interaction. et voila, then 
	//     you don't need complex numbers, but you get small nice
	//     good looking formulas for the renormalized parameters
	// 
	//     can we explain that? yes, we can. if i recalculate chuangs vb matrix,
        //     i simply get the transpose instead of the one he obtained. the
        //     error seems to be in eq. 29 (D21 etc) where he changes to his 
        //     new basis
	//
	// 10. the strain effects of the vb band are the ones of the 6x6 model, 
	//     while the strain effects of the cb band are purely hydrostatic,
	//     (but anisotropic).
	//     for the whole bandgap, we have dEg = a1 * ezz + a2*(exx + eyy)
	//     [phys rev b 54/19 13460 (1996)]
	//     if you check chuangs hamiltonian, you see that u1/u2/u4/u5
	//     have the same strain dependence on their diagonal and u3/u6 
	//     are split off (eq. 34).
	//     so the HH and LH band have (D1 + D3)*ezz + (D2+D4)*(exx+eyy)
	//     as their strain dependence. therefore, the conduction
	//     band deformation can be approximated by
	//     dEc = (a1 - (D1 + D3))*ezz + (a2 - (D2 + D4)) * (exx + eyy)
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
	// get material parameters
	// ------------------------------------------------------
	WurtziteEffectiveMassParams params;	
	params.A1      = hbar_2m0 * this->material->get("hole_effective_mass_A1"); 	
	params.A2      = hbar_2m0 * this->material->get("hole_effective_mass_A2");
	params.A3      = hbar_2m0 * this->material->get("hole_effective_mass_A3");
	params.A4      = hbar_2m0 * this->material->get("hole_effective_mass_A4");
	params.A5      = hbar_2m0 * this->material->get("hole_effective_mass_A5");
	params.A6      = hbar_2m0 * this->material->get("hole_effective_mass_A6");
	params.mc_zz   = this->material->get("electron_effective_mass_c_plane");
	params.mc_xxyy = this->material->get("electron_effective_mass_a_plane");		
				
	const double vb_edge  = this->material->get("valence_band_edge");
	const double cb_edge  = this->material->get("conduction_band_edge");
	const double delta_1  = this->material->get("spin_orbit_split_delta_1");
	const double delta_2  = this->material->get("spin_orbit_split_delta_2");
	const double delta_3  = this->material->get("spin_orbit_split_delta_3");
	double       A7       = 0.0;
	if(Configuration::get_instance()->get("kpmatrix_wurtzite_include_spin_orbit_interaction") == 1.0) {
		A7 = this->material->get("hole_spin_orbit_A7");	
	}		
		
	// ------------------------------------------------------
	// get strain parameters
	// ------------------------------------------------------
	double D1, D2, D3, D4, D5, D6; // valence band
	double ac_zz, ac_xxyy;
		
	if(this->material->is_set("strain_potential_D1") &&
	   this->material->is_set("strain_potential_D2") &&
	   this->material->is_set("strain_potential_D3") &&
	   this->material->is_set("strain_potential_D4") &&
	   this->material->is_set("strain_potential_D5") &&
	   this->material->is_set("strain_potential_D6") &&
	   this->material->is_set("strain_potential_ac_cc") &&
	   this->material->is_set("strain_potential_ac_aa")) {
	   	
		D1 = this->material->get("strain_potential_D1");
		D2 = this->material->get("strain_potential_D2");
		D3 = this->material->get("strain_potential_D3");
		D4 = this->material->get("strain_potential_D4");
		D5 = this->material->get("strain_potential_D5");
		D6 = this->material->get("strain_potential_D6"); 
		ac_zz   = this->material->get("strain_potential_ac_cc");
		ac_xxyy = this->material->get("strain_potential_ac_aa");
		   	
	   	this->strain_dependence_available = true;
	} else {
		this->strain_dependence_available = false;
		if(!surpress_output()) {
			Logger::get_instance()->emit(LOG_WARN, "at least one of the strain potentials is not set. strain dependence is therefore not available!");
		} 	
		D1 = D2 = D3 = D4 = D5 = D6 = ac_zz = ac_xxyy = 0.0;
	}

	// ------------------------------------------------------
	// determine optical momentum matrix elements 
	// ------------------------------------------------------
	double P1, P2;
	bool P1_user_defined = false; 
	bool P2_user_defined = false;
	if(this->material->is_set("optical_momentum_matrix_element_P1")) {
		P1 = this->material->get("optical_momentum_matrix_element_P1");
		P1_user_defined = true;
	} else {
		P1 = get_optical_momentum_matrix_element_P1();	
	}
	if(this->material->is_set("optical_momentum_matrix_element_P2")) {
		P2 = this->material->get("optical_momentum_matrix_element_P2");
		P2_user_defined = true;
	} else {
		P2 = get_optical_momentum_matrix_element_P2();	
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
	// tell the user the parameters before renormalization
	// -------------------------------------------------------
	ostringstream sout; 
	sout.setf(ios::left);
 	sout << "kp 8x8 wurtzite matrix for material " << this->get_material_name() << " uses following parameters:\n"
 		 << "elec. eff. mass c-axis   = " << params.mc_zz << "\n"
 		 << "elec. eff. mass a-plane  = " << params.mc_xxyy << "\n" 		  		 
 		 << "hole effective mass A1   = " << params.A1 / hbar_2m0<< "\n"
 		 << "hole effective mass A2   = " << params.A2 / hbar_2m0 << "\n"
 		 << "hole effective mass A3   = " << params.A3 / hbar_2m0 << "\n"
 		 << "hole effective mass A4   = " << params.A4 / hbar_2m0 << "\n"
 		 << "hole effective mass A5   = " << params.A5 / hbar_2m0 << "\n"
 		 << "hole effective mass A6   = " << params.A6 / hbar_2m0 << "\n"
 		 << "hole spin orbit     A7   = " << setw(15) << A7 << " " << (Configuration::get_instance()->get("kpmatrix_wurtzite_include_spin_orbit_interaction") == 1.0 ? "(enabled)" : "(disabled)") << "\n" 		 
 		 << "spin orbit split delta 1 = " << delta_1 << "\n"
 		 << "spin orbit split delta 2 = " << delta_2 << "\n"
 		 << "spin orbit split delta 3 = " << delta_3 << "\n"
 		 << "conduction band edge     = " << cb_edge << "\n" 		 
 		 << "valence band edge        = " << vb_edge << "\n"
 		 << "vb shift (delta offset)  = " << vb_shift << "\n" 		 
 		 << "optic. momnt. matrix P1  = " << setw(15) << P1 << " (Ep = " << P1*P1/hbar_2m0 << " eV) " << (P1_user_defined ? "(user defined)" : "(calculated from eff.mass)") << "\n"
 		 << "optic. momnt. matrix P2  = " << setw(15) << P2 << " (Ep = " << P2*P2/hbar_2m0 << " eV) " << (P2_user_defined ? "(user defined)" : "(calculated from eff.mass)") << "\n"
 		 << "operator ordering        = ";	
	if(operator_choice == splitting_symmetric) {
		sout << "symmetric splitting\n";
	} else if(operator_choice == splitting_user_defined) {
		sout << "user defined\n";
	} else {
		sout << "automatically setting A- to 0\n";	
	}

	if(this->strain_dependence_available) {
		sout << "strain potential D1      = " << D1 << "\n"
		     << "strain potential D2      = " << D2 << "\n"
		     << "strain potential D3      = " << D3 << "\n"
		     << "strain potential D4      = " << D4 << "\n"
		     << "strain potential D5      = " << D5 << "\n"
		     << "strain potential D6      = " << D6 << "\n"
		     << "strain potential ac_cc = " << ac_xxyy << "\n"
		     << "strain potential ac_aa   = " << ac_zz << "\n";		     		     
	}	 	
	

	
	// ------------------------------------------------------
	// renormalize parameters
	// ------------------------------------------------------
	if(Configuration::get_instance()->get("kpmatrix_wurtzite_8x8_skip_renormalization") != 1.0) {
		params = get_renormalized_parameters(params, P1, P2, cb_edge, vb_edge);
	} else {
		Logger::get_instance()->emit(LOG_WARN, "skipping renormalization of effective mass parameters");
	}
				
	// ------------------------------------------------------
	// determine parameter splitting
	// ------------------------------------------------------
	double A5m, A6m;
	determine_splitted_parameters(params, A5m, A6m);
	params.A6m = A6m;
	params.A5m = A5m;

	// ------------------------------------------------------
	// check non-ellipticity
	// ------------------------------------------------------
	vector<double> eigenvalues;
	determine_ellipticity_eigenvalues(params, eigenvalues);
	const double nonellipticity = get_nonellipticity_ratio(eigenvalues);	

	// ------------------------------------------------------
	// set params to short variable names
	// ------------------------------------------------------
	const double& A1  = params.A1;
	const double& A2  = params.A2;
	const double& A3  = params.A3;
	const double& A4  = params.A4;
	const double& A5  = params.A5;
	const double& A6  = params.A6;		
	const double  A5p = params.A5 - A5m;
	const double  A6p = params.A6 - A6m;
	
	const double& mc_zz   = params.mc_zz;
	const double& mc_xxyy = params.mc_xxyy;
		 		
	// -------------------------------------------------------
	// now tell the user the parameters after the renormalization
	// -------------------------------------------------------
	sout << "derived quantities:\n"
	     << "renormalized mc_eff_xxyy = " << setw(15) << mc_xxyy << " " << ( mc_xxyy <= 0.0 ? "(bad, this will be spurious)" : "(good)") << "\n"
	     << "renormalized mc_eff_zz   = " << setw(15) << mc_zz << " " << ( mc_zz <= 0.0 ? "(bad, this will be spurious)" : "(good)") << "\n" 
		 << "renormalized A1          = " << A1 / hbar_2m0 << "\n"
		 << "renormalized A2          = " << A2 / hbar_2m0 << "\n"
		 << "renormalized A3          = " << A3 / hbar_2m0 << "\n"
		 << "renormalized A4          = " << A4 / hbar_2m0 << "\n"
		 << "renormalized A5          = " << A5 / hbar_2m0 << "\n"
		 << "renormalized A6          = " << A6 / hbar_2m0 << "\n"
		 << "operator ordering: " << "    => A5+ " << A5p / hbar_2m0 << ", A5- => " << A5m / hbar_2m0 
	     << ", A6+ => " << A6p  / hbar_2m0 << ", A6m => " << A6m / hbar_2m0 << "\n"
	     << "nonellipticity           = " << setw(15) << nonellipticity << " " << get_nonellipticity_warning(nonellipticity) << "\n";
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

/*(6,6)*/ this->second_order_strainless[D_DX][D_DX][get_idx(6,6)] =   hbar_2m0 / mc_xxyy;
/*(7,7)*/ this->second_order_strainless[D_DX][D_DX][get_idx(7,7)] =   hbar_2m0 / mc_xxyy;

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

/*(6,6)*/ this->second_order_strainless[D_DY][D_DY][get_idx(6,6)] =   hbar_2m0 / mc_xxyy;
/*(7,7)*/ this->second_order_strainless[D_DY][D_DY][get_idx(7,7)] =   hbar_2m0 / mc_xxyy;
 	
 	// ------------------------------------------------------
 	// build D_DZ D_DZ
 	// ------------------------------------------------------
/*(0,0)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(0,0)] =   A1 + A3;
/*(1,1)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(1,1)] =   A1 + A3;
/*(2,2)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(2,2)] =   A1;
/*(3,3)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(3,3)] =   A1 + A3;
/*(4,4)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(4,4)] =   A1 + A3;
/*(5,5)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(5,5)] =   A1;
 	
/*(6,6)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(6,6)] =   hbar_2m0 / mc_zz;
/*(7,7)*/ this->second_order_strainless[D_DZ][D_DZ][get_idx(7,7)] =   hbar_2m0 / mc_zz;			
			
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
 	// first order D_DX
 	// ------------------------------------------------------
 	// operator ordering assumed to be as used by foreman!
 	// -> basically we make the choice that makes the operator hermitian!
/*(0,6)*/ this->first_order_strainless[op_right][D_DX][get_idx(0,6)] =   i * P2 / constants::sqrt2;
/*(1,6)*/ this->first_order_strainless[op_right][D_DX][get_idx(1,6)] = - i * P2 / constants::sqrt2; 
/*(3,7)*/ this->first_order_strainless[op_right][D_DX][get_idx(3,7)] = - i * P2 / constants::sqrt2; 
/*(4,7)*/ this->first_order_strainless[op_right][D_DX][get_idx(4,7)] =   i * P2 / constants::sqrt2;
/*(6,0)*/ this->first_order_strainless[op_left][D_DX][get_idx(6,0)]  =   i * P2 / constants::sqrt2;
/*(6,1)*/ this->first_order_strainless[op_left][D_DX][get_idx(6,1)]  = - i * P2 / constants::sqrt2; 
/*(7,3)*/ this->first_order_strainless[op_left][D_DX][get_idx(7,3)]  = - i * P2 / constants::sqrt2;
/*(7,4)*/ this->first_order_strainless[op_left][D_DX][get_idx(7,4)]  =   i * P2 / constants::sqrt2;

 	// ------------------------------------------------------
 	// first order D_DY
 	// ------------------------------------------------------
/*(0,6)*/ this->first_order_strainless[op_right][D_DY][get_idx(0,6)] = - P2 / constants::sqrt2;
/*(1,6)*/ this->first_order_strainless[op_right][D_DY][get_idx(1,6)] = - P2 / constants::sqrt2; 
/*(3,7)*/ this->first_order_strainless[op_right][D_DY][get_idx(3,7)] = - P2 / constants::sqrt2; 
/*(4,7)*/ this->first_order_strainless[op_right][D_DY][get_idx(4,7)] = - P2 / constants::sqrt2;
/*(6,0)*/ this->first_order_strainless[op_left][D_DY][get_idx(6,0)]  =   P2 / constants::sqrt2;
/*(6,1)*/ this->first_order_strainless[op_left][D_DY][get_idx(6,1)]  =   P2 / constants::sqrt2; 
/*(7,3)*/ this->first_order_strainless[op_left][D_DY][get_idx(7,3)]  =   P2 / constants::sqrt2;
/*(7,4)*/ this->first_order_strainless[op_left][D_DY][get_idx(7,4)]  =   P2 / constants::sqrt2;

 	// ------------------------------------------------------
 	// first order D_DZ
 	// ------------------------------------------------------
/*(2,6)*/ this->first_order_strainless[op_right][D_DZ][get_idx(2,6)] = - i * P1;
/*(5,7)*/ this->first_order_strainless[op_right][D_DZ][get_idx(5,7)] = - i * P1; 
/*(6,2)*/ this->first_order_strainless[op_left][D_DZ][get_idx(6,2)]  = - i * P1; 
/*(7,5)*/ this->first_order_strainless[op_left][D_DZ][get_idx(7,5)]  = - i * P1;
 	
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
/*(0,2)*/ this->first_order_strainless[op_left][D_DX][get_idx(0,2)]  += - i * A7;
/*(1,2)*/ this->first_order_strainless[op_left][D_DX][get_idx(1,2)]  += - i * A7;
/*(2,0)*/ this->first_order_strainless[op_right][D_DX][get_idx(2,0)] += - i * A7;
/*(2,1)*/ this->first_order_strainless[op_right][D_DX][get_idx(2,1)] += - i * A7;
/*(3,5)*/ this->first_order_strainless[op_left][D_DX][get_idx(3,5)]  += - i * A7;
/*(4,5)*/ this->first_order_strainless[op_left][D_DX][get_idx(4,5)]  += - i * A7;
/*(5,3)*/ this->first_order_strainless[op_right][D_DX][get_idx(5,3)] += - i * A7;
/*(5,4)*/ this->first_order_strainless[op_right][D_DX][get_idx(5,4)] += - i * A7;	

	// ------------------------------------------------------
	// first order strainless D_DY
	// ------------------------------------------------------
/*(0,2)*/ this->first_order_strainless[op_left][D_DY][get_idx(0,2)]  +=   A7;
/*(1,2)*/ this->first_order_strainless[op_left][D_DY][get_idx(1,2)]  += - A7;
/*(2,0)*/ this->first_order_strainless[op_right][D_DY][get_idx(2,0)] += - A7;
/*(2,1)*/ this->first_order_strainless[op_right][D_DY][get_idx(2,1)] +=   A7;
/*(3,5)*/ this->first_order_strainless[op_left][D_DY][get_idx(3,5)]  += - A7;
/*(4,5)*/ this->first_order_strainless[op_left][D_DY][get_idx(4,5)]  +=   A7;
/*(5,3)*/ this->first_order_strainless[op_right][D_DY][get_idx(5,3)] +=   A7;
/*(5,4)*/ this->first_order_strainless[op_right][D_DY][get_idx(5,4)] += - A7;
		
	}
 	
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

/*(6,6)*/ this->zero_order_strainless[get_idx(6,6)] = cb_edge;
/*(7,7)*/ this->zero_order_strainless[get_idx(7,7)] = cb_edge;
		
	
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
/*(6,6)*/ this->zero_order_strain_dependent[EX][EX][get_idx(6,6)] =   ac_xxyy;
/*(7,7)*/ this->zero_order_strain_dependent[EX][EX][get_idx(7,7)] =   ac_xxyy;

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
/*(6,6)*/ this->zero_order_strain_dependent[EY][EY][get_idx(6,6)] =   ac_xxyy;
/*(7,7)*/ this->zero_order_strain_dependent[EY][EY][get_idx(7,7)] =   ac_xxyy;

		// ---------------------------------------------------
		// EZZ
		// ---------------------------------------------------		
/*(0,0)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(0,0)] = D1 + D3 ;
/*(1,1)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(1,1)] = D1 + D3;
/*(2,2)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(2,2)] = D1;
/*(3,3)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(3,3)] = D1 + D3;
/*(4,4)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(4,4)] = D1 + D3;
/*(5,5)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(5,5)] = D1;
/*(6,6)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(6,6)] =   ac_zz;
/*(7,7)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(7,7)] =   ac_zz;

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

		// deprecated: chuangs stuff lead to wrong bands!!!
		// the upper terms are derived from pikus and bir, the lower ones
		// are taken from chuang ...  and they are wrong (c4 symmetry in xy plane)
		if(false) {
		// ---------------------------------------------------
		// EXY
		// ---------------------------------------------------		
/*(0,1)*/ this->zero_order_strain_dependent[EX][EY][get_idx(0,1)] =   i * 2.0 * D5;
/*(1,0)*/ this->zero_order_strain_dependent[EX][EY][get_idx(1,0)] = - i * 2.0 * D5;
/*(3,4)*/ this->zero_order_strain_dependent[EX][EY][get_idx(3,4)] = - i * 2.0 * D5;
/*(4,3)*/ this->zero_order_strain_dependent[EX][EY][get_idx(4,3)] =   i * 2.0 * D5;

		// ---------------------------------------------------
		// EXZ
		// ---------------------------------------------------		
/*(0,2)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(0,2)] = - D6;
/*(1,2)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(1,2)] =   D6;
/*(2,0)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(2,0)] = - D6;
/*(2,1)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(2,1)] =   D6;
/*(3,5)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(3,5)] =   D6;
/*(4,5)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(4,5)] = - D6;
/*(5,3)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(5,3)] = D6;
/*(5,4)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(5,4)] = - D6;


		// ---------------------------------------------------
		// EYZ
		// ---------------------------------------------------		
/*(0,2)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(0,2)] =   i * D6;
/*(1,2)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(1,2)] =   i * D6;
/*(2,0)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(2,0)] = - i * D6;
/*(2,1)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(2,1)] = - i * D6;
/*(3,5)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(3,5)] =   i * D6;
/*(4,5)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(4,5)] =   i * D6;
/*(5,3)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(5,3)] = - i * D6;
/*(5,4)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(5,4)] = - i * D6;
		}
		
	}			
						
}


void KPMatrix8x8Wurtzite::determine_splitted_parameters(const WurtziteEffectiveMassParams& splitted_params, double &A5m, double &A6m) const {

	// --------------------------------------------------
	// symmetrize if somebody is that stupid ...
	// --------------------------------------------------
	if(Configuration::get_instance()->get("kpmatrix_disable_foreman_enable_symmetrization") == 1.0) {
		A5m = splitted_params.A5 / 2.0;
		A6m = splitted_params.A6 / 2.0;
	} else {
		// -------------------------------------------------
		// read splitting from material database
		// -------------------------------------------------
		if(this->material->is_set("hole_effective_mass_A5_minus") && this->material->is_set("hole_effective_mass_A6_minus")) {
			A5m = this->material->get("hole_effective_mass_A5_minus");
			A6m = this->material->get("hole_effective_mass_A6_minus");
		} else {
			// -------------------------------------------------
			// optimal splitting 
			// -------------------------------------------------
			// although we have all the code to determine the optimal
			// splitting, we don't need it: the optimal splitting is
			// A?+ = A? and A?- = 0 ...
			// seems to be something symmetry caused ... 
			A5m = 0.0;
			A6m = 0.0;
		}
	}
	
}

double KPMatrix8x8Wurtzite::get_optical_momentum_matrix_element_P1() const {
	
	const double mz = this->material->get("electron_effective_mass_c_plane");
	const double Eg = this->material->get("conduction_band_edge") - this->material->get("valence_band_edge");
	const double d1 = this->material->get("spin_orbit_split_delta_1");
	const double d2 = this->material->get("spin_orbit_split_delta_2");
	const double d3 = this->material->get("spin_orbit_split_delta_3");
	
	TDKP_ASSERT(Eg > 0, " bandgap < 0???");
	
	double tmp =  constants::hbar_square / (2.0 * constants::m0)
	     * (1.0 / mz - 1.0)	     
	     * (
	     	 (Eg + d1 + d2) * (Eg + 2.0 * d2) 
	     	  - (2.0 * d3 * d3)
	       )
	     / (Eg + 2.0 * d2);
	TDKP_ASSERT(tmp > 0.0, "P1_square < 0.0");	     
	return sqrt(tmp);
}
double KPMatrix8x8Wurtzite::get_optical_momentum_matrix_element_P2() const {

	const double mtr = this->material->get("electron_effective_mass_a_plane");
	const double Eg  = this->material->get("conduction_band_edge") - this->material->get("valence_band_edge");
	const double d1  = this->material->get("spin_orbit_split_delta_1");
	const double d2  = this->material->get("spin_orbit_split_delta_2");
	const double d3  = this->material->get("spin_orbit_split_delta_3");
	
	TDKP_ASSERT(Eg > 0, " bandgap < 0???");
	
	double tmp = constants::hbar_square / (2.0 * constants::m0)
	     	   * (1.0 / mtr - 1.0)
	     	   * Eg
	     	   * (
	     	   		(Eg + d1 + d2) 
	     	   	  * (Eg + 2.0 * d2)
	     	   	  - (2.0 * d3 * d3)
	     	   	 )
	     	   / (
	     	   		(Eg + d1 + d2) 
	     	   	  * (Eg + d2)
	     	   	  - (d3 * d3)
	     	     );
	     	   	 	
	TDKP_ASSERT(tmp > 0.0, "P2_square < 0.0");
	
	return sqrt(tmp);
}

/** renormalize wurzite 6x6 parameters for wurtzite 8x8 model
 * 
 * this is the result of the transformation
 * Hij -> Hij + sum_{v} (HivHvj) / (Ev - Ei)
 * 
 * the renormalization just applies to A1 - A6, and is not applied
 * to the splitting in A5/A6!
 * 
 */ 
WurtziteEffectiveMassParams KPMatrix8x8Wurtzite::get_renormalized_parameters(const WurtziteEffectiveMassParams& params, const double& P1, const double& P2, const double& cb_edge, const double& vb_edge) const {

	WurtziteEffectiveMassParams ret = params;
	
	ret.A1 += P1 * P1 / (cb_edge - vb_edge);
	ret.A2 += 0.0; // does not change
	ret.A3 -= P1 * P1 / (cb_edge - vb_edge);
	ret.A4 += P2 * P2 / (2.0 * (cb_edge - vb_edge));
	ret.A5 += P2 * P2 / (2.0 * (cb_edge - vb_edge));
	ret.A6 += P1 * P2 / (constants::sqrt2 * (cb_edge - vb_edge)); 
	
	double sszz = 1.0 / ret.mc_zz - ((2.0 * constants::m0 / constants::hbar_square) * (P1 * P1 / (cb_edge - vb_edge)));	
	double ssxx	= 1.0 / ret.mc_xxyy - ((2.0 * constants::m0 / constants::hbar_square) * (P2 * P2 / (cb_edge - vb_edge))); 	
	
	TDKP_ASSERT(sszz != 0.0, "bad renormalization! choose different effective mass m_cplane");
	TDKP_ASSERT(ssxx != 0.0, "bad renormalization! choose different effective mass m_aplane");
	
	ret.mc_zz = 1.0 / sszz;
	ret.mc_xxyy = 1.0 / ssxx;
		
	return ret;
	
}

} // end of namespace

