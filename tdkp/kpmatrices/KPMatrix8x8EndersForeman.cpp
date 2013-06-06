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

#include "tdkp/kpmatrices/KPMatrix8x8EndersForeman.h"
#include "tdkp/common/Configuration.h"

namespace tdkp {

// ---------------------------------------------------------
// the enders 8x8 hamiltonian with foreman operator ordering
// ---------------------------------------------------------
const int KPMatrix8x8EndersForeman::sparsity_pattern[88] = {
	0,0,	0,1,	0,2,	0,3,
	1,0,	1,1,	1,2,	1,3,					1,6,	1,7,
	2,0,	2,1,	2,2,	2,3,			2,5,			2,7,
	3,0,	3,1,	3,2,	3,3,			3,5,	3,6,
									4,4,	4,5,	4,6,	4,7,
					5,2,	5,3,	5,4,	5,5,	5,6,	5,7,
			6,1,			6,3,	6,4,	6,5,	6,6,	6,7,
			7,1,	7,2,			7,4,	7,5,	7,6,	7,7
};

inline int KPMatrix8x8EndersForeman::get_idx(int ii, int jj) const {
	for(int ss = 0; ss < 44; ss++) {
		if(sparsity_pattern[2*ss] == ii && sparsity_pattern[2*ss+1] == jj) {
			return ss;	
		}	
	}	
	TDKP_GENERAL_EXCEPTION("invalid idx: " << ii << " and " << jj);
}

/** the diagonal elements in the kp matrix are stored at the following locations */
const int KPMatrix8x8EndersForeman::diagonal_pattern[8] = {0,5,12,19,22,29,36,43};

const int* KPMatrix8x8EndersForeman::get_sparsity_pattern(int &num) const {
	num = 44;
	return &(this->sparsity_pattern[0]);
}

const int* KPMatrix8x8EndersForeman::get_diagonal_pattern(int &num) const {
	num = 8; 
	return &(this->diagonal_pattern[0]);	
}

KPMatrix8x8EndersForeman::KPMatrix8x8EndersForeman(const Material* mat) {
	this->set_material(mat);			
	// matrix space must be allocated from derived class
	// allocate space calls virtual function, which is defined in derived class
	// so it's not available in the default constructor of the base class
	this->allocate_matrix_space();		
}

KPMatrix8x8EndersForeman::KPMatrix8x8EndersForeman() {
	this->allocate_matrix_space();
	this->material = 0;
}

KPMatrix8x8EndersForeman::~KPMatrix8x8EndersForeman() {
	material = 0;
}

bool KPMatrix8x8EndersForeman::check_solution_type_available(KPSolutionType type) const {
	switch(type) {
		case electrons:
			return true;
		case holes:
			return true;
		default:
			TDKP_GENERAL_EXCEPTION("unknown kp solution type");	
	}	
}

void KPMatrix8x8EndersForeman::init_base_and_strain_matrix() {
	
	TDKP_ASSERT(this->material != 0,"this->material != 0");
	
	// -----------------------------------------------------
	// for parameter definition see
	// phys rev b 51 16695 (1995) (enders)
	// and phys rev b 56 R12748 (1997)
	// -----------------------------------------------------	
	double lut1, lut2, lut3, lut1o, lut2o, lut3o, vb_edge, cb_edge, 
	       P,  mc_eff, Ac, hbar_2m0, delta_so, optmat,
	       bandgap;

	double M,L,Nprime, Nm, Np;
	const complex<double> i(0,1);
	hbar_2m0 = constants::hbar_square / (2.0 * constants::m0);

	lut1o    = this->material->get("luttinger_parameter_1");
	lut2o    = this->material->get("luttinger_parameter_2");
	lut3o    = this->material->get("luttinger_parameter_3");
	mc_eff   = this->material->get("electron_effective_mass"); 						
	vb_edge  = this->material->get("valence_band_edge"); 
	delta_so = this->material->get("spin_orbit_splitting");
	optmat   = this->material->get("optical_matrix_element");
	bandgap  = this->material->get("bandgap");
	cb_edge  = vb_edge + bandgap; 
	P        = sqrt(optmat * hbar_2m0);

	// -------------------------------------------------------------
	// adjust luttinger parameters
	// lut parameters are defined for Gamma15 states (hole). as we 
	// now include the cb band into the calculation, we have to subtract
	// the effect of the cb band (Gamma1) as we will handle it directly
	// in the calculation.
	// there is some dispute in the literature what the exact formula is.
	// bahder (phys rev. b. vol 41 1990 11992) uses 
	// lut_i = lut_i - (alpha) * Ep / (3Eg+delta_so), where alpha = (1 for i = 1, 1/2 else)
	// but other e.g. Guy Fishman Phys. Rev. B. Vol. 63, 235302 says
	// one should use the ones obtained by Pidgeon and Brown (1966, Phys. Rev. 146).
	// further foreman uses P&B's, and therefore as following foreman
	// seems to be generally a good idea, i will do that now 
	// -------------------------------------------------------------
	if(Configuration::get_instance()->get("kpmatrix_zincblende_8x8_skip_renormalization") == 0.0) {
		lut1 = lut1o - optmat / (3.0 * bandgap); 
		lut2 = lut2o - optmat / (6.0 * bandgap);
		lut3 = lut3o - optmat / (6.0 * bandgap);
		// foreman eq. (6)
		Ac = hbar_2m0 * (1.0 / mc_eff - optmat / 3.0 * (2.0 / bandgap + 1.0 / (bandgap + delta_so)));
	} else {
		if(!surpress_output()) {
			Logger::get_instance()->emit(LOG_WARN, "user requested to disable parameter renormalization");
		}		
		Ac = hbar_2m0 * (1.0 / mc_eff);
		lut1 = lut1o; 
		lut2 = lut2o;
		lut3 = lut3o;		
	} 

	
	// -------------------------------------------------------------
	// but as you can expect, such a correction can be really critical
	// what may happen is that the effective mass of some band
	// turns into the wrong direction and therefore instead
	// of an electron we get a hole and this will produce spurious modes
	// as the electron would as a hole not see a quantized region but a 
	// barrier ... and there you go with exponential decay.
	// -------------------------------------------------------------
	double optmat_critical_for_L       = bandgap * (lut1o + 4.0 * lut2o);
	double optmat_critical_for_Nprime  = bandgap * 6.0 * lut3o;
	double optmat_critical_for_Ac      = 3.0 / mc_eff / (2.0 / bandgap + 1.0 / (bandgap + delta_so));
	
		                                 
	double bandgap_critical_for_L      = optmat / (lut1o + 4.0 * lut2o);
	double bandgap_critical_for_Nprime = optmat / (6.0 * lut3o);
	double bandgap_critical_for_Ac     = 1.0 / 6.0 * (optmat * mc_eff - 3.0 * delta_so
	                                     + sqrt(optmat*optmat*mc_eff*mc_eff + 18.0 * optmat
	                                           * mc_eff * delta_so + 9.0 * delta_so*delta_so));
	

	// -------------------------------------------------------------
	// foreman transformation
	// see phys rev b. vol 56, nr 20, 12748 eq. (12) and (13)
	// -------------------------------------------------------------			
	L        = - hbar_2m0 * (lut1 + 4.0 * lut2);
	M        = - hbar_2m0 * (lut1 - 2.0 * lut2);
	Nprime   = - hbar_2m0 * 6.0 * lut3;
	Nm       = M - hbar_2m0;
	// ----------------------------------------------------------
	// overwrite Nm if its set in the material file
	// ----------------------------------------------------------
	if(this->material->is_set("foremans_kane_parameter_N_minus")) {
		Nm = hbar_2m0 * this->material->get("foremans_kane_parameter_N_minus");
	}	
	Np       = Nprime - Nm;	 
	
		
	// ---------------------------------------------------------
	// strain potentials
	// ---------------------------------------------------------
	double ac,av,b,d;
	ac = av = b = d = 0.0;
	if(this->material->is_set("strain_potential_ac") &&
 	   this->material->is_set("strain_potential_av") &&
	   this->material->is_set("strain_potential_b")  &&
       this->material->is_set("strain_potential_d")) {
		// ac = hydrostatic strain potential of conduction band
		ac 		 = this->material->get("strain_potential_ac");       	
		// av = hydrostatic strain potential of valence band
		av 		 = this->material->get("strain_potential_av");
		// b,d = shear deformation potential
		b        = this->material->get("strain_potential_b");
		d		 = this->material->get("strain_potential_d");
		this->strain_dependence_available = true;
		
		// --------------------------------------------------------
		// warn on strange potentials
		// --------------------------------------------------------
		if(ac > 0.0 || av < 0.0 || b > 0.0 || d > 0.0) {
			ostringstream pwarn;
			pwarn << "your valence band deformation potentials have possibly a wrong sign. "
			      << "for usual deformation potentials, ac < 0, av > 0, b < 0.0 and d < 0.0 "
			      << "please use a_gap = ac - av for the vb deformation potential.";
			Logger::get_instance()->emit(LOG_WARN, pwarn.str());
		}
				
	} else {
		this->strain_dependence_available = false;
		ac = av = b = d = 0.0e0;
		if(!surpress_output()) {
			string stout("at least one of the strain parameters ac, av, b or d is not defined. strain dependence is therefore not available!");
			Logger::get_instance()->emit(LOG_WARN, stout);
		}
	}
	// -----------------------------------------------------
	// deformation potentials for the enders matrix
	// see enders prb 1995 (4c') and bahder prb 1990 eq. (30)
	// ----------------------------------------------------- 
	const double n = sqrt(3.0) * d;
	const double l = av + 2.0 * b;
	const double m = av - b;
	
	
	// ----------------------------------------------------------
	// analyze non-ellipticity
	// see veprek et al, prb ;-)
	// ----------------------------------------------------------
	double ellipticity_values[] = {M - Nm, M + Nm, L - Np, L+2.0*Np};
	double multiplicity[] = {3.0, 3.0, 2.0, 1.0}; // vielfachheit
	double positive_values = 0.0;
	double negative_values = 0.0;
	for(unsigned int ii = 0; ii < 4; ii++) {
		if(ellipticity_values[ii] < 0.0) {
			negative_values += multiplicity[ii] * ellipticity_values[ii];	
		} else {
			positive_values += multiplicity[ii] * ellipticity_values[ii];
		}	
	}
	string degree_of_nonellipticity_string;
	double degree_of_nonellipticity = negative_values == 0.0 ? 1.0 : tdkp_math::abs(positive_values / negative_values);
	degree_of_nonellipticity_string = get_nonellipticity_warning(degree_of_nonellipticity);
		
	ostringstream sout; 
	sout.setf(ios::left);
 	sout << "kp 8x8 matrix for material " << this->get_material_name() << " uses following parameters:\n";
 	sout << "lut1         = " << lut1o			    << "\n";
 	sout << "lut2         = " << lut2o    			<< "\n" ;
 	sout << "lut3         = " << lut3o    			<< "\n";
 	sout << "mc_eff       = " << mc_eff             << "\n";
 	sout << "bandgap      = " << bandgap            << "\n";
 	sout << "cb_edge      = " << cb_edge            << " (= vb_edge + bandgap)\n";
 	sout << "vb_edge      = " << vb_edge            << "\n";
 	sout << "delta_so     = " << delta_so           << "\n";
 	sout << "opt.mat.el.  = " << optmat             << "\n";
	sout << "energy shift = " << this->energy_shift << "\n";
	
	
	 	
	if(this->strain_dependence_available) {
		sout << "\n"
		     << "strain dependency is available:\n"
		     << "strain ac    = " << ac << "\n"
		     << "strain av    = " << av << "\n"
		     << "strain b     = " << b  << "\n"
		     << "strain d     = " << d  << "\n";
	} 	 	              
 	sout << " derived quantities:\n";
 	sout << "Ac           = " << setw(10) << Ac / hbar_2m0      << (Ac     > 0.0 ? " (good)" : "(unstable, should be > 0)") << "\n";
 	sout << "L            = " << setw(10) << L / hbar_2m0       << (L      < 0.0 ? " (good)" : "(unstable, should be < 0)") << "\n";
 	sout << "M            = " << setw(10) << M / hbar_2m0       << (M      < 0.0 ? " (good)" : "(unstable, should be < 0)") << "\n";
 	sout << "N            = " << setw(10) << Nprime / hbar_2m0  << (Nprime < 0.0 ? " (good)" : "(unstable, should be < 0)") << "\n";
 	sout << "N+           = " << setw(10) << Np / hbar_2m0		<< " ( = N - N-)\n";
 	sout << "N-           = " << setw(10) << Nm / hbar_2m0      << (this->material->is_set("foremans_kane_parameter_N_minus") ? " (read from material file)" : " ( = M - h^2/2m)") << "\n";
  	sout << "P            = " << P                  << "\n"; 
 	sout << "hbar2_2m0    = " << hbar_2m0           << "\n";
 	sout << "non-elliptic = " << setw(10) << degree_of_nonellipticity << degree_of_nonellipticity_string << "\n"; 

	if(!surpress_output()) {
 		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	} 		
 	
 	// -----------------------------------------------------------
 	// the less complicated error message
 	// -----------------------------------------------------------
 	if(Ac <= 0.0 || M >= 0.0 || Nprime >= 0.0 || L >= 0.0) {
 		sout.str("");
 		sout << "you have a problem - your parameters will probably produce spurious\n"
 		     << "modes. the results you will obtain will most probably be wrong and\n"
 		     << "unphysical. usually the reason is due to the unknown splitting of\n"
 		     << "the parameter N into N+ and N-, due to a too big optical matrix\n"
 		     << "or due to a too small bandgap. the critical values for bandgap and\n"
 		     << "matrix element are (where L,Nprime,Ac change sign):\n" 
 		     << "         " << setw(18) << " matrix element "    << " " << setw(18) << " bandgap" << "\n"
 		     << "L      = " << setw(18) << optmat_critical_for_L << " " << setw(18) << bandgap_critical_for_L << "\n"   		     
 		     << "Nprime = " << setw(18) << optmat_critical_for_Nprime << " " << setw(18) << bandgap_critical_for_Nprime << "\n"
 		     << "Ac     = " << setw(18) << optmat_critical_for_Ac << " " << setw(18) << bandgap_critical_for_Ac; 
		Logger::get_instance()->emit(LOG_WARN, sout.str()); 		     
 	}  

	this->initialized = true;


	if(Configuration::get_instance()->get("kpmatrix_disable_foreman_enable_symmetrization") == 1) {
		Logger::get_instance()->emit(LOG_WARN, "you requested to disable the burt/foreman operator ordering and use symmetrization. therefore setting Np = Nm = 0.5 N'");
		Np = Nm = Nprime / 2.0;
	} 
		

	// --------------------------------------------------------------------------
	// diagonal kp terms d/dx d/dx
	// --------------------------------------------------------------------------									
/*(0,0)*/	this->second_order_strainless[D_DX][D_DX][get_idx(0,0)] = Ac;
/*(1,1)*/	this->second_order_strainless[D_DX][D_DX][get_idx(1,1)] = L;
/*(2,2)*/	this->second_order_strainless[D_DX][D_DX][get_idx(2,2)] = M;
/*(3,3)*/	this->second_order_strainless[D_DX][D_DX][get_idx(3,3)] = M;
/*(4,4)*/	this->second_order_strainless[D_DX][D_DX][get_idx(4,4)] = Ac;
/*(5,5)*/	this->second_order_strainless[D_DX][D_DX][get_idx(5,5)] = L;
/*(6,6)*/	this->second_order_strainless[D_DX][D_DX][get_idx(6,6)] = M;
/*(7,7)*/	this->second_order_strainless[D_DX][D_DX][get_idx(7,7)] = M;

	// --------------------------------------------------------------------------
	// diagonal kp terms d/dy d/dy
	// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DY][D_DY][get_idx(0,0)]  = Ac;
/*(1,1)*/	this->second_order_strainless[D_DY][D_DY][get_idx(1,1)]  = M;
/*(2,2)*/	this->second_order_strainless[D_DY][D_DY][get_idx(2,2)]  = L;
/*(3,3)*/	this->second_order_strainless[D_DY][D_DY][get_idx(3,3)]  = M;
/*(4,4)*/	this->second_order_strainless[D_DY][D_DY][get_idx(4,4)]  = Ac;
/*(5,5)*/	this->second_order_strainless[D_DY][D_DY][get_idx(5,5)]  = M;
/*(6,6)*/	this->second_order_strainless[D_DY][D_DY][get_idx(6,6)]  = L;
/*(7,7)*/	this->second_order_strainless[D_DY][D_DY][get_idx(7,7)]  = M;
	
	// --------------------------------------------------------------------------
	// diagonal kp terms d/dz d/dz
	// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(0,0)]   = Ac;
/*(1,1)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(1,1)]   = M;
/*(2,2)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(2,2)]  = M;
/*(3,3)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(3,3)]  = L;
/*(4,4)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(4,4)]  = Ac;
/*(5,5)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(5,5)]  = M;
/*(6,6)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(6,6)]  = M;
/*(7,7)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(7,7)]  = L;	
	
	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dx d/dy
	// --------------------------------------------------------------------------
/*(1,2)*/	this->second_order_strainless[D_DX][D_DY][get_idx(1,2)]   = Np;
/*(2,1)*/	this->second_order_strainless[D_DX][D_DY][get_idx(2,1)]  = Nm;
/*(5,6)*/	this->second_order_strainless[D_DX][D_DY][get_idx(5,6)]  = Np;
/*(6,5)*/	this->second_order_strainless[D_DX][D_DY][get_idx(6,5)]  = Nm;
	
	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dy d/dx
	// --------------------------------------------------------------------------
/*(1,2)*/	this->second_order_strainless[D_DY][D_DX][get_idx(1,2)]   = Nm;
/*(2,1)*/	this->second_order_strainless[D_DY][D_DX][get_idx(2,1)]  = Np;
/*(5,6)*/	this->second_order_strainless[D_DY][D_DX][get_idx(5,6)]  = Nm;
/*(6,5)*/	this->second_order_strainless[D_DY][D_DX][get_idx(6,5)]  = Np;

	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dx d/dz
	// --------------------------------------------------------------------------	
/*(1,3)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(1,3)]   = Np;
/*(3,1)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(3,1)]  = Nm;
/*(5,7)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(5,7)]  = Np;
/*(7,5)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(7,5)]  = Nm;

	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dz d/dx
	// --------------------------------------------------------------------------
/*(1,3)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(1,3)]   = Nm;
/*(3,1)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(3,1)]  = Np;
/*(5,7)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(5,7)]  = Nm;
/*(7,5)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(7,5)]  = Np;

	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dy d/dz
	// --------------------------------------------------------------------------
/*(2,3)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(2,3)]  = Np;
/*(3,2)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(3,2)]  = Nm;
/*(6,7)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(6,7)]  = Np;
/*(7,6)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(7,6)]  = Nm;

	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dz d/dy
	// --------------------------------------------------------------------------
/*(2,3)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(2,3)]  = Nm;
/*(3,2)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(3,2)]  = Np;
/*(6,7)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(6,7)]  = Nm;
/*(7,6)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(7,6)]  = Np;

	// -------------------------------------------------------------------------
	// first order d/dx
	// -------------------------------------------------------------------------
/*(0,1)*/	this->first_order_strainless[op_left][D_DX][get_idx(0,1)]   =   P;
/*(1,0)*/	this->first_order_strainless[op_right][D_DX][get_idx(1,0)]  = - P;
/*(4,5)*/	this->first_order_strainless[op_left][D_DX][get_idx(4,5)]   =   P;
/*(5,4)*/	this->first_order_strainless[op_right][D_DX][get_idx(5,4)]  = - P;

	// -------------------------------------------------------------------------
	// first order d/dy
	// -------------------------------------------------------------------------
/*(0,2)*/	this->first_order_strainless[op_left][D_DY][get_idx(0,2)]   =   P;
/*(2,0)*/	this->first_order_strainless[op_right][D_DY][get_idx(2,0)]  = - P;
/*(4,6)*/	this->first_order_strainless[op_left][D_DY][get_idx(4,6)]   =   P;
/*(6,4)*/	this->first_order_strainless[op_right][D_DY][get_idx(6,4)]  = - P;

	// -------------------------------------------------------------------------
	// first order d/dz
	// -------------------------------------------------------------------------
/*(0,3)*/	this->first_order_strainless[op_left][D_DZ][get_idx(0,3)]   =   P;
/*(3,0)*/	this->first_order_strainless[op_right][D_DZ][get_idx(3,0)]  = - P;
/*(4,7)*/	this->first_order_strainless[op_left][D_DZ][get_idx(4,7)]   =   P;
/*(7,4)*/	this->first_order_strainless[op_right][D_DZ][get_idx(7,4)]  = - P;

	// -------------------------------------------------------------------------
	// zero order terms (bandedge and spin orbit split off stuff)
	// -------------------------------------------------------------------------
/*(0,0)*/	this->zero_order_strainless[get_idx(0,0)]   = cb_edge;
/*(1,1)*/	this->zero_order_strainless[get_idx(1,1)]   = vb_edge - delta_so / 3.0;
/*(1,2)*/	this->zero_order_strainless[get_idx(1,2)]   = - i * delta_so / 3.0;
/*(1,7)*/	this->zero_order_strainless[get_idx(1,7)]   =       delta_so / 3.0;
/*(2,1)*/	this->zero_order_strainless[get_idx(2,1)]  =   i * delta_so / 3.0;
/*(2,2)*/	this->zero_order_strainless[get_idx(2,2)]  = vb_edge - delta_so / 3.0;
/*(2,7)*/	this->zero_order_strainless[get_idx(2,7)]  = - i * delta_so / 3.0;
/*(3,3)*/	this->zero_order_strainless[get_idx(3,3)]  = vb_edge - delta_so / 3.0;
/*(3,5)*/	this->zero_order_strainless[get_idx(3,5)]  =     - delta_so / 3.0;
/*(3,6)*/	this->zero_order_strainless[get_idx(3,6)]  =   i * delta_so / 3.0;
/*(4,4)*/	this->zero_order_strainless[get_idx(4,4)]  = cb_edge;
/*(5,3)*/	this->zero_order_strainless[get_idx(5,3)]  =     - delta_so / 3.0;
/*(5,5)*/	this->zero_order_strainless[get_idx(5,5)]  = vb_edge - delta_so / 3.0;
/*(5,6)*/	this->zero_order_strainless[get_idx(5,6)]  =   i * delta_so / 3.0;
/*(6,3)*/	this->zero_order_strainless[get_idx(6,3)]  = - i * delta_so / 3.0;
/*(6,5)*/	this->zero_order_strainless[get_idx(6,5)]  = - i * delta_so / 3.0;
/*(6,6)*/	this->zero_order_strainless[get_idx(6,6)]  = vb_edge - delta_so / 3.0;
/*(7,1)*/	this->zero_order_strainless[get_idx(7,1)]  =       delta_so / 3.0;
/*(7,2)*/	this->zero_order_strainless[get_idx(7,2)]  =   i * delta_so / 3.0;
/*(7,7)*/	this->zero_order_strainless[get_idx(7,7)]  = vb_edge - delta_so / 3.0;


	// ------------------------------------------------------------------
	// zero order strain terms (deformation potentials)
	// see enders eq. 4c
	// there are only zero order terms, strain dependent effective
	// masses are obtained via k = (1-e)k, 
	// as according to enders, H(k,e) = H((1-e)k) + D(e)
	// and this here is D(e)
	// ------------------------------------------------------------------
	if(this->strain_dependence_available) {
		const short EX = 0; 
		const short EY = 1;
		const short EZ = 2;	
		
		// --------------------------------------------------
		// okay, i agree, its extremely stupid to name the deformation
		// potentials with small caps and the kane params with big caps
		// but anyway, blame enders. im just copying him ;-)
		// --------------------------------------------------
		
		// --------------------------------------------------
		// EXX
		// --------------------------------------------------
/*(0,0)*/ this->zero_order_strain_dependent[EX][EX][get_idx(0,0)] = ac;
/*(1,1)*/ this->zero_order_strain_dependent[EX][EX][get_idx(1,1)] = l;
/*(2,2)*/ this->zero_order_strain_dependent[EX][EX][get_idx(2,2)] = m;
/*(3,3)*/ this->zero_order_strain_dependent[EX][EX][get_idx(3,3)] = m;
/*(4,4)*/ this->zero_order_strain_dependent[EX][EX][get_idx(4,4)] = ac;
/*(5,5)*/ this->zero_order_strain_dependent[EX][EX][get_idx(5,5)] = l;
/*(6,6)*/ this->zero_order_strain_dependent[EX][EX][get_idx(6,6)] = m;
/*(7,7)*/ this->zero_order_strain_dependent[EX][EX][get_idx(7,7)] = m;

		// --------------------------------------------------
		// EYY
		// --------------------------------------------------
/*(0,0)*/ this->zero_order_strain_dependent[EY][EY][get_idx(0,0)]  = ac;
/*(1,1)*/ this->zero_order_strain_dependent[EY][EY][get_idx(1,1)]  = m;
/*(2,2)*/ this->zero_order_strain_dependent[EY][EY][get_idx(2,2)] = l;
/*(3,3)*/ this->zero_order_strain_dependent[EY][EY][get_idx(3,3)] = m;
/*(4,4)*/ this->zero_order_strain_dependent[EY][EY][get_idx(4,4)] = ac;
/*(5,5)*/ this->zero_order_strain_dependent[EY][EY][get_idx(5,5)] = m;
/*(6,6)*/ this->zero_order_strain_dependent[EY][EY][get_idx(6,6)] = l;
/*(7,7)*/ this->zero_order_strain_dependent[EY][EY][get_idx(7,7)] = m;	

		// --------------------------------------------------
		// EZZ
		// --------------------------------------------------
/*(0,0)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(0,0)] = ac;
/*(1,1)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(1,1)] = m;
/*(2,2)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(2,2)] = m;
/*(3,3)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(3,3)] = l;
/*(4,4)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(4,4)] = ac;
/*(5,5)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(5,5)] = m;
/*(6,6)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(6,6)] = m;
/*(7,7)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(7,7)] = l;
	
		// --------------------------------------------------
		// EXY
		// --------------------------------------------------		
/*(1,2)*/ this->zero_order_strain_dependent[EX][EY][get_idx(1,2)] = n;
/*(5,6)*/ this->zero_order_strain_dependent[EX][EY][get_idx(5,6)] = n;

		// --------------------------------------------------
		// EYX
		// --------------------------------------------------		
/*(2,1)*/ this->zero_order_strain_dependent[EY][EX][get_idx(2,1)] = n;
/*(6,5)*/ this->zero_order_strain_dependent[EY][EX][get_idx(6,5)] = n;

		// --------------------------------------------------
		// EXZ
		// --------------------------------------------------		
/*(1,3)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(1,3)] = n;
/*(5,7)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(5,7)] = n;

		// --------------------------------------------------
		// EZX
		// --------------------------------------------------		
/*(3,1)*/ this->zero_order_strain_dependent[EZ][EX][get_idx(3,1)] = n;
/*(7,5)*/ this->zero_order_strain_dependent[EZ][EX][get_idx(7,5)] = n;

		// --------------------------------------------------
		// EYZ
		// --------------------------------------------------		
/*(2,3)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(2,3)] = n;
/*(6,7)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(6,7)] = n;

		// --------------------------------------------------
		// EZY
		// --------------------------------------------------		
/*(3,2)*/ this->zero_order_strain_dependent[EZ][EY][get_idx(3,2)] = n;
/*(7,6)*/ this->zero_order_strain_dependent[EZ][EY][get_idx(7,6)] = n;
	

		// -------------------------------------------------------
		// spin orbit strain dependence
		// there is also a strain dependence of the spin orbit
		// split off terms, given by bahder 1990, (9c)
		// this one here is the implementation of the corresponding terms
		// 
		// the terms are derived in the following way:
		// 1. the enders basis is S, X, Y, Z plus spin 
		// 2. the bahders basis is the angular momentum basis
		// 3. the transformation U(Enders) = Bahders can be built from
		//    the bahders basis defintion (appendix A3)
		// 4. using the transformation matrix, one should obtain 
		//    in bahders basis, Hdelta (4b) as a diagonal matrix
		//    -> you won't get it, there are offdiagonal elements
		//    but if you take U*Conjugate[HSOEnders]*Inverse[U], then
		//    you get the diagonal (and vice versa)
		// 5. therefore to be consistent, the strain dependent terms
		//    calculate here are 
		//    Conjugate[Inverse[U]*<ui|DSO(1) + DSO(2)|uj>*U]
		// 6. ATTENTIONE: there was an erratum in 1992 in bahders
		//    paper where he corrected the strain dependent split off terms   
		// -------------------------------------------------------		
		if(Configuration::get_instance()->get("kpmatrix_include_spin_orbit_strain_dependence") == 1.0) {
			const double dso3 = delta_so / 3.0;
			const short EX = 0; 
			const short EY = 1;
			const short EZ = 2;			
			// --------------------------------------------------
			// EXX
			// --------------------------------------------------						
/*(1,2)*/ this->zero_order_strain_dependent[EX][EX][get_idx(1,2)] += dso3 * (i);
/*(1,7)*/ this->zero_order_strain_dependent[EX][EX][get_idx(1,7)] += dso3 * (-1);
/*(2,1)*/ this->zero_order_strain_dependent[EX][EX][get_idx(2,1)] += dso3 * (-i);
/*(3,5)*/ this->zero_order_strain_dependent[EX][EX][get_idx(3,5)] += dso3 * (1);
/*(5,3)*/ this->zero_order_strain_dependent[EX][EX][get_idx(5,3)] += dso3 * (1);
/*(5,6)*/ this->zero_order_strain_dependent[EX][EX][get_idx(5,6)] += dso3 * (-i);
/*(6,5)*/ this->zero_order_strain_dependent[EX][EX][get_idx(6,5)] += dso3 * (i);
/*(7,1)*/ this->zero_order_strain_dependent[EX][EX][get_idx(7,1)] += dso3 * (-1);

			// --------------------------------------------------
			// EYY
			// --------------------------------------------------						
/*(1,2)*/ this->zero_order_strain_dependent[EY][EY][get_idx(1,2)] += dso3 * (i);
/*(2,1)*/ this->zero_order_strain_dependent[EY][EY][get_idx(2,1)] += dso3 * (-i);
/*(2,7)*/ this->zero_order_strain_dependent[EY][EY][get_idx(2,7)] += dso3 * (i);
/*(3,6)*/ this->zero_order_strain_dependent[EY][EY][get_idx(3,6)] += dso3 * (-i);
/*(5,6)*/ this->zero_order_strain_dependent[EY][EY][get_idx(5,6)] += dso3 * (-i);
/*(6,3)*/ this->zero_order_strain_dependent[EY][EY][get_idx(6,3)] += dso3 * (i);
/*(6,5)*/ this->zero_order_strain_dependent[EY][EY][get_idx(6,5)] += dso3 * (i);
/*(7,2)*/ this->zero_order_strain_dependent[EY][EY][get_idx(7,2)] += dso3 * (-i);


			// --------------------------------------------------
			// EZZ
			// --------------------------------------------------						
/*(1,7)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(1,7)] += dso3 * (-1);
/*(2,7)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(2,7)] += dso3 * (i);
/*(3,5)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(3,5)] += dso3 * (1);
/*(3,6)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(3,6)] += dso3 * (-i);
/*(5,3)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(5,3)] += dso3 * (1);
/*(6,3)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(6,3)] += dso3 * (i);
/*(7,1)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(7,1)] += dso3 * (-1);
/*(7,2)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(7,2)] += dso3 * (-i);


			// --------------------------------------------------
			// EXY
			// --------------------------------------------------						
/*(1,7)*/ this->zero_order_strain_dependent[EX][EY][get_idx(1,7)] += dso3 * (-i);
/*(2,7)*/ this->zero_order_strain_dependent[EX][EY][get_idx(2,7)] += dso3 * (1);
/*(3,5)*/ this->zero_order_strain_dependent[EX][EY][get_idx(3,5)] += dso3 * (i);
/*(3,6)*/ this->zero_order_strain_dependent[EX][EY][get_idx(3,6)] += dso3 * (-1);
/*(5,3)*/ this->zero_order_strain_dependent[EX][EY][get_idx(5,3)] += dso3 * (-i);
/*(6,3)*/ this->zero_order_strain_dependent[EX][EY][get_idx(6,3)] += dso3 * (-1);
/*(7,1)*/ this->zero_order_strain_dependent[EX][EY][get_idx(7,1)] += dso3 * (i);
/*(7,2)*/ this->zero_order_strain_dependent[EX][EY][get_idx(7,2)] += dso3 * (1);


			// --------------------------------------------------
			// EXZ
			// --------------------------------------------------						
/*(1,6)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(1,6)] += dso3 * (-i);
/*(2,3)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(2,3)] += dso3 * (-i);
/*(2,5)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(2,5)] += dso3 * (i);
/*(3,2)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(3,2)] += dso3 * (i);
/*(5,2)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(5,2)] += dso3 * (-i);
/*(6,1)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(6,1)] += dso3 * (i);
/*(6,7)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(6,7)] += dso3 * (i);
/*(7,6)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(7,6)] += dso3 * (-i);


			// --------------------------------------------------
			// EYZ
			// --------------------------------------------------						
/*(1,3)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(1,3)] += dso3 * (-i);
/*(1,6)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(1,6)] += dso3 * (1);
/*(2,5)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(2,5)] += dso3 * (-1);
/*(3,1)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(3,1)] += dso3 * (i);
/*(5,2)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(5,2)] += dso3 * (-1);
/*(5,7)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(5,7)] += dso3 * (i);
/*(6,1)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(6,1)] += dso3 * (1);
/*(7,5)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(7,5)] += dso3 * (-i);



			 
		}		
		
	} // end zero order strain dependent
	
		
}												

} // end of namespace
