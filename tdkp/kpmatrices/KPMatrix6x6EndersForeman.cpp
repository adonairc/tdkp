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

#include "tdkp/kpmatrices/KPMatrix6x6EndersForeman.h"
#include "tdkp/common/Configuration.h"

namespace tdkp {

// ------------------------------------------------------
// kp6x6 matrix derived from the enders 8x8 hamiltonian
// by assuming P -> 0, Eg -> infinity
// using foremans operator ordering.
// ------------------------------------------------------
const int KPMatrix6x6EndersForeman::sparsity_pattern[60] = {
	0,0,	0,1,	0,2,			0,4,	0,5,
	1,0,	1,1,	1,2,	1,3,			1,5,
	2,0,	2,1,	2,2,	2,3,	2,4,	
			3,1,	3,2,	3,3,	3,4,	3,5,
	4,0,			4,2,	4,3,	4,4,	4,5,
	5,0,	5,1,			5,3,	5,4,	5,5					
};

inline int KPMatrix6x6EndersForeman::get_idx(int ii, int jj) const {
	for(int ss = 0; ss < 30; ss++) {
		if(sparsity_pattern[2*ss] == ii && sparsity_pattern[2*ss+1] == jj) {
			return ss;	
		}	
	}	
	TDKP_GENERAL_EXCEPTION("invalid idx: " << ii << " and " << jj);
}

/** the diagonal elements in the kp matrix are stored at the following locations */
const int KPMatrix6x6EndersForeman::diagonal_pattern[6] = {0,6,12,17,23,29};

const int* KPMatrix6x6EndersForeman::get_sparsity_pattern(int &num) const {
	num = 30;
	return &(this->sparsity_pattern[0]);
}

const int* KPMatrix6x6EndersForeman::get_diagonal_pattern(int &num) const {
	num = 6; 
	return &(this->diagonal_pattern[0]);	
}

KPMatrix6x6EndersForeman::KPMatrix6x6EndersForeman(const Material* mat) {
	this->set_material(mat);	
	// let kp matrix assembler in base class know that we dont have first order terms
	this->have_first_order_terms = false;		
	// matrix space must be allocated from derived class
	// allocate space calls virtual function, which is defined in derived class
	// so it's not available in the default constructor of the base class
	this->allocate_matrix_space();		
}

KPMatrix6x6EndersForeman::KPMatrix6x6EndersForeman() {	
	this->allocate_matrix_space();
	this->material = 0;
	// let kp matrix assembler in base class know that we dont have first order terms
	this->have_first_order_terms = false;	
}

KPMatrix6x6EndersForeman::~KPMatrix6x6EndersForeman() {
	material = 0;
}

bool KPMatrix6x6EndersForeman::check_solution_type_available(KPSolutionType type) const {
	switch(type) {
		case electrons:
			return false;
		case holes:
			return true;
		default:
			TDKP_GENERAL_EXCEPTION("unknown kp solution type");	
	}	
}

void KPMatrix6x6EndersForeman::init_base_and_strain_matrix() {
	
	TDKP_ASSERT(this->material != 0,"this->material != 0");
	
	double lut1, lut2, lut3, vb_edge, hbar_2m0, delta_so;
	double M,L,Nprime, Nm, Np;				
	complex<double> i(0,1);
	hbar_2m0 = constants::hbar_square / (2.0 * constants::m0);

	lut1     = hbar_2m0 * this->material->get("luttinger_parameter_1");
	lut2     = hbar_2m0 * this->material->get("luttinger_parameter_2");
	lut3     = hbar_2m0 * this->material->get("luttinger_parameter_3");						
	vb_edge  = this->material->get("valence_band_edge"); 
	delta_so = this->material->get("spin_orbit_splitting"); 

	// -------------------------------------------------------------
	// foreman transformation
	// see phys rev b. vol 56, nr 20, 12748 eq. (12) and (13)
	// -------------------------------------------------------------		
	L        = - (lut1 + 4.0 * lut2);
	M        = - (lut1 - 2.0 * lut2);
	Nprime   = - 6.0 * lut3;
	Nm       = M - hbar_2m0;	
	Np       = Nprime - Nm;

	// ---------------------------------------------------------
	// strain potentials
	// ---------------------------------------------------------
	double av,b,d;
	av = b = d = 0.0;
	if(this->material->is_set("strain_potential_av") &&
	   this->material->is_set("strain_potential_b")  &&
       this->material->is_set("strain_potential_d")) {
		// av = hydrostatic strain potential of valence band
		av 		 = this->material->get("strain_potential_av");
		// b,d = shear deformation potential
		b        = this->material->get("strain_potential_b");
		d		 = this->material->get("strain_potential_d");
		this->strain_dependence_available = true;
		
		// --------------------------------------------------------
		// warn on strange potentials
		// --------------------------------------------------------
		if(av < 0.0 || b > 0.0 || d > 0.0) {
			ostringstream pwarn;
			pwarn << "your valence band deformation potentials have possibly a wrong sign. "
			      << "for usual deformation potentials, av > 0, b < 0.0 and d < 0.0 "
			      << "please use a_gap = ac - av for the vb deformation potential.";
			Logger::get_instance()->emit(LOG_WARN, pwarn.str());
		}
				
	} else {
		this->strain_dependence_available = false;
		av = b = d = 0.0e0;
		if(!surpress_output()) {
			Logger::get_instance()->emit(LOG_WARN, "at least one of the strain parameters av,b or d is not defined. strain dependence is therefore not available!");
		}
	}
	// -----------------------------------------------------
	// deformation potentials for the enders matrix
	// see enders prb 1995 (4c') and bahder prb 1990 eq. (30)
	// ----------------------------------------------------- 
	const double n = sqrt(3.0) * d;
	const double l = av + 2.0 * b;
	const double m = av - b;
		

	ostringstream sout; 
	sout.setf(ios::left);
 	sout << "kp 6x6 matrix for material " << this->get_material_name() << " uses following parameters:\n"
 	     << "lut1         = " << lut1 / hbar_2m0    << "\n"
 	     << "lut2         = " << lut2 / hbar_2m0    << "\n" 
 	     << "lut3         = " << lut3 / hbar_2m0    << "\n"
 	     << "vb_edge      = " << vb_edge  << "\n"
 	     << "delta_so     = " << delta_so << "\n"
 	     << "L            = " << setw(10) << L / hbar_2m0       << (L      < 0.0 ? " (good)" : "(unstable, should be < 0)") << "\n"
 	     << "M            = " << setw(10) << M / hbar_2m0       << (M      < 0.0 ? " (good)" : "(unstable, should be < 0)") << "\n"
		 << "N            = " << setw(10) << Nprime / hbar_2m0  << (Nprime < 0.0 ? " (good)" : "(unstable, should be < 0)") << "\n"
		 << "N+           = " << Np / hbar_2m0      << "\n"
		 << "N-           = " << Nm / hbar_2m0      << "\n" 	     
 	     << "energy shift = " << this->energy_shift; 
 	      	
	if(this->strain_dependence_available) {
		sout << "\n"
		     << "strain dependency is available:\n"
		     << "strain av    = " << av << "\n"
		     << "strain b     = " << b  << "\n"
		     << "strain d     = " << d;			
	}
	if(!surpress_output()) {
 		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	}  

	if(Configuration::get_instance()->get("kpmatrix_disable_foreman_enable_symmetrization") == 1) {
		Logger::get_instance()->emit(LOG_WARN, "you requested to disable the burt/foreman operator ordering and use symmetrization. therefore setting Np = Nm = 0.5 N'");
		Np = Nm = Nprime / 2.0;
	} 
	
	this->initialized = true;

	// --------------------------------------------------------------------------
	// diagonal kp terms d/dx d/dx
	// --------------------------------------------------------------------------									
/*(0,0)*/	this->second_order_strainless[D_DX][D_DX][get_idx(0,0)] = L;
/*(1,1)*/	this->second_order_strainless[D_DX][D_DX][get_idx(1,1)] = M;
/*(2,2)*/	this->second_order_strainless[D_DX][D_DX][get_idx(2,2)] = M;
/*(3,3)*/	this->second_order_strainless[D_DX][D_DX][get_idx(3,3)] = L;
/*(4,4)*/	this->second_order_strainless[D_DX][D_DX][get_idx(4,4)] = M;
/*(5,5)*/	this->second_order_strainless[D_DX][D_DX][get_idx(5,5)] = M;

	// --------------------------------------------------------------------------
	// diagonal kp terms d/dy d/dy
	// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DY][D_DY][get_idx(0,0)] = M;
/*(1,1)*/	this->second_order_strainless[D_DY][D_DY][get_idx(1,1)] = L;
/*(2,2)*/	this->second_order_strainless[D_DY][D_DY][get_idx(2,2)] = M;
/*(3,3)*/	this->second_order_strainless[D_DY][D_DY][get_idx(3,3)] = M;
/*(4,4)*/	this->second_order_strainless[D_DY][D_DY][get_idx(4,4)] = L;
/*(5,5)*/	this->second_order_strainless[D_DY][D_DY][get_idx(5,5)] = M;
	
	// --------------------------------------------------------------------------
	// diagonal kp terms d/dz d/dz
	// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(0,0)] = M;
/*(1,1)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(1,1)] = M;
/*(2,2)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(2,2)] = L;
/*(3,3)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(3,3)] = M;
/*(4,4)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(4,4)] = M;
/*(5,5)*/	this->second_order_strainless[D_DZ][D_DZ][get_idx(5,5)] = L;
	
	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dx d/dy
	// --------------------------------------------------------------------------
/*(0,1)*/	this->second_order_strainless[D_DX][D_DY][get_idx(0,1)]  = Np;
/*(1,0)*/	this->second_order_strainless[D_DX][D_DY][get_idx(1,0)]  = Nm;
/*(3,4)*/	this->second_order_strainless[D_DX][D_DY][get_idx(3,4)] = Np;
/*(4,3)*/	this->second_order_strainless[D_DX][D_DY][get_idx(4,3)] = Nm;
	
	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dy d/dx
	// --------------------------------------------------------------------------
/*(0,1)*/	this->second_order_strainless[D_DY][D_DX][get_idx(0,1)]  = Nm;
/*(1,0)*/	this->second_order_strainless[D_DY][D_DX][get_idx(1,0)]  = Np;
/*(3,4)*/	this->second_order_strainless[D_DY][D_DX][get_idx(3,4)] = Nm;
/*(4,3)*/	this->second_order_strainless[D_DY][D_DX][get_idx(4,3)] = Np;
	
	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dx d/dz
	// --------------------------------------------------------------------------	
/*(0,2)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(0,2)]  = Np;
/*(2,0)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(2,0)]  = Nm;
/*(3,5)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(3,5)] = Np;
/*(5,3)*/	this->second_order_strainless[D_DX][D_DZ][get_idx(5,3)] = Nm;

	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dz d/dx
	// --------------------------------------------------------------------------
/*(0,2)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(0,2)]  = Nm;
/*(2,0)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(2,0)]  = Np;
/*(3,5)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(3,5)] = Nm;
/*(5,3)*/	this->second_order_strainless[D_DZ][D_DX][get_idx(5,3)] = Np;

	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dy d/dz
	// --------------------------------------------------------------------------
/*(1,2)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(1,2)]  = Np;
/*(2,1)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(2,1)]  = Nm;
/*(4,5)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(4,5)] = Np;
/*(5,4)*/	this->second_order_strainless[D_DY][D_DZ][get_idx(5,4)] = Nm;

	// --------------------------------------------------------------------------
	// offdiagonal kp terms d/dz d/dy
	// --------------------------------------------------------------------------
/*(1,2)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(1,2)]  = Nm;
/*(2,1)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(2,1)]  = Np;
/*(4,5)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(4,5)] = Nm;
/*(5,4)*/	this->second_order_strainless[D_DZ][D_DY][get_idx(5,4)] = Np;

	// -------------------------------------------------------------------------
	// zero order terms (bandedge strain etc.)
	// -------------------------------------------------------------------------
/*(0,0)*/	this->zero_order_strainless[get_idx(0,0)]  = vb_edge - delta_so / 3.0;
/*(0,1)*/	this->zero_order_strainless[get_idx(0,1)]  =     - i * delta_so / 3.0;
/*(0,5)*/	this->zero_order_strainless[get_idx(0,5)]  =           delta_so / 3.0;
/*(1,0)*/	this->zero_order_strainless[get_idx(1,0)]  =       i * delta_so / 3.0;
/*(1,1)*/	this->zero_order_strainless[get_idx(1,1)]  = vb_edge - delta_so / 3.0;
/*(1,5)*/	this->zero_order_strainless[get_idx(1,5)]  =     - i * delta_so / 3.0;
/*(2,2)*/	this->zero_order_strainless[get_idx(2,2)] = vb_edge - delta_so / 3.0;
/*(2,3)*/	this->zero_order_strainless[get_idx(2,3)] =         - delta_so / 3.0;
/*(2,4)*/	this->zero_order_strainless[get_idx(2,4)] =       i * delta_so / 3.0;
/*(3,2)*/	this->zero_order_strainless[get_idx(3,2)] =         - delta_so / 3.0;
/*(3,3)*/	this->zero_order_strainless[get_idx(3,3)] = vb_edge - delta_so / 3.0;
/*(3,4)*/	this->zero_order_strainless[get_idx(3,4)] =       i * delta_so / 3.0;
/*(4,2)*/	this->zero_order_strainless[get_idx(4,2)] =     - i * delta_so / 3.0;
/*(4,3)*/	this->zero_order_strainless[get_idx(4,3)] =     - i * delta_so / 3.0;
/*(4,4)*/	this->zero_order_strainless[get_idx(4,4)] = vb_edge - delta_so / 3.0;
/*(5,0)*/	this->zero_order_strainless[get_idx(5,0)] =           delta_so / 3.0;
/*(5,1)*/	this->zero_order_strainless[get_idx(5,1)] =       i * delta_so / 3.0;
/*(5,5)*/	this->zero_order_strainless[get_idx(5,5)] = vb_edge - delta_so / 3.0;

	// -------------------------------------------------------
	// zero order strain dependent terms
	// -------------------------------------------------------
	if(this->strain_dependence_available) {
		const short EX = 0; 
		const short EY = 1;
		const short EZ = 2;
			
			// ---------------------------------------------
			// EXX terms
			// ---------------------------------------------
/*(0,0)*/	this->zero_order_strain_dependent[EX][EX][get_idx(0,0)]  = l;
/*(1,1)*/	this->zero_order_strain_dependent[EX][EX][get_idx(1,1)]  = m;
/*(2,2)*/	this->zero_order_strain_dependent[EX][EX][get_idx(2,2)] = m;

/*(3,3)*/	this->zero_order_strain_dependent[EX][EX][get_idx(3,3)] = l;
/*(4,4)*/	this->zero_order_strain_dependent[EX][EX][get_idx(4,4)] = m;
/*(5,5)*/	this->zero_order_strain_dependent[EX][EX][get_idx(5,5)] = m;

			// ---------------------------------------------
			// EYY terms
			// ---------------------------------------------
/*(0,0)*/	this->zero_order_strain_dependent[EY][EY][get_idx(0,0)]  = m;
/*(1,1)*/	this->zero_order_strain_dependent[EY][EY][get_idx(1,1)]  = l;
/*(2,2)*/	this->zero_order_strain_dependent[EY][EY][get_idx(2,2)] = m;

/*(3,3)*/	this->zero_order_strain_dependent[EY][EY][get_idx(3,3)] = m;
/*(4,4)*/	this->zero_order_strain_dependent[EY][EY][get_idx(4,4)] = l;
/*(5,5)*/	this->zero_order_strain_dependent[EY][EY][get_idx(5,5)] = m;
			// ---------------------------------------------
			// EZZ terms
			// ---------------------------------------------
/*(0,0)*/	this->zero_order_strain_dependent[EZ][EZ][get_idx(0,0)]  = m;
/*(1,1)*/	this->zero_order_strain_dependent[EZ][EZ][get_idx(1,1)]  = m;
/*(2,2)*/	this->zero_order_strain_dependent[EZ][EZ][get_idx(2,2)] = l;

/*(3,3)*/	this->zero_order_strain_dependent[EZ][EZ][get_idx(3,3)] = m;
/*(4,4)*/	this->zero_order_strain_dependent[EZ][EZ][get_idx(4,4)] = m;
/*(5,5)*/	this->zero_order_strain_dependent[EZ][EZ][get_idx(5,5)] = l;
			// ---------------------------------------------
			// EXY terms
			// ---------------------------------------------
/*(0,1)*/	this->zero_order_strain_dependent[EX][EY][get_idx(0,1)] = n;
/*(3,4)*/	this->zero_order_strain_dependent[EX][EY][get_idx(3,4)] = n;
			// ---------------------------------------------
			// EYX terms
			// ---------------------------------------------
/*(1,0)*/	this->zero_order_strain_dependent[EY][EX][get_idx(1,0)] = n;
/*(4,3)*/	this->zero_order_strain_dependent[EY][EX][get_idx(4,3)] = n;
			// ---------------------------------------------
			// EXZ terms
			// ---------------------------------------------
/*(0,2)*/	this->zero_order_strain_dependent[EX][EZ][get_idx(0,2)] = n;
/*(3,5)*/	this->zero_order_strain_dependent[EX][EZ][get_idx(3,5)] = n;

			// ---------------------------------------------
			// EZX terms
			// ---------------------------------------------
/*(2,0)*/	this->zero_order_strain_dependent[EZ][EX][get_idx(2,0)] = n;
/*(5,3)*/	this->zero_order_strain_dependent[EZ][EX][get_idx(5,3)] = n;
			// ---------------------------------------------
			// EYZ terms
			// ---------------------------------------------
/*(1,2)*/	this->zero_order_strain_dependent[EY][EZ][get_idx(1,2)] = n;
/*(4,5)*/	this->zero_order_strain_dependent[EY][EZ][get_idx(4,5)] = n;
			// ---------------------------------------------
			// EZY terms
			// ---------------------------------------------
/*(2,1)*/	this->zero_order_strain_dependent[EZ][EY][get_idx(2,1)] = n;
/*(5,4)*/	this->zero_order_strain_dependent[EZ][EY][get_idx(5,4)] = n;
	
		// -------------------------------------------------------
		// strain dependent spin orbit terms
		// see kp 8x8 model for further explanations!
		// -------------------------------------------------------
		if(Configuration::get_instance()->get("kpmatrix_include_spin_orbit_strain_dependence") == 1.0) {
			const double dso3 = delta_so / 3.0;
			const short EX = 0; 
			const short EY = 1;
			const short EZ = 2;			
			// --------------------------------------------------
			// EXX
			// --------------------------------------------------						
/*(0,1)*/ this->zero_order_strain_dependent[EX][EX][get_idx(0,1)] += dso3 * (i);
/*(0,5)*/ this->zero_order_strain_dependent[EX][EX][get_idx(0,5)] += dso3 * (-1);
/*(1,0)*/ this->zero_order_strain_dependent[EX][EX][get_idx(1,0)] += dso3 * (-i);
/*(2,3)*/ this->zero_order_strain_dependent[EX][EX][get_idx(2,3)] += dso3 * (1);
/*(3,2)*/ this->zero_order_strain_dependent[EX][EX][get_idx(3,2)] += dso3 * (1);
/*(3,4)*/ this->zero_order_strain_dependent[EX][EX][get_idx(3,4)] += dso3 * (-i);
/*(4,3)*/ this->zero_order_strain_dependent[EX][EX][get_idx(4,3)] += dso3 * (i);
/*(5,0)*/ this->zero_order_strain_dependent[EX][EX][get_idx(5,0)] += dso3 * (-1);

			// --------------------------------------------------
			// EYY
			// --------------------------------------------------						
/*(0,1)*/ this->zero_order_strain_dependent[EY][EY][get_idx(0,1)] += dso3 * (i);
/*(1,0)*/ this->zero_order_strain_dependent[EY][EY][get_idx(1,0)] += dso3 * (-i);
/*(1,5)*/ this->zero_order_strain_dependent[EY][EY][get_idx(1,5)] += dso3 * (i);
/*(2,4)*/ this->zero_order_strain_dependent[EY][EY][get_idx(2,4)] += dso3 * (-i);
/*(3,4)*/ this->zero_order_strain_dependent[EY][EY][get_idx(3,4)] += dso3 * (-i);
/*(4,2)*/ this->zero_order_strain_dependent[EY][EY][get_idx(4,2)] += dso3 * (i);
/*(4,3)*/ this->zero_order_strain_dependent[EY][EY][get_idx(4,3)] += dso3 * (i);
/*(5,1)*/ this->zero_order_strain_dependent[EY][EY][get_idx(5,1)] += dso3 * (-i);


			// --------------------------------------------------
			// EZZ
			// --------------------------------------------------						
/*(0,5)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(0,5)] += dso3 * (-1);
/*(1,5)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(1,5)] += dso3 * (i);
/*(2,3)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(2,3)] += dso3 * (1);
/*(2,4)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(2,4)] += dso3 * (-i);
/*(3,2)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(3,2)] += dso3 * (1);
/*(4,2)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(4,2)] += dso3 * (i);
/*(5,0)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(5,0)] += dso3 * (-1);
/*(5,1)*/ this->zero_order_strain_dependent[EZ][EZ][get_idx(5,1)] += dso3 * (-i);


			// --------------------------------------------------
			// EXY
			// --------------------------------------------------						
/*(0,5)*/ this->zero_order_strain_dependent[EX][EY][get_idx(0,5)] += dso3 * (-i);
/*(1,5)*/ this->zero_order_strain_dependent[EX][EY][get_idx(1,5)] += dso3 * (1);
/*(2,3)*/ this->zero_order_strain_dependent[EX][EY][get_idx(2,3)] += dso3 * (i);
/*(2,4)*/ this->zero_order_strain_dependent[EX][EY][get_idx(2,4)] += dso3 * (-1);
/*(3,2)*/ this->zero_order_strain_dependent[EX][EY][get_idx(3,2)] += dso3 * (-i);
/*(4,2)*/ this->zero_order_strain_dependent[EX][EY][get_idx(4,2)] += dso3 * (-1);
/*(5,0)*/ this->zero_order_strain_dependent[EX][EY][get_idx(5,0)] += dso3 * (i);
/*(5,1)*/ this->zero_order_strain_dependent[EX][EY][get_idx(5,1)] += dso3 * (1);


			// --------------------------------------------------
			// EXZ
			// --------------------------------------------------						
/*(0,4)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(0,4)] += dso3 * (-i);
/*(1,2)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(1,2)] += dso3 * (-i);
/*(1,3)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(1,3)] += dso3 * (i);
/*(2,1)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(2,1)] += dso3 * (i);
/*(3,1)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(3,1)] += dso3 * (-i);
/*(4,0)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(4,0)] += dso3 * (i);
/*(4,5)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(4,5)] += dso3 * (i);
/*(5,4)*/ this->zero_order_strain_dependent[EX][EZ][get_idx(5,4)] += dso3 * (-i);


			// --------------------------------------------------
			// EYZ
			// --------------------------------------------------						
/*(0,2)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(0,2)] += dso3 * (-i);
/*(0,4)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(0,4)] += dso3 * (1);
/*(1,3)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(1,3)] += dso3 * (-1);
/*(2,0)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(2,0)] += dso3 * (i);
/*(3,1)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(3,1)] += dso3 * (-1);
/*(3,5)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(3,5)] += dso3 * (i);
/*(4,0)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(4,0)] += dso3 * (1);
/*(5,3)*/ this->zero_order_strain_dependent[EY][EZ][get_idx(5,3)] += dso3 * (-i);

			 
		}		
	
	}			
}												



} // end of namespace
