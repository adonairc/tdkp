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

#include "tdkp/kpmatrices/KPMatrix4x4EndersForeman.h"
#include "tdkp/common/Configuration.h"

namespace tdkp
{


/** sparsity pattern of the kp matrix */
const int KPMatrix4x4EndersForeman::sparsity_pattern[28] = {
	0,0,	0,1,	0,2,
	1,0,	1,1,	1,2,	1,3,
	2,0,	2,1,	2,2,	2,3,
			3,1,	3,2,	3,3
};

/** the diagonal elements in the kp matrix are stored at the following locations */
const int KPMatrix4x4EndersForeman::diagonal_pattern[4] = {0,4,9,13};

const int* KPMatrix4x4EndersForeman::get_sparsity_pattern(int &num) const {
	num = 14;
	return &(this->sparsity_pattern[0]);
}

const int* KPMatrix4x4EndersForeman::get_diagonal_pattern(int &num) const {
	num = 4;
	return &(this->diagonal_pattern[0]);
}

KPMatrix4x4EndersForeman::KPMatrix4x4EndersForeman(const Material* mat) {
	this->set_material(mat);
	// let kp matrix assembler in base class know that we dont have first order terms
	this->have_first_order_terms = false;	
	// matrix space must be allocated from derived class
	// allocate space calls virtual function, which is defined in derived class
	// so it's not available in the default constructor of the base class
	this->allocate_matrix_space();
}

bool KPMatrix4x4EndersForeman::check_solution_type_available(KPSolutionType type) const {
	switch(type) {
		case electrons:
			return false;
		case holes:
			return true;
		default:
			TDKP_GENERAL_EXCEPTION("unknown kp solution type");	
	}	
}

KPMatrix4x4EndersForeman::KPMatrix4x4EndersForeman() {
	// let kp matrix assembler in base class know that we dont have first order terms
	this->have_first_order_terms = false;	
	this->allocate_matrix_space();
	this->material = 0;
}

KPMatrix4x4EndersForeman::~KPMatrix4x4EndersForeman() {
	material = 0;
}

void KPMatrix4x4EndersForeman::init_base_and_strain_matrix() {

	TDKP_ASSERT(this->material != 0,"this->material != 0");

	ostringstream sout;
	double lut1, lut2, lut3, vb_edge, hbar_2m0, An;	
	complex<double> i(0,1);
	hbar_2m0 = constants::hbar_square / (2.0 * constants::m0);

	lut1    = hbar_2m0 * this->material->get("luttinger_parameter_1");
	lut2    = hbar_2m0 * this->material->get("luttinger_parameter_2");
	lut3    = hbar_2m0 * this->material->get("luttinger_parameter_3");
	vb_edge = this->material->get("valence_band_edge");
 	An      = hbar_2m0 + lut1 - 2.0 * lut2 - 3.0 * lut3; // anisotropy factor resulting from proper derivation of kp hamiltonian

 	sout << "kp4x4 matrix for material " << this->get_material_name() << " uses following parameters:\n"
 	     << "lut1         = " << lut1 / hbar_2m0 << "\n"
 	     << "lut2         = " << lut2 / hbar_2m0 << "\n" 
 	     << "lut3         = " << lut3 / hbar_2m0 << "\n"
 	     << "vb_edge      = " << vb_edge << "\n"
 	     << "anisotropy   = " << An / hbar_2m0   << "\n"
 	     << "energy shift = " << this->energy_shift; 
 	   	   	 
 	if(!surpress_output()) {    	
 		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
 	}

	if(Configuration::get_instance()->get("kpmatrix_disable_foreman_enable_symmetrization") == 1) {
		Logger::get_instance()->emit(LOG_WARN, "you requested to disable the burt/foreman operator ordering and use symmetrization. therefore setting the anisotropy to 0");
		An = 0;	
	} 	   	
 	
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
			      << "for usual deformation potentials, av > 0, b < 0.0 and d < 0.0. "
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
	 		
	this->initialized = true;

			// --------------------------------------------------------------------------
			// diagonal kp terms d/dx d/dx
			// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DX][D_DX][0]  = - (lut1 + lut2);
/*(0,2)*/	this->second_order_strainless[D_DX][D_DX][2]  = constants::sqrt3 * lut2;
/*(1,1)*/	this->second_order_strainless[D_DX][D_DX][4]  = - (lut1 - lut2);
/*(1,3)*/	this->second_order_strainless[D_DX][D_DX][6]  = constants::sqrt3 * lut2;
/*(2,0)*/	this->second_order_strainless[D_DX][D_DX][7]  = constants::sqrt3 * lut2;
/*(2,2)*/	this->second_order_strainless[D_DX][D_DX][9]  = - (lut1 - lut2);
/*(3,1)*/	this->second_order_strainless[D_DX][D_DX][11] = constants::sqrt3 * lut2;
/*(3,3)*/	this->second_order_strainless[D_DX][D_DX][13] = - (lut1 + lut2);

			// --------------------------------------------------------------------------
			// diagonal kp terms d/dy d/dy
			// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DY][D_DY][0]  = - (lut1 + lut2);
/*(0,2)*/	this->second_order_strainless[D_DY][D_DY][2]  = - constants::sqrt3 * lut2;
/*(1,1)*/	this->second_order_strainless[D_DY][D_DY][4]  = - (lut1 - lut2);
/*(1,3)*/	this->second_order_strainless[D_DY][D_DY][6]  = - constants::sqrt3 * lut2;
/*(2,0)*/	this->second_order_strainless[D_DY][D_DY][7]  = - constants::sqrt3 * lut2;
/*(2,2)*/	this->second_order_strainless[D_DY][D_DY][9]  = - (lut1 - lut2);
/*(3,1)*/	this->second_order_strainless[D_DY][D_DY][11] = - constants::sqrt3 * lut2;
/*(3,3)*/	this->second_order_strainless[D_DY][D_DY][13] = - (lut1 + lut2);

			// --------------------------------------------------------------------------
			// diagonal kp terms d/dz d/dz
			// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DZ][D_DZ][0]  = - (lut1 - 2.0 * lut2);
/*(1,1)*/	this->second_order_strainless[D_DZ][D_DZ][4]  = - (lut1 + 2.0 * lut2);
/*(2,2)*/	this->second_order_strainless[D_DZ][D_DZ][9]  = - (lut1 + 2.0 * lut2);
/*(3,3)*/	this->second_order_strainless[D_DZ][D_DZ][13] = - (lut1 - 2.0 * lut2);

			// --------------------------------------------------------------------------
			// offdiagonal kp terms kx ky (here comes foreman into the game
			// if the An factor is ommitted, we probably loose coercivity and 
			// get therefore numerical noise instead of converged results
			// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DX][D_DY][0]  = - i * An;
/*(0,2)*/	this->second_order_strainless[D_DX][D_DY][2]  =   i * constants::sqrt3 * lut3;
/*(1,1)*/	this->second_order_strainless[D_DX][D_DY][4]  = - i * An / 3.0;
/*(1,3)*/	this->second_order_strainless[D_DX][D_DY][6]  =   i * constants::sqrt3 * lut3;
/*(2,0)*/	this->second_order_strainless[D_DX][D_DY][7]  = - i * constants::sqrt3 * lut3;
/*(2,2)*/	this->second_order_strainless[D_DX][D_DY][9]  =   i * An / 3.0;
/*(3,1)*/	this->second_order_strainless[D_DX][D_DY][11] = - i * constants::sqrt3 * lut3;
/*(3,3)*/	this->second_order_strainless[D_DX][D_DY][13] =   i * An;

			// --------------------------------------------------------------------------
			// offdiagonal kp terms ky kx 
			// --------------------------------------------------------------------------
/*(0,0)*/	this->second_order_strainless[D_DY][D_DX][0]  =   i * An;
/*(0,2)*/	this->second_order_strainless[D_DY][D_DX][2]  =   i * constants::sqrt3 * lut3;
/*(1,1)*/	this->second_order_strainless[D_DY][D_DX][4]  =   i * An / 3.0;
/*(1,3)*/	this->second_order_strainless[D_DY][D_DX][6]  =   i * constants::sqrt3 * lut3;
/*(2,0)*/	this->second_order_strainless[D_DY][D_DX][7]  = - i * constants::sqrt3 * lut3;
/*(2,2)*/	this->second_order_strainless[D_DY][D_DX][9]  = - i * An / 3.0;
/*(3,1)*/	this->second_order_strainless[D_DY][D_DX][11] = - i * constants::sqrt3 * lut3;
/*(3,3)*/	this->second_order_strainless[D_DY][D_DX][13] = - i * An;

			// --------------------------------------------------------------------------
			// offdiagonal kp terms d/dx d/dz
			// --------------------------------------------------------------------------
/*(0,1)*/	this->second_order_strainless[D_DX][D_DZ][1]  = - (An - 3.0 * lut3) / constants::sqrt3;
/*(1,0)*/	this->second_order_strainless[D_DX][D_DZ][3]  =   (An + 3.0 * lut3) / constants::sqrt3;
/*(1,2)*/	this->second_order_strainless[D_DX][D_DZ][5]  = - 2.0 / 3.0 * An;
/*(2,1)*/	this->second_order_strainless[D_DX][D_DZ][8]  =   2.0 / 3.0 * An;
/*(2,3)*/	this->second_order_strainless[D_DX][D_DZ][10] = - (An + 3.0 * lut3) / constants::sqrt3;
/*(3,2)*/	this->second_order_strainless[D_DX][D_DZ][12] =   (An - 3.0 * lut3) / constants::sqrt3;

			// --------------------------------------------------------------------------
			// offdiagonal kp terms d/dz d/dx
			// --------------------------------------------------------------------------
/*(0,1)*/	this->second_order_strainless[D_DZ][D_DX][1]  =   (An + 3.0 * lut3) / constants::sqrt3;
/*(1,0)*/	this->second_order_strainless[D_DZ][D_DX][3]  = - (An - 3.0 * lut3) / constants::sqrt3;
/*(1,2)*/	this->second_order_strainless[D_DZ][D_DX][5]  =   2.0 / 3.0 * An;
/*(2,1)*/	this->second_order_strainless[D_DZ][D_DX][8]  = - 2.0 / 3.0 * An;
/*(2,3)*/	this->second_order_strainless[D_DZ][D_DX][10] =   (An - 3.0 * lut3) / constants::sqrt3;
/*(3,2)*/	this->second_order_strainless[D_DZ][D_DX][12] = - (An + 3.0 * lut3) / constants::sqrt3;

			// --------------------------------------------------------------------------
			// offdiagonal kp terms d/dy d/dz
			// --------------------------------------------------------------------------
/*(0,1)*/	this->second_order_strainless[D_DY][D_DZ][1]  = - i * (An - 3.0 * lut3) / constants::sqrt3;
/*(1,0)*/	this->second_order_strainless[D_DY][D_DZ][3]  = - i * (An + 3.0 * lut3) / constants::sqrt3;
/*(1,2)*/	this->second_order_strainless[D_DY][D_DZ][5]  = - i * 2.0 / 3.0 * An;
/*(2,1)*/	this->second_order_strainless[D_DY][D_DZ][8]  = - i * 2.0 / 3.0 * An;
/*(2,3)*/	this->second_order_strainless[D_DY][D_DZ][10] = - i * (An + 3.0 * lut3) / constants::sqrt3;
/*(3,2)*/	this->second_order_strainless[D_DY][D_DZ][12] = - i * (An - 3.0 * lut3) / constants::sqrt3;

			// --------------------------------------------------------------------------
			// offdiagonal kp terms d/dz d/dy
			// --------------------------------------------------------------------------
/*(0,1)*/	this->second_order_strainless[D_DZ][D_DY][1]  =   i * (An + 3.0 * lut3) / constants::sqrt3;
/*(1,0)*/	this->second_order_strainless[D_DZ][D_DY][3]  =   i * (An - 3.0 * lut3) / constants::sqrt3;
/*(1,2)*/	this->second_order_strainless[D_DZ][D_DY][5]  =   i * 2.0 / 3.0 * An;
/*(2,1)*/	this->second_order_strainless[D_DZ][D_DY][8]  =   i * 2.0 / 3.0 * An;
/*(2,3)*/	this->second_order_strainless[D_DZ][D_DY][10] =   i * (An - 3.0 * lut3) / constants::sqrt3;
/*(3,2)*/	this->second_order_strainless[D_DZ][D_DY][12] =   i * (An + 3.0 * lut3) / constants::sqrt3;

			// -------------------------------------------------------------------------
			// zero order terms (bandedge strain etc.)
			// -------------------------------------------------------------------------
/*(0,0)*/	this->zero_order_strainless[0]  = vb_edge;
/*(1,1)*/	this->zero_order_strainless[4]  = vb_edge;
/*(2,2)*/	this->zero_order_strainless[9]  = vb_edge;
/*(3,3)*/	this->zero_order_strainless[13] = vb_edge;

	// -------------------------------------------------------
	// create strain dependent terms (only zero order stuff for chuangs matrices)
	// -------------------------------------------------------
	if(this->strain_dependence_available) {

			// -------------------------------------------------------
			// o.k., here is the story:
			// this is the deformation potential matrix transferred
			// from the X,Y,Z basis to the drehimpuls basis using the
			// sebi transformation. basically this means, use
			// the D3 hamiltonian in enders (4c') and apply the unitary 
			// transformation that diagonalized the spin orbit coupling
			// (this is the basis bahder and chuang are usually working)
			// 
			// the same can be read in pikus and bir, p. 310.
			//
			// now if you check the code repository you will notice that 
			// something important changed:
			// the sign of av. 
			// the reason is the following: initially i followed pryor
			// but as it turned out with the whole spurious solution story
			// its better to use the raw enders hamiltonian. 
			// pryor relates his work on bahders 8x8 hamiltonian.
			// now the problem is as follows: bahder published in 1992 an
			// errata to his 1990 hamiltonian, where he corrected 
			// his a' parameter to be -a' (see bahder 1990, eq 30, compare
			// with pikus and bir p 310, and then check bahders errata in 1992)
			// but pryor just to the a' parameter with the wrong sign and
			// declared it to be av.
			//
			// so what happens in his model? if you have compressive strain,
			// then the bandedge change is dEv = av * tr(e) and dEc = ac * tr(e)
			// so the bandedge deformation potential would be ag = ac - av
			//
			// but pryors bandedges change as dEv = -av * tr(e) 
			// so, in order to still have the same ag, he simply says
			// ag = ac + av and therefore takes an ac that is bigger
			// than the usual ag to compensate the wrong shift of his
			// dEv term.
			// so he finally ends up with a more or less correct nett
			// separation. but its not correct.
			// -------------------------------------------------------


			// -------------------------------------------------------
			// exx terms
			// -------------------------------------------------------
/*(0,0)*/	this->zero_order_strain_dependent[0][0][0]  =   av + b / 2.0;
/*(0,2)*/	this->zero_order_strain_dependent[0][0][2]  = - constants::sqrt3 * b / 2.0;
/*(1,1)*/	this->zero_order_strain_dependent[0][0][4]  =   av - b / 2.0;
/*(1,3)*/	this->zero_order_strain_dependent[0][0][6]  = - constants::sqrt3 * b / 2.0;
/*(2,0)*/	this->zero_order_strain_dependent[0][0][7]  = - constants::sqrt3 * b / 2.0;
/*(2,2)*/	this->zero_order_strain_dependent[0][0][9]  =   av - b / 2.0;
/*(3,1)*/	this->zero_order_strain_dependent[0][0][11] = - constants::sqrt3 * b / 2.0;
/*(3,3)*/	this->zero_order_strain_dependent[0][0][13] =   av + b / 2.0;

			// -------------------------------------------------------
			// eyy terms
			// -------------------------------------------------------
/*(0,0)*/	this->zero_order_strain_dependent[1][1][0]  =   av + b / 2.0;
/*(0,2)*/	this->zero_order_strain_dependent[1][1][2]  =   constants::sqrt3 * b / 2.0;
/*(1,1)*/	this->zero_order_strain_dependent[1][1][4]  =   av - b / 2.0;
/*(1,3)*/	this->zero_order_strain_dependent[1][1][6]  =   constants::sqrt3 * b / 2.0;
/*(2,0)*/	this->zero_order_strain_dependent[1][1][7]  =   constants::sqrt3 * b / 2.0;
/*(2,2)*/	this->zero_order_strain_dependent[1][1][9]  =   av - b / 2.0;
/*(3,1)*/	this->zero_order_strain_dependent[1][1][11] =   constants::sqrt3 * b / 2.0;
/*(3,3)*/	this->zero_order_strain_dependent[1][1][13] =   av + b / 2.0;

			// -------------------------------------------------------
			// ezz terms
			// -------------------------------------------------------
/*(0,0)*/	this->zero_order_strain_dependent[2][2][0]  =  av - b;
/*(1,1)*/	this->zero_order_strain_dependent[2][2][4]  =  av + b;
/*(2,2)*/	this->zero_order_strain_dependent[2][2][9]  =  av + b;
/*(3,3)*/	this->zero_order_strain_dependent[2][2][13] =  av - b;

			// -------------------------------------------------------
			// exy terms (ignoring diagonal terms that cancel)
			// -------------------------------------------------------
/*(0,2)*/	this->zero_order_strain_dependent[0][1][2]  = - i * d / 2.0;
/*(1,3)*/	this->zero_order_strain_dependent[0][1][6]  = - i * d / 2.0;
/*(2,0)*/	this->zero_order_strain_dependent[0][1][7]  =   i * d / 2.0;
/*(3,1)*/	this->zero_order_strain_dependent[0][1][11] =   i * d / 2.0;

			// -------------------------------------------------------
			// eyx terms 
			// -------------------------------------------------------
/*(0,2)*/	this->zero_order_strain_dependent[1][0][2]  = - i * d / 2.0;
/*(1,3)*/	this->zero_order_strain_dependent[1][0][6]  = - i * d / 2.0;
/*(2,0)*/	this->zero_order_strain_dependent[1][0][7]  =   i * d / 2.0;
/*(3,1)*/	this->zero_order_strain_dependent[1][0][11] =   i * d / 2.0;

			// -------------------------------------------------------
			// exz terms
			// -------------------------------------------------------
/*(0,1)*/	this->zero_order_strain_dependent[0][2][1]  = - d;
/*(3,2)*/	this->zero_order_strain_dependent[0][2][12] =   d;

			// -------------------------------------------------------
			// ezx terms
			// -------------------------------------------------------
/*(1,0)*/	this->zero_order_strain_dependent[2][0][3]  = - d;
/*(2,3)*/	this->zero_order_strain_dependent[2][0][10] =   d;

			// -------------------------------------------------------
			// eyz
			// -------------------------------------------------------
/*(0,1)*/	this->zero_order_strain_dependent[1][2][1]  = - i * d;
/*(3,2)*/	this->zero_order_strain_dependent[1][2][12] = - i * d;

			// -------------------------------------------------------
			// ezy
			// -------------------------------------------------------
/*(1,0)*/	this->zero_order_strain_dependent[2][1][3]  =   i * d;
/*(2,3)*/	this->zero_order_strain_dependent[2][1][10] =   i * d;
	
	} // end if strain_dependence_available
}





} // end of namespace
