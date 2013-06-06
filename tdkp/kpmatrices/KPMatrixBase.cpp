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

#include "tdkp/kpmatrices/KPMatrixBase.h"
#include "tdkp/common/Configuration.h"

namespace tdkp {

/** base constructor ... */		
KPMatrixBase::KPMatrixBase()
: have_first_order_terms(true), // safe default value
  potential_energy(0.0),
  energy_shift(0.0),
  rotation_matrix(3,3),
  surpress_output_to_user(false)
 {
	this->init();	
}

KPMatrixBase::~KPMatrixBase() {
	this->material = 0;	
}

/** standard constructor
 * 
 * @param material material object to read material properties from
 */
KPMatrixBase::KPMatrixBase(const Material* material_)
: have_first_order_terms(true),
  potential_energy(0.0),
  energy_shift(0.0),
  rotation_matrix(3,3),
  surpress_output_to_user(false)
{
	this->init();
	this->material    	   = material_;
	this->material_changed = true;
	this->material_name    = "unknown";
}

/** init from other kp matrix 
 * 
 * copies rotation / strain tensor / potential etc 
 */ 
void KPMatrixBase::init_from_kpmatrix(const KPMatrixBase& base) {

	// ------------------------------------------------
	// copy relevant params
	// ------------------------------------------------
	this->material = base.material;
	this->material_changed = true;
	this->surpress_output_to_user = base.surpress_output_to_user;
	this->rotated  = base.rotated;
	this->rotation_matrix = base.rotation_matrix;		
	this->energy_shift = base.energy_shift;
	this->calculate();	
	if(base.strain_tensor_set) {
		this->set_strains(base.strain);	
	} else {
		this->strain_tensor_set = false; 	
	}
	this->potential_energy = base.potential_energy;
	this->potential_energy_changed = true;

	 		
}

/** standard class initialization routine */
void KPMatrixBase::init() {
	this->initialized                 = false;
	this->material                    = 0;	
	this->rotated                     = false;
	this->strain_tensor_set           = false;

	this->energy_shift                = 0.0;
	this->potential_energy            = 0.0;
	this->potential_energy_changed    = false;	
	this->strain_tensor_changed       = false;
	this->material_changed            = false;
	this->rotation_changed         	  = false;
	this->strain_dependence_available = false;	
	this->surpress_output_to_user     = false;

	// -------------------------------------------------
	// display ONE warning message when strains have
	// manually been disabled
	// -------------------------------------------------
	static bool warned_on_disabled_strain = false;	
	if(!warned_on_disabled_strain && (
	     Configuration::get_instance()->get("kpmatrix_ignore_second_order_strain_dependence") == 1.0 ||
	     Configuration::get_instance()->get("kpmatrix_ignore_first_order_strain_dependence") == 1.0 ||
	     Configuration::get_instance()->get("kpmatrix_ignore_zero_order_strain_dependence") == 1.0
	   )){
		ostringstream sout;
		sout << "you have manually decided by a configuration directive to "
		     << "disable strain effects in the kp matrix. therefore i will "
		     << "ignore the following strain dependence of the kp matrix: all ";
		int num = 0;
		if(Configuration::get_instance()->get("kpmatrix_ignore_second_order_strain_dependence") == 1.0) {
			sout << "2nd";
			num++;			
		}
		if(Configuration::get_instance()->get("kpmatrix_ignore_first_order_strain_dependence") == 1.0) {
			if(num > 0) {
				sout << ", ";	
			}
			sout << "1st";
		}
		if(Configuration::get_instance()->get("kpmatrix_ignore_zero_order_strain_dependence") == 1.0) {
			if(num > 0) {
				sout << ", ";	
			}
			sout << "0th";
		}   
		sout << " order terms kp terms";
		Logger::get_instance()->emit(LOG_WARN, sout.str());
		warned_on_disabled_strain = true; // only once         	
	}	   	
}

bool KPMatrixBase::ready() const {	
	if(!this->initialized || this->potential_energy_changed || this->material_changed || this->rotation_changed || this->strain_tensor_changed || this->material == 0) {
		return false;	
	} else {
		return true;	
	}
}

/** set rotation of the crystal system
 * 
 * assigns the rotation matrix (given in the usual matrix[rowidx * 3 + colidx] form)
 * to the kp matrix.
 * 
 * you need to call calculate to rotate the matrix to the desired position
 * remember! rotations are NOT cumulative! the whole matrix is recalculated
 * in the unrotated system and then rotated fully
 * 
 */
void KPMatrixBase::set_rotation(const RMatrix<double>& rotation_matrix_) {
	TDKP_ASSERT(rotation_matrix_.rows() == 3 && rotation_matrix_.cols() == 3, "rotation matrix must be a 3x3 matrix");
					
	this->rotation_matrix  = rotation_matrix_;
	this->rotated          = true;
	this->rotation_changed = true;
		
}



/** write kp matrix to file */
void KPMatrixBase::dump(const char* filename) const {
	ofstream fout(filename);
	ostringstream sout;
	if(fout) {
		fout << *this;
		fout.close(); 
		sout << "wrote kp matrix for " << this->get_material_name() << " to file " << filename;
		Logger::get_instance()->emit(LOG_INFO, sout.str()); 		
	} else {		
		sout << "can not open file " << filename << " for writing";
		TDKP_GENERAL_EXCEPTION(sout.str());	
	}	
}


ostream& operator<<(ostream& out, const KPMatrixBase& mat) {
	if(mat.ready()) {
		const char* diffop[] = {"d/dx", "d/dy", "d/dz"};
		//int pairs[] = {0,0, 1,1, 2,2, 1,2, 2,1, 1,3, 3,1, 2,3, 3,2};
		//int aa,bb;
		RMatrix<cplx>*  tmp;
		int width = 12;
		bool have_first_order = false;
		out.precision(3);
		out << " ---------------------------------------------------------\n"
		    << "  kp " 
		    << mat.get_number_of_bands() << "x" << mat.get_number_of_bands() 
		    << " matrix for material " << mat.get_material_name() << "\n\n"
		    << " strainless part second order:\n";
		// for all differential operators
		for(int aa = 0; aa < 3; aa++) {
			for(int bb = 0; bb < 3; bb++) {
				out << " = = = = = = =  " << diffop[aa] << " H(2) " << diffop[bb] << " = = = = = = = \n";
				tmp = mat.get_second_order_matrix(aa,bb);
				// for all matrix entries				
				for(unsigned int ii = 0; ii < tmp->rows(); ii++) {
					for(unsigned int jj = 0; jj < tmp->cols(); jj++) {
						out << setw(width) << (*tmp)(ii,jj) << " ";
					}	
					out << "\n";
				} 	
				out << "\n";
				delete tmp;
			}
			// check if first order terms are nonzero
			if(!have_first_order) {
				// check first order left and right
				for(short zz = 0; zz < 2; zz++) {
					if(zz == 0) {										
						tmp = mat.get_first_order_matrix(KPMatrixBase::op_left, aa);
					} else {
						tmp = mat.get_first_order_matrix(KPMatrixBase::op_right, aa);
					}
					for(unsigned int ii = 0; ii < tmp->rows(); ii++) {
						for(unsigned int jj = 0; jj < tmp->cols(); jj++) {
							if(abs((*tmp)(ii,jj)) != 0.0) {
								have_first_order = true;	
							}
						} 	
					}	
					delete tmp;
				}
			}
		} 		    
		if(have_first_order) {
			out << " strainless part first order:\n";
			for(int aa = 0; aa < 3; aa++) {
				// check first order left and right
				for(short zz = 0; zz < 2; zz++) {
					if(zz == 0) {										
						tmp = mat.get_first_order_matrix(KPMatrixBase::op_left, aa);
						out << " = = = = = = =   H(1) (left) " << diffop[aa] << " = = = = = = = \n";						
					} else {						
						tmp = mat.get_first_order_matrix(KPMatrixBase::op_right, aa);
						out << " = = = = = = =   H(1) (right) " << diffop[aa] << " = = = = = = = \n";						
					}				

					for(unsigned int ii = 0; ii < tmp->rows(); ii++) {
						for(unsigned int jj = 0; jj < tmp->cols(); jj++) {
							out << setw(width) << (*tmp)(ii,jj) << " ";
						}		
						out << "\n";
					}
					out << "\n";	
					delete tmp;
				}											
			}	
		} else {
			out << " strainless part first order is zero ... skipping \n";	
		}
		out << "\n strainless part zero order H(0):\n";
		tmp = mat.get_zero_order_matrix();
		for(unsigned int ii = 0; ii < tmp->rows(); ii++) {
			for(unsigned int jj = 0; jj < tmp->cols(); jj++) {
				out << setw(width) << (*tmp)(ii,jj) << " ";
			}	
			out << "\n";
		}
		out << "\n";			
		delete tmp;    
		
		// output strain dependent part	
		out << "\n straindependent part zero order D(e):\n";
		for(unsigned int ee = 0; ee < 3; ee++) {
			for(unsigned int ff = 0; ff < 3; ff++) {
				out << " = = = = = = =  D(e"<< ee << ff << ") = = = = = = = \n";				
				tmp = mat.get_zero_order_strain_dependent_matrix(ee,ff);
				for(unsigned int ii = 0; ii < tmp->rows(); ii++) {
					for(unsigned int jj = 0; jj < tmp->cols(); jj++) {
						out << setw(width) << (*tmp)(ii,jj) << " ";
					}	
					out << "\n";
				}
				out << "\n";			
				delete tmp;
			}
		}    			
					
	} else {
		out << "kp matrix has not yet been initialized";	
	}	
	return out;
}

/** calculates the rotated second, first and zero order terms
 * 
 * throws exception if material pointer is zero
 */
void KPMatrixBase::calculate() throw(Exception*) {	
		
	TDKP_ASSERT(this->material != 0, "Material must be set before calling calculate");	
	// material or rotation update or not yet initialized -> need full recalculation
	if(!this->initialized || this->material_changed || this->rotation_changed) {		
		this->allocate_matrix_space();
		this->init_base_and_strain_matrix();		
		if(this->rotated) {
			ostringstream sout;
			sout << "rotating kp matrix using the given rotation matrix:\n"
			     << this->rotation_matrix
			     << "\n";		
			if(!surpress_output()) {			     	  
				Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
			}
			// rotate strain matrices to use 'new' strain basis
			this->rotate_strain_matrix(); 
			// then rotate all matrices (also strain dependent) to use the new 'k' basis
			this->rotate_base_matrix();	
		}
		this->material_changed = this->rotation_changed = this->strain_tensor_changed = this->potential_energy_changed = false;	
		this->calculate_zero_order();	
		this->calculate_first_order();
		this->calculate_second_order();		
	// only strain tensor / potential energy changed -> update necessary stuff		
	} else if(this->strain_tensor_changed) {		
		this->calculate_zero_order();
		this->potential_energy_changed = false;			
		this->calculate_second_order();
		if(this->have_first_order_terms) {
			this->calculate_first_order();
		}
		this->strain_tensor_changed = false;
	} 
	if(this->potential_energy_changed) {		
		this->calculate_zero_order();
		this->potential_energy_changed = false;
	}			
}

/** enforce recalculation
 * 
 * calculate does only recalculate the necessary updates. but 
 * if any values e.g. in the material database are changed, the 
 * change won't be detected. therefore one needs then to 
 * enforce a full recalculation of the kp matrix using this
 * function
 */ 
void KPMatrixBase::enforce_recalculation() {
	this->initialized = false;
	this->calculate(); 
}


/** sets the strain tensor but does not recalculate 
 * 
 * call calculate after. recalculation will then only affect strain dependent parts
 * Throws exception if strain dependence is not available (either because not 
 * implemented or if parameters are not set)
 * 
 * @param strain_ StrainTensor with current strain
 */
void KPMatrixBase::set_strains(const StrainTensor& strain_) throw(Exception*) { 
	TDKP_ASSERT(this->initialized, "matrix must be initialized (call calculate()) before strains may be set");
	if(this->strain_dependence_available == false) {
		TDKP_GENERAL_EXCEPTION("strain dependence not available.");
	}
	this->strain = strain_;		
	this->strain_tensor_changed = true;
	this->strain_tensor_set = true;
}

/** set potential energy term that should be added to the bandedges
 * 
 * call calculate after. recalculation will only affect potential energy dependent parts
 * 
 * @param potential_energy_ energy term that should be added to the bandedges
 */
void KPMatrixBase::set_potential(double potential_energy_) {
	this->potential_energy = potential_energy_;
	this->potential_energy_changed = true;
}

/** set energy shift to shift the spectrum up or down to the desired value
 * 
 * purpose: e.g. in solving the equation, one could shift the values to have the desired 
 * ones around 0 which would lead to the biggest largest magnitudes in shift & innode 
 * eigensolver iterations
 */
void KPMatrixBase::set_energy_shift(double eshift) {
	this->energy_shift = eshift;	
	this->potential_energy_changed = true;
}

/** set material */
void KPMatrixBase::set_material(const Material* material_) {
	this->material = material_;
	this->material_changed = true;
}

/** set material name */
void KPMatrixBase::set_material_name(const string& name) {
	this->material_name = name;
}
/** get material name */ 
const string& KPMatrixBase::get_material_name() const {
	return this->material_name;
}

/** allocate space for the kp matrix 
 * 
 * calls get_sparsity_pattern to determine the number of nonzeros in the kp matrix.
 * then, loops trough the kp matrix vectors and calls resize 
 */
void KPMatrixBase::allocate_matrix_space() {
	int num;
	this->get_sparsity_pattern(num);
	for(short ii = 0; ii < 3; ii++) {
		for(short ll = 0; ll < 2; ll++) {
			this->first_order[ll][ii].assign(num, 0.0);				
			this->first_order_strainless[ll][ii].assign(num, 0.0);
		}
		for(short jj = 0; jj < 3; jj++) {
			this->second_order[ii][jj].assign(num, 0.0);
			this->second_order_strainless[ii][jj].assign(num, 0.0);
		}	
	}
	this->zero_order.assign(num, 0.0);			
	this->zero_order_strainless.assign(num, 0.0);			
	for(short ee = 0; ee < 3; ee++) {
		for(short ff = 0; ff < 3; ff++) {			
			this->zero_order_strain_dependent[ee][ff].assign(num, 0.0);			
		}
	}	
}

/** rotates base matrices 
 * 
 * if matrix \f$\mathbf{R}: \mathbf{k} \rightarrow \mathbf{k'}\$  is the rotation matrix, then 
 * the new 'rotated' kp matrices are given by - 
 * for second order: \f$\mathbf{M'}^{(2)}_{ij} = \sum_{kl} r_{ik}r_{jl} \mathbf{M}^{(2)}_{kl} \f$
 * for first order: \f$\mathbf{M'}^{(1)}_{i} = \sum_{j} r_{ij} \mathbf{M}^{(1)}_{j} \f$
 * and zero order of course does not change ...
 * 
 */
void KPMatrixBase::rotate_base_matrix() {
	
	TDKP_ASSERT(this->initialized, "matrix is not initialized");
	
	// allocate space to temporarily store matrices
	int num;
	this->get_sparsity_pattern(num);
	vector<cplx> tmp_2nd[3][3];
	vector<cplx> tmp_1st[2][3];
	
	for(int ii = 0; ii < 3; ii++) {
		for(int jj = 0; jj < 3; jj++) {
			tmp_2nd[ii][jj].assign(num, 0.0);	
		}	
		for(int aa = 0; aa < 2; aa++) {
			tmp_1st[aa][ii].assign(num, 0.0);
		}
	}
	
	// rotate 2nd order stuff
	double r_ik_r_jl;
	for(int ii = 0; ii < 3; ii++) {
		for(int jj = 0; jj < 3; jj++) {
			for(int kk = 0; kk < 3; kk++) {
				for(int ll = 0; ll < 3; ll++) {
					r_ik_r_jl = this->rotation_matrix(ii,kk) * this->rotation_matrix(jj,ll);
					for(int nn = 0; nn < num; nn++) {
						tmp_2nd[ii][jj][nn] += r_ik_r_jl * this->second_order_strainless[kk][ll][nn]; 	
					}					
				}
			}
		}	
	}
	// rotate first order stuff
	if(this->have_first_order_terms) {
		double r_ij;
		for(int ii = 0; ii < 3; ii++) {
			for(int jj = 0; jj < 3; jj++) {
				r_ij = this->rotation_matrix(ii,jj);
				for(int aa = 0; aa < 2; aa++) {
					for(int nn = 0; nn < num; nn++) {
						tmp_1st[aa][ii][nn] += r_ij * this->first_order_strainless[aa][jj][nn];						
					}
				}
			}			
		}
	}
	// copy back
	for(int ii = 0; ii < 3; ii++) {
		// second order
		for(int jj = 0; jj < 3; jj++) {
			for(int nn = 0; nn < num; nn++) {
				this->second_order_strainless[ii][jj][nn] = tmp_2nd[ii][jj][nn];									
			}	
		}			
		// first order
		if(this->have_first_order_terms) {
			for(int aa = 0; aa < 2; aa++) {
				for(int nn = 0; nn < num; nn++) {			
					this->first_order_strainless[aa][ii][nn] = tmp_1st[aa][ii][nn];
				}
			}
		}
	}
	
}

/** rotate strain matrix in order to be written in terms of rotated strains
 * 
 * the strain-dependent part in the kp theory is expressed in terms of the strain tensor
 * \f$\mathbf{e}\f$. but if our system here is rotated, it must be expressed in terms of the 
 * rotated strain tensor \$\mathbf{e'}\$ coming from a preceeding strain calculation. 
 * it will serve then as input to calculate elementwise the strain in the kp matrix
 * 
 * so, let \f$\mathbf{R}\f$ be the rotation matrix with \f$\mathbf{Ru} = \mathbf{u'}\f$
 * and let \f$\mathbf{u}\f$ be the reference configuration and \f$\mathbf{u'}\f$ the rotated.
 * then, the strains transforms like \f$ \mathbf{e'} = \mathbf{ReR}^{T} \f$. 
 * expressing the strain in terms of the rotated system we have
 * \f$\mathbf{e} = \mathbf{R}^{T}\mathbf{e'R} \f$, insert it as the strain, swapping the sums,
 * we get an expression of the strain dependent stuff in terms of the rotated strain.
 *
 * this is important for the zero order terms. for the first and second order
 * terms we forget that. why? 
 * 
 * well: according to enders, prb 1995, the strain tensor e transforms
 * the k vector to (1-e)k and the strained kp hamiltonian is given by
 * 
 *   H(k,e) = H0((1-e)k) + D(e)
 * 
 * as the system H0 is already rotated for k and the strain is given
 * in the same coordinate system as k, i don't have to rotate the
 * second order terms.
 * 
 * the strained first order materices are:
 *   H1^T (1-e)k = ((1-e)^T H)^T k 
 * and the strained second order matrices are
 *   ((1-e)k)^T H2 (1-e)k = k^T (1-e)^T H2 (1-e)k.
 * so the strained kp matrix system is given by
 *   (1-e)^T H (1-e)
 * 
 * and we finally do not need to enter any strain dependent 2nd and 1st 
 * order terms as they can directly be calculated from the unstrained,
 * but rotated matrices.
 */
void KPMatrixBase::rotate_strain_matrix() {
	
	TDKP_ASSERT(this->initialized, "matrix is not initialized");
	
	// allocate space to temporarily store matrices
	int num;
	this->get_sparsity_pattern(num);
	vector<cplx> tmp_0th[3][3];
	for(short ee = 0; ee < 3; ee++) {
		for(short ff = 0; ff < 3; ff++) {	
			tmp_0th[ee][ff].assign(num, 0.0);
		}	
	}
	
	// rotate them
	double r_mk_r_nl;
	for(int kk = 0; kk < 3; kk++) {
		for(int ll = 0; ll < 3; ll++) {
			for(short mm = 0; mm < 3; mm++) {
				for(short nn = 0; nn < 3; nn++) {						
					r_mk_r_nl = this->rotation_matrix(mm,kk) * this->rotation_matrix(nn,ll);					
					// zero order
					for(int dd = 0; dd < num; dd++) {
						tmp_0th[mm][nn][dd] += r_mk_r_nl * this->zero_order_strain_dependent[kk][ll][dd];	
					}
				}
			}
		}	
	}
	// copy back
	for(short ee = 0; ee < 3; ee++) {
		for(short ff = 0; ff < 3; ff++) {
			// zero order
			for(int nn = 0; nn < num; nn++) {
				this->zero_order_strain_dependent[ee][ff][nn] = tmp_0th[ee][ff][nn];	
			}	
		}			
	}
	
}

/** calculate second order kp matrix for a given strain field
 *  
 * how? let (1-e)_ij = a_ij, 
 * as H2s = (1-e)^T H2 (1-e) we have  H2s_kl = sum_mn a_mk hmn anl
 * 
 * this function does not alter the kp matrix!
 */
void KPMatrixBase::build_second_order(const StrainTensor& strain_tensor, short diffop_1, short diffop_2, vector<cplx>& second_order_vec) const {

	// ---------------------------------------------------
	// check if we calculate second order contributions
	// ---------------------------------------------------
	if(Configuration::get_instance()->get("kpmatrix_ignore_second_order_strain_dependence") == 1.0) {
		second_order_vec = this->second_order_strainless[diffop_1][diffop_2];
		return;
	}

	// ----------------------------------------------------
	// reset and initialize second order target tensor
	// ----------------------------------------------------
	int num;	
	this->get_sparsity_pattern(num);	
	second_order_vec.assign(this->second_order_strainless[diffop_1][diffop_2].size(), 0.0);
	
	// -----------------------------------------------------
	// calculate (1-e)
	// -----------------------------------------------------
	double strain_transformation[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
	for(short ee = 0; ee < 3; ee++) {
		for(short ff = 0; ff < 3; ff++) {
			strain_transformation[ee][ff] -= strain_tensor.get(ee,ff);	
		}	
	}	
	// -----------------------------------------------------
	// calculate H2s_kl = amk hmn anl
	// where k is diffop_1 and l is diffop_2
	// -----------------------------------------------------
	double amk_anl;
	for(short mm = 0; mm < 3; mm++) {
		for(short nn = 0; nn < 3; nn++) {
			amk_anl = strain_transformation[mm][diffop_1] 
			        * strain_transformation[nn][diffop_2];
			for(int oo = 0; oo < num; oo++) {
				second_order_vec[oo] += amk_anl * this->second_order_strainless[mm][nn][oo];
			}
		}		
	}
}

/** calculate first order kp matrix for a given strain field 
 * 
 * here, the story is simpler,  H1s = (1-e)^T H1
 * 
 */
void KPMatrixBase::build_first_order(OperatorOrder position, const StrainTensor& strain_tensor, short diffop, vector<cplx>& first_order_vec) const {
	
	TDKP_ASSERT(position == op_left || position == op_right, "position == op_left || position == op_right"); 
	
	// ----------------------------------------------------
	// don't do complicated things if not necessary ;-)
	// ----------------------------------------------------
	if(!this->have_first_order_terms || Configuration::get_instance()->get("kpmatrix_ignore_first_order_strain_dependence") == 1.0) {
		first_order_vec = this->first_order_strainless[position][diffop];
		return;	
	}
	
	// ----------------------------------------------------
	// reset and initialize first order target tensor
	// ----------------------------------------------------
	int num;	
	this->get_sparsity_pattern(num);	
	first_order_vec.assign(this->first_order_strainless[position][diffop].size(), 0.0);
	
	// -----------------------------------------------------
	// calculate (1-e)
	// -----------------------------------------------------
	double strain_transformation[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
	for(short ee = 0; ee < 3; ee++) {
		for(short ff = 0; ff < 3; ff++) {
			strain_transformation[ee][ff] -= strain_tensor.get(ee,ff);	
		}	
	}	
	// ----------------------------------------------------
	// calculate H1s_i = (1-e)_ji Hj
	// ----------------------------------------------------
	for(short jj = 0; jj < 3; jj++) {
		for(short oo = 0; oo < num; oo++) {
			first_order_vec[oo] += strain_transformation[jj][diffop] 
			                     * this->first_order_strainless[position][jj][oo];
		}	
	}
}

void KPMatrixBase::build_zero_order(const StrainTensor& strain_tensor, const double& potential, vector<cplx>& zero_order_vec) const {
	
	int num;		
	// more complicated when strain tensor set
	if(Configuration::get_instance()->get("kpmatrix_ignore_zero_order_strain_dependence") != 1.0) {				
		zero_order_vec = this->zero_order_strainless;		
		this->get_sparsity_pattern(num);
		// eij = eji	
		double eij;
		for(short ee = 0; ee < 3; ee++) {
			// diagonal: (use that eij = eji but dont assume that Kmn = Knm
			eij = strain_tensor.get(ee,ee);
			for(int nn = 0; nn < num; nn++) {
				zero_order_vec[nn] += eij * this->zero_order_strain_dependent[ee][ee][nn];
			}
			// offdiagonal
			for(short ff = ee + 1; ff < 3; ff++) {
				eij = strain_tensor.get(ee,ff);
				for(int nn = 0; nn < num; nn++) {	
					zero_order_vec[nn] += eij * (this->zero_order_strain_dependent[ee][ff][nn] + this->zero_order_strain_dependent[ff][ee][nn]);
				}
			}
		}
	} else {			
		// no strains -> only copy strainless part
		zero_order_vec = this->zero_order_strainless;		
	}
	// add potential energy
	if(potential != 0.0e0 || this->energy_shift != 0.0) {
		this->add_zero_order_energy(potential + this->energy_shift, zero_order_vec);
	}
		
}

void KPMatrixBase::add_zero_order_energy(const double& energy, vector<cplx>& zero_order_vec) const {
	int num;
	const int* diag = this->get_diagonal_pattern(num);	
	for(int nn = 0; nn < num; nn++) {
		zero_order_vec[diag[nn]] += energy;	
	}
}

/** final kp matrix assembly of second order terms
 * 
 *  
 */
void KPMatrixBase::calculate_second_order() {
	for(short ii = 0; ii < 3; ii++) {
		for(short jj = 0; jj < 3; jj++) {
			if(this->strain_tensor_set) {
				this->build_second_order(this->strain, ii, jj, this->second_order[ii][jj]);	
			} else {
				this->second_order[ii][jj] = this->second_order_strainless[ii][jj];	
			}
		}
	}
			
}

/** final kp matrix assembly of first order matrices 
 */
void KPMatrixBase::calculate_first_order() {
	for(short ii = 0; ii < 3; ii++) {
		if(this->strain_tensor_set) {
			this->build_first_order(op_left,  this->strain, ii, this->first_order[op_left][ii]);
			this->build_first_order(op_right, this->strain, ii, this->first_order[op_right][ii]);	
		} else {
			this->first_order[op_left][ii] = this->first_order_strainless[op_left][ii];
			this->first_order[op_right][ii] = this->first_order_strainless[op_right][ii];
		}
	}
	
}

/** calculate final kp matrix zero order terms
 * 
 * write current values to "public readable" data arrays
 *
 */
void KPMatrixBase::calculate_zero_order() {
	
	// more complicated when strain tensor set
	if(this->strain_tensor_set) {
		this->build_zero_order(this->strain, this->potential_energy, this->zero_order);
	} else {
		// no strains -> only copy strainless part
		this->zero_order = this->zero_order_strainless;			
		// add potential energy
		if(this->potential_energy != 0.0e0 || this->energy_shift != 0.0) {
			int num;		
			const int* diag = this->get_diagonal_pattern(num);
			for(int nn = 0; nn < num; nn++) {
				this->zero_order[diag[nn]] += this->potential_energy + this->energy_shift;	
			}
		}
	}	
	
}

const vector<cplx>& KPMatrixBase::get_second_order(short diffop_1, short diffop_2) const {
	TDKP_BOUNDS_ASSERT(diffop_1 >= 0 && diffop_1 < 3 && diffop_2 >= 0 && diffop_2 < 3,"diffop_1 >= 0 && diffop_1 < 3 && diffop_2 >= 0 && diffop_2 < 3");
	TDKP_ASSERT(this->ready(), "kp matrix not ready");
	return this->second_order[diffop_1][diffop_2];
}

const vector<cplx>& KPMatrixBase::get_first_order(OperatorOrder position, short diffop) const {
	TDKP_BOUNDS_ASSERT(diffop >= 0 && diffop < 3, "diffop >= 0 && diffop < 3");
	TDKP_BOUNDS_ASSERT(position == op_left || position == op_right, "position == op_left || position == op_right");
	TDKP_ASSERT(this->ready(), "kp matrix not ready");
	return this->first_order[position][diffop];
}
const vector<cplx>& KPMatrixBase::get_zero_order() const {
	TDKP_ASSERT(this->ready(), "kp matrix not ready");
	return this->zero_order;
}

RMatrix<cplx>* KPMatrixBase::get_second_order_matrix(short diffop_1, short diffop_2) const {
	int num;
	const int* pattern = this->get_sparsity_pattern(num);
	return new RMatrix<cplx>(get_second_order(diffop_1,diffop_2), pattern, num);	
}
RMatrix<cplx>* KPMatrixBase::get_first_order_matrix(OperatorOrder position, short diffop) const {
	int num;
	const int* pattern = this->get_sparsity_pattern(num);
	return new RMatrix<cplx>(get_first_order(position,diffop), pattern, num);		
}
RMatrix<cplx>* KPMatrixBase::get_zero_order_matrix() const {
	int num;
	const int* pattern = this->get_sparsity_pattern(num);	
	return new RMatrix<cplx>(get_zero_order(), pattern, num);			
} 									

RMatrix<cplx>* KPMatrixBase::get_second_order_strainless_matrix(short diffop_1, short diffop_2) const {
	int num;
	const int* pattern = this->get_sparsity_pattern(num);
	return new RMatrix<cplx>(this->second_order_strainless[diffop_1][diffop_2], pattern, num);		
}
RMatrix<cplx>* KPMatrixBase::get_first_order_strainless_matrix(OperatorOrder position, short diffop) const {
	TDKP_ASSERT(position == op_left || position == op_right, "position == op_left || position == op_right");
	int num;	
	const int* pattern = this->get_sparsity_pattern(num);
	return new RMatrix<cplx>(this->first_order_strainless[position][diffop], pattern, num);		
}
RMatrix<cplx>* KPMatrixBase::get_zero_order_strainless_matrix() const {
	int num;
	const int* pattern = this->get_sparsity_pattern(num);
	return new RMatrix<cplx>(this->zero_order_strainless, pattern, num);	
}	
	
RMatrix<cplx>* KPMatrixBase::get_zero_order_strain_dependent_matrix(short strain_idx_1,short strain_idx_2) const {
	int num;
	const int* pattern = this->get_sparsity_pattern(num);
	return new RMatrix<cplx>(this->zero_order_strain_dependent[strain_idx_1][strain_idx_2], pattern, num);	
} 									

string KPMatrixBase::get_nonellipticity_warning(const double& degree_of_nonellipticity) {
		
	string degree_of_nonellipticity_string;

	if(degree_of_nonellipticity == 1.0) {
		degree_of_nonellipticity_string = " (your input values are complete bullshit)";		
	} else if(degree_of_nonellipticity == 0.0) {
		degree_of_nonellipticity_string = " (perfect elliptic)";	
	} else if(degree_of_nonellipticity < 0.05) {
		degree_of_nonellipticity_string = " (usual acceptable)";	
	} else if(degree_of_nonellipticity < 0.1) {
		degree_of_nonellipticity_string = " (may still work)";	
	} else if(degree_of_nonellipticity < 0.2) {
		degree_of_nonellipticity_string = " (looks quite bad)";	
	} else {
		degree_of_nonellipticity_string = " (change your parameters!)";	
	}
	
	return degree_of_nonellipticity_string; 
}  


vector<cplx> KPMatrixBase::evaluate_at(const Vector3D& kvec) const {
	
	cplx i(0.0, 1.0);
	
	vector<cplx> ret(this->get_zero_order().size(), 0.0);
	for(unsigned short ii = 0; ii < 3; ii++) {
		for(unsigned short jj = 0; jj < 3; jj++) {
			const vector<cplx>& rval = this->get_second_order(ii,jj);
			for(unsigned int ss = 0; ss < ret.size(); ss++) {
				ret[ss] += i * kvec(ii) * i * kvec(jj) * (- rval[ss]);
			}	
		}	
		for(unsigned short oo = 0; oo < 2; oo++) {
			const vector<cplx>& rval = this->get_first_order(oo,ii);
			for(unsigned int ss = 0; ss < ret.size(); ss++) {
				ret[ss] += (- i) * kvec(ii) * rval[ss];	
			}
		}
		const vector<cplx>& rval = this->get_zero_order();
		for(unsigned int ss = 0; ss < ret.size(); ss++) {
			ret[ss] += rval[ss];
		}
	}
	return ret;
}

/** determine in-plane rotation symmetry of a kp matrix
 * 
 * this function check a passed rotation matrix (that should initially be rotated to a given rotation system)
 * for rotation symmetries in the yz plane of the rotated system
 * 
 * we only check the final matrices, so the ones that include strain
 * 
 * @param matrix where we are going to perform tests
 * @param initial_rotation initial rotation matrix, so all tests are Rtest * Rinitial
 * @param max_n the maximum Cn symmetry we will test   
 */
void KPMatrixBase::test_Cn_plane_symmetry(unsigned int num_tests, unsigned int max_n) {
	
	const double threshold = 1.0e-10;

	// init random number generator
	timeval* tp = new timeval();
	gettimeofday(tp,NULL);
	srand48(tp->tv_usec);
	delete tp; tp = NULL;
	// create arbitrary direction
	vector<Vector3D> testing_dirs;
	testing_dirs.push_back(Vector3D(0.0, 1.0, 0.0));
	testing_dirs.push_back(Vector3D(0.0, 0.0, 1.0));
	testing_dirs.push_back(Vector3D(0.0, 1.0, 1.0));			
	double length = 2.0;		
	for(unsigned int ii = 0; ii < num_tests; ii++) {
		testing_dirs.push_back(
			Vector3D(
				0,
				length * (0.5 - drand48()),
				length * (0.5 - drand48())
			)
		);
	}
	
	// calculate reference hamiltonians for testing directions
	vector<vector<cplx> > reference_hamiltonians;
	for(unsigned int vv = 0; vv < testing_dirs.size(); vv++) {
		reference_hamiltonians.push_back(this->evaluate_at(testing_dirs[vv]));	
	}

	// --------------------------------------------		
	// for all Cn groups
	// --------------------------------------------	
	vector<bool> symmetry_groups(max_n, false);
	int group_counter = 0;
	ostringstream sout;
	for(unsigned int nn = 2; nn < max_n; nn++) {	
		sout << "checking yz plane with operations in symmetry group C" << nn << ":\n";
		// generating angle
		const double generating_angle = 2.0 * constants::pi / static_cast<double>(nn);
		bool all_rotations_give_equal = true;									
		// for all rotations in the Cn group		
		for(unsigned int ii = 1; ii < nn; ii++) {
			// create current angle
			double phi = generating_angle * static_cast<double>(ii); 
			// create rotation matrix
			RMatrix<double> rotation(3,3);
			rotation(0,0) = 1.0;		
			rotation(1,1) = cos(phi);
			rotation(1,2) = - sin(phi);
			rotation(2,1) = sin(phi);
			rotation(2,2) = cos(phi);			
			int good = 0;
			int bad  = 0;
			// for all testing dirs
			for(unsigned int vv = 0; vv < testing_dirs.size(); vv++) {				
				Vector3D dir = rotation * testing_dirs[vv];
		/*		cout << "testing:\nR = " 
				     << rotation
				     << "\ninitial dir: " 
				     << testing_dirs[vv]
				     << "\nrotated:\n"
				     << dir << "\n";*/
				vector<cplx> res = this->evaluate_at(dir);
				TDKP_ASSERT(res.size() == reference_hamiltonians[vv].size(), "res.size() == reference_hamiltonians[vv].size()"); 
				const vector<cplx>& ref = reference_hamiltonians[vv];
				bool is_good = true;
				for(unsigned int ss = 0; ss < res.size(); ss++) {
					if(tdkp_math::abs(res[ss] - ref[ss]) > threshold) {
						cout << "diff at " << ss << " res: " << res[ss] << " ref: " <<  ref[ss] << " for " << dir << "\n";						
						is_good = false; break;						
					}					
				}
				if(is_good) {
					good++;	
				} else {
					bad++;
				}
			}
			if(bad > 0) {
				all_rotations_give_equal = false;	
			}
			sout << "  r" << ii << ": good = " << good << ", bad = " << bad << "\n";																							
		}
		symmetry_groups[nn] = all_rotations_give_equal;
		if(all_rotations_give_equal) {
			group_counter++;	
		}
							
	}
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());		
	
	// -----------------------------------------------
	// output
	// -----------------------------------------------
	sout.str("");	
	if(group_counter > 0) {
		sout << "the kp matrix has the following symmetry groups:";
		for(unsigned int nn = 2; nn < max_n; nn++) {
			if(symmetry_groups[nn]) {
				sout << " C" << nn; 	
			}
		}
		Logger::get_instance()->emit(LOG_INFO, sout.str());
	} else {
		sout << "could not find any symmetry group for the matrix!";
		Logger::get_instance()->emit(LOG_INFO, sout.str());	
	}	
} 

// --------------------------------------------------------------
// implementation of KPMatrixWurtziteBase
// --------------------------------------------------------------

double KPMatrixWurtziteBase::get_nonellipticity_ratio(const vector<double>& eigenvalues) {
	double pos_val = 0.0;
	double neg_val = 0.0;
	for(unsigned int ii = 0; ii < eigenvalues.size(); ii++) {
		if(eigenvalues[ii] < 0.0) {
			neg_val += eigenvalues[ii]; 	
		} else {
			pos_val += eigenvalues[ii];	
		}
	}	
	if(neg_val == 0.0) {
		return 1.0;	
	}
	return tdkp_math::abs(pos_val / neg_val);
}

/** calculate the ellipticity of the coupled second order differential operator 
 * 
 * damn, if you don't know what i'm talking about, check my paper ;-)
 * */
void KPMatrixWurtziteBase::determine_ellipticity_eigenvalues(const WurtziteEffectiveMassParams& params, vector<double>& eigenvalues) {



	// ------------------------------------------------------
	// set params to short variable names
	// ------------------------------------------------------
	const double& A1  = params.A1;
	const double& A2  = params.A2;
	const double& A3  = params.A3;
	const double& A4  = params.A4;
	const double& A5  = params.A5;
	const double& A5m = params.A5m;
	const double& A6m = params.A6m;	
	const double  A5p = params.A5 - A5m;
	const double  A6p = params.A6 - A6m;
	 			

	RMatrix<cplx> matrix(9,9);
	const complex<double> i(0,1);	
	
	matrix(0,0) =   A2 + A4;
	matrix(0,1) = - A5;
	matrix(0,3) =   i * (A5m - A5p);
	matrix(0,4) = - i * (A5m + A5p);
	matrix(0,8) = - A6p;
	matrix(1,0) = - A5;
	matrix(1,1) =   A2 + A4;
	matrix(1,3) =   i * (A5m + A5p);
	matrix(1,4) = - i * (A5m - A5p);
	matrix(1,8) =   A6p;
	matrix(2,2) =   A2;
	matrix(2,6) = - A6m;
	matrix(2,7) =   A6m;
	matrix(3,0) = - i * (A5m - A5p); 
	matrix(3,1) = - i * (A5m + A5p);
	matrix(3,3) =   A2 + A4;
	matrix(3,4) =   A5;
	matrix(3,8) = - i * A6p;
	matrix(4,0) =   i * (A5m + A5p);
	matrix(4,1) =   i * (A5m - A5p);
	matrix(4,3) =   A5;
	matrix(4,4) =   A2 + A4;
	matrix(4,8) = - i * A6p;
	matrix(5,5) =   A2;
	matrix(5,6) =   i * A6m;
	matrix(5,7) =   i * A6m;
	matrix(6,2) = - A6m;
	matrix(6,5) = - i * A6m;
	matrix(6,6) =   A1 + A3;
	matrix(7,2) =   A6m;
	matrix(7,5) = - i * A6m;
	matrix(7,7) =   A1 + A3;
	matrix(8,0) = - A6p;
	matrix(8,1) =   A6p;
	matrix(8,3) =   i * A6p;
	matrix(8,4) =   i * A6p;
	matrix(8,8) =   A1; 
	  
	TDKP_ASSERT(RMatrix<cplx>::hermitian(matrix), "matrix is not hermitian!");
	vector<cplx> cplx_eigenvalues;
	vector<cplx> eigenvectors;  	
	RMatrix<cplx>::get_eigensystem(matrix, cplx_eigenvalues, eigenvectors);
	TDKP_ASSERT(cplx_eigenvalues.size() == 9, "cplx_eigenvalues.size() == 9");
	eigenvalues.resize(9);
	for(unsigned int ii = 0; ii < cplx_eigenvalues.size(); ii++) {
		eigenvalues[ii] = cplx_eigenvalues[ii].real();
		TDKP_ASSERT(tdkp_math::abs(cplx_eigenvalues[ii].imag()) < 1.0e-12, "eigenvalues are not real!!!!");			
	}
	
}

/** brute force operator optimization
 * 
 * just vary A5m and A6m until ellipticity is reached or  
 * values don't get better ...
 */
 
WurtziteEffectiveMassParams 
KPMatrixWurtziteBase::determine_optimal_splitting(const WurtziteEffectiveMassParams& params)  {
			
				
			
	// -----------------------------------------------
	// start at symmetric splitting
	// -----------------------------------------------
	WurtziteEffectiveMassParams final_params = params;	
	final_params.A5m = final_params.A5 / 2.0;
	final_params.A6m = final_params.A6 / 2.0;
	
	WurtziteEffectiveMassParams tmp_params   = final_params;
	
	double dA5 = tdkp_math::abs(final_params.A5) / 10.0;
	double dA6 = tdkp_math::abs(final_params.A6) / 10.0; 
	vector<double> eigenvalues;
	
	// init first			
	determine_ellipticity_eigenvalues(tmp_params, eigenvalues);
	double current = get_nonellipticity_ratio(eigenvalues);
			
	int max_iterations = 200;
	double stop        = 1.0e-4;
	// iterate
	while(max_iterations-- > 0) {
		
		//cout << "cur: " << setw(15) << current << " A5m: " << setw(15) << final_params.A5m << " A6m: " << setw(15) << final_params.A6m << "\n";
		
		// stop if current is zero
		if(current < stop) {
			break;	
		}		
		// calculate + dA5		
		tmp_params = final_params; tmp_params.A5m += dA5;
		determine_ellipticity_eigenvalues(tmp_params, eigenvalues);
		double dA5_plus = get_nonellipticity_ratio(eigenvalues);
		// calculate - dA5	
		tmp_params = final_params; tmp_params.A5m -= dA5;	
		determine_ellipticity_eigenvalues(tmp_params, eigenvalues);
		double dA5_minus = get_nonellipticity_ratio(eigenvalues);						
		// calculate + dA6		
		tmp_params = final_params; tmp_params.A6m += dA6;
		determine_ellipticity_eigenvalues(tmp_params, eigenvalues);
		double dA6_plus = get_nonellipticity_ratio(eigenvalues);
		// calculate - dA6		
		tmp_params = final_params; tmp_params.A6m -= dA6;
		determine_ellipticity_eigenvalues(tmp_params, eigenvalues);
		double dA6_minus = get_nonellipticity_ratio(eigenvalues);		

		// ---------------------------------------------
		// check if the wobbling is successfull ....
		// ---------------------------------------------
		double best = min(dA5_plus, min(dA5_minus, min(dA6_plus, dA6_minus)));
				 		
		if(best < current) {
			// yep, there is an update
			if(best == dA5_plus) {
				final_params.A5m += dA5;	
			} else if(best == dA5_minus) {
				final_params.A5m -= dA5;
			} else if(best == dA6_plus) {
				final_params.A6m += dA6;
			} else {
				final_params.A6m += dA6;
			}
			current = best;
		} else {
			dA6 /= 2.0;
			dA5 /= 2.0;
		}						
	}	
	if(max_iterations == 0.0) {
		ostringstream sout;
		sout << "could not find an elliptic operator ordering. terminated after 200 iterations with non-ellipticity estimate of " << current;
		Logger::get_instance()->emit(LOG_WARN, sout.str());	
	}
	ostringstream sout;
	sout << "determined optimal parameter splitting:\n"
	     << "A5:              " << setw(5) << final_params.A5 << "\n"
	     << "A5+:             " << setw(5) << final_params.A5 - final_params.A5m << "\n"
	     << "A5-:             " << setw(5) << final_params.A5m << "\n"
	     << "A6:              " << setw(5) << final_params.A6 << "\n"
	     << "A6+:             " << setw(5) << final_params.A6 - final_params.A6m << "\n"
	     << "A6-:             " << setw(5) << final_params.A6m << "\n"
	     << "non-ellitpicity: "	<< setw(5) << current;
		     
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
		      
	return final_params;
}





} // end of namespace
