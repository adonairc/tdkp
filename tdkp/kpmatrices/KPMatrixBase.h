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

#ifndef KPMATRIXBASE_H_
#define KPMATRIXBASE_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/main/Fields.h"

#include <vector>

using namespace std;

namespace tdkp {


/** base class for kp matrices. 
 * 
 * the kp matrices are split into second, first and zero order terms of k
 * they basically model the relation between the bands. 
 * \f$ H = \sum_{i,j} M_{ij}^{(2)}k_{i}k_{j} + \sum_{l}M_{l}^{(1)} k_{l} + M^{0} \f$
 */
class KPMatrixBase {

				
public:

	static const int op_left  = 0;
	static const int op_right = 1;
	typedef int OperatorOrder;

	KPMatrixBase();
	KPMatrixBase(const Material*);
	virtual ~KPMatrixBase();			
	void  init_from_kpmatrix(const KPMatrixBase& base);
	void  set_rotation(const RMatrix<double>& rotation_matrix);
	void  set_material(const Material*);
	void  set_strains(const StrainTensor& tensor) throw(Exception*); 
	void  set_potential(double);
	void  set_energy_shift(double);	
	void  set_material_name(const string& name); 
	const string& get_material_name() const;
	void  calculate() throw(Exception*);
	void  enforce_recalculation();	
	bool  ready() const;
	void  set_output_surpression(bool surpress_output_) { surpress_output_to_user = surpress_output_; }
	bool  surpress_output() const { return surpress_output_to_user; }
	void  dump(const char* filename) const;
	 	
	virtual const int* get_sparsity_pattern(int& length) const = 0;		
	virtual const int* get_diagonal_pattern(int& length) const = 0; 
	virtual int        get_number_of_bands() const = 0;
	virtual bool       check_solution_type_available(KPSolutionType type) const = 0;
	
	// ------------------------------------------------------
	// cache-less assembly
	// ------------------------------------------------------	
	void build_second_order(const StrainTensor& tensor, short diffop_1, short diffop_2, vector<cplx>& copy_to) const;
	void build_first_order (OperatorOrder position, const StrainTensor& tensor, short diffop, vector<cplx>& copy_to) const;	
	void build_zero_order  (const StrainTensor& tensor, const double& potential, vector<cplx>& copy_to) const;
	void add_zero_order_energy(const double& energy, vector<cplx>& zero_order_vec) const;
	
	// ------------------------------------------------------
	// build strained first order momentum contributions
	// ------------------------------------------------------
	// TODO: really create a function for that
			
	const vector<cplx>&   get_second_order(short diffop_1, short diffop_2) const;
	const vector<cplx>&   get_first_order(OperatorOrder position, short diffop) const;
	const vector<cplx>&   get_zero_order() const; 
	
	RMatrix<cplx>* get_second_order_matrix(short diffop_1, short diffop_2) const;
	RMatrix<cplx>* get_first_order_matrix(OperatorOrder position, short diffop) const;
	RMatrix<cplx>* get_zero_order_matrix() const; 	
	
	RMatrix<cplx>* get_second_order_strainless_matrix(short diffop_1, short diffop_2) const;
	RMatrix<cplx>* get_first_order_strainless_matrix(OperatorOrder position, short diffop) const;
	RMatrix<cplx>* get_zero_order_strainless_matrix() const; 		
	
	RMatrix<cplx>* get_zero_order_strain_dependent_matrix(short strain_idx_1,short strain_idx_2) const;
	
	void test_Cn_plane_symmetry(unsigned int num_tests, unsigned int max_n);
	
	vector<complex<double> > evaluate_at(const Vector3D& kvec) const; 
	
protected:
	// no copy constructor
	KPMatrixBase(const KPMatrixBase&) { TDKP_GENERAL_EXCEPTION("don't copy the kp matrix ..."); }
	virtual void init_base_and_strain_matrix() = 0;			
	void allocate_matrix_space(); 
	void init();
	static string get_nonellipticity_warning(const double& degree_of_nonellipticity);
		
	/** whether we have first order terms */
	bool have_first_order_terms;				
	/** current strain tensor */
	StrainTensor   strain;		
	/** indicate whether strain tensor is set */
	bool strain_tensor_set;		
	/** indicate whether strain dependence is available (all parameters are available) */
	bool strain_dependence_available; 
	/** potential energy (e.g. if you do schroedinger-poisson) */
	double potential_energy;	
	/** booleans to indicate that essential variables have changed -> information for update_strategy */

	/** the \f$M_{ij}^{(2)} \f$ 
	 * 
	 * val[ii][jj][ss] 
	 * ii = left diff op
	 * jj = right diff op
	 * ss = sparsity pattern index
	 */
	vector<cplx> second_order[3][3];	
	/** the first order matrix with [left|right][D_DX|D_DY|D_DZ] */	
	vector<cplx> first_order[2][3];
	/** the \f$M^{(0)} \f$ */
	vector<cplx> zero_order;
	/** the \f$M_{ij}^{(2)} without strain \f$ */
	vector<cplx> second_order_strainless[3][3];
	/** the \f$M_{l}^{(1)} without strain \f$ */	
	vector<cplx> first_order_strainless[2][3];
	/** the \f$M^{(0)} without strain \f$ */
	vector<cplx> zero_order_strainless;

	/** the \f$M^{(0)} strain dependent part \f$ */
	vector<cplx> zero_order_strain_dependent[3][3];	
	
	/** general energy shift */
	double energy_shift;	
	/** pointer to material object containing the material properties */
	const Material* material;	
	bool initialized;
		
private:
	void rotate_base_matrix();
	void rotate_strain_matrix();
	void calculate_second_order();
	void calculate_first_order();
	void calculate_zero_order();		
		
	/** the rotation matrix calculated by KPMatrixBase::set_rotation */
	RMatrix<double> rotation_matrix;	
	/** boolean to indicate if rotation has been applied */
	bool   rotated;
	bool   rotation_changed;
	/** boolean to indicate whether strain tensor has changed */
	bool   strain_tensor_changed;	
	bool   potential_energy_changed;		
	bool   material_changed;	
	
	string material_name;
	/** bool whether the output should be surpressed */ 
	bool surpress_output_to_user;
	
};

	ostream& operator<<(ostream& out, const KPMatrixBase& mat);

struct WurtziteEffectiveMassParams {
	WurtziteEffectiveMassParams() { A1 = A2 = A3 = A4 = A5 = A6 = A5m = A6m = mc_xxyy = mc_zz = 0.0; }
	double A1, A2, A3, A4, A5, A5m, A6, A6m, mc_xxyy, mc_zz;
};

/** abstract wurzite class
 *
 * implements some general functions to be used for 6x6 and 8x8 Wurzite models 
 */	
class KPMatrixWurtziteBase : public KPMatrixBase {
public:

	virtual ~KPMatrixWurtziteBase() {}
	
	static void determine_ellipticity_eigenvalues(const WurtziteEffectiveMassParams& params, vector<double>& eigenvalues);		
	static double get_nonellipticity_ratio(const vector<double>& eigenvalues); 
	
	static WurtziteEffectiveMassParams determine_optimal_splitting(const WurtziteEffectiveMassParams& initial);	
	
protected:
	KPMatrixWurtziteBase() {} // not a public class		
	
};	
	

} // end of  namespace

#endif /*KPMATRIXBASE_H_*/
