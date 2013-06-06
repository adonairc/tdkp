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

#ifndef INTRINSICSTRAIN_H_
#define INTRINSICSTRAIN_H_

#include "tdkp/probdefs/LinearProblem.h"

namespace tdkp {

/** calculate strain field due to lattice mismatch */
class IntrinsicStrain : public tdkp::LinearProblem<double>
{
public:
	IntrinsicStrain(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~IntrinsicStrain();

	virtual void set_reference_lattice_constant(double base_lattice);
	void set_strained_axes(short which, bool strained_or_unstrained) { strained_axes[which] = strained_or_unstrained; }	
	StrainField* get_strain_field() const;
	void build_faky_solution();
	// -----------------------------------------------------------
	// for communication between FEMSolver and problem class
	// -----------------------------------------------------------
	virtual const int* get_node_sparsity_pattern(int &num) const;
	virtual void calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* node_internal_indices, int &n) const;   	
	virtual void prepare();
	virtual string get_unique_identifier() const { return string("IntrinsicStrain"); }
	
	void display_solution_info() const {};
	
	void dump();
							
protected:
	void evaluate_strain(const double derivatives[3][3], const RMatrix<double>& intrinsic_strains, StrainTensor& strain_tensor) const;
	void set_sparsity_pattern() throw(Exception*);
	virtual void prepare_initial_intrinsic_strains();

	class region_properties {
		public:
			RMatrix<double> intrinsic_stress_load; // this is not the intrisic stress, but the intrinsic stress load, means: load[ii] = load[ii][jj] * single_1st[jj]
			RMatrix<double> intrinsic_strain;
			double coefficent_matrices[3][3][3][3];
			region_properties() : intrinsic_stress_load(3,3), intrinsic_strain(3,3) {}			
	};			
	
	
	bool   strained_axes[3];
	double reference_lattice_constant;

	int  sparsity_length;
	int* sparsity_copy;				
	vector<region_properties> properties;
	vector<RMatrix<double> >  initial_intrinsic_strains;
	
};

/** anisotropic strain class for the wurzite crystal 
 * 
 * to the present, only strains in principal directions can 
 * be treated. e.g. for a 3D problem, the real space coordinate axes
 * must match the system coordinate axes (up to permutatiions)
 * therefore, one is able to calculate strain profiles in quantum wells
 * in (c) and (m) plane, but not in (r) or (n) plane
 * (see H. Morkoc, Nitride Semiconductors and Devices, Springer, 1999, p 10)
 * 
 * the problem is as follows:
 *    if R denotes the mapping between the coordinate x and the crystal axis
 *    e.g. cx = Rx, then the strain maps as ce = R e R^t and 
 *    e = R^t ce R.
 *    the strain energy is given by 
 *    e_ij C_ijkl e_kl
 *    so, the strain energy in the rotated system is obtained via back-
 *    transforming the rotated strain to the initial system  
 *    (R^t ce R)_ij  Cijkl (R^t ce R)_kl
 * 
 * 	  this one can be rewritten to give a rotated elastic tensor B_ijkl
 *    so that   ce_ij B_ijkl ce_kl gives the strain energy.
 *   
 *    unfortunately, the tensor B_ijkl gets negative values, and values
 *    at e.g. B_1112 or B2223 (coupling e.g. exx and exy)
 * 
 * 	  as the derivation was lengthly already i skipped it for the moment 
 * 
 * 	so if you are interested in (r) or (n) plane wells, you may do the 
 *  following:
 * 		create a 2D structure, where axis 0 is in (a/m)-plane and axis 
 *      1 is c-axis, make a diagonal through the 2D structure and let 
 *      it relax, extract the strain profile and create a 
 *      1D profile out of it ... ;-)
 * 
 *  considering the axes of the wurzite crystal, we use the convention
 *  the axis 0,1 are xy (a/m plane) and axis 2 is c-axis.
 *  this choice can be permuted using the set_axes command  
 * 
 *  
 * 
 */
class IntrinsicStrainWurzite : public IntrinsicStrain {
	
public:	
	IntrinsicStrainWurzite(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual void prepare();
	virtual void set_reference_lattice_constant(double base_lattice);
	virtual void set_reference_lattice_constant(const double& lattice_a, const double& lattice_c);
	void set_axes(short XX, short YY, short ZZ); // we only allow permutations
	
	virtual string get_unique_identifier() const { return string("IntrinsicStrainWZ"); }
	
protected:	
	virtual void prepare_initial_intrinsic_strains();
	
private:

	double reference_lattice_a;
	double reference_lattice_c;

	short DIR_XX, DIR_YY, DIR_ZZ;	
	
				
};


} // end of namespace

#endif /*INTRINSICSTRAIN_H_*/
