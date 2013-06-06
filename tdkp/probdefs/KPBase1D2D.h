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

#ifndef KPBASE1D2D_H_
#define KPBASE1D2D_H_

#include "tdkp/common/Domain.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/main/Bandstructure.h"
#include "tdkp/main/FEMSolverGEVP.h"
#include "tdkp/probdefs/EigenProblem.h"
#include "tdkp/probdefs/KPBase.h"

namespace tdkp {

// ----------------------------------------
// forward declarations for schroedinger PMLs
// ----------------------------------------
class KP4x41D2D;
class KP6x61D2D;
class KP8x81D2D;
class KP6x61D2DWZ;
class KP8x81D2DWZ;

template<class T> class SchroedingerPML;
template<class T> class SchroedingerPML1D2D;

typedef KPBase<SchroedingerProblem<EigenProblem<complex<double>, complex<double>, double> > > KPBase1D2DParent;

class KPBase1D2D : public KPBase<SchroedingerProblem<EigenProblem<complex<double>, complex<double>, double > > > {
		
public:
	KPBase1D2D(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~KPBase1D2D();

	// ----------------------------------------------
	// axis setting
	// ----------------------------------------------	
 	void           set_axes(const Vector3D& k_transversal, const Vector3D& k_x, const Vector3D& k_y);
 	void           set_axes(const Vector3D& k_transversal, const Vector3D& k_q);
			
	// ----------------------------------------------
	// solving
	// ----------------------------------------------	
	virtual void   prepare();
	/** single dispersion solving process */
	virtual void   solve(int num_subbands, double kmin, double kmax, int num_k_values);
	/** solve kp equations for a domain */
	virtual void   solve(int num_subbands, const DomainMaster& domain);
	virtual void   solve(int) { Logger::get_instance()->emit(LOG_ERROR,"you have to use overloaded solve defining the k-space range"); }	
	virtual void   display_solution_info() const;

	// ----------------------------------------------
	// results
	// ----------------------------------------------	
	const BandstructureDomain<complex<double> >& get_bandstructure(int idx = -1) const;
	virtual void delete_solutions();
		
	// ----------------------------------------------
	// energy guessing and control
	// ----------------------------------------------
	virtual void   set_energy_guess(unsigned int kidx, double energy);
	virtual void   remove_energy_guess(unsigned int kidx);
	virtual void   remove_all_energy_guesses();
	void           set_energy_barrier(double energy);
	void           remove_energy_barrier();
	void		   set_upper_energy_limit(const double& upper_limit);		
	virtual void   set_target_energy_offset(double energy_offset);
	void 		   auto_prepare_energy_guess(const DomainMaster& domain); 	
 		
	// ----------------------------------------------
	// intercommunication with FEMSolverGEVP
	// ----------------------------------------------
	const int*     get_node_sparsity_pattern(int &num) const;
	virtual void   calculate_element_matrices(const Element* elem, cplx* lhs, double *rhs, int* node_internal_indices, int &n) const;
	virtual void   add_solution(cplx solution_value, const cplx* solution_vector, int length);	
	virtual EigenProblemType get_problem_type() const;

	
	
protected:
	
	virtual double        get_energy_shift() const;
	virtual KPMatrixBase* get_matrix() const = 0;
	virtual double        get_minimum_bandedges() const;
	virtual bool          use_user_defined_energy_guess() const;

	// ----------------------------------------------
	// functions called during solve process
	// ----------------------------------------------
	void solving_preinform_user(int num_subbands, const DomainMaster& domain) const;
	DomainMaster build_kspace_domain(int num_subbands, double kmin, double kmax, int num_k_values) const;
		
	void   init_base1D2D();
	void   delete_base1D2D();	
	void   add_cache_to_bandstructure_object(Bandstructure<cplx>* band, int kidx);
			
	bool   first_order_terms_exist;
	bool   first_order_ignore_diagonal;
	
	vector<bool>   energy_guess_set;
	vector<double> energy_guess;
	double energy_of_last_solution;
	double energy_offset;
	
	double energy_barrier;
	bool   energy_barrier_set;
	
	double upper_energy_limit;
	bool   upper_energy_limit_set;
	
	/** current transversal k value */
	double          k_transversal;
	/**  current kidx when calculation bandstructre (set by solve) */
	unsigned int    k_idx_current;
	
  	vector<vector<cplx> > solution_cache; 		
	vector<BandstructureDomain<cplx>* > bandstructures;

private:
	static void *parallel_controller_solve(void *arg);

	// ---------------------------------------------
	// friend declarations
	// ---------------------------------------------
	friend class SchroedingerPML1D2D<KP4x41D2D>;
	friend class SchroedingerPML1D2D<KP6x61D2D>;
	friend class SchroedingerPML1D2D<KP8x81D2D>;
	friend class SchroedingerPML1D2D<KP6x61D2DWZ>;
	friend class SchroedingerPML1D2D<KP8x81D2DWZ>;
	
	friend class SchroedingerPML<KP4x41D2D>;
	friend class SchroedingerPML<KP6x61D2D>;
	friend class SchroedingerPML<KP8x81D2D>;
	friend class SchroedingerPML<KP6x61D2DWZ>;
	friend class SchroedingerPML<KP8x81D2DWZ>;


};

}

#endif /*KPBASE1D2D_H_*/
