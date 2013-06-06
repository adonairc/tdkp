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

#ifndef KPBASE3D_H_
#define KPBASE3D_H_

#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/kpmatrices/KPMatrixBase.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/probdefs/EigenProblem3D.h"
#include "tdkp/probdefs/KPBase.h"


namespace tdkp {
	

typedef KPBase<SchroedingerProblem<EigenProblem3D<complex<double>, complex<double> > > > KPBase3DParent;

class KPBase3D : public KPBase<SchroedingerProblem<EigenProblem3D<complex<double>, complex<double> > > > {
		
public:
	KPBase3D(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~KPBase3D();
	const int* get_node_sparsity_pattern(int &num) const;
	virtual void calculate_element_matrices(const Element* elem, cplx* lhs, double *rhs, int* node_internal_indices, int &n) const;   	
	virtual void prepare();
	virtual void solve(int num_solutions);	
	virtual void set_energy_guess(double energy);	
	virtual void drop_energy_guess();	
	void         set_axes(const Vector3D& kx, const Vector3D& ky, const Vector3D& kz);
//	void		 copy_kp_matrices(vector<KPMatrixBase*>& copy_into) const;	
	virtual EigenProblemType get_problem_type() const;
	virtual void add_solution(cplx solution_value, const cplx* solution_vector, int length);
	
	virtual double get_energy_shift() const;
			
protected:	
	virtual KPMatrixBase* get_matrix() const = 0;
	void   init_3D();
	void   delete_3D();	
	bool   first_order_terms_exist;	
	
	bool   energy_guess_set;
	double energy_guess;
						
	
								
};


	
}

#endif /*KPBASE3D_H_*/
