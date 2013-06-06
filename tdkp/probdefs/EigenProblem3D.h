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

#ifndef EIGENPROBLEM3D_H_
#define EIGENPROBLEM3D_H_

#include "tdkp/common/all.h"
#include "tdkp/probdefs/EigenProblem.h"
#include "tdkp/main/Bandstructure.h"

using namespace std;

namespace tdkp {

template<class T, class TMat>
class EigenProblem3D : public tdkp::EigenProblem<T, TMat, double> {
public:

	EigenProblem3D(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~EigenProblem3D();

	// -----------------------------------------------------
	// result access functions
	// -----------------------------------------------------
	virtual  void   display_solution_info() const;
	void 	 delete_solutions();
	int  	 num_solutions() const;
	const T& get_eigenvalue(unsigned int idx) const;
	T& get_eigenvalue(unsigned int idx);

	// -------------------------------------------------------
	// result extraction
	// -------------------------------------------------------
	const BandstructureDomain<complex<double> >& get_bandstructure() const { return solution_container; }

	// -------------------------------------------------------
	// TODO: get rid of these functions
	// -------------------------------------------------------
	EigenSolution<double>* 	   get_probability_object(unsigned int idx) const;
	EigenSolution<T>*          get_solution_object(unsigned int idx) const;
	EigenSolutionSet<double>*  get_probability_objects() const;
	EigenSolutionSet<T>*       get_solution_objects() const;

	// -------------------------------------------------------
	// (intercommunication with FEMSolver during solving)
	// -------------------------------------------------------
	virtual void add_solution(T solution_value, const T* solution_vector, int length);

	/** update bandstructure domain object with calculated solutions */
	void update_bandstructure_container();

protected:
	
	vector<T> solution_value;
private:
	vector<vector<T> >  solution_vector;
	BandstructureDomain<cplx> solution_container;

	const vector<T>& get_raw_eigenvector(unsigned int idx) const;
	double*  get_probability_vector(unsigned int idx) const;


};

template<class T, class TMat>
EigenProblem3D<T, TMat>::~EigenProblem3D() {
	// stl stuff is destroyed automatically
}

template<class T, class TMat>
EigenProblem3D<T, TMat>::EigenProblem3D(const Geometry& geometry_, MaterialDatabase& material_database_)
: EigenProblem<T,TMat,double>(geometry_, material_database_)
{
}

/** add solution from eigenvalue calculation in internal ordering
 *
 * @param value the eigenvalue of the solution
 * @param vector pointer to the eigenvector(will be copied) in the original order
 * @param length length of the eigenvector (used for consistency check
 */
template<class T, class TMat>
void EigenProblem3D<T,TMat>::add_solution(T value, const T* vector, int length) {
	     	 	     
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "EigenProblem3D: added solution " << this->solution_value.size() << " with ev: " << value);	     	 	     
	     	 	     
	this->solution_value.push_back(value);
	// increase solution vector size
	this->solution_vector.resize(this->solution_vector.size() + 1);
	std::vector<T>& target = this->solution_vector.back();
	target.resize(length);
	for(int ii = 0; ii < length; ii++) {
		target[ii] = vector[ii];
	}
	this->normalize_solution(&target[0]); // normalize

}

// ------------------------------------------------------
// result access
// -------------------------------------------------------
template<class T, class TMat>
void EigenProblem3D<T,TMat>::delete_solutions() {
	this->solution_vector.resize(0);
	this->solution_value.resize(0);
	this->update_bandstructure_container();
}

/** @return number of solutions available */
template<class T, class TMat>
int  EigenProblem3D<T,TMat>::num_solutions() const {
	return (int)this->solution_vector.size();
}
/** return the raw eigenvector, basically the data that was set via add_solution
 * @param idx the index of the solution [0, num_solutions() [
 * @return pointer to the eigenvector
 * */
template<class T, class TMat>
const vector<T>& EigenProblem3D<T,TMat>::get_raw_eigenvector(unsigned int idx) const {
	if(this->solution_vector.size() > idx) {
		return this->solution_vector[idx];
	} else {
		TDKP_GENERAL_EXCEPTION("solution index out of bounds");
	}
}

/** returns a new created vector containing the amplitude of the requested band at the nodes
 *
 * the ordering in the vector corresponds to the global node index
 * @param idx solution index [0, num_solutions() [
 * @param band band index [0, get_num_equations_per_node]
 * @return new vector with amplitude values at nodes
 */
 /*
template<class T, class TMat>
T* EigenProblem3D<T,TMat>::get_band_vector(unsigned int idx, unsigned int band) const {
	if(this->solution_vector.size() <= idx) {
		TDKP_GENERAL_EXCEPTION("solution index out of bounds");
	}
	if((signed)band >= this->get_num_equations_per_node()) {
		TDKP_GENERAL_EXCEPTION("band index ouf of bounds");
	}
	typename Geometry::node_iterator it;
	unsigned int num_eq   = this->get_num_equations_per_node();
	unsigned int num_node = this->geometry->get_num_nodes();
	T* res 				  = new T[num_node];
	T* values             = this->solution_vector[idx];
	// loop over nodes and set zero on dirichlet nodes and band value to other
	for(it = this->geometry->nodes_begin(); it != this->geometry->nodes_end(); it++) {
		if((*it)->get_index_internal() == -1) {
			res[(*it)->get_index_global()] = 0.0;
		} else {
			res[(*it)->get_index_global()] = values[(*it)->get_index_internal() * num_eq + band];
		}
	}
	return res;
}*/

/** returns a newly created vector containing the probability of the requested band at the nodes
 *
 * the ordering in the vector corresponds to the global node index
 * @param idx solution index [0, num_solutions() [
 * @param band band index [0, get_num_equations_per_node]
 * @return new vector with amplitude values at nodes
 */
 /*
template<class T, class TMat>
double* EigenProblem3D<T,TMat>::get_band_probability_vector(unsigned int idx, unsigned int band) const {
	T* values			  = this->get_band_vector(idx, band);
	unsigned int num_node = this->geometry->get_num_nodes();
	double* ret			  = new double[num_node];
	double  tmp;
	for(unsigned int ii = 0; ii < num_node; ii++) {
		tmp = abs(values[ii]);
		ret[ii] = tmp * tmp;
	}
	delete[] values;
	return ret;
}*/

/** calculates the probability vector
 *
 * \f$ \node\Phi(\mathbf{r})\node^{2} = \sum_{i = 1}^{\textrm{num bands}} \node f_{i}(\mathbf{r})\node^{2} \f$
 *
 * @param idx
 *
 */
template<class T, class TMat>
double* EigenProblem3D<T, TMat>::get_probability_vector(unsigned int idx) const {
	if(this->solution_vector.size() <= idx) {
		TDKP_GENERAL_EXCEPTION("solution index out of bounds");
	}
	const vector<T>& values   = this->solution_vector[idx];
	double tmp  		      = 0.0;
	int numeq   			  = this->get_num_equations_per_node();
	int numnode 			  = this->geometry.get_num_nodes();
	double* res 			  = new double[numnode];
	int int_id  			  = 0;
	int glb_id  			  = 0;

	Geometry::node_const_iterator it;

	// loop over nodes and store probability
	for(it = this->geometry.nodes_begin(); it != this->geometry.nodes_end(); it++) {
		if((*it)->get_index_internal() == -1) {
			res[(*it)->get_index_global()] = 0.0;
		} else {
			int_id = (*it)->get_index_internal();
			glb_id = (*it)->get_index_global();
			res[glb_id] = 0.0;
			for(int ii = 0; ii < numeq; ii++) {
				tmp = abs(values[int_id * numeq + ii]);
				res[glb_id] += (tmp * tmp);
			}
		}
	}
	return res;
}

template<class T, class TMat>
const T& EigenProblem3D<T,TMat>::get_eigenvalue(unsigned int idx) const {
	if(this->solution_value.size() > idx) {
		return this->solution_value[idx];
	} else {
		TDKP_GENERAL_EXCEPTION("solution index out of bounds");
	}
}

template<class T, class TMat>
T& EigenProblem3D<T,TMat>::get_eigenvalue(unsigned int idx) {
	return const_cast<T&>(
		static_cast<const EigenProblem3D<T,TMat>&>(*this).get_eigenvalue(idx)
	);
}

template<class T, class TMat>
void EigenProblem3D<T,TMat>::display_solution_info() const {
	ostringstream sout;	
	for(unsigned int ii = 0; ii < this->solution_value.size(); ii++) {
		sout.str("");
		sout << "solution [" << ii << "]: " << this->solution_value[ii] << " [eV]";
		Logger::get_instance()->emit(LOG_INFO, sout.str());
	}
}

template<class T, class TMat>
EigenSolution<double>* EigenProblem3D<T,TMat>::get_probability_object(unsigned int idx) const {
	double* prob = this->get_probability_vector(idx);
	EigenSolution<double> * a = new EigenSolution<double>(tdkp_math::only_real(this->get_eigenvalue(idx)), prob, this->geometry.get_num_nodes());
	a->set_identifier(0, string("probability"));
	delete[] prob;
	return a;
}

/** return eigensolution object
 *
 * creates eigensolution object with values on nodes
 */
template<class T, class TMat>
EigenSolution<T>*  EigenProblem3D<T,TMat>::get_solution_object(unsigned int idx) const {

	// ----------------------------------------------------
	// create empty object object, set energy and identifier
	// ----------------------------------------------------
	EigenSolution<T>* tmp = new EigenSolution<T>(this->geometry.get_num_nodes(), this->get_num_equations_per_node());
	tmp->set_energy(this->get_eigenvalue(idx));
	const vector<T>& sol = this->get_raw_eigenvector(idx);
	const int num_eq     = this->get_num_equations_per_node();

	TDKP_ASSERT(sol.size() == this->get_num_equations_per_node() * this->geometry.get_num_nonzero_nodes(), "sol.size() == this->get_num_equations_per_node() * this->geometry->get_num_nonzero_nodes()");

	// ----------------------------------------------------
	// remap from local to global indices
	// ----------------------------------------------------
	// redhat g++4 crashes here ...
#ifndef DEADRAT	
	#pragma omp parallel for default(shared)
#endif	
	for(int ii = 0; ii < (signed)this->geometry.get_num_nodes(); ii++) {
		const Node& node = this->geometry.get_node(ii);
		if(node.get_index_internal() != -1) {
			for(int nn = 0; nn < num_eq; nn++) {
				TDKP_BOUNDS_ASSERT(node.get_index_internal() * num_eq + nn >= 0 && node.get_index_internal() * num_eq + nn < (signed)sol.size(), "node.get_index_internal() * num_eq + nn >= 0 && node.get_index_internal() * num_eq + nn < sol.size()");
				tmp->get_node_value(node.get_index_global(), nn) = sol[node.get_index_internal() * num_eq + nn];
			}
		} else {
			for(int nn = 0; nn < num_eq; nn++) {
				tmp->get_node_value(node.get_index_global(), nn) = 0;
			}
		}
    }

	return tmp;
}

template<class T, class TMat>
EigenSolutionSet<double>*  EigenProblem3D<T,TMat>::get_probability_objects() const {
	EigenSolutionSet<double>* set = new EigenSolutionSet<double>();
	for(int ii = 0; ii < this->num_solutions(); ii++) {
		set->add(this->get_probability_object(ii));
	}
	return set;
}

template<class T, class TMat>
EigenSolutionSet<T>* EigenProblem3D<T,TMat>::get_solution_objects() const {
	EigenSolutionSet<T>* set = new EigenSolutionSet<T>();
	for(int ii = 0; ii < this->num_solutions(); ii++) {
		set->add(this->get_solution_object(ii));
	}
	return set;
}

template<class T, class TMat>
void EigenProblem3D<T,TMat>::update_bandstructure_container() {

	// -------------------------------------------------------
	// reinit bandstructure object anyway
	// -------------------------------------------------------
	unsigned int sol_length = this->geometry.get_num_nodes();

	// -------------------------------------------------------
	// create domain with a singular point of weight 1.0
	// -------------------------------------------------------
	DomainMaster domain(new DomainNodeSingularPoint(1.0));
	this->solution_container.reinit(
		this->get_num_equations_per_node(),
		this->solution_value.size(),
		sol_length,
		domain
	);

	// -------------------------------------------------------
	// pass solutions (copy paste from get_solution_object)
	// TODO: get rid of functions get_solution_object and sets ...
	// -------------------------------------------------------
	for(int oo = 0; oo < (signed)solution_vector.size(); oo++) {
		// ----------------------------------------------------
		// create empty object object, set energy and identifier
		// ----------------------------------------------------
		EigenSolution<cplx>* tmp = new EigenSolution<cplx>(this->geometry.get_num_nodes(), this->get_num_equations_per_node());
		tmp->set_energy(this->get_eigenvalue(oo));
		const vector<T>& sol = this->get_raw_eigenvector(oo);
		const int num_eq     = this->get_num_equations_per_node();

		TDKP_ASSERT(sol.size() == this->get_num_equations_per_node() * this->geometry.get_num_nonzero_nodes(), "sol.size() == this->get_num_equations_per_node() * this->geometry->get_num_nonzero_nodes()");

		// ----------------------------------------------------
		// remap from local to global indices
		// ----------------------------------------------------
#ifndef DEADRAT				
		#pragma omp parallel for default(shared) schedule(static, 5000)
#endif		
		for(int ii = 0; ii < (signed)this->geometry.get_num_nodes(); ii++) {
			const Node& node = this->geometry.get_node(ii);
			if(node.get_index_internal() != -1) {
				for(int nn = 0; nn < num_eq; nn++) {
					TDKP_BOUNDS_ASSERT(node.get_index_internal() * num_eq + nn >= 0 && node.get_index_internal() * num_eq + nn < (signed)sol.size(), "node.get_index_internal() * num_eq + nn >= 0 && node.get_index_internal() * num_eq + nn < sol.size()");
					tmp->get_node_value(node.get_index_global(), nn) = sol[node.get_index_internal() * num_eq + nn];
				}
			} else {
				for(int nn = 0; nn < num_eq; nn++) {
					tmp->get_node_value(node.get_index_global(), nn) = 0;
				}
			}
    	}
    	this->solution_container.add_eigensolution(0, oo, tmp);
	}
}

} // end of namespace

#endif /*EIGENPROBLEM3D_H_*/
