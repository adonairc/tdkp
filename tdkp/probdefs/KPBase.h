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

#ifndef KPBASE_H_
#define KPBASE_H_

#include "tdkp/probdefs/EigenProblem3D.h"
#include "tdkp/kpmatrices/KPMatrixBase.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/common/Configuration.h"
#include <omp.h>

namespace tdkp {

/** common class for 1D - 3D kp problems */
template<class T>
class KPBase : public T {
public: 
	KPBase(const Geometry& geometry_, MaterialDatabase& material_database_);
	virtual ~KPBase();
	StdElementData<double>* get_bandedges() throw(Exception*);	
	virtual void set_solution_type(KPSolutionType type);
	const RMatrix<double>& get_rotation_matrix() const { return rotation_matrix; }
	/** return band edge limits */	
	void get_minmax_edges(double& cb_min, double& cb_max, double& vb_min, double& vb_max);		
	/** get virtual energy shift of the system */
	virtual double get_energy_shift() const = 0;
		
protected:

	/** virtual function that is derived by all subclasses to return the appropriate kp matrix */
	virtual KPMatrixBase* get_matrix() const = 0;
	/** vector with all kp matrices */
	vector<KPMatrixBase*> kp_matrices;
	/** copy of kp matrix sparsity pattern */
	int*     sparsity_copy;
	/** sparsity pattern length */					
	int      sparsity_num;	
	/** k directions, see comment on set_axes functions for 1D2D and 3D for details */
	Vector3D k_directions[3];
	/** system rotation matrix (calculated in set axes) */	
	RMatrix<double> rotation_matrix; 
			
};

// ---------------------------------------------
// KPBase template implementation
// ---------------------------------------------


/** main and the only constructor */
template<class T>
KPBase<T>::KPBase(const Geometry& geometry_, MaterialDatabase& material_database_)
: T(geometry_, material_database_),
  sparsity_copy(0),
  sparsity_num(0),
  rotation_matrix(3,3) 
{
	// -------------------------------------------------
	// default is no rotation
	// -------------------------------------------------
	k_directions[0] = Vector3D(1.0, 0.0, 0.0);
	k_directions[1] = Vector3D(0.0, 1.0, 0.0);
	k_directions[2] = Vector3D(0.0, 0.0, 1.0);
	rotation_matrix(0,0) = rotation_matrix(1,1) = rotation_matrix(2,2) = 1.0;   
}

/** constructor deleteing kp matrices and sparsity patterns */
template<class T>
KPBase<T>::~KPBase() {
	if(this->sparsity_copy) {
		delete[] this->sparsity_copy; this->sparsity_copy = 0;
	}
	for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
		delete this->kp_matrices[ii]; this->kp_matrices[ii] = 0;
	}
	this->kp_matrices.resize(0);	
}


/** set kp solution type (but check against kp matrix if this type is valid) 
 */  
template<class T>
void KPBase<T>::set_solution_type(KPSolutionType type) {
    KPMatrixBase* tmp = this->get_matrix();
    if(tmp->check_solution_type_available(type)) {
        T::set_solution_type(type);
        delete tmp; 
    } else {
        ostringstream sout;
        sout << "you can not solve for " << type << " in the selected kp model";
        delete tmp;
        TDKP_GENERAL_EXCEPTION(sout.str()); 
    }
}

/** get bandedges for each element
 * 
 * particularly interesting when strain is present
 */ 
template<class T>
StdElementData<double>* KPBase<T>::get_bandedges() throw(Exception*) {

    if(this->kp_matrices.size() == 0) {
        TDKP_GENERAL_EXCEPTION("object is not prepared! call .prepare before calculating bandedge.");
    }
    // -------------------------------------------------------------------
    // create std element data and name bandedges
    // -------------------------------------------------------------------
    const int num_equations = this->get_num_equations_per_node(); 
    StdElementData<double>* ret = new StdElementData<double>(this->get_num_equations_per_node(), this->geometry.get_num_elements());
    for(int bb = 0; bb < num_equations; bb++) {
        ostringstream sout;
        sout << "BandEdge" << bb;
        ret->set_identifier(bb, sout.str().c_str());        
    }
    // surpress output (we're going to copy the matrices
    // for a parallel evaluation)
    bool surpressed = kp_matrices.front()->surpress_output();
    for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
    	this->kp_matrices[ii]->set_output_surpression(true);	
    }
       
    // -------------------------------------------------------------------
    // calculate element matrices and extract bandedges
    // -------------------------------------------------------------------	   
    #pragma omp parallel 
    {
    	// create copy of kp matrices for parallel evaluation
    	vector<KPMatrixBase*> copy_kp_matrices(this->kp_matrices.size());
    	for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
    		copy_kp_matrices[ii] = this->get_matrix();
    		copy_kp_matrices[ii]->set_output_surpression(true);
    		copy_kp_matrices[ii]->init_from_kpmatrix(*kp_matrices[ii]);	   
    	}    	
    	
	    KPMatrixBase*  mat       = 0;    
	    RMatrix<cplx>* bg_matrix = 0;
	    vector<cplx>   eigenvalues;
	    vector<cplx>   eigenvectors;      
	    
	    int next, current; next = current = 0;
    	int start_idx, end_idx;
    	if(this->geometry.get_num_elements() > 100000) {
    		get_omp_range(start_idx, end_idx, this->geometry.get_num_elements());
    	} else {
    		start_idx = 0; end_idx = this->geometry.get_num_elements();
    	}
	    const bool progress_bar = Configuration::get_instance()->get("output_bandedges_progress_bar_min_elements") < this->geometry.get_num_elements() && omp_get_thread_num() == 0;    
		if(progress_bar) {    
    		Logger::get_instance()->init_progress_bar("KPBase: calculating the bandedge in elements", end_idx - start_idx);
		}    	
    	
	    for(int ii = start_idx; ii < end_idx; ii++) {
	    	const Element& elem = this->geometry.get_element(ii);
	    	if(elem.enabled() && !elem.get_region().enabled()) {
	    		TDKP_TRACE("stinky element " << elem.get_index_global());	
	    	}
	    		    	
	    	// only for enabled elements
	    	if(elem.enabled()) {    		    		    		    	
		        // get kp matrix and set strains and potential energy
		        mat = copy_kp_matrices[elem.get_region().get_material().get_id()];

		        if(this->strain_field_set()) {
		            mat->set_strains(this->get_strain_field().get(elem.get_index_global()));
		        }
		        if(this->potential_energy_field_set()) {
		        	// find average potential in element
		        	double potential = 0.0;
		        	for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
		        		potential += this->get_potential_energy_field().get_node_value(elem.get_node(nn).get_index_global());
		        	}		        	
		            mat->set_potential(potential / static_cast<double>(elem.get_num_nodes()));	                  
		        }
		        mat->set_energy_shift(0.0);
		        // recalculate and get zero order matrix
		        mat->calculate();
		        bg_matrix = mat->get_zero_order_matrix();   
		                
		        // get eigenvalues of that zero order matrix
		        #pragma omp critical (kpbase_eigensystem)
		        {
		        	RMatrix<cplx>::get_eigensystem(*bg_matrix, eigenvalues, eigenvectors);
		        }      
		        for(int bb = 0; bb < num_equations; bb++) {
		            ret->set_element_value(elem.get_index_global(), bb, eigenvalues[bb].real());
		        }
		        delete bg_matrix; bg_matrix = 0;
	    	}
	        if(progress_bar && next == current) {
	            next = Logger::get_instance()->set_progress_bar(current, end_idx - start_idx);
	        }
	        current++;    	
	    }
	    if(progress_bar) {
    		Logger::get_instance()->end_progress_bar();
    	}   
    	// clean up
    	for(unsigned int ii = 0; ii < copy_kp_matrices.size(); ii++) {
    		delete copy_kp_matrices[ii]; copy_kp_matrices[ii] = 0;	
    	}    	
    } // end omp parallel

	// reset output 
    for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
    	this->kp_matrices[ii]->set_output_surpression(surpressed);	
    }
    TDKP_LOGMSG(LOG_INFO_DEVEL2, "KPBase: bandedges succsessfully calculated");

    return ret;
}

/** return band edge limits 
 * 
 * at return, the arguments contain the band edge extrema. if the values
 * could not be determined (e.g. cb in kp4x4), -1 is assigned to the value
 * 
 * @param cb_min	lowest conduction band edge. all cb states are above it. -1 if not available
 * @param cb_max    upper most conduction band edge 
 * @param vb_min    upper vb edge (so all vb states are below it).
 */
template<class T>
void KPBase<T>::get_minmax_edges(double& cb_min, double& cb_max, double& vb_min, double& vb_max) {

	unsigned int cb_idx = 6;
	unsigned int vb_idx = 5;
	
	TDKP_ASSERT(this->get_num_equations_per_node() == 8, "");	
	StdElementData<double>* bandedges = this->get_bandedges();
	
	cb_min = cb_max = bandedges->get_element_value(0, cb_idx);
	vb_min = vb_max = bandedges->get_element_value(0, vb_idx);
	
	double cb_tmp, vb_tmp;
	for(int ii = 1; ii < bandedges->get_length(); ii++) {
		 cb_tmp = bandedges->get_element_value(ii, cb_idx);
		 vb_tmp = bandedges->get_element_value(ii, vb_idx);
		 if(cb_tmp < cb_min) {
		 	cb_min = cb_tmp;	
		 }				
		 if(cb_tmp > cb_max) {
		 	cb_max = cb_tmp;	
		 }
		 if(vb_tmp > vb_min) {
		 	vb_min = vb_tmp;	
		 }
		 if(vb_tmp < vb_max) {
		 	vb_max = vb_tmp;	
		 }
	}
		
	delete bandedges;

}

} // end of namespace 

#endif /*KPBASE_H_*/
