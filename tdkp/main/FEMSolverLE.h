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

#ifndef FEMSOLVERLE_H_
#define FEMSOLVERLE_H_

#include "tdkp/main/FEMSolver.h"
#include "tdkp/main/SparseMatrixInterface.h"
#include "tdkp/probdefs/LinearProblem.h"
#include "tdkp/solvers/LinearSolver.h"
#include "tdkp/geometry/BoundaryCondition.h"
#include <omp.h>

namespace tdkp {

/** solver class for linear pde equations */
template<class T>
class FEMSolverLE : public tdkp::FEMSolver {

public:	
	typedef SparseMatrixInterface<T>   SMatrix;

	FEMSolverLE(const Geometry& geometry_, LinearProblem<T>& problem_);
	FEMSolverLE(const Geometry& geometry_, LinearProblem<T>& problem_, LinearSolver<T>* solver_);
	virtual ~FEMSolverLE();
	// ---------------------------------------------------
	// simulation setup functions 
	// ---------------------------------------------------
	bool verify_setup() const;	
	void set_be_quiet(bool shut_up) { be_quiet = shut_up; }		
	
	// ---------------------------------------------------
	// assembling 
	// ---------------------------------------------------
	void create_matrix_structures();
	void reset_matrix_to_zero();
	void assemble_system();
	void assemble_rhs_only();
	void solve_system(int num_solutions);
	
	// ---------------------------------------------------
	// access to matrices M in -> out
	// ---------------------------------------------------
	template<class V> 
	void multiply_with_lhs(const vector<V>& in, vector<V>& out) const;
	// ---------------------------------------------------
	// calculate the sum of rhs
	// ---------------------------------------------------
	T    sum_rhs() const; 
	const T* get_rhs() const { return load; } 	
			
protected:

    LinearProblem<T>& problem;    
    LinearSolver<T>*  solver;	
	SMatrix*  stiff;
	T*		  load;
	bool      solver_prepared;
	
	void enforce_partial_boundary_conditions_on_matrix_and_rhs();
	void enforce_partial_boundary_conditions_on_rhs();
	
private:		
	void delete_matrices();
	bool be_quiet;							 			

};

template<class T>
FEMSolverLE<T>::FEMSolverLE(const Geometry& geometry_, LinearProblem<T>& problem_, LinearSolver<T>* solver_)
: FEMSolver(geometry_), 
  problem(problem_),
  solver(solver_),
  stiff(&solver->get_matrix()),
  load(0),
  solver_prepared(false),
  be_quiet(false)
{
	
}
  
template<class T>
FEMSolverLE<T>::FEMSolverLE(const Geometry& geometry_, LinearProblem<T>& problem_) 
: FEMSolver(geometry_),
  problem(problem_),
  solver(0),
  stiff(0),
  load(0),
  solver_prepared(false),
  be_quiet(false) 
{	
}

template<class T>
FEMSolverLE<T>::~FEMSolverLE() {
	this->delete_matrices();
	if(this->solver != 0) {
		delete solver;
		solver = 0;
	}
}

template<class T>
bool FEMSolverLE<T>::verify_setup() const {
	
	if(!this->geometry.verify()) {
		Logger::get_instance()->emit(LOG_ERROR, "geometry is invalid");
		return false;	
	}	

	return true;
}


template<class T>
void FEMSolverLE<T>::create_matrix_structures() {
		
	TimeMeasurements::get_instance().start("matrix graph building");		
		
	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");	
	}
			
    int  	   sparse_num;
    const int* sparse_pat    = this->problem.get_node_sparsity_pattern(sparse_num);    
	int        num_equations = this->problem.get_num_equations_per_node();

	// ------------------------------------------------------------
	// create linear solver and get matrices
	// -----------------------------------------------------------
	if(this->solver == 0) {
		this->solver = LinearSolver<T>::factory(this->geometry.get_num_nonzero_nodes() * num_equations);
		TDKP_POINTER_ASSERT(this->solver);
		this->stiff = &this->solver->get_matrix();
	} else {
		TDKP_ASSERT(this->stiff != 0, "this->stiff != 0");
		TDKP_ASSERT(this->stiff->get_size() == this->geometry.get_num_nonzero_nodes() * num_equations, "this->stiff->get_size() == this->geometry.get_num_nonzero_nodes() * num_equations"); 
		this->stiff->reset();
	}
	solver_prepared = false;
	 	
	int length = num_equations * this->geometry.get_num_nonzero_nodes();
	this->load = new T[length];   
	TDKP_POINTER_ASSERT(this->load);
	for(int ii = 0; ii < length; ii++) {
		this->load[ii] = 0.0;
	}

	if(!stiff || !load) {
		TDKP_GENERAL_EXCEPTION("could not allocate memory for stiffness/mass matrix");	
	}
	
	// -------------------------------------------------------
	// load structure (if available)
	// -------------------------------------------------------
	if(Configuration::get_instance()->get("assembly_use_matrix_structure_cache") == 1.0) {
		if(this->stiff->caching_structure_possible()) {
			ostringstream sout_stiff, sout;
			sout_stiff   << "stiff_" << this->problem.get_unique_identifier() << "_" << this->geometry.get_identifier();
			if(this->stiff->structure_cache_available(sout_stiff.str())) {
				if(omp_get_thread_num() == 0) {
				   	sout << "FEMSolverLE: found structure cache in file. loading them";
				   	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
					TimeMeasurements::get_instance().stop("matrix graph building");				   	
				}
				this->stiff->load_structure(sout_stiff.str());
				return;	
			}	
		}	
	}	
	
	#pragma omp parallel default(shared)
	{
	int        gii,gjj;
	bool       symmetric;

	int processor_first;
	int processor_last;
	
	symmetric = stiff->property_is_set(symmetric_matrix);			
	get_omp_matrix_range(processor_first, processor_last, this->geometry.get_num_nonzero_nodes(), symmetric);

	if(omp_get_num_threads() > 1 && omp_get_thread_num() == 0) {
		ostringstream sout;
		sout << "FEMSolverLE: linewise parallel building of "
		     << (symmetric ? "symmetric" : "nonsymmetric") 
		     << " matrix structure with "
		     << omp_get_num_threads() << " active threads "
		     << "where the average processor takes care of " 
		     << (processor_last - processor_first) * num_equations
		     << " matrix lines.";
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());	
	} 
	if(Logger::get_instance()->get_level() > LOG_INFO && omp_get_thread_num() == 0) {
		ostringstream sout;
		sout << "FEMSolverLE: building " << (symmetric ? "symmetric" : "nonsymmetric") << " matrix structure";
		Logger::get_instance()->init_progress_bar(sout.str(), this->geometry.get_num_elements());
	}
		

		
	// set iterators
	Geometry::element_const_iterator eit;
	Geometry::element_const_iterator end = this->geometry.elements_end();

	int next = 0;
	int cc   = 0;
	
	// loop over elements to build matrix structure
	for(eit = this->geometry.elements_begin(); eit != end; eit++) {
		// only process enabled elements (non-contacts)
		if((*eit)->enabled()) {
			for(unsigned int ii = 0; ii < (*eit)->get_num_nodes(); ii++) {			
				if((*eit)->get_node(ii).get_index_internal() != -1) {
					gii = (*eit)->get_node(ii).get_index_internal();
					// only take the lines of our processor!
					if(processor_first <= gii && gii < processor_last) {
						for(unsigned int jj = 0; jj < (*eit)->get_num_nodes(); jj++) {
							if((*eit)->get_node(jj).get_index_internal() != -1) {
								gjj = (*eit)->get_node(jj).get_index_internal();
								// -----------------------------------						
								// SYMMETRIC 
								// -----------------------------------						
								if(symmetric) {		
									if(gii < gjj) {							
										// loop over sparsity pattern for stiffness matrix
										int last_ss;
										try{
											for(int ss = 0; ss < sparse_num; ss++) {
												last_ss = ss;
												stiff->announce(gii * num_equations + sparse_pat[2*ss], gjj * num_equations + sparse_pat[2*ss + 1]);									
											}
										} catch(Exception* e) {								
											std::cout << "stiff this was element: " << (eit - this->geometry.elements_begin()) / sizeof(Geometry::element_iterator) 
											          << "and the guy we announced was " << gii * num_equations + sparse_pat[2*last_ss] << ", " << gjj * num_equations + sparse_pat[2*last_ss + 1] << "\n";
											throw e;								          
										}						
									} else if (gii == gjj) {
										// loop over sparsity pattern for stiffness matrix
										int lii, ljj, last_ss;
										try{
											for(int ss = 0; ss < sparse_num; ss++) {
												last_ss = ss;
												lii = gii * num_equations + sparse_pat[2*ss];
												ljj = gjj * num_equations + sparse_pat[2*ss + 1];
												if(lii <= ljj) {
													stiff->announce(lii,ljj);									
												}
											}
										} catch(Exception* e) {								
											std::cout << "stiff this was element: " << (eit - this->geometry.elements_begin()) / sizeof(Geometry::element_iterator) 
											          << "and the guy we announced was " << gii * num_equations + sparse_pat[2*last_ss] << ", " << gjj * num_equations + sparse_pat[2*last_ss + 1] << "\n";
											throw e;								          
										}							
									} 
								} else { // if(symmetric)
								// -----------------------------------						
								// END SYMMETRIC / BEGIN NONSYMMETRIC
								// -----------------------------------													
									int last_ss;
									try{
										for(int ss = 0; ss < sparse_num; ss++) {
											last_ss = ss;
											stiff->announce(gii * num_equations + sparse_pat[2*ss], gjj * num_equations + sparse_pat[2*ss + 1]);									
										}
									} catch(Exception* e) {								
										std::cout << "stiff this was element: " << (eit - this->geometry.elements_begin()) / sizeof(Geometry::element_iterator) 
									          << "and the guy we announced was " << gii * num_equations + sparse_pat[2*last_ss] << ", " << gjj * num_equations + sparse_pat[2*last_ss + 1] << "\n";
										throw e;								          
									}				
								}
							}
						}											
					}
				}
			}
		}
		if(next == cc && Logger::get_instance()->get_level() > LOG_INFO && omp_get_thread_num() == 0) {
			next = Logger::get_instance()->set_progress_bar(cc, this->geometry.get_num_elements());	
		}
		cc++;
	}
	if(Logger::get_instance()->get_level() > LOG_INFO && omp_get_thread_num() == 0) {
		Logger::get_instance()->end_progress_bar();
	}
		
	} // end of parallel loop		
		
	// compile structure
	stiff->set_structure();			
	
	TimeMeasurements::get_instance().stop("matrix graph building");
	TimeMeasurements::get_instance().track_memory_usage();
	
	// -------------------------------------------------------
	// save structure (if available and requested)
	// -------------------------------------------------------
	if(omp_get_thread_num() == 0 && Configuration::get_instance()->get("assembly_use_matrix_structure_cache") == 1.0) {
		if(this->stiff->caching_structure_possible()) {
			ostringstream sout_stiff, sout;
			sout_stiff   << "stiff_" << this->problem.get_unique_identifier() << "_" << this->geometry.get_identifier();
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverLE: writing matrix structure to cache");
			this->stiff->save_structure(sout_stiff.str());
		} else {
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverLE: writing matrix structure to cache not supported for current matrix type");
		}
	}				
		
}

template<class T>
void FEMSolverLE<T>::reset_matrix_to_zero() {
	stiff->clear_but_keep_structure();
	const int neq = this->problem.get_num_equations_per_node(); 
	const int num = this->geometry.get_num_nonzero_nodes(); 
	for(int ii = 0; ii < neq * num; ii++) {
		this->load[ii] = (T)0;	
	}
}

template<class T>
void FEMSolverLE<T>::assemble_rhs_only() {

	TimeMeasurements::get_instance().start("rhs assembly");

	// ------------------------------------------------
	// object must be ready
	// ------------------------------------------------
	TDKP_ASSERT(this->problem.is_ready(), "assembling of rhs is only possible if assemble_system was called before");
	if(this->stiff == 0) {
		TDKP_GENERAL_EXCEPTION("matrice structures must be build before assembly.");	
	}
	// ------------------------------------------------
	// tell the problem object that we will only assemble rhs
	// ------------------------------------------------
	this->problem.set_build_rhs_only(true);
	
    
	
	// ------------------------------------------------
	// iterators
	// ------------------------------------------------
	Geometry::element_const_iterator eit;	
	Geometry::element_const_iterator end = this->geometry.elements_end();
	// ------------------------------------------------
	// some space to work on (attention: this is only enough for elements with max Element::max_num_nodes nodes!)
	// ------------------------------------------------
	int num_equations = this->problem.get_num_equations_per_node();
	int   nlhs    = num_equations * num_equations * Element::max_num_nodes * Element::max_num_nodes;
	T*    lhs     = new T[nlhs];
	double* rhs   = new double[nlhs];
	int*  idxs    = new int[Element::max_num_nodes];
	TDKP_POINTER_ASSERT(lhs && rhs && idxs);
	int   n       = 0;
	int   gii;
	int   index_internals[Element::max_num_nodes];
	
	// ------------------------------------------------
	// reset rhs
	// ------------------------------------------------
	const int lsize = this->stiff->get_size();
	for(int ii = 0; ii < lsize; ii++) {
		this->load[ii] = T(0);	
	}	
	
	// -----------------------------------------------------------
	// loop over elements and assemble
	// -----------------------------------------------------------
	for(eit = this->geometry.elements_begin(); eit != end; eit++) {
		
		// only process enabled elements (non-contacts)
		if((*eit)->enabled()) {		
		
			// calculate element matrix				
			this->problem.calculate_element_matrices((*eit), lhs, rhs, idxs, n);		
			
			// get global indices
			for(int ii = 0; ii < n; ii++) {
				index_internals[ii] = (*eit)->get_node(idxs[ii]).get_index_internal();	
			}		
			// for all nodes in element				
			for(int ii = 0; ii < n; ii++) {
				gii = index_internals[ii];				
				// ----------------------------------------
				// assemble rhs
				// ----------------------------------------
				for(int ee = 0; ee < num_equations; ee++) {
					this->load[gii * num_equations + ee] += rhs[ii * num_equations + ee];
				}							
			}
		}
	}

	this->enforce_partial_boundary_conditions_on_rhs();

	delete[] lhs;
	delete[] rhs;
	delete[] idxs;
	lhs = 0; rhs = 0; idxs = 0;
	
	// and back to standard mode
	this->problem.set_build_rhs_only(false);
	
	
	TimeMeasurements::get_instance().stop("rhs assembly");
	TimeMeasurements::get_instance().track_memory_usage();		
		
	
}

template<class T>
void FEMSolverLE<T>::assemble_system() {

	TimeMeasurements::get_instance().start("matrix assembly");

	if(this->stiff == 0) {
		TDKP_GENERAL_EXCEPTION("matrice structures must be build before assembly.");	
	}
	this->problem.prepare();
	solver_prepared = false;
	
	// calculate element matrices
	// get sparsity pattern ( local_nz)
/*
for ii in node_internal_indices
	for jj in node_internal_indices
		global_ii = node_internal_indices[ii]->global_index
		global_jj = node_internal_indices[jj]->global_index
		for nn (every second in size(sparsity_pattern))
			matrix(global_ii + sparsepat(2nn), global_jj + sparsepat(2nn + 1)) += lhs((ii * n + jj) * local_mat_nonsparse + nn)
	*/	

	if(!be_quiet && Logger::get_instance()->get_level() > LOG_INFO) {	
		Logger::get_instance()->init_progress_bar("FEMSolverLE: assembling matrix and rhs", this->geometry.get_num_elements());
	}	

	// ----------------------------------------------
	// shared data
	// ----------------------------------------------
    int  	   sparse_num;    
	const int  num_equations = this->problem.get_num_equations_per_node();    	
	const int* sparse_pat = this->problem.get_node_sparsity_pattern(sparse_num);
	const bool symmetric  = this->stiff->property_is_set(symmetric_matrix);
	// some space to work on (attention: this is only enough for elements with max Element::max_num_nodes nodes!)	
	const int  nlhs = num_equations * num_equations * Element::max_num_nodes * Element::max_num_nodes;	
	
	// ----------------------------------------------
	// parallel section
	// ----------------------------------------------
#pragma omp parallel default(shared)
	{

		// --------------------------------------------------
		// private data
		// --------------------------------------------------
		int        count      = 0;
		int        next       = 0;				
		// iterators
		Geometry::element_const_iterator eit;	
		Geometry::element_const_iterator end = this->geometry.elements_end();
		T*    lhs     = new T[nlhs];
		T*    rhs     = new T[nlhs];
		int   idxs[Element::max_num_nodes];
		TDKP_POINTER_ASSERT(lhs && rhs);
		int   n       = 0;
		int   gii, gjj;	
		int   index_internals[Element::max_num_nodes];
			
		// --------------------------------------------------
		// processor distribution
		// --------------------------------------------------
		int processor_first;
		int processor_last;
		bool my_rows = false;
		int  index_internal;		
		get_omp_matrix_range(processor_first, processor_last, this->geometry.get_num_nonzero_nodes(), symmetric);
		TDKP_ASSERT(processor_first < processor_last, "processor_first (" << processor_first << ") < processor_last (" << processor_last << ")");								
						
		// -----------------------------------------------------------
		// loop over elements and assemble
		// -----------------------------------------------------------
		for(eit = this->geometry.elements_begin(); eit != end; eit++) {
			
			// only process enabled elements (non-contacts)
			if((*eit)->enabled()) {		
			
				// ---------------------------------------------
				// check if element interests the current
				// processor
				// ---------------------------------------------
				my_rows = false;
				for(unsigned int ii = 0; ii < (*eit)->get_num_nodes(); ii++) {
					index_internal = (*eit)->get_node(ii).get_index_internal();
					if(processor_first <= index_internal && index_internal < processor_last) {
						my_rows = true;
						break;	
					}		
				}
				
				// -------------------------------------------
				// if at least one row is covered by the element,
				// we evaluate. if the grid is ordered properly, nodes
				// among an element should usually be on the same processor
				// -------------------------------------------	
				if(my_rows) {								
					// calculate element matrix				
					this->problem.calculate_element_matrices((*eit), lhs, rhs, idxs, n);		
					
					// get global indices
					for(int ii = 0; ii < n; ii++) {
						index_internals[ii] = (*eit)->get_node(idxs[ii]).get_index_internal();	
					}		
					// for all nodes in element				
					for(int ii = 0; ii < n; ii++) {
						gii = index_internals[ii];
						if(processor_first <= gii && gii < processor_last) {
							// ----------------------------------------
							// assemble rhs
							// ----------------------------------------
							for(int ee = 0; ee < num_equations; ee++) {
								this->load[gii * num_equations + ee] += rhs[ii * num_equations + ee];
							}									
							// for all nodes in element (again ;-))
							for(int jj = 0; jj < n; jj++) {	
								gjj = index_internals[jj];
								if(symmetric) {
									if(gii < gjj) {
										// -----------------------------------------------------------------
										// loop over all local matrix elements (for lhs)
										// -----------------------------------------------------------------				
										for(int ss = 0; ss < sparse_num; ss++) {
											stiff->add  (gii * num_equations + sparse_pat[2*ss], gjj * num_equations + sparse_pat[2*ss+1], lhs[(n * ii + jj) * sparse_num + ss]);				
										}
									} else if(gii == gjj) {
										// -----------------------------------------------------------------
										// loop over all local matrix elements (for lhs)
										// -----------------------------------------------------------------				
										int lii, ljj;
										for(int ss = 0; ss < sparse_num; ss++) {
											lii = gii * num_equations + sparse_pat[2*ss];
											ljj = gjj * num_equations + sparse_pat[2*ss+1];
											if(lii <= ljj) {
												stiff->add  (lii,ljj, lhs[(n * ii + jj) * sparse_num + ss]);				
											}
										}
									}
								} else {
									// -----------------------------------------------------------------
									// NONSYMMETRIC matrices
									// -----------------------------------------------------------------
									// loop over all local matrix elements (for lhs)
									// -----------------------------------------------------------------				
									for(int ss = 0; ss < sparse_num; ss++) {
										stiff->add  (gii * num_equations + sparse_pat[2*ss], gjj * num_equations + sparse_pat[2*ss+1], lhs[(n * ii + jj) * sparse_num + ss]);				
									}
								}
							}
						}
					}
				}
			}
			if(!be_quiet && omp_get_thread_num() == 0 && count == next && Logger::get_instance()->get_level() > LOG_INFO) {
				next = Logger::get_instance()->set_progress_bar(count, this->geometry.get_num_elements());
			}				
			count++;
		}
		// kill parallel temporary space
		delete[] lhs;
		delete[] rhs;
		lhs = 0; rhs = 0;		
	} // end parallel
	
	if(!be_quiet && Logger::get_instance()->get_level() > LOG_INFO) {
		Logger::get_instance()->end_progress_bar();
	}
	// enforce partial boundary conditions
	this->enforce_partial_boundary_conditions_on_matrix_and_rhs();
	
	// store matrices to files upon request
	if(Configuration::get_instance()->get("assembly_save_matrices_to_file") != 0.0) {
		stiff->save_to_file("lhs.mat");		
		fstream fout("rhs.mat", ios::out);	
		TDKP_ASSERT(fout, "could not open file rhs.mat for writing load vector");
		fout.precision(15);
		for(unsigned int ii = 0; ii < num_equations * this->geometry.get_num_nonzero_nodes(); ii++) {
			fout << this->load[ii] << "\n";	
		}
		fout.close();
	}

	// check matrices for symmtery upon request
	if(Configuration::get_instance()->get("assembly_check_matrix_for_symmetry") != 0.0) {
		stiff->perform_symmetry_analysis();	
	}
	
	TimeMeasurements::get_instance().stop("matrix assembly");
	TimeMeasurements::get_instance().track_memory_usage();	
				
}

template<class T>
void FEMSolverLE<T>::delete_matrices() {
	if(this->load != 0) {
		delete[] this->load;
		this->load = 0;
	}
	this->stiff = 0;
}
	

/** solve  linear equation
 */
template<class T> 
void FEMSolverLE<T>::solve_system(int num_solutions) {
	
	if(this->stiff == 0 || this->load == 0) {
		TDKP_GENERAL_EXCEPTION("can not solve problem without  having assembled it");	
	}

	T* result = new T[this->stiff->get_size()];
	TDKP_POINTER_ASSERT(result);
	
	TimeMeasurements::get_instance().start("linear equation solving");
	if(!solver_prepared) {	
		solver->prepare();
		solver_prepared = true;
	}
	TimeMeasurements::get_instance().track_memory_usage();
	// ------------------------------------
	// check if rhs is zero
	// ------------------------------------
	bool is_zero = true;
	const int lsize = this->stiff->get_size(); 
	for(int ii = 0; ii < lsize; ii++) {
		if(load[ii] != 0.0) {
			is_zero = false;
			break;
		}		
	}
	// ------------------------------------
	// solve for nonzero rhs
	// ------------------------------------
	if(is_zero) {
		for(int ii = 0; ii < lsize; ii++) {
			result[ii] = 0.0;	
		}
		TDKP_LOGMSG(LOG_WARN, "FEMSolverLE: not solving equation, rhs is zero, so solution is zero.");
	} else {		
		solver->solve_equation(result, this->load);
	}
	TimeMeasurements::get_instance().stop("linear equation solving");
	TimeMeasurements::get_instance().start("result postprocessing");
	this->problem.add_solution(result, this->stiff->get_size());
	TimeMeasurements::get_instance().stop("result postprocessing");	
	
	delete[] result;
	
}


/** multiply vector with lhs matrix
 * 
 * the lhs matrix should be hidden from the problem class side
 * in order to allow simple modification of the matrix and solving
 * routines. therefore we need to provide that interface for e.g. matrix
 * elements to be able to calculate some matrix vector products
 */
template<class T> 
template<class V> 
void FEMSolverLE<T>::multiply_with_lhs(const vector<V>& in, vector<V>& out) const {
	if(stiff) {
		TDKP_ASSERT(in.size() == this->stiff->get_size(), "vector and matrix do not have the same size!");
		out.resize(in.size()); 
		this->stiff->mult_vec(&in[0], &out[0]);
	} else {
		TDKP_GENERAL_EXCEPTION("lhs matrix is undefined and does not exist!");	
	}	
}

template<class T>
T FEMSolverLE<T>::sum_rhs() const {
	T ret(0);
	unsigned int length = this->geometry.get_num_nonzero_nodes();
	for(unsigned int ii = 0; ii < length; ii++) {
		ret += load[ii];			
	}	
	return ret;
}

/** sets row and column of dof to 0, diagonal to 1 and rhs is set to 0) */
template<class T>
void FEMSolverLE<T>::enforce_partial_boundary_conditions_on_matrix_and_rhs() {

	if(geometry.get_boundary_conditions().have_ignore_partially()) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "FEMSolverLE: enforcing partial dirichlet boundary conditions");
		const unsigned int num_equations = this->problem.get_num_equations_per_node();
		const BoundaryCondition& bc      = geometry.get_boundary_conditions();
		unsigned int erased = 0;		
		// for all nodes 
		for(unsigned int ii = 0; ii < geometry.get_num_nodes(); ii++) {
			const Node& node = geometry.get_node(ii);			
			// if node is non-pure dirichlet
			if(node.get_index_internal() != -1) {
				// check all dofs for dirichlet type
				for(unsigned int ee = 0; ee < num_equations; ee++) {  
					if(bc.ignore_node_partially(node,ee)) {
						// statistics
						erased++;
						// set row and column to zero
						stiff->set_row_and_column_to_zero(num_equations * node.get_index_internal() + ee);
						// set diagonal to one
						stiff->set(num_equations * node.get_index_internal() + ee, num_equations * node.get_index_internal() + ee, 1.0);
						// set rhs line to zero
						load[num_equations * node.get_index_internal() + ee] = 0.0;		
					}
				}
			}	
		}
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "FEMSolverLE: enforced " << erased << " dirichlet b.c. on dof level."); 
	}
	
}

/** sets row entries of partial b.c. dof to 0 */
template<class T>
void FEMSolverLE<T>::enforce_partial_boundary_conditions_on_rhs() {

	if(geometry.get_boundary_conditions().have_ignore_partially()) {
		
		const unsigned int num_equations = this->problem.get_num_equations_per_node();
		const BoundaryCondition& bc      = geometry.get_boundary_conditions();
		unsigned int erased = 0;		
		// for all nodes 
		for(unsigned int ii = 0; ii < geometry.get_num_nodes(); ii++) {
			const Node& node = geometry.get_node(ii);			
			// if node is non-pure dirichlet
			if(node.get_index_internal() != -1) {
				// check all dofs for dirichlet type
				for(unsigned int ee = 0; ee < num_equations; ee++) {  
					if(bc.ignore_node_partially(node,ee)) {
						// statistics
						erased++;
						// set row and column to zero
						stiff->set_row_and_column_to_zero(num_equations * node.get_index_internal() + ee);
						// set diagonal to one
						stiff->set(num_equations * node.get_index_internal() + ee, num_equations * node.get_index_internal() + ee, 1.0);
						// set rhs line to zero
						load[num_equations * node.get_index_internal() + ee] = 0.0;		
					}
				}
			}	
		}
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "FEMSolverLE: enforced " << erased << " dirichlet b.c. on dof level.");
	}	
}

} // end of namespace

#endif /*FEMSOLVERLE_H_*/
