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

#ifndef FEMSOLVERGEVP_H_
#define FEMSOLVERGEVP_H_

#include <fstream>
#include <vector>
#include "omp.h"
#include "tdkp/common/all.h"
#include "tdkp/common/Configuration.h"
#include "tdkp/main/FEMSolver.h"
#include "tdkp/main/SparseMatrixInterface.h"
#include "tdkp/probdefs/EigenProblem.h"
#include "tdkp/io/InputParser.h"
#include "tdkp/solvers/EigenSolver.h"
#include "tdkp/geometry/BoundaryCondition.h"

namespace tdkp {
	
				
template<class T, class TMat, class TRHS>				
class FEMSolverGEVP : public FEMSolver {
	
public:

	typedef SparseMatrixInterface<TMat> SMatrix;
	typedef SparseMatrixInterface<TRHS> OMatrix;

	FEMSolverGEVP(const Geometry& geometry_, EigenProblem<T, TMat, TRHS>& problem_);
	virtual ~FEMSolverGEVP();
	// ---------------------------------------------------
	// simulation setup functions 
	// ---------------------------------------------------
	bool verify_setup() const;		
	void set_ordering(Ordering order_) { this->order = order_; }
	
	// ---------------------------------------------------
	// assembling 
	// ---------------------------------------------------
	void create_matrix_structures();
	void reset_matrix_to_zero(); 
	void assemble_system();
	void solve_system(int num_solutions);

	const SMatrix& get_lhs_matrix() const;
	const OMatrix& get_rhs_matrix() const; 
																
			
protected:

	virtual EigenSolver<TMat,TRHS,T>* create_eigensolver_object(unsigned int matrix_size, unsigned int block_size); 

    EigenProblem<T,TMat,TRHS>&    problem; 
	EigenSolver<TMat,TRHS,cplx>*  eigensolver;
	      	
	SMatrix*  stiff;
	OMatrix*  overlap;
	Ordering  order;
	
private:		
								 			
};


	

using namespace tdkp;
using namespace std;


template<class T, class TMat, class TRHS>
FEMSolverGEVP<T,TMat,TRHS>::FEMSolverGEVP(const Geometry& geometry_, EigenProblem<T,TMat,TRHS>& problem_) 
: FEMSolver(geometry_), 
  problem(problem_), 
  eigensolver(0), 
  stiff(0), 
  overlap(0), 
  order(descending) 
{
	
}

template<class T, class TMat, class TRHS>
FEMSolverGEVP<T,TMat,TRHS>::~FEMSolverGEVP() {
	if(this->eigensolver != 0) {
		delete this->eigensolver;
		this->eigensolver = 0;
	}
}

/** return left hand side matrix */
template<class T, class TMat, class TRHS>
const typename FEMSolverGEVP<T,TMat,TRHS>::SMatrix& FEMSolverGEVP<T,TMat,TRHS>::get_lhs_matrix() const {
	TDKP_ASSERT(this->stiff != 0, "stiffnes matrix not set  yet");
	return *stiff;		
}

/** return right hand side matrix */
template<class T, class TMat, class TRHS>
const typename FEMSolverGEVP<T,TMat,TRHS>::OMatrix& FEMSolverGEVP<T,TMat,TRHS>::get_rhs_matrix() const {
	TDKP_ASSERT(this->overlap != 0, "overlap matrix not set yet");
	return *overlap;	
}

template<class T, class TMat, class TRHS>
bool FEMSolverGEVP<T,TMat,TRHS>::verify_setup() const {

	if(!this->geometry.verify()) {
		Logger::get_instance()->emit(LOG_ERROR, "geometry is invalid");
		return false;	
	}	
	return true;
}

template<class T, class TMat, class TRHS>
EigenSolver<TMat,TRHS,T>* FEMSolverGEVP<T,TMat,TRHS>::create_eigensolver_object(unsigned int matrix_size, unsigned int block_size) {
	return EigenSolver<TMat,TRHS,T>::factory(matrix_size, block_size);
}

template<class T, class TMat, class TRHS>
void FEMSolverGEVP<T,TMat,TRHS>::create_matrix_structures() {
		
	TimeMeasurements::get_instance().start("matrix graph building");		

	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");	
	}

	// -------------------------------------------------------------
	// create matrices and set allocate memory
	// -------------------------------------------------------------
	if(this->eigensolver != 0) {
		delete this->eigensolver; this->eigensolver = 0;	
	}
	
    int	   	   sparse_num;
    const int* sparse_pat    = this->problem.get_node_sparsity_pattern(sparse_num);	
	const int  num_equations = this->problem.get_num_equations_per_node();	
		
	this->eigensolver = this->create_eigensolver_object(this->geometry.get_num_nonzero_nodes() * num_equations, num_equations);
	this->stiff       = &(this->eigensolver->get_stiff());
	this->overlap     = &(this->eigensolver->get_overlap());
	
	// setting the block size does not mean that we will also assemble block wise
	// the matrix will later tell me if it likes to have the blocks treated 
	// as blocks or as single entries
	this->stiff->set_block_size(num_equations);
	this->overlap->set_block_size(num_equations);

	// -------------------------------------------------------
	// load structure (if available)
	// -------------------------------------------------------
	if(Configuration::get_instance()->get("assembly_use_matrix_structure_cache") == 1.0) {
		if(this->stiff->caching_structure_possible() && this->overlap->caching_structure_possible()) {
			ostringstream sout_stiff, sout_overlap, sout;
			sout_stiff   << "stiff_" << this->problem.get_unique_identifier() << "_" << this->geometry.get_identifier();
			sout_overlap << "overlap_" << this->problem.get_unique_identifier() << "_" << this->geometry.get_identifier();
			if(this->stiff->structure_cache_available(sout_stiff.str()) && 
			   this->overlap->structure_cache_available(sout_overlap.str())) {
			   	sout << "FEMSolverGEVP: found structure caches in files. loading them";
			   	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
				this->stiff->load_structure(sout_stiff.str());
				this->overlap->load_structure(sout_overlap.str());
				TimeMeasurements::get_instance().stop("matrix graph building");
				return;	
			}	
		}	
	}
	

	// check if rhs matrix will be in reduced form 
	// the blocks on rhs are diagonal with equal values!
	bool reduced_rhs = false;
	if(this->stiff->get_size() / this->overlap->get_size() == (unsigned)num_equations) {
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverGEVP: rhs matrix is treated as having diagonal blocks with equal entries");
		reduced_rhs = true;		
	}	

#pragma omp parallel default(shared)
{
	int        gii,gjj;
	bool       symmetric;
		
	symmetric = stiff->property_is_set(symmetric_matrix);	

	int processor_first;
	int processor_last;
	
	get_omp_matrix_range(processor_first, processor_last, this->geometry.get_num_nonzero_nodes(), symmetric);

	if(omp_get_num_threads() > 1 && omp_get_thread_num() == 0) {
		ostringstream sout;
		sout << "FEMSolverGEVP: linewise parallel building of "
		     << (symmetric ? "symmetric" : "nonsymmetric") 
		     << " matrix structure with "
		     << omp_get_num_threads() << " active threads "
		     << "where the average processor takes care of " 
		     << (processor_last - processor_first) * num_equations
		     << " matrix lines.";
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());				
	}
	      
		
	// set iterators
	Geometry::element_const_iterator eit;
	Geometry::element_const_iterator end = this->geometry.elements_end();

	int next = 0;
	int cc   = 0;
	if(Logger::get_instance()->get_level() > LOG_INFO && omp_get_thread_num() == 0) {
		Logger::get_instance()->init_progress_bar("FEMSolverGEVP: building matrix structure", this->geometry.get_num_elements());	
	}
	
	// loop over elements to build matrix structure
	for(eit = this->geometry.elements_begin(); eit != end; eit++) {
		// only process enabled elements (non-contacts)
		if((*eit)->enabled()) {				
			for(int ii = 0; ii < (signed)(*eit)->get_num_nodes(); ii++) {			
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
									if(!(gii > gjj)) {
										stiff->block_announce(gii,gjj, sparse_pat, sparse_num);																
										if(reduced_rhs) {
											// overlap matrix has small unity diagonal, therefore saving only block structure										
											overlap->announce(gii,gjj);
										} else {
											// announce block diagonal
											for(int cc = 0; cc < num_equations; cc++) {
												overlap->announce(gii * num_equations + cc, gjj * num_equations + cc);	
											}
										}
									} 
								} else { // if(symmetric)
								// -----------------------------------						
								// END SYMMETRIC / BEGIN NONSYMMETRIC
								// -----------------------------------													
									stiff->block_announce(gii,gjj, sparse_pat, sparse_num);							
									if(reduced_rhs) {
										// overlap matrix has small unity diagonal, therefore saving only block structure									
										overlap->announce(gii,gjj);
									} else {
										// announce block diagonal
										for(int cc = 0; cc < num_equations; cc++) {
											overlap->announce(gii * num_equations + cc, gjj * num_equations + cc);	
										}
									}													
								}
							}
						}
					}											
				}
			}
		}
		if(next == cc) {
			if(Logger::get_instance()->get_level() > LOG_INFO && omp_get_thread_num() == 0) {
				next = Logger::get_instance()->set_progress_bar(cc, this->geometry.get_num_elements());
			}	
		}
		cc++;
	}
	if(Logger::get_instance()->get_level() > LOG_INFO && omp_get_thread_num() == 0) {		
		Logger::get_instance()->end_progress_bar();
	}
			
} // end of parallel loop	
	

	TimeMeasurements::get_instance().track_memory_usage();	
	
				
	// compile structure
	stiff->set_structure();			
	overlap->set_structure();

	TimeMeasurements::get_instance().stop("matrix graph building");
	TimeMeasurements::get_instance().track_memory_usage();

	// -------------------------------------------------------
	// save structure (if available and requested)
	// -------------------------------------------------------
	if(Configuration::get_instance()->get("assembly_use_matrix_structure_cache") == 1.0) {
		if(this->stiff->caching_structure_possible() && this->overlap->caching_structure_possible()) {
			ostringstream sout_stiff, sout_overlap, sout;
			sout_stiff   << "stiff_" << this->problem.get_unique_identifier() << "_" << this->geometry.get_identifier();
			sout_overlap << "overlap_" << this->problem.get_unique_identifier() << "_" << this->geometry.get_identifier();
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverGEVP: writing matrix structure to cache");
			this->stiff->save_structure(sout_stiff.str());
			this->overlap->save_structure(sout_overlap.str());
		} else {
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, "FEMSolverGEVP: caching matrix structure is not supported for given matrix type");
		}
	}			
		
}

template<class T, class TMat, class TRHS>
void FEMSolverGEVP<T,TMat,TRHS>::reset_matrix_to_zero() {
	stiff->clear_but_keep_structure();
	overlap->clear_but_keep_structure(); 	
}
	
template<class T, class TMat, class TRHS>
void FEMSolverGEVP<T,TMat,TRHS>::assemble_system() {

	TimeMeasurements::get_instance().start("matrix assembly");

	if(this->overlap == 0 || this->stiff == 0) {
		TDKP_GENERAL_EXCEPTION("matrice structures must be build before assembly.");	
	}

	if(!this->problem.is_ready()) {	
		this->problem.prepare();
		TDKP_ASSERT(this->problem.is_ready(), "class ProblemDefintion must be ready after prepare() is called");
	}
	
    // sparse pat
    int  	   sparse_num;
	const int  num_equations = this->problem.get_num_equations_per_node();
	// (max num node per element is statically set in Element.h)
	// ---------------------------
	// WARNING! DON'T CHANGE nlhs without adjust it in SchroedingerPML!
	// ---------------------------
	const int  nlhs          = num_equations * num_equations * Element::max_num_nodes * Element::max_num_nodes; 
    const int* sparse_pat    = this->problem.get_node_sparsity_pattern(sparse_num);			
	const bool symmetric     = stiff->property_is_set(symmetric_matrix);

	// check if rhs is in reduced size format
	bool reduced_rhs = false;
	if(this->stiff->get_size() / this->overlap->get_size() == (unsigned)num_equations) {
		reduced_rhs = true;		
	}	
	
	// ----------------------------------------------
	// parallel section
	// ----------------------------------------------
#pragma omp parallel default(shared)
	{	
	
		// iterators
		Geometry::element_const_iterator eit; 
		Geometry::element_const_iterator end;
		end = this->geometry.elements_end();
		// -----------------------------------------------------	
		// some space to work on  
		// size controlled via static const int Element::max_num_nodes
		// -----------------------------------------------------
		TMat*   lhs = new TMat[nlhs];
		TRHS*   rhs = new TRHS[nlhs];
		int   	idxs[Element::max_num_nodes]; 
		int   	n  = 0;
		int   	gii, gjj;
		int   	index_internals[Element::max_num_nodes];
		int     current_iteration = 0;
		int     next              = 0;
						
		// --------------------------------------------------
		// processor distribution
		// --------------------------------------------------
		int processor_first;
		int processor_last;
		
		get_omp_matrix_range(processor_first, processor_last, this->geometry.get_num_nonzero_nodes(), symmetric);

		TDKP_ASSERT(processor_first < processor_last, "processor_first (" << processor_first << ") < processor_last (" << processor_last << ")"); 

		if(omp_get_num_threads() > 1 && omp_get_thread_num() == 0) {
			ostringstream sout;
			sout << "FEMSolverGEVP: linewise parallel assembly of matrix with "
			     << omp_get_num_threads() << " active threads "
			     << "where the average processor takes care of " 
			     << processor_last - processor_first
			     << " matrix lines.";			     
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());				
		}					
		if(Logger::get_instance()->get_level() > LOG_INFO && omp_get_thread_num() == 0) {
			Logger::get_instance()->init_progress_bar("FEMSolverGEVP: assembling matrices", this->geometry.get_num_elements());
		}
					
		// -----------------------------------------------------------
		// loop over elements and assemble
		// -----------------------------------------------------------
		bool my_rows = false;
		int  index_internal;		
		
		for(eit = this->geometry.elements_begin(); eit != end; eit++) {
			
			// ---------------------------------------------
			// only process element if its enabled
			// ---------------------------------------------
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
					TDKP_ASSERT(n <= Element::max_num_nodes, "n <= Element::max_num_nodes (fixed size array)"); 		
			
					// get internal indices
					for(int ii = 0; ii < n; ii++) {
						index_internals[ii] = (*eit)->get_node(idxs[ii]).get_index_internal();	
					}		
					// for all nodes in element			 
					for(int ii = 0; ii < n; ii++) {
						gii = index_internals[ii];
						if(processor_first <= gii && gii < processor_last) {
							// for all nodes in element (again ;-))
							for(int jj = 0; jj < n; jj++) {	
								gjj = index_internals[jj];
								if(symmetric) {
									if(gii < gjj) {
										// -----------------------------------------------------------------
										// loop over all local matrix elements (for lhs)
										// -----------------------------------------------------------------				
										for(int ss = 0; ss < sparse_num; ss++) {
											stiff->add(gii * num_equations + sparse_pat[2*ss], gjj * num_equations + sparse_pat[2*ss+1], lhs[(n * ii + jj) * sparse_num + ss]);				
										}
										if(reduced_rhs) {
											// -----------------------------------------------------------------
											// overlap has diagonal local matrix -> rhs has no multiple copies
											// and we only save it once to save memory
											// -----------------------------------------------------------------										
											overlap->add(gii,gjj, rhs[n * ii + jj]);
										} else {
											// -----------------------------------------------------------------
											// overlap is in full form
											// -----------------------------------------------------------------
											for(int cc = 0; cc < num_equations; cc++) {
												overlap->add(gii * num_equations + cc, gjj * num_equations + cc, rhs[n * ii + jj]);
											}																					
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
												stiff->add(lii,ljj, lhs[(n * ii + jj) * sparse_num + ss]);				
											}
										}
										if(reduced_rhs) {
											// -----------------------------------------------------------------
											// overlap has diagonal local matrix -> rhs has no multiple copies
											// and we only save it once to save memory
											// -----------------------------------------------------------------										
											overlap->add(gii,gjj, rhs[n * ii + jj]);
										} else {
											// -----------------------------------------------------------------
											// overlap is in full form
											// -----------------------------------------------------------------
											for(int cc = 0; cc < num_equations; cc++) {
												overlap->add(gii * num_equations + cc, gjj * num_equations + cc, rhs[n * ii + jj]);
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
										stiff->add(gii * num_equations + sparse_pat[2*ss], gjj * num_equations + sparse_pat[2*ss+1], lhs[(n * ii + jj) * sparse_num + ss]);				
									}
									if(reduced_rhs) {
										// -----------------------------------------------------------------
										// overlap has diagonal local matrix -> rhs has no multiple copies
										// and we only save it once to save memory
										// -----------------------------------------------------------------										
										overlap->add(gii,gjj, rhs[n * ii + jj]);
									} else {
										// -----------------------------------------------------------------
										// overlap is in full form
										// -----------------------------------------------------------------
										for(int cc = 0; cc < num_equations; cc++) {
											overlap->add(gii * num_equations + cc, gjj * num_equations + cc, rhs[n * ii + jj]);
										}																					
									}								
								}
							} //  end for(int jj = 0; jj < n; jj++) {	
						} // end if(processor_first <= gii && gii < processor_last) {					
					} // end for(int ii = 0; ii < n; ii++) {
				} // end if(myrows) {
			} // end if(((*eit)->enabled)
			if(omp_get_thread_num() == 0 && Logger::get_instance()->get_level() > LOG_INFO) {
				if(current_iteration > next) {
					next = Logger::get_instance()->set_progress_bar(current_iteration, this->geometry.get_num_elements());
				}	
				current_iteration++;
			}
		} // end loop over elements
		
		// ---------------------------------
		// destroy parallel allocated memory
		// ---------------------------------
		delete[] lhs;
		delete[] rhs;
			
	} // end omp parallel section
	if(Logger::get_instance()->get_level() > LOG_INFO) {
		Logger::get_instance()->end_progress_bar();
	}
	
	// -----------------------------------------------------
	// plug in partial dirichlet boundary conditions
	// -----------------------------------------------------
	if(geometry.get_boundary_conditions().have_ignore_partially()) {
		TDKP_GENERAL_EXCEPTION("sorry, but for eigenproblems, partial enforcing of dirichlet b.c. wasn't implemented as i don't think that this would be reasonable.");	
	}
		
	TimeMeasurements::get_instance().stop("matrix assembly");
	
	// store matrices to files upon request
	if(Configuration::get_instance()->get("assembly_save_matrices_to_file") != 0.0) {
		stiff->save_to_file("lhs.mat");			
		overlap->save_to_file("rhs.mat");
	}

	// check matrices for symmtery upon request
	if(Configuration::get_instance()->get("assembly_check_matrix_for_symmetry") != 0.0) {	
		stiff->perform_symmetry_analysis();
		// should we dump nonsymmetric terms?
		if(Configuration::get_instance()->get("assembly_track_nonhermitian_nodes") != 0.0) {
			const vector<double>& track_nodes = stiff->get_nonhermitian_rows();
			if(track_nodes.size() > 0) {
				TDKP_ASSERT(num_equations * this->geometry.get_num_nonzero_nodes() == track_nodes.size(), "lhs matrix does not have the needed size!");
				StdNodeData<double> nodes_output(num_equations + 1, this->geometry.get_num_nodes()); 	
				for(Geometry::node_const_iterator it = this->geometry.nodes_begin(); it != this->geometry.nodes_end(); it++) {
					if((*it)->get_index_internal() != -1) {    
						double sum = 0;
						for(int jj = 0; jj < num_equations; jj++) {
							nodes_output.set_node_value((*it)->get_index_global(), jj, track_nodes[(*it)->get_index_internal() * num_equations + jj]);
							sum += track_nodes[(*it)->get_index_internal() * num_equations + jj];
						}
						nodes_output.set_node_value((*it)->get_index_global(), num_equations, sum);
					}					 
				}
				InputParser parser;
				Logger::get_instance()->emit(LOG_INFO, "FEMSolverGEVP: wrote data about nonhermitian values at nodes to file named lhs_nonhermitian_nodes.dat"); 
				parser.write_binary(nodes_output, "lhs_nonhermitian_nodes.bin");
			}
			stiff->delete_nonhermitian_rows_data(); 
		}				
	}
	

	
}


		
} // end of namespace


#endif /*ifdef FEMSOLVERGEVP_H_ */
