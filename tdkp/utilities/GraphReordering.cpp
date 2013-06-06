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

#ifndef NOMETIS

#include "tdkp/utilities/GraphReordering.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/main/CSRMatrix.h"
#include "tdkp/io/InputParser.h"

extern "C" {
#include <metis.h>
}
namespace tdkp {

void GraphReordering::calculate_reordering(const Geometry& rgeometry, SparseMatrixProperties matrix_type) {

	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "GraphReordering: calculating reordering");

	const Geometry* geometry = &rgeometry;
	permutations.assign(geometry->get_num_nonzero_nodes(), 0);
	CSRMatrix<double> matrix(geometry->get_num_nonzero_nodes(), matrix_type);
	bool symmetric = matrix.symmetric();
	#pragma omp parallel default(shared)
	{
		int from, until;
		bool my_rows;
		int index_internals[Element::max_num_nodes];
		int nonzero_nodes;
		get_omp_range(from, until, geometry->get_num_nodes());
		for(Geometry::element_const_iterator it = geometry->elements_begin(); 
			it != geometry->elements_end(); it++) {
			// build matrix
			my_rows = false;
			TDKP_ASSERT((*it)->get_num_nodes() <= (unsigned int)Element::max_num_nodes, "(*it)->get_num_nodes() < Element::max_num_nodes (hard limit)");
			nonzero_nodes = 0;  			
			for(unsigned int ii = 0; ii < (*it)->get_num_nodes(); ii++) {
				if((*it)->get_node(ii).get_index_internal() != -1) {
					index_internals[nonzero_nodes] = (*it)->get_node(ii).get_index_internal();
					if(from <= index_internals[nonzero_nodes] && index_internals[nonzero_nodes] < until) {
						my_rows = true;
					}
					nonzero_nodes++;
				}		
			}
			if(my_rows) {
				for(int ii = 0; ii < nonzero_nodes; ii++) {
					if(from <= index_internals[ii] && index_internals[ii] < until) {
						for(int jj = 0; jj < nonzero_nodes; jj++) {
							if(!symmetric || index_internals[ii] <= index_internals[jj]) {
								matrix.announce(index_internals[ii], index_internals[jj]);
							}
						}
					}
				} 					
			}			
		}
	}
	matrix.set_structure();
	{	
		// pass graph to metis 
		int n = matrix.get_size();
		int numflag = 0;
		int options[8] = {1, 1, 2, 1, 0, 1, 0, 5}; // default options
		permutations.assign(n, 0);
		vector<int> iperm(n,0);
		vector<int> prow(n + 1);
		vector<int> icol(matrix.get_num_nonzeros() - n); // without diagonal
		//TDKP_ASSERT(sizeof(idxtype) == sizeof(int), "sizeof(idxtype) == sizeof(int)");
		// extract graph without diagonal (metis needs that)
		int from, until;
		int* matrix_prow = matrix.get_prow();
		int* matrix_icol = matrix.get_icol();
		TDKP_ASSERT(matrix.get_fidx() == 0, "matrix.get_fidx() == 0");
		int current_idx = 0; 
		for(int ii = 0; ii < n; ii++) {
			prow[ii] = current_idx;
			TDKP_BOUNDS_ASSERT(current_idx == matrix_prow[ii] - ii, "current_idx == matrix_prow[ii] - ii");			
			until    = matrix_prow[ii + 1];
			from     = matrix_prow[ii];
			for(int jj = from; jj < until; jj++) {
				if(matrix_icol[jj] != ii) {
					icol[current_idx] = matrix_icol[jj];
					current_idx++;
				}			
			}
		}		
		prow[n] = matrix_prow[n] - n;
		TDKP_ASSERT(current_idx == prow[n], "current_idx == prow[n]"); 
		
		if(Configuration::get_instance()->get("assembly_graph_reordering_algorithm") == 0) {		
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, "GraphReordering: starting METIS_NodeND");
			METIS_NodeND(&n, &prow[0], &icol[0], &numflag,
				options, &iperm[0], &permutations[0]);			 	
		} else {

			// ------------------------------------------------------
			// ------------------------------------------------------
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, "GraphReordering: starting METIS_PartGraphRecursive");
			int num_parts = Configuration::get_instance()->get("assembly_graph_reordering_num_partitions");
			TDKP_ASSERT(num_parts > 0, "set config option assembly_graph_reordering_num_partitions to the desired number of partitions!");			
						
			int size      = matrix.get_size();		
			int wgtflag   = 0; // no weights
			int numflag   = 0; // C ordering
			int options[] = {0,0,0,0,0}; // default options
			int edgecut   = 0; // output
			vector<int> partition;
			partition.resize(size);
			METIS_PartGraphRecursive(
				&size,				// matrix size
				&prow[0],  			// xadjacency (row ptrs)
				&icol[0],			// adjacency (icol indices)
				0,					// vwgt = vertex weights (none)
				0,					// adjwgt edge weights (none)
				&wgtflag,			// weights flag (0 here)
				&numflag,			// C ordering in icol / prow
				&num_parts,			// number of processes == partitions
				options,			// metis options
				&edgecut,			// output
				&partition[0]		// partition
			);
			
			// ----------------------------------------
			// calculate the number of elements per partition
			// ----------------------------------------
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "GraphReordering: determining permuation from graph ordering");
			vector<int> stats(num_parts, 0);
			for(unsigned int ii = 0; ii < partition.size(); ii++) {
				stats[partition[ii]]++;	
			}
			// build offsets
			vector<int> idxs(num_parts,0);
			for(int ii = 1; ii < num_parts; ii++) {
				idxs[ii] = idxs[ii - 1] + stats[ii - 1];
			}
			
			// ----------------------------------------
			// calculate permutation
			// ----------------------------------------
			for(unsigned int ii = 0; ii < partition.size(); ii++) {
				permutations[ii] = idxs[partition[ii]];
				TDKP_BOUNDS_ASSERT(permutations[ii] < (int)partition.size(), "");
				idxs[partition[ii]]++;				
			}
		} 
		// test permutation
		vector<bool> testing(permutations.size(), false);
		for(unsigned int ii = 0; ii < permutations.size(); ii++) {
			TDKP_ASSERT(permutations[ii] >= 0 && permutations[ii] < (signed)permutations.size(), "permutations[ii] >= 0 && permutations[ii] < permutations.size()");
			TDKP_ASSERT(!testing[permutations[ii]], "!testing[permutations[ii]]");
			testing[permutations[ii]] = true; 
		}		 
	}
	
	
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "GraphReodering: calculating reordering finished");

}

void GraphReordering::store_reordering(const char* filename) const {

	if(permutations.size() == 0) {
		TDKP_GENERAL_EXCEPTION("no graph calculated yet!");	
	}
	ofstream fout(filename, ios::binary);
	if(fout) {
		int n = permutations.size();
		fout.write((char*)&n, sizeof(int));
		fout.write((char*)&permutations[0], sizeof(int) * n);
		int check = 7041978;
		fout.write((char*)&check, sizeof(int));
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to file " << filename);	
	} 
		
}
void GraphReordering::load_reordering(const char* filename) {
	ifstream fin(filename, ios::binary);
	if(fin) {
		int n;
		fin.read((char*)&n, sizeof(int));
		TDKP_ASSERT(n > 0, "n > 0");
		permutations.assign(n, 0);
		fin.read((char*)&permutations[0], sizeof(int) * n);
		int check;
		fin.read((char*)&check, sizeof(int));
		fin.close();
		TDKP_ASSERT(check == 7041978, "check == 7041978"); 
		fin.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not read from file " << filename);	
	} 			
}


const vector<int>& GraphReordering::get_reordering() const {
	return permutations;	
}

void GraphReordering::apply_reordering(Geometry& geometry) {
	// ------------------------------------------
	// use node reordering on request
	// ------------------------------------------
	int symmetric = static_cast<int>(Configuration::get_instance()->get("assembly_build_nonsymmetric_matrices")) != 1.0;
	// check if ordering is already available
	ostringstream sout;
	sout << "cache_ordering_" << geometry.get_identifier() << "_" << symmetric 
	     << "_" << geometry.get_num_nonzero_nodes() << ".dat";
	ifstream fin(sout.str().c_str(), ios::binary);
	GraphReordering reordering;
	if(fin) {
		ostringstream mout; 
		mout << "loading and applying reordering from file " << sout.str();
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, mout.str());
		fin.close();
		reordering.load_reordering(sout.str().c_str());
	} else {
		fin.close();
		if(symmetric) {
			TDKP_GENERAL_EXCEPTION("GraphReordering is buggy for symmetric stuff (metis does not work ...)"); 
			reordering.calculate_reordering(geometry, symmetric_matrix);
		} else {
			reordering.calculate_reordering(geometry, nonsymmetric_matrix);
		}		
		//reordering.calculate_reordering(geometry.get_identifier().c_str());
		reordering.store_reordering(sout.str().c_str());
	}
	const vector<int>& reorder = reordering.get_reordering();
	for(Geometry::node_iterator it = geometry.nodes_begin(); it != geometry.nodes_end(); it++) {
		if((*it)->get_index_internal() != -1) {
			(*it)->set_index_internal(reorder[(*it)->get_index_internal()]);	
		}	
	}
}


}
#endif
