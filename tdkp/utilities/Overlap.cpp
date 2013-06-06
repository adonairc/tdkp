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

#include "Overlap.h"

namespace tdkp {

const int Overlap::sparsity_pattern[2] = {0,0};	
	
/** constructor, directly builds fem overlap matrix */
Overlap::Overlap(Geometry& geometry_, MaterialDatabase& matdb_, SparseMatrixProperties matrix_type)
: LinearProblem<double>(geometry_, matdb_),
  fem_assembler(0)
{
	geometry_.set_boundary_conditions(new BCIncludeAll(geometry_));
	fem_assembler = new FEMSolverME<double>(geometry_, *this, matrix_type);	
	fem_assembler->create_matrix_structures();
	fem_assembler->assemble_system();
}

Overlap::~Overlap() {
	delete fem_assembler;
}

/** evaluate overlap (thread safe) */
cplx Overlap::evaluate(const EigenSolution<cplx>& v1, const EigenSolution<cplx>& v2) const {
	
	TDKP_ASSERT(v1.get_basis_size() == v2.get_basis_size(), "v1 and v2 have not the same basis size. therefore the overlap should be 0 anyway. so i refuse to calculate something where the result is already known.");
	TDKP_ASSERT(v1.get_length() == v2.get_length(), "v1 and v2 are not of the same length!");
	TDKP_BOUNDS_ASSERT(v1.get_length() == (signed)geometry.get_num_nonzero_nodes(), "");

	vector<cplx> tmp(v1.get_length());
	vector<cplx> v2a(v1.get_length());		
	const vector<cplx>& v1_data = v1.get_data();
	const vector<cplx>& v2_data = v2.get_data();
	unsigned int basis_size = v1.get_basis_size();		
	// --------------------------------------
	// for every envelope
	// --------------------------------------
	cplx ret(0.0, 0.0);
	for(unsigned int aa = 0; aa < basis_size; aa++) {
		// -------------------------------------
		// copy envelope of bloch function to vec
		// -------------------------------------
		for(unsigned int ii = 0; ii < v2a.size(); ii++) {		
			v2a[ii] = v2_data[ii * basis_size + aa];	
		}
		// -------------------------------------
		// calculate M v2a -> tmp
		// -------------------------------------
		fem_assembler->multiply_with_lhs(v2a, tmp);
		// -------------------------------------
		// dot product
		// -------------------------------------
		for(unsigned int ii = 0; ii < v2a.size(); ii++) {
			ret += conj(v1_data[ii * basis_size + aa]) * tmp[ii];	
		}
	} 

	return ret;
	
}
 	
const int* Overlap::get_node_sparsity_pattern(int& n) const {
	n = 1;
	return this->sparsity_pattern;
} 	
 	
void Overlap::calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* internal_idx, int &n) const {
	
	int    nnode;              /* number of nodes */
		
	n     = 0;
	nnode = (signed)elem->get_num_nodes();
	

	//TDKP_ASSERT(nnode <= Element::max_num_nodes, "nnode <= 4 (working arrays to small ...)");

	// ------------------------------------------------------------------------------
	// build node_internal_indices (nonzero nodes, which we really calculate)
	// ------------------------------------------------------------------------------
	for(int ii = 0; ii < nnode; ii++) {
		if(elem->get_node(ii).get_index_internal() != -1) {
			internal_idx[n++] = ii;
		}
	}

	// ------------------------------------------------------------------------------
	// calculate local mass matrix for nonzero indices
	// ------------------------------------------------------------------------------
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {
			lhs[ii * n + jj] = elem->get_element_integral_0th_order(internal_idx[ii],internal_idx[jj]);
			rhs[ii * n + jj] = 0.0;
		}
	}
}	

/** constructor for quantized systems */
DeCross::DeCross(Geometry& geometry, MaterialDatabase& matdb)
: overlap(0)
{
	overlap = new Overlap(geometry, matdb);
}

/** constructor for bulk */
DeCross::DeCross()
: overlap(0)
{
	
}

/** destructor */
DeCross::~DeCross() {
	// delete object if it exists ... 
	if(overlap != 0) {
		delete overlap; overlap = 0;	
	}	
}

/** calculate overlap between wavefunctions */
cplx DeCross::get_overlap(const EigenSolution<cplx>& lhs, const EigenSolution<cplx>& rhs) const {

	// ------------------------------------------------
	// quantized systems?
	// ------------------------------------------------
	if(overlap) {
		return overlap->evaluate(lhs,rhs);
	} else {
		// ------------------------------------------------
		// bulk		
		// ------------------------------------------------
		TDKP_ASSERT(lhs.get_length() == 1, "thats not a bulk eigensolutions");
		TDKP_ASSERT(rhs.get_length() == 1, "thats not a bulk eigensolutions");
		TDKP_ASSERT(lhs.get_num_data_per_node() == rhs.get_num_data_per_node(), "");
		cplx ret = 0.0;
		for(int ii = 0; ii < lhs.get_num_data_per_node(); ii++) {
			ret += conj(lhs.get_node_value(0,ii))* rhs.get_node_value(0,ii);		
		}
		return ret;			
	}
}

void DeCross::resort(Bandstructure<complex<double> >& bands) const {

	TDKP_ASSERT(bands.get_basis_size() > 1, "resorting does not make sense for effmass bands ...");

	// thats the energy difference where we accept the states to be degenerate
	const double energy_tolerance = Configuration::get_instance()->get("state_resorting_degeneracy_delta_energy_threshold");
	// thats the min non-diagonality between two kidx of a subband. if its below, we guess that it is a accidential crossing and resort 
	const double overlap_tolerance = Configuration::get_instance()->get("state_resorting_overlap_diagonality_threshold");
	// thats the min non-diagonality in the resorting algorithm when we reaccept two states as beeing the right continuation in k-space
	const double accept_overlap_tolerance = Configuration::get_instance()->get("state_resorting_overlap_accept_threshold");
	// min orthogonality between d0d1 upon renormalization
	const double orthogonality_threshold = Configuration::get_instance()->get("state_resorting_minimum_orthogonality_threshold");
	// be quiet or be chatty
	const bool be_chatty = Configuration::get_instance()->get("output_decrossing_details") == 1;
	
	// -------------------------------------------------
	// 1. check for spin degeneracy
	// -------------------------------------------------	
	TDKP_ASSERT(bands.get_number_of_bands() % 2 == 0, "Decross actually expects an even number of spin-degenerate bands, but band structure object has " << bands.get_number_of_bands() << " bands ... so calculate an even number of bands or disable decrossing");
	vector<bool> degenerate(bands.get_number_of_bands() / 2, true);
	bool all_degenerate = true;	
	for(int ii = 0; ii < bands.get_number_of_bands(); ii += 2) {
		for(int kk = 0; kk < bands.get_number_of_k_values(); kk++) {
			if(tdkp_math::abs(bands.get_energy(kk,ii) - bands.get_energy(kk, ii + 1)) > energy_tolerance) {
				degenerate[ii / 2] = all_degenerate = false;
				break;	
			}
		}		
	}
	
	// -------------------------------------------------
	// inform user -> if no band is spin degenerate
	// just stop and return
	// -------------------------------------------------
	if(all_degenerate) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "DeCross: all bands are spin degenerate");	
	} else {		
		bool some_degenerate = false;
		for(int ii = 0; ii < bands.get_number_of_bands() / 2; ii++) {
			if(degenerate[ii]) {
				some_degenerate = true;
				break;	
			}
		}
		if(!some_degenerate) {
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "DeCross: no bands are spin degenerate. nothing to do.");			
		} else {
			// ------------------------------------------------
			// some are degenerate some are not, i am confused and 
			// don't know what to do ... 
			// ------------------------------------------------
			TDKP_LOGMSG(LOG_WARN, "DeCross: problems, some bands seem to be spin degenerate and some bands aren't ... don't know how to handle it.");			
		}
		return;				
	}
				
	// -------------------------------------------------
	// create reference kindices 
	// -------------------------------------------------
	map<unsigned int, unsigned int> kmap;
	kmap[0] = 1;
	//kmap[1] = 0;
	for(int kk = 2; kk < bands.get_number_of_k_values(); kk++) {
		kmap[kk] = kk - 1; 	
	}
	vector<double> decrossing_info(bands.get_number_of_k_values() * ( bands.get_number_of_bands() / 2), 0.0);						
										
	// -------------------------------------------------
	// process each degenerate band combo
	// -------------------------------------------------
	ostringstream sout;	
	for(int ii = 0; ii < bands.get_number_of_bands() / 2; ii++) {
		map<unsigned int, unsigned int>::const_iterator it = kmap.begin();
		while(it != kmap.end()) {
			sout.str("");
			unsigned int kref = it->second;
			unsigned int kidx = it->first;
			bool keep_decrossing = true;
			
			// -------------------------------------------------
			// get eigensolutions
			// -------------------------------------------------
			const EigenSolution<cplx>* r[2];
			EigenSolution<cplx>* d[2];
			for(unsigned int vv = 0; vv < 2; vv++) {
			 	r[vv] = &bands.get_eigensolution(kref, ii * 2 + vv);
			 	d[vv] = &bands.get_eigensolution(kidx, ii * 2 + vv);	
			}
			// -------------------------------------------------
			// build rd^t matrix
			// -------------------------------------------------
			cplx matrix[4];
			sout << "-----------------------------------------------------\n"
			     << " band combo idx " << ii << " at kidx " << kidx << " using kref " << kref << "\n"
			     << " G matrix <w_kref|w_ki>:\n "; 			
			 
			for(unsigned int vv = 0; vv < 2; vv++) {
				sout << "   ";
				for(unsigned int ww = 0; ww < 2; ww++) {
					matrix[vv * 2 + ww] = get_overlap(*r[vv],*d[ww]);
					sout << matrix[vv * 2 + ww] << " \t";  		
				}
				sout << "\n";
			}
			// -------------------------------------------------
			// get determinate of rd^t matrix
			// -------------------------------------------------
			cplx det = tdkp_math::det_2x2_matrix(matrix);							
			sout << " det G = " << det << " (abs: " << tdkp_math::abs(det) << ")\n";
									
			// ---------------------------------------------
			// if det is small, we check some other subbands
			// ---------------------------------------------
			if(tdkp_math::abs(det) <= overlap_tolerance) {
				if(kref > kidx) {
					TDKP_LOGMSG(LOG_INFO_DEVEL2, "detected band crossing of subband combo " << ii \
					<< " between k-points " << kref << " and " << kidx << ". can not fix it for reference k beeing greater than current.");						
				} else {
					TDKP_LOGMSG(LOG_INFO_DEVEL2, "detected band crossing of subband combo " << ii \
					<< " between k-points " << kref << " and " << kidx << ". will try to find the right partner.");
					// ----------------------------------------
					// loop over the rest of bands
					// ----------------------------------------
					EigenSolution<cplx>* other_bands[2];
					bool solved = false;
					for(int jj = ii + 1; jj < bands.get_number_of_bands() / 2; jj++) {
						// -----------------------------------------
						// get wavefunctions
						// -----------------------------------------
						for(unsigned int vv = 0; vv < 2; vv++) {
							other_bands[vv] = &bands.get_eigensolution(kidx, jj * 2 + vv);
						}
						// -----------------------------------------
						// calculate new overlap matrix and its det
						// -----------------------------------------
	 					for(unsigned int vv = 0; vv < 2; vv++) {					
							for(unsigned int ww = 0; ww < 2; ww++) {
								matrix[vv * 2 + ww] = get_overlap(*r[vv],*other_bands[ww]);					
							}
						}
						det = tdkp_math::det_2x2_matrix(matrix);
						TDKP_LOGMSG(LOG_INFO_DEVEL2, "combo " << jj << " has det abs = " << tdkp_math::abs(det)); 
						// --------------------------------------
						// accept if det is not too small
						// --------------------------------------
						if(tdkp_math::abs(det) > accept_overlap_tolerance) {
							TDKP_LOGMSG(LOG_INFO_DEVEL2, "subband combo " << jj << " seems to be the correct pair (|det G| = " \
							<< tdkp_math::abs(det) << ") for subband combo " << ii \
							<< ". will swap the remaining wavefunctions (in k-space).");
							map<unsigned int, unsigned int>::const_iterator swap_kit = it;
							// ----------------------------------------
							// swap for the remaining k points
							// ---------------------------------------
							while(swap_kit != kmap.end()) {
								unsigned int swap_kidx = swap_kit->first;
								// swap ii <--> jj degenerate combo
								for(unsigned int nn = 0; nn < 2; nn++) {
									bands.swap(swap_kidx, ii * 2 + nn, swap_kidx, jj * 2 + nn);	
								}							
								swap_kit++;	
							}
							// ---------------------------------------
							// refill d with new bands
							// ---------------------------------------
							for(unsigned int vv = 0; vv < 2; vv++) {
								d[vv] = &bands.get_eigensolution(kidx, ii * 2 + vv);
							}
							// don't have to update det / matrix as values are still present ...
							solved = true;
							break;
						}
					}
					if(!solved) {
						TDKP_LOGMSG(LOG_WARN, "could not fix crossing in pair " << ii << " at " <<  kidx << " vs " << kref << ". that's a serious problem because it should be possible AROUND the gamma point.");
						for(unsigned int vv = 0; vv < 2; vv++) {				
							for(unsigned int ww = 0; ww < 2; ww++) {
								matrix[vv * 2 + ww] = get_overlap(*r[vv],*d[ww]);						  		
							}
						}
						cplx det = tdkp_math::det_2x2_matrix(matrix);
						if(tdkp_math::abs(det) < 1.0e-10) {
							keep_decrossing = false;	
						} 
					}
				}
			}

			
			decrossing_info[ii * bands.get_number_of_k_values() + kidx] = tdkp_math::abs(det);
			
			// -------------------------------------------------
			// invert matrix = (A^t) (gives new coefficents
			// for the new linear combinations of our wavefunctions
			// these linear combinations should be
			// - diagonal at current k and 
			// - almost diagonal at reference k
			// -------------------------------------------------
			if(keep_decrossing) {
				tdkp_math::invert_2x2_matrix(matrix);
				
				sout << " inv(G):\n";
				for(unsigned int vv = 0; vv < 2; vv++) {
					sout << "   ";
					for(unsigned int ww = 0; ww < 2; ww++) {			
						sout << matrix[vv * 2 + ww] << " \t";  		
					}
					sout << "\n";
				}
				// -------------------------------------------------						
				// build new eigensolution set
				// -------------------------------------------------
				EigenSolution<cplx> n[2];			
				for(unsigned int vv = 0; vv < 2; vv++) {
					EigenSolution<cplx> tmp_sol;
					for(unsigned int ww = 0; ww < 2; ww++) {
						tmp_sol.add(matrix[ww * 2 + vv], *d[ww]);  	
					}
					cplx tmp = get_overlap(tmp_sol, tmp_sol);				
					n[vv].add(1.0 / sqrt(tmp.real()), tmp_sol);
				}			
				*d[0] = n[0]; 
				*d[1] = n[1];
				
				// -------------------------------------------------
				// testing
				// 1. d[0] orthogonal to d[1]
				// -------------------------------------------------
				if(be_chatty) {
					sout << "  testing:\n";
					sout << "  <d[0] |  d[1]> = " << tdkp_math::abs(get_overlap(
														bands.get_eigensolution(kidx, ii * 2),
														bands.get_eigensolution(kidx, ii * 2 + 1)	
													 )) << "\n"
						 << "  <d[0] |  c[0]> = " << tdkp_math::abs(get_overlap(
														bands.get_eigensolution(kidx, ii * 2),
														bands.get_eigensolution(kref, ii * 2)	
													 )) << "\n"
						 << "  <d[0] |  c[1]> = " << tdkp_math::abs(get_overlap(
														bands.get_eigensolution(kidx, ii * 2),
														bands.get_eigensolution(kref, ii * 2 + 1)	
													 )) << "\n"
						 << "  <d[1] |  c[0]> = " << tdkp_math::abs(get_overlap(
														bands.get_eigensolution(kidx, ii * 2 + 1),
														bands.get_eigensolution(kref, ii * 2)	
													 )) << "\n"
						 << "  <d[1] |  c[1]> = " << tdkp_math::abs(get_overlap(
														bands.get_eigensolution(kidx, ii * 2 + 1),
														bands.get_eigensolution(kref, ii * 2 + 1)	
													 )) << "\n"
						 << "  <d[0] |  d[0]> = " << tdkp_math::abs(get_overlap(
														bands.get_eigensolution(kidx, ii * 2),
														bands.get_eigensolution(kidx, ii * 2)	
													 )) << "\n"
						 << "  <d[1] |  d[1]> = " << tdkp_math::abs(get_overlap(
														bands.get_eigensolution(kidx, ii * 2),
														bands.get_eigensolution(kidx, ii * 2)	
													 )) << "\n";
					TDKP_LOGMSG(LOG_INFO, sout.str());												 
				}
				 				
				double d0d1 = tdkp_math::abs(get_overlap(
									bands.get_eigensolution(kidx, ii * 2),
									bands.get_eigensolution(kidx, ii * 2 + 1)	
							  ));
				double d0d0 = tdkp_math::abs(get_overlap(
									bands.get_eigensolution(kidx, ii * 2),
									bands.get_eigensolution(kidx, ii * 2)	
							  ));
				double d1d1 = tdkp_math::abs(get_overlap(
									bands.get_eigensolution(kidx, ii * 2 + 1),
									bands.get_eigensolution(kidx, ii * 2 + 1)	
							  ));
				
				if(kidx == 0) {
					if(d0d1 > orthogonality_threshold || d0d0 - 1 > orthogonality_threshold || d1d1 - 1 > orthogonality_threshold) {
						TDKP_LOGMSG(LOG_WARN, "decrossing algorithm failed at kidx 0! you should probably start at k = 1.0e-4 to circumvent degeneracy at gamma point.");
					}
				} else {			  						  			
					TDKP_ASSERT(d0d1 < orthogonality_threshold, "|<w_ki_s1 | w_ki_s2>| (" << d0d1 << ") < 1.0e-9 (combo " << ii << " at kidx " << kidx << " vs. kref " 	<< kref << "). check bands and increase state_resorting_minimum_orthogonality_threshold.");
					TDKP_ASSERT(d0d0 - 1 < orthogonality_threshold, "|<w_ki_s1 | w_ki_s1>| (" << d0d0 << ")  - 1 < 1.0e-9 (combo " << ii << " at kidx " << kidx << " vs. kref " << kref << "). check bands and increase state_resorting_minimum_orthogonality_threshold.");
					TDKP_ASSERT(d1d1 - 1 < orthogonality_threshold, "|<w_ki_s2 | w_ki_s2>| (" << d1d1 << ")  - 1 < 1.0e-9 (combo " << ii << " at kidx " << kidx << " vs. kref " << kref << "). check bands and increase state_resorting_minimum_orthogonality_threshold.");
				}
			} 																 											 			 
			it++;	
		}
	}
	
	
	// ---------------------------------------------------
	// save decrossing info
	// ---------------------------------------------------	
	if(Configuration::get_instance()->get("output_decrossing_determinates") == 1.0) {
		static int decross_counter = 0;
		const int width = 15;
		ostringstream fname; 
		fname << "decrossing_" << decross_counter++ << ".dat";
		ofstream fout(fname.str().c_str());
		fout << setw(width) << "# kidx";
		for(int ii = 0; ii < bands.get_number_of_bands() / 2; ii++) {
			fout << setw(width) << ii << " ";
		}
		for(int kk = 0; kk < bands.get_number_of_k_values(); kk++) {
			fout << setw(width) << kk;
			for(int ii = 0; ii < bands.get_number_of_bands() / 2; ii++) {
				fout << setw(width) << decrossing_info[ii * bands.get_number_of_k_values() + kk];	
			}
			fout << "\n";
		}
		fout.close();
	}		
}		


}
