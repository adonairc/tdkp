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


#include "tdkp/probdefs/EffectiveMass.h"
#include "tdkp/main/FEMSolverGEVP.h"

namespace tdkp {

template<class PC>
SchroedingerPML<PC>::SchroedingerPML(
	PC& problem
) : EigenProblem<cplx,cplx,cplx>(
		problem.get_geometry(),
		problem.get_material_database()		
	),
	base_problem(problem),
	default_pml_alpha(Configuration::get_instance()->get("pml_default_alpha_real_stretch")),
	default_pml_beta(Configuration::get_instance()->get("pml_default_beta_imag_stretch")),
	default_pml_power(static_cast<unsigned int>(Configuration::get_instance()->get("pml_default_stretch_power")))
{
	this->num_equations_per_node = base_problem.get_num_equations_per_node();
	// ----------------------------------------------
	// init stretching to zero
	// ----------------------------------------------
	stretch[0].assign(get_geometry().get_num_regions(), 0);
	stretch[1].assign(get_geometry().get_num_regions(), 0);
	stretch[2].assign(get_geometry().get_num_regions(), 0);	
	pml_regions.assign(get_geometry().get_num_regions(), false); 
}

template<class PC>
SchroedingerPML3D<PC>::SchroedingerPML3D(
	PC& problem
) : SchroedingerPML<PC>(problem)
{
}

template<class PC>
SchroedingerPML1D2D<PC>::SchroedingerPML1D2D(
	PC& problem
) : SchroedingerPML<PC>(problem)
{
}



/** pass solution to base class */
template<class PC>
void SchroedingerPML<PC>::add_solution(cplx solution_value, const cplx* solution_vector, int length) {
	
	// ------------------------------------------------
	// apply back shift of energies 
	// ------------------------------------------------			
	base_problem.add_solution(solution_value, solution_vector, length); 	
		
}

/** pass sparsity pattern from base class */
template<class PC>
const int* SchroedingerPML<PC>::get_node_sparsity_pattern(int &num) const {
	return base_problem.get_node_sparsity_pattern(num);
}

/** effective mass is slightly different than usual kp calculations */   	
template<>
void SchroedingerPML<EffectiveMass>::calculate_element_matrices(const Element* elem, cplx* lhs, cplx* rhs, int* node_internal_indices, int &n) const;

template<class PC>
void SchroedingerPML<PC>::calculate_element_matrices(const Element* elem, cplx* lhs, cplx* rhs, int* node_internal_indices, int &n) const {
		
	vector<double> crhs(num_equations_per_node*num_equations_per_node*Element::max_num_nodes*Element::max_num_nodes);
	vector<cplx>   clhs(num_equations_per_node*num_equations_per_node*Element::max_num_nodes*Element::max_num_nodes);
	// -----------------------------------------
	// use pml wrapper element if we are in pml region
	// ----------------------------------------- 
	if(elem->enabled() && pml_regions[elem->get_region().get_index_global()]) {
		// -----------------------------------------
		// init wrapper with my element
		// -----------------------------------------
		TDKP_ASSERT(pml_element_wrappers.size() > elem->get_element_unique_type_key(), "");
		TDKP_ASSERT(pml_element_wrappers[elem->get_element_unique_type_key()] != 0, "");		
		ElementPML copy(*pml_element_wrappers[elem->get_element_unique_type_key()], elem);
		for(unsigned int dd = 0; dd < get_geometry().get_dimension(); dd++) {			
			copy.set_pml_stretch(dd,stretch[dd][elem->get_region().get_index_global()]);
		}
		copy.prepare();
		copy.return_real_part();
		base_problem.calculate_element_matrices(&copy, &clhs[0], &crhs[0], node_internal_indices, n);
		// copy real values from decorated class
		for(unsigned int ii = 0; ii < crhs.size(); ii++) {
			rhs[ii] = crhs[ii];
			lhs[ii] = clhs[ii];				
		}
		copy.return_imaginary_part();
		base_problem.calculate_element_matrices(&copy, &clhs[0], &crhs[0], node_internal_indices, n);
		// add imaginary parts
		cplx i(0.0, 1.0);
		for(unsigned int ii = 0; ii < crhs.size(); ii++) {
			rhs[ii] += i * crhs[ii];
			lhs[ii] += i * clhs[ii];				
		}							
	} else {
		// -----------------------------------------
		// no pml, use standard method
		// -----------------------------------------
		base_problem.calculate_element_matrices(elem, &clhs[0], &crhs[0], node_internal_indices, n);
		// copy values from decorated class
		for(unsigned int ii = 0; ii < crhs.size(); ii++) {
			rhs[ii] = crhs[ii];
			lhs[ii] = clhs[ii];				
		}		
	}			
}

template<class PC>
void SchroedingerPML<PC>::prepare() {
	if(!ready) {
		build_coordinate_stretches();
		build_element_wrappers();
	}
	base_problem.prepare();	
	this->ready = true;
}


template<class PC>
void SchroedingerPML3D<PC>::solve(int num_solutions) {
	
	// ----------------------------------------------
	// in any case, we have to prepare our base problem
	// ----------------------------------------------
	this->base_problem.prepare();
	this->base_problem.delete_solutions();
	
	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");	
	}
				
	// ----------------------------------------------
	// ensure we get nonsymmetric matrices
	// ----------------------------------------------
	double tmp_sym = Configuration::get_instance()->get("assembly_build_nonsymmetric_matrices");
	if(tmp_sym != 1.0) {
		TDKP_LOGMSG(LOG_WARN, "SchroedingerPML3D: enforcing use of nonsymmetric matrices");	
	}	
	Configuration::get_instance()->set("assembly_build_nonsymmetric_matrices", 1.0); 		
	FEMSolverGEVP<complex<double>, complex<double>, complex<double> > solver(this->geometry, *this);
	Configuration::get_instance()->set("assembly_build_nonsymmetric_matrices", tmp_sym); 	

	solver.create_matrix_structures();
	
	// ---------------------------------------------
	// PMLS NEED NONSYMMETRIC MATRICES
	// ---------------------------------------------
	TDKP_ASSERT(solver.get_lhs_matrix().property_is_set(nonsymmetric_matrix), "");
	TDKP_ASSERT(solver.get_rhs_matrix().property_is_set(nonsymmetric_matrix), "");	
	
	solver.assemble_system();
	
	if(this->base_problem.get_solution_type() == electrons) {
		solver.set_ordering(ascending);	
	} else {
		solver.set_ordering(descending);			
	}
	
	solver.solve_system(num_solutions);	
	this->base_problem.update_bandstructure_container();
	
}
	
template<class PC>
void SchroedingerPML1D2D<PC>::solve(int num_solutions) {
	TDKP_GENERAL_EXCEPTION("SchroedingerPML1D2D: sorry, but you need to define the kspace for 1D2D kp problems! therefore this function is disabled");
}
template<class PC>
void SchroedingerPML1D2D<PC>::solve(int num_subbands, double kmin, double kmax, int num_k_values) {
	this->solve(num_subbands, this->base_problem.build_kspace_domain(num_subbands, kmin, kmax, num_k_values));
}	

template<class PC>
void SchroedingerPML1D2D<PC>::solve(int num_subbands, const DomainMaster& domain) {
	
	if(this->geometry.get_num_nonzero_nodes() <= 0) {
		TDKP_GENERAL_EXCEPTION("no interior matrices available");
	}
	ostringstream sout;
	this->base_problem.prepare();
    this->base_problem.solving_preinform_user(num_subbands, domain);
      		
	// -------------------------------------------------------------
	// prepare bandstructure object
	// -------------------------------------------------------------
	BandstructureDomain<cplx>* band = new BandstructureDomain<cplx>(this->base_problem.get_num_equations_per_node(), num_subbands, this->base_problem.geometry.get_num_nodes(), domain);
  	this->base_problem.bandstructures.push_back(band);
	
	// -------------------------------------------------------------
	// reset energy of last solution to the place where we expect the 
	// new solutions (bandedges is not a good idea when potential 
	// energy fields are present
	// -------------------------------------------------------------
	this->base_problem.k_idx_current = 0; // nasty, must be 0 before we get the right energy shift (in subsequent calc., the k idx may still be set > 0)
	if(this->base_problem.use_user_defined_energy_guess()) {
		TDKP_ASSERT(this->base_problem.energy_guess.size() > 0, "this->base_problem.energy_guess.size() > 0");
		this->base_problem.energy_of_last_solution = this->base_problem.energy_guess[0];
	} else {
		this->base_problem.energy_of_last_solution = this->base_problem.get_minimum_bandedges();
	}  

	// -------------------------------------------------------------
	// create fem assembler (parallel or non-parallel)
	// -------------------------------------------------------------
	// standard solver object
	double tmp_sym = Configuration::get_instance()->get("assembly_build_nonsymmetric_matrices");
	if(tmp_sym != 1.0) {
		TDKP_LOGMSG(LOG_WARN, "SchroedingerPML1D2D: enforcing use of nonsymmetric matrices");	
	}	
	Configuration::get_instance()->set("assembly_build_nonsymmetric_matrices", 1.0);	
	FEMSolverGEVP<complex<double>, complex<double>, complex<double> >* solver = 0;
	solver = new FEMSolverGEVP<complex<double>, complex<double>, complex<double> >(this->base_problem.geometry, *this);
	Configuration::get_instance()->set("assembly_build_nonsymmetric_matrices", tmp_sym);
		
	// -------------------------------------------------------------
	// create matrix structure
	// -------------------------------------------------------------
	solver->create_matrix_structures();
	
	// ---------------------------------------------
	// PMLS NEED NONSYMMETRIC MATRICES
	// ---------------------------------------------
	TDKP_ASSERT(solver->get_lhs_matrix().property_is_set(nonsymmetric_matrix), "");
	TDKP_ASSERT(solver->get_rhs_matrix().property_is_set(nonsymmetric_matrix), "");		
		
	if(this->base_problem.get_solution_type() == electrons) {
		solver->set_ordering(ascending);
	} else {
		solver->set_ordering(descending);	
	}
	if(Logger::get_instance()->get_level() == LOG_INFO) {
		Logger::get_instance()->init_progress_bar("calculating bandstructure in k space, using points: ", domain.get_number_of_points());
	}	
		
	this->prepare();
	
	// -------------------------------------------------------------
	// tell the matrices to shut up
	// -------------------------------------------------------------
	for(unsigned int ii = 0; ii < this->base_problem.kp_matrices.size(); ii++) {	
		this->base_problem.kp_matrices[ii]->set_output_surpression(true);		
	}
		
	int next    = 0;
	int current = 0;
	// -------------------------------------------------------------
	// loop over domain points
	// -------------------------------------------------------------
  	for(int kk = 0; kk < (signed)domain.get_number_of_points(); kk++) {
		adaptive_omp_threading();  		  		
  		this->base_problem.k_idx_current = kk; // setting current index to store bandstructure at the right place
  		
  		// --------------------------------------------------------
  		// set the length value of the transversal direction
  		// --------------------------------------------------------
  		this->base_problem.k_transversal = domain.get_point(kk).get_coord_abs();
  		
		// --------------------------------------------------------
		// if we have a non-radial well, rotate the kp matrices for
		// the next k point
		// --------------------------------------------------------
		if(this->geometry.get_dimension() == 1 && !domain.radial()) {
			// ----------------------------------------------
			// kz is the free direction in the assembly
			// so, kz' = kz * cos(phi) - ky * sin(phi)
			// and this equals
			//         = kz * dir_x - ky * dir_y
			// where dir is the unit direction of the point 			
			// ----------------------------------------------
			double dir_x = domain.get_point(kk).get_coord(0) / this->base_problem.k_transversal;
			double dir_y = domain.get_point(kk).get_coord(1) / this->base_problem.k_transversal;
			Vector3D new_kz = dir_x * this->base_problem.k_directions[D_DZ] - dir_y * this->base_problem.k_directions[D_DY];
			Vector3D new_ky = dir_x * this->base_problem.k_directions[D_DY] + dir_y * this->base_problem.k_directions[D_DZ];
			RMatrix<double> tmp_rotation(3,3);
			tmp_rotation.set_row(0, this->base_problem.k_directions[D_DX].get_all());
			tmp_rotation.set_row(1, new_ky.get_all());
			tmp_rotation.set_row(2, new_kz.get_all());
			// ----------------------------------------------
			// set rotation to kp matrices
			// ----------------------------------------------
			for(unsigned int ii = 0; ii < this->base_problem.kp_matrices.size(); ii++) {	
				this->base_problem.kp_matrices[ii]->set_rotation(tmp_rotation);				
			}
  			// --------------------------------------------------------
			// information for user
			// --------------------------------------------------------			
			sout.str("");
			sout << "iteration " << kk + 1 << " for kx,ky = (" << domain.get_point(kk).get_coord(0)
			     << ", " << domain.get_point(kk).get_coord(1) <<") "
		     	 << " [1/nm]. expecting kp states at " << - this->base_problem.get_energy_shift() << " [eV] ";
			if(this->base_problem.use_user_defined_energy_guess()) {
				sout << "(user defined) ";	
			}		     
			sout << "for the next cycle.";
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());									   	
		} else {
  			// --------------------------------------------------------
			// information for user
			// --------------------------------------------------------			
			sout.str("");
			sout << "iteration " << kk + 1 << " of " << domain.get_number_of_points() << " for k = " << this->base_problem.k_transversal 
		     	 << " [1/nm]. expecting kp states at " << - this->base_problem.get_energy_shift() << " [eV] ";
			if(this->base_problem.use_user_defined_energy_guess()) {
				sout << "(user defined) ";	
			}		     
			sout << "for the next cycle.";
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());				
		}
		
		// --------------------------------------------------------		
    	// prepare kp matrices for the next iteration 
    	// (give them the new energy shift)
    	// --------------------------------------------------------		    	
		for(unsigned int ii = 0; ii < this->base_problem.kp_matrices.size(); ii++) {						
			this->base_problem.kp_matrices[ii]->set_energy_shift(this->base_problem.get_energy_shift());		
			this->base_problem.kp_matrices[ii]->calculate();
			TDKP_ASSERT(this->base_problem.kp_matrices[ii]->surpress_output(), "somehow the kp matrices where deleted ... this is an inconsistency in prepare and it should happen ... ");			
		}  		
  		
		// --------------------------------------------------------		  		
  		// assemble and solve system  		
  		// --------------------------------------------------------		
    	solver->assemble_system();    	    	    	      	
    	solver->solve_system(num_subbands);
	    	
    	// --------------------------------------------------------		
    	// add solutions from cache to bs object    	
    	// --------------------------------------------------------		
	    this->base_problem.add_cache_to_bandstructure_object(band, kk);
	    			  	
    	if(current == next) {
    		if(Logger::get_instance()->get_level() == LOG_INFO) {
    			next = Logger::get_instance()->set_progress_bar(current, domain.get_number_of_points());
    		}
    	}
    	current++;
   		solver->reset_matrix_to_zero();   		   		
  	}
  	  	  	
  	if(Logger::get_instance()->get_level() == LOG_INFO) {
  		Logger::get_instance()->end_progress_bar();
  	}
  
  	
  	// -------------------------------------------------------------
	// so, the matrices are allowed to talk again ...
	// -------------------------------------------------------------
	for(unsigned int ii = 0; ii < this->base_problem.kp_matrices.size(); ii++) {	
		this->base_problem.kp_matrices[ii]->set_output_surpression(false);
	}
	
	// delete solver object
	delete solver;
	
}

/** pass unique identifer from base class */
template<class PC>
string SchroedingerPML<PC>::get_unique_identifier() const {
	return base_problem.get_unique_identifier();	
}

template<class PC>
EigenProblemType SchroedingerPML<PC>::get_problem_type() const {
	return base_problem.get_problem_type();	
}

template<class PC>
void SchroedingerPML<PC>::display_solution_info() const {
	
	const BandstructureDomain<complex<double> >& ref = base_problem.get_bandstructure();
	// -------------------------------------------------
	// sort imaginary values
	// -------------------------------------------------
	vector<double> imag_values;
	vector<int> indices;
	for(int ii =  0; ii < ref.get_number_of_bands(); ii++) {
		imag_values.push_back(tdkp_math::abs(ref.get_energy(0,ii).imag()));
		indices.push_back(ii);
	}
	tdkp_math::tracked_sort<int,double>(
		indices.begin(),
		indices.end(),
		imag_values.begin(),
		imag_values.end(),
		true
	);
	// -------------------------------------------------
	// output 
	// -------------------------------------------------
	ostringstream sout;
	sout << "SchroedingerPML: we have the following states (in order of stability):\n";
	for(int ii =  0; ii < ref.get_number_of_bands(); ii++) {
		sout << "   state " << indices[ii] 
		     << " at " << ref.get_energy(0,indices[ii]).real() << " eV, tau = "
		     << (constants::hbar / imag_values[ii] / 2.0 * 1e12) << " ps\n";
	}
	TDKP_LOGMSG(LOG_INFO, sout.str());   	
}

template<class PC>
SchroedingerPML<PC>::~SchroedingerPML() {
	this->delete_stretches();
	this->delete_pml_element_wrappers();	
}

template<class PC>
void SchroedingerPML<PC>::set_stretch_parameters(const double& alpha, const double& beta, unsigned int pml_power) {
	default_pml_alpha = alpha;
	default_pml_beta  = beta;
	default_pml_power = pml_power;
	this->delete_stretches();
	this->delete_pml_element_wrappers();	
}
	
	
/** build stretches for all regions names PML */
template<class PC>
void SchroedingerPML<PC>::build_coordinate_stretches() {
	
	this->delete_stretches();
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "SchroedingerPML: checking regions for possible PMLs");
	const string cmp_with("PML_");
	// ------------------------------------------------
	// for all regions
	// ------------------------------------------------	
	for(int rr = 0; rr < (signed)get_geometry().get_num_regions(); rr++) {
		// ------------------------------------------------
		// if region name starts with PML_
		// ------------------------------------------------
		string region_name = get_geometry().get_region(rr).get_name();	
		if(region_name.substr(0,4) == cmp_with) {									
			// ------------------------------------
			// find pml type string
			// ------------------------------------
			int pml_type[] = {0,0,0};
			string pml_type_string = region_name.substr(4);			
			unsigned int ww = 0; 
			while(pml_type_string[ww] != '_' && ww < pml_type_string.size()) {
				ww++;	
			}
			if(ww < pml_type_string.size()) {
				pml_type_string = pml_type_string.substr(ww+1);
			}
			// ------------------------------------
			// determine pml directions
			// ------------------------------------			
			TDKP_ASSERT(pml_type_string.size() % 2 == 0, pml_type_string << " is not a valid PML specifier!");
			unsigned int gg = 0;
			while(gg < pml_type_string.size() / 2) {
				string sub = pml_type_string.substr(gg*2,2);
				int idx = 0;
				int sig = 0;
				switch(sub[0]) {
					case 'p': sig =  1; break;
					case 'm': sig = -1; break;
					default:
						TDKP_GENERAL_EXCEPTION(pml_type_string << " is not a valid PML specifier! could not determine p/m");
				}				
				switch(sub[1]) {
					case 'x': idx = 0; break;
					case 'y': idx = 1; break;
					case 'z': idx = 2; break;
					default:
						TDKP_GENERAL_EXCEPTION(pml_type_string << " is not a valid PML specifier! could not determine x/y/z");					
				}
				pml_type[idx] = sig;
				gg++;						 	
			}  												
			double min_coords[] = {0.0e0, 0.0e0, 0.0e0};
			double max_coords[] = {0.0e0, 0.0e0, 0.0e0};
			bool   first = true;
			// -----------------------------------------------
			// for every element in geometry
			// -----------------------------------------------
			for(Geometry::element_const_iterator eit = get_geometry().elements_begin(); eit != get_geometry().elements_end(); eit++) {
				// -------------------------------------------
				// if element belongs to region
				// -------------------------------------------
				if((*eit)->get_region().get_index_global() == rr) {
					// -----------------------------------------------
					// find minmax coordinates 
					// -----------------------------------------------
					if(first) {					
						for(unsigned int ii = 0; ii < get_geometry().get_dimension(); ii++) {
							max_coords[ii] = min_coords[ii] = (*eit)->get_node(0).get_coord(ii);
						}						
						first = false;	
					} else {
						// for every node
						for(unsigned int vv = 0; vv < (*eit)->get_num_nodes(); vv++) {
							for(unsigned int ii = 0; ii < get_geometry().get_dimension(); ii++) {
								const double& coord = (*eit)->get_node(vv).get_coord(ii);
								if(max_coords[ii] < coord) {
									max_coords[ii] = coord;	
								}
								if(min_coords[ii] > coord) {
									min_coords[ii] = coord;	
								}								
							}
						}	
					}					 
				}
			}
			// --------------------------------------------
			// build pml functions
			// --------------------------------------------
			for(unsigned int ii = 0; ii < get_geometry().get_dimension(); ii++) {
				if(pml_type[ii] != 0) {
					TDKP_ASSERT((signed)stretch[ii].size() > rr && stretch[ii][rr] == 0, ""); // should be empty
					double start, end;
					if(pml_type[ii] > 0) { 
						start = min_coords[ii];
					    end   = max_coords[ii];
					} else {
						start = max_coords[ii];
					    end   = min_coords[ii];					
					}
					stretch[ii][rr] = new PMLStretch(
						start, end, default_pml_alpha, default_pml_beta, default_pml_power
					);
					pml_regions[rr] = true;										
					TDKP_LOGMSG(LOG_INFO_DEVEL1, "SchroedingerPML: " << get_geometry().get_region(rr).get_name() << " is pml in dir " << ii << " ranging from " << start << " to " << end);
				}
			}						
		}
	}
	
}


/** delete coordinate stretches */
template<class PC>
void SchroedingerPML<PC>::delete_stretches() {	
	// -------------------------------------
	// delete coordinate stretches
	// -------------------------------------
	for(unsigned int ii = 0; ii < 3; ii++) {		
		for(unsigned int jj = 0; jj < stretch[ii].size(); jj++) {
			if(stretch[ii][jj] != 0) {
				delete stretch[ii][jj];
				stretch[ii][jj] = 0;	
			}		
		}
	}
	this->ready = false;	
}

/** build element pml wrapper objects (which will do the numerical integration) */ 
template<class PC>
void SchroedingerPML<PC>::build_element_wrappers() {
	Geometry::element_const_iterator eit = get_geometry().elements_begin(); 
	while(eit != get_geometry().elements_end()) {
		
		// --------------------------------------
		// do only for enabled elements
		// --------------------------------------
		if((*eit)->enabled()) {			
			// --------------------------------------
			// resize pml wrappers array if needed
			// --------------------------------------
			if((*eit)->get_element_unique_type_key() >= pml_element_wrappers.size()) {
				pml_element_wrappers.resize((*eit)->get_element_unique_type_key() + 1);		
			}
			// --------------------------------------
			// create new wrapper if necessary
			// --------------------------------------
			if(pml_element_wrappers[(*eit)->get_element_unique_type_key()] == 0) {
				TDKP_LOGMSG(LOG_INFO_DEVEL2, "SchroedingerPML: building wrapper element of type " << (*eit)->get_element_unique_type_key());
				pml_element_wrappers[(*eit)->get_element_unique_type_key()] = new ElementPML((*eit));			
			}
		}
		eit++;
	}
}

template<class PC>
void SchroedingerPML<PC>::delete_pml_element_wrappers() {
	for(unsigned int ii = 0; ii < pml_element_wrappers.size(); ii++) {
		if(pml_element_wrappers[ii] != 0) {
			delete pml_element_wrappers[ii];
			pml_element_wrappers[ii] = 0;	
		}	
	}
	pml_element_wrappers.clear();
}

template<class PC>
void SchroedingerPML<PC>::dump_pml_functions(const char* filename) const {
	
	// -----------------------------------------
	// for every region
	// -----------------------------------------
	for(unsigned int rr = 0; rr < get_geometry().get_num_regions(); rr++) {
		if(pml_regions[rr]) {		
			// ----------------------------------------------
			// open outfile
			// ----------------------------------------------			
			ostringstream fname_dir;
			fname_dir << filename << "_" << get_geometry().get_region(rr).get_name() << ".dat";
			ofstream fout(fname_dir.str().c_str());
			if(!fout) {
				TDKP_GENERAL_EXCEPTION("could not write to file " << fname_dir);
			}
			fout << "# x,[y,z],pml_x r + i[,pml_y r + i,pml_z r + i]\n";
			// ----------------------------------------------
			// for every element
			// ----------------------------------------------
			Geometry::element_const_iterator it = get_geometry().elements_begin();
			while(it != get_geometry().elements_end()) {
				if((*it)->enabled() && (*it)->get_region().get_index_global() == (signed)rr) {
					// mid point
					double cc[] = {0.0e0, 0.0e0, 0.0e0};					
					for(unsigned int nn = 0; nn < (*it)->get_num_nodes(); nn++) {
						const Node& nd = (*it)->get_node(nn);
						for(unsigned int hh = 0; hh < get_geometry().get_dimension(); hh++) {
							cc[hh] += nd.get_coord(hh) / static_cast<double>((*it)->get_num_nodes());
						}		
					}								
					for(unsigned int hh = 0; hh < get_geometry().get_dimension(); hh++) {
						fout << cc[hh] << "  ";	
					}
					fout << "    ";
					// evaluate pml at mid point
					for(unsigned int hh = 0; hh < get_geometry().get_dimension(); hh++) {
						if(stretch[hh][rr] != 0) {
							cplx tmp = stretch[hh][rr]->evaluate(cc[hh]);
							fout << tmp.real() << "  " << tmp.imag() << "           ";  		
						}										
					}
					fout << "\n";
				}
				it++;
			}
			fout.close();							
		}
	}
}


} // end of namespace
