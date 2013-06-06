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

#include "tdkp/clc/CLCSHF.h"
#include "tdkp/solvers/GMRES.h"
#include "tdkp/clc/CLCScreening.h"

namespace tdkp {

/** constructor taking the ingredients for a shf calculation */
CLCSHF1DR::CLCSHF1DR(
	const CLCIngredients1DR& ingredients_
) 
: ingredients(ingredients_),
  delta_k(0.004),
  omega_min(0.0),
  omega_max(0.0),
  omega_num(100)
{
	this->create_dense_domain();
}

CLCSHF1DR::~CLCSHF1DR() { 
	
}

/** set new omega range */
void CLCSHF1DR::set_omega_range(const double& omega_min_, const double& omega_max_, unsigned int omega_num_) {
	TDKP_ASSERT(omega_min_ < omega_max_, "omega_min < omega_max");	
	omega_min = omega_min_;
	omega_max = omega_max_;
	omega_num = omega_num_;	
}
/** set new number of frequency points */
void CLCSHF1DR::set_omega_num(unsigned int omega_num_) {
	omega_num = omega_num_;	
}
void CLCSHF1DR::set_interpolation_density(const double& delta_k_) {
	delta_k = delta_k_;	
	this->create_dense_domain();	
}


void CLCSHF1DR::create_dense_domain() {
	
	// ----------------------------------------------------
	// determine kmax of our integration
	// ----------------------------------------------------
	double radial_kmax = ingredients.get_cb_bands().get_domain().get_last_point().get_coord_abs();
	radial_kmax = min(radial_kmax, ingredients.get_vb_bands().get_domain().get_last_point().get_coord_abs());
	radial_kmax = min(radial_kmax, ingredients.get_matrix_elements().get_domain().get_last_point().get_coord_abs());

	TDKP_ASSERT(radial_kmax > 0.0, "radial_kmax > 0.0"); 

	// ----------------------------------------------------
	// determine number of k points
	// ----------------------------------------------------
	unsigned int num_k_points = static_cast<unsigned int>(ceil(radial_kmax / delta_k));  			
 
 	// ----------------------------------------------------
 	// clean old 
 	// ----------------------------------------------------
 	dense_domain.clean();
 
 	ostringstream sout;
 	sout << "CLC: creating dense radial domain with " << num_k_points 
 	     << " points within [0, " << radial_kmax << "]";
 	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
 
	
	switch(ingredients.get_cb_bands().get_domain().get_dimension()) {
		// --------------------------------------------------
		// quantum dot: exception ...
		// --------------------------------------------------
		case 0:
			TDKP_GENERAL_EXCEPTION("CLCSHF1DR does not work for dots");
			break;
		// --------------------------------------------------
		// quantum wire: create dense trapezoidal line
		// --------------------------------------------------
		case 1:
			create_1D_domain_wire_bands(dense_domain, 0.0, radial_kmax, num_k_points);
			break; 
		// --------------------------------------------------
		// quantum well: create dense radial integration line
		// --------------------------------------------------
		case 2:
			create_2D_domain_radial(dense_domain, 0.0, radial_kmax, num_k_points);
			break;		
		// --------------------------------------------------
		// bulk: create dense 3D radial integration line
		// --------------------------------------------------
		case 3:
			create_3D_domain_radial(dense_domain, Vector3D(1.0, 0.0, 0.0), 0.0, radial_kmax, num_k_points);		
			break;
		default:
			TDKP_GENERAL_EXCEPTION("invalid dimension ... ");
	}
		
}

/** calculate for all transitions that are calculated within the matrix elements */
CLCOpticalResults* CLCSHF1DR::calculate( 
	const double& cb_fermi, 
	const double& vb_fermi, 
	const double& temperature
) const {
	return calculate(
		ingredients.get_matrix_elements().get_num_cb_bands(), 
		ingredients.get_matrix_elements().get_num_vb_bands(),
		cb_fermi,
		vb_fermi,
		temperature
	);
}
		 
/** calculate gain / luminescence for the given number of transitions */
CLCOpticalResults* CLCSHF1DR::calculate(
	unsigned int num_cb_bands, 
	unsigned int num_vb_bands, 
	const double& cb_fermi, 
	const double& vb_fermi, 
	const double& temperature
) const {

	TimeMeasurements::get_instance().start("shf radial evaluation");

	ostringstream sout;
	sout << "CLCSHF 1D-Radial evaluates optics for:\n"
	     << "  " << num_cb_bands << " cb (sub)bands and "<< num_vb_bands << " vb (sub)bands,\n"
	     << "  omega range: " << omega_min * constants::hbar << " -> " 
	     << omega_max * constants::hbar << " [eV],\n"
	     << "  fermi levels of " << cb_fermi << " (cb) and " << vb_fermi << " (vb) at  "
	     << temperature << " Kelvin.\n"
	     << "  Homogeneous broadening is " << ingredients.get_homogeneous_broadening() << ",\n"
	     << "  host optical refractive index is " << ingredients.get_optical_refractive_index() << ",\n"
	     << "  static dielectric screening constants is " << ingredients.get_static_dielectric_constant() << ",\n"
	     << "  quantized volume is " << ingredients.get_geometry_output_division_factor() << ",\n"
	     << "  k-space dimension = " << dense_domain.get_dimension();
	      
	TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());	     
	TDKP_ASSERT(omega_min > 0.0, "omega_min is <= 0.0!");
	TDKP_ASSERT(omega_max > omega_min, "omega_max is smaller than omega_min!");
	
	if(tdkp_math::abs(2.0 * ingredients.get_vb_bands().get_domain().get_last_point().get_coord_abs() 
	   - ingredients.get_coulomb_matrix_elements().get_q_values().back()) > 1.0e-5) {
	   		
		ostringstream sout;
		sout << "CLCSHF: sorry, but you have to ensure that the q-space for the "
		     << "coulomb integrals extends exactly 2*kmax. not less, not more. "
		     << "your kmax = " << ingredients.get_vb_bands().get_domain().get_last_point().get_coord_abs() << ", "
		     << "your qmax = " << ingredients.get_coulomb_matrix_elements().get_q_values().back() << ", "
		     << "|2kmax - qmax| = " << tdkp_math::abs(2.0 * ingredients.get_vb_bands().get_domain().get_last_point().get_coord_abs() - ingredients.get_coulomb_matrix_elements().get_q_values().back());
		TDKP_GENERAL_EXCEPTION(sout.str());		     
	} 		      
		      
	// ------------------------------------------
	// define i
	// ------------------------------------------
	const cplx i(0.0, 1.0);

	// ------------------------------------------
	// which polarizations?
	// ------------------------------------------
	const unsigned int pp_start = static_cast<unsigned int>(Configuration::get_instance()->get("clc_polarization_start"));
	const unsigned int pp_end   = static_cast<unsigned int>(Configuration::get_instance()->get("clc_polarization_end"));
	TDKP_ASSERT(pp_start >= 0, "clc_polarization_start >= 0 failed");
	TDKP_ASSERT(pp_start <= pp_end, "clc_polarization_start <= clc_polarization_end failed");
	TDKP_ASSERT(pp_end <= 2, "clc_polarization_end <= 2 failed ");
	
	// ------------------------------------------
	// check if we have some effmass bands
	// if we have only cb effmass, double num cb bands (to include spin)
	// if we have both cb/vb effmass, then be careful, the
	// summation over spin is already included in the matrix
	// element! (see eg. MomentumOperatorBulkEffectiveMass::eval)
	// ------------------------------------------
	bool cb_effmass = ingredients.get_cb_bands().get_basis_size() == 1 ? true:false;
	bool vb_effmass = ingredients.get_vb_bands().get_basis_size() == 1 ? true:false;
	
	TDKP_ASSERT(!(!cb_effmass && vb_effmass), "!(!cb_effmass && vb_effmass) is not expected as it wouldn't really make sense ..."); 
	
	if(cb_effmass && !vb_effmass) {
		num_cb_bands *= 2;	
	}
	
	// ------------------------------------------
	// set dimension independent prefactor
	// for the suszeptibility
	// ------------------------------------------	
	double prefactor = 1.0 / (ingredients.get_optical_refractive_index() 
	                        * ingredients.get_optical_refractive_index()
	                        * constants::vacuum_permittivity);
	                        
	double spont_prefactor = ingredients.get_optical_refractive_index()
	                       / (
	                       		constants::vacuum_permittivity
	                       	  * constants::hbar
	                       	  * constants::c
	                       	  * constants::c
	                       	  * constants::c
	                       	  * constants::pi
	                       	  * constants::pi
	                         );	                        
	                        
	// -------------------------------------------------------
	// build dimension dependent prefactor
	// -------------------------------------------------------
	double dim_prefactor = 0.0;
	// -------------------------------------------------------
	// the matrix weight factor is used to fix the weight of the radial 
	// domain object. because there the integral over the angles theta and phi
	// is already performed (which gives 2pi (well) or 4pi (bulk))
	// in the coulomb matrices, we actually perform the angular integration
	// numerically and therefore automatically have the angular term in there
	// 
	// further, we also add the dimension dependent terms resulting 
	// form the first change from sum to integral
	// -------------------------------------------------------
	double matrix_weight_factor = 0.0;
	switch(dense_domain.get_dimension()) {
		case 0:
			TDKP_GENERAL_EXCEPTION("the stuff is not suitable for the dot ...");
			break;
		// ----------------------------------------------------
		// -> the radial part is handled inside the DomainMaster, 
		//    so this is only the volume stuff for the sum to integral!
		// ----------------------------------------------------
		case 1:
			TDKP_ASSERT(dense_domain.radial(), "even the quantum wire domain should be radial (-k folded on k)");
			// wire: 
			dim_prefactor = 1.0 / (2.0 * constants::pi * ingredients.get_geometry_output_division_factor());
			// ---------------------------------------------------
			// note, this one is a little bit complicated:
			// the wire has symmetry around k = 0, so E(k) = E(-k),
			// therefore, the summation over k is restricted to k > 0
			// and the coulomb factors are folded togehter
			// Lambda_(ki - kj,ki,kj) + Lambda_(ki + kj,ki,-kj)
			// BUT the domain object has weights such that -k is
			// mapped to +k. therefore the given factor of 2 is
			// already in the k point weight.
			// this is o.k. for the final summation of the gain.
			// but when evaluating the screening exchange shift
			// and solving the self consistent equation for the 
			// polarizations, we have to get rid of that factor
			// as the factor of 2 is included in the N1D(ki,kj)
			// term.
			// the matrix weight factor is the correct place to
			// do this as it is multiplied to all relevant equations
			// ---------------------------------------------------
			// first 2pi -> sum to integral, last 2 is for removing folding of -k onto k
			matrix_weight_factor = 1.0 / (2.0 * constants::pi * 2.0);  
			break;
		case 2:
			// well:
			dim_prefactor = 1.0 / (4.0 * constants::pi * constants::pi * ingredients.get_geometry_output_division_factor());
			// first 4pi^2 -> from sum to integral and last 2pi is for removing angular integral from domain.get_point.get_weight
			matrix_weight_factor = 1.0 / (4.0 * constants::pi * constants::pi * 2.0 * constants::pi);			   
			break;
		case 3:
			// bulk:
			dim_prefactor = 1.0 / (8.0 * constants::pi * constants::pi * constants::pi * ingredients.get_geometry_output_division_factor());
			// first 8pi^3 -> from sum to integral and last 4pi is for removing angular integral from domain.get_point.get_weight
			matrix_weight_factor = 1.0 / (8.0 * constants::pi * constants::pi * constants::pi * 4.0 * constants::pi);			
			break;
		default:
			TDKP_GENERAL_EXCEPTION("invalid dimension");					
	}
	
	// -----------------------------------------
	// add ec^2 / (4.0 pi eb eps0) to matrix factor
	// -----------------------------------------
	matrix_weight_factor *= constants::ec * constants::ec / (4.0 * constants::pi  
	                        		  * ingredients.get_static_dielectric_constant()
	                        		  * constants::vacuum_permittivity);	       
	
	// -----------------------------------------
	// check for included effects
	// -----------------------------------------
	const bool include_self_consistent_polarization = Configuration::get_instance()->get("clc_solve_self_consistent_lambda_k") == 1.0 ? true : false;
	const bool use_pade_approximation = Configuration::get_instance()->get("clc_solve_self_consistent_lambda_k") == 2.0 ? true : false;
	if(!include_self_consistent_polarization && !use_pade_approximation) {
		TDKP_LOGMSG(LOG_WARN, "CLC: you disabled the self-consistent polarization lambda_k");	
	}
	if(use_pade_approximation) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLC: using pade approximation for lambda_k");		
	}
	
	const unsigned int output_omega_idx = Configuration::get_instance()->get("output_clc_lambda_self_consistent");	
	const unsigned int number_of_points = dense_domain.get_number_of_points();
	
	// ------------------------------------------
	// prepare output storage
	// ------------------------------------------
	vector<cplx>   suszeptibility(omega_num * 3, cplx(0.0,0.0));
	vector<double> spont_emission(omega_num * 3, 0.0);
	
	// ------------------------------------------
	// prepare temporary openmp storage 
	// ------------------------------------------	
	unsigned int num_max_threads = omp_get_max_threads();		
	vector<vector<cplx> >   omp_omega_cv_k(omp_get_max_threads()); 
	vector<vector<double> > omp_spont_stats_cv_k(omp_get_max_threads()); 	
	vector<vector<cplx> >   omp_lambda_cv_k(omp_get_max_threads()); 		 
	vector<vector<cplx> >   omp_lambda_cv_k_from_matrix(omp_get_max_threads()); 
	vector<RMatrix<cplx> >  omp_lambda_cv_k_matrix(omp_get_max_threads());
	vector<vector<cplx> >   omp_dipole_mu_k(omp_get_max_threads()); 
	for(unsigned int ii = 0; ii < num_max_threads; ii++) {
		omp_omega_cv_k[ii].resize(number_of_points); 
		omp_spont_stats_cv_k[ii].resize(number_of_points); 	
		omp_lambda_cv_k[ii].resize(number_of_points); 	  
		omp_lambda_cv_k_from_matrix[ii].resize(number_of_points); 
		omp_lambda_cv_k_matrix[ii] = RMatrix<cplx>(number_of_points,number_of_points);
		omp_dipole_mu_k[ii].resize(number_of_points); 		
	}		
		
	// ------------------------------------------
	// prepare temporary storage 
	// ------------------------------------------	
	vector<double> omega_trans(number_of_points);
	vector<cplx>   vb_cb_coulomb_scattering_matrix(number_of_points * number_of_points);
	vector<double> k_point_coords(number_of_points);
	vector<cplx>   momentum_matrix_element(number_of_points);
			             
	// ------------------------------------------
	// precompute k-space coords
	// ------------------------------------------
	for(unsigned int ii = 0; ii < number_of_points; ii++) {		
		k_point_coords[ii] = dense_domain.get_point(ii).get_coord_abs();
	}	             
		    	
	// ------------------------------------------
	// fermi statistics calculator
	// ------------------------------------------
	CLCFermiStats carrier_stats;		
	TimeMeasurements::get_instance().track_memory_usage();

	
	// --------------------------------------------
	// setup linear solver
	// --------------------------------------------	
	GMRES gmres; // quick and dirty adaption of netlibs solver
	unsigned int linear_equation_solver = Configuration::get_instance()->get("clc_linear_equation_solver");


	// -------------------------------------------------------
	// precompute fermi statistics and energies for every band
	// -------------------------------------------------------
	// don't do it twice for effmass bands
	int num_cb_loops = ((cb_effmass  && !vb_effmass) ? num_cb_bands / 2 : num_cb_bands);
	int num_vb_loops = num_vb_bands;
	vector<vector<double> > all_cb_fermi_stats(num_cb_loops);
	vector<vector<double> > all_vb_fermi_stats(num_vb_loops);			
	vector<vector<double> > all_cb_band(num_cb_loops);
	vector<vector<double> > all_vb_band(num_vb_loops);
	// -------------------------------------------------------				
	// cb and vb in same loop (to have bigger and more parallel blocks,
	// for better parallel efficiency
	// -------------------------------------------------------
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "CLC: interpolating bands and precomputing fermi stats");
	//#pragma omp parallel for	
	for(int ii = 0; ii < num_cb_loops + num_vb_loops; ii++) {
		// --------------------------------------
		// cb or vb -> get right arrays
		// --------------------------------------
		const Spline1D* eval_spline      = 0;
		vector<double>* eval_into        = 0;
		vector<double>* eval_fermi_stats = 0;
		KPSolutionType  carrier_type;		
		double          carrier_fermi;		
		if(ii < num_cb_loops) {						
			eval_spline      = &ingredients.get_cb_band(ii);
			eval_into        = &all_cb_band[ii];
			carrier_type     = electrons;
			carrier_fermi    = cb_fermi;
			eval_fermi_stats = &all_cb_fermi_stats[ii];
		} else {
			eval_spline      = &ingredients.get_vb_band(ii - num_cb_loops);
			eval_into        = &all_vb_band[ii - num_cb_loops];
			carrier_type     = holes;	
			carrier_fermi    = vb_fermi;
			eval_fermi_stats = &all_vb_fermi_stats[ii - num_cb_loops];		
		}		
		(*eval_into).resize(number_of_points);
		(*eval_fermi_stats).resize(number_of_points);
		for(int nn = 0; nn < (signed)number_of_points; nn++) {
			(*eval_into)[nn] = (*eval_spline)(k_point_coords[nn]); 
		}			
		carrier_stats.calculate(carrier_type, carrier_fermi, temperature, *eval_into, *eval_fermi_stats);		
	}			

	// ------------------------------------------
	// calculate screening
	// ------------------------------------------	
	CLCScreening screening(ingredients, dense_domain);
	GenericCurve<double,double>* screening_curve = 0; // that's what we pass later 
	if(Configuration::get_instance()->get("clc_use_screening") == 1.0) {
		screening.update(
			((cb_effmass  && !vb_effmass) ? num_cb_bands / 2 : num_cb_bands),
			num_vb_bands,
			cb_fermi, vb_fermi, temperature
		);
		screening_curve = &screening;
		if(Configuration::get_instance()->get("output_clc_screening_epsilon") == 1.0) {
			static int sc_cnt = 0;
			ostringstream fname; fname << "screening" << sc_cnt++ << ".dat";
			ofstream fout(fname.str().c_str());
			for(unsigned int ii = 0; ii < dense_domain.get_number_of_points(); ii++) {
				fout << dense_domain.get_point(ii).get_coord_abs() << "  " 
				     << screening(dense_domain.get_point(ii).get_coord_abs())
				     << "\n"; 	
			}
			fout.close();
		}					
	} else {
		TDKP_LOGMSG(LOG_INFO_DEVEL1, "CLC: screening has been disabled by user request");
		screening_curve = &screening;
	}		
				
	// ------------------------------------------
	// calculate exchange energy ESX
	// ------------------------------------------
	if(Configuration::get_instance()->get("clc_include_exchange_shifts") == 1.0) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLC: screening exchange shifts (at k = 0) are: ");
	} else {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLC: screening exchange shifts are disabled");
	}
	vector<vector<double> > all_cb_exchange_SX(num_cb_loops);
	vector<vector<double> > all_vb_exchange_SX(num_vb_loops);	
	for(int ii = 0; ii < num_cb_loops + num_vb_loops; ii++) {		
		// ------------------------------------------			
		// get right precomputed fermi factors
		// ------------------------------------------
		vector<double>*                 exchange_SX;
		const vector<double>*           fermi_stats;
		const CLCCoulombLambda<double>* coulomb_lambda;
		sout.str("");		
		if(ii < num_cb_loops) {						
			fermi_stats    = &all_cb_fermi_stats[ii];
			exchange_SX    = &all_cb_exchange_SX[ii];
			coulomb_lambda = &ingredients.get_coulomb_matrix_element_cb_cb(ii);
			sout << "  cb " << ii;
		} else {
			fermi_stats    = &all_vb_fermi_stats[ii - num_cb_loops];
			exchange_SX    = &all_vb_exchange_SX[ii - num_cb_loops];
			coulomb_lambda = &ingredients.get_coulomb_matrix_element_vb_vb(ii - num_cb_loops);
			sout << "  vb " << ii - num_cb_loops;		
		}
		(*exchange_SX).assign(number_of_points, 0.0);
		// ------------------------------------------
		// calculate screening contribution matrix
		// ------------------------------------------
		if(Configuration::get_instance()->get("clc_include_exchange_shifts") == 1.0) {		
			#pragma omp parallel for 
			for(int kkii = 0; kkii < (signed)number_of_points; kkii++) {
				for(unsigned int kkjj = 0; kkjj < number_of_points; kkjj++) {					
					(*exchange_SX)[kkii] += coulomb_lambda->get_coulomb_lambda(k_point_coords[kkii], k_point_coords[kkjj], screening_curve)					                      
					                      * dense_domain.get_point(kkjj).get_weight()
					                      * (*fermi_stats)[kkjj]
					                      * matrix_weight_factor; 
				}	
			}
			sout << " = - " << (*exchange_SX)[0] << " [eV], last = - " << (*exchange_SX).back() << " [eV]";
			TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());			
		}	
	}		
	// ------------------------------------------
	// write exchange stuff to file if requested
	// ------------------------------------------
	if(Configuration::get_instance()->get("output_clc_exchange_shifts") == 1.0) {
		static int exchange_counter = 0;
		ostringstream fname;
		fname << "exchange_energy_" << exchange_counter << ".dat";
		ofstream fout(fname.str().c_str());
		if(fout) {
			const int width = 18;
			// --------------------------------------------
			// write header
			// --------------------------------------------
			fout << "# k " << setw(width - 3) << " ";
			for(int kk = -1; kk < (signed)number_of_points; kk++) {
				if(kk >= 0) {
					fout << setw(width) << k_point_coords[kk] << " ";
				}
				for(int ii = 0; ii < num_cb_loops; ii++) {
					if(kk == -1) {
						ostringstream stmp;
						stmp << "cb" << ii;
						fout << setw(width) << stmp.str() << " ";						
					} else {
						fout << setw(width) << all_cb_exchange_SX[ii][kk] << " ";	
					}
				}
				for(int ii = 0; ii < num_vb_loops; ii++) {
					if(kk == -1) {
						ostringstream stmp;
						stmp << "vb" << ii;
						fout << setw(width) << stmp.str() << " ";						
					} else {
						fout << setw(width) << all_vb_exchange_SX[ii][kk] << " ";	
					}						
				}
				fout << "\n";
			}
			fout.close();
		} else {
			TDKP_GENERAL_EXCEPTION("can not write to file " << fname.str());	
		}
		exchange_counter++;	
	}
	
	// ----------------------------------------------
	// calculate coulomb hole self energy 
	// ----------------------------------------------
	vector<double> delta_coulomb_hole_shift(num_vb_bands, 0.0);	
	if(Configuration::get_instance()->get("clc_include_coulomb_hole_self_energy") == 1.0) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLC: coulomb hole self-energies are:");
		for(unsigned int ii = 0; ii < num_vb_bands; ii++) {			
			const CLCRadialCoulombLambda& cmat = ingredients.get_coulomb_matrix_element_vb_vb(ii);
			TDKP_BOUNDS_ASSERT(k_point_coords[0] >= 0.0,"");
			// -------------------------------------------
			// quantum well
			// -------------------------------------------
			if(dense_domain.get_dimension() == 2) { 
				for(unsigned int kk = k_point_coords[0] == 0 ? 1:0; kk < number_of_points; kk++) {
					delta_coulomb_hole_shift[ii] += 
						cmat.get(k_point_coords[kk]) // thats 2*pi *exp(-q ..) 						
						* matrix_weight_factor * 2.0 * constants::pi // 2pi is because matrix weight factor removes angular integration (angular integration is there condensed in get_lambda ...)
				 		* dense_domain.get_point(kk).get_weight() // this is 2pi k dk
				 		/ k_point_coords[kk] // this is to remove k from k dk
				 		* (k_point_coords[kk] / screening(k_point_coords[kk]) - 1.0); // screening returns q*eps_q
				}
				
			// -------------------------------------------
			// quantum wire
			// -------------------------------------------				
			} else if(dense_domain.get_dimension() == 1) {
				for(unsigned int kk = k_point_coords[0] == 0 ? 1:0; kk < number_of_points; kk++) {
					delta_coulomb_hole_shift[ii] +=
						cmat.get(k_point_coords[kk]) // thats 2*K0(-q ..) 						
						* matrix_weight_factor * 2.0 // 2 is because matrix weight factor removes integration for [-inf,0[  
				 		* dense_domain.get_point(kk).get_weight() // this is 2 dk				 		
				 		* (1.0 / screening(k_point_coords[kk]) - 1.0); // screening returns eps_q in the wire case
				}
			// -------------------------------------------
			// bulk material
			// -------------------------------------------				
			} else if(dense_domain.get_dimension() == 3) {
				for(unsigned int kk = k_point_coords[0] == 0 ? 1:0; kk < number_of_points; kk++) {
					delta_coulomb_hole_shift[ii] += 
						cmat.get(k_point_coords[kk]) // thats 4*pi*<g1g2g3g4>  						
						* matrix_weight_factor * 4.0 * constants::pi // 4pi is because matrix weight factor removes angular integration (angular integration is there condensed in get_lambda ...)
				 		* dense_domain.get_point(kk).get_weight() // this is 4pi k^2 dk
				 		/ (k_point_coords[kk] * k_point_coords[kk]) // this is to remove k^2 from k^2 dk
				 		* ((k_point_coords[kk] *k_point_coords[kk]) / screening(k_point_coords[kk]) - 1.0); // screening returns q^2*eps_q				
				 }				
			}	
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLC: valence band " << ii << " = " << delta_coulomb_hole_shift[ii] << " [eV]");					
		}
	} else {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "CLC: disabled coulomb hole self-energy by user request"); 	
	}
	// ------------------------------------------
	// time estimator
	// ------------------------------------------	
	double start_tic   = TimeMeasurements::tic();
	double eta_time    = 0.0;
	double tot_time    = 0.0;			
				
	// ------------------------------------------
	// for cb band
	// ------------------------------------------
	for(unsigned int cc = 0; cc < num_cb_bands; cc++) {
				
		// ------------------------------------------------
		// eff mass band: disperison is same for both spins
		// ------------------------------------------------ 		
		unsigned int cc_idx = (cb_effmass  && !vb_effmass) ? ((cc - (cc % 2)) / 2) : cc;
		// ------------------------------------------------
		// only value index in matrix element differs
		// ------------------------------------------------
		unsigned int vl_idx = (cb_effmass  && !vb_effmass) ? cc % 2:0;		
				
		const vector<double>& cb_band = all_cb_band[cc_idx];
		const vector<double>& cb_fermi_stats = all_cb_fermi_stats[cc_idx];
		const vector<double>& cb_exchange_SX = all_cb_exchange_SX[cc_idx];
		
		
		// ------------------------------------------
		// for every vb band
		// ------------------------------------------
		for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
				
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "evaluating transition for cc = " << cc << ", vv = " << vv << " for " << omega_num << " omega points.\neta time = " << static_cast<int>(eta_time) << " [s], of total = " <<  static_cast<int>(tot_time) << " [s]");				
						
			const vector<double>& vb_band = all_vb_band[vv];
			const vector<double>& vb_fermi_stats = all_vb_fermi_stats[vv];
			const vector<double>& vb_exchange_SX = all_vb_exchange_SX[vv];
			const CLCCoulombLambda<double>& coulomb_Lambda_vb_cb =	
				ingredients.get_coulomb_matrix_element_vb_cb(vv,cc_idx);
							
			// ------------------------------------------
			// precompute coulomb scattering matrix Ni(ki,kj)
			// ------------------------------------------
			if(include_self_consistent_polarization || use_pade_approximation) {				
				TimeMeasurements::get_instance().start("shf coulomb vb-cb matrix");		
				#pragma omp parallel for
				for(int kkii = 0; kkii < (signed)number_of_points; kkii++) {
					for(unsigned int kkjj = 0; kkjj < number_of_points; kkjj++) {
						vb_cb_coulomb_scattering_matrix[kkii * number_of_points + kkjj]
							= coulomb_Lambda_vb_cb.get_coulomb_lambda(k_point_coords[kkii],k_point_coords[kkjj], screening_curve);											
					}
				}
				TimeMeasurements::get_instance().stop("shf coulomb vb-cb matrix");
			}										
										 
			// ------------------------------------------
			// calculate renormalized transition energies
			// ------------------------------------------
			for(int ii = 0; ii < (signed)number_of_points; ii++) {
				omega_trans[ii] = (cb_band[ii] - vb_band[ii] - cb_exchange_SX[ii] - vb_exchange_SX[ii] + delta_coulomb_hole_shift[vv]) / constants::hbar;			 				 		
			}
		
			// ------------------------------------------
			// for every polarization
			// ------------------------------------------
			for(unsigned int pp = pp_start; pp <= pp_end; pp++) {
								
				const CLCLinearCurve& pcv_curve = ingredients.get_momentum_matrix_element(cc_idx, vv, pp, vl_idx);
				
				// -------------------------------------------------
				// extract momentum matrix element at k-space points
				// -------------------------------------------------
				#pragma omp parallel for
				for(int ii = 0; ii < (signed)number_of_points; ii++) {					
					momentum_matrix_element[ii] = pcv_curve(k_point_coords[ii]);	
				}

				#pragma omp parallel 				
				{
					// ------------------------------------------
					// set up thread stuff
					// ------------------------------------------
					int thread_num = omp_get_thread_num();
					int omega_thread_start, omega_thread_end;
					get_omp_range(omega_thread_start, omega_thread_end, omega_num);
					
					vector<cplx>&   omega_cv_k              = omp_omega_cv_k[thread_num]; 
					vector<double>& spont_stats_cv_k        = omp_spont_stats_cv_k[thread_num]; 	
					vector<cplx>&   lambda_cv_k             = omp_lambda_cv_k[thread_num]; 	  
					vector<cplx>&   lambda_cv_k_from_matrix = omp_lambda_cv_k_from_matrix[thread_num]; 
					RMatrix<cplx>&  lambda_cv_k_matrix      = omp_lambda_cv_k_matrix[thread_num];
					vector<cplx>&   dipole_mu_k             = omp_dipole_mu_k[thread_num]; 						

					// ------------------------------------------
					// for every frequency
					// ------------------------------------------
					#pragma omp for schedule(dynamic)
					for(int ww = 0; ww < (signed)omega_num; ww++) {						
						const double omega = omega_min + (omega_max - omega_min) / static_cast<double>(omega_num - 1) * ww;											
						// ------------------------------------------
						// calculate omega_cv_k
						// calculate dipole matrix element
						// and solve for lambda_cv_k (not self consistent)						
						// ------------------------------------------
						for(unsigned int ii = 0; ii < number_of_points; ii++) {
							TDKP_BOUNDS_ASSERT((cb_fermi_stats[ii] + vb_fermi_stats[ii] - 1.0) >= -1.0, "");
							TDKP_BOUNDS_ASSERT((cb_fermi_stats[ii] + vb_fermi_stats[ii] - 1.0) <= 1.0, "");						 
							omega_cv_k[ii] = i 
							               * (cb_fermi_stats[ii] + vb_fermi_stats[ii] - 1.0)
							               / (i * constants::hbar * (omega_trans[ii] - omega) + constants::hbar * ingredients.get_homogeneous_broadening());
							spont_stats_cv_k[ii] = 1.0 / constants::hbar 
							               * (cb_fermi_stats[ii] * vb_fermi_stats[ii]) * ingredients.get_homogeneous_broadening()
							               / ((omega_trans[ii] - omega) * (omega_trans[ii] - omega) + ingredients.get_homogeneous_broadening() * ingredients.get_homogeneous_broadening());
							dipole_mu_k[ii] = constants::ec 
							                / (i * constants::m0 * omega) 
							                * momentum_matrix_element[ii];
							lambda_cv_k[ii] = - dipole_mu_k[ii] * omega_cv_k[ii];							                						               												               					 						               																		 						               
						}					
						
															
						// ------------------------------------------
						// build matrix and solve for lambda_cv_k
						// ------------------------------------------
						if(include_self_consistent_polarization) {						
							//#pragma omp parallel for
							for(int kkii = 0; kkii < (signed)number_of_points; kkii++) {																		
								// coulomb contributions
								for(unsigned int kkjj = 0; kkjj < number_of_points; kkjj++) {
									lambda_cv_k_matrix(kkii,kkjj) = 
										matrix_weight_factor 
									  * omega_cv_k[kkii]
									  * vb_cb_coulomb_scattering_matrix[kkii * number_of_points + kkjj]
									  * dense_domain.get_point(kkjj).get_weight();								
								}
								// diagonal 1
								lambda_cv_k_matrix(kkii,kkii) = lambda_cv_k_matrix(kkii,kkii) + 1.0;
								// lambda without coulomb is a first guess. the guess it very good for gmres and reduces 1/3 of the iterations
								lambda_cv_k_from_matrix[kkii] = lambda_cv_k[kkii]; 						
							}
							// ------------------------------------------------
							// solver self consistency equation using requested solver
							// ------------------------------------------------
							if(linear_equation_solver == 1) {					
								// lapack			
								solve(lambda_cv_k_matrix, lambda_cv_k, lambda_cv_k_from_matrix);
							} else if(linear_equation_solver == 2) {
								// gmres
								gmres.solve(25, lambda_cv_k_matrix, lambda_cv_k, lambda_cv_k_from_matrix);
							} else {
								// ------------------------------------------------
								// both, solve & compare
								// ------------------------------------------------
								vector<cplx> tmp_lambda(number_of_points);
								double tic1 = TimeMeasurements::tic();
								solve(lambda_cv_k_matrix, lambda_cv_k, lambda_cv_k_from_matrix);
								double tic2 = TimeMeasurements::tic();
								gmres.solve(25, lambda_cv_k_matrix, lambda_cv_k, tmp_lambda);
								double tic3 = TimeMeasurements::tic();
								
								double rel = 0.0;
								for(unsigned int ii = 0; ii < number_of_points; ii++) {
									rel += tdkp_math::abs(tmp_lambda[ii] - lambda_cv_k_from_matrix[ii]) / tdkp_math::abs(lambda_cv_k_from_matrix[ii]);								
								}
								TDKP_LOGMSG(LOG_INFO_DEVEL2, "lapack =  " << tic2 - tic1 << " [s], gmres = " << tic3 - tic2 << " [s], avg rel err = " << rel / number_of_points);
							}							 
							// --------------------------------------------
							// debugging / output requests
							// --------------------------------------------
							if(ww == (signed)output_omega_idx) {
								static int counter = 0;
								ostringstream fname;
								// ---------------------------------------
								// write matrix
								// ---------------------------------------
								fname << "matrix_" << counter << ".dat";
								ofstream fout(fname.str().c_str());
								for(unsigned int mii = 0; mii < number_of_points; mii++) {
									for(unsigned int mjj = 0; mjj < number_of_points; mjj++) {
										fout << lambda_cv_k_matrix(mii,mjj).real() << " " << lambda_cv_k_matrix(mii,mjj).imag() << " ";
									}
									fout << "\n";
								}
								fout.close(); 							
								fname.str("");
								fname << "lambda_" << counter << ".dat";
								fout.open(fname.str().c_str());	
								for(unsigned int mii = 0; mii < number_of_points; mii++) {
									fout << k_point_coords[mii] << " "	
									     << lambda_cv_k[mii].real() << " "	
									     << lambda_cv_k[mii].imag() << " " 
									     << lambda_cv_k_from_matrix[mii].real() << " " 
									     << lambda_cv_k_from_matrix[mii].imag() << " "
									     << omega_cv_k[mii].real() << " " 
									     << omega_cv_k[mii].imag() << " " 
									     << "\n";
								}						
								fout.close();
								counter++;
							}	
							lambda_cv_k = lambda_cv_k_from_matrix;						
						} else if(use_pade_approximation) {
							// using pade approximation, means:
							// inserting lambda_k_k from fc. into lambdakk for coulomb terms
							for(int kkii = 0; kkii < (signed)number_of_points; kkii++) {
								lambda_cv_k_from_matrix[kkii] = lambda_cv_k[kkii];														
								// coulomb contributions
								for(unsigned int kkjj = 0; kkjj < number_of_points; kkjj++) {																
									lambda_cv_k_from_matrix[kkii] -= 									
										matrix_weight_factor 
									  * omega_cv_k[kkii]
									  * vb_cb_coulomb_scattering_matrix[kkii * number_of_points + kkjj]
									  * dense_domain.get_point(kkjj).get_weight()
									  * lambda_cv_k[kkjj];
								}
								
							}	
							lambda_cv_k = lambda_cv_k_from_matrix;
						}		
						// ------------------------------------------
						// calculate absorption from lambda_cv_k
						// ------------------------------------------
						for(unsigned int ii = 0; ii < number_of_points; ii++) {						
							suszeptibility[ww * 3 + pp] += conj(dipole_mu_k[ii]) 
														 * lambda_cv_k[ii]
							                             * dense_domain.get_point(ii).get_weight();
							// free carrier spont emission (classical) including
							// bandgap renormalization terms
							// formula is the same as for SLC						                             
							spont_emission[ww * 3 + pp] += real(conj(dipole_mu_k[ii]) * dipole_mu_k[ii])						                             
							                             * omega * omega * omega
							                             * spont_stats_cv_k[ii]
							                             * dense_domain.get_point(ii).get_weight();
						}
					}
				} // end parallel block
				if(Logger::get_instance()->get_level() >= LOG_INFO) {
					double toc = TimeMeasurements::tic() - start_tic;
					unsigned int loop = (cc * num_vb_bands * 3 + vv * 3 + pp) + 1;
					unsigned int num_loops = num_cb_bands * num_vb_bands * 3;
					eta_time = (num_loops - loop) * (toc / loop);
					tot_time = (num_loops) * (toc / loop); 
				}
			}
		}
	}

	// ---------------------------------------------
	// apply constants
	// ---------------------------------------------
	for(unsigned int ii = 0; ii < suszeptibility.size(); ii++) {
		suszeptibility[ii] *= prefactor * dim_prefactor;
		spont_emission[ii] *= spont_prefactor * dim_prefactor;					
	}
		
	TimeMeasurements::get_instance().stop("shf radial evaluation");
	
	return new CLCOpticalResults(
		ingredients,
		omega_min,
		omega_max,
		suszeptibility,
		spont_emission,
		cb_fermi,
		vb_fermi,
		temperature
	);			
						
} 


}
