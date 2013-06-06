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

#include <sstream>

#include "tdkp/probdefs/BulkBandstructureSolver.h"
#include "tdkp/common/Configuration.h"
#include <omp.h>
	
namespace tdkp {
	
extern "C" {
void F77_Name(zheevx)(char* jobz,
	     char* range,
	     char* uplo,
	     int* n,
	     cplx* A,
	     int* lda,
	     double* VL,
	     double* VU,
	     int* IL,
	     int* IU,
	     double* ABSTOL,
	     int* M,
	     double* W,
	     cplx* Z,
	     int* LDZ,
	     cplx* WORK,
	     int* LWORK,
	     double* RWORK,
             int* IWORK,
	     int* IFAIL,
	     int* INFO );
void F77_Name(zheevr)(char* jobz,
	     char* range,
	     char* uplo,
	     int* n,
	     cplx* A,
	     int* lda,
	     double* VL,
	     double* VU,
	     int* IL,
	     int* IU,
	     double* ABSTOL,
	     int* M,
	     double* W,
	     cplx* Z,
	     int* LDZ,
	     int* ISUPPZ, 
	     cplx* WORK,
	     int* LWORK,
	     double* RWORK,
	     int* LRWORK,
         int* IWORK,
         int* LIWORK,	    
	     int* INFO );
void F77_Name(zheev)(char* jobz, char* uplo, int* n, complex<double> *a, int* lda, double* w, complex<double>* work, int* lwork, double* rwork, int* info);	     
	     
}


BulkBandstructureSolver::BulkBandstructureSolver(KPMatrixBase* matrix) 
: solution(0)
{
	this->kp_matrix = matrix;
	this->kp_matrix->calculate();
}

BulkBandstructureSolver::~BulkBandstructureSolver() {
	if(solution != 0) {
		delete solution;
		solution = 0;	
	}
}

/** solve bandstructure in a given direction */
void BulkBandstructureSolver::solve(Vector3D direction, double k_min, double k_max, int num_k_values) {

	TDKP_ASSERT(k_min < k_max, "k_min < k_max"); 
	DomainMaster domain;
	create_3D_domain_radial(domain, direction, k_min, k_max, num_k_values);	
	solve(domain);
	
}

/** solve bandstructure for a given domain */
void BulkBandstructureSolver::solve(DomainMaster& domain) {

	// generate solution space
	int num_bands = this->kp_matrix->get_number_of_bands();
	if(this->solution != 0) {
		delete solution; solution = 0;	
	}
	this->solution = new BandstructureDomain<cplx>(this->kp_matrix->get_number_of_bands(), this->kp_matrix->get_number_of_bands(), 1, domain);

	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "BulkBandstrctureSolver: solving bandstructure for " << domain.get_number_of_points() << " points.");
	

#pragma omp parallel 
	{
	
		// prepare necessary variables	
		RMatrix<cplx>  matrix(num_bands, num_bands);
		vector<cplx>   ev;
		vector<cplx>   evec;
		RMatrix<cplx>* tmp = 0;
	
		double kval;
		cplx i(0.0, 1.0);
		int current = 0; int next = 0; int total = domain.get_number_of_points();
		const int bar_min = 20000; 
		if(omp_get_thread_num() == 0) {
			if(total > bar_min) {
				Logger::get_instance()->init_progress_bar("BulkBandstrctureSolver: calculating bulk bandstructure", total);	
			}
		}
		
		int start, end;
		get_omp_range(start,end, domain.get_number_of_points());
#pragma omp critical (wurstsalat)
{
		TDKP_LOGMSG(LOG_INFO, "thread " << omp_get_thread_num() << " gets points " << start << " to " << end);		
}			
		// -----------------------------------------
		// loop over domain points
		// -----------------------------------------
		for(int kk = start; kk < end; kk++) {
			 
			// ------------------------------------------
			// get current domain point and calculate direction and value
			// -----------------------------------------
			const DomainPoint& point = domain.get_point(kk);		
			Vector3D normalized_direction(point.get_coord(0), point.get_coord(1), point.get_coord(2));
			kval = point.get_coord_abs(); // current k
			if(kval != 0.0) {
				normalized_direction.normalize(); // current direction
			}
			// set very small offset from 0 point (so matrix elements get smother as we lift the degeneracy ...)
			if(kval == 0.0) {								
				kval = Configuration::get_instance()->get("bulksolver_gamma_point_offsetting_value");	
				normalized_direction = Vector3D(1.0, 1.0, 1.0);
				normalized_direction.normalize();	
			}
			 
			// reset matrix
			for(int ii = 0; ii < num_bands; ii++) {
				for(int jj = 0; jj < num_bands; jj++) {
					matrix(ii,jj) = 0.0;
				}
			}
			// fill matrix with second order
			for(int aa = 0; aa < 3; aa++) {
				for(int bb = 0; bb < 3; bb++) {
					tmp = this->kp_matrix->get_second_order_matrix(aa, bb);
					cplx ikaaikbb = kval * kval * (i * normalized_direction(aa) * i * normalized_direction(bb));
					for(int ii = 0; ii < num_bands; ii++) {
						for(int jj = 0; jj < num_bands; jj++) {
							// the minus is because H(2) is the partially integrated H(2) ...
							// means from -d^2/dx^2 -> d/dx d/dx
							matrix(ii,jj) += - (*tmp)(ii,jj) * ikaaikbb;
						}
					}
					delete tmp; tmp = 0;
				}
			}
			// first order
			for(int aa = 0; aa < 3; aa++) {
				for(int gg = 0; gg < 2; gg++) {
					tmp = this->kp_matrix->get_first_order_matrix(gg, aa);
					cplx ikaa = - kval * i * normalized_direction(aa);
					for(int ii = 0; ii < num_bands; ii++) {
						for(int jj = 0; jj < num_bands; jj++) {
							matrix(ii,jj) += (*tmp)(ii,jj) * ikaa;
						}
					}
					delete tmp; tmp = 0;
				}
			}
			// zero order
			tmp = this->kp_matrix->get_zero_order_matrix();
			for(int ii = 0; ii < num_bands; ii++) {
				for(int jj = 0; jj < num_bands; jj++) {
					matrix(ii,jj) += (*tmp)(ii,jj);
				}
			}
			delete tmp; tmp = 0;
			//cout << kk << "\n";
			RMatrix<cplx>::get_eigensystem(matrix, ev, evec);
			//this->solve_for_eigenvalues(matrix, num_bands, ev, evec);
			for(int bb = 0; bb < num_bands; bb++) {
				this->solution->add_eigensolution(kk, bb, new EigenSolution<cplx>(ev[bb], &evec[bb*num_bands], 1, num_bands));
			}
			if(omp_get_thread_num() == 0) {
				if(current == next && total > bar_min) {
					next = Logger::get_instance()->set_progress_bar(current, total);
				}
				current++;
			}
				
		}
		if(omp_get_thread_num() == 0 && total > bar_min) {
			Logger::get_instance()->end_progress_bar();			
		}
	} // end omp parallel region

	//delete[] ev;
	//delete[] evec;
	//delete[] matrix;
}

/** caluclate group velocity
 * 
 * for mister bruggersan: group velocity given as dE/dk = <xT_k dH/dk x_k>
 * rotate kp matrix to desired direction. i always assume it to be the [100] dir
 * 
 * kp dispersion relation is given by
 * Hkp = sum H2ij ki kj + sum H1i ki + H0
 * so, dHkp / dkg = 2 Hgg kg + sum (Hig + Hgi) ki 
 */ 
void BulkBandstructureSolver::calculate_group_velocity(Vector3D direction, double k_min, double k_max, int num_k_values) {
	
	// calculate bandstructure
	direction.normalize(); 
	this->solve(direction, k_min, k_max, num_k_values);
	// get bandstructure
	const Bandstructure<cplx>* bandstructure = &this->get_bandstructure();
	// open outstream for output
	ofstream fout("gruppengesch.dat"); 

	int num_bands = this->kp_matrix->get_number_of_bands();	

	// prepare necessary variables
	cplx* matrix = new cplx[num_bands * num_bands];
	RMatrix<cplx>* tmp       = 0;
	RMatrix<cplx>* tmp_other = 0;

	double kval;
	cplx i(0.0, 1.0);
	// loop over k values
	for(int kk = 0; kk < num_k_values; kk++) {
		// current k
		kval = k_min * (1.0 - (double(kk) / double(num_k_values - 1)))
		     + k_max * (double(kk) / double(num_k_values - 1));
		fout << kval << " \t";
		// loop over directions (D_DX, D_DZ etc)
		for(int dd = 0; dd < 3; dd++) {
			// reset matrix
			for(int ii = 0; ii < num_bands * num_bands; ii++) {
				matrix[ii] = 0.0;
			}
			// fill matrix with second order diagonal terms
			tmp = this->kp_matrix->get_second_order_matrix(dd, dd);
			cplx ikk = - kval * direction(dd);
			for(int ii = 0; ii < num_bands; ii++) {
				for(int jj = 0; jj < num_bands; jj++) {
					// the minus is because H(2) is the partially integrated H(2) ...
					// means from -d^2/dx^2 -> d/dx d/dx
					matrix[ii + jj * num_bands] += - 2.0 * (*tmp)(ii,jj) * ikk;
				}
			}
			delete tmp; tmp = 0;
			// fill matrix with second order offdiagonal terms
			for(int ee = 0; ee < 3; ee++) {
				if(ee != dd) {
					tmp       = this->kp_matrix->get_second_order_matrix(dd, ee);
					tmp_other = this->kp_matrix->get_second_order_matrix(ee, dd);
					cplx ikk = - kval * direction(ee);
					for(int ii = 0; ii < num_bands; ii++) {
						for(int jj = 0; jj < num_bands; jj++) {
							// the minus is because H(2) is the partially integrated H(2) ...
							// means from -d^2/dx^2 -> d/dx d/dx
							matrix[ii + jj * num_bands] += - ((*tmp)(ii,jj) + (*tmp_other)(ii,jj)) * ikk;
						}
					}	
					delete tmp; delete tmp_other;				
				}	
			}
			// first order
			for(unsigned int gg = 0; gg < 2; gg++) { 
				tmp = this->kp_matrix->get_first_order_matrix(gg, dd);
				for(int ii = 0; ii < num_bands; ii++) {
					for(int jj = 0; jj < num_bands; jj++) {
						matrix[ii + jj * num_bands] += (*tmp)(ii,jj) * (-i);
					}
				}
				delete tmp; tmp = 0;	
			}
			// --------------------------------
			// calculate xT matrix x
			// --------------------------------
			vector<cplx> solvec(num_bands);
			vector<cplx> vec_tmp(num_bands);
			cplx res, norm;
			cplx xT_dHdk_x, xTx;
			// for every band
			for(int ii = 0; ii < num_bands; ii++) {
				// get solution vector
				const EigenSolution<cplx>& sol = bandstructure->get_eigensolution(kk,ii);
				for(int jj = 0; jj < num_bands; jj++) {				
					solvec[jj]  = sol.get_node_value(0,jj);
					vec_tmp[jj]	= 0.0;
				}
				// calculate dHdk_x and store to vec_tmp
				for(int cc = 0; cc < num_bands; cc++) {
					for(int dd = 0; dd < num_bands; dd++) {
						vec_tmp[cc] += matrix[cc + dd * num_bands] * solvec[dd];
					}	
				}	
				// caluclate xTx and xT_dHdk_x
				res = norm = 0.0;
				for(int cc = 0; cc < num_bands; cc++) {
					res  += conj(solvec[cc]) * vec_tmp[cc];
					norm += conj(solvec[cc]) * solvec[cc]; 	
				} 
				// write group speed in m/s
				fout << res.real() / constants::hbar * 1.0e-9 << " \t";
			}
		
		}
		fout << "\n";
	}
	delete[] matrix;
	fout.close(); 
}

const BandstructureDomain<cplx>& BulkBandstructureSolver::get_bandstructure() const throw (Exception*){
	return *this->solution;
}


void BulkBandstructureSolver::solve_for_eigenvalues(cplx* matrix, int size, double* eigenvalues, cplx* eigenvectors) {

	// use zheevx or zheevr
	int  info         = 0;
	int  ev_found     = 0;
	vector<cplx> copy;
	for(int ii = 0; ii < size * size; ii++) {
		copy.push_back(matrix[ii]); 	
	}

	string solver_name;

	if(Configuration::get_instance()->get_config_value("eigensolver.inp", "lapack_hermitian_eigensolver") == 1) {	
		// ------------------------------------------------
		// calculate matrix norm and scale entries to one to strengthen convergence behaviour
		// ------------------------------------------------
		double norm  = 0.0; 
		double val;
		int    count = 0;
		for(int ii = 0; ii < size; ii++) {
			for(int jj = 0; jj < size; jj++) {
				val = abs(matrix[ii * size + jj] * conj(matrix[ii * size + jj]));
				norm += val;
				if(val > 0.0) {
					count++;	
				}
			}			 	
		}
		norm = sqrt(norm / count);
		if(norm > 0.0) {
			for(int ii = 0; ii < size; ii++) {
				for(int jj = 0; jj < size; jj++) {
					matrix[ii * size + jj] /= norm;
				}			 	
			}	
		} 
	
		// zheevx
		int     num      = size;
		char 	jobz     = 'V';             // compute eigenvalues and vectors
		char 	range    = 'A';
		char 	uplo     = 'U';
		double 	prec     = 1.0e-14;		
		int 	lwork    = 2*num + 5;
		cplx* 	work     = new cplx[lwork];
		double* rwork    = new double[7 * num];
		int*    iwork    = new int[5 * num];
		int*    ifail    = new int[num];		
		double 	vl       = - 100.0;
		double 	vu       = 100.0;
		int    	iu       = 0;
		int    	il       = 0;
	
		F77_Name(zheevx)(&jobz, &range, &uplo, &num, matrix, &num, &vl, &vu, &il, &iu, &prec, &ev_found, eigenvalues, eigenvectors, &num, work, &lwork, rwork, iwork, ifail, &info);
		
		solver_name = "zheevx";
		
		// scale eigenvalues back to initial norm
		if(norm > 0) {
			for(int ii = 0; ii < size; ii++) {
				eigenvalues[ii] *= norm;	
			}	
		}
		delete[] work;
		delete[] rwork;
		delete[] iwork;
		delete[] ifail;		
		
	} else if(Configuration::get_instance()->get_config_value("eigensolver.inp", "lapack_hermitian_eigensolver") == 2.0) {
		
		// --------------------------------------------------
		// using zheevr
		// --------------------------------------------------	
		int     num      = size;
		char 	jobz     = 'V';             // compute eigenvalues and vectors
		char 	range    = 'A';
		char 	uplo     = 'U';
		double 	prec     = 1.0e-14;		
		int 	lwork    = 5*num + 5;
		cplx* 	work     = new cplx[lwork];
		int     lrwork   = 25*num;
		double* rwork    = new double[lrwork];
		int     liwork   = 11*num;
		int*    iwork    = new int[liwork];		
		double 	vl       = - 100.0;
		double 	vu       = 100.0;
		int    	iu       = 0;
		int    	il       = 0;
		int*    isuppz   = new int[num * 2];
	        	    
		F77_Name(zheevr)(&jobz, &range, &uplo, &num, matrix, &num, &vl, &vu, &il, &iu, 
		        &prec, &ev_found, eigenvalues, eigenvectors, &num, isuppz,
		        work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
	
		solver_name = "zheevr";
	
		delete[] iwork;
		delete[] rwork;
		delete[] work;
	} else {
		
    	char jobz = 'V';
    	char uplo = 'U';
    	int  n    = size;
    	int lda   = size;    	
    	int lwork = 6 * n;
    	complex<double>* work = new complex<double>[lwork];
    	double* rwork         = new double[5*n];
    	int info;		
    	
    	for(int ii = 0; ii < n * n; ii++) {
    		eigenvectors[ii] = matrix[ii];	
    	}
    	F77_Name(zheev)(&jobz, &uplo, &n, eigenvectors, &lda, eigenvalues, work, &lwork, rwork, &info);	
    	    	 
		solver_name = "zheev";
		ev_found = size;    	    	 
    	    	   
    	delete[] work;
    	delete[] rwork;
		
	}	     

	if(info != 0) {
		cout << info << endl;
		for(int ii = 0; ii < size; ii++) {
			cout << "";
			for(int jj = 0; jj < size; jj++) {
				if(jj > 0) {
					cout << " \t";	
				}
				if(abs(copy[ii + jj * size]) == 0) {
					cout << "0";
				} else {
					cout << copy[ii + jj * size].real() << " ";
					if(copy[ii + jj * size].imag() >= 0) {
						cout << "+";
					}
					cout << copy[ii + jj * size].imag() << "i ";
				}
				
			}
			cout << ";\n";	
		}
		cout << "\n";
		ostringstream sout;
		sout << "lapack's " << solver_name << " returned an error during determination of eigenvalues";
		TDKP_GENERAL_EXCEPTION(sout.str());
				
	}
	if(ev_found != size) {
		TDKP_GENERAL_EXCEPTION("lapack's zheev. could not find all eigenvalues which is strange ... ");
	}

	// sort eigenvalues and eigenvectors
	this->sort_results(eigenvalues, eigenvectors, size);

}

/** calculates resulting effective masses from bulk bandstructure */
void BulkBandstructureSolver::analyze_effective_masses() const {
	if(solution == 0) {
		TDKP_GENERAL_EXCEPTION("nonexistent solution set requested");
	}
	BandstructureDomain<cplx>& set = *this->solution;
	const DomainMaster& domain = set.get_domain();
	if(set == 0) {
		Logger::get_instance()->emit(LOG_ERROR, "bandstructure is not of single dispersion type!");
		return;	
	}
	ostringstream sout;
	if(set.get_number_of_k_values() <= 1) {
		TDKP_GENERAL_EXCEPTION("not enough data points calculated!");
	}
	TDKP_ASSERT(domain.radial(), "sorry, but bandstructure must be radial");
	if(domain.get_first_point().get_coord_abs() != 0.0) {
		TDKP_GENERAL_EXCEPTION("bandstructure must start at k = 0");
	}	
	sout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
	     << "  effective masses from bandstructures\n\n";

	double mstar;
	double kk;
	int warnings = 0;
	for(unsigned int dd = 1; dd < domain.get_number_of_points(); dd++) {
		sout << setw(15) << domain.get_point(dd).get_coord_abs() << " ";
	    kk = domain.get_point(dd).get_coord_abs();
		for(int ii = 0; ii < set.get_number_of_bands(); ii++) {
			if(fabs(set.get_energy(dd,ii).imag()) > 1.0e-12) {
				warnings++;
				if(warnings < 10) {
					ostringstream wout; wout << "energy of band " << ii << " for k " << kk << " has nonneglible imaginary parts. this is really fishy.";
					Logger::get_instance()->emit(LOG_WARN, wout.str());
				}
				if(warnings == 9) {
					Logger::get_instance()->emit(LOG_WARN, "will not show any further warnings on imaginary energies");
				}
			}
			mstar = (constants::hbar_square * kk * kk) / (constants::m0 * 2.0 * (set.get_energy(dd,ii).real() - set.get_energy(0,ii).real()));
			sout << setw(15) << mstar << " ";
		}
		sout << "\n";
	}
	Logger::get_instance()->emit(LOG_INFO, sout.str());
}

/** sorts the eigenvalues in ascending order */
void BulkBandstructureSolver::sort_results(double* eigenvalues, cplx* eigenvectors, const int size) const {
	// build index
	vector<int>    indices;
	vector<double> energies;
	vector<cplx>   eigenvectors_copy;
	for(int ii = 0; ii < size; ii++) {
		energies.push_back(eigenvalues[ii]);
		indices.push_back(ii); 	
	}
	// copy eigenvector
	for(int ii = 0; ii < size * size; ii++) {
		eigenvectors_copy.push_back(eigenvectors[ii]);	
	}
	// sort indices (and energies in ascending order)
	tdkp_math::tracked_sort<int,double>(indices.begin(), indices.end(), energies.begin(), energies.end(),true);
	// copy eigenvectors back
	for(int ii = 0; ii < size; ii++) {
		eigenvalues[ii] = energies[ii];
		for(int jj = 0; jj < size; jj++) {
			eigenvectors[ii*size+jj] = eigenvectors_copy[indices[ii]*size + jj];	
		}
	}
	
}

} // end namespace tdkp
