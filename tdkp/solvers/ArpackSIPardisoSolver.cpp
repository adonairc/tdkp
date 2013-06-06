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




#include <iomanip>
#include <iostream>

#include "tdkp/solvers/ArpackSIPardisoSolver.h"
#include "tdkp/common/Configuration.h"

// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING ARPACK
// -----------------------------------------------
#ifdef EIGNSOLV_INCLUDE_ARPACK

using namespace tdkp;

// -------------------------------------------------------------------------
// fortran function definitions
// -------------------------------------------------------------------------
extern "C" {

	// ------------------------------------------------
	// arpack complex driver routine
	// ------------------------------------------------
	int znaupd_(int*    ido,      	// reverse communication flag
				char*   bmat,     	// type of matrix B
				int*    n,        	// dimension of the eigenproblem,
				char*   which,    	// which eigenvalues to calculate
				int*    nev,      	// number of eigenvalues to calcualte
				double* tol,      	// tolerance
				cplx*   resid,    	// complex array of length N containing the residuals
				int*    ncv,      	// number of columns of matrix V
				cplx*   v,      	// final set of arnoldi basis vectors
				int*    ldv,      	// leading dimension of V
				int		iparam[11], // array of length 11 contains control params
				int     ipntr[14],  // array of length 14, used as internal pointer in workd
				cplx*   workd,		// complex array of length 3*N
				cplx*   workl, 		// complex array of length lworkl
				int*    lworkl, 	// min 3*NCV**2 + 5 * NCV
				double* rwork, 		// double array of length ncv
				int*    info);

	// ------------------------------------------------
	// arpack output processing routine to calculate eigenvectors
	// ------------------------------------------------
	int zneupd_(int*    rvec,  		// logical whether to compute eigenvectors
				char*   howmny, 	// form of vectors to compute
				int*    select,     // array of length NCV: which evec to compute
				cplx*  	evals,      // array of length NEV + 1 (on exit, the eigenvalues are stored here)
				cplx*   evecs,      // array of length N * NEV
				int*    ldz,        // leading dimension of array evecs
				cplx*   sigma,      // shift
				cplx*   workev,     // array of length 2 * ncv
				// rest is equal as znaupd
				char*   bmat,
				int*    n,
				char*   which,
				int*    nev,
				double* tol,
				cplx*   resid,
				int*    ncv,
				cplx*   v,
				int*    ldv,
				int		iparam[11],
				int     ipntr[14],
				cplx*   workd,
				cplx*   workl,
				int*    lworkl,
				double* rwork,
				int*    info);

	void cstatn_(); // init statistics to zero
	void dstats_(); // init symmetric double statistics to zero
	/*  nopx   = total number of user OP*x operation
	 *  nbx    = total number of user B*x operation (same as copy when B = I)
	 *  nrorth = total number of reorthogonalization steps taken
	 *  nitref = total number of it. refinement steps in reorthogonalization
	 *  nrstrt = total number of restart steps
	 */
	/*
	c  tcaupd = total time spent in CNAUPD.
	c  tcaup2 = total time spent in CNAUP2.
	c  tcaitr = total time spent in the basic Arnoldi iteration loop,
	c           including iterative refinement in CNAITR.
	c  tceigh = total time spent in computing the Hessenberg eigenvalue
	c           subproblem at each iteration.
	c  tcgets = total time spent in getting the shifts at each iteration.
	c  tcapps = total time spent in applying the shifts at each iteration.
	c  tcconv = total time spent in convergence test at each iteration.
	c
	c==================
	c=== User time  ===
	c==================
	c
	c  tmvopx = total time spent in computing Y = OP * X
	c  tmvbx  = total time spent in computing Y = B * X
	*/
	extern struct {
		int nopx, nbx, nrorth, nitref, nrstrt;
		float tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
              tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
              tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
              tmvopx, tmvbx, tgetv0, titref, trvec;
	} timing_;

	// --------------------------------------------------------------
	// arpack real symmetric driver routine
	// --------------------------------------------------------------
	void dsaupd_(
				int* 	ido,  		// reverse comm
				char* 	bmat, 		// type of inputmatrix: B = 'I' / 'G'
				int*    n,        	// dimension of the eigenproblem,
				char*   which,    	// which eigenvalues to calculate
				int*    nev,      	// number of eigenvalues to calcualte
				double* tol,      	// tolerance
				double* resid,    	// complex array of length N containing the residuals
				int*    ncv,      	// number of columns of matrix V
				double* v,      	// final set of arnoldi basis vectors
				int*    ldv,      	// leading dimension of V
				int		iparam[11], // array of length 11 contains control params
				int     ipntr[14],  // array of length 14, used as internal pointer in workd
				double* workd,		// complex array of length 3*N
				double* workl, 		// complex array of length lworkl
				int*    lworkl, 	// min 3*NCV**2 + 5 * NCV
				int*    info);

	// ------------------------------------------------
	// arpack output processing routine to calculate eigenvectors from symmetric stuff
	// ------------------------------------------------
	int dseupd_(int*    rvec,  		// logical whether to compute eigenvectors
				char*   howmny, 	// form of vectors to compute
				int*    select,     // array of length NCV: which evec to compute
				double* evals,      // array of length NEV + 1 (on exit, the eigenvalues are stored here)
				double* evecs,      // array of length N * NEV
				int*    ldz,        // leading dimension of array evecs
				double* sigma,      // shift
				// rest is equal as znaupd
				char*   bmat,
				int*    n,
				char*   which,
				int*    nev,
				double* tol,
				double* resid,
				int*    ncv,
				double* v,
				int*    ldv,
				int		iparam[11],
				int     ipntr[14],
				double* workd,
				double* workl,
				int*    lworkl,
				int*    info);
}


struct ArpackError {
	int code;
	const char* text;	
};

static ArpackError arpack_errors[] = {
	{1,     "Maximum number of iterations taken. "},
	{2,     "No longer an informational error. Deprecated starting with release 2 of ARPACK."},
	{3,     "No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration. One possibility is to increase the size of NCV relative to NEV."},
	{-1,    "Matrix size N must be positive."},
	{-2,    "Requested number of eigenvalues must be positive"},
	{-3,    "arpack_num_arnoldi_vecs - number of eigenvalues must be >= 1 and less than or equal to N."},
	{-4,    "The maximum number of Arnoldi update iteration must be greater than zero."},
	{-5,    "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'"},
	{-6,    "BMAT must be one of 'I' or 'G'."},
	{-7,    "Length of private work array is not sufficient."},
	{-8,    "Error return from LAPACK eigenvalue calculation"},
	{-9,    "Starting vector is zero."},
	{-10,   "IPARAM(7) must be 1,2,3."},
	{-11,   "IPARAM(7) = 1 and BMAT = 'G' are incompatible."},
	{-12,   "IPARAM(1) must be equal to 0 or 1."},
	{-9999, "Could not build an Arnoldi factorization. User input error highly likely. Please check actual array dimensions and layout. IPARAM(5) returns the size of the current Arnoldi factorization."}
};

int arpack_errors_length = 16; 



ArpackSIPardisoSolver::ArpackSIPardisoSolver(unsigned int size_, unsigned int block_size)
: EigenSolver<cplx,double,cplx>(size_),
  evecs(0),
  evals(0),
  solver(0),
  overlap(0)
 {
	
	TDKP_ASSERT(size_ % block_size == 0, "size must be #nodes * block_size!");
	
	this->solver       = LinearSolver<cplx>::factory(size_);
	if(solver->get_matrix().property_is_set(symmetric_matrix)) {
		this->overlap = new CSRMatrix<double>(size_ / block_size, symmetric_matrix);
	} else {
		this->overlap = new CSRMatrix<double>(size_ / block_size, nonsymmetric_matrix);		
	}
	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "ArpackSIPardisoSolver: is used as the eigensolver");
	
}

ArpackSIPardisoSolver::~ArpackSIPardisoSolver() {
	if(this->evals != 0) {
		delete[] this->evals; this->evals = 0;
	}
	if(this->evecs != 0) {
		delete[] this->evecs; this->evecs = 0;
	}
	if(this->solver) {
  		delete solver; solver = 0;
	}
	if(this->overlap) {
		delete overlap; overlap = 0;	
	}
}


cplx ArpackSIPardisoSolver::eigenvalue(int nn) const throw(Exception*) {
	if(this->evals == 0) {
		TDKP_GENERAL_EXCEPTION("no eigenvalues available");
	}
	TDKP_BOUNDS_ASSERT(nn < this->nev && nn >= 0, "nn < this->nev && nn >= 0");
	return this->evals[nn];
}



const cplx& ArpackSIPardisoSolver::eigenvector(int nn, int vv) const throw(Exception*) {
	if(this->evecs == 0) {
		TDKP_GENERAL_EXCEPTION("no eigenvectors available");
	}
	TDKP_BOUNDS_ASSERT(nn < this->nev && nn >= 0 && vv >= 0 && vv <= this->msize, "nn < this->nev && nn >= 0 && vv >= 0 && v <= msize");
	return this->evecs[nn * this->msize + vv];
}



bool ArpackSIPardisoSolver::find_eigenvectors() throw(Exception*) {

        const Configuration* arpack_controls = Configuration::get_instance();

        if(this->nev == 0 || this->overlap == 0) {
            TDKP_GENERAL_EXCEPTION("eigenvalueproblem not properly assigned");
        }
        TimeMeasurements::get_instance().start("matrix factorization");	     
        solver->prepare();
        TimeMeasurements::get_instance().stop("matrix factorization");
        TimeMeasurements::get_instance().track_memory_usage();

		// does the user want to check the residual during solving of Ax = b? 
		bool test_solver_residual = false;
		if(Configuration::get_instance()->get("arpack_check_solver_residuum") == 1.0) {
			test_solver_residual = true;	
		}

        // my nice progress info output stuff
        const char progress_markers[] = {'-','\\','|','/'};
        int counter = 0;
        int sign    = 0;
        int next    = 0;
        ostringstream progressometer;
        // ----------------------------------------------------------
        // arpack stuff
        // ----------------------------------------------------------
        int  ido        = 0;
        char bmat[1]    = {'G'};
        int  n          = this->get_stiff().get_size();
        char which[2]   = {'L', 'M'};
        int  tmp_nev    = this->nev > 1 ? this->nev:2;
        double tol      = arpack_controls->get("arpack_tolerance"); // 1.0e-10
        cplx* resid     = new cplx[n];
        
        int  ncv        = (int)arpack_controls->get("arpack_min_arnoldi_vecs"); // 36;
        
        // ---------------------------------------------
        // adjust number of arnoldi vecs if needed
        // ---------------------------------------------
        if(ncv < tmp_nev + static_cast<int>(arpack_controls->get("arpack_min_arnoldi_oversize"))) {
        	ncv = tmp_nev + static_cast<int>(arpack_controls->get("arpack_min_arnoldi_oversize"));
        }
        
        cplx* v         = new cplx[n * ncv];
        int  ldv        = n;
        int iparam[11];
        iparam[0]       = 1; // use exact shifts with respect to hessenberg matrix ...
        iparam[2]       = (int)(arpack_controls->get("arpack_max_arnoldi_iterations")) * tmp_nev;  // 150 number of arnoldi update iterations allowed
        iparam[3]       = 1;
        iparam[6]       = 3; // inv(A - sigma B)x = lambda x calculation
        int ipntr[14];
        cplx* workd     = new cplx[3*n];
        int lworkl      = 3*ncv*ncv + 6 * ncv;
        cplx* workl     = new cplx[lworkl];
        double* rwork   = new double[ncv];
        int info        = 0;
                
        // ---------------------------------------------------------
        // OPx -> w arrays
        // ---------------------------------------------------------
        cplx*   multmp     = new cplx[n];
        cplx*   in         = 0;
        cplx*   out        = 0;


        // number of partial differential equations
        int kpn          = (signed)(this->get_stiff().get_size() / this->overlap->get_size());

        cstatn_();
		double z1, z2, z3, z4, z5, z6;
		z1 = z2 = z3 = z4 = z5 = z6 = 0.0e0;
		double tm1, tm2;
		double adaptive_omp_tick = 0.0;
		

				
        while(ido != 99) {
        	tm2 = tic();
        	if(tm2 - adaptive_omp_tick > 3.0) {
//				adaptive_omp_threading();
				adaptive_omp_tick = tm2;
        	}
			tm1 = tic();
			z2 += tm1 - tm2;
            znaupd_(&ido, bmat, &n, which, &tmp_nev, &tol, &resid[0], &ncv,
                    v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
			tm2 = tic();
			z1 += tm2 - tm1;
            // ---------------------------------------------------------
            // reverse communication
            // ---------------------------------------------------------
            if(ido == -1 || ido == 1) {
				
                out = &workd[ipntr[1] - 1];

                // compute Bx -> d
                if(ido == -1) {
                    in = &workd[ipntr[0] - 1];
                /*  if(kpn == 1) {
                        this->overlap->mult_vec(in,multmp);
                    } else { */
                        this->overlap->mult_vec_multiple(in, multmp, kpn);
                    //}
                    in = multmp;
                } else {
                    // Bx -> d already available in workd(ipntr(2))
                    in = &workd[ipntr[2] - 1];
                }
                tm1 = tic();
                z3 += tm1 - tm2;
                // Bx -> d is now in "in"
                // solve now for Aw = in
                solver->solve_equation(out, in, 1);
                tm2 = tic();
                z4 += tm2 - tm1; 
                // --------------------------------------------------
                // test solver
                // --------------------------------------------------
                if(test_solver_residual) {
                	vector<cplx> tst(n,0);
                	solver->get_matrix().mult_vec(out, &tst[0]);
                	double solver_error = 0.0;
                	cplx d;
                	for(int ii = 0; ii < n; ii++) {
                		d = tst[ii] - in[ii];
                		d *= conj(d);
                		solver_error += d.real();	
                	}
                	ostringstream sout;
                	sout << "\r" << setw(3) << timing_.nopx << " -> |Ax - b| = " << sqrt(solver_error) 
                	     << "                                          \n";
					Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());                	     
                }
            } else if(ido == 2) {
                // multiply Bx -> d
                out = &workd[ipntr[1] - 1];
                in = &workd[ipntr[0] - 1];
            /*  if(kpn == 1) {
                    this->overlap->mult_vec(in, out);
                } else { */
                    this->overlap->mult_vec_multiple(in, out, kpn);
                //}
                tm1 = tic();
                z5 += tm1 - tm2;
                tm2 = tm1;
            } else if(ido == 3) {
                TDKP_GENERAL_EXCEPTION("ido wants me to compute shifts.");
            } else if(ido == 99) {
                break;
            } else {
                TDKP_GENERAL_EXCEPTION("ido wants me to do something strange. i don't like ido");
            }
            
            // -------------------------------------------------
            // fancy progress output
            // -------------------------------------------------                       
            if(Logger::get_instance()->get_level() > LOG_INFO) {
	            if(counter++ == next) {
	                progressometer.str("");
	                progressometer << " " << progress_markers[(sign++) % 4]             
	                               << " Bx -> w: "  << std::setw(5) << timing_.nbx 
	                               << " OPx -> w: " << std::setw(5) << timing_.nopx  
	                               << " reosteps: " << std::setw(5) << timing_.nrorth 
	                               << " restarts: " << std::setw(5) << timing_.nrstrt;
	                Logger::get_instance()->set_line_to_stdout(progressometer.str());
	                next = counter + 5;             
	            }	            
            }   
            tm1 = tic();
            z6 += tm1 - tm2;
        }
        if(Logger::get_instance()->get_level() > LOG_INFO) {
        	progressometer << "\n";
        	Logger::get_instance()->set_line_to_stdout(progressometer.str());
        }
		TimeMeasurements::get_instance().track_memory_usage();
		solver->release();

        if(this->evals != 0) {
        	delete[] this->evals;
			this->evals = 0;        			
        }
        if(this->evecs != 0) {
        	delete[] this->evecs;
        	this->evecs = 0;	
        }
        
        ostringstream tout;
        tout << "z1 (znaupd_)       = " << z1 << "\n"
             << "z2 (adaptive_omp)  = " << z2 << "\n"
             << "z3 (Bx -> d/ido 1) = " << z3 << "\n"
             << "z4 (solve eq)      = " << z4 << "\n"
             << "z5 (Bx -> d/ido 2) = " << z5 << "\n"
             << "z6 (progr.meter)   = " << z6 << "\n"
             << "total zi           = " << z1 + z2 + z3 + z4 + z5 + z6 << "\n";  
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, tout.str());             

        int   rvec = 1;
        char  howmny[] = {'A'};
        int*  select = new int[ncv];
        this->evals = new cplx[tmp_nev + 1];
        this->evecs = new cplx[n * tmp_nev];
        int   ldz    = n;
        cplx  sigma(0,0);
        cplx* workev = new cplx[2 * ncv];
        std::ostringstream sout;

        switch(info) {
            case 0:
                sout << "the eigensolver found " << iparam[4] << " of " << this->nev << " eigenvalues";
                if(iparam[4] != tmp_nev) {
                    TDKP_GENERAL_EXCEPTION(sout.str());
                } else if(Logger::get_instance()->get_level() > LOG_INFO) {
                    Logger::get_instance()->emit(LOG_INFO_DEVEL1, sout.str());
                }
                // get ritz eigenvalues/eigenvectors
                zneupd_(&rvec, howmny, select, this->evals, this->evecs, &ldz, &sigma, workev,
                         bmat, &n, which, &tmp_nev, &tol, &resid[0], &ncv,
                        v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

                this->converged_ev = iparam[4] > this->nev ? this->nev : iparam[4];
                // sort solutions
                this->sort_solutions(tmp_nev, this->evecs, this->evals);
                this->orthogonalize_solutions();
                break;

            case 1:
                TDKP_GENERAL_EXCEPTION("eigensolver quit after reaching max iterations");
                break;
            case 3:
                TDKP_GENERAL_EXCEPTION("eigensolver quit because he could not appliy shifts");
                break;
            default:
                int err = -1;
                for(int ii = 0; ii < arpack_errors_length; ii++) {
                    if(arpack_errors[ii].code == info) {
                        err = ii;
                        break;
                    }
                }
                sout.str("");
                sout << "eigensolver arpack returned error code " << info;
                if(err >= 0) {              
                    sout << ": " << arpack_errors[err].text;
                }         
                TDKP_GENERAL_EXCEPTION(sout.str());
                break;
        }

        // -------------------------------------------------------------
        // print stats
        // -------------------------------------------------------------
        if(Configuration::get_instance()->get("output_eigensolver_statistics") == 1) {
            sout.str("");
            sout << timing_.nopx   << std::setw(8) << " = total number of user OP*x operation\n"
                 << timing_.nbx    << std::setw(8) << " = total number of user B*x operation (same as copy when B = I)\n"
                 << timing_.nrorth << std::setw(8) << " = total number of reorthogonalization steps taken\n"
                 << timing_.nitref << std::setw(8) << " = total number of it. refinement steps in reorthogonalization\n"
                 << timing_.nrstrt << std::setw(8) << " = total number of restart steps\n"
                 << timing_.tcaupd << std::setw(8) << " = total time spent in CNAUPD.\n"
                 << timing_.tcaup2 << std::setw(8) << " = total time spent in CNAUP2.\n"
                 << timing_.tcaitr << std::setw(8) << " = total time spent in the basic Arnoldi iteration loop, including iterative refinement in CNAITR.\n"
                 << timing_.tceigh << std::setw(8) << " = total time spent in computing the Hessenberg eigenvalue subproblem at each iteration.\n"
                 << timing_.tcgets << std::setw(8) << " = total time spent in getting the shifts at each iteration.\n"
                 << timing_.tcapps << std::setw(8) << " = total time spent in applying the shifts at each iteration.\n"
                 << timing_.tcconv << std::setw(8) << " = total time spent in convergence test at each iteration.\n"
                 << timing_.tmvopx << std::setw(8) << " = total time spent in computing Y = OP * X\n"
                 << timing_.tmvbx  << std::setw(8) << " = total time spent in computing Y = B * X\n";
            Logger::get_instance()->emit(LOG_INFO, sout.str());
        }
        
    
		if(Configuration::get_instance()->get("arpack_check_eigenvectors_m_orthogonal") == 1.0) {
	        Logger::get_instance()->emit(LOG_INFO_DEVEL2, "testing if dot products of eigenvectors are M orthogonal");  
	        // --------------------------------------------------------
	        // testing orthogonality of eigenvectors
	        // (arpack did not orthogonalize degenerate basis vectors
	        // --------------------------------------------------------
	        vector<cplx> ortho_in(this->get_stiff().get_size());
	        vector<cplx> ortho_out(this->get_stiff().get_size());
	        for(int ii = 0; ii < this->converged_eigenvalues(); ii++) {
	            for(int nn = 0; nn < this->converged_eigenvalues(); nn++) {
	                // calculate  M orthogonal matrix product
	                cplx tmp(0,0);          
	                for(int jj = 0; jj < (signed)this->get_stiff().get_size(); jj++) {
	                    ortho_in[jj] = this->eigenvector(nn,jj);    
	                }
	                this->overlap->mult_vec_multiple(&ortho_in[0], &ortho_out[0], kpn);     
	                for(int jj = 0; jj < (signed)this->get_stiff().get_size(); jj++) {
	                    tmp += conj(this->eigenvector(ii,jj)) * ortho_out[jj];           
	                }
	                if(nn != ii) {
	                    if(abs(tmp) > 1.0e-12) {
	                        ostringstream sout;
	                        sout << "scalar product between " << ii << " and " << nn << " is nonzero: " << tmp;
	                        Logger::get_instance()->emit(LOG_WARN, sout.str());
	                    }
	                }
	            }
	        }   
		}              
        
        delete[] select;
        delete[] workev;

        delete[] resid;
        delete[] v;
        delete[] workd;
        delete[] workl;
        delete[] rwork;
        delete[] multmp;
               
        return true;
}

/** reorthogonalize solutions
 * 
 * there seems to be one problem with arpack: eigenvectors
 * to the same eigenvalue are not orthogonal!
 */
void ArpackSIPardisoSolver::orthogonalize_solutions() {

    int numeq        = this->get_stiff().get_size() / this->overlap->get_size();
    int size         = (signed)this->get_stiff().get_size();
    cplx* nbasis     = new cplx[size * this->nev];
    cplx* tmp        = new cplx[size];
    double threshold = Configuration::get_instance()->get("eigenvector_reorthogonalization_threshold"); 
    cplx norm        = 0;
    
    // some data maybe used to show the user the result
    // of the reorthogonalization procedure
    vector<double> norm_matrix(this->nev * this->nev, -1.0);
    
    for(int nn = 0; nn < this->nev; nn++) {
        // copy new basis vector
        for(int ii = 0; ii < size; ii++) {
            nbasis[nn * size + ii] = this->evecs[nn * size + ii];
        }
        // orthogonalize against existing
        for(int mm = 0; mm < nn; mm++) {
            // Mxj
            this->overlap->mult_vec_multiple(&(nbasis[mm * size]),tmp,numeq);

            norm = 0.0;
            for(int ii = 0; ii < size; ii++) {
                norm += conj(this->evecs[nn * size + ii]) * tmp[ii];
            }
            norm_matrix[nn * this->nev + mm] = abs(norm);
            if(abs(norm) > threshold) {         
                // xi - xj conj((xiMxj))
                for(int ii = 0; ii < size; ii++) {
                    nbasis[nn * size + ii] -= nbasis[mm * size + ii] * conj(norm);
                }
            }
        }
        // orthogonalize resulting basis
        norm = 0.0;
        this->overlap->mult_vec_multiple(&(nbasis[nn * size]),tmp,numeq);

        for(int ii = 0; ii < size; ii++) {
            norm += conj(nbasis[nn * size + ii]) * tmp[ii];
        }
        norm = sqrt(norm);
        // normalize
	    for(int ii = 0; ii < size; ii++) {
	    	nbasis[nn * size + ii] /= norm;
	   	}
       
    }
    if(Configuration::get_instance()->get("eigenvector_skip_rayleigh_quotient_testing") == 0.0) {
	    cplx* tmp2 = new cplx[size];
	    cplx a,b;
	    // test rayleigh quotient
	    for(int nn = 0; nn < this->nev; nn++) {
	        // Mx -> tmp
	    /*  if(numeq == 1) {
	            this->overlap->mult_vec(&(nbasis[nn * size]),tmp);
	        } else {*/
	            // must use mult_vec_multiple as due to the templated abstract interface
	            // i could not use virtual templated functions. therefore i had 
	            // to disallow the usage of arbitrary types in mult_vec and restrict
	            // it to the same type as the matrix is. 
	            // but as overlap is a real valued matrix, but the vectors are complex,
	            // in need this one here ...
	            this->overlap->mult_vec_multiple(&(nbasis[nn * size]),tmp,numeq);
	        //}
	        // Ax -> tmp2
	        this->get_stiff().mult_vec(&(nbasis[nn * size]),tmp2);
	        // x'Ax/x'Mx
	        a = b = 0.0;
	        for(int ii = 0; ii < size; ii++) {
	            a += conj(nbasis[nn * size + ii]) * tmp2[ii];
	            b += conj(nbasis[nn * size + ii]) * tmp[ii];
	        }
	        if(abs(a/b - this->evals[nn]) > 1.0e-10) {
	        	double error = abs(a/b - this->evals[nn]);
	            ostringstream eout;
	            eout << "problem with eigenvector " << nn 
	                 << ": rayleigh quotient after orthogonalization (=" << a/b 
	                 << ") differs from calculated eigenvalue: " << this->evals[nn]
	                 << " by " << error;
	            Logger::get_instance()->emit(LOG_WARN, eout.str());
	            if(error > 1.0e-7) {
	            	TDKP_GENERAL_EXCEPTION(eout.str());
	            }
	        }
	    }
	    delete[] tmp2;
    } else {
    	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "skipping rayleigh quotient testing");	
    }

    // copy back
    for(int ii = 0; ii < size * this->nev; ii++) {
        this->evecs[ii] = nbasis[ii];
    }
    delete[] tmp;
    delete[] nbasis;
    
    // ---------------------------------------------------------
    // write output on request
    // ---------------------------------------------------------
    if(Configuration::get_instance()->get("output_arpack_reorthogonalization") == 1.0) {
        ostringstream sout;
        sout << "The eigensolution obtained from arpack produced the following\n"
             << "norms during the Gram-Schmidt orthogonalization procedure:\n";
        sout.precision(3);
        for(int nn = 0; nn < this->nev; nn++) {
            for(int mm = 0; mm < this->nev; mm++) {
                sout << setw(9);
                if(nn == mm ) {
                    sout << "(self)" << " ";    
                } else if(nn < mm) {
                    sout << " ";
                } else {
                    if(norm_matrix[nn * this->nev + mm] != -1.0) {
                        sout << norm_matrix[nn * this->nev + mm];
                    } else {
                        sout << " ";
                    }
                    sout << "  "; 
                }
            }
            sout << "\n";
        }
        Logger::get_instance()->emit(LOG_INFO, sout.str());
    }
}


ArpackSIRealWrapper::ArpackSIRealWrapper(unsigned int size_, unsigned int block_size_) 
: EigenSolver<double,double,cplx>(size_), // dummy 
  solver(size_,block_size_),
  matrix_wrapper(solver.get_stiff())
{
		
}

ArpackSIRealWrapper::~ArpackSIRealWrapper() {
	
}
void ArpackSIRealWrapper::assign(int nev, EigenProblemType type) {
	solver.assign(nev, type);	
}
void ArpackSIRealWrapper::set_ordering(Ordering order_) {
	solver.set_ordering(order_);	
}
Ordering ArpackSIRealWrapper::get_ordering() const {
	return solver.get_ordering();	
}
bool ArpackSIRealWrapper::find_eigenvectors() throw(Exception*) {
	return solver.find_eigenvectors();	
}
int ArpackSIRealWrapper::converged_eigenvalues() const {
	return solver.converged_eigenvalues();	
}
cplx ArpackSIRealWrapper::eigenvalue(int nn) const throw(Exception*) {
	return solver.eigenvalue(nn);	
}
const cplx& ArpackSIRealWrapper::eigenvector(int nn, int vv) const throw(Exception*) {
	return solver.eigenvector(nn,vv);	
}
	
			
// -----------------------------------------------
// GLOBAL IFDEF FOR EXCLUDING ARPACK
// -----------------------------------------------
#endif
