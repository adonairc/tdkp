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


extern "C" {

/** external jdqz fortran routine */
void F77_Name(jdqz)(
	complex<double>* alpha, complex<double>* beta, complex<double>* eivec, int* wanted, 
	int* n, complex<double>* shift, double* eps, int* kmax, int* jmax, int*  jmin,
    int* method, int*  m, int* l, int* maxnmv, int* maxstep, double* lock, int*  order, 
    int* testspace, complex<double>* work, int* lwork
);
	
}

namespace tdkp {
	
template<class RHS>	
JDQZSolver<RHS>::JDQZSolver(unsigned int size_, unsigned int block_size_)
: EigenSolver<cplx,RHS,cplx>(size_),
  evecs(0),
  evals(0),
  solver(0),
  overlap(0),
  block_size(block_size_)
 {			
	TDKP_ASSERT(size_ % block_size == 0, "size must be #nodes * block_size!");
	TDKP_ASSERT(block_size > 0, "");	
	this->solver = LinearSolver<cplx>::factory(size_);
	if(solver->get_matrix().property_is_set(symmetric_matrix)) {
		this->overlap = new CSRMatrix<RHS>(size_ / block_size, symmetric_matrix);
	} else {
		this->overlap = new CSRMatrix<RHS>(size_ / block_size, nonsymmetric_matrix);		
	}	
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "JDQZSolver: is used as the eigensolver");	
}

template<class RHS>
JDQZSolver<RHS>::~JDQZSolver() {
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

template<>
void JDQZSolver<cplx>::aquire();
template<>
void JDQZSolver<cplx>::release();

template<>
void JDQZSolver<double>::aquire();
template<>
void JDQZSolver<double>::release();

template<class RHS>
cplx JDQZSolver<RHS>::eigenvalue(int nn) const throw(Exception*) {
	if(this->evals == 0) {
		TDKP_GENERAL_EXCEPTION("no eigenvalues available");
	}
	TDKP_BOUNDS_ASSERT(nn < this->nev && nn >= 0, "nn < this->nev && nn >= 0");
	return this->evals[nn];
}

template<class RHS>
const cplx& JDQZSolver<RHS>::eigenvector(int nn, int vv) const throw(Exception*) {
	if(this->evecs == 0) {
		TDKP_GENERAL_EXCEPTION("no eigenvectors available");
	}
	TDKP_BOUNDS_ASSERT(nn < this->nev && nn >= 0 && vv >= 0 && vv <= this->msize, "nn < this->nev && nn >= 0 && vv >= 0 && v <= msize");
	return this->evecs[nn * this->msize + vv];
}

template<class RHS>
bool JDQZSolver<RHS>::find_eigenvectors() throw(Exception*) {

	// --------------------------------------------------
	// jdqz has a search space limit of 50. we stop here if
	// the user requests more than 49 nev's ...
	// --------------------------------------------------
	TDKP_ASSERT(this->nev <= 49, "JDQZSolver: sorry, you want me to calculate " << this->nev << " eigenvalues but jdqz has a fixed limit of 49 eigenvalues (hard coded).");

	// --------------------------------------------------
	// clean up old stuff
	// --------------------------------------------------
    if(this->evals != 0) {
    	delete[] this->evals;
		this->evals = 0;        			
    }
    if(this->evecs != 0) {
    	delete[] this->evecs;
    	this->evecs = 0;	
    }
    
	// --------------------------------------------------
	// aquire global control (for multiplication with A and B etc.)
	// --------------------------------------------------
	this->aquire();			        
    if(this->nev == 0 || this->overlap == 0) {
        TDKP_GENERAL_EXCEPTION("eigenvalueproblem not properly assigned");
    }
    // --------------------------------------------------
    // factorize matrix
    // --------------------------------------------------
    TimeMeasurements::get_instance().track_memory_usage();
    TimeMeasurements::get_instance().start("matrix factorization");	     
    solver->prepare();
    TimeMeasurements::get_instance().stop("matrix factorization");
    TimeMeasurements::get_instance().track_memory_usage();

	// --------------------------------------------------
	// aquire storage space
	// --------------------------------------------------		
	int n       = get_stiff().get_size();
	int kmax    = this->nev + 1;   // num wanted eigensolutions
	int jmax    = max(60, kmax * 3); // max search space
	int jmin    = kmax * 2;  // min search space
	if(jmax > n) {
		jmax = n / 2;
		jmin = n / 4;
		kmax = this->nev + 1;	
	}
	// hard coded max in jmax ...
	if(jmax > 50) {
		jmax = 50;
		jmin = 25;	
	}
	int wanted  = 1;         // yes, i'd like to get 
	int method  = Configuration::get_instance()->get("jdqz_linsolv_method");		
	double eps  = Configuration::get_instance()->get("jdqz_tolerance");    // accuracy
	cplx target = 0.0;       // target (matrixes are shifted to 0 anyway)
	int m       = Configuration::get_instance()->get("jdqz_gmres_search_space_m");        // GMRES search space
	int l       = Configuration::get_instance()->get("jdqz_bicgstab_polynomial");         // BiCGstab GMRES polynomial
	int maxnmv  = Configuration::get_instance()->get("jdqz_linsolv_max_num_matvec");       // max num matrix vector multiplications in gmres/bicgstab 
	int maxstep = Configuration::get_instance()->get("jdqz_max_jd_iterations");      // max num JD iterations
	double lock = Configuration::get_instance()->get("jdqz_eigenvalue_lock");   // eigenvalue tracking
	int order   = Configuration::get_instance()->get("jdqz_eigenvalue_order");      // nearest to target
	int tstspc  = Configuration::get_instance()->get("jdqz_testspace_expansion_method");         // how to expand testspace
	int lwork   = method == 1 ? (4 + 5 + 5*jmax + 3*kmax) : (10 + 6*l + 5*jmax + 3*kmax);
	lwork *= 2;
		
			
	vector<cplx> alpha(jmax);
	vector<cplx> beta(jmax);
	vector<cplx> eivec(kmax * n);
	vector<cplx> work(n * lwork);
	
	F77_Name(jdqz)(
		&alpha[0], &beta[0], &eivec[0], &wanted, &n, &target, &eps, 
		&kmax, &jmax, &jmin, &method, &m, &l, &maxnmv, &maxstep, 
		  &lock, &order, &tstspc, &work[0], &lwork
	);
			     
	TimeMeasurements::get_instance().track_memory_usage();
	solver->release();
          
	// -------------------------------------------------------
	// copy eigenvalues / vectors
	// -------------------------------------------------------
	this->evals = new cplx[kmax];
	this->evecs = new cplx[kmax * get_stiff().get_size()];
	TDKP_ASSERT((signed)get_stiff().get_size() == n, "");
    		
	for(int ii = 0; ii < kmax; ii++) {
		TDKP_ASSERT(tdkp_math::abs(beta[ii]) != 0.0, "JDQZ could not properly solve eigenvalue problem. beta[" << ii << "] is zero!"); 			
		evals[ii] = alpha[ii] / beta[ii];
		for(int jj = 0; jj < n; jj++) {
			evecs[ii * n + jj] = eivec[ii * n + jj];				
		}		
	}
	int tmp_nev = this->nev;
	this->nev = kmax;
	this->sort_solutions(kmax, evecs, evals);
	this->converged_ev = kmax;
	this->orthogonalize_solutions();
	this->converged_ev = this->nev = tmp_nev;              
          
	if(Configuration::get_instance()->get("arpack_check_eigenvectors_m_orthogonal") == 1.0) {
	
        Logger::get_instance()->emit(LOG_INFO_DEVEL2, "testing if dot products of eigenvectors are M orthogonal");  
        // --------------------------------------------------------
        // testing orthogonality of eigenvectors
        // (arpack did not orthogonalize degenerate basis vectors)
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
                this->overlap->mult_vec_multiple(&ortho_in[0], &ortho_out[0], this->get_block_size());     
                for(int jj = 0; jj < (signed)this->get_stiff().get_size(); jj++) {
                    tmp += conj(this->eigenvector(ii,jj)) * ortho_out[jj];           
                }	             
                if(nn != ii) {
                    if(abs(tmp) > 1.0e-12 && sizeof(RHS) == 8) {
                        ostringstream sout;
                        sout << "JDQZSolver: scalar product between " << ii << " and " << nn << " is nonzero: " << tmp;
                        Logger::get_instance()->emit(LOG_WARN, sout.str());
                    }
                }
            }	    
        }   
	}              
                   
	// release global control (for multiplication with A and B etc.)
	this->release();               
           
    return true;
}

/** reorthogonalize solutions
 * 
 * there seems to be one problem with arpack: eigenvectors
 * to the same eigenvalue are not orthogonal!
 */
template<class RHS>
void JDQZSolver<RHS>::orthogonalize_solutions() {

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
            this->overlap->mult_vec_multiple(&(nbasis[nn * size]),tmp,numeq);
	        // Ax -> tmp2
	        this->get_stiff().mult_vec(&(nbasis[nn * size]),tmp2);
	        // x'Ax/x'Mx
	        a = b = 0.0;
	        for(int ii = 0; ii < size; ii++) {
	            a += conj(nbasis[nn * size + ii]) * tmp2[ii];
	            b += conj(nbasis[nn * size + ii]) * tmp[ii];
	        }
	        // -------------------------------------------------
	        // only complain for double (if problem is real)
	        // eigenvectors are NOT orthonormal if PML is present
	        // -------------------------------------------------
	        if(abs(a/b - this->evals[nn]) > 1.0e-10 && sizeof(RHS) == 8) {
	        	double error = abs(a/b - this->evals[nn]);
	            ostringstream eout;
	            eout << "JDQZSolver: problem with eigenvector " << nn 
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
        sout << "JDQZSolver: The eigensolution obtained from arpack produced the following\n"
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
	
	
} // end of namespace
