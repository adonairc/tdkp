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

#include "tdkp/main/CSRMatrix.h"
#include <sstream> 

namespace tdkp
{

template<> cplx CSRMatrix<cplx>::internal_conj(cplx z) const {
	return conj(z);
}
template<> double CSRMatrix<double>::internal_conj(double d) const {
	return d;	
}

template<> bool CSRMatrix<cplx>::symmetric_by_value() const {

	std::ostringstream sout;
	double tolerance = 1.0e-14;
	sout << "checking for hermiticity with tolerance " << tolerance;	
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, sout.str());
	
	if(this->is_symmetric) return true;
	if(this->size < 15000) {
		for(int ii = 0; ii < (signed)this->size; ii++) {			
			for(int jj = ii; jj < (signed)this->size; jj++) {			
				if(tdkp_math::abs(this->get(ii,jj) - conj(this->get(jj,ii))) > tolerance) {
					sout << " difference found at: (" << ii << "," << jj <<") - "
					     << this->get(ii,jj) << " bzw. " << this->get(jj,ii);
					Logger::get_instance()->emit(LOG_INFO, sout.str());
					return false;	
				}
			}
		}	
		return true;	
	} else {
		Logger::get_instance()->emit(LOG_WARN, "not checking whether matrix is symmetric");
		return true;
	}	
}

template<> template<>
bool CSRMatrix<cplx>::symmetry_value_check<cplx>(cplx a, cplx b, double tol) const {
	double rel = (abs(b) > 0.0 ? abs(b) : 1.0);
	if(abs(conj(a) - b) / rel > tol) {
		// just ignore too small values
		if(abs(a) < tol && abs(b) < tol) {
			return true;	
		}
		return false;	
	}
	return true;	
}

template<> template<>
bool CSRMatrix<double>::symmetry_value_check<double>(double a, double b, double tol) const {
	double rel = (b != 0.0 ? b:1.0);
	if(fabs(a - b) / rel > tol) {
		return false;	
	} else {
		return true;	
	}	
}

template<> 
double CSRMatrix<cplx>::symmetry_value_error(cplx a, cplx b) const {
	return abs(conj(a) - b);
}

template<> 
double CSRMatrix<double>::symmetry_value_error(double a, double b) const {
	return fabs(a - b);
}


template<> 
void CSRMatrix<cplx>::mult_vec(const cplx* in, cplx* out) const {

	if(this->prow && this->icol	&& this->nonzeros) {
		int until, from;	
		// --------------------------------------------------
		// distinguish between symmetric and nonsymmetric
		// --------------------------------------------------		
		if(this->is_symmetric) {					
			// clean up out
			for(int ii = 0; ii < (signed)this->size; ii++) {
				out[ii] = 0.0;	
			}							
			// perform symmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - fidx;
				from     = this->prow[ii] - fidx;
				out[ii] += this->nonzeros[from] * in[ii]; // diagonal	
//				TDKP_BOUNDS_ASSERT(from >= 0 && from < this->num_nz && until >= 0 && until <= this->num_nz, "from >= 0 && from < this->size && until >= 0 && until <= this->size");
				// hermitian multiplication
				for(int jj = from + 1; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - fidx)];
					out[(this->icol[jj] - fidx)] += conj(this->nonzeros[jj]) * in[ii];
				}		
			}			
		} else {
			// perform nonsymmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - this->fidx;
				from     = this->prow[ii] - this->fidx;
				out[ii]  = 0.0;			
				for(int jj = from; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - this->fidx)];
				}		
			}			
		}		
	} else {
		TDKP_GENERAL_EXCEPTION("can not multiply with noninitialized matrix");	
	}		
}

template<> template<>
void CSRMatrix<cplx>::mult_vec_multiple(const cplx* in, cplx* out, int kpn) const {
	
	int msize    = (signed)this->size;
	int vec_size = msize * kpn;

	if(this->prow && this->icol	&& this->nonzeros) {
		int until, from, off;
		// --------------------------------------------------
		// distinguish between symmetric and nonsymmetric
		// --------------------------------------------------
		if(this->is_symmetric) {
			// clean up out
			for(int ii = 0; ii < vec_size; ii++) {
				out[ii] = 0;
			}
			// perform symmetric matrix vector product
			for(int ii = 0; ii < msize; ii++) {
				until    = this->prow[ii + 1] - fidx;
				from     = this->prow[ii] - fidx;
				off      = ii * kpn;
				for(int nn = 0; nn < kpn; nn++) {
					out[off + nn] += this->nonzeros[from] * in[off + nn]; // diagonal
				}
				TDKP_BOUNDS_ASSERT(from >= 0 && from < this->num_nz && until >= 0 && until <= this->num_nz, "from >= 0 && from < this->size && until >= 0 && until <= this->size");
				for(int jj = from + 1; jj < until; jj++) {
					for(int nn = 0; nn < kpn; nn++) {
						out[off + nn] += this->nonzeros[jj] * in[(this->icol[jj] - fidx) * kpn + nn];
						out[(this->icol[jj] - fidx) * kpn + nn] += conj(this->nonzeros[jj]) * in[off + nn];
					}
				}
			}
		} else {
			cplx ss[MAX_NUM_EQUATIONS];
			TDKP_ASSERT(kpn <= MAX_NUM_EQUATIONS, "kpn <= MAX_NUM_EQUATIONS");

			// perform nonsymmetric matrix vector product
			int* copy_prow		  = this->prow;
			int* copy_icol		  = this->icol;
			cplx*   copy_nonzeros = this->nonzeros;
#pragma omp parallel for default(none) private(until,from,ss) shared(copy_prow,copy_icol,copy_nonzeros,out,kpn,in) schedule(static,5000)
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = copy_prow[ii + 1] - this->fidx;
				from     = copy_prow[ii] - this->fidx;
				for(int nn = 0; nn < kpn; nn++) {
					ss[nn] = 0.0;
//					out[ii * kpn + nn] = 0.0;
				}
				for(int jj = from; jj < until; jj++) {
					for(int nn = 0; nn < kpn; nn++) {
						ss[nn] += copy_nonzeros[jj] * in[(copy_icol[jj] - this->fidx) * kpn + nn];
//						out[ii * kpn + nn] += this->nonzeros[jj] * in[(this->icol[jj] - this->fidx) * kpn + nn];
					}
				}
				for(int nn = 0; nn < kpn; nn++) {
					out[ii * kpn + nn]  = ss[nn];
				}
			}
		}
	} else {
		TDKP_GENERAL_EXCEPTION("can not multiply with noninitialized matrix");
	}
}

template<> template<>
void CSRMatrix<double>::mult_vec<cplx>(const cplx* in, cplx* out) const {

	if(this->prow && this->icol	&& this->nonzeros) {
		int until, from;
		// --------------------------------------------------
		// distinguish between symmetric and nonsymmetric
		// --------------------------------------------------
		if(this->is_symmetric) {
			// clean up out
			for(int ii = 0; ii < (signed)this->size; ii++) {
				out[ii] = 0.0;
			}
			// perform symmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - fidx;
				from     = this->prow[ii] - fidx;
				out[ii] += this->nonzeros[from] * in[ii]; // diagonal
		//		TDKP_BOUNDS_ASSERT(from >= 0 && from < this->num_nz && until >= 0 && until <= this->num_nz, "from >= 0 && from < this->size && until >= 0 && until <= this->size");
				for(int jj = from + 1; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - fidx)];
					out[(this->icol[jj] - fidx)] += this->nonzeros[jj] * in[ii];
				}
			}
		} else {
			// perform nonsymmetric matrix vector product
			for(int ii = 0; ii < (signed)this->size; ii++) {
				until    = this->prow[ii + 1] - this->fidx;
				from     = this->prow[ii] - this->fidx;
				out[ii]  = 0.0;
				for(int jj = from; jj < until; jj++) {
					out[ii] += this->nonzeros[jj] * in[(this->icol[jj] - this->fidx)];
				}
			}
		}
	} else {
		TDKP_GENERAL_EXCEPTION("can not multiply with noninitialized matrix");
	}
}


} // end namespace
