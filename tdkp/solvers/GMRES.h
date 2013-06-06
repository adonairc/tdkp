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

#ifndef GMRES_H_
#define GMRES_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Vector3D.h"

namespace tdkp {

/** gmres solver class used for microscopic polarization 
 *
 * the coulomb interaction seems to be a small perturbation and therefore
 * calculating the new lambda is quite fast when the iterative method is used
 * over the full matrix inversion.
 * 
 * therefore i simply fetched the gmres version from netlib and quickly adapted it 
 */
class GMRES {
public:
	
	void solve(int nrestart, const RMatrix<cplx>& A,  const vector<cplx>& d, vector<cplx>& x);
		
	class GMRESVector  {
	public:
		GMRESVector();
		explicit GMRESVector(unsigned int n_);
		GMRESVector(const GMRESVector& copy);
		explicit GMRESVector(vector<cplx>& target);
		explicit GMRESVector(const vector<cplx>& source);
		~GMRESVector();
		cplx& operator()(int ii) { return val[ii]; }
		const cplx& operator()(int ii) const { return val[ii]; }
		int size() const { return n; }
		cplx* get_data() { return val; }
		const cplx* get_data() const { return val; }
		void operator*=(const double& rhs);
		void operator*=(const cplx& rhs);
		void operator+=(const GMRESVector& rhs);
		void operator=(const GMRESVector& rhs);
		void assign(const double& x); 
	private:
		int   n;
		bool  my_array;		
		cplx* val;
	};
	
		
private:
	
};

} // end of namespace tdkp

#endif /*GMRES_H_*/
