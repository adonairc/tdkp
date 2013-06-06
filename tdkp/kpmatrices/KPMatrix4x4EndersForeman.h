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

#ifndef KPMATRIX4X4ENDERSFOREMAN_H_
#define KPMATRIX4X4ENDERSFOREMAN_H_

#include "tdkp/kpmatrices/KPMatrixBase.h"

namespace tdkp
{

/** kp 4x4 matrix using |3/2,+/- 3/2> and |3/2, +/- 1/2> as basis functions
 * 
 * the kp matrix was derived by sebastian steiger using foremans modified 
 * kp8x8 enders hamiltonian (to use proper term ordering) and applying 
 * a basis transformation of the |X>, |Y>, |Z> states to the 4x4 states  
 * 
 * see:
 *   Phys. Rev. B, Vol. 51, Nr. 23, p. 16695 (Enders 1995)
 *   Phys. Rev. B, Vol. 56, Nr. 20, p. R 12748 (Foreman 1997)
 *   Phys. Rev. B, Vol. 48, Nr.  7, p. 4964 (Foreman 1993)
 * 
 */

class KPMatrix4x4EndersForeman : public tdkp::KPMatrixBase {
public:
	KPMatrix4x4EndersForeman();	
	KPMatrix4x4EndersForeman(const Material* mat);
	virtual ~KPMatrix4x4EndersForeman();
	const int* get_sparsity_pattern(int& length) const;	
	const int* get_diagonal_pattern(int& length) const;	
	int        get_number_of_bands() const { return 4; }
	bool       check_solution_type_available(KPSolutionType type) const;	
protected:
		
	void init_base_and_strain_matrix();	
	
private:

	static const int sparsity_pattern[28];	
	static const int diagonal_pattern[4];
			
};


} // end namespace

#endif /*KPMATRIX4X4ENDERSFOREMAN_H_*/
