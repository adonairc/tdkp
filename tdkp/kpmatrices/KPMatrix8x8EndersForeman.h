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

#ifndef KPMATRIX8X8ENDERSFOREMAN_H_
#define KPMATRIX8X8ENDERSFOREMAN_H_

#include "tdkp/kpmatrices/KPMatrixBase.h"

namespace tdkp {
	
/** kp8x8 matrix using |S>, |X>, |Y> and |Z> as basis states
 * 
 * corresponds to foremans adaption of the enders 8x8 hamiltonian,
 * assuming an infinite bandgap (and thus decoupling the conduction band)
 */
class KPMatrix8x8EndersForeman : public tdkp::KPMatrixBase {
	
public:
	KPMatrix8x8EndersForeman(const Material* mat);
	KPMatrix8x8EndersForeman();
	virtual ~KPMatrix8x8EndersForeman();
	const int* get_sparsity_pattern(int& length) const;	
	const int* get_diagonal_pattern(int& length) const;	
	int        get_number_of_bands() const { return 8; }	
	bool       check_solution_type_available(KPSolutionType type) const;
		
protected:	
	void init_base_and_strain_matrix();	
	
private:
	inline int get_idx(int ii, int jj) const;
	static const int sparsity_pattern[88];	
	static const int diagonal_pattern[8];	
	
};

} // end of namespace

#endif /*KPMATRIX8X8ENDERSFOREMAN_H_*/
