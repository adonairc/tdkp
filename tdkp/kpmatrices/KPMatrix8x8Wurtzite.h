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

#ifndef KPMATRIX8X8WURTZITE_H_
#define KPMATRIX8X8WURTZITE_H_

#include "tdkp/kpmatrices/KPMatrixBase.h"

namespace tdkp {

class KPMatrix8x8Wurtzite : public tdkp::KPMatrixWurtziteBase
{
public:
	KPMatrix8x8Wurtzite();
	KPMatrix8x8Wurtzite(const Material* mat);	
	virtual ~KPMatrix8x8Wurtzite();

	const int* get_sparsity_pattern(int& length) const;	
	const int* get_diagonal_pattern(int& length) const;	
	int        get_number_of_bands() const { return 8; }	
		
	bool       check_solution_type_available(KPSolutionType type) const;
	double     get_optical_momentum_matrix_element_P1() const;
	double     get_optical_momentum_matrix_element_P2() const;
	
protected:	
	void init_base_and_strain_matrix();
		
	
private:	
	void determine_splitted_parameters(const WurtziteEffectiveMassParams& params, double &A5m, double &A6m) const;
	WurtziteEffectiveMassParams get_renormalized_parameters(const WurtziteEffectiveMassParams& params, const double& P1, const double& P2, const double& cb_edge, const double& vb_edge) const;

	
	inline int get_idx(int ii, int jj) const;
	static const int sparsity_pattern[72];	
	static const int diagonal_pattern[8];		
};

}

#endif /*KPMATRIX8X8WURTZITE_H_*/
